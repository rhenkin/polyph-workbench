module_ltc_explorer_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Long-term conditions"),
    accordion(
      open = FALSE,
      accordion_panel(
        title = "Prevalence table",
        navset_tab(
          nav_panel("Cohort prevalence", dataTableOutput(ns("ltc_freq_table"))),
          nav_panel("Prevalence across demographics",
                    selectizeInput(ns("select_ltcdemog_freq_var"), label = "Select variable:", choices = c("Sex" = "sex","Ethnic group" = "eth_group", "IMD quintile" = "imd_quintile", "PP" = "pp_group", "# LTCs" = "mltc_group")),
                    dataTableOutput(ns("demog_ltc_freq_table")))
        )
      ),
      accordion_panel(
      	"LTC prevalence per drug",
      	card(full_screen = TRUE,
      			 height = "60em",
      		fluidRow(column(6,
      										# uiOutput(ns("drug_dropdown_ui"))

      										selectizeInput(ns("drug_dropdown"), label = "Select 1 or more substance:", choices = NULL, multiple = TRUE,
      																			options = list(
      																				splitOn = ';'
      																			)
      																			),
      										textInput(ns("paste_substances"),
      															label = "Or paste semicolon-separated list:",
      															placeholder = "Paracetamol; Aspirin")
      										),
      						 column(5, numericInput(ns("outcome_age_filter"), label = "Max age at outcome:", value = 100, min = 16, max = 100)),
      						 column(1, textOutput(ns("selected_pats")))),
      		dataTableOutput(ns("ltc_by_presc"))
      	))
      # accordion_panel(
      #   title = "Network",
      #   value = "ltcnet",
      #   visNetworkOutput(ns("ltc_network"), height = "600px")
      # )
    )
  )
}

module_ltc_explorer_server <- function(id, outcome_prescriptions, ltc_data, patient_data, pp_groups_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    ltc_chapters <- fread("chapters.tsv")
    valid_ltcs <- reactive({
      outcome_df <- outcome_prescriptions()
      ltcs <- ltc_data()
      ltcs[patid %in% outcome_df$patid]
    })

    ltc_freq_df <- reactive({
      ltcs <- valid_ltcs()

      term_freq <- ltcs[,.N,.(term)]
      total_patids <- length(unique(ltcs$patid))
      term_freq[, pct_total := signif(N/total_patids, digits = 2)]
      setorder(term_freq,-N)
    })

    output$ltc_freq_table <- renderDataTable({
      ltc_freq_df()
    })

    output$demog_ltc_freq_table <- renderDataTable({
      req(input$select_ltcdemog_freq_var)
      df <- valid_ltcs()
      outcome_df <- outcome_prescriptions()
      demog_var <- input$select_ltcdemog_freq_var

      df <- merge(patient_data(), df, by = "patid")
      df <- merge(df, pp_groups_data(), by = "patid")
      df <- merge(df, unique(outcome_df[,.(patid, mltc_group)]), by = "patid")

      df <- df[term %in% ltc_freq_df()[pct_total>0.005, term]]

      # Calculate statistics
      cat_totals <- calculate_category_totals(df, demog_var)
      df_stats <- calculate_demographic_stats(df, demog_var, ltc_freq_df(), "term")

      # Perform chi-square tests
      chisq_tests <- perform_chisq_tests(df_stats, cat_totals, demog_var, "term")

      # Format and return results
      return(format_results(df_stats, chisq_tests, cat_totals, demog_var, ltc_freq_df(), "term"))

    }, rownames = FALSE)


    observe({
    	prescs <- outcome_prescriptions()
    	updateSelectizeInput(session = session,
    											 inputId = "drug_dropdown",
    											 choices = c("",sort(unique(prescs$substance))),
    											 options = list(delimiter = "; "),
    											 server = TRUE)
    })

    observeEvent(input$paste_substances, {
    	req(outcome_prescriptions())
    	if(input$paste_substances != "") {
    		prescs <- outcome_prescriptions()
    		pasted_items <- trimws(strsplit(input$paste_substances, ";")[[1]])
    		# Filter to only items that exist in your choices
    		valid_items <- pasted_items[pasted_items %in% sort(unique(prescs$substance))]

    		updateSelectizeInput(session, "drug_dropdown", selected = valid_items)
    		updateTextInput(session, "paste_substances", value = "") # Clear the text input
    	}
    })

    age_filtered_prescriptions <- reactive({
    	outcome_prescriptions()[outcome_age <= input$outcome_age_filter*365.25]
    })

    output$selected_pats <- renderText({
    	dt <- age_filtered_prescriptions()
    	paste0("Patients within age range: ", prettyNum(uniqueN(dt$patid), big.mark = ","))
    })

    output$ltc_by_presc <- renderDataTable({
    	req(input$drug_dropdown)
    	input$drug_dropdown
    	ltcs <- valid_ltcs()
    	prescriptions <- age_filtered_prescriptions()
    	ltcs <- ltcs[patid %in% prescriptions$patid]

    	# Patients taking the selected drug
    	patids <- unique(prescriptions[substance %in% input$drug_dropdown, patid])

    	# Patients WITH the drug
    	ltc_freq <- ltcs[patid %in% patids,
    									 list(N_with = uniqueN(patid),
    									 		 Prevalence = round(100*(uniqueN(patid)/length(patids)), digits=2)),
    									 term]

    	# Patients WITHOUT the drug
    	unselected_patids <- ltcs[!patid %in% patids, uniqueN(patid)]
    	not_selected_freq <- ltcs[!patid %in% patids,
    														list(N_without = uniqueN(patid),
    																 Prevalence_Unselected = round(100*(uniqueN(patid)/unselected_patids), digits=2)),
    														term]

    	result <- merge(ltc_freq, not_selected_freq)

    	result <- result[Prevalence >= 0.01]

    	result[, `:=`(
    		total_with_disease = length(patids),
    		total_without_disease = unselected_patids,
    		Prevalence_Ratio = round(Prevalence / Prevalence_Unselected, digits = 2)
    	)]

    	result[is.infinite(Prevalence_Ratio), Prevalence_Ratio := 0]

    	# Calculate 95% CI for the ratio using log method
    	result[, `:=`(
    		p1 = N_with / total_with_disease,
    		p2 = N_without / total_without_disease
    	)]

    	result[, `:=`(
    		log_ratio = log(Prevalence_Ratio),
    		se_log_ratio = sqrt((1/N_with) - (1/total_with_disease) +
    													(1/N_without) - (1/total_without_disease))
    	)]

    	result[, `:=`(
    		CI_lower = round(exp(log_ratio - 1.96 * se_log_ratio), digits = 2),
    		CI_upper = round(exp(log_ratio + 1.96 * se_log_ratio), digits = 2)
    	)]

    	result[,
    				 CI_95 := paste0("(", CI_lower, " - ", CI_upper, ")")
    	]

    	# Clean up intermediate columns if desired
    	result[(CI_lower > 1.0 | CI_upper < 1.0), term := paste0(term, "*")]
    	result[, c("p1", "p2", "log_ratio", "se_log_ratio", "CI_lower", "CI_upper") := NULL]

    	result[order(-Prevalence_Ratio), .(term, Prevalence, Prevalence_Unselected, Prevalence_Ratio, CI_95)]
    }, rownames = FALSE)

    # output$ltc_network <- renderVisNetwork({
    #   outcome_df <- outcome_prescriptions()
    #   ltcs <- ltc_data()
    #   ltcs <- ltcs[patid %in% outcome_df$patid]
    #   # browser()
    #   # df <- merge(ltcs, ltc_chapters, by.x = "term", by.y = "ltc")
    #   # df[,term := body_system]
    #   connection_data <- calculate_network(ltcs, "term", min_prevalence = 0.05, min_co_prevalence = 0.05)
    #   network <- create_network_visualization(connection_data, or_range = c(1,50))
    #   network |>
    #   	visNodes(font = list(
    #   		size = 16,
    #   		face = "arial",
    #   		background = "rgba(255, 255, 255, 0.8)",
    #   		strokeWidth = 0,
    #   		align = "center"
    #   	),
    #   	color = list(background = "#e1e1e1", border = "#e1e1e1"),
    #   	shadow = list(enabled = TRUE, size = 2)) |>
    #   	visEdges(smooth = FALSE,
    #   					 scaling = list(
    #   					 	min = 2,
    #   					 	max = 8
    #   					 ),
    #   					 hoverWidth = 2,
    #   					 color = list(inherit = "both")) |>
    #   	visOptions(highlightNearest = list(enabled = TRUE, algorithm = "all", degree = list(from = 0, to = 1))) |>
    #   	visIgraphLayout()
    # })
  })
}