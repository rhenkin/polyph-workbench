module_prescription_explorer_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Prescriptions"),
    	selectInput(ns("freq_table_bnf"),
    						label = "Select BNF level for analysis:",
    						choices = c("BNF_Chapter", "BNF_Section", "BNF_Paragraph", "BNF_Chemical_Substance")),
      accordion(
        open = FALSE,
        accordion_panel(
          title = "Prevalence table",
          value = "prev_tables",
          navset_tab(
            nav_panel("Cohort prevalence",
            						dataTableOutput(ns("subst_freq_table"))
            					),
            nav_panel("Prevalence across demographics",
                      selectizeInput(ns("select_demog_freq_var"),
                      							 label = "Select variable:",
                      							 choices = c("Sex" = "sex","Ethnic group" = "eth_group", "IMD quintile" = "imd_quintile", "PP" = "pp_group", "# LTCs" = "mltc_group")),
                      dataTableOutput(ns("demog_subs_freq_table")))
            )
        ),
        accordion_panel(
        	title ="Prescription prevalence per LTC",
        	value = "ltc_prev",
        		card(full_screen = TRUE,
        				 height = "60em",
        			fluidRow(column(6, uiOutput(ns("ltc_dropdown_ui"))),
        							 column(5, numericInput(ns("outcome_age_filter"), label = "Outcome age filter:", value = 100, min = 16, max = 100)),
        							 column(1, textOutput(ns("selected_pats")))),
        			dataTableOutput(ns("presc_by_ltc"))
        		)),
        accordion_panel(
        	title = "OR heatmap",
        	value = "or_heatmap_acc",
        	card(
        		full_screen = TRUE,
        		virtualSelectInput(ns("corr_heatmap_strat"), label = "Select group to filter:", choices = NULL,  autoSelectFirstOption = FALSE, search = FALSE, disableOptionGroupCheckbox = TRUE, multiple = FALSE),
        		vegawidgetOutput(ns("corr_heatmap"))
        	)
        ),
        accordion_panel(
        	"Polypharmacy transition per substance",
        	div("Column indicates current number of prescriptions before new prescription is added"),
        	dataTableOutput(ns("transition_table"))
        ),
        accordion_panel(
          "Time to outcome", module_presc_tto_ui(ns("tto_module"))
        ),
        accordion_panel(
        	"Acute prescriptions", module_acute_presc_ui(ns("acute_presc_module"))
        )
    )
  )
}

#' @export
module_prescription_explorer_server <- function(id, outcome_prescriptions, patient_data, ltc_data, selected_outcome, pp_groups_data, acute_outcome_prescriptions) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    #bnf_lookup_simple <- reactiveVal()
    bnf_lookup_simple <- reactive({
    	bnf_simple <- bnf_lookup[
    		,
    		.SD[1, .(BNF_Chemical_Substance, BNF_Subparagraph, BNF_Paragraph, BNF_Section, BNF_Chapter)], BNF_Chemical_Substance
    	] |> unique()
    	bnf_simple
	   })

    subst_pp_df <- reactive({
      df <- outcome_prescriptions()
      validate(need(nrow(df) > 0, "No valid patients found"))
      prepare_prescription_data(df, bnf_lookup_simple(), pp_groups_data(), patient_data(),
                                input$freq_table_bnf)
    })

    heatmap_var_lookup <- reactiveVal()

    observe({
    	req(subst_pp_df())
    	choices <- prepare_heatmap_stratification_choices(subst_pp_df())
    	updateVirtualSelect("corr_heatmap_strat", choices = choices)
    })

    subst_frequency <- reactive({
      calculate_item_frequency(subst_pp_df(), "substance")
    })

    output$subst_freq_table <- renderDataTable({
      df <- copy(subst_frequency())
    	df[,.(substance, N, pct_total = round(pct_total*100, digits = 2))]
    }, rownames = FALSE)

    demog_subs_freq_df <- reactive({
    	req(input$select_demog_freq_var)
    	calculate_demographic_prescription_frequency(
    	  subst_pp_df(),
    	  input$select_demog_freq_var,
    	  subst_frequency(),
    	  min_prevalence = 0.005
    	)
    })

    output$demog_subs_freq_table <- renderDataTable({
    	req(demog_subs_freq_df())
    	df <- demog_subs_freq_df()
    	df[!is.na(padj), pvalue := as.double(sprintf("%.2e", pvalue))]
    	df[!is.na(padj), padj := as.double(sprintf("%.2e", padj))]
    }, rownames = FALSE)

    output$ltc_dropdown_ui <- renderUI({
        ltcs <- ltc_data()
        virtualSelectInput(ns("ltc_dropdown"),
            label = "Select 1 or more LTCs:",
            choices = with(ltc_chapters[ltc %in% unique(ltcs$term)],
            split(ltc, body_system)), multiple = TRUE, search = TRUE)
    })

    age_filtered_prescriptions <- reactive({
    	filter_by_age(subst_pp_df(), input$outcome_age_filter)
    })

    output$selected_pats <- renderText({
    	format_patient_count(age_filtered_prescriptions())
    })

  	output$presc_by_ltc <- renderDataTable({
			req(input$ltc_dropdown)
  		calculate_prescription_prevalence_by_ltc(
  		  ltc_data(),
  		  age_filtered_prescriptions(),
  		  input$ltc_dropdown,
  		  min_prevalence = 0.005
  		)
  	}, rownames = FALSE)


  	output$corr_heatmap <- renderVegawidget({
  		heatmap_data <- prepare_prescription_or_heatmap_data(
  		  subst_pp_df(),
  		  input$corr_heatmap_strat,
  		  min_prevalence = 0.005
  		)
  		or_heatmap(heatmap_data$pairs_dt, heatmap_data$or_matrix) |> as_vegaspec()
  	})

    output$transition_table <- renderDataTable({
    	med_data <- copy(subst_pp_df())
    	setorder(med_data, patid, start_date)
    	level_profiles <- analyze_medications(med_data, mode = "level")
    }, rownames = FALSE)

    module_presc_tto_server(id = "tto_module", subst_pp_df = subst_pp_df, subst_freq = subst_frequency)
    module_acute_presc_server(id = "acute_presc_module", acute_prescriptions_df = acute_outcome_prescriptions)

  })
}