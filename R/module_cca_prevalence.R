module_cca_prevalence_ui <- function(id) {
  ns <- NS(id)

  accordion(
    open = FALSE,
    accordion_panel(
      title = "Long-term conditions",
      value = "ltc_prev_tables",
      navset_tab(
        nav_panel(
          "Group prevalences",
          virtualSelectInput(
            ns("ltc_freq_strat_variable"),
            label = "Select a subset (optional):",
            choices = NULL,
            autoSelectFirstOption = FALSE,
            search = FALSE,
            showValueAsTags = TRUE,
            disableOptionGroupCheckbox = TRUE,
            multiple = FALSE,
            dropboxWrapper = "body"
          ),
          reactableOutput(ns("ltc_freq_table"), height = "50em")
        ),
        nav_panel(
        	title = "Prevalence by prescription filter",
        	card(
        		full_screen = TRUE,
        		height = "60em",
        		div("Select medications to show the prevalence difference only among the patients taking those medications."),
        		div("Table contains conditions with a minimum of 1% prevalence in both cases and controls"),
        		uiOutput(ns("presc_dropdown_ui")),
        		reactableOutput(ns("ltc_by_presc"), height = "50em")
        	)
        )
      )
    ),
    accordion_panel(
      title = "Prescriptions",
      value = "presc_prev_tables",
      navset_tab(
        nav_panel(
          "Group prevalences",
          virtualSelectInput(
            ns("presc_freq_strat_variable"),
            label = "Select a subset (optional):",
            choices = NULL,
            autoSelectFirstOption = FALSE,
            search = FALSE,
            showValueAsTags = TRUE,
            disableOptionGroupCheckbox = TRUE,
            multiple = FALSE,
            dropboxWrapper = "body"
          ),
          reactableOutput(ns("presc_freq_table"), height = "50em")
        ),
        nav_panel(
        	title = "Prevalence by condition filter",
        	card(
        		full_screen = TRUE,
        		height = "60em",
        		div("Select conditions to show the prevalence difference only among the patients that were diagnosed with those conditions."),
        		div("Table contains medications with a minimum of 1% prevalence in both cases and controls"),
        		uiOutput(ns("ltc_dropdown_ui")),
        		reactableOutput(ns("presc_by_ltc"), height = "50em")
        	)
        )
      )
    )
  )
}

module_cca_prevalence_server <- function(id, patient_data_r, prescriptions_r,
                                        ltcs_r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Update stratification choices
    observe({
      req(patient_data_r())
      choices <- create_stratification_choices(patient_data_r())
      updateVirtualSelect("ltc_freq_strat_variable", choices = choices)
      updateVirtualSelect("presc_freq_strat_variable", choices = choices)
    })

    # LTC frequency with optional stratification
    ltc_freq_df <- reactive({
      req(ltcs_r())
      ltcs <- ltcs_r()
      patient_data <- patient_data_r()
      strat_var <- input$ltc_freq_strat_variable

      if (strat_var != "") {
        parts <- strsplit(strat_var, "#")[[1]]
        column_name <- parts[1]
        filter_value <- parts[2]
        selected_patids <- patient_data[get(column_name) == filter_value, patid]
        ltcs <- ltcs[patid %in% selected_patids]
      }

      calculate_frequency_stats(ltcs, "term")
    })


    # Render tables
    output$ltc_freq_table <- renderReactable({
      table_data <- create_prevalence_ratio_table(ltc_freq_df(), "term")
      colnames(table_data) <- c("Condition", "Cases (%)", "Controls (%)", "Case/Control Ratio")
      reactable(table_data[order(-`Case/Control Ratio`)],
      					showPageInfo = FALSE, defaultPageSize = 15)
    })

    # Presc dropdown
    output$presc_dropdown_ui <- renderUI({
    	req(prescriptions_r())
    	prescs <- prescriptions_r()
    	sub_chapter <- bnf_lookup[,.(BNF_Chemical_Substance, BNF_Section)] |> unique()
    	sub_chapter <- sub_chapter[,first(.SD) ,BNF_Chemical_Substance]
    	choices <- prescs[sub_chapter,
    										.(substance, BNF_Section),
    										on = .(substance = BNF_Chemical_Substance)] |> unique()
    	choices <- choices[order(BNF_Section, substance)]
    	with(choices, split(substance, BNF_Section))
    	virtualSelectInput(
    		ns("presc_dropdown"),
    		label = "Select 1 or more substances:",
    		choices = with(choices, split(substance, BNF_Section)),
    		multiple = TRUE,
    		search = TRUE
    	)
    })

    # LTC by prescription table
    output$ltc_by_presc <- renderReactable({
    	req(ltcs_r(), input$presc_dropdown)

    	result <- calculate_prevalence_cca(
    		ltcs_r(),
    		prescriptions_r(),
    		input$presc_dropdown,
    		"substance",
    		"term"
    	)

    	validate(need(!is.null(result), "No results found for filter"))

    	reactable(
    		result,
    		columns = list(
    			p_value = colDef(show = FALSE),  # Hide these columns
    			p_adj = colDef(show = FALSE)
    		),
    		details = function(index) {
    			p_val <- result[index, p_value]
    			p_adj_val <- result[index, p_adj]
    			if (is.na(p_val) & is.na(p_adj_val)) return(NULL)
    			htmltools::div(
    				style = "padding: 16px",
    				htmltools::tags$b("Statistical Testing:"),
    				htmltools::tags$div(
    					style = "margin-top: 8px",
    					sprintf("Raw p-value: %.4f", p_val)
    				),
    				htmltools::tags$div(
    					sprintf("Adjusted p-value: %.4f", p_adj_val)
    				),
    				htmltools::tags$div(
    					style = "margin-top: 8px; font-style: italic; color: #666;",
    					if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
    				)
    			)
    		}
    		,showPageInfo = FALSE, defaultPageSize = 15)
    })

    # Prescription frequency with optional stratification
    presc_freq_df <- reactive({
    	req(prescriptions_r())
    	presc <- prescriptions_r()
    	patient_data <- patient_data_r()
    	strat_var <- input$presc_freq_strat_variable

    	if (strat_var != "") {
    		parts <- strsplit(strat_var, "#")[[1]]
    		column_name <- parts[1]
    		filter_value <- parts[2]
    		selected_patids <- patient_data[get(column_name) == filter_value, patid]
    		presc <- presc[patid %in% selected_patids]
    	}

    	calculate_frequency_stats(presc, "substance")
    })

    output$presc_freq_table <- renderReactable({
      table_data <- create_prevalence_ratio_table(presc_freq_df(), "substance")
      colnames(table_data) <- c("Substance", "Cases (%)", "Controls (%)", "Case/Control Ratio")
      reactable(table_data[order(-`Case/Control Ratio`)],
      					showPageInfo = FALSE, defaultPageSize = 15)
    })

    # LTC dropdown
    output$ltc_dropdown_ui <- renderUI({
      req(ltcs_r())
      ltcs <- ltcs_r()

      virtualSelectInput(
        ns("ltc_dropdown"),
        label = "Select 1 or more LTCs:",
        choices = with(
          ltc_chapters[ltc %in% unique(ltcs$term)],
          split(ltc, body_system)
        ),
        multiple = TRUE,
        search = TRUE
      )
    })

    # Prescription by LTC table
    output$presc_by_ltc <- renderReactable({
      req(ltcs_r(), input$ltc_dropdown)

      result <- calculate_prevalence_cca(
        prescriptions_r(),
        ltcs_r(),
        input$ltc_dropdown,
        "term",
        "substance"
      )

      reactable(
      	result,
      	columns = list(
      		p_value = colDef(show = FALSE),  # Hide these columns
      		p_adj = colDef(show = FALSE)
      	),
      	details = function(index) {
      		p_val <- result[index, p_value]
      		p_adj_val <- result[index, p_adj]

      		htmltools::div(
      			style = "padding: 16px",
      			htmltools::tags$b("Statistical Testing:"),
      			htmltools::tags$div(
      				style = "margin-top: 8px",
      				sprintf("Raw p-value: %.4f", p_val)
      			),
      			htmltools::tags$div(
      				sprintf("Adjusted p-value: %.4f", p_adj_val)
      			),
      			htmltools::tags$div(
      				style = "margin-top: 8px; font-style: italic; color: #666;",
      				if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
      			)
      		)
      	},showPageInfo = FALSE, defaultPageSize = 15)
    })


  })
}