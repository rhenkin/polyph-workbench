module_cca_ui <- function(id) {
	ns <- NS(id)
	nav_panel(
		"Case-control Analysis",
		value = "cca",
		accordion(
			open = NULL,
			id = ns("cca_accordion"),
			accordion_panel(
				title = "Load dataset",
				value = "load_dataset_panel",
				icon = bs_icon("database"),
				selectizeInput(
					ns("dataset_list_selection"),
					"Select dataset:",
					choices = NULL,
					options = list(dropdownParent = 'body')
				),
				actionButton(ns("load_saved_dataset"), "Load saved dataset")
			),
			accordion_panel(
				title = "Analysis",
				value = "analysis_panel",
				icon = bs_icon("bar-chart-line"),
				navset_card_tab(
					nav_panel(
						title = "Overview",
						layout_columns(col_widths = c(2,4,4,2),
							verticalLayout(
								value_box(
									showcase_layout = showcase_left_center(width = "200px", max_height = "150px"),
									height = "150px",
									title = "Cases",
									value = textOutput(ns("value_box_cases")),
									theme = "red",
									showcase = bs_icon("people-fill")
								),
								value_box(
									showcase_layout = showcase_left_center(width = "200px", max_height = "150px"),
									height = "150px",
									title = "Controls",
									value = textOutput(ns("value_box_controls")),
									theme = "blue",
									showcase = bs_icon("people-fill")
								)
							),
							card(vegawidgetOutput(ns(
								"pp_pyramid_plot"
							))),
							card(vegawidgetOutput(ns(
								"mltc_pyramid_plot"
							)))
						),
						layout_columns(
							card(full_screen = TRUE,
									 id = ns("recent_presc_card"),
								card_header("Top 10 Latest prescriptions"),
								card_body(
									vegawidgetOutput(ns("top_recentpresc_bar_plot"))
									)
								),
							card(full_screen = TRUE,
									 vegawidgetOutput(
							ns("topten_ltc_bar_plot")
						)), card(full_screen = TRUE,
										 vegawidgetOutput(
							ns("topten_substance_bar_plot")
						)))
					),
                    
                    
					nav_panel(title = "Advanced", card(
						card_header("Prevalence tables"), card_body(
							accordion(
								open = FALSE,
								accordion_panel(title = "Long-term conditions", value = "ltc_prev_tables", navset_tab(
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
										dataTableOutput(ns("ltc_freq_table"))
									)
								)),
								accordion_panel(title = "Prescriptions", value = "presc_prev_tables", navset_tab(
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
										dataTableOutput(ns("presc_freq_table"))
									)
								)),
								accordion_panel(
									title = "Prescription prevalence per LTC",
									value = "presc_ltc_prev",
									card(
										full_screen = TRUE,
										height = "60em",
										uiOutput(ns("ltc_dropdown_ui")),
										# fluidRow(column(6, uiOutput(ns("ltc_dropdown_ui"))),
										# 				 column(5, numericInput(ns("outcome_age_filter"), label = "Outcome age filter:", value = 100, min = 16, max = 100)),
										# 				 column(1, textOutput(ns("selected_pats")))),
										reactableOutput(ns("presc_by_ltc"), height = "50em")
									)
								),
							)
						)
					)),
					full_screen  = TRUE
				)

			)
		)
	)
}

module_cca_server <- function(id, matched_data) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Reactive values for data storage
    patient_data_r <- reactiveVal()
    prescriptions_r <- reactiveVal()
    ltcs_r <- reactiveVal()
    cases_controls_r <- reactiveVal()
    
    # Update dataset list
    observe({
      updateSelectizeInput(session, "dataset_list_selection", 
                          choices = list.files("studies", "*.rds"))
    })
    
    # Load dataset
    observe({
      req(input$dataset_list_selection)
      
      # Load and prepare data
      dataset <- load_cca_dataset(file.path("studies", input$dataset_list_selection))
      
      # Filter by index date
      filtered <- filter_by_index_date(
        dataset$prescriptions,
        dataset$ltcs,
        dataset$matched_patids,
        lookback_days = 84
      )
      
      # Add burden groups
      patient_data <- add_burden_groups(
        dataset$patient_data,
        filtered$prescriptions,
        filtered$ltcs
      )
      
      # Store in reactive values
      patient_data_r(patient_data)
      prescriptions_r(filtered$prescriptions)
      ltcs_r(filtered$ltcs)
      cases_controls_r(dataset$matched_patids[, .(patid, index_date, substance, group)])
      
      # Update UI
      accordion_panel_close("cca_accordion", "load_dataset_panel")
      accordion_panel_open("cca_accordion", "analysis_panel")
      showNotification(
        paste("Dataset", input$dataset_list_selection, "loaded successfully!"),
        type = "message",
        duration = 3
      )
    }) |> bindEvent(input$load_saved_dataset)
    
    # Value boxes
    output$value_box_cases <- renderText({
      req(patient_data_r())
      prettyNum(nrow(patient_data_r()[treatment == 1]), big.mark = ",")
    })
    
    output$value_box_controls <- renderText({
      req(patient_data_r())
      prettyNum(nrow(patient_data_r()[treatment == 0]), big.mark = ",")
    })
    
    # Pyramid plots
    output$pp_pyramid_plot <- renderVegawidget({
      req(patient_data_r())
      create_burden_pyramid(patient_data_r(), "pp", title = "Polypharmacy burden")
    })
    
    output$mltc_pyramid_plot <- renderVegawidget({
      req(patient_data_r())
      create_burden_pyramid(patient_data_r(), "n_ltc", title = "MLTC burden")
    })
    
    # Top substances plot
    output$top_recentpresc_bar_plot <- renderVegawidget({
      req(cases_controls_r())
      full_screen <- isTruthy(input$recent_presc_card_full_screen)
      create_top_substances_plot(cases_controls_r(), full_screen = full_screen)
    })
    
    # Top conditions plot
    output$topten_ltc_bar_plot <- renderVegawidget({
      req(ltcs_r())
      freq_data <- calculate_frequency_stats(ltcs_r(), "term")
      create_top_conditions_plot(
        freq_data,
        title = "Top 10 conditions most prevalent in cases"
      )
    })
    
    # Top background medications plot
    output$topten_substance_bar_plot <- renderVegawidget({
      req(prescriptions_r())
      freq_data <- calculate_frequency_stats(prescriptions_r(), "substance")
      ratios <- calculate_case_control_ratios(freq_data, "substance", min_case_pct = 5)
      top_items <- ratios[order(-case_pct)][1:10, substance]
      
      term_filtered <- freq_data[substance %in% top_items]
      term_filtered[, `:=`(
        diff = pct[group == "case"] - pct[group == "control"],
        diff_label = sprintf("+%.1f%%", 
          abs(pct[group == "case"] - pct[group == "control"])),
        max_pct = max(pct)
      ), by = substance]
      
      term_filtered[nchar(substance) > 15, 
        substance := paste0(strtrim(substance, 15), "...")
      ]
      
      grouped_bar_plot(
        term_filtered,
        "substance",
        title = "Top 10 Background medications most prevalent in cases"
      ) |> as_vegaspec()
    })
    
    # Prevalence tables sub-module
    module_cca_prevalence_server(
      id = "prevalence",
      patient_data_r = patient_data_r,
      prescriptions_r = prescriptions_r,
      ltcs_r = ltcs_r,
      ltc_chapters = ltc_chapters
    )
  })
}