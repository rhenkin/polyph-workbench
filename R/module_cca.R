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
				radioButtons(
					ns("dataset_source"),
					"Dataset source:",
					choices = c(
						"Use current study from matching (if available)" = "memory",
						"Load saved study from disk" = "disk"
					),
					selected = "memory"
				),
				conditionalPanel(
					condition = "input.dataset_source == 'disk'",
					ns = ns,
					selectizeInput(
						ns("dataset_list_selection"),
						"Select dataset:",
						choices = NULL,
						options = list(dropdownParent = 'body')
					)
				),
				uiOutput(ns("dataset_status")),
				actionButton(ns("load_dataset"), "Load dataset", class = "btn-primary")
			),
			accordion_panel(
				title = "Analysis",
				value = "analysis_panel",
				icon = bs_icon("bar-chart-line"),
				# BNF level selection at the top of analysis panel
				selectInput(ns("cca_bnf_level"),
										label = "Select BNF level for analysis:",
										choices = c("BNF_Chapter", "BNF_Section", "BNF_Paragraph", "BNF_Chemical_Substance"),
										selected = "BNF_Chemical_Substance"),
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
													 card(
													 	card_header("Polypharmacy burden"),
													 	card_body(
													 		vegawidgetOutput(ns("pp_pyramid_plot"))
													 	)
													 ),
													 card(
													 	card_header("MLTC burden"),
													 	card_body(
													 		vegawidgetOutput(ns("mltc_pyramid_plot"))
													 	)
													 )
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
									 card_header("Top 10 Long-term conditions"),
									 card_body(
									 	vegawidgetOutput(
									 		ns("topten_ltc_bar_plot")
									 	)
									 )), card(full_screen = TRUE,
									 				 card_header("Top 10 Background medications"),
									 				 card_body(
									 				 	vegawidgetOutput(
									 				 		ns("topten_substance_bar_plot")
									 				 	))))
					),
					nav_panel(title = "Advanced",
										card(
											card_header("Sociodemographics"), card_body(
												module_cca_sociodemographics_ui(ns("sociodemog"))
											)
										),
										card(
											card_header("Prevalence tables"), card_body(
												module_cca_prevalence_ui(ns("prevalence"))
											)
										),
										card(
											card_header("Co-prescription Analysis"), card_body(
												module_cca_copresc_ui(ns("copresc"))
											)
										)),
					full_screen  = TRUE
				)

			)
		)
	)
}

module_cca_server <- function(id, prepared_study_data_r = NULL) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Reactive values for data storage
		patient_data_r <- reactiveVal()
		prescriptions_r <- reactiveVal()
		prescriptions_aggregated_r <- reactiveVal()  # New: for BNF-aggregated prescriptions
		ltcs_r <- reactiveVal()
		cases_controls_r <- reactiveVal()

		# Update dataset list
		observe({
			updateSelectizeInput(session, "dataset_list_selection",
													 choices = list.files("studies", "*.rds"))
		})

		# Display dataset status
		output$dataset_status <- renderUI({
			if (input$dataset_source == "memory") {
				if (!is.null(prepared_study_data_r) && !is.null(prepared_study_data_r())) {
					study_data <- prepared_study_data_r()
					tags$div(
						class = "alert alert-success",
						style = "margin-top: 10px;",
						tags$strong("Available: "), study_data$metadata$study_name,
						tags$br(),
						sprintf("Cases: %d, Controls: %d",
										study_data$metadata$n_cases,
										study_data$metadata$n_controls_unique)
					)
				} else {
					tags$div(
						class = "alert alert-warning",
						style = "margin-top: 10px;",
						"No study currently available in memory. Create a matched cohort first or load from disk."
					)
				}
			} else {
				if (!is.null(input$dataset_list_selection) && input$dataset_list_selection != "") {
					tags$div(
						class = "alert alert-info",
						style = "margin-top: 10px;",
						tags$strong("Selected: "), input$dataset_list_selection
					)
				} else {
					tags$div(
						class = "alert alert-info",
						style = "margin-top: 10px;",
						"Please select a saved dataset."
					)
				}
			}
		})

		# Load dataset (from either source)
		observe({
			if (input$dataset_source == "memory") {
				# Load from memory
				req(prepared_study_data_r)
				dataset <- prepared_study_data_r()
				req(dataset)

				# Prepare matched_patids structure
				matched_patids <- dataset$matched_patids
				matched_patids[, group := ifelse(treatment == 1, "case", "control")]

				# Filter by index date
				filtered <- filter_by_index_date(
					dataset$all_prescriptions,
					dataset$all_ltc,
					matched_patids,
					lookback_days = 84
				)

				# Add burden groups
				patient_data <- add_burden_groups(
					dataset$all_patient_data,
					filtered$prescriptions,
					filtered$ltcs
				)

				# Store in reactive values
				patient_data_r(patient_data)
				prescriptions_r(filtered$prescriptions)
				ltcs_r(filtered$ltcs)
				cases_controls_r(matched_patids[, .(patid, index_date, substance, group)])

				# Update UI
				accordion_panel_close("cca_accordion", "load_dataset_panel")
				accordion_panel_open("cca_accordion", "analysis_panel")
				showNotification(
					paste("Dataset", dataset$metadata$study_name, "loaded from memory!"),
					type = "message",
					duration = 3
				)

			} else {
				# Load from disk
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
					paste("Dataset", input$dataset_list_selection, "loaded from disk!"),
					type = "message",
					duration = 3
				)
			}
		}) |> bindEvent(input$load_dataset)

		# NEW: Create BNF-aggregated prescriptions reactive
		observe({
			req(prescriptions_r(), input$cca_bnf_level)

			presc <- copy(prescriptions_r())

			# Aggregate substance to selected BNF level
			presc[, substance := bnf_lookup[match(substance, bnf_lookup$BNF_Chemical_Substance),
																			get(input$cca_bnf_level)]]

			# Remove any NA substances that might result from the lookup
			presc <- presc[!is.na(substance)]

			prescriptions_aggregated_r(presc)
		})

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
			create_burden_pyramid(patient_data_r(), "pp")
		})

		output$mltc_pyramid_plot <- renderVegawidget({
			req(patient_data_r())
			create_burden_pyramid(patient_data_r(), "n_ltc")
		})

		# Top substances plot - now uses aggregated prescriptions
		output$top_recentpresc_bar_plot <- renderVegawidget({
			req(cases_controls_r(), input$cca_bnf_level)

			# Aggregate cases_controls to selected BNF level
			cc <- copy(cases_controls_r())
			cc[, substance := bnf_lookup[match(substance, bnf_lookup$BNF_Chemical_Substance),
																	 get(input$cca_bnf_level)]]
			cc <- cc[!is.na(substance)]

			full_screen <- isTruthy(input$recent_presc_card_full_screen)
			create_top_substances_plot(cc, full_screen = full_screen)
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

		# Top background medications plot - now uses aggregated prescriptions
		output$topten_substance_bar_plot <- renderVegawidget({
			req(prescriptions_aggregated_r())
			freq_data <- calculate_frequency_stats(prescriptions_aggregated_r(), "substance")
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

		module_cca_sociodemographics_server(
			id = "sociodemog",
			patient_data_r = patient_data_r
		)

		# Prevalence tables sub-module - now passes aggregated prescriptions
		module_cca_prevalence_server(
			id = "prevalence",
			patient_data_r = patient_data_r,
			prescriptions_r = prescriptions_aggregated_r,  # Changed to use aggregated version
			ltcs_r = ltcs_r
		)

		module_cca_copresc_server(
			id = "copresc",
			prescriptions_r = prescriptions_aggregated_r,
			patient_data_r = patient_data_r,
			bnf_level = input$cca_bnf_level
		)
	})
}