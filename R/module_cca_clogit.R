# R/module_cca_clogit.R

module_cca_clogit_ui <- function(id) {
	ns <- NS(id)

	accordion(
		open = FALSE,
		accordion_panel(
			title = "Conditional Logistic Regression",
			value = "clogit_panel",

			card_body(
				p("Model the association between medications and case status, adjusting for selected long-term conditions and matching strata."),

				h4("Step 1: Select covariates (LTCs)"),

				layout_columns(
					col_widths = c(4, 4, 4),
					numericInput(
						ns("ltc_min_case_prev"),
						"Minimum LTC case prevalence (%):",
						value = 5,
						min = 1,
						max = 50,
						step = 1
					),
					numericInput(
						ns("ltc_min_or"),
						"Minimum LTC OR threshold:",
						value = 1.5,
						min = 1,
						max = 5,
						step = 0.1
					),
					numericInput(
						ns("ltc_max_count"),
						"Maximum number of LTCs:",
						value = 20,
						min = 5,
						max = 50,
						step = 5
					)
				),

				actionButton(ns("filter_ltcs"), "Filter LTCs", class = "btn-secondary"),

				conditionalPanel(
					condition = "output.ltcs_filtered",
					ns = ns,
					hr(),
					p(strong("Selected LTCs for adjustment:")),
					textOutput(ns("ltc_count_text")),
					reactableOutput(ns("filtered_ltcs_table"))
				),

				hr(),

				h4("Step 2: Select medications to model"),

				numericInput(
					ns("med_min_case_prev"),
					"Minimum medication case prevalence (%):",
					value = 2,
					min = 0.5,
					max = 20,
					step = 0.5
				),

				virtualSelectInput(
					ns("selected_meds"),
					"Select medications to model (leave empty for all):",
					choices = NULL,
					multiple = TRUE,
					search = TRUE,
					dropboxWrapper = "body"
				),

				checkboxInput(
					ns("include_interactions"),
					"Model medication interactions (pairwise)",
					value = FALSE
				),

				conditionalPanel(
					condition = "input.include_interactions",
					ns = ns,
					div(
						class = "alert alert-info",
						style = "margin-top: 10px;",
						tags$strong("Note: "),
						"Interaction models will fit separate models for each medication pair. ",
						"Selecting specific medications is strongly recommended to avoid long computation times."
					)
				),

				actionButton(ns("run_clogit"), "Run Models", class = "btn-primary"),

				hr(),

				conditionalPanel(
					condition = "output.results_available",
					ns = ns,
					h4("Results"),
					# Add model description
					div(
						class = "alert alert-light",
						style = "margin-bottom: 15px;",
						tags$strong("Model:"), " Logistic regression adjusting for matching variables and selected LTCs",
						tags$br(),
						uiOutput(ns("model_formula_text"))
					),
					downloadButton(ns("download_results"), "Download Results", class = "btn-sm btn-secondary"),
					br(), br(),
					reactableOutput(ns("clogit_results_table"))
				)
			)
		)
	)
}

module_cca_clogit_server <- function(id, patient_data_r, prescriptions_r, ltcs_r) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Calculate frequency data (reusing existing functions)
		ltc_freq_data <- reactive({
			req(ltcs_r(), patient_data_r())
			freq <- calculate_frequency_stats(ltcs_r(), "term")
			create_prevalence_ratio_table(freq, "term", data_with_group = ltcs_r())
		})

		med_freq_data <- reactive({
			req(prescriptions_r(), patient_data_r())
			freq <- calculate_frequency_stats(prescriptions_r(), "substance")
			create_prevalence_ratio_table(freq, "substance", data_with_group = prescriptions_r())
		})

		# Reactive values
		filtered_ltcs_r <- reactiveVal(NULL)
		clogit_results_r <- reactiveVal(NULL)

		# Filter LTCs based on thresholds
		observeEvent(input$filter_ltcs, {
			req(ltc_freq_data())

			# Get the LTC frequency data (already has ORs calculated)
			ltc_data <- ltc_freq_data()

			# Filter by case prevalence and OR
			filtered <- ltc_data[
				case >= input$ltc_min_case_prev &
					OR >= input$ltc_min_or &
					!is.na(OR)
			]

			# Order by OR descending and take top N
			setorder(filtered, -OR)
			filtered <- filtered[1:min(input$ltc_max_count, nrow(filtered))]

			filtered_ltcs_r(filtered)

			output$ltcs_filtered <- reactive({ TRUE })
			outputOptions(session$output, "ltcs_filtered", suspendWhenHidden = FALSE)

			showNotification(
				sprintf("%d LTCs selected for adjustment", nrow(filtered)),
				type = "message",
				duration = 3
			)
		})

		# Display filtered LTCs count
		output$ltc_count_text <- renderText({
			req(filtered_ltcs_r())
			sprintf("%d LTCs will be included as covariates", nrow(filtered_ltcs_r()))
		})

		# Display filtered LTCs table
		output$filtered_ltcs_table <- renderReactable({
			req(filtered_ltcs_r())

			reactable(
				filtered_ltcs_r()[, .(term, case, control, OR, OR_CI_lower, OR_CI_upper)],
				columns = list(
					term = colDef(name = "LTC", minWidth = 200),
					case = colDef(name = "Case %", format = colFormat(digits = 1)),
					control = colDef(name = "Control %", format = colFormat(digits = 1)),
					OR = colDef(name = "OR", format = colFormat(digits = 2)),
					OR_CI_lower = colDef(name = "CI Lower", format = colFormat(digits = 2)),
					OR_CI_upper = colDef(name = "CI Upper", format = colFormat(digits = 2))
				),
				defaultPageSize = 10,
				compact = TRUE
			)
		})

		# Update medication choices based on prevalence threshold
		observe({
			req(med_freq_data())

			med_data <- med_freq_data()

			# Filter by minimum case prevalence
			eligible_meds <- med_data[case >= input$med_min_case_prev, substance]

			updateVirtualSelect(
				"selected_meds",
				choices = sort(eligible_meds)
			)
		})

		# Dynamic model formula text
		output$model_formula_text <- renderUI({
			if (input$include_interactions) {
				tags$code("treatment ~ med1 + med2 + med1:med2 + LTCs + stratum (reporting interaction term only)")
			} else {
				tags$code("treatment ~ medication + LTCs + stratum")
			}
		})

		# Run conditional logistic regression
		observeEvent(input$run_clogit, {
			req(filtered_ltcs_r(), patient_data_r(), prescriptions_r(), ltcs_r())

			# Determine which medications to model
			if (is.null(input$selected_meds) || length(input$selected_meds) == 0) {
				# Use all eligible medications
				req(med_freq_data())
				med_data <- med_freq_data()
				medications_to_model <- med_data[case >= input$med_min_case_prev, substance]

				showNotification(
					sprintf("No medications selected - running models for all %d eligible medications",
									length(medications_to_model)),
					type = "message",
					duration = 5
				)
			} else {
				medications_to_model <- input$selected_meds
			}

			if (length(medications_to_model) == 0) {
				showNotification(
					"No medications meet the prevalence threshold",
					type = "warning",
					duration = 5
				)
				return()
			}

			# Check if running interaction models
			if (input$include_interactions) {
				# Generate all pairs
				med_pairs <- combn(medications_to_model, 2, simplify = FALSE)

				showNotification(
					sprintf("Running interaction models for %d medication pairs...", length(med_pairs)),
					type = "message",
					duration = NULL,
					id = "clogit_notification"
				)

				tryCatch({
					results <- run_interaction_models(
						med_pairs = med_pairs,
						selected_ltcs = filtered_ltcs_r()$term,
						patient_data = patient_data_r(),
						prescriptions = prescriptions_r(),
						ltcs = ltcs_r()
					)

					clogit_results_r(results)

					output$results_available <- reactive({ TRUE })
					outputOptions(session$output, "results_available", suspendWhenHidden = FALSE)

					removeNotification("clogit_notification")
					showNotification(
						sprintf("Interaction models fitted for %d medication pairs", length(med_pairs)),
						type = "message",
						duration = 3
					)

				}, error = function(e) {
					removeNotification("clogit_notification")
					showNotification(
						paste("Error running interaction models:", e$message),
						type = "error",
						duration = 10
					)
				})
			} else {
				# Original single medication models
				showNotification("Running logistic regression models...",
												 type = "message", duration = NULL, id = "clogit_notification")

				tryCatch({
					results <- run_conditional_logistic_models(
						medications = medications_to_model,
						selected_ltcs = filtered_ltcs_r()$term,
						patient_data = patient_data_r(),
						prescriptions = prescriptions_r(),
						ltcs = ltcs_r()
					)

					clogit_results_r(results)

					output$results_available <- reactive({ TRUE })
					outputOptions(session$output, "results_available", suspendWhenHidden = FALSE)

					removeNotification("clogit_notification")
					showNotification(
						sprintf("Models fitted for %d medications", length(medications_to_model)),
						type = "message",
						duration = 3
					)

				}, error = function(e) {
					removeNotification("clogit_notification")
					showNotification(
						paste("Error running models:", e$message),
						type = "error",
						duration = 10
					)
				})
			}
		})

		# Display results table
		output$clogit_results_table <- renderReactable({
			req(clogit_results_r())

			results <- copy(clogit_results_r())

			# Calculate total cases and controls for percentages
			total_cases <- patient_data_r()[treatment == 1, .N]
			total_controls <- patient_data_r()[treatment == 0, .N]

			# Check if interaction results or single medication results
			if ("med1" %in% names(results)) {
				# Interaction results
				reactable(
					results[, .(med1, med2, interaction_OR, interaction_CI_lower, interaction_CI_upper,
											interaction_p, pct_cases_both, pct_controls_both,
											n_cases_both, n_controls_both, n_ltc_covariates, convergence)],
					columns = list(
						med1 = colDef(name = "Medication 1", minWidth = 150),
						med2 = colDef(name = "Medication 2", minWidth = 150),
						interaction_OR = colDef(name = "Interaction OR", format = colFormat(digits = 3)),
						interaction_CI_lower = colDef(name = "CI Lower", format = colFormat(digits = 3)),
						interaction_CI_upper = colDef(name = "CI Upper", format = colFormat(digits = 3)),
						interaction_p = colDef(name = "P-value", format = colFormat(digits = 4)),
						pct_cases_both = colDef(name = "Cases % (both)", format = colFormat(digits = 2)),
						pct_controls_both = colDef(name = "Controls % (both)", format = colFormat(digits = 2)),
						n_cases_both = colDef(show = FALSE, name = "Cases n (both)", format = colFormat(digits = 0)),
						n_controls_both = colDef(show = FALSE, name = "Controls n (both)", format = colFormat(digits = 0)),
						n_ltc_covariates = colDef(show = FALSE, name = "# LTC Covariates", format = colFormat(digits = 0)),
						convergence = colDef(show = FALSE, name = "Convergence")
					),
					defaultPageSize = 20,
					searchable = TRUE,
					filterable = TRUE,
					defaultSorted = "interaction_p",
					compact = TRUE
				)
			} else {
				# Single medication results (original)
				# Add prevalence columns
				results[, case_prev := round(100 * n_cases / total_cases, 2)]
				results[, control_prev := round(100 * n_controls / total_controls, 2)]

				reactable(
					results[, .(medication, OR, CI_lower, CI_upper, p_value,
											case_prev, control_prev, n_cases, n_controls,
											n_ltc_covariates, convergence)],
					columns = list(
						medication = colDef(name = "Medication", minWidth = 200),
						OR = colDef(name = "OR", format = colFormat(digits = 3)),
						CI_lower = colDef(name = "CI Lower", format = colFormat(digits = 3)),
						CI_upper = colDef(name = "CI Upper", format = colFormat(digits = 3)),
						p_value = colDef(name = "P-value", format = colFormat(digits = 4)),
						case_prev = colDef(name = "Case %", format = colFormat(digits = 2)),
						control_prev = colDef(name = "Control %", format = colFormat(digits = 2)),
						n_cases = colDef(show = FALSE, name = "Cases (n)", format = colFormat(digits = 0)),
						n_controls = colDef(show = FALSE, name = "Controls (n)", format = colFormat(digits = 0)),
						n_ltc_covariates = colDef(show = FALSE, name = "# LTC Covariates", format = colFormat(digits = 0)),
						convergence = colDef(show = FALSE, name = "Convergence")
					),
					defaultPageSize = 20,
					searchable = TRUE,
					filterable = TRUE,
					defaultSorted = "p_value",
					compact = TRUE
				)
			}
		})

		# Download results
		output$download_results <- downloadHandler(
			filename = function() {
				if (input$include_interactions) {
					paste0("clogit_interaction_results_", Sys.Date(), ".csv")
				} else {
					paste0("clogit_results_", Sys.Date(), ".csv")
				}
			},
			content = function(file) {
				fwrite(clogit_results_r(), file)
			}
		)
	})
}