module_cca_logreg_ui <- function(id) {
	ns <- NS(id)

	accordion(
		open = FALSE,
		accordion_panel(
			title = "Logistic Regression",
			value = "logreg_panel",
			icon = bs_icon("calculator"),
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

				navset_tab(
					nav_panel("Background medications model",

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

				virtualSelectInput(
					ns("subgroup_filter"),
					"Filter analysis to specific subgroup (optional):",
					choices = list(
						"All patients" = list("all" = "All patients"),
						"Sex" = list(),
						"Ethnicity" = list(),
						"PP burden" = list(),
						"MLTC burden" = list()
					),
					selected = "all",
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

				actionButton(ns("run_logreg"), "Run Models", class = "btn-primary"),

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
					reactableOutput(ns("logreg_results_table"))
				)
			),
			nav_panel("Prescription × PP interaction model",

								p("Model interactions between the most recent prescription and polypharmacy burden groups, adjusting for selected LTCs."),

								numericInput(
									ns("pp_interaction_min_prev"),
									"Minimum medication case prevalence (%):",
									value = 2,
									min = 0.5,
									max = 20,
									step = 0.5
								),

								virtualSelectInput(
									ns("pp_interaction_meds"),
									"Select medications (leave empty for all):",
									choices = NULL,
									multiple = TRUE,
									search = TRUE,
									dropboxWrapper = "body"
								),

								actionButton(ns("run_pp_interaction"), "Run Model", class = "btn-primary"),

								hr(),

								conditionalPanel(
									condition = "output.pp_interaction_results_available",
									ns = ns,
									h4("Results"),
									div(
										class = "alert alert-light",
										style = "margin-bottom: 15px;",
										tags$strong("Model:"), " Most recent prescription × PP group interaction, adjusting for selected LTCs",
										tags$br(),
										tags$code("treatment ~ medication + pp_group + medication:pp_group + LTCs + stratum"),
										tags$br(),
										tags$small("PP groups: 0-2, 2-5, 5-10, 10+ concurrent medications")
									),
									downloadButton(ns("download_pp_interaction_results"), "Download Results", class = "btn-sm btn-secondary"),
									br(), br(),
									reactableOutput(ns("pp_interaction_results_table"))
								)
			)
			)
			)
		)
	)
}

module_cca_logreg_server <- function(id, patient_data_r, prescriptions_r, ltcs_r, cases_controls_r) {
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


		observe({
			req(patient_data_r())

			patient_data <- patient_data_r()

			subgroup_choices <- list(
				"Sex" = setNames(
					paste0("sex#", unique(patient_data$sex)),
					paste0("Sex: ", unique(patient_data$sex))
				),
				"Ethnicity" = setNames(
					paste0("eth_group#", unique(patient_data$eth_group)),
					paste0("Ethnicity: ", unique(patient_data$eth_group))
				),
				"PP burden" = setNames(
					paste0("pp_group#", unique(patient_data$pp_group)),
					paste0("PP: ", unique(patient_data$pp_group))
				),
				"MLTC burden" = setNames(
					paste0("mltc_group#", unique(patient_data$mltc_group)),
					paste0("MLTC: ", unique(patient_data$mltc_group))
				)
			)

			updateVirtualSelect(
				"subgroup_filter",
				choices = subgroup_choices
			)
		})

		# Reactive values
		filtered_ltcs_r <- reactiveVal(NULL)
		logreg_results_r <- reactiveVal(NULL)

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
					case = colDef(name = "Case %", format = colFormat(digits = 2)),
					control = colDef(name = "Control %", format = colFormat(digits = 2)),
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
		observeEvent(input$run_logreg, {
			req(filtered_ltcs_r(), patient_data_r(), prescriptions_r(), ltcs_r())

			patient_data_filtered <- patient_data_r()
			prescriptions_filtered <- prescriptions_r()
			ltcs_filtered <- ltcs_r()

			if (isTruthy(input$subgroup_filter)) {
				# Parse the subgroup selection (format: "variable#value")
				subgroup_parts <- strsplit(input$subgroup_filter, "#")[[1]]
				subgroup_var <- subgroup_parts[1]
				subgroup_val <- subgroup_parts[2]

				# Filter patient data
				patient_data_filtered <- patient_data_filtered[get(subgroup_var) == subgroup_val]

				# Filter prescriptions and LTCs to matching patients
				filtered_patids <- patient_data_filtered$patid
				prescriptions_filtered <- prescriptions_filtered[patid %in% filtered_patids]
				ltcs_filtered <- ltcs_filtered[patid %in% filtered_patids]

				showNotification(
					sprintf("Filtered to subgroup: %s = %s (%d patients)",
									subgroup_var, subgroup_val, length(filtered_patids)),
					type = "message",
					duration = 3
				)
			}

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

			selected_ltc_terms <- filtered_ltcs_r()$term

			# Check if running interaction models
			if (input$include_interactions) {
				# Generate all pairs
				med_pairs <- combn(medications_to_model, 2, simplify = FALSE)

				showNotification(
					sprintf("Running interaction models for %d medication pairs...", length(med_pairs)),
					type = "message",
					duration = NULL,
					id = "logreg_notification"
				)

				tryCatch({
					results <- run_interaction_models(
						med_pairs = med_pairs,
						selected_ltcs = selected_ltc_terms,
						patient_data = patient_data_filtered,
						prescriptions = prescriptions_filtered,
						ltcs = ltcs_filtered
					)

					logreg_results_r(results)

					output$results_available <- reactive({ TRUE })
					outputOptions(session$output, "results_available", suspendWhenHidden = FALSE)

					removeNotification("logreg_notification")
					showNotification(
						sprintf("Interaction models fitted for %d medication pairs", length(med_pairs)),
						type = "message",
						duration = 3
					)

				}, error = function(e) {
					removeNotification("logreg_notification")
					showNotification(
						paste("Error running interaction models:", e$message),
						type = "error",
						duration = 10
					)
				})
			} else {
				# Original single medication models
				showNotification("Running logistic regression models...",
												 type = "message", duration = NULL, id = "logreg_notification")

				tryCatch({
					results <- run_logistic_models(
						medications = medications_to_model,
						selected_ltcs = selected_ltc_terms,
						patient_data = patient_data_filtered,
						prescriptions = prescriptions_filtered,
						ltcs = ltcs_filtered
					)

					logreg_results_r(results)

					output$results_available <- reactive({ TRUE })
					outputOptions(session$output, "results_available", suspendWhenHidden = FALSE)

					removeNotification("logreg_notification")
					showNotification(
						sprintf("Models fitted for %d medications", length(medications_to_model)),
						type = "message",
						duration = 3
					)

				}, error = function(e) {
					removeNotification("logreg_notification")
					showNotification(
						paste("Error running models:", e$message),
						type = "error",
						duration = 10
					)
				})
			}
		})

		# Display results table
		output$logreg_results_table <- renderReactable({
			req(logreg_results_r())

			results <- copy(logreg_results_r())

			# Calculate total cases and controls for percentages
			total_cases <- patient_data_r()[treatment == 1, .N]
			total_controls <- patient_data_r()[treatment == 0, .N]

			# Check if interaction results or single medication results
			if ("med1" %in% names(results)) {
				# Interaction results
				# Update the reactable output in module_cca_logreg.R
				# Replace the interaction results section with this:

				# Create formatted OR columns
				results[, med1_OR_formatted := sprintf("%.3f (%.3f - %.3f)", med1_OR, med1_CI_lower, med1_CI_upper)]
				results[, med2_OR_formatted := sprintf("%.3f (%.3f - %.3f)", med2_OR, med2_CI_lower, med2_CI_upper)]
				results[, interaction_OR_formatted := sprintf("%.3f (%.3f - %.3f)", interaction_OR, interaction_CI_lower, interaction_CI_upper)]
				results[, combined_OR_formatted := sprintf("%.3f (%.3f - %.3f)", combined_OR, combined_CI_lower, combined_CI_upper)]

				# Interaction results
				reactable(
					results[, .(med1, med2,
											med1_OR_formatted, med1_p,
											med2_OR_formatted, med2_p,
											interaction_OR_formatted, interaction_p,
											combined_OR_formatted,
											pct_cases_both, pct_controls_both,
											n_cases_both, n_controls_both, n_ltc_covariates, convergence)],
					columns = list(
						med1 = colDef(name = "Medication 1", minWidth = 150),
						med2 = colDef(name = "Medication 2", minWidth = 150),

						# Med1 main effect
						med1_OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						med1_p = colDef(name = "P-value", format = colFormat(digits = 4)),

						# Med2 main effect
						med2_OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						med2_p = colDef(name = "P-value", format = colFormat(digits = 4)),

						# Interaction
						interaction_OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						interaction_p = colDef(name = "P-value", format = colFormat(digits = 4)),

						# Combined effect
						combined_OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),

						# Prevalence
						pct_cases_both = colDef(name = "Cases % (both)", format = colFormat(digits = 2)),
						pct_controls_both = colDef(name = "Controls % (both)", format = colFormat(digits = 2)),

						# Hidden columns
						n_cases_both = colDef(show = FALSE, name = "Cases n (both)", format = colFormat(digits = 0)),
						n_controls_both = colDef(show = FALSE, name = "Controls n (both)", format = colFormat(digits = 0)),
						n_ltc_covariates = colDef(show = FALSE, name = "# LTC Covariates", format = colFormat(digits = 0)),
						convergence = colDef(show = FALSE, name = "Convergence")
					),
					columnGroups = list(
						colGroup(name = "Med1 Main Effect (when Med2=0)", columns = c("med1_OR_formatted", "med1_p")),
						colGroup(name = "Med2 Main Effect (when Med1=0)", columns = c("med2_OR_formatted", "med2_p")),
						colGroup(name = "Interaction Effect", columns = c("interaction_OR_formatted", "interaction_p")),
						colGroup(name = "Combined Effect (Med1 + Med2)", columns = c("combined_OR_formatted"))
					),
					defaultPageSize = 20,
					searchable = TRUE,
					filterable = FALSE,
					showPageSizeOptions = TRUE,
					compact = TRUE
				)
			} else {
				# Single medication results (original)
				# Add prevalence columns
				results[, case_prev := round(100 * n_cases / total_cases, 2)]
				results[, control_prev := round(100 * n_controls / total_controls, 2)]
				results[, OR_formatted := sprintf("%.3f (%.3f - %.3f)", OR, CI_lower, CI_upper)]

				reactable(
					results[, .(medication,OR_formatted, p_value,
											case_prev, control_prev)],
					columns = list(
						medication = colDef(name = "Medication", minWidth = 200),
						OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						p_value = colDef(name = "P-value", format = colFormat(digits = 4)),
						case_prev = colDef(name = "Case %", format = colFormat(digits = 2)),
						control_prev = colDef(name = "Control %", format = colFormat(digits = 2))
					),
					defaultPageSize = 20,
					searchable = TRUE,
					filterable = FALSE,
					showPageSizeOptions = TRUE,
					defaultSorted = "medication",
					compact = TRUE
				)
			}
		})

		# Download results
		output$download_results <- downloadHandler(
			filename = function() {
				if (input$include_interactions) {
					paste0("logreg_interaction_results_", Sys.Date(), ".csv")
				} else {
					paste0("logreg_results_", Sys.Date(), ".csv")
				}
			},
			content = function(file) {
				fwrite(logreg_results_r(), file)
			}
		)

		pp_interaction_results_r <- reactiveVal(NULL)

		# Calculate frequency for recent prescriptions from cases_controls_r
		recent_presc_freq_data <- reactive({
			req(cases_controls_r(), patient_data_r())

			# cases_controls_r only has: patid, index_date, substance, group
			# We need to add strata from patient_data_r for the frequency calculations
			recent_presc_with_strata <- merge(
				cases_controls_r(),
				patient_data_r()[, .(patid, strata)],
				by = "patid",
				all.x = TRUE
			)

			freq <- calculate_frequency_stats(recent_presc_with_strata, "substance")
			create_prevalence_ratio_table(freq, "substance", data_with_group = recent_presc_with_strata)
		})

		# Observer to update medication choices for PP interaction model
		observe({
			req(recent_presc_freq_data())

			med_data <- recent_presc_freq_data()
			eligible_meds <- med_data[case >= input$pp_interaction_min_prev, substance]

			updateVirtualSelect(
				"pp_interaction_meds",
				choices = sort(eligible_meds)
			)
		})

		# Run PP interaction model
		observeEvent(input$run_pp_interaction, {
			req(filtered_ltcs_r(), patient_data_r(), cases_controls_r(), ltcs_r())

			# Get selected medications
			if (is.null(input$pp_interaction_meds) || length(input$pp_interaction_meds) == 0) {
				req(recent_presc_freq_data())
				med_data <- recent_presc_freq_data()
				medications_to_model <- med_data[case >= input$pp_interaction_min_prev, substance]

				showNotification(
					sprintf("No medications selected - running models for all %d eligible medications",
									length(medications_to_model)),
					type = "message",
					duration = 5
				)
			} else {
				medications_to_model <- input$pp_interaction_meds
			}

			if (length(medications_to_model) == 0) {
				showNotification(
					"No medications meet the prevalence threshold",
					type = "warning",
					duration = 5
				)
				return()
			}

			selected_ltc_terms <- filtered_ltcs_r()$term

			showNotification("Running PP interaction models...",
											 type = "message", duration = NULL, id = "pp_interaction_notification")

			tryCatch({
				# Run interaction models using cases_controls_r and pp_group
				results <- fit_pp_interaction_models(
					selected_ltcs = selected_ltc_terms,
					medications = medications_to_model,
					patient_data = patient_data_r(),
					recent_prescriptions = cases_controls_r(),
					ltcs = ltcs_r()
				)

				pp_interaction_results_r(results)

				output$pp_interaction_results_available <- reactive({ TRUE })
				outputOptions(output, "pp_interaction_results_available", suspendWhenHidden = FALSE)

				removeNotification("pp_interaction_notification")
				showNotification(
					sprintf("PP interaction models fitted for %d medications", length(medications_to_model)),
					type = "message",
					duration = 3
				)

			}, error = function(e) {
				removeNotification("pp_interaction_notification")
				showNotification(
					paste("Error running PP interaction models:", e$message),
					type = "error",
					duration = 10
				)
			})
		})

		# Results table
		output$pp_interaction_results_table <- renderReactable({
			req(pp_interaction_results_r())

			results <- copy(pp_interaction_results_r())

			reactable(
				results[, .(medication, pp_level,
										combined_OR_formatted, combined_OR,
										main_OR_formatted, main_p,
										interaction_OR_formatted, interaction_p,
										pct_cases, pct_controls)],
				columns = list(
					medication = colDef(name = "Medication", minWidth = 150),
					pp_level = colDef(name = "PP Group", minWidth = 120),
					combined_OR_formatted = colDef(
						name = "Combined OR (95% CI)",
						minWidth = 150,
						style = function(value, index) {
							# Get the actual OR value for conditional formatting
							or_value <- results[index, combined_OR]
							if (is.na(or_value)) return(NULL)
							# Highlight strong effects
							if (or_value > 1.5 || or_value < 0.67) {
								list(background = "#fff3cd", fontWeight = "bold")
							}
						}
					),
					combined_OR = colDef(show = FALSE),  # Hidden column used for sorting/filtering
					main_OR_formatted = colDef(name = "Main OR (95% CI)", minWidth = 150),
					main_p = colDef(
						name = "Main P",
						format = colFormat(digits = 4),
						style = function(value) {
							if (is.na(value)) return(NULL)
							if (value < 0.05) list(fontWeight = "bold")
						}
					),
					interaction_OR_formatted = colDef(name = "Interaction OR (95% CI)", minWidth = 150),
					interaction_p = colDef(
						name = "Interaction P",
						format = colFormat(digits = 4),
						style = function(value) {
							if (is.na(value)) return(NULL)
							if (value < 0.05) list(fontWeight = "bold", color = "#d9534f")
						}
					),
					pct_cases = colDef(name = "Cases %", format = colFormat(digits = 1)),
					pct_controls = colDef(name = "Controls %", format = colFormat(digits = 1))
				),
				filterable = FALSE,
				searchable = FALSE,
				defaultPageSize = 20,
				showPageSizeOptions = TRUE,
				compact = TRUE,
				defaultColDef = colDef(minWidth = 100),
				defaultSorted = list(combined_OR = "desc")
			)
		})

		# Download handler
		output$download_pp_interaction_results <- downloadHandler(
			filename = function() {
				paste0("pp_interaction_results_", Sys.Date(), ".csv")
			},
			content = function(file) {
				results <- pp_interaction_results_r()

				# Select columns for export (include both formatted and raw values)
				export_data <- results[, .(
					medication,
					pp_level,
					combined_OR_formatted,
					combined_OR,
					combined_CI_lower,
					combined_CI_upper,
					main_OR_formatted,
					main_OR,
					main_CI_lower,
					main_CI_upper,
					main_p,
					interaction_OR_formatted,
					interaction_OR,
					interaction_CI_lower,
					interaction_CI_upper,
					interaction_p,
					pct_cases,
					pct_controls,
					n_cases,
					n_controls,
					n_ltc_covariates,
					convergence
				)]

				fwrite(export_data, file)
			}
		)

	})
}