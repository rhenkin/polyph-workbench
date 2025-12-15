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

				# ===== STEP 1: LTC FILTERING (SHARED) =====
				h4("Step 1: Select covariates (LTCs)"),

				layout_columns(
					col_widths = c(4, 4, 4),
					numericInput(
						ns("ltc_min_prev"),
						"Minimum LTC prevalence (%):",
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

				h4("Select additional covariates (optional)"),

				checkboxGroupInput(
					ns("selected_covariates"),
					"Additional covariates:",
					choices = list(
						"Ethnicity group" = "eth_group",
						"IMD Quintiles" = "imd_quintile",
						"Polypharmacy burden group" = "pp_group",
						"MLTC burden group" = "mltc_group"
					),
					selected = NULL  # None selected by default, or select some defaults
				),

				hr(),

				# ===== STEP 2: MEDICATION SELECTION (SHARED) =====
				h4("Step 2: Select medications"),

				# Background medications selection
				h5("Background medications"),
				layout_columns(
					col_widths = c(6, 6),
					numericInput(
						ns("background_med_min_prev"),
						"Minimum background medication prevalence (%):",
						value = 2,
						min = 0.5,
						max = 20,
						step = 0.5
					),
					virtualSelectInput(
						ns("selected_background_meds"),
						"Select background medications:",
						choices = NULL,
						multiple = TRUE,
						search = TRUE,
						dropboxWrapper = "body",
						placeholder = "Select medications..."
					)
				),

				# Recent prescriptions selection
				h5("Most recent prescriptions"),
				layout_columns(
					col_widths = c(6, 6),
					numericInput(
						ns("recent_rx_min_prev"),
						"Minimum recent prescription prevalence (%):",
						value = 2,
						min = 0.5,
						max = 20,
						step = 0.5
					),
					virtualSelectInput(
						ns("selected_recent_rx"),
						"Select recent prescriptions:",
						choices = NULL,
						multiple = TRUE,
						search = TRUE,
						dropboxWrapper = "body",
						placeholder = "Select medications..."
					)
				),

				hr(),

				# ===== STEP 3: MODEL SELECTION =====
				h4("Step 3: Choose model type"),

				radioButtons(
					ns("model_type"),
					"Select model:",
					choices = list(
						"Background medications (main effects)" = "background_main",
						"Background med × background med (pairwise interactions)" = "background_pairwise",
						"Recent prescription × PP burden (interaction)" = "recent_pp",
						"Recent prescription × background med (interaction)" = "recent_background"
					),
					selected = "background_main"
				),

				# Model-specific options
				conditionalPanel(
					condition = "input.model_type == 'background_pairwise'",
					ns = ns,
					numericInput(
						ns("interaction_min_coprescription_prev"),
						"Minimum co-prescription prevalence (%):",
						value = 1,
						min = 0.1,
						max = 10,
						step = 0.1
					),
					div(
						class = "alert alert-info",
						style = "margin-top: 10px;",
						tags$strong("Note: "),
						"You must select specific background medications above. All possible pairs will be tested."
					)
				),

				conditionalPanel(
					condition = "input.model_type == 'recent_background'",
					ns = ns,
					div(
						class = "alert alert-info",
						style = "margin-top: 10px;",
						tags$strong("Note: "),
						"You must select specific medications from both recent prescriptions and background medications above. ",
						"All possible recent × background pairs will be tested."
					)
				),

				# Subgroup filter (optional, for all models)
				virtualSelectInput(
					ns("subgroup_filter"),
					"Filter analysis to specific subgroup (optional):",
					choices = list(
						"Sex" = list(),
						"Ethnicity" = list(),
						"PP burden" = list(),
						"MLTC burden" = list()
					),
					selected = "all",
					search = TRUE,
					hideClearButton = FALSE,
					dropboxWrapper = "body"
				),

				hr(),

				# ===== RUN MODEL =====
				actionButton(ns("run_model"), "Run Model", class = "btn-primary"),

				hr(),

				# ===== RESULTS =====
				conditionalPanel(
					condition = "output.results_available",
					ns = ns,
					h4("Results"),

					# Model description (dynamic based on model type)
					uiOutput(ns("model_description")),

					downloadButton(ns("download_results"), "Download Results", class = "btn-sm btn-secondary"),
					br(), br(),

					# Results tables (conditional display based on model type)
					conditionalPanel(
						condition = "input.model_type == 'background_main'",
						ns = ns,
						reactableOutput(ns("main_effects_table"))
					),

					conditionalPanel(
						condition = "input.model_type == 'background_pairwise'",
						ns = ns,
						reactableOutput(ns("pairwise_table"))
					),

					conditionalPanel(
						condition = "input.model_type == 'recent_pp'",
						ns = ns,
						reactableOutput(ns("pp_interaction_table"))
					),

					conditionalPanel(
						condition = "input.model_type == 'recent_background'",
						ns = ns,
						reactableOutput(ns("recent_background_table"))
					)
				)
			)
		)
	)
}

module_cca_logreg_server <- function(id, patient_data_r, prescriptions_r, ltcs_r, cases_controls_r) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# ===== FREQUENCY DATA (SHARED) =====
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

		recent_presc_freq_data <- reactive({
			req(cases_controls_r(), patient_data_r())
			recent_presc_with_strata <- merge(
				cases_controls_r(),
				patient_data_r()[, .(patid, strata)],
				by = "patid",
				all.x = TRUE
			)
			freq <- calculate_frequency_stats(recent_presc_with_strata, "substance")
			create_prevalence_ratio_table(freq, "substance", data_with_group = recent_presc_with_strata)
		})

		# ===== UPDATE SUBGROUP CHOICES =====
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
					paste0("pp_group#", levels(patient_data$pp_group)),
					paste0("PP: ", levels(patient_data$pp_group))
				),
				"MLTC burden" = setNames(
					paste0("mltc_group#", levels(patient_data$mltc_group)),
					paste0("MLTC: ", levels(patient_data$mltc_group))
				)
			)

			updateVirtualSelect("subgroup_filter", choices = subgroup_choices)
		})

		# ===== REACTIVE VALUES =====
		filtered_ltcs_r <- reactiveVal(NULL)
		model_results_r <- reactiveVal(NULL)

		# ===== STEP 1: FILTER LTCs =====
		observeEvent(input$filter_ltcs, {
			req(ltc_freq_data())

			ltc_data <- ltc_freq_data()
			filtered <- ltc_data[
				case >= input$ltc_min_prev &
					control >= input$ltc_min_prev &
					OR >= input$ltc_min_or &
					!is.na(OR)
			]

			filtered <- filtered[order(-OR)][1:min(input$ltc_max_count, nrow(filtered))]
			filtered_ltcs_r(filtered)

			output$ltcs_filtered <- reactive({ TRUE })
			outputOptions(session$output, "ltcs_filtered", suspendWhenHidden = FALSE)

			showNotification(
				sprintf("%d LTCs selected for adjustment", nrow(filtered)),
				type = "message",
				duration = 3
			)
		})

		output$ltc_count_text <- renderText({
			req(filtered_ltcs_r())
			sprintf("%d LTCs will be included as covariates", nrow(filtered_ltcs_r()))
		})

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

		# ===== STEP 2: UPDATE MEDICATION CHOICES =====
		observe({
			req(med_freq_data())
			med_data <- med_freq_data()
			eligible_meds <- med_data[
				case >= input$background_med_min_prev &
					control >= input$background_med_min_prev,
				substance
			]
			updateVirtualSelect("selected_background_meds", choices = sort(eligible_meds))
		})

		observe({
			req(recent_presc_freq_data())
			recent_data <- recent_presc_freq_data()
			eligible_recent <- recent_data[
				case >= input$recent_rx_min_prev &
					control >= input$recent_rx_min_prev,
				substance
			]
			updateVirtualSelect("selected_recent_rx", choices = sort(eligible_recent))
		})

		# ===== STEP 3: DYNAMIC MODEL DESCRIPTION =====
		output$model_description <- renderUI({
			req(input$model_type)

			descriptions <- list(
				background_main = list(
					text = "Logistic regression adjusting for matching variables and selected LTCs",
					formula = "treatment ~ medication + LTCs + stratum"
				),
				background_pairwise = list(
					text = "Pairwise interaction model adjusting for matching variables and selected LTCs",
					formula = "treatment ~ med1 + med2 + med1:med2 + LTCs + stratum"
				),
				recent_pp = list(
					text = "Most recent prescription × PP group interaction, adjusting for selected LTCs",
					formula = "treatment ~ medication + pp_group + medication:pp_group + LTCs + stratum"
				),
				recent_background = list(
					text = "Recent prescription × background medication interaction, adjusting for selected LTCs",
					formula = "treatment ~ recent_rx + background_med + recent_rx:background_med + LTCs + stratum"
				)
			)

			desc <- descriptions[[input$model_type]]
			div(
				class = "alert alert-light",
				style = "margin-bottom: 15px;",
				tags$strong("Model:"), " ", desc$text,
				tags$br(),
				tags$code(desc$formula)
			)
		})

		# ===== RUN MODEL =====
		observeEvent(input$run_model, {
			req(filtered_ltcs_r(), patient_data_r(), ltcs_r())

			# Validate inputs based on model type
			validation_passed <- switch(
				input$model_type,
				background_main = {
					if (is.null(input$selected_background_meds) || length(input$selected_background_meds) == 0) {
						showNotification(
							"Please select at least one background medication",
							type = "warning",
							duration = 5
						)
						FALSE
					} else TRUE
				},
				background_pairwise = {
					if (is.null(input$selected_background_meds) || length(input$selected_background_meds) < 2) {
						showNotification(
							"Please select at least two background medications for pairwise interactions",
							type = "warning",
							duration = 5
						)
						FALSE
					} else TRUE
				},
				recent_pp = {
					if (is.null(input$selected_recent_rx) || length(input$selected_recent_rx) == 0) {
						showNotification(
							"Please select at least one recent prescription",
							type = "warning",
							duration = 5
						)
						FALSE
					} else TRUE
				},
				recent_background = {
					if (is.null(input$selected_recent_rx) || length(input$selected_recent_rx) == 0 ||
							is.null(input$selected_background_meds) || length(input$selected_background_meds) == 0) {
						showNotification(
							"Please select at least one medication from both recent prescriptions and background medications",
							type = "warning",
							duration = 5
						)
						FALSE
					} else TRUE
				}
			)

			if (!validation_passed) return()

			# Prepare filtered data
			patient_data_filtered <- patient_data_r()
			prescriptions_filtered <- prescriptions_r()
			ltcs_filtered <- ltcs_r()
			case_controls <- cases_controls_r()

			# Apply subgroup filtering if selected
			if (isTruthy(input$subgroup_filter)) {
				subgroup_parts <- strsplit(input$subgroup_filter, "#")[[1]]
				subgroup_var <- subgroup_parts[1]
				subgroup_val <- subgroup_parts[2]

				patient_data_filtered <- patient_data_filtered[get(subgroup_var) == subgroup_val]
				filtered_patids <- patient_data_filtered$patid

				prescriptions_filtered <- prescriptions_filtered[patid %in% filtered_patids]
				ltcs_filtered <- ltcs_filtered[patid %in% filtered_patids]
				case_controls <- case_controls[patid %in% filtered_patids]

				showNotification(
					sprintf("Filtered to subgroup: %s = %s (%d patients)",
									subgroup_var, subgroup_val, length(filtered_patids)),
					type = "message",
					duration = 3
				)
			}

			selected_ltc_terms <- filtered_ltcs_r()$term
			selected_covariates <- input$selected_covariates

			# Run appropriate model based on model_type
			showNotification(
				"Running models...",
				type = "message",
				duration = NULL,
				id = "model_notification"
			)

			tryCatch({
				results <- switch(
					input$model_type,

					# Background medications main effects
					background_main = {
						run_logistic_models(
							medications = input$selected_background_meds,
							selected_ltcs = selected_ltc_terms,
							selected_covariates = selected_covariates,
							patient_data = patient_data_filtered,
							prescriptions = prescriptions_filtered,
							ltcs = ltcs_filtered
						)
					},

					# Background × background pairwise interactions
					background_pairwise = {
						med_pairs <- combn(input$selected_background_meds, 2, simplify = FALSE)

						# Filter by co-prescription prevalence
						med_pairs <- filter_pairs_by_coprescription(
							med_pairs = med_pairs,
							min_coprescription_prev = input$interaction_min_coprescription_prev,
							patient_data = patient_data_filtered,
							prescriptions = prescriptions_filtered
						)

						if (length(med_pairs) == 0) {
							showNotification(
								sprintf("No medication pairs meet the co-prescription prevalence threshold of %.1f%%",
												input$interaction_min_coprescription_prev),
								type = "warning",
								duration = 5
							)
							return()
						}

						run_interaction_models(
							med_pairs = med_pairs,
							selected_ltcs = selected_ltc_terms,
							selected_covariates = selected_covariates,
							patient_data = patient_data_filtered,
							prescriptions = prescriptions_filtered,
							ltcs = ltcs_filtered
						)
					},

					# Recent prescription × PP burden
					recent_pp = {
						selected_covariates_filtered <- setdiff(selected_covariates, "pp_group")

						fit_pp_interaction_models(
							selected_ltcs = selected_ltc_terms,
							selected_covariates = selected_covariates_filtered,
							medications = input$selected_recent_rx,
							patient_data = patient_data_filtered,
							recent_prescriptions = case_controls,
							ltcs = ltcs_filtered
						)
					},

					# Recent prescription × background medication
					recent_background = {
						# Create all pairs of recent × background
						recent_background_pairs <- expand.grid(
							recent = input$selected_recent_rx,
							background = input$selected_background_meds,
							stringsAsFactors = FALSE
						)

						# Run interaction models for each pair
						run_recent_background_interaction_models(
							recent_meds = input$selected_recent_rx,
							background_meds = input$selected_background_meds,
							selected_ltcs = selected_ltc_terms,
							selected_covariates = selected_covariates,
							patient_data = patient_data_filtered,
							recent_prescriptions = case_controls,
							prescriptions = prescriptions_filtered,
							ltcs = ltcs_filtered
						)
					}
				)

				model_results_r(results)

				output$results_available <- reactive({ TRUE })
				outputOptions(output, "results_available", suspendWhenHidden = FALSE)

				removeNotification("model_notification")

				# Show success message with count
				n_results <- switch(
					input$model_type,
					background_main = nrow(results),
					background_pairwise = nrow(results),
					recent_pp = length(unique(results$medication)),
					recent_background = nrow(results)
				)

				showNotification(
					sprintf("Models completed: %d results", n_results),
					type = "message",
					duration = 3
				)

			}, error = function(e) {
				removeNotification("model_notification")
				showNotification(
					paste("Error running models:", e$message),
					type = "error",
					duration = 10
				)
			})
		})

		# ===== RENDER RESULTS TABLES =====

		# Background main effects table
		output$main_effects_table <- renderReactable({
			req(model_results_r())
			req(input$model_type == "background_main")

			results <- copy(model_results_r())
			results[, OR_formatted := sprintf("%.2f (%.2f-%.2f)", OR, CI_lower, CI_upper)]

			reactable(
				results[, .(medication, OR_formatted, p_value, n_cases, n_controls)],
				columns = list(
					medication = colDef(name = "Medication", minWidth = 200),
					OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
					p_value = colDef(name = "P-value", format = colFormat(digits = 4)),
					n_cases = colDef(name = "N Cases"),
					n_controls = colDef(name = "N Controls")
				),
				defaultPageSize = 20,
				searchable = TRUE,
				showPageSizeOptions = TRUE,
				defaultSorted = "medication",
				compact = TRUE
			)
		})

		# Pairwise interactions table
		output$pairwise_table <- renderReactable({
			req(model_results_r())
			req(input$model_type == "background_pairwise")

			results <- copy(model_results_r())
			results[, med1_OR_formatted := sprintf("%.3f (%.3f-%.3f)", med1_OR, med1_CI_lower, med1_CI_upper)]
			results[, med2_OR_formatted := sprintf("%.3f (%.3f-%.3f)", med2_OR, med2_CI_lower, med2_CI_upper)]
			results[, interaction_OR_formatted := sprintf("%.3f (%.3f-%.3f)", interaction_OR, interaction_CI_lower, interaction_CI_upper)]
			results[, combined_OR_formatted := sprintf("%.3f (%.3f-%.3f)", combined_OR, combined_CI_lower, combined_CI_upper)]

			reactable(
				results[, .(med1, med2, med1_OR_formatted, med2_OR_formatted,
										interaction_OR_formatted, combined_OR_formatted,
										interaction_p, pct_cases_both, pct_controls_both)],
				columns = list(
					med1 = colDef(name = "Medication 1", minWidth = 150),
					med2 = colDef(name = "Medication 2", minWidth = 150),
					med1_OR_formatted = colDef(name = "Med1 OR (95% CI)", minWidth = 140),
					med2_OR_formatted = colDef(name = "Med2 OR (95% CI)", minWidth = 140),
					interaction_OR_formatted = colDef(name = "Interaction OR (95% CI)", minWidth = 160),
					combined_OR_formatted = colDef(name = "Combined OR (95% CI)", minWidth = 140),
					interaction_p = colDef(name = "Interaction P", format = colFormat(digits = 4)),
					pct_cases_both = colDef(name = "Cases Both %", format = colFormat(digits = 2)),
					pct_controls_both = colDef(name = "Controls Both %", format = colFormat(digits = 2))
				),
				defaultPageSize = 20,
				searchable = TRUE,
				showPageSizeOptions = TRUE,
				compact = TRUE
			)
		})

		# PP interaction table
		output$pp_interaction_table <- renderReactable({
			req(model_results_r())
			req(input$model_type == "recent_pp")

			results <- copy(model_results_r())

			reactable(
				results[, .(medication, pp_level, main_OR_formatted, interaction_OR_formatted,
										combined_OR_formatted, interaction_p, case_prev, control_prev)],
				columns = list(
					medication = colDef(name = "Medication", minWidth = 150),
					pp_level = colDef(name = "PP Group", minWidth = 100),
					main_OR_formatted = colDef(name = "Main OR (95% CI)", minWidth = 140),
					interaction_OR_formatted = colDef(name = "Interaction OR (95% CI)", minWidth = 160),
					combined_OR_formatted = colDef(name = "Combined OR (95% CI)", minWidth = 140),
					interaction_p = colDef(name = "Interaction P", format = colFormat(digits = 4)),
					case_prev = colDef(name = "Case %", format = colFormat(digits = 2)),
					control_prev = colDef(name = "Control %", format = colFormat(digits = 2))
				),
				defaultPageSize = 20,
				searchable = TRUE,
				showPageSizeOptions = TRUE,
				compact = TRUE,
				groupBy = "medication"
			)
		})

		# Recent × background interaction table
		output$recent_background_table <- renderReactable({
			req(model_results_r())
			req(input$model_type == "recent_background")

			results <- copy(model_results_r())
			results[, recent_OR_formatted := sprintf("%.3f (%.3f-%.3f)", recent_OR, recent_CI_lower, recent_CI_upper)]
			results[, background_OR_formatted := sprintf("%.3f (%.3f-%.3f)", background_OR, background_CI_lower, background_CI_upper)]
			results[, interaction_OR_formatted := sprintf("%.3f (%.3f-%.3f)", interaction_OR, interaction_CI_lower, interaction_CI_upper)]
			results[, combined_OR_formatted := sprintf("%.3f (%.3f-%.3f)", combined_OR, combined_CI_lower, combined_CI_upper)]

			reactable(
				results[, .(recent_med, background_med, recent_OR_formatted, background_OR_formatted,
										interaction_OR_formatted, combined_OR_formatted,
										interaction_p, pct_cases_both, pct_controls_both)],
				columns = list(
					recent_med = colDef(name = "Recent Prescription", minWidth = 150),
					background_med = colDef(name = "Background Medication", minWidth = 150),
					recent_OR_formatted = colDef(name = "Recent OR (95% CI)", minWidth = 140),
					background_OR_formatted = colDef(name = "Background OR (95% CI)", minWidth = 140),
					interaction_OR_formatted = colDef(name = "Interaction OR (95% CI)", minWidth = 160),
					combined_OR_formatted = colDef(name = "Combined OR (95% CI)", minWidth = 140),
					interaction_p = colDef(name = "Interaction P", format = colFormat(digits = 4)),
					pct_cases_both = colDef(name = "Cases Both %", format = colFormat(digits = 2)),
					pct_controls_both = colDef(name = "Controls Both %", format = colFormat(digits = 2))
				),
				defaultPageSize = 20,
				searchable = TRUE,
				showPageSizeOptions = TRUE,
				compact = TRUE
			)
		})

		# ===== DOWNLOAD RESULTS =====
		output$download_results <- downloadHandler(
			filename = function() {
				model_name <- switch(
					input$model_type,
					background_main = "background_main",
					background_pairwise = "background_pairwise",
					recent_pp = "recent_pp_interaction",
					recent_background = "recent_background_interaction"
				)
				paste0("logreg_", model_name, "_", Sys.Date(), ".csv")
			},
			content = function(file) {
				fwrite(model_results_r(), file)
			}
		)

	})
}