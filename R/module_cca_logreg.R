# R/module_cca_logreg.R
# Logistic Regression Module (Refactored)

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
						"Polypharmacy burden (numeric)" = "pp",
						"Polypharmacy burden group" = "pp_group",
						"MLTC burden (numeric)" = "n_ltc",
						"MLTC burden group" = "mltc_group"
					),
					selected = NULL
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
						ns("recent_presc_min_prev"),
						"Minimum recent prescription prevalence (%):",
						value = 2,
						min = 0.5,
						max = 20,
						step = 0.5
					),
					virtualSelectInput(
						ns("selected_recent_presc"),
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
						"Recent prescription (main effects)" = "recent_main",
						"Recent prescription × PP burden (interaction)" = "recent_pp",
						"Recent prescription × background med (interaction)" = "recent_background"
					),
					selected = "background_main"
				),

				# Model-specific options (dynamic)
				uiOutput(ns("model_specific_options")),

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

					# Subgroup info display
					uiOutput(ns("subgroup_info")),

					# Model description (dynamic based on model type)
					uiOutput(ns("model_description")),

					downloadButton(ns("download_results"), "Download Results", class = "btn-sm btn-secondary"),
					br(), br(),

					# Single dynamic results table
					reactableOutput(ns("results_table"))
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
		subgroup_info_r <- reactiveVal(NULL)  # Store subgroup filter info

		# Reset results when model type changes
		observeEvent(input$model_type, {
			model_results_r(NULL)
			subgroup_info_r(NULL)
			output$results_available <- reactive({ FALSE })
			outputOptions(output, "results_available", suspendWhenHidden = FALSE)
		})

		# ===== STEP 1: FILTER LTCs =====
		observeEvent(input$filter_ltcs, {
			req(ltc_freq_data())

			ltc_data <- ltc_freq_data()
			filtered <- ltc_data[
				case >= input$ltc_min_prev &
					control >= input$ltc_min_prev &
					(OR >= input$ltc_min_or | OR < 1/input$ltc_min_or) &
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
				case >= input$recent_presc_min_prev &
					control >= input$recent_presc_min_prev,
				substance
			]
			updateVirtualSelect("selected_recent_presc", choices = sort(eligible_recent))
		})

		# ===== STEP 3: DYNAMIC MODEL OPTIONS =====
		output$model_specific_options <- renderUI({
			req(input$model_type)

			switch(input$model_type,
						 background_pairwise = div(
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

						 recent_main = div(
						 	checkboxInput(
						 		ns("group_recent_meds"),
						 		"Group all selected recent prescriptions into a single indicator variable",
						 		value = FALSE
						 	),
						 	div(
						 		class = "alert alert-info",
						 		style = "margin-top: 5px; margin-bottom: 10px;",
						 		tags$small(
						 			"When checked, a single model will be run testing 'any of the selected medications' ",
						 			"instead of separate models for each medication."
						 		)
						 	)
						 ),

						 recent_pp = div(
						 	checkboxInput(
						 		ns("group_recent_meds"),
						 		"Group all selected recent prescriptions into a single indicator variable",
						 		value = FALSE
						 	),
						 	div(
						 		class = "alert alert-info",
						 		style = "margin-top: 5px; margin-bottom: 10px;",
						 		tags$small(
						 			"When checked, a single model will be run testing 'any of the selected medications' ",
						 			"instead of separate models for each medication."
						 		)
						 	)
						 ),

						 recent_background = div(
						 	checkboxInput(
						 		ns("group_recent_meds"),
						 		"Group all selected recent prescriptions into a single indicator variable",
						 		value = FALSE
						 	),
						 	div(
						 		class = "alert alert-info",
						 		style = "margin-top: 5px; margin-bottom: 10px;",
						 		tags$small(
						 			"When checked, a single model will be run testing 'any of the selected medications' ",
						 			"instead of separate models for each medication."
						 		)
						 	),
						 	div(
						 		class = "alert alert-info",
						 		style = "margin-top: 10px;",
						 		tags$strong("Note: "),
						 		"You must select specific medications from both recent prescriptions and background medications above. ",
						 		"All possible recent × background pairs will be tested."
						 	)
						 ),

						 NULL  # No extra options for background_main
			)
		})

		# ===== DYNAMIC MODEL DESCRIPTION =====
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
				recent_main = list(
					text = "Recent prescription main effects, adjusting for matching variables and selected LTCs",
					formula = "treatment ~ recent_prescription + LTCs + stratum"
				),
				recent_pp = list(
					text = "Most recent prescription × PP group interaction, adjusting for selected LTCs",
					formula = "treatment ~ medication + pp_group + medication:pp_group + LTCs + stratum"
				),
				recent_background = list(
					text = "Recent prescription × background medication interaction, adjusting for selected LTCs",
					formula = "treatment ~ recent_presc + background_med + recent_presc:background_med + LTCs + stratum"
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

		# ===== SUBGROUP INFO DISPLAY =====
		output$subgroup_info <- renderUI({
			info <- subgroup_info_r()

			if (is.null(info)) {
				return(NULL)
			}

			# Create readable variable names
			var_name <- switch(info$variable,
												 sex = "Sex",
												 eth_group = "Ethnicity",
												 pp_group = "Polypharmacy burden",
												 mltc_group = "MLTC burden",
												 info$variable
			)

			div(
				class = "alert alert-info",
				style = "margin-bottom: 15px;",
				tags$strong("Subgroup Analysis:"),
				sprintf(" %s = %s", var_name, info$value),
				tags$br(),
				sprintf("Total patients: %s (%s cases, %s controls)",
								prettyNum(info$n_patients, big.mark = ","),
								prettyNum(info$n_cases, big.mark = ","),
								prettyNum(info$n_controls, big.mark = ","))
			)
		})

		# ===== RUN MODEL =====
		observeEvent(input$run_model, {
			req(filtered_ltcs_r(), patient_data_r(), ltcs_r())

			# Validate inputs
			if (!validate_model_inputs(input$model_type, input)) {
				return()
			}

			# Prepare filtered data
			filtered_data <- prepare_filtered_data(
				input,
				patient_data_r(),
				prescriptions_r(),
				ltcs_r(),
				cases_controls_r()
			)

			selected_ltc_terms <- filtered_ltcs_r()$term

			# Run model
			showNotification(
				"Running models...",
				type = "message",
				duration = NULL,
				id = "model_notification"
			)

			tryCatch({
				results <- execute_model(
					input$model_type,
					input,
					filtered_data,
					selected_ltc_terms
				)

				model_results_r(results)

				# Store subgroup information
				if (isTruthy(input$subgroup_filter) && input$subgroup_filter != "all") {
					subgroup_parts <- strsplit(input$subgroup_filter, "#")[[1]]
					subgroup_info_r(list(
						variable = subgroup_parts[1],
						value = subgroup_parts[2],
						n_patients = nrow(filtered_data$patient_data),
						n_cases = filtered_data$patient_data[treatment == 1, .N],
						n_controls = filtered_data$patient_data[treatment == 0, .N]
					))
				} else {
					subgroup_info_r(NULL)
				}

				output$results_available <- reactive({ TRUE })
				outputOptions(output, "results_available", suspendWhenHidden = FALSE)

				removeNotification("model_notification")

				# Show success message
				showNotification(
					get_success_message(results, input$model_type),
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

		# ===== RENDER RESULTS TABLE (UNIFIED) =====
		output$results_table <- renderReactable({
			req(model_results_r())
			req(input$model_type)

			# Additional safety check: ensure results exist and match model type
			results <- model_results_r()
			if (is.null(results) || nrow(results) == 0) {
				return(NULL)
			}

			render_logreg_results_table(results, input$model_type)
		})

		# ===== DOWNLOAD RESULTS =====
		output$download_results <- downloadHandler(
			filename = function() {
				model_name <- switch(
					input$model_type,
					background_main = "background_main",
					background_pairwise = "background_pairwise",
					recent_main = "recent_main",
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