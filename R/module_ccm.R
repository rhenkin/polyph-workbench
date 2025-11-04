module_ccm_ui <- function(id) {
	ns <- NS(id)
	nav_panel("Case-control Matching",
						fluidRow(
							column(6,
										 textInput(ns("study_name"), "Study name:", value = "my_study"),
										 numericInput(ns("pred_window"), "Prediction window:",
										 						 value = 30, min = 1, max = 365),
										 numericInput(ns("match_ratio"), "Control:Case ratio:",
										 						 value = 4, min = 1, max = 10, step = 1),
										 div("Risk-set matching using sex, binned age at prescription and binned time since multimorbidity"),
										 actionButton(ns("create_cohort"), "Create matched cohort", class = "btn-primary")
							),
							column(6,
										 actionButton(ns("save_study"), "Save study", class = "btn-success"),
										 br(), br(),
										 verbatimTextOutput(ns("save_status"))
							)
						),
						br(),
						conditionalPanel(
							condition = "output.show_results == true",ns = ns,
							card(
							card_header("Matching Results"),
							card_body(
								layout_columns(col_widths = c(2,3,4,4),
									verticalLayout(
										value_box(
											showcase_layout = showcase_left_center(width = "200px", max_height = "150px"),
											height = "150px",
											title = "Cases",
											value = textOutput(ns("vb_cases_n")),
											theme = "red",
											showcase = bs_icon("people-fill")
										),
										value_box(
											showcase_layout = showcase_left_center(width = "200px", max_height = "150px"),
											height = "150px",
											title = "Controls",
											value = textOutput(ns("vb_controls_n")),
											theme = "blue",
											showcase = bs_icon("people-fill")
										)
									),
									card(
										card_header("Summary"),
										card_body(
											vegawidgetOutput(ns("matched_age_dist")),
											vegawidgetOutput(ns("matched_time_dist"))
										)
									)
								)
							)
						)
						)
	)
}

module_ccm_server <- function(id, patient_data, outcome_prescriptions, ltc_data,
															study_dir = "studies") {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Load master risk pool once when module is activated
		master_risk_pool_dataset <- reactive({
			req(file.exists("../data/master_risk_pool_smaller.parquet"))
			showNotification("Opening master risk pool...", duration = 10, type = "message")
			arrow::open_dataset("../data/master_risk_pool_smaller.parquet")
		})


		# Reactive values to store results
		cases_r <- reactiveVal(NULL)
		controls_r <- reactiveVal(NULL)
		prepared_study_data_r <- reactiveVal(NULL)
		# Main workflow: Create matched cohort
		observeEvent(input$create_cohort, {
			req(input$study_name)
			req(outcome_prescriptions())
			req(master_risk_pool_dataset())

			# Show progress
			progress <- Progress$new()
			on.exit(progress$close())
			progress$set(message = "Creating matched cohort", value = 0)

			tryCatch({
				# Call the main workflow

				result <- create_matched_cohort_workflow(
					outcome_prescriptions = outcome_prescriptions(),
					master_risk_pool_dataset = master_risk_pool_dataset(),
					patient_data = patient_data(),
					ltc_data = ltc_data(),
					pred_window = input$pred_window,
					match_ratio = input$match_ratio,
					progress = progress
				)

				# Store results
				cases_r(result$cases)
				controls_r(result$controls)

				# Update UI
				output$show_results <- reactive({ TRUE })
				outputOptions(output, "show_results", suspendWhenHidden = FALSE)

				progress$set(value = 1)
				showNotification("Matched cohort created successfully!",
												 type = "message", duration = 3)

			}, error = function(e) {
				showNotification(paste("Error creating cohort:", e$message),
												 type = "error", duration = 10)
				message("Error: ", e$message)
				print(traceback())
			})
		})

		observe({
			req(!is.null(cases_r()), !is.null(controls_r()))

			cases <- cases_r()
			controls <- controls_r()

			output$vb_cases_n <- renderText({
				prettyNum(uniqueN(cases$patid), big.mark = ",")
			})

			output$vb_controls_n <- renderText({
				prettyNum(uniqueN(controls$patid), big.mark = ",")
			})


			output$matched_age_dist <- renderVegawidget({

				summary_cases <- round(summary(cases$index_age), digits = 2)
				summary_controls <- round(summary(controls$age_at_rx), digits = 2)

				df <- rbind(
					as.data.table(as.list(summary_cases))[, group := "Cases"],
					as.data.table(as.list(summary_controls))[, group := "Controls"]
				)
				colnames(df) <- c("lower", "q1", "median", "mean", "q3", "upper", "group")

				spec <- create_boxplot_spec(df, "Age at prescription")
				spec |> as_vegaspec()

			})

			output$matched_time_dist <- renderVegawidget({

				summary_cases <- round(summary(cases$time_since_first_presc), digits = 2)
				summary_controls <- round(summary(controls$time_since_first_presc), digits = 2)

				df <- rbind(
					as.data.table(as.list(summary_cases))[, group := "Cases"],
					as.data.table(as.list(summary_controls))[, group := "Controls"]
				)
				colnames(df) <- c("lower", "q1", "median", "mean", "q3", "upper", "group")

				spec <- create_boxplot_spec(df,  "Time since first prescription")
				spec |> as_vegaspec()

			})
		}, suspended = TRUE)

		# Prepare study data in memory when matching completes
		observe({
			req(cases_r(), controls_r(), input$study_name)

			tryCatch({
				study_data <- prepare_study_data(
					study_name = input$study_name,
					cases = cases_r(),
					controls = controls_r(),
					gold_patient = gold_patient,
					gold_cp = gold_cp,
					gold_ltc = gold_ltc,
					outcome_prescriptions = outcome_prescriptions()
				)

				prepared_study_data_r(study_data)
				message("Study data prepared in memory: ", input$study_name)

			}, error = function(e) {
				message("Error preparing study data: ", e$message)
			})
		}, suspended = TRUE)

		# Save study
		observeEvent(input$save_study, {
			req(cases_r(), controls_r(), input$study_name)

			progress <- Progress$new()
			on.exit(progress$close())
			progress$set(message = "Saving study...", value = 0)

			tryCatch({
				save_matched_datasets(
					study_name = input$study_name,
					cases = cases_r(),
					controls = controls_r(),
					gold_patient = gold_patient,
					gold_cp = gold_cp,
					gold_ltc = gold_ltc,
					outcome_prescriptions = outcome_prescriptions(),
					study_dir = study_dir
				)

				showNotification(paste("Study", input$study_name, "saved successfully!"),
												 type = "message", duration = 3)

			}, error = function(e) {
				showNotification(paste("Error saving study:", e$message), type = "error", duration = 10)
				message("Error: ", e$message)
			})
		})

		# Save status
		output$save_status <- renderText({
			if (!is.null(cases_r()) && !is.null(controls_r())) {
				"Cohort ready to save"
			} else {
				"Create a matched cohort first"
			}
		})

		return(list(
			cases_r = cases_r,
			controls_r = controls_r,
			prepared_study_data_r = prepared_study_data_r
		))

	})
}
