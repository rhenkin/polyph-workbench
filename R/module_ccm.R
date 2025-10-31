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
		match_summary_r <- reactiveVal(NULL)
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
				# Step 1: Create cases table
				progress$set(detail = "Building cases table...", value = 0.2)

				outcomes <- outcome_prescriptions()[, .(patid, outcome_date = eventdate)] |> unique()

				mrp_dataset <- master_risk_pool_dataset()
				case_patids <- unique(outcomes$patid)

				cases_raw <- mrp_dataset %>%
					dplyr::filter(patid %in% case_patids) %>%
					dplyr::select(patid, prescription_date, substance, stratum_alt, stratum_mm_time, n_ltcs,
												concurrent_cps, age_at_rx, sex, imd_quintile, first_presc_bin, time_since_first_presc) %>%
					dplyr::collect() %>%
					as.data.table()

				cases_raw[ , p_date := as.IDate(prescription_date)]
				cases_raw[, prescription_date := NULL]
				colnames(cases_raw)[length(cases_raw)] <- "prescription_date"
				new_cases <- copy(cases_raw)
				new_cases <- new_cases[outcomes,on = "patid",nomatch =0]

				new_cases <- new_cases[
					prescription_date < outcome_date &
						outcome_date - prescription_date <= input$pred_window
				]
				# Keep only most recent prescriptions
				new_cases <- new_cases[,.SD[which(prescription_date==max(prescription_date))],patid]

				progress$set(detail = "Calculating valid LTCs and CPs...", value = 0.3)

				# Calculate valid LTCs (before index_date, need ≥2)
				ltc <- ltc_data()
				valid_ltcs <- new_cases[ltc,
																				.(patid, term, eventdate = i.eventdate, start_date = x.prescription_date),
																				on = .(patid, prescription_date > eventdate),
														mult="last",
																				nomatch = 0]

				setkey(valid_ltcs, patid, eventdate)
				valid_ltc_patids <- valid_ltcs[, .N, patid][N >= 2, patid]

				valid_ltcs <- valid_ltcs[ltc_chapters, on = .(term = ltc), nomatch = 0]
				setorder(valid_ltcs, patid, eventdate)
				valid_ltcs[, ltc_index := 1:.N, patid]

				cases <- merge(new_cases, valid_ltcs[, .(
					mm_date = .SD[ltc_index==2, eventdate]
				)
				, patid], by = "patid")


				cases <- merge(cases, patient_data()[, .(patid, dob, eth_group)],
											 by = "patid")

				cases[, ':='(
					index_age = as.numeric(prescription_date - dob) / 365.25,
					time_at_risk = round(as.numeric(outcome_date - mm_date)/365.25, digits = 0)
				)]



				cases[, mm_duration_bin := cut(time_at_risk,
																								 breaks = seq(0, 121,
																								 						 2),
																								 include.lowest = TRUE,
																								 right = FALSE)]

				# Add stratification variables
				cases[, ':='(
					age_bin = cut(index_age, breaks = seq(0, 120, 5), include.lowest = TRUE, right = FALSE),
					year = year(prescription_date)
				)]
				#cases[, strata := paste(sex, age_bin, n_ltcs, year, sep = "_")]
				#cases[, strata := paste(sex, mm_duration_bin, age_bin, sep = "_")]
				cases[, strata := stratum_first_presc_bin]
				# Store cases
				cases_r(cases)

				message(sprintf("Cases created: %d patients", nrow(cases)))

				# Step 2: Filter master risk pool to eligible controls
				progress$set(detail = "Filtering master risk pool...", value = 0.5)

				case_patids <- unique(cases$patid)
				strata_needs <- unique(cases$strata)

				# Use Arrow to filter before collecting into memory
				# Only select columns we need for sampling
				eligible_pool <- mrp_dataset %>%
					dplyr::filter(stratum_first_presc_bin %in% strata_needs) %>%
					dplyr::select(patid, prescription_date, substance, sex, age_at_rx, n_ltcs,
								 imd_quintile, age_bin, stratum_first_presc_bin, year, first_presc_bin, time_since_first_presc) %>%
					dplyr::collect() %>%
					as.data.table()

				# Convert Date to IDate for consistency
				eligible_pool[, prescription_date := as.IDate(prescription_date)]

				# Filter out case patients in memory (faster than anti-join in Arrow)
				eligible_pool <- eligible_pool[!patid %in% case_patids]

				message(sprintf("Eligible prescriptions: %d", nrow(eligible_pool)))

				# Step 3: Get most recent prescription per patient per stratum
				progress$set(detail = "Getting unique patients per stratum...", value = 0.6)

				setorder(eligible_pool, patid, prescription_date)
				eligible_unique <- eligible_pool[
					,
					last(.SD),
					by = .(strata = stratum_first_presc_bin, patid)
				]

				message(sprintf("Unique patient-strata combinations: %d", nrow(eligible_unique)))

				# Step 4: Calculate controls needed per stratum
				strata_needs_dt <- cases[, .(n_cases = uniqueN(patid), n_controls_needed = uniqueN(patid) * input$match_ratio), by = strata]
				eligible_unique[strata_needs_dt, n_needed := i.n_controls_needed, on = "strata"]

				# Step 5: Sample within each stratum
				progress$set(detail = "Sampling controls...", value = 0.7)

				all_controls <- eligible_unique[
					!is.na(n_needed),
					{
						n_sample <- min(.N, n_needed[1])
						if(n_sample > 0) .SD[sample(.N, n_sample, replace = FALSE)] else .SD[0]
					},
					by = strata
				]

				all_controls[, n_needed := NULL]
				all_controls[, control_index_date := prescription_date]
				# After the current sampling code, add:
				setorder(all_controls, patid, prescription_date)

				# This removes duplicates
				all_controls <- all_controls[, .SD[sample(.N,1)], by = patid]  # Keep random occurrence per patient

				# Store controls
				controls_r(all_controls)

				message(sprintf("Controls sampled: %d (unique patients: %d)",
												nrow(all_controls), uniqueN(all_controls$patid)))

				# Step 6: Quality checks
				progress$set(detail = "Calculating quality metrics...", value = 0.9)

				match_summary <- merge(
					cases[, .(n_cases = .N), by = strata],
					all_controls[, .(n_controls = .N), by = strata],
					by = "strata",
					all.x = TRUE
				)

				match_summary[is.na(n_controls), n_controls := 0]
				match_summary[, ratio := n_controls / n_cases]

				match_summary_r(match_summary)

				output$show_results <- reactive({ TRUE })
				outputOptions(output, "show_results", suspendWhenHidden = FALSE)

				progress$set(value = 1)
				showNotification("Matched cohort created successfully!", type = "message", duration = 3)

			}, error = function(e) {
				showNotification(paste("Error creating cohort:", e$message), type = "error", duration = 10)
				message("Error: ", e$message)
				print(traceback())
			})
		})


		observe({
			req(cases_r(), controls_r())

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

				spec <- list(
					`$schema` = vega_schema(),
					title = "Age at prescription",
					height = 160,
					# width = width,
					data = list(values = df),
					encoding = list(
						y = list(field = "group", type = "nominal", title = NULL)
					),
					layer = list(
						list(
							mark = list(type = "rule"),
							encoding = list(
								x = list(
									field = "lower",
									type = "quantitative",
									scale = list(zero = FALSE),
									title = NULL
								),
								x2 = list(field = "upper")
							)
						),
						list(
							mark = list(type = "bar", size = 36),
							encoding = list(
								x = list(field = "q1", type = "quantitative"),
								x2 = list(field = "q3"),
								color = list(field = "group", type = "nominal", legend = NULL,
														 scale = list(range = list("#e74c3c", "#2C3E50")))
							)
						),
						list(
							mark = list(type = "tick", color = "white", size = 36),
							encoding = list(
								x = list(field = "median", type = "quantitative")
							)
						),
						list(
							mark = list(type = "text", dy = -25, align = "center"),
							encoding = list(
								x = list(field = "median", type = "quantitative"),
								text = list(field = "median"),
								color = list(value = "black")
							)
						)
					),
					config = list(
						view = list(stroke = NULL),
						font = "Lato, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, \"Helvetica Neue\", Arial, sans-serif, \"Apple Color Emoji\", \"Segoe UI Emoji\", \"Segoe UI Symbol\""
					)
				)
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

				spec <- list(
					`$schema` = vega_schema(),
					title = "Time since first prescription",
					height = 160,
					# width = width,
					data = list(values = df),
					encoding = list(
						y = list(field = "group", type = "nominal", title = NULL)
					),
					layer = list(
						list(
							mark = list(type = "rule"),
							encoding = list(
								x = list(
									field = "lower",
									type = "quantitative",
									scale = list(zero = FALSE),
									title = NULL
								),
								x2 = list(field = "upper")
							)
						),
						list(
							mark = list(type = "bar", size = 36),
							encoding = list(
								x = list(field = "q1", type = "quantitative"),
								x2 = list(field = "q3"),
								color = list(field = "group", type = "nominal", legend = NULL,
														 scale = list(range = list("#e74c3c", "#2C3E50")))
							)
						),
						list(
							mark = list(type = "tick", color = "white", size = 36),
							encoding = list(
								x = list(field = "median", type = "quantitative")
							)
						),
						list(
							mark = list(type = "text", dy = -25, align = "center"),
							encoding = list(
								x = list(field = "median", type = "quantitative"),
								text = list(field = "median"),
								color = list(value = "black")
							)
						)
					),
					config = list(
						view = list(stroke = NULL),
						font = "Lato, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, \"Helvetica Neue\", Arial, sans-serif, \"Apple Color Emoji\", \"Segoe UI Emoji\", \"Segoe UI Symbol\""
					)
				)
				spec |> as_vegaspec()

			})


		})


		# Display cohort summary
		output$cohort_summary <- renderPrint({
			req(cases_r(), controls_r())

			cases <- cases_r()
			controls <- controls_r()

			cat("=== MATCHED COHORT SUMMARY ===\n\n")
			cat(sprintf("Cases: %d\n", uniqueN(cases$patid)))
			cat(sprintf("Controls: %d\n", uniqueN(controls$patid)))
			cat(sprintf("Overall ratio: %.2f:1\n\n", uniqueN(controls$patid) / uniqueN(cases$patid)))

			cat("Age distribution:\n")
			cat(sprintf("  Cases - mean: %.1f years\n", mean(cases$index_age)))
			cat(sprintf("  Controls - mean: %.1f years\n\n", mean(controls$age_at_rx)))

			cat("Sex distribution:\n")
			cat("  Cases:\n")
			print(prop.table(table(cases$sex)))
			cat("  Controls:\n")
			print(prop.table(table(controls$sex)))

			cat("\nDisease burden:\n")
			cat(sprintf("  Cases - mean n_ltc: %.1f\n", mean(cases$n_ltc)))
			cat(sprintf("  Controls - mean n_ltcs: %.1f\n", mean(controls$n_ltcs)))
		})

		# Display quality summary
		output$quality_summary <- renderPrint({
			req(match_summary_r())

			match_summary <- match_summary_r()

			cat("=== MATCHING QUALITY ===\n\n")
			cat("Control:Case Ratio by Stratum:\n")
			print(summary(match_summary$ratio))

			cat("\n\nDistribution of ratios:\n")
			print(table(cut(match_summary$ratio,
											breaks = c(0, 2, 3, 4, 5, Inf),
											labels = c("<2", "2-3", "3-4", "4-5", ">5"))))

			insufficient_strata <- match_summary[ratio < input$match_ratio]
			if(nrow(insufficient_strata) > 0) {
				cat(sprintf("\n%d strata have <%d controls per case\n",
										nrow(insufficient_strata), input$match_ratio))
				cat(sprintf("Total cases affected: %d\n", sum(insufficient_strata$n_cases)))
				cat(sprintf("Proportion of cases: %.1f%%\n",
										sum(insufficient_strata$n_cases) / sum(match_summary$n_cases) * 100))
			}

			duplicate_rate <- 1 - (uniqueN(controls_r()$patid) / nrow(controls_r()))
			cat(sprintf("\nDuplicate control rate: %.3f%%\n", duplicate_rate * 100))
		})

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
		})

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
