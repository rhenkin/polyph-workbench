module_psm_ui_old <- function(id) {
	ns <- NS(id)
	nav_panel("Propensity matching",
						fluidRow(
							column(6,
										 textInput(ns("study_name"), "Study name:", value = "my_study"),
										 actionButton(ns("find_eligible"), "Find eligible controls", class = "btn-primary"),
										 br(), br(),
										 textOutput(ns("number_of_eligible_controls"))
							),
							column(6,
										 selectizeInput(ns("exact_match_ltcs"), "Select LTCs for exact matching:",
										 							 choices = NULL, multiple = TRUE),
										 numericInput(ns("match_ratio"), "Matching ratio:",
										 						 value = 1, min = 1, max = 5, step = 1),
										 actionButton(ns("run_matching"), "Run matching", class = "btn-success"),
										 br(), br(),
										 actionButton(ns("save_matching_results"), "Save results", class = "btn-info")
							)
						),
						br(),
						card(
							card_header("Summary of PSM Results"),
							card_body(
								verbatimTextOutput(ns("matching_results_summary"))
							)
						)
	)
}

module_psm_server_old <- function(id, patient_data, outcome_prescriptions, ltc_data,
														 gold_patient, gold_ltc, gold_cp, gold_acute_presc, study_dir = "studies") {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Reactive values to store data
		eligible_controls_r <- reactiveVal(data.table())
		current_cases_r <- reactiveVal(data.table())
		matching_results_r <- reactiveVal(NULL)

		# Find eligible controls when button is clicked
		observeEvent(input$find_eligible, {
			req(input$study_name)
			browser()

			# Show progress notification
			id <- showNotification("Finding eligible controls...", duration = NULL, closeButton = FALSE)
			on.exit(removeNotification(id), add = TRUE)

			validate(need(nrow(outcome_prescriptions()) > 0, "No outcome prescriptions available"))

			# Prepare cases data from the filtered outcome_prescriptions
			cases_df <- unique(outcome_prescriptions()[, .(patid, eventdate)])
			cases_df[, study_name := input$study_name]

			# Store current cases for later use
			current_cases_r(cases_df)

			message("Finding eligible controls for study ", input$study_name, " at ", Sys.time())

			# Compute eligible controls using the FULL gold data.tables
			tryCatch({
				controls <- get_eligible_controls_extended(
					cases_df = cases_df,
					gold_patient = gold_patient,
					gold_ltc = gold_ltc,
					gold_cp = gold_cp,
					random_seed = 4
				)

				# Store eligible controls
				eligible_controls_r(controls)

				message("Finished finding eligible controls at ", Sys.time())
				showNotification("Eligible controls found successfully!", type = "message", duration = 3)

			}, error = function(e) {
				showNotification(paste("Error finding eligible controls:", e$message), type = "error", duration = 5)
			})
		})

		# Display number of eligible controls
		output$number_of_eligible_controls <- renderText({
			controls <- eligible_controls_r()
			if (nrow(controls) > 0) {
				n_unique_controls <- uniqueN(controls$patid)
				paste("Eligible controls found:", n_unique_controls)
			} else {
				"No eligible controls found yet - click 'Find eligible controls' first"
			}
		})

		# Prepare combined cases and controls data for matching
		cases_controls_r <- reactive({
			validate(need(nrow(outcome_prescriptions()) > 0, "No outcome prescriptions available"))
			validate(need(nrow(eligible_controls_r()) > 0, "No eligible controls found"))

			# Merge cases data
			cases_base <- unique(outcome_prescriptions()[, .(patid, outcome_age)])

			# Add demographics (left join to keep all cases)
			cases_with_demog <- merge(cases_base,
																patient_data()[, .(patid, dob, sex, imd_quintile, eth_group)],
																by = "patid", all.x = TRUE)

			# Add LTC data (left join to keep all cases, even those without LTCs)
			cases_ltcs <- ltc_data()[patid %in% cases_base$patid, .(patid, eventdate, term)]
			cases_ltcs$ltc <- 1
			cases_ltcs <-
				dcast(cases_ltcs, formula = patid ~ term, fill = 0, value.var="ltc")

			cases_merged <- merge(cases_with_demog, cases_ltcs, by = "patid")
			cases_merged[, treatment := 1]

			# Fill any missing LTC columns with 0
			ltc_cols <- setdiff(colnames(cases_merged), c("patid", "outcome_age", "sex", "imd_quintile", "eth_group", "treatment"))
			cases_merged[, (ltc_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = ltc_cols]

			# Add this debug line to check:
			message(sprintf("Cases going into matching: %d", nrow(cases_merged)))

			# Prepare controls data from eligible_controls_r()
			eligible_controls_dt <- copy(eligible_controls_r())

			# Reshape controls LTC data to wide format
			eligible_controls_dt$ltc <- 1
			controls_ltc_wide <- dcast(eligible_controls_dt, patid~term, value.var = "ltc", fill = 0)
			demog <- unique(eligible_controls_dt[,.(patid, index_date, gender, imd_quintile, eth_group, dob)])
			controls_ltc_wide <- merge(controls_ltc_wide, demog)

			# Calculate outcome age for controls
			controls_ltc_wide[, outcome_age := as.numeric(as.IDate(index_date) - as.IDate(dob))]
			controls_ltc_wide[, `:=`(dob = NULL, index_date = NULL)]

			# Rename gender to sex for consistency
			setnames(controls_ltc_wide, "gender", "sex")
			controls_ltc_wide[, treatment := 0]

			# Combine cases and controls
			cases_controls <- rbindlist(list(controls_ltc_wide, cases_merged), fill = TRUE)

			# Clean column names
			colnames(cases_controls) <- make.names(colnames(cases_controls))

			cases_controls <- cases_controls[!is.na(patid) & !is.na(treatment) & !is.na(sex) & !is.na(eth_group) & !is.na(imd_quintile) & !is.na(outcome_age)]

			# Fill NA values with 0 for LTC columns
			ltc_cols <- setdiff(colnames(cases_controls), c("patid", "treatment", "sex", "eth_group", "imd_quintile", "outcome_age"))
			cases_controls[, (ltc_cols) := lapply(.SD, function(x) fifelse(is.na(x), 0, x)), .SDcols = ltc_cols]

			# Convert outcome age to years
			cases_controls[, outcome_age := round(outcome_age / 365.25, digits = 0)]

			cases_controls
		})

		# Update LTC choices for exact matching
		observe({
			req(cases_controls_r())
			cases_controls <- cases_controls_r()
			ltc_cols <- setdiff(colnames(cases_controls), c("patid", "treatment", "sex", "eth_group", "imd_quintile", "outcome_age"))
			updateSelectizeInput(session, "exact_match_ltcs", choices = ltc_cols)
		})

		# Run matching when button is clicked
		observeEvent(input$run_matching, {
			req(cases_controls_r())

			# Show progress notification
			id <- showNotification("Running propensity score matching...", duration = NULL, closeButton = FALSE)
			on.exit(removeNotification(id), add = TRUE)

			cases_controls <- cases_controls_r()

			# Get LTC columns with minimum prevalence
			ltc_cols <- setdiff(colnames(cases_controls), c("patid", "treatment", "sex", "eth_group", "imd_quintile", "outcome_age"))
			min_prev_ltcs <- ltc_cols[colSums(cases_controls[, ..ltc_cols]) > 0.01 * nrow(cases_controls)]

			# Build formula for matching
			#rhs_cols <- c("sex", "eth_group", "imd_quintile", "outcome_age", min_prev_ltcs)
			rhs_cols <- c("outcome_age", min_prev_ltcs)
			if (length(input$exact_match_ltcs) > 0) {
				rhs_cols <- c(rhs_cols, input$exact_match_ltcs)
			}
			rhs <- paste(unique(rhs_cols), collapse = " + ")
			matching_formula <- as.formula(paste("treatment ~ ", rhs))

			tryCatch({
				# Ensure MatchIt is loaded
				if (!requireNamespace("MatchIt", quietly = TRUE)) {
					stop("MatchIt package is required for propensity score matching")
				}

				# Run matching
				m.out <- MatchIt::matchit(
					formula = matching_formula,
					data = cases_controls,
					method = "nearest",
					ratio = input$match_ratio,
					replace = FALSE,
					exact = if (length(input$exact_match_ltcs) > 0) input$exact_match_ltcs else NULL
				)

				# Store results
				matching_results_r(m.out)

				showNotification("Matching completed successfully!", type = "message", duration = 3)

			}, error = function(e) {
				showNotification(paste("Error during matching:", e$message), type = "error", duration = 5)
			})
		})

		# Save matching results - FIXED VERSION
		observeEvent(input$save_matching_results, {
			req(matching_results_r())
			req(input$study_name)

			# Show progress notification
			id <- showNotification("Saving matching results...", duration = NULL, closeButton = FALSE)
			on.exit(removeNotification(id), add = TRUE)

			tryCatch({
				# CRITICAL FIX: Create matched patids with proper eventdate mapping
				matched_data <- MatchIt::match.data(matching_results_r(), data = cases_controls_r())
				matched_patids <- matched_data[, .(patid, treatment)]

				# DEBUG: Check for duplicates in matched_data
				message(sprintf("Matched data rows: %d, Unique patients: %d",
												nrow(matched_data), uniqueN(matched_data$patid)))
				message(sprintf("Cases in matched_data: %d, Controls: %d",
												sum(matched_data$treatment == 1), sum(matched_data$treatment == 0)))

				# FIXED: Handle eventdates separately for cases and controls to avoid duplication
				# For cases: get their actual outcome eventdate (ensure unique per patient)
				cases_in_matched <- matched_patids[treatment == 1, patid]
				cases_with_dates <- current_cases_r()[patid %in% cases_in_matched] |> unique(by = "patid")

				# For controls: use their assigned index_date (ensure unique per patient)
				controls_in_matched <- matched_patids[treatment == 0, patid]
				controls_with_dates <- eligible_controls_r()[patid %in% controls_in_matched, .(patid, eventdate = index_date)] |> unique(by = "patid")

				# FIXED: Ensure no duplicates in matched_patids first
				matched_patids_unique <- unique(matched_patids, by = "patid")

				# FIXED: Combine and merge to avoid duplication
				cases_final <- merge(matched_patids_unique[treatment == 1], cases_with_dates, by = "patid", all.x = TRUE)
				controls_final <- merge(matched_patids_unique[treatment == 0], controls_with_dates, by = "patid", all.x = TRUE)
				controls_final$study_name <- input$study_name

				# Combine the final matched patids
				matched_patids <- rbindlist(list(cases_final, controls_final))

				# DEBUG: Check final matched_patids
				message(sprintf("Final matched_patids - Cases: %d, Controls: %d, Total rows: %d",
												sum(matched_patids$treatment == 1), sum(matched_patids$treatment == 0), nrow(matched_patids)))

				# Verify we have eventdates for all matched patients
				if (any(is.na(matched_patids$eventdate))) {
					warning("Some matched patients missing eventdates")
				}

				save_matched_datasets(
					study_name = input$study_name,
					matched_patids = matched_patids,
					gold_patient = gold_patient,
					gold_cp = gold_cp,
					gold_ltc = gold_ltc,
					cases_df = current_cases_r(),
					outcome_prescriptions = outcome_prescriptions(),
					study_dir = study_dir
				)

				showNotification(paste("Study", input$study_name, "saved successfully!"), type = "message", duration = 3)

			}, error = function(e) {
				showNotification(paste("Error saving results:", e$message), type = "error", duration = 5)
			})
		})

		# Display matching results summary
		output$matching_results_summary <- renderPrint({
			req(matching_results_r())
			summary(matching_results_r(), un = FALSE)
		})

	})
}