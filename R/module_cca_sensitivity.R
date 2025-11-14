#' CCA Sensitivity Analysis Module
#'
#' This module performs post-hoc sensitivity analyses on matched case-control data
#' by stratifying or excluding patients based on LTC presence, without re-matching.
#'
#' @details
#' Provides three types of sensitivity analysis:
#' 1. Stratified analysis - Calculate ORs within subgroups (with/without specific LTCs)
#' 2. Exclusion sensitivity - Recalculate ORs after excluding patients with specific LTCs
#' 3. Comparison summary - Visual comparison of robustness across analyses

module_cca_sensitivity_ui <- function(id) {
	ns <- NS(id)

	accordion(
		open = FALSE,
		accordion_panel(
			title = "Sensitivity Analysis",
			value = "sensitivity_panel",

			p("Test the robustness of findings by analyzing subgroups or excluding patients with specific conditions. All analyses use the same matched cohort - no re-matching required."),

			navset_card_tab(
				# Tab 1: Stratified Analysis
				nav_panel(
					title = "Stratified Analysis",
					icon = icon("layer-group"),

					card_body(
						p("Compare ORs between patients WITH and WITHOUT specific LTCs. Helps identify effect modification or residual confounding."),

						virtualSelectInput(
							ns("stratify_ltcs"),
							label = "Stratify by presence of (max 2 conditions):",
							choices = NULL,
							multiple = TRUE,
							maxValues = 2,
							dropboxWrapper = "body"
						),

						numericInput(
							ns("stratify_min_prev"),
							label = "Minimum drug prevalence (%):",
							value = 2,
							min = 0.5,
							max = 10,
							step = 0.5
						),

						actionButton(ns("run_stratified"), "Run Stratified Analysis",
												 class = "btn-primary", icon = icon("play")),

						hr(),

						h4("Results"),
						p("ORs are calculated separately for patients WITH and WITHOUT the selected condition(s)."),

						conditionalPanel(
							condition = "output.stratified_results_available",
							ns = ns,
							downloadButton(ns("download_stratified"), "Download Results", class = "btn-sm btn-secondary"),
							br(), br()
						),

						reactableOutput(ns("stratified_table"))
					)
				),

				# Tab 2: Exclusion Sensitivity
				nav_panel(
					title = "Exclusion Sensitivity",
					icon = icon("filter"),

					card_body(
						p("Remove patients with specific LTCs and recalculate ORs. If results change substantially, those conditions may be important confounders or effect modifiers."),

						virtualSelectInput(
							ns("exclude_ltcs"),
							label = "Exclude patients with:",
							choices = NULL,
							multiple = TRUE,
							maxValues = 3,
							dropboxWrapper = "body"
						),

						numericInput(
							ns("exclude_min_prev"),
							label = "Minimum drug prevalence (%):",
							value = 2,
							min = 0.5,
							max = 10,
							step = 0.5
						),

						actionButton(ns("run_exclusion"), "Run Exclusion Analysis",
												 class = "btn-primary", icon = icon("play")),

						hr(),

						h4("Results"),
						p("Compares main analysis ORs to ORs after excluding patients with selected conditions."),

						conditionalPanel(
							condition = "output.exclusion_results_available",
							ns = ns,
							downloadButton(ns("download_exclusion"), "Download Results", class = "btn-sm btn-secondary"),
							br(), br()
						),

						reactableOutput(ns("exclusion_table"))
					)
				),

				# Tab 3: Comparison Summary
				nav_panel(
					title = "Robustness Summary",
					icon = icon("chart-line"),

					card_body(
						p("Visual comparison showing how OR estimates change across different analyses."),

						selectInput(
							ns("comparison_drug"),
							label = "Select drug to visualize:",
							choices = NULL
						),

						vegawidgetOutput(ns("forest_plot"), height = "400px"),

						hr(),

						h4("Interpretation Guide"),
						tags$ul(
							tags$li(tags$strong("Robust finding:"), " ORs are similar across all analyses (within 20%)"),
							tags$li(tags$strong("Effect modification:"), " ORs differ substantially between subgroups (e.g., stronger in patients with specific LTC)"),
							tags$li(tags$strong("Potential confounding:"), " OR changes substantially when excluding a condition (suggests that condition was confounding the association)")
						)
					)
				)
			)
		)
	)
}

module_cca_sensitivity_server <- function(id, patient_data_r, prescriptions_aggregated_r, ltcs_r) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Reactive values to store results
		main_or_results <- reactiveVal(NULL)
		stratified_results <- reactiveVal(NULL)
		exclusion_results <- reactiveVal(NULL)

		# Get list of common LTCs for dropdowns
		observe({
			req(ltcs_r(), patient_data_r())

			ltcs <- ltcs_r()
			patient_data <- patient_data_r()

			# Calculate LTC prevalence
			ltc_prevalence <- ltcs[, .(
				n_patients = uniqueN(patid),
				prevalence = uniqueN(patid) / uniqueN(patient_data$patid)
			), by = term]

			# Keep LTCs with >2% prevalence
			common_ltcs <- ltc_prevalence[prevalence > 0.02][order(-prevalence), term]

			# Update both dropdowns
			updateVirtualSelect("stratify_ltcs", choices = common_ltcs)
			updateVirtualSelect("exclude_ltcs", choices = common_ltcs)
		})

		# Calculate main OR results once when data loads
		observe({
			req(prescriptions_aggregated_r(), patient_data_r())

			presc <- prescriptions_aggregated_r()
			patient_data <- patient_data_r()

			# Merge with case/control status
			presc_with_group <- merge(
				presc,
				patient_data[, .(patid, treatment)],
				by = "patid"
			)
			presc_with_group[, group := ifelse(treatment == 1, "case", "control")]

			# Calculate ORs
			main_results <- calculate_prescription_ors(
				presc_with_group,
				min_prevalence = 0.02
			)

			main_or_results(main_results)

			# Update drug dropdown for forest plot
			if (!is.null(main_results) && nrow(main_results) > 0) {
				top_drugs <- main_results[order(-OR)][1:min(20, nrow(main_results)), substance]
				updateSelectInput(session, "comparison_drug", choices = top_drugs)
			}
		})

		# Stratified Analysis
		observeEvent(input$run_stratified, {
			req(input$stratify_ltcs, prescriptions_aggregated_r(), ltcs_r(), patient_data_r())

			showNotification("Running stratified analysis...", type = "message", duration = NULL, id = "stratify_notification")

			tryCatch({
				results <- run_stratified_analysis(
					prescriptions = prescriptions_aggregated_r(),
					ltcs = ltcs_r(),
					patient_data = patient_data_r(),
					stratify_by_ltcs = input$stratify_ltcs,
					min_prevalence = input$stratify_min_prev / 100
				)

				stratified_results(results)
				output$stratified_results_available <- reactive({ TRUE })
				outputOptions(output, "stratified_results_available", suspendWhenHidden = FALSE)

				removeNotification("stratify_notification")
				showNotification("Stratified analysis complete!", type = "message", duration = 3)

			}, error = function(e) {
				removeNotification("stratify_notification")
				showNotification(paste("Error:", e$message), type = "error", duration = 10)
			})
		})

		# Render stratified results table
		output$stratified_table <- renderReactable({
			req(stratified_results())

			results <- stratified_results()

			reactable(
				results,
				columns = list(
					substance = colDef(name = "Drug", minWidth = 150),
					stratum = colDef(name = "Patient Group", minWidth = 120),
					OR = colDef(name = "OR", format = colFormat(digits = 2)),
					CI_lower = colDef(name = "95% CI Lower", format = colFormat(digits = 2)),
					CI_upper = colDef(name = "95% CI Upper", format = colFormat(digits = 2)),
					n_cases = colDef(name = "Cases", format = colFormat(separators = TRUE)),
					n_controls = colDef(name = "Controls", format = colFormat(separators = TRUE))
					#interpretation = colDef(name = "Interpretation", minWidth = 150)
				),
				groupBy = "substance",
				searchable = TRUE,
				defaultPageSize = 20,
				striped = TRUE,
				highlight = TRUE,
				compact = TRUE
			)
		})

		# Exclusion Analysis
		observeEvent(input$run_exclusion, {
			req(input$exclude_ltcs, prescriptions_aggregated_r(), ltcs_r(), patient_data_r(), main_or_results())

			showNotification("Running exclusion analysis...", type = "message", duration = NULL, id = "exclusion_notification")

			tryCatch({
				results <- run_exclusion_analysis(
					prescriptions = prescriptions_aggregated_r(),
					ltcs = ltcs_r(),
					patient_data = patient_data_r(),
					exclude_ltcs = input$exclude_ltcs,
					main_or_results = main_or_results(),
					min_prevalence = input$exclude_min_prev / 100
				)

				exclusion_results(results)
				output$exclusion_results_available <- reactive({ TRUE })
				outputOptions(output, "exclusion_results_available", suspendWhenHidden = FALSE)

				removeNotification("exclusion_notification")
				showNotification("Exclusion analysis complete!", type = "message", duration = 3)

			}, error = function(e) {
				removeNotification("exclusion_notification")
				showNotification(paste("Error:", e$message), type = "error", duration = 10)
			})
		})

		# Render exclusion results table
		output$exclusion_table <- renderReactable({
			req(exclusion_results())

			results <- exclusion_results()

			reactable(
				results,
				columns = list(
					substance = colDef(name = "Drug", minWidth = 150),
					main_OR = colDef(name = "Main OR", format = colFormat(digits = 2)),
					excluded_OR = colDef(name = "Excluded OR", format = colFormat(digits = 2)),
					pct_change = colDef(
						name = "% Change",
						format = colFormat(digits = 1, suffix = "%"),
						style = function(value) {
							if (is.na(value)) return(list(color = "#999"))
							if (abs(value) < 10) list(color = "#27ae60", fontWeight = "bold")
							else if (abs(value) < 25) list(color = "#f39c12", fontWeight = "bold")
							else list(color = "#e74c3c", fontWeight = "bold")
						}
					),
					main_CI = colDef(name = "Main 95% CI", minWidth = 120),
					excluded_CI = colDef(name = "Excluded 95% CI", minWidth = 120),
					n_excluded = colDef(name = "Patients Excluded", format = colFormat(separators = TRUE)),
					robustness = colDef(name = "Robustness", minWidth = 100)
				),
				searchable = TRUE,
				defaultPageSize = 20,
				striped = TRUE,
				highlight = TRUE,
				compact = TRUE,
				defaultSorted = "pct_change",
				defaultSortOrder = "desc"
			)
		})

		# Forest plot for comparison
		output$forest_plot <- renderVegawidget({
			req(input$comparison_drug, main_or_results())

			tryCatch({
				# Collect all results for the selected drug
				drug <- input$comparison_drug
				main_results <- main_or_results()

				# Check if drug exists in main results
				if (!drug %in% main_results$substance) {
					return(NULL)
				}

				# Start with main analysis data
				plot_data <- data.table(
					analysis = "Main Analysis",
					OR = main_results[substance == drug, OR],
					CI_lower = main_results[substance == drug, CI_lower],
					CI_upper = main_results[substance == drug, CI_upper]
				)

				# Add stratified results if available
				if (!is.null(stratified_results())) {
					strat <- stratified_results()[substance == drug]
					if (nrow(strat) > 0) {
						strat_data <- strat[, .(
							analysis = stratum,
							OR = OR,
							CI_lower = CI_lower,
							CI_upper = CI_upper
						)]
						plot_data <- rbindlist(list(plot_data, strat_data), use.names = TRUE)
					}
				}

				# Add exclusion results if available
				if (!is.null(exclusion_results())) {
					excl <- exclusion_results()[substance == drug]
					if (nrow(excl) > 0) {
						# For exclusion analysis, we need to extract CI from the strings
						# Format is "(lower-upper)"
						excl_data <- data.table(
							analysis = "Exclusion Analysis",
							OR = excl$excluded_OR,
							CI_lower = NA_real_,
							CI_upper = NA_real_
						)

						# Parse CI string if needed
						if ("excluded_CI" %in% names(excl)) {
							ci_parts <- gsub("[()]", "", excl$excluded_CI)
							ci_split <- strsplit(ci_parts, "-")[[1]]
							if (length(ci_split) == 2) {
								excl_data$CI_lower <- as.numeric(ci_split[1])
								excl_data$CI_upper <- as.numeric(ci_split[2])
							}
						} else if ("excluded_CI_lower" %in% names(excl)) {
							excl_data$CI_lower <- excl$excluded_CI_lower
							excl_data$CI_upper <- excl$excluded_CI_upper
						}

						plot_data <- rbindlist(list(plot_data, excl_data), use.names = TRUE)
					}
				}

				# Create forest plot using Vega-Lite
				if (nrow(plot_data) > 0) {
					create_forest_plot(plot_data, drug_name = drug)
				} else {
					NULL
				}

			}, error = function(e) {
				message("Error creating forest plot: ", e$message)
				NULL
			})
		})

		# Download handlers
		output$download_stratified <- downloadHandler(
			filename = function() {
				paste0("stratified_analysis_", Sys.Date(), ".csv")
			},
			content = function(file) {
				write.csv(stratified_results(), file, row.names = FALSE)
			}
		)

		output$download_exclusion <- downloadHandler(
			filename = function() {
				paste0("exclusion_analysis_", Sys.Date(), ".csv")
			},
			content = function(file) {
				write.csv(exclusion_results(), file, row.names = FALSE)
			}
		)
	})
}