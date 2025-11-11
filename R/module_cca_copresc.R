#' CCA Co-prescription Analysis Module
#'
#' This module calculates pairwise odds ratios for co-prescription patterns
#' separately in cases and controls, then compares them.
#'
#' @details
#' The module uses utility functions from:
#' - utils_cca_copresc.R: OR calculations
#' - utils_cca_copresc_plots.R: Visualizations

module_cca_copresc_ui <- function(id) {
	ns <- NS(id)

	tagList(
		card(
			card_header("Co-prescription Analysis Settings"),
			layout_columns(
				col_widths = c(4, 4, 4),
				numericInput(
					ns("min_prevalence"),
					"Minimum drug prevalence (%):",
					value = 2,
					min = 0.1,
					max = 10,
					step = 0.1
				),
				numericInput(
					ns("min_coprescription"),
					"Minimum co-prescription (%):",
					value = 1,
					min = 0.1,
					max = 10,
					step = 0.1
				),
				actionButton(
					ns("calculate"),
					"Calculate Co-prescription ORs",
					class = "btn-primary",
					style = "margin-top: 25px;"
				)
			)
		),
		navset_card_tab(
			id = ns("results_tabs"),
			nav_panel(
				title = "Table",
				card_body(
					downloadButton(ns("download_table"), "Download Table"),
					reactableOutput(ns("copresc_table"))
				)
			),
			nav_panel(
				title = "Heatmap",
				card_body(
					vegawidgetOutput(ns("copresc_heatmap"), height = "700px")
				)
			)
		)
	)
}

module_cca_copresc_server <- function(id, prescriptions_r, patient_data_r) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Reactive to store calculated results
		copresc_results <- reactiveVal(NULL)

		# Calculate co-prescription ORs when button clicked
		observeEvent(input$calculate, {
			req(prescriptions_r(), patient_data_r())

			showNotification("Calculating co-prescription ORs...",
											 id = "copresc_calc",
											 duration = NULL,
											 type = "message")

			tryCatch({
				results <- calculate_coprescription_ors(
					prescriptions = prescriptions_r(),
					min_prevalence = input$min_prevalence / 100,
					min_coprescription = input$min_coprescription / 100
				)

				copresc_results(results)

				removeNotification("copresc_calc")

				n_sig <- sum(results$significant == TRUE, na.rm = TRUE)
				showNotification(
					sprintf("Analysis complete: %d drug pairs found, %d with significant differences",
									nrow(results), n_sig),
					type = "message",
					duration = 5
				)
			}, error = function(e) {
				removeNotification("copresc_calc")
				showNotification(
					paste("Error calculating ORs:", e$message),
					type = "error",
					duration = 10
				)
			})
		})

		# Render table
		output$copresc_table <- renderReactable({
			req(copresc_results())

			data <- copresc_results()

			# Check for missing values and handle them
			data <- data[!is.na(case_or) & !is.na(control_or)]

			if (nrow(data) == 0) {
				return(NULL)
			}

			# Format OR columns
			data[, case_or_display := sprintf("%.2f (%.2f-%.2f)",
																				case_or, case_ci_lower, case_ci_upper)]
			data[, control_or_display := sprintf("%.2f (%.2f-%.2f)",
																					 control_or, control_ci_lower, control_ci_upper)]

			# Create display table with abs_or_diff for sorting
			display_data <- data[, .(
				drug1,
				drug2,
				case_or = case_or_display,
				control_or = control_or_display,
				or_diff,
				abs_or_diff = base::abs(or_diff),
				significant
			)]

			reactable(
				display_data,
				columns = list(
					drug1 = colDef(name = "Drug 1", minWidth = 150),
					drug2 = colDef(name = "Drug 2", minWidth = 150),
					case_or = colDef(name = "Case OR (95% CI)", minWidth = 150),
					control_or = colDef(name = "Control OR (95% CI)", minWidth = 150),
					or_diff = colDef(
						name = "OR Difference",
						format = colFormat(digits = 2),
						style = function(value) {
							if (is.na(value)) return(NULL)
							if (!is.finite(value)) return(NULL)
							if (base::abs(value) > 1) {
								color <- if (value > 0) "#e74c3c" else "#3498db"
								list(fontWeight = "bold", color = color)
							}
						}
					),
					abs_or_diff = colDef(show = FALSE),  # Hidden column for sorting
					significant = colDef(
						name = "Significant*",
						cell = function(value) {
							if (is.na(value)) return("")
							if (value) "✓" else ""
						},
						align = "center",
						width = 100
					)
				),
				defaultPageSize = 20,
				searchable = TRUE,
				filterable = TRUE,
				defaultSorted = list(abs_or_diff = "desc"),
				theme = reactableTheme(
					borderColor = "#dfe2e5",
					stripedColor = "#f6f8fa"
				)
			)
		})

		# Render heatmap using vegawidget
		output$copresc_heatmap <- renderVegawidget({
			req(copresc_results())

			data <- copresc_results()

			# Filter out rows with missing values
			data <- data[!is.na(or_diff) & is.finite(or_diff)]

			if (nrow(data) == 0) {
				return(NULL)
			}

			create_copresc_heatmap(data)
		})

		# Download handler
		output$download_table <- downloadHandler(
			filename = function() {
				paste0("copresc_analysis_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(copresc_results())
				data <- copresc_results()
				fwrite(data, file)
			}
		)
	})
}