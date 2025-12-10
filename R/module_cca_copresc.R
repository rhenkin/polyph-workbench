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

	accordion(
		open = FALSE,
		accordion_panel(
			title = "Results",
			value = "copresc_results_panel",
			icon = bs_icon("graph-up"),
			layout_columns(
				col_widths = c(1,1,1),
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
				checkboxInput(ns("or_min_filter"),
											"OR > 1 only?",value =FALSE)
			),
			navset_card_tab(
				id = ns("results_tabs"),
				nav_panel(
					title = "Forest Plot",
					card_body(
						flowLayout(
							selectInput(
								ns("forest_sort"),
								"Sort by:",
								choices = c(
									"Drug 1" = "drug1",
									"Drug 2" = "drug2",
									"Case OR" = "case_or",
									"Control OR" = "control_or",
									"OR Difference" = "or_diff"
								),
								selected = "or_diff"
							),
							checkboxInput(ns("forest_sort_asc"), "Lowest to highest?", value = FALSE)
						),
						vegawidgetOutput(ns("copresc_forest"), height = "800px"),
						layout_column_wrap(width="100px",fixed_width = TRUE,
							div(
								style = "padding-top: 25px;",
								actionButton(ns("prev_page"), "← Previous", class = "btn-secondary")
							),
							div(
								style = "padding-top: 32px; text-align: center;",
								textOutput(ns("page_info"))
							),
							div(
								style = "padding-top: 25px;",
								actionButton(ns("next_page"), "Next →", class = "btn-secondary")
							)
						)
					)
				),
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
						radioButtons(
							ns("heatmap_group"),
							"Show OR for:",
							choices = c("Cases" = "case", "Controls" = "control"),
							selected = "case",
							inline = TRUE
						),
						vegawidgetOutput(ns("copresc_heatmap"), height = "700px")
					)
				)
			)
		)
	)
}

module_cca_copresc_server <- function(id, prescriptions_r, patient_data_r, bnf_level) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Reactive to store calculated results
		#copresc_results <- reactiveVal(NULL)

		# Reactive value for current page
		current_page <- reactiveVal(1)

		# Calculate co-prescription ORs when button clicked
		# observe({
		# 	req(prescriptions_r(), patient_data_r())
		#
		# 	tryCatch({
		# 		results <- calculate_coprescription_ors(
		# 			prescriptions = prescriptions_r(),
		# 			min_prevalence = input$min_prevalence / 100,
		# 			min_coprescription = input$min_coprescription / 100
		# 		)
		#
		# 		copresc_results(results)
		#
		# 		showNotification(
		# 			sprintf("Analysis complete: %d drug pairs found",
		# 							nrow(results)),
		# 			type = "message",
		# 			duration = 5
		# 		)
		# 	}, error = function(e) {
		# 		showNotification(
		# 			paste("Error calculating ORs:", e$message),
		# 			type = "error",
		# 			duration = 10
		# 		)
		# 	})
		# })

		copresc_results <- reactive({
			req(prescriptions_r())
			message("Recalculating co-presc ors")
			calculate_coprescription_ors(
				prescriptions = prescriptions_r(),
				min_prevalence = input$min_prevalence / 100,
				min_coprescription = input$min_coprescription / 100
			)
		})

		display_data <- reactive({
			req(copresc_results())

			data <- copresc_results()

			# Check for missing values and handle them
			data <- data[!is.na(case_or) & !is.na(control_or)]
			if (input$or_min_filter == TRUE) {
				data <- data[case_or > 1 & control_or >= 1]
			}

			if (nrow(data) == 0) {
				return(NULL)
			}

			# Format OR columns
			data[, case_or_display := sprintf("%.2f (%.2f-%.2f)",
																				case_or, case_ci_lower, case_ci_upper)]
			data[case_p < 0.05, case_or_display := paste0(case_or_display, "*")]
			data[, control_or_display := sprintf("%.2f (%.2f-%.2f)",
																					 control_or, control_ci_lower, control_ci_upper)]
			data[control_p < 0.05, control_or_display := paste0(control_or_display, "*")]

			# Create display table with abs_or_diff for sorting and arrow column
			data[, .(
				drug1,
				drug2,
				case_or = case_or_display,
				direction = case_or - control_or,  # Positive = cases higher
				control_or = control_or_display,
				or_diff
			)]

		})

		# Reset to page 1 when sort changes
		observeEvent(input$forest_sort, {
			current_page(1)
		})

		observeEvent(display_data(), {
			current_page(1)
		}, suspended = TRUE)

		# Previous page button
		observeEvent(input$prev_page, {
			if (current_page() > 1) {
				current_page(current_page() - 1)
			}
		})

		# Next page button
		observeEvent(input$next_page, {
			req(forest_data())
			total_pages <- ceiling(nrow(forest_data()) / 25)
			if (current_page() < total_pages) {
				current_page(current_page() + 1)
			}
		})

		# Sorted and filtered data for forest plot
		forest_data <- reactive({
			req(copresc_results())

			data <- copy(copresc_results())

			# Filter for valid data
			data <- data[!is.na(case_or) & !is.na(control_or) &
									 	is.finite(case_or) & is.finite(control_or)]

			if (input$or_min_filter == TRUE) {
				data <- data[case_or > 1 & control_or >= 1]
			}

			if (nrow(data) == 0) {
				return(NULL)
			}

			if (input$forest_sort_asc==TRUE) {
				order_dir <- 1
			} else order_dir <- -1
			setorderv(data, cols = input$forest_sort, order = order_dir)

			return(data)
		})

		# Page info text
		output$page_info <- renderText({
			req(forest_data())
			total_pages <- ceiling(nrow(forest_data()) / 25)
			paste("Page", current_page(), "of", total_pages)
		})

		# Paginated data for current page
		page_data <- reactive({
			req(forest_data())

			data <- forest_data()
			start_row <- (current_page() - 1) * 25 + 1
			end_row <- min(current_page() * 25, nrow(data))

			data <- data[start_row:end_row]

			# Add drug pair label for y-axis
			data[, drug_pair := paste0(drug1, " + ", drug2)]

			return(data)
		})

		# Render table
		output$copresc_table <- renderReactable({

			req(display_data())

			reactable(
				display_data(),
				columns = list(
					drug1 = colDef(name = "Drug 1", minWidth = 150),
					drug2 = colDef(name = "Drug 2", minWidth = 150),
					case_or = colDef(name = "Case OR (95% CI)", minWidth = 150),
					direction = colDef(
						name = "",
						width = 60,
						align = "center",
						cell = function(value) {
							if (is.na(value) || !is.finite(value)) return("")
							arrow <- if (value > 0) "→" else "←"
							color <- if (value > 0) "#e74c3c" else "#3498db"  # Red for cases, blue for controls
							htmltools::tags$span(style = paste0("color: ", color, "; font-size: 20px; font-weight: bold;"), arrow)
						}
					),
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
					)
				),
				defaultPageSize = 20,
				searchable = TRUE,
				fullWidth = FALSE,
				theme = reactableTheme(
					borderColor = "#dfe2e5",
					stripedColor = "#f6f8fa"
				)
			)
		})

		# Render heatmap using vegawidget
		output$copresc_heatmap <- renderVegawidget({
			req(copresc_results(), input$heatmap_group)

			data <- copresc_results()

			# Select the OR column based on user choice
			if (input$heatmap_group == "case") {
				data[, or_value := case_or]
				data[, ci_lower_value := case_ci_lower]
				data[, ci_upper_value := case_ci_upper]
				group_label <- "Cases"
			} else {
				data[, or_value := control_or]
				data[, ci_lower_value := control_ci_lower]
				data[, ci_upper_value := control_ci_upper]
				group_label <- "Controls"
			}

			# Filter out rows with missing values
			data <- data[!is.na(or_value) & is.finite(or_value) &
									 	or_value > 1]

			if (nrow(data) == 0) {
				return(NULL)
			}

			create_copresc_or_heatmap(data, group_label)
		})

		# Render forest plot
		output$copresc_forest <- renderVegawidget({
			req(page_data())

			data <- page_data()

			if (nrow(data) == 0) {
				return(NULL)
			}

			create_copresc_forest_plot(data)
		})

		outputOptions(session$output, "copresc_table", suspendWhenHidden = TRUE)
		outputOptions(session$output, "copresc_forest", suspendWhenHidden = TRUE)
		outputOptions(session$output, "copresc_heatmap", suspendWhenHidden = TRUE)

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