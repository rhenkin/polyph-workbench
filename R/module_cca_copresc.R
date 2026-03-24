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
			title = "Co-prescription Analysis ",
			value = "copresc_results_panel",
			icon = bs_icon("graph-up"),

			# Help text
			div(
				class = "alert alert-info",
				style = "margin-bottom: 15px;",
				HTML("<strong>How to interpret results:</strong> An OR of 2.5
             for Drug A + Drug B means cases are 2.5x more likely to have been prescribed both
             drugs together compared to controls.")
			),

			layout_columns(
				col_widths = c(4, 4, 4),
				numericInput(
					ns("min_prevalence"),
					"Minimum individual drug prevalence (%):",
					value = 2,
					min = 0.1,
					max = 10,
					step = 0.1
				),
				numericInput(
					ns("min_coprescription"),
					"Minimum co-prescription prevalence (%):",
					value = 1,
					min = 0.1,
					max = 10,
					step = 0.1
				),
				checkboxInput(
					ns("or_min_filter"),
					"Show only OR > 1?",
					value = FALSE
				)
			),

			# Substance selection dropdown + calculate button
			layout_columns(
				col_widths = c(8, 4),
				div(
					tags$label("Restrict to specific substances (optional):"),
					uiOutput(ns("substance_filter_ui"))
				),
				div(
					style = "padding-top: 24px;",
					actionButton(
						ns("calculate_btn"),
						"Calculate ORs",
						icon = icon("calculator"),
						class = "btn-primary"
					)
				)
			),

			navset_card_tab(
				id = ns("results_tabs"),

				# Forest Plot Tab
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
									"Odds Ratio" = "or",
									"OR Distance from 1" = "or_diff"
								),
								selected = "or_diff"
							),
							checkboxInput(
								ns("forest_sort_asc"),
								"Ascending order?",
								value = FALSE
							)
						),
						vegawidgetOutput(ns("copresc_forest"), height = "800px"),
						layout_column_wrap(
							width = "100px",
							fixed_width = TRUE,
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

				# Table Tab
				nav_panel(
					title = "Table",
					card_body(
						downloadButton(ns("download_table"), "Download Table"),
						reactableOutput(ns("copresc_table"))
					)
				),

				# Heatmap Tab
				nav_panel(
					title = "Heatmap",
					card_body(
						div(
							class = "alert alert-secondary",
							style = "margin-bottom: 15px;",
							"Heatmap shows case-control odds ratios for all drug pairs.
              Darker colors indicate stronger associations."
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

		# Reactive value for current page
		current_page <- reactiveVal(1)

		# -----------------------------------------------------------------------
		# Substance filter dropdown (same pattern as prevalence module)
		# -----------------------------------------------------------------------

		output$substance_filter_ui <- renderUI({
			req(prescriptions_r(), bnf_level())

			virtualSelectInput(
				ns("substance_filter"),
				label = NULL,
				choices = create_bnf_dropdown_choices(prescriptions_r(), bnf_level()),
				multiple = TRUE,
				search = TRUE,
				showSelectedOptionsFirst = TRUE,
				placeholder = "All substances (leave empty to include all)"
			)
		})

		# -----------------------------------------------------------------------
		# Filtered prescriptions (apply substance selection if provided)
		# -----------------------------------------------------------------------

		filtered_prescriptions <- reactive({
			req(prescriptions_r())
			presc <- prescriptions_r()

			if (isTruthy(input$substance_filter)) {
				presc <- presc[substance %in% input$substance_filter]
			}

			presc
		})

		# -----------------------------------------------------------------------
		# Calculate co-prescription ORs — only fires on button click
		# -----------------------------------------------------------------------

		copresc_results <- reactive({
			req(filtered_prescriptions())
			message("Calculating case-control co-prescription ORs")

			calculate_coprescription_ors(
				prescriptions = filtered_prescriptions(),
				min_prevalence = input$min_prevalence / 100,
				min_coprescription = input$min_coprescription / 100
			)
		}) |> bindEvent(input$calculate_btn)

		# Prepare display data
		display_data <- reactive({
			req(copresc_results())

			data <- copresc_results()

			# Filter out missing ORs
			data <- data[!is.na(or)]

			# Optional filter for OR > 1
			if (input$or_min_filter == TRUE) {
				data <- data[or > 1]
			}

			if (nrow(data) == 0) {
				return(NULL)
			}

			# Format OR with CI and significance marker
			data[, or_display := sprintf("%.2f (%.2f-%.2f)", or, ci_lower, ci_upper)]
			data[p_adjusted < 0.05, or_display := paste0(or_display, "*")]

			# Add absolute OR difference from null (for sorting)
			data[, abs_or_diff := abs(or - 1)]

			return(data)
		})

		# Table output
		output$copresc_table <- renderReactable({
			req(display_data())

			data <- display_data()

			# Select columns for display
			table_data <- data[, .(
				drug1,
				drug2,
				or_display,
				pct_case_both,
				pct_control_both,
				n_case_both,
				n_control_both,
				p_adjusted
			)]

			# Sort based on user selection
			# sort_col <- input$forest_sort
			# sort_asc <- input$forest_sort_asc
			#
			# if (sort_col == "or_diff") {
			# 	if ("abs_or_diff" %in% colnames(data)) {
			# 		if (sort_asc) {
			# 			setorder(table_data, -abs_or_diff)
			# 		} else {
			# 			setorder(table_data, abs_or_diff)
			# 		}
			# 	}
			# } else if (sort_col %in% c("drug1", "drug2")) {
			# 	if (sort_asc) {
			# 		setorderv(table_data, sort_col)
			# 	} else {
			# 		setorderv(table_data, sort_col, -1)
			# 	}
			# } else {
			# 	if (sort_asc) {
			# 		setorder(table_data, or)
			# 	} else {
			# 		setorder(table_data, -or)
			# 	}
			# }

			reactable(
				table_data,
				columns = list(
					drug1 = colDef(name = "Drug 1", minWidth = 150),
					drug2 = colDef(name = "Drug 2", minWidth = 150),
					or_display = colDef(name = "OR (95% CI)", minWidth = 130),
					pct_case_both = colDef(name = "Cases (%)", format = colFormat(digits = 1)),
					pct_control_both = colDef(name = "Controls (%)", format = colFormat(digits = 1)),
					n_case_both = colDef(name = "N Cases"),
					n_control_both = colDef(name = "N Controls"),
					p_adjusted = colDef(name = "p (adj.)", format = colFormat(digits = 3))
				),
				defaultPageSize = 15,
				compact = TRUE,
				highlight = TRUE,
				striped = TRUE
			)
		})

		# Forest plot pagination
		items_per_page <- 20

		sorted_display_data <- reactive({
			req(display_data())
			data <- copy(display_data())

			sort_col <- input$forest_sort
			sort_asc <- input$forest_sort_asc

			if (sort_col == "or_diff") {
				if (sort_asc) setorder(data, abs_or_diff) else setorder(data, -abs_or_diff)
			} else if (sort_col %in% c("drug1", "drug2")) {
				if (sort_asc) setorderv(data, sort_col) else setorderv(data, sort_col, -1)
			} else {
				if (sort_asc) setorder(data, or) else setorder(data, -or)
			}

			data
		})

		total_pages <- reactive({
			req(sorted_display_data())
			ceiling(nrow(sorted_display_data()) / items_per_page)
		})

		observe({
			req(copresc_results())
			current_page(1)
		})

		observeEvent(input$next_page, {
			if (current_page() < total_pages()) current_page(current_page() + 1)
		})

		observeEvent(input$prev_page, {
			if (current_page() > 1) current_page(current_page() - 1)
		})

		output$page_info <- renderText({
			req(total_pages())
			sprintf("Page %d / %d", current_page(), total_pages())
		})

		paged_data <- reactive({
			req(sorted_display_data())
			data <- sorted_display_data()
			start <- (current_page() - 1) * items_per_page + 1
			end   <- min(current_page() * items_per_page, nrow(data))
			page_data <- data[start:end]
			# pair_label required by create_copresc_forest_plot()
			page_data[, pair_label := paste0(drug1, " + ", drug2)]
			page_data
		})

		output$copresc_forest <- renderVegawidget({
			req(paged_data())
			create_copresc_forest_plot(paged_data()) |> as_vegaspec()
		})

		output$copresc_heatmap <- renderVegawidget({
			req(display_data())
			# create_copresc_or_heatmap() expects a column named or_value
			data_symmetric <- rbind(
				display_data()[, .(drug1, drug2, or_value = or)],
				display_data()[, .(drug1 = drug2, drug2 = drug1, or_value = or)]
			)
			create_copresc_or_heatmap(data_symmetric, "Case-Control Comparison") |> as_vegaspec()
		})

		output$download_table <- downloadHandler(
			filename = function() paste0("coprescription_ors_", Sys.Date(), ".csv"),
			content  = function(file) {
				req(display_data())
				fwrite(display_data(), file)
			}
		)
	})
}