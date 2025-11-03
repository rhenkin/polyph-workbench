paired_choices <- c("Ethnic group" = "eth_group", "Sex" = "sex", "IMD quintile" = "imd_quintile", "Outcome age" = "age_group", "Time-to-outcome" = "tto_group", "# LTCs" = "mltc_group")

#' @export
module_overview_ui <- function(id) {
	ns <- NS(id)
	card(
		card_header("Summary"),
		accordion(
			open = FALSE,
			accordion_panel(
				title = "Demographic distribution",
				value = "demodists",
				navset_tab(
					nav_panel("Cohort table", gt_output(ns("demodist_table"))),
					nav_panel("Distribution across polypharmacy levels", gt_output(ns("pp_demodist_table")))
				)
			),
			accordion_panel(
				title = "Overview of polypharmacy distribution for selected outcome",
				value = "pphists",
				# layout_columns(
				#   col_widths = c(6, 6),
				fluidRow(
					vegawidgetOutput(ns("pp_histogram")),
					vegawidgetOutput(ns("tto_histogram")),
					vegawidgetOutput(ns("pp_sex_dist")),
					vegawidgetOutput(ns("pp_age_group_dist")),
					vegawidgetOutput(ns("pp_eth_dist")),
					vegawidgetOutput(ns("pp_imd_dist"))
				)
			),
			accordion_panel(
				title = "Paired distributions of PP",
				value = "pairedpp",
				fluidRow(
					column(width = 2, selectizeInput(ns("paired_x_variable"), "X variable", choices = paired_choices)),
					column(width = 1,
								 actionButton(ns("swap_vars"), label = "↔",
								 						 style = "margin-top: 30px;display: block; margin-left: auto; margin-right: auto;")),
					column(width = 2, selectizeInput(ns("paired_y_variable"), "Y variable", choices = paired_choices, selected = "sex"))
				),
				fluidRow(
					vegawidgetOutput(ns("paired_heatmap"))
				)
			),
			accordion_panel(
				title = "Time-to-outcome versus PP burden from last prescription",
				value = "pptime",
				vegawidgetOutput(ns("tto_curve"))
			)
		)
	)
}

#' @export
module_overview_server <- function(id, outcome_prescriptions, patient_data, outcome_ltc_events, selected_outcome, pp_groups_data) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Core reactive data sources
		polypharmacy_counts <- reactive({
			df <- outcome_prescriptions()
			df <-
				unique(df[, list(
					mltc_group,
					pp = length(unique(substance)),
					first_prescription = min(start_date),
					last_prescription = max(start_date),
					time_to_outcome = as.numeric(eventdate - min(start_date))
				) , .(patid)])
			setkey(df, patid)
		})

		pp_demog_table <- reactive({
			req(outcome_ltc_events(), polypharmacy_counts(), patient_data())
			df1 <- outcome_ltc_events()
			df2 <- polypharmacy_counts()
			df3 <- patient_data()
			df <- merge(df1[df2,], df3, by = "patid", all.x = TRUE, all.y = FALSE)
			df <- df[patid %in% df3$patid]
			# df$age <- df$age_days/365.25
			# breaks <- c(0,44,64,84,Inf)
			# labels <- c("<=44", "45-64", "65-84", "85+")
			# df[, age_group := cut(round(age_days/365.25, digits = 0),
			# 											breaks = breaks,
			# 											labels = labels,
			# 											include.lowest = FALSE,
			# 											ordered_result = TRUE,)]

			df[,age_first_prescription := (first_prescription - dob)]
			df[, first_prescription_age_group := cut(age_first_prescription/365.25,
																							 breaks = breaks,
																							 labels = labels,
																							 include.lowest = TRUE,
																							 ordered_result = TRUE,)]
			tto_labels <- c("<=5", "<=10", "10+")
			tto_breaks <- c(0,5,10,Inf)
			df[, tto_group := cut(time_to_outcome,
														breaks = tto_breaks,
														labels = tto_labels,
														include.lowest = TRUE,
														ordered_result = TRUE,)]
			df <- pp_groups_data()[df]
		})

		data_for_histogram <- reactive({
			req(polypharmacy_counts(), patient_data())
			df <- polypharmacy_counts()
			patient_df <- patient_data()
			merge(df, patient_df, by = "patid", all.x = TRUE, all.y = FALSE)
		})

		# 1. PP Histogram
		output$pp_histogram <- renderVegawidget({
			df <- polypharmacy_counts()
			precalc_histogramPlot(df,
														"pp",
														title = "Burden before diagnosis")
		})

		# 2. TTO Histogram
		output$tto_histogram <- renderVegawidget({
			df <- polypharmacy_counts()
			precalc_histogramPlot(df[,time_to_outcome := time_to_outcome/365.25],
														"time_to_outcome",
														title = "Time from first treatment to outcome")
		})

		# 3. PP by Sex Distribution
		output$pp_sex_dist <- renderVegawidget({
			df <- calculateDensityByGroup_dt(data_for_histogram(), "pp", "sex")
			precalc_facetedViolinPlot(df,
																"pp",
																"sex",
																title = "Burden by sex")
		})

		# 4. PP by Age Group Distribution
		output$pp_age_group_dist <- renderVegawidget({
			validate(need(selected_outcome(), "Select an outcome to view age distribution"))
			df <- polypharmacy_counts()
			df2 <- outcome_ltc_events()
			df2 <- create_age_groups(df2, 5)
			df <- df2[df]
			facetedViolinPlot(df,
												"pp",
												"age_group",
												title = "Burden by age groups",
												w = 100)
		})

		# 5. PP by Ethnicity Distribution
		output$pp_eth_dist <- renderVegawidget({
			df <- calculateDensityByGroup_dt(data_for_histogram(), "pp", "eth_group")
			precalc_facetedViolinPlot(df,
																"pp",
																"eth_group",
																title = "Burden by ethnic group",
																w = 100)
		})

		# 6. PP by IMD Distribution
		output$pp_imd_dist <- renderVegawidget({
			df <- calculateDensityByGroup_dt(data_for_histogram(), "pp", "imd_quintile")
			precalc_facetedViolinPlot(df,
																"pp",
																"imd_quintile",
																title = "Burden by IMD quintile")
		})

		# 7. Demographic distribution tables
		output$demodist_table <- render_gt({
			pp_df <- pp_demog_table()
			to_print <- prepare_cohort_demog_data(pp_df)

			gt(to_print, row_group_as_column = TRUE, groupname_col = "group") |>
				tab_style(
					style = list(
						weight = "bold"
					),
					locations = cells_stub()
				) |>
				fmt_percent(
					columns = "%",
					decimals = 2
				) |>
				cols_nanoplot(columns = "pp_values",
											columns_x_vals = "pp_labels",
											autohide = TRUE,
											new_col_name = "Polypharmacy burden distribution")
		})

		output$pp_demodist_table <- render_gt({
			pp_df <- pp_demog_table()
			final_table <- calc_demog_table_by_pp(pp_df)

			# Create the gt table with the right order of groups
			gt(final_table,
				 rowname_col = "category",
				 groupname_col = "group") |>
				row_group_order("Total")
		})

		# 8. Paired Distributions Section
		observeEvent(input$swap_vars, {
			# Store current values
			x_val <- input$paired_x_variable
			y_val <- input$paired_y_variable

			# Update the inputs with swapped values
			updateSelectizeInput(session, "paired_x_variable", selected = y_val)
			updateSelectizeInput(session, "paired_y_variable", selected = x_val)
		})

		output$paired_heatmap <- renderVegawidget({
			dt <- pp_demog_table()
			pairedBarCharts(dt, input$paired_x_variable, input$paired_y_variable)
		})

		# 9. Association between PP and time-to-outcome
		output$tto_curve <- renderVegawidget({
			df <- outcome_prescriptions()
			#setkey(df, patid)
			results <- df[, list(
				time_to_outcome = as.IDate(eventdate) - max(as.IDate(start_date)),
				pp = .N
			), patid] |> unique()

			summary_stats <- results[order(pp), c(.N,as.list(quantile(time_to_outcome, c(0.25,0.5,0.75)))), pp]
			colnames(summary_stats) <- c("pp", "N", "iqr1", "median", "iqr2")

			tto_line_plot(summary_stats) |> as_vegaspec()
		})

	})
}