#' @export
module_overview_ui <- function(id) {
	ns <- NS(id)
	card(
		card_header("Summary"),
		accordion(
			open = FALSE,
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
				title = "Demographic distribution",
				value = "demodists",
				navset_tab(
					nav_panel("Distribution per polypharmacy levels", gt_output(ns("pp_demodist_table"))),
					nav_panel("Cohort table", gt_output(ns("demodist_table"))),

				)
			),
			accordion_panel(
				title = "Paired distributions of PP",
				value = "pairedpp",
				fluidRow(
					column(width = 2, selectizeInput(ns("paired_x_variable"), "X variable", choices = c("eth_group", "sex", "imd_quintile","age_group", "tto_group"))),
					column(width = 1,
								 actionButton(ns("swap_vars"), label = "↔",
								 						 style = "margin-top: 30px;display: block; margin-left: auto; margin-right: auto;")),
					column(width = 2, selectizeInput(ns("paired_y_variable"), "Y variable", choices = c("eth_group", "sex", "imd_quintile","age_group", "tto_group"), selected = "sex"))
				),
				fluidRow(
					vegawidgetOutput(ns("paired_heatmap"))
				)
			),
			accordion_panel(
				title = "Association between PP and time-to-outcome",
				value = "pptime",
				dataTableOutput(ns("pp_cor_table"))
			)
		)
	)
}

#' @export
module_overview_server <- function(id, outcome_prescriptions, patient_data, outcome_ltc_events, selected_outcome, pp_groups_data) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Helper functions
		create_age_groups <- function(dt, n_groups = 5) {
			dt_copy <- copy(dt)

			# Calculate break points
			breaks <- quantile(dt_copy$age_days/365.25, probs = seq(0, 1, length.out = n_groups + 1))

			# Create labels showing the ranges
			labels <- paste(round(breaks[-length(breaks)]), "to", round(breaks[-1]), "years")

			# Add group column using the range labels
			dt_copy[, age_group := cut(age_days/365.25,
																 breaks = breaks,
																 labels = labels,
																 include.lowest = TRUE,
																 ordered_result = TRUE)]
			setkey(dt_copy, patid)
			return(dt_copy)
		}

		# Core reactive data sources
		polypharmacy_counts <- reactive({
			df <- outcome_prescriptions()
			df <-
				unique(df[, list(
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
			df$age <- df$age_days/365.25
			breaks <- c(0,40,60,80,Inf)
			labels <- c("<=40", "41-60", "61-80", ">80")
			df[, age_group := cut(age_days/365.25,
														breaks = breaks,
														labels = labels,
														include.lowest = TRUE,
														ordered_result = TRUE,)]

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
		pp_frequencies <- function(dt, group_var) {
			# Get all possible pp values
			all_pp <- seq_len(max(dt$pp))

			# Create the frequency table with zeros for missing combinations
			freq_table <- dt[, as.list(table(factor(pp, levels=all_pp))), by=group_var]

			# Convert to the string format you want
			result <- freq_table[, list(
				pp_label = paste0(all_pp, collapse=","),
				pp_values = paste0(unlist(.SD), collapse=",")
			), by = group_var]
			return(result)
		}

		output$demodist_table <- render_gt({
			pp_df <- pp_demog_table()
			total_patids <- length(unique(pp_df$patid))

			sex_df <- pp_df[, .N, sex]
			sex_df[, group := "Sex"]
			sex_df[, pct_total := signif(N/total_patids, digits = 2)]

			sex_df <- merge(sex_df, pp_frequencies(pp_df, "sex"), by = "sex")
			eth_df <- pp_df[, .N, eth_group]
			eth_df[, group := "Ethnic group"]
			eth_df[, pct_total := signif(N/total_patids, digits = 2)]
			eth_df <- merge(eth_df, pp_frequencies(pp_df, "eth_group"), by = "eth_group")

			imd_df <- pp_df[, .N, imd_quintile]
			imd_df[, group := "IMD Quintile"]
			imd_df[, pct_total := signif(N/total_patids, digits = 2)]
			imd_df <- merge(imd_df, pp_frequencies(pp_df, "imd_quintile"), by = "imd_quintile")

			setorder(imd_df, imd_quintile)
			age_df <- pp_df[, .N, age_group]
			age_df[, group := "Age at outcome"]
			age_df[, pct_total := signif(N/total_patids, digits = 2)]
			age_df <- merge(age_df, pp_frequencies(pp_df, "age_group"), by = "age_group")

			setorder(age_df, age_group)

			to_print <-
				rbindlist(
					list(
						sex_df,
						eth_df,
						imd_df,
						age_df
					),
					use.names = FALSE)
			colnames(to_print) <- c("Category", "N", "group", "%", "pp_labels", "pp_values")

			gt(to_print,row_group_as_column = TRUE, groupname_col = "group") |>
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
		output$pp_cor_table <- renderDataTable({
			df <- outcome_prescriptions()
			setkey(df, patid)
			results <- df[, list(
				time_to_outcome = as.IDate(eventdate) - max(as.IDate(start_date)),
				pp = .N
			), patid] |> unique()
			summary_stats <- results[, .(
				n = .N,
				mean_time = mean(time_to_outcome, na.rm = TRUE),
				median_time = median(time_to_outcome, na.rm = TRUE),
				sd_time = sd(time_to_outcome, na.rm = TRUE),
				min_time = min(time_to_outcome, na.rm = TRUE),
				max_time = max(time_to_outcome, na.rm = TRUE)
			), by = pp]
			setorder(summary_stats, pp)
			summary_stats
		})

	})
}