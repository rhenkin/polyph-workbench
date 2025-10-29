module_cca_ui <- function(id) {
	ns <- NS(id)
	nav_panel(
		"Case-control Analysis",
		value = "cca",
		accordion(
			open = NULL,
			id = ns("cca_accordion"),
			accordion_panel(
				title = "Load dataset",
				value = "load_dataset_panel",
				icon = bs_icon("database"),
				selectizeInput(
					ns("dataset_list_selection"),
					"Select dataset:",
					choices = NULL,
					options = list(dropdownParent = 'body')
				),
				actionButton(ns("load_saved_dataset"), "Load saved dataset")
			),
			accordion_panel(
				title = "Analysis",
				value = "analysis_panel",
				icon = bs_icon("bar-chart-line"),
				navset_card_tab(
					nav_panel(
						title = "Overview",
						layout_columns(col_widths = c(2,4,4,2),
							verticalLayout(
								value_box(
									showcase_layout = showcase_left_center(width = "200px", max_height = "150px"),
									height = "150px",
									title = "Cases",
									value = textOutput(ns("value_box_cases")),
									theme = "red",
									showcase = bs_icon("people-fill")
								),
								value_box(
									showcase_layout = showcase_left_center(width = "200px", max_height = "150px"),
									height = "150px",
									title = "Controls",
									value = textOutput(ns("value_box_controls")),
									theme = "blue",
									showcase = bs_icon("people-fill")
								)
							),
							card(vegawidgetOutput(ns(
								"pp_pyramid_plot"
							))),
							card(vegawidgetOutput(ns(
								"mltc_pyramid_plot"
							)))
						),
						layout_columns(
							card(full_screen = TRUE,
									 id = ns("recent_presc_card"),
								card_header("Top 10 Latest prescriptions"),
								card_body(
									vegawidgetOutput(ns("top_recentpresc_bar_plot"))
									)
								),
							card(full_screen = TRUE,
									 vegawidgetOutput(
							ns("topten_ltc_bar_plot")
						)), card(full_screen = TRUE,
										 vegawidgetOutput(
							ns("topten_substance_bar_plot")
						)))
					),
					nav_panel(title = "Advanced", card(
						card_header("Prevalence tables"), card_body(
							accordion(
								open = FALSE,
								accordion_panel(title = "Long-term conditions", value = "ltc_prev_tables", navset_tab(
									nav_panel(
										"Group prevalences",
										virtualSelectInput(
											ns("ltc_freq_strat_variable"),
											label = "Select a subset (optional):",
											choices = NULL,
											autoSelectFirstOption = FALSE,
											search = FALSE,
											showValueAsTags = TRUE,
											disableOptionGroupCheckbox = TRUE,
											multiple = FALSE,
											dropboxWrapper = "body"
										),
										dataTableOutput(ns("ltc_freq_table"))
									)
								)),
								accordion_panel(title = "Prescriptions", value = "presc_prev_tables", navset_tab(
									nav_panel(
										"Group prevalences",
										virtualSelectInput(
											ns("presc_freq_strat_variable"),
											label = "Select a subset (optional):",
											choices = NULL,
											autoSelectFirstOption = FALSE,
											search = FALSE,
											showValueAsTags = TRUE,
											disableOptionGroupCheckbox = TRUE,
											multiple = FALSE,
											dropboxWrapper = "body"
										),
										dataTableOutput(ns("presc_freq_table"))
									)
								)),
								accordion_panel(
									title = "Prescription prevalence per LTC",
									value = "presc_ltc_prev",
									card(
										full_screen = TRUE,
										height = "60em",
										uiOutput(ns("ltc_dropdown_ui")),
										# fluidRow(column(6, uiOutput(ns("ltc_dropdown_ui"))),
										# 				 column(5, numericInput(ns("outcome_age_filter"), label = "Outcome age filter:", value = 100, min = 16, max = 100)),
										# 				 column(1, textOutput(ns("selected_pats")))),
										reactableOutput(ns("presc_by_ltc"), height = "50em")
									)
								),
							)
						)
					)),
					full_screen  = TRUE
				)

			)
		)
	)
}

module_cca_server <- function(id,
									matched_data ,
									gold_patient,
									gold_ltc,
									gold_cp,
									acute_presc_df,
									chapter_menu_data
) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns
		# cases <- matched_data[[1]]
		# controls <- matched_data[[2]]

		cases <- reactiveVal()
		controls <- reactiveVal()

		cases_controls_r <- reactiveVal()


		saved_dataset_r <- reactiveVal()

		observe({
			updateSelectizeInput(session, "dataset_list_selection", choices = list.files("studies", "*.rds"))
		})

		patient_data_r <- reactiveVal()
		prescriptions_r <- reactiveVal()
		ltcs_r <- reactiveVal()

		observe({
			req(input$dataset_list_selection)
			saved_dataset_data <- readRDS(file.path("studies", input$dataset_list_selection))
			saved_dataset_r(saved_dataset_data)

			matched_patids <- saved_dataset_data$matched_patids
			matched_patids[, group := ifelse(treatment==1, "case", "control")]

			patient_data <- saved_dataset_data$all_patient_data
			setkey(patient_data, patid)
			prescriptions <- unique(saved_dataset_data$all_prescriptions)
			setkey(prescriptions, patid, start_date, substance)
			ltcs <- unique(saved_dataset_data$all_ltc)
			setkey(ltcs, patid, eventdate, term)

			umatched_patids <- unique(matched_patids[,.(patid, index_date, group)])
			ltcs <- ltcs[umatched_patids, .(patid, eventdate, age_days, term, group),
									 on = .(patid, eventdate < index_date), nomatch = 0]

			prescriptions <- prescriptions[umatched_patids,
				.(patid, substance = x.substance, index_date = index_date, start_date = x.start_date, stop_date = x.stop_date, duration, group),
									 on = .(patid, start_date <= index_date), nomatch = 0]
			prescriptions <- prescriptions[stop_date >= index_date - 84]

			prescriptions_n <- prescriptions[,list(pp=.N),patid]
			prescriptions_n <- create_value_groups(
				prescriptions_n,
				breaks = c(2, 5, 10),
				right = FALSE,
				value_col = "pp",
				group_col = "pp_group"
			)

			ltcs_n <- ltcs[,list(n_ltc=.N),patid]
			ltcs_n <- create_value_groups(
				ltcs_n,
				breaks = c(2, 5, 10),
				right = FALSE,
				value_col = "n_ltc",
				group_col = "mltc_group"
			)

			patient_data <- patient_data[prescriptions_n]
			patient_data <- patient_data[ltcs_n]

			patient_data_r(patient_data)
			prescriptions_r(prescriptions)
			ltcs_r(ltcs)
			cases_controls_r(matched_patids[, .(patid, index_date, substance, group)])
			accordion_panel_close("cca_accordion", "load_dataset_panel")
			accordion_panel_open("cca_accordion", "analysis_panel")
			showNotification(paste("Dataset", input$dataset_list_selection, "loaded successfully!"),
											 type = "message", duration = 3)
		}) |> bindEvent(input$load_saved_dataset)

		output$top_recentpresc_bar_plot <- renderVegawidget({
			req(cases_controls_r())
			cases_controls <- cases_controls_r()

			cases_controls_n <- cases_controls[,.N,.(group, substance)]
			group_totals <- cases_controls[, .(total = uniqueN(patid)), by = group]
			cases_controls_n[group_totals, pct := round(N/total*100, 2), on = "group"]
			n_sub <- 10
			height <- 250
			if (isTruthy(input$recent_presc_card_full_screen)) {
				n_sub <- 20
				height <- 500
			}

			top_sub <- cases_controls_n[pct >= 1 & group=="case"][order(-pct)][1:n_sub, substance]
			cases_controls_n[, `:=`(max_pct = max(pct),
											  diff_label = sprintf("+%.1f%%", abs(pct[group == "case"] - pct[group == "control"]))
												), substance]
			cases_controls_n <- cases_controls_n[substance %in% top_sub]
			cases_controls_n[nchar(substance)>15, substance := paste0(strtrim(substance, 15), "...")]
			grouped_bar_plot(cases_controls_n, "substance", height = height, width = 280) |> as_vegaspec()
		})

		output$value_box_cases <- renderText({
			req(patient_data_r())
			patient_data <- patient_data_r()
			prettyNum(nrow(patient_data[treatment==1]), big.mark = ",")
		})

		output$value_box_controls <- renderText({
			req(patient_data_r())
			patient_data <- patient_data_r()
			prettyNum(nrow(patient_data[treatment==0]), big.mark = ",")
		})

		output$pp_pyramid_plot <- renderVegawidget({
			req(patient_data_r())
			patient_data <- patient_data_r()
			upper_band <- round(quantile(patient_data$n_ltc, 0.75)*1.5)
			patient_data[, pp_band := ifelse(pp >= upper_band,  paste0(upper_band,"+"), as.character(pp))]
			plot_data <- patient_data[, .N,. (treatment, pp_band)]
			cases_n <- nrow(patient_data[treatment==1])
			controls_n <- nrow(patient_data[treatment==0])
			plot_data[, pct := 0]
			plot_data[treatment==0, pct := N/controls_n]
			plot_data[treatment==1, pct := N/cases_n]
			pyramid_plot(plot_data, "pp_band",
									 rev(c(2:upper_band,paste0(upper_band,"+"))),
									 side_width = 175,
									 title = "Polypharmacy burden") |>
				as_vegaspec()
		})

		output$mltc_pyramid_plot <- renderVegawidget({
			req(patient_data_r())
			patient_data <- patient_data_r()
			upper_band <- round(quantile(patient_data$n_ltc, 0.75)*2)
			patient_data[, mltc_band := ifelse(n_ltc >= upper_band, paste0(upper_band,"+"), as.character(n_ltc))]
			plot_data <- patient_data[, .N,. (treatment, mltc_band)]
			cases_n <- nrow(patient_data[treatment==1])
			controls_n <- nrow(patient_data[treatment==0])
			plot_data[, pct := 0]
			plot_data[treatment==0, pct := N/controls_n]
			plot_data[treatment==1, pct := N/cases_n]
			pyramid_plot(plot_data, "mltc_band",
									 rev(c(2:upper_band,paste0(upper_band,"+"))),
									 title = "MLTC burden") |>
				as_vegaspec()
		})

		# cases_controls <- reactive({
		# 	rbind(cases()[, .(patid, index_date, substance, group = "case")],
		# 				controls()[, .(patid, index_date=control_index_date, substance, group = "control")]
		# 	)
		# })

		observe({
			req(patient_data_r())
			df <- patient_data_r()

			heatmap_vars <- c("sex", "eth_group", "imd_quintile", "pp_group", "mltc_group")
			var_labels <- c("Sex", "Ethnicity", "IMD quintile", "# PP", "# LTC")

			choices <- setNames(
				lapply(seq_along(heatmap_vars), function(i) {
					var <- heatmap_vars[i]
					unique_vals <- sort(unique(df[[var]]))

					# Create named list: display_name = encoded_value
					encoded_list <- setNames(
						paste0(var, "#", unique_vals),  # values (encoded)
						unique_vals                     # names (display labels)
					)
					return(encoded_list)
				}),
				var_labels  # group names
			)

			updateVirtualSelect("ltc_freq_strat_variable", choices = choices)
			updateVirtualSelect("presc_freq_strat_variable", choices = choices)
		})

		ltc_freq_df <- reactive({
			req(ltcs_r())
			ltcs <- ltcs_r()
			patient_data <- patient_data_r()
			strat_var <- input$ltc_freq_strat_variable

			if (strat_var != "") {
				parts <- strsplit(strat_var, "#")[[1]]
				column_name <- parts[1]
				filter_value <- parts[2]
				selected_patids <- patient_data[get(column_name) == filter_value, patid]
				ltcs <- ltcs[patid %in% selected_patids]
			}

			term_freq <- ltcs[, .N, .(group, term)]
			group_totals <- ltcs[, .(total = uniqueN(patid)), by = group]
			term_freq[group_totals, pct := round(N/total*100, 2), on = "group"]
			setorder(term_freq, -pct)
			term_freq
		})

		output$topten_ltc_bar_plot <- renderVegawidget({
			ltc_freq_df <- ltc_freq_df()

			ratios <- ltc_freq_df[, .(
				case_pct = pct[group == "case"],
				control_pct = pct[group == "control"]
			), by = term]

			ratios[, `:=`(ratio = case_pct / control_pct,
										diff = case_pct-control_pct)]

			top_ltc <- ratios[case_pct >= 10][order(-case_pct)][1:10, term]
			# bottom_ltc <- ratios[case_pct >= 10][order(-diff)][1:5, term]

			term_filtered <- ltc_freq_df[term %in% top_ltc]

			term_filtered[, `:=`(
				diff = pct[group == "case"] - pct[group == "control"],
				diff_label = sprintf("+%.1f%%", abs(pct[group == "case"] - pct[group == "control"])),
				max_pct = max(pct)
			), by = term]

			term_filtered[nchar(term)>15, term := paste0(strtrim(term, 15), "...")]

			grouped_bar_plot(term_filtered, "term",
											 # title = "Top 5 conditions with highest difference in each group"
											 title = "Top 10 conditions most prevalent in cases"
											 ) |> as_vegaspec()
		})

		output$ltc_freq_table <- renderDataTable({
			table_data <- ltc_freq_df()
			table_data_wide <- dcast(table_data, term ~ group, value.var = "pct", fill = 0)
			table_data_wide <- table_data_wide[case > 0.5 & control > 0.5]
			table_data_wide[, ratio := round(case/control, digits = 2)]
			colnames(table_data_wide) <- c("Condition", "Cases (%)", "Controls (%)", "Case/Control Ratio")
			datatable(table_data_wide[order(-`Case/Control Ratio`)],
								class = list(stripe = FALSE), rownames = FALSE)
		})

		presc_freq_df <- reactive({
			req(prescriptions_r())
			presc <- prescriptions_r()
			patient_data <- patient_data_r()
			strat_var <- input$presc_freq_strat_variable

			if (strat_var != "") {
				parts <- strsplit(strat_var, "#")[[1]]
				column_name <- parts[1]
				filter_value <- parts[2]
				selected_patids <- patient_data[get(column_name) == filter_value, patid]
				presc <- presc[patid %in% selected_patids]
			}

			substance_freq <- presc[, .N, .(group, substance)]
			group_totals <- presc[, .(total = uniqueN(patid)), by = group]
			substance_freq[group_totals, pct := round(N/total*100, 2), on = "group"]
			setorder(substance_freq, -N)
			substance_freq
		})

		output$topten_substance_bar_plot <- renderVegawidget({
			presc_freq_df <- presc_freq_df()

			ratios <- presc_freq_df[, .(
				case_pct = pct[group == "case"],
				control_pct = pct[group == "control"]
			), by = substance]

			ratios[, ratio := case_pct / control_pct]

			top_ltc <- ratios[case_pct >= 5][order(-case_pct)][1:10, substance]

			term_filtered <- presc_freq_df[substance %in% top_ltc]

			term_filtered[, `:=`(
				diff = pct[group == "case"] - pct[group == "control"],
				diff_label = sprintf("+%.1f%%", abs(pct[group == "case"] - pct[group == "control"])),
				max_pct = max(pct)
			), by = substance]

			term_filtered[nchar(substance)>15, substance := paste0(strtrim(substance, 15), "...")]

			grouped_bar_plot(term_filtered,
											 "substance",
											 title = "Top 10 Background medications most prevalent in cases") |> as_vegaspec()
		})

		output$presc_freq_table <- renderDataTable({
			table_data <- presc_freq_df()
			table_data_wide <- dcast(table_data, substance ~ group, value.var = "pct", fill = 0)
			table_data_wide <- table_data_wide[case > 0.5 & control > 0.5]
			table_data_wide[, ratio := round(case/control, digits = 2)]
			colnames(table_data_wide) <- c("Substance", "Cases (%)", "Controls (%)", "Case/Control Ratio")
			datatable(table_data_wide[order(-`Case/Control Ratio`)],
								class = list(stripe = FALSE), rownames = FALSE)
		})

		output$ltc_dropdown_ui <- renderUI({
			req(ltcs_r())
			ltcs <- ltcs_r()

			virtualSelectInput(ns("ltc_dropdown"),
												 label = "Select 1 or more LTCs:",
												 choices = with(chapter_menu_data[ltc %in% unique(ltcs$term)],
												 							 split(ltc, body_system)),
												 multiple = TRUE,
												 search = TRUE)
		})

		output$presc_by_ltc <- renderReactable({
			req(ltcs_r())
			req(input$ltc_dropdown)
			ltcs <- ltcs_r()
			prescriptions <- prescriptions_r()

			patids <- unique(ltcs[term %in% input$ltc_dropdown, patid])

			# Patients WITH the disease
			presc_freq <- prescriptions[patid %in% patids,
																	list(N_with_disease = uniqueN(patid),
																			 Prevalence = round(100*(uniqueN(patid)/length(patids)), digits=2),
																			 `Median Duration (years)`= round(median(duration/365.2), digits = 2),
																			 `IQR (Q1-Q3)` = paste0("(",
																			 											 round(quantile(duration/365.2, 0.25, na.rm=TRUE), 2), " - ",
																			 											 round(quantile(duration/365.2, 0.75, na.rm=TRUE), 2), ")")),
																	.(group, substance)]


			result <- dcast(presc_freq, substance ~ group, value.var = "Prevalence", fill = 0)
			result <- result[case >= 1 & control >= 1]
			result[, `:=`(
				Prevalence_Ratio = round(case / control, digits = 2)
			)]

			result[is.infinite(Prevalence_Ratio), Prevalence_Ratio := 0]

			# Calculate 95% CI for the ratio using log method
			# result[, `:=`(
			# 	p1 = N_with_disease / total_with_disease,
			# 	p2 = N_without_disease / total_without_disease
			# )]

			# result[, `:=`(
			# 	log_ratio = log(Prevalence_Ratio),
			# 	se_log_ratio = sqrt((1/N_with_disease) - (1/total_with_disease) +
			# 												(1/N_without_disease) - (1/total_without_disease))
			# )]
			#
			# result[, `:=`(
			# 	CI_lower = round(exp(log_ratio - 1.96 * se_log_ratio), digits = 2),
			# 	CI_upper = round(exp(log_ratio + 1.96 * se_log_ratio), digits = 2)
			# )]

			# result[,
			# 			 CI_95 := paste0("(", CI_lower, " - ", CI_upper, ")")
			# ]

			# result[(CI_lower > 1.0 | CI_upper < 1.0), substance := paste0(substance, "*")]
			# Clean up intermediate columns if desired
			# result[, c("p1", "p2", "log_ratio", "se_log_ratio", "CI_lower", "CI_upper") := NULL]

			# result[order(-Prevalence_Ratio), .(substance, Prevalence, Prevalence_Unselected, Prevalence_Ratio, CI_95, `Median Duration (years)`, `Median Duration unselected (years)`, `IQR (Q1-Q3)`, `IQR unsel. (Q1-Q3)`)]
			reactable(result,showPageInfo = FALSE, defaultPageSize = 15)
		})


	})
}