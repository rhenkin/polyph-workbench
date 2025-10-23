module_cca_ui <- function(id) {
	ns <- NS(id)
	nav_panel("Case-control Analysis",
						card(
							card_header("Long-term conditions"),
							card_body(

								selectizeInput(ns("study_list_selection"), "Select study:", choices = NULL,
															 options = list(dropdownParent = 'body')),
								actionButton(ns("load_saved_study"), "Load saved study")
							)
						),
						card(
							card_header("Prevalence tables"),
							card_body(
								accordion(
									open = FALSE,
									accordion_panel(
										title = "Long-term conditions",
										value = "ltc_prev_tables",
										navset_tab(
											nav_panel("Group prevalences",
																dataTableOutput(ns("ltc_freq_table"))
											),
											nav_panel("Stratified prevalences",
																virtualSelectInput(
																	ns("ltc_freq_strat_variable"),
																	label = "Select group to filter:",
																	choices = NULL,
																	autoSelectFirstOption = FALSE,
																	search = FALSE,
																	showValueAsTags = TRUE,
																	disableOptionGroupCheckbox = TRUE,
																	multiple = FALSE, dropboxWrapper = "body"
																),
																)
										)
									)
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
									acute_presc_df
) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns
		# cases <- matched_data[[1]]
		# controls <- matched_data[[2]]

		cases <- reactiveVal()
		controls <- reactiveVal()

		cases_controls <- reactiveVal()


		saved_study_r <- reactiveVal()

		observe({
			updateSelectizeInput(session, "study_list_selection", choices = list.files("studies", "*.rds"))
		})

		patient_data_r <- reactiveVal()
		prescriptions_r <- reactiveVal()
		ltcs_r <- reactiveVal()

		observe({
			req(input$study_list_selection)
			req(input$load_saved_study)
			saved_study_data <- readRDS(file.path("studies", input$study_list_selection))
			saved_study_r(saved_study_data)

			matched_patids <- saved_study_data$matched_patids
			matched_patids[, group := ifelse(treatment==1, "case", "control")]

			patient_data <- saved_study_data$all_patient_data
			setkey(patient_data, patid)
			prescriptions <- unique(saved_study_data$all_prescriptions)
			setkey(prescriptions, patid, start_date, substance)
			ltcs <- unique(saved_study_data$all_ltc)
			setkey(ltcs, patid, eventdate, term)

			umatched_patids <- unique(matched_patids[,.(patid, index_date, group)])
			ltcs <- ltcs[umatched_patids, .(patid, eventdate, age_days, term, group),
									 on = .(patid, eventdate < index_date), nomatch = 0]

			prescriptions <- prescriptions[umatched_patids,
				.(patid, substance = x.substance, index_date = index_date, start_date = x.start_date, stop_date = x.stop_date, duration, group),
									 on = .(patid, start_date <= index_date), nomatch = 0]
			prescriptions <- prescriptions[stop_date >= index_date - 84]

			prescriptions_n <- prescriptions[,.N,patid]
			prescriptions_n <- create_value_groups(
				prescriptions_n,
				breaks = c(2, 5, 10),
				right = FALSE,
				value_col = "N",
				group_col = "pp_group"
			)

			ltcs_n <- ltcs[,.N,patid]
			ltcs_n <- create_value_groups(
				ltcs_n,
				breaks = c(2, 5, 10),
				right = FALSE,
				value_col = "N",
				group_col = "mltc_group"
			)

			patient_data <- patient_data[prescriptions_n]
			patient_data <- patient_data[ltcs_n]

			patient_data_r(patient_data)
			prescriptions_r(prescriptions)
			ltcs_r(ltcs)
			cases_controls(matched_patids[, .(patid, index_date, substance, group)])

			showNotification(paste("Study", input$study_list_selection, "loaded successfully!"),
											 type = "message", duration = 3)
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
		})


		# ltcs_r <- reactive({
		# 	req(cases_controls())
		# 	gold_ltc[cases_controls(), on = .(patid, eventdate <= index_date), nomatch = 0]
		# })

		ltc_freq_df <- reactive({
			ltcs <- ltcs_r()
			term_freq <- ltcs[, .N, .(group, term)]
			group_totals <- ltcs[, .(total = uniqueN(patid)), by = group]
			term_freq[group_totals, pct := round(N/total*100, 2), on = "group"]
			setorder(term_freq, -N)
			term_freq
		})

		output$ltc_freq_table <- renderDataTable({
			table_data <- ltc_freq_df()
			table_data_wide <- dcast(table_data, term ~ group, value.var = "pct", fill = 0)
			table_data_wide[, ratio := round(case/control, digits = 2)]
		}, rownames = FALSE)


	})
}