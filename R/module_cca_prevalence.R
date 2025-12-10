module_cca_prevalence_ui <- function(id) {
	ns <- NS(id)

	accordion(
		open = FALSE,
		accordion_panel(
			title = "Long-term conditions",
			value = "ltc_prev_tables",
			icon = bs_icon("heart-pulse"),
			navset_tab(
				nav_panel(
					"Prevalence table",
					p("Click on a column name to reorder table"),
					virtualSelectInput(
						ns("ltc_freq_strat_variable"),
						label = "Select a subset of patients (optional):",
						choices = NULL,
						autoSelectFirstOption = FALSE,
						search = FALSE,
						showValueAsTags = TRUE,
						disableOptionGroupCheckbox = TRUE,
						multiple = FALSE,
						dropboxWrapper = "body"
					),
					reactableOutput(ns("ltc_freq_table"), height = "50em")
				),
				nav_panel(
					title = "Prevalence by prescription filter",
					card(
						full_screen = TRUE,
						height = "60em",
						div("Select medications to show the prevalence difference only among the patients taking those medications."),
						div("Table contains conditions with a minimum of 1% prevalence in both cases and controls"),
						uiOutput(ns("presc_dropdown_ui")),
						reactableOutput(ns("ltc_by_presc"), height = "50em")
					)
				)
			)
		),
		accordion_panel(
			title = "Prescriptions",
			value = "presc_prev_tables",
			icon = bs_icon("capsule"),
			navset_tab(
				nav_panel(
					"Prevalence table",
					p("Click on a column name to reorder table"),
					virtualSelectInput(
						ns("presc_freq_strat_variable"),
						label = "Select a subset of patients (optional):",
						choices = NULL,
						autoSelectFirstOption = FALSE,
						search = FALSE,
						showValueAsTags = TRUE,
						disableOptionGroupCheckbox = TRUE,
						multiple = FALSE,
						dropboxWrapper = "body"
					),
					reactableOutput(ns("presc_freq_table"), height = "50em")
				),
				nav_panel(
					title = "Prevalence by condition filter",
					card(
						full_screen = TRUE,
						height = "60em",
						div("Select conditions to show the prevalence difference only among the patients that were diagnosed with those conditions."),
						div("Table contains medications with a minimum of 1% prevalence in both cases and controls"),
						uiOutput(ns("ltc_dropdown_ui")),
						reactableOutput(ns("presc_by_ltc"), height = "50em")
					)
				)
			)
		),
		# NEW PANEL: Recent prescriptions (from case_controls_r)
		accordion_panel(
			title = "Recent prescriptions",
			value = "recent_presc_prev_tables",
			icon = bs_icon("capsule-pill"),
			navset_tab(
				nav_panel(
					"Prevalence table",
					p("Click on a column name to reorder table"),
					virtualSelectInput(
						ns("recent_presc_freq_strat_variable"),
						label = "Select a subset of patients (optional):",
						choices = NULL,
						autoSelectFirstOption = FALSE,
						search = FALSE,
						showValueAsTags = TRUE,
						disableOptionGroupCheckbox = TRUE,
						multiple = FALSE,
						dropboxWrapper = "body"
					),
					reactableOutput(ns("recent_presc_freq_table"), height = "50em")
				),
				nav_panel(
					title = "Stratified prevalence",
					card(
						full_screen = TRUE,
						height = "60em",
						div("Filter patients by both long-term conditions AND background medications simultaneously to see prevalence of recent prescriptions in this subset."),
						div("Table contains medications with a minimum of 1% prevalence in both cases and controls"),
						fluidRow(
							column(6,
										 uiOutput(ns("recent_ltc_dropdown_ui"))
							),
							column(6,
										 uiOutput(ns("recent_bg_presc_dropdown_ui"))
							)
						),
						reactableOutput(ns("recent_by_ltc_and_presc"), height = "50em")
					)
				)
			)
		)
	)
}

module_cca_prevalence_server <- function(id, patient_data_r, prescriptions_r,
																				 ltcs_r, cases_controls_r = NULL, bnf_level) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Update stratification choices
		observe({
			req(patient_data_r())
			choices <- create_stratification_choices(patient_data_r())
			updateVirtualSelect("ltc_freq_strat_variable", choices = choices)
			updateVirtualSelect("presc_freq_strat_variable", choices = choices)
			if (!is.null(cases_controls_r)) {
				updateVirtualSelect("recent_presc_freq_strat_variable", choices = choices)
			}
		})

		# LTC frequency with optional stratification
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

			calculate_frequency_stats(ltcs, "term")
		})


		# Render LTC frequency table
		output$ltc_freq_table <- renderReactable({
			# Get the underlying data with group information for OR calculation
			ltcs <- ltcs_r()
			patient_data <- patient_data_r()
			strat_var <- input$ltc_freq_strat_variable

			# Apply same stratification filter
			if (strat_var != "") {
				parts <- strsplit(strat_var, "#")[[1]]
				column_name <- parts[1]
				filter_value <- parts[2]
				selected_patids <- patient_data[get(column_name) == filter_value, patid]
				ltcs <- ltcs[patid %in% selected_patids]
			}

			# Calculate with ORs
			table_data <- create_prevalence_ratio_table(
				ltc_freq_df(),
				"term",
				data_with_group = ltcs
			)

			# Update column names to include OR columns
			or_cols <- c("Condition", "Cases (%)", "Controls (%)",
									 "OR", "OR_CI_lower", "OR_CI_upper")

			table_data[, OR_formatted := sprintf("%.3f (%.3f - %.3f)", OR, OR_CI_lower, OR_CI_upper)]

			if (all(c("OR", "p_value", "p_adj") %in% colnames(table_data))) {
				colnames(table_data) <- c(or_cols, "p_value", "p_adj", "OR_formatted")

				# Create reactable with expandable details for p-values
				reactable(
					table_data,
					columns = list(
						Condition = colDef(
							name = "Condition",
							cell = function(value, index) {
								p_adj_val <- table_data[index, p_adj]
								if (!is.na(p_adj_val) && p_adj_val < 0.05) {
									paste0(value, "*")
								} else {
									value
								}
							}
						),
						OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						`Cases (%)` = colDef(format = colFormat(digits = 2)),
						`Controls (%)` = colDef(format = colFormat(digits = 2)),
						OR = colDef(show = FALSE),
						OR_CI_lower = colDef(show = FALSE),
						OR_CI_upper = colDef(show = FALSE),
						p_value = colDef(show = FALSE),
						p_adj = colDef(show = FALSE)
					),
					details = function(index) {
						p_val <- table_data[index, p_value]
						p_adj_val <- table_data[index, p_adj]
						if (is.na(p_val) & is.na(p_adj_val)) return(NULL)
						htmltools::div(
							style = "padding: 16px",
							htmltools::tags$b("Statistical Testing:"),
							htmltools::tags$div(
								style = "margin-top: 8px",
								sprintf("Raw p-value: %.4f", p_val)
							),
							htmltools::tags$div(
								sprintf("Adjusted p-value: %.4f", p_adj_val)
							),
							htmltools::tags$div(
								style = "margin-top: 8px; font-style: italic; color: #666;",
								if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
							)
						)
					},
					showPageInfo = FALSE,
					compact = TRUE,
					showPageSizeOptions = TRUE,
					defaultPageSize = 15
				)
			} else {
				colnames(table_data) <- c(or_cols)
				reactable(table_data[order(-OR)], showPageInfo = FALSE, defaultPageSize = 15)
			}
		})

		# Prescription frequency with optional stratification
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

			calculate_frequency_stats(presc, "substance")
		})

		output$presc_freq_table <- renderReactable({
			# Get the underlying data with group information for OR calculation
			presc <- prescriptions_r()
			patient_data <- patient_data_r()
			strat_var <- input$presc_freq_strat_variable

			# Apply same stratification filter
			if (strat_var != "") {
				parts <- strsplit(strat_var, "#")[[1]]
				column_name <- parts[1]
				filter_value <- parts[2]
				selected_patids <- patient_data[get(column_name) == filter_value, patid]
				presc <- presc[patid %in% selected_patids]
			}

			# Calculate with ORs
			table_data <- create_prevalence_ratio_table(
				presc_freq_df(),
				"substance",
				data_with_group = presc
			)

			or_cols <- c("Substance", "Cases (%)", "Controls (%)",
									 "OR", "OR_CI_lower", "OR_CI_upper")

			table_data[, OR_formatted := sprintf("%.3f (%.3f - %.3f)", OR, OR_CI_lower, OR_CI_upper)]

			if (all(c("OR", "p_value", "p_adj") %in% colnames(table_data))) {
				colnames(table_data) <- c(or_cols, "p_value", "p_adj", "OR_formatted")

				reactable(
					table_data,
					columns = list(
						Substance = colDef(
							name = "Substance",
							cell = function(value, index) {
								p_adj_val <- table_data[index, p_adj]
								if (!is.na(p_adj_val) && p_adj_val < 0.05) {
									paste0(value, "*")
								} else {
									value
								}
							}
						),
						OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						`Cases (%)` = colDef(format = colFormat(digits = 2)),
						`Controls (%)` = colDef(format = colFormat(digits = 2)),
						OR = colDef(show = FALSE),
						OR_CI_lower = colDef(show = FALSE),
						OR_CI_upper = colDef(show = FALSE),
						p_value = colDef(show = FALSE),
						p_adj = colDef(show = FALSE)
					),
					details = function(index) {
						p_val <- table_data[index, p_value]
						p_adj_val <- table_data[index, p_adj]
						if (is.na(p_val) & is.na(p_adj_val)) return(NULL)
						htmltools::div(
							style = "padding: 16px",
							htmltools::tags$b("Statistical Testing:"),
							htmltools::tags$div(
								style = "margin-top: 8px",
								sprintf("Raw p-value: %.4f", p_val)
							),
							htmltools::tags$div(
								sprintf("Adjusted p-value: %.4f", p_adj_val)
							),
							htmltools::tags$div(
								style = "margin-top: 8px; font-style: italic; color: #666;",
								if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
							)
						)
					},
					showPageInfo = FALSE,
					showPageSizeOptions = TRUE,
					compact = TRUE,
					defaultPageSize = 15
				)
			} else {
				colnames(table_data) <- c(or_cols)
				reactable(table_data[order(-OR)], showPageInfo = FALSE, defaultPageSize = 15)
			}
		})

		# LTC dropdown UI
		output$ltc_dropdown_ui <- renderUI({
			req(ltcs_r())
			ltcs <- ltcs_r()

			# Get LTC chapters for grouping
			sub_chapter <- bnf_lookup[, .(BNF_Chemical_Substance, BNF_Section)] |> unique()
			sub_chapter <- sub_chapter[, first(.SD), BNF_Chemical_Substance]

			# For LTCs, we might not have BNF grouping, so just alphabetical
			ltc_terms <- sort(unique(ltcs$term))

			virtualSelectInput(
				ns("ltc_dropdown"),
				label = "Select 1 or more conditions:",
				choices = ltc_terms,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Prescription dropdown UI
		output$presc_dropdown_ui <- renderUI({
			req(prescriptions_r())

			# Get current BNF level from parent module
			current_level <- bnf_level  # You'll need to pass this as a parameter

			prescs <- prescriptions_r()
			unique_substances <- unique(prescs$substance)

			# Determine grouping strategy based on current level
			if (current_level == "BNF_Chapter") {
				# No grouping needed - just show chapters alphabetically
				choices <- sort(unique_substances)
			} else if (current_level == "BNF_Section") {
				# Group by Chapter, show Sections
				lookup_data <- bnf_lookup[BNF_Section %in% unique_substances,
																	.(BNF_Section, BNF_Chapter)] |> unique()
				lookup_data <- lookup_data[order(BNF_Chapter, BNF_Section)]
				choices <- with(lookup_data, split(BNF_Section, BNF_Chapter))
			} else if (current_level == "BNF_Paragraph") {
				# Group by Section, show Paragraphs
				lookup_data <- bnf_lookup[BNF_Paragraph %in% unique_substances,
																	.(BNF_Paragraph, BNF_Section)] |> unique()
				lookup_data <- lookup_data[order(BNF_Section, BNF_Paragraph)]
				choices <- with(lookup_data, split(BNF_Paragraph, BNF_Section))
			} else {  # BNF_Chemical_Substance
				# Group by Paragraph, show Substances
				lookup_data <- bnf_lookup[BNF_Chemical_Substance %in% unique_substances,
																	.(BNF_Chemical_Substance, BNF_Paragraph)] |> unique()
				lookup_data <- lookup_data[order(BNF_Paragraph, BNF_Chemical_Substance)]
				choices <- with(lookup_data, split(BNF_Chemical_Substance, BNF_Paragraph))
			}

			virtualSelectInput(
				ns("presc_dropdown"),
				label = "Select 1 or more substances:",
				choices = choices,
				multiple = TRUE,
				search = TRUE
			)
		})


		# LTC by prescription table
		output$ltc_by_presc <- renderReactable({
			req(ltcs_r(), input$presc_dropdown)

			result <- calculate_prevalence_cca(
				ltcs_r(),
				prescriptions_r(),
				input$presc_dropdown,
				"substance",
				"term"
			)

			validate(need(!is.null(result), "No results found for filter"))

			result[, OR_formatted := sprintf("%.3f (%.3f - %.3f)", OR, OR_CI_lower, OR_CI_upper)]
			browser()
			reactable(
				result,
				columns = list(
					term = colDef(
						name = "Condition",
						cell = function(value, index) {
							p_adj_val <- result[index, p_adj]
							if (!is.na(p_adj_val) && p_adj_val < 0.05) {
								paste0(value, "*")
							} else {
								value
							}
						}
					),
					OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
					case = colDef(name = "Cases (%)", format = colFormat(digits = 2)),
					control = colDef(name = "Controls (%)", format = colFormat(digits = 2)),
					OR = colDef(show = FALSE),
					OR_CI_lower = colDef(show = FALSE),
					OR_CI_upper = colDef(show = FALSE),
					p_value = colDef(show = FALSE),
					p_adj = colDef(show = FALSE)
				),
				details = function(index) {
					p_val <- result[index, p_value]
					p_adj_val <- result[index, p_adj]
					if (is.na(p_val) & is.na(p_adj_val)) return(NULL)
					htmltools::div(
						style = "padding: 16px",
						htmltools::tags$b("Statistical Testing:"),
						htmltools::tags$div(
							style = "margin-top: 8px",
							sprintf("Raw p-value: %.4f", p_val)
						),
						htmltools::tags$div(
							sprintf("Adjusted p-value: %.4f", p_adj_val)
						),
						htmltools::tags$div(
							style = "margin-top: 8px; font-style: italic; color: #666;",
							if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
						)
					)
				},
				showPageInfo = FALSE,
				showPageSizeOptions = TRUE,
				compact = TRUE,
				defaultPageSize = 15
			)
		})

		# Prescription by LTC table
		output$presc_by_ltc <- renderReactable({
			req(prescriptions_r(), input$ltc_dropdown)

			result <- calculate_prevalence_cca(
				prescriptions_r(),
				ltcs_r(),
				input$ltc_dropdown,
				"term",
				"substance"
			)

			validate(need(!is.null(result), "No results found for filter"))

			result[, OR_formatted := sprintf("%.3f (%.3f - %.3f)", OR, OR_CI_lower, OR_CI_upper)]

			reactable(
				result,
				columns = list(
					substance = colDef(
						name = "Substance",
						cell = function(value, index) {
							p_adj_val <- result[index, p_adj]
							if (!is.na(p_adj_val) && p_adj_val < 0.05) {
								paste0(value, "*")
							} else {
								value
							}
						}
					),
					OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
					case = colDef(name = "Cases (%)", format = colFormat(digits = 2)),
					control = colDef(name = "Controls (%)", format = colFormat(digits = 2)),
					OR = colDef(show = FALSE),
					OR_CI_lower = colDef(show = FALSE),
					OR_CI_upper = colDef(show = FALSE),
					p_value = colDef(show = FALSE),
					p_adj = colDef(show = FALSE)
				),
				details = function(index) {
					p_val <- result[index, p_value]
					p_adj_val <- result[index, p_adj]
					if (is.na(p_val) & is.na(p_adj_val)) return(NULL)
					htmltools::div(
						style = "padding: 16px",
						htmltools::tags$b("Statistical Testing:"),
						htmltools::tags$div(
							style = "margin-top: 8px",
							sprintf("Raw p-value: %.4f", p_val)
						),
						htmltools::tags$div(
							sprintf("Adjusted p-value: %.4f", p_adj_val)
						),
						htmltools::tags$div(
							style = "margin-top: 8px; font-style: italic; color: #666;",
							if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
						)
					)
				},
				showPageInfo = FALSE,
				defaultPageSize = 15
			)
		})

		# ========== NEW SECTION: Recent prescriptions (from case_controls_r) ==========

		# Recent prescription frequency with optional stratification
		recent_presc_freq_df <- reactive({
			req(cases_controls_r())

			# cases_controls_r has: patid, index_date, substance, group
			# We need to add strata for calculate_frequency_stats
			recent_presc_with_strata <- merge(
				cases_controls_r(),
				patient_data_r()[, .(patid, strata)],
				by = "patid"
			)

			strat_var <- input$recent_presc_freq_strat_variable

			if (strat_var != "") {
				parts <- strsplit(strat_var, "#")[[1]]
				column_name <- parts[1]
				filter_value <- parts[2]
				selected_patids <- patient_data_r()[get(column_name) == filter_value, patid]
				recent_presc_with_strata <- recent_presc_with_strata[patid %in% selected_patids]
			}

			calculate_frequency_stats(recent_presc_with_strata, "substance")
		})

		output$recent_presc_freq_table <- renderReactable({
			req(cases_controls_r())

			# Get data with strata for OR calculation
			recent_presc_with_strata <- merge(
				cases_controls_r(),
				patient_data_r()[, .(patid, strata)],
				by = "patid"
			)

			strat_var <- input$recent_presc_freq_strat_variable

			# Apply same stratification filter
			if (strat_var != "") {
				parts <- strsplit(strat_var, "#")[[1]]
				column_name <- parts[1]
				filter_value <- parts[2]
				selected_patids <- patient_data_r()[get(column_name) == filter_value, patid]
				recent_presc_with_strata <- recent_presc_with_strata[patid %in% selected_patids]
			}

			# Calculate with ORs
			table_data <- create_prevalence_ratio_table(
				recent_presc_freq_df(),
				"substance",
				data_with_group = recent_presc_with_strata
			)

			or_cols <- c("Substance", "Cases (%)", "Controls (%)",
									 "OR", "OR_CI_lower", "OR_CI_upper")

			table_data[, OR_formatted := sprintf("%.3f (%.3f - %.3f)", OR, OR_CI_lower, OR_CI_upper)]

			if (all(c("OR", "p_value", "p_adj") %in% colnames(table_data))) {
				colnames(table_data) <- c(or_cols, "p_value", "p_adj", "OR_formatted")

				reactable(
					table_data,
					columns = list(
						Substance = colDef(
							name = "Substance",
							cell = function(value, index) {
								p_adj_val <- table_data[index, p_adj]
								if (!is.na(p_adj_val) && p_adj_val < 0.05) {
									paste0(value, "*")
								} else {
									value
								}
							}
						),
						OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						`Cases (%)` = colDef(format = colFormat(digits = 2)),
						`Controls (%)` = colDef(format = colFormat(digits = 2)),
						OR = colDef(show = FALSE),
						OR_CI_lower = colDef(show = FALSE),
						OR_CI_upper = colDef(show = FALSE),
						p_value = colDef(show = FALSE),
						p_adj = colDef(show = FALSE)
					),
					details = function(index) {
						p_val <- table_data[index, p_value]
						p_adj_val <- table_data[index, p_adj]
						if (is.na(p_val) & is.na(p_adj_val)) return(NULL)
						htmltools::div(
							style = "padding: 16px",
							htmltools::tags$b("Statistical Testing:"),
							htmltools::tags$div(
								style = "margin-top: 8px",
								sprintf("Raw p-value: %.4f", p_val)
							),
							htmltools::tags$div(
								sprintf("Adjusted p-value: %.4f", p_adj_val)
							),
							htmltools::tags$div(
								style = "margin-top: 8px; font-style: italic; color: #666;",
								if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
							)
						)
					},
					showPageInfo = FALSE,
					showPageSizeOptions = TRUE,
					searchable = TRUE,
					compact = TRUE,
					defaultPageSize = 15
				)
			} else {
				colnames(table_data) <- c(or_cols)
				reactable(table_data[order(-OR)], showPageInfo = FALSE, defaultPageSize = 15)
			}
		})

		# LTC dropdown for recent prescriptions filter
		output$recent_ltc_dropdown_ui <- renderUI({
			req(ltcs_r())
			ltc_terms <- sort(unique(ltcs_r()$term))

			virtualSelectInput(
				ns("recent_ltc_dropdown"),
				label = "Filter by conditions:",
				choices = ltc_terms,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Background prescription dropdown for recent prescriptions filter
		output$recent_bg_presc_dropdown_ui <- renderUI({
			req(prescriptions_r())

			current_level <- bnf_level

			prescs <- prescriptions_r()
			unique_substances <- unique(prescs$substance)

			if (current_level == "BNF_Chapter") {
				choices <- sort(unique_substances)
			} else if (current_level == "BNF_Section") {
				lookup_data <- bnf_lookup[BNF_Section %in% unique_substances,
																	.(BNF_Section, BNF_Chapter)] |> unique()
				lookup_data <- lookup_data[order(BNF_Chapter, BNF_Section)]
				choices <- with(lookup_data, split(BNF_Section, BNF_Chapter))
			} else if (current_level == "BNF_Paragraph") {
				lookup_data <- bnf_lookup[BNF_Paragraph %in% unique_substances,
																	.(BNF_Paragraph, BNF_Section)] |> unique()
				lookup_data <- lookup_data[order(BNF_Section, BNF_Paragraph)]
				choices <- with(lookup_data, split(BNF_Paragraph, BNF_Section))
			} else {
				lookup_data <- bnf_lookup[BNF_Chemical_Substance %in% unique_substances,
																	.(BNF_Chemical_Substance, BNF_Paragraph)] |> unique()
				lookup_data <- lookup_data[order(BNF_Paragraph, BNF_Chemical_Substance)]
				choices <- with(lookup_data, split(BNF_Chemical_Substance, BNF_Paragraph))
			}

			virtualSelectInput(
				ns("recent_bg_presc_dropdown"),
				label = "Filter by background medications:",
				choices = choices,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Recent prescriptions by BOTH LTC and background medication filter
		output$recent_by_ltc_and_presc <- renderReactable({
			req(cases_controls_r())

			# Need at least one filter selected
			has_ltc_filter <- !is.null(input$recent_ltc_dropdown) && length(input$recent_ltc_dropdown) > 0
			has_presc_filter <- !is.null(input$recent_bg_presc_dropdown) && length(input$recent_bg_presc_dropdown) > 0

			validate(need(has_ltc_filter || has_presc_filter,
										"Please select at least one condition or background medication"))

			# Start with all patients
			filtered_patids <- unique(patient_data_r()$patid)

			# Filter by LTCs if selected
			if (has_ltc_filter) {
				ltc_patids <- unique(ltcs_r()[term %in% input$recent_ltc_dropdown, patid])
				filtered_patids <- bit64::as.integer64(intersect(as.character(filtered_patids), as.character(ltc_patids)))
			}

			# Filter by background prescriptions if selected
			if (has_presc_filter) {
				presc_patids <- unique(prescriptions_r()[substance %in% input$recent_bg_presc_dropdown, patid])
				filtered_patids <- bit64::as.integer64(intersect(as.character(filtered_patids), as.character(presc_patids)))
			}

			validate(need(length(filtered_patids) > 0,
										"No patients match the selected filters"))

			# Filter recent prescriptions to these patients and add strata
			recent_presc_filtered <- cases_controls_r()[patid %in% filtered_patids]
			recent_presc_filtered <- merge(
				recent_presc_filtered,
				patient_data_r()[, .(patid, strata)],
				by = "patid"
			)

			# Calculate frequencies
			freq <- calculate_frequency_stats(recent_presc_filtered, "substance")

			# Calculate prevalence table with ORs
			result <- create_prevalence_ratio_table(
				freq,
				"substance",
				min_pct = 1,
				data_with_group = recent_presc_filtered
			)

			validate(need(!is.null(result) && nrow(result) > 0,
										"No results found for the selected filters"))

			or_cols <- c("Substance", "Cases (%)", "Controls (%)",
									 "OR", "OR_CI_lower", "OR_CI_upper")

			result[, OR_formatted := sprintf("%.3f (%.3f - %.3f)", OR, OR_CI_lower, OR_CI_upper)]

			if (all(c("OR", "p_value", "p_adj") %in% colnames(result))) {
				colnames(result) <- c(or_cols, "p_value", "p_adj", "OR_formatted")

				reactable(
					result,
					columns = list(
						Substance = colDef(
							name = "Substance",
							cell = function(value, index) {
								p_adj_val <- result[index, p_adj]
								if (!is.na(p_adj_val) && p_adj_val < 0.05) {
									paste0(value, "*")
								} else {
									value
								}
							}
						),
						OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
						`Cases (%)` = colDef(format = colFormat(digits = 2)),
						`Controls (%)` = colDef(format = colFormat(digits = 2)),
						OR = colDef(show = FALSE),
						OR_CI_lower = colDef(show = FALSE),
						OR_CI_upper = colDef(show = FALSE),
						p_value = colDef(show = FALSE),
						p_adj = colDef(show = FALSE)
					),
					details = function(index) {
						p_val <- result[index, p_value]
						p_adj_val <- result[index, p_adj]
						if (is.na(p_val) & is.na(p_adj_val)) return(NULL)
						htmltools::div(
							style = "padding: 16px",
							htmltools::tags$b("Statistical Testing:"),
							htmltools::tags$div(
								style = "margin-top: 8px",
								sprintf("Raw p-value: %.4f", p_val)
							),
							htmltools::tags$div(
								sprintf("Adjusted p-value: %.4f", p_adj_val)
							),
							htmltools::tags$div(
								style = "margin-top: 8px; font-style: italic; color: #666;",
								if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
							)
						)
					},
					showPageInfo = FALSE,
					searchable = TRUE,
					showPageSizeOptions = TRUE,
					defaultPageSize = 15,
					compact = TRUE
				)
			} else {
				colnames(result) <- c(or_cols)
				setorder(result, -OR)
				reactable(result, showPageInfo = FALSE, defaultPageSize = 15)
			}
		})

	})
}