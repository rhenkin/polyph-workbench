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
						hideClearButton = FALSE,
						disableOptionGroupCheckbox = TRUE,
						multiple = FALSE,
						dropboxWrapper = "body"
					),
					downloadButton(ns("download_ltc_freq"), "Download Table", class = "btn-sm btn-secondary"),
					br(), br(),
					reactableOutput(ns("ltc_freq_table"), height = "50em")
				),
				nav_panel(
					title = "Prevalence by prescription filter",
					card(
						full_screen = TRUE,
						height = "60em",
						div("Select medications to show the prevalence difference only among the patients taking those medications."),
						div("Table contains conditions with a minimum of 1% prevalence in both cases and controls"),
						fluidRow(
							column(6,
										 div(strong("Filter by background (continuous) prescriptions:")),
										 uiOutput(ns("presc_dropdown_ui"))
							),
							column(6,
										 div(strong("Filter by recent prescriptions (optional):")),
										 uiOutput(ns("ltc_by_presc_recent_dropdown_ui"))
							)
						),
						downloadButton(ns("download_ltc_by_presc"), "Download Table", class = "btn-sm btn-secondary"),
						br(), br(),
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
						disableOptionGroupCheckbox = TRUE,
						multiple = FALSE,
						hideClearButton = FALSE,
						dropboxWrapper = "body"
					),
					downloadButton(ns("download_presc_freq"), "Download Table", class = "btn-sm btn-secondary"),
					br(), br(),
					reactableOutput(ns("presc_freq_table"), height = "50em")
				),
				nav_panel(
					title = "Prevalence by condition filter",
					card(
						full_screen = TRUE,
						height = "60em",
						div("Select conditions to show the prevalence difference only among the patients that were diagnosed with those conditions."),
						div("Table contains medications with a minimum of 1% prevalence in both cases and controls"),
						fluidRow(
							column(6,
										 div(strong("Filter by LTCs:")),
										 uiOutput(ns("ltc_dropdown_ui"))
							),
							column(6,
										 div(strong("Filter by recent prescriptions (optional):")),
										 uiOutput(ns("presc_by_ltc_recent_dropdown_ui"))
							)
						),
						downloadButton(ns("download_presc_by_ltc"), "Download Table", class = "btn-sm btn-secondary"),
						br(), br(),
						reactableOutput(ns("presc_by_ltc"), height = "50em")
					)
				)
			)
		),
		accordion_panel(
			title = "Recent prescriptions",
			value = "recent_presc_prev_tables",
			icon = bs_icon("prescription2"),
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
						disableOptionGroupCheckbox = TRUE,
						hideClearButton = FALSE,
						dropboxWrapper = "body"
					),
					downloadButton(ns("download_recent_presc_freq"), "Download Table", class = "btn-sm btn-secondary"),
					br(), br(),
					reactableOutput(ns("recent_presc_freq_table"), height = "50em")
				),
				nav_panel(
					title = "Prevalence by condition and background medication filter",
					card(
						full_screen = TRUE,
						height = "60em",
						div("Select conditions and background medications to show the prevalence difference only among patients with those characteristics."),
						div("Table contains medications with a minimum of 1% prevalence in both cases and controls"),
						fluidRow(
							column(6,
										 uiOutput(ns("recent_ltc_dropdown_ui"))
							),
							column(6,
										 uiOutput(ns("recent_bg_presc_dropdown_ui"))
							)
						),
						downloadButton(ns("download_recent_by_filters"), "Download Table", class = "btn-sm btn-secondary"),
						br(), br(),
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

		# Recent prescription dropdown for LTC by prescription table
		output$ltc_by_presc_recent_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(cases_controls_r())

			unique_substances <- sort(unique(cases_controls_r()$substance))

			virtualSelectInput(
				ns("ltc_by_presc_recent_dropdown"),
				label = "Select recent prescriptions:",
				choices = unique_substances,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Recent prescription dropdown for prescription by LTC table
		output$presc_by_ltc_recent_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(cases_controls_r())

			unique_substances <- sort(unique(cases_controls_r()$substance))

			virtualSelectInput(
				ns("presc_by_ltc_recent_dropdown"),
				label = "Select recent prescriptions:",
				choices = unique_substances,
				multiple = TRUE,
				search = TRUE
			)
		})

		# ============================================================================
		# LTC PREVALENCE TABLES
		# ============================================================================

		# Store rendered table data
		ltc_freq_table_data <- reactiveVal(NULL)
		ltc_by_presc_table_data <- reactiveVal(NULL)
		presc_freq_table_data <- reactiveVal(NULL)
		presc_by_ltc_table_data <- reactiveVal(NULL)
		recent_presc_freq_table_data <- reactiveVal(NULL)
		recent_by_filters_table_data <- reactiveVal(NULL)

		# LTC frequency with optional stratification
		ltc_freq_df <- reactive({
			req(ltcs_r())
			ltcs_filtered <- apply_patient_stratification(
				ltcs_r(),
				input$ltc_freq_strat_variable,
				patient_data_r()
			)
			calculate_frequency_stats(ltcs_filtered, "term")
		})

		# Render LTC frequency table
		output$ltc_freq_table <- renderReactable({
			req(ltcs_r())

			# Apply stratification and calculate with ORs
			ltcs_filtered <- apply_patient_stratification(
				ltcs_r(),
				input$ltc_freq_strat_variable,
				patient_data_r()
			)

			result_table <- create_prevalence_ratio_table(
				ltc_freq_df(),
				"term",
				data_with_group = ltcs_filtered
			)

			ltc_freq_table_data(result_table)

			# Add formatted OR column
			result_table <- add_or_formatted_column(result_table)

			# Rename columns with standard names
			col_names <- get_prevalence_column_names("LTC")
			colnames(result_table) <- c(col_names, "p_value", "p_adj", "OR_formatted")

			# Build column definitions
			columns <- list(
				LTC = colDef(
					name = "LTC",
					minWidth = 200,
					cell = create_significance_cell(result_table)
				),
				`Cases (%)` = colDef(format = colFormat(digits = 2)),
				`Controls (%)` = colDef(format = colFormat(digits = 2))
			)
			columns <- c(columns, get_or_reactable_columns())

			# Render with standard config
			do.call(reactable, c(
				list(
					data = result_table,
					columns = columns,
					details = function(index) create_pvalue_details(index, result_table)
				),
				get_standard_reactable_config()
			))
		})

		# ============================================================================
		# LTC BY PRESCRIPTION FILTER
		# ============================================================================

		# Dynamic prescription dropdown
		output$presc_dropdown_ui <- renderUI({
			req(prescriptions_r())
			current_level <- bnf_level
			unique_substances <- unique(prescriptions_r()$substance)

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
			} else {  # BNF_Chemical_Substance
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
				search = TRUE,
				showSelectedOptionsFirst = TRUE
			)
		})

		# LTC by prescription table

		output$ltc_by_presc <- renderReactable({
			req(ltcs_r())

			# Need at least one filter selected
			req(isTruthy(input$presc_dropdown) | isTruthy(input$ltc_by_presc_recent_dropdown))

			# Start with all patients
			filtered_patids <- unique(ltcs_r()$patid)

			# Filter by background prescriptions if selected
			if (isTruthy(input$presc_dropdown)) {
				bg_presc_patids <- unique(prescriptions_r()[substance %in% input$presc_dropdown, patid])
				filtered_patids <- bit64::as.integer64(intersect(as.character(filtered_patids), as.character(bg_presc_patids)))
			}

			# Filter by recent prescriptions if selected
			if (!is.null(cases_controls_r) && isTruthy(input$ltc_by_presc_recent_dropdown)) {
				recent_presc_patids <- unique(cases_controls_r()[substance %in% input$ltc_by_presc_recent_dropdown, patid])
				filtered_patids <- bit64::as.integer64(intersect(as.character(filtered_patids), as.character(recent_presc_patids)))
			}

			# Filter LTCs to these patients
			ltcs_filtered <- ltcs_r()[patid %in% filtered_patids]

			# Calculate frequency stats (this function needs group, strata, and term columns)
			freq_stats <- calculate_frequency_stats(ltcs_filtered, "term")

			# Create prevalence ratio table
			result_table <- create_prevalence_ratio_table(
				freq_stats,
				"term",
				data_with_group = ltcs_filtered
			)

			validate(need(nrow(result_table) > 0, "No results found for filter"))

			ltc_by_presc_table_data(result_table)

			# Add formatted OR column
			result_table <- add_or_formatted_column(result_table)

			# Rename columns with standard names
			col_names <- get_prevalence_column_names("LTC")
			colnames(result_table) <- c(col_names, "p_value", "p_adj", "OR_formatted")

			# Build column definitions
			columns <- list(
				LTC = colDef(
					name = "LTC",
					minWidth = 200,
					cell = create_significance_cell(result_table)
				),
				`Cases (%)` = colDef(format = colFormat(digits = 2)),
				`Controls (%)` = colDef(format = colFormat(digits = 2))
			)
			columns <- c(columns, get_or_reactable_columns())

			# Render with standard config
			do.call(reactable, c(
				list(
					data = result_table,
					columns = columns,
					details = function(index) create_pvalue_details(index, result_table)
				),
				get_standard_reactable_config()
			))
		})

		# ============================================================================
		# PRESCRIPTION PREVALENCE TABLES
		# ============================================================================

		# Prescription frequency with optional stratification
		presc_freq_df <- reactive({
			req(prescriptions_r())
			presc_filtered <- apply_patient_stratification(
				prescriptions_r(),
				input$presc_freq_strat_variable,
				patient_data_r()
			)
			calculate_frequency_stats(presc_filtered, "substance")
		})

		# Render prescription frequency table
		output$presc_freq_table <- renderReactable({
			req(prescriptions_r())

			# Apply stratification and calculate with ORs
			presc_filtered <- apply_patient_stratification(
				prescriptions_r(),
				input$presc_freq_strat_variable,
				patient_data_r()
			)

			result_table <- create_prevalence_ratio_table(
				presc_freq_df(),
				"substance",
				data_with_group = presc_filtered
			)

			presc_freq_table_data(result_table)

			# Add formatted OR column
			result_table <- add_or_formatted_column(result_table)

			# Rename columns with standard names
			col_names <- get_prevalence_column_names("Substance")
			colnames(result_table) <- c(col_names, "p_value", "p_adj", "OR_formatted")

			# Build column definitions
			columns <- list(
				Substance = colDef(
					name = "Substance",
					minWidth = 200,
					cell = create_significance_cell(result_table)
				),
				`Cases (%)` = colDef(format = colFormat(digits = 2)),
				`Controls (%)` = colDef(format = colFormat(digits = 2))
			)
			columns <- c(columns, get_or_reactable_columns())

			# Render with standard config
			do.call(reactable, c(
				list(
					data = result_table,
					columns = columns,
					details = function(index) create_pvalue_details(index, result_table)
				),
				get_standard_reactable_config()
			))
		})

		# ============================================================================
		# PRESCRIPTION BY LTC FILTER
		# ============================================================================

		# Dynamic LTC dropdown
		output$ltc_dropdown_ui <- renderUI({
			req(ltcs_r())
			virtualSelectInput(
				ns("ltc_dropdown"),
				label = "Select 1 or more conditions:",
				choices = sort(unique(ltcs_r()$term)),
				multiple = TRUE,
				search = TRUE
			)
		})

		# Prescription by LTC table
		output$presc_by_ltc <- renderReactable({
			req(prescriptions_r())

			# Need at least one filter selected
			req(isTruthy(input$ltc_dropdown) | isTruthy(input$presc_by_ltc_recent_dropdown))

			# Start with all patients
			filtered_patids <- unique(prescriptions_r()$patid)

			# Filter by LTCs if selected
			if (isTruthy(input$ltc_dropdown)) {
				ltc_patids <- unique(ltcs_r()[term %in% input$ltc_dropdown, patid])
				filtered_patids <- bit64::as.integer64(intersect(as.character(filtered_patids), as.character(ltc_patids)))
			}

			# Filter by recent prescriptions if selected
			if (!is.null(cases_controls_r) && isTruthy(input$presc_by_ltc_recent_dropdown)) {
				recent_presc_patids <- unique(cases_controls_r()[substance %in% input$presc_by_ltc_recent_dropdown, patid])
				filtered_patids <- bit64::as.integer64(intersect(as.character(filtered_patids), as.character(recent_presc_patids)))
			}

			# Filter prescriptions to these patients
			presc_filtered <- prescriptions_r()[patid %in% filtered_patids]

			# Calculate frequency stats (this function needs group, strata, and substance columns)
			freq_stats <- calculate_frequency_stats(presc_filtered, "substance")

			# Create prevalence ratio table
			result_table <- create_prevalence_ratio_table(
				freq_stats,
				"substance",
				data_with_group = presc_filtered
			)

			validate(need(nrow(result_table) > 0, "No results found for filter"))

			presc_by_ltc_table_data(result_table)

			# Add formatted OR column
			result_table <- add_or_formatted_column(result_table)

			# Rename columns with standard names
			col_names <- get_prevalence_column_names("Substance")
			colnames(result_table) <- c(col_names, "p_value", "p_adj", "OR_formatted")

			# Build column definitions
			columns <- list(
				Substance = colDef(
					name = "Substance",
					minWidth = 200,
					cell = create_significance_cell(result_table)
				),
				`Cases (%)` = colDef(format = colFormat(digits = 2)),
				`Controls (%)` = colDef(format = colFormat(digits = 2))
			)
			columns <- c(columns, get_or_reactable_columns())

			# Render with standard config
			do.call(reactable, c(
				list(
					data = result_table,
					columns = columns,
					details = function(index) create_pvalue_details(index, result_table)
				),
				get_standard_reactable_config()
			))
		})
		# ============================================================================
		# RECENT PRESCRIPTIONS (from cases_controls_r)
		# ============================================================================

		if (!is.null(cases_controls_r)) {

			# Recent prescription frequency with optional stratification
			recent_presc_freq_df <- reactive({
				req(cases_controls_r())

				# Add strata for calculate_frequency_stats
				recent_presc_with_strata <- merge(
					cases_controls_r(),
					patient_data_r()[, .(patid, strata)],
					by = "patid"
				)

				recent_presc_filtered <- apply_patient_stratification(
					recent_presc_with_strata,
					input$recent_presc_freq_strat_variable,
					patient_data_r()
				)

				calculate_frequency_stats(recent_presc_filtered, "substance")
			})

			# Render recent prescription frequency table
			output$recent_presc_freq_table <- renderReactable({
				req(cases_controls_r())

				# Get data with strata for OR calculation
				recent_presc_with_strata <- merge(
					cases_controls_r(),
					patient_data_r()[, .(patid, strata)],
					by = "patid"
				)

				recent_presc_filtered <- apply_patient_stratification(
					recent_presc_with_strata,
					input$recent_presc_freq_strat_variable,
					patient_data_r()
				)

				result_table <- create_prevalence_ratio_table(
					recent_presc_freq_df(),
					"substance",
					data_with_group = recent_presc_filtered
				)

				recent_presc_freq_table_data(result_table)

				# Add formatted OR column
				result_table <- add_or_formatted_column(result_table)

				# Rename columns with standard names
				col_names <- get_prevalence_column_names("Substance")
				colnames(result_table) <- c(col_names, "p_value", "p_adj", "OR_formatted")

				# Build column definitions
				columns <- list(
					Substance = colDef(
						name = "Substance",
						minWidth = 200,
						cell = create_significance_cell(result_table)
					),
					`Cases (%)` = colDef(format = colFormat(digits = 2)),
					`Controls (%)` = colDef(format = colFormat(digits = 2))
				)
				columns <- c(columns, get_or_reactable_columns())

				# Render with standard config
				do.call(reactable, c(
					list(
						data = result_table,
						columns = columns,
						details = function(index) create_pvalue_details(index, result_table)
					),
					get_standard_reactable_config()
				))
			})

			# Dynamic LTC and background prescription dropdowns for filtering
			output$recent_ltc_dropdown_ui <- renderUI({
				req(ltcs_r())
				virtualSelectInput(
					ns("recent_ltc_dropdown"),
					label = "Filter by LTCs:",
					choices = sort(unique(ltcs_r()$term)),
					multiple = TRUE,
					search = TRUE
				)
			})

			output$recent_bg_presc_dropdown_ui <- renderUI({
				req(prescriptions_r())

				unique_substances <- sort(unique(prescriptions_r()$substance))

				if (bnf_level == "BNF_Chapter") {
					choices <- sort(unique_substances)
				} else if (bnf_level == "BNF_Section") {
					lookup_data <- bnf_lookup[BNF_Section %in% unique_substances,
																		.(BNF_Section, BNF_Chapter)] |> unique()
					lookup_data <- lookup_data[order(BNF_Chapter, BNF_Section)]
					choices <- with(lookup_data, split(BNF_Section, BNF_Chapter))
				} else if (bnf_level == "BNF_Paragraph") {
					lookup_data <- bnf_lookup[BNF_Paragraph %in% unique_substances,
																		.(BNF_Paragraph, BNF_Section)] |> unique()
					lookup_data <- lookup_data[order(BNF_Section, BNF_Paragraph)]
					choices <- with(lookup_data, split(BNF_Paragraph, BNF_Section))
				} else {  # BNF_Chemical_Substance
					lookup_data <- bnf_lookup[BNF_Chemical_Substance %in% unique_substances,
																		.(BNF_Chemical_Substance, BNF_Paragraph)] |> unique()
					lookup_data <- lookup_data[order(BNF_Paragraph, BNF_Chemical_Substance)]
					choices <- with(lookup_data, split(BNF_Chemical_Substance, BNF_Paragraph))
				}

				virtualSelectInput(
					ns("recent_bg_presc_dropdown"),
					label = "Filter by background prescriptions:",
					choices = choices,
					multiple = TRUE,
					search = TRUE
				)
			})

			# Recent prescriptions by LTC and background prescription filter
			output$recent_by_ltc_and_presc <- renderReactable({
				req(cases_controls_r())
				req(isTruthy(input$recent_ltc_dropdown) | isTruthy(input$recent_bg_presc_dropdown))

				# Filter patients by LTC
				if (isTruthy(input$recent_ltc_dropdown)) {
					ltc_patids <- unique(ltcs_r()[term %in% input$recent_ltc_dropdown, patid])
				} else ltc_patids <- unique(ltcs_r()$patid)


				# Filter patients by background prescription
				if (isTruthy(input$recent_bg_presc_dropdown)) {
					bg_presc_patids <- unique(prescriptions_r()[substance %in% input$recent_bg_presc_dropdown, patid])
				} else bg_presc_patids <- unique(prescriptions_r()$patid)


				# Get intersection
				filtered_patids <- bit64::as.integer64(intersect(as.character(ltc_patids), as.character(bg_presc_patids)))

				# Filter recent prescriptions to these patients
				recent_presc_filtered <- merge(
					cases_controls_r()[patid %in% filtered_patids],
					patient_data_r()[, .(patid, strata)],
					by = "patid"
				)

				# Calculate frequencies
				freq <- calculate_frequency_stats(recent_presc_filtered, "substance")

				# Calculate prevalence table with ORs
				result_table <- create_prevalence_ratio_table(
					freq,
					"substance",
					min_pct = 1,
					data_with_group = recent_presc_filtered
				)

				validate(need(!is.null(result_table) && nrow(result_table) > 0,
											"No results found for the selected filters"))

				recent_by_filters_table_data(result_table)

				# Add formatted OR column
				result_table <- add_or_formatted_column(result_table)

				# Rename columns with standard names
				col_names <- get_prevalence_column_names("Substance")
				colnames(result_table) <- c(col_names, "p_value", "p_adj", "OR_formatted")

				# Build column definitions
				columns <- list(
					Substance = colDef(
						name = "Substance",
						minWidth = 200,
						cell = create_significance_cell(result_table)
					),
					`Cases (%)` = colDef(format = colFormat(digits = 2)),
					`Controls (%)` = colDef(format = colFormat(digits = 2))
				)
				columns <- c(columns, get_or_reactable_columns())

				# Render with standard config
				do.call(reactable, c(
					list(
						data = result_table,
						columns = columns,
						details = function(index) create_pvalue_details(index, result_table)
					),
					get_standard_reactable_config()
				))
			})
		}

		# ============================================================================
		# DOWNLOAD HANDLERS
		# ============================================================================

		output$download_ltc_freq <- downloadHandler(
			filename = function() {
				strat_suffix <- if (!is.null(input$ltc_freq_strat_variable) &&
														input$ltc_freq_strat_variable != "all") {
					paste0("_", gsub(" ", "_", input$ltc_freq_strat_variable))
				} else {
					""
				}
				paste0("ltc_prevalence", strat_suffix, "_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(ltc_freq_table_data())
				fwrite(ltc_freq_table_data(), file)
			}
		)

		output$download_ltc_by_presc <- downloadHandler(
			filename = function() {
				paste0("ltc_by_prescription_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(ltc_by_presc_table_data())
				fwrite(ltc_by_presc_table_data(), file)
			}
		)

		output$download_presc_freq <- downloadHandler(
			filename = function() {
				strat_suffix <- if (!is.null(input$presc_freq_strat_variable) &&
														input$presc_freq_strat_variable != "all") {
					paste0("_", gsub(" ", "_", input$presc_freq_strat_variable))
				} else {
					""
				}
				paste0("prescription_prevalence", strat_suffix, "_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(presc_freq_table_data())
				fwrite(presc_freq_table_data(), file)
			}
		)

		output$download_presc_by_ltc <- downloadHandler(
			filename = function() {
				paste0("prescription_by_ltc_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(presc_by_ltc_table_data())
				fwrite(presc_by_ltc_table_data(), file)
			}
		)

		if (!is.null(cases_controls_r)) {
			output$download_recent_presc_freq <- downloadHandler(
				filename = function() {
					strat_suffix <- if (!is.null(input$recent_presc_freq_strat_variable) &&
															input$recent_presc_freq_strat_variable != "all") {
						paste0("_", gsub(" ", "_", input$recent_presc_freq_strat_variable))
					} else {
						""
					}
					paste0("recent_prescription_prevalence", strat_suffix, "_", Sys.Date(), ".csv")
				},
				content = function(file) {
					req(recent_presc_freq_table_data())
					fwrite(recent_presc_freq_table_data(), file)
				}
			)

			output$download_recent_by_filters <- downloadHandler(
				filename = function() {
					paste0("recent_prescription_by_filters_", Sys.Date(), ".csv")
				},
				content = function(file) {
					req(recent_by_filters_table_data())
					fwrite(recent_by_filters_table_data(), file)
				}
			)
		}
	})
}