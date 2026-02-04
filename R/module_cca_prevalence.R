
module_cca_prevalence_ui <- function(id) {
	ns <- NS(id)

	accordion(
		open = FALSE,
		accordion_panel(
			title = "Long-term conditions",
			value = "ltc_prev_tables",
			icon = bs_icon("heart-pulse"),
			card(
				full_screen = TRUE,
				height = "60em",
				p("Click on a column name to reorder table"),
				div("Table contains conditions with a minimum of 1% prevalence in both cases and controls"),

				# Filter section
				fluidRow(
					column(3,
								 div(strong("Filter by patient subgroup (optional):")),
								 virtualSelectInput(
								 	ns("ltc_freq_strat_variable"),
								 	label = NULL,
								 	choices = NULL,
								 	autoSelectFirstOption = FALSE,
								 	search = FALSE,
								 	hideClearButton = FALSE,
								 	disableOptionGroupCheckbox = TRUE,
								 	multiple = FALSE,
								 	dropboxWrapper = "body"
								 )
					),
					column(3,
								 div(strong("Filter by LTCs (optional):")),
								 uiOutput(ns("ltc_filter_dropdown_ui"))
					),
					column(3,
								 # FIXED: Changed from "background prescriptions" to be consistent
								 div(strong("Filter by background prescriptions (optional):")),
								 uiOutput(ns("presc_dropdown_ui"))
					),
					column(3,
								 div(strong("Filter by recent prescriptions (optional):")),
								 uiOutput(ns("ltc_by_presc_recent_dropdown_ui"))
					)
				),

				# Patient count display
				uiOutput(ns("ltc_patient_count")),
				br(),

				downloadButton(ns("download_ltc_freq"), "Download Table", class = "btn-sm btn-secondary"),
				br(), br(),
				reactableOutput(ns("ltc_freq_table"), height = "50em")
			)
		),
		accordion_panel(
			title = "Prescriptions",
			value = "presc_prev_tables",
			icon = bs_icon("capsule"),
			card(
				full_screen = TRUE,
				height = "60em",
				p("Click on a column name to reorder table"),
				div("Table contains medications with a minimum of 1% prevalence in both cases and controls"),

				# Filter section
				fluidRow(
					column(3,
								 div(strong("Filter by patient subgroup (optional):")),
								 virtualSelectInput(
								 	ns("presc_freq_strat_variable"),
								 	label = NULL,
								 	choices = NULL,
								 	autoSelectFirstOption = FALSE,
								 	search = FALSE,
								 	disableOptionGroupCheckbox = TRUE,
								 	multiple = FALSE,
								 	hideClearButton = FALSE,
								 	dropboxWrapper = "body"
								 )
					),
					column(3,
								 # FIXED: Changed from "prescriptions" to "background prescriptions"
								 div(strong("Filter by background prescriptions (optional):")),
								 uiOutput(ns("presc_filter_dropdown_ui"))
					),
					column(3,
								 div(strong("Filter by LTCs (optional):")),
								 uiOutput(ns("ltc_dropdown_ui"))
					),
					column(3,
								 div(strong("Filter by recent prescriptions (optional):")),
								 uiOutput(ns("presc_by_ltc_recent_dropdown_ui"))
					)
				),

				# Patient count display
				uiOutput(ns("presc_patient_count")),
				br(),

				downloadButton(ns("download_presc_freq"), "Download Table", class = "btn-sm btn-secondary"),
				br(), br(),
				reactableOutput(ns("presc_freq_table"), height = "50em")
			)
		),
		accordion_panel(
			title = "Recent prescriptions",
			value = "recent_presc_prev_tables",
			icon = bs_icon("prescription2"),
			card(
				full_screen = TRUE,
				height = "60em",
				p("Click on a column name to reorder table"),
				div("Table contains medications with a minimum of 1% prevalence in both cases and controls"),

				# Filter section
				fluidRow(
					column(3,
								 div(strong("Filter by patient subgroup (optional):")),
								 virtualSelectInput(
								 	ns("recent_presc_freq_strat_variable"),
								 	label = NULL,
								 	choices = NULL,
								 	autoSelectFirstOption = FALSE,
								 	search = FALSE,
								 	disableOptionGroupCheckbox = TRUE,
								 	hideClearButton = FALSE,
								 	dropboxWrapper = "body"
								 )
					),
					column(3,
								 div(strong("Filter by LTCs (optional):")),
								 uiOutput(ns("recent_ltc_dropdown_ui"))
					),
					column(3,
								 div(strong("Filter by background prescriptions (optional):")),
								 uiOutput(ns("recent_bg_presc_dropdown_ui"))
					),
					column(3,
								 # Empty column to maintain layout
								 div()
					)
				),

				# Patient count display
				uiOutput(ns("recent_presc_patient_count")),
				br(),

				downloadButton(ns("download_recent_presc_freq"), "Download Table", class = "btn-sm btn-secondary"),
				br(), br(),
				reactableOutput(ns("recent_presc_freq_table"), height = "50em")
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

		# ============================================================================
		# DYNAMIC DROPDOWNS FOR FILTERS
		# ============================================================================

		# LTC filter dropdown (for LTC table - to see co-morbidities)
		output$ltc_filter_dropdown_ui <- renderUI({
			req(ltcs_r())
			unique_ltcs <- sort(unique(ltcs_r()$term))

			virtualSelectInput(
				ns("ltc_filter_dropdown"),
				label = NULL,
				choices = unique_ltcs,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Background prescription dropdown (for LTC table)
		output$presc_dropdown_ui <- renderUI({
			req(prescriptions_r(), bnf_level())  # Added bnf_level dependency

			# Use bnf_level directly (it's reactive from parent)
			current_level <- bnf_level()
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
				label = NULL,
				choices = choices,
				multiple = TRUE,
				search = TRUE,
				showSelectedOptionsFirst = TRUE
			)
		})


		# Prescription filter dropdown (for prescription table - to see medications prescribed together)
		output$presc_filter_dropdown_ui <- renderUI({
			req(prescriptions_r(), bnf_level())  # Added bnf_level dependency

			current_level <- bnf_level()
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
				ns("presc_filter_dropdown"),
				label = NULL,
				choices = choices,
				multiple = TRUE,
				search = TRUE,
				showSelectedOptionsFirst = TRUE
			)
		})

		# LTC dropdown (for prescription tables)
		output$ltc_dropdown_ui <- renderUI({
			req(ltcs_r())
			unique_ltcs <- sort(unique(ltcs_r()$term))

			virtualSelectInput(
				ns("ltc_dropdown"),
				label = NULL,
				choices = unique_ltcs,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Recent prescription dropdown for LTC table
		output$ltc_by_presc_recent_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(cases_controls_r(), bnf_level())

			current_level <- bnf_level()
			unique_substances <- sort(unique(cases_controls_r()$substance))


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
				ns("ltc_by_presc_recent_dropdown"),
				label = NULL,
				choices = choices,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Recent prescription dropdown for prescription table
		output$presc_by_ltc_recent_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(cases_controls_r(), bnf_level())
			current_level <- bnf_level()
			unique_substances <- sort(unique(cases_controls_r()$substance))

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
				ns("presc_by_ltc_recent_dropdown"),
				label = NULL,
				choices = choices,
				multiple = TRUE,
				search = TRUE
			)
		})

		# LTC dropdown for recent prescription table
		output$recent_ltc_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(ltcs_r())

			unique_ltcs <- sort(unique(ltcs_r()$term))

			virtualSelectInput(
				ns("recent_ltc_dropdown"),
				label = NULL,
				choices = unique_ltcs,
				multiple = TRUE,
				search = TRUE
			)
		})

		# Background prescription dropdown for recent prescription table
		output$recent_bg_presc_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(prescriptions_r(), bnf_level())  # Added bnf_level dependency

			current_level <- bnf_level()
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
				ns("recent_bg_presc_dropdown"),
				label = NULL,
				choices = choices,
				multiple = TRUE,
				search = TRUE,
				showSelectedOptionsFirst = TRUE
			)
		})

		# ============================================================================
		# HELPER FUNCTION FOR PATIENT COUNTS
		# ============================================================================

		get_filtered_patient_counts <- function(data, strat_var, filter1_patids = NULL, filter2_patids = NULL, filter3_patids = NULL, patient_data) {
			# Get total counts for each group (before filtering)
			total_case_count <- uniqueN(data[group == "case", patid])
			total_control_count <- uniqueN(data[group == "control", patid])

			# Start with all patients
			filtered_data <- data

			# Apply stratification first
			if (!is.null(strat_var) && strat_var != "") {
				filtered_data <- apply_patient_stratification(filtered_data, strat_var, patient_data)
			}

			# Then apply additional filters
			all_patids <- unique(filtered_data$patid)

			if (!is.null(filter1_patids)) {
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter1_patids)))
			}

			if (!is.null(filter2_patids)) {
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter2_patids)))
			}

			if (!is.null(filter3_patids)) {
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter3_patids)))
			}

			# Filter to final patient set
			filtered_data <- filtered_data[patid %in% all_patids]

			# Count by group (after filtering)
			case_count <- uniqueN(filtered_data[group == "case", patid])
			control_count <- uniqueN(filtered_data[group == "control", patid])
			total_count <- case_count + control_count

			# Calculate percentages relative to total in each group
			case_pct <- if (total_case_count > 0) round(case_count / total_case_count * 100, 1) else 0
			control_pct <- if (total_control_count > 0) round(control_count / total_control_count * 100, 1) else 0

			list(
				case_n = case_count,
				case_pct = case_pct,
				control_n = control_count,
				control_pct = control_pct,
				total_n = total_count
			)
		}


		# ============================================================================
		# LTC PREVALENCE TABLES
		# ============================================================================

		# Store rendered table data
		ltc_freq_table_data <- reactiveVal(NULL)
		presc_freq_table_data <- reactiveVal(NULL)
		recent_presc_freq_table_data <- reactiveVal(NULL)

		# Patient count for LTC table
		output$ltc_patient_count <- renderUI({
			req(ltcs_r(), patient_data_r())

			# Get filter patient IDs
			ltc_filter_patids <- NULL
			bg_presc_patids <- NULL
			recent_presc_patids <- NULL

			if (isTruthy(input$ltc_filter_dropdown)) {
				ltc_filter_patids <- unique(ltcs_r()[term %in% input$ltc_filter_dropdown, patid])
			}

			if (isTruthy(input$presc_dropdown)) {
				bg_presc_patids <- unique(prescriptions_r()[substance %in% input$presc_dropdown, patid])
			}

			if (!is.null(cases_controls_r) && isTruthy(input$ltc_by_presc_recent_dropdown)) {
				recent_presc_patids <- unique(cases_controls_r()[substance %in% input$ltc_by_presc_recent_dropdown, patid])
			}

			counts <- get_filtered_patient_counts(
				ltcs_r(),
				input$ltc_freq_strat_variable,
				ltc_filter_patids,
				bg_presc_patids,
				recent_presc_patids,
				patient_data_r()
			)

			div(
				style = "padding: 10px; background-color: #f8f9fa; border-radius: 4px; margin-top: 10px;",
				strong("Filtered patients: "),
				sprintf("Cases: %s (%s%%) | Controls: %s (%s%%) | Total: %s",
								prettyNum(counts$case_n, big.mark = ","),
								counts$case_pct,
								prettyNum(counts$control_n, big.mark = ","),
								counts$control_pct,
								prettyNum(counts$total_n, big.mark = ","))
			)
		})

		# LTC frequency table
		output$ltc_freq_table <- renderReactable({
			req(ltcs_r(), patient_data_r())

			# Apply stratification first
			ltcs_filtered <- ltcs_r()
			if (!is.null(input$ltc_freq_strat_variable) && input$ltc_freq_strat_variable != "") {
				ltcs_filtered <- apply_patient_stratification(
					ltcs_filtered,
					input$ltc_freq_strat_variable,
					patient_data_r()
				)
			}

			# Then apply additional filters
			all_patids <- unique(ltcs_filtered$patid)

			if (isTruthy(input$ltc_filter_dropdown)) {
				ltc_filter_patids <- unique(ltcs_r()[term %in% input$ltc_filter_dropdown, patid])
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(ltc_filter_patids)))
			}

			if (isTruthy(input$presc_dropdown)) {
				bg_presc_patids <- unique(prescriptions_r()[substance %in% input$presc_dropdown, patid])
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(bg_presc_patids)))
			}

			if (!is.null(cases_controls_r) && isTruthy(input$ltc_by_presc_recent_dropdown)) {
				recent_presc_patids <- unique(cases_controls_r()[substance %in% input$ltc_by_presc_recent_dropdown, patid])
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(recent_presc_patids)))
			}

			ltcs_filtered <- ltcs_filtered[patid %in% all_patids]

			# Calculate frequency stats
			freq_stats <- calculate_frequency_stats(ltcs_filtered, "term")
			result_table <- create_prevalence_ratio_table(
				freq_stats,
				"term",
				data_with_group = ltcs_filtered
			)

			validate(need(nrow(result_table) > 0, "No results found for selected filters"))

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
		# PRESCRIPTION PREVALENCE TABLES
		# ============================================================================

		# Patient count for prescription table
		output$presc_patient_count <- renderUI({
			req(prescriptions_r(), patient_data_r())

			# Get filter patient IDs
			presc_filter_patids <- NULL
			ltc_patids <- NULL
			recent_presc_patids <- NULL

			if (isTruthy(input$presc_filter_dropdown)) {
				presc_filter_patids <- unique(prescriptions_r()[substance %in% input$presc_filter_dropdown, patid])
			}

			if (isTruthy(input$ltc_dropdown)) {
				ltc_patids <- unique(ltcs_r()[term %in% input$ltc_dropdown, patid])
			}

			if (!is.null(cases_controls_r) && isTruthy(input$presc_by_ltc_recent_dropdown)) {
				recent_presc_patids <- unique(cases_controls_r()[substance %in% input$presc_by_ltc_recent_dropdown, patid])
			}

			counts <- get_filtered_patient_counts(
				prescriptions_r(),
				input$presc_freq_strat_variable,
				presc_filter_patids,
				ltc_patids,
				recent_presc_patids,
				patient_data_r()
			)

			div(
				style = "padding: 10px; background-color: #f8f9fa; border-radius: 4px; margin-top: 10px;",
				strong("Filtered patients: "),
				sprintf("Cases: %s (%s%%) | Controls: %s (%s%%) | Total: %s",
								prettyNum(counts$case_n, big.mark = ","),
								counts$case_pct,
								prettyNum(counts$control_n, big.mark = ","),
								counts$control_pct,
								prettyNum(counts$total_n, big.mark = ","))
			)
		})

		# Prescription frequency table
		output$presc_freq_table <- renderReactable({
			req(prescriptions_r(), patient_data_r())

			# Apply stratification first
			presc_filtered <- prescriptions_r()
			if (!is.null(input$presc_freq_strat_variable) && input$presc_freq_strat_variable != "") {
				presc_filtered <- apply_patient_stratification(
					presc_filtered,
					input$presc_freq_strat_variable,
					patient_data_r()
				)
			}

			# Then apply additional filters
			all_patids <- unique(presc_filtered$patid)

			if (isTruthy(input$ltc_dropdown)) {
				ltc_patids <- unique(ltcs_r()[term %in% input$ltc_dropdown, patid])
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(ltc_patids)))
			}

			if (!is.null(cases_controls_r) && isTruthy(input$presc_by_ltc_recent_dropdown)) {
				recent_presc_patids <- unique(cases_controls_r()[substance %in% input$presc_by_ltc_recent_dropdown, patid])
				all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(recent_presc_patids)))
			}

			presc_filtered <- presc_filtered[patid %in% all_patids]

			# Calculate frequency stats
			freq_stats <- calculate_frequency_stats(presc_filtered, "substance")
			result_table <- create_prevalence_ratio_table(
				freq_stats,
				"substance",
				data_with_group = presc_filtered
			)

			validate(need(nrow(result_table) > 0, "No results found for selected filters"))

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
		# RECENT PRESCRIPTIONS (from cases_controls_r)
		# ============================================================================

		if (!is.null(cases_controls_r)) {

			# Patient count for recent prescription table
			output$recent_presc_patient_count <- renderUI({
				req(cases_controls_r(), patient_data_r())

				# Get filter patient IDs
				ltc_patids <- NULL
				bg_presc_patids <- NULL

				if (isTruthy(input$recent_ltc_dropdown)) {
					ltc_patids <- unique(ltcs_r()[term %in% input$recent_ltc_dropdown, patid])
				}

				if (isTruthy(input$recent_bg_presc_dropdown)) {
					bg_presc_patids <- unique(prescriptions_r()[substance %in% input$recent_bg_presc_dropdown, patid])
				}

				# Add strata for consistency with other data
				recent_presc_with_strata <- merge(
					cases_controls_r(),
					patient_data_r()[, .(patid, strata)],
					by = "patid"
				)

				# FIXED: Only pass 2 filter parameters instead of 3
				# Removed the third parameter which was incorrectly trying to filter recent prescriptions by recent prescriptions
				counts <- get_filtered_patient_counts(
					recent_presc_with_strata,
					input$recent_presc_freq_strat_variable,
					ltc_patids,
					bg_presc_patids,
					NULL,  # No third filter needed
					patient_data_r()
				)

				div(
					style = "padding: 10px; background-color: #f8f9fa; border-radius: 4px; margin-top: 10px;",
					strong("Filtered patients: "),
					sprintf("Cases: %s (%s%%) | Controls: %s (%s%%) | Total: %s",
									prettyNum(counts$case_n, big.mark = ","),
									counts$case_pct,
									prettyNum(counts$control_n, big.mark = ","),
									counts$control_pct,
									prettyNum(counts$total_n, big.mark = ","))
				)
			})
			# Recent prescription frequency table
			output$recent_presc_freq_table <- renderReactable({
				req(cases_controls_r(), patient_data_r())

				# Add strata for calculate_frequency_stats
				recent_presc_with_strata <- merge(
					cases_controls_r(),
					patient_data_r()[, .(patid, strata)],
					by = "patid"
				)

				# Apply stratification first
				recent_presc_filtered <- recent_presc_with_strata
				if (!is.null(input$recent_presc_freq_strat_variable) && input$recent_presc_freq_strat_variable != "") {
					recent_presc_filtered <- apply_patient_stratification(
						recent_presc_filtered,
						input$recent_presc_freq_strat_variable,
						patient_data_r()
					)
				}

				# Then apply additional filters
				all_patids <- unique(recent_presc_filtered$patid)

				if (isTruthy(input$recent_ltc_dropdown)) {
					ltc_patids <- unique(ltcs_r()[term %in% input$recent_ltc_dropdown, patid])
					all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(ltc_patids)))
				}

				if (isTruthy(input$recent_bg_presc_dropdown)) {
					bg_presc_patids <- unique(prescriptions_r()[substance %in% input$recent_bg_presc_dropdown, patid])
					all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(bg_presc_patids)))
				}

				recent_presc_filtered <- recent_presc_filtered[patid %in% all_patids]

				# Calculate frequency stats
				freq_stats <- calculate_frequency_stats(recent_presc_filtered, "substance")
				result_table <- create_prevalence_ratio_table(
					freq_stats,
					"substance",
					data_with_group = recent_presc_filtered
				)

				validate(need(nrow(result_table) > 0, "No results found for selected filters"))

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
		}

		# ============================================================================
		# DOWNLOAD HANDLERS
		# ============================================================================

		output$download_ltc_freq <- downloadHandler(
			filename = function() {
				paste0("ltc_prevalence_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(ltc_freq_table_data())
				fwrite(ltc_freq_table_data(), file)
			}
		)

		output$download_presc_freq <- downloadHandler(
			filename = function() {
				paste0("prescription_prevalence_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(presc_freq_table_data())
				fwrite(presc_freq_table_data(), file)
			}
		)

		output$download_recent_presc_freq <- downloadHandler(
			filename = function() {
				paste0("recent_prescription_prevalence_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(recent_presc_freq_table_data())
				fwrite(recent_presc_freq_table_data(), file)
			}
		)

	})
}