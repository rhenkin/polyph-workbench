
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

		vs_input <- function(id, choices) {
			virtualSelectInput(
				id,
				label = NULL,
				choices = choices,
				multiple = TRUE,
				search = TRUE,
				showSelectedOptionsFirst = TRUE
			)
		}

		# LTC filter dropdown (for LTC table - to see co-morbidities)
		output$ltc_filter_dropdown_ui <- renderUI({
			req(ltcs_r())

			vs_input(ns("ltc_filter_dropdown"),
							 create_ltc_dropdown_choices(ltcs_r()))
		})

		# Background prescription dropdown (for LTC table)
		output$presc_dropdown_ui <- renderUI({
			req(prescriptions_r(), bnf_level())  # Added bnf_level dependency

			vs_input(ns("presc_dropdown"),
							 create_bnf_dropdown_choices(prescriptions_r(), bnf_level()))
		})


		# Prescription filter dropdown (for prescription table - to see medications prescribed together)
		output$presc_filter_dropdown_ui <- renderUI({
			req(prescriptions_r(), bnf_level())  # Added bnf_level dependency

			vs_input(ns("presc_filter_dropdown"),
							 create_bnf_dropdown_choices(prescriptions_r(), bnf_level()))
		})

		# LTC dropdown (for prescription tables)
		output$ltc_dropdown_ui <- renderUI({
			req(ltcs_r())

			vs_input(ns("ltc_dropdown"),
							 create_ltc_dropdown_choices(ltcs_r()))
		})

		# Recent prescription dropdown for LTC table
		output$ltc_by_presc_recent_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(cases_controls_r(), bnf_level())

			vs_input(ns("ltc_by_presc_recent_dropdown"),
							 create_bnf_dropdown_choices(cases_controls_r(), bnf_level()))

		})

		# Recent prescription dropdown for prescription table
		output$presc_by_ltc_recent_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(cases_controls_r(), bnf_level())

			vs_input(ns("presc_by_ltc_recent_dropdown"),
							 create_bnf_dropdown_choices(cases_controls_r(), bnf_level()))
		})

		# LTC dropdown for recent prescription table
		output$recent_ltc_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(ltcs_r())

			vs_input(ns("recent_ltc_dropdown"),
							 create_ltc_dropdown_choices(ltcs_r()))
		})

		# Background prescription dropdown for recent prescription table
		output$recent_bg_presc_dropdown_ui <- renderUI({
			if (is.null(cases_controls_r)) return(NULL)
			req(prescriptions_r(), bnf_level())  # Added bnf_level dependency

			vs_input(ns("recent_bg_presc_dropdown"),
							 create_bnf_dropdown_choices(prescriptions_r(), bnf_level()))
		})


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

			render_patient_count_ui(counts)
		})

		strat_filter_ltc_data <- reactive({
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

			ltcs_filtered

		})


		ltc_freq_patids <- reactive({
			ltcs_filtered <- strat_filter_ltc_data()

			ltc_filter_patids <- if (isTruthy(input$ltc_filter_dropdown)) {
				unique(ltcs_r()[term %in% input$ltc_filter_dropdown, patid])
			}
			bg_presc_patids <- if (isTruthy(input$presc_dropdown)) {
				unique(prescriptions_r()[substance %in% input$presc_dropdown, patid])
			}
			recent_presc_patids <- if (!is.null(cases_controls_r) && isTruthy(input$ltc_by_presc_recent_dropdown)) {
				unique(cases_controls_r()[substance %in% input$ltc_by_presc_recent_dropdown, patid])
			}

			get_filtered_patids(ltcs_filtered, ltc_filter_patids, bg_presc_patids, recent_presc_patids)
		})

		ltc_freq_data <- reactive({
			req(ltcs_r(), patient_data_r(), length(ltc_freq_patids())>0)

			calculate_filtered_prevalence_table(
				strat_filter_ltc_data(),
				ltc_freq_patids(),
				"term"
			)
		}) |> bindCache(ltc_freq_patids())

		# LTC frequency table
		output$ltc_freq_table <- renderReactable({
			render_prevalence_reactable(ltc_freq_data(), "LTC")
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

			render_patient_count_ui(counts)
		})

		strat_filter_presc_data <- reactive({
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

			presc_filtered
		})

		presc_freq_patids <- reactive({
			presc_filtered <- strat_filter_presc_data()

			presc_filter_patids <- if (isTruthy(input$presc_filter_dropdown)) {
				unique(prescriptions_r()[substance %in% input$presc_filter_dropdown, patid])
			}
			ltc_patids <- if (isTruthy(input$ltc_dropdown)) {
				unique(ltcs_r()[term %in% input$ltc_dropdown, patid])
			}
			recent_presc_patids <- if (!is.null(cases_controls_r) && isTruthy(input$presc_by_ltc_recent_dropdown)) {
				unique(cases_controls_r()[substance %in% input$presc_by_ltc_recent_dropdown, patid])
			}

			get_filtered_patids(presc_filtered, presc_filter_patids, ltc_patids, recent_presc_patids)
		})

		presc_freq_data <- reactive({
			req(prescriptions_r(), patient_data_r(), length(presc_freq_patids())>0)

			calculate_filtered_prevalence_table(
				strat_filter_presc_data(),
				presc_freq_patids(),
				"substance"
			)

		}) |> bindCache(presc_freq_patids())

		# Prescription frequency table
		output$presc_freq_table <- renderReactable({
			render_prevalence_reactable(presc_freq_data(), "Substance")
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

				counts <- get_filtered_patient_counts(
					recent_presc_with_strata,
					input$recent_presc_freq_strat_variable,
					ltc_patids,
					bg_presc_patids,
					NULL,
					patient_data_r()
				)

				render_patient_count_ui(counts)
			})

			strat_filter_recent_presc_data <- reactive({
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

				recent_presc_filtered
			})

			recent_presc_freq_patids <- reactive({
				recent_presc_filtered <- strat_filter_recent_presc_data()

				ltc_patids <- if (isTruthy(input$recent_ltc_dropdown)) {
					unique(ltcs_r()[term %in% input$recent_ltc_dropdown, patid])
				}
				bg_presc_patids <- if (isTruthy(input$recent_bg_presc_dropdown)) {
					unique(prescriptions_r()[substance %in% input$recent_bg_presc_dropdown, patid])
				}

				get_filtered_patids(recent_presc_filtered, ltc_patids, bg_presc_patids, NULL)
			})

			recent_presc_freq_data <- reactive({
				req(cases_controls_r(), patient_data_r(), length(recent_presc_freq_patids())>0)

				calculate_filtered_prevalence_table(
					strat_filter_recent_presc_data(),
					recent_presc_freq_patids(),
					"substance"
				)
			}) |> bindCache(recent_presc_freq_patids(), cases_controls_r())

			# Recent prescription frequency table
			output$recent_presc_freq_table <- renderReactable({
				render_prevalence_reactable(recent_presc_freq_data(), "Substance")
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
				req(ltc_freq_data())
				to_save <- merge(ltc_freq_data(), ltc_chapters, by.x = "term", by.y = "ltc")
				fwrite(to_save, file)
			}
		)

		output$download_presc_freq <- downloadHandler(
			filename = function() {
				paste0("prescription_prevalence_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(presc_freq_data())
				fwrite(presc_freq_data(), file)
			}
		)

		output$download_recent_presc_freq <- downloadHandler(
			filename = function() {
				paste0("recent_prescription_prevalence_", Sys.Date(), ".csv")
			},
			content = function(file) {
				req(recent_presc_freq_data())
				fwrite(recent_presc_freq_data(), file)
			}
		)

	})
}