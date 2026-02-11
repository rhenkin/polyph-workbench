# R/module_cca_temporal.R

#' Temporal Distribution Analysis UI
#'
#' @param id Module namespace ID
#' @return Shiny UI element
module_cca_temporal_ui <- function(id) {
	ns <- NS(id)

	accordion(
		open = FALSE,
		accordion_panel(
			title = "Time-to analysis",
			value = "ltc_prev_tables",
			icon = bs_icon("heart-pulse"),
			card(
				full_screen = TRUE,

			# Risk window info banner
			uiOutput(ns("risk_window_banner")),

			# Filter controls
			card(
				card_header("Filters", class = "bg-light"),
				fluidRow(
					column(4,
								 virtualSelectInput(
								 	ns("strat_variable"),
								 	label = "Patient subgroup (optional):",
								 	choices = NULL,
								 	autoSelectFirstOption = FALSE,
								 	search = FALSE,
								 	multiple = FALSE,
								 	hideClearButton = FALSE,
								 	dropboxWrapper = "body"
								 )
					),
					column(4,
								 div(strong("Filter by LTCs (optional):")),
								 uiOutput(ns("ltc_filter_ui"))
					),
					column(4,
								 div(strong("Filter by prescription (optional):")),
								 uiOutput(ns("substance_filter_ui"))
					)
				),

				# Patient count
				fluidRow(
					column(12,
								 uiOutput(ns("patient_count_display"))
					)
				)
			),

			# Main visualization tabs
			accordion(
				id = "temporal_accs",

				# Tab 1: Overall distribution histogram
				accordion_panel(
					title = "Distribution",
					icon = bs_icon("bar-chart"),
					card_body(
						p("Overall distribution of events across the risk window.
	            Early events (within first week) are more likely medication-induced."),
						vegawidgetOutput(ns("temporal_histogram"), height = "400px"),
						br(),
						reactableOutput(ns("summary_stats_table"))
					)
				),

				# Tab 2: By substance
				accordion_panel(
					title = "By Substance",
					icon = bs_icon("capsule"),
					card_body(
						p("Temporal patterns for individual substances."),

						fluidRow(
							column(6,
										 radioButtons(
										 	ns("substance_view_type"),
										 	"View type:",
										 	choices = c("Table" = "table", "Heatmap" = "heatmap"),
										 	inline = TRUE
										 )
							),
							column(6,
										 conditionalPanel(
										 	condition = "input.substance_view_type == 'heatmap'",
										 	ns = ns,
										 	numericInput(
										 		ns("heatmap_top_n"),
										 		"Number of substances to show:",
										 		value = 20,
										 		min = 5,
										 		max = 50,
										 		step = 5
										 	)
										 )
							)
						),

						conditionalPanel(
							condition = "input.substance_view_type == 'table'",
							ns = ns,
							downloadButton(ns("download_substance_table"),
														 "Download Table",
														 class = "btn-sm btn-secondary"),
							br(), br(),
							reactableOutput(ns("substance_table"), height = "600px")
						),

						conditionalPanel(
							condition = "input.substance_view_type == 'heatmap'",
							ns = ns,
							vegawidgetOutput(ns("substance_heatmap"), height = "600px")
						)
					)
				),

				# Tab 3: Stratified analysis
				accordion_panel(
					title = "Stratified",
					icon = bs_icon("grid-3x3"),
					card_body(
						p("Compare temporal patterns across patient subgroups."),

						fluidRow(
							column(6,
										 selectInput(
										 	ns("stratify_by"),
										 	"Stratify by:",
										 	choices = c(
										 		"Polypharmacy level" = "pp_group",
										 		"Sex" = "sex",
										 		#"Age group" = "age_group",
										 		"IMD quintile" = "imd_quintile",
										 		"Number of LTCs" = "mltc_group"
										 	)
										 )
							),
							column(6,
										 radioButtons(
										 	ns("strat_plot_type"),
										 	"Plot type:",
										 	choices = c("Faceted bars" = "facet", "Stacked bars" = "stacked"),
										 	inline = TRUE
										 )
							)
						),

						vegawidgetOutput(ns("stratified_plot"), height = "500px")
					)
				),

				# Tab 4: Early vs Late comparison
				accordion_panel(
					title = "Early vs Late",
					icon = bs_icon("clock-history"),
					card_body(
						p("Compare substance distribution splitting at median."),

						fluidRow(
							column(6,
										 card(
										 	card_header("Early Events", class = "bg-danger text-white"),
										 	uiOutput(ns("early_events_summary"))
										 )
							),
							column(6,
										 card(
										 	card_header("Late Events", class = "bg-warning"),
										 	uiOutput(ns("late_events_summary"))
										 )
							)
						),

						br(),

						h5("Top substances comparison"),
						fluidRow(
							column(6,
										 h6("Early Events"),
										 reactableOutput(ns("early_substances_table"))
							),
							column(6,
										 h6("Late Events"),
										 reactableOutput(ns("late_substances_table"))
							)
						)
			)
		)
		)
)
)
	)
}

#' Temporal Distribution Analysis Server
#'
#' @param id Module namespace ID
#' @param prescriptions_r Reactive returning prescriptions data.table
#' @param patient_data_r Reactive returning patient data.table
#' @param ltcs_r Reactive returning LTCs data.table
#' @param metadata_r Reactive returning metadata list with pred_window
#' @param bnf_level Character, the BNF aggregation level
module_cca_temporal_server <- function(id, prescriptions_r, patient_data_r,
																			 ltcs_r, metadata_r, bnf_level = "BNF_Chemical_Substance") {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# ========================================================================
		# 1. REACTIVE DATA PREPARATION
		# ========================================================================

		# Get risk window from metadata
		risk_window_r <- reactive({
			req(metadata_r())
			metadata <- metadata_r()
			# Extract pred_window from metadata
			pred_window <- metadata$flow_diagram$pred_window

			if (is.null(pred_window) || is.na(pred_window)) {
				# Fallback: try to get from metadata directly
				pred_window <- metadata$pred_window
			}

			if (is.null(pred_window) || is.na(pred_window)) {
				warning("pred_window not found in metadata. Using default of 30 days.")
				return(30)
			}

			return(pred_window)
		})

		# Prepare base temporal data (cases only)
		temporal_data_base <- reactive({
			req(prescriptions_r(), risk_window_r())
			presc <- copy(prescriptions_r())

			# Filter to cases only - controls don't have temporal relationship
			presc <- presc[group == "case"]

			# Must have required columns - CCA uses index_date and start_date
			req(all(c("patid", "index_date", "outcome_date") %in% names(presc)))

			# Add temporal grouping
			presc <- add_temporal_groups(presc, risk_window_r())

			# Merge with patient data for stratification variables
			if (!is.null(patient_data_r())) {
				patient_info <- patient_data_r()[group == "case",
																				 .(patid, sex, eth_group, imd_quintile, pp_group, mltc_group)]
				presc <- merge(presc, patient_info, by = "patid", all.x = TRUE)
			}

			return(presc)
		})

		# Apply all filters
		temporal_data_filtered <- reactive({
			req(temporal_data_base())

			dt <- copy(temporal_data_base())

			# Filter by stratification variable
			if (!is.null(input$strat_variable) && input$strat_variable != "") {
				dt <- filter_by_stratification(dt, input$strat_variable)
			}

			# Filter by LTCs
			if (!is.null(input$ltc_filter) && length(input$ltc_filter) > 0) {
				req(ltcs_r())
				ltc_data <- ltcs_r()[group == "case"]

				# Get patients with ALL selected LTCs
				valid_patids <- ltc_data[term %in% input$ltc_filter, unique(patid)]

				dt <- dt[patid %in% valid_patids]
			}

			# Filter by substances
			if (!is.null(input$substance_filter) && length(input$substance_filter) > 0) {
				dt <- dt[substance %in% input$substance_filter]
			}

			return(dt)
		})

		# ========================================================================
		# 2. UPDATE FILTER CHOICES
		# ========================================================================

		observe({
			req(patient_data_r())
			choices <- create_stratification_choices(patient_data_r())
			updateVirtualSelect("strat_variable", choices = choices)
		})

		output$ltc_filter_ui <- renderUI({
			req(ltcs_r())

			virtualSelectInput(
				ns("ltc_filter"),
				label = NULL,
				choices = create_ltc_dropdown_choices(ltcs_r()),
				multiple = TRUE,
				search = TRUE,
				showSelectedOptionsFirst = TRUE,
				dropboxWrapper = "body"
			)
		})

		output$substance_filter_ui <- renderUI({

			req(temporal_data_base(), bnf_level())

			virtualSelectInput(
				ns("substance_filter"),
				label = NULL,
				choices = create_bnf_dropdown_choices(temporal_data_base(), bnf_level()),
				multiple = TRUE,
				search = TRUE,
				showSelectedOptionsFirst = TRUE,
				dropboxWrapper = "body"
			)
		})

		# ========================================================================
		# 3. DISPLAY OUTPUTS
		# ========================================================================

		# Risk window banner
		output$risk_window_banner <- renderUI({
			req(risk_window_r())

			window <- risk_window_r()
			groups <- create_day_groups(window)

			div(
				class = "alert alert-info",
				style = "margin-bottom: 20px;",
				icon("info-circle"), " ",
				strong("Risk window: "), sprintf("%d days", window), " | ",
				strong("Grouping: "), paste(groups$labels, collapse = " • ")
			)
		})

		# Patient count display
		output$patient_count_display <- renderUI({
			req(temporal_data_filtered())

			n_total <- uniqueN(temporal_data_base()$patid)
			n_filtered <- uniqueN(temporal_data_filtered()$patid)

			div(
				style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px; margin-top: 10px;",
				strong("Showing: "),
				span(prettyNum(n_filtered, big.mark = ","), style = "font-size: 1.2em; color: #0066cc;"),
				" cases",
				if (n_filtered < n_total) {
					span(" (filtered from ", prettyNum(n_total, big.mark = ","), " total)")
				}
			)
		})

		# Tab 1: Histogram
		output$temporal_histogram <- renderVegawidget({
			req(temporal_data_filtered(), risk_window_r())

			create_temporal_histogram(temporal_data_filtered(), risk_window_r())
		})

		output$summary_stats_table <- renderReactable({
			req(temporal_data_filtered())
			create_temporal_summary_stats(temporal_data_filtered())
		})

		# Tab 2: Substance analysis
		substance_table_data <- reactive({
			req(temporal_data_filtered())
			create_temporal_substance_table(temporal_data_filtered())
		})

		output$substance_table <- renderReactable({
			req(substance_table_data())
			render_temporal_substance_table(substance_table_data())
		})

		output$substance_heatmap <- renderVegawidget({
			req(temporal_data_filtered(), input$heatmap_top_n)
			create_temporal_heatmap(temporal_data_filtered(), input$heatmap_top_n)
		})

		output$download_substance_table <- downloadHandler(
			filename = function() {
				paste0("temporal_by_substance_", Sys.Date(), ".csv")
			},
			content = function(file) {
				fwrite(substance_table_data(), file)
			}
		)

		# Tab 3: Stratified analysis
		output$stratified_plot <- renderVegawidget({
			req(temporal_data_filtered(), input$stratify_by, input$strat_plot_type)

			if (input$strat_plot_type == "facet") {
				create_stratified_temporal_plot(
					temporal_data_filtered(),
					input$stratify_by,
					risk_window_r()
				)
			} else {
				create_stacked_temporal_plot(
					temporal_data_filtered(),
					input$stratify_by,
					risk_window_r()
				)
			}
		})

		# Tab 4: Early vs Late comparison
		early_late_split <- reactive({
			req(temporal_data_filtered())
			dt <- temporal_data_filtered()
			median_day <- median(dt$days_to_event, na.rm = TRUE)

			list(
				early = dt[days_to_event <= median_day],
				late = dt[days_to_event > median_day],
				median_day = median_day
			)
		})

		output$early_events_summary <- renderUI({
			req(early_late_split())

			early <- early_late_split()$early
			median_day <- early_late_split()$median_day

			tagList(
				p(strong("Events within first ", round(median_day), " days")),
				p("N = ", prettyNum(uniqueN(early$patid), big.mark = ",")),
				p("Median days to event: ", round(median(early$days_to_event, na.rm = TRUE), 1)),
				p("Median polypharmacy: ", round(median(as.numeric(gsub("[^0-9]", "", early$pp_group)), na.rm = TRUE), 1))
			)
		})

		output$late_events_summary <- renderUI({
			req(early_late_split())

			late <- early_late_split()$late
			median_day <- early_late_split()$median_day

			tagList(
				p(strong("Events after ", round(median_day), " days")),
				p("N = ", prettyNum(uniqueN(late$patid), big.mark = ",")),
				p("Median days to event: ", round(median(late$days_to_event, na.rm = TRUE), 1)),
				p("Median polypharmacy: ", round(median(as.numeric(gsub("[^0-9]", "", late$pp_group)), na.rm = TRUE), 1))
			)
		})

		output$early_substances_table <- renderReactable({
			req(early_late_split())
			create_top_substances_table(early_late_split()$early, top_n = 10)
		})

		output$late_substances_table <- renderReactable({
			req(early_late_split())
			create_top_substances_table(early_late_split()$late, top_n = 10)
		})
	})
}