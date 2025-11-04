library(bslib)
library(data.table)
library(DT)
library(MASS)
library(shiny)
library(vegawidget)
library(htmltools)
library(gt)
library(reactable)
library(shinyjs)
library(shinyWidgets)
library(ggplot2)
library(ppcor)
library(bsicons)

#setDTthreads(6)

set.seed(42)

source("global.R")

#' @export
ui <-
	page_fluid(
		useBusyIndicators(spinners = FALSE),
		useShinyjs(),  # Initialize shinyjs
		tags$head(
			# CSS to style the checkmark
			tags$style("
			.vscomp-dropbox-container  {z-index:99999 !important;}

      .action-complete {
        color: #28a745;
        margin-left: 10px;
        display: none;
      }

       /* Landing modal styling */
      .landing-box {
        border: 2px solid #e0e0e0;
        border-radius: 8px;
        padding: 30px;
        text-align: center;
        cursor: pointer;
        transition: all 0.3s ease;
        background-color: #ffffff;
        height: 100%;
        display: flex;
        flex-direction: column;
        justify-content: center;
        align-items: center;
        min-height: 250px;
      }

      .landing-box:hover {
        border-color: #2C3E50;
        box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        transform: translateY(-4px);
      }

      .landing-box-icon {
        font-size: 4rem;
        margin-bottom: 20px;
        color: #2C3E50;
      }

      .landing-box-title {
        font-size: 1.5rem;
        font-weight: bold;
        margin-bottom: 10px;
        color: #2C3E50;
      }

      .landing-box-description {
        font-size: 1rem;
        color: #666;
        line-height: 1.5;
      }

      .landing-modal-header {
        text-align: center;
        margin-bottom: 30px;
      }

      .landing-modal-title {
        font-size: 2.5rem;
        font-weight: bold;
        color: #2C3E50;
        margin-bottom: 10px;
      }

      .landing-modal-subtitle {
        font-size: 1.2rem;
        color: #666;
      }
    	")
		),
		theme = bs_theme(preset = "flatly", font_scale = 0.8),
		layout_sidebar(
			# Left sidebar (collapsible)
			sidebar = sidebar(
				title = "Patient Filters",
				id = "left_menu",
				position = "left",
				open = "open",
				bg = "#f7f7f7",
				collapsible = TRUE,
				width = "20%",  # Roughly equivalent to your 1fr
				left_card_ui("patient_filter_module")
			),

			# Nested layout for right sidebar
			layout_sidebar(
				sidebar = sidebar(
					title = "Drug Filters",
					id = "right_menu",
					position = "right",
					open = "open",
					bg = "#f7f7f7",
					collapsible = TRUE,
					width = "20%",  # Roughly equivalent to your 1fr
					right_card_ui("drug_filter_module")
				),

				# Main content (will take remaining space, ~4fr equivalent)
				middle_card_ui("polyph_module")
			),
			class = "p-0"
		)
	)

#' @export
server <- function(input, output, session) {
	# Reactive value to track if modal has been shown
	modal_shown <- reactiveVal(FALSE)

	# Show landing modal on startup
	observe(priority = 9999, x = {
		if (!modal_shown()) {
			showModal(
				modalDialog(
					size = "xl",
					easyClose = FALSE,
					footer = NULL,
					div(
						class = "landing-modal-header",
						div(class = "landing-modal-title", "Polypharmacy Workbench"),
						div(class = "landing-modal-subtitle",
								"Choose your analysis pathway to explore patterns in prescribing and multimorbidity")
					),
					layout_columns(
						col_widths = c(6, 6),
						# Box 1: Outcome Explorer
						actionLink(
							inputId = "nav_to_outcome_explorer",
							label = div(
								class = "landing-box",
								bs_icon("search", class = "landing-box-icon"),
								div(class = "landing-box-title", "Outcome Explorer"),
								div(class = "landing-box-description",
										"Explore prescribing patterns in patients with specific outcomes and create studies by matching to controls.")
							)
						),
						# Box 3: Case-Control Analysis
						actionLink(
							inputId = "nav_to_analysis",
							label = div(
								class = "landing-box",
								bs_icon("bar-chart-line", class = "landing-box-icon"),
								div(class = "landing-box-title", "Case-Control Analysis"),
								div(class = "landing-box-description",
										"Analyse previously saved matched case-control datasets. Compare prescribing patterns, disease prevalence, and polypharmacy burden between cases and controls.")
							)
						)
					)
				)
			)
			modal_shown(TRUE)
		}
	})

	# Navigation handlers
	observeEvent(input$nav_to_outcome_explorer, {
		removeModal()
		nav_select("polyph_module-middle_tab", "outcome_explorer", session)
	})

	observeEvent(input$nav_to_analysis, {
		removeModal()
		nav_select("polyph_module-middle_tab", "cca", session)
	})

	# Get patient data directly from data.table
	# patient_df <- gold_patient[, .(patid, sex = gender, eth_group, imd_quintile)]

	# Get outcome list directly from data.table
	outcome_list <- sort(unique(gold_outcomes$outcome))

	stored_queries <- reactiveValues()

	patientFilter_r <- reactiveValues()

	# Create a default reactive for analyzed_patient_count (before middle_card is initialized)
	analyzed_patient_count <- reactiveVal(0)

	# Run left side menu code
	# From module: patientFilter_r
	# To module: patient numbers and list out of outcomes
	left_card_server(
		id = "patient_filter_module",
		patientFilter_r = patientFilter_r,
		outcome_list = outcome_list,
		analyzed_patient_count = analyzed_patient_count
	)

	patient_data <- reactive({
		req(input$nav_to_outcome_explorer, patientFilter_r$outcome)
		message("Querying patient data ", Sys.time())

		# Use data.table function instead of database query
		res <- buildPatientQuery(gold_patient,gold_outcomes,gold_ltc,patientFilter_r$selected_ltcs,
														 patientFilter_r$input_list$eth_group,
														 patientFilter_r$input_list$sex,
														 patientFilter_r$input_list$imd_quintile,
														 patientFilter_r$outcome)

		# Store query info (now just the parameters used)
		stored_queries$patient <- list(
			terms = patientFilter_r$selected_ltcs,
			eth_group = patientFilter_r$input_list$eth_group,
			sex = patientFilter_r$input_list$sex,
			imd_quintile = patientFilter_r$input_list$imd_quintile,
			outcome = patientFilter_r$outcome
		)

		# Filter out NA imd_quintile values
		res <- res[imd_quintile != "NA"]
		setkey(res, patid)
		res
	}) |> bindCache(
		patientFilter_r$selected_ltcs,
		patientFilter_r$input_list$eth_group,
		patientFilter_r$input_list$sex,
		patientFilter_r$input_list$imd_quintile,
		patientFilter_r$outcome
	)

	outcome_data <- reactive({
		req(patientFilter_r$outcome)
		message("Querying outcome data ", Sys.time())

		# Use data.table function instead of database query
		# res <- buildOutcomeQuery(gold_outcomes, patientFilter_r$outcome)
		res <- buildOutcomeLtcOptionsQuery(gold_patient, gold_ltc, gold_outcomes,
																			 patientFilter_r$outcome,
																			 patientFilter_r$selected_ltcs,
																			 patientFilter_r$input_list$eth_group,
																			 patientFilter_r$input_list$sex,
																			 patientFilter_r$input_list$imd_quintile)
		# Store query info
		stored_queries$outcome <- list(outcome = patientFilter_r$outcome)

		setkey(res, patid, eventdate)
		res
	}) |> bindCache(
		patientFilter_r$outcome
	)

	prescription_data <- reactive({
		req(patientFilter_r$outcome)
		message("Querying prescription data ", Sys.time())

		# Use data.table function instead of database query
		res <- buildPrescriptionQuery(gold_cp,gold_patient,gold_ltc,gold_outcomes,
																	patientFilter_r$selected_ltcs,
																	NULL,
																	patientFilter_r$input_list$eth_group,
																	patientFilter_r$input_list$sex,
																	patientFilter_r$input_list$imd_quintile,
																	patientFilter_r$outcome)

		# Store query info
		stored_queries$prescription <- list(
			terms = patientFilter_r$selected_ltcs,
			eth_group = patientFilter_r$input_list$eth_group,
			sex = patientFilter_r$input_list$sex,
			imd_quintile = patientFilter_r$input_list$imd_quintile,
			outcome = patientFilter_r$outcome
		)

		setkey(res, patid, start_date)
		res
	}) |> bindCache(
		patientFilter_r$selected_ltcs,
		patientFilter_r$input_list$eth_group,
		patientFilter_r$input_list$sex,
		patientFilter_r$input_list$imd_quintile,
		patientFilter_r$outcome
	)

	acute_presc_df <- reactive({
		req(patientFilter_r$outcome)
		message("Querying acute prescription data ", Sys.time())

		# Use data.table function instead of database query
		res <- buildAcutePrescriptionQuery(gold_acute_presc,gold_patient,gold_outcomes,gold_ltc,
																			 patientFilter_r$selected_ltcs,
																			 NULL,
																			 patientFilter_r$input_list$eth_group,
																			 patientFilter_r$input_list$sex,
																			 patientFilter_r$input_list$imd_quintile,
																			 patientFilter_r$outcome)

		setkey(res, patid, start_date)
		res
	}) |> bindCache(
		patientFilter_r$selected_ltcs,
		patientFilter_r$input_list$eth_group,
		patientFilter_r$input_list$sex,
		patientFilter_r$input_list$imd_quintile,
		patientFilter_r$outcome
	)

	right_card_inputs <- right_card_server(id = "drug_filter_module")
	filtered_prescription_data <- reactive({
		req(prescription_data())
		filtered_substances <-
			get_filtered_substances(right_card_inputs$included, right_card_inputs$excluded, bnf_lookup)
		tmp_df <- prescription_data()[substance %in% filtered_substances]
		tmp_df[duration >= right_card_inputs$minimum_duration()]
	})

	ltc_data <- reactive({
		req(patientFilter_r$outcome)
		message("Querying LTC data ", Sys.time())

		# Use data.table function instead of database query
		res <- buildLtcQuery(gold_patient,gold_outcomes,gold_ltc,
												 patientFilter_r$selected_ltcs,
												 NULL,
												 patientFilter_r$input_list$eth_group,
												 patientFilter_r$input_list$sex,
												 patientFilter_r$input_list$imd_quintile,
												 patientFilter_r$outcome)
		# Store query info
		stored_queries$ltc <- list(
			terms = patientFilter_r$selected_ltcs,
			eth_group = patientFilter_r$input_list$eth_group,
			sex = patientFilter_r$input_list$sex,
			imd_quintile = patientFilter_r$input_list$imd_quintile,
			outcome = patientFilter_r$outcome
		)

		setkey(res, patid, age_days)
		res
	}) |> bindCache(
		patientFilter_r$selected_ltcs,
		patientFilter_r$input_list$eth_group,
		patientFilter_r$input_list$sex,
		patientFilter_r$input_list$imd_quintile,
		patientFilter_r$outcome
	)

	# Middle module
	middle_card_output <- middle_card_server(
		id = "polyph_module",
		prescription_data = filtered_prescription_data,
		ltc_data = ltc_data,
		patient_data = patient_data,
		outcome_data = outcome_data,
        acute_presc_df = acute_presc_df,
		min_nltc = patientFilter_r$input_list$min_nltc,
		stored_queries = stored_queries,
		polypharmacy_threshold = right_card_inputs$polypharmacy_threshold,
		earliest_treatment_end = right_card_inputs$earliest_treatment_end
	)

	# Update analyzed patient count when outcome_prescriptions changes
	observe({
		req(middle_card_output$outcome_prescriptions())
		analyzed_patient_count(uniqueN(middle_card_output$outcome_prescriptions()$patid))
	})

	observe({
		req(middle_card_output$seltab())
		sidebar_toggle(
			id = "left_menu",
			open = middle_card_output$seltab() == "outcome_explorer"
		)
		sidebar_toggle(
			id = "right_menu",
			open = middle_card_output$seltab() == "outcome_explorer"
		)
	})
}

shiny::shinyApp(ui, server)
