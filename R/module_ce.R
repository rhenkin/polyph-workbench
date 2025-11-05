# Case explorer module
# First tab in the middle card
# Calls sub-modules: overview, prescription explorer, LTC explorer
# Shows number of patients with both outcome and prescription data
module_ce_ui <- function(id) {
	ns <- NS(id)
	nav_panel("Outcome explorer", value = "outcome_explorer",
							navset_tab(
								nav_panel(
									title = "Overview",
									conditionalPanel(condition = "output.data_loaded == null", ns = ns,
																	 div("Select an outcome, patient filters and click 'Start analysis' to begin",
																	 		class = "mt-4 p-5 h4 border rounded")),
									conditionalPanel(condition = "output.data_loaded == true", ns = ns,
										layout_columns(
											col_widths = c(3,3,3,3),
											value_box(
												showcase_layout = showcase_top_right(width = "200px", max_height = "150px"),
												height = "150px",
												title = "Cases",
												value = textOutput(ns("value_box_cases")),
												theme = "red",
												showcase = bs_icon("people-fill")
											),
											value_box(
												showcase_layout = showcase_top_right(width = "200px", max_height = "150px"),
												height = "150px",
												title = "Median age",
												value = textOutput(ns("value_box_median_age")),
												theme = "red",
												showcase = bs_icon("person-lines-fill"),
												p("years")
											),
											value_box(
												showcase_layout = showcase_top_right(width = "200px", max_height = "150px"),
												height = "150px",
												title = "Median multimorbidity",
												value = textOutput(ns("value_box_median_ltc")),
												theme = "red",
												showcase = bs_icon("heart-pulse-fill"),
												p("long-term conditions")
											),
											value_box(
												showcase_layout = showcase_top_right(width = "200px", max_height = "150px"),
												height = "150px",
												title = "Median polypharmacy",
												value = textOutput(ns("value_box_median_pp")),
												theme = "red",
												showcase = bs_icon("capsule-pill"),
												p("concurrent medications")
											)
										),
										layout_columns(
											col_widths = c(3,4,4),
											card(
												card_header("Polypharmacy Distribution"),
												vegawidgetOutput(ns("pp_distribution_plot")) #, width = 300)
											),
											card(
												card_header("Top 10 Long-term Conditions"),
												vegawidgetOutput(ns("top_ltcs_plot")) # , width = 300)
											),
											card(
												card_header("Top 10 Substances"),
												vegawidgetOutput(ns("top_substances_plot"))
											)
										)
									)
								),
								nav_panel(
									title = "Advanced",
									module_overview_ui(ns("overview_module")),
									module_prescription_explorer_ui(ns("prescription_explorer_module")),
									module_ltc_explorer_ui(ns("ltc_explorer_module"))
								)
							)
	)
}

module_ce_server <- function(id,
      selected_outcome,
      outcome_prescriptions,
      patient_data,
      outcome_data,
      ltc_data,
      pp_groups_data,
      acute_outcome_prescriptions) {
  moduleServer(
    id,
    function(input, output, session) {
      ns <- session$ns

      output$data_loaded <- reactive({
      	req(outcome_prescriptions())
      	!is.null(selected_outcome())
      })
      outputOptions(output, "data_loaded", suspendWhenHidden = FALSE)

      # Overview tab outputs
      output$value_box_cases <- renderText({
      	req(outcome_prescriptions())
      	prettyNum(uniqueN(outcome_prescriptions()$patid), big.mark = ",")
      })

      output$value_box_median_age <- renderText({
      	req(outcome_prescriptions())
      	prettyNum(median(round(outcome_prescriptions()$outcome_age/365.25, digits = 0)), big.mark = ",")
      })

      output$value_box_median_ltc <- renderText({
      	req(outcome_prescriptions())
      	prettyNum(median(outcome_prescriptions()$n_ltc), big.mark = ",")
      })

      output$value_box_median_pp <- renderText({
      	req(pp_groups_data())
      	prettyNum(median(pp_groups_data()$pp), big.mark = ",")
      })

      output$pp_distribution_plot <- renderVegawidget({
      	req(outcome_prescriptions(), pp_groups_data())
      	create_pp_distribution_plot(
      		pp_groups_data(),
      		width = "container",
      		title = "Polypharmacy Distribution"
      	)
      })

      output$top_ltcs_plot <- renderVegawidget({
      	req(ltc_data(), outcome_prescriptions())
      	create_top_ltcs_plot(
      		ltc_data(),
      		outcome_prescriptions(),
      		n_top = 10,
      		height = 400,
      		width = "container",
      		title = "Top 10 Long-term Conditions"
      	)
      })

      output$top_substances_plot <- renderVegawidget({
      	req(outcome_prescriptions())
      	create_top_substances_plot_ce(
      		outcome_prescriptions(),
      		n_top = 10,
      		height = 400,
      		width = "container",
      		title = "Top 10 Substances"
      	)
      })

    module_overview_server(
      id = "overview_module",
      outcome_prescriptions = outcome_prescriptions,
      patient_data = patient_data,
      outcome_ltc_events = outcome_data,
      selected_outcome = selected_outcome,
      pp_groups_data = pp_groups_data
    )

    module_prescription_explorer_server(
      id = "prescription_explorer_module",
      outcome_prescriptions = outcome_prescriptions,
      patient_data = patient_data,
      ltc_data = ltc_data,
      selected_outcome = selected_outcome,
      pp_groups_data = pp_groups_data,
      acute_outcome_prescriptions = acute_outcome_prescriptions
    )

    module_ltc_explorer_server(
      id = "ltc_explorer_module",
      outcome_prescriptions = outcome_prescriptions,
      ltc_data = ltc_data,
      patient_data = patient_data,
      pp_groups_data = pp_groups_data
    )
    }
  )
}