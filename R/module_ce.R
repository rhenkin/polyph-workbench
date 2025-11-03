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
									layout_columns(
										col_widths = c(2,5,5),
										value_box(
											showcase_layout = showcase_left_center(width = "200px", max_height = "150px"),
											height = "150px",
											title = "Cases",
											value = textOutput(ns("value_box_cases")),
											theme = "red",
											showcase = bs_icon("people-fill")
										),
										card(
											card_header("Polypharmacy Distribution"),
											vegawidgetOutput(ns("pp_distribution_plot"), width = 300)
										),
										card(
											card_header("Top 10 Long-term Conditions"),
											vegawidgetOutput(ns("top_ltcs_plot"), , width = 300)
										)),

										layout_columns(
											col_widths = c(4, 8),
											card(
												card_header("Top 10 Substances"),
												vegawidgetOutput(ns("top_substances_plot"))
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

      # Overview tab outputs
      output$value_box_cases <- renderText({
      	req(outcome_prescriptions())
      	prettyNum(uniqueN(outcome_prescriptions()$patid), big.mark = ",")
      })

      output$pp_distribution_plot <- renderVegawidget({
      	req(outcome_prescriptions(), pp_groups_data())
      	create_pp_distribution_plot(
      		pp_groups_data(),
      		title = "Polypharmacy Distribution"
      	)
      })

      output$top_ltcs_plot <- renderVegawidget({
      	req(ltc_data(), outcome_prescriptions())
      	create_top_ltcs_plot(
      		ltc_data(),
      		outcome_prescriptions(),
      		n_top = 10,
      		title = "Top 10 Long-term Conditions"
      	)
      })

      output$top_substances_plot <- renderVegawidget({
      	req(outcome_prescriptions())
      	create_top_substances_plot_ce(
      		outcome_prescriptions(),
      		n_top = 10,
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