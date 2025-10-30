# Case explorer module
# First tab in the middle card
# Calls sub-modules: overview, prescription explorer, LTC explorer
# Shows number of patients with both outcome and prescription data
module_ce_ui <- function(id) {
    ns <- NS(id)    
    nav_panel("Outcome explorer", value = "outcome_explorer",
      card(
        fillable = FALSE,
        module_overview_ui(ns("overview_module")),
        module_prescription_explorer_ui(ns("prescription_explorer_module")),
        module_ltc_explorer_ui(ns("ltc_explorer_module"))
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