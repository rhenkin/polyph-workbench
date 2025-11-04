#' @export
middle_card_ui <- function(id) {
  ns <- NS(id)
  navset_pill(
    id = ns("middle_tab"),
    module_ce_ui(ns("ce_module")),
    module_ccm_ui(ns("ccm_module")),
    module_cca_ui(ns("cca_module"))
  )
}

#' @export
middle_card_server <- function(id, prescription_data, ltc_data, patient_data,
                              outcome_data, acute_presc_df, min_nltc,
                              stored_queries, polypharmacy_threshold, earliest_treatment_end) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    # Data layer - reactive data calculations
    queried_terms <- reactive({ stored_queries$ltc$terms })

    outcome_prescriptions <- reactive({
    	req(prescription_data(), outcome_data(), ltc_data())
      calculate_outcome_prescriptions(
        presc_df = prescription_data(),
        outcome_df = outcome_data(),
        ltcs = ltc_data(),
        polypharmacy_threshold = polypharmacy_threshold(),
        earliest_treatment_end = earliest_treatment_end(),
        min_nltc = min_nltc,
        queried_terms = queried_terms()
      )
    })

    acute_outcome_prescriptions <- reactive({
      calculate_acute_outcome_prescriptions(
        outcome_df = outcome_data(),
        acute_df = acute_presc_df(),
        outcome_prescriptions = outcome_prescriptions(),
        earliest_treatment_end = earliest_treatment_end()
      )
    })

    pp_groups_data <- reactive({
      req(outcome_prescriptions())
      calculate_pp_groups(outcome_prescriptions())
    })

    selected_outcome <- reactive({
      req(outcome_data())
      unique(outcome_data()$term)
    })

    # Module servers
    module_ce_server(
      id = "ce_module",
      selected_outcome = selected_outcome,
      outcome_prescriptions = outcome_prescriptions,
      patient_data = patient_data,
      outcome_data = outcome_data,
      ltc_data = ltc_data,
      pp_groups_data = pp_groups_data,
      acute_outcome_prescriptions = acute_outcome_prescriptions
    )

    matched_data <- module_ccm_server(
      id = "ccm_module",
      patient_data = patient_data,
      outcome_prescriptions = outcome_prescriptions,
      ltc_data = ltc_data,
      study_dir = "studies"
    )

    module_cca_server(
      id = "cca_module",
      prepared_study_data_r = matched_data$prepared_study_data_r
    )

    return(list(
      seltab = reactive(input$middle_tab),
      outcome_prescriptions = outcome_prescriptions
    ))
  })
}