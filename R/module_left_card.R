#' @export
left_card_ui <- function(id) {
  ns <- NS(id)
  # card(card_header("Patient filtering"),
  list(
       textOutput(ns("selectedPatientInfo")),
       actionButton(ns("trigger_update"),label = "Update filter"),
       selectizeInput(ns("select_outcome"), label = "Select an outcome:", choices = NULL),
       div(dataTableOutput(ns("chapter_table"), height = "300px"), style = "height: 400px"),
       #checkboxInput(ns("mm_checkbox"), "2 LTCs before outcome?", value = TRUE),
  		 numericInput(ns("nltc_checkbox"), "Number of LTCs before outcome", value = 2, step = 1, min = 0),
       checkboxGroupInput(ns("sex_checkbox"), "Sex:", choices = c("Male", "Female"), inline = TRUE),
       checkboxGroupInput(ns("eth_checkbox"),
                          "Ethnicity:",
                          choices = c("White", "Black or Black British","Chinese or Other Group",
                                      "Unknown", "Asian or Asian British", "Mixed")),
       checkboxGroupInput(ns("imd_checkbox"), "IMD Quintile:", choices = 1:5, inline = TRUE)
  )
}

#' @export
left_card_server <- function(id, chapter_menu_data, patientFilter_r, selected_patient_number, total_patient_number, outcome_list) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    updateSelectizeInput(session, "select_outcome", choices = outcome_list)

    output$chapter_table <- renderDataTable({
      chapter_menu_data
    }, rownames = FALSE,
    options = list(
      order = list(list(1, 'asc')),
      dom = "<f<t>>",
      info = FALSE,
      paging = FALSE,
      scrollY = "300px"
    ))

    output$selectedPatientInfo <-  renderText({
      paste0("Filtered eligible patients: ", selected_patient_number(), "/", total_patient_number)
    })

    observeEvent(input$trigger_update,
                 {
                   selected_ltc_rowids <- input$chapter_table_rows_selected
                   selected_ltcs <- chapter_menu_data[selected_ltc_rowids,]$ltc
                   patientFilter_r$outcome <- input$select_outcome
                   patientFilter_r$selected_ltcs <- selected_ltcs
                   patientFilter_r$input_list <- list("min_nltc" = input$nltc_checkbox,
                                                      "sex" = input$sex_checkbox,
                                                      "eth_group" = input$eth_checkbox,
                                                      "imd_quintile" = input$imd_checkbox)
                 })

  })
}