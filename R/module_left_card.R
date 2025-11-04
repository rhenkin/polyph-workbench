#' @export
left_card_ui <- function(id) {
  ns <- NS(id)
  # card(card_header("Patient filtering"),
  list(
  		 div("Choose outcome and patient criteria and click \"Start analysis\" to begin"),
       textOutput(ns("selectedPatientInfo")),
       actionButton(ns("trigger_update"),label = "Start analysis"),
       selectizeInput(ns("select_outcome"), label = "Select an outcome:", choices = NULL),
       div("Select one or more conditions to define inclusion criteria:"),
  		 div(dataTableOutput(ns("chapter_table"), height = "300px"), style = "height: 400px"),
       #checkboxInput(ns("mm_checkbox"), "2 LTCs before outcome?", value = TRUE),
  		 numericInput(ns("nltc_checkbox"), "Number of conditions before outcome:", value = 2, step = 1, min = 0),
       div("Select one or more demographic categories below to define inclusion criteria:"),
  		 checkboxGroupInput(ns("sex_checkbox"), "Sex:", choices = c("Male", "Female"), inline = TRUE),
       checkboxGroupInput(ns("eth_checkbox"),
                          "Ethnicity:",
                          choices = c("White", "Black or Black British","Chinese or Other Group",
                                      "Unknown", "Asian or Asian British", "Mixed")),
       checkboxGroupInput(ns("imd_checkbox"), "IMD Quintile:", choices = 1:5, inline = TRUE)
  )
}

#' @export
left_card_server <- function(id, patientFilter_r, outcome_list, analyzed_patient_count) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    updateSelectizeInput(session, "select_outcome", choices = outcome_list)

    output$chapter_table <- renderDataTable({
      datatable(ltc_chapters,
      	class = list(stripe = FALSE),
		    rownames = FALSE,
		    options = list(
		      order = list(list(1, 'asc')),
		      dom = "<f<t>>",
		      info = FALSE,
		      paging = FALSE,
		      scrollY = "300px"
    )) })
	outputOptions(session$output, "chapter_table", suspendWhenHidden = TRUE)

    output$selectedPatientInfo <-  renderText({
      # req(analyzed_patient_count() > 0)
      paste0("Selected patients: ", prettyNum(analyzed_patient_count(), big.mark = ","))
    }) |> bindEvent(analyzed_patient_count(), ignoreInit = TRUE)

    observeEvent(input$trigger_update,
                 {
                   selected_ltc_rowids <- input$chapter_table_rows_selected
                   selected_ltcs <- ltc_chapters[selected_ltc_rowids,]$ltc
                   patientFilter_r$outcome <- input$select_outcome
                   patientFilter_r$selected_ltcs <- selected_ltcs
                   patientFilter_r$input_list <- list("min_nltc" = input$nltc_checkbox,
                                                      "sex" = input$sex_checkbox,
                                                      "eth_group" = input$eth_checkbox,
                                                      "imd_quintile" = input$imd_checkbox)
                 })

  })
}