module_presc_tto_ui <- function(id) {
  ns <- NS(id)
  card(
    selectizeInput(ns("select_substance"), label = "Select substance:", choices = character(0),
                   options = list(dropdownParent = "body")),
    textOutput(ns("survres"))
  )
}

module_presc_tto_server <- function(id, subst_pp_df) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    observe({
      req(subst_pp_df())
      updateSelectizeInput(
        session,
        "select_substance",
        choices = sort(unique(subst_pp_df()$BNF_Chemical_Substance)),
        selected = character(0),
        server = TRUE
      )
    })

    selected_subst_df <- reactive({
      req(input$select_substance)
      browser()
      subst_pp_df()[BNF_Chemical_Substance==input$select_substance]
    })
    observe({
      req(input$select_substance)
      browser()
      df[order(patid,-start_date)][, .SD[1], by = patid][,mean(duration),BNF_Chemical_Substance][order(-V1)]

    })

    #ltcs <- ltc_data()
  })
}