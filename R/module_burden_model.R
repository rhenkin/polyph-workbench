#' @export
module_burden_model_ui <- function(id) {
  ns <- NS(id)
  card(#card_header("Drug association with polypharmacy burden"),
    #checkboxInput(ns("control_demog"), "Control for demographics?"),
    vegawidgetOutput(ns("drug_pp_plot")),
    downloadButton(ns("download_irr_table"), "Download table")
  )
}

#' @export
module_burden_model_server <- function(id, subst_pp_df, patient_data, subst_frequency) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    top_50_substance <- reactive({
      validate(need(nrow(subst_pp_df()) > 0, "Waiting for data to load"))
      
      subst_pp_df()[substance %in% subst_frequency()[1:50, substance]]
    })
    
    # data_prep <- reactive({
    #   prepare_substance_data(top_50_substance(), patient_data())
    # })
    # results <- reactiveVal()
    # observe({
    #   req(data_prep(), top_50_substance())
    #   
    #   prep_data <- data_prep()
    #   model_dt <- prep_data$model_dt
    #   substances <- prep_data$substances
    #   curr_results <- results()
    #   # Process each substance and update the plot
    #   for (substance in substances) {
    #     result <- calculate_single_substance_irr(substance, model_dt)
    #     new_results <- rbindlist(list(curr_results, result))
    #     results(new_results)
    #   }
    # })
    # vw_shiny_set_data(outputId = "drug_pp_plot", 
    #                               name = "source",
    #                               value = results())
    # 
    
    irr_data <- reactive({
      calculate_substance_irr(top_50_substance(), patient_data() , TRUE) #input$control_demog)})
    })

    output$drug_pp_plot <- renderVegawidget({
      dot_plot(irr_data()[pct_sample > 0.01])
    })
    
    output$download_irr_table <- downloadHandler(
      filename = function() {
        "incidence_table.csv"
      },
      content = function(file) {
        utils::write.csv(irr_data()[pct_sample > 0.01], file, row.names = FALSE)
      }, contentType = "text/csv"
    )
    
  })
}