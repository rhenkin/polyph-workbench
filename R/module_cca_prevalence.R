module_cca_prevalence_ui <- function(id) {
  ns <- NS(id)
  
  accordion(
    open = FALSE,
    accordion_panel(
      title = "Long-term conditions",
      value = "ltc_prev_tables",
      navset_tab(
        nav_panel(
          "Group prevalences",
          virtualSelectInput(
            ns("ltc_freq_strat_variable"),
            label = "Select a subset (optional):",
            choices = NULL,
            autoSelectFirstOption = FALSE,
            search = FALSE,
            showValueAsTags = TRUE,
            disableOptionGroupCheckbox = TRUE,
            multiple = FALSE,
            dropboxWrapper = "body"
          ),
          dataTableOutput(ns("ltc_freq_table"))
        )
      )
    ),
    accordion_panel(
      title = "Prescriptions",
      value = "presc_prev_tables",
      navset_tab(
        nav_panel(
          "Group prevalences",
          virtualSelectInput(
            ns("presc_freq_strat_variable"),
            label = "Select a subset (optional):",
            choices = NULL,
            autoSelectFirstOption = FALSE,
            search = FALSE,
            showValueAsTags = TRUE,
            disableOptionGroupCheckbox = TRUE,
            multiple = FALSE,
            dropboxWrapper = "body"
          ),
          dataTableOutput(ns("presc_freq_table"))
        )
      )
    ),
    accordion_panel(
      title = "Prescription prevalence per LTC",
      value = "presc_ltc_prev",
      card(
        full_screen = TRUE,
        height = "60em",
        uiOutput(ns("ltc_dropdown_ui")),
        reactableOutput(ns("presc_by_ltc"), height = "50em")
      )
    )
  )
}

module_cca_prevalence_server <- function(id, patient_data_r, prescriptions_r, 
                                        ltcs_r) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    
    # Update stratification choices
    observe({
      req(patient_data_r())
      choices <- create_stratification_choices(patient_data_r())
      updateVirtualSelect("ltc_freq_strat_variable", choices = choices)
      updateVirtualSelect("presc_freq_strat_variable", choices = choices)
    })
    
    # LTC frequency with optional stratification
    ltc_freq_df <- reactive({
      req(ltcs_r())
      ltcs <- ltcs_r()
      patient_data <- patient_data_r()
      strat_var <- input$ltc_freq_strat_variable
      
      if (strat_var != "") {
        parts <- strsplit(strat_var, "#")[[1]]
        column_name <- parts[1]
        filter_value <- parts[2]
        selected_patids <- patient_data[get(column_name) == filter_value, patid]
        ltcs <- ltcs[patid %in% selected_patids]
      }
      
      calculate_frequency_stats(ltcs, "term")
    })
    
    # Prescription frequency with optional stratification
    presc_freq_df <- reactive({
      req(prescriptions_r())
      presc <- prescriptions_r()
      patient_data <- patient_data_r()
      strat_var <- input$presc_freq_strat_variable
      
      if (strat_var != "") {
        parts <- strsplit(strat_var, "#")[[1]]
        column_name <- parts[1]
        filter_value <- parts[2]
        selected_patids <- patient_data[get(column_name) == filter_value, patid]
        presc <- presc[patid %in% selected_patids]
      }
      
      calculate_frequency_stats(presc, "substance")
    })
    
    # Render tables
    output$ltc_freq_table <- renderDataTable({
      table_data <- create_prevalence_ratio_table(ltc_freq_df(), "term")
      colnames(table_data) <- c("Condition", "Cases (%)", "Controls (%)", "Case/Control Ratio")
      datatable(table_data[order(-`Case/Control Ratio`)],
                class = list(stripe = FALSE), rownames = FALSE)
    })
    
    output$presc_freq_table <- renderDataTable({
      table_data <- create_prevalence_ratio_table(presc_freq_df(), "substance")
      colnames(table_data) <- c("Substance", "Cases (%)", "Controls (%)", "Case/Control Ratio")
      datatable(table_data[order(-`Case/Control Ratio`)],
                class = list(stripe = FALSE), rownames = FALSE)
    })
    
    # LTC dropdown
    output$ltc_dropdown_ui <- renderUI({
      req(ltcs_r())
      ltcs <- ltcs_r()
      
      virtualSelectInput(
        ns("ltc_dropdown"),
        label = "Select 1 or more LTCs:",
        choices = with(
          ltc_chapters[ltc %in% unique(ltcs$term)],
          split(ltc, body_system)
        ),
        multiple = TRUE,
        search = TRUE
      )
    })
    
    # Prescription by LTC table
    output$presc_by_ltc <- renderReactable({
      req(ltcs_r(), input$ltc_dropdown)
      
      result <- calculate_presc_by_ltc(
        prescriptions_r(),
        ltcs_r(),
        input$ltc_dropdown
      )
      
      reactable(result, showPageInfo = FALSE, defaultPageSize = 15)
    })
  })
}