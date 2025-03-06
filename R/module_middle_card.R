#' @export
middle_card_ui <- function(id) {
  ns <- NS(id)
  card(card_header("Polypharmacy viewer"),
       card_body(
         fillable = FALSE,
         #selectizeInput(ns("select_ltc_outcome"), "Select outcome:", choices = NULL, options = list(dropdownParent = "body")),
         div(style="min-height: 2rem",
         		fluidRow(
         			column(8, textOutput(ns("npats_outcome_prescriptions"))),
         			column(4, downloadButton(ns("download_code"), "Export Code")))),
         module_overview_ui(ns("overview_module")),
         module_outcome_explorer_ui(ns("outcome_explorer_module")),
         module_ltc_explorer_ui(ns("ltc_explorer_module")),
       )
  )
}

#' @export
middle_card_server <- function(id, prescription_data, ltc_data, patient_data, outcome_data, bnf_lookup, multimorbid_check, stored_queries,
															 polypharmacy_threshold, earliest_treatment_end) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    # observe({
    #   req(selected_ltc_list())
    #   updateSelectizeInput(session, "select_ltc_outcome", choices = c("", selected_ltc_list()))
    # })

    output$download_code <- downloadHandler(
    	filename = function() {
    		paste0("data_processing_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".R")
    	},
    	contentType = "text/r-script",
    	content = function(file) {
    		# Use the stored queries directly
    		export_sql_commands(
    			stored_queries$prescription,
    			stored_queries$ltc,
    			stored_queries$patient,
    			stored_queries$outcome,
    			file
    		)
    	}
    )

    outcome_ltc_events <- reactive({
      req(input$select_ltc_outcome)
      # query_obj <- buildLtcDatesQuery(dbQuoteString(db_pool, selected_outcome()))
      # poolWithTransaction(db_pool, function(con) {
      #   res <- dbGetQuery(con, query_obj)
      #   data.table$as.data.table(res)
      # })
      df <- ltc_data()
      df[term==input$select_ltc_outcome]
    })
    # Current method
    # outcome_prescriptions <- reactive({
    #   # Cache the reactive values to avoid multiple calls
    #   presc_df <- prescription_data()
    #   outcome_df <- outcome_data()
    #   ltcs <- ltc_data()
    #
    #   validate(need(nrow(ltcs) > 0, "No valid patients found"))
    #
    #   message("Computing outcome prescriptions")
    #
    #   # Pre-filter prescription data if possible
    #   # Assuming most outcomes won't be more than 30 days after stop_date
    #   # This reduces the data.table join size
    #   presc_df <- presc_df[stop_date >= (min(outcome_df$eventdate) - 30)]
    #
    #   # Perform the join with pre-computed indices
    #   diag_presc <- outcome_df[presc_df,
    #                            .(
    #                              patid,
    #                              eventdate = x.eventdate,
    #                              outcome_age = age_days,
    #                              substance,
    #                              start_date = i.start_date,
    #                              stop_date = i.stop_date,
    #                              duration
    #                            ),
    #                            on = .(patid, eventdate >= start_date),
    #                            nomatch = 0
    #   ][, .SD[length(unique(substance)) > 1],
    #                                        by = .(patid)]
    #   setkey(diag_presc,patid)
    #
    #   # Pre-filter LTC data if possible to reduce join size
    #   ltcs <- ltcs[patid %in% diag_presc$patid]
    #
    #   # Join with LTC data and filter for multiple conditions
    #   nltcs <- ltcs[diag_presc,
    #                .(patid, eventdate, outcome_age, age_days, term),
    #                on = .(patid, age_days < outcome_age),
    #                nomatch = 0][, .(n_ltc = length(unique(term))),
    #                             by = .(patid)]
    #   setkey(nltcs, patid)
    #   diag_presc[nltcs][n_ltc>1]
    # })

    outcome_prescriptions <- reactive({
      # Cache the reactive values and validate early
      presc_df <- prescription_data()
      outcome_df <- outcome_data()
      ltcs <- ltc_data()

      validate(need(nrow(ltcs) > 0, "No valid patients found"))

      message("Computing outcome prescriptions ", Sys.time())

      # Pre-compute date range limits to minimize data
      #browser()
      # min_outcome_date <- pmin(outcome_df$eventdate)
      # date_buffer <- 30
      # min_required_date <- min_outcome_date - date_buffer

      # More aggressive initial filtering of prescription data
      #presc_df <- presc_df[stop_date >= min_required_date &
      #                       patid %in% outcome_df$patid]

      presc_df <- presc_df[patid %in% outcome_df$patid]
      # Set keys before joining to improve performance
      # setkey(presc_df, patid, start_date)
      # setkey(outcome_df, patid, eventdate)

      # Two-stage join approach for better memory efficiency
      # First get relevant patient-dates
      potential_matches <- outcome_df[, .(patid, eventdate, age_days)]
      setkey(potential_matches, patid)

      # Then do the main join with pre-filtered data
      diag_presc <- potential_matches[presc_df,
                                      .(
                                        patid,
                                        eventdate = x.eventdate,
                                        outcome_age = age_days,
                                        substance,
                                        start_date,
                                        stop_date,
                                        duration
                                      ),
                                      on = .(patid, eventdate > start_date),
                                      nomatch = 0
      ]
      diag_presc <- diag_presc[stop_date >= eventdate - earliest_treatment_end()]
      # Use more efficient grouping operation
      multi_substance <- diag_presc[, if(uniqueN(substance) >= polypharmacy_threshold()) .SD,
                                    by = patid]

      # Early exit if no matching patients
      if (nrow(multi_substance) == 0) return(data.table())

      # Pre-filter LTC data more aggressively
      ltcs <- ltcs[patid %in% multi_substance$patid]
      # setkey(ltcs, patid, age_days)

      # Optimize the LTC join with pre-filtering
      nltcs <- ltcs[multi_substance,
                    .(patid, eventdate, outcome_age, age_days, term),
                    on = .(patid, age_days < outcome_age),
                    nomatch = 0
      ][, .(n_ltc = uniqueN(term)), by = patid]

      # Final filtering
      setkey(nltcs, patid)
      message("Finished ", Sys.time())
      # Keep only those patients that are already multimorbid before outcome
      # Maybe this could be an option on the interface?
      min_ltc <- 0
      if (multimorbid_check) min_ltc <- 1
      multi_substance[nltcs[n_ltc > min_ltc]]
    })

    output$npats_outcome_prescriptions <- renderText({
      paste0("Number of patients with polypharmacy within 30 days of outcome: ", length(unique(outcome_prescriptions()$patid)))
    })
    outputOptions(output, "npats_outcome_prescriptions", priority = 20)

    selected_outcome <- reactive({
      unique(outcome_data()$term)
    })

    pp_groups_data <- reactive({
      df <- outcome_prescriptions()
      pat_data <- patient_data()
      df <- df[patid %in% pat_data$patid, list(pp = length(unique(substance))), .(patid)]
      setkey(df, patid)
      dt <- create_value_groups(df,
                          breaks = c(2,4,6,8,11),
                          value_col = "pp",
                          group_col = "pp_group")
      dt
    })

    module_overview_server("overview_module", outcome_prescriptions, patient_data, outcome_data, selected_outcome, pp_groups_data)
    module_outcome_explorer_server("outcome_explorer_module", outcome_prescriptions, patient_data, selected_outcome, bnf_lookup, pp_groups_data)
    module_ltc_explorer_server("ltc_explorer_module", outcome_prescriptions, ltc_data, patient_data)
  })
}