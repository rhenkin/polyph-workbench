#' @export
middle_card_ui <- function(id) {
  ns <- NS(id)
  navset_pill(
  	nav_panel("Polypharmacy viewer",
  	card(
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
  	),
  	module_psm_ui(ns("psm_module")),
  	module_cca_ui(ns("cca_module"))
  )
}

#' @export
middle_card_server <- function(id, db_pool, prescription_data, ltc_data, patient_data, outcome_data, bnf_lookup, chapter_menu_data, min_nltc, stored_queries,
															 polypharmacy_threshold, earliest_treatment_end, acute_presc_df, gold_patient, gold_ltc, gold_cp) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

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

    queried_terms <- reactive({ stored_queries$ltc$terms })

    outcome_ltc_events <- reactive({
      req(input$select_ltc_outcome)
      df <- ltc_data()
      df[term==input$select_ltc_outcome]
    })

    outcome_prescriptions <- reactive({
      # Cache the reactive values and validate early
      presc_df <- prescription_data()
      outcome_df <- outcome_data()
      ltcs <- ltc_data()

      validate(need(nrow(ltcs) > 0, "No valid patients found"))

      message("Computing outcome prescriptions ", Sys.time())
      presc_df <- presc_df[patid %in% outcome_df$patid]

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

			diag_presc <- diag_presc[, list(start_date = first(start_date),
      												 		 stop_date = last(stop_date),
      												 		 duration = sum(duration)),
      												 .(patid,substance,eventdate,outcome_age)]

			diag_presc[stop_date >= eventdate, duration := duration - (stop_date-eventdate)]

      setkey(diag_presc, patid, start_date)

      # Calculate polypharmacy count - for each unique start date, all drugs get the cumulative count
      diag_presc[, drugs_count := 1:.N, by = patid]

      # Handle drugs with same start date - they should all get the highest count for that date
      diag_presc[, polyph_number := max(drugs_count), by = .(patid, start_date)]

      # Use more efficient grouping operation
      diag_presc[, N := .N, patid]
      multi_substance <- diag_presc[N >= polypharmacy_threshold()]

      # Early exit if no matching patients
      if (nrow(multi_substance) == 0) return(data.table())

      # Pre-filter LTC data more aggressively
      ltcs <- ltcs[patid %in% multi_substance$patid]

      # Optimize the LTC join with pre-filtering
      unique_patients <- unique(multi_substance[, .(patid, eventdate, outcome_age)])


      first_ltc <- ltcs[unique_patients,
      									.(patid, eventdate = i.eventdate, outcome_age = i.outcome_age, age_days = x.age_days, term = x.term),
      									on = .(patid, age_days < outcome_age),
      									nomatch = 0
      ]
      if (!is.null(queried_terms())) {
	      pats_with_terms <- first_ltc[,all(queried_terms() %in% term), patid][V1==TRUE,patid]
	      first_ltc <- first_ltc[patid %in% pats_with_terms]
      }
      nltcs <- first_ltc[, list(n_ltc = .N), patid]

      # Final filtering
      setkey(nltcs, patid)
      message("Finished ", Sys.time())

      final_df <- multi_substance[nltcs[n_ltc >= min_nltc]]
      final_df <- create_value_groups(final_df, breaks = c(2,5,10), right =FALSE, value_col = "n_ltc", group_col = "mltc_group")
      return(final_df)
    })

    acute_outcome_prescriptions <- reactive({
    	outcome_df <- outcome_data()
    	acute_df <- acute_presc_df()
    	valid_patids <- unique(outcome_prescriptions()$patid)
    	acute_df <- acute_df[patid %in% valid_patids]
    	merged <- outcome_df[acute_df, .(patid, eventdate=x.eventdate, substance,start_date,outcome_age=age_days), on = .(patid, eventdate > start_date), nomatch=0]
    	merged <- merged[start_date >= eventdate - earliest_treatment_end()]
    	setkey(merged, patid)
			merged
    })

    output$npats_outcome_prescriptions <- renderText({
      paste0("Number of patients with polypharmacy within 30 days of outcome: ", length(unique(outcome_prescriptions()$patid)))
    })
    outputOptions(output, "npats_outcome_prescriptions", priority = 20)

    selected_outcome <- reactive({
    	req(outcome_data())
      unique(outcome_data()$term)
    })

    pp_groups_data <- reactive({
    	req(outcome_prescriptions())
      df <- outcome_prescriptions()
      pat_data <- patient_data()
      grouped_df <- df[,list(pp=.N),patid]
      dt <- create_value_groups(grouped_df,
                          breaks = c(2, 5, 10),
      										right = FALSE,
                          value_col = "pp",
                          group_col = "pp_group")
      dt
    })

    module_overview_server("overview_module", outcome_prescriptions, patient_data, outcome_data, selected_outcome, pp_groups_data)
    module_outcome_explorer_server("outcome_explorer_module", outcome_prescriptions, patient_data, ltc_data, selected_outcome, bnf_lookup, chapter_menu_data, pp_groups_data, acute_outcome_prescriptions)
    module_ltc_explorer_server("ltc_explorer_module", outcome_prescriptions, ltc_data, patient_data, pp_groups_data)
    matched_data <- module_psm_server(
    	id = "psm_module",
    	patient_data = patient_data,
    	outcome_prescriptions = outcome_prescriptions,
    	ltc_data = ltc_data,
    	gold_patient = gold_patient,
    	gold_ltc = gold_ltc,
    	gold_cp = gold_cp,
    	acute_presc_df,
    	study_dir = "studies"  # optional, defaults to "studies"
    )
    module_cca_server(id = "cca_module",
    									matched_data = matched_data,
    									gold_patient = gold_patient,
    									gold_ltc = gold_ltc,
    									gold_cp = gold_cp,
    									acute_presc_df
    									)
  })
}