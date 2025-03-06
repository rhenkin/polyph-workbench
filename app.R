dyn.load('/share/apps/rocky9/spack/apps/linux-rocky9-x86_64_v4/gcc-12.2.0/postgresql/15.2-mf5fwdl/lib/libpq.so.5')
library(bslib)
library(data.table)
library(DBI)
library(DT)
library(MASS)
library(pool)
library(RPostgres)
library(shiny)
library(vegawidget)
library(htmltools)
library(igraph)
library(visNetwork)
library(gt)
library(reactable)

#setDTthreads(6)

# Create a pool object for PostgreSQL
db_pool <- dbPool(
  drv = Postgres(),
  dbname = "iphs_finer_cprd",
  host = "db1",
  port = 5432,
  user = "iphs_finer_cprd",
  password = "ies4phahph3A"
)

bnf_lookup <- fread("new_bnf_lkp.csv")

#' @export
ui <-
  page_fluid(
  	useBusyIndicators(),
    theme = bs_theme(preset = "lumen"),
    layout_column_wrap(
      width = NULL, fill = TRUE,
      style = css(grid_template_columns = "1fr 3fr 1fr"),
      left_card_ui("patient_filter_module"),
      middle_card_ui("polyph_module"),
      right_card_ui("drug_filter_module")
    )
  )

#' @export
server <- function(input, output, session) {
    patient_df <- poolWithTransaction(db_pool, function(con) {
      res <- dbGetQuery(con, sprintf("SELECT patid,gender as sex,eth_group,imd_quintile FROM gold_cp_patient"))
      #as.numeric(res[1,1])
      as.data.table(res)
    })
    total_patient_number <- nrow(patient_df)
    selected_patient_number <- reactiveVal(0)

    outcome_list <- poolWithTransaction(db_pool, function(con) {
      res <- as.data.table(
        dbGetQuery(con, "SELECT DISTINCT outcome FROM outcomes ORDER BY outcome ASC")
      )
      res[,1]
    })

    stored_queries <- reactiveValues()

    patientFilter_r <- reactiveValues()
    # Run left side menu code
    # From module: patientFilter_r
    # To module: patient numbers and list out of outcomes
    left_card_server(
      "patient_filter_module",
      patientFilter_r,
      selected_patient_number,
      total_patient_number,
      outcome_list
    )

    patient_data <- reactive({
      req(patientFilter_r$outcome)
      message("Querying patient data ", Sys.time())
      query_obj <- buildPatientQuery(patientFilter_r$selected_ltcs,
                                     patientFilter_r$input_list$eth_group,
                                     patientFilter_r$input_list$sex,
                                     patientFilter_r$input_list$imd_quintile,
                                     patientFilter_r$outcome)

      stored_queries$patient <- query_obj

      poolWithTransaction(db_pool, function(con) {
        res <- as.data.table(
          dbGetQuery(con, query_obj$query, params = unlist(query_obj$params))
        )
        res[,dob := as.IDate(dob)]
        res$patid <- bit64::as.integer64(res$patid)
        setkey(res,patid)
        res[imd_quintile != "NA"]
      })
    }) |> bindCache(
      patientFilter_r$selected_ltcs,
      patientFilter_r$input_list$eth_group,
      patientFilter_r$input_list$sex,
      patientFilter_r$input_list$imd_quintile,
      patientFilter_r$outcome
    )

    selected_patient_number <- reactive({
      nrow(patient_data())
    })

    outcome_data <- reactive({
      req(patientFilter_r$outcome)
      message("Querying outcome data ", Sys.time())
      query_obj <- buildOutcomeQuery(patientFilter_r$outcome)
      stored_queries$outcome <- query_obj
      poolWithTransaction(db_pool, function(con) {
        res <- as.data.table(
          dbGetQuery(con, query_obj$query, params = unlist(query_obj$params))
        )
        res[, eventdate := as.IDate(eventdate)]
        setkey(res,patid,eventdate)
        res
      })
    }) |> bindCache(
      patientFilter_r$outcome
    )

    prescription_data <- reactive({
      req(patientFilter_r$outcome)
      message("Querying prescription data ", Sys.time())
      query_obj <- buildPrescriptionQuery(patientFilter_r$selected_ltcs,
                                          NULL,
                                          patientFilter_r$input_list$eth_group,
                                          patientFilter_r$input_list$sex,
                                          patientFilter_r$input_list$imd_quintile,
                                          patientFilter_r$outcome)
      stored_queries$prescription <- query_obj
      poolWithTransaction(db_pool, function(con) {
        res <- as.data.table(
          dbGetQuery(con, query_obj$query, params = unlist(query_obj$params))
        )
        res[, start_date := as.IDate(start_date)]
        res[, stop_date := as.IDate(stop_date)]
        setkey(res,patid,start_date)
        res
      })
    }) |> bindCache(
      patientFilter_r$selected_ltcs,
      patientFilter_r$input_list$eth_group,
      patientFilter_r$input_list$sex,
      patientFilter_r$input_list$imd_quintile,
      patientFilter_r$outcome
    )

    right_card_inputs <- right_card_server(id = "drug_filter_module", bnf_lookup)
    filtered_prescription_data <- reactive({
      req(prescription_data())
      filtered_substances <-
        get_filtered_substances(right_card_inputs$included, right_card_inputs$excluded, bnf_lookup)
      tmp_df <- prescription_data()[substance %in% filtered_substances]
      tmp_df[duration >= right_card_inputs$minimum_duration()]
    })

    ltc_data <- reactive({
      req(patientFilter_r$outcome)
      message("Querying LTC data ", Sys.time())
      query_obj <- buildLtcQuery(patientFilter_r$selected_ltcs,
                                 NULL,
                                 patientFilter_r$input_list$eth_group,
                                 patientFilter_r$input_list$sex,
                                 patientFilter_r$input_list$imd_quintile,
                                 patientFilter_r$outcome)
      stored_queries$ltc <- query_obj
      poolWithTransaction(db_pool, function(con) {
        res <- as.data.table(
          dbGetQuery(con, query_obj$query, params = unlist(query_obj$params))
        )
        res[,eventdate := as.IDate(eventdate)]
        res$patid <- bit64::as.integer64(res$patid)
        setkey(res, patid, age_days)
        res
      })
    }) |> bindCache(
      patientFilter_r$selected_ltcs,
      patientFilter_r$input_list$eth_group,
      patientFilter_r$input_list$sex,
      patientFilter_r$input_list$imd_quintile,
      patientFilter_r$outcome
    )

    # Middle module
    middle_card_server(id = "polyph_module",
                       filtered_prescription_data,
                       ltc_data,
                       patient_data,
                       outcome_data,
                       bnf_lookup,
                       patientFilter_r$input_list$multmorbid,
    									 stored_queries,
    									 right_card_inputs$polypharmacy_threshold,
    									 right_card_inputs$earliest_treatment_end)

    # Close the pool when the session ends
    session$onSessionEnded(function() {
      poolClose(db_pool)
    })
}

shiny::shinyApp(ui, server)
