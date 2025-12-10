#' Long-term Conditions Explorer UI Module
#'
#' Creates the user interface for exploring long-term conditions (LTCs) prevalence
#' within a cohort, including overall prevalence, demographic breakdowns, and
#' prevalence by prescribed drug.
#'
#' @param id Character string. The module namespace ID.
#'
#' @return A Shiny UI card containing accordion panels for:
#'   - LTC prevalence tables (overall and by demographics)
#'   - LTC prevalence analysis by drug
#'
#' @details
#' The UI provides two main sections:
#' 1. Prevalence tables: Shows LTC frequency in the cohort and across demographic groups
#' 2. LTC prevalence per drug: Allows filtering by specific substances to analyze
#'    LTC prevalence among patients prescribed those drugs
module_ltc_explorer_ui <- function(id) {
    ns <- NS(id)
    card(
        card_header("Long-term conditions"),
        accordion(
            open = FALSE,
            # First accordion panel: LTC prevalence tables
            accordion_panel(
                title = "Prevalence table",
                icon = bs_icon("table"),
                navset_tab(
                    # Tab 1: Overall cases prevalence
                    nav_panel("Prevalence", dataTableOutput(ns("ltc_freq_table"))),
                    # Tab 2: Prevalence stratified by demographic variables
                    nav_panel(
                        "Prevalence across demographics",
                        selectizeInput(ns("select_ltcdemog_freq_var"), label = "Select variable:", choices = c("Sex" = "sex", "Ethnic group" = "eth_group", "IMD quintile" = "imd_quintile", "PP" = "pp_group", "# LTCs" = "mltc_group")),
                        dataTableOutput(ns("demog_ltc_freq_table"))
                    )
                )
            ),
            # Second accordion panel: LTC prevalence by drug
            accordion_panel(
                "LTC prevalence per drug",
                icon = bs_icon("prescription2"),
                card(
                    full_screen = TRUE,
                    height = "60em",
                    fluidRow(
                        column(
                            6,
                            # Dropdown for selecting multiple substances
                            selectizeInput(ns("drug_dropdown"),
                                label = "Select 1 or more substance:", choices = NULL, multiple = TRUE,
                                options = list(
                                    splitOn = ";"
                                )
                            ),
                            # Text input for pasting a list of substances
                            textInput(ns("paste_substances"),
                                label = "Or paste semicolon-separated list:",
                                placeholder = "Paracetamol; Aspirin"
                            )
                        ),
                        # Age filter for patients
                        column(5, numericInput(ns("outcome_age_filter"), label = "Max age at outcome:", value = 100, min = 16, max = 100)),
                        # Display count of patients within age range
                        column(1, textOutput(ns("selected_pats")))
                    ),
                    # Table showing LTC prevalence for selected drugs
                    dataTableOutput(ns("ltc_by_presc"))
                )
            )
        )
    )
}

#' Long-term Conditions Explorer Server Module
#'
#' Server logic for the LTC explorer module. Handles data filtering, calculation
#' of LTC prevalence statistics, and rendering of interactive tables.
#'
#' @param id Character string. The module namespace ID.
#' @param outcome_prescriptions Reactive expression returning a data.table of
#'   prescription records with at least columns: patid, substance, outcome_age
#' @param ltc_data Reactive expression returning a data.table of long-term
#'   condition records with at least column: patid
#' @param patient_data Reactive expression returning a data.table of patient
#'   demographic data
#' @param pp_groups_data Reactive expression returning a data.table of polypharmacy
#'   group classifications
#'
#' @return Server logic (no return value)
#'
#' @details
#' The server module performs the following operations:
#' - Filters LTC data to only patients in the outcome cohort
#' - Calculates overall and demographic-specific LTC prevalence
#' - Updates drug selection dropdowns dynamically based on prescription data
#' - Filters prescriptions by age and calculates LTC prevalence by drug
module_ltc_explorer_server <- function(id, outcome_prescriptions, ltc_data, patient_data, pp_groups_data) {
    moduleServer(id, function(input, output, session) {
        ns <- session$ns

        # Reactive: Filter LTC data to only include patients in the outcome cohort
        valid_ltcs <- reactive({
            outcome_df <- outcome_prescriptions()
            ltcs <- ltc_data()
            ltcs[patid %in% outcome_df$patid]
        })

        # Reactive: Calculate overall LTC frequency for the cohort
        ltc_freq_df <- reactive({
            ltcs <- valid_ltcs()
            cases_freq <- calculate_ltc_frequency(ltcs)
            overall_freq <- calculate_ltc_frequency(gold_ltc)
            result <- merge(
            	cases_freq[, .(term, N, pct_cases = pct_total)],
            	overall_freq[, .(term, pct_overall = pct_total)],
            	by = "term",
            	all.x = TRUE
            )
						result[, N := prettyNum(N, big.mark = ",")]
            setnames(result, c("pct_cases", "pct_overall"),
            				 c("Cases (%)", "Overall (%)"))

            result

        })

        # Output: Render the overall LTC frequency table
        output$ltc_freq_table <- renderDataTable({
            ltc_freq_df()
        }, rownames = FALSE)

        # Output: Render LTC frequency table stratified by demographic variables
        output$demog_ltc_freq_table <- renderDataTable(
            {
                req(input$select_ltcdemog_freq_var)
                df <- valid_ltcs()
                outcome_df <- outcome_prescriptions()
                demog_var <- input$select_ltcdemog_freq_var

                # Calculate LTC prevalence broken down by the selected demographic variable
                calculate_demographic_ltc_frequency(
                    ltcs = df,
                    patient_data = patient_data(),
                    pp_groups_data = pp_groups_data(),
                    outcome_df = outcome_df,
                    demog_var = demog_var,
                    ltc_freq_df = ltc_freq_df()
                )
            },
            rownames = FALSE
        )

        # Observer: Dynamically update drug dropdown with available substances from prescriptions
        observe({
            prescs <- outcome_prescriptions()
            updateSelectizeInput(
                session = session,
                inputId = "drug_dropdown",
                choices = c("", sort(unique(prescs$substance))),
                options = list(delimiter = "; "),
                server = TRUE
            )
        })

        # Observer: Handle pasting of semicolon-separated substance list
        # Validates pasted substances and updates the drug dropdown accordingly
        observeEvent(input$paste_substances, {
            req(outcome_prescriptions())
            if (input$paste_substances != "") {
                prescs <- outcome_prescriptions()
                # Parse the pasted list (semicolon-separated)
                pasted_items <- trimws(strsplit(input$paste_substances, ";")[[1]])
                # Filter to only valid substances that exist in the prescription data
                valid_items <- pasted_items[pasted_items %in% sort(unique(prescs$substance))]

                # Update dropdown with validated substances
                updateSelectizeInput(session, "drug_dropdown", selected = valid_items)
                # Clear the text input after processing
                updateTextInput(session, "paste_substances", value = "")
            }
        })

        # Reactive: Filter prescriptions by maximum age at outcome
        # Age is converted from years to days for comparison (365.25 days per year)
        age_filtered_prescriptions <- reactive({
            filter_by_age(outcome_prescriptions(), input$outcome_age_filter)
        })

        # Output: Display count of unique patients within the age filter range
        output$selected_pats <- renderText({
            format_patient_count(age_filtered_prescriptions())
        })

        # Output: Render table showing LTC prevalence for selected drug(s)
        output$ltc_by_presc <- renderDataTable(
            {
                req(input$drug_dropdown)
                input$drug_dropdown
                ltcs <- valid_ltcs()
                prescriptions <- age_filtered_prescriptions()
                # Further filter LTCs to patients in the age-filtered prescription cohort
                ltcs <- ltcs[patid %in% prescriptions$patid]

                # Calculate LTC prevalence statistics for the selected drug(s)
                calculate_ltc_prevalence_by_drug(
                    ltcs = ltcs,
                    prescriptions = prescriptions,
                    selected_drugs = input$drug_dropdown
                )
            },
            rownames = FALSE
        )
    })
}
