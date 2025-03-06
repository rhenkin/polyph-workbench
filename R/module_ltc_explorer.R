module_ltc_explorer_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Long-term conditions"),
    accordion(
      open = FALSE,
      accordion_panel(
        title = "Frequency table",
        navset_tab(
          nav_panel("Cohort frequency", dataTableOutput(ns("ltc_freq_table"))),
          nav_panel("Frequency across demographics",
                    selectizeInput(ns("select_ltcdemog_freq_var"), label = "Select variable:", choices = c("sex", "eth_group",  "imd_quintile", "pp_group")),
                    dataTableOutput(ns("demog_ltc_freq_table")))
        )
      )
      # accordion_panel(
      #   title = "Network",
      #   value = "ltcnet",
      #   visNetworkOutput(ns("ltc_network"), height = "600px")
      # )
    )
  )
}

module_ltc_explorer_server <- function(id, outcome_prescriptions, ltc_data, patient_data) {
  moduleServer(id, function(input, output, session) {
    ns <- NS(id)
    ltc_chapters <- fread("chapters.tsv")
    valid_ltcs <- reactive({
      outcome_df <- outcome_prescriptions()
      ltcs <- ltc_data()
      ltcs[patid %in% outcome_df$patid]
    })

    ltc_freq_df <- reactive({
      ltcs <- valid_ltcs()
      term_freq <- ltcs[,.N,.(term)]
      total_patids <- length(unique(ltcs$patid))
      term_freq[, pct_total := signif(N/total_patids, digits = 2)]
      setorder(term_freq,-N)
    })

    output$ltc_freq_table <- renderDataTable({
      ltc_freq_df()
    })

    output$demog_ltc_freq_table <- renderDataTable({
      req(input$select_ltcdemog_freq_var)
      df <- valid_ltcs()
      demog_var <- input$select_ltcdemog_freq_var
      #df <- patient_data()[df]
      df <- merge(patient_data(), df, by = "patid")
      cat_totals <- df[, .(Total = uniqueN(patid)), by = demog_var]
      # df_stats <- df[, .(N = uniqueN(patid)), by = c(demog_var, "term")][
      #   cat_totals,
      #   on = demog_var
      # ][
      #   , pct := sprintf("%.2f", round(N/Total, digit=2))
      # ]

      df_stats <- df[, .(
        N_category = uniqueN(patid)  # Count of patients in each category
      ), by = c(demog_var, "term")]
      df_stats <- merge(df_stats,
                        ltc_freq_df()[,.(term, N)],
                        by = "term")
      df_stats[, pct := signif(N_category/N, digit = 2)]
      # Now calculate percentage using N from subst_frequency
      df_counts <- copy(df_stats)
      cohort_props <- cat_totals[, .(prop = Total/sum(Total))]

      chisq_tests <- df_counts[, {
        obs <- c(N_category)  # Make sure we have a vector
        if(length(obs) == length(cat_totals[[demog_var]])) {  # Check we have both male and female
          #exp <- sum(N_category) * cohort_props$prop
          if(sum(N_category) >= 10) {
            test <- suppressWarnings(chisq.test(obs, p = cohort_props$prop))
            pval <- test$p.value
          } else {
            pval <- NA_real_
          }
        } else {
          pval <- NA_real_
        }
        .(pvalue = pval)
      }, by=term]

      chisq_tests[, padj := p.adjust(pvalue, method="bonferroni")]


      # Format as percentage
      result <- dcast(df_stats,
                      term ~ get(demog_var),
                      value.var = "pct",
                      fill = 0)
      result[chisq_tests,`:=`(pvalue = pvalue, padj = i.padj), on="term"]
      # Rename columns to include totals
      old_names <- as.character(unique(df_stats[[demog_var]]))
      new_names <- paste0(old_names, " (", signif(cat_totals$Total/sum(cat_totals$Total), digit = 2), ")")
      #new_names <- paste0(old_names, " (n=", cat_totals$Total, ")")

      setnames(result,
               old = old_names,
               new = new_names)
      result[!is.na(padj), pvalue := sprintf("%.3f", pvalue)]
      result[!is.na(padj), padj := sprintf("%.3f", padj)]
      result <- merge(result, ltc_freq_df()[,.(term, N, pct_total)])
      result[order(-N)]
    }, rownames = FALSE)

    # output$ltc_network <- renderVisNetwork({
    #   outcome_df <- outcome_prescriptions()
    #   ltcs <- ltc_data()
    #   ltcs <- ltcs[patid %in% outcome_df$patid]
    #   # browser()
    #   # df <- merge(ltcs, ltc_chapters, by.x = "term", by.y = "ltc")
    #   # df[,term := body_system]
    #   network_dt <- calculate_disease_network(ltcs, filter_method = "prevalence")
    #   g <- create_igraph_network(network_dt[above_threshold==TRUE])
    #   network_vis <- prepare_visnetwork(g)
    # })
  })
}