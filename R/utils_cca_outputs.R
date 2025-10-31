#' Create pyramid plot for burden comparison
#' @param patient_data data.table with treatment, burden column
#' @param burden_col Name of burden column (pp or n_ltc)
#' @param title Plot title
#' @return vegaspec object
create_burden_pyramid <- function(patient_data, burden_col = "pp", title = "Burden") {
  upper_band <- round(quantile(patient_data[[burden_col]], 0.75) * 1.5)
  band_col <- paste0(burden_col, "_band")

  patient_data[, (band_col) := ifelse(
    get(burden_col) >= upper_band,
    paste0(upper_band, "+"),
    as.character(get(burden_col))
  )]

  plot_data <- patient_data[, .N, by = c("treatment", band_col)]
  cases_n <- nrow(patient_data[treatment == 1])
  controls_n <- nrow(patient_data[treatment == 0])

  plot_data[, pct := 0]
  plot_data[treatment == 0, pct := N / controls_n]
  plot_data[treatment == 1, pct := N / cases_n]

  pyramid_plot(
    plot_data,
    band_col,
    rev(c(2:upper_band, paste0(upper_band, "+"))),
    side_width = 175,
    title = title
  ) |> as_vegaspec()
}

#' Create top substances bar plot
#' @param cases_controls data.table with group and substance
#' @param n_top Number of top substances to show
#' @param height Plot height
#' @param full_screen Whether plot is in full screen mode
#' @return vegaspec object
create_top_substances_plot <- function(cases_controls, n_top = 10, height = 250,
                                      full_screen = FALSE) {
  if (full_screen) {
    n_top <- 20
    height <- 500
  }

  cases_controls_n <- cases_controls[, .N, .(group, substance)]
  group_totals <- cases_controls[, .(total = uniqueN(patid)), by = group]
  cases_controls_n[group_totals, pct := round(N / total * 100, 2), on = "group"]

  top_sub <- cases_controls_n[pct >= 1 & group == "case"][order(-pct)][1:n_top, substance]

  cases_controls_n[, `:=`(
    max_pct = max(pct),
    diff_label = sprintf("+%.1f%%", abs(pct[group == "case"] - pct[group == "control"]))
  ), substance]

  cases_controls_n <- cases_controls_n[substance %in% top_sub]
  cases_controls_n[nchar(substance) > 15,
    substance := paste0(strtrim(substance, 15), "...")
  ]

  grouped_bar_plot(cases_controls_n, "substance", height = height, width = 280) |>
    as_vegaspec()
}

#' Create top conditions bar plot for case-control comparison
#' @param freq_data Frequency data from calculate_frequency_stats
#' @param n_top Number of top items to show
#' @param title Plot title
#' @return vegaspec object
create_top_conditions_plot <- function(freq_data, n_top = 10, title = NULL) {
  ratios <- calculate_case_control_ratios(freq_data, "term", min_case_pct = 10)
  top_items <- ratios[order(-case_pct)][1:n_top, term]

  term_filtered <- freq_data[term %in% top_items]
  term_filtered[, `:=`(
    diff = pct[group == "case"] - pct[group == "control"],
    diff_label = sprintf("+%.1f%%", abs(pct[group == "case"] - pct[group == "control"])),
    max_pct = max(pct)
  ), by = term]

  term_filtered[nchar(term) > 15, term := paste0(strtrim(term, 15), "...")]

  grouped_bar_plot(term_filtered, "term", title = title) |> as_vegaspec()
}

#' Create wide prevalence table with ratios
#' @param freq_data Frequency data with group, item, and pct columns
#' @param item_col Name of the item column
#' @param min_pct Minimum prevalence threshold
#' @return data.table in wide format with ratios
create_prevalence_ratio_table <- function(freq_data, item_col, min_pct = 0.5) {
  table_data_wide <- dcast(freq_data,
    as.formula(paste(item_col, "~ group")),
    value.var = "pct",
    fill = 0
  )

  table_data_wide <- table_data_wide[case > min_pct & control > min_pct]
  table_data_wide[, ratio := round(case / control, digits = 2)]

  table_data_wide
}

#' Format and prepare stratification choices for UI
#' @param patient_data data.table with demographic variables
#' @return Named list suitable for virtualSelect choices
create_stratification_choices <- function(patient_data) {
  heatmap_vars <- c("sex", "eth_group", "imd_quintile", "pp_group", "mltc_group")
  var_labels <- c("Sex", "Ethnicity", "IMD quintile", "# PP", "# LTC")

  choices <- setNames(
    lapply(seq_along(heatmap_vars), function(i) {
      var <- heatmap_vars[i]
      unique_vals <- sort(unique(patient_data[[var]]))

      encoded_list <- setNames(
        paste0(var, "#", unique_vals),
        unique_vals
      )
      encoded_list
    }),
    var_labels
  )

  choices
}

#' Calculate prescription prevalence by LTC
#' @param prescriptions data.table of prescriptions
#' @param ltcs data.table of LTCs
#' @param selected_ltcs Character vector of selected LTC terms
#' @return data.table with prevalence by group
calculate_presc_by_ltc_cca <- function(prescriptions, ltcs, selected_ltcs) {
  patids <- unique(ltcs[term %in% selected_ltcs, patid])

  presc_freq <- prescriptions[patid %in% patids,
    list(
      N_with_disease = uniqueN(patid),
      Prevalence = round(100 * (uniqueN(patid) / length(patids)), digits = 2),
      `Median Duration (years)` = round(median(duration / 365.2), digits = 2),
      `IQR (Q1-Q3)` = paste0(
        "(",
        round(quantile(duration / 365.2, 0.25, na.rm = TRUE), 2), " - ",
        round(quantile(duration / 365.2, 0.75, na.rm = TRUE), 2), ")"
      )
    ),
    .(group, substance)
  ]

  result <- dcast(presc_freq, substance ~ group, value.var = "Prevalence", fill = 0)
  result <- result[case >= 1 & control >= 1]
  result[, Prevalence_Ratio := round(case / control, digits = 2)]
  result[is.infinite(Prevalence_Ratio), Prevalence_Ratio := 0]

  result
}