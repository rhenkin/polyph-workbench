#' Generic function to process demographic group data
#'
#' @param dt Data table containing demographic and pp data
#' @param group_var Character string of the column name to group by
#' @param group_label Character string for the group label in output
#' @param total_n Total number of patients for percentage calculation
#' @param pp_var Character string of the polypharmacy variable to group by (default "pp_group")
#' @return Data table with columns: group, category, pp_var, N, pct (or pct_total for cohort table)
process_demog_group <- function(dt, group_var, group_label, total_n = NULL, pp_var = NULL) {
  if (is.null(pp_var)) {
    # For cohort table - calculate overall percentages
    result <- dt[, .N, by = group_var]
    result[, group := group_label]
    result[, pct_total := signif(N/total_n, digits = 2)]
    setnames(result, old = group_var, new = "category")
  } else {
    # For PP distribution table - calculate percentages within PP groups
    result <- dt[, .N, by = c(group_var, pp_var)]
    result[, group := group_label]
    result[, pct := signif(N/sum(N), digits = 2), by = pp_var]
    setnames(result, old = group_var, new = "category")
  }
  return(result)
}

calc_demog_table_by_pp <- function(pp_df) {
  pp_totals <- pp_df[, .(N = .N), by = pp_group]
  total_patients <- pp_df[, .N]

  # Add row percentages for totals
  pp_totals[, pct_of_total := signif(N/total_patients, digits = 2)]
  pp_totals[, value := sprintf("%d [%.1f%%]", N, pct_of_total * 100)]
  pp_totals[, group := "Total"]
  pp_totals[, category := ""]

  # Define demographic groups to process
  demog_groups <- list(
    list(var = "sex", label = "Sex"),
    list(var = "eth_group", label = "Ethnic group"),
    list(var = "imd_quintile", label = "IMD Quintile"),
    list(var = "age_group", label = "Age at outcome"),
    list(var = "mltc_group", label = "# LTCs"),
    list(var = "tto_group", label = "Time to outcome since first prescription (years)"),
    list(var = "first_prescription_age_group", label = "First prescription age group")
  )

  # Process IMD ordering before processing
  pp_df[, imd_quintile := factor(imd_quintile, c(1,2,3,4,5, NA), ordered = TRUE)]

  # Process all demographic groups using the generic function
  demog_dfs <- lapply(demog_groups, function(g) {
    df <- process_demog_group(pp_df, g$var, g$label, pp_var = "pp_group")
    # Apply ordering for IMD if needed
    if (g$var == "imd_quintile") {
      setorder(df, category)
    }
    return(df)
  })

  # Combine all dataframes including the totals as a group
  to_print <- suppressWarnings(rbindlist(
    c(list(pp_totals), demog_dfs),
    use.names = TRUE,
    fill = TRUE
  ))

  # For the demographic categories (not totals), create the N (%) string
  to_print[group != "Total", value := sprintf("%d (%.1f%%)", N, pct * 100)]

  # Now reshape wide
  final_table <- dcast(to_print, group + category ~ pp_group, value.var = "value", fill = 0)
  final_table
}

#' Prepare data for cohort demographic distribution table
#'
#' @param pp_df Data table with patient polypharmacy and demographic data
#' @return Data table ready for gt table formatting with columns: group, category, N, pct_total, pp_labels, pp_values
prepare_cohort_demog_data <- function(pp_df) {
  total_patids <- length(unique(pp_df$patid))

  # Define demographic groups to process
  demog_groups <- list(
    list(var = "sex", label = "Sex"),
    list(var = "eth_group", label = "Ethnic group"),
    list(var = "imd_quintile", label = "IMD Quintile"),
    list(var = "age_group", label = "Age at outcome"),
    list(var = "mltc_group", label = "# LTCs")
  )

  # Helper function to calculate PP frequencies
  pp_frequencies <- function(dt, group_var) {
    all_pp <- seq_len(max(dt$pp))
    freq_table <- dt[, as.list(table(factor(pp, levels=all_pp))), by=group_var]
    result <- freq_table[, list(
      pp_label = paste0(all_pp, collapse=","),
      pp_values = paste0(unlist(.SD), collapse=",")
    ), by = group_var]
    setnames(result, old = group_var, new = "category")
    return(result)
  }

  # Process all demographic groups
  demog_dfs <- lapply(demog_groups, function(g) {
    df <- process_demog_group(pp_df, g$var, g$label, total_n = total_patids)
    df <- merge(df, pp_frequencies(pp_df, g$var), by = "category")

    # Apply ordering for specific groups
    if (g$var == "imd_quintile") {
      setorder(df, category)
    } else if (g$var == "age_group") {
      setorder(df, category)
    }
    return(df)
  })

  # Combine all demographic groups
  to_print <- rbindlist(demog_dfs, use.names = FALSE)
  colnames(to_print) <- c("Category", "N", "group", "%", "pp_labels", "pp_values")

  return(to_print)
}
