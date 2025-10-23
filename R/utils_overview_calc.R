calc_demog_table_by_pp <- function(pp_df) {
  pp_totals <- pp_df[, .(N = .N), by = pp_group]
  total_patients <- pp_df[, .N]

  # Add row percentages for totals
  pp_totals[, pct_of_total := signif(N/total_patients, digits = 2)]
  pp_totals[, value := sprintf("%d [%.1f%%]", N, pct_of_total * 100)]
  pp_totals[, group := "Total"]
  pp_totals[, category := ""]

  # Sex calculations
  sex_df <- pp_df[, .N, .(sex, pp_group)]
  sex_df[, group := "Sex"]
  sex_df[, pct := signif(N/sum(N), digits = 2), by = pp_group]
  sex_df[, category := sex]

  # Ethnic group calculations
  eth_df <- pp_df[, .N, .(eth_group, pp_group)]
  eth_df[, group := "Ethnic group"]
  eth_df[, pct := signif(N/sum(N), digits = 2), by = pp_group]
  eth_df[, category := eth_group]

  # IMD calculations
  pp_df[, imd_quintile := factor(imd_quintile, c(1,2,3,4,5, NA), ordered = TRUE)]
  imd_df <- pp_df[, .N, .(imd_quintile, pp_group)]
  imd_df[, group := "IMD Quintile"]
  imd_df[, pct := signif(N/sum(N), digits = 2), by = pp_group]
  imd_df[, category := imd_quintile]
  setorder(imd_df, imd_quintile)

  # Age calculations
  age_df <- pp_df[, .N, .(age_group, pp_group)]
  age_df[, group := "Age at outcome"]
  age_df[, pct := signif(N/sum(N), digits = 2), by = pp_group]
  age_df[, category := age_group]

  mltc_df <- pp_df[, .N, .(mltc_group, pp_group)]
  mltc_df[, group := "# LTCs"]
	mltc_df[, pct := signif(N/sum(N), digits = 2), by = pp_group]
  mltc_df[, category := mltc_group]

  tto_df <- pp_df[, .N, .(tto_group, pp_group)]
  tto_df[, group := "Time to outcome since first prescription (years)"]
  tto_df[, pct := signif(N/sum(N), digits = 2), by = pp_group]
  tto_df[, category := tto_group]

  fp_age_df <- pp_df[, .N, .(first_prescription_age_group, pp_group)]
  fp_age_df[, group := "First prescription age group"]
  fp_age_df[, pct := signif(N/sum(N), digits = 2), by = pp_group]
  fp_age_df[, category := first_prescription_age_group]

  # Combine all dataframes including the totals as a group
  to_print <- suppressWarnings(rbindlist(
    list(pp_totals, sex_df, eth_df, imd_df, age_df, mltc_df, tto_df, fp_age_df),
    use.names = TRUE,
    fill = TRUE
  ))

  # For the demographic categories (not totals), create the N (%) string
  to_print[group != "Total", value := sprintf("%d (%.1f%%)", N, pct * 100)]

  # Now reshape wide
  final_table <- dcast(to_print, group + category ~ pp_group, value.var = "value", fill = 0)
  final_table
}
