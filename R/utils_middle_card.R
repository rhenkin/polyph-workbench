#' Calculate outcome prescriptions
#' @export
calculate_outcome_prescriptions <- function(presc_df, outcome_df, ltcs, polypharmacy_threshold, 
                                           earliest_treatment_end, min_nltc, queried_terms = NULL) {
  validate(need(nrow(ltcs) > 0, "No valid patients found"))
  
  message("Computing outcome prescriptions ", Sys.time())
  presc_df <- presc_df[patid %in% outcome_df$patid]
  
  # Two-stage join approach for better memory efficiency
  potential_matches <- outcome_df[, .(patid, eventdate, age_days)]
  setkey(potential_matches, patid)
  
  diag_presc <- potential_matches[presc_df,
                                  .(patid, eventdate = x.eventdate, outcome_age = age_days,
                                    substance, start_date, stop_date, duration),
                                  on = .(patid, eventdate > start_date),
                                  nomatch = 0]
  diag_presc <- diag_presc[stop_date >= eventdate - earliest_treatment_end]
  
  diag_presc <- diag_presc[, list(start_date = first(start_date),
                                  stop_date = last(stop_date),
                                  duration = sum(duration)),
                           .(patid, substance, eventdate, outcome_age)]
  
  diag_presc[stop_date >= eventdate, duration := duration - (stop_date - eventdate)]
  
  setkey(diag_presc, patid, start_date)
  
  # Calculate polypharmacy count
  diag_presc[, drugs_count := 1:.N, by = patid]
  diag_presc[, polyph_number := max(drugs_count), by = .(patid, start_date)]
  diag_presc[, N := .N, patid]
  multi_substance <- diag_presc[N >= polypharmacy_threshold]
  
  if (nrow(multi_substance) == 0) return(data.table())
  
  # Filter LTC data
  ltcs <- ltcs[patid %in% multi_substance$patid]
  unique_patients <- unique(multi_substance[, .(patid, eventdate, outcome_age)])
  
  first_ltc <- ltcs[unique_patients,
                    .(patid, eventdate = i.eventdate, outcome_age = i.outcome_age, 
                      age_days = x.age_days, term = x.term),
                    on = .(patid, age_days < outcome_age),
                    nomatch = 0]
  
  if (!is.null(queried_terms)) {
    pats_with_terms <- first_ltc[, all(queried_terms %in% term), patid][V1 == TRUE, patid]
    first_ltc <- first_ltc[patid %in% pats_with_terms]
  }
  
  nltcs <- first_ltc[, list(n_ltc = .N), patid]
  setkey(nltcs, patid)
  
  message("Finished ", Sys.time())
  
  final_df <- multi_substance[nltcs[n_ltc >= min_nltc]]
  final_df <- create_value_groups(final_df, breaks = c(2, 5, 10), right = FALSE, 
                                  value_col = "n_ltc", group_col = "mltc_group")
  return(final_df)
}

#' Calculate acute outcome prescriptions
#' @export
calculate_acute_outcome_prescriptions <- function(outcome_df, acute_df, outcome_prescriptions, 
                                                 earliest_treatment_end) {
  valid_patids <- unique(outcome_prescriptions$patid)
  acute_df <- acute_df[patid %in% valid_patids]
  
  merged <- outcome_df[acute_df, 
                       .(patid, eventdate = x.eventdate, substance, start_date, 
                         outcome_age = age_days), 
                       on = .(patid, eventdate > start_date), 
                       nomatch = 0]
  merged <- merged[start_date >= eventdate - earliest_treatment_end]
  setkey(merged, patid)
  
  return(merged)
}

#' Calculate polypharmacy groups
#' @export
calculate_pp_groups <- function(outcome_prescriptions) {
  grouped_df <- outcome_prescriptions[, list(pp = .N), patid]
  dt <- create_value_groups(grouped_df, breaks = c(2, 5, 10), right = FALSE,
                           value_col = "pp", group_col = "pp_group")
  return(dt)
}

#' Prepare substance data with BNF mapping
#' @export
prepare_substance_data <- function(outcome_prescriptions, bnf_lookup, bnf_level, 
                                  pp_groups_data, patient_data) {
  df <- copy(outcome_prescriptions)
  validate(need(nrow(df) > 0, "No valid patients found"))
  
  setkey(df, patid)
  df[, substance := bnf_lookup[match(df$substance, bnf_lookup$BNF_Chemical_Substance), 
                               get(bnf_level)]]
  df <- merge(df, pp_groups_data, by = "patid")
  df <- merge(df, patient_data, by = "patid")
  
  return(df)
}

#' Calculate substance frequency
#' @export
calculate_substance_frequency <- function(subst_pp_df) {
  unique_patid <- subst_pp_df[, .(substance, patid)] |> unique()
  subst_freq <- unique_patid[, .N, .(substance)]
  total_patids <- uniqueN(subst_pp_df$patid)
  subst_freq[, pct_total := N / total_patids]
  setorder(subst_freq, -N)
  
  return(subst_freq)
}

#' Calculate prescriptions by LTC
#' @export
calculate_presc_by_ltc <- function(prescriptions, ltcs, selected_ltcs, age_filter) {
  patids <- unique(ltcs[term %in% selected_ltcs, patid])
  prescriptions <- prescriptions[outcome_age <= age_filter * 365.25]
  
  # Patients WITH the disease
  presc_freq <- prescriptions[patid %in% patids,
                              list(N_with_disease = uniqueN(patid),
                                   Prevalence = round(100 * (uniqueN(patid) / length(patids)), 
                                                     digits = 2),
                                   `Median Duration (years)` = round(median(duration / 365.2), 
                                                                     digits = 2),
                                   `IQR (Q1-Q3)` = paste0("(",
                                                         round(quantile(duration / 365.2, 0.25, na.rm = TRUE), 2), 
                                                         " - ",
                                                         round(quantile(duration / 365.2, 0.75, na.rm = TRUE), 2), 
                                                         ")")),
                              substance]
  
  # Patients WITHOUT the disease
  unselected_patids <- prescriptions[!patid %in% patids, uniqueN(patid)]
  not_selected_freq <- prescriptions[!patid %in% patids,
                                    list(N_without_disease = uniqueN(patid),
                                         Prevalence_Unselected = round(100 * (uniqueN(patid) / unselected_patids), 
                                                                       digits = 2),
                                         `Median Duration unselected (years)` = round(median(duration / 365.2), 
                                                                                      digits = 2),
                                         `IQR unsel. (Q1-Q3)` = paste0("(",
                                                                       round(quantile(duration / 365.2, 0.25, na.rm = TRUE), 2), 
                                                                       " - ",
                                                                       round(quantile(duration / 365.2, 0.75, na.rm = TRUE), 2), 
                                                                       ")")),
                                    substance]
  
  result <- merge(presc_freq, not_selected_freq)
  result <- result[Prevalence >= 0.005]
  
  # Calculate prevalence ratio and confidence intervals
  result[, `:=`(
    total_with_disease = length(patids),
    total_without_disease = unselected_patids,
    Prevalence_Ratio = round(Prevalence / Prevalence_Unselected, digits = 2)
  )]
  
  result[is.infinite(Prevalence_Ratio), Prevalence_Ratio := 0]
  
  result[, `:=`(
    p1 = N_with_disease / total_with_disease,
    p2 = N_without_disease / total_without_disease
  )]
  
  result[, `:=`(
    log_ratio = log(Prevalence_Ratio),
    se_log_ratio = sqrt((1 / N_with_disease) - (1 / total_with_disease) +
                       (1 / N_without_disease) - (1 / total_without_disease))
  )]
  
  result[, `:=`(
    CI_lower = round(exp(log_ratio - 1.96 * se_log_ratio), digits = 2),
    CI_upper = round(exp(log_ratio + 1.96 * se_log_ratio), digits = 2)
  )]
  
  result[, CI_95 := paste0("(", CI_lower, " - ", CI_upper, ")")]
  result[(CI_lower > 1.0 | CI_upper < 1.0), substance := paste0(substance, "*")]
  
  result[, c("p1", "p2", "log_ratio", "se_log_ratio", "CI_lower", "CI_upper") := NULL]
  
  result[order(-Prevalence_Ratio), 
         .(substance, Prevalence, Prevalence_Unselected, Prevalence_Ratio, CI_95, 
           `Median Duration (years)`, `Median Duration unselected (years)`, 
           `IQR (Q1-Q3)`, `IQR unsel. (Q1-Q3)`)]
}