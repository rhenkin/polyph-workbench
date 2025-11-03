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
  final_df[, age_years := round(outcome_age/365.25, digits = 0)]
  final_df <- create_value_groups(final_df, breaks = c(0,45,65,85), right = FALSE,
  																value_col = "age_years", group_col = "age_group",
  																label_fmt = "%g-%g",
  																first_prefix = "<=",
  																last_suffix = "+"
  																)

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