# LTC Data Transformation Utilities
# This file contains reusable functions for processing long-term condition (LTC) data
# Note: Generic frequency calculation functions are available in utils_calculations.R

#' Calculate LTC frequency
#'
#' @param ltcs data.table containing LTC data with columns: patid, term
#' @return data.table with columns: term, N, pct_total (sorted by frequency descending)
calculate_ltc_frequency <- function(ltcs) {
    result <- calculate_item_frequency(ltcs, "term")
    result[, pct_total := signif(pct_total, digits = 2)]
    return(result)
}

#' Calculate demographic LTC frequency with statistical tests
#'
#' @param ltcs data.table containing LTC data
#' @param patient_data data.table containing patient demographic data
#' @param pp_groups_data data.table containing PP group assignments
#' @param outcome_df data.table containing outcome prescriptions with mltc_group
#' @param demog_var character string specifying the demographic variable to analyze
#' @param ltc_freq_df data.table containing overall LTC frequencies (from calculate_ltc_frequency)
#' @param min_prevalence numeric threshold for including LTCs (default: 0.005)
#' @return formatted data.table with demographic breakdowns and chi-square test results
calculate_demographic_ltc_frequency <- function(ltcs, patient_data, pp_groups_data,
                                                outcome_df, demog_var, ltc_freq_df,
                                                min_prevalence = 0.005) {
    # Merge all relevant data
    df <- merge(patient_data, ltcs, by = "patid")
    df <- merge(df, pp_groups_data, by = "patid")
    df <- merge(df, unique(outcome_df[, .(patid, mltc_group)]), by = "patid")

    # Use generic demographic frequency calculation
    return(calculate_demographic_frequency(
        data = df,
        item_col = "term",
        demog_var = demog_var,
        frequency_data = ltc_freq_df,
        min_prevalence = min_prevalence
    ))
}

#' Calculate LTC prevalence by prescription drug(s)
#'
#' @param ltcs data.table containing LTC data filtered to relevant patients
#' @param prescriptions data.table containing prescription data (age-filtered)
#' @param selected_drugs character vector of drug substance names to analyze
#' @param min_prevalence numeric threshold for including results (default: 0.01)
#' @return data.table with prevalence ratios and confidence intervals
calculate_ltc_prevalence_by_drug <- function(ltcs, prescriptions, selected_drugs,
                                             min_prevalence = 0.01) {
    # Patients taking the selected drug(s)
    patids <- unique(prescriptions[substance %in% selected_drugs, patid])

    # Use generic prevalence ratio calculation
    result <- calculate_prevalence_ratio(
        data = ltcs,
        item_col = "term",
        selected_patids = patids,
        min_prevalence = min_prevalence,
        duration_col = NULL
    )

    # Return selected columns in desired order
    return(result[, .(term, Prevalence, Prevalence_Unselected, Prevalence_Ratio, CI_95)])
}
