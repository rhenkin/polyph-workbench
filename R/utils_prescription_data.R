# Prescription Data Transformation Utilities
# This file contains reusable functions for processing prescription data

#' Prepare prescription data with BNF lookup and patient data
#'
#' @param outcome_prescriptions data.table containing prescription data
#' @param bnf_lookup data.table containing BNF classification hierarchy
#' @param pp_groups_data data.table containing PP group assignments
#' @param patient_data data.table containing patient demographic data
#' @param bnf_level character string specifying BNF level (e.g., "BNF_Chapter", "BNF_Section")
#' @return data.table with merged prescription, patient, and PP group data
prepare_prescription_data <- function(outcome_prescriptions, bnf_lookup, pp_groups_data,
                                     patient_data, bnf_level) {
  df <- copy(outcome_prescriptions)
  setkey(df, patid)

  # Map substance to selected BNF level
  df[, substance := bnf_lookup[match(df$substance, bnf_lookup$BNF_Chemical_Substance),
                                get(bnf_level)]]

  # Merge with PP groups and patient data
  df <- merge(df, pp_groups_data, by = "patid")
  df <- merge(df, patient_data, by = "patid")

  return(df)
}

#' Calculate demographic prescription frequency with statistical tests
#'
#' @param prescription_df data.table containing prescription data merged with demographics
#' @param demog_var character string specifying the demographic variable to analyze
#' @param prescription_freq_df data.table containing overall prescription frequencies
#' @param min_prevalence numeric threshold for including prescriptions (default: 0.005)
#' @return formatted data.table with demographic breakdowns and chi-square test results
calculate_demographic_prescription_frequency <- function(prescription_df, demog_var,
                                                         prescription_freq_df,
                                                         min_prevalence = 0.005) {
  return(calculate_demographic_frequency(
    data = prescription_df,
    item_col = "substance",
    demog_var = demog_var,
    frequency_data = prescription_freq_df,
    min_prevalence = min_prevalence
  ))
}

#' Calculate prescription prevalence by LTC(s)
#'
#' @param ltc_data data.table containing LTC data
#' @param prescriptions data.table containing prescription data (age-filtered)
#' @param selected_ltcs character vector of LTC terms to analyze
#' @param min_prevalence numeric threshold for including results (default: 0.005)
#' @return data.table with prevalence ratios, confidence intervals, and duration statistics
calculate_prescription_prevalence_by_ltc <- function(ltc_data, prescriptions,
                                                     selected_ltcs,
                                                     min_prevalence = 0.005) {
  # Patients WITH the selected LTC(s)
  patids <- unique(ltc_data[term %in% selected_ltcs, patid])

  # Use generic prevalence ratio calculation with duration column
  result <- calculate_prevalence_ratio(
    data = prescriptions,
    item_col = "substance",
    selected_patids = patids,
    min_prevalence = min_prevalence,
    duration_col = "duration"
  )

  # Return selected columns in desired order
  return(result[, .(substance, Prevalence, Prevalence_Unselected, Prevalence_Ratio, CI_95,
                    `Median Duration (years)`, `Median Duration unselected (years)`,
                    `IQR (Q1-Q3)`, `IQR unsel. (Q1-Q3)`)])
}

#' Prepare heatmap stratification choices
#'
#' @param prescription_df data.table containing prescription data with demographic variables
#' @return named list of stratification choices suitable for virtualSelectInput
prepare_heatmap_stratification_choices <- function(prescription_df) {
  heatmap_vars <- c("sex", "eth_group", "imd_quintile", "pp_group", "mltc_group")
  var_labels <- c("Sex", "Ethnicity", "IMD quintile", "# PP", "# LTC")

  choices <- setNames(
    lapply(seq_along(heatmap_vars), function(i) {
      var <- heatmap_vars[i]
      unique_vals <- sort(unique(prescription_df[[var]]))

      # Create named list: display_name = encoded_value
      encoded_list <- setNames(
        paste0(var, "#", unique_vals),  # values (encoded)
        unique_vals                      # names (display labels)
      )
      return(encoded_list)
    }),
    var_labels  # group names
  )

  return(choices)
}

#' Prepare prescription OR heatmap data
#'
#' @param prescription_df data.table containing prescription data
#' @param stratification_filter character string in format "variable#value" (empty string for no filter)
#' @param min_prevalence numeric threshold for including prescriptions (default: 0.005)
#' @return data.table with OR calculations ready for plotting
prepare_prescription_or_heatmap_data <- function(prescription_df, stratification_filter = "",
                                                 min_prevalence = 0.005) {
  df <- copy(prescription_df)
  df$yes <- 1

  # Apply stratification filter if specified
  if (stratification_filter != "") {
    parts <- strsplit(stratification_filter, "#")[[1]]
    column_name <- parts[1]
    filter_value <- parts[2]
    df <- df[get(column_name) == filter_value]
  }

  # Create patient-substance matrix
  mat <- dcast(df, patid ~ substance, value.var = "yes",
               fun.aggregate = is.numeric, fill = 0)

  # Calculate prevalence and filter to common prescriptions
  prev <- as.list(mat[, lapply(.SD[, -1], mean)])
  keep_conditions <- names(prev[which(prev > min_prevalence)])

  # Calculate ORs
  ors <- calc_all_ors_vectorized(as.matrix(mat[, ..keep_conditions]))

  # Create pairs data.table
  pairs_dt <- data.table(
    expand.grid(
      drug1 = keep_conditions,
      drug2 = keep_conditions,
      stringsAsFactors = FALSE
    )
  )
  pairs_dt[, or := as.vector(ors$or)]
  pairs_dt[, ci_lower := as.vector(ors$ci_lower)]
  pairs_dt[, ci_upper := as.vector(ors$ci_upper)]

  # Filter to significant ORs > 1
  pairs_dt <- pairs_dt[!is.na(or) & or > 1 & (ci_lower > 1.0 | ci_upper < 1.0)]

  return(list(pairs_dt = pairs_dt, or_matrix = ors$or))
}
