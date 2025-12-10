#' Load and prepare case-control study dataset
#' @param dataset_path Path to the saved dataset file
#' @return List containing patient_data, prescriptions, ltcs, and cases_controls
load_cca_dataset <- function(dataset_path) {
  saved_dataset_data <- readRDS(dataset_path)

  matched_patids <- saved_dataset_data$matched_patids
  matched_patids[, group := ifelse(treatment == 1, "case", "control")]

  patient_data <- saved_dataset_data$all_patient_data
  setkey(patient_data, patid)

  prescriptions <- unique(saved_dataset_data$all_prescriptions)
  setkey(prescriptions, patid, start_date, substance)

  ltcs <- unique(saved_dataset_data$all_ltc)
  setkey(ltcs, patid, eventdate, term)

  list(
    matched_patids = matched_patids,
    patient_data = patient_data,
    prescriptions = prescriptions,
    ltcs = ltcs
  )
}

#' Filter prescriptions and LTCs by index date
#' @param prescriptions data.table of prescriptions
#' @param ltcs data.table of long-term conditions
#' @param matched_patids data.table with patid, index_date, group
#' @param lookback_days Number of days before index_date to include prescriptions
#' @return List with filtered prescriptions and ltcs
filter_by_index_date <- function(prescriptions, ltcs, matched_patids, lookback_days = 84) {
  umatched_patids <- unique(matched_patids[, .(patid, index_date, group, strata)])

  ltcs_filtered <- ltcs[umatched_patids,
    .(patid, eventdate, age_days, term, group, strata),
    on = .(patid, strata, eventdate < index_date),
    nomatch = 0
  ]

  presc_filtered <- prescriptions[umatched_patids,
    .(patid, substance = x.substance, index_date = index_date,
      start_date = x.start_date, stop_date = x.stop_date, duration, group, strata),
    on = .(patid, strata, start_date <= index_date),
    nomatch = 0
  ]

  presc_filtered <- presc_filtered[stop_date >= index_date - lookback_days]

  list(prescriptions = presc_filtered, ltcs = ltcs_filtered)
}

#' Add polypharmacy and MLTC groups to patient data
#' @param patient_data data.table of patient demographics
#' @param prescriptions data.table of prescriptions
#' @param ltcs data.table of long-term conditions
#' @return Updated patient_data with pp_group and mltc_group columns
add_burden_groups <- function(patient_data, prescriptions, ltcs) {
  prescriptions_n <- prescriptions[, list(pp = .N), .(patid, strata)]
  prescriptions_n <- create_value_groups(
    prescriptions_n,
    breaks = c(2, 5, 10),
    right = FALSE,
    value_col = "pp",
    group_col = "pp_group"
  )

  ltcs_n <- ltcs[, list(n_ltc = .N), .(patid, strata)]
  ltcs_n <- create_value_groups(
    ltcs_n,
    breaks = c(2, 5, 10),
    right = FALSE,
    value_col = "n_ltc",
    group_col = "mltc_group"
  )

  patient_data <- merge(patient_data, prescriptions_n, by = c("patid", "strata"))
  patient_data <- merge(patient_data, ltcs_n, by = c("patid", "strata"))

  patient_data
}

#' Calculate frequency statistics for conditions or substances
#' @param data data.table with group and term/substance columns
#' @param item_col Name of the column containing items (term or substance)
#' @return data.table with N, pct by group and item
calculate_frequency_stats <- function(data, item_col) {
  # freq <- data[, .(N = uniqueN(patid)), by = c("group", item_col)]
  # group_totals <- data[, .(total = uniqueN(patid)), by = group]
  freq <- data[, .(N = uniqueN(.SD,by = c("patid","strata"))), by = c("group", item_col)]
  group_totals <- data[,.(total = uniqueN(.SD,by = c("patid","strata"))),group]
  freq[group_totals, pct := round(N / total * 100, 2), on = "group"]
  setorder(freq, -pct)
  freq
}

#' Calculate prevalence ratios between cases and controls
#' @param freq_data data.table from calculate_frequency_stats
#' @param item_col Name of the item column
#' @param min_case_pct Minimum case prevalence to include
#' @return data.table with ratios and differences
calculate_case_control_ratios <- function(freq_data, item_col, min_case_pct = 1) {
  ratios <- freq_data[, .(
    case_pct = pct[group == "case"],
    control_pct = pct[group == "control"]
  ), by = item_col]

  ratios[, `:=`(
    ratio = case_pct / control_pct,
    diff = case_pct - control_pct
  )]

  ratios[case_pct >= min_case_pct]
}

#' Apply patient stratification filter
#'
#' Filters data based on stratification variable in format "column#value"
#'
#' @param data data.table to filter
#' @param strat_var Stratification variable string in format "column#value", or empty string
#' @param patient_data data.table with patient information
#' @return Filtered data.table
apply_patient_stratification <- function(data, strat_var, patient_data) {
	if (strat_var == "" || is.null(strat_var)) return(data)

	parts <- strsplit(strat_var, "#")[[1]]
	column_name <- parts[1]
	filter_value <- parts[2]
	selected_patids <- patient_data[get(column_name) == filter_value, patid]

	data[patid %in% selected_patids]
}