#' Main workflow for creating matched case-control cohort
#'
#' @param outcome_prescriptions data.table of outcome prescriptions
#' @param master_risk_pool_dataset Arrow dataset
#' @param patient_data data.table of patient demographics
#' @param ltc_data data.table of LTC records
#' @param pred_window numeric prediction window in days
#' @param match_ratio numeric control:case ratio
#' @return list(cases, controls, match_summary)
create_matched_cohort_workflow <- function(outcome_prescriptions,
																					 master_risk_pool_dataset,
																					 patient_data,
																					 ltc_data,
																					 pred_window = 30,
																					 match_ratio = 4,
																					 progress,
																					 patient_filters = NULL) {

	# Step 1: Build cases table
	progress$set(message = "Filtering prescription risk table", value = 0.2, detail = "Building cases table")
	cases <- build_cases_table(
		outcome_prescriptions = outcome_prescriptions,
		master_risk_pool = master_risk_pool_dataset,
		ltc_data = ltc_data,
		patient_data = patient_data,
		pred_window = pred_window
	)

	# Step 2: Get eligible control pool
	progress$set(message = "Filtering prescription risk table", value = 0.5, detail= "Finding required controls")
	eligible_pool <- filter_eligible_control_pool(
		cases = cases,
		master_risk_pool = master_risk_pool_dataset,
		ltc_data = ltc_data,
		patient_filters = patient_filters
	)

	# Step 3: Sample controls
	progress$set(message = "Sampling controls...", value = 0.7)
	controls <- sample_controls_by_strata(
		eligible_pool = eligible_pool,
		cases = cases,
		match_ratio = match_ratio
	)

	return(list(
		cases = cases,
		controls = controls
	))
}

#' Build cases table with stratification variables
#'
#' @param outcome_prescriptions data.table
#' @param master_risk_pool Arrow dataset
#' @param ltc_data data.table
#' @param patient_data data.table
#' @param pred_window numeric days
#' @return data.table of cases with strata
build_cases_table <- function(outcome_prescriptions,
															master_risk_pool,
															ltc_data,
															patient_data,
															pred_window) {

	# Get case patient IDs and outcome dates
	outcomes <- outcome_prescriptions[, .(patid, outcome_date = eventdate)] |> unique()
	case_patids <- unique(outcomes$patid)

	# Pull cases from master risk pool
	cases_raw <- extract_cases_from_mrp(
		master_risk_pool = master_risk_pool,
		case_patids = case_patids
	)

	# Join with outcomes and filter by prediction window
	new_cases <- filter_cases_by_pred_window(
		cases_raw = cases_raw,
		outcomes = outcomes,
		pred_window = pred_window
	)

	# Validate LTCs and calculate multimorbidity date
	new_cases <- add_multimorbidity_info(
		cases = new_cases,
		ltc_data = ltc_data
	)

	# Add demographics and calculate age
	new_cases <- add_case_demographics(
		cases = new_cases,
		patient_data = patient_data
	)

	# Create stratification variables
	new_cases <- create_stratification_variables(new_cases)

	return(new_cases)
}

#' Extract case prescriptions from master risk pool
extract_cases_from_mrp <- function(master_risk_pool, case_patids) {
	cases_raw <- master_risk_pool %>%
		dplyr::filter(patid %in% case_patids) %>%
		dplyr::select(patid, prescription_date, substance, stratum_first_presc_bin,
									n_ltcs, concurrent_cps, age_at_rx,
									sex, imd_quintile, first_presc_bin, time_since_first_presc) %>%
		dplyr::collect() %>%
		as.data.table()

	# Date format conversion
	cases_raw[, p_date := as.IDate(prescription_date)]
	cases_raw[, prescription_date := NULL]
	setnames(cases_raw, "p_date", "prescription_date")

	return(cases_raw)
}

#' Filter cases to only prescriptions within prediction window before outcome
filter_cases_by_pred_window <- function(cases_raw, outcomes, pred_window) {
	new_cases <- cases_raw[outcomes, on = "patid", nomatch = 0]

	new_cases <- new_cases[
		prescription_date < outcome_date &
			outcome_date - prescription_date <= pred_window
	]

	# Keep only most recent prescription per patient
	new_cases <- new_cases[, .SD[which(prescription_date == max(prescription_date))], patid]

	return(new_cases)
}

#' Add multimorbidity information (≥2 LTCs before prescription)
add_multimorbidity_info <- function(cases, ltc_data) {
	# Get valid LTCs before prescription date
	valid_ltcs <- cases[ltc_data,
											.(patid, term, eventdate = i.eventdate, start_date = x.prescription_date),
											on = .(patid, prescription_date > eventdate),
											mult = "last",
											nomatch = 0
	]

	setkey(valid_ltcs, patid, eventdate)

	# Filter to patients with ≥2 LTCs
	valid_ltc_patids <- valid_ltcs[, .N, patid][N >= 2, patid]
	valid_ltcs <- valid_ltcs[patid %in% valid_ltc_patids]

	# Find the date of the 2nd LTC (multimorbidity date)
	setorder(valid_ltcs, patid, eventdate)
	valid_ltcs[, ltc_index := 1:.N, patid]
	mm_dates <- valid_ltcs[ltc_index == 2, .(patid, mm_date = eventdate)]

	# Merge back to cases
	cases <- merge(cases, mm_dates, by = "patid")

	return(cases)
}

#' Add patient demographics and calculate ages
add_case_demographics <- function(cases, patient_data) {
	cases <- merge(cases,
								 patient_data[, .(patid, dob, eth_group)],
								 by = "patid")

	cases[, ':='(
		index_age = as.numeric(prescription_date - dob) / 365.25,
		time_at_risk = round(as.numeric(outcome_date - mm_date) / 365.25, digits = 0)
	)]

	return(cases)
}

#' Create stratification variables for matching
create_stratification_variables <- function(cases) {
	cases[, mm_duration_bin := cut(
		time_at_risk,
		breaks = seq(0, 121, 2),
		include.lowest = TRUE,
		right = FALSE
	)]

	cases[, ':='(
		age_bin = cut(index_age, breaks = seq(0, 120, 5),
									include.lowest = TRUE, right = FALSE),
		year = year(prescription_date)
	)]

	# Primary stratification variable
	cases[, strata := stratum_first_presc_bin]

	return(cases)
}

#' Filter master risk pool to eligible controls
#'
#' @param cases data.table of cases with strata
#' @param master_risk_pool Arrow dataset
#' @return data.table of eligible control pool
filter_eligible_control_pool <- function(cases, master_risk_pool, ltc_data, patient_filters = NULL) {

	case_patids <- unique(cases$patid)
	strata_needs <- unique(cases$strata)

	# Use Arrow to filter before collecting
	eligible_pool <- master_risk_pool |>
		dplyr::filter(!patid %in% case_patids) |>
		dplyr::filter(stratum_first_presc_bin %in% strata_needs)

	# Apply patient filters if provided
	if (!is.null(patient_filters)) {
		input_list <- patient_filters$input_list

		if (!is.null(patient_filters$selected_ltcs)) {

			patids_with_ltc <- ltc_data[term %in% patient_filters$selected_ltcs, patid]
			eligible_pool <- eligible_pool |>
				dplyr::filter(patid %in% patids_with_ltc)
			message(sprintf("Filtering controls by LTC: %s",
											paste(patient_filters$selected_ltcs, collapse = ", ")))
		}

		# Filter by minimum number of LTCs
		if (!is.null(input_list$min_nltc) && input_list$min_nltc > 0) {
			eligible_pool <- eligible_pool |>
				dplyr::filter(n_ltcs >= input_list$min_nltc)
			message(sprintf("Filtering controls by min LTCs: %d", input_list$min_nltc))
		}

		# Filter by sex (if specified)
		if (!is.null(input_list$sex) && length(input_list$sex) > 0) {
			eligible_pool <- eligible_pool |>
				dplyr::filter(sex %in% input_list$sex)
			message(sprintf("Filtering controls by sex: %s", paste(input_list$sex, collapse = ", ")))
		}

		# Filter by ethnicity (if specified)
		if (!is.null(input_list$eth_group) && length(input_list$eth_group) > 0) {
			# Need to check if eth_group exists in master_risk_pool
			if ("eth_group" %in% names(master_risk_pool)) {
				eligible_pool <- eligible_pool |>
					dplyr::filter(eth_group %in% input_list$eth_group)
				message(sprintf("Filtering controls by ethnicity: %s",
												paste(input_list$eth_group, collapse = ", ")))
			} else {
				warning("eth_group not available in master_risk_pool - skipping ethnicity filter")
			}
		}

		# Filter by IMD quintile (if specified)
		if (!is.null(patient_filters$imd_quintile) && length(patient_filters$imd_quintile) > 0) {
			# Convert to numeric for comparison
			imd_values <- as.numeric(patient_filters$imd_quintile)
			eligible_pool <- eligible_pool |>
				dplyr::filter(imd_quintile %in% imd_values)
			message(sprintf("Filtering controls by IMD quintile: %s",
											paste(patient_filters$imd_quintile, collapse = ", ")))
		}
	}

	# Now collect and convert to data.table
	eligible_pool <- eligible_pool |>
		dplyr::select(patid, prescription_date, substance, sex, age_at_rx,
									n_ltcs, imd_quintile, age_bin, stratum_first_presc_bin,
									year, first_presc_bin, time_since_first_presc) |>
		dplyr::collect() |>
		as.data.table()

	# Convert date format
	eligible_pool[, prescription_date := as.IDate(prescription_date)]

	# Exclude case patients
	# eligible_pool <- eligible_pool[!patid %in% case_patids]


	return(eligible_pool)
}

#' Sample controls within strata
#'
#' @param eligible_pool data.table of eligible controls
#' @param cases data.table of cases
#' @param match_ratio numeric desired control:case ratio
#' @return data.table of sampled controls
sample_controls_by_strata <- function(eligible_pool, cases, match_ratio) {

	# Get unique patient per stratum (most recent prescription)
	setorder(eligible_pool, patid, prescription_date)
	eligible_unique <- eligible_pool[
		,
		last(.SD),
		by = .(strata = stratum_first_presc_bin, patid)
	]

	message(sprintf("Unique patient-strata combinations: %d", nrow(eligible_unique)))

	# Calculate controls needed per stratum
	strata_needs_dt <- cases[
		,
		.(n_cases = uniqueN(patid), n_controls_needed = uniqueN(patid) * match_ratio),
		by = strata
	]

	eligible_unique[strata_needs_dt, n_needed := i.n_controls_needed, on = "strata"]

	# Sample within each stratum
	all_controls <- eligible_unique[
		!is.na(n_needed),
		{
			n_sample <- min(.N, n_needed[1])
			if(n_sample > 0) .SD[sample(.N, n_sample, replace = FALSE)] else .SD[0]
		},
		by = strata
	]

	all_controls[, n_needed := NULL]
	all_controls[, control_index_date := prescription_date]

	# Remove duplicate patients (keep random occurrence)
	setorder(all_controls, patid, prescription_date)
	all_controls <- all_controls[, .SD[sample(.N, 1)], by = patid]

	message(sprintf("Controls sampled: %d (unique patients: %d)",
									nrow(all_controls), uniqueN(all_controls$patid)))

	return(all_controls)
}

