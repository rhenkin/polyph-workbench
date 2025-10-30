#' Prepare study data structure from matching results (without saving)
#' @param study_name character string
#' @param cases data.table with case information
#' @param controls data.table with control information
#' @param gold_patient data.table with patient demographic data
#' @param gold_cp data.table with prescription data
#' @param gold_ltc data.table with long-term condition data
#' @param outcome_prescriptions data.table with outcome prescription data
#' @return list with all study data components
prepare_study_data <- function(study_name, cases, controls, gold_patient, gold_cp, gold_ltc,
																outcome_prescriptions) {
	# Validate inputs
	if (!all(c("patid", "outcome_date") %in% colnames(cases))) {
		stop("cases must have columns: patid, index_date, outcome_date")
	}
	if (!all(c("patid", "control_index_date") %in% colnames(controls))) {
		stop("controls must have columns: patid, control_index_date")
	}

	# Get patient IDs
	case_patids <- unique(cases$patid)
	control_patids <- unique(controls$patid)

	n_cases <- length(case_patids)
	n_controls <- uniqueN(controls$patid)

	message(sprintf("Preparing study data: %d cases, %d unique controls (%d total control observations)",
									n_cases, n_controls, nrow(controls)))

	# ========================================================================
	# Create matched prescriptions
	# ========================================================================

	# Cases: use the filtered outcome_prescriptions
	cases_presc <- unique(cases[, .(patid, index_date = outcome_date, eventdate = outcome_date, outcome_age = index_age * 365.25)])
	cases_presc <- cases_presc[gold_cp, on = "patid", nomatch = 0][start_date <= index_date & stop_date >= index_date - 84]
	cases_presc[, treatment := 1]

	# Controls: get prescriptions active at their index date
	control_index_lookup <- unique(controls[, .(patid, index_date = control_index_date)])

	# Get all CP prescriptions for controls that were active at index date
	controls_presc_raw <- control_index_lookup[gold_cp, on = "patid", nomatch = 0][
		start_date <= index_date & stop_date >= index_date - 84  # 84 days lookback
	]

	# Add patient demographics to calculate outcome_age
	controls_dob <- gold_patient[patid %in% control_patids, .(patid, dob)]
	controls_presc <- merge(controls_presc_raw, controls_dob, by = "patid")

	# Create structure matching outcome_prescriptions
	controls_presc[, `:=`(
		eventdate = index_date,
		outcome_age = as.numeric(as.IDate(index_date) - as.IDate(dob)),
		treatment = 0
	)]
	controls_presc[, `:=`(dob = NULL, index_date = NULL)]

	# Combine prescription data
	all_prescriptions <- rbindlist(list(
		cases_presc[, .(patid, eventdate, outcome_age, substance, start_date, stop_date, duration, treatment)],
		controls_presc[, .(patid, eventdate, outcome_age, substance, start_date, stop_date, duration, treatment)]
	), fill = TRUE)
	all_prescriptions[, study_name := study_name]
    setkey(all_prescriptions, patid)

	# ========================================================================
	# Create matched patient data
	# ========================================================================

	cases_patient_data <- gold_patient[patid %in% case_patids]
	cases_patient_data[, treatment := 1]

	controls_patient_data <- gold_patient[patid %in% control_patids]
	controls_patient_data[, treatment := 0]

	all_patient_data <- rbindlist(list(cases_patient_data, controls_patient_data))
	all_patient_data[, study_name := study_name]
	setkey(all_patient_data, patid)

	# ========================================================================
	# Create matched LTC data
	# ========================================================================

	# Cases: get all LTC data before their outcome
	cases_index <- cases[, .(patid, index_date = outcome_date)]
	cases_ltc <- cases_index[gold_ltc, on = "patid", nomatch = 0][
		eventdate < index_date
	]
	cases_ltc[, treatment := 1]
	cases_ltc[, index_date := NULL]

	# Controls: get LTC data before their control index date
	controls_ltc <- control_index_lookup[gold_ltc, on = "patid", nomatch = 0][
		eventdate < index_date
	]
	controls_ltc[, treatment := 0]
	controls_ltc[, index_date := NULL]

	# Combine LTC data
	all_ltc <- rbindlist(list(
		cases_ltc[, .(patid, eventdate, age_days, term, treatment)],
		controls_ltc[, .(patid, eventdate, age_days, term, treatment)]
	), fill = TRUE)
	all_ltc[, study_name := study_name]
    setkey(all_ltc, patid, eventdate)

	# ========================================================================
	# Create matched_patids for backwards compatibility
	# ========================================================================

	# Cases
	matched_cases <- cases[, .(
		patid,
		index_date = outcome_date,
		substance,
		treatment = 1,
		study_name = study_name
	)]

	# Controls
	matched_controls <- controls[, .(
		patid,
		index_date = control_index_date,
		substance,
		treatment = 0,
		study_name = study_name
	)]

	matched_patids <- rbindlist(list(matched_cases, matched_controls), use.names = TRUE)
    setkey(matched_patids, patid)

	# ========================================================================
	# Create study metadata
	# ========================================================================

	metadata <- list(
		study_name = study_name,
		n_cases = n_cases,
		n_controls_unique = n_controls,
		n_controls_total = nrow(controls),
		matching_ratio_unique = n_controls / n_cases,
		matching_ratio_total = nrow(controls) / n_cases,
		case_date_range = range(cases$outcome_date),
		control_date_range = range(controls$control_index_date),
		created_date = Sys.Date(),
		matching_method = "stratified_sampling",
		strata_variable = "sex_age_bin_n_ltc_year"
	)

	# ========================================================================
	# Prepare final structure
	# ========================================================================
	colnames(all_patient_data)[3] <- "sex"
	study_data <- list(
		all_prescriptions = all_prescriptions,
		all_patient_data = all_patient_data,
		all_ltc = all_ltc,
		matched_patids = matched_patids,
		cases = cases,
		controls = controls,
		metadata = metadata
	)

	message(sprintf("Study data prepared - Cases: %d, Unique controls: %d (Total: %d)",
									n_cases, n_controls, nrow(controls)))

	return(study_data)
}

save_matched_datasets <- function(study_name, cases, controls, gold_patient, gold_cp, gold_ltc,
																	outcome_prescriptions, study_dir = "studies") {
	# Create directory if it doesn't exist
	if (!dir.exists(study_dir)) {
		dir.create(study_dir, recursive = TRUE)
	}

	# Prepare the study data using the helper function
	study_data <- prepare_study_data(
		study_name = study_name,
		cases = cases,
		controls = controls,
		gold_patient = gold_patient,
		gold_cp = gold_cp,
		gold_ltc = gold_ltc,
		outcome_prescriptions = outcome_prescriptions
	)

	# Save to disk
	saveRDS(study_data, file.path(study_dir, paste0(study_name, "_study_data.rds")))

	message("Study '", study_name, "' saved successfully to ", study_dir)
	message(sprintf("Matching ratio - Unique: %.2f:1, Total: %.2f:1",
									study_data$metadata$matching_ratio_unique,
									study_data$metadata$matching_ratio_total))

	invisible(study_data)
}

compute_eligible_controls_core_new <- function(cases_df, gold_patient, gold_ltc, random_seed = 4) {
	# Set random seed for reproducibility
	set.seed(random_seed)

	# Get case patient IDs and age range
	case_patids <- unique(cases_df$patid)
	cases_with_age <- merge(cases_df, gold_patient[, .(patid, dob)], by = "patid")
	cases_with_age[, outcome_age := as.numeric(as.IDate(eventdate) - as.IDate(dob))]

	min_case_age <- min(cases_with_age$outcome_age)
	max_case_age <- max(cases_with_age$outcome_age)

	# Get all potential controls - vectorized filtering
	all_controls <- gold_patient[!patid %in% case_patids]

	# Calculate age range dates for all controls at once
	all_controls[, `:=`(
		min_age_date = as.IDate(dob) + min_case_age,
		max_age_date = as.IDate(dob) + max_case_age
	)]

	# Vectorized eligibility check
	eligible_controls <- all_controls[
		# Patient was alive during potential age range
		(is.na(dod) | as.IDate(dod) > min_age_date) &
			(is.na(tod) | as.IDate(tod) > min_age_date) &
			# Age range overlaps with study period
			min_age_date <= as.IDate("2020-12-31") &
			max_age_date >= as.IDate("2005-01-01")
	]

	# Assign random ages and calculate index dates - vectorized
	set.seed(random_seed)
	eligible_controls[, assigned_age := sample(min_case_age:max_case_age, .N, replace = TRUE)]
	eligible_controls[, index_date := as.IDate(dob) + assigned_age]

	# Final validation - vectorized
	valid_controls <- eligible_controls[
		(is.na(dod) | index_date < as.IDate(dod)) &
			(is.na(tod) | index_date < as.IDate(tod)) &
			index_date >= as.IDate("2005-01-01") &
			index_date <= as.IDate("2020-12-31"),
		.(patid, index_date)
	]

	# Rest of the function remains the same...
	setkey(valid_controls, patid)
	setkey(gold_ltc, patid)

	controls_ltc <- gold_ltc[valid_controls, on = "patid", nomatch = 0][
		eventdate < index_date
	]

	ltc_counts <- controls_ltc[, .(n_ltc_terms = uniqueN(term)), by = .(patid, index_date)]
	controls_with_valid_ltc <- ltc_counts[n_ltc_terms >= 2, .(patid, index_date)]

	controls_with_terms <- controls_ltc[controls_with_valid_ltc, on = c("patid", "index_date"), nomatch = 0][
		, .(patid, index_date, term)
	]

	return(controls_with_terms)
}

compute_eligible_controls_core <- function(cases_df, gold_patient, gold_ltc, gold_cp, random_seed = 4) {
	# Set random seed for reproducibility
	set.seed(random_seed)

	# Get all case dates for random assignment to controls
	case_dates <- unique(cases_df$eventdate)
	date_count <- length(case_dates)

	if (date_count == 0) {
		stop("No case dates found in cases_df")
	}

	# Get case patient IDs
	case_patids <- unique(cases_df$patid)

	# Get all potential controls (patients not in cases)
	all_controls <- gold_patient[!patid %in% case_patids, .(patid)]

	# Assign random case dates to controls
	potential_controls <- all_controls[, .(
		patid,
		index_date = sample(case_dates, .N, replace = TRUE)
	)]

	# Validate controls based on death dates and study period
	# Join with patient data to get death dates
	valid_controls <- merge(potential_controls, gold_patient, by = "patid")[
		# Control must be alive at index date (tod = transfer out date, dod = date of death)
		(is.na(tod) | index_date < as.IDate(tod)) &
			(is.na(dod) | index_date < as.IDate(dod)) &
			# Index date must be within study period
			index_date <= as.IDate("2020-12-31"),
		.(patid, index_date)
	]

	# Find controls with at least 2 distinct LTC terms before index date
	# First get all LTC records for valid controls before their index date
	# Use a more careful join to avoid cartesian products
	setkey(valid_controls, patid)
	setkey(gold_ltc, patid)

	controls_ltc <- gold_ltc[valid_controls, on = "patid", nomatch = 0][
		eventdate < index_date
	]

	# Count distinct LTC terms per patient
	ltc_counts <- controls_ltc[, .(n_ltc_terms = uniqueN(term)), by = .(patid, index_date)]

	# Keep only controls with 2+ distinct LTC terms
	controls_with_valid_ltc <- ltc_counts[n_ltc_terms >= 2, .(patid, index_date)]

	controls_cps <- gold_cp[valid_controls, on = "patid", nomatch = 0][
		stop_date >= index_date-84
	]

	controls_with_valid_ltc <- controls_with_valid_ltc[patid %in% unique(controls_cps$patid)]

	# Get all LTC terms for these valid controls
	# Use a more efficient approach - join back with the already filtered controls_ltc
	controls_with_terms <- controls_ltc[controls_with_valid_ltc, on = c("patid", "index_date"), nomatch = 0][
		, .(patid, index_date, term)
	]

	return(controls_with_terms)
}

#' Get eligible controls with extended demographic information
#'
#' @param cases_df data.table with columns: patid, eventdate, study_name
#' @param gold_patient data.table with patient demographic data
#' @param gold_ltc data.table with long-term condition data
#' @param random_seed integer for reproducible results
#' @return data.table with patid, index_date, term, gender, imd_quintile, eth_group, dob
get_eligible_controls_extended <- function(cases_df, gold_patient, gold_ltc, gold_cp, random_seed = 4) {
	# Get the core eligible controls
	controls_core <- compute_eligible_controls_core(cases_df, gold_patient, gold_ltc, gold_cp, random_seed)

	# Join with patient demographic data
	controls_extended <- controls_core[gold_patient, on = "patid", nomatch = 0][
		, .(patid, index_date, term, gender, imd_quintile, eth_group, dob)
	]

	return(controls_extended)
}

#' Save study results to files
#'
#' @param study_name character string
#' @param cases_df data.table with case information
#' @param matched_patids data.table with matched patient IDs and treatment indicator
#' @param study_dir directory to save files (default: "studies")
save_study_results <- function(study_name, cases_df, matched_patids, study_dir = "studies") {
	# Create directory if it doesn't exist
	if (!dir.exists(study_dir)) {
		dir.create(study_dir, recursive = TRUE)
	}

	# Save study metadata
	metadata <- list(
		study_name = study_name,
		n_cases = nrow(cases_df),
		n_controls = nrow(matched_patids[treatment == 0]),
		case_dates = unique(cases_df$eventdate),
		created_date = Sys.Date()
	)

	saveRDS(metadata, file.path(study_dir, paste0(study_name, "_metadata.rds")))

	# Save matched patient IDs
	fwrite(matched_patids, file.path(study_dir, paste0(study_name, "_matched_patids.csv")))

	message("Study '", study_name, "' saved to ", study_dir)
}

#' Load study results from files
#'
#' @param study_name character string
#' @param study_dir directory containing study files (default: "studies")
#' @return list with metadata and matched_patids
load_study_results <- function(study_name, study_dir = "studies") {
	metadata_file <- file.path(study_dir, paste0(study_name, "_metadata.rds"))
	patids_file <- file.path(study_dir, paste0(study_name, "_matched_patids.csv"))

	if (!file.exists(metadata_file) || !file.exists(patids_file)) {
		stop("Study '", study_name, "' not found in ", study_dir)
	}

	metadata <- readRDS(metadata_file)
	matched_patids <- fread(patids_file)

	return(list(
		metadata = metadata,
		matched_patids = matched_patids
	))
}

#' List available studies
#'
#' @param study_dir directory containing study files (default: "studies")
#' @return character vector of study names
list_studies <- function(study_dir = "studies") {
	if (!dir.exists(study_dir)) {
		return(character(0))
	}

	metadata_files <- list.files(study_dir, pattern = "_metadata.rds$", full.names = FALSE)
	study_names <- gsub("_metadata.rds$", "", metadata_files)

	return(study_names)
}

#' Reconstruct matched datasets from saved study and original data.tables
#'
#' @param study_name character string
#' @param gold_patient data.table with patient demographic data
#' @param gold_cp data.table with prescription data
#' @param gold_ltc data.table with long-term condition data
#' @param study_dir directory containing study files (default: "studies")
#' @return list with matched prescription, patient, and LTC data
reconstruct_matched_datasets <- function(study_name, gold_patient, gold_cp, gold_ltc, study_dir = "studies") {
	# Load study results
	study_results <- load_study_results(study_name, study_dir)
	matched_patids <- study_results$matched_patids

	# Get case and control patient IDs
	case_patids <- matched_patids[treatment == 1, patid]
	control_patids <- matched_patids[treatment == 0, patid]
	control_index_dates <- matched_patids[treatment == 0, .(patid, index_date = eventdate)]

	# Reconstruct prescription data
	# Cases: get their actual outcome prescriptions
	cases_presc <- gold_cp[patid %in% case_patids]
	cases_presc[, treatment := 1]

	# Controls: get prescriptions before their assigned index date
	controls_presc <- control_index_dates[gold_cp, on = "patid", nomatch = 0][
		start_date < index_date & stop_date >= index_date - 84  # 84 days lookback
	]
	controls_presc[, treatment := 0]

	# Combine prescription data
	matched_prescriptions <- rbindlist(list(
		cases_presc[, .(patid, substance, start_date, stop_date, duration, treatment)],
		controls_presc[, .(patid, substance, start_date, stop_date, duration, treatment)]
	), fill = TRUE)

	# Reconstruct patient data
	matched_patients <- gold_patient[patid %in% c(case_patids, control_patids)]
	matched_patients[, treatment := fifelse(patid %in% case_patids, 1, 0)]

	# Reconstruct LTC data
	# Cases: get all their LTC data
	cases_ltc <- gold_ltc[patid %in% case_patids]
	cases_ltc[, treatment := 1]

	# Controls: get LTC data before their index date
	controls_ltc <- control_index_dates[gold_ltc, on = "patid", nomatch = 0][
		eventdate < index_date
	]
	controls_ltc[, treatment := 0]

	# Combine LTC data
	matched_ltc <- rbindlist(list(
		cases_ltc[, .(patid, eventdate, age_days, term, treatment)],
		controls_ltc[, .(patid, eventdate, age_days, term, treatment)]
	), fill = TRUE)

	return(list(
		prescriptions = matched_prescriptions,
		patients = matched_patients,
		ltc = matched_ltc,
		metadata = study_results$metadata
	))
}

save_matched_datasets_old <- function(study_name, matched_patids, gold_patient, gold_cp, gold_ltc,
																	cases_df, outcome_prescriptions, study_dir = "studies") {
	# Create directory if it doesn't exist
	if (!dir.exists(study_dir)) {
		dir.create(study_dir, recursive = TRUE)
	}

	# CRITICAL: Validate matched_patids structure
	if (!all(c("patid", "treatment", "eventdate") %in% colnames(matched_patids))) {
		stop("matched_patids must have columns: patid, treatment, eventdate")
	}

	# Get case and control patient IDs with proper verification
	case_patids <- matched_patids[treatment == 1, patid]
	control_patids <- matched_patids[treatment == 0, patid]

	# FIXED: Calculate expected cases correctly (unique patients, not prescription rows)
	n_cases_expected <- uniqueN(cases_df$patid)  # This should be ~83,975
	n_cases_actual <- length(case_patids)
	n_controls_actual <- length(control_patids)

	message(sprintf("Expected cases (unique patients): %d, Actual cases: %d, Controls: %d",
									n_cases_expected, n_cases_actual, n_controls_actual))

	# Get control index dates properly
	control_index_dates <- matched_patids[treatment == 0, .(patid, index_date = eventdate)]

	# Create matched prescriptions
	# Cases: use the filtered outcome_prescriptions that already have correct eventdate and outcome_age
	cases_presc <- outcome_prescriptions[patid %in% case_patids]
	cases_presc[, treatment := 1]

	# Controls: create equivalent structure to outcome_prescriptions
	# First get prescriptions before their index date
	controls_presc_raw <- control_index_dates[gold_cp, on = "patid", nomatch = 0][
		start_date < index_date & stop_date >= index_date - 84  # 84 days lookback
	]

	# Add patient demographics to calculate outcome_age
	controls_dob <- gold_patient[patid %in% control_patids, .(patid, dob)]
	controls_presc <- merge(controls_presc_raw, controls_dob, by = "patid")

	# Create structure matching outcome_prescriptions
	controls_presc[, `:=`(
		eventdate = index_date,
		outcome_age = as.numeric(as.IDate(index_date) - as.IDate(dob)),
		treatment = 0
	)]
	controls_presc[, `:=`(dob = NULL, index_date = NULL)]  # Clean up

	# Combine prescription data with consistent structure
	all_prescriptions <- rbindlist(list(
		cases_presc[, .(patid, eventdate, outcome_age, substance, start_date, stop_date, duration, treatment)],
		controls_presc[, .(patid, eventdate, outcome_age, substance, start_date, stop_date, duration, treatment)]
	), fill = TRUE)
	all_prescriptions[, study_name := study_name]

	# Create matched patient data - FIXED to use only matched patients
	cases_patient_data <- gold_patient[patid %in% case_patids]
	cases_patient_data[, treatment := 1]

	controls_patient_data <- gold_patient[patid %in% control_patids]
	controls_patient_data[, treatment := 0]

	all_patient_data <- rbindlist(list(cases_patient_data, controls_patient_data))
	all_patient_data[, study_name := study_name]

	# CRITICAL VERIFICATION: Check patient counts
	final_cases <- all_patient_data[treatment == 1, .N]
	final_controls <- all_patient_data[treatment == 0, .N]

	message(sprintf("Final dataset - Cases: %d, Controls: %d", final_cases, final_controls))

	# Create matched LTC data
	# Cases: get all their LTC data
	cases_ltc <- gold_ltc[patid %in% case_patids]
	cases_ltc[, treatment := 1]

	# Controls: get LTC data before their index date
	controls_ltc <- control_index_dates[gold_ltc, on = "patid", nomatch = 0][
		eventdate < index_date
	]
	controls_ltc[, treatment := 0]

	# Combine LTC data
	all_ltc <- rbindlist(list(
		cases_ltc[, .(patid, eventdate, age_days, term, treatment)],
		controls_ltc[, .(patid, eventdate, age_days, term, treatment)]
	), fill = TRUE)
	all_ltc[, study_name := study_name]

	# Create study metadata with verification
	metadata <- list(
		study_name = study_name,
		n_cases = final_cases,
		n_controls = final_controls,
		n_cases_expected = n_cases_expected,
		case_dates = unique(cases_df$eventdate),
		created_date = Sys.Date(),
		matching_ratio = final_controls / final_cases
	)

	# Save everything as a single RDS file
	study_data <- list(
		all_prescriptions = all_prescriptions,
		all_patient_data = all_patient_data,
		all_ltc = all_ltc,
		matched_patids = matched_patids,
		metadata = metadata
	)

	saveRDS(study_data, file.path(study_dir, paste0(study_name, "_study_data.rds")))

	message("Study '", study_name, "' saved as RDS file to ", study_dir)
	message("Final verification - Cases: ", final_cases, ", Controls: ", final_controls)
}