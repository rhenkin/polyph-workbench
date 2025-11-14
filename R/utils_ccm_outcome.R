#' Outcome-Based Case-Control Matching Workflow
#'
#' This workflow creates matched case-control cohorts where cases are defined
#' by outcome events alone, without requiring a new prescription within a
#' prediction window. Controls are matched on demographics and multimorbidity
#' status, then polypharmacy is measured as an exposure variable.

#' Main workflow for outcome-based matching
#'
#' @param outcome_events data.table with patid and eventdate (outcome date)
#' @param patient_data data.table of patient demographics (cases only)
#' @param ltc_data data.table of LTC records (cases only)
#' @param match_ratio numeric control:case ratio
#' @param index_offset numeric days before outcome to set as index date (default 0)
#' @param patient_filters list of patient filter criteria
#' @return list(cases, controls)
#' @note Uses global variables for controls: gold_patient, gold_cp, gold_ltc
create_outcome_based_cohort_workflow <- function(outcome_events,
																								 patient_data,
																								 ltc_data,
																								 match_ratio = 4,
																								 index_offset = 0,
																								 progress,
																								 patient_filters = NULL) {

	# Step 1: Build cases table with index dates (uses filtered reactive data)
	progress$set(message = "Building cases...", value = 0.2)
	cases <- build_outcome_based_cases(
		outcome_events = outcome_events,
		patient_data = patient_data,
		ltc_data = ltc_data,
		index_offset = index_offset,
		patient_filters = patient_filters
	)

	# Step 2: Find eligible control pool (uses global gold_* variables)
	progress$set(message = "Finding eligible controls...", value = 0.5)
	eligible_controls <- find_outcome_based_eligible_controls(
		cases = cases,
		patient_filters = patient_filters
	)

	# Step 3: Sample controls by strata
	progress$set(message = "Sampling controls...", value = 0.7)
	controls <- sample_outcome_based_controls(
		eligible_controls = eligible_controls,
		cases = cases,
		match_ratio = match_ratio
	)

	# Step 4: Measure polypharmacy for both cases and controls using global gold_cp
	progress$set(message = "Measuring polypharmacy...", value = 0.9)
	cases <- measure_polypharmacy_at_index(cases, gold_cp)
	controls <- measure_polypharmacy_at_index(controls, gold_cp)

	return(list(
		cases = cases,
		controls = controls
	))
}

#' Build cases table for outcome-based matching
#'
#' @param outcome_events data.table with patid, eventdate
#' @param patient_data data.table with patient demographics (cases only)
#' @param ltc_data data.table with LTC records (cases only)
#' @param index_offset numeric days before outcome for index date
#' @param patient_filters list of filter criteria
#' @return data.table of cases with index_date and stratification variables
build_outcome_based_cases <- function(outcome_events,
																			patient_data,
																			ltc_data,
																			index_offset = 0,
																			patient_filters = NULL) {

	# Start with unique outcome events
	cases <- unique(outcome_events[, .(patid, outcome_date = eventdate)])

	# Set index date (outcome date or offset before)
	cases[, index_date := outcome_date - index_offset]

	# Add multimorbidity information (uses filtered ltc_data for cases)
	cases <- add_outcome_based_multimorbidity(
		cases = cases,
		ltc_data = ltc_data,
		patient_filters = patient_filters
	)

	# Add demographics from filtered patient_data (cases only)
	cases <- merge(cases,
								 patient_data[, .(patid, dob, sex, eth_group, imd_quintile)],
								 by = "patid")

	# Calculate age at index
	cases[, index_age := as.numeric(index_date - dob) / 365.25]

	# Calculate time at risk (time from multimorbidity onset to index)
	cases[, time_at_risk := round(as.numeric(index_date - mm_date) / 365.25, digits = 0)]

	# Create stratification variables
	cases <- create_outcome_based_strata(cases)

	return(cases)
}

#' Add multimorbidity information for outcome-based cases
#'
#' @param cases data.table with patid and index_date
#' @param ltc_data data.table with LTC records
#' @param patient_filters list with selected_ltcs if applicable
#' @return data.table with mm_date added
add_outcome_based_multimorbidity <- function(cases, ltc_data, patient_filters = NULL) {

	# Get LTCs before index date
	valid_ltcs <- cases[ltc_data,
											.(patid, term, eventdate = i.eventdate, index_date = x.index_date),
											on = .(patid, index_date > eventdate),
											nomatch = 0]

	setkey(valid_ltcs, patid, eventdate)

	# Filter by selected LTCs if specified
	if (!is.null(patient_filters$selected_ltcs)) {
		ltc_counts <- valid_ltcs[term %in% patient_filters$selected_ltcs,
														 .(n_selected_ltcs = uniqueN(term)),
														 by = patid]

		n_required <- length(patient_filters$selected_ltcs)
		patids_with_all_ltcs <- ltc_counts[n_selected_ltcs == n_required, patid]

		valid_ltcs <- valid_ltcs[patid %in% patids_with_all_ltcs]
		cases <- cases[patid %in% patids_with_all_ltcs]
	}

	# Filter to patients with ≥2 LTCs
	valid_ltc_patids <- valid_ltcs[, .N, patid][N >= 2, patid]
	valid_ltcs <- valid_ltcs[patid %in% valid_ltc_patids]

	# Find date of 2nd LTC (multimorbidity date)
	setorder(valid_ltcs, patid, eventdate)
	valid_ltcs[, ltc_index := 1:.N, patid]
	mm_dates <- valid_ltcs[ltc_index == 2, .(patid, mm_date = eventdate)]

	# Merge back to cases
	cases <- merge(cases, mm_dates, by = "patid")

	return(cases)
}

#' Create stratification variables for outcome-based matching
#'
#' @param cases data.table with demographic and temporal variables
#' @return data.table with strata variable added
create_outcome_based_strata <- function(cases) {

	# Bin time at risk (time since multimorbidity)
	cases[, mm_duration_bin := cut(
		time_at_risk,
		breaks = seq(0, 121, 2),
		include.lowest = TRUE,
		right = FALSE
	)]

	# Bin age
	cases[, age_bin := cut(
		index_age,
		breaks = seq(0, 120, 5),
		include.lowest = TRUE,
		right = FALSE
	)]

	# Extract year
	cases[, year := year(index_date)]

	# Create composite strata variable
	cases[, strata := paste(sex, age_bin, mm_duration_bin, year, sep = "_")]

	return(cases)
}

#' Find eligible controls for outcome-based matching
#'
#' @param cases data.table of cases with strata
#' @param patient_filters list of filter criteria
#' @return data.table of eligible controls with assigned index dates
#' @note Uses global variables: gold_patient, gold_ltc (for controls)
find_outcome_based_eligible_controls <- function(cases,
																								 patient_filters = NULL) {

	case_patids <- unique(cases$patid)

	# Get all potential controls (not cases) from gold_patient (global)
	potential_controls <- gold_patient[!patid %in% case_patids]

	# Apply basic patient filters
	if (!is.null(patient_filters)) {
		input_list <- patient_filters$input_list

		if (!is.null(input_list$sex) && length(input_list$sex) > 0) {
			potential_controls <- potential_controls[sex %in% input_list$sex]
		}

		if (!is.null(input_list$eth_group) && length(input_list$eth_group) > 0) {
			potential_controls <- potential_controls[eth_group %in% input_list$eth_group]
		}

		if (!is.null(patient_filters$imd_quintile) && length(patient_filters$imd_quintile) > 0) {
			imd_values <- as.numeric(patient_filters$imd_quintile)
			potential_controls <- potential_controls[imd_quintile %in% imd_values]
		}
	}

	# Get unique strata information from cases
	strata_info <- unique(cases[, .(strata, index_date, sex, age_bin, mm_duration_bin, year)])

	# Prepare control candidates with necessary columns
	control_candidates <- potential_controls[, .(patid, dob, sex, eth_group, imd_quintile, tod, dod)]

	# Filter: Must be alive and in study at each strata's index date
	# Calculate age for ALL possible index dates first (pre-compute age bins)
	setkey(control_candidates, patid)
	setkey(strata_info, sex)

	# Process each unique strata efficiently using data.table
	# Group strata by sex to reduce redundant filtering
	strata_by_sex <- split(strata_info, by = "sex")

	eligible_pool_list <- lapply(names(strata_by_sex), function(sex_val) {
		strata_subset <- strata_by_sex[[sex_val]]
		controls_subset <- control_candidates[sex == sex_val]

		if (nrow(controls_subset) == 0) return(NULL)

		# For each strata in this sex group, find matching controls
		strata_results <- lapply(1:nrow(strata_subset), function(i) {
			strata_row <- strata_subset[i]

			# Filter controls for this specific strata
			matched <- controls_subset[
				# Must be alive and in study at this index date
				(is.na(dod) | as.IDate(dod) > strata_row$index_date) &
					(is.na(tod) | as.IDate(tod) > strata_row$index_date)
			]

			if (nrow(matched) == 0) return(NULL)

			# Calculate age at this specific index date
			matched[, ':='(
				index_date = strata_row$index_date,
				strata = strata_row$strata
			)]
			matched[, index_age := as.numeric(index_date - dob) / 365.25]
			matched[, age_bin := cut(
				index_age,
				breaks = seq(0, 120, 5),
				include.lowest = TRUE,
				right = FALSE
			)]

			# Filter to matching age bin
			matched <- matched[age_bin == strata_row$age_bin]

			if (nrow(matched) == 0) return(NULL)

			return(matched)
		})

		# Combine results for this sex
		rbindlist(strata_results, fill = TRUE)
	})

	# Combine all results
	eligible_pool <- rbindlist(eligible_pool_list, fill = TRUE)

	if (nrow(eligible_pool) == 0) {
		message("No eligible controls found after demographic matching")
		return(data.table())
	}

	message(sprintf("Potential control-strata combinations before MM validation: %d", nrow(eligible_pool)))

	# Now validate multimorbidity status and duration using gold_ltc (global)
	eligible_pool <- validate_outcome_based_multimorbidity(
		eligible_pool = eligible_pool,
		patient_filters = patient_filters
	)

	message(sprintf("Eligible control-strata combinations after MM validation: %d", nrow(eligible_pool)))

	return(eligible_pool)
}

#' Validate multimorbidity status for eligible controls
#'
#' @param eligible_pool data.table of potential controls with index_date
#' @param patient_filters list with selected_ltcs if applicable
#' @return data.table with only controls who had ≥2 LTCs before index_date
#' @note Uses global variable: gold_ltc (for controls)
validate_outcome_based_multimorbidity <- function(eligible_pool,
																									patient_filters = NULL) {

	# Get LTCs before each control's index date from gold_ltc (global)
	# Use efficient join
	setkey(eligible_pool, patid, index_date)
	setkey(gold_ltc, patid, eventdate)

	valid_ltcs <- eligible_pool[gold_ltc,
															.(patid, index_date = x.index_date, strata = x.strata,
																term = i.term, eventdate = i.eventdate),
															on = .(patid, index_date > eventdate),
															nomatch = 0]

	setkey(valid_ltcs, patid, index_date, strata, eventdate)

	# Filter by selected LTCs if specified
	if (!is.null(patient_filters$selected_ltcs)) {
		# Count how many DISTINCT selected LTCs each control-strata combo has
		ltc_counts <- valid_ltcs[term %in% patient_filters$selected_ltcs,
														 .(n_selected_ltcs = uniqueN(term)),
														 by = .(patid, index_date, strata)]

		n_required <- length(patient_filters$selected_ltcs)
		valid_combinations <- ltc_counts[n_selected_ltcs == n_required,
																		 .(patid, index_date, strata)]

		# Filter to only those with all selected LTCs
		valid_ltcs <- valid_ltcs[valid_combinations,
														 on = .(patid, index_date, strata),
														 nomatch = 0]
	}

	# Filter to control-strata combinations with ≥2 LTCs
	ltc_count_by_control <- valid_ltcs[, .(n_ltcs = uniqueN(term)),
																		 by = .(patid, index_date, strata)]
	valid_controls <- ltc_count_by_control[n_ltcs >= 2, .(patid, index_date, strata)]

	# Find multimorbidity date (2nd LTC) for each control-strata combination
	setorder(valid_ltcs, patid, index_date, strata, eventdate)
	valid_ltcs[, ltc_index := 1:.N, by = .(patid, index_date, strata)]
	mm_dates <- valid_ltcs[ltc_index == 2, .(patid, index_date, strata, mm_date = eventdate)]

	# Calculate time at risk
	mm_dates[, time_at_risk := round(as.numeric(index_date - mm_date) / 365.25, digits = 0)]

	# Calculate mm_duration_bin
	mm_dates[, mm_duration_bin := cut(
		time_at_risk,
		breaks = seq(0, 121, 2),
		include.lowest = TRUE,
		right = FALSE
	)]

	# Extract mm_duration_bin from strata to verify match
	mm_dates[, mm_duration_bin_from_strata := gsub("^[^_]+_[^_]+_([^_]+)_.*$", "\\1", strata)]

	# Keep only controls with matching mm_duration_bin
	valid_controls_with_mm <- mm_dates[
		as.character(mm_duration_bin) == mm_duration_bin_from_strata,
		.(patid, index_date, strata, mm_date, time_at_risk)
	]

	# Merge back to eligible pool
	eligible_pool <- eligible_pool[valid_controls_with_mm,
																 on = .(patid, index_date, strata),
																 nomatch = 0]

	return(eligible_pool)
}

#' Sample controls within strata
#'
#' @param eligible_controls data.table of eligible controls
#' @param cases data.table of cases
#' @param match_ratio numeric desired control:case ratio
#' @return data.table of sampled controls
sample_outcome_based_controls <- function(eligible_controls, cases, match_ratio) {

	# Calculate controls needed per stratum
	strata_needs <- cases[, .(n_cases = uniqueN(patid),
														n_controls_needed = uniqueN(patid) * match_ratio),
												by = strata]

	eligible_controls[strata_needs, n_needed := i.n_controls_needed, on = "strata"]

	# Sample within each stratum
	sampled_controls <- eligible_controls[
		!is.na(n_needed),
		{
			n_sample <- min(.N, n_needed[1])
			if (n_sample > 0) .SD[sample(.N, n_sample, replace = FALSE)] else .SD[0]
		},
		by = strata
	]

	sampled_controls[, n_needed := NULL]

	# Remove duplicate patients (can only be control once)
	setorder(sampled_controls, patid, index_date)
	sampled_controls <- sampled_controls[, .SD[sample(.N, 1)], by = patid]

	message(sprintf("Controls sampled: %d (unique patients: %d)",
									nrow(sampled_controls), uniqueN(sampled_controls$patid)))

	return(sampled_controls)
}

#' Measure polypharmacy at index date
#'
#' @param cohort data.table with patid and index_date (cases or controls)
#' @param gold_cp data.table of prescription records
#' @param lookback_days numeric days before index to consider prescriptions active
#' @return data.table with polypharmacy count added
measure_polypharmacy_at_index <- function(cohort, gold_cp, lookback_days = 84) {

	# Find prescriptions active at index date
	active_prescriptions <- cohort[gold_cp,
																 .(patid, index_date = x.index_date,
																 	substance = i.substance,
																 	start_date = i.start_date,
																 	stop_date = i.stop_date),
																 on = .(patid,
																 			 index_date >= start_date,
																 			 index_date <= stop_date + lookback_days),
																 nomatch = 0]

	# Count unique substances per patient
	pp_counts <- active_prescriptions[, .(polypharmacy_count = uniqueN(substance)),
																		by = .(patid, index_date)]

	# Merge back to cohort (0 for those with no prescriptions)
	cohort <- merge(cohort, pp_counts, by = c("patid", "index_date"), all.x = TRUE)
	cohort[is.na(polypharmacy_count), polypharmacy_count := 0]

	return(cohort)
}