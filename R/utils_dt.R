buildWhereConditions <- function(dt,
																 terms = NULL,
																 start_date = NULL,
																 eth_group_f = NULL,
																 sex_f = NULL,
																 imd_quintile_f = NULL,
																 outcome_f = NULL) {

	# Start with the full data.table
	result <- dt

	# Apply filters based on provided parameters
	# if (!is.null(terms)) {
	# 	result <- result[term %in% terms]
	# }

	if (!is.null(start_date)) {
		result <- result[eventdate >= as.Date(start_date)]
	}

	if (!is.null(eth_group_f)) {
		result <- result[eth_group %in% eth_group_f]
	}

	if (!is.null(sex_f)) {
		result <- result[gender %in% sex_f]
	}

	if (!is.null(imd_quintile_f)) {
		result <- result[imd_quintile %in% imd_quintile_f]
	}

	if (!is.null(outcome_f)) {
		result <- result[outcome == outcome_f]
	}

	return(result)
}

buildAcutePrescriptionQuery <- function(gold_acute_presc,gold_patient,gold_outcomes,gold_ltc,terms = NULL,
																				start_date = NULL,
																				eth_group = NULL,
																				sex = NULL,
																				imd_quintile = NULL,
																				outcome = NULL) {

	# Start with gold_outcomes as base
	base_dt <- copy(gold_outcomes)

	# Join with gold_ltc if terms specified
	if (!is.null(terms)) {
		subset_ltc <- gold_ltc[term %in% terms]
		base_dt <- base_dt[subset_ltc, on = "patid", nomatch = 0]
	}

	# Join with gold_patient if demographic filters specified
	if (!is.null(eth_group) || !is.null(sex) || !is.null(imd_quintile)) {
		base_dt <- base_dt[gold_patient, on = "patid", nomatch = 0]
	}

	# Apply filtering conditions
	filtered_dt <- buildWhereConditions(base_dt, terms, start_date, eth_group, sex, imd_quintile, outcome)

	# Get distinct patient IDs
	patient_ids <- unique(filtered_dt$patid)

	# Filter gold_acute_presc for these patients
	result <- gold_acute_presc[patid %in% patient_ids, .(patid, substance, start_date, stop_date)]

	return(result)
}

buildPrescriptionQuery <- function(gold_cp,gold_patient,gold_ltc,gold_outcomes,terms = NULL,
																	 start_date = NULL,
																	 eth_group = NULL,
																	 sex = NULL,
																	 imd_quintile = NULL,
																	 outcome = NULL) {
	# Start with gold_outcomes as base
	base_dt <- copy(gold_outcomes)

	# Join with gold_ltc if terms specified
	if (!is.null(terms)) {
		subset_ltc <- gold_ltc[term %in% terms]
		base_dt <- base_dt[subset_ltc, on = "patid", nomatch = 0]
	}

	# Join with gold_patient if demographic filters specified
	if (!is.null(eth_group) || !is.null(sex) || !is.null(imd_quintile)) {
		base_dt <- base_dt[gold_patient, on = "patid", nomatch = 0]
	}

	# Apply filtering conditions
	filtered_dt <- buildWhereConditions(base_dt, terms, start_date, eth_group, sex, imd_quintile, outcome)

	# Get distinct patient IDs
	patient_ids <- unique(filtered_dt$patid)

	# Filter gold_cp for these patients
	result <- gold_cp[patid %in% patient_ids, .(patid, substance, start_date, stop_date, duration)]

	return(result)
}

#' @export
buildLtcQuery <- function(gold_patient,gold_outcomes,gold_ltc,terms = NULL,
													start_date = NULL,
													eth_group = NULL,
													sex = NULL,
													imd_quintile = NULL,
													outcome_f = NULL) {
	# Start with gold_outcomes and join with gold_ltc
	subset_ltc <- copy(gold_ltc)
	if (!is.null(terms)) {
		subset_patids <- gold_ltc[term %in% terms, patid]
		subset_ltc <- subset_ltc[patid %in% subset_patids]
	}
	base_dt <- gold_outcomes[outcome==outcome_f][subset_ltc, on = "patid", nomatch = 0]

	# Filter where LTC eventdate is <= outcome eventdate
	base_dt <- base_dt[i.eventdate < eventdate]

	# Join with gold_patient if demographic filters specified
	if (!is.null(eth_group) || !is.null(sex) || !is.null(imd_quintile)) {
		base_dt <- base_dt[gold_patient, on = "patid", nomatch = 0]
	}

	# Apply filtering conditions
	filtered_dt <- buildWhereConditions(base_dt, terms, start_date, eth_group, sex, imd_quintile)

	# Select and rename columns as needed
	result <- filtered_dt[order(patid,i.eventdate), .(patid, eventdate=i.eventdate,age_days=i.age_days , term)]

	return(result)
}

#' @export
buildPatientQuery <- function(gold_patient,gold_outcomes,gold_ltc,terms = NULL,
															eth_group = NULL,
															sex = NULL,
															imd_quintile = NULL,
															outcome = NULL) {

	# Start with gold_patient joined with gold_outcomes
	base_dt <- gold_patient[gold_outcomes, on = "patid", nomatch = 0]

	# Join with gold_ltc if terms specified
	if (!is.null(terms)) {
		subset_ltc <- gold_ltc[term %in% terms]
		base_dt <- base_dt[subset_ltc, on = "patid", nomatch = 0]
	}

	# Apply filtering conditions (excluding start_date as it's not relevant here)
	filtered_dt <- buildWhereConditions(base_dt, terms, NULL, eth_group, sex, imd_quintile, outcome)

	# Select distinct patients with required columns
	result <- unique(filtered_dt[, .(patid, dob, sex = gender, eth_group, imd_quintile)])

	return(result)
}

buildOutcomeQuery <- function(gold_outcomes, outcome_f) {
	result <- gold_outcomes[outcome == outcome_f, .(patid, eventdate, age_days, term = outcome)]
	return(result)
}

#' @export
buildOutcomeLtcOptionsQuery <- function(gold_patient,
																				gold_ltc,
																				gold_outcomes,
																				outcome_f,
																				terms = NULL,
																				start_date = NULL,
																				eth_group = NULL,
																				sex = NULL,
																				imd_quintile = NULL) {

	# Start with gold_patient joined with gold_ltc
	subset_ltc <- copy(gold_ltc)
	if (!is.null(terms)) {
		subset_ltc <- gold_ltc[term %in% terms]
	}
	base_dt <- gold_patient[subset_ltc, on = "patid", nomatch = 0]

	# Apply filtering conditions (excluding outcome as it's not relevant here)
	filtered_dt <- buildWhereConditions(base_dt, terms, start_date, eth_group, sex, imd_quintile, NULL)
	outcomes <- copy(gold_outcomes)
	outcomes <- outcomes[outcome == outcome_f, .(patid, eventdate, age_days, term = outcome)]
	outcomes <- outcomes[patid %in% filtered_dt$patid]

	return(outcomes)
}

#' @export
buildLtcDatesQuery <- function(ltc_name) {
	# Filter gold_ltc for the specific term and patients in gold_cp
	cp_patients <- unique(gold_cp$patid)
	result <- gold_ltc[term == ltc_name & patid %in% cp_patients, .(patid, eventdate, age_days)]
	return(result)
}