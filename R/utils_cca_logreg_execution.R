# R/utils_cca_logreg_execution.R
# Model execution and validation utilities

#' Validate model inputs based on model type
#'
#' @param model_type character string indicating model type
#' @param inputs list of input values from Shiny
#' @return logical TRUE if valid, FALSE otherwise (with notification)
validate_model_inputs <- function(model_type, inputs) {
	switch(model_type,
				 background_main = {
				 	if (is.null(inputs$selected_background_meds) || length(inputs$selected_background_meds) == 0) {
				 		showNotification(
				 			"Please select at least one background medication",
				 			type = "warning",
				 			duration = 5
				 		)
				 		return(FALSE)
				 	}
				 	TRUE
				 },

				 background_pairwise = {
				 	if (is.null(inputs$selected_background_meds) || length(inputs$selected_background_meds) < 2) {
				 		showNotification(
				 			"Please select at least two background medications for pairwise interactions",
				 			type = "warning",
				 			duration = 5
				 		)
				 		return(FALSE)
				 	}
				 	TRUE
				 },

				 recent_pp = {
				 	if (is.null(inputs$selected_recent_presc) || length(inputs$selected_recent_presc) == 0) {
				 		showNotification(
				 			"Please select at least one recent prescription",
				 			type = "warning",
				 			duration = 5
				 		)
				 		return(FALSE)
				 	}
				 	TRUE
				 },

				 recent_main = {
				 	if (is.null(inputs$selected_recent_presc) || length(inputs$selected_recent_presc) == 0) {
				 		showNotification(
				 			"Please select at least one recent prescription",
				 			type = "warning",
				 			duration = 5
				 		)
				 		return(FALSE)
				 	}
				 	TRUE
				 },

				 recent_background = {
				 	if (is.null(inputs$selected_recent_presc) || length(inputs$selected_recent_presc) == 0 ||
				 			is.null(inputs$selected_background_meds) || length(inputs$selected_background_meds) == 0) {
				 		showNotification(
				 			"Please select at least one medication from both recent prescriptions and background medications",
				 			type = "warning",
				 			duration = 5
				 		)
				 		return(FALSE)
				 	}
				 	TRUE
				 }
	)
}

#' Execute recent prescription main effects model
execute_recent_main <- function(medications, group_medications, selected_ltcs,
																selected_covariates, data) {
	fit_recent_main_models(
		selected_ltcs = selected_ltcs,
		selected_covariates = selected_covariates,
		medications = medications,
		patient_data = data$patient_data,
		recent_prescriptions = data$cases_controls,
		ltcs = data$ltcs,
		group_medications = group_medications
	)
}

#' Prepare filtered data for modeling
#'
#' @param inputs list of input values from Shiny
#' @param patient_data data.table of patient data
#' @param prescriptions data.table of prescriptions
#' @param ltcs data.table of LTCs
#' @param cases_controls data.table of cases and controls
#' @return list with filtered datasets
prepare_filtered_data <- function(inputs, patient_data, prescriptions, ltcs, cases_controls) {

	# Start with full datasets
	patient_data_filtered <- copy(patient_data)
	prescriptions_filtered <- copy(prescriptions)
	ltcs_filtered <- copy(ltcs)
	case_controls <- copy(cases_controls)

	# Apply subgroup filtering if selected
	if (isTruthy(inputs$subgroup_filter) && inputs$subgroup_filter != "all") {
		subgroup_parts <- strsplit(inputs$subgroup_filter, "#")[[1]]
		subgroup_var <- subgroup_parts[1]
		subgroup_val <- subgroup_parts[2]

		patient_data_filtered <- patient_data_filtered[get(subgroup_var) == subgroup_val]
		filtered_patids <- patient_data_filtered$patid

		prescriptions_filtered <- prescriptions_filtered[patid %in% filtered_patids]
		ltcs_filtered <- ltcs_filtered[patid %in% filtered_patids]
		case_controls <- case_controls[patid %in% filtered_patids]
	}

	list(
		patient_data = patient_data_filtered,
		prescriptions = prescriptions_filtered,
		ltcs = ltcs_filtered,
		cases_controls = case_controls
	)
}

#' Execute appropriate model based on model type
#'
#' @param model_type character string indicating model type
#' @param inputs list of input values from Shiny
#' @param data list of filtered datasets
#' @param selected_ltc_terms character vector of selected LTC terms
#' @return data.table with model results
execute_model <- function(model_type, inputs, data, selected_ltc_terms) {

	selected_covariates <- inputs$selected_covariates

	switch(model_type,

				 background_main = {
				 	execute_background_main(
				 		medications = inputs$selected_background_meds,
				 		selected_ltcs = selected_ltc_terms,
				 		selected_covariates = selected_covariates,
				 		data = data
				 	)
				 },

				 background_pairwise = {
				 	execute_background_pairwise(
				 		medications = inputs$selected_background_meds,
				 		min_coprescription_prev = inputs$interaction_min_coprescription_prev,
				 		selected_ltcs = selected_ltc_terms,
				 		selected_covariates = selected_covariates,
				 		data = data
				 	)
				 },

				 recent_pp = {
				 	execute_recent_pp(
				 		medications = inputs$selected_recent_presc,
				 		group_medications = inputs$group_recent_meds,
				 		selected_ltcs = selected_ltc_terms,
				 		selected_covariates = selected_covariates,
				 		data = data
				 	)
				 },

				 recent_main = {
				 	execute_recent_main(
				 		medications = inputs$selected_recent_presc,
				 		group_medications = inputs$group_recent_meds,
				 		selected_ltcs = selected_ltc_terms,
				 		selected_covariates = selected_covariates,
				 		data = data
				 	)
				 },

				 recent_background = {
				 	execute_recent_background(
				 		recent_meds = inputs$selected_recent_presc,
				 		background_meds = inputs$selected_background_meds,
				 		group_recent = inputs$group_recent_meds,
				 		selected_ltcs = selected_ltc_terms,
				 		selected_covariates = selected_covariates,
				 		data = data
				 	)
				 }
	)
}

#' Execute background main effects model
execute_background_main <- function(medications, selected_ltcs, selected_covariates, data) {
	run_logistic_models(
		medications = medications,
		selected_ltcs = selected_ltcs,
		selected_covariates = selected_covariates,
		patient_data = data$patient_data,
		prescriptions = data$prescriptions,
		ltcs = data$ltcs
	)
}

#' Execute background pairwise interactions model
execute_background_pairwise <- function(medications, min_coprescription_prev,
																				selected_ltcs, selected_covariates, data) {

	# Filter to medication pairs with sufficient co-prescription
	med_pairs <- create_medication_pairs(medications)

	# Calculate co-prescription prevalence
	presc_wide <- dcast(
		data$prescriptions[substance %in% medications],
		patid ~ substance,
		fun.aggregate = length,
		value.var = "substance"
	)

	# Convert to binary
	for (col in medications) {
		if (col %in% names(presc_wide)) {
			presc_wide[, (col) := as.integer(get(col) > 0)]
		}
	}

	# Merge with patient data to get treatment status
	presc_with_treatment <- merge(
		presc_wide,
		data$patient_data[, .(patid, treatment)],
		by = "patid"
	)

	# Filter pairs by co-prescription prevalence
	filtered_pairs <- list()
	for (pair in med_pairs) {
		med1 <- pair[1]
		med2 <- pair[2]

		if (med1 %in% names(presc_with_treatment) && med2 %in% names(presc_with_treatment)) {
			# Calculate prevalence
			case_prev <- presc_with_treatment[treatment == 1, sum(get(med1) == 1 & get(med2) == 1) / .N * 100]
			control_prev <- presc_with_treatment[treatment == 0, sum(get(med1) == 1 & get(med2) == 1) / .N * 100]

			if (case_prev >= min_coprescription_prev && control_prev >= min_coprescription_prev) {
				filtered_pairs[[length(filtered_pairs) + 1]] <- pair
			}
		}
	}

	if (length(filtered_pairs) == 0) {
		showNotification(
			"No medication pairs meet the minimum co-prescription prevalence threshold",
			type = "warning",
			duration = 5
		)
		return(data.table())
	}

	# Run interaction models
	run_interaction_models(
		med_pairs = filtered_pairs,
		selected_ltcs = selected_ltcs,
		selected_covariates = selected_covariates,
		patient_data = data$patient_data,
		prescriptions = data$prescriptions,
		ltcs = data$ltcs
	)
}

#' Execute recent prescription × PP burden model
execute_recent_pp <- function(medications, group_medications, selected_ltcs,
															selected_covariates, data) {
	fit_pp_interaction_models(
		selected_ltcs = selected_ltcs,
		selected_covariates = selected_covariates,
		medications = medications,
		patient_data = data$patient_data,
		recent_prescriptions = data$cases_controls,
		ltcs = data$ltcs,
		group_medications = group_medications
	)
}

#' Execute recent × background interaction model
execute_recent_background <- function(recent_meds, background_meds, group_recent,
																			selected_ltcs, selected_covariates, data) {
	run_recent_background_interaction_models(
		recent_meds = recent_meds,
		background_meds = background_meds,
		selected_ltcs = selected_ltcs,
		selected_covariates = selected_covariates,
		patient_data = data$patient_data,
		recent_prescriptions = data$cases_controls,
		prescriptions = data$prescriptions,
		ltcs = data$ltcs,
		group_recent = group_recent
	)
}

#' Create all unique medication pairs
#'
#' @param medications character vector of medication names
#' @return list of character vectors, each containing a pair
create_medication_pairs <- function(medications) {
	n <- length(medications)
	pairs <- list()

	if (n < 2) return(pairs)

	for (i in 1:(n-1)) {
		for (j in (i+1):n) {
			pairs[[length(pairs) + 1]] <- c(medications[i], medications[j])
		}
	}

	return(pairs)
}

#' Get success message for model completion
#'
#' @param results data.table with model results
#' @param model_type character string indicating model type
#' @return character string with success message
get_success_message <- function(results, model_type) {
	n_results <- switch(model_type,
											background_main = nrow(results),
											background_pairwise = nrow(results),
											recent_pp = length(unique(results$medication)),
											recent_main = length(unique(results$medication)),
											recent_background = nrow(results)
	)

	sprintf("Models completed: %d results", n_results)
}