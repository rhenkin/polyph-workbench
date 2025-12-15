# R/utils_cca_logreg.R

#' Run logistic regression models for multiple medications
#'
#' @param medications character vector of medication names
#' @param selected_ltcs character vector of LTC terms to include as covariates
#' @param patient_data data.table with patid, treatment, and demographic info
#' @param prescriptions data.table with patid, substance, group
#' @param ltcs data.table with patid, term, group
#' @return data.table with model results
run_logistic_models <- function(medications, selected_ltcs, selected_covariates,
																						patient_data, prescriptions, ltcs) {

	# Strip asterisks from names (added by significance testing)
	medications <- gsub("\\*$", "", medications)
	selected_ltcs <- gsub("\\*$", "", selected_ltcs)

	# Prepare the base dataset
	model_data <- prepare_logreg_data(
		patient_data = patient_data,
		prescriptions = prescriptions,
		ltcs = ltcs,
		medications = medications,
		selected_ltcs = selected_ltcs
	)

	# Fit models for each medication
	results_list <- lapply(medications, function(med) {
		fit_single_logreg_model(
			model_data = model_data,
			medication = med,
			selected_ltcs = selected_ltcs,
			selected_covariates = selected_covariates
		)
	})

	# Combine results
	results <- rbindlist(results_list, fill = TRUE)

	return(results)
}

#' Run interaction models for medication pairs
#'
#' @param med_pairs list of character vectors, each with 2 medication names
#' @param selected_ltcs character vector of LTC terms to include as covariates
#' @param patient_data data.table with patid, treatment, and demographic info
#' @param prescriptions data.table with patid, substance, group
#' @param ltcs data.table with patid, term, group
#' @return data.table with interaction model results
run_interaction_models <- function(med_pairs, selected_ltcs, selected_covariates,
																	 patient_data, prescriptions, ltcs) {

	# Strip asterisks from names
	selected_ltcs <- gsub("\\*$", "", selected_ltcs)

	# Get all unique medications from pairs
	all_meds <- unique(unlist(med_pairs))
	all_meds <- gsub("\\*$", "", all_meds)

	# Prepare the base dataset with all medications
	model_data <- prepare_logreg_data(
		patient_data = patient_data,
		prescriptions = prescriptions,
		ltcs = ltcs,
		medications = all_meds,
		selected_ltcs = selected_ltcs
	)

	# Fit interaction model for each pair
	results_list <- lapply(med_pairs, function(pair) {
		med1 <- gsub("\\*$", "", pair[1])
		med2 <- gsub("\\*$", "", pair[2])

		fit_interaction_model(
			model_data = model_data,
			med1 = med1,
			med2 = med2,
			selected_ltcs = selected_ltcs,
			selected_covariates = selected_covariates
		)
	})

	# Combine results
	results <- rbindlist(results_list, fill = TRUE)

	return(results)
}

#' Fit interaction model for a medication pair (WITH MAIN EFFECTS AND COMBINED EFFECT)
#'
#' @param model_data data.table prepared for modeling
#' @param med1 character string of first medication name
#' @param med2 character string of second medication name
#' @param selected_ltcs character vector of LTC terms
#' @return data.table with interaction model results INCLUDING main effects and combined effect
fit_interaction_model <- function(model_data, med1, med2, selected_ltcs, selected_covariates) {

	# Get column names
	med1_col <- paste0("med_", make.names(med1))
	med2_col <- paste0("med_", make.names(med2))
	ltc_cols <- paste0("ltc_", make.names(selected_ltcs))

	# Check if medication columns exist
	if (!med1_col %in% names(model_data) || !med2_col %in% names(model_data)) {
		return(data.table(
			med1 = med1,
			med2 = med2,
			med1_OR = NA_real_,
			med1_CI_lower = NA_real_,
			med1_CI_upper = NA_real_,
			med1_p = NA_real_,
			med2_OR = NA_real_,
			med2_CI_lower = NA_real_,
			med2_CI_upper = NA_real_,
			med2_p = NA_real_,
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
			combined_OR = NA_real_,
			combined_CI_lower = NA_real_,
			combined_CI_upper = NA_real_,
			pct_cases_both = NA_real_,
			pct_controls_both = NA_real_,
			n_cases_both = NA_integer_,
			n_controls_both = NA_integer_,
			n_ltc_covariates = length(selected_ltcs),
			convergence = "medication_not_found"
		))
	}

	# Check for variation in both medications
	if (model_data[, uniqueN(get(med1_col))] < 2 || model_data[, uniqueN(get(med2_col))] < 2) {
		n_cases_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 1, .N]
		n_controls_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 0, .N]
		total_cases <- model_data[treatment == 1, .N]
		total_controls <- model_data[treatment == 0, .N]
		pct_cases_both <- round(100 * n_cases_both / total_cases, 2)
		pct_controls_both <- round(100 * n_controls_both / total_controls, 2)

		return(data.table(
			med1 = med1,
			med2 = med2,
			med1_OR = NA_real_,
			med1_CI_lower = NA_real_,
			med1_CI_upper = NA_real_,
			med1_p = NA_real_,
			med2_OR = NA_real_,
			med2_CI_lower = NA_real_,
			med2_CI_upper = NA_real_,
			med2_p = NA_real_,
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
			combined_OR = NA_real_,
			combined_CI_lower = NA_real_,
			combined_CI_upper = NA_real_,
			pct_cases_both = pct_cases_both,
			pct_controls_both = pct_controls_both,
			n_cases_both = n_cases_both,
			n_controls_both = n_controls_both,
			n_ltc_covariates = length(selected_ltcs),
			convergence = "no_variation_in_medication"
		))
	}

	# Keep only LTC columns that exist and have variation
	existing_ltc_cols <- ltc_cols[ltc_cols %in% names(model_data)]

	ltc_cols_with_variation <- character(0)
	for (col in existing_ltc_cols) {
		if (model_data[, uniqueN(get(col))] > 1) {
			ltc_cols_with_variation <- c(ltc_cols_with_variation, col)
		}
	}

	# Build formula with interaction
	covariate_terms <- c(med1_col, med2_col, ltc_cols_with_variation)
	# Add selected covariates if any
	if (!is.null(selected_covariates) && length(selected_covariates) > 0) {
		covariate_terms <- c(covariate_terms, selected_covariates)
	}
	interaction_term <- paste0(med1_col, ":", med2_col)
	model_data$group <- factor(model_data$strata)
	formula_str <- paste0("treatment ~ ", paste(covariate_terms, collapse = " + "),
												" + ", interaction_term, " + factor(group)")

	# Fit model
	tryCatch({
		model <- glm(
			formula = as.formula(formula_str),
			data = model_data,
			family = binomial(link = "logit")
		)

		# Extract results for ALL terms
		coef_summary <- summary(model)$coefficients

		# EXTRACT MED1 MAIN EFFECT
		if (med1_col %in% rownames(coef_summary)) {
			med1_row <- coef_summary[med1_col, , drop = FALSE]
			med1_or <- exp(med1_row[1, "Estimate"])
			med1_se <- med1_row[1, "Std. Error"]
			med1_ci_lower <- exp(med1_row[1, "Estimate"] - 1.96 * med1_se)
			med1_ci_upper <- exp(med1_row[1, "Estimate"] + 1.96 * med1_se)
			med1_p <- med1_row[1, "Pr(>|z|)"]
		} else {
			med1_or <- NA_real_
			med1_ci_lower <- NA_real_
			med1_ci_upper <- NA_real_
			med1_p <- NA_real_
		}

		# EXTRACT MED2 MAIN EFFECT
		if (med2_col %in% rownames(coef_summary)) {
			med2_row <- coef_summary[med2_col, , drop = FALSE]
			med2_or <- exp(med2_row[1, "Estimate"])
			med2_se <- med2_row[1, "Std. Error"]
			med2_ci_lower <- exp(med2_row[1, "Estimate"] - 1.96 * med2_se)
			med2_ci_upper <- exp(med2_row[1, "Estimate"] + 1.96 * med2_se)
			med2_p <- med2_row[1, "Pr(>|z|)"]
		} else {
			med2_or <- NA_real_
			med2_ci_lower <- NA_real_
			med2_ci_upper <- NA_real_
			med2_p <- NA_real_
		}

		# EXTRACT INTERACTION TERM
		interaction_row_name <- paste0(med1_col, ":", med2_col)
		if (!interaction_row_name %in% rownames(coef_summary)) {
			interaction_row_name <- paste0(med2_col, ":", med1_col)
		}

		if (!interaction_row_name %in% rownames(coef_summary)) {
			n_cases_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 1, .N]
			n_controls_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 0, .N]
			total_cases <- model_data[treatment == 1, .N]
			total_controls <- model_data[treatment == 0, .N]
			pct_cases_both <- round(100 * n_cases_both / total_cases, 2)
			pct_controls_both <- round(100 * n_controls_both / total_controls, 2)

			return(data.table(
				med1 = med1,
				med2 = med2,
				med1_OR = med1_or,
				med1_CI_lower = med1_ci_lower,
				med1_CI_upper = med1_ci_upper,
				med1_p = med1_p,
				med2_OR = med2_or,
				med2_CI_lower = med2_ci_lower,
				med2_CI_upper = med2_ci_upper,
				med2_p = med2_p,
				interaction_OR = NA_real_,
				interaction_CI_lower = NA_real_,
				interaction_CI_upper = NA_real_,
				interaction_p = NA_real_,
				combined_OR = NA_real_,
				combined_CI_lower = NA_real_,
				combined_CI_upper = NA_real_,
				pct_cases_both = pct_cases_both,
				pct_controls_both = pct_controls_both,
				n_cases_both = n_cases_both,
				n_controls_both = n_controls_both,
				n_ltc_covariates = length(ltc_cols_with_variation),
				convergence = "interaction_term_not_found"
			))
		}

		interaction_row <- coef_summary[interaction_row_name, , drop = FALSE]

		# Calculate OR and CI for interaction
		int_or <- exp(interaction_row[1, "Estimate"])
		int_se <- interaction_row[1, "Std. Error"]
		int_ci_lower <- exp(interaction_row[1, "Estimate"] - 1.96 * int_se)
		int_ci_upper <- exp(interaction_row[1, "Estimate"] + 1.96 * int_se)
		int_p <- interaction_row[1, "Pr(>|z|)"]

		# CALCULATE COMBINED EFFECT (med1 * med2 * interaction)
		combined_or <- med1_or * med2_or * int_or

		# Calculate combined CI using delta method approximation
		# Variance on log scale: Var(log(combined)) = Var(log(med1)) + Var(log(med2)) + Var(log(int))
		combined_log_se <- sqrt(med1_se^2 + med2_se^2 + int_se^2)
		combined_log_or <- log(combined_or)
		combined_ci_lower <- exp(combined_log_or - 1.96 * combined_log_se)
		combined_ci_upper <- exp(combined_log_or + 1.96 * combined_log_se)

		# Count cases and controls with BOTH medications
		n_cases_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 1, .N]
		n_controls_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 0, .N]

		total_cases <- model_data[treatment == 1, .N]
		total_controls <- model_data[treatment == 0, .N]

		pct_cases_both <- round(100 * n_cases_both / total_cases, 2)
		pct_controls_both <- round(100 * n_controls_both / total_controls, 2)

		data.table(
			med1 = med1,
			med2 = med2,
			med1_OR = med1_or,
			med1_CI_lower = med1_ci_lower,
			med1_CI_upper = med1_ci_upper,
			med1_p = med1_p,
			med2_OR = med2_or,
			med2_CI_lower = med2_ci_lower,
			med2_CI_upper = med2_ci_upper,
			med2_p = med2_p,
			interaction_OR = int_or,
			interaction_CI_lower = int_ci_lower,
			interaction_CI_upper = int_ci_upper,
			interaction_p = int_p,
			combined_OR = combined_or,
			combined_CI_lower = combined_ci_lower,
			combined_CI_upper = combined_ci_upper,
			pct_cases_both = pct_cases_both,
			pct_controls_both = pct_controls_both,
			n_cases_both = n_cases_both,
			n_controls_both = n_controls_both,
			n_ltc_covariates = length(ltc_cols_with_variation),
			convergence = "success"
		)

	}, error = function(e) {
		n_cases_both_err <- tryCatch(
			model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 1, .N],
			error = function(e2) NA_integer_
		)
		n_controls_both_err <- tryCatch(
			model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 0, .N],
			error = function(e2) NA_integer_
		)

		total_cases <- tryCatch(model_data[treatment == 1, .N], error = function(e2) NA_integer_)
		total_controls <- tryCatch(model_data[treatment == 0, .N], error = function(e2) NA_integer_)

		pct_cases_both_err <- if (!is.na(n_cases_both_err) && !is.na(total_cases)) {
			round(100 * n_cases_both_err / total_cases, 2)
		} else {
			NA_real_
		}

		pct_controls_both_err <- if (!is.na(n_controls_both_err) && !is.na(total_controls)) {
			round(100 * n_controls_both_err / total_controls, 2)
		} else {
			NA_real_
		}

		data.table(
			med1 = med1,
			med2 = med2,
			med1_OR = NA_real_,
			med1_CI_lower = NA_real_,
			med1_CI_upper = NA_real_,
			med1_p = NA_real_,
			med2_OR = NA_real_,
			med2_CI_lower = NA_real_,
			med2_CI_upper = NA_real_,
			med2_p = NA_real_,
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
			combined_OR = NA_real_,
			combined_CI_lower = NA_real_,
			combined_CI_upper = NA_real_,
			pct_cases_both = pct_cases_both_err,
			pct_controls_both = pct_controls_both_err,
			n_cases_both = n_cases_both_err,
			n_controls_both = n_controls_both_err,
			n_ltc_covariates = length(ltc_cols_with_variation),
			convergence = paste("error:", substr(e$message, 1, 50))
		)
	})
}

#' Prepare data for logistic regression
#'
#' @param patient_data data.table with patient info (has treatment, sex, pp_group, mltc_group)
#' @param prescriptions data.table with medications
#' @param ltcs data.table with long-term conditions
#' @param medications character vector of medications to include
#' @param selected_ltcs character vector of LTCs to include
#' @return data.table in wide format ready for modeling
prepare_logreg_data <- function(patient_data, prescriptions, ltcs,
																medications, selected_ltcs) {

	# Start with patient data
	base_data <- copy(patient_data[, .(patid, treatment, eth_group, imd_quintile, sex, pp_group, mltc_group, strata)])

	# Create binary indicators for medications
	med_indicators <- create_medication_indicators(prescriptions, medications)

	# Create binary indicators for LTCs
	ltc_indicators <- create_ltc_indicators_logreg(ltcs, selected_ltcs)

	# Merge everything
	model_data <- merge(base_data, med_indicators, by = "patid", all.x = TRUE)
	model_data <- merge(model_data, ltc_indicators, by = "patid", all.x = TRUE)

	# Replace NA with 0 for indicators
	med_cols <- paste0("med_", make.names(medications))
	ltc_cols <- paste0("ltc_", make.names(selected_ltcs))

	for (col in c(med_cols, ltc_cols)) {
		if (col %in% names(model_data)) {
			model_data[is.na(get(col)), (col) := 0]
		}
	}

	# Convert factors to ensure proper contrasts
	model_data[, sex := factor(sex)]
	model_data[, eth_group := factor(eth_group)]
	model_data[, pp_group := factor(pp_group)]
	model_data[, mltc_group := factor(mltc_group)]

	return(model_data)
}

#' Create binary medication indicators
#'
#' @param prescriptions data.table with patid and substance
#' @param medications character vector of medications
#' @return data.table with patid and binary indicators
create_medication_indicators <- function(prescriptions, medications) {

	# Filter to selected medications
	presc_subset <- prescriptions[substance %in% medications, .(patid, substance)]

	if (nrow(presc_subset) == 0) {
		# Return empty data.table with just patid
		return(data.table(patid = unique(prescriptions$patid)))
	}

	# Create binary indicators
	presc_subset[, value := 1]

	# Cast to wide format
	med_wide <- dcast(
		presc_subset,
		patid ~ substance,
		value.var = "value",
		fun.aggregate = function(x) as.integer(length(x) > 0),
		fill = 0
	)

	# Rename columns to valid R names with prefix
	old_names <- names(med_wide)[-1]

	if (length(old_names) > 0) {
		new_names <- paste0("med_", make.names(old_names))
		setnames(med_wide, old = old_names, new = new_names)
	}

	return(med_wide)
}

#' Create binary LTC indicators for logistic regression
#'
#' @param ltcs data.table with patid and term
#' @param selected_ltcs character vector of LTC terms
#' @return data.table with patid and binary indicators
create_ltc_indicators_logreg <- function(ltcs, selected_ltcs) {

	# Filter to selected LTCs
	ltcs_subset <- ltcs[term %in% selected_ltcs, .(patid, term)]

	if (nrow(ltcs_subset) == 0) {
		# Return empty data.table with just patid
		return(data.table(patid = unique(ltcs$patid)))
	}

	# Create binary indicators
	ltcs_subset[, value := 1]

	# Cast to wide format
	ltc_wide <- dcast(
		ltcs_subset,
		patid ~ term,
		value.var = "value",
		fun.aggregate = function(x) as.integer(length(x) > 0),
		fill = 0
	)

	# Rename columns to valid R names with prefix
	old_names <- names(ltc_wide)[-1]

	if (length(old_names) > 0) {
		new_names <- paste0("ltc_", make.names(old_names))
		setnames(ltc_wide, old = old_names, new = new_names)
	}

	return(ltc_wide)
}

#' Fit logistic regression for a single medication
#'
#' @param model_data data.table prepared for modeling
#' @param medication character string of medication name
#' @param selected_ltcs character vector of LTC terms
#' @return data.table with model results (single row)
fit_single_logreg_model <- function(model_data, medication, selected_ltcs, selected_covariates) {

	# Get column names
	med_col <- paste0("med_", make.names(medication))
	ltc_cols <- paste0("ltc_", make.names(selected_ltcs))

	# Check if medication column exists
	if (!med_col %in% names(model_data)) {
		return(data.table(
			medication = medication,
			OR = NA_real_,
			CI_lower = NA_real_,
			CI_upper = NA_real_,
			p_value = NA_real_,
			n_cases = NA_integer_,
			n_controls = NA_integer_,
			n_ltc_covariates = length(selected_ltcs),
			convergence = "medication_not_found"
		))
	}

	# Check for variation in medication
	if (model_data[, uniqueN(get(med_col))] < 2) {
		return(data.table(
			medication = medication,
			OR = NA_real_,
			CI_lower = NA_real_,
			CI_upper = NA_real_,
			p_value = NA_real_,
			n_cases = model_data[get(med_col) == 1 & treatment == 1, .N],
			n_controls = model_data[get(med_col) == 1 & treatment == 0, .N],
			n_ltc_covariates = length(selected_ltcs),
			convergence = "no_variation_in_medication"
		))
	}

	# Keep only LTC columns that exist and have variation
	existing_ltc_cols <- ltc_cols[ltc_cols %in% names(model_data)]

	ltc_cols_with_variation <- character(0)
	for (col in existing_ltc_cols) {
		if (model_data[, uniqueN(get(col))] > 1) {
			ltc_cols_with_variation <- c(ltc_cols_with_variation, col)
		}
	}

	# Build formula: treatment ~ medication + ltcs + sex + pp_group + mltc_group
	covariate_terms <- c(med_col, ltc_cols_with_variation)
	# Add selected covariates if any
	if (!is.null(selected_covariates) && length(selected_covariates) > 0) {
		covariate_terms <- c(covariate_terms, selected_covariates)
	}
	model_data$group <- factor(model_data$strata)
	formula_str <- paste0("treatment ~ ", paste(covariate_terms, collapse = " + "), "+ factor(group)")

	# Fit model
	tryCatch({
		model <- glm(
			formula = as.formula(formula_str),
			data = model_data,
			family = binomial(link = "logit")
		)

		# Extract results for the medication
		coef_summary <- summary(model)$coefficients
		medication_row <- coef_summary[med_col, , drop = FALSE]

		# Calculate OR and CI
		or <- exp(medication_row[1, "Estimate"])
		se <- medication_row[1, "Std. Error"]
		ci_lower <- exp(medication_row[1, "Estimate"] - 1.96 * se)
		ci_upper <- exp(medication_row[1, "Estimate"] + 1.96 * se)
		p_value <- medication_row[1, "Pr(>|z|)"]

		# Count cases and controls with medication
		n_cases <- model_data[get(med_col) == 1 & treatment == 1, .N]
		n_controls <- model_data[get(med_col) == 1 & treatment == 0, .N]

		data.table(
			medication = medication,
			OR = or,
			CI_lower = ci_lower,
			CI_upper = ci_upper,
			p_value = p_value,
			n_cases = n_cases,
			n_controls = n_controls,
			n_ltc_covariates = length(ltc_cols_with_variation),
			convergence = "success"
		)

	}, error = function(e) {
		n_cases_err <- tryCatch(
			model_data[get(med_col) == 1 & treatment == 1, .N],
			error = function(e2) NA_integer_
		)
		n_controls_err <- tryCatch(
			model_data[get(med_col) == 1 & treatment == 0, .N],
			error = function(e2) NA_integer_
		)

		data.table(
			medication = medication,
			OR = NA_real_,
			CI_lower = NA_real_,
			CI_upper = NA_real_,
			p_value = NA_real_,
			n_cases = n_cases_err,
			n_controls = n_controls_err,
			n_ltc_covariates = length(ltc_cols_with_variation),
			convergence = paste("error:", substr(e$message, 1, 50))
		)
	})
}

#' Fit PP interaction models using most recent prescriptions
#'
#' @param selected_ltcs character vector of LTC terms to adjust for
#' @param medications character vector of medications to model
#' @param patient_data data.table with patient info including pp_group
#' @param recent_prescriptions data.table from cases_controls_r with most recent prescriptions
#' @param ltcs data.table with LTCs
#' @return data.table with interaction model results
fit_pp_interaction_models <- function(selected_ltcs, selected_covariates, medications,
																			patient_data, recent_prescriptions, ltcs) {

	# Strip asterisks from LTC names
	selected_ltcs <- gsub("\\*$", "", selected_ltcs)

	# Prepare base model data with recent prescriptions
	model_data <- prepare_pp_interaction_data(
		patient_data = patient_data,
		recent_prescriptions = recent_prescriptions,
		ltcs = ltcs,
		medications = medications,
		selected_ltcs = selected_ltcs
	)

	# Fit model for each medication
	results_list <- lapply(medications, function(med) {
		fit_single_pp_interaction_model(
			model_data = model_data,
			medication = med,
			selected_ltcs = selected_ltcs,
			selected_covariates = selected_covariates
		)
	})

	# Combine results and remove NULL entries
	results <- rbindlist(results_list[!sapply(results_list, is.null)], fill = TRUE)

	return(results)
}

#' Prepare data for PP interaction models
#'
#' @param patient_data data.table with patient characteristics
#' @param recent_prescriptions data.table with most recent prescriptions (cases_controls_r)
#' @param ltcs data.table with LTCs
#' @param medications character vector of medications
#' @param selected_ltcs character vector of LTCs to include
#' @return data.table in wide format ready for modeling
prepare_pp_interaction_data <- function(patient_data, recent_prescriptions, ltcs,
																				medications, selected_ltcs) {

	# Start with patient data - including pp_group
	base_data <- copy(patient_data[, .(patid, treatment, strata, eth_group, imd_quintile, sex, pp, pp_group)])

	# CRITICAL: Ensure pp_group is an UNORDERED factor
	# This prevents R from using polynomial contrasts (L, Q, C)
	base_data[, pp_group := factor(pp_group, ordered = FALSE)]

	# Filter recent_prescriptions to only medications we're modeling
	recent_presc_filtered <- copy(recent_prescriptions[substance %in% medications, .(patid, substance)])

	# Remove duplicates first
	recent_presc_filtered <- unique(recent_presc_filtered)

	# Create wide format with medication indicators
	med_wide <- dcast(
		recent_presc_filtered,
		patid ~ substance,
		fun.aggregate = length,
		value.var = "substance"
	)

	# Convert counts to binary indicators (0/1)
	med_cols <- setdiff(names(med_wide), "patid")
	for (col in med_cols) {
		med_wide[, (col) := as.integer(get(col) > 0)]
	}

	# Rename medication columns with med_ prefix
	setnames(med_wide, med_cols, paste0("med_", make.names(med_cols)))

	# Merge medications with base data
	base_data <- merge(base_data, med_wide, by = "patid", all.x = TRUE)

	# Fill missing medication indicators with 0
	med_indicator_cols <- grep("^med_", names(base_data), value = TRUE)
	for (col in med_indicator_cols) {
		base_data[is.na(get(col)), (col) := 0L]
	}

	# Create LTC indicators - filter to selected LTCs only
	ltc_filtered <- copy(ltcs[term %in% selected_ltcs, .(patid, term)])
	ltc_filtered <- unique(ltc_filtered)

	# Create wide format with LTC indicators
	ltc_wide <- dcast(
		ltc_filtered,
		patid ~ term,
		fun.aggregate = length,
		value.var = "term"
	)

	# Convert counts to binary indicators
	ltc_cols <- setdiff(names(ltc_wide), "patid")
	for (col in ltc_cols) {
		ltc_wide[, (col) := as.integer(get(col) > 0)]
	}

	# Rename LTC columns with ltc_ prefix
	setnames(ltc_wide, ltc_cols, paste0("ltc_", make.names(ltc_cols)))

	# Merge LTCs with base data
	base_data <- merge(base_data, ltc_wide, by = "patid", all.x = TRUE)

	# Fill missing LTC indicators with 0
	ltc_indicator_cols <- grep("^ltc_", names(base_data), value = TRUE)
	for (col in ltc_indicator_cols) {
		base_data[is.na(get(col)), (col) := 0L]
	}

	return(base_data)
}

#' Fit single PP interaction model
#'
#' @param model_data data.table prepared for modeling
#' @param medication character medication name
#' @param selected_ltcs character vector of LTC terms
#' @return data.table with model results
fit_single_pp_interaction_model <- function(model_data, medication, selected_ltcs, selected_covariates) {

	med_col <- paste0("med_", make.names(medication))

	# Check medication exists and has variation
	if (!med_col %in% names(model_data)) {
		return(NULL)
	}

	if (model_data[, uniqueN(get(med_col))] < 2) {
		return(NULL)
	}

	# Get LTC columns with variation
	ltc_cols <- paste0("ltc_", make.names(selected_ltcs))
	ltc_cols <- ltc_cols[ltc_cols %in% names(model_data)]
	ltc_cols_with_variation <- character(0)

	for (col in ltc_cols) {
		if (model_data[, uniqueN(get(col))] > 1) {
			ltc_cols_with_variation <- c(ltc_cols_with_variation, col)
		}
	}

	# Build formula with interaction
	interaction_term <- paste0(med_col, ":pp_group")
	covariate_terms <- c(med_col, "pp_group", ltc_cols_with_variation)
	# Add selected covariates if any
	if (!is.null(selected_covariates) && length(selected_covariates) > 0) {
		covariate_terms <- c(covariate_terms, selected_covariates)
	}
	model_data$group <- factor(model_data$strata)

	formula_str <- paste0("treatment ~ ",
												paste(covariate_terms, collapse = " + "),
												" + ", interaction_term,
												" + factor(group)")

	# Fit model
	tryCatch({
		model <- glm(
			formula = as.formula(formula_str),
			data = model_data,
			family = binomial(link = "logit")
		)

		# Extract coefficients
		coef_summary <- summary(model)$coefficients

		# Get main medication effect
		main_effect <- extract_coefficient(coef_summary, med_col)

		# Get interaction effects for PP groups
		interaction_results <- extract_pp_group_interactions(coef_summary, med_col)

		# Calculate combined effects
		result <- data.table(
			medication = medication,
			pp_level = interaction_results$pp_level,
			main_OR = main_effect$OR,
			main_CI_lower = main_effect$CI_lower,
			main_CI_upper = main_effect$CI_upper,
			main_p = main_effect$p,
			interaction_OR = interaction_results$OR,
			interaction_CI_lower = interaction_results$CI_lower,
			interaction_CI_upper = interaction_results$CI_upper,
			interaction_p = interaction_results$p
		)

		# Calculate combined OR and CIs
		result[, combined_OR := main_OR * interaction_OR]
		result[, combined_CI_lower := main_CI_lower * interaction_CI_lower]
		result[, combined_CI_upper := main_CI_upper * interaction_CI_upper]

		# Create formatted strings for display
		result[, main_OR_formatted := sprintf("%.2f (%.2f-%.2f)", main_OR, main_CI_lower, main_CI_upper)]
		result[, interaction_OR_formatted := sprintf("%.2f (%.2f-%.2f)", interaction_OR, interaction_CI_lower, interaction_CI_upper)]
		result[, combined_OR_formatted := sprintf("%.2f (%.2f-%.2f)", combined_OR, combined_CI_lower, combined_CI_upper)]

		result[, n_ltc_covariates := length(ltc_cols_with_variation)]
		result[, convergence := "success"]

		# Add prevalence information
		result <- add_prevalence_info(result, model_data, med_col)

		return(result)

	}, error = function(e) {
		data.table(
			medication = medication,
			pp_level = NA_character_,
			main_OR = NA_real_,
			main_CI_lower = NA_real_,
			main_CI_upper = NA_real_,
			main_p = NA_real_,
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
			combined_OR = NA_real_,
			combined_CI_lower = NA_real_,
			combined_CI_upper = NA_real_,
			main_OR_formatted = NA_character_,
			interaction_OR_formatted = NA_character_,
			combined_OR_formatted = NA_character_,
			n_ltc_covariates = length(ltc_cols_with_variation),
			convergence = paste("error:", substr(e$message, 1, 50))
		)
	})
}

#' Extract single coefficient from model summary
#'
#' @param coef_summary matrix from summary(model)$coefficients
#' @param term_name character name of coefficient
#' @return list with OR, CI_lower, CI_upper, p
extract_coefficient <- function(coef_summary, term_name) {
	if (!term_name %in% rownames(coef_summary)) {
		return(list(OR = NA_real_, CI_lower = NA_real_, CI_upper = NA_real_, p = NA_real_))
	}

	row <- coef_summary[term_name, , drop = FALSE]
	estimate <- row[1, "Estimate"]
	se <- row[1, "Std. Error"]

	list(
		OR = exp(estimate),
		CI_lower = exp(estimate - 1.96 * se),
		CI_upper = exp(estimate + 1.96 * se),
		p = row[1, "Pr(>|z|)"]
	)
}

#' Extract PP group interaction coefficients
#'
#' @param coef_summary matrix from summary(model)$coefficients
#' @param med_col character medication column name
#' @return data.table with pp_level, OR, CI_lower, CI_upper, p
extract_pp_group_interactions <- function(coef_summary, med_col) {
	# Find all interaction terms matching pattern med_col:pp_group
	interaction_pattern <- paste0("^", med_col, ":pp_group")
	interaction_rows <- grep(interaction_pattern, rownames(coef_summary), value = TRUE)

	# Extract results for each non-reference group
	if (length(interaction_rows) > 0) {
		results <- lapply(interaction_rows, function(term) {
			row <- coef_summary[term, , drop = FALSE]
			estimate <- row[1, "Estimate"]
			se <- row[1, "Std. Error"]

			# Extract PP level from term name
			# e.g., "med_Opioid.Analgesics:pp_group[2,5)" -> "[2,5)"
			pp_level <- gsub(paste0(".*:pp_group"), "", term)

			data.table(
				pp_level = pp_level,
				OR = exp(estimate),
				CI_lower = exp(estimate - 1.96 * se),
				CI_upper = exp(estimate + 1.96 * se),
				p = row[1, "Pr(>|z|)"]
			)
		})

		results_dt <- rbindlist(results)
	} else {
		# If no interaction terms found, might be no variation
		results_dt <- data.table(
			pp_level = "no_variation",
			OR = NA_real_,
			CI_lower = NA_real_,
			CI_upper = NA_real_,
			p = NA_real_
		)
	}

	# Add reference group (first level of pp_group, typically [0,2))
	# The reference group has OR = 1.0 by definition
	reference_row <- data.table(
		pp_level = "reference",
		OR = 1.0,
		CI_lower = 1.0,
		CI_upper = 1.0,
		p = NA_real_
	)

	# Combine reference and other levels
	final_results <- rbind(reference_row, results_dt)

	return(final_results)
}

#' Add prevalence information to results
#'
#' @param result data.table with model results
#' @param model_data data.table used for modeling
#' @param med_col character medication column name
#' @return data.table with added prevalence columns
add_prevalence_info <- function(result, model_data, med_col) {

	total_cases <- model_data[treatment == 1, .N]
	total_controls <- model_data[treatment == 0, .N]

	n_cases <- model_data[get(med_col) == 1 & treatment == 1, .N]
	n_controls <- model_data[get(med_col) == 1 & treatment == 0, .N]

	result[, `:=`(
		pct_cases = round(100 * n_cases / total_cases, 2),
		pct_controls = round(100 * n_controls / total_controls, 2),
		n_cases = n_cases,
		n_controls = n_controls
	)]

	return(result)
}

#' Filter medication pairs by co-prescription prevalence
#'
#' @param med_pairs list of medication pairs
#' @param min_coprescription_prev minimum co-prescription prevalence threshold
#' @param patient_data data.table with patid and treatment
#' @param prescriptions data.table with patid and substance
#' @return list of medication pairs meeting threshold
filter_pairs_by_coprescription <- function(med_pairs, min_coprescription_prev,
																					 patient_data, prescriptions) {

	# Calculate total cases and controls
	total_cases <- patient_data[treatment == 1, .N]
	total_controls <- patient_data[treatment == 0, .N]

	# Filter pairs
	filtered_pairs <- Filter(function(pair) {
		med1 <- gsub("\\*$", "", pair[1])
		med2 <- gsub("\\*$", "", pair[2])

		# Get patients with med1
		patids_med1 <- prescriptions[substance == med1, unique(patid)]
		# Get patients with med2
		patids_med2 <- prescriptions[substance == med2, unique(patid)]
		# Get patients with both
		patids_both <- bit64::as.integer64(intersect(as.character(patids_med1), as.character(patids_med2)))

		# Calculate prevalence in cases and controls
		n_cases_both <- patient_data[patid %in% patids_both & treatment == 1, .N]
		n_controls_both <- patient_data[patid %in% patids_both & treatment == 0, .N]

		pct_cases_both <- 100 * n_cases_both / total_cases
		pct_controls_both <- 100 * n_controls_both / total_controls

		# Keep pair if EITHER cases or controls meet threshold
		return(pct_cases_both >= min_coprescription_prev | pct_controls_both >= min_coprescription_prev)

	}, med_pairs)

	return(filtered_pairs)
}

# R/utils_cca_logreg.R - Add these functions to the existing file

#' Run recent prescription × background medication interaction models
#'
#' @param recent_meds character vector of recent prescription names
#' @param background_meds character vector of background medication names
#' @param selected_ltcs character vector of LTC terms to include as covariates
#' @param patient_data data.table with patid, treatment, and demographic info
#' @param recent_prescriptions data.table from cases_controls_r with most recent prescriptions
#' @param prescriptions data.table with background prescriptions (patid, substance, group)
#' @param ltcs data.table with patid, term, group
#' @return data.table with interaction model results
run_recent_background_interaction_models <- function(recent_meds, background_meds,
																										 selected_ltcs, selected_covariates, patient_data,
																										 recent_prescriptions, prescriptions, ltcs) {

	# Strip asterisks from names
	recent_meds <- gsub("\\*$", "", recent_meds)
	background_meds <- gsub("\\*$", "", background_meds)
	selected_ltcs <- gsub("\\*$", "", selected_ltcs)

	# Prepare the dataset with both recent and background medications
	model_data <- prepare_recent_background_data(
		patient_data = patient_data,
		recent_prescriptions = recent_prescriptions,
		prescriptions = prescriptions,
		ltcs = ltcs,
		recent_meds = recent_meds,
		background_meds = background_meds,
		selected_ltcs = selected_ltcs
	)

	# Create all pairs of recent × background
	all_pairs <- expand.grid(
		recent = recent_meds,
		background = background_meds,
		stringsAsFactors = FALSE
	)

	# Fit interaction model for each pair
	results_list <- lapply(seq_len(nrow(all_pairs)), function(i) {
		fit_recent_background_interaction_model(
			model_data = model_data,
			recent_med = all_pairs$recent[i],
			background_med = all_pairs$background[i],
			selected_ltcs = selected_ltcs,
			selected_covariates = selected_covariates
		)
	})

	# Combine results and remove NULL entries
	results <- rbindlist(results_list[!sapply(results_list, is.null)], fill = TRUE)

	return(results)
}

#' Prepare data for recent × background interaction models
#'
#' @param patient_data data.table with patient characteristics
#' @param recent_prescriptions data.table with most recent prescriptions (cases_controls_r)
#' @param prescriptions data.table with background prescriptions
#' @param ltcs data.table with LTCs
#' @param recent_meds character vector of recent medications
#' @param background_meds character vector of background medications
#' @param selected_ltcs character vector of LTCs to include
#' @return data.table in wide format ready for modeling
prepare_recent_background_data <- function(patient_data, recent_prescriptions, prescriptions,
																					 ltcs, recent_meds, background_meds, selected_ltcs) {

	# Start with patient data
	base_data <- copy(patient_data[, .(patid, treatment, strata, eth_group, imd_quintile, sex, pp_group, mltc_group)])

	# Create indicators for RECENT prescriptions
	recent_filtered <- copy(recent_prescriptions[substance %in% recent_meds, .(patid, substance)])
	recent_filtered <- unique(recent_filtered)

	if (nrow(recent_filtered) > 0) {
		recent_wide <- dcast(
			recent_filtered,
			patid ~ substance,
			fun.aggregate = length,
			value.var = "substance"
		)

		# Convert to binary and rename with recent_ prefix
		recent_cols <- setdiff(names(recent_wide), "patid")
		for (col in recent_cols) {
			recent_wide[, (col) := as.integer(get(col) > 0)]
		}
		setnames(recent_wide, recent_cols, paste0("recent_", make.names(recent_cols)))

		# Merge with base data
		base_data <- merge(base_data, recent_wide, by = "patid", all.x = TRUE)

		# Fill missing with 0
		recent_indicator_cols <- grep("^recent_", names(base_data), value = TRUE)
		for (col in recent_indicator_cols) {
			base_data[is.na(get(col)), (col) := 0L]
		}
	}

	# Create indicators for BACKGROUND medications
	background_filtered <- copy(prescriptions[substance %in% background_meds, .(patid, substance)])
	background_filtered <- unique(background_filtered)

	if (nrow(background_filtered) > 0) {
		background_wide <- dcast(
			background_filtered,
			patid ~ substance,
			fun.aggregate = length,
			value.var = "substance"
		)

		# Convert to binary and rename with background_ prefix
		background_cols <- setdiff(names(background_wide), "patid")
		for (col in background_cols) {
			background_wide[, (col) := as.integer(get(col) > 0)]
		}
		setnames(background_wide, background_cols, paste0("background_", make.names(background_cols)))

		# Merge with base data
		base_data <- merge(base_data, background_wide, by = "patid", all.x = TRUE)

		# Fill missing with 0
		background_indicator_cols <- grep("^background_", names(base_data), value = TRUE)
		for (col in background_indicator_cols) {
			base_data[is.na(get(col)), (col) := 0L]
		}
	}

	# Create LTC indicators
	ltc_filtered <- copy(ltcs[term %in% selected_ltcs, .(patid, term)])
	ltc_filtered <- unique(ltc_filtered)

	if (nrow(ltc_filtered) > 0) {
		ltc_wide <- dcast(
			ltc_filtered,
			patid ~ term,
			fun.aggregate = length,
			value.var = "term"
		)

		ltc_cols <- setdiff(names(ltc_wide), "patid")
		for (col in ltc_cols) {
			ltc_wide[, (col) := as.integer(get(col) > 0)]
		}
		setnames(ltc_wide, ltc_cols, paste0("ltc_", make.names(ltc_cols)))

		base_data <- merge(base_data, ltc_wide, by = "patid", all.x = TRUE)

		ltc_indicator_cols <- grep("^ltc_", names(base_data), value = TRUE)
		for (col in ltc_indicator_cols) {
			base_data[is.na(get(col)), (col) := 0L]
		}
	}

	return(base_data)
}

#' Fit interaction model for recent × background medication pair
#'
#' @param model_data data.table prepared for modeling
#' @param recent_med character string of recent medication name
#' @param background_med character string of background medication name
#' @param selected_ltcs character vector of LTC terms
#' @return data.table with interaction model results
fit_recent_background_interaction_model <- function(model_data, recent_med,
																										background_med, selected_ltcs, selected_covariates) {

	# Get column names
	recent_col <- paste0("recent_", make.names(recent_med))
	background_col <- paste0("background_", make.names(background_med))
	ltc_cols <- paste0("ltc_", make.names(selected_ltcs))

	# Check both medications exist
	if (!recent_col %in% names(model_data) || !background_col %in% names(model_data)) {
		return(NULL)
	}

	# Check for variation in both medications
	if (model_data[, uniqueN(get(recent_col))] < 2 ||
			model_data[, uniqueN(get(background_col))] < 2) {
		return(NULL)
	}

	# Calculate prevalence of both medications together
	n_cases_both <- model_data[get(recent_col) == 1 & get(background_col) == 1 & treatment == 1, .N]
	n_controls_both <- model_data[get(recent_col) == 1 & get(background_col) == 1 & treatment == 0, .N]
	total_cases <- model_data[treatment == 1, .N]
	total_controls <- model_data[treatment == 0, .N]

	pct_cases_both <- tryCatch(
		(n_cases_both / total_cases) * 100,
		error = function(e) NA_real_
	)
	pct_controls_both <- tryCatch(
		(n_controls_both / total_controls) * 100,
		error = function(e) NA_real_
	)

	# Keep only LTC columns that exist and have variation
	existing_ltc_cols <- ltc_cols[ltc_cols %in% names(model_data)]
	ltc_cols_with_variation <- character(0)
	for (col in existing_ltc_cols) {
		if (model_data[, uniqueN(get(col))] > 1) {
			ltc_cols_with_variation <- c(ltc_cols_with_variation, col)
		}
	}

	# Build formula with interaction
	covariate_terms <- c(recent_col, background_col, ltc_cols_with_variation)
	# Add selected covariates if any
	if (!is.null(selected_covariates) && length(selected_covariates) > 0) {
		covariate_terms <- c(covariate_terms, selected_covariates)
	}
	interaction_term <- paste0(recent_col, ":", background_col)
	model_data$group <- factor(model_data$strata)
	formula_str <- paste0("treatment ~ ",
												paste(covariate_terms, collapse = " + "),
												" + ", interaction_term,
												" + factor(group)")
	# Fit model
	tryCatch({
		model <- glm(
			formula = as.formula(formula_str),
			data = model_data,
			family = binomial(link = "logit")
		)

		coef_summary <- summary(model)$coefficients

		message("Coefficients: ", paste(rownames(coef_summary), collapse = ", "))

		# Extract RECENT medication main effect
		if (recent_col %in% rownames(coef_summary)) {
			recent_row <- coef_summary[recent_col, , drop = FALSE]
			recent_or <- exp(recent_row[1, "Estimate"])
			recent_se <- recent_row[1, "Std. Error"]
			recent_ci_lower <- exp(recent_row[1, "Estimate"] - 1.96 * recent_se)
			recent_ci_upper <- exp(recent_row[1, "Estimate"] + 1.96 * recent_se)
			recent_p <- recent_row[1, "Pr(>|z|)"]
		} else {
			recent_or <- NA_real_
			recent_ci_lower <- NA_real_
			recent_ci_upper <- NA_real_
			recent_p <- NA_real_
		}

		# Extract BACKGROUND medication main effect
		if (background_col %in% rownames(coef_summary)) {
			background_row <- coef_summary[background_col, , drop = FALSE]
			background_or <- exp(background_row[1, "Estimate"])
			background_se <- background_row[1, "Std. Error"]
			background_ci_lower <- exp(background_row[1, "Estimate"] - 1.96 * background_se)
			background_ci_upper <- exp(background_row[1, "Estimate"] + 1.96 * background_se)
			background_p <- background_row[1, "Pr(>|z|)"]
		} else {
			background_or <- NA_real_
			background_ci_lower <- NA_real_
			background_ci_upper <- NA_real_
			background_p <- NA_real_
		}

		# Extract INTERACTION effect
		interaction_name <- paste0(recent_col, ":", background_col)
		if (interaction_name %in% rownames(coef_summary)) {
			interaction_row <- coef_summary[interaction_name, , drop = FALSE]
			interaction_or <- exp(interaction_row[1, "Estimate"])
			interaction_se <- interaction_row[1, "Std. Error"]
			interaction_ci_lower <- exp(interaction_row[1, "Estimate"] - 1.96 * interaction_se)
			interaction_ci_upper <- exp(interaction_row[1, "Estimate"] + 1.96 * interaction_se)
			interaction_p <- interaction_row[1, "Pr(>|z|)"]
		} else {
			interaction_or <- NA_real_
			interaction_ci_lower <- NA_real_
			interaction_ci_upper <- NA_real_
			interaction_p <- NA_real_
		}

		# Calculate COMBINED effect (product of main effects and interaction)
		combined_or <- recent_or * background_or * interaction_or
		combined_ci_lower <- recent_ci_lower * background_ci_lower * interaction_ci_lower
		combined_ci_upper <- recent_ci_upper * background_ci_upper * interaction_ci_upper

		# Return results
		data.table(
			recent_med = recent_med,
			background_med = background_med,
			recent_OR = recent_or,
			recent_CI_lower = recent_ci_lower,
			recent_CI_upper = recent_ci_upper,
			recent_p = recent_p,
			background_OR = background_or,
			background_CI_lower = background_ci_lower,
			background_CI_upper = background_ci_upper,
			background_p = background_p,
			interaction_OR = interaction_or,
			interaction_CI_lower = interaction_ci_lower,
			interaction_CI_upper = interaction_ci_upper,
			interaction_p = interaction_p,
			combined_OR = combined_or,
			combined_CI_lower = combined_ci_lower,
			combined_CI_upper = combined_ci_upper,
			pct_cases_both = pct_cases_both,
			pct_controls_both = pct_controls_both,
			n_cases_both = n_cases_both,
			n_controls_both = n_controls_both,
			n_ltc_covariates = length(ltc_cols_with_variation),
			convergence = "success"
		)

	}, error = function(e) {
		data.table(
			recent_med = recent_med,
			background_med = background_med,
			recent_OR = NA_real_,
			recent_CI_lower = NA_real_,
			recent_CI_upper = NA_real_,
			recent_p = NA_real_,
			background_OR = NA_real_,
			background_CI_lower = NA_real_,
			background_CI_upper = NA_real_,
			background_p = NA_real_,
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
			combined_OR = NA_real_,
			combined_CI_lower = NA_real_,
			combined_CI_upper = NA_real_,
			pct_cases_both = pct_cases_both,
			pct_controls_both = pct_controls_both,
			n_cases_both = n_cases_both,
			n_controls_both = n_controls_both,
			n_ltc_covariates = length(selected_ltcs),
			convergence = paste("error:", substr(e$message, 1, 50))
		)
	})
}