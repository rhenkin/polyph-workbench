# R/utils_cca_clogit.R

#' Run logistic regression models for multiple medications
#'
#' @param medications character vector of medication names
#' @param selected_ltcs character vector of LTC terms to include as covariates
#' @param patient_data data.table with patid, treatment, and demographic info
#' @param prescriptions data.table with patid, substance, group
#' @param ltcs data.table with patid, term, group
#' @return data.table with model results
run_conditional_logistic_models <- function(medications, selected_ltcs,
																						patient_data, prescriptions, ltcs) {

	# Strip asterisks from names (added by significance testing)
	medications <- gsub("\\*$", "", medications)
	selected_ltcs <- gsub("\\*$", "", selected_ltcs)

	# Prepare the base dataset
	model_data <- prepare_clogit_data(
		patient_data = patient_data,
		prescriptions = prescriptions,
		ltcs = ltcs,
		medications = medications,
		selected_ltcs = selected_ltcs
	)

	# Fit models for each medication
	results_list <- lapply(medications, function(med) {
		fit_single_clogit_model(
			model_data = model_data,
			medication = med,
			selected_ltcs = selected_ltcs
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
run_interaction_models <- function(med_pairs, selected_ltcs,
																	 patient_data, prescriptions, ltcs) {

	# Strip asterisks from names
	selected_ltcs <- gsub("\\*$", "", selected_ltcs)

	# Get all unique medications from pairs
	all_meds <- unique(unlist(med_pairs))
	all_meds <- gsub("\\*$", "", all_meds)

	# Prepare the base dataset with all medications
	model_data <- prepare_clogit_data(
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
			selected_ltcs = selected_ltcs
		)
	})

	# Combine results
	results <- rbindlist(results_list, fill = TRUE)

	return(results)
}

#' Fit interaction model for a medication pair
#'
#' @param model_data data.table prepared for modeling
#' @param med1 character string of first medication name
#' @param med2 character string of second medication name
#' @param selected_ltcs character vector of LTC terms
#' @return data.table with interaction model results (single row)
fit_interaction_model <- function(model_data, med1, med2, selected_ltcs) {

	# Get column names
	med1_col <- paste0("med_", make.names(med1))
	med2_col <- paste0("med_", make.names(med2))
	ltc_cols <- paste0("ltc_", make.names(selected_ltcs))

	# Check if medication columns exist
	if (!med1_col %in% names(model_data) || !med2_col %in% names(model_data)) {
		return(data.table(
			med1 = med1,
			med2 = med2,
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
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
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
			pct_cases_both = pct_cases_both,
			pct_controls_both = pct_controls_both,
			n_cases_both = n_cases_both,
			n_controls_both = n_controls_both,
			n_ltc_covariates = length(selected_ltcs),
			convergence = "no_variation_in_medications"
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

	# Build formula with interaction: treatment ~ med1 + med2 + med1:med2 + ltcs + strata
	interaction_term <- paste0(med1_col, ":", med2_col)
	covariate_terms <- c(med1_col, med2_col, ltc_cols_with_variation)
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

		# Extract results for the INTERACTION term only
		coef_summary <- summary(model)$coefficients

		# The interaction term in the coefficient table
		interaction_row_name <- paste0(med1_col, ":", med2_col)

		if (!interaction_row_name %in% rownames(coef_summary)) {
			# Try alternative format
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
				interaction_OR = NA_real_,
				interaction_CI_lower = NA_real_,
				interaction_CI_upper = NA_real_,
				interaction_p = NA_real_,
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
		or <- exp(interaction_row[1, "Estimate"])
		se <- interaction_row[1, "Std. Error"]
		ci_lower <- exp(interaction_row[1, "Estimate"] - 1.96 * se)
		ci_upper <- exp(interaction_row[1, "Estimate"] + 1.96 * se)
		p_value <- interaction_row[1, "Pr(>|z|)"]

		# Count cases and controls with BOTH medications and calculate percentages
		n_cases_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 1, .N]
		n_controls_both <- model_data[get(med1_col) == 1 & get(med2_col) == 1 & treatment == 0, .N]

		total_cases <- model_data[treatment == 1, .N]
		total_controls <- model_data[treatment == 0, .N]

		pct_cases_both <- round(100 * n_cases_both / total_cases, 2)
		pct_controls_both <- round(100 * n_controls_both / total_controls, 2)

		data.table(
			med1 = med1,
			med2 = med2,
			interaction_OR = or,
			interaction_CI_lower = ci_lower,
			interaction_CI_upper = ci_upper,
			interaction_p = p_value,
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
			interaction_OR = NA_real_,
			interaction_CI_lower = NA_real_,
			interaction_CI_upper = NA_real_,
			interaction_p = NA_real_,
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
prepare_clogit_data <- function(patient_data, prescriptions, ltcs,
																medications, selected_ltcs) {

	# Start with patient data
	base_data <- copy(patient_data[, .(patid, treatment, sex, pp_group, mltc_group, strata)])

	# Create binary indicators for medications
	med_indicators <- create_medication_indicators(prescriptions, medications)

	# Create binary indicators for LTCs
	ltc_indicators <- create_ltc_indicators_clogit(ltcs, selected_ltcs)

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
create_ltc_indicators_clogit <- function(ltcs, selected_ltcs) {

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
fit_single_clogit_model <- function(model_data, medication, selected_ltcs) {

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