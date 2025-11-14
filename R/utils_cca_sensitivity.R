#' Utility Functions for CCA Sensitivity Analysis
#'
#' Functions to perform post-hoc sensitivity analyses on matched case-control data

#' Calculate prescription ORs for case-control data
#'
#' @param presc_with_group data.table with columns: patid, substance, group (case/control)
#' @param min_prevalence minimum prevalence threshold (default 0.02)
#' @return data.table with OR, CI, and prevalence statistics
#'
#' @details This is a wrapper around the existing create_prevalence_ratio_table
#' function from utils_cca_outputs.R to maintain consistency across the tool.
calculate_prescription_ors <- function(presc_with_group, min_prevalence = 0.02) {

	# Calculate frequency stats
	freq_data <- presc_with_group[, .N, by = .(substance, group)]

	# Get group totals
	group_totals <- presc_with_group[, .(total = uniqueN(patid)), by = group]

	# Add percentages
	freq_data[group_totals, pct := round(N / total * 100, 2), on = "group"]

	# Use existing create_prevalence_ratio_table function
	# This reuses the vectorized OR calculation already in the project
	result <- create_prevalence_ratio_table(
		freq_data = freq_data,
		item_col = "substance",
		min_pct = min_prevalence * 100,  # Convert proportion to percentage
		data_with_group = presc_with_group,
		test_method = "chisq",
		p_adjust_method = "BH"
	)

	if (is.null(result) || nrow(result) == 0) {
		return(data.table())
	}

	# Standardize column names for consistency
	# create_prevalence_ratio_table returns: substance, case, control, OR, OR_CI_lower, OR_CI_upper, p_value, p_adj
	expected_cols <- c("substance", "case", "control", "OR", "OR_CI_lower", "OR_CI_upper", "p_value", "p_adj")

	# Rename CI columns to match expected names
	if ("OR_CI_lower" %in% names(result)) {
		setnames(result, "OR_CI_lower", "CI_lower")
	}
	if ("OR_CI_upper" %in% names(result)) {
		setnames(result, "OR_CI_upper", "CI_upper")
	}

	# Rename case/control prevalence columns
	if ("case" %in% names(result) && "control" %in% names(result)) {
		setnames(result,
						 old = c("case", "control"),
						 new = c("case_prev", "control_prev"))
	}

	# Calculate sample sizes if not already present
	n_case_total <- uniqueN(presc_with_group[group == "case", patid])
	n_control_total <- uniqueN(presc_with_group[group == "control", patid])

	result[, n_cases := round(case_prev * n_case_total / 100)]
	result[, n_controls := round(control_prev * n_control_total / 100)]

	setorder(result, -OR)

	return(result)
}

#' Create LTC indicator matrix for patients
#'
#' @param ltcs data.table with patid and term columns
#' @param patids vector of patient IDs to include
#' @return data.table with patid and binary indicators for each LTC
create_ltc_indicators <- function(ltcs, patids) {

	# Filter to relevant patients
	ltcs_subset <- ltcs[patid %in% patids]

	# Create binary indicators
	ltc_matrix <- dcast(
		ltcs_subset,
		patid ~ term,
		fun.aggregate = length,
		value.var = "term"
	)

	# Convert to 0/1
	ltc_cols <- names(ltc_matrix)[-1]
	ltc_matrix[, (ltc_cols) := lapply(.SD, function(x) as.integer(x > 0)), .SDcols = ltc_cols]

	return(ltc_matrix)
}

#' Run stratified analysis by LTC presence
#'
#' @param prescriptions data.table of prescriptions
#' @param ltcs data.table of LTCs
#' @param patient_data data.table with treatment variable
#' @param stratify_by_ltcs character vector of LTC terms to stratify by
#' @param min_prevalence minimum drug prevalence threshold
#' @return data.table with stratified OR results
run_stratified_analysis <- function(prescriptions, ltcs, patient_data,
																		stratify_by_ltcs, min_prevalence = 0.02) {

	# Create LTC indicators
	all_patids <- unique(patient_data$patid)
	ltc_indicators <- create_ltc_indicators(ltcs, all_patids)

	# Merge prescriptions with patient data and LTC indicators
	presc_full <- merge(prescriptions, patient_data[, .(patid, treatment)], by = "patid")
	presc_full <- merge(presc_full, ltc_indicators, by = "patid", all.x = TRUE)
	presc_full[, group := ifelse(treatment == 1, "case", "control")]

	# For each LTC, calculate ORs for with/without
	all_results <- list()

	for (ltc in stratify_by_ltcs) {

		# Check if LTC column exists
		if (!ltc %in% names(presc_full)) {
			warning(paste("LTC", ltc, "not found in data"))
			next
		}

		# Replace NA with 0 for the LTC indicator
		if (ltc %in% names(presc_full)) {
			presc_full[is.na(get(ltc)), (ltc) := 0]
		}

		# Stratum: WITH the LTC
		with_ltc <- presc_full[get(ltc) == 1]
		if (uniqueN(with_ltc$patid) > 50) {  # Need minimum sample size
			or_with <- calculate_prescription_ors(with_ltc, min_prevalence)
			if (nrow(or_with) > 0) {
				or_with[, stratum := paste0("With ", ltc)]
				all_results[[paste0(ltc, "_with")]] <- or_with
			}
		}

		# Stratum: WITHOUT the LTC
		without_ltc <- presc_full[get(ltc) == 0]
		if (uniqueN(without_ltc$patid) > 50) {
			or_without <- calculate_prescription_ors(without_ltc, min_prevalence)
			if (nrow(or_without) > 0) {
				or_without[, stratum := paste0("Without ", ltc)]
				all_results[[paste0(ltc, "_without")]] <- or_without
			}
		}
	}

	# Combine all results
	if (length(all_results) == 0) {
		stop("No valid strata found. Try different LTCs or lower sample size requirements.")
	}

	combined <- rbindlist(all_results, use.names = TRUE, fill = TRUE)

	# Add interpretation
	combined[, interpretation := fcase(
		abs(OR - 1) < 0.3, "Minimal effect",
		OR > 1.5, "Elevated risk",
		OR < 0.67, "Reduced risk",
		default = "Moderate effect"
	)]

	# Order by substance then stratum
	setorder(combined, substance, stratum)

	return(combined)
}

#' Run exclusion sensitivity analysis
#'
#' @param prescriptions data.table of prescriptions
#' @param ltcs data.table of LTCs
#' @param patient_data data.table with treatment variable
#' @param exclude_ltcs character vector of LTC terms to exclude
#' @param main_or_results data.table of main analysis ORs
#' @param min_prevalence minimum drug prevalence threshold
#' @return data.table comparing main vs exclusion ORs
run_exclusion_analysis <- function(prescriptions, ltcs, patient_data,
																	 exclude_ltcs, main_or_results, min_prevalence = 0.02) {

	# Find patients with any of the excluded LTCs
	patients_to_exclude <- unique(ltcs[term %in% exclude_ltcs, patid])

	message(sprintf("Excluding %d patients with selected LTCs", length(patients_to_exclude)))

	# Filter prescriptions and patient data
	presc_filtered <- prescriptions[!patid %in% patients_to_exclude]
	patient_data_filtered <- patient_data[!patid %in% patients_to_exclude]

	# Calculate ORs on filtered data
	presc_with_group <- merge(
		presc_filtered,
		patient_data_filtered[, .(patid, treatment)],
		by = "patid"
	)
	presc_with_group[, group := ifelse(treatment == 1, "case", "control")]

	or_excluded <- calculate_prescription_ors(presc_with_group, min_prevalence)

	if (nrow(or_excluded) == 0) {
		stop("No substances met the minimum prevalence threshold after exclusion")
	}

	# Merge with main results for comparison
	# Ensure we have the right column names
	comparison <- merge(
		main_or_results[, .(substance,
												main_OR = OR,
												main_CI_lower = CI_lower,
												main_CI_upper = CI_upper,
												main_n_cases = n_cases,
												main_n_controls = n_controls)],
		or_excluded[, .(substance,
										excluded_OR = OR,
										excluded_CI_lower = CI_lower,
										excluded_CI_upper = CI_upper,
										excluded_n_cases = n_cases,
										excluded_n_controls = n_controls)],
		by = "substance",
		all = FALSE  # Only substances in both analyses
	)

	if (nrow(comparison) == 0) {
		stop("No overlapping substances between main and exclusion analyses. Try lowering the prevalence threshold.")
	}

	# Calculate change metrics
	comparison[, `:=`(
		pct_change = ((excluded_OR - main_OR) / main_OR) * 100,
		n_excluded = main_n_cases + main_n_controls - excluded_n_cases - excluded_n_controls
	)]

	# Format CIs as strings
	comparison[, main_CI := sprintf("(%.2f-%.2f)", main_CI_lower, main_CI_upper)]
	comparison[, excluded_CI := sprintf("(%.2f-%.2f)", excluded_CI_lower, excluded_CI_upper)]

	# Add robustness assessment
	comparison[, robustness := fcase(
		abs(pct_change) < 10, "✓ Robust",
		abs(pct_change) < 25, "⚠ Moderate change",
		default = "✗ Substantial change"
	)]

	# Select final columns
	result <- comparison[, .(
		substance,
		main_OR,
		excluded_OR,
		pct_change,
		main_CI,
		excluded_CI,
		n_excluded,
		robustness
	)]

	# Sort by absolute value of pct_change (create temp column for sorting)
	result[, abs_pct_change := abs(pct_change)]
	setorder(result, -abs_pct_change)
	result[, abs_pct_change := NULL]  # Remove temp column

	return(result)
}

#' Create forest plot comparing ORs across analyses (Vega-Lite version)
#'
#' @param plot_data data.table with columns: analysis, OR, CI_lower, CI_upper
#' @param drug_name character string for plot title
#' @return vegaspec object
#'
#' @details Creates a forest plot using Vega-Lite to match the tool's existing
#' visualization style. Shows OR point estimates with confidence intervals for
#' each analysis scenario.
create_forest_plot <- function(plot_data, drug_name) {

	# Filter out invalid values
	plot_data <- plot_data[!is.na(OR) & is.finite(OR) &
												 	!is.na(CI_lower) & is.finite(CI_lower) &
												 	!is.na(CI_upper) & is.finite(CI_upper)]

	if (nrow(plot_data) == 0) {
		return(NULL)
	}

	# Get analysis order (reverse for bottom-to-top display)
	analysis_order <- rev(plot_data$analysis)

	# Calculate x-axis range
	min_ci <- min(plot_data$CI_lower, na.rm = TRUE)
	max_ci <- max(plot_data$CI_upper, na.rm = TRUE)

	# Use log scale if range is wide
	use_log_scale <- (max_ci / min_ci) > 5

	# Create Vega-Lite spec matching your existing style
	spec <- list(
		`$schema` = vegawidget::vega_schema(),
		data = list(values = plot_data),
		width = 600,
		height = max(200, 60 * nrow(plot_data)),
		title = list(
			text = paste("Sensitivity Analysis:", drug_name),
			subtitle = "Odds Ratio estimates across different analyses",
			fontSize = 16,
			fontWeight = "bold"
		),
		encoding = list(
			y = list(
				field = "analysis",
				type = "nominal",
				title = NULL,
				sort = analysis_order,
				axis = list(
					labelFontSize = 12,
					labelLimit = 300
				)
			)
		),
		layer = list(
			# Reference line at OR = 1
			list(
				mark = list(
					type = "rule",
					strokeDash = c(4, 4),
					color = "gray",
					strokeWidth = 1.5
				),
				encoding = list(
					x = list(datum = 1),
					y = list()
				)
			),
			# Confidence interval lines
			list(
				mark = list(
					type = "rule",
					size = 3,
					color = "#3498db"
				),
				encoding = list(
					x = list(
						field = "CI_lower",
						type = "quantitative",
						title = "Odds Ratio",
						scale = if (use_log_scale) {
							list(
								type = "log",
								domain = c(max(0.1, min_ci * 0.9), max_ci * 1.1)
							)
						} else {
							list(
								domain = c(0, max_ci * 1.1)
							)
						},
						axis = list(
							grid = TRUE,
							tickCount = 8
						)
					),
					x2 = list(
						field = "CI_upper"
					)
				)
			),
			# OR point estimates
			list(
				mark = list(
					type = "point",
					filled = TRUE,
					size = 100,
					color = "#3498db"
				),
				encoding = list(
					x = list(
						field = "OR",
						type = "quantitative"
					),
					tooltip = list(
						list(field = "analysis", type = "nominal", title = "Analysis"),
						list(field = "OR", type = "quantitative", title = "OR", format = ".2f"),
						list(field = "CI_lower", type = "quantitative", title = "95% CI Lower", format = ".2f"),
						list(field = "CI_upper", type = "quantitative", title = "95% CI Upper", format = ".2f")
					)
				)
			)
		),
		config = list(
			view = list(stroke = NULL),
			font = "Lato, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"
		)
	)

	return(vegawidget::as_vegaspec(spec))
}