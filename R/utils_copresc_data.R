#' Utility functions for co-prescription analysis calculations
#' Pairwise odds ratio calculations for case-control comparisons
#'
#' @note Depends on calc_all_ors_vectorized() from utils_stats.R

#' Calculate pairwise co-prescription ORs for cases and controls separately
#'
#' @param prescriptions data.table with columns: patid, substance, group
#' @param min_prevalence Minimum prevalence threshold (as proportion)
#' @param min_coprescription Minimum co-prescription threshold (as proportion)
#' @return data.table with pairwise ORs for cases and controls
calculate_coprescription_ors <- function(prescriptions,
																				 min_prevalence = 0.02,
																				 min_coprescription = 0.01) {

	# Separate cases and controls
	case_presc <- prescriptions[group == "case"]
	control_presc <- prescriptions[group == "control"]

	# Calculate ORs for each group
	case_ors <- calculate_group_pairwise_ors(
		case_presc,
		min_prevalence,
		min_coprescription
	)

	control_ors <- calculate_group_pairwise_ors(
		control_presc,
		min_prevalence,
		min_coprescription
	)

	# Merge results
	setnames(case_ors,
					 c("or", "ci_lower", "ci_upper", "p_value"),
					 c("case_or", "case_ci_lower", "case_ci_upper", "case_p"))

	setnames(control_ors,
					 c("or", "ci_lower", "ci_upper", "p_value"),
					 c("control_or", "control_ci_lower", "control_ci_upper", "control_p"))

	# Merge on drug pairs
	merged <- merge(
		case_ors[, .(drug1, drug2, case_or, case_ci_lower, case_ci_upper, case_p)],
		control_ors[, .(drug1, drug2, control_or, control_ci_lower, control_ci_upper, control_p)],
		by = c("drug1", "drug2"),
		all = TRUE
	)

	# Calculate OR difference
	merged[, or_diff := case_or - control_or]

	# Test for significant difference using Z-test
	# SE of difference = sqrt(SE1^2 + SE2^2)
	# SE of log(OR) from CI: SE = (log(upper) - log(lower)) / (2 * 1.96)
	# merged[, case_se := (log(case_ci_upper) - log(case_ci_lower)) / (2 * 1.96)]
	# merged[, control_se := (log(control_ci_upper) - log(control_ci_lower)) / (2 * 1.96)]
	# merged[, diff_se := sqrt(case_se^2 + control_se^2)]
	# merged[, z_stat := (log(case_or) - log(control_or)) / diff_se]
	# merged[, p_diff := 2 * pnorm(-base::abs(z_stat))]
	# merged[, p_diff_adjusted := p.adjust(p_diff, method = "fdr")]
	# merged[, significant := p_diff_adjusted < 0.05]

	# Clean up temporary columns
	# merged[, c("case_se", "control_se", "diff_se", "z_stat") := NULL]

	# Sort by absolute OR difference (using order, not setorder with abs)
	# merged[, abs_or_diff := base::abs(or_diff)]
	# setorder(merged, -abs_or_diff)
	# merged[, abs_or_diff := NULL]

	return(merged)
}


#' Calculate pairwise ORs for a single group (cases or controls)
#'
#' @param presc_data data.table with patid and substance columns
#' @param min_prevalence Minimum prevalence threshold
#' @param min_coprescription Minimum co-prescription threshold
#' @return data.table with drug pairs and their ORs
calculate_group_pairwise_ors <- function(presc_data,
																				 min_prevalence = 0.02,
																				 min_coprescription = 0.01) {

	# Get unique patients and substances
	n_patients <- uniqueN(presc_data$patid)

	# Filter by prevalence
	# drug_counts <- presc_data[, .(n = uniqueN(patid)), by = substance]
	drug_counts <- presc_data[, .(n = .N), by = substance]
	drug_counts[, prevalence := n / n_patients]
	common_drugs <- drug_counts[prevalence >= min_prevalence, substance]

	if (length(common_drugs) < 2) {
		return(data.table(
			drug1 = character(),
			drug2 = character(),
			co_count = integer(),
			or = numeric(),
			ci_lower = numeric(),
			ci_upper = numeric(),
			p_value = numeric()
		))
	}

	presc_filtered <- presc_data[substance %in% common_drugs]

	# Create patient-drug matrix
	presc_filtered[, has_drug := 1L]
	drug_matrix <- dcast(
		presc_filtered,
		patid ~ substance,
		value.var = "has_drug",
		fun.aggregate = function(x) as.integer(length(x) > 0),
		fill = 0L
	)

	# Remove patid column and convert to matrix
	drug_matrix[, patid := NULL]
	mat <- as.matrix(drug_matrix)
	drugs <- colnames(mat)
	n_drugs <- length(drugs)

	# Get OR results WITH contingency tables
	or_results <- calc_all_ors_vectorized(mat, return_contingency = TRUE)
	# Extract matrices
	or_matrix <- or_results$or
	ci_lower_matrix <- or_results$ci_lower
	ci_upper_matrix <- or_results$ci_upper
	contingency <- or_results$contingency

	results_list <- list()
	for (i in seq_len(n_drugs - 1)) {
		for (j in (i + 1):n_drugs) {

			a <- or_results$a[i, j]
			b <- or_results$b[i, j]
			c <- or_results$c[i, j]
			d <- or_results$d[i, j]

			contingency_table <- matrix(c(a, b, c, d), nrow = 2)
			p_value <- fisher.test(contingency_table)$p.value

			results_list[[length(results_list) + 1]] <- data.table(
				drug1 = drugs[i],
				drug2 = drugs[j],
				co_count = a,  # This is the co-occurrence count
				or = or_matrix[i, j],
				ci_lower = ci_lower_matrix[i, j],
				ci_upper = ci_upper_matrix[i, j],
				p_value = p_value
			)
		}
	}

	if (length(results_list) == 0) {
		return(data.table(
			drug1 = character(),
			drug2 = character(),
			co_count = integer(),
			or = numeric(),
			ci_lower = numeric(),
			ci_upper = numeric(),
			p_value = numeric()
		))
	}

	results <- rbindlist(results_list)

	# Adjust p-values
	results[, p_adjusted := p.adjust(p_value, method = "fdr")]

	return(results)
}