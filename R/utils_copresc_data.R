#' Calculate pairwise co-prescription ORs comparing cases vs controls
#'
#' This function calculates the odds ratio of having BOTH medications
#' in cases compared to controls (proper case-control comparison).
#'
#' @param prescriptions data.table with columns: patid, substance, group
#' @param min_prevalence Minimum prevalence threshold (as proportion)
#' @param min_coprescription Minimum co-prescription threshold (as proportion)
#' @return data.table with case-control ORs for each drug pair
calculate_coprescription_ors <- function(prescriptions,
																				 min_prevalence = 0.02,
																				 min_coprescription = 0.01) {

	# Get unique patients by group
	case_patids <- unique(prescriptions[group == "case", patid])
	control_patids <- unique(prescriptions[group == "control", patid])
	n_case <- length(case_patids)
	n_control <- length(control_patids)

	# Filter drugs by individual prevalence in either cases or controls
	drug_counts <- prescriptions[, .(
		n_case = uniqueN(patid[group == "case"]),
		n_control = uniqueN(patid[group == "control"])
	), by = substance]

	drug_counts[, `:=`(
		prev_case = n_case / n_case,
		prev_control = n_control / n_control
	)]

	# Keep drugs meeting prevalence threshold in either group
	common_drugs <- drug_counts[
		prev_case >= min_prevalence | prev_control >= min_prevalence,
		substance
	]

	if (length(common_drugs) < 2) {
		return(data.table(
			drug1 = character(),
			drug2 = character(),
			or = numeric(),
			ci_lower = numeric(),
			ci_upper = numeric(),
			p_value = numeric(),
			n_case_both = integer(),
			n_control_both = integer(),
			pct_case_both = numeric(),
			pct_control_both = numeric()
		))
	}

	# Filter to common drugs
	presc_filtered <- prescriptions[substance %in% common_drugs]

	# Create separate patient-drug matrices for cases and controls
	# CASES
	case_presc <- presc_filtered[group == "case"]
	case_presc[, has_drug := 1L]
	case_matrix <- dcast(
		case_presc,
		patid ~ substance,
		value.var = "has_drug",
		fun.aggregate = function(x) as.integer(length(x) > 0),
		fill = 0L
	)
	case_matrix[, patid := NULL]

	# CONTROLS
	control_presc <- presc_filtered[group == "control"]
	control_presc[, has_drug := 1L]
	control_matrix <- dcast(
		control_presc,
		patid ~ substance,
		value.var = "has_drug",
		fun.aggregate = function(x) as.integer(length(x) > 0),
		fill = 0L
	)
	control_matrix[, patid := NULL]

	# Get drugs that actually appear in both matrices
	drugs <- sort(common_drugs)
	case_drugs <- names(case_matrix)
	control_drugs <- names(control_matrix)

	# Only keep drugs that appear in BOTH matrices
	drugs <- intersect(drugs, intersect(case_drugs, control_drugs))
	n_drugs <- length(drugs)

	if (n_drugs < 2) {
		return(data.table(
			drug1 = character(),
			drug2 = character(),
			or = numeric(),
			ci_lower = numeric(),
			ci_upper = numeric(),
			p_value = numeric(),
			n_case_both = integer(),
			n_control_both = integer(),
			pct_case_both = numeric(),
			pct_control_both = numeric()
		))
	}

	# Reorder/align columns (now safe since drugs are in both)
	case_mat <- as.matrix(case_matrix[, ..drugs])
	control_mat <- as.matrix(control_matrix[, ..drugs])

	# Calculate co-occurrence matrices for each group
	# This gives us: for each pair (i,j), how many patients have BOTH drugs
	case_coprescription <- crossprod(case_mat)  # drugs x drugs matrix
	control_coprescription <- crossprod(control_mat)

	# Now calculate case-control OR for each pair
	# For each drug pair (i,j):
	# a = cases with both drugs
	# b = cases without both drugs
	# c = controls with both drugs
	# d = controls without both drugs

	results_list <- list()

	for (i in seq_len(n_drugs - 1)) {
		for (j in (i + 1):n_drugs) {

			# Number with BOTH drugs in each group
			a <- case_coprescription[i, j]
			c <- control_coprescription[i, j]

			# Number without BOTH drugs in each group
			b <- n_case - a
			d <- n_control - c

			# Calculate percentages
			pct_case_both <- round(100 * a / n_case, 2)
			pct_control_both <- round(100 * c / n_control, 2)

			# Skip if co-prescription prevalence too low
			if (pct_case_both < (min_coprescription * 100) &&
					pct_control_both < (min_coprescription * 100)) {
				next
			}

			# Calculate OR and CI (only if all cells > 0)
			if (a > 0 && b > 0 && c > 0 && d > 0) {
				or <- (a * d) / (b * c)

				# Standard error on log scale
				se_log_or <- sqrt(1/a + 1/b + 1/c + 1/d)
				log_or <- log(or)

				ci_lower <- exp(log_or - 1.96 * se_log_or)
				ci_upper <- exp(log_or + 1.96 * se_log_or)

				# Fisher's exact test p-value
				contingency_table <- matrix(c(a, b, c, d), nrow = 2)
				p_value <- fisher.test(contingency_table)$p.value

			} else {
				or <- NA_real_
				ci_lower <- NA_real_
				ci_upper <- NA_real_
				p_value <- NA_real_
			}

			results_list[[length(results_list) + 1]] <- data.table(
				drug1 = drugs[i],
				drug2 = drugs[j],
				or = or,
				ci_lower = ci_lower,
				ci_upper = ci_upper,
				p_value = p_value,
				n_case_both = as.integer(a),
				n_control_both = as.integer(c),
				pct_case_both = pct_case_both,
				pct_control_both = pct_control_both
			)
		}
	}

	if (length(results_list) == 0) {
		return(data.table(
			drug1 = character(),
			drug2 = character(),
			or = numeric(),
			ci_lower = numeric(),
			ci_upper = numeric(),
			p_value = numeric(),
			n_case_both = integer(),
			n_control_both = integer(),
			pct_case_both = numeric(),
			pct_control_both = numeric()
		))
	}

	results <- rbindlist(results_list)

	# Adjust p-values for multiple testing
	results[, p_adjusted := p.adjust(p_value, method = "fdr")]

	# Sort by OR (descending)
	setorder(results, -or, na.last = TRUE)

	return(results)
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