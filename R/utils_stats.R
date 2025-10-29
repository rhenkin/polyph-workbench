perform_chisq_tests <- function(df_counts, cat_totals, demog_var, condition_col) {
	cohort_props <- cat_totals[, .(prop = Total/sum(Total))]

	# Use by= with the condition column variable
	chisq_tests <- df_counts[, {
		obs <- c(N_category)
		pval <- NA_real_
		if (length(obs) == length(cat_totals[[demog_var]]) && sum(N_category) >= 10) {
			test <- suppressWarnings(chisq.test(obs, p = cohort_props$prop))
			pval <- test$p.value
		}
		.(pvalue = pval)
	}, by = condition_col]

	chisq_tests[, padj := p.adjust(pvalue, method = "bonferroni")]
	return(chisq_tests)
}

calc_all_ors_vectorized <- function(drug_matrix) {
	n_drugs <- ncol(drug_matrix)
	n_patients <- nrow(drug_matrix)
	drug_names <- colnames(drug_matrix)

	# Pre-calculate everything we need
	drug_sums <- colSums(drug_matrix)
	cross_products <- crossprod(drug_matrix)  # t(drug_matrix) %*% drug_matrix

	# CORRECTED: Create matrices for all combinations
	drug_sums_i <- matrix(drug_sums, n_drugs, n_drugs, byrow = TRUE)   # Row sums repeated (drug i)
	drug_sums_j <- matrix(drug_sums, n_drugs, n_drugs, byrow = FALSE)  # Column sums repeated (drug j)

	# Calculate a, b, c, d for ALL pairs at once
	a <- cross_products                           # both drugs
	b <- drug_sums_i - a                         # drug i only (now correct)
	c <- drug_sums_j - a                         # drug j only (now correct)
	d <- n_patients - a - b - c                  # neither

	# Vectorized OR calculation (avoid division by zero)
	valid <- (a > 0) & (b > 0) & (c > 0) & (d > 0)

	or_matrix <- matrix(NA, n_drugs, n_drugs)
	ci_lower_matrix <- matrix(NA, n_drugs, n_drugs)
	ci_upper_matrix <- matrix(NA, n_drugs, n_drugs)

	# Only calculate where valid
	or_matrix[valid] <- (a * d / (b * c))[valid]

	# FIXED: Only calculate CI where valid
	log_or_valid <- log(or_matrix[valid])
	se_log_or_valid <- sqrt((1/a + 1/b + 1/c + 1/d)[valid])

	ci_lower_matrix[valid] <- exp(log_or_valid - 1.96 * se_log_or_valid)
	ci_upper_matrix[valid] <- exp(log_or_valid + 1.96 * se_log_or_valid)

	# Set diagonal to NA (drug with itself isn't meaningful)
	diag(or_matrix) <- NA
	diag(ci_lower_matrix) <- NA
	diag(ci_upper_matrix) <- NA

	rownames(or_matrix) <- colnames(or_matrix) <- drug_names
	rownames(ci_lower_matrix) <- colnames(ci_lower_matrix) <- drug_names
	rownames(ci_upper_matrix) <- colnames(ci_upper_matrix) <- drug_names

	return(list(or = or_matrix, ci_lower = ci_lower_matrix, ci_upper = ci_upper_matrix))
}