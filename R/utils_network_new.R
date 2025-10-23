# Function to fit the log-linear model for a substance pair
fit_substance_pair_model <- function(
		analysis_data,
		substance_counts,
		s1,
		s2,
		confounder = "",
		stratify = FALSE) {  # Keep stratify parameter for internal function logic

	# Get prevalence values
	prev1 <- substance_counts[substance == s1, prevalence][[1]]
	prev2 <- substance_counts[substance == s2, prevalence][[1]]

	# Calculate co-prevalence
	co_count_val <- sum(analysis_data[[s1]] & analysis_data[[s2]])
	co_prev <- co_count_val / nrow(analysis_data)

	# If not stratifying, fit a single model
	if (!stratify || confounder == "") {
		contingency_data <- analysis_data[, .N, by = c(s1, s2)]
		formula_str <- sprintf("N ~ %s + %s + %s:%s", s1, s2, s1, s2)

		# Fit the log-linear model
		model <- tryCatch(
			glm(formula = as.formula(formula_str),
					family = poisson(link = "log"),
					data = contingency_data),
			error = function(e) NULL
		)

		if (is.null(model)) {
			return(NULL)
		}

		# Extract coefficient and p-value for the interaction term
		coef_summary <- summary(model)$coefficients
		interaction_term <- paste0(s1, ":", s2)
		alt_interaction_term <- paste0(s2, ":", s1)

		coef_row <- grep(paste0("^", interaction_term, "$"), rownames(coef_summary))
		if (length(coef_row) == 0) {
			coef_row <- grep(paste0("^", alt_interaction_term, "$"), rownames(coef_summary))
		}

		if (length(coef_row) == 0) {
			return(NULL)
		}

		beta <- coef_summary[coef_row, "Estimate"]
		se <- coef_summary[coef_row, "Std. Error"]
		p_value <- coef_summary[coef_row, "Pr(>|z|)"]

		# Calculate confidence intervals
		exp_ci <- suppressMessages(exp(confint(model)))
		ci_lower <- exp_ci[nrow(exp_ci), 1]
		ci_upper <- exp_ci[nrow(exp_ci), 2]

		result <- data.table(
			substance1 = s1,
			substance2 = s2,
			prevalence1 = prev1,
			prevalence2 = prev2,
			co_prevalence = co_prev,
			co_count = co_count_val,
			beta = beta,
			se = se,
			odds_ratio = exp(beta),
			ci_lower = ci_lower,
			ci_upper = ci_upper,
			p_value = p_value,
			model_type = "single"
		)

		return(result)
	} else {
		# For stratified models, we'll fit a model for each level of the confounder
		# and then combine them using meta-analysis

		# Get unique values of the confounder
		strata <- unique(analysis_data[[confounder]])

		# Fit model for each stratum
		stratum_results <- lapply(strata, function(stratum_val) {
			# Subset data for this stratum
			stratum_data <- analysis_data[get(confounder) == stratum_val]

			if (nrow(stratum_data) < 10) {
				return(NULL)  # Skip strata with too few observations
			}

			# Create contingency table for this stratum
			contingency_data <- stratum_data[, .N, by = c(s1, s2)]

			# If the contingency table is empty or has zero counts, skip
			if (nrow(contingency_data) == 0 || any(contingency_data$N == 0)) {
				return(NULL)
			}

			formula_str <- sprintf("N ~ %s + %s + %s:%s", s1, s2, s1, s2)

			model <- tryCatch(
				glm(formula = as.formula(formula_str),
						family = poisson(link = "log"),
						data = contingency_data),
				error = function(e) NULL
			)

			if (is.null(model)) {
				return(NULL)
			}

			# Extract coefficient and standard error for the interaction term
			coef_summary <- summary(model)$coefficients
			interaction_term <- paste0(s1, ":", s2)
			alt_interaction_term <- paste0(s2, ":", s1)

			coef_row <- grep(paste0("^", interaction_term, "$"), rownames(coef_summary))
			if (length(coef_row) == 0) {
				coef_row <- grep(paste0("^", alt_interaction_term, "$"), rownames(coef_summary))
			}

			if (length(coef_row) == 0) {
				return(NULL)
			}

			beta <- coef_summary[coef_row, "Estimate"]
			se <- coef_summary[coef_row, "Std. Error"]
			p_value <- coef_summary[coef_row, "Pr(>|z|)"]

			data.table(
				stratum = stratum_val,
				beta = beta,
				se = se,
				p_value = p_value,
				weight = nrow(stratum_data)
			)
		})

		# Combine non-null results
		stratum_results <- rbindlist(stratum_results[!sapply(stratum_results, is.null)])

		if (nrow(stratum_results) == 0) {
			return(NULL)  # No valid strata models
		}

		# Perform meta-analysis if we have at least 2 strata
		if (nrow(stratum_results) >= 2) {
			# Use the rma function from metafor for random-effects meta-analysis
			meta_res <- tryCatch(
				metafor::rma(yi = beta, sei = se, weights = weight, data = stratum_results),
				error = function(e) NULL
			)

			if (!is.null(meta_res)) {
				# Extract the meta-analysis results
				meta_beta <- meta_res$b[1]  # Overall effect
				meta_se <- meta_res$se  # Standard error of overall effect
				meta_p <- meta_res$pval  # P-value
				meta_ci_lower <- exp(meta_res$ci.lb)  # Lower CI
				meta_ci_upper <- exp(meta_res$ci.ub)  # Upper CI
				meta_odds_ratio <- exp(meta_beta)  # Odds ratio

				result <- data.table(
					substance1 = s1,
					substance2 = s2,
					prevalence1 = prev1,
					prevalence2 = prev2,
					co_prevalence = co_prev,
					co_count = co_count_val,
					beta = meta_beta,
					se = meta_se,
					odds_ratio = meta_odds_ratio,
					ci_lower = meta_ci_lower,
					ci_upper = meta_ci_upper,
					p_value = meta_p,
					model_type = "meta-analysis",
					n_strata = nrow(stratum_results),
					i2 = meta_res$I2  # Heterogeneity measure
				)

				return(result)
			}
		}

		# If meta-analysis failed or we only have one stratum, use a weighted average instead
		total_weight <- sum(stratum_results$weight)
		weighted_beta <- sum(stratum_results$beta * stratum_results$weight) / total_weight

		# Approximate the standard error for the weighted average
		# This is a simplified approach - in a real application you might want a more precise formula
		weighted_se <- sqrt(sum((stratum_results$se^2) * (stratum_results$weight/total_weight)^2))

		weighted_p <- 2 * pnorm(-abs(weighted_beta / weighted_se))  # Two-tailed p-value
		weighted_ci_lower <- exp(weighted_beta - 1.96 * weighted_se)
		weighted_ci_upper <- exp(weighted_beta + 1.96 * weighted_se)

		result <- data.table(
			substance1 = s1,
			substance2 = s2,
			prevalence1 = prev1,
			prevalence2 = prev2,
			co_prevalence = co_prev,
			co_count = co_count_val,
			beta = weighted_beta,
			se = weighted_se,
			odds_ratio = exp(weighted_beta),
			ci_lower = weighted_ci_lower,
			ci_upper = weighted_ci_upper,
			p_value = weighted_p,
			model_type = "weighted-average",
			n_strata = nrow(stratum_results)
		)

		return(result)
	}
}

# Main function for substance network analysis
calculate_substance_net_meta <- function(dt,
																		confounder = "",
																		min_prevalence = 0.01,
																		min_co_prevalence = 0.02,
																		p_adjust_method = "fdr") {

	# Create wide format data with one column per substance (binary indicators)
	# First, get unique patients and their confounder value if applicable
	if (confounder != "") {
		patients <- unique(dt[, c("patid", confounder), with = FALSE])
	} else {
		patients <- unique(dt[, .(patid)])
	}

	dt[, substance := make.names(substance)]
	# Calculate substance prevalence
	substance_counts <- dt[, .N, by = substance]
	substance_counts[, prevalence := N / length(unique(dt$patid))]

	# Filter substances by minimum prevalence
	valid_substances <- substance_counts[prevalence >= min_prevalence, substance]
	message(sprintf("Found %d substances with prevalence >= %.1f%%",
									length(valid_substances), min_prevalence * 100))

	# Only keep rows with valid substances
	dt_filtered <- dt[substance %in% valid_substances]

	# Create a binary matrix of patient-substance combinations
	substance_matrix <- dcast(dt_filtered, patid ~ substance,
														fun.aggregate = length, value.var = "substance")

	# Replace counts > 0 with 1 to create binary indicators
	substance_cols <- setdiff(names(substance_matrix), "patid")
	for (col in substance_cols) {
		substance_matrix[, (col) := as.integer(get(col) > 0)]
	}

	# Merge with patient data to get the confounder
	if (confounder != "") {
		analysis_data <- merge(substance_matrix, patients, by = "patid")
	} else {
		analysis_data <- substance_matrix
	}

	# Calculate co-prevalence between all pairs of substances
	n_patients <- nrow(analysis_data)
	substance_cols <- setdiff(names(substance_matrix), "patid")

	# Create matrix representation for faster calculations
	substance_matrix_only <- as.matrix(substance_matrix[, ..substance_cols])

	# Calculate co-occurrence matrix directly
	co_occurrence <- t(substance_matrix_only) %*% substance_matrix_only

	# Convert to co-prevalence
	co_prevalence_matrix <- co_occurrence / n_patients

	# Convert to data.table format
	substance_pairs <- data.table(expand.grid(
		substance1 = substance_cols,
		substance2 = substance_cols,
		stringsAsFactors = FALSE
	))

	# Only keep pairs where substance1 < substance2 (to avoid redundancy)
	substance_pairs <- substance_pairs[substance1 < substance2]

	# Get co-prevalence values from the matrix
	co_prevalence_data <- substance_pairs[, .(
		substance1 = substance1,
		substance2 = substance2,
		co_count = unlist(mapply(function(s1, s2) co_occurrence[s1, s2], substance1, substance2)),
		co_prevalence = unlist(mapply(function(s1, s2) co_prevalence_matrix[s1, s2], substance1, substance2))
	)]

	# Filter pairs by minimum co-prevalence
	valid_pairs <- co_prevalence_data[co_prevalence >= min_co_prevalence]
	message(sprintf("Found %d substance pairs with co-prevalence >= %.1f%%",
									nrow(valid_pairs), min_co_prevalence * 100))

	# Run models for each valid pair
	model_results <- lapply(1:nrow(valid_pairs), function(i) {
		s1 <- valid_pairs$substance1[i]
		s2 <- valid_pairs$substance2[i]

		# Call the separate model fitting function
		fit_substance_pair_model(
			analysis_data = analysis_data,
			substance_counts = substance_counts,
			s1 = s1,
			s2 = s2,
			confounder = confounder,
			stratify = confounder != ""  # Automatically stratify if confounder is present
		)
	})

	model_results <- rbindlist(model_results[!sapply(model_results, is.null)])

	# Adjust p-values for multiple comparisons
	if (nrow(model_results) > 0) {
		model_results[, expected_prevalence := prevalence1 * prevalence2]
		model_results[, prevalence_ratio := co_prevalence / expected_prevalence]
		model_results[, p_adjusted := p.adjust(p_value, method = p_adjust_method)]
		# Add significance flags based on adjusted p-values
		model_results[, significance := ifelse(p_adjusted < 0.05, "significant", "not significant")]
	}

	# Create nodes for visNetwork
	nodes <- data.table(
		id = unique(c(model_results$substance1, model_results$substance2)),
		label = unique(c(model_results$substance1, model_results$substance2)),
		title = unique(c(model_results$substance1, model_results$substance2))
	)

	# Add prevalence information to nodes
	nodes[substance_counts, prevalence := i.prevalence, on = .(id = substance)]

	# Create edges with undirected relationships
	edges <- data.table(
		from = model_results$substance1,
		to = model_results$substance2,
		value = abs(model_results$beta), # Line width based on effect size
		odds_ratio = model_results$odds_ratio,
		title = sprintf("OR: %.2f (Co-prev: %.1f%%, 95%% CI: %.2f-%.2f, p=%.3f, adj-p=%.3f)",
										model_results$odds_ratio,
										model_results$co_prevalence * 100,
										model_results$ci_lower,
										model_results$ci_upper,
										model_results$p_value,
										model_results$p_adjusted),
		color = ifelse(model_results$significance == "significant",
									 ifelse(model_results$odds_ratio > 1, "#cc0000", "#00008a"),
									 "gray"),
		# Add model type to edge attributes
		model_type = model_results$model_type
	)

	# Return results
	return(list(
		nodes = nodes,
		edges = edges,
		models = model_results,
		substance_counts = substance_counts,
		co_prevalence_data = co_prevalence_data
	))
}