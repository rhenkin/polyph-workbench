calculate_substance_net <- function(dt,
																		confounder = "",
																		min_prevalence = 0.01,
																		min_co_prevalence = 0.02,
																		p_adjust_method = "fdr") {

	# Create wide format data with one column per substance (binary indicators)
	# First, get unique patients and their outcome_age
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

	# Merge with patient data to get outcome_age
	if (confounder != "") {
		analysis_data <- merge(substance_matrix, patients, by = "patid")
	} else {
		analysis_data <- substance_matrix
	}


	# Discretize outcome_age into categories
	# This is needed for log-linear models to handle continuous covariates
	# analysis_data[, age_cat := cut(outcome_age/365.25,
	# 															 breaks = quantile(outcome_age/365.25, probs = seq(0, 1, 0.25), na.rm = TRUE),
	# 															 include.lowest = TRUE)]

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

	# Run log-linear models for each valid pair
	# We'll use the glm function with family=poisson, which fits log-linear models
	model_results <- lapply(1:nrow(valid_pairs), function(i) {
		s1 <- valid_pairs$substance1[i]
		s2 <- valid_pairs$substance2[i]
		# Make sure we're using the correct syntax for data.table subsetting
		co_prev <- valid_pairs[substance1 == s1 & substance2 == s2, co_prevalence][[1]]
		co_count_val <- valid_pairs[substance1 == s1 & substance2 == s2, co_count][[1]]
		prev1 <- substance_counts[substance == s1, prevalence][[1]]
		prev2 <- substance_counts[substance == s2, prevalence][[1]]
		# Create a contingency table for this pair with age as a stratification factor
		# This creates a data frame with counts for each combination of s1, s2, and age_cat
		if (confounder != "") {
			contingency_data <- analysis_data[, .N, by = c(s1, s2, confounder)]
			formula_str <- sprintf("N ~ %s + %s + %s + %s:%s", s1, s2, confounder, s1, s2)
		} else {
			contingency_data <- analysis_data[, .N, by = c(s1, s2)]
			formula_str <- sprintf("N ~ %s + %s + %s:%s", s1, s2, s1, s2)
		}


		# If the contingency table is empty or has zero counts, skip this pair
		if (nrow(contingency_data) == 0 || any(contingency_data$N == 0)) {
			return(NULL)
		}

		# Log-linear model formula with age_cat as a control variable
		# The : represents the interaction term we're interested in


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

		# The interaction term will be named s1:s2 or s2:s1
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
		exp_ci <- suppressMessages(exp(confint(model)))
		p_value <- coef_summary[coef_row, "Pr(>|z|)"]
		odds_ratio <- exp(beta)  # In log-linear models, exp(coefficient) is interpreted as odds ratio
		# Return as data.table
		result <- data.table(
			substance1 = s1,
			substance2 = s2,
			prevalence1 = prev1,
			prevalence2 = prev2,
			co_prevalence = co_prev,
			co_count = co_count_val,
			beta = beta,
			odds_ratio = exp(beta),
			ci_lower = exp_ci[nrow(exp_ci), 1],
			ci_upper = exp_ci[nrow(exp_ci), 2],
			p_value = p_value
		)
		return(result)
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
		#group = "substance"
	)

	# Add prevalence information to nodes
	# for (i in 1:nrow(nodes)) {
	# 	s <- nodes$id[i]
	# 	prev <- substance_counts[substance == s, prevalence]
	# 	nodes[i, prevalence := prev]
	# }
	nodes[substance_counts, prevalence := i.prevalence, on = .(id = substance)]

	# Create edges with undirected relationships
	edges <- data.table(
		from = model_results$substance1,
		to = model_results$substance2,
		value = abs(model_results$beta), # Line width based on effect size,
		odds_ratio = model_results$odds_ratio,
		title = sprintf("OR: %.2f (Co-prev: %.1f%%, 95%% CI: %2.f-%2.f p=%.3f, adj-p=%.3f)",
										model_results$odds_ratio,
										model_results$co_prevalence * 100,
										model_results$ci_lower,
										model_results$ci_upper,
										model_results$p_value,
										model_results$p_adjusted),
		color = ifelse(model_results$significance == "significant",
									 ifelse(model_results$odds_ratio > 1, "#cc0000", "#00008a"),
									 "gray")
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

create_substance_network <- function(results,
																				or_range = NULL,
																				min_edge_value = NULL,
																				max_edges = NULL,
																				node_size_by = "prevalence",
																				show_only_significant = TRUE,
																				height = "600px",
																				width = "100%") {
	# Get nodes and edges from results
	nodes <- as.data.frame(results$nodes)
	edges <- as.data.frame(results$edges)
	# Filter to show only significant relationships if requested
	if (show_only_significant) {
		significant_edges <- edges[edges$color != "gray", ]
		# Only keep nodes that have at least one significant connection
		connected_nodes <- unique(c(significant_edges$from, significant_edges$to))
		nodes <- nodes[nodes$id %in% connected_nodes, ]
		edges <- significant_edges
	}

	edges <- edges[(edges$odds_ratio >= or_range[[1]]) & (edges$odds_ratio <= or_range[[2]]), ]
	connected_nodes <- unique(c(edges$from, edges$to))
	nodes <- nodes[nodes$id %in% connected_nodes, ]
	# Filter by minimum edge value if specified
	if (!is.null(min_edge_value)) {
		edges <- edges[edges$value >= min_edge_value, ]
		# Update nodes to only include those with remaining connections
		connected_nodes <- unique(c(edges$from, edges$to))
		nodes <- nodes[nodes$id %in% connected_nodes, ]
	}
	# Limit number of edges if specified
	if (!is.null(max_edges) && nrow(edges) > max_edges) {
		# Sort by absolute value (strength of association) and keep top edges
		edges <- edges[order(-edges$value), ][1:max_edges, ]
		# Update nodes to only include those with remaining connections
		connected_nodes <- unique(c(edges$from, edges$to))
		nodes <- nodes[nodes$id %in% connected_nodes, ]
	}
	# Scale node size based on selected metric
	if (node_size_by == "prevalence") {
		# Scale node size by prevalence, with a reasonable range
		nodes$value <- nodes$prevalence * 100  # Convert to percentage for better scaling
		# Add prevalence percentage to node labels
		nodes$title <- paste0(nodes$label, "<br>Prevalence: ",
													round(nodes$prevalence * 100, 1), "%")
	} else if (node_size_by == "connections") {
		# Count connections for each node
		from_counts <- table(edges$from)
		to_counts <- table(edges$to)
		all_counts <- table(c(names(from_counts), names(to_counts)))
		# Add connection count to node size
		nodes$value <- all_counts[nodes$id]
		# Add connection count to node labels
		nodes$title <- paste0(nodes$label, "<br>Connections: ", all_counts[nodes$id])
	}
	# Create the visNetwork visualization
	network <- visNetwork(nodes, edges, height = height, width = width)
	return(network)
}