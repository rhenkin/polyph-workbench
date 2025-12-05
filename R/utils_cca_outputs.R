#' Create pyramid plot for burden comparison
#' @param patient_data data.table with treatment, burden column
#' @param burden_col Name of burden column (pp or n_ltc)
#' @param title Plot title
#' @return vegaspec object
create_burden_pyramid <- function(patient_data, burden_col = "pp", title = "Burden", ...) {
  upper_band <- round(quantile(patient_data[[burden_col]], 0.75) * 1.75)

  band_col <- paste0(burden_col, "_band")

  patient_data[, (band_col) := ifelse(
    get(burden_col) >= upper_band,
    paste0(upper_band, "+"),
    as.character(get(burden_col))
  )]

  plot_data <- patient_data[, .N, by = c("treatment", band_col)]
  cases_n <- nrow(patient_data[treatment == 1])
  controls_n <- nrow(patient_data[treatment == 0])

  plot_data[, pct := 0]
  plot_data[treatment == 0, pct := N / controls_n]
  plot_data[treatment == 1, pct := N / cases_n]

  pyramid_plot(
    plot_data,
    band_col,
    rev(c(2:upper_band, paste0(upper_band, "+"))),
    side_width = 175,
    title = title, ...
  ) |> as_vegaspec()
}

#' Create top substances bar plot
#' @param cases_controls data.table with group and substance
#' @param n_top Number of top substances to show
#' @param height Plot height
#' @param full_screen Whether plot is in full screen mode
#' @return vegaspec object
create_top_substances_plot <- function(cases_controls, n_top = 10, height = 250,
                                      full_screen = FALSE) {
  if (full_screen) {
    n_top <- 20
    height <- 500
  }

  cases_controls_n <- cases_controls[, .N, .(group, substance)]
  group_totals <- cases_controls[, .(total = .N), by = group]
  cases_controls_n[group_totals, pct := round(N / total * 100, 2), on = "group"]

  top_sub <- cases_controls_n[pct >= 1 & group == "case"][order(-pct)][1:n_top, substance]

  cases_controls_n[, `:=`(
    max_pct = max(pct),
    diff_label = sprintf("+%.1f%%", abs(pct[group == "case"] - pct[group == "control"]))
  ), substance]

  cases_controls_n <- cases_controls_n[substance %in% top_sub]
  cases_controls_n[nchar(substance) > 15,
    substance := paste0(strtrim(substance, 15), "...")
  ]

  grouped_bar_plot(cases_controls_n, "substance", height = height, width = 280) |>
    as_vegaspec()
}

#' Create top conditions bar plot for case-control comparison
#' @param freq_data Frequency data from calculate_frequency_stats
#' @param n_top Number of top items to show
#' @param title Plot title
#' @return vegaspec object
create_top_conditions_plot <- function(freq_data, n_top = 10, title = NULL) {
  ratios <- calculate_case_control_ratios(freq_data, "term", min_case_pct = 10)
  top_items <- ratios[order(-case_pct)][1:n_top, term]

  term_filtered <- freq_data[term %in% top_items]
  term_filtered[, `:=`(
    diff = pct[group == "case"] - pct[group == "control"],
    diff_label = sprintf("+%.1f%%", abs(pct[group == "case"] - pct[group == "control"])),
    max_pct = max(pct)
  ), by = term]

  term_filtered[nchar(term) > 15, term := paste0(strtrim(term, 15), "...")]

  grouped_bar_plot(term_filtered, "term", title = title) |> as_vegaspec()
}

#' Create wide prevalence table with ratios
#' @param freq_data Frequency data with group, item, and pct columns
#' @param item_col Name of the item column
#' @param min_pct Minimum prevalence threshold
#' @return data.table in wide format with ratios
create_prevalence_ratio_table_old <- function(freq_data, item_col, min_pct = 0.5) {
  table_data_wide <- dcast(freq_data,
    as.formula(paste(item_col, "~ group")),
    value.var = "pct",
    fill = 0
  )

  table_data_wide <- table_data_wide[case > min_pct & control > min_pct]
  table_data_wide[, ratio := round(case / control, digits = 2)]

  table_data_wide
}

#' Enhanced Create wide prevalence table with ratios and odds ratios
#'
#' This is an enhanced version of create_prevalence_ratio_table that adds
#' vectorized odds ratio calculations using the fast method from the OR heatmaps.
#'
#' @param freq_data Frequency data with group, item, and pct columns
#' @param item_col Name of the item column
#' @param min_pct Minimum prevalence threshold
#' @param data_with_group Optional: Original data.table with patid, group, and item_col for OR calculation
#' @param calculate_or Logical: whether to calculate odds ratios (requires data_with_group)
#' @param test_method Test method for p-values: "chisq" or "fisher"
#' @param p_adjust_method P-value adjustment method (default: "BH")
#' @return data.table in wide format with ratios and optionally ORs with CIs and p-values
#'
#' @details
#' This function enhances the basic prevalence ratio calculation with:
#' - Vectorized odds ratio calculation (same method as OR heatmaps)
#' - Confidence intervals for odds ratios
#' - Statistical significance testing (chi-square or Fisher's exact)
#' - Multiple testing correction
#' - Significance markers (asterisks for p_adj < 0.05)
#'
#' @examples
#' # Basic usage with just prevalence ratios (old behavior)
#' result <- create_prevalence_ratio_table(freq_data, "term", min_pct = 0.5)
#'
#' # Enhanced usage with odds ratios
#' result <- create_prevalence_ratio_table(
#'   freq_data, "term", min_pct = 0.5,
#'   data_with_group = ltcs_data,  # Must include patid, group, and term columns
#'   calculate_or = TRUE
#' )
create_prevalence_ratio_table <- function(freq_data, item_col, min_pct = 0.5,
																					data_with_group,
																					test_method = "chisq",
																					p_adjust_method = "BH") {

	# Validate inputs
	if (missing(data_with_group) || is.null(data_with_group)) {
		stop("data_with_group is required for OR calculation")
	}

	# Create wide format table with prevalences
	table_data_wide <- dcast(freq_data,
													 as.formula(paste(item_col, "~ group")),
													 value.var = "pct",
													 fill = 0)

	# Filter by minimum prevalence
	table_data_wide <- table_data_wide[case > min_pct & control > min_pct]

	if (nrow(table_data_wide) == 0) {
		return(NULL)
	}

	# === ODDS RATIO CALCULATION (Vectorized Matrix Method) ===

	# Get case and control patient IDs
	case_patids <- unique(data_with_group[group == "case", patid])
	control_patids <- unique(data_with_group[group == "control", patid])
	n_case <- length(case_patids)
	n_control <- length(control_patids)

	# Get items (substances or terms)
	items <- table_data_wide[[item_col]]

	# Build case matrix - vectorized patient x item matrix
	case_data <- data_with_group[patid %in% case_patids & get(item_col) %in% items]
	case_matrix <- dcast(case_data, patid ~ get(item_col),
											 fun.aggregate = function(x) as.integer(length(x) > 0),
											 value.var = "patid")
	case_matrix[, patid := NULL]
	case_matrix <- as.matrix(case_matrix)

	# Build control matrix - vectorized patient x item matrix
	control_data <- data_with_group[patid %in% control_patids & get(item_col) %in% items]
	control_matrix <- dcast(control_data, patid ~ get(item_col),
													fun.aggregate = function(x) as.integer(length(x) > 0),
													value.var = "patid")
	control_matrix[, patid := NULL]
	control_matrix <- as.matrix(control_matrix)

	# Vectorized calculation of 2x2 contingency table cells for ALL items at once
	item_case <- as.numeric(colSums(case_matrix))          # a: has item AND is case
	item_control <- as.numeric(colSums(control_matrix))    # b: has item AND is control
	no_item_case <- as.numeric(n_case) - item_case         # c: no item AND is case
	no_item_control <- as.numeric(n_control) - item_control # d: no item AND is control

	# Vectorized OR calculation: OR = (a * d) / (b * c)
	# Only calculate where all cells are > 0 to avoid division by zero
	valid <- (item_case > 0) & (item_control > 0) & (no_item_case > 0) & (no_item_control > 0)
	or <- rep(NA_real_, length(items))
	or[valid] <- (item_case[valid] * no_item_control[valid]) / (item_control[valid] * no_item_case[valid])

	# Vectorized CI calculation using log method
	ci_lower <- rep(NA_real_, length(items))
	ci_upper <- rep(NA_real_, length(items))

	if (any(valid)) {
		log_or <- log(or[valid])
		se_log_or <- sqrt(1/item_case[valid] + 1/item_control[valid] +
												1/no_item_case[valid] + 1/no_item_control[valid])
		ci_lower[valid] <- exp(log_or - 1.96 * se_log_or)
		ci_upper[valid] <- exp(log_or + 1.96 * se_log_or)
	}

	# Add OR and CI columns to result
	table_data_wide[, OR := round(or, 2)]
	table_data_wide[, OR_CI_lower := round(ci_lower, 2)]
	table_data_wide[, OR_CI_upper := round(ci_upper, 2)]

	# === STATISTICAL TESTING ===

	if (test_method == "fisher") {
		# Fisher test needs to be row-by-row (slower but exact for small samples)
		p_values <- sapply(seq_along(items), function(i) {
			contingency <- matrix(c(item_case[i], item_control[i],
															no_item_case[i], no_item_control[i]), nrow = 2)
			fisher.test(contingency)$p.value
		})
		table_data_wide[, p_value := p_values]

	} else if (test_method == "chisq") {
		# Fully vectorized chi-square test - VERY FAST!
		n_total <- n_case + n_control
		row1_total <- item_case + item_control          # total with item
		row2_total <- no_item_case + no_item_control    # total without item
		col1_total <- item_case + no_item_case          # = n_case
		col2_total <- item_control + no_item_control    # = n_control

		# Chi-square statistic: X^2 = n(ad - bc)^2 / [(a+b)(c+d)(a+c)(b+d)]
		chisq_stat <- (n_total * (item_case * no_item_control - item_control * no_item_case)^2) /
			(row1_total * row2_total * col1_total * col2_total)

		# P-value from chi-square distribution with df=1
		p_values <- pchisq(chisq_stat, df = 1, lower.tail = FALSE)
		table_data_wide[, p_value := p_values]
	}

	# Adjust for multiple testing (vectorized)
	table_data_wide[, p_adj := p.adjust(p_value, method = p_adjust_method)]
	table_data_wide[, p_adj := round(p_adj, 4)]

	# Add asterisk to significant items (p_adj < 0.05)
	table_data_wide[p_adj < 0.05, (item_col) := paste0(get(item_col), "*")]

	return(table_data_wide)
}

#' Format and prepare stratification choices for UI
#' @param patient_data data.table with demographic variables
#' @return Named list suitable for virtualSelect choices
create_stratification_choices <- function(patient_data) {
  heatmap_vars <- c("sex", "eth_group", "imd_quintile", "pp_group", "mltc_group")
  var_labels <- c("Sex", "Ethnicity", "IMD quintile", "# PP", "# LTC")

  choices <- setNames(
    lapply(seq_along(heatmap_vars), function(i) {
      var <- heatmap_vars[i]
      unique_vals <- sort(unique(patient_data[[var]]))

      encoded_list <- setNames(
        paste0(var, "#", unique_vals),
        unique_vals
      )
      encoded_list
    }),
    var_labels
  )

  choices
}

calculate_prevalence_cca <- function(data1, data2, selected_values,
																		 filter_col, group_col,
																		 test_method = "chisq",
																		 p_adjust_method = "BH") {

	# Step 1: Filter to patients with selected values in data2
	patids <- unique(data2[get(filter_col) %in% selected_values, patid])

	if (length(patids) == 0) {
		return(NULL)
	}

	# Step 2: Filter data1 to these patients
	data1_filtered <- data1[patid %in% patids]

	# Step 3: Calculate frequencies for these filtered patients
	# This creates the freq_data that create_prevalence_ratio_table expects
	freq <- data1_filtered[, .N, by = c("group", group_col)]
	# group_totals <- data1_filtered[, .(total = uniqueN(patid)), by = group]
	group_totals <- data1_filtered[, .(total = .N), by = group]
	freq[group_totals, pct := round(N / total * 100, 2), on = "group"]

	# Step 4: Reuse create_prevalence_ratio_table to get ORs and statistics
	# Note: min_pct = 1 means we need at least 1% prevalence in BOTH case and control
	result <- create_prevalence_ratio_table(
		freq_data = freq,
		item_col = group_col,
		min_pct = 1,  # Minimum 1% prevalence in both groups
		data_with_group = data1_filtered,  # Pass filtered data for OR calculation
		test_method = test_method,
		p_adjust_method = p_adjust_method
	)

	return(result)
}

# calculate_prevalence_cca <- function(data1, data2, selected_values,
# 																		 filter_col, group_col,
# 																		 test_significance = TRUE,
# 																		 p_adjust_method = "BH",
# 																		 test_method = "chisq",
# 																		 use_odds_ratio = TRUE) {
# 	patids <- unique(data2[get(filter_col) %in% selected_values, patid])
#
# 	presc_freq <- data1[patid %in% patids,
# 											list(
# 												N_with_disease = uniqueN(patid),
# 												Prevalence = round(100 * (uniqueN(patid) / length(patids)), digits = 2)
# 											),
# 											by = c("group", group_col)
# 	]
#
# 	result <- dcast(presc_freq,
# 									as.formula(paste(group_col, "~ group")),
# 									value.var = "Prevalence",
# 									fill = 0)
# 	if (nrow(result) == 0) return(NULL)
# 	result <- result[case >= 1 & control >= 1]
#
# 	if (use_odds_ratio && test_significance) {
# 		# Get case and control patient IDs
# 		case_patids <- unique(data1[patid %in% patids & group == "case", patid])
# 		control_patids <- unique(data1[patid %in% patids & group == "control", patid])
# 		n_case <- length(case_patids)
# 		n_control <- length(control_patids)
#
# 		# Get items
# 		items <- result[[group_col]]
#
# 		# Build case matrix
# 		case_data <- data1[patid %in% case_patids & get(group_col) %in% items]
# 		case_matrix <- dcast(case_data, patid ~ get(group_col),
# 												 fun.aggregate = function(x) as.integer(length(x) > 0),
# 												 value.var = "patid")
# 		case_matrix[, patid := NULL]
# 		case_matrix <- as.matrix(case_matrix)
#
# 		# Build control matrix
# 		control_data <- data1[patid %in% control_patids & get(group_col) %in% items]
# 		control_matrix <- dcast(control_data, patid ~ get(group_col),
# 														fun.aggregate = function(x) as.integer(length(x) > 0),
# 														value.var = "patid")
# 		control_matrix[, patid := NULL]
# 		control_matrix <- as.matrix(control_matrix)
#
# 		# Vectorized calculation of contingency table cells (convert to numeric to avoid overflow)
# 		item_case <- as.numeric(colSums(case_matrix))
# 		item_control <- as.numeric(colSums(control_matrix))
# 		no_item_case <- as.numeric(n_case) - item_case
# 		no_item_control <- as.numeric(n_control) - item_control
#
# 		# Vectorized OR calculation
# 		valid <- (item_case > 0) & (item_control > 0) & (no_item_case > 0) & (no_item_control > 0)
# 		or <- rep(NA_real_, length(items))
# 		or[valid] <- (item_case[valid] * no_item_control[valid]) / (item_control[valid] * no_item_case[valid])
#
# 		# Vectorized CI calculation
# 		ci_lower <- rep(NA_real_, length(items))
# 		ci_upper <- rep(NA_real_, length(items))
#
# 		if (any(valid)) {
# 			log_or <- log(or[valid])
# 			se_log_or <- sqrt(1/item_case[valid] + 1/item_control[valid] +
# 													1/no_item_case[valid] + 1/no_item_control[valid])
# 			ci_lower[valid] <- exp(log_or - 1.96 * se_log_or)
# 			ci_upper[valid] <- exp(log_or + 1.96 * se_log_or)
# 		}
#
# 		# Add to result
# 		result[, OR := round(or, 2)]
# 		result[, OR_CI_lower := round(ci_lower, 2)]
# 		result[, OR_CI_upper := round(ci_upper, 2)]
#
# 		if (test_method == "fisher") {
# 			# Fisher test still needs to be done row-by-row
# 			p_values <- sapply(seq_along(items), function(i) {
# 				contingency <- matrix(c(item_case[i], item_control[i],
# 																no_item_case[i], no_item_control[i]), nrow = 2)
# 				fisher.test(contingency)$p.value
# 			})
# 			result[, p_value := p_values]
#
# 		} else if (test_method == "chisq") {
# 			# Fully vectorized chi-square (correct formula)
# 			n_total <- n_case + n_control
# 			row1_total <- item_case + item_control
# 			row2_total <- no_item_case + no_item_control
# 			col1_total <- item_case + no_item_case  # = n_case
# 			col2_total <- item_control + no_item_control  # = n_control
#
# 			# Correct chi-square formula
# 			chisq_stat <- (n_total * (item_case * no_item_control - item_control * no_item_case)^2) /
# 				(row1_total * row2_total * col1_total * col2_total)
#
# 			p_values <- pchisq(chisq_stat, df = 1, lower.tail = FALSE)
# 			result[, p_value := p_values]
# 		}
#
# 		# Adjust for multiple testing (vectorized)
# 		result[, p_adj := p.adjust(p_value, method = p_adjust_method)]
#
# 		# Add asterisk to significant items (vectorized)
# 		result[p_adj < 0.05, (group_col) := paste0(get(group_col), "*")]
#
# 		result[, p_adj := round(p_adj, 4)]
#
#
# 	} else {
# 		# Fall back to simple prevalence ratio
# 		result[, Prevalence_Ratio := round(case / control, digits = 2)]
# 		result[is.infinite(Prevalence_Ratio), Prevalence_Ratio := 0]
# 	}
#
# 	result
# }