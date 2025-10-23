# Helper function to format results
format_results <- function(df_stats, chisq_tests, cat_totals, demog_var, frequency_data, condition_col) {
	# Create formula for dcast dynamically
	formula_str <- paste(condition_col, "~", demog_var)

	# Create wide format results
	result <- dcast(df_stats,
									formula = as.formula(formula_str),
									value.var = "pct",
									fill = 0,
									fun.aggregate = max)

	# Add p-values
	result[chisq_tests, `:=`(pvalue = pvalue, padj = i.padj), on = condition_col]

	# Rename columns with proportions
	old_names <- as.character(unique(df_stats[[demog_var]]))
	# Make sure cat_totals order matches old_names order
	cat_totals_ordered <- cat_totals[match(old_names, cat_totals[[demog_var]]),]
	new_names <- paste0(
		old_names,
		" (",
		signif(cat_totals_ordered$Total/sum(cat_totals_ordered$Total), digit = 2),
		")"
	)
	setnames(result, old = old_names, new = new_names)

	# Add frequency information
	merge_cols <- c(condition_col, "N", "pct_total")
	result <- merge(result, frequency_data[, ..merge_cols], by = condition_col)

	# Format p-values and pct_total for better display
	result[, pvalue := ifelse(is.na(pvalue), NA_real_,
														ifelse(pvalue < 0.001, 0, round(pvalue, 3)))]
	result[, padj := ifelse(is.na(padj), NA_real_,
													ifelse(padj < 0.001, 0, round(padj, 3)))]
	result[, pct_total := round(pct_total, 2)]

	return(result[order(-N, na.last = TRUE)])
}