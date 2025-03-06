# Helper function to format results
format_results <- function(df_stats, chisq_tests, cat_totals, demog_var, subst_frequency) {
	# Create wide format results
	result <- dcast(
		df_stats,
		substance ~ get(demog_var),
		value.var = "pct",
		fill = 0,
		fun.aggregate = max
	)

	# Add p-values
	result[chisq_tests, `:=`(pvalue = pvalue, padj = i.padj), on = "substance"]

	# Rename columns with proportions
	old_names <- as.character(unique(df_stats[[demog_var]]))
	new_names <- paste0(
		old_names,
		" (",
		signif(cat_totals$Total/sum(cat_totals$Total), digit = 2),
		")"
	)
	setnames(result, old = old_names, new = new_names)

	# Add frequency information
	result <- merge(
		result,
		subst_frequency[, .(substance, N, pct_total)],
		by = "substance"
	)

	return(result[order(-N, na.last = TRUE)])
}