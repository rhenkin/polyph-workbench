# Helper function to perform chi-square tests
perform_chisq_tests <- function(df_counts, cat_totals, demog_var) {
	cohort_props <- cat_totals[, .(prop = Total/sum(Total))]

	chisq_tests <- df_counts[, {
		obs <- c(N_category)
		pval <- NA_real_

		if (length(obs) == length(cat_totals[[demog_var]]) && sum(N_category) >= 10) {
			test <- suppressWarnings(chisq.test(obs, p = cohort_props$prop))
			pval <- test$p.value
		}

		.(pvalue = pval)
	}, by = substance]

	chisq_tests[, padj := p.adjust(pvalue, method = "fdr")]
	return(chisq_tests)
}
