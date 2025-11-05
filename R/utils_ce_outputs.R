#' Utility functions for Case Explorer (Outcome Explorer) outputs
#' These create visualizations for cases only (not case-control comparisons)

#' Create simple horizontal bar plot for single group
#' @param df data.table with item column and pct column
#' @param y_var Name of the y-axis variable (item column)
#' @param height Plot height
#' @param width Plot width
#' @param title Plot title
#' @return vega spec list
simple_bar_plot <- function(df, y_var, height = 250, width = 300, title = NULL, sort_values= NULL) {
	spec <- list(
		`$schema` = vega_schema(),
		title = title,
		height = height,
		width = width,
		data = list(values = df),
		mark = list(
			type = "bar",
			color = "#e74c3c"  # Red color matching cases in CCA
		),
		encoding = list(
			x = list(
				field = "pct",
				type = "quantitative",
				title = "Prevalence (%)",
				scale = list(domainMax = max(df$pct) * 1.05)
			),
			y = list(
				field = y_var,
				type = "nominal",
				axis = list(grid = FALSE, title = NULL),
				sort = sort_values
			)
		),
		config = list(
			view = list(stroke = NULL),
			font = "Lato, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, \"Helvetica Neue\", Arial, sans-serif, \"Apple Color Emoji\", \"Segoe UI Emoji\", \"Segoe UI Symbol\"",
			axis = list(labelFontSize = 12)
		)
	)
	spec
}

#' Create PP distribution bar plot for cases
#' @param outcome_prescriptions data.table with outcome prescriptions
#' @param pp_groups_data data.table with PP groups
#' @param title Plot title
#' @return vegaspec object
create_pp_distribution_plot <- function(pp_groups_data, height = 250, width = 280, title = "Polypharmacy Distribution") {
	pp_dist <- pp_groups_data[, .N, by = pp_group]
	total <- sum(pp_dist$N)
	pp_dist[, pct := round(N / total * 100, 2)]
	pp_dist$pp_group <- factor(	pp_dist$pp_group, labels = c("2-4", "5-9", "10+"))
	simple_bar_plot(
		pp_dist,
		"pp_group",
		height = 250,
		width = width,
		title = title,
		sort_values = list("2-4", "5-9", "10+")
	) |>
		as_vegaspec()
}

#' Create top 10 LTCs bar plot for cases
#' @param ltc_data data.table with LTC data
#' @param outcome_prescriptions data.table with outcome prescriptions (to filter patients)
#' @param n_top Number of top items to show
#' @param title Plot title
#' @return vegaspec object
create_top_ltcs_plot <- function(ltc_data, outcome_prescriptions, n_top = 10, height = 250, width = 280, title = "Top 10 Long-term Conditions") {
	# Filter to patients in outcome_prescriptions
	valid_patids <- unique(outcome_prescriptions$patid)
	ltc_filtered <- ltc_data[patid %in% valid_patids]

	# Calculate frequency
	ltc_freq <- ltc_filtered[, .(N = uniqueN(patid)), by = term]
	total_patients <- uniqueN(ltc_filtered$patid)
	ltc_freq[, pct := round(N / total_patients * 100, 2)]

	# Get top N
	setorder(ltc_freq, -pct)
	top_ltcs <- ltc_freq[1:min(n_top, nrow(ltc_freq))]

	# Truncate long names
	top_ltcs[nchar(term) > 25, term := paste0(strtrim(term, 25), "...")]

	simple_bar_plot(top_ltcs, "term", height = height, width = width, title = title) |>
		as_vegaspec()
}

#' Create top 10 substances bar plot for cases
#' @param outcome_prescriptions data.table with outcome prescriptions
#' @param n_top Number of top items to show
#' @param title Plot title
#' @return vegaspec object
create_top_substances_plot_ce <- function(outcome_prescriptions, n_top = 10, height = 250, width = 280, title = "Top 10 Substances") {
	# Calculate frequency
	subst_freq <- outcome_prescriptions[, .(N = uniqueN(patid)), by = substance]
	total_patients <- uniqueN(outcome_prescriptions$patid)
	subst_freq[, pct := round(N / total_patients * 100, 2)]

	# Get top N
	setorder(subst_freq, -pct)
	top_substances <- subst_freq[1:min(n_top, nrow(subst_freq))]

	# Truncate long names
	top_substances[nchar(substance) > 25, substance := paste0(strtrim(substance, 25), "...")]

	simple_bar_plot(top_substances, "substance", height = height, width = width, title = title) |>
		as_vegaspec()
}