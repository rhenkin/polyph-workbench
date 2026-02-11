#' Create heatmap of top substances by temporal pattern
create_temporal_heatmap <- function(dt, top_n = 20) {

	# Get top substances by total count
	top_substances <- dt[, .N, by = substance][order(-N)][1:min(top_n, .N), substance]

	dt_filtered <- dt[substance %in% top_substances]

	# Calculate percentages within each substance
	substance_day <- dt_filtered[, .N, by = .(substance, day_group)]
	substance_day[, total := sum(N), by = substance]
	substance_day[, percent := 100 * N / total]

	spec <- list(
		`$schema` = vega_schema(),
		width = 600,
		height = 400,
		data = list(values = substance_day),
		mark = list(type = "rect", tooltip = TRUE),
		encoding = list(
			x = list(
				field = "day_group",
				type = "ordinal",
				title = "Days from prescription to event",
				axis = list(labelAngle = -45)
			),
			y = list(
				field = "substance",
				type = "nominal",
				title = "Substance",
				sort = list(field = "total", order = "descending")
			),
			color = list(
				field = "percent",
				type = "quantitative",
				title = "% within substance",
				scale = list(scheme = "blues")
			),
			tooltip = list(
				list(field = "substance", title = "Substance"),
				list(field = "day_group", title = "Period"),
				list(field = "N", title = "Cases"),
				list(field = "percent", title = "% of substance", format = ".1f")
			)
		)
	)

	return(as_vegaspec(spec))
}

#' Create histogram of days to event
#' @param dt data.table with days_to_event
#' @param pred_window numeric
#' @return vegaspec
create_temporal_histogram <- function(dt, pred_window) {

	# Count by day group
	summary_dt <- dt[, .(N = .N), by = day_group]
	summary_dt[, percent := round(100 * N / sum(N), 1)]

	# Get grouping info for proper binning
	groups <- create_day_groups(pred_window)

	spec <- list(
		`$schema` = vega_schema(),
		width = 600,
		height = 300,
		data = list(values = summary_dt),
		mark = list(
			type = "bar",
			tooltip = TRUE,
			color = "#4682B4"
		),
		encoding = list(
			x = list(
				field = "day_group",
				type = "ordinal",
				title = "Days from prescription to event",
				axis = list(labelAngle = -45)
			),
			y = list(
				field = "N",
				type = "quantitative",
				title = "Number of cases"
			),
			tooltip = list(
				list(field = "day_group", title = "Period"),
				list(field = "N", title = "Cases"),
				list(field = "percent", title = "Percent (%)")
			)
		)
	)

	return(as_vegaspec(spec))
}


#' Create faceted histogram by stratification variable
create_stratified_temporal_plot <- function(dt, strat_var, pred_window) {

	# Summary by day group and strat variable
	summary_dt <- dt[, .N, by = c("day_group", strat_var)]
	setnames(summary_dt, strat_var, "strat_value")

	# Calculate percentages within each stratum
	summary_dt[, total := sum(N), by = strat_value]
	summary_dt[, percent := round(100 * N / total, 1)]

	spec <- list(
		`$schema` = vega_schema(),
		width = 150,
		height = 150,
		data = list(values = summary_dt),
		mark = list(type = "bar", tooltip = TRUE),
		encoding = list(
			x = list(
				field = "day_group",
				type = "ordinal",
				title = "Days to event",
				axis = list(labelAngle = -45, labelLimit = 80)
			),
			y = list(
				field = "percent",
				type = "quantitative",
				title = "% of cases"
			),
			color = list(
				field = "strat_value",
				type = "nominal",
				title = strat_var,
				legend = list(orient = "bottom")
			),
			column = list(
				field = "strat_value",
				type = "nominal",
				title = NULL
			),
			tooltip = list(
				list(field = "strat_value", title = strat_var),
				list(field = "day_group", title = "Period"),
				list(field = "N", title = "Cases"),
				list(field = "percent", title = "%", format = ".1f")
			)
		)
	)

	return(as_vegaspec(spec))
}


#' Create stacked bar plot for temporal stratified analysis
#'
#' @param dt data.table with day_group and stratification variable
#' @param strat_var Name of stratification variable column
#' @param pred_window Risk window size
#' @return vegaspec object
#' @export
create_stacked_temporal_plot <- function(dt, strat_var, pred_window) {

	# Summary by day group and strat variable
	summary_dt <- dt[, .N, by = c("day_group", strat_var)]
	setnames(summary_dt, strat_var, "strat_value")

	# Calculate percentages within each day group
	summary_dt[, total_day := sum(N), by = day_group]
	summary_dt[, percent := round(100 * N / total_day, 1)]

	spec <- list(
		`$schema` = vega_schema(),
		width = 600,
		height = 400,
		data = list(values = summary_dt),
		mark = list(type = "bar", tooltip = TRUE),
		encoding = list(
			x = list(
				field = "day_group",
				type = "ordinal",
				title = "Days from prescription to event",
				axis = list(labelAngle = -45)
			),
			y = list(
				field = "N",
				type = "quantitative",
				title = "Number of cases",
				stack = TRUE
			),
			color = list(
				field = "strat_value",
				type = "nominal",
				title = strat_var,
				legend = list(orient = "right")
			),
			tooltip = list(
				list(field = "day_group", title = "Period"),
				list(field = "strat_value", title = strat_var),
				list(field = "N", title = "Cases"),
				list(field = "percent", title = "% of period", format = ".1f")
			)
		)
	)

	return(as_vegaspec(spec))
}