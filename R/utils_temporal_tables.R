# R/utils_temporal_tables.R

#' Create temporal substance table showing distribution across day groups
#'
#' @param dt data.table with substance, day_group, and patient identifiers
#' @return data.table in wide format with percentages only
#' @export
create_temporal_substance_table <- function(dt) {

	# Count by substance and day group
	substance_day <- dt[, .N, by = .(substance, day_group)]

	# Pivot to wide format
	substance_wide <- dcast(substance_day,
													substance ~ day_group,
													value.var = "N",
													fill = 0)

	# Calculate total across all day groups
	day_cols <- setdiff(names(substance_wide), "substance")
	substance_wide[, Total_N := rowSums(.SD), .SDcols = day_cols]

	# Calculate percentages for each day group (% of substance's total in each day)
	for (col in day_cols) {
		pct_col <- paste0(col, "_pct")
		substance_wide[, (pct_col) := round(100 * get(col) / Total_N, 1)]
	}

	# Sort by total count (most common substances first)
	setorder(substance_wide, -Total_N)

	# Create final display columns: substance, Total_N, then percentages only
	pct_cols <- paste0(day_cols, "_pct")
	final_cols <- c("substance", "Total_N", pct_cols)

	# Rename percentage columns to remove _pct suffix for display
	result <- substance_wide[, ..final_cols]

	# Rename columns for cleaner display
	for (i in seq_along(day_cols)) {
		setnames(result, paste0(day_cols[i], "_pct"), day_cols[i])
	}

	result <- result[Total_N>10]

	return(result)
}


#' Render temporal substance table with reactable (percentage-only version)
#'
#' @param dt data.table from create_temporal_substance_table
#' @return reactable object
#' @export
render_temporal_substance_table <- function(dt) {

	# Get day group columns (everything except substance and Total_N)
	day_cols <- setdiff(names(dt), c("substance", "Total_N"))

	# Create column definitions
	col_defs <- list(
		substance = colDef(
			name = "Substance",
			minWidth = 200,
			sticky = "left",
			style = list(fontWeight = 600, borderRight = "2px solid #dee2e6")
		),
		Total_N = colDef(
			name = "Total Cases",
			minWidth = 100,
			align = "center",
			style = list(fontWeight = 600, backgroundColor = "#f8f9fa")
		)
	)

	# Add columns for each day group (showing percentages)
	for (col in day_cols) {
		col_defs[[col]] <- colDef(
			name = col,
			minWidth = 100,
			align = "center",
			cell = function(value) {
				if (value == 0) {
					span(style = "color: #ccc;", "—")
				} else {
					# Color code by percentage
					color <- if (value >= 30) {
						"#d32f2f"  # Red for high %
					} else if (value >= 20) {
						"#f57c00"  # Orange for medium-high %
					} else if (value >= 10) {
						"#fbc02d"  # Yellow for medium %
					} else {
						"#666"     # Gray for low %
					}

					span(
						style = sprintf("font-weight: 600; color: %s;", color),
						sprintf("%.1f%%", value)
					)
				}
			}
			# style = function(value) {
			# 	# Background color gradient based on percentage
			# 	if (value == 0) {
			# 		list(backgroundColor = "#fff")
			# 	} else {
			# 		# Light gradient from white to blue
			# 		opacity <- min(value / 50, 1)  # Cap at 50% for full opacity
			# 		list(backgroundColor = sprintf("rgba(33, 150, 243, %.2f)", opacity * 0.2))
			# 	}
			# }
		)
	}

	reactable(
		dt,
		columns = col_defs,
		defaultPageSize = 20,
		searchable = TRUE,
		highlight = TRUE,
		bordered = TRUE,
		striped = FALSE,  # Turn off striping since we're using background colors
		compact = TRUE,
		defaultSorted = "Total_N",
		defaultSortOrder = "desc",
		theme = reactableTheme(
			headerStyle = list(
				background = "#f8f9fa",
				borderBottom = "2px solid #dee2e6",
				fontWeight = 600
			)
		),
		# Add a summary footer showing that each row sums to 100%
		defaultColDef = colDef(
			footerStyle = list(fontWeight = 600)
		)
	)
}

#' Create summary statistics table for temporal distribution
#'
#' @param dt data.table with day_group column
#' @return reactable object
#' @export
create_temporal_summary_stats <- function(dt) {

	# Calculate statistics by day group
	summary_stats <- dt[, .(
		N = .N,
		`% of Total` = round(.N / nrow(dt) * 100, 1),
		`Mean Days` = round(mean(days_to_event), 1),
		`Median Days` = round(median(days_to_event), 1),
		`Min Days` = min(days_to_event),
		`Max Days` = max(days_to_event)
	), by = day_group]

	# Add overall row
	overall <- data.table(
		day_group = "Overall",
		N = nrow(dt),
		`% of Total` = 100.0,
		`Mean Days` = round(mean(dt$days_to_event), 1),
		`Median Days` = round(median(dt$days_to_event), 1),
		`Min Days` = min(dt$days_to_event),
		`Max Days` = max(dt$days_to_event)
	)

	summary_stats <- rbindlist(list(summary_stats, overall), fill = TRUE)

	setorder(summary_stats, day_group)

	reactable(
		summary_stats,
		columns = list(
			day_group = colDef(name = "Period", minWidth = 120),
			N = colDef(name = "Cases", format = colFormat(separators = TRUE)),
			`% of Total` = colDef(name = "% of Total", format = colFormat(digits = 1)),
			`Mean Days` = colDef(name = "Mean Days", format = colFormat(digits = 1)),
			`Median Days` = colDef(name = "Median Days", format = colFormat(digits = 1)),
			`Min Days` = colDef(name = "Min", align = "right"),
			`Max Days` = colDef(name = "Max", align = "right")
		),
		highlight = TRUE,
		bordered = TRUE,
		compact = TRUE,
		defaultPageSize = 20,
		style = list(fontSize = "0.9em")
	)
}

#' Create top substances table for early vs late comparison
#'
#' @param dt data.table with substance column
#' @param top_n Number of top substances to show
#' @return reactable object
#' @export
create_top_substances_table <- function(dt, top_n = 10) {

	# Calculate substance frequencies
	substance_freq <- dt[, .(
		N = .N,
		`%` = round(.N / nrow(dt) * 100, 1),
		`Mean Days` = round(mean(days_to_event), 1),
		`Median Days` = round(median(days_to_event), 1)
	), by = substance]

	setorder(substance_freq, -N)
	substance_freq <- substance_freq[1:min(top_n, nrow(substance_freq))]

	reactable(
		substance_freq,
		columns = list(
			substance = colDef(name = "Substance", minWidth = 180),
			N = colDef(name = "Cases", format = colFormat(separators = TRUE)),
			`%` = colDef(name = "%", format = colFormat(digits = 1)),
			`Mean Days` = colDef(name = "Mean Days", format = colFormat(digits = 1)),
			`Median Days` = colDef(name = "Median Days", format = colFormat(digits = 1))
		),
		highlight = TRUE,
		bordered = TRUE,
		compact = TRUE,
		defaultPageSize = top_n,
		pagination = FALSE,
		style = list(fontSize = "0.85em")
	)
}