# R/utils_temporal_helpers.R

# R/utils_temporal_grouping.R

#' Create adaptive day grouping based on risk window
#' @param pred_window numeric, the prediction window in days
#' @return list with breaks, labels, and bin_width
create_day_groups <- function(pred_window) {

	if (pred_window <= 7) {
		# Daily bins for short windows
		breaks <- 0:pred_window
		labels <- as.character(1:pred_window)
		bin_width <- 1

	} else if (pred_window <= 14) {
		# Every 2 days
		breaks <- seq(0, pred_window, by = 2)
		labels <- paste0(breaks[-length(breaks)] + 1, "-", breaks[-1])
		bin_width <- 2

	} else if (pred_window <= 30) {
		# Weekly bins
		breaks <- c(0, 1, 7, 14, 21, min(30, pred_window))
		if (pred_window > 30) breaks <- c(breaks, pred_window)
		labels <- c("Day 1", "2-7 days", "8-14 days", "15-21 days", "22-30 days")
		if (pred_window > 30) labels <- c(labels, paste0("31-", pred_window, " days"))
		bin_width <- "variable"

	} else if (pred_window <= 90) {
		# Bi-weekly then monthly
		breaks <- c(0, 7, 14, 30, 60, min(90, pred_window))
		if (pred_window > 90) breaks <- c(breaks, pred_window)
		labels <- c("1-7 days", "8-14 days", "15-30 days", "31-60 days", "61-90 days")
		if (pred_window > 90) labels <- c(labels, paste0("91-", pred_window, " days"))
		bin_width <- "variable"

	} else {
		# Monthly bins for very long windows
		breaks <- c(0, 7, 30, 60, 90, 120, 180, min(365, pred_window))
		if (pred_window > 365) breaks <- c(breaks, pred_window)
		labels <- c("1-7 days", "8-30 days", "31-60 days", "61-90 days",
								"91-120 days", "121-180 days", "181-365 days")
		if (pred_window > 365) labels <- c(labels, paste0("365+ days"))
		bin_width <- "variable"
	}

	return(list(
		breaks = breaks,
		labels = labels,
		bin_width = bin_width,
		pred_window = pred_window
	))
}

#' Add day group variable to prescription data
#' @param dt data.table with eventdate/index_date and outcome_date columns
#' @param pred_window numeric, the prediction window
#' @return data.table with days_to_event and day_group columns
#' @export
add_temporal_groups <- function(dt, pred_window) {
	# Calculate days to event
	# Handle both CCA data formats: eventdate (from saved studies) or index_date (from filter_by_index_date)
	if ("eventdate" %in% names(dt)) {
		dt[, days_to_event := as.numeric(outcome_date - eventdate)]
	} else if ("index_date" %in% names(dt)) {
		dt[, days_to_event := as.numeric(outcome_date - index_date)]
	} else {
		stop("Data must have either 'eventdate' or 'index_date' column along with 'outcome_date'. ",
				 "Available columns: ", paste(names(dt), collapse = ", "))
	}

	# Get grouping structure
	groups <- create_day_groups(pred_window)

	# Add grouped variable
	dt[, day_group := cut(days_to_event,
												breaks = groups$breaks,
												labels = groups$labels,
												include.lowest = TRUE,
												right = TRUE)]

	return(dt)
}
#' Filter data by stratification variable
#'
#' Uses the encoded format "column#value" from virtualSelect
#'
#' @param data data.table to filter
#' @param strat_var Stratification string in format "column#value"
#' @return Filtered data.table
#' @export
filter_by_stratification <- function(data, strat_var) {

	if (is.null(strat_var) || strat_var == "") {
		return(data)
	}

	# Parse the encoded stratification variable
	parts <- strsplit(strat_var, "#")[[1]]

	if (length(parts) != 2) {
		warning("Invalid stratification variable format: ", strat_var)
		return(data)
	}

	column_name <- parts[1]
	filter_value <- parts[2]

	# Check column exists
	if (!column_name %in% names(data)) {
		warning(sprintf("Column '%s' not found in data", column_name))
		return(data)
	}

	# Apply filter
	data[get(column_name) == filter_value]
}

#' Filter patients by LTCs (require ALL selected LTCs)
#'
#' @param ltc_data data.table with patid and term columns
#' @param selected_ltcs Character vector of LTC terms
#' @return Vector of patient IDs
#' @export
filter_patients_by_ltcs <- function(ltc_data, selected_ltcs) {

	if (is.null(selected_ltcs) || length(selected_ltcs) == 0) {
		return(unique(ltc_data$patid))
	}

	# Count how many of the selected LTCs each patient has
	ltc_counts <- ltc_data[term %in% selected_ltcs,
												 .(n_selected = uniqueN(term)),
												 by = patid]

	# Keep only patients with ALL selected LTCs
	n_required <- length(selected_ltcs)
	valid_patids <- ltc_counts[n_selected == n_required, patid]

	return(valid_patids)
}
#'
#' #' Create dropdown choices for LTCs
#' #'
#' #' @param ltc_data data.table with term column
#' #' @param top_n Number of most common LTCs to include
#' #' @return Character vector of LTC choices
#' #' @export
#' create_ltc_dropdown_choices <- function(ltc_data, top_n = 50) {
#' 	ltc_counts <- ltc_data[, .N, by = term][order(-N)]
#' 	ltc_counts[1:min(top_n, nrow(ltc_counts)), term]
#' }
#'
#' #' Create dropdown choices for BNF substances
#' #'
#' #' @param prescription_data data.table with substance column
#' #' @param bnf_level Current BNF aggregation level
#' #' @param top_n Number of most common substances to include
#' #' @return Character vector of substance choices
#' #' @export
#' create_bnf_dropdown_choices <- function(prescription_data, bnf_level, top_n = 100) {
#' 	# Get the appropriate column based on BNF level
#' 	substance_col <- if (bnf_level == "BNF_Chemical_Substance") {
#' 		"substance"
#' 	} else {
#' 		bnf_level
#' 	}
#'
#' 	if (!substance_col %in% names(prescription_data)) {
#' 		substance_col <- "substance"  # Fallback
#' 	}
#'
#' 	substance_counts <- prescription_data[, .N, by = get(substance_col)][order(-N)]
#' 	setnames(substance_counts, "get", "substance")
#' 	substance_counts[1:min(top_n, nrow(substance_counts)), substance]
#' }

#' Calculate filtered patient counts for display
#'
#' @param data data.table with group column
#' @param strat_var Stratification variable
#' @param filter_patids List of patient ID vectors from different filters
#' @return List with case and control counts
#' @export
get_filtered_patient_counts_temporal <- function(data, strat_var = NULL,
																				filter_patids_1 = NULL,
																				filter_patids_2 = NULL,
																				filter_patids_3 = NULL) {

	dt <- copy(data)

	# Apply stratification
	if (!is.null(strat_var) && strat_var != "" && !is.null(dt)) {
		dt <- filter_by_stratification(dt, strat_var)
	}

	# Apply patient ID filters (intersection)
	all_filters <- list(filter_patids_1, filter_patids_2, filter_patids_3)
	all_filters <- all_filters[!sapply(all_filters, is.null)]

	if (length(all_filters) > 0) {
		# Intersect all filter lists
		valid_patids <- Reduce(intersect, all_filters)
		dt <- dt[patid %in% valid_patids]
	}

	# Count unique patients by group
	counts <- dt[, .(n = uniqueN(patid)), by = group]

	list(
		n_cases = counts[group == "case", n],
		n_controls = counts[group == "control", n],
		n_total = uniqueN(dt$patid)
	)
}

#' Render patient count display UI
#'
#' @param counts List from get_filtered_patient_counts
#' @return Shiny UI element
#' @export
render_patient_count_ui_temp <- function(counts) {

	if (is.null(counts$n_cases)) counts$n_cases <- 0
	if (is.null(counts$n_controls)) counts$n_controls <- 0

	div(
		style = "padding: 10px; background-color: #f8f9fa; border-radius: 5px; margin: 10px 0;",
		div(
			style = "display: flex; justify-content: space-around; text-align: center;",
			div(
				div(style = "font-size: 0.9em; color: #666;", "Cases"),
				div(style = "font-size: 1.4em; font-weight: 600; color: #dc3545;",
						prettyNum(counts$n_cases, big.mark = ","))
			),
			div(
				div(style = "font-size: 0.9em; color: #666;", "Controls"),
				div(style = "font-size: 1.4em; font-weight: 600; color: #28a745;",
						prettyNum(counts$n_controls, big.mark = ","))
			),
			div(
				div(style = "font-size: 0.9em; color: #666;", "Total"),
				div(style = "font-size: 1.4em; font-weight: 600; color: #007bff;",
						prettyNum(counts$n_total, big.mark = ","))
			)
		)
	)
}