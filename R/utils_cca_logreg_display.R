# R/utils_cca_logreg_display.R
# Display and formatting utilities for logistic regression results

#' Render logistic regression results table
#'
#' Single unified function to render results for all model types
#'
#' @param results data.table with model results
#' @param model_type character string indicating model type
#' @return reactable object
render_logreg_results_table <- function(results, model_type) {
	if (is.null(results) || nrow(results) == 0) {
		return(NULL)
	}

	# Format OR columns based on model type
	results <- format_or_columns(results, model_type)

	# Get model-specific configuration
	config <- get_table_config(model_type)

	# Create reactable with configuration
	reactable(
		results[, config$columns, with = FALSE],
		columns = config$column_defs,
		defaultPageSize = 20,
		searchable = TRUE,
		showPageSizeOptions = TRUE,
		compact = TRUE,
		defaultSorted = config$default_sort
	)
}

#' Format OR columns for display
#'
#' @param results data.table with model results
#' @param model_type character string indicating model type
#' @return data.table with formatted OR columns
format_or_columns <- function(results, model_type) {
	results <- copy(results)

	switch(model_type,
				 background_main = {
				 	results[, OR_formatted := sprintf("%.2f (%.2f-%.2f)", OR, CI_lower, CI_upper)]
				 	# Calculate percentages for display
				 	if ("n_cases" %in% names(results) && "n_controls" %in% names(results) &&
				 			"total_cases" %in% names(results) && "total_controls" %in% names(results)) {
				 		results[, pct_cases := round(100 * n_cases / total_cases, 2)]
				 		results[, pct_controls := round(100 * n_controls / total_controls, 2)]
				 	}
				 },
				 background_pairwise = {
				 	results[, med1_OR_formatted := sprintf("%.3f (%.3f-%.3f)", med1_OR, med1_CI_lower, med1_CI_upper)]
				 	results[, med2_OR_formatted := sprintf("%.3f (%.3f-%.3f)", med2_OR, med2_CI_lower, med2_CI_upper)]
				 	results[, interaction_OR_formatted := sprintf("%.3f (%.3f-%.3f)", interaction_OR, interaction_CI_lower, interaction_CI_upper)]
				 	results[, combined_OR_formatted := sprintf("%.3f (%.3f-%.3f)", combined_OR, combined_CI_lower, combined_CI_upper)]
				 },
				 recent_pp = {
				 	# Already formatted in the modeling function
				 	NULL
				 },

				 recent_main = {
				 	# Already formatted in the modeling function
				 	NULL
				 },
				 recent_background = {
				 	results[, recent_OR_formatted := sprintf("%.3f (%.3f-%.3f)", recent_OR, recent_CI_lower, recent_CI_upper)]
				 	results[, background_OR_formatted := sprintf("%.3f (%.3f-%.3f)", background_OR, background_CI_lower, background_CI_upper)]
				 	results[, interaction_OR_formatted := sprintf("%.3f (%.3f-%.3f)", interaction_OR, interaction_CI_lower, interaction_CI_upper)]
				 	results[, combined_OR_formatted := sprintf("%.3f (%.3f-%.3f)", combined_OR, combined_CI_lower, combined_CI_upper)]
				 }
	)

	return(results)
}

#' Get table configuration for model type
#'
#' @param model_type character string indicating model type
#' @return list with columns and column_defs
get_table_config <- function(model_type) {
	switch(model_type,
				 background_main = list(
				 	columns = c("medication", "OR_formatted", "p_value", "pct_cases", "pct_controls"),
				 	column_defs = list(
				 		medication = colDef(name = "Medication", minWidth = 200),
				 		OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
				 		p_value = colDef(name = "P-value", format = colFormat(digits = 4)),
				 		pct_cases = colDef(name = "Cases %", format = colFormat(digits = 2)),
				 		pct_controls = colDef(name = "Controls %", format = colFormat(digits = 2))
				 	),
				 	default_sort = "medication"
				 ),

				 background_pairwise = list(
				 	columns = c("med1", "med2", "med1_OR_formatted", "med2_OR_formatted",
				 							"interaction_OR_formatted", "combined_OR_formatted",
				 							"interaction_p", "pct_cases_both", "pct_controls_both"),
				 	column_defs = list(
				 		med1 = colDef(name = "Medication 1", minWidth = 150),
				 		med2 = colDef(name = "Medication 2", minWidth = 150),
				 		med1_OR_formatted = colDef(name = "Med1 OR (95% CI)", minWidth = 140),
				 		med2_OR_formatted = colDef(name = "Med2 OR (95% CI)", minWidth = 140),
				 		interaction_OR_formatted = colDef(name = "Interaction OR (95% CI)", minWidth = 160),
				 		combined_OR_formatted = colDef(name = "Combined OR (95% CI)", minWidth = 140),
				 		interaction_p = colDef(name = "Interaction P", format = colFormat(digits = 4)),
				 		pct_cases_both = colDef(name = "Cases Both %", format = colFormat(digits = 2)),
				 		pct_controls_both = colDef(name = "Controls Both %", format = colFormat(digits = 2))
				 	),
				 	default_sort = NULL
				 ),

				 recent_pp = list(
				 	columns = c("medication", "pp_level", "main_OR_formatted", "interaction_OR_formatted",
				 							"combined_OR_formatted", "interaction_p", "case_prev", "control_prev"),
				 	column_defs = list(
				 		medication = colDef(name = "Medication", minWidth = 150),
				 		pp_level = colDef(name = "PP Group", minWidth = 100),
				 		main_OR_formatted = colDef(name = "Main OR (95% CI)", minWidth = 140),
				 		interaction_OR_formatted = colDef(name = "Interaction OR (95% CI)", minWidth = 160),
				 		combined_OR_formatted = colDef(name = "Combined OR (95% CI)", minWidth = 140),
				 		interaction_p = colDef(name = "Interaction P", format = colFormat(digits = 4)),
				 		case_prev = colDef(name = "Case %", format = colFormat(digits = 2)),
				 		control_prev = colDef(name = "Control %", format = colFormat(digits = 2))
				 	),
				 	default_sort = NULL
				 ),

				 recent_main = list(
				 	columns = c("medication", "OR_formatted", "p_value", "pct_cases", "pct_controls"),
				 	column_defs = list(
				 		medication = colDef(name = "Medication", minWidth = 200),
				 		OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
				 		p_value = colDef(name = "P-value", format = colFormat(digits = 4)),
				 		pct_cases = colDef(name = "Cases %", format = colFormat(digits = 2)),
				 		pct_controls = colDef(name = "Controls %", format = colFormat(digits = 2))
				 	),
				 	default_sort = "medication"
				 ),

				 recent_background = list(
				 	columns = c("recent_med", "background_med", "recent_OR_formatted", "background_OR_formatted",
				 							"interaction_OR_formatted", "combined_OR_formatted",
				 							"interaction_p", "pct_cases_both", "pct_controls_both"),
				 	column_defs = list(
				 		recent_med = colDef(name = "Recent Prescription", minWidth = 150),
				 		background_med = colDef(name = "Background Medication", minWidth = 150),
				 		recent_OR_formatted = colDef(name = "Recent OR (95% CI)", minWidth = 140),
				 		background_OR_formatted = colDef(name = "Background OR (95% CI)", minWidth = 140),
				 		interaction_OR_formatted = colDef(name = "Interaction OR (95% CI)", minWidth = 160),
				 		combined_OR_formatted = colDef(name = "Combined OR (95% CI)", minWidth = 140),
				 		interaction_p = colDef(name = "Interaction P", format = colFormat(digits = 4)),
				 		pct_cases_both = colDef(name = "Cases Both %", format = colFormat(digits = 2)),
				 		pct_controls_both = colDef(name = "Controls Both %", format = colFormat(digits = 2))
				 	),
				 	default_sort = NULL
				 )
	)
}