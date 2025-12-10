#' Reactable Utility Functions for CCA Module
#'
#' These functions provide reusable components for rendering prevalence tables
#' with odds ratios, confidence intervals, and statistical significance.

#' Get standard OR column definitions for reactable
#'
#' @return Named list of colDef objects for OR-related columns
get_or_reactable_columns <- function() {
	list(
		OR_formatted = colDef(name = "OR (95% CI)", minWidth = 140),
		OR = colDef(show = FALSE),
		OR_CI_lower = colDef(show = FALSE),
		OR_CI_upper = colDef(show = FALSE),
		p_value = colDef(show = FALSE),
		p_adj = colDef(show = FALSE)
	)
}

#' Create significance cell renderer for item column
#'
#' Adds asterisk to items with significant adjusted p-values
#'
#' @param result_table data.table with p_adj column (required for closure)
#' @return Function suitable for reactable cell parameter
create_significance_cell <- function(result_table) {
	function(value, index) {
		p_adj_val <- result_table[index, p_adj]
		if (!is.na(p_adj_val) && p_adj_val < 0.05) {
			paste0(value, "*")
		} else {
			value
		}
	}
}

#' Create p-value details panel for reactable
#'
#' @param index Row index
#' @param result_table data.table with p_value and p_adj columns
#' @return htmltools div with formatted p-value information, or NULL if no p-values
create_pvalue_details <- function(index, result_table) {
	p_val <- result_table[index, p_value]
	p_adj_val <- result_table[index, p_adj]
	if (is.na(p_val) & is.na(p_adj_val)) return(NULL)

	htmltools::div(
		style = "padding: 16px",
		htmltools::tags$b("Statistical Testing:"),
		htmltools::tags$div(
			style = "margin-top: 8px",
			sprintf("Raw p-value: %.4f", p_val)
		),
		htmltools::tags$div(
			sprintf("Adjusted p-value: %.4f", p_adj_val)
		),
		htmltools::tags$div(
			style = "margin-top: 8px; font-style: italic; color: #666;",
			if (p_adj_val < 0.05) "Statistically significant (p < 0.05)" else "Not significant"
		)
	)
}

#' Get standard reactable configuration
#'
#' @return Named list of standard reactable parameters
get_standard_reactable_config <- function() {
	list(
		showPageInfo = FALSE,
		searchable = TRUE,
		showPageSizeOptions = TRUE,
		defaultPageSize = 15,
		compact = TRUE
	)
}

#' Get standard column name mapping for prevalence tables
#'
#' @param item_name Display name for the item column (e.g., "LTC", "Substance")
#' @return Named character vector for column renaming
get_prevalence_column_names <- function(item_name = "Item") {
	c(item_name, "Cases (%)", "Controls (%)", "OR", "OR_CI_lower", "OR_CI_upper")
}