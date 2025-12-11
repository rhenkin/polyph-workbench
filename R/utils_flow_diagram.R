#' Generate CONSORT-style flow diagram text from study metadata
#'
#' @param metadata Study metadata list containing flow_diagram data
#' @return Character string with formatted flow diagram text
#' @export
generate_flow_diagram_text <- function(metadata) {

	if (is.null(metadata$flow_diagram)) {
		stop("Flow diagram data not found in metadata. Was this study created with an older version?")
	}

	fd <- metadata$flow_diagram

	# Build the text output
	lines <- c()

	# Header
	lines <- c(lines, "========================================")
	lines <- c(lines, "FLOW DIAGRAM")
	lines <- c(lines, paste("Study:", metadata$study_name))
	lines <- c(lines, paste("Generated:", fd$generated_date))
	lines <- c(lines, "========================================")
	lines <- c(lines, "")

	# Database starting point
	if (!is.null(fd$n_total_database) && !is.na(fd$n_total_database)) {
		lines <- c(lines, sprintf("Total database population: n=%s",
															format(fd$n_total_database, big.mark = ",")))
		lines <- c(lines, "")
	}

	# Initial cohort
	if (!is.null(fd$n_initial_cohort) && !is.na(fd$n_initial_cohort)) {
		lines <- c(lines, sprintf("Initial cohort with outcome and meeting eligibility criteria¹: n=%s",
															format(fd$n_initial_cohort, big.mark = ",")))

		if (!is.na(fd$n_excluded_initial)) {
			lines <- c(lines, sprintf("  └─ Excluded: n=%s",
																format(fd$n_excluded_initial, big.mark = ",")))
		}
		lines <- c(lines, "")
	}

	# Cases section
	lines <- c(lines, "CASES:")
	lines <- c(lines, "─────────")

	pred_window_text <- if (!is.null(fd$pred_window)) {
		sprintf("within %d-day prediction window", fd$pred_window)
	} else {
		"within prediction window"
	}

	lines <- c(lines, sprintf("Outcome prescriptions %s: n=%s",
														pred_window_text,
														format(fd$n_cases_after_pred_window, big.mark = ",")))

	if (!is.na(fd$n_excluded_pred_window)) {
		lines <- c(lines, sprintf("  └─ Excluded (outside prediction window): n=%s",
															format(fd$n_excluded_pred_window, big.mark = ",")))
	}

	lines <- c(lines, "")
	lines <- c(lines, sprintf("Final cases: n=%s",
														format(metadata$n_cases, big.mark = ",")))
	lines <- c(lines, "")
	lines <- c(lines, "")

	# Controls section
	lines <- c(lines, "CONTROLS:")
	lines <- c(lines, "─────────")

	if (!is.null(fd$n_eligible_control_pool) && !is.na(fd$n_eligible_control_pool)) {
		lines <- c(lines, sprintf("Eligible control pool² (patient-strata combinations): n=%s",
															format(fd$n_eligible_control_pool, big.mark = ",")))
		lines <- c(lines, "")
	}

	match_ratio_text <- if (!is.null(fd$match_ratio)) {
		sprintf("(ratio %d:1)", fd$match_ratio)
	} else {
		""
	}

	lines <- c(lines, sprintf("Sampled for matching %s: n=%s",
														match_ratio_text,
														format(fd$n_controls_sampled_total, big.mark = ",")))
	lines <- c(lines, "")
	lines <- c(lines, sprintf("Final controls: n=%s (unique patients: n=%s)",
														format(metadata$n_controls_total, big.mark = ","),
														format(metadata$n_controls_unique, big.mark = ",")))
	lines <- c(lines, "")
	lines <- c(lines, "")

	# Final cohort summary
	lines <- c(lines, "FINAL MATCHED COHORT:")
	lines <- c(lines, "─────────────────────")
	lines <- c(lines, sprintf("Cases: n=%s", format(metadata$n_cases, big.mark = ",")))
	lines <- c(lines, sprintf("Controls (total): n=%s", format(metadata$n_controls_total, big.mark = ",")))
	lines <- c(lines, sprintf("Controls (unique): n=%s", format(metadata$n_controls_unique, big.mark = ",")))
	lines <- c(lines, sprintf("Matching ratio (total): %.2f:1", metadata$matching_ratio_total))
	lines <- c(lines, sprintf("Matching ratio (unique): %.2f:1", metadata$matching_ratio_unique))
	lines <- c(lines, "")
	lines <- c(lines, "")

	# Footnotes
	lines <- c(lines, "NOTES:")
	lines <- c(lines, "──────")
	lines <- c(lines, sprintf("¹ Eligibility criteria: %s", fd$eligibility_criteria))
	lines <- c(lines, "² Non-cases with matching strata and meeting same eligibility criteria")
	lines <- c(lines, "")
	lines <- c(lines, sprintf("Matching method: %s", metadata$matching_method))
	lines <- c(lines, sprintf("Stratification variable: %s", metadata$strata_variable))
	lines <- c(lines, sprintf("Study created: %s", metadata$created_date))

	# Combine all lines
	flow_text <- paste(lines, collapse = "\n")

	return(flow_text)
}


#' Preview flow diagram in console (for testing)
#'
#' @param metadata Study metadata
#' @export
preview_flow_diagram <- function(metadata) {
	flow_text <- generate_flow_diagram_text(metadata)
	cat(flow_text)
	invisible(flow_text)
}