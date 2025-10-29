#' Calculate frequency of items (generic version)
#'
#' @param data data.table containing data with columns: patid and the item column
#' @param item_col character string specifying the column name for items (e.g., "term", "substance")
#' @return data.table with columns: item_col, N, pct_total (sorted by frequency descending)
calculate_item_frequency <- function(data, item_col) {
  # Get unique patient-item combinations
  cols_to_select <- c("patid", item_col)
  unique_patid <- unique(data[, ..cols_to_select])

  # Calculate frequency using proper by syntax with character variable
  item_freq <- unique_patid[, .N, by = c(item_col)]
  total_patids <- uniqueN(data$patid)
  item_freq[, pct_total := N / total_patids]
  setorderv(item_freq, "N", order = -1)

  return(item_freq)
}

#' Calculate prevalence ratio with confidence intervals (generic version)
#'
#' @param data data.table containing the full dataset with patid and item columns
#' @param item_col character string specifying the column name for items (e.g., "term", "substance")
#' @param selected_patids character vector of patient IDs in the "with" group
#' @param min_prevalence numeric threshold for including results (default: 0.01)
#' @param duration_col optional character string specifying duration column name for duration statistics
#' @return data.table with prevalence ratios, confidence intervals, and optional duration statistics
calculate_prevalence_ratio <- function(data, item_col, selected_patids,
                                       min_prevalence = 0.01, duration_col = NULL) {
  # Patients WITH the selected condition
  with_freq <- data[
    patid %in% selected_patids,
    {
      result <- list(
        N_with = uniqueN(patid),
        Prevalence = round(100 * (uniqueN(patid) / length(selected_patids)), digits = 2)
      )
      # Add duration statistics if duration column is provided
      if (!is.null(duration_col) && duration_col %in% names(.SD)) {
        result[[paste0("Median Duration (years)")]] <- round(median(get(duration_col) / 365.2), digits = 2)
        result[["IQR (Q1-Q3)"]] <- paste0(
          "(",
          round(quantile(get(duration_col) / 365.2, 0.25, na.rm = TRUE), 2), " - ",
          round(quantile(get(duration_col) / 365.2, 0.75, na.rm = TRUE), 2), ")"
        )
      }
      result
    },
    by = item_col
  ]

  # Patients WITHOUT the selected condition
  unselected_patids <- data[!patid %in% selected_patids, uniqueN(patid)]
  without_freq <- data[
    !patid %in% selected_patids,
    {
      result <- list(
        N_without = uniqueN(patid),
        Prevalence_Unselected = round(100 * (uniqueN(patid) / unselected_patids), digits = 2)
      )
      # Add duration statistics for unselected group if duration column is provided
      if (!is.null(duration_col) && duration_col %in% names(.SD)) {
        result[["Median Duration unselected (years)"]] <- round(median(get(duration_col) / 365.2), digits = 2)
        result[["IQR unsel. (Q1-Q3)"]] <- paste0(
          "(",
          round(quantile(get(duration_col) / 365.2, 0.25, na.rm = TRUE), 2), " - ",
          round(quantile(get(duration_col) / 365.2, 0.75, na.rm = TRUE), 2), ")"
        )
      }
      result
    },
    by = item_col
  ]

  # Merge results
  result <- merge(with_freq, without_freq, by = item_col)

  # Filter by minimum prevalence
  result <- result[Prevalence >= min_prevalence]

  # Calculate prevalence ratio
  result[, `:=`(
    total_with = length(selected_patids),
    total_without = unselected_patids,
    Prevalence_Ratio = round(Prevalence / Prevalence_Unselected, digits = 2)
  )]

  result[is.infinite(Prevalence_Ratio), Prevalence_Ratio := 0]

  # Calculate 95% CI for the ratio using log method
  result[, `:=`(
    p1 = N_with / total_with,
    p2 = N_without / total_without
  )]

  result[, `:=`(
    log_ratio = log(Prevalence_Ratio),
    se_log_ratio = sqrt(
      (1 / N_with) - (1 / total_with) +
      (1 / N_without) - (1 / total_without)
    )
  )]

  result[, `:=`(
    CI_lower = round(exp(log_ratio - 1.96 * se_log_ratio), digits = 2),
    CI_upper = round(exp(log_ratio + 1.96 * se_log_ratio), digits = 2)
  )]

  result[, CI_95 := paste0("(", CI_lower, " - ", CI_upper, ")")]

  # Mark significant results (CI doesn't include 1.0)
  result[(CI_lower > 1.0 | CI_upper < 1.0), (item_col) := paste0(get(item_col), "*")]

  # Clean up intermediate columns
  result[, c("p1", "p2", "log_ratio", "se_log_ratio", "CI_lower", "CI_upper",
             "total_with", "total_without") := NULL]

  # Sort by prevalence ratio
  setorderv(result, "Prevalence_Ratio", order = -1)

  return(result)
}

#' Filter data by age threshold
#'
#' @param data data.table containing age column
#' @param max_age numeric maximum age in years
#' @param age_col character string specifying the age column name (default: "outcome_age")
#' @return filtered data.table
filter_by_age <- function(data, max_age, age_col = "outcome_age") {
  data[get(age_col) <= max_age * 365.25]
}

#' Format patient count with thousands separator
#'
#' @param data data.table to count unique patients from
#' @return character string with formatted count
format_patient_count <- function(data) {
  paste0("Patients within age range: ", prettyNum(uniqueN(data$patid), big.mark = ","))
}

#' Create groups from continuous numeric values
#' @param dt data.table containing value column
#' @param breaks numeric vector of break points
#' @param value_col name of column to group (default: "value")
#' @param group_col name of output grouped column (default: "group")
#' @param label_fmt string format for labels (default: "[%g-%g)")
#' @param right logical, intervals closed on right? (default: TRUE)
#' @param label_suffix optional suffix for labels
#' @return data.table with new group column
create_value_groups <- function(dt,
                                breaks,
                                value_col = "value",
                                group_col = "group",
                                label_fmt = "[%g-%g]",
                                right = TRUE,
                                label_suffix = NULL) {
	dt_copy <- copy(dt)
  if(max(breaks) < max(dt_copy[[value_col]], na.rm=TRUE)) {
    breaks <- c(breaks, max(dt_copy[[value_col]], na.rm=TRUE))
  }
  labels <- sprintf(label_fmt, breaks[-length(breaks)], breaks[-1]-1)
  if(!is.null(label_suffix)) labels <- paste0(labels, label_suffix)
  labels[length(labels)] <- sub("\\)", "]", labels[length(labels)])

  dt_copy[, (group_col) := cut(get(value_col),
                               breaks = breaks,
                               labels = labels,
                               right = right,
                               include.lowest = TRUE,
                               ordered = TRUE)][,.("patid", group_col)]
  setkey(dt_copy, patid)
  return(dt_copy)
}

create_age_groups <- function(dt, n_groups = 5) {
	dt_copy <- copy(dt)

	# Calculate break points
	breaks <- quantile(dt_copy$age_days/365.25, probs = seq(0, 1, length.out = n_groups + 1))

	# Create labels showing the ranges
	labels <- paste(round(breaks[-length(breaks)]), "to", round(breaks[-1]), "years")

	# Add group column using the range labels
	dt_copy[, age_group := cut(age_days/365.25,
														 breaks = breaks,
														 labels = labels,
														 include.lowest = TRUE,
														 ordered_result = TRUE)]
	setkey(dt_copy, patid)
	return(dt_copy)
}

# Helper function to calculate category totals
calculate_category_totals <- function(df, demog_var) {
	df[, .(Total = uniqueN(patid)), by = demog_var]
}

# Helper function to calculate demographic statistics
calculate_demographic_stats <- function(df, demog_var, frequency_data, condition_col) {
	# Calculate counts by demographic category and substance
	df_stats <- df[, .(
		N_category = uniqueN(patid)
	), by = c(demog_var, condition_col)]

	# Merge with substance frequency data
	merge_cols <- c(condition_col, "N")
	df_stats <- merge(df_stats, frequency_data[, ..merge_cols], by = condition_col)

	# Calculate percentages
	df_stats[, pct := signif(N_category/N, digit = 2)]

	return(df_stats)
}

#' Generic demographic frequency calculation with statistical tests
#'
#' @param data data.table containing the full dataset (must include patid and item_col)
#' @param item_col character string specifying the item column name (e.g., "term", "substance")
#' @param demog_var character string specifying the demographic variable to analyze
#' @param frequency_data data.table containing overall frequencies (from calculate_item_frequency)
#' @param min_prevalence numeric threshold for including items (default: 0.005)
#' @return formatted data.table with demographic breakdowns and chi-square test results
calculate_demographic_frequency <- function(data, item_col, demog_var, frequency_data,
                                           min_prevalence = 0.005) {
  df <- copy(data)

  # Filter to common items only
  df <- df[get(item_col) %in% frequency_data[pct_total > min_prevalence, get(item_col)]]

  # Calculate statistics
  cat_totals <- calculate_category_totals(df, demog_var)
  df_stats <- calculate_demographic_stats(df, demog_var, frequency_data, item_col)

  # Perform chi-square tests
  chisq_tests <- perform_chisq_tests(df_stats, cat_totals, demog_var, item_col)

  # Format and return results
  return(format_results(df_stats, chisq_tests, cat_totals, demog_var,
                       frequency_data, item_col))
}

poly_level_category <- function(level) {
	dplyr::case_when(
		level == 0 ~ "[0]",
		level == 1 ~ "[1]",
		level == 2 ~ "[2]",
		level == 3 ~ "[3]",
		level == 4 ~ "[4]",
		level >= 5 & level <= 10 ~ "[5-10]",
		level >= 11 ~ "[11+]",
		TRUE ~ NA_character_
	)
}

analyze_medications <- function(input_med_data, mode = "transition") {
	# Define the correct order for polypharmacy categories
	level_order <- c("[0]", "[1]", "[2]", "[3]", "[4]", "[5-10]", "[11+]")

	med_data <- data.table::copy(input_med_data)
	setkey(med_data, patid, eventdate)

	# Allocate the result columns directly
	med_data[, med_seq := seq_len(.N), by = patid]
	med_data[, poly_level_before := med_seq - 1]
	med_data[, poly_level_after := med_seq]

	# Apply the categorization function once to the entire columns
	med_data[, level_cat_before := poly_level_category(poly_level_before)]
	med_data[, level_cat_after := poly_level_category(poly_level_after)]

	# Calculate transitions
	med_data[, is_transition := level_cat_before != level_cat_after]
	med_data[, days_to_outcome := eventdate - start_date]

	# Select only the columns we need for the result
	med_data <- med_data[, .(
		patid, substance, med_seq,
		poly_level_before, poly_level_after,
		level_cat_before, level_cat_after,
		is_transition, days_to_outcome
	)]
	patient_meds <- med_data

	# Filter for frequent medications
	frequent_meds <- patient_meds[, .(
		count = .N,
		prevalence = .N / nrow(patient_meds) * 100
	), by = substance][prevalence >= 1, substance]

	# Choose analysis based on mode
	if(mode == "transition") {
		substance_data <- patient_meds[substance %in% frequent_meds, .(
			count = .N
		), by = .(substance, level_cat_before, level_cat_after)]

		# Order for transitions
		result <- dcast(
			substance_data,
			substance ~ paste(level_cat_before, "→", level_cat_after),
			value.var = "count",
			fill = 0
		)
	} else if(mode == "level") {
		# For level mode, we're interested in the polypharmacy level BEFORE
		# this substance was added (the context in which it was prescribed)
		substance_data <- patient_meds[substance %in% frequent_meds, .(
			count = .N
		), by = .(substance, level_cat_before)]

		# Order categories correctly
		substance_data[, level_cat_before := factor(level_cat_before, levels = level_order)]

		result <- dcast(
			substance_data,
			substance ~ level_cat_before,
			value.var = "count",
			fill = 0
		)
	}

	# Calculate total
	level_cols <- intersect(names(result), c(level_order, paste(level_order, "→", level_order, sep="")))
	result[, total := rowSums(.SD), .SDcols = level_cols]

	# Calculate percentages
	result[, (level_cols) := lapply(.SD, function(x) round((x / total) * 100, 1)),
				 .SDcols = level_cols]

	# Return ordered results
	return(result[order(-total)])
}