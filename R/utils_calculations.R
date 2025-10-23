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

create_pp_groups <- function(dt) {
  dt_copy <- copy(dt)

  # Create more intuitive breaks that still roughly follow log scale
  breaks <- c(2, 4, 6, 8, 11, max(dt_copy$pp, na.rm = TRUE))

  dt_copy[, pp_group := findInterval(pp, breaks, rightmost.closed = TRUE)]

  # Create labels
  labels <- sprintf("[%g-%g)", breaks[-length(breaks)], breaks[-1])
  labels[length(labels)] <- sub("\\)", "]", labels[length(labels)])

  dt_copy[, pp_group := factor(labels[pp_group],
                               levels = labels,  # This preserves the order
                               ordered = TRUE)]

  setkey(dt_copy, patid)
  return(dt_copy)
}

calculate_age_standardized_medians <- function(dt,
                                               value_col = "pp",
                                               group_col = "eth_group",
                                               age_col = "age_group",
                                               bootstrap_num = NULL) {

  # Calculate the total population size
  total_n <- dt[, .N]

  # Calculate the age distribution in the total population
  age_dist <- dt[, .(weight = .N/total_n), by = age_col]

  # Calculate medians for each combination of group and age
  group_age_medians <- dt[, .(
    median_value = as.numeric(median(get(value_col), na.rm = TRUE)),
    n = .N
  ), by = c(group_col, age_col)]

  # Merge with age distribution weights
  group_age_medians <- merge(group_age_medians, age_dist, by = age_col)

  # Calculate weighted medians for each group
  standardized_results <- group_age_medians[, .(
    age_standardized_median = sum(median_value * weight),
    total_n = sum(n)
  ), by = group_col]

  # Add crude (non-standardized) medians for comparison
  crude_medians <- dt[, .(
    crude_median = as.numeric(median(get(value_col), na.rm = TRUE))
  ), by = group_col]

  # Merge standardized and crude results
  final_results <- merge(standardized_results, crude_medians, by = group_col)

  # Calculate CIs only if bootstrap_num is provided
  if (!is.null(bootstrap_num)) {
    boot_ci <- function(dt, group) {
      boot_results <- vector("numeric", bootstrap_num)

      for(i in 1:bootstrap_num) {
        # Bootstrap sample within each age group
        boot_sample <- dt[, .(
          boot_median = as.numeric(median(sample(get(value_col),
                                      size = .N,
                                      replace = TRUE),
                               na.rm = TRUE))
        ), by = c(group_col, age_col)]

        # Merge with weights and calculate standardized median
        boot_sample <- merge(boot_sample, age_dist, by = age_col)
        boot_results[i] <- boot_sample[get(group_col) == group,
                                       sum(boot_median * weight)]
      }

      # Calculate 95% CI
      quantile(boot_results, probs = c(0.025, 0.975))
    }

    # Add confidence intervals for each group
    final_results[, c("ci_lower", "ci_upper") := {
      ci <- boot_ci(dt, get(group_col))
      .(ci[1], ci[2])
    }, by = group_col]
  }

  # Order results by standardized median
  setorder(final_results, -age_standardized_median)

  return(final_results)
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


# poly_level_category <- function(level) {
# 	dplyr::case_when(
# 		level == 0 ~ "[0]",
# 		level == 1 ~ "[1]",
# 		level >= 2 & level <= 3 ~ "[2-3]",
# 		level >= 4 & level <= 5 ~ "[4-5]",
# 		level >= 6 & level <= 7 ~ "[6-7]",
# 		level >= 8 & level <= 10 ~ "[8-10]",
# 		level >= 11 ~ "[11-29]",
# 		TRUE ~ NA_character_
# 	)
# }

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
	#level_order <- c("[0]", "[1]", "[2-3]", "[4-5]", "[6-7]", "[8-10]", "[11-29]")
	level_order <- c("[0]", "[1]", "[2]", "[3]", "[4]", "[5-10]", "[11+]")

	# Process medication data by patient
	# patient_meds <- med_data[, {
	# 	# Create sequence number for each medication
	# 	med_seq <- 1:.N
	#
	# 	# Calculate the actual polypharmacy level BEFORE adding each medication
	# 	# This represents the number of medications BEFORE this substance was added
	# 	poly_level_before <- med_seq - 1
	#
	# 	# Categories
	# 	level_cat_before <- poly_level_category(poly_level_before)
	# 	level_cat_after <- poly_level_category(med_seq)
	#
	# 	# Flag transitions
	# 	is_transition <- level_cat_before != level_cat_after
	#
	# 	list(
	# 		patid = patid,
	# 		substance = substance,
	# 		med_seq = med_seq,
	# 		poly_level_before = poly_level_before,
	# 		poly_level_after = med_seq,
	# 		level_cat_before = level_cat_before,
	# 		level_cat_after = level_cat_after,
	# 		is_transition = is_transition,
	# 		days_to_outcome = eventdate-start_date
	# 	)
	# }, by = patid]
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

	# # Alternatively, if you don't want to modify the original data:
	# patient_meds <- copy(med_data)[, `:=`(
	# 	med_seq = seq_len(.N),
	# 	poly_level_before = seq_len(.N) - 1,
	# 	poly_level_after = seq_len(.N)
	# ), by = patid]
	#
	# patient_meds[, `:=`(
	# 	level_cat_before = poly_level_category(poly_level_before),
	# 	level_cat_after = poly_level_category(poly_level_after),
	# 	is_transition = NULL,
	# 	days_to_outcome = eventdate - start_date
	# )]
	#
	# patient_meds[, is_transition := level_cat_before != level_cat_after]


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