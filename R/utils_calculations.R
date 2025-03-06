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
calculate_demographic_stats <- function(df, demog_var, subst_frequency) {
	# Calculate counts by demographic category and substance
	df_stats <- df[, .(
		N_category = uniqueN(patid)
	), by = c(demog_var, "substance")]

	# Merge with substance frequency data
	df_stats <- merge(
		df_stats,
		subst_frequency[, .(substance, N)],
		by = "substance"
	)

	# Calculate percentages
	df_stats[, pct := signif(N_category/N, digit = 2)]

	return(df_stats)
}