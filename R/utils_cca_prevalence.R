# UI helper
create_ltc_dropdown_choices <- function(ltcs_data) {
	with(ltc_chapters[ltc %in% unique(ltcs_data$term)],
			 split(ltc, body_system))
}

create_bnf_dropdown_choices <- function(prescriptions_data, bnf_level) {
	unique_substances <- unique(prescriptions_data$substance)

	if (bnf_level == "BNF_Chapter") {
		return(sort(unique_substances))
	} else if (bnf_level == "BNF_Section") {
		lookup_data <- bnf_lookup[BNF_Section %in% unique_substances,
															.(BNF_Section, BNF_Chapter)] |> unique()
		lookup_data <- lookup_data[order(BNF_Chapter, BNF_Section)]
		return(with(lookup_data, split(BNF_Section, BNF_Chapter)))
	} else if (bnf_level == "BNF_Paragraph") {
		lookup_data <- bnf_lookup[BNF_Paragraph %in% unique_substances,
															.(BNF_Paragraph, BNF_Section)] |> unique()
		lookup_data <- lookup_data[order(BNF_Section, BNF_Paragraph)]
		return(with(lookup_data, split(BNF_Paragraph, BNF_Section)))
	} else {  # BNF_Chemical_Substance
		lookup_data <- bnf_lookup[BNF_Chemical_Substance %in% unique_substances,
															.(BNF_Chemical_Substance, BNF_Paragraph)] |> unique()
		lookup_data <- lookup_data[order(BNF_Paragraph, BNF_Chemical_Substance)]
		return(with(lookup_data, split(BNF_Chemical_Substance, BNF_Paragraph)))
	}
}

# UI
render_patient_count_ui <- function(counts) {
	div(
		style = "padding: 10px; background-color: #f8f9fa; border-radius: 4px; margin-top: 10px;",
		strong("Filtered patients: "),
		sprintf("Cases: %s (%s%%) | Controls: %s (%s%%) | Total: %s",
						prettyNum(counts$case_n, big.mark = ","),
						counts$case_pct,
						prettyNum(counts$control_n, big.mark = ","),
						counts$control_pct,
						prettyNum(counts$total_n, big.mark = ","))
	)
}

render_prevalence_reactable <- function(result_table, item_column_name) {
	validate(need(nrow(result_table) > 0, "No results found for selected filters"))

	# Add formatted OR column
	result_table <- add_or_formatted_column(result_table)

	# Rename columns with standard names
	col_names <- get_prevalence_column_names(item_column_name)
	colnames(result_table) <- c(col_names, "p_value", "p_adj", "OR_formatted")

	# Build column definitions
	columns <- list(
		colDef(
			name = item_column_name,
			minWidth = 200,
			cell = create_significance_cell(result_table)
		),
		`Cases (%)` = colDef(format = colFormat(digits = 2)),
		`Controls (%)` = colDef(format = colFormat(digits = 2))
	)
	names(columns)[1] <- item_column_name
	columns <- c(columns, get_or_reactable_columns())

	# Render with standard config
	do.call(reactable, c(
		list(
			data = result_table,
			columns = columns,
			details = function(index) create_pvalue_details(index, result_table)
		),
		get_standard_reactable_config()
	))
}

# Calculations
calculate_filtered_prevalence_table <- function(filtered_data, all_patids, id_column) {
	# Filter to final patient set
	filtered_data <- filtered_data[patid %in% all_patids]

	# Calculate frequency stats
	freq_stats <- calculate_frequency_stats(filtered_data, id_column)

	# Create prevalence ratio table
	create_prevalence_ratio_table(freq_stats, id_column, data_with_group = filtered_data)
}

get_filtered_patient_counts <- function(data, strat_var, filter1_patids = NULL, filter2_patids = NULL, filter3_patids = NULL, patient_data) {
	# Get total counts for each group (before filtering)
	total_case_count <- uniqueN(data[group == "case", patid])
	total_control_count <- uniqueN(data[group == "control", patid])

	# Start with all patients
	filtered_data <- data

	# Apply stratification first
	if (!is.null(strat_var) && strat_var != "") {
		filtered_data <- apply_patient_stratification(filtered_data, strat_var, patient_data)
	}

	# Then apply additional filters
	all_patids <- unique(filtered_data$patid)

	if (!is.null(filter1_patids)) {
		all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter1_patids)))
	}

	if (!is.null(filter2_patids)) {
		all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter2_patids)))
	}

	if (!is.null(filter3_patids)) {
		all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter3_patids)))
	}

	# Filter to final patient set
	filtered_data <- filtered_data[patid %in% all_patids]

	# Count by group (after filtering)
	case_count <- uniqueN(filtered_data[group == "case", patid])
	control_count <- uniqueN(filtered_data[group == "control", patid])
	total_count <- case_count + control_count

	# Calculate percentages relative to total in each group
	case_pct <- if (total_case_count > 0) round(case_count / total_case_count * 100, 1) else 0
	control_pct <- if (total_control_count > 0) round(control_count / total_control_count * 100, 1) else 0

	list(
		case_n = case_count,
		case_pct = case_pct,
		control_n = control_count,
		control_pct = control_pct,
		total_n = total_count
	)
}

get_filtered_patids <- function(base_data, filter1_patids = NULL, filter2_patids = NULL, filter3_patids = NULL) {
	all_patids <- unique(base_data$patid)

	if (!is.null(filter1_patids)) {
		all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter1_patids)))
	}

	if (!is.null(filter2_patids)) {
		all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter2_patids)))
	}

	if (!is.null(filter3_patids)) {
		all_patids <- bit64::as.integer64(intersect(as.character(all_patids), as.character(filter3_patids)))
	}

	all_patids[order(all_patids)]
}
