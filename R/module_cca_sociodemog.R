module_cca_sociodemographics_ui <- function(id) {
	ns <- NS(id)
	accordion(
		open = FALSE,
		accordion_panel(
			title = "Table",
			value = "demog_tables",
			gt::gt_output(ns("sociodemog_table"))
		)
	)
}

module_cca_sociodemographics_server <- function(id, patient_data_r) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		# Calculate sociodemographic statistics
		sociodemog_data <- reactive({
			req(patient_data_r())
			browser()
			patient_df <- patient_data_r()

			# Calculate totals
			cases_total <- patient_df[treatment == 1, .N]
			controls_total <- patient_df[treatment == 0, .N]

			# Define demographic groups to process
			demog_groups <- list(
				list(var = "sex", label = "Sex"),
				list(var = "eth_group", label = "Ethnicity"),
				list(var = "imd_quintile", label = "IMD Quintile"),
				list(var = "mltc_group", label = "# LTCs"),
				list(var = "pp_group", label = "# PP")
			)

			# Process each demographic group
			results_list <- lapply(demog_groups, function(g) {
				# Calculate for cases (treatment == 1)
				cases_stats <- patient_df[treatment == 1, .N, by = c(g$var)]
				cases_stats[, group := g$label]
				cases_stats[, cases_pct := round(N / cases_total * 100, 1)]
				cases_stats[, cases_value := sprintf("%d (%.1f%%)", N, cases_pct)]
				setnames(cases_stats, old = g$var, new = "category")

				# Calculate for controls (treatment == 0)
				controls_stats <- patient_df[treatment == 0, .N, by = c(g$var)]
				controls_stats[, controls_pct := round(N / controls_total * 100, 1)]
				controls_stats[, controls_value := sprintf("%d (%.1f%%)", N, controls_pct)]
				setnames(controls_stats, old = g$var, new = "category")

				# Merge cases and controls
				merged <- merge(
					cases_stats[, .(group, category, cases_value)],
					controls_stats[, .(category, controls_value)],
					by = "category",
					all = TRUE
				)

				# Handle missing values (replace with "0 (0.0%)")
				merged[is.na(cases_value), cases_value := "0 (0.0%)"]
				merged[is.na(controls_value), controls_value := "0 (0.0%)"]

				# Order IMD quintile if applicable
				if (g$var == "imd_quintile") {
					merged[, category := factor(category, levels = c("1", "2", "3", "4", "5"))]
					setorder(merged, category)
				}

				return(merged)
			})

			# Combine all demographic groups
			final_table <- rbindlist(results_list, use.names = TRUE)

			return(final_table)
		})

		# Render table
		output$sociodemog_table <- gt::render_gt({
			req(sociodemog_data())

			gt::gt(
				sociodemog_data(),
				rowname_col = "category",
				groupname_col = "group"
			) |>
				gt::cols_label(
					cases_value = "Cases",
					controls_value = "Controls"
				) |>
				gt::tab_style(
					style = list(
						gt::cell_text(weight = "bold")
					),
					locations = gt::cells_row_groups()
				)
		})
	})
}