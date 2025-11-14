#' Visualization functions for co-prescription analysis
#' Creates heatmaps using vegawidget (matching CE module style)

#' Create OR heatmap for co-prescription analysis (single group)
#'
#' @param copresc_results data.table with drug1, drug2, or_value, significant columns
#' @param group_label Label for the group (e.g., "Cases" or "Controls")
#' @return vegaspec object
create_copresc_or_heatmap <- function(copresc_results, group_label = "Cases") {

	# Create long format data for heatmap
	long_df <- copresc_results[, .(
		drug1,
		drug2,
		or_value
		# significant,
		# label = ifelse(significant, "*", "")
	)]

	# Add reverse pairs to make symmetric
	long_df_reverse <- copresc_results[, .(
		drug1 = drug2,
		drug2 = drug1,
		or_value
		# significant,
		# label = ifelse(significant, "*", "")
	)]

	# Combine both directions
	long_df <- rbind(long_df, long_df_reverse)

	# Get unique drugs and create ordering
	all_drugs <- unique(c(long_df$drug1, long_df$drug2))

	# Create matrix for clustering
	or_matrix <- matrix(1, length(all_drugs), length(all_drugs),
											dimnames = list(all_drugs, all_drugs))

	for (i in seq_len(nrow(long_df))) {
		d1 <- long_df$drug1[i]
		d2 <- long_df$drug2[i]
		val <- long_df$or_value[i]

		or_matrix[d1, d2] <- val
	}

	# Replace NA with 1 for clustering
	or_matrix[is.na(or_matrix)] <- 1

	# Hierarchical clustering
	if (nrow(or_matrix) > 2) {
		dist_matrix <- dist(or_matrix, method = "euclidean")
		hclust_result <- hclust(dist_matrix, method = "ward.D2")
		drug_order <- rownames(or_matrix)[hclust_result$order]
	} else {
		drug_order <- all_drugs
	}

	# Calculate max OR for scale
	long_df[, or_log := log(or_value)]
	max_or <- max(long_df$or_value, na.rm = TRUE)

	# Create Vega-Lite spec
	spec <- list(
		`$schema` = vegawidget::vega_schema(),
		data = list(values = long_df),
		hconcat = list(
			# Main heatmap
			list(
				title = paste0("Co-prescription ORs: ", group_label),
				width = 600,
				height = 600,
				transform = list(list(filter = list(param = "brush"))),
				layer = list(
					# Heatmap rectangles
					list(
						mark = list(
							type = "rect",
							stroke = "white",
							strokeWidth = 1
						),
						encoding = list(
							x = list(
								field = "drug1",
								type = "nominal",
								title = NULL,
								sort = drug_order,
								axis = list(
									labelAngle = -45,
									labelFontSize = list(
										expr = "brush && brush.drug1 ? max(16, min(8, 18 - length(brush.drug1) * 0.5)) : 12"
									)
								)
							),
							y = list(
								field = "drug2",
								type = "nominal",
								title = NULL,
								sort = drug_order,
								axis = list(
									labelFontSize = list(
										expr = "brush && brush.drug2 ? max(16, min(8, 18 - length(brush.drug2) * 0.5)) : 12"
									)
								)
							),
							color = list(
								field = "or_log",
								type = "quantitative",
								scale = list(
									scheme = "goldred",
									reverse = FALSE,
									domain = c(0, max(long_df$or_log))
								),
								legend = list(labelExpr = "round(exp(datum.value))")
							),
							tooltip = list(
								list(field = "drug1", type = "nominal", title = "Drug 1"),
								list(field = "drug2", type = "nominal", title = "Drug 2"),
								list(field = "or_value", type = "quantitative", title = "OR", format = ".2f")
								# list(field = "significant", type = "nominal", title = "Significant")
							)
						)
					)
					# Significance markers
					# list(
					# 	mark = list(
					# 		type = "text",
					# 		fontSize = 20,
					# 		fontWeight = "bold"
					# 	),
					# 	encoding = list(
					# 		x = list(field = "drug1", type = "nominal", sort = drug_order),
					# 		y = list(field = "drug2", type = "nominal", sort = drug_order),
					# 		text = list(field = "label", type = "nominal"),
					# 		color = list(value = "white"),
					# 		opacity = list(
					# 			condition = list(
					# 				test = "datum.significant == true",
					# 				value = 1
					# 			),
					# 			value = 0
					# 		)
					# 	)
					# )
				)
			),
			# Brush overview (zoom control)
			list(
				title = "Drag to Zoom",
				width = 150,
				height = 150,
				params = list(list(
					name = "brush",
					select = list(type = "interval", encodings = c("x", "y"))
				)),
				mark = list(
					type = "rect",
					stroke = "white",
					strokeWidth = 0.5
				),
				encoding = list(
					x = list(
						field = "drug1",
						type = "nominal",
						sort = drug_order,
						axis = list(
							title = NULL,
							labels = FALSE,
							ticks = FALSE
						)
					),
					y = list(
						field = "drug2",
						type = "nominal",
						sort = drug_order,
						axis = list(
							title = NULL,
							labels = FALSE,
							ticks = FALSE
						)
					),
					color = list(
						field = "or_log",
						type = "quantitative",
						scale = list(
							scheme = "goldred",
							reverse = FALSE,
							domain = c(0, max(long_df$or_log))
						),
						legend = NULL
					)
				)
			)
		),
		resolve = list(scale = list(color = "shared")),
		config = list(
			view = list(stroke = NULL),
			font = "Lato, -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif"
		)
	)

	return(vegawidget::as_vegaspec(spec))
}

#' Create Forest Plot for Co-prescription Analysis
#'
#' Creates a forest plot showing case and control ORs side-by-side for each drug pair
#' Uses rule marks for CI lines and point marks for OR estimates
#'
#' @param data data.table with columns: drug_pair, drug1, drug2,
#'   case_or, case_ci_lower, case_ci_upper,
#'   control_or, control_ci_lower, control_ci_upper
#' @return vegaspec object
create_copresc_forest_plot <- function(data) {

	# Prepare data in long format for grouped display
	# Each drug pair will have two rows: one for case, one for control
	data_long <- rbindlist(list(
		data[, .(
			drug_pair = drug_pair,
			drug1 = drug1,
			drug2 = drug2,
			group = "case",
			or = case_or,
			ci_lower = case_ci_lower,
			ci_upper = case_ci_upper
		)],
		data[, .(
			drug_pair = drug_pair,
			drug1 = drug1,
			drug2 = drug2,
			group = "control",
			or = control_or,
			ci_lower = control_ci_lower,
			ci_upper = control_ci_upper
		)]
	))

	# Filter out invalid values
	data_long <- data_long[!is.na(or) & is.finite(or) &
												 	!is.na(ci_lower) & is.finite(ci_lower) &
												 	!is.na(ci_upper) & is.finite(ci_upper)]

	if (nrow(data_long) == 0) {
		return(NULL)
	}

	# Get unique drug pairs in order for y-axis
	drug_pair_order <- unique(data$drug_pair)

	# Calculate max CI for x-axis scale
	max_ci <- max(data_long$ci_upper, na.rm = TRUE)
	x_max <- min(max_ci * 1.1, 20)  # Cap at 20 for readability

	# Create Vega-Lite spec
	spec <- list(
		`$schema` = vegawidget::vega_schema(),
		data = list(values = data_long),
		width = 600,
		height = 25 * length(drug_pair_order),
		title = "Co-prescription Odds Ratios: Cases vs Controls",
		encoding = list(
			y = list(
				field = "drug_pair",
				type = "nominal",
				title = NULL,
				sort = drug_pair_order,
				axis = list(
					labelFontSize = 11,
					labelLimit = 250,
					minExtent = 250
				)
			),
			yOffset = list(
				field = "group",
				type = "nominal",
				scale = list(
					domain = c("case", "control"),
					range = c(0, 20)  # Offset to show case/control on same row
				)
			)
		),
		layer = list(
			#Reference line at OR = 1
			list(
				mark = list(
					type = "rule",
					strokeDash = c(3, 3),
					color = "gray"
				),
				encoding = list(
					x = list(datum = 1),
					y = list()
				)
			),
			# CI lines (rule mark with x and x2)
			list(
				mark = list(
					type = "rule",
					size = 2
				),
				encoding = list(
					x = list(
						field = "ci_lower",
						type = "quantitative",
						title = "Odds Ratio",
						# scale = list(
						# 	domain = c(0, x_max),
						# 	type = "log",
						# 	nice = FALSE
						# ),
						axis = list(
							grid = TRUE,
							tickCount = 10
						)
					),
					x2 = list(
						field = "ci_upper"
					),
					color = list(
						field = "group",
						type = "nominal",
						scale = list(
							domain = c("case", "control"),
							range = c("#e74c3c", "#3498db")  # Red for cases, blue for controls
						),
						legend = list(
							title = "Group",
							orient = "top"
						)
					)
				)
			),
			# OR point estimates
			list(
				mark = list(
					type = "point",
					filled = TRUE,
					size = 80
				),
				encoding = list(
					x = list(
						field = "or",
						type = "quantitative"
					),
					color = list(
						field = "group",
						type = "nominal"
					),
					tooltip = list(
						list(field = "drug_pair", type = "nominal", title = "Drug Pair"),
						list(field = "group", type = "nominal", title = "Group"),
						list(field = "or", type = "quantitative", title = "OR", format = ".2f"),
						list(field = "ci_lower", type = "quantitative", title = "CI Lower", format = ".2f"),
						list(field = "ci_upper", type = "quantitative", title = "CI Upper", format = ".2f")
					)
				)
			)
		)
	)

	vegawidget::as_vegaspec(spec)
}