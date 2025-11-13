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