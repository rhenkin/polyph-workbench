#' Visualization functions for co-prescription analysis
#' Creates heatmaps using vegawidget (matching CE module style)

#' Create OR difference heatmap for co-prescription analysis
#'
#' @param copresc_results data.table with drug1, drug2, or_diff, significant columns
#' @return vegaspec object
create_copresc_heatmap <- function(copresc_results) {

	# Create long format data for heatmap
	long_df <- copresc_results[, .(
		drug1,
		drug2,
		or_diff,
		significant,
		label = ifelse(significant, "*", "")
	)]

	# Get unique drugs and create ordering
	all_drugs <- unique(c(long_df$drug1, long_df$drug2))

	# Create matrix for clustering
	or_diff_matrix <- matrix(0, length(all_drugs), length(all_drugs),
													 dimnames = list(all_drugs, all_drugs))

	for (i in seq_len(nrow(long_df))) {
		d1 <- long_df$drug1[i]
		d2 <- long_df$drug2[i]
		val <- long_df$or_diff[i]

		or_diff_matrix[d1, d2] <- val
		or_diff_matrix[d2, d1] <- val
	}

	# Replace NA with 0 for clustering
	or_diff_matrix[is.na(or_diff_matrix)] <- 0

	# Hierarchical clustering
	if (nrow(or_diff_matrix) > 2) {
		dist_matrix <- dist(or_diff_matrix, method = "euclidean")
		hclust_result <- hclust(dist_matrix, method = "ward.D2")
		drug_order <- rownames(or_diff_matrix)[hclust_result$order]
	} else {
		drug_order <- all_drugs
	}

	# Calculate max absolute value for symmetric color scale
	max_abs_val <- max(base::abs(long_df$or_diff), na.rm = TRUE)

	# Create Vega-Lite spec
	spec <- list(
		`$schema` = vegawidget::vega_schema(),
		data = list(values = long_df),
		hconcat = list(
			# Main heatmap
			list(
				title = "Co-prescription OR Differences: Cases vs Controls",
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
								field = "or_diff",
								type = "quantitative",
								scale = list(
									scheme = "redblue",
									reverse = TRUE,
									domain = c(-max_abs_val, max_abs_val),
									domainMid = 0
								),
								legend = list(
									title = "OR Difference (Case - Control)",
									orient = "right"
								)
							),
							tooltip = list(
								list(field = "drug1", type = "nominal", title = "Drug 1"),
								list(field = "drug2", type = "nominal", title = "Drug 2"),
								list(field = "or_diff", type = "quantitative", title = "OR Difference", format = ".2f"),
								list(field = "significant", type = "nominal", title = "Significant (p<0.05)")
							)
						)
					),
					# Significance markers
					list(
						mark = list(
							type = "text",
							fontSize = 20,
							fontWeight = "bold"
						),
						encoding = list(
							x = list(field = "drug1", type = "nominal", sort = drug_order),
							y = list(field = "drug2", type = "nominal", sort = drug_order),
							text = list(field = "label", type = "nominal"),
							color = list(value = "black"),
							opacity = list(
								condition = list(
									test = "datum.significant == true",
									value = 1
								),
								value = 0
							)
						)
					)
				)
			),
			# Brush overview (zoom control)
			list(
				title = "Brush to Zoom",
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
						field = "or_diff",
						type = "quantitative",
						scale = list(
							scheme = "redblue",
							reverse = TRUE,
							domain = c(-max_abs_val, max_abs_val),
							domainMid = 0
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