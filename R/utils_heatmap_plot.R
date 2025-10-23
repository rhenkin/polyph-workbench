cor_heatmap <- function(pcor_matrix) {

	df_long <- data.table(
		Source = rep(rownames(pcor_matrix), ncol(pcor_matrix)),
		Target = rep(colnames(pcor_matrix), each = nrow(pcor_matrix)),
		Value = as.vector(pcor_matrix)
	)
	df_long <- df_long[Value != 0]

	hc <- hclust(dist(abs(pcor_matrix)))
	cluster_order <- hc$order
	var_order <- rownames(pcor_matrix)[cluster_order]

	list(
		`$schema` = vegawidget::vega_schema(),
		data = list(values = df_long),
		hconcat = list(
			list(
				title = "Partial correlation heatmap",
				width = 400,
				height = 400,
				transform = list(list(filter = list(param = "brush"))),
				mark = list(
					type = "rect",
					stroke = "white",
					strokeWidth = 1
				),
				encoding = list(
					x = list(
						field = "Source",
						type = "nominal",
						title = NULL,
						sort = var_order ,
						axis = list(labelAngle = -45, labelFontSize = 12)
					),
					y = list(
						field = "Target",
						type = "nominal",
						title = NULL,
						sort = var_order ,
						axis = list(labelFontSize = 12)
					),
					color = list(
						field = "Value",
						type = "quantitative",
						title = "Correlation",
						scale = list(scheme = "plasma", #reverse = FALSE,
												 domain = c(0, max(df_long$Value)))
					),
					tooltip = list(
						list(field = "Source", type = "nominal"),
						list(field = "Target", type = "nominal"),
						list(
							field = "Value",
							type = "quantitative",
							format = ".2f"
						)
					)
				)
			),
			list(
				title = "Brush to Zoom",
				width = 120,
				height = 120,
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
						field = "Source",
						type = "nominal",
						sort = var_order ,
						axis = list(
							title = NULL,
							labels = FALSE,
							ticks = FALSE
						)
					),
					y = list(
						field = "Target",
						type = "nominal",
						sort = var_order ,
						axis = list(
							title = NULL,
							labels = FALSE,
							ticks = FALSE
						)
					),
					color = list(
						field = "Value",
						type = "quantitative",
						scale = list(scheme = "plasma", #reverse = FALSE,
												 domain = c(0, max(df_long$Value)))
					#	legend = NULL
					)
				)
			)
		),
		resolve = list(scale = list(color = "shared"))
	)
}

or_heatmap <- function(long_or_df, or_mat) {

	# or_matrix <- dcast(long_or_df, drug1 ~ drug2, value.var = "or", fun.aggregate = length, fill = 1)
	#
	# # Set row names and remove drug1 column
	# or_matrix_clean <- as.matrix(or_matrix[, -1])
	# rownames(or_matrix_clean) <- or_matrix$drug1
	#
	# # Handle missing values (set to 1 for clustering, since OR=1 means no association)
	or_mat[is.na(or_mat)] <- 1

	# Hierarchical clustering
	dist_matrix <- dist(or_mat, method = "euclidean")
	hclust_result <- hclust(dist_matrix, method = "ward.D2")

	# Get the reordered drug names
	var_order <- rownames(or_mat)[hclust_result$order]

	long_or_df[, or_log := log(or)]
	#long_or_df[, or_ci := ]
	# hc <- hclust(dist(abs(pcor_matrix)))
	# cluster_order <- hc$order
	# var_order <- rownames(pcor_matrix)[cluster_order]
	#var_order <- sort(unique(long_or_df$drug1))
	list(
		`$schema` = vegawidget::vega_schema(),
		data = list(values = long_or_df),
		hconcat = list(
			list(
				title = "OR heatmap",
				width = 600,
				height = 600,
				transform = list(list(filter = list(param = "brush"))),
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
						sort = var_order ,
						axis = list(labelAngle = -45, labelFontSize = list(expr = "brush && brush.drug1 ? max(16, min(8, 18 - length(brush.drug1) * 0.5)) : 12"))
					),
					y = list(
						field = "drug2",
						type = "nominal",
						title = NULL,
						sort = var_order ,
						axis = list(labelFontSize = list(expr = "brush && brush.drug2 ? max(16, min(8, 18 - length(brush.drug2) * 0.5)) : 12"))
					),
					color = list(
						field = "or_log",
						type = "quantitative",
						title = "Odds ratio",
						# scale = list(type = "quantile", scheme = list(name = "goldred", count = 5))
						scale = list(scheme = "goldred", reverse = FALSE,
						 						 domain = c(0, max(long_or_df$or_log))),
						legend = list(labelExpr = "round(exp(datum.value))")
					),
					tooltip = list(
						list(field = "drug1", type = "nominal"),
						list(field = "drug2", type = "nominal"),
						list(
							field = "or",
							type = "quantitative",
							format = ".2f"
						)
					)
				)
			),
			list(
				title = "Brush to Zoom",
				width = 160,
				height = 160,
				params = list(list(
					name = "brush",
					select = list(type = "interval", encodings = c("x", "y"),
												mark = list(fill = "#ccc", stroke = "#333", fillOpacity = 0.5))
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
						sort = var_order ,
						axis = list(
							title = NULL,
							labels = FALSE,
							ticks = FALSE
						)
					),
					y = list(
						field = "drug2",
						type = "nominal",
						sort = var_order ,
						axis = list(
							title = NULL,
							labels = FALSE,
							ticks = FALSE
						)
					),
					color = list(
						field = "or_log",
						type = "quantitative",
						#scale = list(type = "quantile", scheme = list(name = "goldred", count = 5))
						scale = list(scheme = "goldred", reverse = FALSE,
												 domain = c(0, max(long_or_df$or_log)))
						#	legend = NULL
					)
				)
			)
		),
		resolve = list(scale = list(color = "shared"))
	)
}

