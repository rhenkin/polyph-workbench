create_boxplot_spec <- function(df, title) {
	spec <- list(
		`$schema` = vega_schema(),
		title = title,
		height = 160,
		# width = width,
		data = list(values = df),
		encoding = list(
			y = list(field = "group", type = "nominal", title = NULL)
		),
		layer = list(
			list(
				mark = list(type = "rule"),
				encoding = list(
					x = list(
						field = "lower",
						type = "quantitative",
						scale = list(zero = FALSE),
						title = NULL
					),
					x2 = list(field = "upper")
				)
			),
			list(
				mark = list(type = "bar", size = 36),
				encoding = list(
					x = list(field = "q1", type = "quantitative"),
					x2 = list(field = "q3"),
					color = list(field = "group", type = "nominal", legend = NULL,
											 scale = list(range = list("#e74c3c", "#2C3E50")))
				)
			),
			list(
				mark = list(type = "tick", color = "white", size = 36),
				encoding = list(
					x = list(field = "median", type = "quantitative")
				)
			),
			list(
				mark = list(type = "text", dy = -25, align = "center"),
				encoding = list(
					x = list(field = "median", type = "quantitative"),
					text = list(field = "median"),
					color = list(value = "black")
				)
			)
		),
		config = list(
			view = list(stroke = NULL),
			font = "Lato, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, \"Helvetica Neue\", Arial, sans-serif, \"Apple Color Emoji\", \"Segoe UI Emoji\", \"Segoe UI Symbol\""
		)
	)
	spec
}