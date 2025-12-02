
grouped_bar_plot <- function(df, y_var, height = 250, width = 300, title = NULL) {



	spec <-  list(
		`$schema` = vega_schema(),
		title = title,
		height = height,
		width = width,
		data = list(values = df),
		layer = list(
			# Bar layer
			list(
				mark = list(type = "bar", tooltip = TRUE),
				encoding = list(
					x = list(
						field = "pct",
						type = "quantitative",
						title = "Prevalence (%)",
						scale = list(domainMax = max(df$pct)*1.05)
					),
					y = list(
						field = y_var,
						type = "ordinal",
						axis = list(grid = FALSE, title = NULL),
						sort = list(field = "max_pct", op = "max", order = "descending")
					),
					yOffset = list(field = "group"),
					color = list(
						field = "group",
						scale = list(range = list("#e74c3c", "#2C3E50")),
						legend = NULL
					)
				)
			),
			# Text label layer
			list(
				mark = list(type = "text", align = "left", dx = 5),
				encoding = list(
					x = list(
						field = "max_pct",
						type = "quantitative"
					),
					y = list(
						field = y_var,
						type = "ordinal",
						sort = list(field = "max_pct", op = "max", order = "descending")
					),
					text = list(field = "diff_label"),
					color = list(value = "black")
				),
				transform = list(
					list(filter = "datum.group == 'case'")
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

pyramid_plot <- function(df, y_var, y_sort, height = 250, side_width = 150, title = NULL) {

# Create the Vega-Lite specification as a list
	#y_sort <- rev(c(2:19,"20+"))
	max_domain <- max(df$pct)+0.01
spec <- list(
	`$schema` = vega_schema(),
	title = title,
	data = list(values = df),
	spacing = 0,
	hconcat = list(
		# Female side
		list(
			height = height,
			width = side_width,
			transform = list(
				list(filter = list(field = "treatment", equal = 1))
			),
			title = "Cases",
			mark = "bar",
			encoding = list(
				y = list(
					field = y_var,
					type = "nominal",
					axis = NULL,
					sort = y_sort
				),
				x = list(
					aggregate = "sum",
					field = "pct",
					title = "N",
					scale = list(domain = c(0,max_domain)),
					axis = list(values = seq(from=0.04, to=0.2, by =0.02),
											format = ".0%"),
					sort = "descending"
				),
				color = list(
					field = "treatment",
					scale = list(range = list("#2C3E50", "#e74c3c")),
					legend = NULL
				)
			)
		),
		# Middle age labels
		list(
			height = height,
			width = 20,
			view = list(stroke = NULL),
			mark = list(
				type = "text",
				align = "center"
			),
			encoding = list(
				y = list(field = y_var, type = "nominal", axis = NULL,
								 sort = y_sort),
				text = list(field = y_var, type = "nominal")
			)
		),
		# Male side
		list(
			height = height,
			width = side_width,
			transform = list(
				list(filter = list(field = "treatment", equal = 0))
			),
			title = "Controls",
			mark = "bar",
			encoding = list(
				y = list(
					field = y_var,
					type = "nominal",
					title = NULL,
					axis = NULL,
					sort = y_sort
				),
				x = list(
					aggregate = "sum",
					field = "pct",
					title = "N",
					scale = list(domain = c(0,max_domain)),
					axis = list(values = seq(from=0.04, to=0.2, by =0.02),
											format = ".0%")
				),
				color = list(
					field = "treatment",
					legend = NULL
				)
			)
		)
	),
	config = list(
		view = list(stroke = NULL),
		font = "Lato, -apple-system, BlinkMacSystemFont, \"Segoe UI\", Roboto, \"Helvetica Neue\", Arial, sans-serif, \"Apple Color Emoji\", \"Segoe UI Emoji\", \"Segoe UI Symbol\""
	)
)

return(spec)

}