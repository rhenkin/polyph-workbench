histogram_calc <- function(dt, x_field, bin_step = 5) {
  # Get values from the column
  x_values <- dt[[x_field]]
  # Calculate bin boundaries
  # min_val <- max(2,floor(min(x_values, na.rm = TRUE) - (bin_step/2)))  # Subtract half bin_step
  # max_val <- ceiling(max(x_values, na.rm = TRUE) + (bin_step/2))  # Add half bin_step
  # breaks <- seq(from = min_val, to = max_val, by = bin_step)

  # Use hist() function to get counts, but don't plot
  h <- hist(x_values, plot = FALSE)

  # Create data.frame with results
  data.frame(
    bin_interval = paste0("[", paste(h$breaks[-length(h$breaks)], h$breaks[-1], sep="-"), ")"),
    count = h$counts
  )
}

precalc_histogramPlot <- function(df, x_field, label = NULL, w = 200, h = 300, title = NULL, bin_step = 5) {
  if (is.null(label)) label <- x_field

  # Calculate histogram
  hist_data <- histogram_calc(df, x_field, bin_step)
  list(
    `$schema` = vega_schema(),
    data = list(values = hist_data),
    title = title,
    mark = "bar",
    width = w,
    encoding = list(
      x = list(
        field = "bin_interval",
        type = "ordinal",
        title = label,
        sort = hist_data$bin_interval
        # scale = list(
        #   # Ensure x-axis extends slightly beyond the data
        #   domain = c(min(hist_data$bin_start), max(hist_data$bin_start) + bin_step)
        # )
      ),
      y = list(
        field = "count",
        type = "quantitative",
        title = "Count"
      ),
      width = list(value = bin_step)  # Make bars touch each other
    )
  ) |> as_vegaspec()
}

horizontalBarPlot <- function(df, y_field, label = NULL, w = 200, h = 300,  title = NULL) {

  if (is.null(label)) label <- y_field

  list(
    `$schema` = vega_schema(),
    data = list(values = df),
    title = title,
    mark = "bar",
    width = w,
    height = h,
    encoding = list(
      x = list(field = "N", type = "quantitative"),
      y = list(field = y_field, type = "ordinal", title = label)
    )
  ) |> as_vegaspec()

}


#' @export
histogramPlot <- function(df, x_field, label = NULL, w = 200, h = 300,  title = NULL) {

  if (is.null(label)) label <- x_field

  list(
    `$schema` = vega_schema(),
    data = list(values = df),
    title = title,
    mark = "bar",
    width = w,
    #height = h,
    encoding = list(
      x = list(bin = list(step = 5), field = x_field, title = label),
      y = list(aggregate = "count")
    )
  ) |> as_vegaspec()

}

#' @export
facetedHistogramPlot <- function(df, x_field, facet_var, label = NULL, w = 300, h = 300,  title = NULL) {

  if (is.null(label)) label <- x_field

  list(
    `$schema` = vega_schema(),
    data = list(values = df),
    facet = list(row = list(field = facet_var)),
    title = title,
    spec = list(
      mark = "bar",
      encoding = list(
        x = list(bin = list(step = 5), field = x_field, title = label),
        y = list(aggregate = "count")
      )
    )
  ) |> as_vegaspec()
}

calculateDensityByGroup_dt <- function(dt, x_field, facet_var, n = 512, from = NULL, to = NULL) {
  # Ensure input is data.table
  if (!is.data.table(dt)) dt <- as.data.table(dt)

  # If from/to not specified, determine range from data
  if (is.null(from)) from <- dt[, min(get(x_field), na.rm = TRUE)]
  if (is.null(to)) to <- dt[, max(get(x_field), na.rm = TRUE)]

  # Function to calculate density for each group
  calc_density <- function(x) {
    dens <- density(x, n = n, from = from, to = to)
    list(x = dens$x, density = dens$y)
  }

  # Calculate density by group using data.table
  result_dt <- dt[,
                  calc_density(get(x_field)),
                  by = get(facet_var)
  ]

  # Rename columns to match expected format
  setnames(result_dt,
           old = c("get", "x", "density"),
           new = c(facet_var, x_field, "density"))

  return(result_dt)
}

precalc_facetedViolinPlot <- function(df, x_field, facet_var, label = NULL, w = 80, title = NULL) {
  if (is.null(label)) label <- x_field

  # Use the pre-calculated density values directly
  list(
    `$schema` = vega_schema(),
    data = list(values = df),  # df should now contain pre-calculated density values
    width = w,
    mark = list(type = "area", orient = "horizontal"),
    title = title,
    encoding = list(
      y = list(field = x_field, type = "quantitative"),
      x = list(field = "density", type = "quantitative", title = FALSE,
               axis = list(grid = FALSE, labels = FALSE, ticks = TRUE),
               stack = "center", impute = NULL),
      column = list(field = facet_var, type = "nominal", title = NULL,
                    sort = sort(unique(df[[facet_var]])),
                    header = list(labelOrient = "bottom", labelPadding = 0,
                                  titleOrient = "bottom")),
      color = list(field = facet_var, type = "nominal", legend = NULL)
    ),
    config = list(
      view = list(stroke = NULL),
      facet = list(spacing = 0)
    )
  ) |> as_vegaspec()
}


#' @export
facetedViolinPlot <- function(df, x_field, facet_var, label = NULL, w = 80, title = NULL) {

  if (is.null(label)) label <- x_field

  list(
    `$schema` = vega_schema(),
    data = list(values = df),
    width = w,
    transform = list(
      list(density = x_field, as = c(x_field, "density"), groupby = list(facet_var), extent = c(0, 50))
    ),
    mark = list(type = "area", orient = "horizontal"),
    title = title,
    encoding = list(
      y = list(field = x_field, type = "quantitative"),
      x = list(field = "density", type = "quantitative", title = FALSE, axis = list(grid = FALSE, labels = FALSE, ticks = TRUE), stack = "center", impute = NULL),
      column = list(field = facet_var, type = "nominal", title = NULL, sort = sort(unique(df[[facet_var]])), header = list(labelOrient = "bottom", labelPadding = 0, titleOrient = "bottom")),
      color = list(field = facet_var, type = "nominal", legend = NULL)
    ),
    config = list(
      #view = list(continuousWidth = 300, continuousHeight = 300, stroke = NULL),
      view = list(stroke = NULL),
      facet = list(spacing = 0)
    )
  ) |> as_vegaspec()

}


pairedBarCharts <- function(df, x, y) {
  plot_data <- df[, .(count = .N),
                  by = c(x, y, "pp")
  ][, proportion := count/sum(count),
    by = c(x, y)]

  # Create the Vega-Lite spec with preprocessed data
  vega_spec <- list(
    `$schema` = "https://vega.github.io/schema/vega-lite/v5.json",
    data = list(
      values = plot_data
    ),
    facet = list(
      row = list(field = x, type = "nominal", sort = levels(plot_data[[x]])),
      column = list(field = y, type = "nominal", sort = levels(plot_data[[y]]))
    ),
    spec = list(
      width = 150,
      height = 150,
      mark = list(
        type = "bar",
        tooltip = TRUE
      ),
      encoding = list(
        x = list(
          field = "pp",
          type = "quantitative",
          title = "PP Value",
          scale = list(domain = list(2, 30)),
          axis = list(tickMinStep = 1)
        ),
        y = list(
          field = "proportion",
          type = "quantitative",
          title = NA,
          axis = list(format = ".1%")
        ),
        tooltip = list(
          list(field = "pp", type = "quantitative", title = "PP Value"),
          list(field = "count", type = "quantitative", title = "Count"),
          list(field = "proportion", type = "quantitative", title = "Proportion", format = ".1%")
        )
      )
    ),
    config = list(
      view = list(stroke = NULL),
      facet = list(spacing = 8),
      header = list(
        labelFontSize = 12
      ),
      axis = list(
        labelFont = "sans-serif",
        titleFont = "sans-serif",
        labelFontSize = 12,
        titleFontSize = 14
      )
    )
  ) |> as_vegaspec()
}