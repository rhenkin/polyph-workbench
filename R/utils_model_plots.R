#' @export
dot_plot <- function(df, add_color = FALSE) {
  spec <- list(
    `$schema` = vega_schema(),
    width = "container",
    height = 600,
    title = list(
      text = "Association between being prescribed one of top 50 drugs and high polypharmacy burden",
      subtitle = paste0("Adjusted by Age at outcome and number of LTCs / Ordered by IRR / Total number of patients in models: ", df[1,n_total]) #, " / % indicates positive cases")
    ),
    data = list(
      values = df
      #name = "source"
    ),
    layer = list(
      # Confidence interval lines
      list(
        mark = "rule",
        encoding = list(
          x = list(
            field = "lower_CI",
            type = "quantitative"
          ),
          x2 = list(
            field = "upper_CI"
          ),
          y = list(
            field = "substance",
            type = "nominal",
            sort = "-IRR"
          )
        )
      ),
      
      # Point estimates
      list(
        mark = list(
          type = "point",
          filled = TRUE,
          size = 100
        ),
        encoding = list(
          x = list(
            field = "IRR",
            type = "quantitative",
            title = "Incidence Rate Ratio (95% CI)",
            axis = list(labelFontSize = 12)
          ),
          y = list(
            field = "substance",
            type = "nominal",
            sort = "-IRR",
            axis = list(labelFontSize = 12)
          ),
          tooltip = list(
            list(field = "substance", type = "nominal", title = "BNF code"),
            list(field = "IRR", type = "quantitative", format = ".2f", title = "IRR"),
            list(field = "n_exposed", type = "quantitative", title = "N exposed")
          )
        )
      ),
      
      # Text annotations for n_exposed
      list(
        mark = list(
          type = "text",
          align = "left",
          dx = 5,
          fontSize = 11
        ),
        encoding = list(
          x = list(
            field = "upper_CI",
            type = "quantitative"
          ),
          y = list(
            field = "substance",
            type = "nominal",
            sort = "-IRR"
          ),
          text = list(
            field = "n_exposed",
            type = "nominal"
          )
        )
      ),
      
      # Reference line at IRR = 1
      list(
        mark = "rule",
        encoding = list(
          x = list(datum = 1),
          color = list(value = "#000000"),
          size = list(value = 1)
        )
      )
    ),
    
    config = list(
      axis = list(
        grid = TRUE,
        gridColor = "#DCDCDC"
      )
    )
  )
  if (add_color) {
    spec$layer[[1]]$color <- list(
      field = "significance",
      type = "nominal"
    )
    spec$layer[[2]]$color <-  list(
      field = "significance",
      type = "nominal",
      title = "FDR-adjusted p-value",
      scale = list(
        domain = list("Significant", "Non-significant"),
        values = list("#ff8800", "#0088ff")
      )
    )
    spec$layer[[2]]$tooltip[[4]] <-  list(field = "p_value", type = "quantitative", format = ".3f", title = "p-value")
  }
  spec |> as_vegaspec()
}