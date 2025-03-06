#' @export
module_burden_group_ui <- function(id) {
  ns <- NS(id)
  card(card_header("Drug association with burden groups"),
       dataTableOutput(ns("substance_table")))
}

#' @export
module_burden_group_server <- function(id, subst_pp_df) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns
    output$substance_table <- renderDataTable({
      df <- subst_pp_df()
      # tmp <- df[,.(patid,pp_group)] |> unique()
      # print(table(tmp$pp_group))
      df <- df[, .N, .(substance, pp_group)]
      df <-
        dcast(
          df,
          substance ~ pp_group,
          value.var = "N",
          fill = 0
        )
      row_totals <- rowSums(df[, -1])
      df <- df[order(row_totals, decreasing = T),]
    }, rownames = FALSE)
  })
}