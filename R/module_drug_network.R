module_drug_network_ui <- function(id) {
	ns <- NS(id)
	accordion_panel(
		"Co-occurrence networks",
		accordion(
			accordion_panel(
				title = "Network",
				p(
					paste("Network built using log-linear models between pairs of drugs.",
					"If checked, connections are displayed if FDR-adjusted p-value is lower than 0.05. Meta-analysis is performed if a control variable is selected.",
					"Red connections indicate OR > 1 and blue connections indicate OR < 1. Prevalence filtering is applied per stratum.")
				),
				div(
					style = "display: flex",
					numericInput(ns("subst_min_prev"), label = "Minimum prevalence for drug:", min = 0, value = 0.01, max = 1, step = 0.01),
					numericInput(ns("pair_min_prev"), label = "Minimum co-occurrence for pairs:", min = 0, value = 0.02, max = 1, step = 0.01),
					selectizeInput(ns("network_confounder"), label = "Control variable:", choices = c("", "age_group", "pp_group", "imd_quintile", "eth_group", "sex")),
					checkboxInput(ns("network_signif_only"), label = "adj. p-value < 0.05", value = TRUE),
					sliderInput(ns("network_or_filter"), label = "OR range", min = -50, max = 50, value = c(1, 50))
				),
				div(
					style = "display: flex",
					selectizeInput(ns("network_strat_variable"), label = "Stratifying variable:", choices = c("sex", "eth_group",  "imd_quintile", "age_group", "pp_group"), selected = "",
												 options=list(
												 	allowEmptyOption=TRUE,
												 	showEmptyOptionInDropdown=TRUE,
												 	emptyOptionLabel="No stratification"
												 )),
					selectizeInput(ns("network_strat_value"), "Select stratum:", choices = NULL, selected = "",
												 options=list(
												 	allowEmptyOption=TRUE,
												 	showEmptyOptionInDropdown=TRUE,
												 	emptyOptionLabel="No stratification"
												 )),
					selectizeInput(ns("network_agg_level"), "Drug grouping:",  choices = c("BNF_Chemical_Substance", "BNF_Paragraph", "BNF_Section", "BNF_Chapter"))
				),
				visNetworkOutput(ns("substance_network_graph"), height = "600px")
			),
			accordion_panel(
				title = "Table",
				dataTableOutput(ns("substance_network_table"))
			)
		)
	)
}

module_drug_network_server <- function(id, subst_pp_df, bnf_lookup) {
		moduleServer(id, function(input, output, session) {
			ns <- NS(id)


			base_network_data <- reactive({
				req(subst_pp_df())
				df <- subst_pp_df()[,.(patid, sex, eth_group, imd_quintile, outcome_age, substance, pp_group)]
				df[, outcome_age := round(outcome_age/365.25)]
				df <- create_value_groups(df, quantile(df$outcome_age, probs = seq(0, 1, length.out = 6)), "outcome_age", "age_group")
				df[,substance := bnf_lookup[match(substance,bnf_lookup$BNF_Chemical_Substance),input$network_agg_level,with=FALSE]]
				df
			})

			observe({
				req(!is.null(input$network_strat_variable))

				selectize_options <- list(
					allowEmptyOption = TRUE,
					showEmptyOptionInDropdown = TRUE,
					emptyOptionLabel = "No stratification"
				)

				if (input$network_strat_variable != "") {
					select_values <- sort(unique(base_network_data()[[input$network_strat_variable]]))
					updateSelectizeInput(session, "network_strat_value", choices = select_values, options = selectize_options)
				} else {
					updateSelectizeInput(session, "network_strat_value", choices = c(""))
				}
			}) |>
				bindEvent(input$network_strat_variable)

			substance_network_models <- reactive({
				req(base_network_data())
				df <- base_network_data()
				# Only apply filtering if both variable AND value are properly selected
				net_strat_var <- isolate({ input$network_strat_variable})
				if (!is.null(input$network_strat_value) && input$network_strat_value != "") {
					# Use tryCatch to handle any errors during filtering
					tryCatch({
						df <- df[get(net_strat_var) == input$network_strat_value]
					}, error = function(e) {
						# Log the error but keep df as is
						warning("Error in filtering network data: ", e$message)
					})
				}
				#calculate_substance_net(df, input$network_confounder, min_prevalence = input$subst_min_prev, min_co_prevalence = input$pair_min_prev)
				calculate_substance_net_meta(df, input$network_confounder, min_prevalence = input$subst_min_prev, min_co_prevalence = input$pair_min_prev)
			})

			output$substance_network_graph <- renderVisNetwork({
				net_data <- substance_network_models()
				signif_flag <- input$network_signif_only
				if (signif_flag) {
					validate(need(nrow(net_data$models[significance == "significant"]) > 0, "No network identified in current selection"))
				}
				network <- create_substance_network(net_data, or_range = input$network_or_filter, show_only_significant = signif_flag)
				network |>
					visNodes(font = list(
						size = 16,
						face = "arial",
						background = "rgba(255, 255, 255, 0.8)",
						strokeWidth = 0,
						align = "center"
					),
					color = list(background = "#e1e1e1", border = "#e1e1e1"),
					shadow = list(enabled = TRUE, size = 2)) |>
					visEdges(smooth = FALSE,
									 scaling = list(
									 	min = 2,
									 	max = 8
									 ),
									 hoverWidth = 2,
									 color = list(inherit = "both")) |>
					visOptions(highlightNearest = list(enabled = TRUE, algorithm = "all", degree = list(from = 0, to = 1))) |>
					visIgraphLayout() |>
					# visLayout(randomSeed = 123) |>
					# visPhysics(solver = "forceAtlas2Based",
					# 					 forceAtlas2Based = list(gravitationalConstant = -100),
					# 					 stabilization = FALSE) |>
					visInteraction(navigationButtons = TRUE)
			})

			output$substance_network_table <- renderDataTable({
				net_data <- substance_network_models()$models
				net_data <- net_data[,. (substance1,substance2,co_count,odds_ratio,ci_lower,ci_upper, p_value, p_adjusted)]
				datatable(net_data, rownames = FALSE) |> formatRound(c("odds_ratio","ci_lower", "ci_upper", "p_value", "p_adjusted"), digits = 4)
			})

		})
}