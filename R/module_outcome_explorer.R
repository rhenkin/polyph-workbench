module_outcome_explorer_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Prescriptions"),
      accordion(
        open = FALSE,
        accordion_panel(
          title = "Frequency table",
          navset_tab(
            nav_panel("Cohort frequency", dataTableOutput(ns("subst_freq_table"))),
            nav_panel("Frequency across demographics",
                      selectizeInput(ns("select_demog_freq_var"),
                      							 label = "Select variable:",
                      							 choices = c("sex", "eth_group",  "imd_quintile", "pp_group")),
                      dataTableOutput(ns("demog_subs_freq_table"))),
            nav_panel("Age-adjustment across demographics",
            					div(
            						"Odds ratios from logistic regression models for substances with prevalence >= 1% and unadjusted significance in the previous table, showing associations with patient characteristics (adjusted for age at outcome). Values >1 indicate higher prescription odds, <1 indicate lower odds. * denotes statistical significance (FDR-adjusted p<0.05)."
            					),
                      fluidRow(
                      	column(6,selectizeInput(ns("select_demog_model_var"), label = "Select variable:", choices = c("sex", "eth_group",  "imd_quintile", "pp_group"), options = list(onInitialize = I('function() { this.setValue(""); }'))),
            						column(6,actionButton(ns("run_subst_model"), "Run model")))
            						),
                      dataTableOutput(ns("freq_model_table")))
            )
        ),
        # accordion_panel(
        #   "Co-occurrence network", visNetworkOutput(ns("substance_network"), height = "600px")
        # ),
        # accordion_panel(
        #   "Drug-burden association", module_burden_model_ui(ns("model_module"))
        # ),
        # accordion_panel(
        #   "Drug-burden groups", module_burden_group_ui(ns("bgroup_module"))
        # ),
        accordion_panel(
          "Time to outcome", module_presc_tto_ui(ns("tto_module"))
        ),
        accordion_panel(
        	"Non-repeat prescriptions", module_nr_presc_ui(ns("nrpresc_module"))
        )
    )
  )
}

#' @export
module_outcome_explorer_server <- function(id, outcome_prescriptions, patient_data, selected_outcome, bnf_lookup, pp_groups_data) {
  moduleServer(id, function(input, output, session) {
    ns <- NS(id)

    subst_pp_df <- reactive({
      df <- outcome_prescriptions()
      validate(need(nrow(df) > 0, "No valid patients found"))
      # df <- merge(df, bnf_lookup[,.(BNF_Presentation_Code, BNF_Chemical_Substance)],
      #             by.x = "substance",
      #             by.y = "BNF_Presentation_Code", all.x = TRUE)
      setkey(df, patid)
      df <- merge(df, pp_groups_data(), by="patid")
      df
    })

    subst_frequency <- reactive({
      df <- subst_pp_df()
      unique_patid <- df[,.(substance,patid)] |> unique()
      subst_freq <- unique_patid[,.N,.(substance)]
      total_patids <- length(unique(df$patid))
      subst_freq[, pct_total := signif(N/total_patids, digits = 2)]
      setorder(subst_freq,-N)
    })

    output$subst_freq_table <- renderDataTable({
      subst_frequency()[,.(substance, N, pct_total)]
    }, rownames = FALSE)


#     output$substance_network <- renderVisNetwork({
#     	selected_subst <- unique(subst_frequency()[pct_total>0.01,substance])
#     	df <- outcome_prescriptions()[substance %in% selected_subst]
# 			browser()
# 			total_patients <- as.numeric(uniqueN(df$patid))
# 			subst_freq <- df[, .(
# 				frequency = as.numeric(.N),
# 				prevalence = as.numeric(.N)/uniqueN(df$patid)
# 			), by = substance]
#
# 			# Calculate co-occurrences
# 			cooccurrence <- df[, .(substance1 = substance, patid)]
# 			setkey(cooccurrence, patid)
#
# 			# Step 1: Calculate statistics for each pair
# 			pairs <- cooccurrence[cooccurrence,
# 														.(substance1 = i.substance1, substance2 = x.substance1),
# 														on = "patid",
# 														allow.cartesian = TRUE
# 			][substance1 < substance2
# 			][, .(cij = as.numeric(.N)), by = .(substance1, substance2)]
#
# 			# Add individual disease counts and prevalence
# 			result <- pairs[subst_freq, on = .(substance1 = substance)
# 			][, `:=`(
# 				ci = as.numeric(frequency),
# 				prev1 = prevalence
# 			)][subst_freq, on = .(substance2 = substance)
# 			][, `:=`(
# 				cj = as.numeric(frequency),
# 				prev2 = prevalence
# 			)]
#
# 			# Calculate expected co-occurrence under independence
# 			result[, expected_cij := total_patients * prev1 * prev2]
#
# 			# Calculate relative risk and lift
# 			result[, `:=`(
# 				relative_risk = (cij/total_patients) / (prev1 * prev2),
# 				lift = cij / expected_cij
# 			)]
#
# 			# Calculate phi coefficient and t-value using numeric types
# 			result[, `:=`(
# 				numerator = cij * total_patients - ci * cj,
# 				denominator = sqrt(as.numeric(ci) * cj * (total_patients - ci) * (total_patients - cj))
# 			)]
#
# 			# Complete the calculations
# 			result[, `:=` (
# 						phi = numerator/denominator,
# 				    expected = total_patients * ci * cj,
# 				    ochiai = cij / sqrt(ci * cj),
# 				    observed_to_expected = cij / (total_patients * ci * cj),
# 				    relative_difference = (cij - (total_patients * cij * cj)) / (total_patients * ci * cj)
# 				  )]
#
# 			# Clean up intermediate columns
# 			result[, c("numerator", "denominator") := NULL]
# 		result <- result[!is.na(phi)]
# 			edges_df <- data.frame(
# 				from = result$substance1,
# 				to = result$substance2,
# 				weight = result$phi+1
# 			)
#
# 			# Create igraph object
# 			g <- graph_from_data_frame(
# 				d = edges_df,
# 				directed = FALSE
# 			)
#
# 			# Calculate node metrics
# 			V(g)$degree <- degree(g)
# 			V(g)$betweenness <- betweenness(g)
# 			V(g)$community <- as.numeric(membership(cluster_louvain(g)))
#
# 			nodes <- data.frame(
# 				id = as.character(V(g)$name),
# 				label = as.character(V(g)$name),
# 				value = as.numeric(V(g)$degree),
# 				group = as.numeric(V(g)$community),
# 				title = sprintf("Degree: %d<br>Betweenness: %.2f",
# 												as.numeric(V(g)$degree),
# 												as.numeric(V(g)$betweenness)),
# 				stringsAsFactors = FALSE
# 			)
#
# 			# Prepare edges with explicit weight handling
# 			edge_df <- data.frame(
# 				from = as.character(get.edgelist(g)[,1]),
# 				to = as.character(get.edgelist(g)[,2])
# 			)
#
# 			# Add weights if they exist
# 			if ("weight" %in% edge.attributes(g)) {
# 				edge_df$value <- E(g)$weight
# 				edge_df$title <- sprintf("Weight: %.2f", E(g)$weight)
# 			} else {
# 				edge_df$value <- 1
# 				edge_df$title <- "Weight: 1"
# 			}
#
# 			# Create visualization
# 			visNetwork(nodes, edge_df) %>%
# 				visIgraphLayout() %>%
# 				# visPhysics(solver = "forceAtlas2Based",
# 				#            forceAtlas2Based = list(gravitationalConstant = -100)) %>%
# 				visOptions(highlightNearest = TRUE,
# 									 nodesIdSelection = TRUE) %>%
# 				visLayout(randomSeed = 123) %>%
# 				visInteraction(navigationButtons = TRUE)  %>%
# 				visEdges(smooth = FALSE)  %>%
# 				visNodes(scaling = list(min = 5, max = 30), font = list(size = 30))
#
#     })

    demog_subs_freq_df <- reactive({
    	req(input$select_demog_freq_var)

    	# Prepare data
    	df <- merge(patient_data(), subst_pp_df(), by = "patid")
    	demog_var <- input$select_demog_freq_var

    	# Calculate statistics
    	cat_totals <- calculate_category_totals(df, demog_var)
    	df_stats <- calculate_demographic_stats(df, demog_var, subst_frequency())

    	# Perform chi-square tests
    	chisq_tests <- perform_chisq_tests(df_stats, cat_totals, demog_var)

    	# Format and return results
    	format_results(df_stats, chisq_tests, cat_totals, demog_var, subst_frequency())
    })

    output$demog_subs_freq_table <- renderDataTable({
    	req(demog_subs_freq_df())
    	df <- demog_subs_freq_df()
    	df[!is.na(padj), pvalue := as.double(sprintf("%.3f", pvalue))]
    	df[!is.na(padj), padj := as.double(sprintf("%.3f", padj))]
    }, rownames = FALSE)

    substance_model_data <- eventReactive(input$run_subst_model, {
    req(demog_subs_freq_df())
    df <- demog_subs_freq_df()
    calculate_substance_model(subst_pp_df()[substance %in% df[pvalue<0.05 & pct_total>0.01, substance][1:10]], patient_data(),
    													c(input$select_demog_model_var, "outcome_age"))
    })
    output$freq_model_table <- renderDataTable({
			substance_model_data()
    })

    # module_burden_model_server("model_module", subst_pp_df, patient_data, subst_frequency)
    # module_burden_group_server("bgroup_module", subst_pp_df)
    # module_drug_survival_server("survival_module", subst_pp_df)
  })
}