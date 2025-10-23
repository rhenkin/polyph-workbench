module_outcome_explorer_ui <- function(id) {
  ns <- NS(id)
  card(
    card_header("Prescriptions"),
    	selectInput(ns("freq_table_bnf"),
    						label = "Select BNF level for analysis:",
    						choices = c("BNF_Chapter", "BNF_Section", "BNF_Paragraph", "BNF_Chemical_Substance")),
      accordion(
        open = FALSE,
        accordion_panel(
          title = "Prevalence table",
          value = "prev_tables",
          navset_tab(
            nav_panel("Cohort prevalence",
            						dataTableOutput(ns("subst_freq_table"))
            					),
            nav_panel("Prevalence across demographics",
                      selectizeInput(ns("select_demog_freq_var"),
                      							 label = "Select variable:",
                      							 choices = c("Sex" = "sex","Ethnic group" = "eth_group", "IMD quintile" = "imd_quintile", "PP" = "pp_group", "# LTCs" = "mltc_group")),
                      dataTableOutput(ns("demog_subs_freq_table"))),
            nav_panel("Age-adjustment across demographics",
            					div(
            						"Odds ratios from logistic regression models for substances with prevalence >= 0.5% and unadjusted significance in the previous table, showing associations with patient characteristics (adjusted for age at outcome). Values >1 indicate higher prescription odds, <1 indicate lower odds. * denotes statistical significance (FDR-adjusted p<0.05)."
            					),
                      fluidRow(
                      	column(6,selectizeInput(ns("select_demog_model_var"), label = "Select variable:", choices = c("Sex" = "sex","Ethnic group" = "eth_group", "IMD quintile" = "imd_quintile", "PP" = "pp_group", "# LTCs" = "mltc_group"), options = list(onInitialize = I('function() { this.setValue(""); }'))),
            						column(6,actionButton(ns("run_subst_model"), "Run model")))
            						),
                      dataTableOutput(ns("freq_model_table")))
            )
        ),
        accordion_panel(
        	title ="Prescription prevalence per LTC",
        	value = "ltc_prev",
        		card(full_screen = TRUE,
        				 height = "60em",
        			fluidRow(column(6, uiOutput(ns("ltc_dropdown_ui"))),
        							 column(5, numericInput(ns("outcome_age_filter"), label = "Outcome age filter:", value = 100, min = 16, max = 100)),
        							 column(1, textOutput(ns("selected_pats")))),
        			dataTableOutput(ns("presc_by_ltc"))
        		)),
        # module_drug_network_ui(ns("drug_network_module")),
        accordion_panel(
        	title = "OR heatmap",
        	value = "or_heatmap_acc",
        	card(
        		full_screen = TRUE,
        		# selectInput(ns("heatmap_bnf"),
        		# 						label = "Select BNF level:",
        		# 						choices = c("BNF_Chapter", "BNF_Section", "BNF_Paragraph", "BNF_Chemical_Substance")),
        		virtualSelectInput(ns("corr_heatmap_strat"), label = "Select group to filter:", choices = NULL,  autoSelectFirstOption = FALSE, search = FALSE, disableOptionGroupCheckbox = TRUE, multiple = FALSE),
        		vegawidgetOutput(ns("corr_heatmap"))
        	)
        ),
        accordion_panel(
        	"Polypharmacy transition per substance",
        	div("Column indicates current number of prescriptions before new prescription is added"),
        	dataTableOutput(ns("transition_table"))
        ),

        # accordion_panel(
        #   "Drug-burden association", module_burden_model_ui(ns("model_module"))
        # ),
        # accordion_panel(
        #   "Drug-burden groups", module_burden_group_ui(ns("bgroup_module"))
        # ),
        accordion_panel(
          "Time to outcome", module_presc_tto_ui(ns("tto_module"))
        )
        # accordion_panel(
        # 	"Acute prescriptions", module_acute_presc_ui(ns("acute_presc_module"))
        # )
    )
  )
}

#' @export
module_outcome_explorer_server <- function(id, outcome_prescriptions, patient_data, ltc_data, selected_outcome, bnf_lookup, chapter_menu_data, pp_groups_data, acute_outcome_prescriptions) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    bnf_lookup <- bnf_lookup[,
    												 .SD[1,.(BNF_Chemical_Substance, BNF_Subparagraph, BNF_Paragraph, BNF_Section, BNF_Chapter)], BNF_Chemical_Substance] |> unique()

    subst_pp_df <- reactive({
      df <- copy(outcome_prescriptions())
      validate(need(nrow(df) > 0, "No valid patients found"))
      setkey(df, patid)
      df[, substance := bnf_lookup[match(df$substance, bnf_lookup$BNF_Chemical_Substance), get(input$freq_table_bnf)]]
      df <- merge(df, pp_groups_data(), by="patid")
      df <- merge(df, patient_data(), by="patid")
      df
    })

    heatmap_var_lookup <- reactiveVal()

    observe({
    	req(subst_pp_df())
    	df <- subst_pp_df()
    	heatmap_vars <- c("sex", "eth_group", "imd_quintile", "pp_group", "mltc_group")
    	var_labels <- c("Sex", "Ethnicity", "IMD quintile", "# PP", "# LTC")

    	choices <- setNames(
    		lapply(seq_along(heatmap_vars), function(i) {
    			var <- heatmap_vars[i]
    			unique_vals <- sort(unique(df[[var]]))

    			# Create named list: display_name = encoded_value
    			encoded_list <- setNames(
    				paste0(var, "#", unique_vals),  # values (encoded)
    				unique_vals                     # names (display labels)
    			)
    			return(encoded_list)
    		}),
    		var_labels  # group names
    	)

    	updateVirtualSelect("corr_heatmap_strat", choices = choices)
    })

    subst_frequency <- reactive({
      df <- copy(subst_pp_df())

      unique_patid <- df[,.(substance,patid)] |> unique()
      subst_freq <- unique_patid[,.N,.(substance)]
      total_patids <- uniqueN(df$patid)
      subst_freq[, pct_total := N/total_patids] #round(100*N/total_patids, digits = 2)]
      setorder(subst_freq, -N)
    })

    output$subst_freq_table <- renderDataTable({
      df <- copy(subst_frequency())
    	df[,.(substance, N, pct_total = round(pct_total*100, digits = 2))]
    }, rownames = FALSE)

    demog_subs_freq_df <- reactive({
    	req(input$select_demog_freq_var)

    	# Prepare data
    	#df <- merge(patient_data(), subst_pp_df(), by = "patid")
    	df <- copy(subst_pp_df())
    	#df[, substance := bnf_lookup[match(df$substance, bnf_lookup$BNF_Chemical_Substance), get(input$freq_table_bnf)]]
    	demog_var <- input$select_demog_freq_var

    	df <- df[substance %in% subst_frequency()[pct_total>0.005, substance]]

    	# Calculate statistics
    	cat_totals <- calculate_category_totals(df, demog_var)
    	df_stats <- calculate_demographic_stats(df, demog_var, subst_frequency(), "substance")

    	# Perform chi-square tests
    	chisq_tests <- perform_chisq_tests(df_stats, cat_totals, demog_var, "substance")

    	# Format and return results
    	format_results(df_stats, chisq_tests, cat_totals, demog_var, subst_frequency(), "substance")
    })

    output$demog_subs_freq_table <- renderDataTable({
    	req(demog_subs_freq_df())
    	df <- demog_subs_freq_df()
    	df[!is.na(padj), pvalue := as.double(sprintf("%.2e", pvalue))]
    	df[!is.na(padj), padj := as.double(sprintf("%.2e", padj))]
    }, rownames = FALSE)

    substance_model_data <- eventReactive(input$run_subst_model, {
	    req(demog_subs_freq_df())
	    df <- demog_subs_freq_df()
	    data_for_model <- copy(subst_pp_df())
	    #data_for_model[, substance := bnf_lookup[match(data_for_model$substance, bnf_lookup$BNF_Chemical_Substance), get(input$freq_table_bnf)]]
	    calculate_substance_model(data_for_model[substance %in% df[pvalue<0.05 & pct_total>0.005, substance][1:20]], patient_data(),
	    													c(input$select_demog_model_var, "outcome_age"))
    })
    output$freq_model_table <- renderDataTable({
			substance_model_data()
    }, rownames = FALSE)

    output$ltc_dropdown_ui <- renderUI({
    	ltcs <- ltc_data()
    	virtualSelectInput(ns("ltc_dropdown"),label = "Select 1 or more LTCs:", choices = with(chapter_menu_data[ltc %in% unique(ltcs$term)], split(ltc, body_system)), multiple = TRUE, search = TRUE)
    })

    age_filtered_prescriptions <- reactive({
    	subst_pp_df()[outcome_age <= input$outcome_age_filter*365.25]
    })

    output$selected_pats <- renderText({
    	dt <- age_filtered_prescriptions()
    	paste0("Patients within age range: ", prettyNum(uniqueN(dt$patid), big.mark = ","))
    })

  	output$presc_by_ltc <- renderDataTable({
			req(input$ltc_dropdown)
  		ltcs <- ltc_data()
  		patids <- unique(ltcs[term %in% input$ltc_dropdown, patid])
  		prescriptions <- age_filtered_prescriptions()

  		# Patients WITH the disease
  		presc_freq <- prescriptions[patid %in% patids,
  																list(N_with_disease = uniqueN(patid),
  																		 Prevalence = round(100*(uniqueN(patid)/length(patids)), digits=2),
  																		 `Median Duration (years)`= round(median(duration/365.2), digits = 2),
  																		 `IQR (Q1-Q3)` = paste0("(",
  																		 											 round(quantile(duration/365.2, 0.25, na.rm=TRUE), 2), " - ",
  																		 											 round(quantile(duration/365.2, 0.75, na.rm=TRUE), 2), ")")),
  																substance]

  		# Patients WITHOUT the disease
  		unselected_patids <- prescriptions[!patid %in% patids, uniqueN(patid)]
  		not_selected_freq <- prescriptions[!patid %in% patids,
  																			 list(N_without_disease = uniqueN(patid),
  																			 		 Prevalence_Unselected = round(100*(uniqueN(patid)/unselected_patids), digits=2),
  																			 		 `Median Duration unselected (years)`= round(median(duration/365.2), digits = 2),
  																			 		 `IQR unsel. (Q1-Q3)` = paste0("(",
  																			 		 											round(quantile(duration/365.2, 0.25, na.rm=TRUE), 2), " - ",
  																			 		 											round(quantile(duration/365.2, 0.75, na.rm=TRUE), 2), ")")),
  																			 substance]

  		result <- merge(presc_freq, not_selected_freq)

  		result <- result[Prevalence >= 0.005]

  		result[, `:=`(
  			total_with_disease = length(patids),
  			total_without_disease = unselected_patids,
  			Prevalence_Ratio = round(Prevalence / Prevalence_Unselected, digits = 2)
  		)]

  		result[is.infinite(Prevalence_Ratio), Prevalence_Ratio := 0]

  		# Calculate 95% CI for the ratio using log method
  		result[, `:=`(
  			p1 = N_with_disease / total_with_disease,
  			p2 = N_without_disease / total_without_disease
  		)]

  		result[, `:=`(
  			log_ratio = log(Prevalence_Ratio),
  			se_log_ratio = sqrt((1/N_with_disease) - (1/total_with_disease) +
  														(1/N_without_disease) - (1/total_without_disease))
  		)]

  		result[, `:=`(
  			CI_lower = round(exp(log_ratio - 1.96 * se_log_ratio), digits = 2),
  			CI_upper = round(exp(log_ratio + 1.96 * se_log_ratio), digits = 2)
  		)]

  		result[,
  			CI_95 := paste0("(", CI_lower, " - ", CI_upper, ")")
  		]

  		result[(CI_lower > 1.0 | CI_upper < 1.0), substance := paste0(substance, "*")]
  		# Clean up intermediate columns if desired
  		result[, c("p1", "p2", "log_ratio", "se_log_ratio", "CI_lower", "CI_upper") := NULL]

  		result[order(-Prevalence_Ratio), .(substance, Prevalence, Prevalence_Unselected, Prevalence_Ratio, CI_95, `Median Duration (years)`, `Median Duration unselected (years)`, `IQR (Q1-Q3)`, `IQR unsel. (Q1-Q3)`)]

  	}, rownames = FALSE)


  	output$corr_heatmap <- renderVegawidget({
  		df <- copy(subst_pp_df())
  		#df[, substance := bnf_lookup[match(df$substance, bnf_lookup$BNF_Chemical_Substance), get(input$heatmap_bnf)]]
  		df$yes <- 1

  		selected <- input$corr_heatmap_strat
  		if (selected != "") {
  			parts <- strsplit(selected, "#")[[1]]
  			column_name <- parts[1]
  			filter_value <- parts[2]

  			# Filter your data.table
  			df <- df[get(column_name) == filter_value]
  		}

  		mat <- dcast(df, patid~substance, value.var = "yes", fun.aggregate = is.numeric, fill = 0)
  		prev <- as.list(mat[, lapply(.SD[, -1], mean)])
  		keep_conditions <- names(prev[which(prev > 0.005)])





  		#browser()

  		# pairs_dt <- data.table(expand.grid(drug1 = keep_conditions, drug2 = keep_conditions,
  		# 																	 stringsAsFactors = FALSE))
  		# pairs_dt <- pairs_dt[drug1 < drug2]
  		#
  		# results <- calc_or_vec(pairs_dt$drug1, pairs_dt$drug2, mat)
  		#
  		# # Results will be a matrix where each column is (or, ci_lower, ci_upper)
  		# pairs_dt[, `:=`(
  		# 	or = results[1, ],
  		# 	ci_lower = results[2, ],
  		# 	ci_upper = results[3, ]
  		# )]
  		#
  		# pairs_dt <- pairs_dt[!is.na(or) & or>1 & (ci_lower > 1.0 | ci_upper < 1.0)]
  		# pairs_symmetric <- rbind(
  		# 	pairs_dt,
  		# 	pairs_dt[, .(drug1 = drug2, drug2 = drug1, or, ci_lower, ci_upper)]
  		# )
  		ors <- calc_all_ors_vectorized(as.matrix(mat[,..keep_conditions]))
  		pairs_dt <- data.table(expand.grid(drug1 = keep_conditions, drug2 = keep_conditions,
  		 																	 stringsAsFactors = FALSE))
  		pairs_dt[, or := as.vector(ors$or)]
  		pairs_dt[, ci_lower := as.vector(ors$ci_lower)]
  		pairs_dt[, ci_upper := as.vector(ors$ci_upper)]
  		pairs_dt <- pairs_dt[!is.na(or) & or>1 & (ci_lower > 1.0 | ci_upper < 1.0)]

  		# pairs_symmetric <- rbind(
  		# 	pairs_dt,
  		# 	pairs_dt[, .(drug1 = drug2, drug2 = drug1, or, ci_lower, ci_upper)]
  		# )
  		or_heatmap(pairs_dt, ors$or) |> as_vegaspec()

			# Partial corrleation if needed
			# pcor_result <- pcor(mat[,..keep_conditions])
			# pcor_matrix <- pcor_result$estimate
			# p_values <- pcor_result$p.value
			# rownames(pcor_matrix) <- colnames(pcor_matrix) <- keep_conditions
			# rownames(p_values) <- colnames(p_values) <- keep_conditions
			#
			# p_adjusted <- p.adjust(as.vector(p_values), method = "bonferroni")
			# dim(p_adjusted) <- dim(p_values)  # Restore matrix structure
			# dimnames(p_adjusted) <- dimnames(p_values)
			#
			# # Apply standard significance threshold to FDR-corrected p-values
			# significance_threshold <- 0.05
			# significant_positive <- p_adjusted < significance_threshold & pcor_matrix > 0
			# cat("Keeping", sum(significant_positive, na.rm = TRUE), "out of", length(p_adjusted), "correlations\n")
			# pcor_matrix[p_adjusted >= significance_threshold | pcor_matrix <= 0] <- 0
			#
			# # Ensure diagonal is zero
			# diag(pcor_matrix) <- 0
			# cor_heatmap(pcor_matrix) |> as_vegaspec()

  	})

    output$transition_table <- renderDataTable({
    	med_data <- copy(subst_pp_df())
    	#med_data[, substance := bnf_lookup[match(med_data$substance, bnf_lookup$BNF_Chemical_Substance), get(input$freq_table_bnf)]]
    	setorder(med_data, patid, start_date)


    	level_profiles <- analyze_medications(med_data, mode = "level")

    	# Combined trajectory and transitions in one step
    	# transitions <- med_data[, {
    	# 	# Create sequence number for each medication
    	# 	med_seq <- 1:.N
    	# 	# Previous level starts at 0, then takes previous med_seq values
    	# 	prev_level <- c(0, med_seq[-.N])
    	#
    	# 	# Calculate categories directly
    	# 	prev_level_cat <- poly_level_category(prev_level)
    	# 	curr_level_cat <- poly_level_category(med_seq)
    	#
    	# 	# Flag transitions between categories
    	# 	is_transition <- prev_level_cat != curr_level_cat
    	#
    	# 	list(
    	# 		patid = patid,
    	# 		substance = substance,
    	# 		med_seq = med_seq,
    	# 		poly_level = med_seq,
    	# 		poly_level_cat = curr_level_cat,
    	# 		prev_poly_level = prev_level,
    	# 		prev_poly_level_cat = prev_level_cat,
    	# 		is_transition = is_transition,
    	# 		days_to_outcome = eventdate-start_date
    	# 	)
    	# }, by = patid]
    	#
    	# # Calculate medication prevalence and filter in one step
    	# frequent_meds <- transitions[, .(
    	# 	count = .N,
    	# 	prevalence = .N / nrow(transitions) * 100
    	# ), by = substance][prevalence >= 1, substance]
    	#
    	# # Transition analysis for frequent medications only
    	# transition_substances <- transitions[is_transition & substance %in% frequent_meds, .(
    	# 	count = .N
    	# ), by = .(substance, prev_poly_level_cat, poly_level_cat)]
    	#
    	# # Create the substance profiles
    	# substance_profiles <- dcast(
    	# 	transition_substances,
    	# 	substance ~ paste(prev_poly_level_cat, "→", poly_level_cat),
    	# 	value.var = "count",
    	# 	fill = 0
    	# )
    	#
    	# # Calculate total for each substance
    	# substance_profiles[, total := rowSums(.SD), .SDcols = patterns("\\[")]
    	#
    	# # Convert to proportions
    	# prop_cols <- names(substance_profiles)[grep("\\[", names(substance_profiles))]
    	# substance_profiles[, (prop_cols) := lapply(.SD, function(x) round((x / total) * 100, digits = 2)),
    	# 									 .SDcols = prop_cols]
    	#
    	#
    	# substance_profiles[order(-total)]
    }, rownames = FALSE)

    # module_burden_model_server("model_module", subst_pp_df, patient_data, subst_frequency)
    # module_burden_group_server("bgroup_module", subst_pp_df)
    # module_drug_network_server("drug_network_module", subst_pp_df, bnf_lookup)
    module_presc_tto_server("tto_module", subst_pp_df, subst_frequency)
    #module_acute_presc_server("acute_presc_module", acute_outcome_prescriptions)

  })
}