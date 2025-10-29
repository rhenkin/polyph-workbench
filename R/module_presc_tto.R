module_presc_tto_ui <- function(id) {
  ns <- NS(id)
  # card(full_screen = TRUE,
  # 		 min_height = 600,
  	navset_tab(
  		nav_panel(title = "TTO prevalence",
  							div("Prevalence (in %) of timings from prescripion start to outcome. Only the last prescription before outcome is counted. Time is measured in days"),
  							dataTableOutput(ns("tto_table"))
  		),
  		nav_panel(title = "Stratified prevalence",
  							div("Prevalence (in %) of timings from prescripion start to outcome. Only last prescription before outcome is counted. Time is measured in days"),
  							fluidRow(
  								column(2, selectizeInput(ns("strat_tto_variable"), "Select a variable", choices = c("Ethnic group" = "eth_group", "Sex" = "sex", "IMD quintile" = "imd_quintile", "PP group" = "pp_group", "# LTCs" = "mltc_group")),),
  								column(2, selectizeInput(ns("strat_tto_values"), "Select a subgroup:", choices = NULL))
  							),
  							dataTableOutput(ns("strat_tto_table"))
  							)
  	)
  # )
}

module_presc_tto_server <- function(id, subst_pp_df, subst_freq) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

    last_presc_info <- function(df) {
    	df[,list(tto = eventdate - start_date)]
    }

    calculate_analysis_data <- function(df) {    	

    	drug_analysis <- df[df[, .I[which.max(start_date)], by = patid]$V1][
    		, `:=`(tto = eventdate - start_date)][
    			, `:=`(
    				within_1m = tto <= 30,
    				within_3m = tto %between% c(31, 90),
    				within_6m = tto %between% c(91, 180),
    				within_1y = tto %between% c(181, 365),
    				after_1y = tto > 365
    			)][
    				, .(
    					N = .N,
    					`< 30d` = round(sum(within_1m, na.rm = TRUE) / .N * 100, 2),
    					`< 3m` = round(sum(within_3m, na.rm = TRUE) / .N * 100, 2),
    					`< 6m` = round(sum(within_6m, na.rm = TRUE) / .N * 100, 2),
    					`< 1y` = round(sum(within_1y, na.rm = TRUE) / .N * 100, 2),
    					`1+y` = round(sum(after_1y, na.rm = TRUE) / .N * 100, 2),
    					`Median TTO` = as.double(median(tto, na.rm = TRUE)),
    					`Mean TTO` = round(as.double(mean(tto, na.rm = TRUE)), 2)
    				), by = substance]

    	return(drug_analysis)
    }

    output$tto_table <- renderDataTable({
    	req(subst_pp_df())
    	df <- copy(subst_pp_df())
    	subset_df <- df[substance %in% subst_freq()[pct_total>0.005, substance]]
    	setkey(subset_df, patid, substance)

    	calculate_analysis_data(subset_df)

    }, rownames = FALSE)

    observe({
    	req(input$strat_tto_variable)
    	df <- copy(subst_pp_df())
    	updateSelectizeInput(session, "strat_tto_values", choices = sort(unique(df[[input$strat_tto_variable]])))
    })

    output$strat_tto_table <- renderDataTable({
    	req(subst_pp_df())
    	req(input$strat_tto_values)
    	df <- copy(subst_pp_df())
    	subset_df <- df[substance %in% subst_freq()[pct_total>0.005, substance]]
    	subset_df <- subset_df[get(input$strat_tto_variable)==input$strat_tto_values]

    	calculate_analysis_data(subset_df)

    }, rownames = FALSE)

  })
}