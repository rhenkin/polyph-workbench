module_acute_presc_ui <- function(id) {
	ns <- NS(id)
	card(
		dataTableOutput(ns("acute_subst_freq_table"))
	)
}

module_acute_presc_server <-  function(id, acute_prescriptions_df) {
	moduleServer(id, function(input, output, session) {
		ns <- NS(id)

		subst_freq_df <- reactive({
			req(acute_prescriptions_df())
			df <- acute_prescriptions_df()
			df[,list(N = uniqueN(patid), "Mean time-to-outcome (days)" = mean(eventdate-start_date)),substance][order(-`Mean time-to-outcome (days)`)]
		})

		output$acute_subst_freq_table <- renderDataTable({
			req(subst_freq_df())
			subst_freq_df()
		}, rownames = FALSE)

	})
}