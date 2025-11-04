bnf_table_ui <- function(id) {
	ns <- NS(id)
	actionButton(ns("show_bnf_table"), "BNF Table")
}

bnf_table_server <- function(id, bnf_lookup) {
	moduleServer(id, function(input, output, session) {
		ns <- session$ns

		bnf_to_display <- reactiveVal()
		observe({
			unique_bnf_table <- unique(bnf_lookup[BNF_Chemical_Substance != "",.(BNF_Chapter, BNF_Section, BNF_Paragraph, BNF_Chemical_Substance)])
			bnf_to_display(unique_bnf_table)
		}, suspended = TRUE)


		observeEvent(input$show_bnf_table, {

			# Get initial chapters
			chapters <- unique(bnf_to_display()$BNF_Chapter)
			chapters <- sort(chapters[!is.na(chapters)])

			showModal(modalDialog(
				size = "xl", # Made larger to accommodate the list boxes
				title = "BNF Table",
				card(
					card_header("BNF Browser"),
					card_body(
						fluidRow(
							# Three HTML select inputs side by side
							column(
								3,
								h5("BNF Chapter"),
								selectInput(
									inputId = ns("bnf_chapter"),
									label = NULL,
									choices = chapters,
									selected = NULL,
									multiple = TRUE,
									selectize = FALSE, # This makes it a classic HTML select
									size = 20 # Number of visible options
								)
							),
							column(
								3,
								div(
									style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
									h5("BNF Section", style = "margin: 0;"),
									actionButton(
										ns("select_all_sections"),
										"Select All",
										size = "sm",
										style = "padding: 2px 8px; font-size: 0.75rem;"
									)
								),
								selectInput(
									inputId = ns("bnf_section"),
									label = NULL,
									choices = NULL,
									selected = NULL,
									multiple = TRUE,
									selectize = FALSE,
									size = 20
								)
							),
							column(
								3,
								div(
									style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
									h5("BNF Paragraph", style = "margin: 0;"),
									actionButton(
										ns("select_all_paragraphs"),
										"Select All",
										size = "sm",
										style = "padding: 2px 8px; font-size: 0.75rem;"
									)
								),
								selectInput(
									inputId = ns("bnf_paragraph"),
									label = NULL,
									choices = NULL,
									selected = NULL,
									multiple = TRUE,
									selectize = FALSE,
									size = 20
								)
							),
							# Chemical substances list
							column(
								3,
								div(
									style = "display: flex; justify-content: space-between; align-items: center; margin-bottom: 5px;",
									h5("Chemical Substances", style = "margin: 0;"),
									actionButton(
										ns("copy_substances"),
										"Copy All",
										size = "sm",
										style = "padding: 2px 8px; font-size: 0.75rem;"
									)
								),
								uiOutput(ns("substances_display"))
							)
						),

						# Custom CSS for better styling
						tags$style(
							HTML(paste0("
								#", ns("bnf_chapter"), ", #", ns("bnf_section"), ", #", ns("bnf_paragraph"), " {
									height: 300px !important;
									font-size: 0.875rem;
									border: 1px solid #dee2e6;
								}

								#", ns("bnf_chapter"), " option:checked, #", ns("bnf_section"), " option:checked, #", ns("bnf_paragraph"), " option:checked {
									background-color: #0d6efd;
									color: white;
								}

								/* Improve spacing */
								.form-group {
									margin-bottom: 0;
								}
							"))
						),
						tags$script(HTML("
    Shiny.addCustomMessageHandler('copyToClipboard', function(message) {
        navigator.clipboard.writeText(message.text).then(function() {
            console.log('Copied to clipboard successfully');
        }).catch(function(err) {
            console.error('Failed to copy to clipboard: ', err);
            // Fallback for older browsers
            var textArea = document.createElement('textarea');
            textArea.value = message.text;
            document.body.appendChild(textArea);
            textArea.select();
            document.execCommand('copy');
            document.body.removeChild(textArea);
        });
    });
"))
					)
				),
				easyClose = TRUE,
				footer = NULL
			))
		}, suspended = TRUE)

		# Select all sections button
		observeEvent(input$select_all_sections, {
			if (!is.null(input$bnf_section)) {
				# Get all available choices for sections
				current_choices <- names(input$bnf_section)
				if (is.null(current_choices)) {
					# If no names, get from the select input choices
					session_data <- session$userData
					if (exists("current_section_choices", envir = session_data)) {
						current_choices <- session_data$current_section_choices
					}
				}

				# Get choices from the current selectInput
				filtered_data <- bnf_to_display()[BNF_Chapter %in% input$bnf_chapter]
				sections <- unique(filtered_data$BNF_Section)
				sections <- sort(sections[!is.na(sections)])

				updateSelectInput(
					session,
					"bnf_section",
					selected = sections
				)
			}
		})

		# Select all paragraphs button
		observeEvent(input$select_all_paragraphs, {
			if (!is.null(input$bnf_chapter) && !is.null(input$bnf_section) &&
					length(input$bnf_chapter) > 0 && length(input$bnf_section) > 0) {

				filtered_data <- bnf_to_display()[
					BNF_Chapter %in% input$bnf_chapter &
						BNF_Section %in% input$bnf_section
				]
				paragraphs <- unique(filtered_data$BNF_Paragraph)
				paragraphs <- sort(paragraphs[!is.na(paragraphs)])

				updateSelectInput(
					session,
					"bnf_paragraph",
					selected = paragraphs
				)
			}
		})

		# Update sections when chapters are selected
		observeEvent(input$bnf_chapter, {
			if (!is.null(input$bnf_chapter) && length(input$bnf_chapter) > 0) {
				filtered_data <- bnf_to_display()[BNF_Chapter %in% input$bnf_chapter]
				sections <- unique(filtered_data$BNF_Section)
				sections <- sort(sections[!is.na(sections)])

				updateSelectInput(
					session,
					"bnf_section",
					choices = sections,
					selected = NULL
				)
			} else {
				updateSelectInput(session, "bnf_section", choices = NULL, selected = NULL)
			}

			# Clear downstream selections
			updateSelectInput(session, "bnf_paragraph", choices = NULL, selected = NULL)
		})

		# Update paragraphs when sections are selected
		observeEvent(input$bnf_section, {
			if (!is.null(input$bnf_chapter) && !is.null(input$bnf_section) &&
					length(input$bnf_chapter) > 0 && length(input$bnf_section) > 0) {

				filtered_data <- bnf_to_display()[
					BNF_Chapter %in% input$bnf_chapter &
						BNF_Section %in% input$bnf_section
				]
				paragraphs <- unique(filtered_data$BNF_Paragraph)
				paragraphs <- sort(paragraphs[!is.na(paragraphs)])

				updateSelectInput(
					session,
					"bnf_paragraph",
					choices = paragraphs,
					selected = NULL
				)
			} else {
				updateSelectInput(session, "bnf_paragraph", choices = NULL, selected = NULL)
			}
		}, suspended = TRUE)

		# Render substances display
		output$substances_display <- renderUI({
			req(input$bnf_chapter, input$bnf_section, input$bnf_paragraph)
			if (!is.null(input$bnf_chapter) && !is.null(input$bnf_section) && !is.null(input$bnf_paragraph) &&
					length(input$bnf_chapter) > 0 && length(input$bnf_section) > 0 && length(input$bnf_paragraph) > 0) {

				filtered_data <- bnf_to_display()[
					BNF_Chapter %in% input$bnf_chapter &
						BNF_Section %in% input$bnf_section &
						BNF_Paragraph %in% input$bnf_paragraph
				]
				substances <- unique(filtered_data$BNF_Chemical_Substance)
				substances <- sort(substances[!is.na(substances)])

				if (length(substances) > 0) {
					list_items <- lapply(substances, function(substance) {
						div(
							class = "mb-1 p-2",
							style = "border-bottom: 1px solid #dee2e6; cursor: default;",
							substance
						)
					})

					div(
						style = "height: inherit; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 0.375rem; padding: 0.5rem; background-color: #f8f9fa; font-size: 0.875rem;",
						list_items
					)
				} else {
					div(
						style = "height: 300px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 0.375rem; padding: 0.5rem; background-color: #f8f9fa; font-size: 0.875rem;",
						div(
							class = "text-muted",
							style = "font-style: italic;",
							"No substances found for this selection"
						)
					)
				}
			} else {
				div(
					style = "height: 300px; overflow-y: auto; border: 1px solid #dee2e6; border-radius: 0.375rem; padding: 0.5rem; background-color: #f8f9fa; font-size: 0.875rem;",
					div(
						class = "text-muted",
						style = "font-style: italic;",
						"Select items above to see substances..."
					)
				)
			}
		})


		# Copy substances to clipboard
		observeEvent(input$copy_substances, {
			if (!is.null(input$bnf_chapter) && !is.null(input$bnf_section) && !is.null(input$bnf_paragraph) &&
					length(input$bnf_chapter) > 0 && length(input$bnf_section) > 0 && length(input$bnf_paragraph) > 0) {

				filtered_data <- bnf_to_display()[
					BNF_Chapter %in% input$bnf_chapter &
						BNF_Section %in% input$bnf_section &
						BNF_Paragraph %in% input$bnf_paragraph
				]
				substances <- unique(filtered_data$BNF_Chemical_Substance)
				substances <- sort(substances[!is.na(substances)])

				if (length(substances) > 0) {
					# Create text to copy (one substance per line)
					copy_text <- paste(substances, collapse = ";")

					# Use JavaScript to copy to clipboard
					session$sendCustomMessage(
						type = "copyToClipboard",
						message = list(text = copy_text)
					)

					# Show confirmation
					showNotification("Chemical substances copied to clipboard!", type = "message", duration = 2)
				}
			}
		})

		# Add JavaScript to handle the custom message
		session$onSessionEnded(function() {
			insertUI(
				selector = "head",
				ui = tags$script(HTML("
					Shiny.addCustomMessageHandler('updateSubstances', function(message) {
						var element = document.querySelector(message.selector);
						if (element) {
							element.innerHTML = message.content;
						}
					});
				")),
				immediate = TRUE
			)
		})
	})
}