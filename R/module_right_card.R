right_card_ui <- function(id) {
  ns <- NS(id)
  # card(card_header("Drug filtering"),
  list(
  	div("Prescription filters can be defined before or after starting an analysis."),
       tagList(
         tags$style(HTML(sprintf("
      #%s .excluded-item, #%s .included-item {
        display: flex;
        justify-content: space-between;
        align-items: center;
        padding: 5px;
        margin: 2px 0;
        background-color: #f8f9fa;
        border-radius: 4px;
      }
      #%s .delete-btn {
        color: red;
        cursor: pointer;
        padding: 0 5px;
      }
      #%s .excluded-list, #%s .included-list {
        max-width: 400px;
      }
      #%s .included-item {
        background-color: #e3f2fd;
      }
      #%s .select-container {
        display: flex;
        align-items: last baseline;
        gap: 10px;
        margin-bottom: 10px;
      }
      #%s .action-btn {
        padding: 0px 10px;
        font-size: 18px;
        cursor: pointer;
      }
    ", ns("bnf-module"), ns("bnf-module"), ns("bnf-module"),
                                 ns("bnf-module"), ns("bnf-module"), ns("bnf-module"),
                                 ns("bnf-module"), ns("bnf-module")
         ))),
         div(
           id = ns("bnf-module"),
           numericInput(ns("polyph_threshold"),
           						 "Polypharmacy threshold:",
           						 value = 2, min = 2, max = 10),

           numericInput(ns("minimum_cp_duration"),
                        "Minimum treatment duration:",
                        value = 0),

           numericInput(ns("earliest_treatment_end"),
           						 "Earliest treatment end until outcome (days):",
           						 value = 84),

           # Chapter selection
           bnf_table_ui(ns("bnf_table_module")),
           div("Choose BNF levels for inclusion or exclusion"),
           div(class = "select-container",
               actionButton(ns("add_chapter"), "+", class = "action-btn"),
               selectInput(ns("chapter"),
                           "Chapter",
                           choices = sort(unique(bnf_lookup$BNF_Chapter))),
               actionButton(ns("exclude_chapter"), "-", class = "action-btn")
           ),

           # Section selection
           div(class = "select-container",
               actionButton(ns("add_section"), "+", class = "action-btn"),
               selectInput(ns("section"),
                           "Section",
                           choices = NULL),
               actionButton(ns("exclude_section"), "-", class = "action-btn")
           ),

           # Paragraph selection
           div(class = "select-container",
               actionButton(ns("add_paragraph"), "+", class = "action-btn"),
               selectInput(ns("paragraph"),
                           "Paragraph",
                           choices = NULL),
               actionButton(ns("exclude_paragraph"), "-", class = "action-btn")
           ),

           # Subparagraph selection
           div(class = "select-container",
               actionButton(ns("add_subparagraph"), "+", class = "action-btn"),
               selectInput(ns("subparagraph"),
                           "Subparagraph",
                           choices = NULL),
               actionButton(ns("exclude_subparagraph"), "-", class = "action-btn")
           ),

           hr(),

           span("Included BNFs:"),
           div(id = ns("included-list"), class = "included-list"),

           hr(),

           span("Excluded BNFs:"),
           div(id = ns("excluded-list"), class = "excluded-list")
         )
       )
  )
}

right_card_server <- function(id) {
  moduleServer(id, function(input, output, session) {
    ns <- session$ns

  #   lookup_for_modal <- unique(bnf_lookup[,.(BNF_Chapter, BNF_Section, BNF_Paragraph, BNF_Chemical_Substance)])
  #
  #   observeEvent(input$show_bnf_table, {
  #   	showModal(modalDialog(size = "l",
  #   		title = "BNF Table",
  # 			dataTableOutput(ns("bnf_table")),
  #   		easyClose = TRUE,
  #   		footer = NULL
  #   	))
  #   })
  #
  #   output$bnf_table <- renderDT({
  #   	datatable(lookup_for_modal,
  #   						rownames = FALSE, filter = "top")
  #   })

    # Reactive values to store included and excluded BNFs
    included_bnfs <- reactiveVal(list())
    excluded_bnfs <- reactiveVal(default_excluded_bnf)

    # Update Chapter choices on initialization
    # observe({
    #   chapters <- sort(unique(bnf_lookup$BNF_Chapter))
    #   updateSelectInput(session, "chapter", choices = chapters)
    # })

    # Update Section based on Chapter
    observe({
      req(input$chapter)
      sections <- sort(unique(bnf_lookup$BNF_Section[bnf_lookup$BNF_Chapter == input$chapter]))
      updateSelectInput(session, "section", choices = sections)
    })

    # Update Paragraph based on Section
    observe({
      req(input$section)
      paragraphs <- sort(unique(bnf_lookup$BNF_Paragraph[bnf_lookup$BNF_Section == input$section]))
      updateSelectInput(session, "paragraph", choices = paragraphs)
    })

    # Update Subparagraph based on Paragraph
    observe({
      req(input$paragraph)
      subparagraphs <- sort(unique(bnf_lookup$BNF_Subparagraph[bnf_lookup$BNF_Paragraph == input$paragraph]))
      updateSelectInput(session, "subparagraph", choices = subparagraphs)
    })

    # Helper function to add items to a list
    add_to_list <- function(category, value, current_list) {
      if (!is.null(value)) {
        current <- current_list()
        new_item <- list(
          category = category,
          value = value
        )
        # Add only if unique
        exists <- FALSE
        for (existing in current) {
          if (identical(new_item, existing)) {
            exists <- TRUE
            break
          }
        }
        if (!exists) {
          current[[length(current) + 1]] <- new_item
          current_list(current)
        }
      }
    }

    # Add handlers for each level
    observeEvent(input$add_chapter, {
      add_to_list("Chapter", input$chapter, included_bnfs)
    })

    observeEvent(input$add_section, {
      add_to_list("Section", input$section, included_bnfs)
    })

    observeEvent(input$add_paragraph, {
      add_to_list("Paragraph", input$paragraph, included_bnfs)
    })

    observeEvent(input$add_subparagraph, {
      add_to_list("Subparagraph", input$subparagraph, included_bnfs)
    })

    # Exclude handlers for each level
    observeEvent(input$exclude_chapter, {
      add_to_list("Chapter", input$chapter, excluded_bnfs)
    })

    observeEvent(input$exclude_section, {
      add_to_list("Section", input$section, excluded_bnfs)
    })

    observeEvent(input$exclude_paragraph, {
      add_to_list("Paragraph", input$paragraph, excluded_bnfs)
    })

    observeEvent(input$exclude_subparagraph, {
      add_to_list("Subparagraph", input$subparagraph, excluded_bnfs)
    })

    # Helper function to update UI list (unchanged from original)
    update_ui_list <- function(items, list_id, item_class) {
      removeUI(
        selector = paste0("#", session$ns(list_id), " .", item_class),
        multiple = TRUE,
        immediate = TRUE
      )

      if (length(items) > 0) {
        html <- lapply(seq_along(items), function(i) {
          item <- items[[i]]
          div(
            class = item_class,
            span(paste0(item$category, ": ", item$value)),
            span(
              class = "delete-btn",
              onclick = sprintf(
                "Shiny.setInputValue('%s', %d, {priority: 'event'})",
                session$ns(paste0("delete_", sub("-item", "", item_class), "_item")),
                i
              ),
              "×"
            )
          )
        })

        insertUI(
          selector = paste0("#", session$ns(list_id)),
          where = "beforeEnd",
          ui = html,
          immediate = TRUE
        )
      }
    }

    # Update inclusion list UI
    observe({
      update_ui_list(included_bnfs(), "included-list", "included-item")
    })

    # Update exclusion list UI
    observe({
      update_ui_list(excluded_bnfs(), "excluded-list", "excluded-item")
    })

    # Handle deletion of included items
    observeEvent(input$delete_included_item, {
      index <- input$delete_included_item
      current <- included_bnfs()
      if (index <= length(current)) {
        current <- current[-index]
        included_bnfs(current)
      }
    })

    # Handle deletion of excluded items
    observeEvent(input$delete_excluded_item, {
      index <- input$delete_excluded_item
      current <- excluded_bnfs()
      if (index <= length(current)) {
        current <- current[-index]
        excluded_bnfs(current)
      }
    })

    bnf_table_server(id = "bnf_table_module", bnf_lookup = bnf_lookup)

    # Return both reactive values
    return(list(
      included = included_bnfs,
      excluded = excluded_bnfs,
      minimum_duration = reactive({ input$minimum_cp_duration }),
      polypharmacy_threshold = reactive({ as.numeric(input$polyph_threshold) }),
      earliest_treatment_end = reactive({ as.numeric(input$earliest_treatment_end) })
    ))
  })
}