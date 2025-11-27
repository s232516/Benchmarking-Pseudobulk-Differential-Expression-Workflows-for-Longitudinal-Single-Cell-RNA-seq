mod_roc_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      sidebar = bslib::accordion(
        id = ns("roc_accordion"),

        bslib::accordion_panel(
          title = "ROC / PR settings",
          icon = shiny::icon("sliders"),

          shiny::selectInput(
            ns("analysis"),
            "Select simulation dataset:",
            choices = NULL
          ),

          shiny::selectInput(
            ns("method"),
            "Select DE method:",
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select formula", shiny::icon("circle-info")),
              "Select which DE model formula to evaluate."
            ),
            choices = NULL
          )
        )
      ),


      bslib::layout_columns(
        col_widths = 12,
        bslib::navset_card_tab(


          bslib::nav_panel(
            title = "ROC Curve",
            icon = shiny::icon("chart-line"),

            div(
              style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",

              # ROC plot card
              bslib::card(
                style = "margin-bottom: 1rem;",
                bslib::card_header(
                  div(
                    style="display:flex; justify-content:space-between; align-items:center;",
                    "ROC Curve",
                    shiny::downloadButton(ns("download_roc"), "")
                  )
                ),
                bslib::card_body(
                  shiny::plotOutput(ns("roc_plot"), height = "550px")
                )
              ),

              # ROC table card
              bslib::card(
                bslib::card_header("ROC AUC Table"),
                bslib::card_body(
                  div(
                    style = "width: 100%",
                    DT::DTOutput(ns("roc_table"))
                  )
                )
              )
            )
          ),


          bslib::nav_panel(
            title = "PR Curve",
            icon = shiny::icon("chart-area"),

            div(
              style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",

              # PR plot card
              bslib::card(
                style = "margin-bottom: 1rem;",
                bslib::card_header(
                  div(
                    style="display:flex; justify-content:space-between; align-items:center;",
                    "PR Curve",
                    shiny::downloadButton(ns("download_pr"), "")
                  )
                ),
                bslib::card_body(
                  shiny::plotOutput(ns("pr_plot"), height = "550px")
                )
              ),

              # PR table card
              bslib::card(
                bslib::card_header("PR Curve Table"),
                bslib::card_body(
                  div(
                    style = "width: 100%",
                    DT::DTOutput(ns("pr_table"))
                  )
                )
              )
            )
          )
        )
      )
    )
  )
}


#' mod_roc Server Functions
#'
#' @noRd
mod_roc_server <- function(input, output, session, data) {
  ns <- session$ns


  observe({
    shiny::updateSelectInput(
      session,
      "analysis",
      choices = unique(data$analysis),
      selected = unique(data$analysis)[1]
    )
  })

  observeEvent(input$analysis, {
    req(input$analysis)
    df <- data |> dplyr::filter(analysis == input$analysis)

    shiny::updateSelectInput(
      session, "method",
      choices = unique(df$method),
      selected = unique(df$method)[1]
    )
  })

  observeEvent(c(input$analysis,input$method), {
    req(input$analysis, input$method)
    df <- data |> dplyr::filter(
      analysis == input$analysis,
      method == input$method
    )
    shiny::updateSelectInput(
      session, "formula",
      choices = unique(df$formula),
      selected = unique(df$formula)[1]
    )
  })


  rocpr_data <- reactive({
    req(input$analysis, input$method, input$formula)

    data |>
      dplyr::filter(
        analysis == input$analysis,
        method == input$method,
        formula == input$formula,
        !is.na(trueDE)
      ) |>
      dplyr::mutate(score = -log10(padj))
  })

  # ----------------------
  # SHARED EVALMOD OBJECT
  # ----------------------
  rocpr_obj <- reactive({
    df <- rocpr_data()

    # prevent evalmod from running on empty data
    validate(
      need(!is.null(df) && nrow(df) > 0,
           "Loading...")
    )

    precrec::evalmod(
      scores = df$score,
      labels = df$trueDE
    )
  })

  # ----------------------
  # ROC PLOT
  # ----------------------
  output$roc_plot <- shiny::renderPlot({
    req(rocpr_obj())

    plt <- suppressWarnings(ggplot2::autoplot(rocpr_obj(), curvetype = "ROC")) +
      ggplot2::labs(
        title = paste("ROC Curve -", input$method),
        subtitle = paste("Analysis:", input$analysis)
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(aspect.ratio = 1)

    plt
  })

  # ----------------------
  # ROC TABLE
  # ----------------------
  output$roc_table <- DT::renderDT({
    req(rocpr_obj())

    auc <- precrec::auc(rocpr_obj())
    DT::datatable(
      auc,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  # ----------------------
  # PR PLOT
  # ----------------------
  output$pr_plot <- shiny::renderPlot({
    req(rocpr_obj())

    plt <- suppressWarnings(ggplot2::autoplot(rocpr_obj(), curvetype = "PR")) +
      ggplot2::labs(
        title = paste("PR Curve -", input$method),
        subtitle = paste("Analysis:", input$analysis)
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(aspect.ratio = 1)

    plt
  })

  # ----------------------
  # PR TABLE  (same style as volcano)
  # ----------------------
  output$pr_table <- DT::renderDT({
    req(rocpr_obj())
    auc <- precrec::auc(rocpr_obj())
    DT::datatable(
      auc,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  # ----------------------
  # DOWNLOAD HANDLERS
  # ----------------------
  output$download_roc <- shiny::downloadHandler(
    filename = function() {
      paste0("ROC_", input$analysis, "_", input$method, ".png")
    },
    content = function(file) {
      ggsave(file, output$roc_plot(), width = 7, height = 6)
    }
  )

  output$download_pr <- shiny::downloadHandler(
    filename = function() {
      paste0("PR_", input$analysis, "_", input$method, ".png")
    },
    content = function(file) {
      ggsave(file, output$pr_plot(), width = 7, height = 6)
    }
  )
}
