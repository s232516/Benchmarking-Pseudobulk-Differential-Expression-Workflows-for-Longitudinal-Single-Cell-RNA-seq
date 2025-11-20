#' mod_roc UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
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
          title = "ROC settings",
          icon = shiny::icon("chart-line"),

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

          shiny::checkboxInput(
            ns("smooth"),
            "Smooth ROC curve",
            value = FALSE
          )
        )
      ),

      bslib::layout_columns(
        col_widths = 12,

        bslib::card(
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

        bslib::card(
          bslib::card_header("ROC AUC Table"),
          bslib::card_body(
            DT::DTOutput(ns("roc_table"))
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

    # populate analysis dropdown
    observe({
      shiny::updateSelectInput(
        session,
        "analysis",
        choices = unique(data$analysis),
        selected = unique(data$analysis)[1]
      )
    })

    # populate method dropdown based on analysis
    observeEvent(input$analysis, {
      df <- data |>
        dplyr::filter(analysis == input$analysis)

      shiny::updateSelectInput(
        session,
        "method",
        choices = unique(df$method),
        selected = unique(df$method)[1]
      )
    })

    # reactive data for the selected method & analysis
    roc_data <- reactive({
      req(input$analysis, input$method)

      df <- data |>
        dplyr::filter(
          analysis == input$analysis,
          method == input$method
        ) |>
        dplyr::filter(!is.na(trueDE)) |>
        dplyr::mutate(score = -log10(.data$padj))

      df
    })

    # compute ROC using precrec
    roc_obj <- reactive({
      df <- roc_data()
      precrec::evalmod(
        scores = df$score,
        labels = df$trueDE
      )
    })

    # plot ROC
    output$roc_plot <- shiny::renderPlot({
      req(roc_obj())
      plt <- suppressWarnings(
        ggplot2:::autoplot(roc_obj(), curvetype = "ROC")
      ) +
        ggplot2::labs(
          title = paste("ROC Curve -", input$method),
          subtitle = paste("Analysis:", input$analysis)
        ) +
        ggplot2::theme_bw(base_size = 16)  +
        ggplot2::theme(
          aspect.ratio = 1,
          text = ggplot2::element_text(size = 18)
        )


      if (input$smooth) {
        plt <- plt + ggplot2::geom_smooth(se = FALSE)
      }

      plt
    })

    # AUC table
    output$roc_table <- DT::renderDT({
      auc <- precrec::auc(roc_obj())
      DT::datatable(auc)
    })

    # download handler
    output$download_roc <- shiny::downloadHandler(
      filename = function() {
        paste0("ROC_", input$analysis, "_", input$method, ".png")
      },
      content = function(file) {
        ggsave(file, plot = output$roc_plot(), width = 7, height = 6)
      }
    )
}
