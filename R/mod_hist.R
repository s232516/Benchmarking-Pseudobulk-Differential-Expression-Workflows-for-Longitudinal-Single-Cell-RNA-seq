#' mod_volcano UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_hist_logfc_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      # Sidebar controls
      sidebar = bslib::accordion(
        id = ns("tx_accordion"),
        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),
          id = ns("tx_analysis"),
          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select dataset",
                                       shiny::icon("circle-info")),
              "Select which is the dataset that you want to use for
                dimensionality reduction."
            ),
            choices = NULL
          ),
          shiny::selectInput(
            inputId = ns("de_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select dataset",
                shiny::icon("circle-info")
              ),
              "Select which is the dataset that you want to use for
              dimensionality reduction."
            ),
            choices = NULL
          ),
          shiny::selectInput(
            inputId = ns("formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select formula",
                shiny::icon("circle-info")
              ),
              "Select which is the formula that you want to use for
              dimensionality reduction."
            ),
            choices = NULL
        ),

        shiny::sliderInput(
          ns("bins"),
          "Histogram bins:",
          min = 5, max = 80, value = 40
        )
      )
      ),

      # MAIN AREA: scrollable container with plot (fixed height) + table below
      bslib::layout_columns(
        col_widths = 12,
        div(
          style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",
          # volcano plot card
          bslib::card(
            style = "margin-bottom: 1rem;",
            bslib::card_header(
              div(
                style = "display: flex; justify-content: space-between;
                align-items: center;",
                "Histogram of log2 Fold Changes",
                shiny::downloadButton(
                  class = "clickable-button-light",
                  outputId = ns("download_plot"),
                  label = ""
                )
              )
            ),
            bslib::card_body(
              div(
                style = "width: 100%; height: 550px; overflow: hidden;",
                shiny::plotOutput(ns("hist_logfc"), height = "550px",
                                  width = "100%")
              )
            )
          )
        )
      )
    )
  )
}

#' mod_volcano Server Functions
#'
#' @noRd
mod_hist_logfc_server <- function(input, output, session, data) {
  ns <- session$ns

  observeEvent(input$data, {
    df <- data |>
      dplyr::filter(analysis == input$data)
    choices <- unique(df$method)
    shiny::updateSelectInput(
      session,
      "de_method",
      choices = choices,
      selected = choices[1]
    )
  })

  observe({
    shiny::updateSelectInput(
      session,
      "data",
      choices = unique(data$analysis),
      selected = unique(data$analysis)[1]
    )
  })

  observeEvent(c(input$data,input$de_method), {
    df <- data |>
      dplyr::filter(analysis == input$data,
                    method == input$de_method)
    choices <- unique(df$formula)
    shiny::updateSelectInput(
      session,
      "formula",
      choices = choices,
      selected = choices[1]
    )
  })

  df <- reactive({
    req( input$de_method, input$data)
    df <- data |>
      dplyr::filter(method == input$de_method) |>
      dplyr::filter(analysis == input$data) |>
      dplyr::filter(formula == input$formula)
  })

    output$hist_logfc <- shiny::renderPlot({
      req(df())

      df <- df()

      ggplot2::ggplot(df, ggplot2::aes(x = log2FoldChange, fill = trueDE)) +
        ggplot2::geom_histogram(color = "black", alpha = 0.7, bins = input$bins,
                                position = "identity") +
        ggplot2::scale_fill_manual(
          values = c(`TRUE` = "red", `FALSE` = "steelblue")
        ) +
        ggplot2::labs(
          title = "Histogram of log2 Fold Changes",
          x = "Log2 Fold Change",
          y = "Count",
          fill = "True DE"
        )  +
        ggplot2::scale_x_continuous(
          breaks = scales::breaks_pretty(n = 10)
        ) +
        ggplot2::theme_bw(base_size = 16) +
        ggplot2::theme(
          aspect.ratio = 1,
          axis.title.x = ggplot2::element_text(size = 30),
          axis.title.y = ggplot2::element_text(size = 30),
          axis.text.x  = ggplot2::element_text(size = 20),
          axis.text.y  = ggplot2::element_text(size = 20)
        )
    })

    file_name <- reactive({
      req(input$data, input$de_method, input$formula)
      paste0(
        "hist_log2FC__",
        input$data, "__",
        gsub("\\s+", "_", input$de_method),
        ".png"
      )
    })



    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        file_name()
      },
      contentType = "image/png",
      content = function(file) {
        req(df(), input$bins)

        df_local <- df()

        plot1 <- ggplot2::ggplot(df_local, ggplot2::aes(x = log2FoldChange, fill = trueDE)) +
          ggplot2::geom_histogram(
            color = "black", alpha = 0.7, bins = input$bins,
            position = "identity"
          ) +
          ggplot2::scale_fill_manual(values = c(`TRUE` = "red", `FALSE` = "steelblue")) +
          ggplot2::labs(
            title = "Histogram of log2 Fold Changes",
            x = "Log2 Fold Change",
            y = "Count",
            fill = "True DE"
          ) +
          ggplot2::scale_x_continuous(breaks = scales::breaks_pretty(n = 6)) +
          ggplot2::theme_bw(base_size = 16) +
          ggplot2::theme(
            aspect.ratio = 1,
            text = ggplot2::element_text(size = 18),
            axis.title.x = ggplot2::element_text(size = 30),
            axis.title.y = ggplot2::element_text(size = 30),
            axis.text.x  = ggplot2::element_text(size = 20),
            axis.text.y  = ggplot2::element_text(size = 20)
          )

        ggplot2::ggsave(file, plot = plot1, width = 8, height = 6, dpi = 300)
      }
    )

}

