#' mod_goi UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_goi_ui <- function(id) {
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
            inputId = ns("gene"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select gene for analysis:",
                shiny::icon("circle-info")
              ),
              "Select which is the gene that you want to use."
            ),
            choices = NULL
          )
        )
      ),

      # MAIN AREA: scrollable container with plot (fixed height) + table below
      bslib::layout_columns(
        col_widths = 12,
        div(
          style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",
          # goi plot card
          bslib::card(
            style = "margin-bottom: 1rem;",
            bslib::card_header(
              div(
                style = "display: flex; justify-content: space-between;
                align-items: center;",
                "Gene  of Interest Plot",
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
                shiny::plotOutput(ns("goi_plot"),
                  height = "550px",
                  width = "100%"
                )
              )
            )
          )
        )
      )
    )
  )
}


#' mod_goi Server Functions
#'
#' @noRd
mod_goi_server <- function(input, output, session, data) {
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

  observe({
    df <- data |>
      dplyr::filter(method == input$de_method) |>
      dplyr::filter(analysis == input$data)
    shiny::updateSelectInput(
      session,
      "gene",
      choices = sort(unique(df$gene)),
      selected = unique(df$gene)[1]
    )
  })


  goi_data <- reactive({
    req(input$de_method, input$data)
    df <- data |>
      dplyr::filter(method == input$de_method) |>
      dplyr::filter(analysis == input$data) |>
      dplyr::filter(gene == input$gene)
    return(df)
  })

  output$goi_plot <- shiny::renderPlot({
    req(goi_data())
    if (nrow(goi_data()) != 0) {
      df <- as.data.frame(goi_data())
      ggplot2::ggplot(df, ggplot2::aes(
        x = analysis,
        y = log2FoldChange,
        group = gene,
        color = gene
      )) +
        ggplot2::geom_line(size = 1.1) +
        ggplot2::geom_point(size = 2) +
        ggplot2::theme_minimal(base_size = 16) +
        ggplot2::labs(
          title = "Effect Size Change by Gene",
          x = "Condition / Comparison",
          y = "Log2 Fold Change",
          color = "Gene"
        ) +
        ggplot2::theme(
          axis.text.x = ggplot2::element_text(angle = 45, hjust = 1),
          legend.position = "right",
          aspect.ratio = 1,
          text = ggplot2::element_text(size = 18)
        )
    }
  })

  #
  #   output$download_plot <- shiny::downloadHandler(
  #     filename = function() {
  #       paste("goi_plot", ".png", sep = "")
  #     },
  #     content = function(file) {
  #       df <- as.data.frame(goi_data())
  #
  #       color_order <- c(
  #         "Upregulated" = "#b2182b", "Downregulated" = "#2166ac",
  #         "Not significant" = "#c0c0c0"
  #       )
  #
  #       goi_plot(
  #         data = df, x = "log2FoldChange", y = "pvalue",
  #         color_col = "significance", colors = color_order,
  #         labels = "gene", y_axis_title = "P-value",
  #         x_axis_title = "Effect Size", hline = -log10(input$pvalue),
  #         vline = input$logfc
  #       ) + ggplot2::theme(
  #         aspect.ratio = 1,
  #         text = ggplot2::element_text(size = 18)
  #       )
  #     }
  #   )
}
