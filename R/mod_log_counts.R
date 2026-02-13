#' mod_volcano UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_logfc_count_ui <- function(id) {
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
              style = "display = flex; justify-content: space-between;
            align-items: center;",
              shiny::plotOutput(
                outputId = ns("hist_logfc"),
                height = "550px",
                width = "800px"
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
mod_logfc_count_server <- function(input, output, session, data, sim) {
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

  observeEvent(c(input$data, input$de_method), {
    df <- data |>
      dplyr::filter(
        analysis == input$data,
        method == input$de_method
      )
    choices <- unique(df$formula)
    shiny::updateSelectInput(
      session,
      "formula",
      choices = choices,
      selected = choices[1]
    )
  })

  df <- reactive({
    req(input$data, input$de_method)
    df <- data |>
      dplyr::filter(method == input$de_method) |>
      dplyr::filter(analysis == input$data) |>
      dplyr::filter(formula == input$formula) |>
      dplyr::filter(trueDE)
    df
  })

  simu <- reactive({
    req(input$data)

    counts <- sim[[input$data]]@assays@data$counts

    # FORCE 2D matrix
    counts <- counts[, , drop = FALSE]

    df <- data.frame(
      gene = ifelse(
        grepl("^Gene", rownames(counts)),
        paste0("Gene_", readr::parse_number(rownames(counts))),
        rownames(counts)
      ),
      counts = Matrix::rowSums(counts),
      row.names = NULL
    )
    print(df)
    df
  })


  final_df <- reactive({
    req(df(), simu())
    df <- df()
    sim <- simu()
    merged <- df |>
      dplyr::left_join(sim)
    merged
  })

  output$hist_logfc <- shiny::renderPlot({
    req(final_df())
    df <- final_df()

    ggplot2::ggplot(df, ggplot2::aes(x = counts, y = log2FoldChange)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        title = "Histogram of log2 Fold Changes (True DE)",
        x = "Count",
        y = "Log2 Fold Change"
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      )
  })
}
