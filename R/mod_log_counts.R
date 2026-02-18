mod_logfc_count_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      sidebar = bslib::accordion(
        id = ns("tx_accordion"),
        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),
          id = ns("tx_analysis"),

          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select dataset", shiny::icon("circle-info")),
              "Select which dataset you want to use."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("de_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select method", shiny::icon("circle-info")),
              "Select which differential expression method you want to use."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select formula", shiny::icon("circle-info")),
              "Select which design formula you want to use."
            ),
            choices = NULL
          )
        )
      ),

      bslib::layout_columns(
        col_widths = 12,
        div(
          style = "overflow-y:auto; padding-right:10px; padding-left:10px;",

          bslib::card(
            style = "margin-bottom: 1rem;",
            bslib::card_header(
              div(
                style = "display:flex; justify-content:space-between; align-items:center;",
                "Counts vs log2 Fold Change (True DE)",
                shiny::downloadButton(
                  outputId = ns("download_plot"),
                  label    = "",
                  class    = "clickable-button-light"
                )
              )
            ),
            bslib::card_body(
              style = "display:flex; justify-content:space-between; align-items:center;",
              shiny::plotOutput(
                outputId = ns("hist_logfc"),
                height   = "550px",
                width    = "100%"
              )
            )
          )
        )
      )
    )
  )
}

mod_logfc_count_server <- function(input, output, session, data, sim) {

  observe({
    shiny::updateSelectInput(
      session, "data",
      choices  = unique(data$analysis),
      selected = unique(data$analysis)[1]
    )
  })

  observeEvent(input$data, {
    df <- data |>
      dplyr::filter(analysis == input$data)

    choices <- unique(df$method)

    shiny::updateSelectInput(
      session, "de_method",
      choices  = choices,
      selected = choices[1]
    )
  }, ignoreInit = TRUE)

  observeEvent(c(input$data, input$de_method), {
    df <- data |>
      dplyr::filter(
        analysis == input$data,
        method   == input$de_method
      )

    choices <- unique(df$formula)

    shiny::updateSelectInput(
      session, "formula",
      choices  = choices,
      selected = choices[1]
    )
  }, ignoreInit = TRUE)

  df_de <- reactive({
    req(input$data, input$de_method, input$formula)
    data |>
      dplyr::filter(
        analysis == input$data,
        method   == input$de_method,
        formula  == input$formula,
        trueDE
      )
  })

  df_counts <- reactive({
    req(input$data)

    sce <- sim[[input$data]]
    req(!is.null(sce))

    counts <- SummarizedExperiment::assay(sce, "counts")
    counts <- counts[, , drop = FALSE]

    dplyr::tibble(
      gene   = rownames(counts),
      counts = Matrix::rowSums(counts)
    ) |>
      dplyr::mutate(
        gene = dplyr::if_else(
          grepl("^Gene", .data$gene),
          paste0("Gene_", readr::parse_number(.data$gene)),
          .data$gene
        )
      )
  })

  final_df <- reactive({
    req(df_de(), df_counts())
    dplyr::left_join(df_de(), df_counts(), by = "gene")
  })

  plot_obj <- reactive({
    req(final_df())
    df <- final_df()

    ggplot2::ggplot(df, ggplot2::aes(x = counts, y = log2FoldChange)) +
      ggplot2::geom_point(size = 3, alpha = 0.8) +
      ggplot2::labs(
        title = "Counts vs log2FC (True DE)",
        x = "Total counts per gene",
        y = "Estimated log2 fold-change"
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      )
  })

  output$hist_logfc <- shiny::renderPlot({
    plot_obj()
  })

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      req(input$data, input$de_method, input$formula)
      paste0(
        "counts_vs_log2FC__",
        input$data, "__",
        gsub("\\s+", "_", input$de_method), "__",
        gsub("[^A-Za-z0-9]+", "_", input$formula),
        ".png"
      )
    },
    contentType = "image/png",
    content = function(file) {
      ggplot2::ggsave(file, plot = plot_obj(), width = 8, height = 6, dpi = 300)
    }
  )
}
