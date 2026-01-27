#' mod_volcano UI Function
#'
#' @description A shiny Module.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_volcano_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      # ---- Sidebar ----
      sidebar = bslib::accordion(
        id = ns("tx_accordion"),

        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),

          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select dataset",
                shiny::icon("circle-info")
              ),
              "Select which dataset to analyse."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("de_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select DE method",
                shiny::icon("circle-info")
              ),
              "Select the differential expression method."
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
              "Select the model formula."
            ),
            choices = NULL
          ),

          shiny::checkboxInput(
            inputId = ns("highlight_truede"),
            label = "Highlight true DE genes",
            value = FALSE
          ),

          shiny::radioButtons(
            inputId = ns("pvalue"),
            label = "Adjusted p-value threshold",
            choices = c(0.1, 0.05, 0.01, 0.005, 0.001),
            selected = 0.05
          ),

          shiny::sliderInput(
            inputId = ns("logfc"),
            label = "Log2 fold-change threshold",
            min = 0,
            max = 4,
            value = 0,
            step = 0.01
          )
        )
      ),

      # ---- Main area ----
      bslib::layout_columns(
        col_widths = 12,

        div(
          style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",

          # ---- Volcano plot ----
          bslib::card(
            style = "margin-bottom: 1rem;",

            bslib::card_header(
              div(
                style = "display:flex; justify-content:space-between; align-items:center;",
                "Volcano Plot",
                shiny::downloadButton(
                  outputId = ns("download_plot"),
                  label = "",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("volcano_plot"),
                height = "550px",
                width = "100%"
              )
            )
          ),

          # ---- Table ----
          bslib::card(
            bslib::card_header("Volcano Plot Data Table"),
            bslib::card_body(
              DT::DTOutput(ns("volcano_table"))
            )
          )
        )
      )
    )
  )
}
#' mod_volcano Server Function
#'
#' @noRd
mod_volcano_server <- function(input, output, session, data) {

  # ---- Update dataset selector ----
  observe({
    shiny::updateSelectInput(
      session,
      "data",
      choices = unique(data$analysis),
      selected = unique(data$analysis)[1]
    )
  })

  observeEvent(input$data, {
    df <- data |>
      dplyr::filter(analysis == input$data)

    shiny::updateSelectInput(
      session,
      "de_method",
      choices = unique(df$method),
      selected = unique(df$method)[1]
    )
  })

  observeEvent(c(input$data, input$de_method), {
    df <- data |>
      dplyr::filter(
        analysis == input$data,
        method == input$de_method
      )

    shiny::updateSelectInput(
      session,
      "formula",
      choices = unique(df$formula),
      selected = unique(df$formula)[1]
    )
  })

  # ---- Update logFC slider dynamically ----
  observe({
    shiny::updateSliderInput(
      session,
      "logfc",
      min = 0,
      max = round(
        max(abs(data$log2FoldChange), na.rm = TRUE)
      ),
      value = 0,
      step = 0.01
    )
  })

  # ---- Volcano data ----
  volcano_data <- reactive({
    req(input$data, input$de_method, input$formula)

    data |>
      dplyr::filter(
        analysis == input$data,
        method == input$de_method,
        formula == input$formula
      ) |>
      dplyr::mutate(
        neg_log10_pvalue = -log10(pvalue),
        is_true_de = trueDE == 1,
        significance = dplyr::case_when(
          padj <= as.numeric(input$pvalue) &
            log2FoldChange > input$logfc ~ "Upregulated",
          padj <= as.numeric(input$pvalue) &
            log2FoldChange < -input$logfc ~ "Downregulated",
          TRUE ~ "Not significant"
        )
      )
  })

  # ---- Volcano plot object (single source of truth) ----
  volcano_plot_obj <- reactive({
    df <- as.data.frame(volcano_data())

    validate(
      need(nrow(df) > 0, "No data available")
    )

    max_es <- ceiling(max(abs(df$log2FoldChange), na.rm = TRUE))
    x_lims <- c(-max_es, max_es)

    if (input$highlight_truede) {
      df$significance <- ifelse(df$is_true_de, "true DE", df$significance)

      counts <- df |>
        dplyr::count(significance) |>
        tibble::deframe()

      colors <- c(
        "true DE" = "#C9980B",
        "Upregulated" = "#b2182b",
        "Downregulated" = "#2166ac",
        "Not significant" = "#c0c0c0"
      )

      labels <- paste0(names(counts), " (", counts, ")")
      names(labels) <- names(counts)

    } else {
      counts <- df |>
        dplyr::count(significance) |>
        tibble::deframe()

      colors <- c(
        "Upregulated" = "#b2182b",
        "Downregulated" = "#2166ac",
        "Not significant" = "#c0c0c0"
      )

      labels <- paste0(names(counts), " (", counts, ")")
      names(labels) <- names(counts)

    }

    # smallest adjusted p-value among non-significant genes
    padj_line <- df |>
      dplyr::filter(
       significance == "Not significant") |>
      dplyr::select(padj,gene)
    print(padj_line)
    padj_line = min(padj_line$padj, na.rm = TRUE)
    print(padj_line)


    # handle edge case: all genes significant
    if (length(padj_line) == 0 || is.na(padj_line)) {
      padj_line <- NULL
    }


    volcano_plot(
      data = df,
      x = "log2FoldChange",
      y = "padj",
      color_col = "significance",
      colors = colors,
      color_labels = labels,
      labels = "gene",
      x_axis_title = "Effect size",
      y_axis_title = "P-value Adjusted",
      plot_title = paste("Volcano Plot -", input$de_method),
      x_lims = x_lims,
      hline = padj_line,
      vline = input$logfc,
      label_size = 5
    ) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      )
  })

  # ---- Render plot ----
  output$volcano_plot <- shiny::renderPlot({
    volcano_plot_obj()
  })

  # ---- Table ----
  output$volcano_table <- DT::renderDT({
    DT::datatable(
      volcano_data(),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  # ---- Desktop/plots directory ----
  plots_dir <- file.path(
    path.expand("~"),
    "Desktop",
    "plots"
  )

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }

  # ---- Download ----
  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "volcano_",
        input$data, "_",
        input$de_method, "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {

      desktop_file <- file.path(
        plots_dir,
        paste0(
          "volcano_",
          input$data, "_",
          input$de_method, "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = volcano_plot_obj(),
        width = 7,
        height = 6,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = volcano_plot_obj(),
        width = 7,
        height = 7,
        dpi = 300
      )
    }
  )
}

