#' mod_de_heatmap UI Function
#'
#' @description A shiny Module for a DE heatmap (blue = down, red = up),
#' styled to match the volcano module.
#'
#' @param id Internal parameter for {shiny}.
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_de_heatmap_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      # ---- Sidebar ----
      sidebar = bslib::accordion(
        id = ns("hm_accordion"),

        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),

          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select dataset", shiny::icon("circle-info")),
              "Select which dataset to analyse."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("de_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select DE method", shiny::icon("circle-info")),
              "Select the differential expression method."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select formula", shiny::icon("circle-info")),
              "Select the model formula."
            ),
            choices = NULL
          ),

          # same pattern as volcano
          shiny::selectInput(
            inputId = ns("norm"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Normalization", shiny::icon("circle-info")),
              "Filter results by whether log-normalization was used (only applies when the 'norm' column exists)."
            ),
            choices = c("All" = "all"),
            selected = "all"
          ),

          shiny::radioButtons(
            inputId = ns("padj_thr"),
            label = "Adjusted p-value threshold",
            choices = c(0.1, 0.05, 0.01, 0.005, 0.001),
            selected = 0.01
          ),

          shiny::sliderInput(
            inputId = ns("min_abs_effect"),
            label = "Min |effect| threshold",
            min = 0,
            max = 4,
            value = 0,
            step = 0.01
          ),

          shiny::numericInput(
            inputId = ns("top_n"),
            label = "Max genes shown (top by padj)",
            value = 30,
            min = 5,
            max = 500,
            step = 5
          )
        )
      ),

      # ---- Main area ----
      bslib::layout_columns(
        col_widths = 12,

        div(
          style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",

          bslib::card(
            style = "margin-bottom: 1rem;",

            bslib::card_header(
              div(
                style = "display:flex; justify-content:space-between; align-items:center;",
                "DE Heatmap",
                shiny::downloadButton(
                  outputId = ns("download_heatmap"),
                  label = "",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("heatmap_plot"),
                height = "650px",
                width = "100%"
              )
            )
          ),

          bslib::card(
            bslib::card_header("Heatmap Gene Table"),
            bslib::card_body(
              DT::DTOutput(ns("heatmap_table"))
            )
          )
        )
      )
    )
  )
}


#' mod_de_heatmap Server Function
#'
#' @noRd
mod_de_heatmap_server <- function(input, output, session, data) {

  # ---- Update dataset selector ----
  observe({
    shiny::updateSelectInput(
      session,
      "data",
      choices = unique(data$analysis),
      selected = unique(data$analysis)[1]
    )
  })

  # ---- Update method selector ----
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

  # ---- Update formula + norm + effect choices ----
  observeEvent(c(input$data, input$de_method), {
    df <- data |>
      dplyr::filter(analysis == input$data, method == input$de_method)

    shiny::updateSelectInput(
      session,
      "formula",
      choices = unique(df$formula),
      selected = unique(df$formula)[1]
    )

    # norm choices (same logic as volcano)
    if ("norm" %in% names(df) && any(!is.na(df$norm))) {
      avail <- sort(unique(df$norm[!is.na(df$norm)]))
      choices <- c("All" = "all")
      if (TRUE %in% avail)  choices <- c(choices, "Normalized" = "TRUE")
      if (FALSE %in% avail) choices <- c(choices, "Not normalized" = "FALSE")

      shiny::updateSelectInput(
        session,
        "norm",
        choices = choices,
        selected = "all"
      )
    } else {
      shiny::updateSelectInput(
        session,
        "norm",
        choices = c("All" = "all"),
        selected = "all"
      )
    }

  })

  # ---- Update effect threshold slider range dynamically ----
  observeEvent(c(input$data, input$de_method, input$formula, input$norm), {
    req(input$data, input$de_method, input$formula, input$norm)

    df <- data |>
      dplyr::filter(
        analysis == input$data,
        method == input$de_method,
        formula == input$formula
      )

    if (input$norm != "all" && "norm" %in% names(df)) {
      df <- df |> dplyr::filter(norm == as.logical(input$norm))
    }


      mx <- suppressWarnings(max(abs(df$log2FoldChange), na.rm = TRUE))
      if (!is.finite(mx) || is.na(mx)) mx <- 4
      shiny::updateSliderInput(
        session,
        "min_abs_effect",
        min = 0,
        max = ceiling(mx),
        value = 0,
        step = 0.01
      )

  })

  # ---- Data used for heatmap ----
  heatmap_data <- shiny::reactive({
    req(input$data, input$de_method, input$formula, input$norm,
        input$padj_thr, input$min_abs_effect, input$top_n)

    df <- data |>
      dplyr::filter(
        analysis == input$data,
        method == input$de_method,
        formula == input$formula
      )

    if (input$norm != "all" && "norm" %in% names(df)) {
      df <- df |> dplyr::filter(norm == as.logical(input$norm))
    }



    df |>
      dplyr::filter(
        !is.na(padj),
        padj <= as.numeric(input$padj_thr),
        !is.na(.data$log2FoldChange),
        abs(.data$log2FoldChange) >= input$min_abs_effect
      ) |>
      dplyr::arrange(padj) |>
      dplyr::slice_head(n = input$top_n) |>
      dplyr::mutate(
        effect = .data$log2FoldChange,
        direction = dplyr::case_when(
          effect > 0 ~ "Upregulated",
          effect < 0 ~ "Downregulated",
          TRUE ~ "Zero"
        )
      )
  })

  # ---- Heatmap plot object (single source of truth) ----
  heatmap_plot_obj <- shiny::reactive({
    df <- heatmap_data()
    validate(shiny::need(nrow(df) > 0, "No genes pass the current thresholds."))

    # Order genes (cluster or by padj)
    gene_order <- df$gene
    mat <- matrix(df$effect, ncol = 1)
    rownames(mat) <- df$gene
    colnames(mat) <- "log2FoldChange"


    gene_order <- df$gene


    plot_df <- df |>
      dplyr::mutate(
        gene = factor(gene), # top at top visually
        metric = .data$log2FoldChange
      )

    max_abs <- max(abs(plot_df$effect), na.rm = TRUE)
    if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

    ggplot2::ggplot(plot_df, ggplot2::aes(x = .data$metric, y = .data$gene, fill = .data$metric)) +
      ggplot2::geom_tile(color = "white", linewidth = 0.3) +
      ggplot2::scale_fill_gradient2(
        low = "#2166ac",
        mid = "white",
        high = "#b2182b",
        midpoint = 0,
        limits = c(-max_abs, max_abs),
        name = "Effect"
      ) +
      ggplot2::labs(
        x = NULL,
        y = NULL,
        title = paste("DE heatmap -", input$de_method),
        subtitle = paste0(
          "padj ≤ ", input$padj_thr,
          ", |effect| ≥ ", input$min_abs_effect,
          ", top ", input$top_n, " genes"
        )
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18),
        axis.text.y = ggplot2::element_text(size = 10),
        panel.grid = ggplot2::element_blank()
      )
  })

  # ---- Render plot ----
  output$heatmap_plot <- shiny::renderPlot({
    heatmap_plot_obj()
  })

  # ---- Table ----
  output$heatmap_table <- DT::renderDT({
    DT::datatable(
      heatmap_data() |>
        dplyr::select(gene, padj, effect, direction, dplyr::any_of(c("trueDE", "true_logFC"))),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  # ---- Desktop/plots directory ----
  plots_dir <- file.path(path.expand("~"), "Desktop", "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # ---- Download ----
  output$download_heatmap <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "heatmap_",
        input$data, "_",
        gsub("[^A-Za-z0-9_]+", "_", input$de_method), "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {

      desktop_file <- file.path(
        plots_dir,
        paste0(
          "heatmap_",
          input$data, "_",
          gsub("[^A-Za-z0-9_]+", "_", input$de_method), "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = heatmap_plot_obj(),
        width = 7,
        height = 7,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = heatmap_plot_obj(),
        width = 7,
        height = 7,
        dpi = 300
      )
    }
  )
}
