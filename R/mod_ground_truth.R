mod_gt_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",
      sidebar = bslib::accordion(
        id = ns("gt_accordion"),
        bslib::accordion_panel(
          title = "GT settings",
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
          ),
          shiny::radioButtons(
            inputId = ns("conf_level"),
            label = "Confidence level of the intervals",
            choices = c(
              "80%" = 0.80,
              "90%" = 0.90,
              "95%" = 0.95,
              "97.5%" = 0.975,
              "99%" = 0.99
            ),
            selected = 0.95
          ),
          shiny::selectInput(
            inputId = ns("gene"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select genes", shiny::icon("circle-info")),
              "Select which genes you want to display."
            ),
            choices = NULL,
            multiple = TRUE
          ),
        )
      ),
      bslib::layout_columns(
        col_widths = 12,
        bslib::nav_panel(
          title = "Ground Truth",
          icon = shiny::icon("chart-line"),
          div(
            style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",
            # ROC plot card
            bslib::card(
              style = "margin-bottom: 1rem;",
              bslib::card_header(
                div(
                  style = "display:flex; justify-content:space-between; align-items:center;",
                  "Ground Truth",
                  shiny::downloadButton(ns("download_gt"), "")
                )
              ),
              bslib::card_body(
                shiny::plotOutput(ns("gt_plot"), height = "550px")
              )
            ),

            # ROC table card
            bslib::card(
              bslib::card_header("Grund Truth Table"),
              bslib::card_body(
                div(
                  style = "width: 100%",
                  DT::DTOutput(ns("gt_table"))
                )
              )
            )
          )
        )
      )
    )
  )
}


#' mod_gt Server Functions
#'
#' @noRd
mod_gt_server <- function(input, output, session, data) {
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

  observeEvent(c(input$analysis, input$method), {
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

  observeEvent(c(input$analysis, input$method, input$formula), {
    req(input$analysis, input$method, input$formula)
    df <- data |> dplyr::filter(
      analysis == input$analysis,
      method == input$method,
      formula == input$formula
    )
    selected <- df |>
      dplyr::mutate(abs_log2FC = abs(log2FoldChange)) |>
      dplyr::arrange(dplyr::desc(abs_log2FC)) |>
      dplyr::slice_head(n = 20)

    shiny::updateSelectInput(
      session, "gene",
      choices = unique(df$gene),
      selected = unique(selected$gene)
    )
  })


  gt_data <- reactive({
    req(input$analysis, input$method, input$formula)

    data |>
      dplyr::filter(
        analysis == input$analysis,
        method == input$method,
        formula == input$formula,
        gene %in% input$gene
      ) |>
      dplyr::mutate(score = -log10(padj)) |>
      dplyr::mutate(
        lower = log2FoldChange - qnorm(1 -
          (1 - as.numeric(input$conf_level)) / 2) * lfcSE,
        upper = log2FoldChange + qnorm(1 -
          (1 - as.numeric(input$conf_level)) / 2) * lfcSE
      )
  })

  output$gt_plot <- shiny::renderPlot({
    req(gt_data())
    df <- gt_data() |>
      dplyr::arrange(log2FoldChange) |>
      dplyr::mutate(gene = factor(gene, levels = gene))


    ggplot2::ggplot(
      df,
      ggplot2::aes(
        y = gene
      )
    ) +
      # error bars for CIs
      ggplot2::geom_errorbar(
        ggplot2::aes(xmin = lower, xmax = upper),
        width = 0.2,
        linewidth = 1.2,
        color = "steelblue"
      ) +
      ggplot2::geom_point(
        ggplot2::aes(x = log2FoldChange, color = trueDE),
        size = 3
      ) +
      ggplot2::scale_color_manual(values = c(
        `TRUE` = "#C9980B",
        `FALSE` = "steelblue")
        ) +
      ggplot2::geom_point(
        ggplot2::aes(x = true_logFC),
        color = "#C9980B", size = 2
      ) +
      ggplot2::labs(
        title = paste("", input$gene, "in", input$analysis),
        x     = "Log2 Fold Change",
        y     = "Genes"
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        axis.text.x = ggplot2::element_text(angle = 30, hjust = 1),
        legend.position = "right",
        aspect.ratio = 1
      )
  })

  # ----------------------
  # GT TABLE
  # ----------------------
  output$gt_table <- DT::renderDT({
    req(gt_data())

    dt <- gt_data() |>
      dplyr::select(
        gene,
        log2FoldChange,
        lfcSE,
        lower,
        upper,
        true_logFC,
        trueDE,
        pvalue,
        padj
      ) |>
      dplyr::arrange(log2FoldChange)

    DT::datatable(
      dt,
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

  output$download_gt <- shiny::downloadHandler(
    filename = function() {
      # Suggest the full path to Windows
      file.path("C:/Users/nerea/OneDrive/Escritorio/imagenes",
                paste0("ground_truth.png"))
    },
    content = function(file) {
      req(gt_data())
      df <- gt_data() |>
        dplyr::arrange(log2FoldChange) |>
        dplyr::mutate(gene = factor(gene, levels = gene))

      ggplot2::ggsave(
        filename = file,
        plot = ggplot2::ggplot(
          df,
          ggplot2::aes(y = gene)
        ) +
          ggplot2::geom_errorbar(
            ggplot2::aes(xmin = lower, xmax = upper),
            width = 0.2,
            linewidth = 1.2,
            color = "#C9980B"
          ) +
          ggplot2::geom_point(
            ggplot2::aes(x = log2FoldChange),
            color = "steelblue", size = 3
          ) +
          ggplot2::geom_point(
            ggplot2::aes(x = true_logFC),
            color = "", size = 2
          ) +
          ggplot2::labs(
            title = paste("", input$gene, "in", input$analysis),
            x     = "Log2 Fold Change",
            y     = "Genes"
          ) +
          ggplot2::theme_bw(base_size = 16),
        width = 8,
        height = 8,
        dpi = 300
      )
    }
  )
}
