#' mod_roc UI Function
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
mod_roc_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      # ---- Sidebar ----
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
              trigger = shiny::tagList(
                "Select formula",
                shiny::icon("circle-info")
              ),
              "Select which DE model formula to evaluate."
            ),
            choices = NULL
          )
        )
      ),

      # ---- Main panel ----
      bslib::layout_columns(
        col_widths = 12,

        bslib::navset_card_tab(

          # ===== ROC =====
          bslib::nav_panel(
            title = "ROC Curve",
            icon = shiny::icon("chart-line"),

            bslib::card(
              style = "margin-bottom: 1rem;",
              bslib::card_header(
                div(
                  style = "display:flex; justify-content:space-between; align-items:center;",
                  "ROC Curve",
                  shiny::downloadButton(ns("download_roc"), "Download")
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
          ),

          # ===== PR =====
          bslib::nav_panel(
            title = "PR Curve",
            icon = shiny::icon("chart-area"),

            bslib::card(
              style = "margin-bottom: 1rem;",
              bslib::card_header(
                div(
                  style = "display:flex; justify-content:space-between; align-items:center;",
                  "PR Curve",
                  shiny::downloadButton(ns("download_pr"), "Download")
                )
              ),
              bslib::card_body(
                shiny::plotOutput(ns("pr_plot"), height = "550px")
              )
            ),

            bslib::card(
              bslib::card_header("PR Curve Table"),
              bslib::card_body(
                DT::DTOutput(ns("pr_table"))
              )
            )
          )
        )
      )
    )
  )
}

#' mod_roc Server Function
#'
#' @noRd
mod_roc_server <- function(input, output, session, data) {

  # ---- Update selectors ----
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
      session,
      "method",
      choices = unique(df$method),
      selected = unique(df$method)[1]
    )
  })

  observeEvent(c(input$analysis, input$method), {
    req(input$analysis, input$method)

    df <- data |>
      dplyr::filter(
        analysis == input$analysis,
        method == input$method
      )

    shiny::updateSelectInput(
      session,
      "formula",
      choices = unique(df$formula),
      selected = unique(df$formula)[1]
    )
  })

  # ---- Filtered ROC/PR data ----
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

  # ---- evalmod object ----
  rocpr_obj <- reactive({
    df <- rocpr_data()

    validate(
      need(nrow(df) > 0, "Loading...")
    )

    precrec::evalmod(
      scores = df$score,
      labels = df$trueDE
    )
  })

  # ---- FDR = 0.05 operating point ----
  roc_fdr05_point <- reactive({
    df <- rocpr_data()

    df <- df |>
      dplyr::arrange(padj) |>
      dplyr::mutate(
        discoveries = dplyr::row_number(),
        false_pos = cumsum(trueDE == 0),
        true_pos  = cumsum(trueDE == 1),
        fdr = false_pos / discoveries
      )

    df_fdr <- df |> dplyr::filter(fdr <= 0.01)

    validate(
      need(nrow(df_fdr) > 0, "No FDR ≤ 0.01 point available")
    )

    pt <- df_fdr |> dplyr::slice_tail(n = 1)

    n_pos <- sum(df$trueDE == 1)
    n_neg <- sum(df$trueDE == 0)

    tibble::tibble(
      FPR = pt$false_pos / n_neg,
      TPR = pt$true_pos  / n_pos,
      padj = pt$padj
    )
  })

  # ---- ROC plot ----
  roc_plot_obj <- reactive({
    req(rocpr_obj())

    p <- suppressWarnings(
      ggplot2::autoplot(rocpr_obj(), curvetype = "ROC")
    ) +
      ggplot2::labs(
        title = paste("ROC Curve -", input$method),
        subtitle = paste("Analysis:", input$analysis)
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(aspect.ratio = 1)

    if (!is.null(roc_fdr05_point())) {
      pt <- roc_fdr05_point()

      p <- p +
        ggplot2::geom_point(
          data = pt,
          ggplot2::aes(x = FPR, y = TPR),
          color = "#C9980B",
          inherit.aes = FALSE,
          size = 3
        ) +
        ggplot2::geom_text(
          data = pt,
          ggplot2::aes(
            x = FPR,
            y = TPR,
            label = "FDR = 0.01"
          ),
          hjust = -0.1,
          vjust = -0.5,
          inherit.aes = FALSE,
          size = 5,
          color = "#C9980B"
        )
    }

    p
  })

  # ---- PR plot ----
  pr_plot_obj <- reactive({
    req(rocpr_obj())

    suppressWarnings(
      ggplot2::autoplot(rocpr_obj(), curvetype = "PR")
    ) +
      ggplot2::labs(
        title = paste("PR Curve -", input$method),
        subtitle = paste("Analysis:", input$analysis)
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(aspect.ratio = 1)
  })

  # ---- Render plots ----
  output$roc_plot <- shiny::renderPlot({
    roc_plot_obj()
  })

  output$pr_plot <- shiny::renderPlot({
    pr_plot_obj()
  })

  # ---- Tables ----
  output$roc_table <- DT::renderDT({
    DT::datatable(
      precrec::auc(rocpr_obj()),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })

  output$pr_table <- DT::renderDT({
    DT::datatable(
      precrec::auc(rocpr_obj()),
      options = list(
        pageLength = 10,
        scrollX = TRUE
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

  # ---- Download ROC ----
  output$download_roc <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "ROC_", input$analysis, "_", input$method, "_", Sys.Date(), ".png"
      )
    },
    content = function(file) {

      desktop_file <- file.path(
        plots_dir,
        paste0(
          "ROC_", input$analysis, "_", input$method, "_", Sys.Date(), ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = roc_plot_obj(),
        width = 7,
        height = 6,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = roc_plot_obj(),
        width = 7,
        height = 6,
        dpi = 300
      )
    }
  )

  # ---- Download PR ----
  output$download_pr <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "PR_", input$analysis, "_", input$method, "_", Sys.Date(), ".png"
      )
    },
    content = function(file) {

      desktop_file <- file.path(
        plots_dir,
        paste0(
          "PR_", input$analysis, "_", input$method, "_", Sys.Date(), ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = pr_plot_obj(),
        width = 7,
        height = 7,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = pr_plot_obj(),
        width = 7,
        height = 7,
        dpi = 300
      )
    }
  )
}
