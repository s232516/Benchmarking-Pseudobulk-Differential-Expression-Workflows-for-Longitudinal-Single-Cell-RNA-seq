#' mod_counts_hist UI Function
#'
#' @description Histogram of raw and log-normalized counts.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
mod_counts_hist_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      # ---- Sidebar ----
      sidebar = bslib::accordion(
        id = ns("counts_accordion"),

        bslib::accordion_panel(
          title = "Counts histogram settings",
          icon = shiny::icon("sliders"),

          shiny::selectInput(
            ns("dataset"),
            "Select simulation:",
            choices = NULL
          ),

          shiny::radioButtons(
            ns("scale"),
            "Histogram scale:",
            choices = c(
              "Raw counts" = "raw",
              "Log-normalized counts" = "lognorm"
            ),
            selected = "raw"
          ),

          shiny::sliderInput(
            ns("bins"),
            "Number of bins:",
            min = 10, max = 200, value = 50
          )
        )
      ),

      # ---- Main area ----
      bslib::layout_columns(
        col_widths = 12,

        bslib::nav_panel(
          title = "Counts distribution",
          icon = shiny::icon("chart-bar"),

          bslib::card(
            bslib::card_header(
              div(
                style = "display:flex; justify-content:space-between; align-items:center;",
                "Counts histogram",
                shiny::downloadButton(
                  outputId = ns("download_counts_hist"),
                  label = "",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                ns("counts_hist_plot"),
                height = "550px"
              )
            )
          )
        )
      )
    )
  )
}
#' mod_counts_hist Server Function
#'
#' @noRd
mod_counts_hist_server <- function(input, output, session, data) {

  # ---- Populate dataset selector ----
  observe({
    shiny::updateSelectInput(
      session,
      "dataset",
      choices = names(data),
      selected = names(data)[1]
    )
  })

  # ---- Extract counts matrix safely ----
  counts_mat <- reactive({
    req(input$dataset)
    SummarizedExperiment::assay(data[[input$dataset]], "counts")
  })

  # ---- Raw counts (flattened) ----
  raw_counts <- reactive({
    req(counts_mat())
    as.vector(counts_mat())
  })

  # ---- Log-normalized counts (exact formula) ----
  lognorm_counts <- reactive({
    req(counts_mat())

    norm_mat <- apply(counts_mat(), 2, function(x) {
      log2((1e6 * x / sum(x)) + 1)
    })

    as.vector(norm_mat)
  })

  # ---- Choose values to plot ----
  plot_values <- reactive({
    req(input$scale)

    if (identical(input$scale, "raw")) {
      raw_counts()
    } else {
      lognorm_counts()
    }
  })

  # ---- Plot object ----
  counts_hist_plot_obj <- reactive({
    req(input$dataset, input$scale, input$bins)

    vals <- plot_values()
    validate(need(!is.null(vals) && length(vals) > 0, "No data available."))

    vals <- as.numeric(vals)

    # ---- Exclude zeros ----
    vals <- vals[vals > 0]
    validate(need(length(vals) > 0, "No non-zero values to display."))

    df_plot <- data.frame(value = vals)

    x_lab <- if (identical(input$scale, "raw")) {
      "Raw counts"
    } else {
      expression(log[2]~"(CPM + 1)")
    }

    title_txt <- paste(
      if (identical(input$scale, "raw"))
        "Raw counts distribution"
      else
        "Log-normalized counts distribution",
      "\u2013",
      input$dataset
    )

    ggplot2::ggplot(df_plot, ggplot2::aes(x = .data$value)) +
      ggplot2::geom_histogram(
        bins = input$bins,
        fill = "grey80",
        color = "grey30"
      ) +
      ggplot2::labs(
        x = x_lab,
        y = "Frequency",
        title = title_txt
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        text = ggplot2::element_text(size = 18),
        aspect.ratio = 1
      )
  })

  # ---- Render ----
  output$counts_hist_plot <- shiny::renderPlot({
    counts_hist_plot_obj()
  })

  # ---- Desktop/plots directory ----
  plots_dir <- file.path(path.expand("~"), "Desktop", "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # ---- Download handler ----
  output$download_counts_hist <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "counts_hist_",
        input$dataset, "_",
        input$scale, "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {

      desktop_file <- file.path(
        plots_dir,
        paste0(
          "counts_hist_",
          input$dataset, "_",
          input$scale, "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = counts_hist_plot_obj(),
        width = 7,
        height = 7,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = counts_hist_plot_obj(),
        width = 7,
        height = 7,
        dpi = 300
      )
    }
  )
}
