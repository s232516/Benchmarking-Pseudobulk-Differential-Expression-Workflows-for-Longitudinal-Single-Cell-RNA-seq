#' mod_ci_hist UI Function
#'
#' @description A shiny Module for CI coverage histograms.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_ci_hist_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      sidebar = bslib::accordion(
        id = ns("ci_accordion"),

        bslib::accordion_panel(
          title = "CI coverage settings",
          icon = shiny::icon("sliders"),

          shiny::selectInput(
            ns("method"),
            "Select DE method:",
            choices = NULL
          ),

          shiny::selectInput(
            ns("formula"),
            "Select model formula:",
            choices = NULL
          ),

          shiny::selectInput(
            ns("metric"),
            "Coverage metric:",
            choices = c(
              "All genes"        = "prop_ci_total",
              "Upregulated"      = "prop_ci_up",
              "Downregulated"    = "prop_ci_down",
              "Not significant"  = "prop_ci_not_sig"
            ),
            selected = "prop_ci_total"
          ),

          shiny::sliderInput(
            ns("bins"),
            "Histogram bins:",
            min = 5, max = 80, value = 40
          )
        )
      ),

      bslib::layout_columns(
        col_widths = 12,

        bslib::nav_panel(
          title = "CI Coverage Histogram",
          icon = shiny::icon("chart-bar"),

          bslib::card(
            bslib::card_header(
              div(
                style = "display:flex; justify-content:space-between; align-items:center;",
                "RescueSim CI coverage",
                shiny::downloadButton(
                  outputId = ns("download_ci_hist"),
                  label = "Download",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(ns("ci_hist_plot"), height = "550px")
            )
          )
        )
      )
    )
  )
}

#' mod_ci_hist Server Functions
#'
#' @noRd
mod_ci_hist_server <- function(input, output, session, data) {

  # ---- Only rescueSim datasets ----
  rescue_data <- reactive({
    data |>
      dplyr::filter(grepl("^rescueSim", analysis))
  })

  # ---- Populate method selector ----
  observe({
    df <- rescue_data()
    req(nrow(df) > 0)

    shiny::updateSelectInput(
      session,
      "method",
      choices = unique(df$method),
      selected = unique(df$method)[1]
    )
  })

  # ---- Populate formula selector ----
  observeEvent(input$method, {
    req(input$method)

    df <- rescue_data() |>
      dplyr::filter(method == input$method)

    shiny::updateSelectInput(
      session,
      "formula",
      choices = unique(df$formula),
      selected = unique(df$formula)[1]
    )
  })

  # ---- Filtered dataset used for plotting ----
  plot_data <- reactive({
    req(input$method, input$formula, input$metric)

    rescue_data() |>
      dplyr::filter(
        method == input$method,
        formula == input$formula
      ) |>
      dplyr::select(
        analysis,
        value = !!rlang::sym(input$metric)
      ) |>
      dplyr::filter(!is.na(value))
  })

  # ---- Plot object (single source of truth) ----
  ci_hist_plot_obj <- reactive({
    df <- plot_data()
    validate(need(nrow(df) > 0, "No data available for this selection."))

    # exact x values present in the data
    x_vals <- sort(unique(df$value))

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(
        bins = input$bins,
        fill = "grey80",
        color = "grey30"
      ) +
      ggplot2::geom_vline(
        xintercept = x_vals,
        color = "#C9980B",
        linewidth = 1,
        alpha = 0.8
      ) +
      ggplot2::scale_x_continuous(
        breaks = x_vals,
        labels = round(x_vals, 3)
      ) +
      ggplot2::labs(
        x = "CI coverage proportion",
        y = "Number of rescueSim datasets",
        title = paste("CI coverage distribution –", input$metric)
      ) +
      ggplot2::theme_bw(base_size = 16)  +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18),
        axis.text.x = ggplot2::element_text(
          angle = 45,
          hjust = 1
        )
      )
  })

  # ---- Render ----
  output$ci_hist_plot <- shiny::renderPlot({
    ci_hist_plot_obj()
  })

  # ---- Desktop/plots directory ----
  plots_dir <- file.path(path.expand("~"), "Desktop", "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # ---- Download ----
  output$download_ci_hist <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "ci_hist_",
        input$method, "_",
        gsub("[^A-Za-z0-9_]+", "_", input$formula), "_",
        input$metric, "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {

      desktop_file <- file.path(
        plots_dir,
        paste0(
          "ci_hist_",
          input$method, "_",
          gsub("[^A-Za-z0-9_]+", "_", input$formula), "_",
          input$metric, "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = ci_hist_plot_obj(),
        width = 8,
        height = 6,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = ci_hist_plot_obj(),
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  )
}
