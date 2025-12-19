#' mod_volcano UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
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
            bslib::card_header("RescueSim CI coverage"),
            bslib::card_body(
              shiny::plotOutput(ns("ci_hist_plot"), height = "550px")
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
mod_ci_hist_server <- function(input, output, session, data) {

  ns <- session$ns

  # Only rescueSim datasets
  rescue_data <- reactive({
    data |>
      dplyr::filter(grepl("^rescueSim", analysis))
  })

  # Populate method selector
  observe({
    shiny::updateSelectInput(
      session,
      "method",
      choices = unique(rescue_data()$method),
      selected = unique(rescue_data()$method)[1]
    )
  })

  # Populate formula selector
  observeEvent(input$method, {
    df <- rescue_data() |>
      dplyr::filter(method == input$method)

    shiny::updateSelectInput(
      session,
      "formula",
      choices = unique(df$formula),
      selected = unique(df$formula)[1]
    )
  })

  # Filtered dataset used for plotting
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

  output$ci_hist_plot <- shiny::renderPlot({
    req(plot_data())

    df <- plot_data()

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(
        bins = input$bins,
        fill = "grey80",
        color = "grey30"
      ) +
      ggplot2::geom_vline(
        ggplot2::aes(xintercept = value),
        color = "#C9980B",
        linewidth = 1,
        alpha = 0.8
      ) + ggplot2::scale_x_continuous(
        breaks = scales::breaks_pretty(n = 10)
      ) +
      ggplot2::labs(
        x = "CI coverage proportion",
        y = "Number of rescueSim datasets",
        title = paste(
          "CI coverage distribution –",
          input$metric
        )
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18))
  })
}


