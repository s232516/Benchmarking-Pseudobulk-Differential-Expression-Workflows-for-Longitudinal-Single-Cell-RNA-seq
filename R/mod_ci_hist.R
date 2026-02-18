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
            ns("analysis_group"),
            "Select dataset:",
            choices = c(
              "splatPop"        = "splatpop",
              "rescueSim"       = "rescueSim",
              "rescueSim (norm)" = "rescueSim_norm"
            ),
            selected = "rescueSim"
          ),

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
                shiny::textOutput(ns("ci_hist_header"), inline = TRUE),
                shiny::downloadButton(
                  outputId = ns("download_ci_hist"),
                  label = "Download",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                ns("ci_hist_plot"),
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

mod_ci_hist_server <- function(input, output, session, data) {

  # --- Map selection -> filter rule ---
  filtered_data <- shiny::reactive({
    shiny::req(input$analysis_group)

    d <- data

    if (input$analysis_group == "splatpop") {
      d <- dplyr::filter(d, grepl("^splatpop", analysis))
    } else if (input$analysis_group == "rescueSim") {
      # exclude norm
      d <- dplyr::filter(d, grepl("^rescueSim", analysis) & !grepl("\\(norm\\)", analysis))
    } else if (input$analysis_group == "rescueSim_norm") {
      d <- dplyr::filter(d, grepl("^rescueSim", analysis) & grepl("\\norm\\", analysis))
    }
    print(d)

    d
  })

  # --- Header text (dynamic) ---
  output$ci_hist_header <- shiny::renderText({
    shiny::req(input$analysis_group)
    switch(
      input$analysis_group,
      "splatpop" = "splatPop CI coverage",
      "rescueSim" = "rescueSim CI coverage",
      "rescueSim_norm" = "rescueSim (norm) CI coverage"
    )
  })

  # --- Populate method selector ---
  shiny::observe({
    df <- filtered_data()
    shiny::req(nrow(df) > 0)

    shiny::updateSelectInput(
      session,
      "method",
      choices = unique(df$method),
      selected = unique(df$method)[1]
    )
  })

  # --- Populate formula selector ---
  shiny::observeEvent(input$method, {
    shiny::req(input$method)

    df <- filtered_data() |>
      dplyr::filter(method == input$method)

    shiny::updateSelectInput(
      session,
      "formula",
      choices = unique(df$formula),
      selected = unique(df$formula)[1]
    )
  })

  # --- Filtered dataset used for plotting ---
  plot_data <- shiny::reactive({
    shiny::req(input$method, input$formula, input$metric)
    print('hola')
    print(filtered_data())

    filtered_data() |>
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

  # --- Plot object ---
  ci_hist_plot_obj <- shiny::reactive({
    df <- plot_data()
    shiny::validate(shiny::need(nrow(df) > 0, "No data available for this selection."))

    x_vals <- sort(unique(df$value))

    n_label <- switch(
      input$analysis_group,
      "splatpop" = "Number of splatPop datasets",
      "rescueSim" = "Number of rescueSim datasets",
      "rescueSim_norm" = "Number of rescueSim (norm) datasets"
    )

    title_prefix <- switch(
      input$analysis_group,
      "splatpop" = "CI coverage distribution – splatPop",
      "rescueSim" = "CI coverage distribution – rescueSim",
      "rescueSim_norm" = "CI coverage distribution – rescueSim (norm)"
    )

    ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(
        bins = input$bins,
        fill = "grey80",
        color = "grey30"
      ) +
      ggplot2::scale_x_continuous(
        breaks = x_vals,
        labels = round(x_vals, 3)
      ) +
      ggplot2::labs(
        x = "CI coverage proportion",
        y = n_label,
        title = paste(title_prefix, "(", input$metric, ")")
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18),
        axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)
      )
  })

  # --- Render ---
  output$ci_hist_plot <- shiny::renderPlot({
    ci_hist_plot_obj()
  })

  # --- Desktop/plots directory ---
  plots_dir <- file.path(path.expand("~"), "Desktop", "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)

  # --- Download ---
  output$download_ci_hist <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "ci_hist_",
        input$analysis_group, "_",
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
          input$analysis_group, "_",
          input$method, "_",
          gsub("[^A-Za-z0-9_]+", "_", input$formula), "_",
          input$metric, "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(desktop_file, ci_hist_plot_obj(), width = 8, height = 6, dpi = 300)
      ggplot2::ggsave(file,         ci_hist_plot_obj(), width = 8, height = 6, dpi = 300)
    }
  )
}
