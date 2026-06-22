mod_volcano_ui <- function(id) {
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

          shiny::radioButtons(
            inputId = ns("sim_source"),
            label = "Simulation set:",
            choices = c(
              "Fair benchmark" = "fair",
              "Original simulations" = "old"
            ),
            selected = "fair"
          ),

          shiny::conditionalPanel(
            condition = "input.sim_source == 'old'",
            ns = ns,

            shiny::selectInput(
              inputId = ns("data"),
              label = bslib::tooltip(
                trigger = shiny::tagList("Select dataset",
                                         shiny::icon("circle-info")),
                "Select which original simulation dataset to analyse."
              ),
              choices = NULL
            )
          ),

          shiny::conditionalPanel(
            condition = "input.sim_source == 'fair'",
            ns = ns,

            shiny::selectInput(
              inputId = ns("scenario"),
              label = "Scenario:",
              choices = NULL
            ),

            shiny::selectInput(
              inputId = ns("n_subj"),
              label = "Subjects per group:",
              choices = NULL
            ),

            shiny::selectInput(
              inputId = ns("seed"),
              label = "Seed:",
              choices = NULL
            )
          ),

          shiny::selectInput(
            inputId = ns("de_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select DE method",
                                       shiny::icon("circle-info")),
              "Select the differential expression method."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select formula",
                                       shiny::icon("circle-info")),
              "Select the model formula."
            ),
            choices = NULL
          ),

          # NEW: normalization selector
          shiny::selectInput(
            inputId = ns("norm"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Normalization", shiny::icon("circle-info")),
              "Filter results by whether log-normalization was used (only applies when the 'norm' column exists)."
            ),
            choices = c("All" = "all"),   # server will update this
            selected = "all"
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

      bslib::layout_columns(
        col_widths = 12,

        div(
          style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",

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


mod_volcano_server <- function(input, output, session, data) {

   # Choose dataset column

  dataset_col <- if ("analysis_join" %in% colnames(data)) {
    "analysis_join"
  } else {
    "analysis"
  }
  # Base cleaned data

  base_data <- shiny::reactive({
    df <- as.data.frame(data)

    n <- nrow(df)

    # Fix scenario.x / scenario.y problem
    if (!"scenario" %in% names(df)) {
      scenario_x <- if ("scenario.x" %in% names(df)) as.character(df$scenario.x)
      else rep(NA_character_, n)
      scenario_y <- if ("scenario.y" %in% names(df)) as.character(df$scenario.y)
      else rep(NA_character_, n)

      df$scenario <- dplyr::coalesce(scenario_x, scenario_y)
    }

    # Fix nSubjsPerGroup.x / nSubjsPerGroup.y problem
    if (!"nSubjsPerGroup" %in% names(df)) {
      n_x <- if ("nSubjsPerGroup.x" %in% names(df))
        as.numeric(df$nSubjsPerGroup.x) else rep(NA_real_, n)
      n_y <- if ("nSubjsPerGroup.y" %in% names(df))
        as.numeric(df$nSubjsPerGroup.y) else rep(NA_real_, n)

      df$nSubjsPerGroup <- dplyr::coalesce(n_x, n_y)
    }

    # Fix seed.x / seed.y problem
    if (!"seed" %in% names(df)) {
      seed_x <- if ("seed.x" %in% names(df)) as.numeric(df$seed.x)
      else rep(NA_real_, n)
      seed_y <- if ("seed.y" %in% names(df)) as.numeric(df$seed.y)
      else rep(NA_real_, n)

      df$seed <- dplyr::coalesce(seed_x, seed_y)
    }

    if ("trueDE" %in% names(df)) {
      df$trueDE <- as.logical(df$trueDE)
    }

    df
  })


  # Update simulation selectors

  observe({
    df <- base_data()

    old_datasets <- df |>
      dplyr::filter(!grepl("^fair_", .data[[dataset_col]])) |>
      dplyr::pull(.data[[dataset_col]]) |>
      unique() |>
      sort()

    fair_df <- df |>
      dplyr::filter(grepl("^fair_", .data[[dataset_col]]))

    scenarios <- fair_df |>
      dplyr::filter(!is.na(scenario)) |>
      dplyr::pull(scenario) |>
      unique() |>
      sort()

    n_subj <- fair_df |>
      dplyr::filter(!is.na(nSubjsPerGroup)) |>
      dplyr::pull(nSubjsPerGroup) |>
      unique() |>
      sort()

    seeds <- fair_df |>
      dplyr::filter(!is.na(seed)) |>
      dplyr::pull(seed) |>
      unique() |>
      sort()

    shiny::updateSelectInput(
      session,
      "data",
      choices = old_datasets,
      selected = if (length(old_datasets) > 0) old_datasets[1] else NULL
    )

    shiny::updateSelectInput(
      session,
      "scenario",
      choices = scenarios,
      selected = if (length(scenarios) > 0) scenarios[1] else NULL
    )

    shiny::updateSelectInput(
      session,
      "n_subj",
      choices = as.character(n_subj),
      selected = if (length(n_subj) > 0) as.character(n_subj[1]) else NULL
    )

    shiny::updateSelectInput(
      session,
      "seed",
      choices = as.character(seeds),
      selected = if (length(seeds) > 0) as.character(seeds[1]) else NULL
    )
  })

  # ---------------------------------------------------------------------------
  # Selected simulation data
  # ---------------------------------------------------------------------------

  selected_simulation_data <- shiny::reactive({
    df <- base_data()

    if (is.null(input$sim_source) || input$sim_source == "fair") {
      req(input$scenario, input$n_subj, input$seed)

      df |>
        dplyr::filter(
          grepl("^fair_", .data[[dataset_col]]),
          scenario == input$scenario,
          nSubjsPerGroup == as.numeric(input$n_subj),
          seed == as.numeric(input$seed)
        )
    } else {
      req(input$data)

      df |>
        dplyr::filter(.data[[dataset_col]] == input$data)
    }
  })

  selected_simulation_label <- shiny::reactive({
    df <- selected_simulation_data()

    if (nrow(df) == 0) {
      return("No simulation selected")
    }

    unique(df[[dataset_col]])[1]
  })

  # ---------------------------------------------------------------------------
  # Update method selector
  # ---------------------------------------------------------------------------

  observeEvent(
    {
      input$sim_source
      input$data
      input$scenario
      input$n_subj
      input$seed
    },
    {
      df <- selected_simulation_data()

      methods <- df |>
        dplyr::filter(!is.na(method)) |>
        dplyr::pull(method) |>
        unique() |>
        sort()

      shiny::updateSelectInput(
        session,
        "de_method",
        choices = methods,
        selected = if (length(methods) > 0) methods[1] else NULL
      )
    },
    ignoreInit = FALSE
  )

  # ---------------------------------------------------------------------------
  # Update formula and normalization selectors
  # ---------------------------------------------------------------------------

  observeEvent(
    {
      input$de_method
    },
    {
      req(input$de_method)

      df <- selected_simulation_data() |>
        dplyr::filter(method == input$de_method)

      formulas <- df |>
        dplyr::filter(!is.na(formula)) |>
        dplyr::pull(formula) |>
        unique() |>
        sort()

      shiny::updateSelectInput(
        session,
        "formula",
        choices = formulas,
        selected = if (length(formulas) > 0) formulas[1] else NULL
      )

      if ("norm" %in% names(df) && any(!is.na(df$norm))) {
        avail <- sort(unique(df$norm[!is.na(df$norm)]))

        choices <- c("All" = "all")

        if (TRUE %in% avail) {
          choices <- c(choices, "Normalized" = "TRUE")
        }

        if (FALSE %in% avail) {
          choices <- c(choices, "Not normalized" = "FALSE")
        }

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
    },
    ignoreInit = FALSE
  )

  # ---------------------------------------------------------------------------
  # Update logFC slider
  # ---------------------------------------------------------------------------

  observe({
    df <- selected_simulation_data()

    max_lfc <- max(abs(df$log2FoldChange), na.rm = TRUE)

    if (!is.finite(max_lfc)) {
      max_lfc <- 4
    }

    shiny::updateSliderInput(
      session,
      "logfc",
      min = 0,
      max = ceiling(max_lfc),
      value = 0,
      step = 0.01
    )
  })

  # ---------------------------------------------------------------------------
  # Volcano data
  # ---------------------------------------------------------------------------

  volcano_data <- shiny::reactive({
    req(input$de_method, input$formula, input$norm)

    df <- selected_simulation_data() |>
      dplyr::filter(
        method == input$de_method,
        formula == input$formula
      )

    if (input$norm != "all" && "norm" %in% names(df)) {
      df <- df |>
        dplyr::filter(norm == as.logical(input$norm))
    }

    df |>
      dplyr::mutate(
        pvalue_safe = dplyr::if_else(
          is.na(pvalue),
          NA_real_,
          pmax(pvalue, .Machine$double.xmin)
        ),
        padj_safe = dplyr::if_else(
          is.na(padj),
          NA_real_,
          pmax(padj, .Machine$double.xmin)
        ),
        neg_log10_pvalue = -log10(pvalue_safe),
        neg_log10_padj = -log10(padj_safe),
        is_true_de = dplyr::if_else(is.na(trueDE), FALSE, trueDE),
        significance = dplyr::case_when(
          !is.na(padj) &
            padj <= as.numeric(input$pvalue) &
            log2FoldChange > input$logfc ~ "Upregulated",

          !is.na(padj) &
            padj <= as.numeric(input$pvalue) &
            log2FoldChange < -input$logfc ~ "Downregulated",

          TRUE ~ "Not significant"
        )
      )
  })

  # ---------------------------------------------------------------------------
  # Volcano plot
  # ---------------------------------------------------------------------------

  volcano_plot_obj <- shiny::reactive({
    df <- as.data.frame(volcano_data())

    validate(
      shiny::need(nrow(df) > 0, "No data available for this selection.")
    )

    max_es <- ceiling(max(abs(df$log2FoldChange), na.rm = TRUE))

    if (!is.finite(max_es) || max_es == 0) {
      max_es <- 1
    }

    x_lims <- c(-max_es, max_es)

    if (isTRUE(input$highlight_truede)) {
      df$significance <- ifelse(
        df$is_true_de,
        "true DE",
        df$significance
      )

      counts <- df |>
        dplyr::count(significance) |>
        tibble::deframe()

      colors <- c(
        "true DE" = "#C9980B",
        "Upregulated" = "#b2182b",
        "Downregulated" = "#2166ac",
        "Not significant" = "#c0c0c0"
      )
    } else {
      counts <- df |>
        dplyr::count(significance) |>
        tibble::deframe()

      colors <- c(
        "Upregulated" = "#b2182b",
        "Downregulated" = "#2166ac",
        "Not significant" = "#c0c0c0"
      )
    }

    labels <- paste0(names(counts), " (", counts, ")")
    names(labels) <- names(counts)

    volcano_plot(
      data = df,
      x = "log2FoldChange",
      y = "pvalue_safe",
      color_col = "significance",
      colors = colors,
      color_labels = labels,
      labels = "gene",
      x_axis_title = "Effect size",
      y_axis_title = "P-value (-log10 scale)",
      plot_title = paste(
        "Volcano Plot -",
        selected_simulation_label(),
        "-",
        input$de_method
      ),
      x_lims = x_lims,
      vline = input$logfc,
      label_size = 5
    ) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      )
  })

  # ---------------------------------------------------------------------------
  # Render plot and table
  # ---------------------------------------------------------------------------

  output$volcano_plot <- shiny::renderPlot({
    volcano_plot_obj()
  })

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

  # ---------------------------------------------------------------------------
  # Download
  # ---------------------------------------------------------------------------

  plots_dir <- file.path(path.expand("~"), "Desktop", "plots")

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "volcano_",
        selected_simulation_label(),
        "_",
        input$de_method,
        "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0(
          "volcano_",
          selected_simulation_label(),
          "_",
          input$de_method,
          "_",
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
