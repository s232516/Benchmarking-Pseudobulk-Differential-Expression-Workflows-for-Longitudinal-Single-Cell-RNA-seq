mod_ci_hist_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 300,
      open = "always",

      sidebar = bslib::accordion(
        id = ns("ci_accordion"),

        bslib::accordion_panel(
          title = "CI coverage settings",
          icon = shiny::icon("sliders"),

          shiny::radioButtons(
            inputId = ns("sim_group"),
            label = "Simulation group:",
            choices = c(
              "All simulations" = "all",
              "Original simulations" = "old",
              "New simulations" = "new"
            ),
            selected = "new"
          ),

          shiny::conditionalPanel(
            condition = "input.sim_group == 'new'",
            ns = ns,

            shiny::selectInput(
              inputId = ns("scenario_filter"),
              label = "Scenario:",
              choices = c("All" = "all"),
              selected = "all"
            ),

            shiny::selectInput(
              inputId = ns("n_subj_filter"),
              label = "Subjects per group:",
              choices = c("All" = "all"),
              selected = "all"
            ),

            shiny::selectInput(
              inputId = ns("seed_filter"),
              label = "Seed:",
              choices = c("All" = "all"),
              selected = "all"
            )
          ),

          shiny::selectInput(
            inputId = ns("method"),
            label = "DE method:",
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("formula"),
            label = "Model formula:",
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("norm"),
            label = "Normalization:",
            choices = c("All" = "all"),
            selected = "all"
          ),

          shiny::selectInput(
            inputId = ns("metric"),
            label = "Coverage metric:",
            choices = c(
              "All genes" = "prop_ci_total",
              "Significant genes" = "prop_ci_sig",
              "Not significant genes" = "prop_ci_not_sig"
            ),
            selected = "prop_ci_total"
          ),

          shiny::sliderInput(
            inputId = ns("bins"),
            label = "Histogram bins:",
            min = 5,
            max = 80,
            value = 30
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
                shiny::textOutput(ns("ci_hist_header"), inline = TRUE),
                shiny::downloadButton(
                  outputId = ns("download_ci_hist"),
                  label = "",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("ci_hist_plot"),
                height = "550px",
                width = "100%"
              )
            )
          ),

          bslib::card(
            bslib::card_header("Coverage values"),
            bslib::card_body(
              DT::DTOutput(ns("ci_table"))
            )
          )
        )
      )
    )
  )
}


mod_ci_hist_server <- function(input, output, session, data) {

  safe_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
    gsub("_+", "_", x)
  }

  divide_safe <- function(num, den) {
    dplyr::if_else(
      !is.na(den) & den > 0,
      as.numeric(num) / as.numeric(den),
      NA_real_
    )
  }

  scenario_order <- c(
    "global_null",
    "time_only_null",
    "weak",
    "medium",
    "strong"
  )

  base_data <- shiny::reactive({
    df <- as.data.frame(data)
    n <- nrow(df)

    # -------------------------------------------------------------------------
    # Standardise identifiers
    # -------------------------------------------------------------------------

    if (!"analysis" %in% names(df)) {
      if ("analysis_join" %in% names(df)) {
        df$analysis <- df$analysis_join
      } else if ("sim" %in% names(df)) {
        df$analysis <- df$sim
      } else {
        df$analysis <- NA_character_
      }
    }

    if (!"analysis_join" %in% names(df)) {
      df$analysis_join <- df$analysis
    }

    if (!"sim" %in% names(df)) {
      df$sim <- df$analysis
    }

    if (!"norm" %in% names(df)) {
      df$norm <- NA
    }

    df <- df |>
      dplyr::mutate(
        analysis = as.character(analysis),
        analysis_join = as.character(analysis_join),
        sim = as.character(sim),
        simulation_id = dplyr::coalesce(
          as.character(analysis_join),
          as.character(analysis),
          as.character(sim)
        )
      )

    # -------------------------------------------------------------------------
    # Standardise scenario / n / seed
    # -------------------------------------------------------------------------

    if (!"scenario" %in% names(df)) {
      scenario_x <- if ("scenario.x" %in% names(df)) {
        as.character(df$scenario.x)
      } else {
        rep(NA_character_, n)
      }

      scenario_y <- if ("scenario.y" %in% names(df)) {
        as.character(df$scenario.y)
      } else {
        rep(NA_character_, n)
      }

      df$scenario <- dplyr::coalesce(scenario_x, scenario_y)
    }

    if (!"nSubjsPerGroup" %in% names(df)) {
      n_x <- if ("nSubjsPerGroup.x" %in% names(df)) {
        as.numeric(df$nSubjsPerGroup.x)
      } else {
        rep(NA_real_, n)
      }

      n_y <- if ("nSubjsPerGroup.y" %in% names(df)) {
        as.numeric(df$nSubjsPerGroup.y)
      } else {
        rep(NA_real_, n)
      }

      df$nSubjsPerGroup <- dplyr::coalesce(n_x, n_y)
    }

    if (!"seed" %in% names(df)) {
      seed_x <- if ("seed.x" %in% names(df)) {
        as.numeric(df$seed.x)
      } else {
        rep(NA_real_, n)
      }

      seed_y <- if ("seed.y" %in% names(df)) {
        as.numeric(df$seed.y)
      } else {
        rep(NA_real_, n)
      }

      df$seed <- dplyr::coalesce(seed_x, seed_y)
    }

    # Parse fair names if scenario/n/seed are missing.
    fair_name <- dplyr::coalesce(
      as.character(df$analysis_join),
      as.character(df$analysis),
      as.character(df$sim)
    )

    parsed <- regexec(
      "^fair_(.+)_n([0-9]+)_seed([0-9]+).*",
      fair_name
    )

    parsed <- regmatches(fair_name, parsed)

    parsed_scenario <- rep(NA_character_, n)
    parsed_n <- rep(NA_real_, n)
    parsed_seed <- rep(NA_real_, n)

    for (i in seq_along(parsed)) {
      if (length(parsed[[i]]) == 4) {
        parsed_scenario[i] <- parsed[[i]][2]
        parsed_n[i] <- as.numeric(parsed[[i]][3])
        parsed_seed[i] <- as.numeric(parsed[[i]][4])
      }
    }

    df$scenario <- dplyr::coalesce(
      as.character(df$scenario),
      parsed_scenario
    )

    df$nSubjsPerGroup <- dplyr::coalesce(
      as.numeric(df$nSubjsPerGroup),
      parsed_n
    )

    df$seed <- dplyr::coalesce(
      as.numeric(df$seed),
      parsed_seed
    )

    # New simulations are the fair/scenario simulations.
    # Old simulations are the original rescueSim-style ones.
    df <- df |>
      dplyr::mutate(
        .is_new = !is.na(scenario) | grepl("^fair_", simulation_id),
        .is_old = !.is_new
      )

    # -------------------------------------------------------------------------
    # Build coverage metrics
    # -------------------------------------------------------------------------

    if (!"prop_ci_total" %in% names(df)) {
      if (all(c("in_ci_total", "total_genes") %in% names(df))) {
        df$prop_ci_total <- divide_safe(df$in_ci_total, df$total_genes)
      } else {
        df$prop_ci_total <- NA_real_
      }
    }

    if (!"prop_ci_up" %in% names(df)) {
      if (all(c("in_ci_up", "total_up") %in% names(df))) {
        df$prop_ci_up <- divide_safe(df$in_ci_up, df$total_up)
      } else {
        df$prop_ci_up <- NA_real_
      }
    }

    if (!"prop_ci_down" %in% names(df)) {
      if (all(c("in_ci_down", "total_down") %in% names(df))) {
        df$prop_ci_down <- divide_safe(df$in_ci_down, df$total_down)
      } else {
        df$prop_ci_down <- NA_real_
      }
    }

    if (!"prop_ci_not_sig" %in% names(df)) {
      if (all(c("in_ci_not_sig", "total_not_sig") %in% names(df))) {
        df$prop_ci_not_sig <- divide_safe(
          df$in_ci_not_sig,
          df$total_not_sig
        )
      } else {
        df$prop_ci_not_sig <- NA_real_
      }
    }

    # Significant = up + down together.
    # Significant = up + down together.
    # This creates both the count columns and the proportion column.
    if (all(c("in_ci_up", "in_ci_down", "total_up", "total_down") %in% names(df))) {
      df$total_sig <- as.numeric(df$total_up) + as.numeric(df$total_down)
      df$in_ci_sig <- as.numeric(df$in_ci_up) + as.numeric(df$in_ci_down)

      df$prop_ci_sig <- divide_safe(df$in_ci_sig, df$total_sig)
    } else {
      if (!"total_sig" %in% names(df)) {
        df$total_sig <- NA_real_
      }

      if (!"in_ci_sig" %in% names(df)) {
        df$in_ci_sig <- NA_real_
      }

      if (!"prop_ci_sig" %in% names(df)) {
        if (all(c("prop_ci_up", "prop_ci_down") %in% names(df))) {
          df$prop_ci_sig <- rowMeans(
            df[, c("prop_ci_up", "prop_ci_down")],
            na.rm = TRUE
          )

          df$prop_ci_sig[is.nan(df$prop_ci_sig)] <- NA_real_
        } else {
          df$prop_ci_sig <- NA_real_
        }
      }
    }

    df
  })

  # ---------------------------------------------------------------------------
  # Update scenario / n / seed filters for new simulations
  # ---------------------------------------------------------------------------

  shiny::observe({
    df <- base_data() |>
      dplyr::filter(.is_new)

    scenarios <- df |>
      dplyr::filter(!is.na(scenario)) |>
      dplyr::pull(scenario) |>
      unique()

    scenarios <- c(
      intersect(scenario_order, scenarios),
      setdiff(sort(scenarios), scenario_order)
    )

    current_scenario <- input$scenario_filter

    selected_scenario <- if (
      !is.null(current_scenario) &&
      current_scenario %in% c("all", scenarios)
    ) {
      current_scenario
    } else if ("weak" %in% scenarios) {
      "weak"
    } else {
      "all"
    }

    shiny::updateSelectInput(
      session,
      "scenario_filter",
      choices = c("All" = "all", scenarios),
      selected = selected_scenario
    )

    df_n <- df

    if (!is.null(selected_scenario) && selected_scenario != "all") {
      df_n <- df_n |>
        dplyr::filter(scenario == selected_scenario)
    }

    n_subj <- df_n |>
      dplyr::filter(!is.na(nSubjsPerGroup)) |>
      dplyr::pull(nSubjsPerGroup) |>
      unique() |>
      sort()

    current_n <- input$n_subj_filter

    selected_n <- if (
      !is.null(current_n) &&
      current_n %in% c("all", as.character(n_subj))
    ) {
      current_n
    } else {
      "all"
    }

    shiny::updateSelectInput(
      session,
      "n_subj_filter",
      choices = c("All" = "all", as.character(n_subj)),
      selected = selected_n
    )

    df_seed <- df_n

    if (!is.null(selected_n) && selected_n != "all") {
      df_seed <- df_seed |>
        dplyr::filter(nSubjsPerGroup == as.numeric(selected_n))
    }

    seeds <- df_seed |>
      dplyr::filter(!is.na(seed)) |>
      dplyr::pull(seed) |>
      unique() |>
      sort()

    current_seed <- input$seed_filter

    selected_seed <- if (
      !is.null(current_seed) &&
      current_seed %in% c("all", as.character(seeds))
    ) {
      current_seed
    } else {
      "all"
    }

    shiny::updateSelectInput(
      session,
      "seed_filter",
      choices = c("All" = "all", as.character(seeds)),
      selected = selected_seed
    )
  })

  # ---------------------------------------------------------------------------
  # Filter data by simulation group and new-simulation filters
  # ---------------------------------------------------------------------------

  group_data <- shiny::reactive({
    shiny::req(input$sim_group)

    df <- base_data()

    if (input$sim_group == "old") {
      df <- df |>
        dplyr::filter(.is_old)
    }

    if (input$sim_group == "new") {
      df <- df |>
        dplyr::filter(.is_new)

      if (!is.null(input$scenario_filter) && input$scenario_filter != "all") {
        df <- df |>
          dplyr::filter(scenario == input$scenario_filter)
      }

      if (!is.null(input$n_subj_filter) && input$n_subj_filter != "all") {
        df <- df |>
          dplyr::filter(nSubjsPerGroup == as.numeric(input$n_subj_filter))
      }

      if (!is.null(input$seed_filter) && input$seed_filter != "all") {
        df <- df |>
          dplyr::filter(seed == as.numeric(input$seed_filter))
      }
    }

    df
  })

  output$ci_hist_header <- shiny::renderText({
    if (is.null(input$sim_group)) {
      return("CI coverage")
    }

    if (input$sim_group == "new") {
      return(
        paste0(
          "New simulations CI coverage",
          " | scenario: ", input$scenario_filter,
          " | n: ", input$n_subj_filter,
          " | seed: ", input$seed_filter
        )
      )
    }

    switch(
      input$sim_group,
      "old" = "Original simulations CI coverage",
      "all" = "All simulations CI coverage",
      "CI coverage"
    )
  })

  # ---------------------------------------------------------------------------
  # Update method selector
  # ---------------------------------------------------------------------------

  shiny::observeEvent(
    {
      input$sim_group
      input$scenario_filter
      input$n_subj_filter
      input$seed_filter
    },
    {
      df <- group_data()

      methods <- df |>
        dplyr::filter(!is.na(method)) |>
        dplyr::pull(method) |>
        unique() |>
        sort()

      shiny::updateSelectInput(
        session,
        "method",
        choices = methods,
        selected = if (length(methods) > 0) methods[1] else NULL
      )
    },
    ignoreInit = FALSE
  )

  # ---------------------------------------------------------------------------
  # Update formula and norm selectors
  # ---------------------------------------------------------------------------

  shiny::observeEvent(
    {
      input$method
      input$sim_group
      input$scenario_filter
      input$n_subj_filter
      input$seed_filter
    },
    {
      shiny::req(input$method)

      df <- group_data() |>
        dplyr::filter(method == input$method)

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
        avail <- sort(unique(as.character(df$norm[!is.na(df$norm)])))

        choices <- c("All" = "all")

        if ("TRUE" %in% avail) {
          choices <- c(choices, "Normalized" = "TRUE")
        }

        if ("FALSE" %in% avail) {
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
  # Plot/table data
  # ---------------------------------------------------------------------------

  plot_data <- shiny::reactive({
    shiny::req(input$method, input$formula, input$metric)

    df <- group_data() |>
      dplyr::filter(
        method == input$method,
        formula == input$formula
      )

    if (!is.null(input$norm) && input$norm != "all" && "norm" %in% names(df)) {
      df <- df |>
        dplyr::filter(as.character(norm) == input$norm)
    }

    shiny::validate(
      shiny::need(
        input$metric %in% names(df),
        paste("Metric not found:", input$metric)
      )
    )

    out <- df |>
      dplyr::mutate(
        value = .data[[input$metric]]
      ) |>
      dplyr::filter(!is.na(value)) |>
      dplyr::select(
        simulation_id,
        analysis,
        analysis_join,
        scenario,
        nSubjsPerGroup,
        seed,
        method,
        formula,
        norm,
        dplyr::any_of(c(
          "conf_level",
          "total_genes",
          "in_ci_total",
          "total_up",
          "in_ci_up",
          "total_down",
          "in_ci_down",
          "total_sig",
          "in_ci_sig",
          "total_not_sig",
          "in_ci_not_sig"
        )),
        prop_ci_total,
        prop_ci_sig,
        prop_ci_not_sig,
        value
      )

    shiny::validate(
      shiny::need(
        nrow(out) > 0,
        paste0(
          "No CI coverage values found for metric '",
          input$metric,
          "'. Check that ci_data contains this metric and the selected filters are not empty."
        )
      )
    )

    out
  })

  # ---------------------------------------------------------------------------
  # Plot
  # ---------------------------------------------------------------------------

  ci_hist_plot_obj <- shiny::reactive({
    df <- plot_data()

    metric_label <- switch(
      input$metric,
      "prop_ci_total" = "All genes",
      "prop_ci_sig" = "Significant genes",
      "prop_ci_not_sig" = "Not significant genes",
      input$metric
    )

    subtitle_text <- paste0(
      "Method: ", input$method,
      " | Formula: ", input$formula
    )

    if (input$sim_group == "new") {
      subtitle_text <- paste0(
        subtitle_text,
        " | scenario: ", input$scenario_filter,
        " | n: ", input$n_subj_filter,
        " | seed: ", input$seed_filter
      )
    }

    conf_line <- if ("conf_level" %in% names(df)) {
      suppressWarnings(unique(df$conf_level))
    } else {
      numeric(0)
    }

    conf_line <- conf_line[is.finite(conf_line)]

    p <- ggplot2::ggplot(df, ggplot2::aes(x = value)) +
      ggplot2::geom_histogram(
        bins = input$bins,
        fill = "grey80",
        color = "grey30"
      ) +
      ggplot2::scale_x_continuous(
        breaks = scales::breaks_extended(6),
        labels = scales::percent_format(accuracy = 1)
      ) +
      ggplot2::labs(
        x = "CI coverage proportion",
        y = "Number of simulation datasets",
        title = paste("CI coverage distribution -", metric_label),
        subtitle = subtitle_text
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18),
        axis.title.x = ggplot2::element_text(size = 20),
        axis.title.y = ggplot2::element_text(size = 20),
        axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14)
      )

    if (length(conf_line) == 1) {
      p <- p +
        ggplot2::geom_vline(
          xintercept = conf_line,
          linetype = "dashed",
          color = "grey40"
        )
    }

    p
  })

  output$ci_hist_plot <- shiny::renderPlot({
    ci_hist_plot_obj()
  })

  output$ci_table <- DT::renderDT({
    table_df <- plot_data() |>
      dplyr::mutate(
        prop_ci_total = round(prop_ci_total, 4),
        prop_ci_sig = round(prop_ci_sig, 4),
        prop_ci_not_sig = round(prop_ci_not_sig, 4),
        value = round(value, 4)
      ) |>
      dplyr::arrange(
        scenario,
        nSubjsPerGroup,
        seed,
        simulation_id
      )

    DT::datatable(
      table_df,
      options = list(
        pageLength = 15,
        scrollX = TRUE
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

  output$download_ci_hist <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "ci_hist_",
        input$sim_group, "_",
        safe_filename(input$scenario_filter), "_",
        safe_filename(input$n_subj_filter), "_",
        safe_filename(input$seed_filter), "_",
        safe_filename(input$method), "_",
        safe_filename(input$formula), "_",
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
          input$sim_group, "_",
          safe_filename(input$scenario_filter), "_",
          safe_filename(input$n_subj_filter), "_",
          safe_filename(input$seed_filter), "_",
          safe_filename(input$method), "_",
          safe_filename(input$formula), "_",
          input$metric, "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = ci_hist_plot_obj(),
        width = 8,
        height = 7,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = ci_hist_plot_obj(),
        width = 8,
        height = 7,
        dpi = 300
      )
    }
  )
}
