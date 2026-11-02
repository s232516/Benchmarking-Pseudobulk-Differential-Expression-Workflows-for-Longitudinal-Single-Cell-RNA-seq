#' mod_tpr_line UI Function
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
mod_tpr_line_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::card(
    bslib::layout_sidebar(
      fillable = TRUE,

      sidebar = bslib::accordion(
        id = ns("tpr_accordion"),

        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),

          shiny::selectInput(
            inputId = ns("de_method"),
            label = "DE method:",
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("formula"),
            label = "Formula:",
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("norm"),
            label = "Normalization:",
            choices = c("All" = "all"),
            selected = "all"
          ),

          shiny::selectizeInput(
            inputId = ns("scenario"),
            label = "Scenario:",
            choices = NULL,
            multiple = TRUE,
            options = list(
              plugins = list("remove_button")
            )
          ),

          shiny::selectInput(
            inputId = ns("seed"),
            label = "Seed:",
            choices = c("All" = "all"),
            selected = "all"
          ),

          shiny::numericInput(
            inputId = ns("alpha"),
            label = "Adjusted p-value threshold:",
            value = 0.01,
            min = 0,
            max = 1,
            step = 0.001
          ),

          shiny::checkboxInput(
            inputId = ns("show_errorbars"),
            label = "Show seed variability",
            value = TRUE
          )
        )
      ),

      bslib::layout_columns(
        col_widths = 12,

        bslib::card(
          fill = FALSE,
          fill_screen = TRUE,

          bslib::card_header(
            div(
              style = "display: flex; justify-content: space-between; align-items: center;",
              "TPR by sample size",
              shiny::downloadButton(
                outputId = ns("download_plot"),
                label = "",
                class = "clickable-button-light"
              )
            )
          ),

          bslib::card_body(
            div(
              style = "width: 800px; height: 550px; overflow: auto;",
              shiny::plotOutput(
                outputId = ns("tpr_line_plot"),
                height = "550px",
                width = "800px"
              )
            )
          )
        ),

        bslib::card(
          bslib::card_header("TPR values"),
          bslib::card_body(
            DT::DTOutput(ns("tpr_table"))
          )
        )
      )
    )
  )
}


#' mod_tpr_line Server Function
#'
#' @noRd
mod_tpr_line_server <- function(input, output, session, data) {

  safe_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
    gsub("_+", "_", x)
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

    scenario_main <- if ("scenario" %in% names(df)) {
      as.character(df$scenario)
    } else {
      rep(NA_character_, n)
    }

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

    df$scenario <- dplyr::coalesce(
      scenario_main,
      scenario_x,
      scenario_y
    )

    n_main <- if ("nSubjsPerGroup" %in% names(df)) {
      as.numeric(df$nSubjsPerGroup)
    } else {
      rep(NA_real_, n)
    }

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

    df$nSubjsPerGroup <- dplyr::coalesce(
      n_main,
      n_x,
      n_y
    )

    seed_main <- if ("seed" %in% names(df)) {
      as.numeric(df$seed)
    } else {
      rep(NA_real_, n)
    }

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

    df$seed <- dplyr::coalesce(
      seed_main,
      seed_x,
      seed_y
    )

    # -------------------------------------------------------------------------
    # Parse fair names if metadata is missing
    # -------------------------------------------------------------------------

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

    df |>
      dplyr::mutate(
        trueDE = as.logical(trueDE)
      )
  })

  # ---------------------------------------------------------------------------
  # Only use new/fair simulations
  # ---------------------------------------------------------------------------

  new_data <- shiny::reactive({
    df <- base_data() |>
      dplyr::filter(
        !is.na(scenario),
        !is.na(nSubjsPerGroup),
        !is.na(seed)
      )

    shiny::validate(
      shiny::need(
        nrow(df) > 0,
        "No new/fair simulations found. Check that the data has scenario, nSubjsPerGroup and seed."
      )
    )

    df
  })

  # ---------------------------------------------------------------------------
  # Update method selector
  # ---------------------------------------------------------------------------

  shiny::observe({
    df <- new_data()

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
  })

  # ---------------------------------------------------------------------------
  # Update formula / norm / scenario / seed selectors
  # ---------------------------------------------------------------------------

  shiny::observeEvent(
    input$de_method,
    {
      shiny::req(input$de_method)

      df <- new_data() |>
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

      scenarios <- df |>
        dplyr::filter(!is.na(scenario)) |>
        dplyr::pull(scenario) |>
        unique()

      scenarios <- c(
        intersect(scenario_order, scenarios),
        setdiff(sort(scenarios), scenario_order)
      )

      shiny::updateSelectizeInput(
        session,
        "scenario",
        choices = scenarios,
        selected = scenarios,
        server = TRUE
      )

      seeds <- df |>
        dplyr::filter(!is.na(seed)) |>
        dplyr::pull(seed) |>
        unique() |>
        sort()

      shiny::updateSelectInput(
        session,
        "seed",
        choices = c("All" = "all", as.character(seeds)),
        selected = "all"
      )
    },
    ignoreInit = FALSE
  )

  # ---------------------------------------------------------------------------
  # Filtered data
  # ---------------------------------------------------------------------------

  filtered_data <- shiny::reactive({
    shiny::req(input$de_method, input$formula, input$norm)

    df <- new_data() |>
      dplyr::filter(
        method == input$de_method,
        formula == input$formula
      )

    if (!is.null(input$norm) && input$norm != "all") {
      df <- df |>
        dplyr::filter(as.character(norm) == input$norm)
    }

    if (!is.null(input$scenario) && length(input$scenario) > 0) {
      df <- df |>
        dplyr::filter(scenario %in% input$scenario)
    }

    if (!is.null(input$seed) && input$seed != "all") {
      df <- df |>
        dplyr::filter(seed == as.numeric(input$seed))
    }

    df
  })

  # ---------------------------------------------------------------------------
  # TPR per seed
  # ---------------------------------------------------------------------------

  tpr_by_seed <- shiny::reactive({
    df <- filtered_data()
    alpha <- if (is.null(input$alpha)) 0.01 else input$alpha

    shiny::validate(
      shiny::need("trueDE" %in% names(df), "Column trueDE not found."),
      shiny::need("padj" %in% names(df), "Column padj not found.")
    )

    df |>
      dplyr::filter(
        !is.na(trueDE),
        !is.na(padj)
      ) |>
      dplyr::mutate(
        called = padj <= alpha
      ) |>
      dplyr::group_by(
        scenario,
        nSubjsPerGroup,
        seed,
        method,
        formula,
        norm
      ) |>
      dplyr::summarise(
        alpha = alpha,
        total_genes = dplyr::n(),

        actual_DE = sum(trueDE, na.rm = TRUE),
        actual_not_DE = sum(!trueDE, na.rm = TRUE),

        TP = sum(trueDE & called, na.rm = TRUE),
        FP = sum(!trueDE & called, na.rm = TRUE),
        TN = sum(!trueDE & !called, na.rm = TRUE),
        FN = sum(trueDE & !called, na.rm = TRUE),

        TPR = dplyr::if_else(
          actual_DE > 0,
          TP / actual_DE,
          NA_real_
        ),

        .groups = "drop"
      ) |>
      dplyr::filter(!is.na(TPR))
  })

  # ---------------------------------------------------------------------------
  # Mean TPR over seeds
  # ---------------------------------------------------------------------------

  tpr_summary <- shiny::reactive({
    df <- tpr_by_seed()

    shiny::validate(
      shiny::need(
        nrow(df) > 0,
        "No TPR values available for this selection."
      )
    )

    df |>
      dplyr::group_by(
        scenario,
        nSubjsPerGroup
      ) |>
      dplyr::summarise(
        mean_TPR = mean(TPR, na.rm = TRUE),
        sd_TPR = stats::sd(TPR, na.rm = TRUE),
        n_seeds = dplyr::n_distinct(seed),
        se_TPR = dplyr::if_else(
          n_seeds > 1,
          sd_TPR / sqrt(n_seeds),
          NA_real_
        ),
        lower = pmax(0, mean_TPR - se_TPR),
        upper = pmin(1, mean_TPR + se_TPR),
        .groups = "drop"
      ) |>
      dplyr::arrange(
        scenario,
        nSubjsPerGroup
      )
  })

  # ---------------------------------------------------------------------------
  # Plot object
  # ---------------------------------------------------------------------------

  tpr_line_plot_obj <- shiny::reactive({
    df <- tpr_summary()

    p <- ggplot2::ggplot(
      df,
      ggplot2::aes(
        x = nSubjsPerGroup,
        y = mean_TPR,
        color = scenario,
        group = scenario
      )
    ) +
      ggplot2::geom_line(linewidth = 1) +
      ggplot2::geom_point(size = 3) +
      ggplot2::scale_x_continuous(
        breaks = sort(unique(df$nSubjsPerGroup))
      ) +
      ggplot2::scale_y_continuous(
        limits = c(0, 1),
        breaks = seq(0, 1, by = 0.2)
      ) +
      ggplot2::labs(
        title = "True positive rate by sample size",
        subtitle = paste0(
          "Method: ", input$de_method,
          " | Formula: ", input$formula,
          " | alpha: ", input$alpha
        ),
        x = "Subjects per group",
        y = "True positive rate",
        color = "Scenario"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18),
        axis.title.x = ggplot2::element_text(size = 18),
        axis.title.y = ggplot2::element_text(size = 18),
        axis.text.x = ggplot2::element_text(size = 14),
        axis.text.y = ggplot2::element_text(size = 14),
        legend.title = ggplot2::element_text(size = 14),
        legend.text = ggplot2::element_text(size = 12)
      )

    if (isTRUE(input$show_errorbars)) {
      p <- p +
        ggplot2::geom_errorbar(
          ggplot2::aes(
            ymin = lower,
            ymax = upper
          ),
          width = 0.15,
          alpha = 0.6
        )
    }

    p
  })

  # ---------------------------------------------------------------------------
  # Render plot
  # ---------------------------------------------------------------------------

  output$tpr_line_plot <- shiny::renderPlot({
    tpr_line_plot_obj()
  })

  # ---------------------------------------------------------------------------
  # Render table
  # ---------------------------------------------------------------------------

  output$tpr_table <- DT::renderDT({
    table_df <- tpr_summary() |>
      dplyr::mutate(
        mean_TPR = round(mean_TPR, 4),
        sd_TPR = round(sd_TPR, 4),
        se_TPR = round(se_TPR, 4),
        lower = round(lower, 4),
        upper = round(upper, 4)
      )

    DT::datatable(
      table_df,
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  # ---------------------------------------------------------------------------
  # Download handler
  # ---------------------------------------------------------------------------

  plots_dir <- file.path(
    path.expand("~"),
    "Desktop",
    "plots"
  )

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "tpr_line_",
        safe_filename(input$de_method),
        "_",
        safe_filename(input$formula),
        "_",
        Sys.Date(),
        ".png"
      )
    },

    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0(
          "tpr_line_",
          safe_filename(input$de_method),
          "_",
          safe_filename(input$formula),
          "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = tpr_line_plot_obj(),
        width = 8,
        height = 7,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = tpr_line_plot_obj(),
        width = 8,
        height = 7,
        dpi = 300
      )
    }
  )
}
