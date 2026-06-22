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

      sidebar = bslib::accordion(
        id = ns("roc_accordion"),

        bslib::accordion_panel(
          title = "ROC / PR settings",
          icon = shiny::icon("sliders"),

          shiny::radioButtons(
            inputId = ns("compare_mode"),
            label = "Compare:",
            choices = c(
              "DE workflows for one simulation" = "methods",
              "Simulations for one DE workflow" = "simulations"
            ),
            selected = "methods"
          ),

          shiny::conditionalPanel(
            condition = "input.compare_mode == 'methods'",
            ns = ns,

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
                inputId = ns("analysis"),
                label = "Original simulation:",
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
            )
          ),

          shiny::conditionalPanel(
            condition = "input.compare_mode == 'simulations'",
            ns = ns,

            shiny::radioButtons(
              inputId = ns("sim_compare_source"),
              label = "Simulation set:",
              choices = c(
                "Fair benchmark" = "fair",
                "Original simulations" = "old"
              ),
              selected = "fair"
            ),

            shiny::conditionalPanel(
              condition = "input.sim_compare_source == 'fair'",
              ns = ns,

              shiny::selectInput(
                inputId = ns("scenario_filter"),
                label = "Scenario:",
                choices = NULL
              ),

              shiny::selectInput(
                inputId = ns("n_subj_filter"),
                label = "Subjects per group:",
                choices = NULL
              ),

              shiny::selectInput(
                inputId = ns("seed_filter"),
                label = "Seed:",
                choices = NULL
              )
            ),

            shiny::conditionalPanel(
              condition = "input.sim_compare_source == 'old'",
              ns = ns,

              shiny::selectizeInput(
                inputId = ns("old_analyses_multi"),
                label = "Original simulations:",
                choices = NULL,
                multiple = TRUE,
                options = list(
                  plugins = list("remove_button")
                )
              )
            ),

            shiny::selectInput(
              inputId = ns("method"),
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
            )
          ),

          shiny::numericInput(
            inputId = ns("alpha"),
            label = "FDR point:",
            value = 0.01,
            min = 0,
            max = 1,
            step = 0.001
          ),

          shiny::checkboxInput(
            inputId = ns("zoom_fdr_region"),
            label = "Zoom around FDR points",
            value = FALSE
          ),

          shiny::conditionalPanel(
            condition = "input.zoom_fdr_region == true",
            ns = ns,

            shiny::numericInput(
              inputId = ns("roc_xmax"),
              label = "ROC max FPR:",
              value = 0.05,
              min = 0.001,
              max = 1,
              step = 0.005
            ),

            shiny::numericInput(
              inputId = ns("pr_xmax"),
              label = "PR max recall:",
              value = 0.20,
              min = 0.001,
              max = 1,
              step = 0.01
            )
          ),

          shiny::checkboxInput(
            inputId = ns("show_legend"),
            label = "Show legend",
            value = TRUE
          )
        )
      ),

      bslib::layout_columns(
        col_widths = 12,

        div(
          style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",

          bslib::navset_tab(
            id = ns("roc_pr_tabs"),

            bslib::nav_panel(
              title = "ROC Curve",
              icon = shiny::icon("chart-line"),

              bslib::card(
                style = "margin-bottom: 1rem;",

                bslib::card_header(
                  div(
                    style = "display:flex; justify-content:space-between; align-items:center;",
                    "ROC Curve",
                    shiny::downloadButton(
                      outputId = ns("download_roc"),
                      label = "",
                      class = "clickable-button-light"
                    )
                  )
                ),

                bslib::card_body(
                  shiny::plotOutput(
                    outputId = ns("roc_plot"),
                    height = "550px",
                    width = "100%"
                  )
                )
              ),

              bslib::card(
                bslib::card_header("ROC AUC Table"),
                bslib::card_body(
                  DT::DTOutput(ns("roc_table"))
                )
              )
            ),

            bslib::nav_panel(
              title = "PR Curve",
              icon = shiny::icon("chart-area"),

              bslib::card(
                style = "margin-bottom: 1rem;",

                bslib::card_header(
                  div(
                    style = "display:flex; justify-content:space-between; align-items:center;",
                    "PR Curve",
                    shiny::downloadButton(
                      outputId = ns("download_pr"),
                      label = "",
                      class = "clickable-button-light"
                    )
                  )
                ),

                bslib::card_body(
                  shiny::plotOutput(
                    outputId = ns("pr_plot"),
                    height = "550px",
                    width = "100%"
                  )
                )
              ),

              bslib::card(
                bslib::card_header("PR Curve Table"),
                bslib::card_body(
                  DT::DTOutput(ns("pr_table"))
                )
              )
            )
          ),

          bslib::card(
            style = "margin-top: 1rem;",
            bslib::card_header("FDR point values"),
            bslib::card_body(
              DT::DTOutput(ns("fdr_points_table"))
            )
          ),

          bslib::card(
            style = "margin-top: 1rem;",
            bslib::card_header("Confusion matrix at selected FDR"),
            bslib::card_body(
              DT::DTOutput(ns("confusion_table"))
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

  safe_filename <- function(x) {
    x <- gsub("[^A-Za-z0-9_\\-]+", "_", x)
    gsub("_+", "_", x)
  }

  wrap_labels <- function(x, width = 28) {
    vapply(
      x,
      function(z) paste(strwrap(z, width = width), collapse = "\n"),
      character(1)
    )
  }

  empty_point_tbl <- function() {
    tibble::tibble(
      curve_id = character(),
      FPR = numeric(),
      TPR = numeric(),
      recall = numeric(),
      precision = numeric(),
      real_fdr = numeric(),
      n_called = numeric(),
      tp = numeric(),
      fp = numeric()
    )
  }

  scenario_order <- c(
    "global_null",
    "time_only_null",
    "weak",
    "medium",
    "strong"
  )

  dataset_col <- if ("analysis_join" %in% colnames(data)) {
    "analysis_join"
  } else {
    "analysis"
  }

  base_data <- shiny::reactive({
    df <- as.data.frame(data)
    n <- nrow(df)

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

    df |>
      dplyr::mutate(
        trueDE = as.logical(trueDE),
        pvalue_safe = dplyr::if_else(
          is.na(pvalue),
          NA_real_,
          pmax(pvalue, .Machine$double.xmin)
        ),
        score = -log10(pvalue_safe)
      )
  })

  shiny::observe({
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
      unique()

    scenarios <- c(
      intersect(scenario_order, scenarios),
      setdiff(sort(scenarios), scenario_order)
    )

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

    default_scenario <- if ("weak" %in% scenarios) {
      "weak"
    } else if (length(scenarios) > 0) {
      scenarios[1]
    } else {
      NULL
    }

    shiny::updateSelectInput(
      session,
      "analysis",
      choices = old_datasets,
      selected = if (length(old_datasets) > 0) old_datasets[1] else NULL
    )

    shiny::updateSelectizeInput(
      session,
      "old_analyses_multi",
      choices = old_datasets,
      selected = head(old_datasets, 5),
      server = TRUE
    )

    shiny::updateSelectInput(
      session,
      "scenario",
      choices = scenarios,
      selected = default_scenario
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

    shiny::updateSelectInput(
      session,
      "scenario_filter",
      choices = c("All" = "all", scenarios),
      selected = default_scenario
    )

    shiny::updateSelectInput(
      session,
      "n_subj_filter",
      choices = c("All" = "all", as.character(n_subj)),
      selected = "all"
    )

    shiny::updateSelectInput(
      session,
      "seed_filter",
      choices = c("All" = "all", as.character(seeds)),
      selected = "all"
    )
  })

  selected_one_sim_data <- shiny::reactive({
    df <- base_data()

    if (is.null(input$sim_source) || input$sim_source == "fair") {
      shiny::req(input$scenario, input$n_subj, input$seed)

      df |>
        dplyr::filter(
          grepl("^fair_", .data[[dataset_col]]),
          scenario == input$scenario,
          nSubjsPerGroup == as.numeric(input$n_subj),
          seed == as.numeric(input$seed)
        )
    } else {
      shiny::req(input$analysis)

      df |>
        dplyr::filter(.data[[dataset_col]] == input$analysis)
    }
  })

  selected_one_sim_label <- shiny::reactive({
    df <- selected_one_sim_data()

    if (nrow(df) == 0) {
      return("no_simulation")
    }

    unique(df[[dataset_col]])[1]
  })

  selected_many_sims_data <- shiny::reactive({
    df <- base_data()

    if (is.null(input$sim_compare_source) || input$sim_compare_source == "fair") {
      out <- df |>
        dplyr::filter(grepl("^fair_", .data[[dataset_col]]))

      if (!is.null(input$scenario_filter) && input$scenario_filter != "all") {
        out <- out |>
          dplyr::filter(scenario == input$scenario_filter)
      }

      if (!is.null(input$n_subj_filter) && input$n_subj_filter != "all") {
        out <- out |>
          dplyr::filter(nSubjsPerGroup == as.numeric(input$n_subj_filter))
      }

      if (!is.null(input$seed_filter) && input$seed_filter != "all") {
        out <- out |>
          dplyr::filter(seed == as.numeric(input$seed_filter))
      }

      out
    } else {
      shiny::req(input$old_analyses_multi)

      df |>
        dplyr::filter(
          !grepl("^fair_", .data[[dataset_col]]),
          .data[[dataset_col]] %in% input$old_analyses_multi
        )
    }
  })

  shiny::observeEvent(
    {
      input$compare_mode
      input$sim_compare_source
      input$scenario_filter
      input$n_subj_filter
      input$seed_filter
      input$old_analyses_multi
    },
    {
      df <- selected_many_sims_data()

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

  shiny::observeEvent(
    {
      input$method
      input$sim_compare_source
      input$scenario_filter
      input$n_subj_filter
      input$seed_filter
      input$old_analyses_multi
    },
    {
      shiny::req(input$method)

      df <- selected_many_sims_data() |>
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

  rocpr_data <- shiny::reactive({
    if (is.null(input$compare_mode) || input$compare_mode == "methods") {
      df <- selected_one_sim_data() |>
        dplyr::mutate(curve_id = method)
    } else {
      shiny::req(input$method, input$formula, input$norm)

      df <- selected_many_sims_data() |>
        dplyr::filter(
          method == input$method,
          formula == input$formula
        )

      if (input$norm != "all" && "norm" %in% names(df)) {
        df <- df |>
          dplyr::filter(norm == as.logical(input$norm))
      }

      df <- df |>
        dplyr::mutate(
          curve_id = dplyr::if_else(
            grepl("^fair_", .data[[dataset_col]]),
            paste0(
              scenario,
              " | n=",
              nSubjsPerGroup,
              " | seed=",
              seed
            ),
            .data[[dataset_col]]
          )
        )
    }

    df |>
      dplyr::filter(
        !is.na(trueDE),
        !is.na(score),
        is.finite(score)
      ) |>
      dplyr::group_by(curve_id) |>
      dplyr::filter(dplyr::n_distinct(trueDE) == 2) |>
      dplyr::ungroup()
  })

  calc_curves_one <- function(df_one, curve_name) {
    df_one <- df_one |>
      dplyr::arrange(dplyr::desc(score))

    n_pos <- sum(df_one$trueDE, na.rm = TRUE)
    n_neg <- sum(!df_one$trueDE, na.rm = TRUE)

    if (n_pos == 0 || n_neg == 0) {
      return(NULL)
    }

    tp <- cumsum(df_one$trueDE)
    fp <- cumsum(!df_one$trueDE)

    tpr <- tp / n_pos
    fpr <- fp / n_neg

    precision <- tp / pmax(tp + fp, 1)
    recall <- tpr

    roc_df <- tibble::tibble(
      curve_id = curve_name,
      curvetype = "ROC",
      x = c(0, fpr),
      y = c(0, tpr)
    )

    pr_df <- tibble::tibble(
      curve_id = curve_name,
      curvetype = "PRC",
      x = c(0, recall),
      y = c(1, precision)
    )

    dplyr::bind_rows(roc_df, pr_df)
  }

  rocpr_curves <- shiny::reactive({
    df <- rocpr_data()

    shiny::validate(
      shiny::need(
        nrow(df) > 0,
        "No valid ROC/PR data. Null-only simulations cannot produce ROC/PR curves."
      )
    )

    split_df <- split(df, df$curve_id)

    curves <- lapply(names(split_df), function(nm) {
      calc_curves_one(split_df[[nm]], nm)
    })

    curves <- dplyr::bind_rows(curves)

    shiny::validate(
      shiny::need(
        nrow(curves) > 0,
        "No valid curves available for this selection."
      )
    )

    curves
  })

  auc_trapezoid <- function(x, y) {
    ok <- is.finite(x) & is.finite(y)
    x <- x[ok]
    y <- y[ok]

    if (length(x) < 2) {
      return(NA_real_)
    }

    ord <- order(x, y)
    x <- x[ord]
    y <- y[ord]

    sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  }

  auc_table <- shiny::reactive({
    rocpr_curves() |>
      dplyr::group_by(curve_id, curvetype) |>
      dplyr::summarise(
        auc = auc_trapezoid(x, y),
        .groups = "drop"
      ) |>
      dplyr::arrange(curvetype, dplyr::desc(auc))
  })

  roc_fdr_point <- shiny::reactive({
    df <- rocpr_data()
    alpha <- if (is.null(input$alpha)) 0.01 else input$alpha

    df |>
      dplyr::mutate(
        methodDE = !is.na(padj) & padj <= alpha,
        truepos = trueDE & methodDE,
        falsepos = !trueDE & methodDE
      ) |>
      dplyr::group_by(curve_id) |>
      dplyr::summarise(
        n_pos = sum(trueDE, na.rm = TRUE),
        n_neg = sum(!trueDE, na.rm = TRUE),
        tp = sum(truepos, na.rm = TRUE),
        fp = sum(falsepos, na.rm = TRUE),
        n_called = tp + fp,
        observed_fdr = dplyr::if_else(n_called > 0, fp / n_called, NA_real_),
        TPR = dplyr::if_else(n_pos > 0, tp / n_pos, NA_real_),
        FPR = dplyr::if_else(n_neg > 0, fp / n_neg, NA_real_),
        .groups = "drop"
      ) |>
      dplyr::filter(!is.na(TPR), !is.na(FPR))
  })

  pr_fdr_point <- shiny::reactive({
    df <- rocpr_data()
    alpha <- if (is.null(input$alpha)) 0.01 else input$alpha

    df |>
      dplyr::mutate(
        methodDE = !is.na(padj) & padj <= alpha,
        truepos = trueDE & methodDE,
        falsepos = !trueDE & methodDE
      ) |>
      dplyr::group_by(curve_id) |>
      dplyr::summarise(
        n_pos = sum(trueDE, na.rm = TRUE),
        tp = sum(truepos, na.rm = TRUE),
        fp = sum(falsepos, na.rm = TRUE),
        n_called = tp + fp,
        observed_fdr = dplyr::if_else(n_called > 0, fp / n_called, NA_real_),
        recall = dplyr::if_else(n_pos > 0, tp / n_pos, NA_real_),
        precision = dplyr::if_else(n_called > 0, tp / n_called, NA_real_),
        .groups = "drop"
      ) |>
      dplyr::filter(!is.na(recall), !is.na(precision))
  })

  real_fdr_points <- shiny::reactive({
    df <- rocpr_data()
    alpha <- if (is.null(input$alpha)) 0.01 else input$alpha

    split_df <- split(df, df$curve_id)

    out <- lapply(names(split_df), function(nm) {
      df_one <- split_df[[nm]] |>
        dplyr::arrange(dplyr::desc(score))

      n_pos <- sum(df_one$trueDE, na.rm = TRUE)
      n_neg <- sum(!df_one$trueDE, na.rm = TRUE)

      if (n_pos == 0 || n_neg == 0) {
        return(NULL)
      }

      df_one <- df_one |>
        dplyr::mutate(
          tp = cumsum(trueDE),
          fp = cumsum(!trueDE),
          n_called = tp + fp,
          real_fdr = dplyr::if_else(
            n_called > 0,
            fp / n_called,
            NA_real_
          ),
          TPR = tp / n_pos,
          FPR = fp / n_neg,
          recall = TPR,
          precision = dplyr::if_else(
            n_called > 0,
            tp / n_called,
            NA_real_
          )
        ) |>
        dplyr::filter(
          n_called > 0,
          is.finite(real_fdr)
        )

      if (nrow(df_one) == 0) {
        return(NULL)
      }

      below <- df_one |>
        dplyr::filter(real_fdr <= alpha)

      if (nrow(below) > 0) {
        chosen <- below |>
          dplyr::slice_max(
            order_by = n_called,
            n = 1,
            with_ties = FALSE
          )
      } else {
        chosen <- df_one |>
          dplyr::mutate(distance = abs(real_fdr - alpha)) |>
          dplyr::slice_min(
            order_by = distance,
            n = 1,
            with_ties = FALSE
          )
      }

      chosen |>
        dplyr::transmute(
          curve_id = nm,
          FPR = FPR,
          TPR = TPR,
          recall = recall,
          precision = precision,
          real_fdr = real_fdr,
          n_called = n_called,
          tp = tp,
          fp = fp
        )
    })

    out <- dplyr::bind_rows(out)

    if (nrow(out) == 0) {
      return(empty_point_tbl())
    }

    out
  })

  fdr_points_table <- shiny::reactive({
    alpha <- if (is.null(input$alpha)) 0.01 else input$alpha

    nominal_roc <- roc_fdr_point()

    nominal_pr <- pr_fdr_point() |>
      dplyr::select(
        curve_id,
        recall,
        precision
      )

    nominal <- nominal_roc |>
      dplyr::left_join(
        nominal_pr,
        by = "curve_id"
      ) |>
      dplyr::mutate(
        point_type = paste0("Nominal padj <= ", alpha),
        empirical_fdr = observed_fdr
      ) |>
      dplyr::select(
        point_type,
        curve_id,
        n_called,
        tp,
        fp,
        empirical_fdr,
        FPR,
        TPR,
        recall,
        precision
      )

    real <- real_fdr_points() |>
      dplyr::mutate(
        point_type = paste0("Empirical FDR <= ", alpha),
        empirical_fdr = real_fdr
      ) |>
      dplyr::select(
        point_type,
        curve_id,
        n_called,
        tp,
        fp,
        empirical_fdr,
        FPR,
        TPR,
        recall,
        precision
      )

    dplyr::bind_rows(nominal, real) |>
      dplyr::mutate(
        dplyr::across(
          c(empirical_fdr, FPR, TPR, recall, precision),
          ~ round(.x, 5)
        )
      ) |>
      dplyr::arrange(curve_id, point_type)
  })

  confusion_table <- shiny::reactive({
    df <- rocpr_data()
    alpha <- if (is.null(input$alpha)) 0.01 else input$alpha

    # Choose simulation column safely
    sim_col <- intersect(
      c(dataset_col, "analysis_join", "analysis", "sim"),
      names(df)
    )

    if (length(sim_col) > 0) {
      df$simulation <- as.character(df[[sim_col[1]]])
    } else {
      df$simulation <- as.character(df$curve_id)
    }

    df <- df |>
      dplyr::mutate(
        truth = as.logical(trueDE),
        called = !is.na(padj) & padj <= alpha
      )

    group_cols <- c(
      "simulation",
      "scenario",
      "nSubjsPerGroup",
      "seed",
      "method",
      "formula",
      "norm"
    )

    group_cols <- intersect(group_cols, names(df))

    df |>
      dplyr::group_by(dplyr::across(dplyr::all_of(group_cols))) |>
      dplyr::summarise(
        alpha = alpha,
        total_genes = dplyr::n(),

        TP = sum(truth & called, na.rm = TRUE),
        TN = sum(!truth & !called, na.rm = TRUE),
        FP = sum(!truth & called, na.rm = TRUE),
        FN = sum(truth & !called, na.rm = TRUE),

        .groups = "drop"
      ) |>
      dplyr::arrange(
        simulation,
        method,
        formula,
        norm
      )
  })

  roc_plot_obj <- shiny::reactive({
    curves <- rocpr_curves() |>
      dplyr::filter(curvetype == "ROC")

    pt <- roc_fdr_point()
    real_pt <- real_fdr_points()

    roc_xmax <- if (
      is.null(input$roc_xmax) ||
      !is.finite(input$roc_xmax) ||
      input$roc_xmax <= 0
    ) {
      0.05
    } else {
      input$roc_xmax
    }

    roc_x_lims <- if (isTRUE(input$zoom_fdr_region)) {
      c(0, roc_xmax)
    } else {
      c(0, 1)
    }

    p <- ggplot2::ggplot(
      curves,
      ggplot2::aes(x = x, y = y, colour = curve_id, group = curve_id)
    ) +
      ggplot2::geom_step(linewidth = 0.9, alpha = 0.9) +
      ggplot2::geom_abline(
        slope = 1,
        intercept = 0,
        linetype = "dashed",
        colour = "grey50"
      ) +
      ggplot2::scale_colour_discrete(
        labels = function(x) wrap_labels(x, width = 28)
      ) +
      ggplot2::labs(
        x = "False positive rate",
        y = "True positive rate",
        colour = NULL
      ) +
      ggplot2::coord_cartesian(
        xlim = roc_x_lims,
        ylim = c(0, 1),
        expand = FALSE
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18),
        axis.title.x = ggplot2::element_text(size = 18),
        axis.title.y = ggplot2::element_text(size = 18),
        axis.text.x = ggplot2::element_text(size = 13),
        axis.text.y = ggplot2::element_text(size = 13),
        legend.position = if (isTRUE(input$show_legend)) "right" else "none",
        legend.text = ggplot2::element_text(size = 9),
        legend.key.width = grid::unit(0.7, "cm"),
        legend.key.height = grid::unit(0.45, "cm"),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          ncol = 1,
          byrow = TRUE
        )
      )

    if (nrow(pt) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = pt,
          ggplot2::aes(x = FPR, y = TPR, colour = curve_id),
          inherit.aes = FALSE,
          size = 2.8
        )
    }

    if (nrow(real_pt) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = real_pt,
          ggplot2::aes(x = FPR, y = TPR, colour = curve_id),
          inherit.aes = FALSE,
          shape = 4,
          size = 4,
          stroke = 1.3
        )
    }

    p
  })

  pr_plot_obj <- shiny::reactive({
    curves <- rocpr_curves() |>
      dplyr::filter(curvetype == "PRC")

    pt <- pr_fdr_point()
    real_pt <- real_fdr_points()

    pr_xmax <- if (
      is.null(input$pr_xmax) ||
      !is.finite(input$pr_xmax) ||
      input$pr_xmax <= 0
    ) {
      0.20
    } else {
      input$pr_xmax
    }

    pr_x_lims <- if (isTRUE(input$zoom_fdr_region)) {
      c(0, pr_xmax)
    } else {
      c(0, 1)
    }

    p <- ggplot2::ggplot(
      curves,
      ggplot2::aes(x = x, y = y, colour = curve_id, group = curve_id)
    ) +
      ggplot2::geom_step(linewidth = 0.9, alpha = 0.9) +
      ggplot2::scale_colour_discrete(
        labels = function(x) wrap_labels(x, width = 28)
      ) +
      ggplot2::labs(
        x = "Recall",
        y = "Precision",
        colour = NULL
      ) +
      ggplot2::coord_cartesian(
        xlim = pr_x_lims,
        ylim = c(0, 1),
        expand = FALSE
      ) +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18),
        axis.title.x = ggplot2::element_text(size = 18),
        axis.title.y = ggplot2::element_text(size = 18),
        axis.text.x = ggplot2::element_text(size = 13),
        axis.text.y = ggplot2::element_text(size = 13),
        legend.position = if (isTRUE(input$show_legend)) "right" else "none",
        legend.text = ggplot2::element_text(size = 9),
        legend.key.width = grid::unit(0.7, "cm"),
        legend.key.height = grid::unit(0.45, "cm"),
        plot.margin = ggplot2::margin(10, 10, 10, 10)
      ) +
      ggplot2::guides(
        colour = ggplot2::guide_legend(
          ncol = 1,
          byrow = TRUE
        )
      )

    if (nrow(pt) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = pt,
          ggplot2::aes(x = recall, y = precision, colour = curve_id),
          inherit.aes = FALSE,
          size = 2.8
        )
    }

    if (nrow(real_pt) > 0) {
      p <- p +
        ggplot2::geom_point(
          data = real_pt,
          ggplot2::aes(x = recall, y = precision, colour = curve_id),
          inherit.aes = FALSE,
          shape = 4,
          size = 4,
          stroke = 1.3
        )
    }

    p
  })

  output$roc_plot <- shiny::renderPlot({
    roc_plot_obj()
  })

  output$pr_plot <- shiny::renderPlot({
    pr_plot_obj()
  })

  output$roc_table <- DT::renderDT({
    DT::datatable(
      auc_table() |>
        dplyr::filter(curvetype == "ROC"),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })

  output$pr_table <- DT::renderDT({
    DT::datatable(
      auc_table() |>
        dplyr::filter(curvetype == "PRC"),
      options = list(
        pageLength = 10,
        scrollX = TRUE
      )
    )
  })

  output$fdr_points_table <- DT::renderDT({
    DT::datatable(
      fdr_points_table(),
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  output$confusion_table <- DT::renderDT({
    DT::datatable(
      confusion_table(),
      options = list(
        pageLength = 20,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  plots_dir <- file.path(
    path.expand("~"),
    "Desktop",
    "plots"
  )

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }

  current_label <- shiny::reactive({
    if (is.null(input$compare_mode) || input$compare_mode == "methods") {
      paste0("methods_", selected_one_sim_label())
    } else {
      paste0("simulations_", input$method)
    }
  })

  output$download_roc <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "ROC_",
        safe_filename(current_label()),
        "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0(
          "ROC_",
          safe_filename(current_label()),
          "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = roc_plot_obj(),
        width = 8,
        height = 6,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = roc_plot_obj(),
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  )

  output$download_pr <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "PR_",
        safe_filename(current_label()),
        "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0(
          "PR_",
          safe_filename(current_label()),
          "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = pr_plot_obj(),
        width = 8,
        height = 6,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = pr_plot_obj(),
        width = 8,
        height = 6,
        dpi = 300
      )
    }
  )
}
