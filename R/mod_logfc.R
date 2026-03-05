mod_es_compare_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 340,
      open = "always",

      sidebar = bslib::accordion(
        id = ns("tx_accordion"),

        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),

          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select dataset", shiny::icon("circle-info")),
              "Select which dataset/analysis to compare within."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("norm"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Normalization", shiny::icon("circle-info")),
              "Filter results by whether log-normalization was used (only applies when the 'norm' column exists)."
            ),
            choices = c("All" = "all"),
            selected = "all"
          ),

          shiny::hr(),

          shiny::selectInput(
            inputId = ns("x_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList("X-axis method", shiny::icon("circle-info")),
              "Choose the method for X-axis effect sizes."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("x_formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList("X-axis formula", shiny::icon("circle-info")),
              "Choose the formula for the selected X-axis method."
            ),
            choices = NULL
          ),

          shiny::hr(),

          shiny::selectInput(
            inputId = ns("y_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Y-axis method", shiny::icon("circle-info")),
              "Choose the method for Y-axis effect sizes."
            ),
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("y_formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Y-axis formula", shiny::icon("circle-info")),
              "Choose the formula for the selected Y-axis method."
            ),
            choices = NULL
          ),

          shiny::hr(),

          shiny::checkboxInput(
            inputId = ns("highlight_truede"),
            label = "Highlight true DE genes (yellow)",
            value = TRUE
          ),

          shiny::sliderInput(
            inputId = ns("logfc_thr"),
            label = "Absolute logFC threshold (optional filter)",
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
                "Effect Size Scatter",
                shiny::downloadButton(
                  outputId = ns("download_plot"),
                  label = "",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("scatter_plot"),
                height = "550px",
                width = "100%"
              )
            )
          ),

          bslib::card(
            bslib::card_header("Scatter Plot Data Table"),
            bslib::card_body(
              DT::DTOutput(ns("scatter_table"))
            )
          )
        )
      )
    )
  )
}
mod_es_compare_server <- function(input, output, session, data) {

  methods_with_true <- function(methods_vec) {
    unique(c(sort(methods_vec), "True logFC"))
  }

  has_formula_for_method <- function(df, method_name) {
    if (!"formula" %in% names(df)) return(FALSE)
    if (is.null(method_name) || is.na(method_name)) return(FALSE)
    if (!"method" %in% names(df)) return(FALSE)
    if (!method_name %in% df$method) return(FALSE)

    vals <- df$formula[df$method == method_name]
    vals <- vals[!is.na(vals)]
    vals <- as.character(vals)
    vals <- vals[nzchar(vals)]
    length(unique(vals)) > 0
  }

  get_true_logfc_tbl <- function(df) {
    if ("true_logFC" %in% names(df)) {
      return(df |>
               dplyr::select(gene, true_logfc = true_logFC) |>
               dplyr::distinct())
    }

    df |>
      dplyr::select(gene) |>
      dplyr::distinct() |>
      dplyr::mutate(true_logfc = NA_real_)
  }

  base_filtered <- shiny::reactive({
    req(input$data, input$norm)

    df <- data |>
      dplyr::filter(analysis == input$data)

    if (input$norm != "all" && "norm" %in% names(df)) {
      df <- df |>
        dplyr::filter(norm == as.logical(input$norm))
    }
    df
  })

  # ---- dataset selector ----
  observe({
    shiny::updateSelectInput(
      session,
      "data",
      choices = unique(data$analysis),
      selected = unique(data$analysis)[1]
    )
  })

  # ---- norm selector ----
  observeEvent(input$data, {
    df <- data |>
      dplyr::filter(analysis == input$data)

    if ("norm" %in% names(df) && any(!is.na(df$norm))) {
      avail <- sort(unique(df$norm[!is.na(df$norm)]))
      choices <- c("All" = "all")
      if (TRUE %in% avail)  choices <- c(choices, "Normalized" = "TRUE")
      if (FALSE %in% avail) choices <- c(choices, "Not normalized" = "FALSE")

      shiny::updateSelectInput(session, "norm", choices = choices, selected = "all")
    } else {
      shiny::updateSelectInput(session, "norm", choices = c("All" = "all"), selected = "all")
    }
  })

  # ---- methods ----
  observeEvent(c(input$data, input$norm), {
    df <- base_filtered()
    method_choices <- methods_with_true(unique(df$method))

    shiny::updateSelectInput(session, "x_method", choices = method_choices, selected = method_choices[1])
    shiny::updateSelectInput(
      session, "y_method",
      choices = method_choices,
      selected = if (length(method_choices) > 1) method_choices[2] else method_choices[1]
    )
  })

  # ---- formulas (axis-specific) ----
  observeEvent(c(input$data, input$norm, input$x_method), {
    req(input$x_method)
    df <- base_filtered()

    if (identical(input$x_method, "True logFC") || !has_formula_for_method(df, input$x_method)) {
      shiny::updateSelectInput(session, "x_formula", choices = c("N/A" = "na"), selected = "na")
    } else {
      choices <- sort(unique(df$formula[df$method == input$x_method]))
      choices <- choices[!is.na(choices)]
      choices <- as.character(choices)
      choices <- choices[nzchar(choices)]
      shiny::updateSelectInput(session, "x_formula", choices = choices, selected = choices[1])
    }
  })

  observeEvent(c(input$data, input$norm, input$y_method), {
    req(input$y_method)
    df <- base_filtered()

    if (identical(input$y_method, "True logFC") || !has_formula_for_method(df, input$y_method)) {
      shiny::updateSelectInput(session, "y_formula", choices = c("N/A" = "na"), selected = "na")
    } else {
      choices <- sort(unique(df$formula[df$method == input$y_method]))
      choices <- choices[!is.na(choices)]
      choices <- as.character(choices)
      choices <- choices[nzchar(choices)]
      shiny::updateSelectInput(session, "y_formula", choices = choices, selected = choices[1])
    }
  })

  # ---- logFC range ----
  observe({
    shiny::updateSliderInput(
      session,
      "logfc_thr",
      min = 0,
      max = round(max(abs(data$log2FoldChange), na.rm = TRUE)),
      value = 0,
      step = 0.01
    )
  })

  # ----------------------------
  # Compare data (FIXED trueDE handling)
  # ----------------------------
  compare_data <- shiny::reactive({
    req(input$x_method, input$y_method)

    df <- base_filtered()

    x_needs_formula <- (!identical(input$x_method, "True logFC")) && has_formula_for_method(df, input$x_method)
    y_needs_formula <- (!identical(input$y_method, "True logFC")) && has_formula_for_method(df, input$y_method)

    if (x_needs_formula) req(input$x_formula)
    if (y_needs_formula) req(input$y_formula)

    # X table
    if (identical(input$x_method, "True logFC")) {
      x_tbl <- get_true_logfc_tbl(df) |>
        dplyr::transmute(gene, logfc_x = true_logfc) |>
        dplyr::distinct()
    } else {
      x_df <- df |>
        dplyr::filter(method == input$x_method)

      if (x_needs_formula) {
        x_df <- x_df |>
          dplyr::filter(formula == input$x_formula)
      }

      x_tbl <- x_df |>
        dplyr::select(
          gene,
          logfc_x = log2FoldChange,
          trueDE  = dplyr::any_of("trueDE")
        ) |>
        dplyr::distinct()
    }

    # Y table
    if (identical(input$y_method, "True logFC")) {
      y_tbl <- get_true_logfc_tbl(df) |>
        dplyr::transmute(gene, logfc_y = true_logfc) |>
        dplyr::distinct()
    } else {
      y_df <- df |>
        dplyr::filter(method == input$y_method)

      if (y_needs_formula) {
        y_df <- y_df |>
          dplyr::filter(formula == input$y_formula)
      }

      y_tbl <- y_df |>
        dplyr::select(
          gene,
          logfc_y = log2FoldChange,
          trueDE  = dplyr::any_of("trueDE")
        ) |>
        dplyr::distinct()
    }

    out <- dplyr::full_join(x_tbl, y_tbl, by = "gene")

    # unify trueDE if present from either side
    if ("trueDE.x" %in% names(out) || "trueDE.y" %in% names(out)) {
      out <- out |>
        dplyr::mutate(trueDE = dplyr::coalesce(.data$`trueDE.x`, .data$`trueDE.y`)) |>
        dplyr::select(-dplyr::any_of(c("trueDE.x", "trueDE.y")))
    }

    # ---- FIX: guarantee trueDE exists ----
    if (!"trueDE" %in% names(out)) {
      out$trueDE <- NA_real_
    }

    # optional abs(logFC) filter
    if (is.numeric(input$logfc_thr) && input$logfc_thr > 0) {
      out <- out |>
        dplyr::filter(
          (is.na(logfc_x) | abs(logfc_x) >= input$logfc_thr) |
            (is.na(logfc_y) | abs(logfc_y) >= input$logfc_thr)
        )
    }

    # ---- Highlight logic ----
    out <- out |>
      dplyr::mutate(
        is_true_de = !is.na(trueDE) & trueDE == 1
      )

    # OPTIONAL BEHAVIOR:
    # If you want "when an axis method is trueDE, treat ALL genes as true", uncomment:
    # if (input$x_method %in% c("trueDE","TrueDE") || input$y_method %in% c("trueDE","TrueDE")) {
    #   out$is_true_de <- TRUE
    # }

    out
  })

  # ----------------------------
  # Plot
  # ----------------------------
  scatter_plot_obj <- shiny::reactive({
    df <- as.data.frame(compare_data())
    validate(shiny::need(nrow(df) > 0, "No data available after filtering."))

    df <- df |>
      dplyr::filter(!is.na(logfc_x), !is.na(logfc_y))

    validate(shiny::need(nrow(df) > 0, "No complete (X,Y) pairs to plot."))

    df <- df |>
      dplyr::mutate(
        color_group = dplyr::if_else(is_true_de & isTRUE(input$highlight_truede), "true DE", "other")
      )

    df_base <- base_filtered()
    x_needs_formula <- (!identical(input$x_method, "True logFC")) && has_formula_for_method(df_base, input$x_method)
    y_needs_formula <- (!identical(input$y_method, "True logFC")) && has_formula_for_method(df_base, input$y_method)

    x_label <- if (identical(input$x_method, "True logFC")) {
      "Effect size (log2FC): True logFC"
    } else {
      paste0("Effect size (log2FC): ", input$x_method)
    }

    y_label <- if (identical(input$y_method, "True logFC")) {
      "Effect size (log2FC): True logFC"
    }  else {
      paste0("Effect size (log2FC): ", input$y_method)
    }

    ggplot2::ggplot(df, ggplot2::aes(x = logfc_x, y = logfc_y)) +
      ggplot2::geom_hline(yintercept = 0, linetype = "dashed", alpha = 0.4) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed", alpha = 0.4) +
      ggplot2::geom_point(ggplot2::aes(color = color_group), alpha = 0.85,
                          size = 2) +
      ggplot2::scale_color_manual(
        values = c("true DE" = "#C9980B", "other" = "#c0c0c0"),
        breaks = c("true DE", "other"),
        labels = c("True DE", "Other")
      ) +
      ggplot2::labs(
        x = x_label,
        y = y_label,
        color = NULL
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, color = "black",
                           linetype = "dashed", size = 1) +
      ggplot2::coord_equal() +
      ggplot2::theme_bw(base_size = 16) +
      ggplot2::theme(
                     aspect.ratio = 1,
                     axis.title.x = ggplot2::element_text(size = 25),
                     axis.title.y = ggplot2::element_text(size = 25),
                     axis.text.x  = ggplot2::element_text(size = 20),
                     axis.text.y  = ggplot2::element_text(size = 20),
                     legend.position = "none"
                     )
  })

  output$scatter_plot <- shiny::renderPlot({
    scatter_plot_obj()
  })

  output$scatter_table <- DT::renderDT({
    DT::datatable(
      compare_data(),
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
      )
    )
  })

  # ----------------------------
  # Download
  # ----------------------------
  plots_dir <- file.path(path.expand("~"), "Desktop", "plots")
  if (!dir.exists(plots_dir)) dir.create(plots_dir, recursive = TRUE)
  make_fname <- function() {
    x_part <- if (identical(input$x_method, "True logFC")) {
      "truth"
    } else if (!is.null(input$x_formula) && input$x_formula != "na") {
      paste0(gsub("\\s+", "_", input$x_method), "_", gsub("\\s+", "_", input$x_formula))
    } else {
      gsub("\\s+", "_", input$x_method)
    }

    y_part <- if (identical(input$y_method, "True logFC")) {
      "truth"
    } else if (!is.null(input$y_formula) && input$y_formula != "na") {
      paste0(gsub("\\s+", "_", input$y_method), "_", gsub("\\s+", "_", input$y_formula))
    } else {
      gsub("\\s+", "_", input$y_method)
    }

    paste0("lfc_", input$data, "_", x_part, "_vs_", y_part, ".png")
  }

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      req(input$data, input$x_method, input$y_method)
      make_fname()
    },
    content = function(file) {
      req(scatter_plot_obj())

      desktop_file <- file.path(plots_dir, make_fname())

      ggplot2::ggsave(desktop_file, plot = scatter_plot_obj(), width = 7, height = 7, dpi = 300)
      ggplot2::ggsave(file,         plot = scatter_plot_obj(), width = 7, height = 7, dpi = 300)
    }
  )
}
