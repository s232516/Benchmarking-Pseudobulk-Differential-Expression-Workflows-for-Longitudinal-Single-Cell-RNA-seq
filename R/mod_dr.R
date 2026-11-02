#' mod_dr UI Function
#'
#' @description A shiny Module for dimensionality reduction.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
mod_dr_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::card(
    bslib::layout_sidebar(
      fillable = TRUE,

      sidebar = bslib::accordion(
        id = ns("dr_accordion"),

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
                trigger = shiny::tagList(
                  "Select dataset",
                  shiny::icon("circle-info")
                ),
                "Select the original simulation dataset."
              ),
              choices = character(0)
            )
          ),

          shiny::conditionalPanel(
            condition = "input.sim_source == 'fair'",
            ns = ns,

            shiny::selectInput(
              inputId = ns("scenario"),
              label = "Scenario:",
              choices = character(0)
            ),

            shiny::selectInput(
              inputId = ns("n_subj"),
              label = "Subjects per group:",
              choices = character(0)
            ),

            shiny::selectInput(
              inputId = ns("seed"),
              label = "Seed:",
              choices = character(0)
            )
          ),

          shiny::selectInput(
            inputId = ns("colour_by"),
            label = "Colour by:",
            choices = character(0)
          )
        )
      ),

      bslib::layout_columns(
        col_widths = 12,

        bslib::navset_card_tab(

          bslib::nav_panel(
            title = "PCA",
            icon = shiny::icon("chart-simple"),

            bslib::card_header(
              div(
                style = "display: flex; justify-content: flex-end;",
                shiny::downloadButton(
                  outputId = ns("pca_download_plot"),
                  label = "Download",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("pca_plot"),
                height = "550px",
                width = "800px"
              )
            )
          ),

          bslib::nav_panel(
            title = "UMAP",
            icon = shiny::icon("layer-group"),

            bslib::card_header(
              div(
                style = "display: flex; justify-content: flex-end;",
                shiny::downloadButton(
                  outputId = ns("umap_download_plot"),
                  label = "Download",
                  class = "clickable-button-light"
                )
              )
            ),

            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("umap_plot"),
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
#' mod_dr Server Function
#'
#' @noRd
mod_dr_server <- function(input, output, session, data) {

  # ---------------------------------------------------------------------------
  # Helpers
  # ---------------------------------------------------------------------------

  safe_filename <- function(x) {
    x |>
      gsub("[^A-Za-z0-9_\\-]+", "_", x = _) |>
      gsub("_+", "_", x = _)
  }

  parse_sim_metadata <- function(sim_names) {
    out <- lapply(sim_names, function(nm) {
      is_fair <- grepl("^fair_", nm)

      scenario <- NA_character_
      nSubjsPerGroup <- NA_real_
      seed <- NA_real_

      if (is_fair) {
        m <- regexec("^fair_(.+)_n([0-9]+)_seed([0-9]+)$", nm)
        r <- regmatches(nm, m)[[1]]

        if (length(r) == 4) {
          scenario <- r[2]
          nSubjsPerGroup <- as.numeric(r[3])
          seed <- as.numeric(r[4])
        }
      }

      data.frame(
        sim_name = nm,
        sim_source = ifelse(is_fair, "fair", "old"),
        scenario = scenario,
        nSubjsPerGroup = nSubjsPerGroup,
        seed = seed,
        stringsAsFactors = FALSE
      )
    })

    dplyr::bind_rows(out)
  }

  scenario_order <- c(
    "global_null",
    "time_only_null",
    "weak",
    "medium",
    "strong"
  )

  # ---------------------------------------------------------------------------
  # Simulation metadata
  # ---------------------------------------------------------------------------

  sim_meta <- shiny::reactive({
    parse_sim_metadata(names(data))
  })

  # ---------------------------------------------------------------------------
  # Update old/fair selectors
  # ---------------------------------------------------------------------------

  observe({
    meta <- sim_meta()

    old_datasets <- meta |>
      dplyr::filter(sim_source == "old") |>
      dplyr::pull(sim_name) |>
      sort()

    fair_meta <- meta |>
      dplyr::filter(sim_source == "fair")

    scenarios <- fair_meta |>
      dplyr::filter(!is.na(scenario)) |>
      dplyr::pull(scenario) |>
      unique()

    scenarios <- c(
      intersect(scenario_order, scenarios),
      setdiff(sort(scenarios), scenario_order)
    )

    shiny::updateSelectInput(
      session = session,
      inputId = "data",
      choices = old_datasets,
      selected = if (length(old_datasets) > 0) old_datasets[1] else NULL
    )

    shiny::updateSelectInput(
      session = session,
      inputId = "scenario",
      choices = scenarios,
      selected = if (length(scenarios) > 0) scenarios[1] else NULL
    )
  })

  observeEvent(
    input$scenario,
    {
      req(input$scenario)

      meta <- sim_meta() |>
        dplyr::filter(
          sim_source == "fair",
          scenario == input$scenario
        )

      n_subj <- meta |>
        dplyr::filter(!is.na(nSubjsPerGroup)) |>
        dplyr::pull(nSubjsPerGroup) |>
        unique() |>
        sort()

      shiny::updateSelectInput(
        session = session,
        inputId = "n_subj",
        choices = as.character(n_subj),
        selected = if (length(n_subj) > 0) as.character(n_subj[1]) else NULL
      )
    },
    ignoreInit = FALSE
  )

  observeEvent(
    {
      input$scenario
      input$n_subj
    },
    {
      req(input$scenario, input$n_subj)

      meta <- sim_meta() |>
        dplyr::filter(
          sim_source == "fair",
          scenario == input$scenario,
          nSubjsPerGroup == as.numeric(input$n_subj)
        )

      seeds <- meta |>
        dplyr::filter(!is.na(seed)) |>
        dplyr::pull(seed) |>
        unique() |>
        sort()

      shiny::updateSelectInput(
        session = session,
        inputId = "seed",
        choices = as.character(seeds),
        selected = if (length(seeds) > 0) as.character(seeds[1]) else NULL
      )
    },
    ignoreInit = FALSE
  )

  # ---------------------------------------------------------------------------
  # Selected simulation
  # ---------------------------------------------------------------------------

  selected_simulation_name <- shiny::reactive({
    meta <- sim_meta()

    if (is.null(input$sim_source) || input$sim_source == "fair") {
      req(input$scenario, input$n_subj, input$seed)

      selected <- meta |>
        dplyr::filter(
          sim_source == "fair",
          scenario == input$scenario,
          nSubjsPerGroup == as.numeric(input$n_subj),
          seed == as.numeric(input$seed)
        ) |>
        dplyr::pull(sim_name)

      validate(
        shiny::need(length(selected) > 0, "No fair simulation found for this selection.")
      )

      selected[1]
    } else {
      req(input$data)
      input$data
    }
  })

  selected_sce <- shiny::reactive({
    nm <- selected_simulation_name()

    validate(
      shiny::need(nm %in% names(data), paste("Simulation not found:", nm))
    )

    data[[nm]]
  })

  # ---------------------------------------------------------------------------
  # Prepare SCE for plotting
  # ---------------------------------------------------------------------------

  plot_sce <- shiny::reactive({
    sce <- selected_sce()

    # Add logcounts if missing
    if (!"logcounts" %in% SummarizedExperiment::assayNames(sce)) {
      if ("counts" %in% SummarizedExperiment::assayNames(sce)) {
        sce <- scater::logNormCounts(sce)
      }
    }

    # Add PCA if missing
    if (!"PCA" %in% SingleCellExperiment::reducedDimNames(sce)) {
      sce <- scater::runPCA(
        sce,
        ncomponents = 10
      )
    }

    # Add UMAP if missing
    if (!"UMAP" %in% SingleCellExperiment::reducedDimNames(sce)) {
      sce <- scater::runUMAP(
        sce,
        dimred = "PCA"
      )
    }

    sce
  })

  # ---------------------------------------------------------------------------
  # Update colour_by choices
  # ---------------------------------------------------------------------------

  observeEvent(
    selected_simulation_name(),
    {
      sce <- selected_sce()

      cd_cols <- colnames(
        as.data.frame(SummarizedExperiment::colData(sce))
      )

      preferred <- c(
        "Condition",
        "group",
        "time",
        "Sample",
        "sampleID",
        "subjectID",
        "CellID"
      )

      choices <- c(
        intersect(preferred, cd_cols),
        setdiff(cd_cols, preferred)
      )

      selected <- dplyr::case_when(
        "Condition" %in% choices ~ "Condition",
        "group" %in% choices ~ "group",
        "time" %in% choices ~ "time",
        length(choices) > 0 ~ choices[1],
        TRUE ~ NA_character_
      )

      shiny::updateSelectInput(
        session = session,
        inputId = "colour_by",
        choices = choices,
        selected = selected
      )
    },
    ignoreInit = FALSE
  )

  # ---------------------------------------------------------------------------
  # Plot objects
  # ---------------------------------------------------------------------------

  pca_plot_obj <- shiny::reactive({
    req(input$colour_by)

    sce <- plot_sce()

    validate(
      shiny::need(
        input$colour_by %in% colnames(as.data.frame(SummarizedExperiment::colData(sce))),
        "Selected colour variable is not available in this simulation."
      )
    )

    scater::plotPCA(
      sce,
      colour_by = input$colour_by,
      point_size = 3
    ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      ) +
      ggplot2::labs(
        title = paste("PCA Plot -", selected_simulation_name())
      )
  })

  umap_plot_obj <- shiny::reactive({
    req(input$colour_by)

    sce <- plot_sce()

    validate(
      shiny::need(
        input$colour_by %in% colnames(as.data.frame(SummarizedExperiment::colData(sce))),
        "Selected colour variable is not available in this simulation."
      )
    )

    scater::plotUMAP(
      sce,
      colour_by = input$colour_by,
      point_size = 3
    ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      ) +
      ggplot2::labs(
        title = paste("UMAP Plot -", selected_simulation_name())
      )
  })

  # ---------------------------------------------------------------------------
  # Render plots
  # ---------------------------------------------------------------------------

  output$pca_plot <- shiny::renderPlot({
    pca_plot_obj()
  })

  output$umap_plot <- shiny::renderPlot({
    umap_plot_obj()
  })

  # ---------------------------------------------------------------------------
  # Downloads
  # ---------------------------------------------------------------------------

  plots_dir <- file.path(
    path.expand("~"),
    "Desktop",
    "plots"
  )

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }

  output$pca_download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "PCA_",
        safe_filename(selected_simulation_name()),
        "_",
        input$colour_by,
        "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0(
          "PCA_",
          safe_filename(selected_simulation_name()),
          "_",
          input$colour_by,
          "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = pca_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = pca_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )
    }
  )

  output$umap_download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0(
        "UMAP_",
        safe_filename(selected_simulation_name()),
        "_",
        input$colour_by,
        "_",
        Sys.Date(),
        ".png"
      )
    },
    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0(
          "UMAP_",
          safe_filename(selected_simulation_name()),
          "_",
          input$colour_by,
          "_",
          Sys.Date(),
          ".png"
        )
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = umap_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = umap_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )
    }
  )
}
