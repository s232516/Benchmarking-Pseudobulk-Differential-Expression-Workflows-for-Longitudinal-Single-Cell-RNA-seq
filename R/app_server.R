#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {

  # ---------------------------------------------------------------------------
  # Pins board
  # ---------------------------------------------------------------------------

  board_path <- file.path(
    "C:/Users/nerea/OneDrive/Documentos/MsC",
    "nrgn-data"
  )

  message("Working directory: ", getwd())
  message("Pins board path: ", normalizePath(board_path, mustWork = FALSE))

  board <- pins::board_local(board_path)

  available_pins <- pins::pin_list(board)

  message("Available pins: ", paste(available_pins, collapse = ", "))

  # ---------------------------------------------------------------------------
  # Helper: safely read optional pins
  # ---------------------------------------------------------------------------

  read_pin_or_null <- function(board, name) {
    message("Trying to read pin: ", name)

    available <- pins::pin_list(board)

    if (!name %in% available) {
      message("Pin is not listed in this board: ", name)
      return(NULL)
    }

    tryCatch(
      {
        x <- pins::pin_read(board, name)

        message(
          "Successfully read pin: ", name,
          " | class: ", paste(class(x), collapse = ", "),
          " | length: ", length(x)
        )

        x
      },
      error = function(e) {
        message("Error while reading pin: ", name)
        message("Actual error message: ", conditionMessage(e))
        NULL
      }
    )
  }

  # ---------------------------------------------------------------------------
  # Helper: combine old and fair objects
  # ---------------------------------------------------------------------------

  combine_old_and_fair <- function(old, fair) {
    if (is.null(fair)) {
      return(old)
    }

    if (is.data.frame(old) && is.data.frame(fair)) {
      return(dplyr::bind_rows(old, fair))
    }

    if (is.list(old) && is.list(fair)) {
      return(c(old, fair))
    }

    c(old, fair)
  }

  # ---------------------------------------------------------------------------
  # Simulations
  # ---------------------------------------------------------------------------

  old_simulation <- read_pin_or_null(board, "simulations_final")
  old_ground_truth <- read_pin_or_null(board, "ground_truth_final")

  fair_simulation <- read_pin_or_null(board, "fair_rescue_simulations")
  fair_ground_truth <- read_pin_or_null(board, "fair_rescue_ground_truth")

  validate(
    shiny::need(!is.null(old_simulation), "Could not load simulations_final."),
    shiny::need(!is.null(old_ground_truth), "Could not load ground_truth_final.")
  )

  simulation <- combine_old_and_fair(
    old = old_simulation,
    fair = fair_simulation
  )

  ground_truth <- combine_old_and_fair(
    old = old_ground_truth,
    fair = fair_ground_truth
  )

  message(
    "Final simulation object | class: ",
    paste(class(simulation), collapse = ", "),
    " | length: ",
    length(simulation)
  )

  message(
    "Final ground truth object | class: ",
    paste(class(ground_truth), collapse = ", "),
    " | length: ",
    length(ground_truth)
  )

  # ---------------------------------------------------------------------------
  # Method comparison data
  # ---------------------------------------------------------------------------

  compare <- read_pin_or_null(board, "compare_data_paper_plus_fair")

  if (is.null(compare)) {
    message("Using old compare_data_paper pin instead.")
    compare <- read_pin_or_null(board, "compare_data_paper")
  }

  validate(
    shiny::need(!is.null(compare), "Could not load compare data.")
  )

  if ("trueDE" %in% colnames(compare)) {
    compare <- compare |>
      dplyr::mutate(
        trueDE = as.logical(trueDE)
      )
  }

  if ("analysis_join" %in% colnames(compare)) {
    compare <- compare |>
      dplyr::mutate(
        analysis_display = analysis_join
      )
  } else if ("analysis" %in% colnames(compare)) {
    compare <- compare |>
      dplyr::mutate(
        analysis_display = analysis
      )
  }

  message(
    "Compare data loaded | rows: ",
    nrow(compare),
    " | columns: ",
    ncol(compare)
  )

  # ---------------------------------------------------------------------------
  # CI data
  # ---------------------------------------------------------------------------

  ci_data <- read_pin_or_null(board, "rescueSim_ci_paper_plus_fair")

  if (is.null(ci_data)) {
    message("Using old rescueSim_ci_paper pin instead.")
    ci_data <- read_pin_or_null(board, "rescueSim_ci_paper")
  }

  validate(
    shiny::need(!is.null(ci_data), "Could not load CI data.")
  )

  message(
    "CI data loaded | class: ",
    paste(class(ci_data), collapse = ", "),
    " | rows: ",
    ifelse(is.data.frame(ci_data), nrow(ci_data), NA),
    " | length: ",
    length(ci_data)
  )

  # ---------------------------------------------------------------------------
  # Modules
  # ---------------------------------------------------------------------------

  callModule(
    module = mod_dr_server,
    id = "dim_red",
    data = simulation
  )

  callModule(
    module = mod_volcano_server,
    id = "volcano",
    data = compare
  )

  callModule(
    module = mod_roc_server,
    id = "roc",
    data = compare
  )

  callModule(
    module = mod_ci_hist_server,
    id = "ci_module",
    data = ci_data
  )

  # callModule(mod_gt_server, "gt", data = compare)

  # callModule(mod_hist_logfc_server, "hist_module", data = compare)

  # callModule(mod_goi_server, "goi", sce_list = simulation)

  # callModule(
  #   mod_logfc_count_server,
  #   "lc_module",
  #   data = compare,
  #   sim = simulation
  # )

  # callModule(mod_counts_hist_server, "ch_module", data = simulation)

  # callModule(
  #   mod_filter_summary_server,
  #   "fs_module",
  #   compare_data = compare,
  #   ground_truth = ground_truth
  # )

  # callModule(mod_es_compare_server, "fc", data = compare)

  # callModule(mod_de_heatmap_server, "de_heatmap", data = compare)
}
