#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # load fixed data

  board <- pins::board_local("nrgn-data")

  method_data  <- pins::pin_read(board, "method_data_final")
  simulation <- pins::pin_read(board, "simulations")
  compare <- pins::pin_read(board, "compare_data_final")
  ci_data <- pins::pin_read(board, "rescueSim_ci_final")
  ground_truth <- pins::pin_read(board, "ground_truth")

  # callModule(mod_overview_server, "overview", sim = sim_pop_cond)
  callModule(mod_dr_server, "dim_red", data = simulation)
  callModule(mod_volcano_server, "volcano", data = compare)
  callModule(mod_roc_server, "roc", data = compare)
  callModule(mod_gt_server, "gt", data = compare)
  callModule(mod_hist_logfc_server, "hist_module", data = compare)
  callModule(mod_ci_hist_server, "ci_module", data = ci_data)
  callModule(mod_goi_server, "goi", sce_list = simulation)
  callModule(mod_logfc_count_server, "lc_module", data = compare,
             sim = simulation)
  callModule(mod_counts_hist_server, "ch_module", data = simulation)
  callModule(mod_filter_summary_server, "fs_module", compare_data = compare,
  ground_truth = ground_truth)
}
