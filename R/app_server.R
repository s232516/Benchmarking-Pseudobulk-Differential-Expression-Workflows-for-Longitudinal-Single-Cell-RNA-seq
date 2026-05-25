#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # load fixed data

  board <- pins::board_local("nrgn-data")
  simulation <- pins::pin_read(board, "simulations_final")
  compare <- pins::pin_read(board, "compare_data_paper")
  # for normalized counts comment this line.
  compare <- compare |>
    dplyr::filter(
      .data$norm == TRUE |
        .data$method %in% c("DESeq2 pseudobulk", "edgeR pseudobulk",
                            "limma-voom pseudobulk")
    )
  ci_data <- pins::pin_read(board, "rescueSim_ci_paper")
  ground_truth <- pins::pin_read(board, "ground_truth_final")

  # callModule(mod_overview_server, "overview", sim = sim_pop_cond)
  callModule(mod_dr_server, "dim_red", data = simulation)
  callModule(mod_volcano_server, "volcano", data = compare) # done
  callModule(mod_roc_server, "roc", data = compare) # norm
  callModule(mod_gt_server, "gt", data = compare)
  callModule(mod_hist_logfc_server, "hist_module", data = compare)
  callModule(mod_ci_hist_server, "ci_module", data = ci_data)
  callModule(mod_goi_server, "goi", sce_list = simulation) # done
  callModule(mod_logfc_count_server, "lc_module", data = compare,
             sim = simulation)
  callModule(mod_counts_hist_server, "ch_module", data = simulation)# done
  callModule(mod_filter_summary_server, "fs_module", compare_data = compare,
  ground_truth = ground_truth)
  callModule(mod_es_compare_server, "fc", data = compare)
  # callModule(mod_de_heatmap_server, "de_heatmap", data = compare)
}
