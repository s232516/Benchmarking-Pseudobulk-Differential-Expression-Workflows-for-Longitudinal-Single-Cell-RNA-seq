#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  tagList(
    golem_add_external_resources(),
    shinyjs::useShinyjs(),
    bslib::page_navbar(
      title = 'Single Cell Methods',
      id = "sc",
      bg = "#001965",
      inverse = "auto",
      theme = bslib::bs_theme(
        bg = "#FFFFFF",
        fg = "#001965",
        primary = "#99a3c1"
      ) |>
        bslib::bs_add_rules(
          list()
        ),
      fillable = TRUE,
      # bslib::nav_panel('Overview', mod_overview_ui('overview')),
      bslib::nav_panel('DR Plots', mod_dr_ui('dim_red')),
      bslib::nav_panel('Volcano Plots', mod_volcano_ui('volcano')),
      bslib::nav_panel('ROC Plots', mod_roc_ui('roc')),
      bslib::nav_panel('TPR lineplot', mod_tpr_line_ui("tpr_line")),
      # bslib::nav_panel('Ground Truth Plot', mod_gt_ui('gt')),
      # bslib::nav_panel('Gene of Interest', mod_goi_ui('goi')),
      # bslib::nav_panel('Histogram', mod_hist_logfc_ui('hist_module')),
      bslib::nav_panel('CI Coverage Histogram', mod_ci_hist_ui('ci_module')),
      bslib::nav_panel('CI summary plot',
                       mod_ci_summary_ui("ci_summary_module"))

      # bslib::nav_panel('Log2 Fold Change', mod_logfc_count_ui('lc_module')),
      # bslib::nav_panel('Counts Histogram', mod_counts_hist_ui('ch_module')),
      # bslib::nav_panel('Filter Summary', mod_filter_summary_ui('fs_module')),
      # bslib::nav_panel('LFC',  mod_es_compare_ui('fc'))
      # # bslib::nav_panel('Heatmap', mod_de_heatmap_ui("de_heatmap"))
    )
  )
}

#' Add external Resources to the Application
#'
#' This function is internally used to add external
#' resources inside the Shiny application.
#'
#' @import shiny
#' @importFrom golem add_resource_path activate_js favicon bundle_resources
#' @noRd
golem_add_external_resources <- function() {
  add_resource_path(
    "www",
    app_sys("app/www")
  )

  tags$head(
    # Link your CSS
    tags$link(rel = "stylesheet", type = "text/css", href = "www/custom.css")
  )
}
