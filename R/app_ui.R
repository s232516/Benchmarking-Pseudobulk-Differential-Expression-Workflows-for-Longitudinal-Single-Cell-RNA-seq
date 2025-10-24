#' The application User-Interface
#'
#' @param request Internal parameter for `{shiny}`.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_ui <- function(request) {
  bslib::page_navbar(
    title = 'Single Cell Methods',
    id = "sc",
    bg = "#001965",
    bslib::nav_panel('Overview', mod_overview_ui('overview')),
    bslib::nav_panel('DR Plots', mod_dr_ui('dim_red')),
    bslib::nav_panel('Volcano Plots', mod_volcano_ui('volcano')),
    bslib::nav_panel('Error Bars', mod_errorbars_ui('error_bars'))
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
    favicon(),
    bundle_resources(
      path = app_sys("app/www"),
      app_title = "MsC"
    )
    # Add here other external resources
    # for example, you can add shinyalert::useShinyalert()
  )
}
