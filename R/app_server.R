#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # load fixed data

  board <- pins::board_local("nrgn-data")

  method_data  <- pins::pin_read(board, "method_data")
  simulation <- pins::pin_read(board, "simulations")
  compare <- pins::pin_read(board, "compare_data")

  # callModule(mod_overview_server, "overview", sim = sim_pop_cond)
  callModule(mod_dr_server, "dim_red", data = simulation)
  callModule(mod_volcano_server, "volcano", data = method_data)
  callModule(mod_roc_server, "roc", data = compare)

  # callModule(mod_errorbars_server, "error_bars",
  #   data = method_data,
  #   sim = sim_pop_cond
  # )
  # callModule(mod_goi_server, "goi", data = method_data)
}
