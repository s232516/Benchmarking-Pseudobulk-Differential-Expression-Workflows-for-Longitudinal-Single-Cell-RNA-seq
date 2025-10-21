#' The application server-side
#'
#' @param input,output,session Internal parameters for {shiny}.
#'     DO NOT REMOVE.
#' @import shiny
#' @noRd
app_server <- function(input, output, session) {
  # load fixed data
  data_path <- file.path('data','simData.rda')
  if (!file.exists(data_path)) stop('data/simData.rda not found. Run R/data_prep.R first.')
  load(data_path) # provides sim.pop.cond, mast_res, deseq_res


  # make a small data list for modules
  shared <- reactiveValues(
    sim = sim.pop.cond,
    mast = mast_res,
    deseq = deseq_res
  )


  callModule(mod_overview_server, 'overview_1', shared = shared)
  callModule(mod_pca_server, 'pca_1', shared = shared, colorBy = reactive(input$colorBy))
  callModule(mod_umap_server, 'umap_1', shared = shared, colorBy = reactive(input$colorBy))
  callModule(mod_volcano_server, 'volcano_1', shared = shared)
}


# To run interactively without golem packaging
# if (interactive()) {
#   shiny::shinyApp(ui = app_ui, server = app_server)
# }
