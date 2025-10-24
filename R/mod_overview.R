#' mod_overview UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_overview_ui <- function(id) {
  ns <- NS(id)
  tagList(
    verbatimTextOutput(ns('summary'))
  )
}

#' mod_overview Server Functions
#'
#' @noRd
mod_overview_server <- function(input, output, session, sim) {
  ns <- session$ns
  output$summary <- renderPrint({
    cat('Number of genes:', nrow(sim), '
')
    cat('Number of cells:', ncol(sim), '
')
    cat('Condition table:
')
    print(table(SummarizedExperiment::colData(sim)$Condition))
    counts_mat <- SingleCellExperiment::counts(sim)
    cat(sprintf('Proportion zeros in counts: %.3f
', mean(counts_mat == 0)))
  })

}
