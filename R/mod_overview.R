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
    verbatimTextOutput(ns('summary')),
    downloadButton(ns('download_counts'), 'Download normalized counts (CSV)')
  )
}

#' mod_overview Server Functions
#'
#' @noRd
mod_overview_server <- function(input, output, session, shared) {
  ns <- session$ns
  output$summary <- renderPrint({
    sim <- shared$sim
    cat('Number of genes:', nrow(sim), '
')
    cat('Number of cells:', ncol(sim), '
')
    cat('Condition table:
')
    print(table(colData(sim)$Condition))
    counts_mat <- counts(sim)
    cat(sprintf('Proportion zeros in counts: %.3f
', mean(counts_mat == 0)))
  })
  output$download_counts <- downloadHandler(
    filename = function() 'normalized_counts_fixed.csv',
    content = function(file) {
      sim <- shared$sim
      nm <- assay(sim, 'logcounts')
      write.csv(as.data.frame(nm), file)
    }
  )
}
