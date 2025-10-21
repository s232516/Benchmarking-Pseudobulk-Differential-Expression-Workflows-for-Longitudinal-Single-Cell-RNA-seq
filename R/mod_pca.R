#' mod_pca UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_pca_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotlyOutput(ns('pcaPlot'), height = '600px')
  )
}

#' mod_pca Server Functions
#'
#' @noRd
mod_pca_server <- function(input, output, session, shared, colorBy) {
  output$pcaPlot <- plotly::renderPlotly({
    sim <- shared$sim
    pca <- scater::reducedDim(sim, 'PCA')
    df <- as.data.frame(pca[, 1:3])
    df$Sample <- SummarizedExperiment::colData(sim)$Sample
    df$Condition <- SummarizedExperiment::colData(sim)$Condition
    colToUse <- ifelse(colorBy() == 'Sample', df$Sample, df$Condition)

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = PC1,
      y = PC2,
      color = colToUse,
      text = paste0('Sample: ', Sample, '<br>Condition: ', Condition)
    )) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggplot2::labs(color = colorBy()) +
      ggplot2::theme_minimal()

    plotly::ggplotly(p, tooltip = 'text')
  })
}

