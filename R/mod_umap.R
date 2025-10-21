#' mod_umap UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
library(shiny)
library(ggplot2)
library(plotly)


mod_umap_ui <- function(id) {
  ns <- NS(id)
  tagList(
    plotlyOutput(ns('umapPlot'), height = '600px')
  )
}

#' mod_umap Server Functions
#'
#' @noRd
mod_umap_server <- function(input, output, session, shared, colorBy) {
  output$umapPlot <- plotly::renderPlotly({
    sim <- shared$sim
    umap <- reducedDim(sim, 'UMAP')
    df <- as.data.frame(umap[,1:2])
    colnames(df) <- c('UMAP1','UMAP2')
    df$Sample <- colData(sim)$Sample
    df$Condition <- colData(sim)$Condition
    colToUse <- ifelse(colorBy() == 'Sample', df$Sample, df$Condition)
    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = UMAP1,
      y = UMAP2,
      color = colToUse,
      text = paste0('Sample: ', Sample, '<br>Condition: ', Condition)
    )) +
      ggplot2::geom_point(alpha = 0.8, size = 2) +
      ggplot2::labs(color = colorBy()) +
      ggplot2::theme_minimal()

    plotly::ggplotly(p, tooltip = 'text')

  })
}
