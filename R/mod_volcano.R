#' mod_volcano UI Function
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


mod_volcano_ui <- function(id) {
  ns <- NS(id)
  tagList(
    sliderInput(ns('pCut'), 'Adjusted p-value cutoff', min = 0, max = 0.1, value = 0.05, step = 0.005),
    numericInput(ns('lfcCut'), 'Log2 fold-change cutoff', value = 1, step = 0.1),
    plotlyOutput(ns('volcanoPlot'), height = '600px')
  )
}

#' mod_volcano Server Functions
#'
#' @noRd
mod_volcano_server <- function(input, output, session, shared) {
  ns <- session$ns

  volcano_data <- reactive({
    df <- switch(input$de_method,
                 'MAST' = shared$mast,
                 'DESeq2' = shared$deseq)
    req(df)
    df %>% mutate(significant = padj <= input$pCut & abs(log2FoldChange) >= input$lfcCut)
  })

  output$volcanoPlot <- plotly::renderPlotly({
    df <- volcano_data()
    p <- volcano_plot(df, x = 'log2FoldChange', y = 'padj', labels = 'gene',
                      hline = input$pCut, vline = input$lfcCut)
    ggplotly(p, tooltip = 'text')
  })

  output$volcanoTable <- DT::renderDT({
    df <- volcano_data()
    datatable(df, options = list(pageLength = 10))
  })

  output$download_table <- downloadHandler(
    filename = function() paste0(input$de_method, '_volcano_table.csv'),
    content = function(file) {
      df <- volcano_data()
      write.csv(df, file, row.names = FALSE)
    }
  )
}
