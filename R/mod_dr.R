#' mod_dr UI Function
#'
#' @description A shiny Module for dimensionality reduction.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_dr_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::layout_sidebar(
      fillable = TRUE,
      sidebar = bslib::accordion(
          id = ns("dr_accordion"),
          bslib::accordion_panel(
            title = "Analysis",
            icon = shiny::icon("bars"),
            id = ns("dr_analysis"),
            shiny::selectInput(
              inputId = ns("data"),
              label = bslib::tooltip(
                trigger = shiny::tagList("Select dataset",
                                         shiny::icon("circle-info")),
                "Select which is the dataset that you want to use for
                dimensionality reduction."
              ),
              choices = character(0)
            )
          )
        ),
      bslib::layout_columns(
        col_widths = 12,
        bslib::navset_card_tab(
          id = ns("dr_tabs"),
          bslib::nav_panel(
            title = "PCA",
            id = ns("pca"),
            shiny::plotOutput(ns("pca_plot"), height = "600px")
          ),
          bslib::nav_panel(
            title = "UMAP",
            icon = shiny::icon("layer-group"),
            id = ns("umap"),
            shiny::plotOutput(ns("umap_plot"), height = "600px")
          )
        )
      )
    )
  )
}

#' mod_dr Server Functions
#'
#' @noRd
mod_dr_server <- function(input, output, session, data) {

  output$pca_plot <- shiny::renderPlot({

    scater::plotPCA(
      data,
      colour_by = "Sample",
      shape_by = "Condition",
      point_size = 3
    ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 12)
      ) +
      ggplot2::labs(title = "PCA Plot")
  })

  output$umap_plot <- shiny::renderPlot({

    # Step 4: plot UMAP
    scater::plotUMAP(
      data,
      colour_by = "Condition",
      point_size = 3
    ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        text = ggplot2::element_text(size = 12)
      ) +
      ggplot2::labs(title = "UMAP Plot")
  })
}
