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
            bslib::nav_panel(
              title = "PCA",
              icon = shiny::icon("chart-simple"),
                bslib::card_header(
                  div(
                    style = "display: flex; justify-content: flex-end;
            align-items: center;",
                    shiny::downloadButton(
                      class = "clickable-button-light",
                      outputId = ns("pca_download_plot"),
                      label = ""
                    )
                  )),
                bslib::card_body(
                  style = "display = flex; justify-content: space-between;
            align-items: center;",
                  div(
                    style = "width: 800px; height:800px; overflow: auto;",
                    shiny::plotOutput(ns("pca_plot"), height = "550px",
                                      width = "800px"
                    )
                  )
                )
            ),
            bslib::nav_panel(
              title = "UMAP",
              icon = shiny::icon("layer-group"),
              bslib::card_header(
                div(
                  style = "display: flex; justify-content: flex-end;
            align-items: center;",
                  shiny::downloadButton(
                    class = "clickable-button-light",
                    outputId = ns("umap_download_plot"),
                    label = ""
                  )
                )),
              bslib::card_body(
                style = "display = flex; justify-content: space-between;
            align-items: center;",
                div(
                  style = "width: 800px; height: 550px; overflow: auto;",
                  shiny::plotOutput(ns("umap_plot"), height = "550px",
                                    width = "800px"
                  )
                )
              )
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
  observe({
    shiny::updateSelectInput(
      session,
      "data",
      choices = names(data),
      selected = names(data)[1]
    )
  })

  dr_df <- reactive({
    data[[input$data]]
  })

  output$pca_plot <- shiny::renderPlot({
    df <- dr_df()

    scater::plotPCA(
      df,
      colour_by = "Condition",
      point_size = 3
    ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      ) +
      ggplot2::labs(title = "PCA Plot")
  })

  output$umap_plot <- shiny::renderPlot({
    df <- dr_df()

    # Step 4: plot UMAP
    scater::plotUMAP(
      df,
      colour_by = "Condition",
      point_size = 3
    ) +
      ggplot2::theme_bw() +
      ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      ) +
      ggplot2::labs(title = "UMAP Plot")
  })
}
