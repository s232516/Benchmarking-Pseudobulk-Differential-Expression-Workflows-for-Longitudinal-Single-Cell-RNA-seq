#' mod_dr UI Function
#'
#' @description A shiny Module for dimensionality reduction.
#'
#' @param id Internal parameter for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_dr_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::card(
    bslib::layout_sidebar(
      fillable = TRUE,

      # ---- Sidebar ----
      sidebar = bslib::accordion(
        id = ns("dr_accordion"),
        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),
          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select dataset",
                shiny::icon("circle-info")
              ),
              "Select the dataset used for dimensionality reduction."
            ),
            choices = character(0)
          )
        )
      ),

      # ---- Main panel ----
      bslib::layout_columns(
        col_widths = 12,
        bslib::navset_card_tab(
          # ===== PCA =====
          bslib::nav_panel(
            title = "PCA",
            icon = shiny::icon("chart-simple"),
            bslib::card_header(
              div(
                style = "display: flex; justify-content: flex-end;",
                shiny::downloadButton(
                  outputId = ns("pca_download_plot"),
                  label = "Download",
                  class = "clickable-button-light"
                )
              )
            ),
            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("pca_plot"),
                height = "550px",
                width = "800px"
              )
            )
          ),

          # ===== UMAP =====
          bslib::nav_panel(
            title = "UMAP",
            icon = shiny::icon("layer-group"),
            bslib::card_header(
              div(
                style = "display: flex; justify-content: flex-end;",
                shiny::downloadButton(
                  outputId = ns("umap_download_plot"),
                  label = "Download",
                  class = "clickable-button-light"
                )
              )
            ),
            bslib::card_body(
              shiny::plotOutput(
                outputId = ns("umap_plot"),
                height = "550px",
                width = "800px"
              )
            )
          )
        )
      )
    )
  )
}

#' mod_dr Server Function
#'
#' @noRd
mod_dr_server <- function(input, output, session, data) {
  # ---- Update dataset selector ----
  observe({
    shiny::updateSelectInput(
      session = session,
      inputId = "data",
      choices = names(data),
      selected = names(data)[1]
    )
  })

  # ---- Reactive dataset ----
  dr_df <- reactive({
    req(input$data)
    data[[input$data]]
  })

  # ---- PCA plot object ----
  pca_plot_obj <- reactive({
    req(dr_df())

    scater::plotPCA(
      dr_df(),
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

  # ---- UMAP plot object ----
  umap_plot_obj <- reactive({
    req(dr_df())

    scater::plotUMAP(
      dr_df(),
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

  # ---- Render plots ----
  output$pca_plot <- shiny::renderPlot({
    pca_plot_obj()
  })

  output$umap_plot <- shiny::renderPlot({
    umap_plot_obj()
  })

  # ---- Desktop/plots directory ----
  plots_dir <- file.path(
    path.expand("~"),
    "Desktop",
    "plots"
  )

  if (!dir.exists(plots_dir)) {
    dir.create(plots_dir, recursive = TRUE)
  }

  # ---- PCA download ----
  output$pca_download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0("PCA_", input$data, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0("PCA_", input$data, "_", Sys.Date(), ".png")
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = pca_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = pca_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )
    }
  )

  # ---- UMAP download ----
  output$umap_download_plot <- shiny::downloadHandler(
    filename = function() {
      paste0("UMAP_", input$data, "_", Sys.Date(), ".png")
    },
    content = function(file) {
      desktop_file <- file.path(
        plots_dir,
        paste0("UMAP_", input$data, "_", Sys.Date(), ".png")
      )

      ggplot2::ggsave(
        filename = desktop_file,
        plot = umap_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )

      ggplot2::ggsave(
        filename = file,
        plot = umap_plot_obj(),
        width = 8,
        height = 8,
        dpi = 300
      )
    }
  )
}
