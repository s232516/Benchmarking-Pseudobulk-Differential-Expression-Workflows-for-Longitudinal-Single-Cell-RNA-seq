#' mod_errorbar UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd

mod_errorbars_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::layout_sidebar(
      fillable = TRUE,
      sidebar = bslib::accordion(
        id = ns("eb_accordion"),
        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),
          id = ns("eb_analysis"),
          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select dataset",
                                       shiny::icon("circle-info")),
              "Select which is the dataset that you want to use for
              dimensionality reduction."
            ),
            choices = NULL
          ),
          shiny::selectInput(
            inputId = ns("de_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList("Select dataset",
                                       shiny::icon("circle-info")),
              "Select which is the dataset that you want to use for
              dimensionality reduction."
            ),
            choices = NULL
          )
        )
      ),
      bslib::layout_columns(
        col_widths = 12,
        bslib::card(
          title = "Error Bars Plot",
          icon = shiny::icon("volcano"),
          id = ns("eb_card"),
          # Create a header for downloading the plot
          bslib::card_header(
            shiny::downloadButton(
              outputId = ns("download_plot"),
              class = "btn-primary"
            )
          ),
          shiny::plotOutput(ns("eb_plot"), height = "600px")
        )
      )
    )
  )
}

#' mod_erorbar Server Functions
#' @noRd
mod_errorbars_server <- function(input, output, session, data, sim) {
  ns <- session$ns
  observe({
    shiny::updateSelectInput(
      session,
      'de_method',
      choices = unique(data$method),
      selected = unique(data$method)[1]
    )
  })

  observe({
    shiny::updateSelectInput(
      session,
      'data',
      choices = unique(data$analysis),
      selected =  unique(data$analysis)[1]
    )
  })

  # UPDATE IF I HAVE MORE DATA

  real_data <- reactive({
    expr <- as.matrix(SingleCellExperiment::logcounts(sim))

    condition <- SingleCellExperiment::colData(sim)$Condition

    true_effects <- data.frame(
      gene = rownames(expr),
      TrueEffect = apply(expr, 1, function(x) {
        mean(x[condition == "Treatment"]) - mean(x[condition == "Placebo"])
      })
    )
  })

  eb_data <- reactive({
    df <- data |>
      dplyr::filter(method == input$de_method) |>
      dplyr::filter(analysis == input$data)
  })

  errorbar_data <- reactive({
    final <- eb_data() |>
    dplyr::mutate(
      EstimatedEffect = log2FoldChange,
      SE = abs(log2FoldChange / stat)
    )  |>
      dplyr::select(gene, EstimatedEffect, SE)
  })

  merged_data <- reactive({
    df <- errorbar_data() |>
      dplyr::left_join(real_data(),
                       by = "gene")
  })


  # Render Plot
  output$eb_plot <- shiny::renderPlot({
    df <-  merged_data()

    p <- ggplot2::ggplot(df, ggplot2::aes(
      x = TrueEffect,
      y = EstimatedEffect
    )) +
      ggplot2::geom_point(alpha = 0.8) +
      ggplot2::geom_errorbar(
        ggplot2::aes(ymin = EstimatedEffect - 1.96 * SE,
                     ymax = EstimatedEffect + 1.96 * SE),
        width = 0,
        alpha = 0.6
      ) +
      ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
      ggplot2::labs(
        title = "True vs Estimated Effect Sizes with 95% CI",
        x = "True Simulated Effect",
        y = "Estimated Effect"
      ) +
      ggplot2::theme_bw()
    p
  })
}
