#' mod_volcano UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
mod_volcano_ui <- function(id) {
  ns <- shiny::NS(id)
  bslib::card(
    bslib::layout_sidebar(
      fillable = TRUE,
      sidebar = bslib::accordion(
        id = ns("tx_accordion"),
        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),
          id = ns("tx_analysis"),
          shiny::selectInput(
            inputId = ns("data"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select dataset",
                shiny::icon("circle-info")
              ),
              "Select which is the dataset that you want to use for
              dimensionality reduction."
            ),
            choices = NULL
          ),
          shiny::selectInput(
            inputId = ns("de_method"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select dataset",
                shiny::icon("circle-info")
              ),
              "Select which is the dataset that you want to use for
              dimensionality reduction."
            ),
            choices = NULL
          ),
          shiny::radioButtons(
            inputId = ns("pvalue"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "P-value",
                shiny::icon("circle-info")
              ),
              "Select the p-value to use in the analysis."
            ),
            choices = c(0.1, 0.05, 0.01, 0.005, 0.001),
            selected = 0.05
          ),
          shiny::sliderInput(
            inputId = ns("logfc"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Log2 Fold Change",
                shiny::icon("circle-info")
              ),
              "Select the log2 fold change to use in the analysis."
            ),
            min = 0,
            max = 4,
            value = 0
          )
        )
      ),
      bslib::layout_columns(
        col_widths = 12,
        bslib::card(
          fill = FALSE,
          fill_screen = TRUE,
          # Create a header for downloading the plot
          bslib::card_header(
            div(
            style = "display: flex; justify-content: space-between;
            align-items: center;",
            "Volcano Plot",
            shiny::downloadButton(
              class = "clickable-button-light",
              outputId = ns("download_plot"),
              label = ""
            )
          )),
          bslib::card_body(
            style = "display = flex; justify-content: space-between;
            align-items: center;",
            div(
              style = "width: 800px; height: 550px; overflow: auto;",
              shiny::plotOutput(ns("volcano_plot"), height = "550px",
                                width = "800px"
              )
            )
          )
        )
      )
    )
  )
}
#' mod_volcano Server Functions
#'
#' @noRd
mod_volcano_server <- function(input, output, session, data) {
  ns <- session$ns

  observe({
    shiny::updateSelectInput(
      session,
      "de_method",
      choices = unique(data$method),
      selected = unique(data$method)[1]
    )
  })

  observe({
    shiny::updateSelectInput(
      session,
      "data",
      choices = unique(data$analysis),
      selected = unique(data$analysis)[1]
    )
  })


  shiny::updateSliderInput(
    session,
    "logfc",
    min = 0,
    max = round(max(abs(data$log2FoldChange),
      na.rm = TRUE
    )),
    value = 0,
    step = 0.01
  )
  # Bonferroni Correction

  volcano_data <- reactive({
    df <- data |>
      dplyr::filter(method == input$de_method) |>
      dplyr::filter(analysis == input$data) |>
      dplyr::mutate(neg_log10_pvalue = -log10(.data$pvalue)) |>
      dplyr::arrange(.data$pvalue) |>
      dplyr::mutate(pvalue_adj = dplyr::case_when(
        (.data$padj <= as.numeric(input$pvalue)) ~ .data$pvalue,
        TRUE ~ 0
      )) |>
      dplyr::mutate(significance = dplyr::case_when(
        (.data$padj <= as.numeric(input$pvalue)) & (.data$log2FoldChange >
                                                      input$logfc)
        ~ "Upregulated",
        (.data$padj <= as.numeric(input$pvalue)) & (.data$log2FoldChange <
                                                      -input$logfc)
        ~ "Downregulated",
        TRUE ~ "Not significant"
      ))
  })

  output$volcano_plot <- shiny::renderPlot({
    req(volcano_data())
    df <- as.data.frame(volcano_data())
    max_es <- ceiling(max(abs(df$log2FoldChange),
      na.rm = TRUE
    ))
    max_pvalue <- max(df$pvalue_adj)
    y_thres <- NULL
    if (max_pvalue > 0) {
      y_thres <- max_pvalue
    }

    x_lims <- c(-max_es, max_es)


    color_order <- c(
      "Upregulated" = "#b2182b", "Downregulated" = "#2166ac",
      "Not significant" = "#c0c0c0"
    )

    volcano_plot(
      data = df,
      x = "log2FoldChange",
      y = "pvalue",
      color_col = "significance",
      colors = color_order,
      labels = "gene",
      y_axis_title = "P-value",
      x_axis_title = "Effect Size",
      plot_title = paste("Volcano Plot -", input$de_method),
      x_lims = x_lims,
      hline = y_thres,
      vline = input$logfc,
      label_size = 5
    ) + ggplot2::theme(
      aspect.ratio = 1,
      text = ggplot2::element_text(size = 18)
    )
  })

  output$volcanoTable <- DT::renderDT({
    df <- volcano_data()
    datatable(df, options = list(pageLength = 10))
  })

  output$download_plot <- shiny::downloadHandler(
    filename = function() {
      paste("volcano_plot", ".png", sep = "")
    },
    content = function(file) {
      df <- as.data.frame(volcano_data())

      color_order <- c(
        "Upregulated" = "#b2182b", "Downregulated" = "#2166ac",
        "Not significant" = "#c0c0c0"
      )

      volcano_plot(
        data = df, x = "log2FoldChange", y = "pvalue",
        color_col = "significance", colors = color_order,
        labels = "gene", y_axis_title = "P-value",
        x_axis_title = "Effect Size", hline = -log10(input$pvalue),
        vline = input$logfc
      ) + ggplot2::theme(
        aspect.ratio = 1,
        text = ggplot2::element_text(size = 18)
      )
    }
  )
}
