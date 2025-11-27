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

  bslib::page_fillable(
    bslib::layout_sidebar(
      fillable = TRUE,
      width = 280,
      open = "always",

      # Sidebar controls
      sidebar = bslib::accordion(
        id = ns("tx_accordion"),
        bslib::accordion_panel(
          title = "Analysis",
          icon = shiny::icon("bars"),
          id = ns("tx_analysis"),
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
            inputId = ns("formula"),
            label = bslib::tooltip(
              trigger = shiny::tagList(
                "Select formula",
                shiny::icon("circle-info")
              ),
              "Select which is the formula that you want to use for
              dimensionality reduction."
            ),
            choices = NULL
          ),
          shiny::checkboxInput(
            inputId = ns("highlight_truede"),
            label = "Highlight true DE genes",
            value = FALSE
          ),
          shiny::radioButtons(
            inputId = ns("pvalue"),
            label = "P-value",
            choices = c(0.1, 0.05, 0.01, 0.005, 0.001),
            selected = 0.05
          ),
          shiny::sliderInput(
            inputId = ns("logfc"),
            label = "Log2 Fold Change",
            min = 0,
            max = 4,
            value = 0
          )
        )
      ),

      # MAIN AREA: scrollable container with plot (fixed height) + table below
      bslib::layout_columns(
        col_widths = 12,
        div(
          style = "overflow-y: auto; padding-right: 10px; padding-left: 10px;",
          # volcano plot card
          bslib::card(
            style = "margin-bottom: 1rem;",
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
              )
            ),
            bslib::card_body(
              div(
                style = "width: 100%; height: 550px; overflow: hidden;",
                shiny::plotOutput(ns("volcano_plot"), height = "550px",
                                  width = "100%")
              )
            )
          ),

          # table card below (scrollable to it)
          bslib::card(
            bslib::card_header("Volcano Plot Data Table"),
            bslib::card_body(
              div(
                style = "width: 100%",
                DT::DTOutput(ns("volcano_table"))
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

  observeEvent(input$data, {
    df <- data |>
      dplyr::filter(analysis == input$data)
    choices <- unique(df$method)
    shiny::updateSelectInput(
      session,
      "de_method",
      choices = choices,
      selected = choices[1]
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

  observeEvent(c(input$data,input$de_method), {
    df <- data |>
      dplyr::filter(analysis == input$data,
                    method == input$de_method)
    choices <- unique(df$formula)
    shiny::updateSelectInput(
      session,
      "formula",
      choices = choices,
      selected = choices[1]
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
    req( input$de_method, input$data)
    df <- data |>
      dplyr::filter(method == input$de_method) |>
      dplyr::filter(analysis == input$data) |>
      dplyr::filter(formula == input$formula) |>
      dplyr::mutate(neg_log10_pvalue = -log10(.data$pvalue),
                    is_true_de = ifelse(trueDE == 1, TRUE, FALSE)) |>
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
    if(nrow(volcano_data()) != 0){
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

    # If highlight is active, override significance colors for true DE genes
    if (input$highlight_truede) {
      df$significance <- ifelse(df$is_true_de, "true DE", df$significance)
      counts <- df |>
        dplyr::count(significance) |>
        tibble::deframe()

      # build legend labels with counts
      legend_labels <- c(
        "true DE" = paste0("true DE (", counts["true DE"], ")"),
        "Upregulated" = paste0("Upregulated (", counts["Upregulated"], ")"),
        "Downregulated" = paste0("Downregulated (", counts["Downregulated"], ")"),
        "Not significant" = paste0("Not significant (", counts["Not significant"], ")")
      )

      color_order <- c(
        "true DE" = "#C9980B",       # gold
        "Upregulated" = "#b2182b",
        "Downregulated" = "#2166ac",
        "Not significant" = "#c0c0c0"
      )
    } else {
      color_order <- c(
        "Upregulated" = "#b2182b",
        "Downregulated" = "#2166ac",
        "Not significant" = "#c0c0c0"
      )

      counts <- df |>
        dplyr::count(significance) |>
        tibble::deframe()

      # build legend labels with counts
      legend_labels <- c(
        "Upregulated" = paste0("Upregulated (", counts["Upregulated"], ")"),
        "Downregulated" = paste0("Downregulated (", counts["Downregulated"], ")"),
        "Not significant" = paste0("Not significant (", counts["Not significant"], ")")
      )


    }


    volcano_plot(
      data = df,
      x = "log2FoldChange",
      y = "pvalue",
      color_col = "significance",
      colors = color_order,
      color_labels = legend_labels,
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
    }
  })

  output$volcano_table <- DT::renderDT({
    df <- volcano_data()

    DT::datatable(
      df,
      options = list(
        pageLength = 10,
        scrollX = TRUE,
        searching = TRUE,
        ordering = TRUE
        )
    )
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
