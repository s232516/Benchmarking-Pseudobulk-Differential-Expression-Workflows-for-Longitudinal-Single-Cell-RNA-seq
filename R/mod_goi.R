#' mod_goi UI Function
#'
#' @description A shiny Module.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd
#'
#' @importFrom shiny NS tagList
#' mod_goi UI Function
#'
#' @description Gene-of-Interest LFC trajectories across timepoints and treatment arms.
#'
#' @param id,input,output,session Internal parameters for {shiny}.
#'
#' @noRd

#' mod_goi_sce UI Function
#'
#' @description Gene-of-Interest trajectories over time from SCE objects.
#'
#' @param id Shiny module id
#'
#' @noRd
mod_goi_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_sidebar(
      width = 280,
      open = "always",
      sidebar = bslib::accordion(
        id = ns("goi_accordion"),
        bslib::accordion_panel(
          title = "GOI settings",
          icon  = shiny::icon("dna"),

          shiny::selectInput(
            inputId = ns("analysis"),
            label   = "Select dataset:",
            choices = NULL
          ),

          shiny::selectInput(
            inputId = ns("gene"),
            label   = "Select gene:",
            choices = NULL
          )
        )
      ),
      bslib::layout_columns(
        col_widths = 12,

        bslib::card(
          style = "margin-bottom: 1rem;",
          bslib::card_header(
            div(
              style = "display:flex; justify-content:space-between; align-items:center;",
              "Gene trajectory over time",
              shiny::downloadButton(
                outputId = ns("download_plot"),
                label    = "",
                class    = "clickable-button-light"
              )
            )
          ),
          bslib::card_body(
            style = "display = flex; justify-content: space-between;
            align-items: center;",
            shiny::plotOutput(
              outputId = ns("goi_plot"),
              height = "550px",
              width = "800px"
            )
          )
        )
      )
    )
  )
}





#' mod_goi_sce Server Functions
#'
#' @param id Shiny module id
#' @param sce_list named list of SCE objects, e.g. list(jan = sce_jan, rescuesim = sce_res)
#'
#' @noRd
mod_goi_server <- function(input, output, session, sce_list) {
    ns <- session$ns

    observe({
      available <- intersect(names(sce_list), c("jan", "rescueSim"))
      shiny::updateSelectInput(
        session,
        "analysis",
        choices = available,
        selected = available[1]
      )
    })


    observeEvent(input$analysis, {
      req(input$analysis)
      sce <- sce_list[[input$analysis]]
      gene_names <- rownames(sce)  # change if you use rowData column instead
      shiny::updateSelectInput(
        session,
        "gene",
        choices = sort(gene_names),
        selected = gene_names[1]
      )
    })

    goi_summary <- reactive({
      req(input$analysis, input$gene)

      sce <- sce_list[[input$analysis]]

      # ---- column names: adjust here if needed ----
      cd <- as.data.frame(SummarizedExperiment::colData(sce))

      stopifnot(all(c("subjectID", "group", "time") %in% names(cd)))

      # extract expression for this gene across cells
      expr_vec <- as.numeric(SingleCellExperiment::logcounts(sce)[input$gene, ])

      cell_df <- dplyr::tibble(
        expr      = expr_vec,
        subjectID = cd$subjectID,
        group     = cd$group,
        time      = cd$time
      )

      # average over cells → one value per subject × group × time
      subj_df <- cell_df |>
        dplyr::group_by(subjectID, group, time) |>
        dplyr::summarise(
          expr_mean = mean(expr, na.rm = TRUE),
          .groups   = "drop"
        )

      # average over subjects within each group × time + compute CI
      summary_df <- subj_df |>
        dplyr::group_by(group, time) |>
        dplyr::summarise(
          mean_expr = mean(expr_mean, na.rm = TRUE),
          sd_expr   = stats::sd(expr_mean, na.rm = TRUE),
          n_subj    = dplyr::n(),
          se_expr   = sd_expr / sqrt(n_subj),
          ci_low    = mean_expr - 1.96 * se_expr,
          ci_high   = mean_expr + 1.96 * se_expr,
          .groups   = "drop"
        )

      # ensure time is ordered (if it’s something like Baseline, EoT, etc.)
      summary_df$time <- factor(summary_df$time, levels = unique(summary_df$time))

      summary_df
    })

    output$goi_plot <- shiny::renderPlot({
      df<- goi_summary()
      # prevent evalmod from running on empty data
      validate(
        need(!is.null(df) && nrow(df) > 0,
             "Loading...")
      )
      df <- goi_summary()
      validate(need(nrow(df) > 0, "No data available for this gene."))
      pastel_cols <- c(
        "Treatment" = "#FFB3BA",
        "Placebo"   = "#BAE1FF"
      )


      ggplot2::ggplot(
        df,
        ggplot2::aes(
          x     = time,
          y     = mean_expr,
          group = group,
          color = group,
          fill  = group
        )
      ) +
        # error bars for CIs
        ggplot2::geom_errorbar(
          ggplot2::aes(x = time, ymin = ci_low, ymax = ci_high),
          width    = 0.2,
          linewidth = 1.2
        ) +
        ggplot2::geom_line(size = 1.2) +
        ggplot2::geom_point(size = 3) +
        ggplot2::labs(
          title = paste("Trajectory for", input$gene, "in", input$analysis),
          x     = "Timepoint",
          y     = "Mean log-expression",
          color = "Treatment",
          fill  = "Treatment"
        ) +
        ggplot2::scale_fill_manual(values = pastel_cols) +
        ggplot2::theme_bw(base_size = 16) +
        ggplot2::theme(
          axis.text.x     = ggplot2::element_text(angle = 30, hjust = 1),
          legend.position = "right",
          aspect.ratio    = 1
        )
    })

    output$download_plot <- shiny::downloadHandler(
      filename = function() {
        paste0("GOI_", input$gene, "_", input$analysis, ".png")
      },
      content = function(file) {
        p <- output$goi_plot()  # uses the same plot
        ggplot2::ggsave(file, p, width = 8, height = 6)
      }
    )
}
