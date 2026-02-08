#' mod_filter_summary UI
#'
#' @param id shiny id
#' @noRd
mod_filter_summary_ui <- function(id) {
  ns <- shiny::NS(id)

  bslib::page_fillable(
    bslib::layout_columns(
      col_widths = 12,

      bslib::card(
        bslib::card_header("Filtered genes – DE composition"),
        bslib::card_body(
          DT::DTOutput(ns("filter_summary_table"))
        )
      )
    )
  )
}
#' mod_filter_summary server
#'
#' @noRd

#' mod_filter_summary Server Function
#'
#' @noRd
mod_filter_summary_server <- function(input, output, session,
                                      compare_data,
                                      ground_truth) {

  # ---- Safety checks ----
  stopifnot(is.list(ground_truth))
  stopifnot(all(c("analysis", "method", "formula", "gene", "trueDE") %in%
                  colnames(compare_data)))

  # ---- Core computation ----
  filter_summary <- reactive({

    compare_data |>
      dplyr::group_by(analysis, method, formula) |>
      dplyr::group_modify(function(df, keys) {

        sim <- keys$analysis

        # Ground truth for this simulation
        gt <- ground_truth[[sim]]
        validate(
          need(!is.null(gt),
               paste("Ground truth not found for", sim))
        )

        total_genes <- nrow(gt)

        # Genes kept after filtering
        kept_genes <- unique(df$gene)
        n_kept <- length(kept_genes)

        # Genes filtered out
        filtered_genes <- setdiff(gt$gene, kept_genes)
        n_filtered <- length(filtered_genes)

        # True DE genes in ground truth
        trueDE_genes <- gt$gene[gt$trueDE]

        # Filtered genes that were truly DE
        filtered_trueDE <- intersect(filtered_genes, trueDE_genes)
        n_filtered_trueDE <- length(filtered_trueDE)

        tibble::tibble(
          n_total_genes = total_genes,
          n_kept_genes = n_kept,
          n_filtered_genes = n_filtered,
          n_filtered_trueDE = n_filtered_trueDE,
          prop_filtered_trueDE =
            ifelse(n_filtered > 0,
                   n_filtered_trueDE / n_filtered,
                   NA_real_)
        )
      }) |>
      dplyr::ungroup()
  })

  # ---- Render table ----
  output$filter_summary_table <- DT::renderDT({

    df <- filter_summary()

    DT::datatable(
      df,
      rownames = FALSE,
      options = list(
        pageLength = 10,
        autoWidth = TRUE,
        scrollX = TRUE
      )
    )
  })
}

