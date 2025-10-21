#' helpers
#'
#' @description A volcano plot function
#'
#' @return The return value, if any, from executing the utility.
#'
#' @noRd


volcano_plot <- function(data,
                         x = "effect_size",
                         y = "pvalue",
                         color_col = NULL,
                         colors = NULL,
                         color_labels = NULL,
                         color_legend_title = NULL,
                         num_labels = 10,
                         custom_labels = NULL,
                         labels = NULL,
                         label_col = NULL, # deprecated
                         label_size = 3,
                         point_size = 3,
                         hline = 0.05,
                         vline = 0,
                         plot_title = NULL,
                         plot_subtitle = NULL,
                         plot_caption = NULL,
                         x_axis_title = x,
                         y_axis_title = y,
                         x_lims = NULL,
                         y_breaks = NULL) {

  # If all values are NA in y column, warn the user

  if (all(is.na(data[[y]]))) {
    warning(paste("All values in", y, "column are NA"))
  }

  # Transform y values to avoid underflow

  data[[y]] <- ifelse(
    data[[y]] < 1e-300,
    yes = 1e-300,
    no = data[[y]]
  )

  # Set x_lim

  if (is.null(x_lims)) {
    max_x <- ceiling(max(abs(data[[x]]), na.rm = TRUE))
    x_lims <- c(-max_x, max_x)
  }

  # Set y axis ticks

  if (is.null(y_breaks)) {

    # Set number of y axis ticks
    n_ticks <- 5
    y_breaks <-
      seq(
        from = max(-log10(data[[y]]), na.rm = TRUE),
        to = 0,
        length.out = n_ticks
      )[1:(n_ticks - 1)] |>
      ceiling()

    # Add scientific notation
    y_breaks <- c(as.numeric(paste0("1e-", y_breaks)), 1)

  }

  if (is.null(color_col)) {

    p <- data |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data[[x]],
          y = .data[[y]]
        )
      )

  } else {

    p <- data |>
      ggplot2::ggplot(
        ggplot2::aes(
          x = .data[[x]],
          y = .data[[y]],
          color = .data[[color_col]]
        )
      )
  }

  # Add horizontal/vertical lines

  if (!is.null(hline)) {

    p <- p +
      ggplot2::geom_hline(
        yintercept = hline,
        lty = 2,
        col = "grey"
      )
  }

  if (!is.null(vline)) {
    p <- p +
      ggplot2::geom_vline(
        xintercept = c(-vline, vline),
        lty = 2,
        col = "grey"
      )
  }

  # Add points layer

  p <- p +
    ggplot2::geom_point(size = point_size, alpha = 0.8)



  # Add label layer

  if (!is.null(labels) || !is.null(label_col)) {

    # Determine effective label column, prioritizing the 'labels' argument

    effective_label_col <- ifelse(!is.null(labels), labels, label_col)



    # Select the lowest num_labels based on p-value

    label_data <- data |>
      dplyr::arrange(.data[[y]]) |>
      dplyr::slice_head(n = num_labels)



    # Incorporate custom labels if provided

    if (!is.null(custom_labels)) {

      label_data <- label_data |>
        dplyr::bind_rows(
          data |>
            dplyr::filter(.data[[effective_label_col]] %in% custom_labels)
        ) |>
        dplyr::filter(!is.na(.data[[effective_label_col]])) |>
        dplyr::distinct()
    }



    # Create base aesthetic mapping

    label_aes <- ggplot2::aes(label = .data[[effective_label_col]])

    # Add color aesthetic only if color_col exists

    if (!is.null(color_col)) {
      label_aes$colour <- rlang::expr(.data[[color_col]])
    }
    p <- p +
      ggrepel::geom_label_repel(
        data = label_data,
        mapping = label_aes,
        size = label_size,
        hjust = 0.5,
        vjust = 0.5,
        max.overlaps = 50,
        min.segment.length = grid::unit(0, "lines"),
        alpha = 0.9,
        show.legend = FALSE
      )
  }

  # Transform axes

  p <- p +
    ggplot2::scale_y_continuous(
      breaks = y_breaks,
      expand = c(0.05, 0.05),
      trans = neglog_trans(base = 10)
    ) +
    ggplot2::scale_x_continuous(
      limits = x_lims,
      breaks = scales::breaks_extended(5),
    )



  # Add labels

  p <- p +

    ggplot2::labs(
      x = x_axis_title,
      y = y_axis_title,
      color = color_legend_title,
      title = plot_title,
      subtitle = plot_subtitle,
      caption = plot_caption

    )



  # Control colors

  if (!is.null(colors)) {
    if (!is.null(color_labels)) {

      p <- p +
        ggplot2::scale_color_manual(
          values = colors,
          labels = color_labels
        )

    } else {

      p <- p +
        ggplot2::scale_color_manual(
          values = colors
        )
    }
  }

  # Theme adjustments

  p <- p +
    ggplot2::theme_bw() +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank())

  p
}

# Scale y neg log10 helper functions

neglog_breaks <- function(n = 5, base = 10) {
  function(x) {
    rng <- rev(-log(range(x, na.rm = TRUE), base = base))
    min <- floor(rng[1])
    max <- ceiling(rng[2])
    if (max == min) {
      return(base^min)
    }
    by <- floor((max - min) / n) + 1
    base^(seq(min, max, by = by))
  }
}

neglog_trans <- function(base = exp(1)) {
  scales::trans_new("neglog",
                    transform = function(x) -log(x, base),
                    inverse = function(x) base^(-x),
                    breaks = neglog_breaks(base = base),
                    domain = c(1e-400, Inf)

  )

}
