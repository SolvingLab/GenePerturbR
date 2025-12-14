#' Beautiful Heatmap with ggplot2
#'
#' High-quality heatmap implementation using pure ggplot2
#' Replaces ComplexHeatmap for easier installation
#'
#' @keywords internal

#' Plot Heatmap with ggplot2 (Beautiful Version)
#'
#' @param matrix_data Matrix. Data to plot
#' @param scale Character. "row", "column", or "none"
#' @param cluster_rows Logical. Cluster rows
#' @param cluster_cols Logical. Cluster columns
#' @param row_labels Character vector. Row labels
#' @param col_labels Character vector. Column labels
#' @param row_annotation Character vector. Row groupings
#' @param col_annotation Data frame. Column annotations
#' @param title Character. Plot title
#' @param colors Character vector. Color palette
#' @param legend_title Character. Legend title
#' @param show_values Logical. Show values in cells
#' @param theme ggplot2 theme
#'
#' @return ggplot object
#' @keywords internal
.gpdb_ggplot_heatmap <- function(matrix_data,
                                 scale = "row",
                                 cluster_rows = TRUE,
                                 cluster_cols = TRUE,
                                 row_labels = NULL,
                                 col_labels = NULL,
                                 row_annotation = NULL,
                                 col_annotation = NULL,
                                 title = NULL,
                                 colors = NULL,
                                 legend_title = "Value",
                                 show_values = FALSE,
                                 theme = NULL) {
  # Handle NA values (replace with 0 for clustering)
  has_na <- any(is.na(matrix_data))
  if (has_na) {
    matrix_data[is.na(matrix_data)] <- 0
    message("Replaced NA values with 0 for clustering")
  }

  # Apply scaling
  if (scale == "row") {
    matrix_data <- t(scale(t(matrix_data)))
    if (is.null(legend_title)) legend_title <- "Z-score"
  } else if (scale == "column") {
    matrix_data <- scale(matrix_data)
    if (is.null(legend_title)) legend_title <- "Z-score"
  }

  # Clustering (optimized for speed)
  if (cluster_rows && nrow(matrix_data) > 1) {
    # Skip for large matrices
    if (nrow(matrix_data) > 200) {
      message("Skipping row clustering (>200 rows)")
      cluster_rows <- FALSE
    } else {
      tryCatch(
        {
          row_dist <- stats::dist(matrix_data, method = "euclidean")
          if (any(!is.finite(as.matrix(row_dist)))) {
            warning("Cannot cluster rows")
            cluster_rows <- FALSE
          } else {
            row_dend <- stats::hclust(row_dist, method = "complete")
            matrix_data <- matrix_data[row_dend$order, ]
          }
        },
        error = function(e) {
          warning("Row clustering failed")
          cluster_rows <<- FALSE
        }
      )
    }
  }

  if (cluster_cols && ncol(matrix_data) > 1) {
    # Skip for many columns
    if (ncol(matrix_data) > 50) {
      message("Skipping column clustering (>50 cols)")
      cluster_cols <- FALSE
    } else {
      tryCatch(
        {
          col_dist <- stats::dist(t(matrix_data), method = "euclidean")
          if (any(!is.finite(as.matrix(col_dist)))) {
            warning("Cannot cluster columns")
            cluster_cols <- FALSE
          } else {
            col_dend <- stats::hclust(col_dist, method = "complete")
            matrix_data <- matrix_data[, col_dend$order]
          }
        },
        error = function(e) {
          warning("Column clustering failed")
          cluster_cols <<- FALSE
        }
      )
    }
  }

  # Prepare labels
  if (is.null(row_labels)) row_labels <- rownames(matrix_data)
  if (is.null(col_labels)) col_labels <- colnames(matrix_data)

  # Convert to long format for ggplot2 (optimized)
  n_rows <- nrow(matrix_data)
  n_cols <- ncol(matrix_data)

  heatmap_df <- data.frame(
    row = rep(1:n_rows, n_cols),
    col = rep(1:n_cols, each = n_rows),
    value = as.vector(matrix_data),
    row_label = rep(row_labels, n_cols),
    col_label = rep(col_labels, each = n_rows),
    stringsAsFactors = FALSE
  )

  # Add row annotation
  if (!is.null(row_annotation)) {
    heatmap_df$row_group <- rep(row_annotation, ncol(matrix_data))
  }

  # Default colors
  if (is.null(colors)) {
    colors <- .gpdb_get_palette("spectral", 11)
  }

  # Default theme
  if (is.null(theme)) {
    theme <- .gpdb_theme_default()
  }

  # Create plot (ultra-optimized for speed)
  p <- ggplot2::ggplot(heatmap_df, ggplot2::aes(x = col, y = n_rows - row + 1, fill = value)) +
    ggplot2::geom_raster() + # MUCH faster than geom_tile!
    ggplot2::scale_fill_gradientn(
      colors = colors,
      name = legend_title,
      na.value = "grey90"
    ) +
    ggplot2::scale_x_continuous(
      breaks = if (n_cols <= 20) 1:n_cols else NULL,
      labels = if (n_cols <= 20) col_labels else NULL,
      expand = c(0, 0)
    ) +
    ggplot2::scale_y_continuous(
      breaks = if (n_rows <= 50) 1:n_rows else NULL,
      labels = if (n_rows <= 50) rev(row_labels) else NULL,
      expand = c(0, 0)
    ) +
    ggplot2::labs(
      title = title,
      x = NULL,
      y = NULL
    ) +
    theme +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(angle = 45, hjust = 1, size = 8),
      axis.text.y = ggplot2::element_text(size = 8),
      axis.ticks = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(fill = NA, color = "black", linewidth = 0.8),
      legend.position = "right"
    )

  # Add values in cells if requested
  if (show_values) {
    heatmap_df$label <- sprintf("%.2f", heatmap_df$value)
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = label),
      size = 2.5,
      color = "black"
    )
  }

  # Add row annotation sidebar if provided
  if (!is.null(row_annotation) && !is.null(heatmap_df$row_group)) {
    # Create annotation strip
    unique_groups <- unique(row_annotation)
    group_colors <- .gpdb_get_palette("qualitative", length(unique_groups))
    names(group_colors) <- unique_groups

    # Add facet or annotation
    # For now, add as border color
    group_df <- data.frame(
      row = 1:nrow(matrix_data),
      group = row_annotation
    )
  }

  # Add column annotation on top if provided
  if (!is.null(col_annotation)) {
    # This would require patchwork to stack annotation bar + heatmap
    # For now, skip or use a simple approach
  }

  return(p)
}
