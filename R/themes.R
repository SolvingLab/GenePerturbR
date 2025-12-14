#' Theme and Color System for GenePerturbR
#'
#' Unified theme and color palette system for all visualizations
#'
#' @keywords internal

#' GenePerturbR Default Theme
#'
#' @param base_size Base font size
#' @param base_family Base font family
#'
#' @return ggplot2 theme
#' @keywords internal
.gpdb_theme_default <- function(base_size = 12, base_family = "") {
  ggplot2::theme_classic(base_size = base_size, base_family = base_family) +
    ggplot2::theme(
      # Plot title
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = base_size + 2,
        face = "bold",
        color = "black"
      ),
      # Axis
      axis.text = ggplot2::element_text(size = base_size - 2, color = "black"),
      axis.title = ggplot2::element_text(size = base_size, face = "bold", color = "black"),
      axis.line = ggplot2::element_line(color = "black", linewidth = 0.6),
      # Legend
      legend.position = "bottom",
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = base_size - 1, color = "black"),
      # Panel
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white", color = NA)
    )
}

#' GenePerturbR Color Palettes
#'
#' @keywords internal
.gpdb_colors <- function() {
  list(
    # 11-color spectral gradient (for heatmaps)
    spectral = c(
      "#5E4FA2", "#3288BD", "#66C2A5", "#ABDDA4", "#E6F598",
      "#FFFFBF", "#FEE08B", "#FDAE61", "#F46D43", "#D53E4F", "#9E0142"
    ),

    # DEG 3-class (volcano plots)
    deg = c(
      down = "#3288BD",
      nosig = "#E3E3E3",
      up = "#D53E4F"
    ),

    # Sample groups
    group = c(
      control = "#3288BD",
      treatment = "#D53E4F",
      treat = "#D53E4F"
    ),

    # Gradient for bar plots (6-color)
    gradient_blue_green = c("#66C2A5", "#41A98E", "#3288BD", "#3a6ea5", "#5E4FA2", "#4A3A8C"),

    # Qualitative palette (9-color, for categories)
    qualitative = c(
      "#41A98E", "#ED6355", "#EFA63A", "#3a6ea5",
      "#06d6a0", "#C5DC6C", "#71C8E1", "#d4a373", "#B5739D"
    )
  )
}

#' Get Color Palette
#'
#' @param palette_name Character. Name of palette
#' @param n Integer. Number of colors needed (will interpolate if needed)
#'
#' @return Character vector of colors
#' @keywords internal
.gpdb_get_palette <- function(palette_name = "spectral", n = NULL) {
  palettes <- .gpdb_colors()

  if (!palette_name %in% names(palettes)) {
    palette_name <- "spectral"
  }

  colors <- palettes[[palette_name]]

  # Interpolate if n is specified and different from palette length
  if (!is.null(n) && n != length(colors)) {
    colors <- grDevices::colorRampPalette(colors)(n)
  }

  return(colors)
}
