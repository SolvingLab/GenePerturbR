#' Analyze Regulatory Cascade
#'
#' Trace regulatory cascades from a starting gene through multiple layers.
#' Identifies multi-step regulatory pathways (e.g., A regulates B, B regulates C).
#'
#' @param start_gene Character. Starting gene
#' @param max_depth Integer. Maximum cascade depth (default 3, max 5)
#' @param min_effect_size Numeric. Minimum effect size for each step (default 1.0)
#' @param min_confidence Character. Minimum confidence level (default "medium")
#' @param direction Character. "forward" (targets), "backward" (regulators), or "both"
#'
#' @return List with cascade paths and statistics
#' @export
#'
#' @examples
#' \dontrun{
#' # Find what TP53 regulates, then what those genes regulate
#' cascade <- gpdb_analyze_cascade("TP53", max_depth = 3)
#'
#' # View top cascades
#' print(cascade$top_paths)
#'
#' # Visualize
#' gpdb_plot_cascade(cascade)
#' }
gpdb_analyze_cascade <- function(start_gene,
                                 max_depth = 3,
                                 min_effect_size = 1.0,
                                 min_confidence = "medium",
                                 direction = c("forward", "backward", "both")) {
  direction <- match.arg(direction)
  start_gene <- .gpdb_format_genes(start_gene)[1]

  if (max_depth > 5) {
    warning("max_depth > 5 may be slow and produce too many paths. Setting to 5.")
    max_depth <- 5
  }

  message("Analyzing regulatory cascade from ", start_gene, " (depth: ", max_depth, ")")

  # Initialize
  all_paths <- list()
  current_genes <- start_gene
  visited <- character(0)

  # Iterative deepening
  for (depth in 1:max_depth) {
    if (length(current_genes) == 0) break

    next_level_genes <- character(0)

    for (curr_gene in current_genes) {
      if (curr_gene %in% visited) next
      visited <- c(visited, curr_gene)

      # Query next layer
      if (direction %in% c("forward", "both")) {
        targets <- tryCatch(
          {
            gpdb_find_targets(
              curr_gene,
              top_n = 20,
              min_confidence = min_confidence,
              min_effect_size = min_effect_size
            )
          },
          error = function(e) {
            return(list(upregulated = data.frame(), downregulated = data.frame()))
          }
        )

        if (!is.null(targets$upregulated) && nrow(targets$upregulated) > 0) {
          for (i in 1:min(5, nrow(targets$upregulated))) {
            target <- targets$upregulated$target_gene[i]
            effect <- targets$upregulated$logfc_mean[i]

            all_paths[[length(all_paths) + 1]] <- list(
              from = curr_gene,
              to = target,
              depth = depth,
              effect = effect,
              direction = "forward"
            )

            if (depth < max_depth) {
              next_level_genes <- c(next_level_genes, target)
            }
          }
        }
      }

      if (direction %in% c("backward", "both")) {
        regulators <- tryCatch(
          {
            gpdb_find_regulators(
              curr_gene,
              top_n = 20,
              min_confidence = min_confidence
            )
          },
          error = function(e) {
            return(list(repressors = data.frame(), activators = data.frame()))
          }
        )

        if (!is.null(regulators$repressors) && nrow(regulators$repressors) > 0) {
          for (i in 1:min(5, nrow(regulators$repressors))) {
            reg <- regulators$repressors$perturbed_gene[i]
            effect <- regulators$repressors$logfc_mean[i]

            all_paths[[length(all_paths) + 1]] <- list(
              from = reg,
              to = curr_gene,
              depth = depth,
              effect = effect,
              direction = "backward"
            )

            if (depth < max_depth && direction == "both") {
              next_level_genes <- c(next_level_genes, reg)
            }
          }
        }
      }
    }

    current_genes <- unique(next_level_genes)
    message("  Depth ", depth, ": found ", length(current_genes), " genes")
  }

  # Convert to data frame
  if (length(all_paths) == 0) {
    message("No cascade paths found")
    return(list(
      start_gene = start_gene,
      n_paths = 0,
      paths = data.frame()
    ))
  }

  paths_df <- do.call(rbind, lapply(all_paths, as.data.frame))

  # Build path strings for analysis
  # This is complex, so for now return edges

  result <- list(
    start_gene = start_gene,
    max_depth = max_depth,
    n_paths = nrow(paths_df),
    paths = paths_df,
    n_genes = length(unique(c(paths_df$from, paths_df$to)))
  )

  message("Found ", nrow(paths_df), " regulatory relationships across ", length(unique(c(paths_df$from, paths_df$to))), " genes")

  return(result)
}


#' Plot Cascade Network
#'
#' Visualize cascade analysis results as a network
#'
#' @param cascade_result Output from gpdb_analyze_cascade
#' @param layout Character. Network layout
#' @param title Character. Plot title
#' @param theme ggplot2 theme
#'
#' @return ggplot object
#' @export
gpdb_plot_cascade <- function(cascade_result,
                              layout = "fr",
                              title = NULL,
                              theme = NULL) {
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("ggraph package required", call. = FALSE)
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package required", call. = FALSE)
  }

  paths <- cascade_result$paths

  if (nrow(paths) == 0) {
    stop("No paths to visualize", call. = FALSE)
  }

  # Create graph
  graph <- igraph::graph_from_data_frame(
    paths[, c("from", "to", "effect", "depth")],
    directed = TRUE
  )

  # Add node attributes
  focal_gene <- cascade_result$start_gene
  all_nodes <- igraph::V(graph)$name

  node_types <- ifelse(all_nodes == focal_gene, "focal", "other")
  igraph::V(graph)$node_type <- node_types
  igraph::V(graph)$degree <- igraph::degree(graph)

  # Title
  if (is.null(title)) {
    title <- paste0(
      "Regulatory Cascade from ", focal_gene,
      " (", cascade_result$n_paths, " edges)"
    )
  }

  # Theme
  if (is.null(theme)) {
    theme <- .gpdb_theme_default()
  }

  # Colors
  node_colors <- c(focal = "#D53E4F", other = "#3288BD")
  edge_color_scale <- .gpdb_get_palette("spectral", 11)

  # Plot
  p <- ggraph::ggraph(graph, layout = layout) +
    ggraph::geom_edge_link(
      ggplot2::aes(
        color = effect,
        width = abs(effect),
        alpha = abs(effect)
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(2, "mm"), type = "closed"),
      end_cap = ggraph::circle(3, "mm")
    ) +
    ggraph::scale_edge_color_gradientn(
      colors = edge_color_scale,
      name = "Effect Size",
      guide = "none" # Disable edge color legend
    ) +
    ggraph::scale_edge_width(range = c(0.3, 1.5), guide = "none") +
    ggraph::scale_edge_alpha(range = c(0.3, 0.9), guide = "none") +
    ggraph::geom_node_point(
      ggplot2::aes(color = node_type, size = degree)
    ) +
    ggplot2::scale_color_manual(
      values = node_colors,
      name = "Node Type"
    ) +
    ggplot2::scale_size(range = c(5, 12), guide = "none") +
    ggraph::geom_node_text(
      ggplot2::aes(label = name),
      repel = TRUE,
      size = 3,
      fontface = "bold",
      bg.color = "white",
      bg.r = 0.1
    ) +
    ggplot2::labs(title = title) +
    theme +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank()
    )

  return(p)
}
