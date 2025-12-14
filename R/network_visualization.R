#' Plot Gene Regulatory Network
#'
#' Visualize regulatory relationships as a network graph.
#' Shows which genes regulate/are regulated by a focal gene.
#'
#' @param gene Character. Central gene of interest
#' @param top_regulators Integer. Number of top regulators to show (default 10)
#' @param top_targets Integer. Number of top targets to show (default 10)
#' @param min_confidence Character. Minimum evidence level ("high", "medium", "low")
#' @param layout Character. Network layout: "fr" (force-directed), "circle", "star"
#' @param node_size Numeric. Size of nodes (default 8)
#' @param edge_width Numeric. Width of edges (default 1)
#' @param show_labels Logical. Show gene labels (default TRUE)
#' @param title Character. Plot title
#' @param theme ggplot2 theme object (NULL uses default)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Network of MYC regulators and targets
#' p <- gpdb_plot_network("MYC", top_regulators = 10, top_targets = 10)
#' print(p)
#'
#' # Star layout
#' p2 <- gpdb_plot_network("TP53", layout = "star")
#' }
gpdb_plot_network <- function(gene,
                              top_regulators = 10,
                              top_targets = 10,
                              min_confidence = "medium",
                              layout = c("fr", "circle", "star"),
                              node_size = 8,
                              edge_width = 1,
                              show_labels = TRUE,
                              title = NULL,
                              theme = NULL) {
  if (!requireNamespace("ggraph", quietly = TRUE)) {
    stop("ggraph package required. Install with: install.packages('ggraph')",
      call. = FALSE
    )
  }

  if (!requireNamespace("igraph", quietly = TRUE)) {
    stop("igraph package required. Install with: install.packages('igraph')",
      call. = FALSE
    )
  }

  layout <- match.arg(layout)
  gene <- .gpdb_format_genes(gene)[1]

  # Get regulators (genes that regulate our focal gene)
  regulators <- gpdb_find_regulators(
    gene,
    top_n = top_regulators,
    min_confidence = min_confidence
  )

  # Get targets (genes regulated by our focal gene)
  targets <- gpdb_find_targets(
    gene,
    top_n = top_targets,
    min_confidence = min_confidence
  )

  # Build edge list
  edges <- data.frame()

  # Regulators -> focal gene
  if (!is.null(regulators$repressors) && nrow(regulators$repressors) > 0) {
    rep_edges <- data.frame(
      from = regulators$repressors$perturbed_gene,
      to = gene,
      weight = abs(regulators$repressors$logfc_mean),
      type = "repressor",
      effect = "negative"
    )
    edges <- rbind(edges, head(rep_edges, top_regulators))
  }

  if (!is.null(regulators$activators) && nrow(regulators$activators) > 0) {
    act_edges <- data.frame(
      from = regulators$activators$perturbed_gene,
      to = gene,
      weight = abs(regulators$activators$logfc_mean),
      type = "activator",
      effect = "positive"
    )
    edges <- rbind(edges, head(act_edges, top_regulators))
  }

  # Focal gene -> targets
  if (!is.null(targets$upregulated) && nrow(targets$upregulated) > 0) {
    up_edges <- data.frame(
      from = gene,
      to = targets$upregulated$target_gene,
      weight = abs(targets$upregulated$logfc_mean),
      type = "upregulated",
      effect = "positive"
    )
    edges <- rbind(edges, head(up_edges, top_targets))
  }

  if (!is.null(targets$downregulated) && nrow(targets$downregulated) > 0) {
    down_edges <- data.frame(
      from = gene,
      to = targets$downregulated$target_gene,
      weight = abs(targets$downregulated$logfc_mean),
      type = "downregulated",
      effect = "negative"
    )
    edges <- rbind(edges, head(down_edges, top_targets))
  }

  if (nrow(edges) == 0) {
    stop("No regulatory relationships found for ", gene, call. = FALSE)
  }

  # Create graph
  graph <- igraph::graph_from_data_frame(edges, directed = TRUE)

  # Add node attributes
  all_nodes <- unique(c(edges$from, edges$to))
  node_types <- ifelse(all_nodes == gene, "focal",
    ifelse(all_nodes %in% edges$from[edges$to == gene], "regulator", "target")
  )

  igraph::V(graph)$node_type <- node_types[match(igraph::V(graph)$name, all_nodes)]
  igraph::V(graph)$degree <- igraph::degree(graph)

  # Default title
  if (is.null(title)) {
    title <- paste0("Regulatory Network of ", gene)
  }

  # Default theme
  if (is.null(theme)) {
    theme <- .gpdb_theme_default()
  }

  # Colors
  node_colors <- c(
    focal = "#D53E4F", # Red for focal gene
    regulator = "#3288BD", # Blue for regulators
    target = "#66C2A5" # Green for targets
  )

  edge_colors <- c(
    positive = "#41A98E", # Green for activation
    negative = "#ED6355" # Red for repression
  )

  # Create network plot with ggraph
  p <- ggraph::ggraph(graph, layout = layout) +
    # Edges
    ggraph::geom_edge_link(
      ggplot2::aes(
        color = effect,
        width = weight,
        alpha = weight
      ),
      arrow = ggplot2::arrow(length = ggplot2::unit(3, "mm"), type = "closed"),
      end_cap = ggraph::circle(5, "mm")
    ) +
    ggraph::scale_edge_color_manual(
      values = edge_colors,
      name = "Effect",
      guide = "none"
    ) +
    ggraph::scale_edge_width(
      range = c(0.5, 2) * edge_width,
      guide = "none"
    ) +
    ggraph::scale_edge_alpha(
      range = c(0.3, 0.9),
      guide = "none"
    ) +
    # Nodes
    ggraph::geom_node_point(
      ggplot2::aes(
        color = node_type,
        size = degree
      )
    ) +
    ggplot2::scale_color_manual(
      values = node_colors,
      name = "Node Type",
      labels = c(
        focal = "Focal Gene",
        regulator = "Regulator",
        target = "Target"
      )
    ) +
    ggplot2::scale_size(
      range = c(node_size * 0.8, node_size * 1.5),
      guide = "none"
    ) +
    # Labels
    {
      if (show_labels) {
        ggraph::geom_node_text(
          ggplot2::aes(label = name),
          repel = TRUE,
          size = 3.5,
          fontface = "bold",
          bg.color = "white",
          bg.r = 0.1
        )
      }
    } +
    ggplot2::labs(title = title) +
    theme +
    ggplot2::theme(
      axis.text = ggplot2::element_blank(),
      axis.title = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      axis.line = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank()
    ) +
    ggplot2::coord_equal()

  message("Network: ", nrow(edges), " edges, ", length(all_nodes), " nodes")

  return(p)
}
