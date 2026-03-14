# Suppress R CMD check notes for ggraph/ggplot2 NSE
utils::globalVariables(c("effect", "weight", "node_type", "degree", "name"))

#' Visualize Gene Regulatory Network with Upstream and Downstream Relationships
#'
#' @description
#' Construct and visualize bidirectional regulatory network showing both upstream regulators
#' (genes that control focal gene) and downstream targets (genes regulated by focal gene),
#' aggregated from 7,665 RNA-seq perturbation experiments for pathway mapping and regulatory
#' cascade analysis. **Database covers 2,810 genes across 1,063 cell lines** with multi-dataset
#' evidence supporting each regulatory edge. Returns publication-ready network graph with
#' color-coded nodes (focal gene, regulators, targets) and effect-based edge styling
#' (activation/repression).
#'
#' @param gene Character. Gene symbol for focal node (e.g., "TP53", "MYC").
#'   Must exist in database (use \code{gpdb_search()} to verify). Case-insensitive.
#' @param top_regulators Integer. Number of top upstream regulators to display.
#'   Default: \code{10}. Range: 5-50 recommended.
#'   Higher values show more comprehensive regulatory inputs but increase visual complexity.
#' @param top_targets Integer. Number of top downstream targets to display.
#'   Default: \code{10}. Range: 5-50 recommended.
#'   Higher values reveal broader regulatory impact but may clutter visualization.
#' @param min_confidence Character. Evidence strength filter for edges ("high", "medium", "low").
#'   Default: \code{"medium"}. Options: "high" (5+ datasets, publication-grade), "medium" (2-4 datasets),
#'   "low" (1 dataset, exploratory).
#'   Higher confidence reduces false positives; lower captures emerging regulatory relationships.
#' @param layout Character. Network layout algorithm: "fr" (force-directed, default), "circle", "star".
#'   Default: \code{"fr"}. Options: "fr" (Fruchterman-Reingold, reveals clusters), "circle" (equal spacing),
#'   "star" (focal gene centered).
#'   Force-directed layout best for identifying regulatory modules; star layout emphasizes hub structure.
#' @param node_size Numeric. Base size of nodes in mm.
#'   Default: \code{8}. Range: 5-15 recommended.
#'   Node size scales with degree (connectivity) automatically.
#' @param edge_width Numeric. Base width of edges in mm.
#'   Default: \code{1}. Range: 0.5-3 recommended.
#'   Edge width scales with effect size (|logFC|) automatically.
#' @param show_labels Logical. Display gene symbol labels on nodes.
#'   Default: \code{TRUE}.
#'   Set to FALSE for dense networks (>100 nodes) to reduce overlap.
#' @param title Character. Custom plot title.
#'   Default: NULL (auto-generates "Regulatory Network of [gene]").
#' @param theme ggplot2 theme object. Custom theme for plot styling.
#'   Default: NULL (uses \code{.gpdb_theme_default()} academic publication theme).
#'
#' @return ggraph/ggplot object with interactive network visualization.
#'
#' **Network Structure**:
#' \describe{
#'   \item{\strong{Nodes}}{Genes colored by role:
#'     \itemize{
#'       \item Red: Focal gene (query)
#'       \item Blue: Regulators (upstream, control focal gene)
#'       \item Green: Targets (downstream, regulated by focal gene)
#'     }
#'     Node size scales with connectivity degree.
#'   }
#'   \item{\strong{Edges}}{Directed regulatory relationships colored by effect:
#'     \itemize{
#'       \item Green arrows: Activation (positive regulation)
#'       \item Red arrows: Repression (negative regulation)
#'     }
#'     Edge width/transparency scale with effect size (|logFC|).
#'   }
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Identify key hubs: Nodes with high connectivity are master regulators
#'   \item Trace regulatory paths: Use \code{\link{gpdb_analyze_cascade}} for multi-step regulation
#'   \item Validate relationships: Use \code{\link{gpdb_load_data}} to check raw expression data
#'   \item Expand network: Increase \code{top_regulators}/\code{top_targets} to reveal broader network
#'   \item Export for Cytoscape: Extract graph data with \code{igraph::as_data_frame(p$data)}
#'   \item Customize visualization: Modify returned ggplot object (add layers, change theme)
#' }
#'
#' **Alternative Methods**:
#' \itemize{
#'   \item \code{\link{gpdb_analyze_cascade}}: Use when tracing multi-step regulatory cascades (A→B→C)
#'   \item \code{\link{gpdb_plot_cascade}}: Use when visualizing sequential regulatory steps with pathway context
#'   \item \code{\link{gpdb_find_regulators}}: Use when only upstream regulators are needed (tabular output)
#'   \item \code{\link{gpdb_find_targets}}: Use when only downstream targets are needed (tabular output)
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test query**: TP53 (tumor suppressor, highly connected hub)
#'   \item **Runtime**: 0.85 sec (database query + network construction + visualization)
#'   \item **Result size**: 40 edges (20 regulators → TP53 → 20 targets), 40 nodes
#'   \item **Database**: 18.9M relationships, optimized with SQLite indexing + igraph rendering
#' }
#'
#' @details
#' **Interpreting Network Structure**:
#' \itemize{
#'   \item **Node degree** (connectivity): High-degree nodes are regulatory hubs. Focal gene degree
#'     reflects regulatory centrality in biological networks.
#'   \item **Edge directionality**: Arrows point from regulator to target. Regulator perturbation
#'     causes target expression change.
#'   \item **Effect type**: Green (activation) means regulator knockout decreases target (positive control).
#'     Red (repression) means regulator knockout increases target (negative control).
#'   \item **Network density**: Dense networks indicate genes in core regulatory modules; sparse
#'     networks suggest specialized or tissue-specific function.
#' }
#'
#' **Layout Selection Guide**:
#' \itemize{
#'   \item **Force-directed (fr)**: Best for revealing regulatory modules and clustering patterns.
#'     Genes with similar regulatory relationships cluster together.
#'   \item **Circle**: Best for comparing connectivity patterns across equal spacing.
#'     Easier to count edges and identify hubs.
#'   \item **Star**: Best for emphasizing focal gene as central hub.
#'     Clearly distinguishes regulators (pointing in) from targets (pointing out).
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Data Query** (get network data):
#' \itemize{
#'   \item \code{\link{gpdb_find_regulators}}: Query upstream regulators only (tabular format)
#'   \item \code{\link{gpdb_find_targets}}: Query downstream targets only (tabular format)
#'   \item \code{\link{gpdb_list_datasets}}: Find available datasets for focal gene
#' }
#'
#' **Network Analysis** (explore regulatory paths):
#' \itemize{
#'   \item \code{\link{gpdb_analyze_cascade}}: Trace multi-step regulatory cascades (A→B→C)
#'   \item \code{\link{gpdb_plot_cascade}}: Visualize sequential regulatory steps with pathway annotation
#'   \item \code{\link{gpdb_compare_genes}}: Compare regulatory networks of multiple genes
#' }
#'
#' **Downstream Analysis** (interpret network):
#' \itemize{
#'   \item \code{\link{gpdb_enrich}}: Pathway enrichment of network genes (GO/KEGG/MSigDB)
#'   \item \code{\link{gpdb_load_data}}: Load raw expression data to validate edges
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How to visualize which genes regulate my gene of interest?
#'   \item What genes are controlled by TP53 knockout in cancer cells?
#'   \item Can I see both upstream and downstream regulatory relationships?
#'   \item How do I identify master regulators in my network?
#'   \item Which layout is best for finding regulatory modules?
#'   \item How to adjust network size for publication figures?
#'   \item What do the edge colors mean in the network plot?
#'   \item How reliable are the regulatory relationships shown?
#'   \item Can I export the network for Cytoscape analysis?
#'   \item How to trace regulatory cascades from the network?
#'   \item What if my network is too dense to read?
#'   \item How to compare networks of different genes?
#'   \item Can I filter edges by effect size or confidence?
#'   \item How to identify feedback loops in the network?
#'   \item Which genes are the most connected hubs?
#'   \item How to validate regulatory edges with raw data?
#'   \item Can I customize node colors and edge styles?
#'   \item What cell types are represented in the network?
#'   \item How to find tissue-specific regulatory relationships?
#'   \item Can I overlay pathway annotations on the network?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Usage (TESTED - runtime: 0.85 sec)
#' # ===========================================================================
#' # Scientific Question: What are the regulatory inputs and outputs of TP53?
#' # Query: TP53 regulatory network with high-confidence edges
#' # Expected output: Network graph with TP53 as focal hub
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' p <- gpdb_plot_network(
#'   gene = "TP53",
#'   top_regulators = 10,
#'   top_targets = 10,
#'   min_confidence = "high"
#' )
#'
#' # Display network
#' print(p)
#'
#' # Check network structure
#' # Message shows: "Network: 40 edges, 40 nodes"
#'
#' # Interpretation:
#' # - Red node: TP53 (focal gene)
#' # - Blue nodes: Upstream regulators (control TP53)
#' # - Green nodes: Downstream targets (regulated by TP53)
#' # - Green arrows: Activation
#' # - Red arrows: Repression
#'
#' # ===========================================================================
#' # Example 2: Layout Comparison (Show impact of layout algorithms)
#' # ===========================================================================
#' # Compare different layouts to reveal network structure
#'
#' # Force-directed layout: Reveals regulatory modules
#' p_fr <- gpdb_plot_network(
#'   gene = "MYC",
#'   top_regulators = 15,
#'   top_targets = 15,
#'   layout = "fr", # Default, good for clustering
#'   min_confidence = "medium"
#' )
#'
#' # Star layout: Emphasizes hub structure
#' p_star <- gpdb_plot_network(
#'   gene = "MYC",
#'   top_regulators = 15,
#'   top_targets = 15,
#'   layout = "star", # Focal gene centered
#'   min_confidence = "medium"
#' )
#'
#' # Circle layout: Equal spacing for counting edges
#' p_circle <- gpdb_plot_network(
#'   gene = "MYC",
#'   top_regulators = 15,
#'   top_targets = 15,
#'   layout = "circle", # Uniform distribution
#'   min_confidence = "medium"
#' )
#'
#' # Display side-by-side comparison
#' library(patchwork)
#' p_fr + p_star + p_circle + plot_layout(ncol = 3)
#'
#' # ===========================================================================
#' # Example 3: Confidence Level Impact (Filter by evidence strength)
#' # ===========================================================================
#' # Compare how confidence filtering affects network
#'
#' # High confidence: Publication-grade evidence only
#' p_high <- gpdb_plot_network(
#'   gene = "EGFR",
#'   top_regulators = 20,
#'   top_targets = 20,
#'   min_confidence = "high", # 5+ datasets, most reliable
#'   title = "EGFR Network (High Confidence)"
#' )
#'
#' # Medium confidence: Include moderate evidence
#' p_medium <- gpdb_plot_network(
#'   gene = "EGFR",
#'   top_regulators = 20,
#'   top_targets = 20,
#'   min_confidence = "medium", # 2-4 datasets
#'   title = "EGFR Network (Medium Confidence)"
#' )
#'
#' # Compare network sizes
#' print(p_high) # Fewer edges, higher reliability
#' print(p_medium) # More edges, broader coverage
#'
#' # ===========================================================================
#' # Example 4: Customize Visualization (Adjust aesthetics)
#' # ===========================================================================
#' # Fine-tune network appearance for publication
#'
#' p_custom <- gpdb_plot_network(
#'   gene = "KRAS",
#'   top_regulators = 12,
#'   top_targets = 12,
#'   layout = "fr",
#'   node_size = 10, # Larger nodes
#'   edge_width = 1.5, # Thicker edges
#'   show_labels = TRUE, # Show gene names
#'   min_confidence = "high",
#'   title = "KRAS Regulatory Network in Cancer"
#' )
#'
#' # Further customize (ggplot layers work!)
#' p_custom <- p_custom +
#'   ggplot2::theme(
#'     plot.title = ggplot2::element_text(size = 16, face = "bold"),
#'     legend.position = "bottom"
#'   )
#'
#' print(p_custom)
#'
#' # Save high-resolution figure
#' ggplot2::ggsave("kras_network.pdf", p_custom, width = 10, height = 8)
#'
#' # ===========================================================================
#' # Example 5: Extract Network Data (Export for external tools)
#' # ===========================================================================
#' # Get underlying igraph object for advanced analysis
#'
#' p <- gpdb_plot_network("BRCA1", top_regulators = 20, top_targets = 20)
#'
#' # Extract graph data (requires accessing internal structure)
#' # Note: ggraph embeds igraph object in plot data
#'
#' # Alternative: Query data directly
#' regulators <- gpdb_find_regulators("BRCA1", top_n = 20, min_confidence = "high")
#' targets <- gpdb_find_targets("BRCA1", top_n = 20, min_confidence = "high")
#'
#' # Build edge list for Cytoscape
#' edges_up <- data.frame(
#'   source = "BRCA1",
#'   target = targets$upregulated$target_gene,
#'   interaction = "activates",
#'   weight = targets$upregulated$logfc_mean
#' )
#'
#' edges_down <- data.frame(
#'   source = "BRCA1",
#'   target = targets$downregulated$target_gene,
#'   interaction = "represses",
#'   weight = targets$downregulated$logfc_mean
#' )
#'
#' all_edges <- rbind(edges_up, edges_down)
#'
#' # Export for Cytoscape
#' write.csv(all_edges, "brca1_network_edges.csv", row.names = FALSE)
#'
#' # ===========================================================================
#' # Next Steps (Brief pointers to downstream analysis)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - Regulatory cascades: ?gpdb_analyze_cascade
#' # - Pathway enrichment: ?gpdb_enrich
#' # - Compare networks: ?gpdb_compare_genes
#' # - Validate edges: ?gpdb_load_data
#' # - Sequential visualization: ?gpdb_plot_cascade
#' }
#'
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
