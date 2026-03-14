#' Visualize Differential Expression with Volcano Plot
#'
#' Create publication-ready volcano plot showing gene expression changes from a
#' single perturbation experiment. Displays log2 fold change vs -log10(adjusted p-value)
#' with automatic gene labeling and classification. **Visualizes DEG results from one dataset**
#' (17,000+ genes per experiment) for identifying most significantly changed genes.
#'
#' @param dataset_id Character. Dataset ID from GPSAdb (e.g., "D28200").
#'   Default: None (required). Use \code{gpdb_list_datasets()} to find available datasets.
#'   Each dataset represents one perturbation experiment (gene knockout/knockdown).
#' @param padj_cutoff Numeric. Adjusted p-value significance threshold.
#'   Default: \code{0.05}. Range: 0-1. Genes below this threshold are considered significant.
#'   Lower values (0.01) give higher confidence; higher values (0.1) include more candidates.
#' @param logfc_cutoff Numeric. Log2 fold change threshold for biological significance.
#'   Default: \code{1} (2-fold change). Range: 0-5. Higher values focus on strong effects.
#'   Common: 1 (2-fold), 1.5 (2.8-fold), 2 (4-fold).
#' @param nlabel Integer. Number of top significant genes to auto-label.
#'   Default: \code{10}. Range: 0-50. Set to 0 to disable auto-labeling.
#'   Genes are ranked by adjusted p-value.
#' @param highlight_genes Character vector. Specific genes to label (overrides nlabel).
#'   Default: \code{NULL} (use auto-labeling). Provide gene symbols to highlight genes of interest.
#'   Example: \code{c("TP53", "MYC", "KRAS")}.
#' @param colors Character vector of length 3. Custom colors for Down/NoSig/Up genes.
#'   Default: \code{NULL} (uses publication-friendly blue/grey/red palette).
#'   Order: c(Down, NoSig, Up). Example: \code{c("#3288BD", "#E3E3E3", "#D53E4F")}.
#' @param point.maxsize Numeric. Maximum point size for visualizing effect magnitude.
#'   Default: \code{4}. Range: 1-10. Points are sized by absolute log2 fold change.
#' @param point.alpha Numeric. Transparency of data points.
#'   Default: \code{0.8}. Range: 0-1. Lower values reduce overplotting for dense data.
#' @param intercept.width Numeric. Line width for significance thresholds.
#'   Default: \code{0.65}. Range: 0.3-1.5. Dashed lines mark padj and logFC cutoffs.
#' @param label.size Numeric. Font size for gene labels.
#'   Default: \code{3.5}. Range: 2-6. Larger values for presentation slides.
#' @param label.bg Character. Background color for gene labels.
#'   Default: \code{"white"}. Improves readability over dense point clouds.
#' @param label.bg.r Numeric. Corner radius for label background.
#'   Default: \code{0.1}. Range: 0-0.3. Higher values give rounded corners.
#' @param legend.position Character. Legend placement.
#'   Default: \code{"bottom"}. Options: "top", "bottom", "left", "right", "none".
#' @param title Character. Plot title.
#'   Default: \code{NULL} (auto-generated: "Gene in CellLine (Method)").
#'   Provide custom title to override.
#' @param theme.plot ggplot2 theme object. Custom theme for plot styling.
#'   Default: \code{NULL} (uses \code{.gpdb_theme_default()}, publication-ready academic theme).
#'   Can provide custom ggplot2 theme.
#'
#' @return ggplot object showing differential expression volcano plot.
#'
#' **Plot Elements**:
#' \describe{
#'   \item{X-axis}{Log2 fold change (positive = upregulated, negative = downregulated)}
#'   \item{Y-axis}{-Log10(adjusted p-value). Higher values = more significant}
#'   \item{Points}{Each point = one gene. Color shows regulation status (Up/Down/NoSig)}
#'   \item{Size}{Point size reflects effect magnitude (absolute log2FC)}
#'   \item{Dashed lines}{Significance thresholds (vertical = logFC, horizontal = padj)}
#'   \item{Labels}{Top significant genes or user-specified genes}
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Save plot: \code{ggplot2::ggsave("volcano.pdf", p, width = 8, height = 6)}
#'   \item Customize colors: Modify \code{colors} parameter with custom palette
#'   \item Extract DEG data: Use \code{gpdb_load_deg(dataset_id)} to get underlying data
#'   \item Compare datasets: Plot multiple datasets for different cell lines/conditions
#'   \item Enrichment analysis: Feed significant genes to \code{gpdb_enrich()}
#'   \item Network visualization: Use \code{gpdb_plot_network()} for regulatory relationships
#' }
#'
#' **Alternative Visualizations**:
#' \itemize{
#'   \item \code{\link{gpdb_plot_heatmap_deg}}: Heatmap of DEG across multiple datasets.
#'     Use when: Comparing expression patterns across experiments.
#'   \item \code{\link{gpdb_plot_heatmap_expr}}: Expression heatmap with hierarchical clustering.
#'     Use when: Exploring sample-to-sample relationships.
#'   \item \code{\link{gpdb_plot_enrichment}}: Enrichment results visualization.
#'     Use when: Interpreting pathway/GO term enrichment from DEG.
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test dataset**: D28200 (TP53 in PANC-1 cells, siRNA)
#'   \item **Runtime**: 0.61 sec (load DEG data + ggplot rendering)
#'   \item **Data size**: 17,095 genes (264 up, 95 down at default thresholds)
#'   \item **Result**: ggplot object with 5 layers (points, lines, labels)
#' }
#'
#' @details
#' **Interpreting Volcano Plots**:
#' \itemize{
#'   \item **Quadrants**: Upper corners = significant and strong effect. Top center = significant but weak effect.
#'     Sides = strong effect but not significant. Center = no change.
#'   \item **Point color**: Blue (downregulated), Grey (not significant), Red (upregulated)
#'   \item **Point size**: Larger points have stronger effect sizes (higher |logFC|)
#'   \item **Threshold lines**: Dashed lines mark your chosen significance cutoffs
#'   \item **Gene labels**: Automatically selects most significant genes or shows your highlights
#' }
#'
#' **Publication Theme**:
#' \itemize{
#'   \item Clean white background with minimal gridlines
#'   \item Colorblind-friendly default palette (blue-grey-red)
#'   \item Bold axis labels and title for readability
#'   \item Returns ggplot object for further customization
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Data Access** (get underlying data):
#' \itemize{
#'   \item \code{\link{gpdb_list_datasets}}: Find available datasets for a gene
#'   \item \code{\link{gpdb_get_info}}: Get detailed dataset metadata
#'   \item \code{\link{gpdb_load_deg}}: Load DEG results (input for this plot)
#' }
#'
#' **Downstream Analysis** (use plot results):
#' \itemize{
#'   \item \code{\link{gpdb_enrich}}: Pathway enrichment on significant genes
#'   \item \code{\link{gpdb_plot_network}}: Visualize gene regulatory network
#'   \item \code{\link{gpdb_find_targets}}: Query aggregated targets across datasets
#' }
#'
#' **Related Visualizations** (other plot types):
#' \itemize{
#'   \item \code{\link{gpdb_plot_heatmap_deg}}: DEG heatmap across multiple datasets
#'   \item \code{\link{gpdb_plot_heatmap_expr}}: Expression heatmap with clustering
#'   \item \code{\link{gpdb_plot_enrichment}}: Enrichment dotplot/barplot
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I visualize differential expression results from my perturbation experiment?
#'   \item What genes are most significantly changed after TP53 knockout?
#'   \item How do I create a publication-ready volcano plot?
#'   \item Can I highlight specific genes of interest in the volcano plot?
#'   \item How do I adjust significance thresholds for volcano plots?
#'   \item What do the colors and sizes mean in a volcano plot?
#'   \item How do I compare DEG patterns across different datasets visually?
#'   \item Which genes should I focus on for downstream validation?
#'   \item How do I customize volcano plot colors for my presentation?
#'   \item What's the difference between adjusted p-value and fold change thresholds?
#'   \item How do I export high-resolution volcano plots for publication?
#'   \item Can I plot multiple datasets side-by-side for comparison?
#'   \item How do I interpret volcano plot quadrants biologically?
#'   \item What threshold values are recommended for drug target discovery?
#'   \item How do I label more genes or disable auto-labeling?
#'   \item Can I change the theme for dark background presentations?
#'   \item How do I identify outlier genes with extreme fold changes?
#'   \item What's the best way to show subtle but significant changes?
#'   \item How do I overlay pathway-specific genes on volcano plot?
#'   \item Can I use volcano plots for comparing different perturbation methods?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Usage (TESTED - runtime: 0.61 sec)
#' # ===========================================================================
#' # Scientific Question: What genes are significantly changed by TP53 perturbation?
#' # Dataset: TP53 siRNA in PANC-1 cells
#' # Expected output: Volcano plot with 264 up, 95 down (at default thresholds)
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' # Find available datasets
#' datasets <- gpdb_list_datasets(gene = "TP53")
#' cat("Found", nrow(datasets), "datasets for TP53\n")
#'
#' # Create basic volcano plot
#' p1 <- gpdb_plot_volcano(dataset_id = datasets$dataset_id[1])
#'
#' # Plot is ready for viewing or saving
#' print(p1)
#'
#' # Save to file
#' # ggplot2::ggsave("volcano_tp53.pdf", p1, width = 8, height = 6)
#'
#' # ===========================================================================
#' # Example 2: Customize Significance Thresholds
#' # ===========================================================================
#' # Compare different stringency levels
#'
#' # Stringent thresholds (high confidence)
#' p2_stringent <- gpdb_plot_volcano(
#'   dataset_id = datasets$dataset_id[1],
#'   padj_cutoff = 0.01, # More significant
#'   logfc_cutoff = 2, # 4-fold change
#'   nlabel = 20
#' )
#'
#' # Standard thresholds (balanced)
#' p2_standard <- gpdb_plot_volcano(
#'   dataset_id = datasets$dataset_id[1],
#'   padj_cutoff = 0.05, # Standard
#'   logfc_cutoff = 1, # 2-fold change
#'   nlabel = 10
#' )
#'
#' # Lenient thresholds (exploratory)
#' p2_lenient <- gpdb_plot_volcano(
#'   dataset_id = datasets$dataset_id[1],
#'   padj_cutoff = 0.1, # More permissive
#'   logfc_cutoff = 0.5, # 1.4-fold change
#'   nlabel = 5
#' )
#'
#' # ===========================================================================
#' # Example 3: Highlight Specific Genes
#' # ===========================================================================
#' # Focus on genes of interest (e.g., cell cycle genes)
#'
#' # Highlight candidate drug targets
#' cell_cycle_genes <- c("CDK1", "CDK2", "CDK4", "CCND1", "CCNE1", "MYC")
#'
#' p3 <- gpdb_plot_volcano(
#'   dataset_id = datasets$dataset_id[1],
#'   highlight_genes = cell_cycle_genes,
#'   point.maxsize = 5,
#'   label.size = 4
#' )
#'
#' # Verify which genes are present
#' deg_data <- gpdb_load_deg(datasets$dataset_id[1])
#' found_genes <- cell_cycle_genes[cell_cycle_genes %in% deg_data$gene]
#' cat("Highlighted genes found:", paste(found_genes, collapse = ", "), "\n")
#'
#' # ===========================================================================
#' # Example 4: Customize Visual Appearance
#' # ===========================================================================
#' # Publication-ready customization
#'
#' # Custom color palette (colorblind-friendly)
#' custom_colors <- c(
#'   "#0072B2", # Down (blue)
#'   "#999999", # NoSig (grey)
#'   "#D55E00" # Up (orange)
#' )
#'
#' p4 <- gpdb_plot_volcano(
#'   dataset_id = datasets$dataset_id[1],
#'   colors = custom_colors,
#'   point.maxsize = 6,
#'   point.alpha = 0.6,
#'   legend.position = "right",
#'   title = "TP53 Perturbation in Pancreatic Cancer Cells"
#' )
#'
#' # Further customization with ggplot2
#' library(ggplot2)
#' p4_custom <- p4 +
#'   theme(plot.title = element_text(size = 16, face = "bold")) +
#'   labs(caption = "Data source: GPSAdb (D28200)")
#'
#' # ===========================================================================
#' # Next Steps (Workflow Integration)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - Load DEG data: ?gpdb_load_deg
#' # - Enrichment analysis: ?gpdb_enrich
#' # - Network visualization: ?gpdb_plot_network
#' # - Compare across datasets: ?gpdb_plot_heatmap_deg
#' # - Find aggregated targets: ?gpdb_find_targets
#' }
#'
gpdb_plot_volcano <- function(dataset_id,
                              padj_cutoff = 0.05,
                              logfc_cutoff = 1,
                              nlabel = 10,
                              highlight_genes = NULL,
                              colors = NULL,
                              point.maxsize = 4,
                              point.alpha = 0.8,
                              intercept.width = 0.65,
                              label.size = 3.5,
                              label.bg = "white",
                              label.bg.r = 0.1,
                              legend.position = "bottom",
                              title = NULL,
                              theme.plot = NULL) {
  if (!requireNamespace("ggrepel", quietly = TRUE)) {
    stop("ggrepel package required. Install with: install.packages('ggrepel')",
      call. = FALSE
    )
  }

  .gpdb_validate_dataset_id(dataset_id)

  # Load DEG data
  deg_data <- gpdb_load_deg(dataset_id)

  # Get dataset info for title
  if (is.null(title)) {
    info <- gpdb_get_info(dataset_id)
    title <- paste0(info$gene, " in ", info$cell_line, " (", info$method, ")")
  }

  # Filter valid data
  deg_data <- deg_data[!is.na(deg_data$logFC) & !is.na(deg_data$adj.P.Val), ]

  if (nrow(deg_data) == 0) {
    stop("No valid data for plotting", call. = FALSE)
  }

  if (!"gene" %in% names(deg_data)) {
    deg_data$gene <- rownames(deg_data)
  }

  # Classify genes into Up/Down/NoSig
  deg_data$Type <- ifelse(
    deg_data$adj.P.Val < padj_cutoff & deg_data$logFC > logfc_cutoff, "Up",
    ifelse(
      deg_data$adj.P.Val < padj_cutoff & deg_data$logFC < -logfc_cutoff, "Down",
      "NoSig"
    )
  )

  deg_data$Type <- factor(deg_data$Type, levels = c("Down", "NoSig", "Up"))

  # Calculate -log10(P)
  deg_data$neg_log10_p <- -log10(deg_data$adj.P.Val)

  # Select genes to label
  if (!is.null(highlight_genes)) {
    # User-specified genes
    deg_data$label <- ifelse(deg_data$gene %in% highlight_genes, deg_data$gene, NA)
  } else {
    # Auto-select top genes by significance
    deg_sig <- deg_data[deg_data$Type != "NoSig", ]
    if (nrow(deg_sig) > 0) {
      deg_sig <- deg_sig[order(deg_sig$adj.P.Val), ]
      top_genes <- head(deg_sig$gene, nlabel)
      deg_data$label <- ifelse(deg_data$gene %in% top_genes, deg_data$gene, NA)
    } else {
      deg_data$label <- NA
    }
  }

  # Prepare labeled data
  deg_labeled <- deg_data[!is.na(deg_data$label), ]

  # Use default colors if not provided
  if (is.null(colors)) {
    deg_colors <- .gpdb_colors()$deg
    colors <- c(deg_colors["down"], deg_colors["nosig"], deg_colors["up"])
  }
  names(colors) <- c("Down", "NoSig", "Up")

  # Use default theme if not provided
  if (is.null(theme.plot)) {
    theme.plot <- .gpdb_theme_default()
  }

  # Create base plot
  p <- ggplot2::ggplot(deg_data, ggplot2::aes(x = logFC, y = neg_log10_p)) +
    ggplot2::geom_point(
      ggplot2::aes(color = Type, size = abs(logFC)),
      alpha = point.alpha
    ) +
    ggplot2::geom_vline(
      xintercept = c(-logfc_cutoff, logfc_cutoff),
      linetype = "dashed",
      color = "grey30",
      linewidth = intercept.width
    ) +
    ggplot2::geom_hline(
      yintercept = -log10(padj_cutoff),
      linetype = "dashed",
      color = "grey30",
      linewidth = intercept.width
    ) +
    ggplot2::scale_color_manual(
      values = colors,
      labels = c("Down" = "Downregulated", "NoSig" = "Not Significant", "Up" = "Upregulated")
    ) +
    ggplot2::scale_size_area(max_size = point.maxsize) +
    ggplot2::labs(
      x = bquote(~ Log[2] ~ "(Fold Change)"),
      y = bquote(~ -Log[10] ~ "(Adjusted P-value)"),
      title = title,
      color = NULL
    ) +
    ggplot2::xlim(-max(abs(deg_data$logFC)) * 1.05, max(abs(deg_data$logFC)) * 1.05) +
    theme.plot +
    ggplot2::theme(
      axis.text = ggplot2::element_text(size = 10, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 12, colour = "black", face = "bold"),
      axis.title.y = ggplot2::element_text(size = 12, colour = "black", face = "bold"),
      panel.grid.major = ggplot2::element_line(color = "grey90", linewidth = 0.3),
      panel.grid.minor = ggplot2::element_blank(),
      panel.background = ggplot2::element_rect(fill = "white"),
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14, colour = "black", face = "bold"),
      legend.position = legend.position,
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11, colour = "black")
    ) +
    ggplot2::guides(
      size = "none",
      color = ggplot2::guide_legend(
        order = 0,
        override.aes = list(size = point.maxsize, alpha = 1)
      )
    )

  # Add gene labels with white background (plotthis2 style)
  if (nrow(deg_labeled) > 0) {
    # First add circle markers for labeled genes
    p <- p + ggplot2::geom_point(
      data = deg_labeled,
      ggplot2::aes(x = logFC, y = neg_log10_p),
      size = point.maxsize,
      shape = 21,
      stroke = 1,
      color = "black",
      fill = NA
    )

    # Add labels with white background using bg.color and bg.r parameters
    p <- p + ggrepel::geom_text_repel(
      data = deg_labeled,
      ggplot2::aes(x = logFC, y = neg_log10_p, label = label),
      size = label.size,
      color = "black",
      bg.color = label.bg, # Background color
      bg.r = label.bg.r, # Background radius
      box.padding = 0.5,
      point.padding = 0.3,
      segment.color = "grey40",
      max.overlaps = 100,
      min.segment.length = 0,
      seed = 8525
    )
  }

  # Calculate and report summary
  n_up <- sum(deg_data$Type == "Up", na.rm = TRUE)
  n_down <- sum(deg_data$Type == "Down", na.rm = TRUE)
  n_total <- nrow(deg_data)

  message("Volcano plot: ", n_up, " up, ", n_down, " down (of ", n_total, " total genes)")
  if (nrow(deg_labeled) > 0) {
    message("Labeled ", nrow(deg_labeled), " genes")
  }

  return(p)
}


#' Compare Differential Expression Across Multiple Datasets with Heatmap
#'
#' Create publication-ready heatmap comparing gene expression changes (log2FC) across
#' multiple perturbation experiments. Each column represents one dataset, each row represents
#' one gene. **Visualizes cross-dataset DEG patterns** for identifying consistent perturbation
#' effects across cell lines, methods, or conditions.
#'
#' @param dataset_ids Character vector. Dataset IDs from GPSAdb to compare (minimum 2 required).
#'   Default: None (required). Use \code{gpdb_list_datasets()} to find datasets for same gene.
#'   Example: \code{c("D28200", "D28199", "D27994")} for comparing TP53 perturbations.
#' @param genes Character vector. Specific genes to visualize in heatmap.
#'   Default: \code{NULL} (auto-select top DEGs). Provide gene symbols to focus on genes of interest.
#'   Example: \code{c("TP53", "MDM2", "CDKN1A", "BAX")}. Case-insensitive.
#' @param top_n Integer. Number of top DEGs to show when genes not specified.
#'   Default: \code{50}. Range: 10-200. Auto-selects most significantly changed genes across datasets.
#'   Higher values show more genes but may reduce readability for large heatmaps.
#' @param cluster_rows Logical. Perform hierarchical clustering on genes (rows).
#'   Default: \code{TRUE}. Groups genes with similar expression patterns across datasets.
#'   Set \code{FALSE} to preserve gene order or when genes have specific ordering.
#' @param cluster_cols Logical. Perform hierarchical clustering on datasets (columns).
#'   Default: \code{TRUE}. Groups datasets with similar perturbation effects.
#'   Set \code{FALSE} to preserve dataset order (e.g., time series).
#' @param scale Character. Data scaling method for visualization.
#'   Default: \code{"row"} (z-score per gene). Options: "row", "column", "none".
#'   "row" scaling highlights relative changes across datasets; "none" shows raw log2FC values.
#' @param colors Character vector. Custom color palette for heatmap gradient.
#'   Default: \code{NULL} (uses 11-color spectral palette: blue-white-red).
#'   Provide vector of colors for custom gradient. Diverging palettes recommended for log2FC.
#' @param show_values Logical. Display numeric values in heatmap cells.
#'   Default: \code{FALSE}. Set \code{TRUE} for small heatmaps (<20 genes) to show exact values.
#'   Not recommended for large heatmaps due to overcrowding.
#' @param title Character. Custom plot title.
#'   Default: \code{NULL} (auto-generated: "Gene Expression (N datasets)").
#'   Provide custom title to override.
#' @param theme ggplot2 theme object. Custom theme for plot styling.
#'   Default: \code{NULL} (uses \code{.gpdb_theme_default()}, publication-ready academic theme).
#'   Can provide custom ggplot2 theme.
#'
#' @return ggplot object showing DEG heatmap for cross-dataset comparison.
#'
#' **Plot Elements**:
#' \describe{
#'   \item{Rows}{Genes (auto-selected top DEGs or user-specified)}
#'   \item{Columns}{Datasets (each experiment/perturbation)}
#'   \item{Color}{Expression change (blue = downregulated, red = upregulated, white = no change)}
#'   \item{Scaling}{Default z-score normalization per gene for cross-dataset comparability}
#'   \item{Clustering}{Hierarchical clustering groups similar patterns (optional)}
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Save plot: \code{ggplot2::ggsave("heatmap_deg.pdf", p, width = 10, height = 8)}
#'   \item Extract gene clusters: Identify genes with consistent up/down regulation
#'   \item Focus on specific genes: Re-plot with \code{genes} parameter for targets
#'   \item Compare conditions: Use to compare same gene across cell lines or methods
#'   \item Enrichment analysis: Feed consistent DEGs to \code{gpdb_enrich()}
#'   \item Investigate outliers: Identify datasets with unique expression patterns
#'   \item Single dataset detail: Use \code{gpdb_plot_volcano()} for individual dataset exploration
#' }
#'
#' **Alternative Visualizations**:
#' \itemize{
#'   \item \code{\link{gpdb_plot_heatmap_expr}}: Sample-level expression heatmap for ONE dataset.
#'     Use when: Need to see individual sample variation within experiment.
#'   \item \code{\link{gpdb_plot_volcano}}: Volcano plot for ONE dataset.
#'     Use when: Want to identify significant DEGs in single experiment.
#'   \item \code{\link{gpdb_plot_network}}: Gene regulatory network visualization.
#'     Use when: Exploring gene-gene relationships and regulatory interactions.
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test 1 - Auto-selected genes**: 5 datasets (TP53 in different cell lines), top 30 DEGs
#'   \item **Runtime**: 0.55 sec (load DEGs + clustering + ggplot rendering)
#'   \item **Result size**: 30 genes × 5 datasets heatmap
#'   \item **Test 2 - Specific genes**: 5 datasets, 6 specified genes (TP53 pathway)
#'   \item **Runtime**: 0.29 sec (faster with fewer genes)
#'   \item **Result**: ggplot object with 1 layer (geom_tile)
#' }
#'
#' @details
#' **Interpreting Cross-Dataset Heatmaps**:
#' \itemize{
#'   \item **Consistent patterns**: Genes showing same direction (red or blue) across datasets indicate
#'     robust perturbation effects. These are high-confidence targets.
#'   \item **Dataset-specific effects**: Genes showing different patterns may indicate context-dependent
#'     regulation (cell-type specific, method-dependent).
#'   \item **Color intensity**: Stronger colors (darker red/blue) indicate larger fold changes after z-score scaling.
#'   \item **Clustering**: Dendrogram groups reveal gene modules and dataset similarities.
#'   \item **Missing values**: White/grey cells indicate gene not detected or not significant in that dataset.
#' }
#'
#' **When to Use This Heatmap**:
#' \itemize{
#'   \item Compare same gene perturbation across multiple cell lines
#'   \item Identify conserved vs cell-type-specific targets
#'   \item Validate findings across different experimental methods (CRISPR, siRNA, shRNA)
#'   \item Find gene signatures consistent across perturbation contexts
#'   \item Prioritize targets with reproducible effects for validation
#' }
#'
#' **Scaling Methods**:
#' \itemize{
#'   \item **"row" (default)**: Z-score per gene. Best for comparing relative changes across datasets.
#'     Highlights which datasets have strongest effect for each gene.
#'   \item **"column"**: Z-score per dataset. Best for comparing which genes change most within each dataset.
#'   \item **"none"**: Raw log2FC values. Best when absolute fold changes are important.
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Data Access** (get input data):
#' \itemize{
#'   \item \code{\link{gpdb_list_datasets}}: Find available datasets for a gene
#'   \item \code{\link{gpdb_get_info}}: Get detailed metadata for datasets
#'   \item \code{\link{gpdb_load_batch}}: Load multiple DEG datasets efficiently
#' }
#'
#' **Downstream Analysis** (use heatmap results):
#' \itemize{
#'   \item \code{\link{gpdb_enrich}}: Pathway enrichment on consistent DEGs
#'   \item \code{\link{gpdb_find_targets}}: Query aggregated targets across datasets
#'   \item \code{\link{gpdb_compare_genes}}: Statistical comparison of perturbation effects
#' }
#'
#' **Related Visualizations** (complementary plots):
#' \itemize{
#'   \item \code{\link{gpdb_plot_volcano}}: Volcano plot for single dataset detail
#'   \item \code{\link{gpdb_plot_heatmap_expr}}: Sample-level expression heatmap
#'   \item \code{\link{gpdb_plot_network}}: Gene regulatory network
#'   \item \code{\link{gpdb_plot_enrichment}}: Enrichment results visualization
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I compare gene expression changes across multiple experiments?
#'   \item Which genes are consistently affected by TP53 perturbation in different cell lines?
#'   \item How do I visualize DEG patterns across multiple datasets?
#'   \item Can I create a heatmap comparing CRISPR vs siRNA knockdown effects?
#'   \item How do I identify cell-type-specific vs universal targets?
#'   \item What's the best way to visualize perturbation reproducibility?
#'   \item How do I cluster datasets by similar expression patterns?
#'   \item Can I show specific genes of interest across experiments?
#'   \item How do I compare the same perturbation in different tissues?
#'   \item What do the colors mean in cross-dataset heatmaps?
#'   \item How do I interpret z-score scaling in heatmaps?
#'   \item Should I use row or column scaling for my comparison?
#'   \item How do I identify genes with context-dependent regulation?
#'   \item Can I disable clustering to keep my dataset order?
#'   \item How many datasets should I include in heatmap comparison?
#'   \item What's the difference between DEG heatmap and expression heatmap?
#'   \item How do I export high-resolution heatmaps for publication?
#'   \item Can I customize the color palette for my heatmap?
#'   \item How do I handle missing values in cross-dataset comparisons?
#'   \item What size heatmap is optimal for readability?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Cross-Dataset Comparison (TESTED - runtime: 0.55 sec)
#' # ===========================================================================
#' # Scientific Question: Which genes are consistently affected by TP53 perturbation?
#' # Compare: TP53 knockdown/knockout across 5 different cell lines
#' # Expected output: Heatmap showing 30 top DEGs × 5 datasets
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' # Find available TP53 perturbation datasets
#' datasets <- gpdb_list_datasets(gene = "TP53")
#' cat("Found", nrow(datasets), "datasets for TP53\n")
#'
#' # Select first 5 datasets for comparison
#' dataset_ids <- head(datasets$dataset_id, 5)
#'
#' # Create cross-dataset DEG heatmap
#' p1 <- gpdb_plot_heatmap_deg(
#'   dataset_ids = dataset_ids,
#'   top_n = 30, # Show top 30 DEGs
#'   cluster_rows = TRUE, # Cluster genes
#'   cluster_cols = TRUE # Cluster datasets
#' )
#'
#' # Plot is ready
#' print(p1)
#'
#' # Save to file
#' # ggplot2::ggsave("tp53_heatmap.pdf", p1, width = 10, height = 8)
#'
#' # ===========================================================================
#' # Example 2: Focus on Specific Genes (Pathway-Centric)
#' # ===========================================================================
#' # Question: How do TP53 target genes respond across different contexts?
#' # Focus on known TP53 pathway genes
#'
#' # Define TP53 pathway genes
#' tp53_targets <- c(
#'   "CDKN1A", # p21, cell cycle arrest
#'   "MDM2", # Negative regulator
#'   "BAX", # Apoptosis
#'   "PUMA", # Apoptosis (BBC3)
#'   "NOXA", # Apoptosis (PMAIP1)
#'   "GADD45A", # DNA repair
#'   "TP53INP1" # Cell cycle
#' )
#'
#' p2 <- gpdb_plot_heatmap_deg(
#'   dataset_ids = dataset_ids,
#'   genes = tp53_targets, # Specific genes only
#'   cluster_rows = TRUE, # Group by expression pattern
#'   scale = "row" # Z-score per gene
#' )
#'
#' # Identify consistent vs variable targets
#' cat("Heatmap shows which TP53 targets are universally vs context-specifically regulated\n")
#'
#' # ===========================================================================
#' # Example 3: Compare Different Perturbation Methods
#' # ===========================================================================
#' # Question: Are CRISPR and siRNA knockdowns comparable?
#'
#' # Get datasets with different methods
#' datasets_all <- gpdb_list_datasets(gene = "TP53")
#'
#' # Select mix of methods (example)
#' # Note: In real analysis, filter by 'method' column in datasets
#' crispr_ids <- datasets_all$dataset_id[grep("CRISPR", datasets_all$method, ignore.case = TRUE)]
#' sirna_ids <- datasets_all$dataset_id[grep("siRNA", datasets_all$method, ignore.case = TRUE)]
#'
#' # Mix methods for comparison
#' mixed_ids <- c(head(crispr_ids, 3), head(sirna_ids, 3))
#'
#' if (length(mixed_ids) >= 2) {
#'   p3 <- gpdb_plot_heatmap_deg(
#'     dataset_ids = mixed_ids,
#'     top_n = 40,
#'     cluster_cols = TRUE, # Will group by method similarity
#'     title = "TP53 Perturbation: CRISPR vs siRNA"
#'   )
#'
#'   cat("Check if datasets cluster by method (CRISPR vs siRNA)\n")
#' }
#'
#' # ===========================================================================
#' # Example 4: Customize Appearance
#' # ===========================================================================
#' # Publication-ready customization
#'
#' # Custom diverging color palette (RdBu - colorblind friendly)
#' library(ggplot2)
#'
#' p4 <- gpdb_plot_heatmap_deg(
#'   dataset_ids = dataset_ids,
#'   top_n = 25,
#'   scale = "row",
#'   cluster_rows = TRUE,
#'   cluster_cols = FALSE, # Keep dataset order
#'   title = "TP53 Target Gene Expression Across Cell Lines"
#' )
#'
#' # Further customization with ggplot2
#' p4_custom <- p4 +
#'   theme(
#'     axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
#'     axis.text.y = element_text(size = 8),
#'     plot.title = element_text(size = 14, face = "bold")
#'   ) +
#'   labs(caption = "Data source: GPSAdb")
#'
#' # ===========================================================================
#' # Example 5: Compare Same Gene Across Tissues
#' # ===========================================================================
#' # Question: Is TP53 effect tissue-specific?
#'
#' # Get dataset info to filter by tissue
#' all_tp53 <- gpdb_list_datasets(gene = "TP53")
#'
#' # Example: Compare specific tissues (if available)
#' # Note: Filter by tissue_type or cell_line in real analysis
#' tissue_ids <- head(all_tp53$dataset_id, 6)
#'
#' p5 <- gpdb_plot_heatmap_deg(
#'   dataset_ids = tissue_ids,
#'   top_n = 35,
#'   cluster_cols = TRUE, # Cluster by tissue similarity
#'   scale = "row",
#'   title = "TP53 Effects: Tissue Comparison"
#' )
#'
#' cat("Look for tissue-specific gene signatures in clustering\n")
#'
#' # ===========================================================================
#' # Next Steps (Workflow Integration)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - Find datasets: ?gpdb_list_datasets
#' # - Load batch DEGs: ?gpdb_load_batch
#' # - Enrich consistent genes: ?gpdb_enrich
#' # - Single dataset detail: ?gpdb_plot_volcano
#' # - Network analysis: ?gpdb_plot_network
#' }
#'
gpdb_plot_heatmap_deg <- function(dataset_ids,
                                  genes = NULL,
                                  top_n = 50,
                                  cluster_rows = TRUE,
                                  cluster_cols = TRUE,
                                  scale = "row",
                                  colors = NULL,
                                  show_values = FALSE,
                                  title = NULL,
                                  theme = NULL) {
  if (length(dataset_ids) < 2) {
    stop("Need at least 2 datasets for heatmap comparison", call. = FALSE)
  }

  # Load DEG data
  deg_list <- gpdb_load_batch(dataset_ids, type = "deg")

  # Build matrix
  if (!is.null(genes)) {
    genes <- .gpdb_format_genes(genes)
    logfc_matrix <- matrix(NA, nrow = length(genes), ncol = length(deg_list))
    rownames(logfc_matrix) <- genes
    colnames(logfc_matrix) <- names(deg_list)

    for (i in seq_along(deg_list)) {
      deg <- deg_list[[i]]
      for (g in genes) {
        idx <- which(toupper(deg$gene) == g)
        if (length(idx) > 0) {
          logfc_matrix[g, i] <- deg$logFC[idx[1]]
        }
      }
    }
  } else {
    # Top DEGs
    all_genes <- unique(unlist(lapply(deg_list, function(x) {
      x <- x[!is.na(x$adj.P.Val) & x$adj.P.Val < 0.05, ]
      head(x$gene[order(abs(x$logFC), decreasing = TRUE)], top_n)
    })))

    all_genes <- head(all_genes[all_genes != ""], top_n)

    logfc_matrix <- matrix(NA, nrow = length(all_genes), ncol = length(deg_list))
    rownames(logfc_matrix) <- all_genes
    colnames(logfc_matrix) <- names(deg_list)

    for (i in seq_along(deg_list)) {
      deg <- deg_list[[i]]
      for (g in all_genes) {
        idx <- which(deg$gene == g)
        if (length(idx) > 0) {
          logfc_matrix[g, i] <- deg$logFC[idx[1]]
        }
      }
    }
  }

  # Remove all-NA rows
  logfc_matrix <- logfc_matrix[rowSums(!is.na(logfc_matrix)) > 0, , drop = FALSE]

  if (nrow(logfc_matrix) == 0) {
    stop("No overlapping genes found", call. = FALSE)
  }

  if (is.null(title)) {
    title <- paste0("Gene Expression (", ncol(logfc_matrix), " datasets)")
  }

  # Use default colors if not provided
  if (is.null(colors)) {
    colors <- .gpdb_get_palette("spectral", 11)
  }

  # Use default theme if not provided
  if (is.null(theme)) {
    theme <- .gpdb_theme_default()
  }

  # Use ggplot2 heatmap
  legend_title <- if (scale == "row") "Z-score" else "log2FC"

  ht <- .gpdb_ggplot_heatmap(
    logfc_matrix,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    row_labels = rownames(logfc_matrix),
    col_labels = colnames(logfc_matrix),
    title = title,
    colors = colors,
    legend_title = legend_title,
    show_values = show_values,
    theme = theme
  )

  message("Heatmap: ", nrow(logfc_matrix), " genes × ", ncol(logfc_matrix), " datasets")

  return(ht)
}


#' Visualize Sample-Level Expression with Heatmap for Single Dataset
#'
#' Create publication-ready heatmap showing raw expression values across samples within
#' a single perturbation experiment. Each column represents one sample (control or treatment),
#' each row represents one gene. **Visualizes sample-to-sample variation** for top upregulated
#' and downregulated genes with automatic sample group annotations.
#'
#' @param dataset_id Character. Dataset ID from GPSAdb to visualize.
#'   Default: None (required). Use \code{gpdb_list_datasets()} to find available datasets.
#'   Example: \code{"D28200"} for TP53 siRNA in PANC-1 cells.
#' @param top_up Integer. Number of top upregulated genes to show.
#'   Default: \code{20}. Range: 5-50. Genes ranked by log2 fold change (highest first).
#'   Higher values show more upregulated genes but may reduce readability.
#' @param top_down Integer. Number of top downregulated genes to show.
#'   Default: \code{20}. Range: 5-50. Genes ranked by log2 fold change (lowest first).
#'   Total genes displayed = top_up + top_down.
#' @param scale Character. Data scaling method for visualization.
#'   Default: \code{"row"} (z-score per gene). Options: "row", "column", "none".
#'   "row" scaling highlights relative expression across samples; "none" shows raw normalized expression.
#' @param cluster_rows Logical. Perform hierarchical clustering on genes (rows).
#'   Default: \code{TRUE}. Groups genes with similar expression patterns across samples.
#'   Set \code{FALSE} to preserve up/down gene ordering.
#' @param cluster_cols Logical. Perform hierarchical clustering on samples (columns).
#'   Default: \code{TRUE}. Groups samples by expression similarity (may separate control vs treatment).
#'   Set \code{FALSE} to preserve sample order from metadata.
#' @param show_gene_names Logical. Display gene names on y-axis.
#'   Default: \code{TRUE}. Set \code{FALSE} for cleaner appearance when showing many genes (>40).
#' @param colors Character vector. Custom color palette for heatmap gradient.
#'   Default: \code{NULL} (uses 11-color spectral palette: blue-white-red).
#'   Provide vector of colors for custom gradient.
#' @param show_values Logical. Display numeric values in heatmap cells.
#'   Default: \code{FALSE}. Set \code{TRUE} for small heatmaps (<20 genes) to show exact values.
#'   Not recommended for large heatmaps due to overcrowding.
#' @param title Character. Custom plot title.
#'   Default: \code{NULL} (auto-generated: "Gene in CellLine (Method)").
#'   Provide custom title to override.
#' @param theme ggplot2 theme object. Custom theme for plot styling.
#'   Default: \code{NULL} (uses \code{.gpdb_theme_default()}, publication-ready academic theme).
#'   Can provide custom ggplot2 theme.
#'
#' @return ggplot object showing sample-level expression heatmap for single dataset.
#'
#' **Plot Elements**:
#' \describe{
#'   \item{Rows}{Genes (top upregulated + top downregulated from DEG analysis)}
#'   \item{Columns}{Samples (individual replicates from control and treatment groups)}
#'   \item{Color}{Expression level (blue = low, white = medium, red = high)}
#'   \item{Scaling}{Default z-score normalization per gene for visualizing relative changes}
#'   \item{Clustering}{Hierarchical clustering reveals sample groupings and gene modules}
#'   \item{Annotations}{Row annotation shows Up/Down status; column annotation shows sample groups}
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Save plot: \code{ggplot2::ggsave("heatmap_expr.pdf", p, width = 10, height = 12)}
#'   \item Check sample clustering: Verify control and treatment samples cluster separately
#'   \item Identify outliers: Spot samples with unusual expression patterns
#'   \item Extract gene groups: Find co-expressed gene modules from row clustering
#'   \item Validate targets: Confirm top DEG genes show clear treatment effects
#'   \item Quality control: Assess technical variation within groups
#'   \item Compare with DEG: Cross-reference with \code{gpdb_plot_volcano()} for same dataset
#' }
#'
#' **Alternative Visualizations**:
#' \itemize{
#'   \item \code{\link{gpdb_plot_volcano}}: Volcano plot for identifying significant DEGs in same dataset.
#'     Use when: Want to see fold change vs significance for all genes.
#'   \item \code{\link{gpdb_plot_heatmap_deg}}: DEG heatmap across multiple datasets.
#'     Use when: Comparing expression changes (log2FC) across experiments, not within-sample variation.
#'   \item \code{\link{gpdb_plot_network}}: Gene regulatory network visualization.
#'     Use when: Exploring gene-gene regulatory relationships.
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test dataset**: D28200 (TP53 siRNA in PANC-1 cells)
#'   \item **Runtime**: 1.04 sec (load expression + DEG + clustering + ggplot rendering)
#'   \item **Data size**: 29,586 genes total, 6 samples (3 control + 3 treatment)
#'   \item **Result size**: 39 genes × 6 samples heatmap (20 up + 19 down)
#'   \item **Result**: ggplot object with 1 layer (geom_tile)
#' }
#'
#' @details
#' **Interpreting Sample-Level Heatmaps**:
#' \itemize{
#'   \item **Sample clustering**: Ideally, control samples cluster together and treatment samples cluster together,
#'     indicating clear perturbation effects. Mixed clustering may indicate weak effects or technical issues.
#'   \item **Gene expression blocks**: Distinct color blocks indicate genes coordinately regulated across samples.
#'     Upregulated genes should show high expression (red) in treatment, low (blue) in control.
#'   \item **Replicate consistency**: Samples within same group should have similar color patterns.
#'     Large variation suggests technical noise or biological heterogeneity.
#'   \item **Row annotations**: Left sidebar shows whether gene is Up or Down regulated based on DEG analysis.
#'   \item **Column annotations**: Top bar shows sample groups (control vs treatment) if metadata available.
#' }
#'
#' **When to Use This Heatmap**:
#' \itemize{
#'   \item Quality control: Check if samples cluster by treatment group
#'   \item Validate DEG results: Visually confirm top genes show expected patterns
#'   \item Identify outlier samples: Detect technical failures or biological outliers
#'   \item Explore co-expression: Find genes with coordinated expression across samples
#'   \item Assess treatment effect: Visualize magnitude and consistency of perturbation
#'   \item Compare replicates: Evaluate within-group variation
#' }
#'
#' **Scaling Methods**:
#' \itemize{
#'   \item **"row" (default)**: Z-score per gene. Best for comparing relative expression across samples.
#'     Highlights which samples have highest/lowest expression for each gene.
#'   \item **"column"**: Z-score per sample. Best for comparing which genes are highest expressed within each sample.
#'   \item **"none"**: Raw log2(CPM+1) values. Best when absolute expression levels are important.
#' }
#'
#' **Difference from DEG Heatmap**:
#' \itemize{
#'   \item **gpdb_plot_heatmap_expr**: Shows RAW EXPRESSION across SAMPLES within ONE dataset.
#'     Visualizes sample-to-sample variation and replicate consistency.
#'   \item **gpdb_plot_heatmap_deg**: Shows LOG2 FOLD CHANGE across DATASETS.
#'     Compares perturbation effects across multiple experiments.
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Data Access** (get input data):
#' \itemize{
#'   \item \code{\link{gpdb_list_datasets}}: Find available datasets for a gene
#'   \item \code{\link{gpdb_get_info}}: Get detailed metadata for dataset
#'   \item \code{\link{gpdb_load_data}}: Load raw expression matrix and metadata
#'   \item \code{\link{gpdb_load_deg}}: Load DEG results used for gene selection
#' }
#'
#' **Complementary Analysis** (same dataset):
#' \itemize{
#'   \item \code{\link{gpdb_plot_volcano}}: Volcano plot for ALL genes in same dataset
#'   \item \code{\link{gpdb_enrich}}: Pathway enrichment on top DEGs
#'   \item \code{\link{gpdb_find_targets}}: Query aggregated targets across all datasets
#' }
#'
#' **Related Visualizations** (different perspectives):
#' \itemize{
#'   \item \code{\link{gpdb_plot_heatmap_deg}}: Cross-dataset DEG comparison heatmap
#'   \item \code{\link{gpdb_plot_network}}: Gene regulatory network
#'   \item \code{\link{gpdb_plot_enrichment}}: Enrichment results visualization
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I visualize gene expression across samples in my experiment?
#'   \item How do I check if control and treatment samples cluster separately?
#'   \item Can I see which genes show consistent expression patterns across replicates?
#'   \item How do I identify outlier samples in my perturbation experiment?
#'   \item What's the difference between expression heatmap and DEG heatmap?
#'   \item How do I visualize top upregulated and downregulated genes?
#'   \item Can I create a heatmap showing sample-to-sample variation?
#'   \item How do I assess replicate consistency in RNA-seq data?
#'   \item What do the colors mean in sample-level expression heatmaps?
#'   \item How do I interpret hierarchical clustering of samples?
#'   \item Should I use row or column scaling for expression heatmaps?
#'   \item How many genes should I include in the heatmap?
#'   \item Can I disable gene name labels for cleaner visualization?
#'   \item How do I check quality control of my perturbation experiment?
#'   \item What indicates a successful perturbation in expression heatmap?
#'   \item How do I identify co-expressed gene modules?
#'   \item Can I customize which genes to show in the heatmap?
#'   \item How do I export high-resolution heatmaps for publication?
#'   \item What's the best visualization for within-experiment variation?
#'   \item How do I validate DEG results visually?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Sample-Level Expression Heatmap (TESTED - runtime: 1.04 sec)
#' # ===========================================================================
#' # Scientific Question: Do control and treatment samples cluster separately?
#' # Dataset: TP53 siRNA in PANC-1 cells (3 control + 3 treatment samples)
#' # Expected output: Heatmap showing 39 genes (20 up + 19 down) × 6 samples
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' # Find available datasets
#' datasets <- gpdb_list_datasets(gene = "TP53")
#' cat("Found", nrow(datasets), "datasets for TP53\n")
#'
#' # Select first dataset
#' dataset_id <- datasets$dataset_id[1]
#'
#' # Get dataset info
#' info <- gpdb_get_info(dataset_id)
#' cat("Dataset:", info$gene, "in", info$cell_line, "\n")
#' cat(
#'   "Samples:", info$n_samples, "(", info$n_control, "control +",
#'   info$n_treat, "treatment)\n"
#' )
#'
#' # Create sample-level expression heatmap
#' p1 <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 20, # Top 20 upregulated genes
#'   top_down = 20, # Top 20 downregulated genes
#'   scale = "row", # Z-score per gene
#'   cluster_rows = TRUE, # Cluster genes
#'   cluster_cols = TRUE # Cluster samples
#' )
#'
#' # Plot is ready
#' print(p1)
#'
#' # Check: Do treatment samples cluster together?
#' cat("Check column dendrogram to verify sample grouping\n")
#'
#' # Save to file
#' # ggplot2::ggsave("expression_heatmap.pdf", p1, width = 10, height = 12)
#'
#' # ===========================================================================
#' # Example 2: Customize Gene Number for Focused View
#' # ===========================================================================
#' # Focus on fewer genes for cleaner visualization
#'
#' # Small focused heatmap (top 10 up + down)
#' p2_small <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 10, # Fewer genes
#'   top_down = 10,
#'   scale = "row"
#' )
#'
#' # Large comprehensive heatmap (top 30 up + down)
#' p2_large <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 30, # More genes
#'   top_down = 30,
#'   scale = "row",
#'   show_gene_names = FALSE # Hide labels for many genes
#' )
#'
#' # ===========================================================================
#' # Example 3: Disable Clustering to Preserve Ordering
#' # ===========================================================================
#' # Keep genes ordered by fold change (up at top, down at bottom)
#'
#' p3 <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 15,
#'   top_down = 15,
#'   cluster_rows = FALSE, # Preserve up/down ordering
#'   cluster_cols = FALSE, # Preserve sample order from metadata
#'   scale = "row"
#' )
#'
#' cat("Genes ordered by fold change: upregulated (top) → downregulated (bottom)\n")
#'
#' # ===========================================================================
#' # Example 4: Compare Different Scaling Methods
#' # ===========================================================================
#' # Understand impact of scaling on visualization
#'
#' # Row scaling (default) - compare across samples
#' p4_row <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 15,
#'   top_down = 15,
#'   scale = "row", # Z-score per gene
#'   title = "Row Scaling: Z-score per Gene"
#' )
#'
#' # No scaling - raw normalized values
#' p4_none <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 15,
#'   top_down = 15,
#'   scale = "none", # Raw log2(CPM+1)
#'   title = "No Scaling: Raw Expression"
#' )
#'
#' # Row scaling emphasizes relative changes; no scaling shows absolute levels
#'
#' # ===========================================================================
#' # Example 5: Quality Control Check
#' # ===========================================================================
#' # Identify potential issues with samples or perturbation
#'
#' # Create heatmap
#' p5 <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 20,
#'   top_down = 20,
#'   cluster_cols = TRUE # Enable sample clustering
#' )
#'
#' # Interpretation guide:
#' cat("\n=== Quality Control Checklist ===\n")
#' cat("✓ Control samples cluster together?\n")
#' cat("✓ Treatment samples cluster together?\n")
#' cat("✓ Clear separation between groups?\n")
#' cat("✓ Upregulated genes show high expression in treatment?\n")
#' cat("✓ Downregulated genes show low expression in treatment?\n")
#' cat("✓ Consistent patterns across replicates?\n")
#'
#' # ===========================================================================
#' # Example 6: Combine with Volcano Plot for Validation
#' # ===========================================================================
#' # Cross-validate DEG results between volcano and expression heatmap
#'
#' # Load DEG results to identify top genes
#' deg_data <- gpdb_load_deg(dataset_id)
#' top_genes <- head(deg_data$gene[order(deg_data$adj.P.Val)], 10)
#' cat("Top 10 significant genes:", paste(top_genes, collapse = ", "), "\n")
#'
#' # Create volcano plot
#' p6_volcano <- gpdb_plot_volcano(
#'   dataset_id = dataset_id,
#'   highlight_genes = top_genes
#' )
#'
#' # Create expression heatmap
#' p6_heatmap <- gpdb_plot_heatmap_expr(
#'   dataset_id = dataset_id,
#'   top_up = 15,
#'   top_down = 15
#' )
#'
#' # Compare: Genes significant in volcano should show clear patterns in heatmap
#'
#' # ===========================================================================
#' # Next Steps (Workflow Integration)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - Load expression data: ?gpdb_load_data
#' # - Volcano plot: ?gpdb_plot_volcano
#' # - Cross-dataset comparison: ?gpdb_plot_heatmap_deg
#' # - Enrichment analysis: ?gpdb_enrich
#' # - Network visualization: ?gpdb_plot_network
#' }
#'
gpdb_plot_heatmap_expr <- function(dataset_id,
                                   top_up = 20,
                                   top_down = 20,
                                   scale = "row",
                                   cluster_rows = TRUE,
                                   cluster_cols = TRUE,
                                   show_gene_names = TRUE,
                                   colors = NULL,
                                   show_values = FALSE,
                                   title = NULL,
                                   theme = NULL) {
  .gpdb_validate_dataset_id(dataset_id)

  # Load data (expression already has gene_name as rownames after preprocessing)
  message("Loading expression data and DEG results...")
  data_obj <- gpdb_load_data(dataset_id, normalize = TRUE)
  deg_data <- gpdb_load_deg(dataset_id, filter = FALSE)

  # Get dataset info for title
  info <- data_obj$info
  if (is.null(title)) {
    title <- paste0(info$gene, " in ", info$cell_line, " (", info$method, ")")
  }

  # Select top genes from DEG
  deg_data <- deg_data[!is.na(deg_data$logFC) & !is.na(deg_data$adj.P.Val), ]
  deg_data <- deg_data[order(deg_data$logFC, decreasing = TRUE), ]

  # Use gene column from DEG (fallback to gene_id if gene is empty)
  if ("gene" %in% names(deg_data)) {
    deg_genes <- ifelse(
      !is.na(deg_data$gene) & deg_data$gene != "",
      deg_data$gene,
      deg_data$gene_id
    )
  } else {
    deg_genes <- deg_data$gene_id
  }

  deg_data$gene_symbol <- deg_genes
  deg_data <- deg_data[order(deg_data$logFC, decreasing = TRUE), ]

  # Select top genes
  top_up_genes <- head(deg_data$gene_symbol[deg_data$logFC > 0], top_up)
  top_down_genes <- tail(deg_data$gene_symbol[deg_data$logFC < 0], top_down)

  selected_genes <- c(top_up_genes, top_down_genes)
  selected_genes <- selected_genes[!is.na(selected_genes) & selected_genes != ""]

  if (length(selected_genes) == 0) {
    stop("No valid genes found in DEG data", call. = FALSE)
  }

  # Extract expression matrix (rownames are already gene symbols)
  expr_matrix <- data_obj$expression
  expr_gene_names <- rownames(expr_matrix)

  # Match genes (case-insensitive)
  expr_genes_upper <- toupper(expr_gene_names)
  selected_genes_upper <- toupper(selected_genes)

  gene_idx <- which(expr_genes_upper %in% selected_genes_upper)

  if (length(gene_idx) == 0) {
    stop("No selected genes found in expression matrix", call. = FALSE)
  }

  # Filter expression matrix
  expr_matrix <- expr_matrix[gene_idx, , drop = FALSE]
  matched_genes <- rownames(expr_matrix)

  message("Selected ", length(gene_idx), " genes for heatmap")

  # Convert to numeric matrix
  expr_matrix <- as.matrix(expr_matrix)

  # Remove any non-numeric columns
  numeric_cols <- apply(expr_matrix, 2, function(x) is.numeric(x))
  expr_matrix <- expr_matrix[, numeric_cols, drop = FALSE]

  if (ncol(expr_matrix) == 0) {
    stop("No numeric expression data found", call. = FALSE)
  }

  # Apply scaling if requested
  if (scale == "row") {
    expr_matrix <- t(scale(t(expr_matrix)))
  } else if (scale == "column") {
    expr_matrix <- scale(expr_matrix)
  }

  # Determine which genes are up vs down
  matched_genes_upper <- toupper(matched_genes)
  top_up_upper <- toupper(top_up_genes)
  top_down_upper <- toupper(top_down_genes)

  gene_labels <- ifelse(matched_genes_upper %in% top_up_upper,
    "Upregulated",
    "Downregulated"
  )

  # Use default colors if not provided
  if (is.null(colors)) {
    colors <- .gpdb_get_palette("spectral", 11)
  }

  # Use default theme if not provided
  if (is.null(theme)) {
    theme <- .gpdb_theme_default()
  }

  # Prepare column annotation (sample groups)
  col_annotation <- NULL
  metadata <- data_obj$metadata
  if (!is.null(metadata) && "group" %in% names(metadata)) {
    col_annotation <- data.frame(Group = metadata$group)
    message("Added sample group annotations")
  }

  # Create beautiful ggplot2 heatmap
  legend_title <- if (scale == "row") "Z-score" else "log2(CPM+1)"

  ht <- .gpdb_ggplot_heatmap(
    expr_matrix,
    scale = scale,
    cluster_rows = cluster_rows,
    cluster_cols = cluster_cols,
    row_labels = if (show_gene_names) matched_genes else NULL,
    col_labels = colnames(expr_matrix),
    row_annotation = gene_labels,
    col_annotation = col_annotation,
    title = title,
    colors = colors,
    legend_title = legend_title,
    show_values = show_values,
    theme = theme
  )

  n_up_matched <- sum(gene_labels == "Upregulated")
  n_down_matched <- sum(gene_labels == "Downregulated")

  message(
    "Heatmap: ", n_up_matched, " up + ",
    n_down_matched, " down = ",
    nrow(expr_matrix), " genes × ", ncol(expr_matrix), " samples"
  )

  return(ht)
}


#' Compare Dataset Distribution Across Experimental Contexts
#'
#' Create publication-ready bar plot visualizing how datasets for a gene are distributed
#' across different experimental contexts (tissues, cell lines, or perturbation methods).
#' **Summarizes data availability** for exploring context-specific vs universal gene functions
#' and identifying well-studied vs understudied experimental systems.
#'
#' @param gene Character. Gene symbol to query.
#'   Default: None (required). Gene name will be formatted automatically (case-insensitive).
#'   Example: \code{"TP53"} or \code{"tp53"}. Use \code{gpdb_search()} to verify gene availability.
#' @param stratify_by Character. Stratification variable for grouping datasets.
#'   Default: \code{"tissue"} (first option). Options: "tissue", "cell_line", "method".
#'   Determines how datasets are categorized in bar plot. Each choice reveals different aspects:
#'   - "tissue": Shows tissue-type coverage (e.g., Lung, Liver, Brain)
#'   - "cell_line": Shows specific cell line usage (e.g., HeLa, A549, MCF7)
#'   - "method": Shows perturbation method distribution (e.g., CRISPR, siRNA, shRNA)
#' @param colors Character vector. Custom color palette for bars.
#'   Default: \code{NULL} (uses gradient blue-green palette, auto-scaled to number of categories).
#'   Provide vector of colors for custom palette. Colors interpolated if fewer than categories.
#' @param title Character. Custom plot title.
#'   Default: \code{NULL} (auto-generated: "Gene Datasets by Stratification").
#'   Provide custom title to override. Title is centered and bolded.
#' @param x.text.angle Numeric. Rotation angle for x-axis labels in degrees.
#'   Default: \code{45}. Range: 0-90. Use 45-90 for long category names (e.g., cell lines).
#'   Set 0 for horizontal labels if category names are short.
#' @param bar.alpha Numeric. Transparency of bars.
#'   Default: \code{0.9}. Range: 0-1. Lower values create semi-transparent bars.
#'   Use 0.7-0.9 for overlapping plots or layered visualizations.
#' @param bar.width Numeric. Width of bars relative to category spacing.
#'   Default: \code{0.7}. Range: 0.3-1.0. Higher values create wider bars with less spacing.
#'   Use 0.5-0.7 for cleaner appearance with many categories.
#' @param add.text Logical. Display count labels on top of bars.
#'   Default: \code{TRUE}. Shows exact dataset count above each bar.
#'   Set \code{FALSE} for cleaner appearance when counts are not critical.
#' @param theme.plot ggplot2 theme object. Custom theme for plot styling.
#'   Default: \code{NULL} (uses \code{.gpdb_theme_default()}, publication-ready academic theme).
#'   Can provide custom ggplot2 theme.
#'
#' @return ggplot object showing dataset distribution bar plot.
#'
#' **Plot Elements**:
#' \describe{
#'   \item{X-axis}{Categories (tissues, cell lines, or methods) sorted by dataset count (ascending)}
#'   \item{Y-axis}{Number of datasets in each category}
#'   \item{Bars}{Each bar represents one category, colored with gradient (left = fewer, right = more)}
#'   \item{Labels}{Dataset counts displayed on top of bars (if add.text = TRUE)}
#'   \item{Sorting}{Categories ordered from fewest to most datasets for easy comparison}
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Save plot: \code{ggplot2::ggsave("comparison.pdf", p, width = 10, height = 6)}
#'   \item Identify well-studied contexts: Focus on categories with most datasets (right side)
#'   \item Find understudied contexts: Target categories with few datasets (left side) for novelty
#'   \item Select datasets for analysis: Use top categories to choose representative experiments
#'   \item Compare across stratifications: Create all three plots (tissue/cell_line/method) for comprehensive view
#'   \item Filter datasets: Use \code{gpdb_list_datasets()} with context filters based on plot insights
#'   \item Plan experiments: Identify gaps in experimental coverage for new studies
#' }
#'
#' **Alternative Visualizations**:
#' \itemize{
#'   \item \code{\link{gpdb_plot_heatmap_deg}}: Heatmap comparing expression changes across selected datasets.
#'     Use when: Want to see actual perturbation effects after identifying datasets.
#'   \item \code{\link{gpdb_plot_network}}: Gene regulatory network for the gene.
#'     Use when: Interested in gene-gene relationships rather than dataset distribution.
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test gene**: TP53 (71 datasets across multiple contexts)
#'   \item **Test 1 - By tissue**: 0.44 sec, 23 tissue types
#'   \item **Test 2 - By cell_line**: 0.15 sec, 39 cell lines
#'   \item **Test 3 - By method**: 0.04 sec, 4 perturbation methods
#'   \item **Result**: ggplot object with bar plot, sorted categories, gradient coloring
#' }
#'
#' @details
#' **Interpreting Comparison Plots**:
#' \itemize{
#'   \item **Left to right ordering**: Bars sorted by count (ascending). Left = least studied, Right = most studied.
#'   \item **Color gradient**: Visual encoding of dataset density. Darker/brighter colors (right) = more data available.
#'   \item **Height = availability**: Taller bars indicate more experimental replicates and robustness.
#'   \item **Category gaps**: Missing categories indicate unexplored experimental contexts.
#'   \item **Count labels**: Exact numbers help prioritize dataset selection for analysis.
#' }
#'
#' **When to Use Each Stratification**:
#' \itemize{
#'   \item **tissue**: Answer "Which organs/tissues have been studied for this gene?"
#'     Use for: Understanding tissue-specific functions, selecting physiologically relevant contexts.
#'   \item **cell_line**: Answer "Which cell lines are commonly used for this gene?"
#'     Use for: Choosing appropriate cell models, comparing immortalized vs primary cells.
#'   \item **method**: Answer "Which perturbation methods are available for this gene?"
#'     Use for: Assessing data consistency across techniques (CRISPR vs siRNA), method selection.
#' }
#'
#' **Use Cases**:
#' \itemize{
#'   \item **Data availability assessment**: Quickly understand how much data exists before analysis
#'   \item **Context selection**: Choose well-studied contexts for robust findings or understudied for novelty
#'   \item **Gap identification**: Discover unexplored tissue types or cell lines for new experiments
#'   \item **Method comparison**: Check if gene has been studied with multiple perturbation techniques
#'   \item **Experimental design**: Identify which cell lines or tissues are commonly used
#'   \item **Publication planning**: Understand field coverage before starting new projects
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Dataset Discovery** (find datasets):
#' \itemize{
#'   \item \code{\link{gpdb_list_datasets}}: List all datasets for a gene with detailed metadata
#'   \item \code{\link{gpdb_get_info}}: Get detailed information for specific dataset
#'   \item \code{\link{gpdb_search}}: Search for genes in database
#' }
#'
#' **Downstream Analysis** (analyze selected datasets):
#' \itemize{
#'   \item \code{\link{gpdb_plot_heatmap_deg}}: Compare expression changes across datasets
#'   \item \code{\link{gpdb_plot_volcano}}: Volcano plot for individual dataset
#'   \item \code{\link{gpdb_find_targets}}: Query aggregated targets across all datasets
#' }
#'
#' **Related Visualizations** (other summaries):
#' \itemize{
#'   \item \code{\link{gpdb_plot_network}}: Gene regulatory network
#'   \item \code{\link{gpdb_plot_enrichment}}: Enrichment results visualization
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I see which tissues have been studied for my gene?
#'   \item Which cell lines are commonly used for TP53 experiments?
#'   \item How many datasets are available for each perturbation method?
#'   \item What's the distribution of datasets across experimental contexts?
#'   \item How do I identify well-studied vs understudied tissue types?
#'   \item Can I visualize data availability before starting analysis?
#'   \item Which experimental systems have the most datasets?
#'   \item How do I choose appropriate cell lines for my analysis?
#'   \item Are there enough CRISPR datasets for my gene of interest?
#'   \item What tissue types are missing from the database?
#'   \item How do I compare dataset coverage between different genes?
#'   \item Can I see which methods are most commonly used?
#'   \item How do I assess data availability for experimental design?
#'   \item What's the best way to visualize dataset distribution?
#'   \item Which contexts should I prioritize for robust findings?
#'   \item How do I identify gaps in experimental coverage?
#'   \item Can I customize the colors and appearance of the bar plot?
#'   \item How do I export high-resolution comparison plots?
#'   \item What does the ordering of bars tell me about data availability?
#'   \item How do I interpret dataset counts for analysis planning?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Tissue Distribution (TESTED - runtime: 0.44 sec)
#' # ===========================================================================
#' # Scientific Question: Which tissues have been studied for TP53?
#' # Use Case: Select tissue-relevant contexts for analysis
#' # Expected output: Bar plot showing 23 tissue types with TP53 data
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' # Visualize TP53 dataset distribution by tissue
#' p1_tissue <- gpdb_plot_comparison(
#'   gene = "TP53",
#'   stratify_by = "tissue"
#' )
#'
#' # Plot is ready
#' print(p1_tissue)
#'
#' # Interpretation: Taller bars (right) = more studied tissues
#' # Look for tissues relevant to your research question
#'
#' # Save to file
#' # ggplot2::ggsave("tp53_tissue_distribution.pdf", p1_tissue, width = 12, height = 6)
#'
#' # ===========================================================================
#' # Example 2: Cell Line Distribution
#' # ===========================================================================
#' # Question: Which cell lines are commonly used for TP53 studies?
#'
#' p2_cellline <- gpdb_plot_comparison(
#'   gene = "TP53",
#'   stratify_by = "cell_line"
#' )
#'
#' print(p2_cellline)
#'
#' # Interpretation: Identifies most frequently used cell models
#' # Right side = well-established cell line models
#' # Left side = less common but potentially novel models
#'
#' # ===========================================================================
#' # Example 3: Method Distribution
#' # ===========================================================================
#' # Question: What perturbation methods are available for TP53?
#'
#' p3_method <- gpdb_plot_comparison(
#'   gene = "TP53",
#'   stratify_by = "method"
#' )
#'
#' print(p3_method)
#'
#' # Interpretation: Shows method diversity and availability
#' # Compare CRISPR vs siRNA vs shRNA data availability
#' # Multiple methods = cross-validation opportunity
#'
#' # ===========================================================================
#' # Example 4: Customize Appearance
#' # ===========================================================================
#' # Publication-ready customization
#'
#' library(ggplot2)
#'
#' # Custom color palette (warm colors)
#' custom_colors <- c("#FEE5D9", "#FCAE91", "#FB6A4A", "#DE2D26", "#A50F15")
#'
#' p4_custom <- gpdb_plot_comparison(
#'   gene = "TP53",
#'   stratify_by = "tissue",
#'   colors = custom_colors, # Custom gradient
#'   x.text.angle = 60, # Steeper angle
#'   bar.alpha = 0.85, # Slightly transparent
#'   bar.width = 0.6, # Narrower bars
#'   title = "TP53 Perturbation Studies: Tissue Coverage"
#' )
#'
#' # Further customization with ggplot2
#' p4_final <- p4_custom +
#'   labs(
#'     subtitle = "Based on GPSAdb (71 datasets)",
#'     caption = "Data source: Guo et al. (2022) Nucleic Acids Res"
#'   ) +
#'   theme(
#'     plot.subtitle = element_text(hjust = 0.5, size = 11, color = "grey30")
#'   )
#'
#' # ===========================================================================
#' # Example 5: Compare Multiple Genes
#' # ===========================================================================
#' # Question: How does data availability compare between genes?
#'
#' # Compare TP53 vs MYC tissue coverage
#' p5_tp53 <- gpdb_plot_comparison("TP53", stratify_by = "tissue")
#' p5_myc <- gpdb_plot_comparison("MYC", stratify_by = "tissue")
#'
#' # Print side-by-side comparison
#' # library(patchwork)  # For combining plots
#' # p5_tp53 + p5_myc
#'
#' cat("Compare which gene has broader tissue coverage\n")
#'
#' # ===========================================================================
#' # Example 6: Workflow Integration - From Comparison to Analysis
#' # ===========================================================================
#' # Complete workflow: Explore → Select → Analyze
#'
#' # Step 1: Visualize tissue distribution
#' p6 <- gpdb_plot_comparison("TP53", stratify_by = "tissue")
#' print(p6)
#'
#' # Step 2: List datasets to see details
#' datasets <- gpdb_list_datasets(gene = "TP53")
#' cat("Total datasets:", nrow(datasets), "\n")
#'
#' # Step 3: Filter by high-coverage tissue (example: if "Lung" has most data)
#' # lung_datasets <- datasets[datasets$tissue == "Lung", ]
#' # cat("Lung datasets:", nrow(lung_datasets), "\n")
#'
#' # Step 4: Select representative datasets for analysis
#' # selected_ids <- head(lung_datasets$dataset_id, 5)
#'
#' # Step 5: Compare expression patterns
#' # heatmap <- gpdb_plot_heatmap_deg(selected_ids, top_n = 30)
#'
#' cat("Workflow: Comparison plot → Dataset selection → Expression analysis\n")
#'
#' # ===========================================================================
#' # Example 7: Identify Research Gaps
#' # ===========================================================================
#' # Find understudied contexts for novel research
#'
#' # Get raw data for analysis
#' datasets_all <- gpdb_list_datasets(gene = "TP53")
#'
#' # Summarize tissue coverage
#' tissue_summary <- table(datasets_all$tissue)
#' tissue_sorted <- sort(tissue_summary)
#'
#' cat("\n=== Research Gap Analysis ===\n")
#' cat("Understudied tissues (≤2 datasets):\n")
#' print(tissue_sorted[tissue_sorted <= 2])
#'
#' cat("\n Well-studied tissues (≥5 datasets):\n")
#' print(tissue_sorted[tissue_sorted >= 5])
#'
#' # Visualize
#' p7 <- gpdb_plot_comparison("TP53", stratify_by = "tissue")
#' cat("\nLeft bars = potential for novel discoveries\n")
#' cat("Right bars = robust data for validation\n")
#'
#' # ===========================================================================
#' # Next Steps (Workflow Integration)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - List datasets: ?gpdb_list_datasets
#' # - Get dataset info: ?gpdb_get_info
#' # - Cross-dataset heatmap: ?gpdb_plot_heatmap_deg
#' # - Find targets: ?gpdb_find_targets
#' # - Network analysis: ?gpdb_plot_network
#' }
#'
gpdb_plot_comparison <- function(gene,
                                 stratify_by = c("tissue", "cell_line", "method"),
                                 colors = NULL,
                                 title = NULL,
                                 x.text.angle = 45,
                                 bar.alpha = 0.9,
                                 bar.width = 0.7,
                                 add.text = TRUE,
                                 theme.plot = NULL) {
  gene <- .gpdb_format_genes(gene)[1]
  stratify_by <- match.arg(stratify_by)

  # Use default theme if not provided
  if (is.null(theme.plot)) {
    theme.plot <- .gpdb_theme_default()
  }

  # Get datasets
  datasets <- gpdb_list_datasets(gene = gene)

  if (nrow(datasets) == 0) {
    stop("No datasets found for gene: ", gene, call. = FALSE)
  }

  # Aggregate by stratification variable
  if (stratify_by == "tissue") {
    plot_data <- datasets |>
      dplyr::group_by(tissue) |>
      dplyr::summarise(n_datasets = dplyr::n(), .groups = "drop")
    x_var <- "tissue"
    x_label <- "Tissue"
  } else if (stratify_by == "cell_line") {
    plot_data <- datasets |>
      dplyr::group_by(cell_line) |>
      dplyr::summarise(n_datasets = dplyr::n(), .groups = "drop")
    x_var <- "cell_line"
    x_label <- "Cell Line"
  } else {
    plot_data <- datasets |>
      dplyr::group_by(method) |>
      dplyr::summarise(n_datasets = dplyr::n(), .groups = "drop")
    x_var <- "method"
    x_label <- "Method"
  }

  # Sort by n_datasets (ascending) - from left to right, fewer to more
  plot_data <- plot_data[order(plot_data$n_datasets), ]
  plot_data[[x_var]] <- factor(plot_data[[x_var]], levels = plot_data[[x_var]])

  if (is.null(title)) {
    title <- paste0(
      gene, " Datasets by ",
      tools::toTitleCase(gsub("_", " ", stratify_by))
    )
  }

  # Create color gradient if not provided
  if (is.null(colors)) {
    # Use default gradient palette
    colors <- .gpdb_get_palette("gradient_blue_green", nrow(plot_data))
  } else if (length(colors) < nrow(plot_data)) {
    colors <- grDevices::colorRampPalette(colors)(nrow(plot_data))
  }

  # Create plot
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = .data[[x_var]], y = .data[["n_datasets"]], fill = .data[[x_var]])) +
    ggplot2::geom_bar(
      stat = "identity",
      alpha = bar.alpha,
      width = bar.width,
      color = NA # No border
    ) +
    ggplot2::scale_fill_manual(values = colors) +
    ggplot2::labs(
      title = title,
      x = x_label,
      y = "Number of Datasets",
      fill = NULL
    ) +
    theme.plot +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5,
        size = 14,
        face = "bold",
        colour = "black"
      ),
      axis.text.y = ggplot2::element_text(size = 10, colour = "black"),
      axis.text.x = ggplot2::element_text(
        size = 11,
        colour = "black",
        angle = x.text.angle,
        hjust = ifelse(x.text.angle > 0, 1, 0.5),
        vjust = ifelse(x.text.angle > 0, 1, 0.5)
      ),
      axis.title.x = ggplot2::element_text(
        size = 12,
        colour = "black",
        face = "bold"
      ),
      axis.title.y = ggplot2::element_text(
        size = 12,
        colour = "black",
        face = "bold"
      ),
      legend.position = "none" # Remove legend since colors encode the same info as x-axis
    )

  # Add count labels on top of bars
  if (add.text) {
    p <- p + ggplot2::geom_text(
      ggplot2::aes(label = n_datasets),
      vjust = -0.5,
      size = 3.5,
      fontface = "bold",
      color = "black"
    )
  }

  message("Comparison plot: ", nrow(plot_data), " categories (sorted by count)")

  return(p)
}
