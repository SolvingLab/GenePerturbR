#' Plot Volcano Plot for Single Dataset
#'
#' Create volcano plot showing differential expression from ONE dataset.
#' Displays log2 fold change vs -log10(adjusted p-value).
#'
#' @param dataset_id Character. Dataset ID (e.g., "D10001"). Required.
#' @param padj_cutoff Numeric. Significance threshold (default 0.05)
#' @param logfc_cutoff Numeric. Fold change threshold in log2 scale (default 1 = 2-fold)
#' @param nlabel Integer. Number of top genes to label (default 10)
#' @param highlight_genes Character vector. Specific genes to highlight (optional)
#' @param colors Three colors for Up/NoSig/Down genes (default: blue/grey/red)
#' @param point.maxsize Maximum size of points (default 4)
#' @param point.alpha Alpha transparency of points (default 0.8)
#' @param intercept.width Width of threshold lines (default 0.65)
#' @param label.size Size of gene labels (default 3.5)
#' @param label.bg Background color for labels (default "white")
#' @param label.bg.r Background corner radius for labels (default 0.1)
#' @param legend.position Position of legend (default "bottom")
#' @param title Character. Plot title (auto-generated if NULL)
#' @param theme.plot ggplot2 theme (default theme_bw)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic volcano plot
#' p1 <- gpdb_plot_volcano("D10001")
#'
#' # Highlight specific genes
#' p2 <- gpdb_plot_volcano("D10001", highlight_genes = c("MYC", "TP53", "KRAS"))
#'
#' # Customize appearance
#' p3 <- gpdb_plot_volcano("D10001",
#'   padj_cutoff = 0.01,
#'   logfc_cutoff = 2,
#'   nlabel = 20
#' )
#' }
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


#' Plot DEG Heatmap Across Multiple Datasets
#'
#' Compare differential expression (logFC) across multiple datasets.
#' Each column represents one dataset, each row represents one gene.
#' This function is designed for cross-dataset comparison of perturbation effects.
#'
#' @param dataset_ids Character vector. Dataset IDs to compare (2+ required)
#' @param genes Character vector. Specific genes to show (optional)
#' @param top_n Integer. If genes not specified, show top N DEGs (default 50)
#' @param cluster_rows Logical. Cluster genes (default TRUE)
#' @param cluster_cols Logical. Cluster datasets (default TRUE)
#' @param scale Character. "row" (z-score), "column", or "none" (default "row")
#' @param colors Character vector. Color palette (NULL uses default spectral)
#' @param show_values Logical. Show values in cells (default FALSE)
#' @param title Character. Plot title
#' @param theme ggplot2 theme (NULL uses default)
#'
#' @return ggplot object
#' @export
#'
#' @seealso \code{\link{gpdb_plot_heatmap_expr}} for single-dataset sample-level heatmaps
#'
#' @examples
#' \dontrun{
#' # Compare TP53 knockdown effects across datasets
#' datasets <- gpdb_list_datasets(gene = "TP53")
#' gpdb_plot_heatmap_deg(head(datasets$dataset_id, 5), top_n = 30)
#' }
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


#' Plot Expression Heatmap for Single Dataset
#'
#' Visualize sample-level gene expression for one dataset.
#' Each column represents one sample, each row represents one gene.
#' Shows top upregulated and downregulated genes with sample group annotations.
#'
#' @param dataset_id Character. Dataset ID to visualize
#' @param top_up Integer. Number of top upregulated genes (default 20)
#' @param top_down Integer. Number of top downregulated genes (default 20)
#' @param scale Character. "row" (z-score), "column", or "none" (default "row")
#' @param cluster_rows Logical. Cluster genes (default TRUE)
#' @param cluster_cols Logical. Cluster samples (default TRUE)
#' @param show_gene_names Logical. Show gene names (default TRUE)
#' @param colors Character vector. Color palette (NULL uses default)
#' @param show_values Logical. Show values in cells (default FALSE)
#' @param title Character. Plot title (auto-generated if NULL)
#' @param theme ggplot2 theme (NULL uses default)
#'
#' @return ggplot object
#' @export
#'
#' @seealso \code{\link{gpdb_plot_heatmap_deg}} for multi-dataset DEG comparison
#'
#' @examples
#' \dontrun{
#' # Plot sample-level expression heatmap
#' gpdb_plot_heatmap_expr("D10001", top_up = 20, top_down = 20)
#'
#' # Customize number of genes
#' gpdb_plot_heatmap_expr("D10001", top_up = 15, top_down = 15)
#' }
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


#' Plot Comparison Across Datasets
#'
#' Create bar plot showing dataset distribution by tissue, cell line, or method.
#' Bars are sorted by count (ascending) and colored in gradient.
#'
#' @param gene Character. Gene symbol
#' @param stratify_by Character. "tissue", "cell_line", or "method"
#' @param colors Color palette. Default uses a beautiful gradient
#' @param title Character. Plot title (centered)
#' @param x.text.angle Angle of x-axis text (default 45)
#' @param bar.alpha Alpha transparency of bars (default 0.9)
#' @param bar.width Width of bars (default 0.7)
#' @param add.text Logical. Add count labels on bars (default TRUE)
#' @param theme.plot ggplot2 theme (default theme_classic)
#'
#' @return ggplot object
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare TP53 datasets by tissue
#' p1 <- gpdb_plot_comparison("TP53", stratify_by = "tissue")
#'
#' # Compare by cell line
#' p2 <- gpdb_plot_comparison("TP53", stratify_by = "cell_line")
#'
#' # Compare by method
#' p3 <- gpdb_plot_comparison("TP53", stratify_by = "method")
#' }
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


# ================================================================================
# Backward Compatibility Aliases
# ================================================================================

#' @rdname gpdb_plot_heatmap_deg
#' @export
gpdb_plot_heatmap <- gpdb_plot_heatmap_deg

#' @rdname gpdb_plot_heatmap_expr
#' @export
gpdb_plot_heatmap_single <- gpdb_plot_heatmap_expr
