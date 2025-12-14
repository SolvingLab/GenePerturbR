#' Over-Representation Analysis for Gene Perturbation Results
#'
#' @description
#' Perform enrichment analysis on upregulated and downregulated genes separately.
#' Supports multiple enrichment databases and returns both significant and non-significant results.
#'
#' @param targets Either:
#'   \itemize{
#'     \item A gpdb_find_targets result (list with $upregulated and $downregulated)
#'     \item A named list with $up and $down gene vectors
#'     \item A character vector of gene symbols (will analyze as single list)
#'   }
#' @param enrich.type Character. Database to use:
#'   "GO", "KEGG", "Wiki", "Reactome", "MsigDB", "Mesh", "HgDisease", "Enrichrdb", "Self"
#' @param organism Character. "hs" (human), "mm" (mouse), etc. Default "hs"
#' @param filter.pcg Logical. Filter to protein-coding genes only (default TRUE).
#'   Recommended to enable for accurate enrichment as most databases (GO/KEGG) primarily annotate PCGs.
#' @param min.genes Integer. Minimum genes required after filtering (default 10, warns if below)
#' @param split.by Character. How to split genes:
#'   \itemize{
#'     \item "auto": Automatically detect structure (default)
#'     \item "direction": Separate up/down regulated genes
#'     \item "none": Analyze all genes together
#'   }
#' @param top_up Integer. Number of top upregulated genes to analyze (default: all)
#' @param top_down Integer. Number of top downregulated genes to analyze (default: all)
#' @param p.cutoff Numeric. P-value cutoff for significance (default 0.05)
#' @param q.cutoff Numeric. Q-value cutoff for significance (default 0.05)
#' @param background.genes Character vector. Background genes (default: NULL, uses database genes)
#' @param GO.ont Character. GO ontology: "bp", "cc", "mf", or "all" (default "bp")
#' @param KEGG.category Character. KEGG category (default "pathway")
#' @param Msigdb.category Character. MSigDB category (default "H" for Hallmark)
#' @param return.all Logical. Return all pathways or only significant ones (default TRUE)
#' @param ... Additional arguments passed to underlying enrichment functions
#'
#' @return A list containing:
#'   \itemize{
#'     \item upregulated: ORA results for upregulated genes (list with $All and $Sig)
#'     \item downregulated: ORA results for downregulated genes (list with $All and $Sig)
#'     \item params: Parameters used for the analysis
#'   }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # From gpdb_find_targets result
#' targets <- gpdb_find_targets("TP53", top_n = 100)
#' enrich_res <- gpdb_enrich(targets, enrich.type = "GO")
#'
#' # From gene list
#' genes <- c("MYC", "TP53", "KRAS", "EGFR")
#' enrich_res <- gpdb_enrich(genes, enrich.type = "KEGG", split.by = "none")
#'
#' # Custom gene lists
#' gene_lists <- list(
#'   up = c("MYC", "CCND1", "CDK4"),
#'   down = c("TP53", "RB1", "CDKN1A")
#' )
#' enrich_res <- gpdb_enrich(gene_lists, enrich.type = "GO")
#' }
gpdb_enrich <- function(targets,
                        enrich.type = c("GO", "KEGG", "Wiki", "Reactome", "MsigDB", "Mesh", "HgDisease", "Enrichrdb", "Self"),
                        organism = "hs",
                        filter.pcg = TRUE,
                        min.genes = 10,
                        split.by = c("auto", "direction", "none"),
                        top_up = NULL,
                        top_down = NULL,
                        p.cutoff = 0.05,
                        q.cutoff = 0.05,
                        background.genes = NULL,
                        GO.ont = "bp",
                        KEGG.category = "pathway",
                        Msigdb.category = "H",
                        return.all = TRUE,
                        ...) {
  # Match arguments
  enrich.type <- match.arg(enrich.type)
  split.by <- match.arg(split.by)

  # Check BioEnricher availability
  if (!requireNamespace("BioEnricher", quietly = TRUE)) {
    stop("BioEnricher package required. Install from: https://github.com/ZaoQuliu/BioEnricher",
      call. = FALSE
    )
  }

  # Parse input and extract gene lists
  gene_lists <- .parse_enrich_input(targets, split.by, top_up, top_down)

  # Filter to protein-coding genes (recommended for accurate enrichment)
  if (filter.pcg) {
    gene_lists <- .filter_protein_coding(gene_lists, organism, min.genes)
  }

  # Store parameters
  params <- list(
    enrich.type = enrich.type,
    organism = organism,
    filter.pcg = filter.pcg,
    split.by = split.by,
    p.cutoff = p.cutoff,
    q.cutoff = q.cutoff,
    GO.ont = GO.ont,
    KEGG.category = KEGG.category,
    Msigdb.category = Msigdb.category,
    n_genes_up = if (!is.null(gene_lists$up)) length(gene_lists$up) else 0,
    n_genes_down = if (!is.null(gene_lists$down)) length(gene_lists$down) else 0
  )

  # Perform enrichment for each direction
  result <- list(params = params)

  if (!is.null(gene_lists$up) && length(gene_lists$up) > 0) {
    message("\n=== Enriching upregulated genes (n=", length(gene_lists$up), ") ===")
    result$upregulated <- BioEnricher::ORA(
      genes = gene_lists$up,
      enrich.type = enrich.type,
      organism = organism,
      p.cutoff = p.cutoff,
      q.cutoff = q.cutoff,
      background.genes = background.genes,
      GO.ont = GO.ont,
      KEGG.category = KEGG.category,
      Msigdb.category = Msigdb.category,
      ...
    )
  } else {
    result$upregulated <- list(All = data.frame(), Sig = data.frame())
  }

  if (!is.null(gene_lists$down) && length(gene_lists$down) > 0) {
    message("\n=== Enriching downregulated genes (n=", length(gene_lists$down), ") ===")
    result$downregulated <- BioEnricher::ORA(
      genes = gene_lists$down,
      enrich.type = enrich.type,
      organism = organism,
      p.cutoff = p.cutoff,
      q.cutoff = q.cutoff,
      background.genes = background.genes,
      GO.ont = GO.ont,
      KEGG.category = KEGG.category,
      Msigdb.category = Msigdb.category,
      ...
    )
  } else {
    result$downregulated <- list(All = data.frame(), Sig = data.frame())
  }

  # Summary
  n_up_sig <- nrow(result$upregulated$Sig)
  n_down_sig <- nrow(result$downregulated$Sig)
  n_up_all <- nrow(result$upregulated$All)
  n_down_all <- nrow(result$downregulated$All)

  message("\n=== Enrichment Summary ===")
  message("Upregulated:   ", n_up_sig, " significant / ", n_up_all, " total pathways")
  message("Downregulated: ", n_down_sig, " significant / ", n_down_all, " total pathways")

  class(result) <- c("gpdb_enrichment", "list")
  return(result)
}


#' Plot Paired Enrichment Dotplot
#'
#' @description
#' Create side-by-side dotplots showing enrichment results for upregulated
#' and downregulated genes. Inspired by GSEA paired plots.
#'
#' @param enrich_result A gpdb_enrichment object from gpdb_enrich()
#' @param show.term.num Integer. Number of top pathways to show per side (default 15)
#' @param use.all Logical. Use all results or only significant ones (default FALSE, use significant)
#' @param x Character. X-axis variable: "FoldEnrich", "pvalue", "p.adjust", "Count" (default "FoldEnrich")
#' @param color.by Character. Color variable: "p.adjust", "pvalue", "Count", "FoldEnrich" (default "p.adjust")
#' @param size.by Character. Size variable: "Count", "pvalue", "p.adjust", "FoldEnrich" (default "Count")
#' @param colors Character vector. Color palette (default: BrBG diverging palette)
#' @param size.range Numeric vector of length 2. Size range for points (default c(3, 8))
#' @param title Character. Main title (default: auto-generated)
#' @param up.title Character. Title for upregulated panel (default "Upregulated Genes")
#' @param down.title Character. Title for downregulated panel (default "Downregulated Genes")
#' @param legend.position Character. Legend position: "right", "bottom", "none" (default "right")
#' @param theme ggplot2 theme (default theme_bw)
#'
#' @return A patchwork object combining two dotplots
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' targets <- gpdb_find_targets("TP53", top_n = 100)
#' enrich_res <- gpdb_enrich(targets, enrich.type = "GO")
#' gpdb_plot_enrichment(enrich_res)
#'
#' # Customize appearance
#' gpdb_plot_enrichment(
#'   enrich_res,
#'   show.term.num = 20,
#'   x = "pvalue",
#'   color.by = "FoldEnrich"
#' )
#' }
gpdb_plot_enrichment <- function(enrich_result,
                                 show.term.num = 15,
                                 use.all = FALSE,
                                 x = "FoldEnrich",
                                 color.by = "p.adjust",
                                 size.by = "Count",
                                 colors = c(
                                   "#003c30", "#01665e", "#35978f", "#80cdc1", "#c7eae5",
                                   "#f6e8c3", "#dfc27d", "#bf812d", "#8c510a", "#543005"
                                 ),
                                 size.range = c(3, 8),
                                 title = NULL,
                                 up.title = "Upregulated Genes",
                                 down.title = "Downregulated Genes",
                                 legend.position = "right",
                                 theme = ggplot2::theme_bw(base_rect_size = 1.5)) {
  # Validate input
  if (!inherits(enrich_result, "gpdb_enrichment")) {
    stop("enrich_result must be a gpdb_enrichment object from gpdb_enrich()", call. = FALSE)
  }

  # Check required packages
  if (!requireNamespace("patchwork", quietly = TRUE)) {
    stop("patchwork package required. Install with: install.packages('patchwork')", call. = FALSE)
  }

  # Extract data
  up_data <- if (use.all) enrich_result$upregulated$All else enrich_result$upregulated$Sig
  down_data <- if (use.all) enrich_result$downregulated$All else enrich_result$downregulated$Sig

  # Check if we have data
  has_up <- !is.null(up_data) && nrow(up_data) > 0
  has_down <- !is.null(down_data) && nrow(down_data) > 0

  if (!has_up && !has_down) {
    stop("No enrichment results to plot. Try use.all = TRUE to show non-significant results.",
      call. = FALSE
    )
  }

  # Create plots
  plot_up <- if (has_up) {
    .create_enrich_dotplot(
      up_data,
      show.term.num = show.term.num,
      x = x,
      color.by = color.by,
      size.by = size.by,
      colors = colors,
      size.range = size.range,
      title = up.title,
      legend.position = legend.position,
      theme = theme
    )
  } else {
    .create_empty_plot(up.title)
  }

  plot_down <- if (has_down) {
    .create_enrich_dotplot(
      down_data,
      show.term.num = show.term.num,
      x = x,
      color.by = color.by,
      size.by = size.by,
      colors = rev(colors), # Reverse colors for down-regulated
      size.range = size.range,
      title = down.title,
      legend.position = legend.position,
      theme = theme
    )
  } else {
    .create_empty_plot(down.title)
  }

  # Create main title
  if (is.null(title)) {
    title <- paste0(enrich_result$params$enrich.type, " Enrichment Analysis")
    if (enrich_result$params$enrich.type == "GO") {
      title <- paste0(title, " (", toupper(enrich_result$params$GO.ont), ")")
    }
  }

  # Build caption
  caption_parts <- c(
    paste0("Database: ", enrich_result$params$enrich.type),
    paste0("Cutoff: p<", enrich_result$params$p.cutoff, ", q<", enrich_result$params$q.cutoff)
  )
  caption_text <- paste(caption_parts, collapse = " | ")

  # Combine plots
  combined <- patchwork::wrap_plots(plot_up, plot_down, ncol = 2) +
    patchwork::plot_annotation(
      title = title,
      caption = caption_text,
      theme = ggplot2::theme(
        plot.title = ggplot2::element_text(hjust = 0.5, size = 14, face = "bold", colour = "black"),
        plot.caption = ggplot2::element_text(hjust = 0, size = 9, colour = "gray40")
      )
    )

  # Calculate dimensions
  n_pathways <- max(
    min(nrow(up_data), show.term.num),
    min(nrow(down_data), show.term.num),
    5,
    na.rm = TRUE
  )

  # Get max pathway label length for width calculation
  all_labels <- c(
    if (has_up) as.character(head(up_data$Description, show.term.num)) else character(0),
    if (has_down) as.character(head(down_data$Description, show.term.num)) else character(0)
  )
  max_label_len <- if (length(all_labels) > 0) max(nchar(all_labels), na.rm = TRUE) else 50

  width <- max(14, 8 + max_label_len * 0.08)
  height <- max(6, 3 + n_pathways * 0.25)

  attr(combined, "width") <- width
  attr(combined, "height") <- height

  return(combined)
}


# ==============================================================================
# Internal Helper Functions
# ==============================================================================

#' Parse input for enrichment analysis
#' @keywords internal
.parse_enrich_input <- function(targets, split.by, top_up, top_down) {
  # Case 1: gpdb_find_targets result
  if (is.list(targets) && all(c("upregulated", "downregulated") %in% names(targets))) {
    up_genes <- if (!is.null(targets$upregulated)) {
      genes <- targets$upregulated$target_gene
      if (!is.null(top_up)) genes <- head(genes, top_up)
      genes
    } else {
      NULL
    }

    down_genes <- if (!is.null(targets$downregulated)) {
      genes <- targets$downregulated$target_gene
      if (!is.null(top_down)) genes <- head(genes, top_down)
      genes
    } else {
      NULL
    }

    return(list(up = up_genes, down = down_genes))
  }

  # Case 2: Named list with $up and $down
  if (is.list(targets) && all(c("up", "down") %in% names(targets))) {
    up_genes <- targets$up
    down_genes <- targets$down
    if (!is.null(top_up)) up_genes <- head(up_genes, top_up)
    if (!is.null(top_down)) down_genes <- head(down_genes, top_down)
    return(list(up = up_genes, down = down_genes))
  }

  # Case 3: Character vector
  if (is.character(targets)) {
    if (split.by == "none") {
      # Analyze all together as "up"
      return(list(up = targets, down = NULL))
    } else {
      warning("Cannot split character vector by direction. Use split.by='none' or provide structured input.",
        call. = FALSE
      )
      return(list(up = targets, down = NULL))
    }
  }

  # Case 4: Data frame (assume has logFC column)
  if (is.data.frame(targets)) {
    if ("logFC" %in% names(targets) || "logfc" %in% names(targets)) {
      logfc_col <- if ("logFC" %in% names(targets)) "logFC" else "logfc"
      gene_col <- if ("target_gene" %in% names(targets)) {
        "target_gene"
      } else if ("gene" %in% names(targets)) {
        "gene"
      } else {
        stop("Cannot find gene column in data frame", call. = FALSE)
      }

      up_genes <- targets[targets[[logfc_col]] > 0, gene_col]
      down_genes <- targets[targets[[logfc_col]] < 0, gene_col]

      if (!is.null(top_up)) up_genes <- head(up_genes, top_up)
      if (!is.null(top_down)) down_genes <- head(down_genes, top_down)

      return(list(up = up_genes, down = down_genes))
    }
  }

  stop("Cannot parse input. Provide gpdb_find_targets result, named list, or character vector.",
    call. = FALSE
  )
}


#' Create enrichment dotplot (internal)
#' @keywords internal
.create_enrich_dotplot <- function(res,
                                   show.term.num,
                                   x,
                                   color.by,
                                   size.by,
                                   colors,
                                   size.range,
                                   title,
                                   legend.position,
                                   theme) {
  # Prepare x-axis
  x_col <- x
  x_lab <- x

  if (x == "pvalue") {
    res$plot_x <- -log10(res$pvalue)
    x_col <- "plot_x"
    x_lab <- bquote(~ -Log[10] ~ italic("P-value"))
  } else if (x == "p.adjust") {
    res$plot_x <- -log10(res$p.adjust)
    x_col <- "plot_x"
    x_lab <- bquote(~ -Log[10] ~ "FDR")
  } else if (x == "FoldEnrich") {
    x_lab <- "Fold Enrichment"
  }

  # Prepare color
  color_col <- color.by
  color_title <- color.by

  if (color.by == "pvalue") {
    res$plot_color <- -log10(res$pvalue)
    color_col <- "plot_color"
    color_title <- bquote(~ -Log[10] ~ italic("P-value"))
  } else if (color.by == "p.adjust") {
    res$plot_color <- -log10(res$p.adjust)
    color_col <- "plot_color"
    color_title <- bquote(~ -Log[10] ~ "FDR")
  }

  # Prepare size
  size_col <- size.by
  size_title <- size.by

  if (size.by == "pvalue") {
    res$plot_size <- -log10(res$pvalue)
    size_col <- "plot_size"
    size_title <- bquote(~ -Log[10] ~ italic("P-value"))
  } else if (size.by == "p.adjust") {
    res$plot_size <- -log10(res$p.adjust)
    size_col <- "plot_size"
    size_title <- bquote(~ -Log[10] ~ "FDR")
  }

  # Select top pathways
  show.term.num <- min(show.term.num, nrow(res))
  if (show.term.num == 0) {
    return(.create_empty_plot(title))
  }

  plot_data <- res %>%
    dplyr::arrange(pvalue) %>%
    dplyr::slice(1:show.term.num) %>%
    dplyr::arrange(dplyr::desc(.data[[x_col]])) %>%
    dplyr::mutate(Description = factor(Description, levels = rev(Description)))

  # Create plot
  p <- ggplot2::ggplot(
    plot_data,
    ggplot2::aes(
      x = .data[[x_col]],
      y = Description,
      fill = .data[[color_col]],
      size = .data[[size_col]]
    )
  ) +
    ggplot2::geom_point(shape = 21, color = "black") +
    theme +
    ggplot2::labs(
      fill = color_title,
      size = size_title,
      x = x_lab,
      y = NULL,
      title = title
    ) +
    ggplot2::scale_fill_gradientn(colours = colors) +
    ggplot2::scale_size(range = size.range) +
    ggplot2::scale_y_discrete(
      labels = function(x) Hmisc::capitalize(tolower(x)),
      position = "right"
    ) +
    ggplot2::theme(
      axis.text.x = ggplot2::element_text(size = 10, colour = "black"),
      axis.title.x = ggplot2::element_text(size = 13, colour = "black", face = "bold"),
      axis.text.y = ggplot2::element_text(size = 13, colour = "black"),
      axis.title.y = ggplot2::element_blank(),
      panel.background = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 11, colour = "black"),
      legend.title = ggplot2::element_text(size = 13, colour = "black", face = "bold"),
      legend.background = ggplot2::element_blank(),
      legend.key = ggplot2::element_blank(),
      legend.position = legend.position,
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, colour = "black", face = "bold")
    )

  return(p)
}


#' Create empty plot with message
#' @keywords internal
.create_empty_plot <- function(title) {
  ggplot2::ggplot() +
    ggplot2::annotate(
      "text",
      x = 0.5, y = 0.5,
      label = "No significant pathways",
      size = 5,
      color = "gray50"
    ) +
    ggplot2::labs(title = title) +
    ggplot2::theme_void() +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 12, face = "bold")
    )
}


#' Filter to protein-coding genes
#' @keywords internal
.filter_protein_coding <- function(gene_lists, organism, min.genes = 10) {
  if (!requireNamespace("BioEnricher", quietly = TRUE)) {
    message("BioEnricher not available, skipping PCG filtering")
    return(gene_lists)
  }

  message("Filtering protein-coding genes (set filter.pcg=FALSE to disable)...")

  # Filter upregulated genes
  if (!is.null(gene_lists$up) && length(gene_lists$up) > 0) {
    n_before <- length(gene_lists$up)
    gene_lists$up <- BioEnricher::pickPCG(gene_lists$up, org = organism)
    n_after <- length(gene_lists$up)
    pct <- round(n_after / n_before * 100, 1)

    message(
      "  Upregulated: ", n_before, " → ", n_after,
      " genes (", pct, "% PCG)"
    )

    if (n_after < min.genes) {
      warning(
        "Only ", n_after, " protein-coding genes in upregulated list. ",
        "Consider increasing top_up or set filter.pcg = FALSE",
        call. = FALSE
      )
    }
  }

  # Filter downregulated genes
  if (!is.null(gene_lists$down) && length(gene_lists$down) > 0) {
    n_before <- length(gene_lists$down)
    gene_lists$down <- BioEnricher::pickPCG(gene_lists$down, org = organism)
    n_after <- length(gene_lists$down)
    pct <- round(n_after / n_before * 100, 1)

    message(
      "  Downregulated: ", n_before, " → ", n_after,
      " genes (", pct, "% PCG)"
    )

    if (n_after < min.genes) {
      warning(
        "Only ", n_after, " protein-coding genes in downregulated list. ",
        "Consider increasing top_down or set filter.pcg = FALSE",
        call. = FALSE
      )
    }
  }

  return(gene_lists)
}
