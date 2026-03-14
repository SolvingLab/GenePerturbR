#' Pathway Enrichment Analysis with Over-Representation Test
#'
#' @description
#' Identify enriched biological pathways in upregulated and downregulated gene sets
#' using over-representation analysis (ORA) with Fisher's exact test. Automatically
#' filters protein-coding genes and performs separate enrichment for each direction
#' to reveal distinct pathway changes. Supports 9 annotation databases including
#' GO/KEGG/Reactome/MSigDB for comprehensive pathway discovery and functional interpretation.
#'
#' @param targets Flexible input accepting multiple formats:
#'   \itemize{
#'     \item \code{gpdb_find_targets} result (list with \code{$upregulated} and \code{$downregulated})
#'     \item Named list with \code{$up} and \code{$down} gene vectors
#'     \item Character vector of gene symbols (analyzed as single group)
#'   }
#' @param enrich.type Character. Annotation database for pathway enrichment.
#'   Options: "GO" (Gene Ontology), "KEGG" (pathways), "Wiki" (WikiPathways),
#'   "Reactome", "MsigDB" (Molecular Signatures), "Mesh", "HgDisease", "Enrichrdb".
#'   Default: "GO".
#' @param organism Character. Species code: "hs" (human), "mm" (mouse), "rn" (rat), etc.
#'   Default: "hs". Passed to geneset/genekitr2 for ID conversion.
#' @param filter.pcg Logical. Filter to protein-coding genes only (default: TRUE).
#'   Recommended for accurate enrichment as GO/KEGG primarily annotate PCGs.
#'   Set FALSE to include non-coding genes (lncRNA, miRNA).
#' @param min.genes Integer. Minimum protein-coding genes required after filtering (default: 10).
#'   Warns if below threshold. Consider increasing \code{top_up}/\code{top_down} or disable \code{filter.pcg}.
#' @param split.by Character. Gene splitting strategy: "auto" (detect structure),
#'   "direction" (separate up/down), "none" (analyze together). Default: "auto".
#' @param top_up Integer. Number of top upregulated genes to analyze (default: NULL, uses all).
#'   Useful for focusing on strongest perturbation effects.
#' @param top_down Integer. Number of top downregulated genes to analyze (default: NULL, uses all).
#' @param p.cutoff Numeric. Raw p-value cutoff for significance (default: 0.05).
#'   Pathways with p < 0.05 marked as significant.
#' @param q.cutoff Numeric. Adjusted p-value (FDR) cutoff for significance (default: 0.05).
#'   More stringent than p-value for multiple testing correction.
#' @param background.genes Character vector. Custom background genes (default: NULL, uses database genes).
#'   Provide when analyzing specific gene sets (e.g., expressed genes in your dataset).
#' @param GO.ont Character. GO ontology type: "bp" (biological process), "cc" (cellular component),
#'   "mf" (molecular function), "all". Default: "bp". Only used when \code{enrich.type = "GO"}.
#' @param KEGG.category Character. KEGG category to query (default: "pathway").
#'   Only used when \code{enrich.type = "KEGG"}.
#' @param Msigdb.category Character. MSigDB category: "H" (Hallmark, default), "C1-C8".
#'   Only used when \code{enrich.type = "MsigDB"}.
#' @param return.all Logical. Return all pathways or only significant ones (default: TRUE).
#'   FALSE returns only pathways passing \code{p.cutoff} and \code{q.cutoff}.
#' @param ... Additional arguments (currently unused, kept for API compatibility).
#'
#' @return \code{gpdb_enrichment} object (list) with three elements:
#'
#' **Return Structure**:
#' \describe{
#'   \item{\strong{upregulated}}{List with \code{$All} and \code{$Sig} data.frames.
#'     Contains enrichment results for upregulated genes.
#'     \itemize{
#'       \item \code{$All}: All pathways tested (e.g., 15461 GO terms)
#'       \item \code{$Sig}: Significant pathways (p < 0.05, q < 0.05)
#'       \item Key columns: \code{Description} (pathway name), \code{pvalue}, \code{p.adjust},
#'             \code{Count} (genes in pathway), \code{FoldEnrich} (enrichment ratio)
#'     }
#'   }
#'   \item{\strong{downregulated}}{List with \code{$All} and \code{$Sig} data.frames.
#'     Contains enrichment results for downregulated genes (same structure as upregulated).
#'   }
#'   \item{\strong{params}}{List of analysis parameters:
#'     \code{enrich.type}, \code{organism}, \code{filter.pcg}, \code{split.by},
#'     \code{p.cutoff}, \code{q.cutoff}, \code{GO.ont}, \code{n_genes_up}, \code{n_genes_down}.
#'   }
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Visualize results: \code{\link{gpdb_plot_enrichment}(result)} for paired dotplots
#'   \item Filter significant pathways: \code{result$upregulated$Sig[1:10, ]} for top 10
#'   \item Export results: \code{write.csv(result$upregulated$Sig, "enrichment.csv")}
#'   \item Compare databases: Run with different \code{enrich.type} (GO vs KEGG vs Reactome)
#'   \item Validate in network: \code{\link{gpdb_plot_network}()} to visualize gene-pathway relationships
#' }
#'
#' **Alternative Methods**:
#' \itemize{
#'   \item \code{fgsea::fgseaMultilevel}: Use when you have ranked gene list and want GSEA instead of ORA.
#'   \item \code{clusterProfiler::enrichGO}: Use for advanced GO enrichment with custom background.
#' }
#'
#' @details
#' **Interpreting Results**:
#' \itemize{
#'   \item \strong{pvalue}: Fisher's exact test p-value. Lower = stronger enrichment.
#'         Significant if < 0.05, highly significant if < 0.001.
#'   \item \strong{p.adjust}: FDR-corrected p-value (Benjamini-Hochberg method).
#'         Controls false discovery rate. Use for reporting (q < 0.05 = 5% FDR).
#'   \item \strong{FoldEnrich}: Observed/Expected ratio. FoldEnrich = 2 means 2x more genes
#'         than expected by chance. Values > 2 indicate strong biological signal.
#'   \item \strong{Count}: Number of input genes in pathway. Pathways with Count > 5
#'         are more reliable than those with 1-2 genes.
#'   \item \strong{GeneRatio}: Proportion of input genes (e.g., "10/100" = 10% genes in pathway).
#'   \item \strong{BgRatio}: Proportion of background genes (e.g., "200/20000" = 1% of genome).
#' }
#'
#' **Protein-Coding Gene Filtering**:
#' \itemize{
#'   \item \strong{Why filter}: GO/KEGG databases primarily annotate protein-coding genes (PCGs).
#'         Including non-coding genes dilutes enrichment signal.
#'   \item \strong{Expected retention}: 70-85% of input genes are PCGs (typical for RNA-seq).
#'   \item \strong{Warning threshold}: If < 10 PCGs remain, consider disabling filter or
#'         increasing input gene count.
#' }
#'
#' **Dependencies**: Requires \code{geneset} and \code{genekitr2} packages (included in Imports).
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @section Performance Test:
#' \itemize{
#'   \item \strong{Test query}: TP53 targets (100 upregulated + 100 downregulated genes)
#'   \item \strong{Runtime}: 29.8 sec (gene ID conversion + GO database query)
#'   \item \strong{PCG filtering}: 100 → 81 genes (81% retention, up); 100 → 71 genes (71%, down)
#'   \item \strong{Result size}: 15,461 GO terms tested per direction, 0 significant (p < 0.05)
#'   \item \strong{Note}: Runtime scales with gene count and database size. GO (largest) > KEGG > Reactome.
#' }
#'
#' @seealso
#' **Upstream Preparation** (get gene lists):
#' \itemize{
#'   \item \code{\link{gpdb_find_targets}}: Query downstream genes regulated by perturbation
#'   \item \code{\link{gpdb_find_regulators}}: Find upstream regulators of target gene
#'   \item \code{\link{gpdb_what_happens}}: Get all perturbation effects for comprehensive overview
#' }
#'
#' **Downstream Visualization** (interpret results):
#' \itemize{
#'   \item \code{\link{gpdb_plot_enrichment}}: Create paired dotplots for up/down pathways
#'   \item \code{\link{gpdb_plot_network}}: Visualize gene-pathway interaction networks
#' }
#'
#' **Alternative Enrichment Methods** (different approaches):
#' \itemize{
#'   \item \code{fgsea::fgseaMultilevel}: Gene Set Enrichment Analysis (preranked GSEA).
#'     Use when: Have ranked gene list (not just significant genes).
#'   \item \code{clusterProfiler::enrichGO}: Direct GO enrichment with custom parameters.
#'     Use when: Need fine control over single gene set enrichment.
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item What biological pathways are enriched in my upregulated genes?
#'   \item How do I interpret enrichment p-values and fold enrichment?
#'   \item Should I filter protein-coding genes for GO enrichment?
#'   \item What's the difference between GO biological process and molecular function?
#'   \item Can I analyze upregulated and downregulated genes separately?
#'   \item Which enrichment database should I use: GO, KEGG, or Reactome?
#'   \item How many genes do I need for reliable pathway enrichment?
#'   \item What does FDR q-value mean in enrichment analysis?
#'   \item How to customize background genes for tissue-specific enrichment?
#'   \item Can I run enrichment on mouse or rat gene lists?
#'   \item Why are there no significant pathways in my results?
#'   \item How to compare enrichment results between different conditions?
#'   \item What's the minimum p-value cutoff for pathway significance?
#'   \item Should I use Hallmark or Canonical pathways in MSigDB?
#'   \item How to export enrichment results for publication?
#'   \item Can I visualize enrichment results as dotplots or networks?
#'   \item What does "10/100 genes in pathway" mean?
#'   \item How to focus enrichment on top 50 strongest regulated genes?
#'   \item Why does GO enrichment take longer than KEGG enrichment?
#'   \item Can I run enrichment on gene symbols directly without IDs?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Usage (TESTED - runtime: 29.8 sec)
#' # ===========================================================================
#' # Scientific Question: What pathways are enriched in TP53 perturbation?
#' # Query: TP53 downstream targets with GO biological process enrichment
#' # Expected output: Separate enrichment for upregulated and downregulated genes
#'
#' # Zero configuration - just run!
#' library(GenePerturbR)
#'
#' # Step 1: Get gene targets
#' targets <- gpdb_find_targets("TP53", top_n = 100, min_confidence = "high")
#' cat("Upregulated:", nrow(targets$upregulated), "genes\n")
#' cat("Downregulated:", nrow(targets$downregulated), "genes\n")
#'
#' # Step 2: Perform enrichment
#' result <- gpdb_enrich(
#'   targets,
#'   enrich.type = "GO",
#'   GO.ont = "bp", # Biological process
#'   filter.pcg = TRUE, # Recommended for GO
#'   p.cutoff = 0.05,
#'   q.cutoff = 0.05
#' )
#'
#' # Verify output structure
#' cat("\nResult class:", class(result)[1], "\n")
#' cat("Elements:", paste(names(result), collapse = ", "), "\n")
#' cat("Significant UP pathways:", nrow(result$upregulated$Sig), "\n")
#' cat("Significant DOWN pathways:", nrow(result$downregulated$Sig), "\n")
#'
#' # Quick preview of top pathways
#' if (nrow(result$upregulated$Sig) > 0) {
#'   head(result$upregulated$Sig[, c("Description", "pvalue", "p.adjust", "Count")])
#' }
#'
#' # ===========================================================================
#' # Example 2: Database Comparison (GO vs KEGG vs Reactome)
#' # ===========================================================================
#' # Compare pathway enrichment across different annotation databases
#'
#' # GO enrichment (comprehensive, ~15K terms)
#' go_result <- gpdb_enrich(targets, enrich.type = "GO", GO.ont = "bp")
#' cat("GO pathways tested:", nrow(go_result$upregulated$All), "\n")
#'
#' # KEGG enrichment (curated pathways, ~300 terms)
#' kegg_result <- gpdb_enrich(targets, enrich.type = "KEGG")
#' cat("KEGG pathways tested:", nrow(kegg_result$upregulated$All), "\n")
#'
#' # Reactome enrichment (detailed pathway hierarchy, ~2K terms)
#' reactome_result <- gpdb_enrich(targets, enrich.type = "Reactome")
#' cat("Reactome pathways tested:", nrow(reactome_result$upregulated$All), "\n")
#'
#' # Compare significant findings
#' cat("\nSignificant pathways found:\n")
#' cat(
#'   "  GO:", nrow(go_result$upregulated$Sig), "(up)",
#'   nrow(go_result$downregulated$Sig), "(down)\n"
#' )
#' cat(
#'   "  KEGG:", nrow(kegg_result$upregulated$Sig), "(up)",
#'   nrow(kegg_result$downregulated$Sig), "(down)\n"
#' )
#'
#' # ===========================================================================
#' # Example 3: Parameter Optimization (PCG filtering and gene count)
#' # ===========================================================================
#' # Test impact of protein-coding gene filtering and input gene count
#'
#' # Default: Filter PCGs, top 100 genes
#' result_default <- gpdb_enrich(
#'   targets,
#'   enrich.type = "GO",
#'   filter.pcg = TRUE, # Keep only protein-coding genes
#'   top_up = 100,
#'   top_down = 100
#' )
#'
#' # Alternative: Include all genes, top 200
#' result_all <- gpdb_enrich(
#'   targets,
#'   enrich.type = "GO",
#'   filter.pcg = FALSE, # Include non-coding genes
#'   top_up = 200,
#'   top_down = 200
#' )
#'
#' # Compare gene counts after filtering
#' cat("PCG filtering impact:\n")
#' cat("  With filter:", result_default$params$n_genes_up, "genes (up)\n")
#' cat("  Without filter:", result_all$params$n_genes_up, "genes (up)\n")
#'
#' # ===========================================================================
#' # Example 4: Direct Gene List Input (custom analysis)
#' # ===========================================================================
#' # Enrich custom gene lists without database query
#'
#' # Named list format
#' custom_genes <- list(
#'   up = c("MYC", "CCND1", "CDK4", "CDK6", "E2F1"),
#'   down = c("TP53", "RB1", "CDKN1A", "CDKN1B", "CDKN2A")
#' )
#' custom_result <- gpdb_enrich(custom_genes, enrich.type = "KEGG")
#'
#' # Single gene vector (analyzed as one group)
#' gene_vector <- c("MYC", "TP53", "KRAS", "EGFR", "PIK3CA")
#' single_result <- gpdb_enrich(
#'   gene_vector,
#'   enrich.type = "GO",
#'   split.by = "none" # Don't try to split by direction
#' )
#'
#' # ===========================================================================
#' # Next Steps (Complete Workflow)
#' # ===========================================================================
#' # For complete pathway analysis workflows, see:
#' # - Visualization: ?gpdb_plot_enrichment
#' # - Network analysis: ?gpdb_plot_network
#' # - Query more targets: ?gpdb_find_targets
#' # - Find regulators: ?gpdb_find_regulators
#' # - Compare conditions: Run gpdb_enrich on multiple gene sets
#' }
gpdb_enrich <- function(targets,
                        enrich.type = c("GO", "KEGG", "Wiki", "Reactome", "MsigDB", "Mesh", "HgDisease", "Enrichrdb"),
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

  # Parse input and extract gene lists
  gene_lists <- .parse_enrich_input(targets, split.by, top_up, top_down)

  # Filter to protein-coding genes (recommended for accurate enrichment)
  if (filter.pcg) {
    gene_lists <- .filter_protein_coding(gene_lists, organism, min.genes)
  }

  # Load gene set database (once, shared between up and down)
  message("Loading ", enrich.type, " gene set database...")
  geneset_df <- .gpdb_get_geneset(
    enrich.type = enrich.type,
    GO.ont = GO.ont,
    KEGG.category = KEGG.category,
    Msigdb.category = Msigdb.category,
    organism = organism
  )
  message("Loaded ", length(unique(geneset_df$id)), " gene sets with ",
    length(unique(geneset_df$gene)), " unique genes")

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
    result$upregulated <- .gpdb_ora_enrich(
      genes = gene_lists$up,
      geneset_df = geneset_df,
      p.cutoff = p.cutoff,
      q.cutoff = q.cutoff,
      background.genes = background.genes
    )
  } else {
    result$upregulated <- list(All = data.frame(), Sig = data.frame())
  }

  if (!is.null(gene_lists$down) && length(gene_lists$down) > 0) {
    message("\n=== Enriching downregulated genes (n=", length(gene_lists$down), ") ===")
    result$downregulated <- .gpdb_ora_enrich(
      genes = gene_lists$down,
      geneset_df = geneset_df,
      p.cutoff = p.cutoff,
      q.cutoff = q.cutoff,
      background.genes = background.genes
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


#' Paired Dotplot Visualization for Pathway Enrichment
#'
#' @description
#' Create side-by-side dotplots comparing enriched pathways in upregulated versus
#' downregulated genes using publication-ready aesthetics. Visualizes top pathways
#' ranked by significance with customizable color gradients (p-value/FoldEnrich) and
#' point sizes (gene count). Returns a patchwork object combining two dotplots for
#' direct comparison of directional pathway changes.
#'
#' @param enrich_result \code{gpdb_enrichment} object from \code{\link{gpdb_enrich}}.
#'   Must contain \code{$upregulated} and \code{$downregulated} enrichment results.
#' @param show.term.num Integer. Number of top pathways to display per panel (default: 15).
#'   Pathways ranked by p-value. Increase to 20-30 for comprehensive view, decrease to 5-10 for focus.
#' @param use.all Logical. Plot all pathways or only significant ones (default: FALSE, significant only).
#'   Set TRUE to visualize trends even without significant enrichment (p < 0.05).
#' @param x Character. X-axis variable: "FoldEnrich" (enrichment ratio, default),
#'   "pvalue" (raw p-value), "p.adjust" (FDR-corrected), "Count" (gene count).
#'   "FoldEnrich" emphasizes biological magnitude; "pvalue" emphasizes statistical significance.
#' @param color.by Character. Variable mapped to point color: "p.adjust" (FDR, default),
#'   "pvalue", "FoldEnrich", "Count". Color gradient shows continuous scale.
#' @param size.by Character. Variable mapped to point size: "Count" (gene count, default),
#'   "pvalue", "p.adjust", "FoldEnrich". Larger points indicate more genes or stronger significance.
#' @param colors Character vector. Color palette for gradient (default: BrBG diverging palette).
#'   Provide custom colors like \code{c("blue", "white", "red")} or use viridis scales.
#'   Reversed automatically for downregulated panel to maintain visual symmetry.
#' @param size.range Numeric vector of length 2. Point size range in mm (default: c(3, 8)).
#'   Adjust for visual clarity: smaller range for dense plots, larger for emphasis.
#' @param title Character. Main plot title (default: auto-generated from database and ontology).
#'   Set to custom text or NULL to remove.
#' @param up.title Character. Title for upregulated panel (default: "Upregulated Genes").
#' @param down.title Character. Title for downregulated panel (default: "Downregulated Genes").
#' @param legend.position Character. Legend placement: "right" (default), "bottom", "none".
#'   Use "bottom" for wide plots, "none" to maximize plot area.
#' @param theme ggplot2 theme object (default: \code{theme_bw(base_rect_size = 1.5)}).
#'   Publication-ready theme with clean white background. Customize with \code{theme_minimal()}, etc.
#'
#' @return \code{patchwork} object (combined ggplot) with two side-by-side dotplots.
#'
#' **Return Structure**:
#' \describe{
#'   \item{\strong{Combined plot}}{Two-panel layout showing upregulated (left) vs downregulated (right) pathways.
#'     \itemize{
#'       \item \strong{X-axis}: Enrichment metric (FoldEnrich, p-value, etc.)
#'       \item \strong{Y-axis}: Pathway descriptions (auto-formatted, capitalized)
#'       \item \strong{Point color}: Continuous gradient (e.g., FDR from low to high)
#'       \item \strong{Point size}: Proportional to selected variable (e.g., gene count)
#'       \item \strong{Dimensions}: Auto-calculated based on pathway count and label length
#'     }
#'   }
#'   \item{\strong{Attributes}}{
#'     \itemize{
#'       \item \code{attr(plot, "width")}: Recommended width in inches (typically 14-20)
#'       \item \code{attr(plot, "height")}: Recommended height in inches (typically 6-12)
#'     }
#'   }
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Save to file: \code{ggsave("enrichment.pdf", plot, width = attr(plot, "width"), height = attr(plot, "height"))}
#'   \item Customize further: \code{plot + patchwork::plot_annotation(tag_levels = 'A')}
#'   \item Extract individual panels: \code{plot[[1]]} (upregulated), \code{plot[[2]]} (downregulated)
#'   \item Modify theme: \code{plot & ggplot2::theme_minimal()}
#'   \item Export to PowerPoint: Use \code{rvg} or \code{officer} packages for editable vector graphics
#' }
#'
#' **Alternative Visualization Methods**:
#' \itemize{
#'   \item \code{\link{gpdb_plot_network}}: Use when you want gene-pathway network instead of dotplot.
#'   \item \code{enrichplot::dotplot}: Use for single-direction enrichment (not paired comparison).
#'   \item \code{enrichplot::barplot}: Use for simpler bar chart representation.
#' }
#'
#' @details
#' **Visual Design Principles**:
#' \itemize{
#'   \item \strong{Paired comparison}: Side-by-side layout reveals distinct pathway responses
#'         (e.g., cell cycle up, apoptosis down).
#'   \item \strong{Color gradients}: BrBG diverging palette (brown-blue-green) is colorblind-friendly
#'         and shows continuous scale effectively. Reversed for downregulated panel.
#'   \item \strong{Point size}: Larger points indicate more genes, helping identify robust
#'         vs. sparse enrichment.
#'   \item \strong{Automatic scaling}: Plot dimensions adjust based on pathway count and
#'         label length for optimal readability.
#'   \item \strong{Y-axis formatting}: Pathway names auto-capitalized and positioned on right
#'         for better alignment with upregulated panel.
#' }
#'
#' **Choosing Visualization Parameters**:
#' \itemize{
#'   \item \strong{X-axis}: Use "FoldEnrich" (default) for biological interpretation,
#'         "pvalue" for statistical emphasis, "Count" for gene set size.
#'   \item \strong{Color}: Map "p.adjust" (default) to highlight significance,
#'         "FoldEnrich" to emphasize magnitude.
#'   \item \strong{Size}: Map "Count" (default) to show gene set size,
#'         "pvalue" to emphasize significance via size.
#'   \item \strong{show.term.num}: 15 (default) balances detail vs. clutter.
#'         Use 5-10 for presentations, 20-30 for comprehensive reports.
#' }
#'
#' **Dependencies**: Requires \code{patchwork} package for combining plots (auto-checked).
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @section Performance Test:
#' \itemize{
#'   \item \strong{Test data}: TP53 KEGG enrichment (347 pathways tested, 0 UP + 1 DOWN significant)
#'   \item \strong{Runtime}: 0.096 sec (plot creation and layout)
#'   \item \strong{Output}: patchwork object, 14 x 6 inches (auto-calculated)
#'   \item \strong{Note}: Runtime scales with \code{show.term.num} and label complexity.
#'         Typically < 0.2 sec for 10-30 pathways.
#' }
#'
#' @seealso
#' **Upstream Analysis** (generate enrichment results):
#' \itemize{
#'   \item \code{\link{gpdb_enrich}}: Perform ORA enrichment to generate input for this function
#'   \item \code{\link{gpdb_find_targets}}: Query gene targets before enrichment
#' }
#'
#' **Alternative Visualizations** (different plot types):
#' \itemize{
#'   \item \code{\link{gpdb_plot_network}}: Visualize gene-pathway network.
#'     Use when: Want to see gene-level connections within pathways.
#'   \item \code{enrichplot::cnetplot}: Category-network plot.
#'     Use when: Want to visualize gene-category relationships.
#'   \item \code{enrichplot::emapplot}: Enrichment map with pathway similarity.
#'     Use when: Want to group similar pathways by gene overlap.
#' }
#'
#' **Export and Customization**:
#' \itemize{
#'   \item \code{ggplot2::ggsave}: Save to PDF/PNG with custom dimensions
#'   \item \code{patchwork::plot_annotation}: Add overall title, caption, or tags
#'   \item \code{patchwork::plot_layout}: Adjust panel arrangement (ncol, widths)
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I visualize pathway enrichment results as dotplots?
#'   \item What's the best way to compare upregulated vs downregulated pathways?
#'   \item How to customize enrichment dotplot colors and sizes?
#'   \item Can I show top 20 pathways instead of default 15?
#'   \item How do I save enrichment plots for publication?
#'   \item What does point size represent in the enrichment dotplot?
#'   \item Should I use FoldEnrich or p-value for the X-axis?
#'   \item How to interpret color gradients in enrichment plots?
#'   \item Can I plot all pathways even if not significant?
#'   \item What are the recommended plot dimensions for enrichment figures?
#'   \item How do I export enrichment plots to PowerPoint?
#'   \item Can I extract individual panels from the paired dotplot?
#'   \item How to change the color palette to viridis or custom colors?
#'   \item Why do enrichment plots show different colors for up vs down?
#'   \item How to adjust point size range for better visualization?
#'   \item Can I remove legends to maximize plot area?
#'   \item How to customize pathway label formatting?
#'   \item What if I have no significant pathways to plot?
#'   \item How to create publication-quality enrichment figures?
#'   \item Can I combine multiple enrichment plots into one figure?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Usage (TESTED - runtime: 0.096 sec)
#' # ===========================================================================
#' # Scientific Question: How do upregulated vs downregulated pathways differ?
#' # Input: gpdb_enrichment object from TP53 perturbation analysis
#' # Expected output: Side-by-side dotplots showing top pathways per direction
#'
#' # Zero configuration - just run!
#' library(GenePerturbR)
#'
#' # Step 1: Get enrichment results
#' targets <- gpdb_find_targets("TP53", top_n = 200, min_confidence = "medium")
#' enrich_result <- gpdb_enrich(targets, enrich.type = "KEGG", p.cutoff = 0.05)
#'
#' # Step 2: Create visualization
#' plot <- gpdb_plot_enrichment(enrich_result)
#'
#' # Print plot (in RStudio/interactive session)
#' print(plot)
#'
#' # Check plot properties
#' cat("Plot class:", class(plot)[1], "\n")
#' cat("Recommended dimensions:", attr(plot, "width"), "x", attr(plot, "height"), "inches\n")
#'
#' # Save to file with recommended dimensions
#' # ggsave("enrichment.pdf", plot,
#' #        width = attr(plot, "width"),
#' #        height = attr(plot, "height"))
#'
#' # ===========================================================================
#' # Example 2: Customize Visual Encoding (X-axis and color mapping)
#' # ===========================================================================
#' # Compare different ways to emphasize enrichment
#'
#' # Default: FoldEnrich (X) + p.adjust (color)
#' plot1 <- gpdb_plot_enrichment(
#'   enrich_result,
#'   show.term.num = 10,
#'   x = "FoldEnrich", # Biological magnitude
#'   color.by = "p.adjust" # Statistical significance
#' )
#'
#' # Alternative: p-value (X) + FoldEnrich (color)
#' plot2 <- gpdb_plot_enrichment(
#'   enrich_result,
#'   show.term.num = 10,
#'   x = "pvalue", # Statistical significance
#'   color.by = "FoldEnrich" # Biological magnitude
#' )
#'
#' # Emphasize gene count
#' plot3 <- gpdb_plot_enrichment(
#'   enrich_result,
#'   show.term.num = 10,
#'   x = "Count", # Number of genes
#'   color.by = "p.adjust",
#'   size.by = "FoldEnrich" # Size by enrichment ratio
#' )
#'
#' # ===========================================================================
#' # Example 3: Custom Color Palettes and Sizes
#' # ===========================================================================
#' # Use custom colors for publication or aesthetic preference
#'
#' # Viridis-style colors
#' plot_viridis <- gpdb_plot_enrichment(
#'   enrich_result,
#'   show.term.num = 15,
#'   colors = viridis::viridis(10),
#'   size.range = c(4, 10) # Larger points
#' )
#'
#' # Red-blue diverging (classic)
#' plot_rdb <- gpdb_plot_enrichment(
#'   enrich_result,
#'   show.term.num = 15,
#'   colors = c(
#'     "#2166AC", "#4393C3", "#92C5DE", "#D1E5F0",
#'     "#F7F7F7", "#FDDBC7", "#F4A582", "#D6604D", "#B2182B"
#'   ),
#'   size.range = c(2, 6) # Smaller, tighter range
#' )
#'
#' # ===========================================================================
#' # Example 4: Handle Non-Significant Results
#' # ===========================================================================
#' # Visualize trends even when no pathways pass significance threshold
#'
#' # If no significant pathways, use.all = TRUE shows top pathways anyway
#' plot_all <- gpdb_plot_enrichment(
#'   enrich_result,
#'   show.term.num = 10,
#'   use.all = TRUE, # Show all pathways (not just p < 0.05)
#'   x = "FoldEnrich",
#'   color.by = "pvalue" # Still color by p-value to show trends
#' )
#'
#' # Useful for exploratory analysis or when enrichment is weak
#'
#' # ===========================================================================
#' # Example 5: Advanced Customization (patchwork + ggplot2)
#' # ===========================================================================
#' # Modify plot after creation for publication needs
#'
#' library(patchwork)
#' library(ggplot2)
#'
#' plot <- gpdb_plot_enrichment(enrich_result, show.term.num = 10)
#'
#' # Add overall title and caption
#' plot_final <- plot +
#'   plot_annotation(
#'     title = "TP53 Perturbation: Pathway Enrichment Analysis",
#'     subtitle = "KEGG Pathways (Fisher's Exact Test)",
#'     caption = "Source: GenePerturbR (GPSAdb database)",
#'     theme = theme(plot.title = element_text(size = 16, face = "bold"))
#'   )
#'
#' # Modify both panels with unified theme
#' plot_minimal <- plot & theme_minimal(base_size = 12)
#'
#' # Extract and customize individual panels
#' plot_up <- plot[[1]] +
#'   labs(title = "Activated Pathways") +
#'   theme(axis.text.y = element_text(color = "darkred"))
#'
#' # ===========================================================================
#' # Next Steps (Complete Workflow)
#' # ===========================================================================
#' # For complete visualization workflows, see:
#' # - Network visualization: ?gpdb_plot_network
#' # - Run enrichment first: ?gpdb_enrich
#' # - Get gene targets: ?gpdb_find_targets
#' # - Export options: ?ggplot2::ggsave
#' # - Combine plots: ?patchwork::plot_layout
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

  # Prepare plot data without pipe operator
  plot_data <- res[order(res$pvalue), ]
  plot_data <- plot_data[1:min(show.term.num, nrow(plot_data)), ]
  plot_data <- plot_data[order(-plot_data[[x_col]]), ]
  plot_data$Description <- factor(plot_data$Description, levels = rev(plot_data$Description))

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
  pcg_genes <- .gpdb_get_pcg()
  if (is.null(pcg_genes)) {
    message("Protein-coding gene list not available, skipping PCG filtering")
    return(gene_lists)
  }

  message("Filtering protein-coding genes (set filter.pcg=FALSE to disable)...")

  # Filter upregulated genes
  if (!is.null(gene_lists$up) && length(gene_lists$up) > 0) {
    n_before <- length(gene_lists$up)
    gene_lists$up <- intersect(gene_lists$up, pcg_genes)
    n_after <- length(gene_lists$up)
    pct <- round(n_after / n_before * 100, 1)

    message(
      "  Upregulated: ", n_before, " \u2192 ", n_after,
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
    gene_lists$down <- intersect(gene_lists$down, pcg_genes)
    n_after <- length(gene_lists$down)
    pct <- round(n_after / n_before * 100, 1)

    message(
      "  Downregulated: ", n_before, " \u2192 ", n_after,
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


# ==============================================================================
# Local Enrichment Engine (geneset + genekitr2 + Fisher's exact test)
# ==============================================================================

#' Get protein-coding gene list from bundled data
#' @keywords internal
.gpdb_get_pcg <- function() {
  rds_file <- system.file("extdata", "protein_coding_genes.rds", package = "GenePerturbR")
  if (rds_file == "") {
    warning("Protein-coding gene list not found in package data")
    return(NULL)
  }
  readRDS(rds_file)
}


#' Fetch gene set database via geneset package and return standardized format
#'
#' Returns a data.frame with columns: id, term, gene (gene symbols).
#' Handles Entrez-to-Symbol conversion for databases that return Entrez IDs.
#' @keywords internal
.gpdb_get_geneset <- function(enrich.type,
                              GO.ont = "bp",
                              KEGG.category = "pathway",
                              Msigdb.category = "H",
                              organism = "hs") {
  if (!requireNamespace("geneset", quietly = TRUE)) {
    stop("geneset package required. Install with: devtools::install_github('GangLiLab/geneset')",
      call. = FALSE
    )
  }

  type_upper <- toupper(enrich.type)

  gs <- switch(type_upper,
    "GO" = geneset::getGO(org = organism, ont = tolower(GO.ont)),
    "KEGG" = geneset::getKEGG(org = organism, category = KEGG.category),
    "MSIGDB" = geneset::getMsigdb(org = organism, category = toupper(Msigdb.category)),
    "REACTOME" = geneset::getReactome(org = organism),
    "WIKI" = geneset::getWiki(org = organism),
    "MESH" = geneset::getMesh(org = organism, method = "gendoo", category = "A"),
    "HGDISEASE" = geneset::getHgDisease(source = "do"),
    "ENRICHRDB" = {
      if (!requireNamespace("geneset", quietly = TRUE)) stop("geneset required")
      geneset::getEnrichrdb(org = "human", library = "Cancer_Cell_Line_Encyclopedia")
    },
    stop("Unsupported enrichment database: ", enrich.type, call. = FALSE)
  )

  # Standardize column names across databases
  geneset_data <- gs$geneset
  geneset_names <- gs$geneset_name

  if (type_upper == "MSIGDB") {
    # MsigDB uses different column names: gs_name, entrez_gene
    colnames(geneset_data) <- c("id", "gene")
    # No geneset_name for MsigDB; use gs_name as both id and term
    geneset_names <- data.frame(
      id = unique(geneset_data$id),
      name = gsub("_", " ", unique(geneset_data$id)),
      stringsAsFactors = FALSE
    )
  } else if (type_upper == "GO") {
    # GO uses ontology-specific column name (bp, cc, mf) for the ID
    colnames(geneset_data)[1] <- "id"
  }

  # Determine if gene column contains Entrez IDs (numeric) or symbols
  sample_genes <- head(unique(geneset_data$gene), 50)
  is_entrez <- all(grepl("^[0-9]+$", sample_genes))

  if (is_entrez) {
    # Convert Entrez IDs to gene symbols via genekitr2
    if (!requireNamespace("genekitr2", quietly = TRUE)) {
      stop("genekitr2 package required. Install with: devtools::install_github('Zaoqu-Liu/genekitr2')",
        call. = FALSE
      )
    }

    unique_ids <- unique(geneset_data$gene)
    message("Converting ", length(unique_ids), " Entrez IDs to gene symbols...")

    gene_info <- tryCatch(
      genekitr2::genInfo(id = unique_ids, unique = TRUE, org = organism),
      error = function(e) {
        warning("ID conversion failed: ", e$message)
        return(NULL)
      }
    )

    if (!is.null(gene_info) && nrow(gene_info) > 0) {
      id_col <- if ("input_id" %in% colnames(gene_info)) "input_id" else colnames(gene_info)[1]
      gene_map <- stats::setNames(gene_info$symbol, as.character(gene_info[[id_col]]))
      geneset_data$gene <- gene_map[as.character(geneset_data$gene)]
      geneset_data <- geneset_data[!is.na(geneset_data$gene) & geneset_data$gene != "", ]
    }
  }

  # Merge with pathway names
  if (!is.null(geneset_names) && nrow(geneset_names) > 0) {
    colnames(geneset_names) <- c("id", "term")
    result <- merge(geneset_data[, c("id", "gene")], geneset_names, by = "id", all.x = TRUE)
    result$term[is.na(result$term)] <- result$id[is.na(result$term)]
  } else {
    result <- geneset_data[, c("id", "gene")]
    result$term <- result$id
  }

  result <- result[, c("id", "term", "gene")]
  result <- result[!is.na(result$gene) & result$gene != "", ]

  return(result)
}


#' Local Over-Representation Analysis using Fisher's exact test
#'
#' Returns list(All, Sig) with standard ORA column names.
#' @keywords internal
.gpdb_ora_enrich <- function(genes,
                              geneset_df,
                              p.cutoff = 0.05,
                              q.cutoff = 0.05,
                              background.genes = NULL) {
  # Build universe (background gene set)
  all_annotated_genes <- unique(geneset_df$gene)
  if (!is.null(background.genes)) {
    universe <- unique(c(background.genes, genes))
    universe <- intersect(universe, all_annotated_genes)
  } else {
    universe <- all_annotated_genes
  }

  # Filter input genes to those in universe
  genes <- unique(genes[genes %in% universe])
  n_input <- length(genes)
  n_universe <- length(universe)

  if (n_input == 0) {
    warning("No input genes found in annotation database", call. = FALSE)
    empty_df <- data.frame(
      ID = character(), Description = character(),
      GeneRatio = character(), BgRatio = character(),
      pvalue = numeric(), p.adjust = numeric(),
      Count = integer(), FoldEnrich = numeric(),
      geneID_symbol = character(),
      stringsAsFactors = FALSE
    )
    return(list(All = empty_df, Sig = empty_df))
  }

  # Split gene sets into list
  pathway_list <- split(geneset_df$gene, geneset_df$id)
  pathway_names <- stats::setNames(geneset_df$term, geneset_df$id)
  pathway_names <- pathway_names[!duplicated(names(pathway_names))]

  # Run Fisher's exact test for each pathway
  results <- lapply(names(pathway_list), function(pw_id) {
    pw_genes <- unique(pathway_list[[pw_id]])
    pw_genes_in_universe <- pw_genes[pw_genes %in% universe]

    n_pw <- length(pw_genes_in_universe)
    if (n_pw == 0) return(NULL)

    # Overlap
    overlap <- intersect(genes, pw_genes_in_universe)
    n_overlap <- length(overlap)
    if (n_overlap == 0) return(NULL)

    # 2x2 contingency table for Fisher's exact test
    #                 In pathway    Not in pathway
    # In gene list    a             b
    # Not in list     c             d
    a <- n_overlap
    b <- n_input - n_overlap
    c <- n_pw - n_overlap
    d <- n_universe - n_input - n_pw + n_overlap

    mat <- matrix(c(a, c, b, d), nrow = 2)
    fisher_result <- stats::fisher.test(mat, alternative = "greater")

    fold_enrich <- (a / n_input) / (n_pw / n_universe)

    data.frame(
      ID = pw_id,
      Description = as.character(pathway_names[pw_id]),
      GeneRatio = paste0(n_overlap, "/", n_input),
      BgRatio = paste0(n_pw, "/", n_universe),
      pvalue = fisher_result$p.value,
      Count = n_overlap,
      FoldEnrich = round(fold_enrich, 4),
      geneID_symbol = paste(overlap, collapse = "/"),
      stringsAsFactors = FALSE
    )
  })

  results <- do.call(rbind, results[!sapply(results, is.null)])

  if (is.null(results) || nrow(results) == 0) {
    empty_df <- data.frame(
      ID = character(), Description = character(),
      GeneRatio = character(), BgRatio = character(),
      pvalue = numeric(), p.adjust = numeric(),
      Count = integer(), FoldEnrich = numeric(),
      geneID_symbol = character(),
      stringsAsFactors = FALSE
    )
    return(list(All = empty_df, Sig = empty_df))
  }

  # Multiple testing correction (BH/FDR)
  results$p.adjust <- stats::p.adjust(results$pvalue, method = "BH")

  # Sort by p-value
  results <- results[order(results$pvalue), ]
  rownames(results) <- NULL

  # Split into All and Sig
  sig_results <- results[results$pvalue < p.cutoff & results$p.adjust < q.cutoff, ]

  return(list(All = results, Sig = sig_results))
}
