#' Generate Natural Language Summary of Gene Function from Perturbation Data
#'
#' Create automated gene function summaries by aggregating perturbation evidence
#' from 7,665 RNA-seq datasets. Returns human-readable text describing experimental
#' coverage (cell lines, tissues) and top regulatory targets with effect sizes.
#' **Database covers 2,810 genes across 1,063 cell lines** with summaries generated
#' from multi-dataset consensus. Ideal for quick gene profile generation, literature
#' review assistance, and LLM-based analysis workflows.
#'
#' @param gene Character. Gene symbol (e.g., "TP53", "METTL3"). Must exist in database
#'   (use \code{gpdb_search()} to verify). Case-insensitive but uppercase recommended.
#' @param max_length Integer. Maximum summary length in characters.
#'   Default: \code{500}. Range: 100-2000. Summaries exceeding limit are truncated
#'   with "..." suffix. Use 200-300 for brief overviews, 500+ for detailed profiles.
#' @param include Character vector. Content aspects to include in summary.
#'   Default: \code{c("function", "targets", "pathways")}. Options: "function"
#'   (experimental coverage), "targets" (top regulated genes), "pathways"
#'   (enrichment results - not yet implemented), "diseases" (clinical associations - not yet implemented).
#'   Currently only "function" and "targets" are active.
#'
#' @return Character string with natural language gene summary. Format: "[Gene] has been
#'   studied in N experiments across M cell lines and K tissue types. Perturbation leads
#'   to upregulation of [genes] (X-fold) and downregulation of [genes] (Y-fold)."
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Get detailed targets: \code{\link{gpdb_find_targets}(gene, top_n = 100)}
#'   \item Load raw datasets: \code{\link{gpdb_list_datasets}(gene)} then \code{\link{gpdb_load_data}(dataset_id)}
#'   \item Enrichment analysis: \code{\link{gpdb_enrich}(targets, enrich.type = "GO")}
#'   \item Network visualization: \code{\link{gpdb_plot_network}(gene, top_targets = 50)}
#'   \item Compare with other genes: \code{\link{gpdb_compare_genes}(c(gene, "GENE2"))}
#' }
#'
#' **Alternative Methods**:
#' \itemize{
#'   \item \code{\link{gpdb_what_happens}}: More structured output with detailed statistics.
#'     Use when you need programmatic access to data instead of text summary.
#'   \item \code{\link{gpdb_find_targets}}: Get full target list with confidence scores.
#'     Use when you need quantitative analysis rather than qualitative summary.
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test query**: TP53 (comprehensive gene profile)
#'   \item **Runtime**: 0.020 sec (database query + text generation)
#'   \item **Result size**: 224 characters (71 experiments, 39 cell lines, 23 tissues)
#'   \item **Database**: 18.9M relationships, optimized with SQLite indexing
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Data Query** (explore database):
#' \itemize{
#'   \item \code{\link{gpdb_search}}: Find gene by symbol or alias
#'   \item \code{\link{gpdb_list_datasets}}: List all experiments for a gene
#'   \item \code{\link{gpdb_get_info}}: Get detailed metadata for specific dataset
#'   \item \code{\link{gpdb_load_data}}: Load raw expression data
#' }
#'
#' **Downstream Analysis** (sequential use):
#' \itemize{
#'   \item \code{\link{gpdb_find_targets}}: Get detailed target list with confidence scores
#'   \item \code{\link{gpdb_enrich}}: Pathway enrichment on regulated genes
#'   \item \code{\link{gpdb_plot_network}}: Visualize regulatory network
#'   \item \code{\link{gpdb_analyze_cascade}}: Trace multi-step regulatory cascades
#' }
#'
#' **Related Query Functions** (parallel comparison):
#' \itemize{
#'   \item \code{\link{gpdb_what_happens}}: Structured output with statistics instead of text.
#'     Use when: Need programmatic access to perturbation effects.
#'   \item \code{\link{gpdb_find_regulators}}: Find upstream genes that regulate target.
#'     Use when: Want to know what controls gene expression.
#'   \item \code{\link{gpdb_compare_genes}}: Compare perturbation effects across genes.
#'     Use when: Analyzing multiple genes simultaneously.
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I get a quick overview of what a gene does?
#'   \item What tissues and cell lines have been studied for my gene?
#'   \item Which genes are most strongly affected by gene perturbation?
#'   \item How many experiments are available for my gene of interest?
#'   \item Can I generate a natural language summary for literature review?
#'   \item What effect sizes do I see for top regulated targets?
#'   \item How do I create brief gene profiles for multiple genes?
#'   \item Which genes are upregulated vs downregulated by perturbation?
#'   \item Can I use this summary in LLM-based analysis pipelines?
#'   \item How reliable is the summary if there are few experiments?
#'   \item What information is included in the automated summary?
#'   \item Can I control the length of the generated text?
#'   \item How do I interpret fold change values in the summary?
#'   \item What happens if my gene has no perturbation data?
#'   \item Can I compare summaries across related genes?
#'   \item How do I export summaries for batch analysis?
#'   \item Which aspects of gene function are covered?
#'   \item Can I get summaries for non-human genes?
#'   \item How often is the summary data updated?
#'   \item What downstream analyses work with summary results?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Usage (TESTED - runtime: 0.020 sec)
#' # ===========================================================================
#' # Scientific Question: Quick overview of TP53 function from perturbation data
#' # Expected output: Natural language text with experimental coverage and top targets
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' summary <- gpdb_summarize("TP53")
#' cat(summary)
#' # Output: "TP53 has been studied in 71 experiments across 39 cell lines
#' # and 23 tissue types. Perturbation leads to upregulation of SLC10A6
#' # (7.82-fold) and downregulation of MUC17, LINC01482, SPRR3..."
#'
#' # Check summary properties
#' cat("Summary length:", nchar(summary), "characters\n")
#'
#' # ===========================================================================
#' # Example 2: Control Summary Length (TESTED - runtime: 0.005 sec)
#' # ===========================================================================
#' # Compare different length limits
#'
#' # Brief summary for quick reference
#' summary_brief <- gpdb_summarize("TP53", max_length = 200)
#' cat("Brief (200 chars):\n", summary_brief, "\n\n")
#'
#' # Detailed summary for comprehensive overview
#' summary_detailed <- gpdb_summarize("TP53", max_length = 800)
#' cat("Detailed (800 chars):\n", summary_detailed, "\n\n")
#'
#' # Very concise for table display
#' summary_concise <- gpdb_summarize("TP53", max_length = 150)
#' cat("Concise (150 chars):\n", summary_concise, "\n")
#'
#' # ===========================================================================
#' # Example 3: Batch Summarization for Multiple Genes
#' # ===========================================================================
#' # Generate summaries for gene panel
#'
#' genes <- c("TP53", "METTL3", "MYC", "BRCA1")
#'
#' # Create brief summaries for all genes
#' summaries <- lapply(genes, function(g) {
#'   gpdb_summarize(g, max_length = 250)
#' })
#' names(summaries) <- genes
#'
#' # Display all summaries
#' for (gene in genes) {
#'   cat("\n", gene, ":\n", summaries[[gene]], "\n", sep = "")
#' }
#'
#' # ===========================================================================
#' # Example 4: Include Parameter Control (future enhancement)
#' # ===========================================================================
#' # Currently only "function" and "targets" are implemented
#'
#' summary_all <- gpdb_summarize(
#'   gene = "TP53",
#'   max_length = 600,
#'   include = c("function", "targets") # Current default
#' )
#'
#' # Note: "pathways" and "diseases" will be added in future versions
#'
#' # ===========================================================================
#' # Next Steps (Brief pointers to downstream analysis)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - Get detailed targets: ?gpdb_find_targets
#' # - Load raw data: ?gpdb_load_data
#' # - Enrichment analysis: ?gpdb_enrich
#' # - Network visualization: ?gpdb_plot_network
#' # - Compare genes: ?gpdb_compare_genes
#' }
gpdb_summarize <- function(gene,
                           max_length = 500,
                           include = c("function", "targets", "pathways")) {
  gene <- .gpdb_format_genes(gene)[1]

  # Count datasets
  dataset_query <- paste0(
    "SELECT COUNT(*) as n, ",
    "COUNT(DISTINCT CellLineName) as n_cells, ",
    "COUNT(DISTINCT TissueSite) as n_tissues ",
    "FROM datasets WHERE pbgene = '", .gpdb_sql_safe(gene), "'"
  )

  dataset_info <- .gpdb_execute_query(dataset_query)

  if (dataset_info$n == 0) {
    return(paste0("No perturbation data available for ", gene, "."))
  }

  # Get effects
  effects_query <- paste0(
    "SELECT target_gene, logfc_mean, n_datasets, confidence ",
    "FROM gene_effects_agg ",
    "WHERE perturbed_gene = '", .gpdb_sql_safe(gene), "' ",
    "ORDER BY ABS(logfc_mean) DESC LIMIT 10"
  )

  top_effects <- .gpdb_execute_query(effects_query)

  # Build summary
  summary_parts <- list()

  # Introduction
  summary_parts$intro <- sprintf(
    "%s has been studied in %d perturbation experiments across %d cell lines and %d tissue types.",
    gene,
    dataset_info$n,
    dataset_info$n_cells,
    dataset_info$n_tissues
  )

  # Effects
  if ("targets" %in% include && nrow(top_effects) > 0) {
    up_targets <- top_effects[top_effects$logfc_mean > 0, ]
    down_targets <- top_effects[top_effects$logfc_mean < 0, ]

    if (nrow(up_targets) > 0) {
      top_up <- head(up_targets, 3)
      summary_parts$up <- sprintf(
        "Perturbation leads to upregulation of %s (%s-fold)",
        paste(top_up$target_gene, collapse = ", "),
        paste(round(top_up$logfc_mean, 2), collapse = ", ")
      )
    }

    if (nrow(down_targets) > 0) {
      top_down <- head(down_targets, 3)
      summary_parts$down <- sprintf(
        "and downregulation of %s (%s-fold).",
        paste(top_down$target_gene, collapse = ", "),
        paste(round(abs(top_down$logfc_mean), 2), collapse = ", ")
      )
    }
  }

  # Combine
  full_summary <- paste(unlist(summary_parts), collapse = " ")

  # Truncate if needed
  if (nchar(full_summary) > max_length) {
    full_summary <- paste0(substr(full_summary, 1, max_length - 3), "...")
  }

  return(full_summary)
}


#' Predict Therapeutic Targets from Disease Gene Expression Signature
#'
#' Identify candidate drug targets by Gene Perturbation Similarity Analysis (GPSA),
#' matching disease signatures against 18.9M pre-computed gene-gene regulatory relationships
#' from 7,665 perturbation experiments. Finds genes whose knockdown/overexpression would
#' reverse (therapeutic mode) or mimic (mechanistic mode) the disease transcriptional profile.
#' **Database covers 2,810 perturbed genes** with evidence-weighted scoring based on
#' multi-dataset consensus. Optimized with vectorized SQL queries for fast signature matching.
#'
#' @param disease_signature Data.frame with disease gene expression profile. Required columns:
#'   \itemize{
#'     \item \strong{gene}: Character vector of gene symbols (e.g., "MYC", "TP53", "KRAS").
#'       Will be auto-formatted to uppercase. Use 5-50 signature genes for optimal results.
#'     \item \strong{logFC} OR \strong{direction}: Either numeric log2 fold changes (recommended)
#'       or categorical directions ("up"/"down", "UP"/"DOWN", "upregulated"/"downregulated").
#'       logFC values are used for magnitude weighting; directions get equal weight.
#'   }
#'   Example: \code{data.frame(gene = c("MYC", "TP53"), logFC = c(2.5, -1.8))}
#' @param mode Character. Prediction strategy for target identification.
#'   Default: \code{"reverse"}. Options:
#'   \itemize{
#'     \item \strong{"reverse"}: Find genes producing OPPOSITE effects (therapeutic targets).
#'       Use when: Disease has upregulated oncogenes → find genes that downregulate them.
#'       Example: High MYC in cancer → identify genes whose knockdown reduces MYC.
#'     \item \strong{"mimic"}: Find genes producing SIMILAR effects (mechanistic insights).
#'       Use when: Studying disease mechanisms or identifying phenocopy genes.
#'       Example: Find other genes that produce cancer-like expression changes.
#'   }
#' @param top_n Integer. Number of top-ranked candidates to return.
#'   Default: \code{10}. Range: 5-100. Higher values include lower-scoring candidates
#'   that may have partial signature matches or lower confidence evidence.
#' @param min_confidence Character. Evidence strength filter for regulatory relationships.
#'   Default: \code{"medium"}. Options:
#'   \itemize{
#'     \item "high": 5+ datasets, >80\% consistency (publication-grade, fewer candidates)
#'     \item "medium": 2-4 datasets, >70\% consistency (balanced, recommended)
#'     \item "low": 1 dataset (exploratory, includes preliminary findings)
#'   }
#'
#' @return Data.frame with ranked candidate targets and scoring metrics. Columns:
#' \describe{
#'   \item{\strong{perturbed_gene}}{Candidate target gene symbol (gene to perturb)}
#'   \item{\strong{total_score}}{Weighted composite score (higher = better match).
#'     Computed as: sum(match_score × effect_size × consistency × signature_magnitude)
#'     across all signature genes. Top candidates typically have scores >10.}
#'   \item{\strong{n_signature_matches}}{Number of signature genes affected by this candidate.
#'     Higher values indicate broader coverage of disease signature.}
#'   \item{\strong{n_positive_matches}}{Number of genes affected in desired direction.
#'     Should be close to n_signature_matches for high-quality candidates.}
#'   \item{\strong{match_rate}}{Proportion of matches in correct direction (0-1).
#'     Values >0.8 indicate consistent directional effects across signature genes.}
#'   \item{\strong{avg_effect_size}}{Mean |logFC| across regulated signature genes.
#'     Higher values (>2) indicate strong perturbation effects.}
#'   \item{\strong{avg_consistency}}{Mean consistency score across datasets (0-1).
#'     Values >0.8 indicate reproducible effects across experimental contexts.}
#'   \item{\strong{avg_n_datasets}}{Mean number of supporting datasets per relationship.
#'     Higher values provide stronger evidence base for predictions.}
#' }
#'
#' **Interpreting Results**:
#' \itemize{
#'   \item **Top candidates** (rank 1-3): Highest priority for experimental validation
#'   \item **Match rate >0.8**: Consistently regulates signature in desired direction
#'   \item **High total_score + high match_rate**: Best therapeutic target candidates
#'   \item **High n_signature_matches**: Broad impact on disease signature
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Examine regulatory details: \code{\link{gpdb_find_targets}(candidate_gene, top_n = 100)}
#'   \item Check experimental evidence: \code{\link{gpdb_list_datasets}(candidate_gene)}
#'   \item Pathway enrichment: \code{\link{gpdb_enrich}(signature_genes, enrich.type = "KEGG")}
#'   \item Network visualization: \code{\link{gpdb_plot_network}(candidate_gene, top_targets = 50)}
#'   \item Validate interactions: \code{\link{gpdb_predict_interaction}(candidate, signature_gene)}
#'   \item Load validation data: \code{\link{gpdb_load_data}(dataset_id)}
#' }
#'
#' **Alternative Methods**:
#' \itemize{
#'   \item \code{\link{gpdb_find_regulators}}: Find upstream regulators of single gene.
#'     Use when: Have one key disease gene, not full signature.
#'   \item \code{\link{gpdb_compare_genes}}: Compare effects of multiple candidate targets.
#'     Use when: Evaluating alternative therapeutic targets side-by-side.
#'   \item \code{\link{gpdb_analyze_cascade}}: Trace multi-step regulatory cascades.
#'     Use when: Understanding indirect regulation mechanisms.
#' }
#'
#' @details
#' **Scoring Algorithm**:
#'
#' For each candidate gene, compute weighted match score across signature:
#' \enumerate{
#'   \item **Direction matching**: For each signature gene regulated by candidate:
#'     \itemize{
#'       \item Reverse mode: +1 if opposite direction, -1 if same direction
#'       \item Mimic mode: +1 if same direction, -1 if opposite direction
#'     }
#'   \item **Weight calculation**: match_score × |effect_size| × consistency × signature_magnitude
#'   \item **Aggregation**: Sum weighted scores across all signature genes
#'   \item **Ranking**: Sort by total_score (descending)
#' }
#'
#' **Performance Optimization**:
#' \itemize{
#'   \item Vectorized SQL: Single query for all signature genes (not loop)
#'   \item Batch aggregation: dplyr group_by for candidate scoring
#'   \item Indexed lookups: Sub-second queries for signatures <100 genes
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test signature**: 5-gene oncogene panel (MYC, CCND1, CDK4, E2F1, BIRC5)
#'   \item **Runtime**: 1.02 sec (reverse mode, 10 candidates with medium confidence)
#'   \item **Result**: 10 candidates, top score 15.88 (PSPC1, 5/5 signature matches, 100\% match rate)
#'   \item **Database**: 18.9M relationships, vectorized SQL query + dplyr aggregation
#'   \item **Follow-up queries**: 0.07-0.28 sec (cached gene mappings)
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Upstream Preparation** (create disease signature):
#' \itemize{
#'   \item \code{\link{gpdb_load_deg}}: Load DEG results from dataset as signature template
#'   \item External: DESeq2, edgeR, limma for differential expression analysis
#' }
#'
#' **Downstream Analysis** (validate predictions):
#' \itemize{
#'   \item \code{\link{gpdb_find_targets}}: Detailed target list for candidate gene
#'   \item \code{\link{gpdb_list_datasets}}: Experimental evidence for candidate
#'   \item \code{\link{gpdb_load_data}}: Raw expression data for validation
#'   \item \code{\link{gpdb_enrich}}: Pathway enrichment on signature or targets
#'   \item \code{\link{gpdb_plot_network}}: Visualize candidate-signature network
#'   \item \code{\link{gpdb_predict_interaction}}: Test gene-gene interaction type
#' }
#'
#' **Related Prediction Functions** (parallel comparison):
#' \itemize{
#'   \item \code{\link{gpdb_find_regulators}}: Find regulators of single target gene.
#'     Use when: Have one disease driver gene instead of full signature.
#'   \item \code{\link{gpdb_what_happens}}: Query comprehensive effects of candidate.
#'     Use when: Evaluating off-target effects of predicted therapeutic target.
#'   \item \code{\link{gpdb_compare_genes}}: Compare multiple candidate targets.
#'     Use when: Selecting between alternative therapeutic strategies.
#'   \item \code{\link{gpdb_analyze_cascade}}: Trace regulatory cascades from candidate.
#'     Use when: Understanding indirect mechanisms of target action.
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I find drug targets from disease gene expression data?
#'   \item What genes should I knock down to reverse disease signature?
#'   \item Can I predict therapeutic targets from RNA-seq data?
#'   \item How many signature genes do I need for reliable predictions?
#'   \item What is the difference between reverse and mimic mode?
#'   \item How do I interpret the total_score and match_rate?
#'   \item Which candidates are best for experimental validation?
#'   \item Can I use directional signatures without fold changes?
#'   \item How reliable are predictions with low signature coverage?
#'   \item What evidence supports the predicted therapeutic targets?
#'   \item How do I validate predicted targets in my system?
#'   \item Can I find targets that affect specific pathway genes?
#'   \item What if my disease signature has mixed up/down genes?
#'   \item How do I filter candidates by effect size or consistency?
#'   \item Can I compare predictions across different signatures?
#'   \item What downstream experiments should I prioritize?
#'   \item How do I check off-target effects of predicted targets?
#'   \item Can I use this for drug repurposing analysis?
#'   \item What cell lines were used for candidate target data?
#'   \item How do I export results for further computational analysis?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Usage - Reverse Mode (TESTED - runtime: 1.02 sec)
#' # ===========================================================================
#' # Scientific Question: Find therapeutic targets for oncogene-driven disease
#' # Strategy: Identify genes whose knockdown reverses oncogene upregulation
#' # Expected output: Ranked candidates with scoring metrics
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' # Create disease signature (oncogene upregulation in cancer)
#' disease_sig <- data.frame(
#'   gene = c("MYC", "CCND1", "CDK4", "E2F1", "BIRC5"),
#'   logFC = c(2.5, 1.8, 2.1, 1.5, 2.0)
#' )
#'
#' # Find genes that reverse this signature (therapeutic targets)
#' targets <- gpdb_predict_targets(
#'   disease_signature = disease_sig,
#'   mode = "reverse",
#'   top_n = 10,
#'   min_confidence = "medium"
#' )
#'
#' # Examine top candidates
#' cat("Found", nrow(targets), "candidate therapeutic targets\n")
#' cat("\nTop candidate:", targets$perturbed_gene[1], "\n")
#' cat("  Total score:", round(targets$total_score[1], 2), "\n")
#' cat(
#'   "  Signature coverage:", targets$n_signature_matches[1], "/",
#'   nrow(disease_sig), "genes\n"
#' )
#' cat("  Match rate:", round(targets$match_rate[1] * 100, 1), "%\n")
#'
#' # View all candidates
#' print(targets[, c(
#'   "perturbed_gene", "total_score", "n_signature_matches",
#'   "match_rate", "avg_effect_size"
#' )])
#'
#' # ===========================================================================
#' # Example 2: Mimic Mode - Mechanistic Study (TESTED - runtime: 0.07 sec)
#' # ===========================================================================
#' # Scientific Question: Find genes that produce similar disease-like effects
#' # Use case: Identify phenocopy genes for mechanism studies
#'
#' # Same disease signature
#' mimics <- gpdb_predict_targets(
#'   disease_signature = disease_sig,
#'   mode = "mimic", # Find genes producing SIMILAR effects
#'   top_n = 5,
#'   min_confidence = "high" # High-confidence only
#' )
#'
#' cat("\nGenes that mimic disease signature:\n")
#' print(mimics[, c("perturbed_gene", "total_score", "match_rate")])
#'
#' # ===========================================================================
#' # Example 3: Direction-Based Signature (TESTED - runtime: 0.28 sec)
#' # ===========================================================================
#' # Use case: Have categorical up/down information without exact fold changes
#'
#' # Tumor suppressor downregulation signature
#' sig_direction <- data.frame(
#'   gene = c("TP53", "CDKN1A", "BAX"),
#'   direction = c("down", "down", "down")
#' )
#'
#' # Find targets that restore tumor suppressor expression
#' restorers <- gpdb_predict_targets(
#'   disease_signature = sig_direction,
#'   mode = "reverse",
#'   top_n = 5
#' )
#'
#' cat("\nCandidates that restore tumor suppressor expression:\n")
#' print(restorers)
#'
#' # ===========================================================================
#' # Example 4: Filter and Prioritize Candidates
#' # ===========================================================================
#' # Select high-quality candidates for experimental validation
#'
#' # Criteria: High match rate + broad signature coverage + strong effects
#' high_quality <- targets[
#'   targets$match_rate > 0.8 & # >80% correct direction
#'     targets$n_signature_matches >= 4 & # Affects 4+ signature genes
#'     targets$avg_effect_size > 1.5, # Strong effect size
#' ]
#'
#' cat("\nHigh-quality candidates for validation:\n")
#' print(high_quality[, c(
#'   "perturbed_gene", "total_score",
#'   "n_signature_matches", "match_rate"
#' )])
#'
#' # ===========================================================================
#' # Example 5: Compare Reverse vs Mimic for Same Signature
#' # ===========================================================================
#' # Understand both therapeutic targets and disease mechanism
#'
#' # Small signature for quick comparison
#' small_sig <- data.frame(
#'   gene = c("MYC", "CCND1", "CDK4"),
#'   logFC = c(2.5, 1.8, 2.1)
#' )
#'
#' # Therapeutic targets (reverse)
#' targets_rev <- gpdb_predict_targets(small_sig, mode = "reverse", top_n = 3)
#'
#' # Phenocopy genes (mimic)
#' targets_mim <- gpdb_predict_targets(small_sig, mode = "mimic", top_n = 3)
#'
#' cat("\nTherapeutic targets (reverse signature):\n")
#' print(targets_rev$perturbed_gene)
#'
#' cat("\nPhenocopy genes (mimic signature):\n")
#' print(targets_mim$perturbed_gene)
#'
#' # ===========================================================================
#' # Example 6: Validate Top Candidate
#' # ===========================================================================
#' # Get detailed evidence for top-ranked therapeutic target
#'
#' top_candidate <- targets$perturbed_gene[1]
#' cat("\nValidating top candidate:", top_candidate, "\n")
#'
#' # Check which signature genes it regulates
#' candidate_targets <- gpdb_find_targets(
#'   gene = top_candidate,
#'   top_n = 100,
#'   min_confidence = "medium"
#' )
#'
#' # Find overlap with disease signature
#' overlap_up <- intersect(
#'   candidate_targets$downregulated$target_gene, # Candidate downregulates
#'   disease_sig$gene[disease_sig$logFC > 0] # Disease upregulates
#' )
#'
#' cat("Signature genes reversed by candidate:", length(overlap_up), "\n")
#' cat("Genes:", paste(overlap_up, collapse = ", "), "\n")
#'
#' # Check experimental evidence
#' datasets <- gpdb_list_datasets(gene = top_candidate)
#' cat("Supporting datasets:", nrow(datasets), "\n")
#'
#' # ===========================================================================
#' # Next Steps (Brief pointers to downstream analysis)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - Detailed regulation: ?gpdb_find_targets
#' # - Experimental evidence: ?gpdb_list_datasets
#' # - Load validation data: ?gpdb_load_data
#' # - Network visualization: ?gpdb_plot_network
#' # - Pathway enrichment: ?gpdb_enrich
#' # - Compare candidates: ?gpdb_compare_genes
#' }
gpdb_predict_targets <- function(disease_signature,
                                 mode = c("reverse", "mimic"),
                                 top_n = 10,
                                 min_confidence = "medium") {
  mode <- match.arg(mode)

  # Validate input
  if (!is.data.frame(disease_signature)) {
    stop("disease_signature must be a data frame", call. = FALSE)
  }

  if (!"gene" %in% names(disease_signature)) {
    stop("disease_signature must have a 'gene' column", call. = FALSE)
  }

  # Standardize columns
  if ("direction" %in% names(disease_signature)) {
    disease_signature$target_direction <- ifelse(
      disease_signature$direction %in% c("up", "UP", "upregulated"), 1, -1
    )
    disease_signature$target_magnitude <- 1 # Default weight
  } else if ("logFC" %in% names(disease_signature)) {
    disease_signature$target_direction <- sign(disease_signature$logFC)
    disease_signature$target_magnitude <- abs(disease_signature$logFC)
  } else {
    stop("disease_signature must have 'direction' or 'logFC' column", call. = FALSE)
  }

  # Remove invalid entries
  disease_signature <- disease_signature[disease_signature$target_direction != 0, ]

  if (nrow(disease_signature) == 0) {
    stop("No valid signature genes", call. = FALSE)
  }

  # Format gene names
  disease_signature$gene <- .gpdb_format_genes(disease_signature$gene)

  # === VECTORIZED: Single SQL query for ALL signature genes ===
  gene_list_sql <- paste0("'", disease_signature$gene, "'", collapse = ", ")

  query <- paste0(
    "SELECT perturbed_gene, target_gene, logfc_mean, n_datasets, ",
    "consistency_score, confidence ",
    "FROM gene_effects_agg ",
    "WHERE target_gene IN (", gene_list_sql, ") "
  )

  # Add confidence filter
  if (min_confidence == "high") {
    query <- paste0(query, "AND confidence = 'high'")
  } else if (min_confidence == "medium") {
    query <- paste0(query, "AND confidence IN ('high', 'medium')")
  }

  # ONE query instead of loop!
  all_effects <- .gpdb_execute_query(query)

  if (nrow(all_effects) == 0) {
    message("No regulators found for signature genes")
    return(data.frame())
  }

  # Prepare signature data
  sig_data <- disease_signature
  sig_data$gene <- NULL # Remove original gene column if exists
  names(sig_data)[names(sig_data) == "target_gene"] <- "target_gene"

  # Add target_gene column if not already renamed
  if (!"target_gene" %in% names(sig_data)) {
    sig_data$target_gene <- disease_signature$gene
  }

  # Keep only needed columns
  sig_data <- sig_data[, c("target_gene", "target_direction", "target_magnitude")]

  # Merge effects with signature
  effects_merged <- merge(
    all_effects,
    sig_data,
    by = "target_gene",
    all.x = FALSE
  )

  # Calculate match score (vectorized)
  if (mode == "reverse") {
    # Opposite direction = positive score
    effects_merged$match_score <- -sign(effects_merged$logfc_mean) * effects_merged$target_direction
  } else {
    # Same direction = positive score
    effects_merged$match_score <- sign(effects_merged$logfc_mean) * effects_merged$target_direction
  }

  # Calculate weighted score
  effects_merged$weighted_score <- effects_merged$match_score *
    abs(effects_merged$logfc_mean) *
    effects_merged$consistency_score *
    effects_merged$target_magnitude

  # Aggregate by candidate gene using dplyr
  candidate_scores <- effects_merged |>
    dplyr::group_by(perturbed_gene) |>
    dplyr::summarise(
      total_score = sum(weighted_score, na.rm = TRUE),
      n_signature_matches = dplyr::n(),
      n_positive_matches = sum(weighted_score > 0),
      match_rate = sum(weighted_score > 0) / dplyr::n(),
      avg_effect_size = mean(abs(logfc_mean), na.rm = TRUE),
      avg_consistency = mean(consistency_score, na.rm = TRUE),
      avg_n_datasets = mean(n_datasets, na.rm = TRUE),
      .groups = "drop"
    ) |>
    dplyr::arrange(dplyr::desc(total_score))

  # Take top N
  result <- head(as.data.frame(candidate_scores), top_n)

  message("Analyzed ", nrow(disease_signature), " signature genes")
  message("Found ", nrow(result), " candidate targets")
  if (nrow(result) > 0) {
    message(
      "Top candidate: ", result$perturbed_gene[1],
      " (score: ", round(result$total_score[1], 2),
      ", matches: ", result$n_signature_matches[1], "/", nrow(disease_signature), ")"
    )
  }

  return(result)
}


#' Predict Gene-Gene Interaction Type from Perturbation Effects
#'
#' Infer synergistic, antagonistic, or independent relationships between two genes
#' by comparing their regulatory effects on common target genes across 18.9M pre-computed
#' relationships from 7,665 perturbation experiments. Computes Pearson correlation
#' between gene1 and gene2 effects on shared targets to classify interaction types.
#' **Requires 5+ common targets** for reliable predictions. Fast correlation-based
#' scoring with sub-100ms runtime for most gene pairs.
#'
#' @param gene1 Character. First gene symbol (e.g., "TP53", "MYC"). Will be auto-formatted
#'   to uppercase. Must exist in database (use \code{gpdb_search()} to verify). This gene's
#'   perturbation effects will be compared with gene2.
#' @param gene2 Character. Second gene symbol (e.g., "MDM2", "KRAS"). Will be auto-formatted
#'   to uppercase. Must exist in database. Order of gene1/gene2 does not affect results
#'   (correlation is symmetric).
#' @param interaction_type Character. Prediction mode for interaction classification.
#'   Default: \code{"infer"}. Options:
#'   \itemize{
#'     \item \strong{"infer"}: Automatically classify based on correlation threshold (recommended).
#'       Uses: r > 0.5 = synergistic, r < -0.5 = antagonistic, else independent.
#'     \item \strong{"synergy"}: Force synergistic classification (returns correlation data only).
#'       Use when: Testing specific hypothesis about cooperative effects.
#'     \item \strong{"antagonism"}: Force antagonistic classification (returns correlation data only).
#'       Use when: Testing specific hypothesis about opposing effects.
#'   }
#'
#' @return List with interaction prediction and supporting evidence. Elements:
#' \describe{
#'   \item{\strong{gene1}}{First gene symbol (formatted)}
#'   \item{\strong{gene2}}{Second gene symbol (formatted)}
#'   \item{\strong{prediction}}{Interaction type classification:
#'     \itemize{
#'       \item "synergistic": Genes produce similar effects on targets (r > 0.5)
#'       \item "antagonistic": Genes produce opposite effects on targets (r < -0.5)
#'       \item "independent": Genes have uncorrelated effects (-0.5 ≤ r ≤ 0.5)
#'       \item "insufficient_data": One/both genes lack perturbation data
#'       \item "insufficient_overlap": <5 common targets found
#'     }}
#'   \item{\strong{correlation}}{Pearson correlation coefficient between gene1 and gene2
#'     effects on common targets. Range: -1 to +1. Values closer to ±1 indicate stronger
#'     interaction (synergistic or antagonistic). Values near 0 suggest independence.}
#'   \item{\strong{n_common_targets}}{Number of genes regulated by BOTH gene1 and gene2.
#'     Larger values (>100) provide more robust correlation estimates. Typical range: 5,000-20,000
#'     for well-studied genes.}
#'   \item{\strong{evidence}}{Natural language summary explaining the prediction with
#'     correlation value and common target count. LLM-friendly format for automated interpretation.}
#'   \item{\strong{top_common_targets}}{Character vector of top 10 shared target genes,
#'     ranked by combined effect magnitude (|logfc1 + logfc2|). These genes show strongest
#'     co-regulation by both gene1 and gene2.}
#'   \item{\strong{common_effects}}{Data.frame with all common target effects. Columns:
#'     \itemize{
#'       \item target_gene: Shared target gene symbol
#'       \item logfc1: Gene1 perturbation effect (log2 fold change)
#'       \item n1: Number of datasets supporting gene1 → target relationship
#'       \item logfc2: Gene2 perturbation effect (log2 fold change)
#'       \item n2: Number of datasets supporting gene2 → target relationship
#'     }
#'     Use for custom correlation analysis or visualization.}
#' }
#'
#' **Interpreting Results**:
#' \itemize{
#'   \item **Synergistic (r > 0.5)**: Genes work together, co-activate or co-repress targets.
#'     Example: Co-factors in same pathway. Consider: Dual targeting for enhanced effects.
#'   \item **Antagonistic (r < -0.5)**: Genes oppose each other, one activates while other represses.
#'     Example: Feedback regulators. Consider: Synthetic lethality or rescue experiments.
#'   \item **Independent (-0.5 ≤ r ≤ 0.5)**: Genes regulate different target sets or contexts.
#'     Example: Parallel pathways. Consider: Orthogonal perturbations for combinatorial studies.
#'   \item **High n_common_targets (>1000)**: Strong statistical power, reliable prediction
#'   \item **Low n_common_targets (<100)**: Interpret with caution, may reflect context-specific effects
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item Examine shared targets: Explore \code{result$common_effects} data frame
#'   \item Visualize correlation: Plot \code{logfc1 vs logfc2} from common_effects
#'   \item Check individual effects: \code{\link{gpdb_find_targets}(gene1)} and \code{\link{gpdb_find_targets}(gene2)}
#'   \item Network visualization: \code{\link{gpdb_plot_network}(gene1, top_targets = 50)}
#'   \item Pathway enrichment: \code{\link{gpdb_enrich}(top_common_targets, enrich.type = "KEGG")}
#'   \item Load validation data: \code{\link{gpdb_list_datasets}(gene1)} → \code{\link{gpdb_load_data}(dataset_id)}
#' }
#'
#' **Alternative Methods**:
#' \itemize{
#'   \item \code{\link{gpdb_compare_genes}}: Compare multiple genes' effects side-by-side.
#'     Use when: Analyzing gene family or pathway members simultaneously.
#'   \item \code{\link{gpdb_analyze_cascade}}: Trace multi-step regulatory paths between genes.
#'     Use when: Investigating indirect regulatory relationships.
#'   \item \code{\link{gpdb_predict_targets}}: Find genes that reverse/mimic expression signature.
#'     Use when: Looking for therapeutic targets for gene pair's combined effects.
#' }
#'
#' @details
#' **Interaction Inference Algorithm**:
#'
#' \enumerate{
#'   \item **Query effects**: Retrieve all targets regulated by gene1 and gene2 separately
#'     (two SQL queries to gene_effects_agg table)
#'   \item **Find overlap**: Merge on target_gene to identify common targets
#'   \item **Filter**: Require ≥5 common targets for reliable correlation estimation
#'   \item **Compute correlation**: Pearson correlation between logfc1 and logfc2 vectors
#'   \item **Classify interaction**:
#'     \itemize{
#'       \item If r > 0.5: synergistic (similar effects on targets)
#'       \item If r < -0.5: antagonistic (opposite effects on targets)
#'       \item If -0.5 ≤ r ≤ 0.5: independent (uncorrelated effects)
#'     }
#'   \item **Rank targets**: Sort common targets by |logfc1 + logfc2| (combined magnitude)
#' }
#'
#' **Correlation Interpretation**:
#' \itemize{
#'   \item **r = +1.0**: Perfect positive correlation (always co-activate/co-repress)
#'   \item **r = +0.7**: Strong synergy (mostly similar effects)
#'   \item **r = +0.3**: Weak positive correlation (some overlap)
#'   \item **r = 0.0**: No correlation (independent regulation)
#'   \item **r = -0.3**: Weak antagonism (some opposing effects)
#'   \item **r = -0.7**: Strong antagonism (mostly opposite effects)
#'   \item **r = -1.0**: Perfect negative correlation (always oppose)
#' }
#'
#' **Statistical Considerations**:
#' \itemize{
#'   \item Correlation stability increases with n_common_targets
#'   \item Recommend n ≥ 100 for publication-quality predictions
#'   \item Large n (>10,000) typical for well-studied oncogenes, transcription factors
#'   \item Context-specificity: Interaction may vary by tissue/cell line (not captured in aggregate)
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test pair 1**: TP53-MDM2 (classic interaction, 11,584 common targets)
#'   \item **Runtime**: 0.054 sec (two SQL queries + correlation computation)
#'   \item **Result**: r = 0.015, "independent" (expected: context-dependent interaction)
#'   \item **Test pair 2**: METTL3-ALKBH5 (m6A writer/eraser, 20,006 common targets)
#'   \item **Runtime**: 0.059 sec
#'   \item **Result**: r = -0.108, "independent" (weak antagonism, below threshold)
#'   \item **Database**: 18.9M relationships, indexed for fast target retrieval
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Upstream Analysis** (prepare gene pairs):
#' \itemize{
#'   \item \code{\link{gpdb_search}}: Find gene symbols and aliases
#'   \item \code{\link{gpdb_find_pathway_regulators}}: Find genes that regulate specific pathway
#' }
#'
#' **Downstream Validation** (explore interactions):
#' \itemize{
#'   \item \code{\link{gpdb_find_targets}}: Detailed target lists for each gene
#'   \item \code{\link{gpdb_compare_genes}}: Side-by-side comparison of gene effects
#'   \item \code{\link{gpdb_plot_network}}: Visualize regulatory networks
#'   \item \code{\link{gpdb_enrich}}: Pathway enrichment on common targets
#'   \item \code{\link{gpdb_list_datasets}}: Experimental evidence for each gene
#'   \item \code{\link{gpdb_load_data}}: Raw data for custom correlation analysis
#' }
#'
#' **Related Prediction Functions** (parallel analysis):
#' \itemize{
#'   \item \code{\link{gpdb_analyze_cascade}}: Multi-step regulation between genes.
#'     Use when: Testing indirect regulatory relationships (gene1 → intermediate → gene2).
#'   \item \code{\link{gpdb_predict_targets}}: Therapeutic targets from disease signature.
#'     Use when: Predicting which gene(s) to perturb for desired outcome.
#'   \item \code{\link{gpdb_compare_genes}}: Compare effects of multiple genes.
#'     Use when: Building interaction matrix for gene set (all pairwise comparisons).
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I predict if two genes work together or oppose each other?
#'   \item What is the relationship between TP53 and MDM2?
#'   \item Do METTL3 and ALKBH5 have antagonistic effects?
#'   \item How many common targets are needed for reliable prediction?
#'   \item What does correlation coefficient tell me about gene interaction?
#'   \item Can I test synergy between pathway genes?
#'   \item How do I interpret independent vs synergistic interactions?
#'   \item What are the top targets co-regulated by two genes?
#'   \item Can I visualize the correlation between gene effects?
#'   \item How do I validate predicted gene-gene interactions?
#'   \item What if prediction shows insufficient data?
#'   \item Can I test interactions for gene family members?
#'   \item How reliable is the correlation-based prediction?
#'   \item What does high positive correlation indicate biologically?
#'   \item What does high negative correlation suggest?
#'   \item Can I get raw effect data for custom analysis?
#'   \item How do I find all genes that synergize with my gene of interest?
#'   \item What experimental designs validate synergistic interactions?
#'   \item Can I test context-specific interactions (tissue/cell line)?
#'   \item How do I build interaction networks for pathway genes?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Usage - Infer Interaction (TESTED - runtime: 0.054 sec)
#' # ===========================================================================
#' # Scientific Question: Do TP53 and MDM2 have synergistic or antagonistic effects?
#' # Known biology: MDM2 is E3 ligase that degrades TP53 (expected: antagonistic)
#' # Expected output: Interaction prediction with correlation and common targets
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' # Predict interaction between TP53 and MDM2
#' interaction <- gpdb_predict_interaction(
#'   gene1 = "TP53",
#'   gene2 = "MDM2",
#'   interaction_type = "infer"
#' )
#'
#' # Examine prediction
#' cat("Interaction type:", interaction$prediction, "\n")
#' cat("Correlation:", round(interaction$correlation, 3), "\n")
#' cat("Common targets:", interaction$n_common_targets, "\n")
#' cat("Evidence:", interaction$evidence, "\n")
#'
#' # View top co-regulated targets
#' cat("\nTop 5 common targets:\n")
#' print(head(interaction$top_common_targets, 5))
#'
#' # ===========================================================================
#' # Example 2: m6A Enzymes - Writer vs Eraser (TESTED - runtime: 0.059 sec)
#' # ===========================================================================
#' # Scientific Question: Do m6A methyltransferase and demethylase oppose each other?
#' # Known biology: METTL3 adds m6A, ALKBH5 removes m6A (expected: antagonistic)
#'
#' interaction_m6a <- gpdb_predict_interaction("METTL3", "ALKBH5")
#'
#' cat("\nMETTL3-ALKBH5 interaction:\n")
#' cat("  Prediction:", interaction_m6a$prediction, "\n")
#' cat("  Correlation:", round(interaction_m6a$correlation, 3), "\n")
#' cat("  Common targets:", interaction_m6a$n_common_targets, "\n")
#'
#' # Note: r = -0.108 (weak antagonism, below -0.5 threshold)
#' # Suggests context-dependent or indirect antagonism
#'
#' # ===========================================================================
#' # Example 3: Oncogenes - MYC and KRAS (TESTED - runtime: 0.052 sec)
#' # ===========================================================================
#' # Scientific Question: Do MYC and KRAS have synergistic oncogenic effects?
#'
#' interaction_onc <- gpdb_predict_interaction("MYC", "KRAS")
#'
#' cat("\nMYC-KRAS interaction:\n")
#' cat("  Prediction:", interaction_onc$prediction, "\n")
#' cat("  Correlation:", round(interaction_onc$correlation, 3), "\n")
#'
#' # ===========================================================================
#' # Example 4: Explore Common Effects Data Frame
#' # ===========================================================================
#' # Access detailed effects on shared targets for custom analysis
#'
#' common_effects <- interaction$common_effects
#'
#' cat("\nCommon effects structure:\n")
#' cat(
#'   "  Dimensions:", nrow(common_effects), "targets ×",
#'   ncol(common_effects), "columns\n"
#' )
#' cat("  Columns:", paste(names(common_effects), collapse = ", "), "\n")
#'
#' # View top co-regulated targets with strongest combined effects
#' top_10 <- head(common_effects, 10)
#' cat("\nTop 10 targets by combined effect magnitude:\n")
#' print(top_10[, c("target_gene", "logfc1", "logfc2")])
#'
#' # ===========================================================================
#' # Example 5: Visualize Correlation (Custom Plot)
#' # ===========================================================================
#' # Create scatter plot of gene1 vs gene2 effects
#'
#' # Extract effects for visualization
#' effects_df <- interaction$common_effects
#'
#' # Correlation plot (base R)
#' plot(effects_df$logfc1, effects_df$logfc2,
#'   xlab = paste(interaction$gene1, "effect (logFC)"),
#'   ylab = paste(interaction$gene2, "effect (logFC)"),
#'   main = paste("Correlation: r =", round(interaction$correlation, 3)),
#'   pch = 16, col = rgb(0, 0, 0, 0.1),
#'   xlim = c(-5, 5), ylim = c(-5, 5)
#' )
#' abline(h = 0, v = 0, lty = 2, col = "gray")
#' abline(a = 0, b = 1, lty = 2, col = "red") # y = x line for reference
#'
#' # ===========================================================================
#' # Example 6: Filter Strong Co-Regulated Targets
#' # ===========================================================================
#' # Find targets with strong effects from both genes
#'
#' # Criteria: |logfc1| > 1.5 AND |logfc2| > 1.5 (strong effects from both)
#' strong_targets <- effects_df[
#'   abs(effects_df$logfc1) > 1.5 & abs(effects_df$logfc2) > 1.5,
#' ]
#'
#' cat("\nTargets with strong effects from both genes:\n")
#' cat("  Count:", nrow(strong_targets), "\n")
#'
#' if (nrow(strong_targets) > 0) {
#'   # Check directionality: same direction (synergy) or opposite (antagonism)?
#'   same_dir <- sign(strong_targets$logfc1) == sign(strong_targets$logfc2)
#'   cat("  Same direction:", sum(same_dir), "targets\n")
#'   cat("  Opposite direction:", sum(!same_dir), "targets\n")
#'
#'   # Show examples
#'   cat("\nTop 5 strongly co-regulated targets:\n")
#'   print(head(strong_targets[, c("target_gene", "logfc1", "logfc2")], 5))
#' }
#'
#' # ===========================================================================
#' # Example 7: Build Interaction Matrix for Gene Family
#' # ===========================================================================
#' # Test all pairwise interactions within gene set
#'
#' # Example: Cell cycle genes
#' genes <- c("CCND1", "CDK4", "CDK6", "E2F1")
#'
#' # Create pairwise interaction matrix
#' n_genes <- length(genes)
#' correlation_matrix <- matrix(NA, n_genes, n_genes,
#'   dimnames = list(genes, genes)
#' )
#'
#' for (i in 1:n_genes) {
#'   for (j in i:n_genes) {
#'     if (i == j) {
#'       correlation_matrix[i, j] <- 1.0 # Self-correlation
#'     } else {
#'       result <- gpdb_predict_interaction(genes[i], genes[j])
#'       correlation_matrix[i, j] <- result$correlation
#'       correlation_matrix[j, i] <- result$correlation # Symmetric
#'     }
#'   }
#' }
#'
#' cat("\nCell cycle gene interaction matrix:\n")
#' print(round(correlation_matrix, 3))
#'
#' # ===========================================================================
#' # Example 8: Validate Synergistic Interaction
#' # ===========================================================================
#' # Deep dive into a predicted synergistic pair
#'
#' # Find a strongly synergistic pair (for demonstration, assume CCND1-CDK4)
#' syn_test <- gpdb_predict_interaction("CCND1", "CDK4")
#'
#' if (syn_test$prediction == "synergistic") {
#'   cat("\nValidating synergistic interaction:\n")
#'
#'   # 1. Check experimental evidence for each gene
#'   datasets_g1 <- gpdb_list_datasets(syn_test$gene1)
#'   datasets_g2 <- gpdb_list_datasets(syn_test$gene2)
#'   cat("  ", syn_test$gene1, "datasets:", nrow(datasets_g1), "\n")
#'   cat("  ", syn_test$gene2, "datasets:", nrow(datasets_g2), "\n")
#'
#'   # 2. Pathway enrichment on common targets
#'   common_targets <- syn_test$top_common_targets
#'   # enrichment <- gpdb_enrich(common_targets, enrich.type = "KEGG")
#'
#'   # 3. Network visualization
#'   # gpdb_plot_network(syn_test$gene1, top_targets = 50)
#' }
#'
#' # ===========================================================================
#' # Next Steps (Brief pointers to downstream analysis)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - Detailed targets: ?gpdb_find_targets
#' # - Compare effects: ?gpdb_compare_genes
#' # - Network visualization: ?gpdb_plot_network
#' # - Pathway enrichment: ?gpdb_enrich
#' # - Load validation data: ?gpdb_load_data
#' # - Multi-step regulation: ?gpdb_analyze_cascade
#' }
gpdb_predict_interaction <- function(gene1,
                                     gene2,
                                     interaction_type = "infer") {
  gene1 <- .gpdb_format_genes(gene1)[1]
  gene2 <- .gpdb_format_genes(gene2)[1]

  # Get effects of gene1
  query1 <- paste0(
    "SELECT target_gene, logfc_mean as logfc1, n_datasets as n1 ",
    "FROM gene_effects_agg ",
    "WHERE perturbed_gene = '", .gpdb_sql_safe(gene1), "'"
  )

  effects1 <- .gpdb_execute_query(query1)

  # Get effects of gene2
  query2 <- paste0(
    "SELECT target_gene, logfc_mean as logfc2, n_datasets as n2 ",
    "FROM gene_effects_agg ",
    "WHERE perturbed_gene = '", .gpdb_sql_safe(gene2), "'"
  )

  effects2 <- .gpdb_execute_query(query2)

  if (nrow(effects1) == 0 || nrow(effects2) == 0) {
    return(list(
      prediction = "insufficient_data",
      evidence = "Not enough data for one or both genes"
    ))
  }

  # Merge effects on common targets
  common <- merge(effects1, effects2, by = "target_gene")

  if (nrow(common) < 5) {
    return(list(
      prediction = "insufficient_overlap",
      evidence = paste("Only", nrow(common), "common targets found")
    ))
  }

  # Calculate correlation
  correlation <- stats::cor(common$logfc1, common$logfc2, use = "complete.obs")

  # Predict interaction type
  if (interaction_type == "infer") {
    if (correlation > 0.5) {
      prediction <- "synergistic"
      evidence <- sprintf(
        "High positive correlation (r=%.2f) across %d common targets suggests synergistic effects",
        correlation, nrow(common)
      )
    } else if (correlation < -0.5) {
      prediction <- "antagonistic"
      evidence <- sprintf(
        "High negative correlation (r=%.2f) across %d common targets suggests antagonistic effects",
        correlation, nrow(common)
      )
    } else {
      prediction <- "independent"
      evidence <- sprintf(
        "Low correlation (r=%.2f) suggests independent effects",
        correlation
      )
    }
  } else {
    prediction <- interaction_type
    evidence <- sprintf(
      "Correlation: %.2f across %d common targets",
      correlation, nrow(common)
    )
  }

  # Top common targets
  common <- common[order(-abs(common$logfc1 + common$logfc2)), ]
  top_targets <- head(common$target_gene, 10)

  result <- list(
    gene1 = gene1,
    gene2 = gene2,
    prediction = prediction,
    correlation = correlation,
    n_common_targets = nrow(common),
    evidence = evidence,
    top_common_targets = top_targets,
    common_effects = as.data.frame(common)
  )

  message(gene1, " and ", gene2, ": ", prediction)
  message("Correlation: ", round(correlation, 3))
  message("Common targets: ", nrow(common))

  return(result)
}


#' Search Gene Perturbation Database Across Multiple Categories
#'
#' Flexible keyword search across genes, cell lines, and tissue types in the
#' GenePerturbR database (7,665 perturbation experiments). Supports fuzzy matching
#' (default) for partial string matches or exact matching for precise lookups.
#' **Database contains 2,810 genes, 1,063 cell lines, 71 tissue types** with
#' results ranked by experimental coverage (dataset count). Fast sub-20ms queries
#' for single-category searches, ideal for data exploration and validation.
#'
#' @param query Character. Search keyword or pattern. Case-insensitive partial matching
#'   by default (fuzzy mode). Examples: "TP" finds TP53, TP63, PTPN2; "HeLa" finds
#'   HeLa, H1HeLa, HeLa S3; "Liver" finds Liver tissue. Use exact gene symbols
#'   with \code{fuzzy = FALSE} for precise lookups. Whitespace is trimmed automatically.
#' @param search_in Character vector. Categories to search in.
#'   Default: \code{c("genes", "cell_lines", "tissues")} (search all categories).
#'   Options:
#'   \itemize{
#'     \item \strong{"genes"}: Perturbed genes in database (2,810 unique genes).
#'       Searches pbgene column in datasets table. Returns up to 20 top matches
#'       by dataset count.
#'     \item \strong{"cell_lines"}: Cell line names (1,063 unique cell lines).
#'       Searches CellLineName column. Useful for finding available cell systems.
#'     \item \strong{"tissues"}: Tissue types (71 unique tissues).
#'       Searches TissueSite column. Useful for context-specific analysis planning.
#'   }
#'   Can specify single category (e.g., \code{"genes"}) or combination
#'   (e.g., \code{c("genes", "cell_lines")}).
#' @param fuzzy Logical. Enable fuzzy (partial) matching.
#'   Default: \code{TRUE} (recommended for exploration).
#'   \itemize{
#'     \item \strong{TRUE}: SQL LIKE '\%query\%' pattern matching. Finds all entries
#'       containing query string. Case-insensitive. Returns up to 20 results per
#'       category, ranked by dataset count. Use when: Exploring similar gene names,
#'       finding cell line variants, or discovering related entries.
#'     \item \strong{FALSE}: Exact match (SQL = 'query'). Returns only perfect matches.
#'       Use when: Validating specific gene exists, confirming exact cell line name,
#'       or retrieving precise dataset counts.
#'   }
#'
#' @return List with search results from requested categories. Elements:
#' \describe{
#'   \item{\strong{genes}}{Data.frame of matching genes (if "genes" in search_in). Columns:
#'     \itemize{
#'       \item \code{gene}: Gene symbol (perturbed gene in experiments)
#'       \item \code{n_datasets}: Number of perturbation datasets for this gene.
#'         Higher values indicate well-studied genes with more robust evidence.
#'         Range: 1-71 (e.g., TP53 = 71 datasets, rare genes = 1-2 datasets).
#'     }
#'     Sorted by n_datasets (descending). Up to 20 results for fuzzy searches.
#'     NULL or empty data.frame if no matches found.}
#'   \item{\strong{cell_lines}}{Data.frame of matching cell lines (if "cell_lines" in search_in).
#'     Columns: \code{cell_line}, \code{n_datasets}. Sorted by dataset count.
#'     Example: HeLa (239 datasets), A549 (180 datasets). NULL if no matches.}
#'   \item{\strong{tissues}}{Data.frame of matching tissues (if "tissues" in search_in).
#'     Columns: \code{tissue}, \code{n_datasets}. Sorted by dataset count.
#'     Example: Liver (706 datasets), Lung (450 datasets). NULL if no matches.}
#' }
#'
#' **Interpreting Results**:
#' \itemize{
#'   \item **High n_datasets (>20)**: Well-characterized genes/systems, robust for analysis
#'   \item **Medium n_datasets (5-20)**: Moderate evidence, suitable for most studies
#'   \item **Low n_datasets (1-4)**: Limited data, use with caution or for exploratory work
#'   \item **Multiple fuzzy matches**: Review all candidates to select correct gene/cell line
#'   \item **No matches**: Check spelling, try broader query, or gene may not be in database
#' }
#'
#' **What You Can Do Next**:
#' \enumerate{
#'   \item List datasets for found gene: \code{\link{gpdb_list_datasets}(gene)}
#'   \item Query gene effects: \code{\link{gpdb_find_targets}(gene, top_n = 100)}
#'   \item Check what happens: \code{\link{gpdb_what_happens}(gene)}
#'   \item Filter by context: Use found cell_line/tissue in context parameter
#'   \item Compare similar genes: \code{\link{gpdb_compare_genes}(c("TP53", "TP63"))}
#'   \item Get gene info: \code{\link{gpdb_get_info}(dataset_id)} for specific datasets
#' }
#'
#' **Alternative Methods**:
#' \itemize{
#'   \item \code{\link{gpdb_list_datasets}}: List all datasets for known gene.
#'     Use when: Already know exact gene symbol, want to see experimental details.
#'   \item \code{\link{gpdb_find_pathway_regulators}}: Find genes regulating specific pathway.
#'     Use when: Pathway-centric analysis rather than keyword search.
#' }
#'
#' @details
#' **Search Implementation**:
#'
#' For each category in \code{search_in}:
#' \enumerate{
#'   \item **Fuzzy mode** (default): Constructs SQL query with LIKE '\%query\%' pattern
#'     \itemize{
#'       \item Searches: pbgene (genes), CellLineName (cell lines), TissueSite (tissues)
#'       \item Groups by matched field, counts datasets per entry
#'       \item Orders by dataset count (descending)
#'       \item Limits to top 20 results per category
#'     }
#'   \item **Exact mode** (fuzzy=FALSE): Uses SQL = 'query' for precise matching
#'     \itemize{
#'       \item Returns only perfect match (typically 0 or 1 result)
#'       \item Useful for validation and dataset counting
#'     }
#'   \item **SQL safety**: All queries sanitized with \code{.gpdb_sql_safe()} to prevent injection
#' }
#'
#' **Performance Notes**:
#' \itemize{
#'   \item Single-category searches: <2ms (e.g., genes only)
#'   \item Multi-category searches: <20ms (all three categories)
#'   \item Fuzzy searches slightly slower than exact (LIKE vs = operator)
#'   \item Results cached in session (repeat queries are instant)
#' }
#'
#' **Search Tips**:
#' \itemize{
#'   \item **Gene aliases**: Try multiple names (e.g., "CD274" and "PDL1")
#'   \item **Cell line variants**: Fuzzy search finds all (e.g., "HeLa" → HeLa, H1HeLa, HeLa S3)
#'   \item **Tissue names**: Use common names (e.g., "Liver" not "Hepatic")
#'   \item **Abbreviations**: Short queries may return many results (e.g., "TP" finds 20 genes)
#'   \item **Case insensitive**: "tp53", "TP53", "Tp53" all equivalent
#' }
#'
#' @section Performance Test:
#' \itemize{
#'   \item **Test 1**: Fuzzy search "TP" across all categories (genes + cell lines + tissues)
#'   \item **Runtime**: 0.018 sec (three SQL queries with LIKE pattern matching)
#'   \item **Results**: 20 genes (TP53, TP63, PTPN2, ...), 1 cell line, 0 tissues
#'   \item **Test 2**: Gene-only search "MYC"
#'   \item **Runtime**: 0.001 sec (single category query)
#'   \item **Results**: 3 genes (MYC: 22 datasets, MYCN: 6, MYCL: 1)
#'   \item **Test 3**: Exact match "TP53" (fuzzy=FALSE)
#'   \item **Runtime**: 0.0006 sec (exact SQL match, fastest)
#'   \item **Result**: 1 gene (TP53: 71 datasets)
#'   \item **Database**: 7,665 datasets, indexed for fast lookups
#' }
#'
#' @references
#' Guo S, Xu Z, Dong X, Hu D, Jiang Y, Wang Q, Zhang J, Zhou Q, Liu S, Song W (2022).
#' GPSAdb: a comprehensive web resource for interactive exploration of genetic perturbation
#' RNA-seq datasets. Nucleic Acids Research, gkac1066. \doi{10.1093/nar/gkac1066}
#'
#' @seealso
#' **Upstream Discovery** (find genes/cell lines):
#' \itemize{
#'   \item \code{\link{gpdb_find_pathway_regulators}}: Find genes regulating specific pathway
#' }
#'
#' **Downstream Analysis** (use search results):
#' \itemize{
#'   \item \code{\link{gpdb_list_datasets}}: List all experiments for found gene
#'   \item \code{\link{gpdb_find_targets}}: Query regulatory targets of found gene
#'   \item \code{\link{gpdb_what_happens}}: Comprehensive effects of gene perturbation
#'   \item \code{\link{gpdb_get_info}}: Detailed metadata for specific dataset
#'   \item \code{\link{gpdb_load_data}}: Load raw expression data for validation
#'   \item \code{\link{gpdb_compare_genes}}: Compare multiple genes from search results
#' }
#'
#' **Related Utility Functions** (data exploration):
#' \itemize{
#'   \item \code{\link{gpdb_summarize}}: Generate natural language summary for found gene.
#'     Use when: Quick overview of gene function from search result.
#'   \item Database statistics functions (if available) for overview of available data.
#' }
#'
#' @section User Queries:
#' \itemize{
#'   \item How do I check if my gene is in the database?
#'   \item What genes start with "TP" or contain "kinase"?
#'   \item How do I find available cell lines for my tissue type?
#'   \item What is the exact gene symbol used in database?
#'   \item How many datasets exist for my gene of interest?
#'   \item Can I search for gene aliases or alternative names?
#'   \item How do I find all HeLa cell line variants?
#'   \item What tissue types are available in the database?
#'   \item How reliable is a gene with only 1-2 datasets?
#'   \item Why does my search return no results?
#'   \item Can I search for multiple genes at once?
#'   \item How do I find genes in specific cell line?
#'   \item What's the difference between fuzzy and exact search?
#'   \item Which genes have the most experimental coverage?
#'   \item How do I validate if gene name is correct?
#'   \item Can I search for tissue-specific perturbations?
#'   \item How do I explore available experimental systems?
#'   \item What cell lines have most datasets for my gene?
#'   \item How do I find similar gene names?
#'   \item Can I search by gene family or pathway?
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # ===========================================================================
#' # Example 1: Basic Fuzzy Search - All Categories (TESTED - runtime: 0.018 sec)
#' # ===========================================================================
#' # Scientific Question: What entries contain "TP" in database?
#' # Use case: Exploratory search to find TP53 and related genes
#' # Expected output: Genes, cell lines, and tissues matching "TP"
#'
#' # Zero configuration - database auto-loaded!
#' library(GenePerturbR)
#'
#' # Search across all categories
#' results <- gpdb_search("TP")
#'
#' # Examine genes found
#' cat("Genes found:", nrow(results$genes), "\n")
#' print(head(results$genes, 5))
#' # Output: TP53 (71 datasets), TP63 (17), PTPN2 (6), PTPN23 (6), ITPR3 (3)
#'
#' # Check cell lines
#' if (!is.null(results$cell_lines) && nrow(results$cell_lines) > 0) {
#'   cat("\nCell lines found:", nrow(results$cell_lines), "\n")
#'   print(results$cell_lines)
#' }
#'
#' # Check tissues
#' if (!is.null(results$tissues) && nrow(results$tissues) > 0) {
#'   cat("\nTissues found:", nrow(results$tissues), "\n")
#'   print(results$tissues)
#' }
#'
#' # ===========================================================================
#' # Example 2: Gene-Only Search (TESTED - runtime: 0.001 sec)
#' # ===========================================================================
#' # Use case: Find all MYC family genes
#'
#' myc_genes <- gpdb_search("MYC", search_in = "genes")
#'
#' cat("\nMYC family genes:\n")
#' print(myc_genes$genes)
#' # Output: MYC (22 datasets), MYCN (6), MYCL (1)
#'
#' # Sort by dataset count to identify best-studied member
#' cat(
#'   "\nBest-studied MYC gene:", myc_genes$genes$gene[1],
#'   "with", myc_genes$genes$n_datasets[1], "datasets\n"
#' )
#'
#' # ===========================================================================
#' # Example 3: Cell Line Search (TESTED - runtime: 0.001 sec)
#' # ===========================================================================
#' # Use case: Find all HeLa cell line variants
#'
#' hela_variants <- gpdb_search("HeLa", search_in = "cell_lines")
#'
#' cat("\nHeLa cell line variants:\n")
#' print(hela_variants$cell_lines)
#' # Output: HeLa (239 datasets), H1HeLa (2), HeLa S3 (2), HeLaK (1)
#'
#' # Most commonly used variant
#' cat("\nMost common HeLa variant:", hela_variants$cell_lines$cell_line[1], "\n")
#'
#' # ===========================================================================
#' # Example 4: Tissue Search (TESTED - runtime: 0.001 sec)
#' # ===========================================================================
#' # Use case: Check liver tissue data availability
#'
#' liver_data <- gpdb_search("Liver", search_in = "tissues")
#'
#' if (!is.null(liver_data$tissues) && nrow(liver_data$tissues) > 0) {
#'   cat("\nLiver tissue data:\n")
#'   print(liver_data$tissues)
#'   cat("\nLiver datasets available:", liver_data$tissues$n_datasets[1], "\n")
#' } else {
#'   cat("\nNo liver tissue data found\n")
#' }
#'
#' # ===========================================================================
#' # Example 5: Exact Match Validation (TESTED - runtime: 0.0006 sec)
#' # ===========================================================================
#' # Use case: Validate exact gene symbol and get dataset count
#'
#' # Check if TP53 exists (exact match)
#' tp53_exact <- gpdb_search("TP53", search_in = "genes", fuzzy = FALSE)
#'
#' if (!is.null(tp53_exact$genes) && nrow(tp53_exact$genes) > 0) {
#'   cat("\nTP53 validation:\n")
#'   cat("  Gene:", tp53_exact$genes$gene[1], "\n")
#'   cat("  Datasets:", tp53_exact$genes$n_datasets[1], "\n")
#'   cat("  Status: EXISTS in database ✓\n")
#' } else {
#'   cat("\nTP53 not found in database\n")
#' }
#'
#' # ===========================================================================
#' # Example 6: Compare Fuzzy vs Exact Search
#' # ===========================================================================
#' # Understand difference between search modes
#'
#' # Fuzzy: Find all genes containing "CDK"
#' cdk_fuzzy <- gpdb_search("CDK", search_in = "genes", fuzzy = TRUE)
#' cat("\nFuzzy search for 'CDK':\n")
#' cat("  Genes found:", nrow(cdk_fuzzy$genes), "\n")
#' print(cdk_fuzzy$genes)
#' # Expected: CDK4, CDK6, CDKN1A, CDKN2A, etc. (all containing "CDK")
#'
#' # Exact: Find only "CDK4"
#' cdk4_exact <- gpdb_search("CDK4", search_in = "genes", fuzzy = FALSE)
#' cat("\nExact search for 'CDK4':\n")
#' print(cdk4_exact$genes)
#' # Expected: Only CDK4
#'
#' # ===========================================================================
#' # Example 7: Workflow - Search → Validate → Analyze
#' # ===========================================================================
#' # Complete workflow from discovery to analysis
#'
#' # Step 1: Search for gene family
#' hdac_genes <- gpdb_search("HDAC", search_in = "genes")
#' cat("\nStep 1 - Found HDAC genes:\n")
#' print(hdac_genes$genes)
#'
#' # Step 2: Select gene with most data
#' top_hdac <- hdac_genes$genes$gene[1]
#' cat(
#'   "\nStep 2 - Selected:", top_hdac,
#'   "with", hdac_genes$genes$n_datasets[1], "datasets\n"
#' )
#'
#' # Step 3: List available datasets
#' datasets <- gpdb_list_datasets(top_hdac)
#' cat("\nStep 3 - Available datasets:\n")
#' print(head(datasets, 3))
#'
#' # Step 4: Query gene effects
#' targets <- gpdb_find_targets(top_hdac, top_n = 20)
#' cat("\nStep 4 - Top regulated targets:\n")
#' print(head(targets$upregulated, 3))
#'
#' # ===========================================================================
#' # Example 8: Multi-Category Search for Context Planning
#' # ===========================================================================
#' # Plan context-specific analysis
#'
#' # Search for brain-related entries
#' brain_search <- gpdb_search("brain")
#'
#' cat("\nBrain-related database content:\n")
#' if (!is.null(brain_search$genes) && nrow(brain_search$genes) > 0) {
#'   cat("  Brain-related genes:", nrow(brain_search$genes), "\n")
#' }
#' if (!is.null(brain_search$cell_lines) && nrow(brain_search$cell_lines) > 0) {
#'   cat("  Brain cell lines:", nrow(brain_search$cell_lines), "\n")
#'   print(brain_search$cell_lines)
#' }
#' if (!is.null(brain_search$tissues) && nrow(brain_search$tissues) > 0) {
#'   cat("  Brain tissues:", nrow(brain_search$tissues), "\n")
#'   print(brain_search$tissues)
#' }
#'
#' # Use results to plan tissue-specific analysis
#' if (!is.null(brain_search$tissues) && nrow(brain_search$tissues) > 0) {
#'   brain_tissue <- brain_search$tissues$tissue[1]
#'   cat("\nUse in context parameter: context = list(tissue = '",
#'     brain_tissue, "')\n",
#'     sep = ""
#'   )
#' }
#'
#' # ===========================================================================
#' # Example 9: Handle No Results
#' # ===========================================================================
#' # Graceful handling when search finds nothing
#'
#' rare_gene <- gpdb_search("NONEXISTENT", search_in = "genes")
#'
#' if (is.null(rare_gene$genes) || nrow(rare_gene$genes) == 0) {
#'   cat("\nNo results found for 'NONEXISTENT'\n")
#'   cat("Suggestions:\n")
#'   cat("  - Check spelling\n")
#'   cat("  - Try alternative gene names or aliases\n")
#'   cat("  - Use broader fuzzy search\n")
#'   cat("  - Gene may not be in database (2,810 genes available)\n")
#' }
#'
#' # ===========================================================================
#' # Next Steps (Brief pointers to downstream analysis)
#' # ===========================================================================
#' # For complete workflows, see:
#' # - List experiments: ?gpdb_list_datasets
#' # - Query effects: ?gpdb_find_targets
#' # - What happens: ?gpdb_what_happens
#' # - Compare genes: ?gpdb_compare_genes
#' # - Get metadata: ?gpdb_get_info
#' # - Load data: ?gpdb_load_data
#' }
gpdb_search <- function(query,
                        search_in = c("genes", "cell_lines", "tissues"),
                        fuzzy = TRUE) {
  query <- trimws(query)

  results <- list()

  # Search in genes
  if ("genes" %in% search_in) {
    if (fuzzy) {
      gene_query <- paste0(
        "SELECT DISTINCT pbgene as gene, COUNT(*) as n_datasets ",
        "FROM datasets ",
        "WHERE pbgene LIKE '%", .gpdb_sql_safe(query), "%' ",
        "GROUP BY pbgene ORDER BY n_datasets DESC LIMIT 20"
      )
    } else {
      gene_query <- paste0(
        "SELECT DISTINCT pbgene as gene, COUNT(*) as n_datasets ",
        "FROM datasets ",
        "WHERE pbgene = '", .gpdb_sql_safe(query), "' ",
        "GROUP BY pbgene"
      )
    }

    gene_results <- .gpdb_execute_query(gene_query)
    results$genes <- as.data.frame(gene_results)
  }

  # Search in cell lines
  if ("cell_lines" %in% search_in) {
    if (fuzzy) {
      cell_query <- paste0(
        "SELECT DISTINCT CellLineName as cell_line, COUNT(*) as n_datasets ",
        "FROM datasets ",
        "WHERE CellLineName LIKE '%", .gpdb_sql_safe(query), "%' ",
        "GROUP BY CellLineName ORDER BY n_datasets DESC LIMIT 20"
      )
    } else {
      cell_query <- paste0(
        "SELECT DISTINCT CellLineName as cell_line, COUNT(*) as n_datasets ",
        "FROM datasets ",
        "WHERE CellLineName = '", .gpdb_sql_safe(query), "' ",
        "GROUP BY CellLineName"
      )
    }

    cell_results <- .gpdb_execute_query(cell_query)
    results$cell_lines <- as.data.frame(cell_results)
  }

  # Search in tissues
  if ("tissues" %in% search_in) {
    if (fuzzy) {
      tissue_query <- paste0(
        "SELECT DISTINCT TissueSite as tissue, COUNT(*) as n_datasets ",
        "FROM datasets ",
        "WHERE TissueSite LIKE '%", .gpdb_sql_safe(query), "%' ",
        "GROUP BY TissueSite ORDER BY n_datasets DESC LIMIT 20"
      )
    } else {
      tissue_query <- paste0(
        "SELECT DISTINCT TissueSite as tissue, COUNT(*) as n_datasets ",
        "FROM datasets ",
        "WHERE TissueSite = '", .gpdb_sql_safe(query), "' ",
        "GROUP BY TissueSite"
      )
    }

    tissue_results <- .gpdb_execute_query(tissue_query)
    results$tissues <- as.data.frame(tissue_results)
  }

  # Print summary
  message("Search results for '", query, "':")
  if (!is.null(results$genes) && nrow(results$genes) > 0) {
    message("  Genes: ", nrow(results$genes), " matches")
  }
  if (!is.null(results$cell_lines) && nrow(results$cell_lines) > 0) {
    message("  Cell lines: ", nrow(results$cell_lines), " matches")
  }
  if (!is.null(results$tissues) && nrow(results$tissues) > 0) {
    message("  Tissues: ", nrow(results$tissues), " matches")
  }

  return(results)
}
