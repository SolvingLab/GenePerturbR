#' Summarize Gene Function
#'
#' Generate comprehensive natural language summary of gene function based on perturbation data
#'
#' @param gene Character. Gene symbol
#' @param max_length Integer. Maximum length of summary in characters (default 500)
#' @param include Character vector. Aspects to include: "function", "targets", "pathways", "diseases"
#'
#' @return Character string with gene summary
#' @export
#'
#' @examples
#' \dontrun{
#' # Get comprehensive summary
#' summary <- gpdb_summarize("METTL3")
#' cat(summary)
#'
#' # Brief summary
#' summary_brief <- gpdb_summarize("TP53", max_length = 200)
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


#' Predict Drug Targets Based on Disease Signature
#'
#' Uses Gene Perturbation Similarity Analysis (GPSA) to identify candidate therapeutic targets.
#' The algorithm finds genes whose knockdown would reverse (or mimic) a disease gene expression signature.
#'
#' @param disease_signature Data frame with two required columns:
#'   \itemize{
#'     \item gene: Gene symbols
#'     \item logFC or direction: Either numeric log fold changes or "up"/"down" labels
#'   }
#' @param mode Character. Analysis mode:
#'   \itemize{
#'     \item "reverse": Find genes that produce OPPOSITE effects (for therapy, default)
#'       Example: If disease has high MYC, find genes whose knockdown reduces MYC
#'     \item "mimic": Find genes that produce SIMILAR effects (for mechanism study)
#'   }
#' @param top_n Integer. Number of top candidate targets to return (default 10)
#' @param min_confidence Character. Minimum evidence quality:
#'   \itemize{
#'     \item "high": Only use relationships from 5+ datasets (most reliable)
#'     \item "medium": Include relationships from 2+ datasets (default)
#'     \item "low": Include all data (exploratory)
#'   }
#'
#' @details
#' Scoring Logic:
#' For each candidate gene, compute a weighted score based on how well it matches the signature:
#' \itemize{
#'   \item Match score: +1 if desired direction, -1 if opposite
#'   \item Weight factors: effect size × consistency × signature magnitude
#'   \item Total score: Sum across all signature genes
#' }
#'
#' Higher scores indicate better candidates that affect more signature genes in the desired direction.
#'
#' @return Data frame with candidate targets ranked by priority
#' @export
#'
#' @examples
#' \dontrun{
#' # Create a disease signature
#' disease_sig <- data.frame(
#'   gene = c("MYC", "TP53", "KRAS"),
#'   logFC = c(2.5, -1.8, 3.2)
#' )
#'
#' # Find genes that reverse this signature
#' targets <- gpdb_predict_targets(disease_sig, mode = "reverse")
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


#' Predict Gene Interaction
#'
#' Predict whether two genes have synergistic or antagonistic effects
#'
#' @param gene1 Character. First gene symbol
#' @param gene2 Character. Second gene symbol
#' @param interaction_type Character. "infer", "synergy", or "antagonism"
#'
#' @return List with interaction prediction and evidence
#' @export
#'
#' @examples
#' \dontrun{
#' # Predict interaction between TP53 and MDM2
#' interaction <- gpdb_predict_interaction("TP53", "MDM2")
#'
#' # Check for synergy between m6A writer and eraser
#' interaction <- gpdb_predict_interaction("METTL3", "ALKBH5")
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


#' Search Database
#'
#' Flexible search across genes, cell lines, and datasets
#'
#' @param query Character. Search query (gene name, cell line, or keyword)
#' @param search_in Character vector. Where to search: "genes", "cell_lines", "tissues"
#' @param fuzzy Logical. Whether to use fuzzy matching (default TRUE)
#'
#' @return List with search results from each category
#' @export
#'
#' @examples
#' \dontrun{
#' # Search for anything related to "liver"
#' results <- gpdb_search("liver")
#'
#' # Search only in genes
#' results <- gpdb_search("TP", search_in = "genes")
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
