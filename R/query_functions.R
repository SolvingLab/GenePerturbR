#' What Happens When a Gene is Knocked Out or Down
#'
#' Query the effects of gene perturbation aggregated across all available datasets.
#' Returns comprehensive information about which genes are affected, their effect sizes,
#' and the reliability of these findings.
#'
#' @param gene Character. Gene symbol to query (e.g., "TP53", "METTL3")
#' @param context List. Optional filters to restrict analysis:
#'   \itemize{
#'     \item cell_line: Specific cell line (e.g., "K-562", "HeLa")
#'     \item tissue: Specific tissue type (e.g., "Liver", "Lung")
#'   }
#' @param aggregate Logical. Whether to return aggregated multi-dataset results (default TRUE)
#' @param top_n Integer. Number of top target genes to return (default 50)
#'
#' @return List containing:
#'   \itemize{
#'     \item summary: Natural language text summary
#'     \item top_upregulated: Genes increased by perturbation
#'     \item top_downregulated: Genes decreased by perturbation
#'     \item all_effects: Complete results
#'     \item stats: Statistics including n_datasets, consistency scores, etc.
#'   }
#'
#' @details
#' Results are aggregated from multiple independent experiments. Key metrics:
#' \itemize{
#'   \item logfc_mean: Average log2 fold change across datasets
#'   \item n_datasets: Number of independent experiments supporting this relationship
#'   \item consistency_score: Fraction of datasets showing same direction (0.5-1.0)
#'   \item confidence: Evidence quality classification (see below)
#' }
#'
#' Confidence levels indicate evidence strength:
#' \itemize{
#'   \item "high": 5+ datasets (reliable, reproducible finding)
#'   \item "medium": 2-4 datasets (moderate evidence)
#'   \item "low": 1 dataset (preliminary, needs validation)
#' }
#'
#' @return List with summary text, top targets, and statistics
#' @export
#'
#' @examples
#' \dontrun{
#' # What happens when TP53 is knocked out?
#' result <- gpdb_what_happens("TP53")
#' cat(result$summary)
#'
#' # In a specific cell line
#' result <- gpdb_what_happens("METTL3", context = list(cell_line = "K-562"))
#' }
gpdb_what_happens <- function(gene,
                              context = NULL,
                              aggregate = TRUE,
                              top_n = 50) {
  gene <- .gpdb_format_genes(gene)[1]

  # Build query using SQL builder
  query <- .gpdb_build_query(
    table = "gene_effects_agg",
    select = "*",
    filters = list(
      perturbed_gene = gene,
      cell_line = context$cell_line,
      tissue = context$tissue
    ),
    order_by = "ABS(logfc_mean) DESC"
  )

  effects <- .gpdb_execute_query(query)

  if (nrow(effects) == 0) {
    message("No data found for gene: ", gene)
    return(list(
      summary = paste("No perturbation data available for", gene),
      data = data.frame(),
      stats = list(n_datasets = 0)
    ))
  }

  # Get dataset count
  dataset_query <- paste0(
    "SELECT COUNT(*) as n FROM datasets WHERE pbgene = '",
    .gpdb_sql_safe(gene), "'"
  )
  n_datasets <- .gpdb_execute_query(dataset_query)$n

  # Separate up and down regulated
  up_targets <- effects[effects$logfc_mean > 0, ]
  down_targets <- effects[effects$logfc_mean < 0, ]

  # Get top targets
  top_up <- head(up_targets[order(-up_targets$logfc_mean), ], top_n)
  top_down <- head(down_targets[order(down_targets$logfc_mean), ], top_n)

  # Generate summary text
  summary_text <- .gpdb_generate_summary(
    gene, n_datasets, top_up, top_down, NULL
  )

  # Statistics
  stats <- list(
    gene = gene,
    n_datasets = n_datasets,
    n_targets_total = nrow(effects),
    n_upregulated = nrow(up_targets),
    n_downregulated = nrow(down_targets),
    n_high_confidence = sum(effects$confidence == "high"),
    avg_effect_size = mean(abs(effects$logfc_mean), na.rm = TRUE)
  )

  result <- list(
    summary = summary_text,
    top_upregulated = as.data.frame(top_up),
    top_downregulated = as.data.frame(top_down),
    all_effects = as.data.frame(effects),
    stats = stats
  )

  message("Found ", n_datasets, " datasets for ", gene)
  message(
    "Total targets: ", nrow(effects),
    " (", nrow(up_targets), " up, ", nrow(down_targets), " down)"
  )

  return(result)
}


#' Find Regulators of Target Gene
#'
#' Identify genes that regulate the expression of a target gene when they are knocked out/down.
#' For example, if gene A knockdown causes gene B to increase, then A is a repressor of B.
#'
#' @param target_gene Character. Target gene symbol (e.g., "MYC", "TP53")
#' @param direction Character. Regulation direction:
#'   \itemize{
#'     \item "up": Find genes whose knockdown INCREASES target (i.e., repressors of target)
#'     \item "down": Find genes whose knockdown DECREASES target (i.e., activators of target)
#'     \item "any": Find both types (default)
#'   }
#' @param top_n Integer. Number of top regulators to return (default 50)
#' @param min_datasets Integer. Minimum number of independent datasets required (default 2).
#'   Higher values give more reliable results but fewer candidates.
#' @param min_confidence Character. Filter by evidence strength:
#'   \itemize{
#'     \item "high": n_datasets >= 5 (strong evidence, highly reliable)
#'     \item "medium": n_datasets >= 2 (moderate evidence, default)
#'     \item "low": n_datasets = 1 (single study, exploratory)
#'     \item "any": Include all levels
#'   }
#' @param return_separate Logical. If TRUE and direction="any", return activators and repressors
#'   as separate lists for clarity (default TRUE)
#'
#' @return Data frame or list with regulatory relationships
#' @export
#'
#' @examples
#' \dontrun{
#' # What regulates MYC? (returns separate lists)
#' regulators <- gpdb_find_regulators("MYC", top_n = 20)
#' regulators$repressors # Genes whose knockdown increases MYC
#' regulators$activators # Genes whose knockdown decreases MYC
#'
#' # High confidence only
#' hc_regs <- gpdb_find_regulators("MYC", min_confidence = "high", top_n = 10)
#' }
gpdb_find_regulators <- function(target_gene,
                                 direction = c("any", "up", "down"),
                                 top_n = 50,
                                 min_datasets = 2,
                                 min_confidence = c("medium", "any", "high", "low"),
                                 return_separate = TRUE) {
  target_gene <- .gpdb_format_genes(target_gene)[1]
  direction <- match.arg(direction)
  min_confidence <- match.arg(min_confidence)

  # Common query parameters
  select_cols <- paste(
    "perturbed_gene, target_gene, logfc_mean, logfc_sd,",
    "n_datasets, consistency_score, effect_size, confidence,",
    "n_tissues, n_celllines"
  )

  base_filters <- list(
    target_gene = target_gene,
    min_datasets = min_datasets,
    min_confidence = if (min_confidence != "any") min_confidence else NULL
  )

  # Direction-specific queries using SQL builder
  if (direction == "up") {
    # Repressors: genes whose knockdown INCREASES target
    query <- .gpdb_build_query(
      select = select_cols,
      filters = c(base_filters, list(direction = "up")),
      order_by = "logfc_mean DESC",
      limit = top_n
    )

    regulators <- .gpdb_execute_query(query)
    regulators$regulation_type <- "repressor"

    message("Found ", nrow(regulators), " repressors of ", target_gene)
    return(as.data.frame(regulators))
  } else if (direction == "down") {
    # Activators: genes whose knockdown DECREASES target
    query <- .gpdb_build_query(
      select = select_cols,
      filters = c(base_filters, list(direction = "down")),
      order_by = "logfc_mean ASC",
      limit = top_n
    )

    regulators <- .gpdb_execute_query(query)
    regulators$regulation_type <- "activator"

    message("Found ", nrow(regulators), " activators of ", target_gene)
    return(as.data.frame(regulators))
  } else {
    # direction = "any"
    if (return_separate) {
      # Query both directions
      query_up <- .gpdb_build_query(
        select = select_cols,
        filters = c(base_filters, list(direction = "up")),
        order_by = "logfc_mean DESC",
        limit = top_n
      )

      query_down <- .gpdb_build_query(
        select = select_cols,
        filters = c(base_filters, list(direction = "down")),
        order_by = "logfc_mean ASC",
        limit = top_n
      )

      repressors <- .gpdb_execute_query(query_up)
      repressors$regulation_type <- "repressor"

      activators <- .gpdb_execute_query(query_down)
      activators$regulation_type <- "activator"

      message(
        "Found ", nrow(repressors), " repressors and ",
        nrow(activators), " activators of ", target_gene
      )

      return(list(
        repressors = as.data.frame(repressors),
        activators = as.data.frame(activators),
        summary = paste0(
          target_gene, " has ", nrow(repressors),
          " repressors and ", nrow(activators), " activators"
        )
      ))
    } else {
      # Return combined
      query <- .gpdb_build_query(
        select = select_cols,
        filters = base_filters,
        order_by = "ABS(logfc_mean) DESC",
        limit = if (!is.null(top_n)) top_n * 2 else NULL
      )

      regulators <- .gpdb_execute_query(query)
      regulators$regulation_type <- ifelse(regulators$logfc_mean > 0,
        "repressor", "activator"
      )

      message("Found ", nrow(regulators), " regulators of ", target_gene)
      return(as.data.frame(regulators))
    }
  }
}


#' Find Targets of Gene
#'
#' Identify which genes are affected when a gene is knocked out/down.
#' Returns genes that are significantly changed (padj < 0.05) across multiple datasets.
#'
#' @param gene Character. Gene symbol to perturb (e.g., "TP53", "METTL3")
#' @param direction Character. Which targets to return:
#'   \itemize{
#'     \item "both": Return both upregulated and downregulated targets (default)
#'     \item "up": Only genes INCREASED by knockdown
#'     \item "down": Only genes DECREASED by knockdown
#'   }
#' @param top_n Integer. Number of top targets to return per direction (default 50)
#' @param min_confidence Character. Minimum evidence strength required:
#'   \itemize{
#'     \item "high": Only relationships supported by 5+ datasets (most reliable)
#'     \item "medium": Include relationships with 2-4 datasets (default)
#'     \item "low": Include single-dataset observations (exploratory)
#'   }
#' @param min_effect_size Numeric. Minimum absolute log2 fold change (default 0.5)
#'
#' @return Data frame or list with target genes and statistics
#' @export
#'
#' @examples
#' \dontrun{
#' # What does METTL3 regulate? (returns separate lists)
#' targets <- gpdb_find_targets("METTL3", top_n = 30)
#' targets$upregulated # Genes upregulated by METTL3 knockdown
#' targets$downregulated # Genes downregulated by METTL3 knockdown
#'
#' # Only strong effects
#' strong_targets <- gpdb_find_targets("TP53", min_effect_size = 1.5)
#' }
gpdb_find_targets <- function(gene,
                              direction = c("both", "up", "down"),
                              top_n = 50,
                              min_confidence = "medium",
                              min_effect_size = 0.5) {
  gene <- .gpdb_format_genes(gene)[1]
  direction <- match.arg(direction)

  # Common query parameters
  select_cols <- paste(
    "target_gene, logfc_mean, logfc_sd, n_datasets,",
    "consistency_score, effect_size, confidence, tissues, celllines"
  )

  base_filters <- list(
    perturbed_gene = gene,
    min_effect_size = min_effect_size,
    min_confidence = if (min_confidence != "low") min_confidence else NULL
  )

  # Direction-specific queries using SQL builder
  if (direction == "up") {
    query <- .gpdb_build_query(
      select = select_cols,
      filters = c(base_filters, list(direction = "up")),
      order_by = "logfc_mean DESC",
      limit = top_n
    )

    targets <- .gpdb_execute_query(query)
    targets$direction <- "up"

    message("Found ", nrow(targets), " upregulated targets of ", gene)
    return(as.data.frame(targets))
  } else if (direction == "down") {
    query <- .gpdb_build_query(
      select = select_cols,
      filters = c(base_filters, list(direction = "down")),
      order_by = "logfc_mean ASC",
      limit = top_n
    )

    targets <- .gpdb_execute_query(query)
    targets$direction <- "down"

    message("Found ", nrow(targets), " downregulated targets of ", gene)
    return(as.data.frame(targets))
  } else {
    # Both directions
    query_up <- .gpdb_build_query(
      select = select_cols,
      filters = c(base_filters, list(direction = "up")),
      order_by = "logfc_mean DESC",
      limit = top_n
    )

    query_down <- .gpdb_build_query(
      select = select_cols,
      filters = c(base_filters, list(direction = "down")),
      order_by = "logfc_mean ASC",
      limit = top_n
    )

    up_targets <- .gpdb_execute_query(query_up)
    up_targets$direction <- "up"

    down_targets <- .gpdb_execute_query(query_down)
    down_targets$direction <- "down"

    message(
      "Found ", nrow(up_targets), " upregulated and ",
      nrow(down_targets), " downregulated targets of ", gene
    )

    return(list(
      upregulated = as.data.frame(up_targets),
      downregulated = as.data.frame(down_targets),
      summary = paste0(
        gene, " regulates ", nrow(up_targets), " upregulated and ",
        nrow(down_targets), " downregulated targets"
      )
    ))
  }
}


#' Compare Multiple Genes
#'
#' Compare the effects of multiple gene perturbations
#'
#' @param genes Character vector. Gene symbols to compare
#' @param metric Character. Comparison metric: "overlap", "correlation", "difference"
#' @param context List. Optional context filters
#'
#' @return List with comparison results and visualization data
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare TP53 and RB1
#' comparison <- gpdb_compare_genes(c("TP53", "RB1"))
#'
#' # Compare m6A writers
#' m6a_comp <- gpdb_compare_genes(c("METTL3", "METTL14"))
#' }
gpdb_compare_genes <- function(genes,
                               metric = c("overlap", "correlation", "difference"),
                               context = NULL) {
  genes <- .gpdb_format_genes(genes)
  metric <- match.arg(metric)

  if (length(genes) < 2) {
    stop("At least 2 genes required for comparison", call. = FALSE)
  }

  # OPTIMIZED: Single SQL query with IN clause (10x faster for multiple genes)
  gene_list_sql <- paste0("'", sapply(genes, .gpdb_sql_safe), "'", collapse = ", ")
  query <- paste0(
    "SELECT perturbed_gene, target_gene, logfc_mean, n_datasets, confidence ",
    "FROM gene_effects_agg ",
    "WHERE perturbed_gene IN (", gene_list_sql, ")"
  )

  all_effects_df <- .gpdb_execute_query(query)

  if (nrow(all_effects_df) == 0) {
    stop("No data found for any of the genes", call. = FALSE)
  }

  # Split by gene
  all_effects <- split(all_effects_df, all_effects_df$perturbed_gene)

  # Check which genes have no data
  missing_genes <- setdiff(genes, names(all_effects))
  if (length(missing_genes) > 0) {
    warning("No data found for: ", paste(missing_genes, collapse = ", "), call. = FALSE)
  }

  if (length(all_effects) < 2) {
    stop("Insufficient data for comparison (need at least 2 genes with data)", call. = FALSE)
  }

  # Calculate overlaps
  all_targets <- lapply(all_effects, function(x) x$target_gene)

  # Common targets
  common_targets <- Reduce(intersect, all_targets)

  # Unique targets for each gene
  unique_targets <- list()
  for (i in seq_along(all_targets)) {
    others <- setdiff(seq_along(all_targets), i)
    unique_targets[[genes[i]]] <- setdiff(
      all_targets[[i]],
      unlist(all_targets[others])
    )
  }

  # Build result
  result <- list(
    genes = genes,
    n_common = length(common_targets),
    common_targets = common_targets,
    unique_targets = unique_targets,
    all_effects = all_effects,
    overlap_matrix = NULL # TODO: implement pairwise overlap matrix
  )

  message("Comparison of ", length(genes), " genes:")
  message("  Common targets: ", length(common_targets))
  for (i in seq_along(unique_targets)) {
    message(
      "  ", names(unique_targets)[i], " unique: ",
      length(unique_targets[[i]])
    )
  }

  return(result)
}


#' Compare Contexts
#'
#' Compare gene perturbation effects across different contexts (tissues/cell lines)
#'
#' @param gene Character. Gene symbol
#' @param contexts List. Named list of context specifications
#' @param metric Character. Comparison metric
#'
#' @return List with context-specific results
#' @export
#'
#' @examples
#' \dontrun{
#' # Compare TP53 effects in different tissues
#' result <- gpdb_compare_contexts(
#'   "TP53",
#'   contexts = list(
#'     liver = list(tissue = "Liver"),
#'     lung = list(tissue = "Lung")
#'   )
#' )
#' }
gpdb_compare_contexts <- function(gene,
                                  contexts,
                                  metric = "specificity") {
  gene <- .gpdb_format_genes(gene)[1]

  # Get effects for each context
  context_effects <- list()

  for (context_name in names(contexts)) {
    ctx <- contexts[[context_name]]
    result <- gpdb_what_happens(gene, context = ctx, aggregate = TRUE)
    context_effects[[context_name]] <- result$all_effects
  }

  # Calculate context-specific and common targets
  all_targets <- lapply(context_effects, function(x) x$target_gene)
  common_targets <- Reduce(intersect, all_targets)

  specific_targets <- list()
  for (ctx_name in names(all_targets)) {
    others_idx <- setdiff(seq_along(all_targets), which(names(all_targets) == ctx_name))
    specific_targets[[ctx_name]] <- setdiff(
      all_targets[[ctx_name]],
      unlist(all_targets[others_idx])
    )
  }

  result <- list(
    gene = gene,
    contexts = names(contexts),
    common_targets = common_targets,
    specific_targets = specific_targets,
    all_effects = context_effects
  )

  message("Context comparison for ", gene, ":")
  message("  Common targets: ", length(common_targets))
  for (ctx in names(specific_targets)) {
    message("  ", ctx, " specific: ", length(specific_targets[[ctx]]))
  }

  return(result)
}


#' Find Genes Affecting Pathway
#'
#' Identify genes whose perturbation affects a specific pathway
#'
#' @param pathway Character. Pathway name (partial match supported)
#' @param effect Character. Effect type: "activate", "suppress", or "any"
#' @param evidence_level Character. Evidence level filter
#'
#' @return Data frame with genes affecting the pathway
#' @export
#'
#' @examples
#' \dontrun{
#' # Which genes affect cell cycle?
#' cc_genes <- gpdb_pathway_genes("cell cycle")
#'
#' # Genes that activate apoptosis
#' apop_genes <- gpdb_pathway_genes("apoptosis", effect = "activate")
#' }
gpdb_pathway_genes <- function(pathway,
                               effect = c("any", "activate", "suppress"),
                               evidence_level = "strong") {
  effect <- match.arg(effect)
  con <- .gpdb_get_connection()

  # TODO: This requires pathway_summary table from data preparation
  # For now, return a placeholder

  message("Pathway analysis not yet implemented")
  message("Pathway: ", pathway)

  return(data.frame(
    message = "Feature coming soon after pathway enrichment integration"
  ))
}


#' Analyze Gene Family
#'
#' Analyze and compare members of a gene family
#'
#' @param family_name Character. Gene family name or custom list
#' @param members Character vector. Optional custom gene list
#' @param interaction_type Character. Type of interaction to analyze
#'
#' @return List with family analysis results
#' @export
#'
#' @examples
#' \dontrun{
#' # Analyze m6A RNA modification genes
#' m6a_analysis <- gpdb_gene_family("m6A",
#'   members = c("METTL3", "METTL14", "ALKBH5", "FTO")
#' )
#' }
gpdb_gene_family <- function(family_name,
                             members = NULL,
                             interaction_type = "all") {
  # Predefined gene families
  families <- list(
    "m6A" = c(
      "METTL3", "METTL14", "WTAP", "ALKBH5", "FTO",
      "YTHDF1", "YTHDF2", "YTHDF3", "IGF2BP1", "IGF2BP3"
    ),
    "histone_methylation" = c("EZH2", "KMT2A", "KMT2D", "KDM1A", "KDM6A", "KDM5A"),
    "splicing" = c("SF3B1", "SF3B4", "U2AF1", "U2AF2", "SRSF1", "SRSF2")
  )

  # Get gene list
  if (is.null(members)) {
    if (family_name %in% names(families)) {
      members <- families[[family_name]]
    } else {
      stop("Unknown family: ", family_name, ". Please provide 'members' argument.",
        call. = FALSE
      )
    }
  }

  members <- .gpdb_format_genes(members)

  message("Analyzing gene family: ", family_name)
  message("Members: ", paste(members, collapse = ", "))

  # Get effects for each member
  member_effects <- list()
  for (gene in members) {
    result <- gpdb_what_happens(gene, aggregate = TRUE)
    if (!is.null(result$all_effects) && nrow(result$all_effects) > 0) {
      member_effects[[gene]] <- result
    }
  }

  # Find common and specific targets
  all_targets <- lapply(member_effects, function(x) x$all_effects$target_gene)

  if (length(all_targets) >= 2) {
    common_targets <- Reduce(intersect, all_targets)
  } else {
    common_targets <- character(0)
  }

  result <- list(
    family = family_name,
    members = members,
    n_members_with_data = length(member_effects),
    member_effects = member_effects,
    common_targets = common_targets,
    n_common = length(common_targets)
  )

  message("Found data for ", length(member_effects), "/", length(members), " members")
  message("Common targets: ", length(common_targets))

  return(result)
}
