#' List Datasets
#'
#' @param gene Character. Gene symbol
#' @param cell_line Character. Cell line name
#' @param tissue Character. Tissue type
#' @param method Character. Perturbation method
#' @param min_quality Numeric. Not implemented in current version
#' @param quality_tier Character. Not implemented in current version
#'
#' @return Data frame with dataset information
#' @export
gpdb_list_datasets <- function(gene = NULL, cell_line = NULL, tissue = NULL,
                               method = NULL, min_quality = NULL, quality_tier = NULL) {
  query <- "SELECT Dataset as dataset_id, pbgene as gene, pb_ENSEMBL as ensembl_id,
            gene_biotype, method, nSample as n_samples, CellLineName as cell_line,
            TissueSite as tissue, Datasource, accession
            FROM datasets WHERE 1=1"

  if (!is.null(gene)) query <- paste0(query, " AND pbgene = '", .gpdb_sql_safe(gene), "'")
  if (!is.null(cell_line)) query <- paste0(query, " AND CellLineName = '", .gpdb_sql_safe(cell_line), "'")
  if (!is.null(tissue)) query <- paste0(query, " AND TissueSite = '", .gpdb_sql_safe(tissue), "'")
  if (!is.null(method)) query <- paste0(query, " AND method = '", .gpdb_sql_safe(method), "'")

  result <- .gpdb_execute_query(paste0(query, " ORDER BY Dataset DESC"))

  if (nrow(result) == 0) {
    message("No datasets found")
    return(data.frame())
  }

  message("Found ", nrow(result), " datasets")
  return(as.data.frame(result))
}


#' Get Dataset Information
#'
#' @param dataset_id Character. Dataset ID
#' @return List with metadata
#' @export
gpdb_get_info <- function(dataset_id) {
  .gpdb_validate_dataset_id(dataset_id)

  query <- paste0(
    "SELECT Dataset as dataset_id, pbgene as gene, pb_ENSEMBL as ensembl_id, ",
    "gene_biotype, method, nSample as n_samples, CellLineName as cell_line, ",
    "TissueSite as tissue, Datasource, accession ",
    "FROM datasets WHERE Dataset = '", .gpdb_sql_safe(dataset_id), "'"
  )

  result <- .gpdb_execute_query(query)
  if (nrow(result) == 0) stop("Dataset not found: ", dataset_id, call. = FALSE)

  return(as.list(result[1, ]))
}


#' Load Complete Dataset (Expression + Metadata)
#'
#' Load expression matrix AND sample metadata together.
#' Expression matrix has gene symbols as rownames (pre-processed).
#'
#' @param dataset_id Character. Dataset ID
#' @param normalize Logical. Apply log2(CPM+1) normalization (default FALSE, keep raw counts)
#'
#' @return List containing:
#'   \itemize{
#'     \item expression: Data frame with gene symbols as rownames (genes x samples)
#'     \item metadata: Sample metadata with 'group' column (treatment/control)
#'     \item info: Dataset information with helpful tips
#'   }
#'
#' @details
#' Expression matrix structure:
#' \itemize{
#'   \item rownames = gene symbols (e.g., "TP53", "MYC")
#'   \item columns = samples (e.g., "GSM123456")
#'   \item values = raw counts (or log2(CPM+1) if normalized)
#' }
#'
#' Metadata structure:
#' \itemize{
#'   \item 'group' column = "treatment" or "control"
#'   \item Other columns may include batch info, cell type, etc.
#' }
#'
#' @export
#'
#' @examples
#' \dontrun{
#' # Load dataset
#' data <- gpdb_load_data("D10001")
#'
#' # Check structure
#' head(rownames(data$expression)) # gene symbols
#' head(data$metadata$group) # treatment/control
#' print(data$info) # dataset information
#'
#' # Use for DESeq2 analysis
#' library(DESeq2)
#' dds <- DESeqDataSetFromMatrix(
#'   countData = data$expression,
#'   colData = data$metadata,
#'   design = ~group
#' )
#' dds <- DESeq(dds)
#' }
gpdb_load_data <- function(dataset_id, normalize = FALSE) {
  .gpdb_validate_dataset_id(dataset_id)

  processed_dir <- .gpdb_get_processed_dir()

  # Load expression matrix
  expr_path <- file.path(processed_dir, "expression", paste0(dataset_id, ".qs"))
  if (!file.exists(expr_path)) {
    stop("Expression data not found for: ", dataset_id, call. = FALSE)
  }
  expr_data <- .gpdb_load_qs(expr_path)

  # Load sample metadata
  meta_path <- file.path(processed_dir, "metadata", paste0(dataset_id, ".qs"))
  if (file.exists(meta_path)) {
    meta_data <- .gpdb_load_qs(meta_path)
  } else {
    warning("Metadata not found for ", dataset_id, call. = FALSE)
    meta_data <- NULL
  }

  # Get dataset info
  info <- gpdb_get_info(dataset_id)

  # Convert gene_id to gene_name (使用缓存机制，快速！)
  expr_data <- .gpdb_convert_to_gene_name(expr_data)

  # Normalize if requested
  if (normalize) {
    expr_matrix <- as.matrix(expr_data)
    lib_size <- colSums(expr_matrix)
    cpm <- sweep(expr_matrix, 2, lib_size / 1e6, "/")
    expr_data <- as.data.frame(log2(cpm + 1))
    rownames(expr_data) <- rownames(expr_matrix)
    message("Applied log2(CPM+1) normalization")
  }

  # Add helpful info
  info$expression_format <- "Data frame with gene symbols as rownames"
  info$n_genes <- nrow(expr_data)
  info$n_samples <- ncol(expr_data)

  if (!is.null(meta_data) && "group" %in% names(meta_data)) {
    groups <- table(meta_data$group)
    info$sample_groups <- paste(names(groups), "=", groups, collapse = ", ")
    info$metadata_key_column <- "group (treatment vs control)"
  }

  result <- list(
    expression = as.data.frame(expr_data),
    metadata = if (!is.null(meta_data)) as.data.frame(meta_data) else NULL,
    info = info
  )

  message("Loaded dataset ", dataset_id, ": ", info$gene, " in ", info$cell_line)
  message(
    "  Expression: ", nrow(result$expression), " genes (rownames) × ",
    ncol(result$expression), " samples"
  )
  if (!is.null(result$metadata)) {
    message("  Metadata: ", nrow(result$metadata), " samples")
    if ("group" %in% names(result$metadata)) {
      message("  Groups: ", info$sample_groups)
    }
  }

  return(result)
}


#' Load DEG Results Only
#'
#' Load pre-computed differential expression results.
#' For loading raw data + metadata together, use gpdb_load_data().
#'
#' @param dataset_id Character. Dataset ID
#' @param filter Logical. Keep only significant genes
#' @param padj_cutoff Numeric. P-value cutoff
#' @param logfc_cutoff Numeric. LogFC cutoff
#'
#' @return Data frame with DEG results
#' @export
gpdb_load_deg <- function(dataset_id, filter = FALSE,
                          padj_cutoff = 0.05, logfc_cutoff = 0) {
  .gpdb_validate_dataset_id(dataset_id)

  file_path <- file.path(.gpdb_get_processed_dir(), "deg", paste0(dataset_id, ".qs"))
  if (!file.exists(file_path)) stop("DEG not found for: ", dataset_id, call. = FALSE)

  deg_data <- .gpdb_load_qs(file_path)

  if (filter) {
    deg_data <- deg_data[
      !is.na(deg_data$adj.P.Val) &
        deg_data$adj.P.Val < padj_cutoff &
        abs(deg_data$logFC) > logfc_cutoff,
    ]
    if (nrow(deg_data) > 0) {
      message("Filtered to ", nrow(deg_data), " significant genes")
    }
  }

  return(as.data.frame(deg_data))
}


#' Load Sample Metadata Only
#'
#' Load sample-level metadata.
#' For loading expression + metadata together, use gpdb_load_data().
#'
#' @param dataset_id Character. Dataset ID
#' @return Data frame
#' @export
gpdb_load_metadata <- function(dataset_id) {
  .gpdb_validate_dataset_id(dataset_id)

  file_path <- file.path(.gpdb_get_processed_dir(), "metadata", paste0(dataset_id, ".qs"))
  if (!file.exists(file_path)) stop("Metadata not found for: ", dataset_id, call. = FALSE)

  return(as.data.frame(.gpdb_load_qs(file_path)))
}


#' Batch Load Datasets
#'
#' @param dataset_ids Character vector. Dataset IDs
#' @param type Character. "data" (expression+metadata), "deg", or "metadata"
#' @param show_progress Logical. Show progress bar (default TRUE)
#' @param ... Additional arguments
#' @return Named list
#' @export
gpdb_load_batch <- function(dataset_ids,
                            type = c("data", "deg", "metadata"),
                            show_progress = TRUE,
                            ...) {
  type <- match.arg(type)

  # Batch query metadata first (faster!)
  if (type == "data" || type == "metadata") {
    metadata_info <- .gpdb_batch_query_metadata(dataset_ids)
    message("Pre-loaded metadata for ", nrow(metadata_info), " datasets")
  }

  load_fn <- switch(type,
    "data" = gpdb_load_data,
    "deg" = gpdb_load_deg,
    "metadata" = gpdb_load_metadata
  )

  result <- list()
  message("Loading ", length(dataset_ids), " datasets...")

  pb <- .gpdb_progress_bar(length(dataset_ids), "Loading", show_progress)

  for (i in seq_along(dataset_ids)) {
    tryCatch(
      {
        result[[dataset_ids[i]]] <- load_fn(dataset_ids[i], ...)
      },
      error = function(e) {
        warning("Failed: ", dataset_ids[i], call. = FALSE)
      }
    )
    .gpdb_update_progress(pb, i)
  }

  .gpdb_close_progress(pb)
  message("Loaded ", length(result), "/", length(dataset_ids))

  return(result)
}
