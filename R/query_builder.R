#' SQL Query Builder for GenePerturbR
#'
#' Unified SQL query construction system to eliminate code duplication
#' and improve query performance
#'
#' @keywords internal

#' Build SQL Query
#'
#' @param table Character. Table name
#' @param select Character. Columns to select
#' @param filters List. Named list of filter conditions
#' @param order_by Character. ORDER BY clause
#' @param limit Integer. LIMIT clause
#'
#' @return Character. SQL query string
#' @keywords internal
.gpdb_build_query <- function(table = "gene_effects_agg",
                              select = "*",
                              filters = list(),
                              order_by = NULL,
                              limit = NULL) {
  query <- paste0("SELECT ", select, " FROM ", table, " WHERE 1=1")

  # Perturbed gene filter
  if (!is.null(filters$perturbed_gene)) {
    query <- paste0(query, " AND perturbed_gene = '", .gpdb_sql_safe(filters$perturbed_gene), "'")
  }

  # Target gene filter
  if (!is.null(filters$target_gene)) {
    query <- paste0(query, " AND target_gene = '", .gpdb_sql_safe(filters$target_gene), "'")
  }

  # Target gene IN list
  if (!is.null(filters$target_gene_in)) {
    gene_list <- paste0("'", sapply(filters$target_gene_in, .gpdb_sql_safe), "'", collapse = ", ")
    query <- paste0(query, " AND target_gene IN (", gene_list, ")")
  }

  # Minimum datasets
  if (!is.null(filters$min_datasets)) {
    query <- paste0(query, " AND n_datasets >= ", filters$min_datasets)
  }

  # Confidence filter
  if (!is.null(filters$min_confidence)) {
    if (filters$min_confidence == "high") {
      query <- paste0(query, " AND confidence = 'high'")
    } else if (filters$min_confidence == "medium") {
      query <- paste0(query, " AND confidence IN ('high', 'medium')")
    }
  }

  # Effect size filter
  if (!is.null(filters$min_effect_size)) {
    query <- paste0(query, " AND ABS(logfc_mean) >= ", filters$min_effect_size)
  }

  # Direction filter
  if (!is.null(filters$direction)) {
    if (filters$direction == "up") {
      query <- paste0(query, " AND logfc_mean > 0")
    } else if (filters$direction == "down") {
      query <- paste0(query, " AND logfc_mean < 0")
    }
  }

  # Context filters
  if (!is.null(filters$cell_line)) {
    query <- paste0(query, " AND celllines LIKE '%", .gpdb_sql_safe(filters$cell_line), "%'")
  }

  if (!is.null(filters$tissue)) {
    query <- paste0(query, " AND tissues LIKE '%", .gpdb_sql_safe(filters$tissue), "%'")
  }

  # Add ORDER BY
  if (!is.null(order_by)) {
    query <- paste0(query, " ORDER BY ", order_by)
  }

  # Add LIMIT
  if (!is.null(limit)) {
    query <- paste0(query, " LIMIT ", limit)
  }

  return(query)
}

#' Execute Query with Error Handling
#'
#' @param query Character. SQL query
#' @param error_msg Character. Custom error message
#'
#' @return Data frame. Query results
#' @keywords internal
.gpdb_execute_query <- function(query, error_msg = "Database query failed") {
  con <- .gpdb_get_connection()

  tryCatch(
    {
      result <- DBI::dbGetQuery(con, query)
      return(as.data.frame(result))
    },
    error = function(e) {
      stop(error_msg, ": ", e$message, call. = FALSE)
    }
  )
}

#' Query Dataset Metadata in Batch
#'
#' @param dataset_ids Character vector. Dataset IDs
#'
#' @return Data frame. Metadata for all datasets
#' @keywords internal
.gpdb_batch_query_metadata <- function(dataset_ids) {
  con <- .gpdb_get_connection()

  ids_sql <- paste0("'", sapply(dataset_ids, .gpdb_sql_safe), "'", collapse = ", ")
  query <- paste0(
    "SELECT Dataset as dataset_id, pbgene as gene, pb_ENSEMBL as ensembl_id, ",
    "gene_biotype, method, nSample as n_samples, CellLineName as cell_line, ",
    "TissueSite as tissue, Datasource, accession ",
    "FROM datasets WHERE Dataset IN (", ids_sql, ")"
  )

  DBI::dbGetQuery(con, query)
}
