#' Get database connection
#' @keywords internal
.gpdb_get_connection <- function() {
  db_base <- Sys.getenv("APIKIT_DB_PATH")

  if (db_base == "") {
    stop("APIKIT_DB_PATH not set. Configure: Sys.setenv(APIKIT_DB_PATH='/path')",
      call. = FALSE
    )
  }

  db_path <- file.path(db_base, "gpsadb", "processed", "gpsadb.db")

  if (!file.exists(db_path)) {
    stop("Database not found: ", db_path, "\n",
      "Please run data preparation script first.",
      call. = FALSE
    )
  }

  # Reuse connection if available and valid
  if (!is.null(.gpdb_env$con) &&
    !is.null(.gpdb_env$db_path) &&
    .gpdb_env$db_path == db_path) {
    # Test if connection is still valid
    test <- try(DBI::dbGetQuery(.gpdb_env$con, "SELECT 1"), silent = TRUE)
    if (!inherits(test, "try-error")) {
      return(.gpdb_env$con)
    }
  }

  # Create new connection
  .gpdb_env$con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
  .gpdb_env$db_path <- db_path

  return(.gpdb_env$con)
}

#' Get processed data directory
#' @keywords internal
.gpdb_get_processed_dir <- function() {
  db_base <- Sys.getenv("APIKIT_DB_PATH")

  if (db_base == "") {
    stop("APIKIT_DB_PATH not set. Configure: Sys.setenv(APIKIT_DB_PATH='/path')",
      call. = FALSE
    )
  }

  processed_dir <- file.path(db_base, "gpsadb", "processed")

  if (!dir.exists(processed_dir)) {
    stop("Processed data directory not found: ", processed_dir,
      call. = FALSE
    )
  }

  return(processed_dir)
}

#' Load qs file
#' @keywords internal
.gpdb_load_qs <- function(file_path) {
  if (!file.exists(file_path)) {
    stop("File not found: ", file_path, call. = FALSE)
  }

  qs::qread(file_path)
}

#' Format gene names
#' @keywords internal
.gpdb_format_genes <- function(genes) {
  # Remove NA, empty strings, and trim whitespace
  genes <- trimws(genes)
  genes <- genes[!is.na(genes) & genes != ""]

  # Convert to uppercase for consistency
  toupper(genes)
}

#' Validate dataset ID
#' @keywords internal
.gpdb_validate_dataset_id <- function(dataset_id) {
  if (!grepl("^D\\d+$", dataset_id)) {
    stop("Invalid dataset ID format: ", dataset_id,
      ". Expected format: D10001",
      call. = FALSE
    )
  }

  TRUE
}

#' SQL-safe string
#' @keywords internal
.gpdb_sql_safe <- function(x) {
  gsub("'", "''", as.character(x))
}

#' Generate text summary
#' @keywords internal
.gpdb_generate_summary <- function(gene, n_datasets, top_up, top_down, pathways) {
  summary <- paste0(
    "Gene: ", gene, "\n",
    "Evidence: ", n_datasets, " perturbation datasets\n\n"
  )

  if (!is.null(top_up) && nrow(top_up) > 0) {
    summary <- paste0(summary, "Top upregulated targets:\n")
    for (i in 1:min(5, nrow(top_up))) {
      summary <- paste0(
        summary, "  - ", top_up$target_gene[i],
        " (", round(top_up$logfc_mean[i], 2), "-fold, ",
        top_up$n_datasets[i], " datasets)\n"
      )
    }
    summary <- paste0(summary, "\n")
  }

  if (!is.null(top_down) && nrow(top_down) > 0) {
    summary <- paste0(summary, "Top downregulated targets:\n")
    for (i in 1:min(5, nrow(top_down))) {
      summary <- paste0(
        summary, "  - ", top_down$target_gene[i],
        " (", round(top_down$logfc_mean[i], 2), "-fold, ",
        top_down$n_datasets[i], " datasets)\n"
      )
    }
    summary <- paste0(summary, "\n")
  }

  if (!is.null(pathways) && nrow(pathways) > 0) {
    summary <- paste0(
      summary, "Enriched pathways: ",
      paste(head(pathways$pathway, 3), collapse = ", "),
      "\n"
    )
  }

  return(summary)
}
