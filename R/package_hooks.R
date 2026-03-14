NULL

# Package environment to store database connection
.gpdb_env <- new.env(parent = emptyenv())

#' @keywords internal
.onLoad <- function(libname, pkgname) {
  # Initialize connection pool
  .gpdb_env$con <- NULL
  .gpdb_env$db_path <- NULL
}

#' @keywords internal
.onUnload <- function(libpath) {
  # Close database connection on unload
  if (!is.null(.gpdb_env$con)) {
    try(DBI::dbDisconnect(.gpdb_env$con), silent = TRUE)
  }
}

#' @keywords internal
.onAttach <- function(libname, pkgname) {
  # Check if APIKIT_DB_PATH is set
  db_base <- Sys.getenv("APIKIT_DB_PATH")

  if (db_base == "") {
    packageStartupMessage(
      "GenePerturbR: APIKIT_DB_PATH not set.\n",
      "Configure in R: Sys.setenv(APIKIT_DB_PATH='/path/to/API_DB')\n",
      "Or add to ~/.Renviron: APIKIT_DB_PATH=/path/to/API_DB"
    )
  } else {
    db_path <- file.path(db_base, "gpsadb", "processed", "gpsadb.db")
    if (!file.exists(db_path)) {
      packageStartupMessage(
        "GenePerturbR: Database not found at: ", db_path, "\n",
        "Current APIKIT_DB_PATH: ", db_base, "\n",
        "Please check path or run: source(system.file('data-raw/prepare_data.R', package='GenePerturbR'))"
      )
    } else {
      n_pairs <- tryCatch(
        {
          con <- DBI::dbConnect(RSQLite::SQLite(), db_path)
          n <- DBI::dbGetQuery(con, "SELECT COUNT(*) as n FROM gene_effects_agg")$n
          DBI::dbDisconnect(con)
          format(n, big.mark = ",")
        },
        error = function(e) "unknown"
      )

      packageStartupMessage(
        "GenePerturbR loaded successfully\n",
        "Database: ", basename(db_path), " (", n_pairs, " gene pairs)"
      )
    }
  }
}
