#' Arrow Parity Test Functions
#'
#' This file provides test functions to verify that .qs and .arrow files
#' produce identical data, ensuring data integrity during migration.
#'
#' @author MicrobiomeAnalyst Team

#' Test parity between qs and arrow file outputs
#'
#' Compares data from a .qs file with its corresponding .arrow file
#' to verify that the shadow_save() function preserves data integrity.
#'
#' @param qs_path Path to the .qs file
#' @return TRUE if parity verified, stops with error otherwise
#' @export
test_arrow_parity <- function(qs_path) {
  if (!file.exists(qs_path)) {
    stop("QS file not found: ", qs_path)
  }

  arrow_path <- sub("\\.qs$", ".arrow", qs_path)
  if (!file.exists(arrow_path)) {
    stop("Arrow file not found: ", arrow_path)
  }

  # Read both formats
  qs_data <- qs::qread(qs_path)
  arrow_data <- arrow::read_feather(arrow_path)

  # Convert qs data to comparable format
  if (is.matrix(qs_data)) {
    qs_data <- as.data.frame(qs_data)
  }

  # Dimension checks
  if (nrow(qs_data) != nrow(arrow_data)) {
    stop(sprintf("Row count mismatch: qs=%d, arrow=%d", nrow(qs_data), nrow(arrow_data)))
  }
  if (ncol(qs_data) != ncol(arrow_data)) {
    stop(sprintf("Column count mismatch: qs=%d, arrow=%d", ncol(qs_data), ncol(arrow_data)))
  }

  # Column name check
  if (!all(names(qs_data) == names(arrow_data))) {
    stop("Column names mismatch")
  }

  # Value comparison for numeric columns
  for (col in names(qs_data)) {
    if (is.numeric(qs_data[[col]])) {
      diff <- max(abs(qs_data[[col]] - arrow_data[[col]]), na.rm = TRUE)
      if (diff > 1e-10) {
        stop(sprintf("Values mismatch in column '%s': max diff = %e", col, diff))
      }
    }
  }

  message("PARITY VERIFIED: ", basename(qs_path))
  return(TRUE)
}

#' Run parity tests on all .qs files in a directory
#'
#' Iterates through all .qs files in the specified directory and
#' runs parity tests against their corresponding .arrow files.
#'
#' @param dir_path Directory containing .qs files
#' @return Invisible logical vector of test results
#' @export
run_all_parity_tests <- function(dir_path) {
  qs_files <- list.files(dir_path, pattern = "\\.qs$", full.names = TRUE)
  results <- sapply(qs_files, function(f) {
    tryCatch({
      test_arrow_parity(f)
    }, error = function(e) {
      message("FAILED: ", basename(f), " - ", e$message)
      FALSE
    })
  })

  passed <- sum(results)
  total <- length(results)
  message(sprintf("\n=== PARITY SUMMARY: %d/%d passed ===", passed, total))
  return(invisible(results))
}
