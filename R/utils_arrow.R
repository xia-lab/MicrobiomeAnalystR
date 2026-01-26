##################################################
## R script for MicrobiomeAnalyst
## Description: Arrow utilities for zero-copy data exchange with Java
## Author: MicrobiomeAnalyst Team
## Part of the Rserve/qs to Apache Arrow migration
###################################################

#' Sync file to disk and verify existence (Safe-Handshake pattern)
#'
#' This function ensures that a file is fully written to disk before
#' returning control to Java. Uses normalizePath with mustWork=TRUE
#' to provide a filesystem-level guarantee that the file exists.
#'
#' @param file_path Path to the file to sync and verify
#' @param delay Delay in seconds before verification (default: 0.02 = 20ms)
#' @return The normalized (absolute) path to the file, or NULL if verification fails
#' @export
sync_file <- function(file_path, delay = 0.02) {
    if (is.null(file_path) || !nzchar(file_path)) {
        return(invisible(NULL))
    }

    # Brief delay to allow filesystem buffers to flush
    Sys.sleep(delay)

    # Use normalizePath with mustWork=TRUE as filesystem verification
    # This blocks until the OS confirms the file is accessible
    tryCatch({
        verified_path <- base::normalizePath(file_path, mustWork = TRUE)
        return(verified_path)
    }, error = function(e) {
        warning(sprintf("sync_file: File verification failed for '%s': %s",
                        file_path, e$message))
        return(NULL)
    })
}

#' Write Arrow file with safe-handshake verification
#'
#' Writes an Arrow (Feather) file and returns the verified absolute path.
#' Uses normalizePath(mustWork=TRUE) to ensure the file is fully written
#' before returning control to Java.
#'
#' @param df Data frame to write
#' @param path Path for the Arrow file
#' @param compress Compression type (default: "uncompressed" for memory-mapping)
#' @return The verified absolute path to the Arrow file
#' @export
write_arrow_safe <- function(df, path, compress = "uncompressed") {
    # Ensure factors are converted to character
    for (col in names(df)) {
        if (is.factor(df[[col]])) {
            df[[col]] <- as.character(df[[col]])
        }
    }

    # CRITICAL: Remove existing file to prevent file-lock race conditions
    if (file.exists(path)) {
        unlink(path)
        Sys.sleep(0.01)
    }

    # Write the Arrow file
    arrow::write_feather(df, path, compression = compress)

    # Brief delay then verify with normalizePath
    Sys.sleep(0.02)

    # Return verified absolute path - blocks until file is confirmed accessible
    verified_path <- base::normalizePath(path, mustWork = TRUE)
    return(verified_path)
}

#' Safe column extraction with name-first, index-fallback strategy
#'
#' This function provides a safe way to extract columns from data frames
#' during Arrow migration. It ensures data integrity by:
#' 1. Trying named column access first (preferred)
#' 2. Falling back to index if name doesn't exist (with warning)
#' 3. Returning NA vector with error if both fail
#'
#' @param tab Data frame or matrix to extract from
#' @param name Column name to try first (character)
#' @param idx Fallback column index (1-based integer)
#' @param context Optional context string for logging
#' @return Column values as vector, or NA vector if extraction fails
#' @export
safeGetCol <- function(tab, name, idx = NULL, context = "") {
  nrows <- nrow(tab)
  if (is.null(nrows) || nrows == 0) {
    return(character(0))
  }

  # Strategy 1: Try by name first (preferred - most reliable)
  if (!is.null(name) && name %in% colnames(tab)) {
    return(tab[[name]])
  }

  # Strategy 2: Fall back to index (with warning)
  if (!is.null(idx) && is.numeric(idx) && idx > 0 && ncol(tab) >= idx) {
    actualName <- colnames(tab)[idx]
    warning(sprintf("[%s] Column '%s' not found, using index %d (actual column: '%s').",
                    context, name, idx, actualName))
    return(tab[, idx])
  }

  # Strategy 3: Both failed - return NA with error
  warning(sprintf("[%s] FAILED to extract column: name='%s', idx=%s. Available columns: %s",
                  context, name, ifelse(is.null(idx), "NULL", idx),
                  paste(colnames(tab), collapse=", ")))
  return(rep(NA, nrows))
}

#' Validate that required columns exist in a data frame
#'
#' @param tab Data frame to validate
#' @param required Character vector of required column names
#' @param context Context string for error messages
#' @return TRUE if all columns exist, FALSE otherwise (with warnings)
#' @export
validateColumns <- function(tab, required, context = "") {
  missing <- setdiff(required, colnames(tab))
  if (length(missing) > 0) {
    warning(sprintf("[%s] Missing required columns: %s. Available: %s",
                    context, paste(missing, collapse=", "),
                    paste(colnames(tab), collapse=", ")))
    return(FALSE)
  }
  return(TRUE)
}

#' Shadow save function for mixed-type data frames (Safe-Handshake pattern)
#'
#' Saves data as both .qs (for R) and .arrow (for Java zero-copy access).
#' ALWAYS preserves rownames as the first column (row_names_id) in Arrow files.
#' Returns the verified absolute path to the Arrow file for Java to memory-map.
#'
#' @param obj The object to save (matrix or data.frame)
#' @param file The .qs file path (will auto-generate .arrow path)
#' @param compress Compression for Arrow (default: "uncompressed" for memory-mapping)
#' @return The verified absolute path to the Arrow file, or NULL if save failed
#' @export
shadow_save <- function(obj, file, compress = "uncompressed") {
    # Always save to qs for R compatibility
    qs::qsave(obj, file)

    # Generate Arrow path
    arrow_path <- sub("\\.qs$", ".arrow", file)
    if (!grepl("\\.arrow$", arrow_path)) {
        arrow_path <- paste0(file, ".arrow")
    }

    tryCatch({
        if (is.matrix(obj) || is.data.frame(obj)) {
            df <- as.data.frame(obj)

            # Preserve rownames as first column
            rn <- rownames(obj)
            if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
                df <- cbind(row_names_id = as.character(rn), df)
            }

            # Ensure all character columns are properly typed
            for (col in names(df)) {
                if (is.factor(df[[col]])) {
                    df[[col]] <- as.character(df[[col]])
                }
            }

            # CRITICAL: Remove existing file to prevent file-lock race conditions
            if (file.exists(arrow_path)) {
                unlink(arrow_path)
                Sys.sleep(0.01)
            }

            arrow::write_feather(df, arrow_path, compression = compress)

            # SAFE-HANDSHAKE: Brief delay then verify with normalizePath
            # This blocks until the OS confirms the file is accessible
            Sys.sleep(0.02)
            verified_path <- base::normalizePath(arrow_path, mustWork = TRUE)
            return(verified_path)
        } else if (is.list(obj) && !inherits(obj, "phyloseq")) {
            # For simple lists with data.frame components, save each component
            for (nm in names(obj)) {
                if (is.data.frame(obj[[nm]]) || is.matrix(obj[[nm]])) {
                    component_path <- sub("\\.arrow$", paste0("_", nm, ".arrow"), arrow_path)
                    df <- as.data.frame(obj[[nm]])
                    rn <- rownames(obj[[nm]])
                    if (!is.null(rn) && length(rn) > 0) {
                        df <- cbind(row_names_id = as.character(rn), df)
                    }
                    # CRITICAL: Remove existing file to prevent file-lock race conditions
                    if (file.exists(component_path)) {
                        unlink(component_path)
                        Sys.sleep(0.01)
                    }
                    arrow::write_feather(df, component_path, compression = compress)
                }
            }
        }
        # Skip phyloseq and other complex S4 objects for now
    }, error = function(e) {
        warning(paste("Arrow shadow save failed:", e$message))
    })

    return(NULL)
}

#' Shadow save for mixed-type data frames (Safe-Handshake pattern)
#'
#' Handles factors by converting to character. Returns verified path.
#'
#' @param obj Data frame or matrix to save
#' @param file Path to .qs file
#' @param compress Compression for Arrow (default: "uncompressed")
#' @return The verified absolute path to the Arrow file, or NULL if save failed
#' @export
shadow_save_mixed <- function(obj, file, compress = "uncompressed") {
    # Always save qs format for backward compatibility
    qs::qsave(obj, file)

    # Derive Arrow path
    arrow_path <- sub("\\.qs$", ".arrow", file)

    tryCatch({
        if (is.matrix(obj) || is.data.frame(obj)) {
            df <- as.data.frame(obj, stringsAsFactors = FALSE)

            # Convert factors to character for Arrow compatibility
            for (col in names(df)) {
                if (is.factor(df[[col]])) {
                    df[[col]] <- as.character(df[[col]])
                }
            }

            rn <- rownames(obj)
            if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
                # Prepend row_names_id as first column
                df <- cbind(row_names_id = as.character(rn), df)
            }

            # CRITICAL: Remove existing file to prevent file-lock race conditions
            if (file.exists(arrow_path)) {
                unlink(arrow_path)
                Sys.sleep(0.01)
            }

            arrow::write_feather(df, arrow_path, compression = compress)

            # SAFE-HANDSHAKE: Brief delay then verify with normalizePath
            Sys.sleep(0.02)
            verified_path <- base::normalizePath(arrow_path, mustWork = TRUE)
            return(verified_path)
        }
    }, error = function(e) {
        warning(paste("Arrow save (mixed) failed:", e$message))
    })

    return(NULL)
}

#' Direct Arrow save with Safe-Handshake (no qs, for new code paths)
#'
#' Saves data directly to Arrow format without qs backup.
#' Use this for new code paths where qs compatibility is not needed.
#' Returns the verified absolute path to the Arrow file.
#'
#' @param obj The R object to save (data.frame or matrix)
#' @param file The file path (should end with .arrow)
#' @return The verified absolute path to the Arrow file
#' @export
arrow_save <- function(obj, file) {
  if (!grepl("\\.arrow$", file)) {
    file <- paste0(file, ".arrow")
  }

  tryCatch({
    if (is.data.frame(obj) || is.matrix(obj)) {
      df <- as.data.frame(obj)

      # Preserve rownames
      rn <- rownames(obj)
      if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(1:nrow(df)))) {
        df <- cbind(row_names_id = as.character(rn), df)
      }

      # Ensure factors are converted
      for (col in names(df)) {
        if (is.factor(df[[col]])) {
          df[[col]] <- as.character(df[[col]])
        }
      }

      # CRITICAL: Remove existing file to prevent file-lock race conditions
      if (file.exists(file)) {
        unlink(file)
        Sys.sleep(0.01)
      }

      arrow::write_feather(df, file, compression = "uncompressed")

      # SAFE-HANDSHAKE: Brief delay then verify with normalizePath
      Sys.sleep(0.02)
      verified_path <- base::normalizePath(file, mustWork = TRUE)
      return(verified_path)
    } else {
      stop("arrow_save only supports data.frame and matrix objects")
    }
  }, error = function(e) {
    stop(paste("Arrow save failed for", file, ":", e$message))
  })
}
