##################################################
## R script for MicrobiomeAnalyst
## Description: Arrow utilities for zero-copy data exchange with Java
## Author: MicrobiomeAnalyst Team
## Part of the Rserve/qs to Apache Arrow migration
##
## Key Functions:
##   shadow_save()  - Saves .qs + .arrow (for backward compatibility)
##   arrow_save()   - Saves .arrow only (for new code paths)
##
## Safe-Handshake Pattern:
##   1. unlink() existing file to prevent file-lock conflicts
##   2. write_feather() to create new Arrow file
##   3. normalizePath(mustWork=TRUE) to verify file is accessible
###################################################

# ==================== Internal Helper ====================

#' Internal: Write Arrow file with safe-handshake
#'
#' Core function that handles the actual Arrow write with all safety measures.
#' Used internally by shadow_save() and arrow_save().
#'
#' @param df Data frame to write (already prepared with row_names_id if needed)
#' @param arrow_path Path for the Arrow file
#' @return The verified absolute path, or NULL on failure
#' @keywords internal
.write_arrow_internal <- function(df, arrow_path) {
    # Convert factors to character for Arrow compatibility
    for (col in names(df)) {
        if (is.factor(df[[col]])) {
            df[[col]] <- as.character(df[[col]])
        }
    }

    # CRITICAL: Remove existing file to prevent file-lock race conditions
    # (Java may still have the old file memory-mapped)
    if (file.exists(arrow_path)) {
        unlink(arrow_path)
        Sys.sleep(0.01)
    }

    # Write Arrow file
    arrow::write_feather(df, arrow_path, compression = "uncompressed")

    # SAFE-HANDSHAKE: Brief delay then verify with normalizePath
    # This blocks until the OS confirms the file is accessible
    Sys.sleep(0.02)
    verified_path <- base::normalizePath(arrow_path, mustWork = TRUE)
    return(verified_path)
}

#' Internal: Prepare data frame with row_names_id column
#'
#' @param obj Matrix or data.frame
#' @return Data frame with row_names_id as first column (if rownames exist)
#' @keywords internal
.prepare_df_with_rownames <- function(obj) {
    df <- as.data.frame(obj, stringsAsFactors = FALSE)
    rn <- rownames(obj)
    # Only add row_names_id if rownames are meaningful (not just 1:n)
    if (!is.null(rn) && length(rn) > 0 && !all(rn == as.character(seq_len(nrow(df))))) {
        df <- cbind(row_names_id = as.character(rn), df)
    }
    return(df)
}

# ==================== Public API ====================

#' Shadow save: Saves both .qs (for R) and .arrow (for Java)
#'
#' Primary function for saving data during migration. Maintains backward
#' compatibility with existing R code while enabling zero-copy Java access.
#'
#' ALWAYS preserves rownames as the first column (row_names_id) in Arrow files.
#'
#' @param obj The object to save (matrix, data.frame, or list of these)
#' @param file The .qs file path (Arrow path is auto-derived)
#' @return The verified absolute path to the Arrow file, or NULL if not applicable
#' @export
#'
#' @examples
#' \dontrun{
#' # Basic usage
#' shadow_save(result_matrix, "results.qs")
#'
#' # Java then reads from "results.arrow" via:
#' # DataTable table = sb.getOrCreateDataTable("results", arrowPath);
#' }
shadow_save <- function(obj, file) {
    # Always save to qs for R compatibility
    qs::qsave(obj, file)

    # Generate Arrow path from qs path
    arrow_path <- sub("\\.qs$", ".arrow", file)
    if (!grepl("\\.arrow$", arrow_path)) {
        arrow_path <- paste0(file, ".arrow")
    }

    tryCatch({
        if (is.matrix(obj) || is.data.frame(obj)) {
            # Single matrix/data.frame
            df <- .prepare_df_with_rownames(obj)
            return(.write_arrow_internal(df, arrow_path))

        } else if (is.list(obj) && !inherits(obj, "phyloseq") && !is.null(names(obj))) {
            # Named list with data.frame/matrix components
            # Save each component as separate Arrow file
            for (nm in names(obj)) {
                if (is.data.frame(obj[[nm]]) || is.matrix(obj[[nm]])) {
                    component_path <- sub("\\.arrow$", paste0("_", nm, ".arrow"), arrow_path)
                    df <- .prepare_df_with_rownames(obj[[nm]])
                    .write_arrow_internal(df, component_path)
                }
            }
        }
        # Skip phyloseq and other complex S4 objects (qs handles them)
    }, error = function(e) {
        warning(paste("Arrow shadow save failed:", e$message))
    })

    return(NULL)
}

#' Arrow save: Saves directly to .arrow (no .qs backup)
#'
#' Use this for new code paths where qs compatibility is not needed,
#' or when you only need the Arrow file (e.g., result matrices for Java display).
#'
#' @param obj The R object to save (data.frame or matrix)
#' @param file The file path (will add .arrow extension if missing)
#' @return The verified absolute path to the Arrow file
#' @export
#'
#' @examples
#' \dontrun{
#' arrow_save(confusion_matrix, "rf_confusion_mat.arrow")
#' arrow_save(result_mat, "res_mat.arrow")
#' }
arrow_save <- function(obj, file) {
    # Ensure .arrow extension
    if (!grepl("\\.arrow$", file)) {
        file <- paste0(file, ".arrow")
    }

    if (!is.data.frame(obj) && !is.matrix(obj)) {
        stop("arrow_save only supports data.frame and matrix objects")
    }

    tryCatch({
        df <- .prepare_df_with_rownames(obj)
        return(.write_arrow_internal(df, file))
    }, error = function(e) {
        stop(paste("Arrow save failed for", file, ":", e$message))
    })
}

# ==================== Utility Functions ====================

#' Safe column extraction with name-first, index-fallback strategy
#'
#' Provides safe column extraction during Arrow migration by:
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
                    paste(colnames(tab), collapse = ", ")))
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
                        context, paste(missing, collapse = ", "),
                        paste(colnames(tab), collapse = ", ")))
        return(FALSE)
    }
    return(TRUE)
}
