#' Arrow Shadow Save Utility Functions
#'
#' This file provides helper functions for the Arrow migration.
#' The shadow_save() function writes data in both .qs and .arrow formats
#' to enable gradual migration from qs to Apache Arrow.
#'
#' @author MicrobiomeAnalyst Team

#' Shadow save: writes both .qs and .arrow formats
#'
#' This function saves data in the original .qs format and additionally
#' creates an .arrow (Feather) file for Java-side reading. The .qs file
#' maintains backward compatibility while the .arrow file enables
#' direct Java reading without RServe.
#'
#' @param obj The R object to save (data.frame, matrix, or list)
#' @param file The file path (should end with .qs)
#' @return Invisible NULL. Called for side effects.
#' @export
shadow_save <- function(obj, file) {
  # Original qs save (keep for backward compatibility)
  qs::qsave(obj, file)

  # Arrow shadow save
  arrow_path <- sub("\\.qs$", ".arrow", file)
  if (!grepl("\\.arrow$", arrow_path)) {
    arrow_path <- paste0(file, ".arrow")
  }

  # Convert to data.frame if possible, then save as feather
  tryCatch({
    if (is.data.frame(obj) || is.matrix(obj)) {
      df <- as.data.frame(obj)
      arrow::write_feather(df, arrow_path, compression = "uncompressed")
    } else if (is.list(obj) && !inherits(obj, "phyloseq")) {
      # For simple lists with data.frame components, save each component
      for (nm in names(obj)) {
        if (is.data.frame(obj[[nm]]) || is.matrix(obj[[nm]])) {
          component_path <- sub("\\.arrow$", paste0("_", nm, ".arrow"), arrow_path)
          arrow::write_feather(as.data.frame(obj[[nm]]), component_path, compression = "uncompressed")
        }
      }
    }
    # Skip phyloseq and other complex S4 objects for now
  }, error = function(e) {
    warning(paste("Arrow shadow save failed for", file, ":", e$message))
  })

  invisible(NULL)
}
