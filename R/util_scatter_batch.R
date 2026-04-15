##################################################
## Batch Encasing Computation — ellipsoid only
##################################################

#' Compute ellipsoid encasing for multiple groups in single subprocess call
#' @param filenm Output JSON filename
#' @param type Encasing type (always "ellipse")
#' @param groups_json JSON string with group info
#' @param level Confidence level (default 0.95)
#' @param omics Omics type (default "NA")
#' @return JSON filename
#' @export
ComputeEncasingBatch <- function(filenm, type, groups_json, level = 0.95, omics = "NA") {
  tryCatch({
    level <- as.numeric(level)

    if (!file.exists("pos.xyz.qs")) {
      sink(filenm); cat("{}"); sink()
      return(filenm)
    }
    pos.xyz <- qs::qread("pos.xyz.qs")

    groups_list <- RJSONIO::fromJSON(groups_json)
    if (is.data.frame(groups_list)) {
      groups_list <- split(groups_list, seq_len(nrow(groups_list)))
    }

    # Parse groups and extract coords in master
    group_data <- lapply(seq_along(groups_list), function(i) {
      g <- groups_list[[i]]
      if (is.character(g)) { gn <- unname(g["group"]); ids <- unname(g["ids"]) }
      else if (is.data.frame(g)) { gn <- g$group[1]; ids <- g$ids[1] }
      else { gn <- g$group; ids <- g$ids }
      nms <- strsplit(ids, "; ")[[1]]
      inx <- rownames(pos.xyz) %in% nms
      list(group = gn, coords = as.matrix(pos.xyz[inx, c(1:3)]))
    })

    # Compute all ellipsoids in single subprocess
    result_list <- rsclient_isolated_exec(
      func_body = function(input_data) {
        Sys.setenv(RGL_USE_NULL = TRUE)
        lapply(input_data$groups, function(g) {
          coords <- g$coords
          if (nrow(coords) < 4) return(list(group = g$group, mesh = list(), error = "Insufficient points"))
          tryCatch({
            pos <- cov(coords, y = NULL, use = "everything")
            center <- colMeans(coords)
            t_val <- sqrt(qchisq(input_data$level, 3))
            mesh <- list()
            mesh[[1]] <- rgl::ellipse3d(x = as.matrix(pos), centre = center, t = t_val)
            list(group = g$group, mesh = mesh, error = NULL)
          }, error = function(e) {
            list(group = g$group, mesh = list(), error = e$message)
          })
        })
      },
      input_data = list(groups = group_data, level = level),
      packages = c("rgl", "qs"),
      timeout = 120,
      output_type = "qs"
    )

    if (!is.list(result_list) || isFALSE(result_list$success)) {
      sink(filenm); cat("{}"); sink()
    } else {
      sink(filenm); cat(RJSONIO::toJSON(result_list)); sink()
    }
  }, error = function(e) {
    message("[ComputeEncasingBatch] ", e$message)
    sink(filenm); cat("{}"); sink()
  })
  return(filenm)
}
