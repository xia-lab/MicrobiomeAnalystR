##################################################
## R script for MicrobiomeAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## Y. Lu   
## J. Xia, jeff.xia@mcgill.ca
###################################################
#procrustes or diablo
my.reduce.dimension <- function(mbSetObj, reductionOpt= "procrustes", method="globalscore", dimn=10,analysisVar, diabloPar=0.2){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  message("[my.reduce.dimension] reductionOpt=", reductionOpt, " method=", method)
  if(method == ""){
    method="globalscore"
  }
  dimn = as.numeric(dimn);

  d.list = vector("list",length = 2)
  names(d.list) <- c("mic","met")
  omics.type = vector();
  featureNms <- vector();
  combined.res <- list()
  
  omics.type <- c("microbiome","metabolomics")
  
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
  }
  
  d.list[["mic"]] = list()
  # Use normalized/auto-scaled data for integration
  if(micDataType=='ko'){
    d.list[["mic"]][["data.proc"]] = list(OTU = current.proc$mic$data.proc)
  } else {
    # Wrap as named list matching the analyzed taxonomy level
    analyzed_lvl <- names(which(!sapply(phyloseq_objs$res_deAnal, is.null)))
    if(length(analyzed_lvl) == 0) analyzed_lvl <- names(phyloseq_objs$count_tables)
    d.list[["mic"]][["data.proc"]] = setNames(list(current.proc$mic$data.proc), analyzed_lvl[1])
  }
  if(micDataType=='ko'){
    d.list[["mic"]][["comp.res"]]  = current.proc$mic$res_deAnal[,c(3,4,1)]
    names(d.list[["mic"]][["comp.res"]])[3] = "T.Stats"
    d.list[["mic"]][["comp.res"]] = list(OTU=d.list[["mic"]][["comp.res"]])
    d.list[["mic"]][["enrich.nms"]] = list(OTU=rownames(current.proc$mic$res_deAnal))
    
  }else{
    d.list[["mic"]][["comp.res"]] = lapply(phyloseq_objs$res_deAnal, function(x){names(x)[1] ="T.Stats"; return(x[,c(3,4,1), drop=FALSE])})
    d.list[["mic"]][["enrich.nms"]] = lapply(phyloseq_objs$res_deAnal ,function(x) rownames(x))
  }
  d.list[["mic"]][["meta"]] = data.frame(mbSetObj$dataSet$sample_data)
  
  d.list[["met"]] = list()
  d.list[["met"]][["data.proc"]] = if(!is.null(current.proc$met$data.norm)) current.proc$met$data.norm else current.proc$met$data.proc
  met_de <- current.proc$met$res_deAnal
  if (is.null(met_de) || !is.data.frame(met_de)) met_de <- data.frame(matrix(0, nrow=1, ncol=3))
  ncols <- min(3, ncol(met_de))
  d.list[["met"]][["comp.res"]] = met_de[, 1:ncols, drop=FALSE]
  d.list[["met"]][["enrich.nms"]] = rownames(current.proc$met$res_deAnal)
  d.list[["met"]][["meta"]] = data.frame(mbSetObj$dataSet$sample_data)
  
  newmeta <- rbind(  d.list[["mic"]][["meta"]],d.list[["met"]][["meta"]])
  comp.res1 = lapply( d.list[["mic"]][["comp.res"]],function(x) rbind(x, d.list[["met"]][["comp.res"]]) )
  
  enrich.nms1 = lapply(d.list[["mic"]][["enrich.nms"]],function(x) c(x, d.list[["met"]][["enrich.nms"]]) )
  comp.res.inx1 =lapply(d.list[["mic"]][["comp.res"]],function(x) c(rep(1,nrow(x)), rep(2,nrow( d.list[["met"]][["comp.res"]])) ));
  featureNms <- lapply( d.list[["mic"]][["data.proc"]] ,function(x) c(rownames(x), rownames(d.list[["met"]][["data.proc"]]) ));
  
  combined.res$comp.res = comp.res1
  combined.res$enrich_ids = enrich.nms1
  combined.res$comp.res.inx = comp.res.inx1
  combined.res$meta = newmeta
  if(reductionOpt == "diablo"){
    meta.df <- NULL
    if(!is.null(current.proc$meta_para$sample_data) && is.data.frame(current.proc$meta_para$sample_data)){
      meta.df <- current.proc$meta_para$sample_data
    } else if(!is.null(mbSetObj$dataSet$sample_data) && is.data.frame(mbSetObj$dataSet$sample_data)){
      meta.df <- mbSetObj$dataSet$sample_data
    }
    if(is.null(meta.df) || ncol(meta.df) < 1){
      stop("[DIABLO] Sample metadata is missing or empty.")
    }

    resolve_meta_name <- function(nm, choices){
      if(length(nm) == 0 || is.null(nm) || is.na(nm)) return(NA_character_)
      nm <- as.character(nm)[1]
      if(nm %in% choices) return(nm)
      choices_trim <- trimws(as.character(choices))
      hit <- which(choices_trim == trimws(nm))
      if(length(hit) > 0) return(choices[hit[1]])
      hit <- which(tolower(choices_trim) == tolower(trimws(nm)))
      if(length(hit) > 0) return(choices[hit[1]])
      return(NA_character_)
    }

    analysisVar.orig <- as.character(analysisVar)[1]
    analysisVar <- resolve_meta_name(analysisVar.orig, colnames(meta.df))
    if(is.na(analysisVar) && !is.null(mbSetObj$dataSet$sample_data) && is.data.frame(mbSetObj$dataSet$sample_data)){
      meta.df2 <- mbSetObj$dataSet$sample_data
      alt <- resolve_meta_name(analysisVar.orig, colnames(meta.df2))
      if(!is.na(alt)){
        meta.df <- meta.df2
        analysisVar <- alt
      }
    }
    if(is.na(analysisVar)){
      stop(paste0("[DIABLO] Analysis variable '", analysisVar.orig, "' is not available in sample metadata."))
    }
    if(analysisVar != analysisVar.orig){
      message("[DIABLO] Using metadata variable '", analysisVar, "' (resolved from '", analysisVar.orig, "').")
    }

    meta.type.nms <- names(mbSetObj$dataSet$meta.types)
    meta.type.key <- resolve_meta_name(analysisVar, meta.type.nms)
    diablo.meta.type <- if(!is.na(meta.type.key)) mbSetObj$dataSet$meta.types[meta.type.key] else NA
    diablo.meta.type <- if(length(diablo.meta.type) > 0) as.character(diablo.meta.type[[1]]) else NA_character_
    if(is.na(diablo.meta.type) || !(diablo.meta.type %in% c("disc", "cont"))){
      meta.var.raw <- meta.df[, analysisVar, drop = TRUE]
      if(is.numeric(meta.var.raw)){
        diablo.meta.type <- "cont"
      } else {
        meta.num <- suppressWarnings(as.numeric(as.character(meta.var.raw)))
        if(sum(!is.na(meta.num)) == length(meta.var.raw) && length(unique(meta.num)) > 5){
          diablo.meta.type <- "cont"
        } else {
          diablo.meta.type <- "disc"
        }
      }
      message("[DIABLO] Metadata type for '", analysisVar, "' was missing/invalid; inferred as '", diablo.meta.type, "'.")
    }

    diabloPar <- as.numeric(diabloPar); #default diabloPar was 0.2
    diablo.res <- list()
    dats <- vector("list",length=length(d.list$mic$data.proc))
    names(dats) <- names(d.list$mic$data.proc)

    res <- pos.xyz <- pos.xyz2 <- misc <- loading.pos.xyz <- loadingNames <- vector("list",length=length(d.list$mic$data.proc))
    names(res) <- names(pos.xyz) <- names(pos.xyz2) <- names(misc) <-
      names(loading.pos.xyz) <- names(loadingNames) <- names(d.list$mic$data.proc)

    for(l in 1:length(dats)){

      dats[[l]][["mic"]] <- d.list$mic$data.proc[[l]]
      dats[[l]][["met"]] <- d.list$met$data.proc
      dats[[l]]= lapply(dats[[l]], function(x){
        x <- data.frame(x, stringsAsFactors = FALSE);
        x <- t(x);
      })

      # Prepare Y and remove invalid samples for this analysis variable
      if(isTRUE(diablo.meta.type == "disc")){
        Y_factor <- as.factor(meta.df[, analysisVar, drop = TRUE])
        keep.inx <- !is.na(Y_factor) & as.character(Y_factor) != ""
        Y_factor <- droplevels(Y_factor[keep.inx])
      } else {
        meta.var <- meta.df[, analysisVar, drop = TRUE]
        Y_numeric <- suppressWarnings(as.numeric(as.character(meta.var)))
        keep.inx <- is.finite(Y_numeric)
        Y_numeric <- Y_numeric[keep.inx]
      }

      dats[[l]] <- lapply(dats[[l]], function(x) x[keep.inx, , drop = FALSE])

      n_samples <- nrow(dats[[l]]$mic)
      if(n_samples < 2){
        stop(paste0("[DIABLO] Not enough valid samples for '", analysisVar, "' after filtering missing/non-numeric values."))
      }

      max_ncomp_by_samples <- n_samples - 1  # Maximum components = n_samples - 1

      # Calculate appropriate ncomp considering sample size and class size
      if(isTRUE(diablo.meta.type == "disc")){
        if(length(unique(Y_factor)) < 2){
          stop(paste0("[DIABLO] Analysis variable '", analysisVar, "' needs at least 2 groups after filtering."))
        }
        min_class_size <- min(table(Y_factor))
        max_ncomp_by_classes <- min_class_size - 1  # Each class needs enough samples
        max_ncomp <- max(1, min(max_ncomp_by_samples, max_ncomp_by_classes, ncol(dats[[l]]$mic), dimn))
        Y <- Y_factor
      } else {
        max_ncomp <- max(1, min(max_ncomp_by_samples, ncol(dats[[l]]$mic), dimn))
        Y <- matrix(Y_numeric, ncol = 1)
        rownames(Y) <- rownames(dats[[l]]$mic)
      }

      ncomp <- max_ncomp

      design = matrix(diabloPar, ncol = length(dats[[l]]), nrow = length(dats[[l]]),
                      dimnames = list(names(dats[[l]]), names(dats[[l]])));
      diag(design) = 0;

      # mixOmics quarantined — run in RSclient, matching original callr pattern exactly
      work_dir <- getwd()
      # Image names set by Java via diablo.img.names global
      img.names <- if (exists("diablo.img.names")) diablo.img.names else list(
        pca="diagnostic_pca_diablo_0.png", loading="diagnostic_loading_0.png",
        diag="diagnostic_components_diablo_0.png", circos="diagnostic_circos_diablo_0.png")

      qs::qsave(list(dats = dats[[l]], Y = Y, ncomp = ncomp, design = design,
                      meta_type = diablo.meta.type,
                      pca_img = file.path(work_dir, img.names$pca),
                      loading_img = file.path(work_dir, img.names$loading),
                      diag_img = file.path(work_dir, img.names$diag),
                      circos_img = file.path(work_dir, img.names$circos)), "diablo_input.qs")
      diablo_result <- run_func_via_rsclient(
        func = function(work_dir) {
          setwd(work_dir)
          suppressPackageStartupMessages({
            require(mixOmics); require(Cairo); require(grid)
            require(gridExtra); require(gridGraphics); require(cowplot)
          })
          input <- qs::qread("diablo_input.qs")
          message("[DIABLO] Starting model fitting with ", input$ncomp, " components...")
          start_time <- Sys.time()
          if (input$meta_type == "disc") {
            model <- mixOmics::block.splsda(X = input$dats, Y = input$Y,
                                            ncomp = input$ncomp, design = input$design, near.zero.var = TRUE)
          } else {
            model <- mixOmics::block.spls(X = input$dats, Y = input$Y,
                                          ncomp = input$ncomp, design = input$design, mode = "regression", near.zero.var = TRUE)
          }
          message("[DIABLO] Model fitting completed in ", round(difftime(Sys.time(), start_time, units="secs"), 1), " seconds")
          # Extract plain data
          pos.xyz <- model$variates[[1]]
          pos.xyz2 <- model$variates[[2]]
          loadings <- lapply(model$loadings, function(x) {
            pos <- as.data.frame(x); rn <- rownames(x)
            for (i in 1:ncol(pos)) {
              col_min <- min(pos[, i], na.rm = TRUE); col_max <- max(pos[, i], na.rm = TRUE)
              if (col_max - col_min > 0) pos[, i] <- (pos[, i] - col_min) / (col_max - col_min) - 0.5
            }
            rownames(pos) <- rn; pos
          })
          loading.pos.xyz <- rbind(loadings[[1]], loadings[[2]])
          var.vec <- if ("prop_expl_var" %in% names(model)) model$prop_expl_var
                     else if ("explained_variance" %in% names(model)) model$explained_variance
                     else list(mic = 0, met = 0)

          # Save diablo.res.qs for scatter viewer
          diablo.res <- list(dim.res = setNames(list(model), names(input$dats)[1]))
          qs::qsave(diablo.res, "diablo.res.qs")
          # Save raw model separately for perf() — avoids double qs corruption
          qs::qsave(model, "diablo_model.qs")

          # Generate ALL diagnostic plots while model is in memory
          # Helper to safely generate a plot
          .safe_plot <- function(img_file, width, height, plot_fn) {
            tryCatch({
              Cairo(file = img_file, width = width, height = height, type = "png", bg = "white", unit = "in", dpi = 96)
              plot_fn()
              dev.off()
            }, error = function(e) {
              tryCatch(dev.off(), error = function(e2) NULL)
              message("[DIABLO plot] ", e$message)
            })
          }

          # 1. plotDiablo (PCA diagnostic) — only plot components that exist
          message("[DIABLO] Generating PCA diagnostic plots...")
          .safe_plot(input$pca_img, 8, 15, function() {
            ncomp_max <- min(3, model$ncomp[1])
            fig.list <- list()
            for (nc in 1:ncomp_max) {
              local_nc <- nc
              fig.list[[nc]] <- tryCatch(
                as_grob(function() { plotDiablo(model, ncomp = local_nc) }),
                error = function(e) { message("[plotDiablo ncomp=", local_nc, "] ", e$message); NULL }
              )
            }
            fig.list <- fig.list[!sapply(fig.list, is.null)]
            if (length(fig.list) > 0) {
              grid.arrange(grobs = fig.list, nrow = length(fig.list))
            } else {
              plot.new(); text(0.5, 0.5, "plotDiablo failed for all components", cex = 1.2)
            }
          })

          # 2. Loading plots — custom robust rendering (avoids mixOmics::plotLoadings margin issues)
          message("[DIABLO] Generating loading plots...")
          .safe_plot(input$loading_img, 10, 18, function() {  # Reduced from 13×24 to 10×18 for faster rendering
            ncomp.raw <- suppressWarnings(as.integer(model$ncomp))
            ncomp.max <- if(length(ncomp.raw) == 0) NA_integer_ else suppressWarnings(min(ncomp.raw, na.rm = TRUE))
            if(!is.finite(ncomp.max) || is.na(ncomp.max) || ncomp.max < 1){
              ncomp.max <- suppressWarnings(min(vapply(model$loadings, function(v) ncol(as.matrix(v)), integer(1)), na.rm = TRUE))
            }
            if(!is.finite(ncomp.max) || is.na(ncomp.max) || ncomp.max < 1){
              ncomp.max <- 1L
            }
            ncomp.plot <- min(3L, as.integer(ncomp.max))

            op <- par(no.readonly = TRUE)
            on.exit(par(op), add = TRUE)
            par(mfrow = c(ncomp.plot, 1), mar = c(5, 4, 3, 1))

            for (nc in seq_len(ncomp.plot)) {
              tryCatch({
                top.rows <- lapply(names(model$loadings), function(block.nm){
                  mat <- as.matrix(model$loadings[[block.nm]])
                  if(ncol(mat) < nc){
                    return(NULL)
                  }
                  vals <- mat[, nc, drop = TRUE]
                  keep <- is.finite(vals)
                  vals <- vals[keep]
                  if(length(vals) == 0){
                    return(NULL)
                  }
                  ord <- order(abs(vals), decreasing = TRUE)
                  ord <- ord[seq_len(min(10, length(ord)))]
                  data.frame(
                    feature = names(vals)[ord],
                    loading = as.numeric(vals[ord]),
                    block = block.nm,
                    stringsAsFactors = FALSE
                  )
                })

                top.df <- do.call(rbind, top.rows)
                if(is.null(top.df) || nrow(top.df) == 0){
                  stop(paste0("No finite loadings available for component ", nc, "."))
                }

                top.df$feature <- paste0(top.df$block, ": ", top.df$feature)
                top.df <- top.df[order(abs(top.df$loading), decreasing = TRUE), , drop = FALSE]
                top.df$feature <- factor(top.df$feature, levels = rev(unique(top.df$feature)))

                p <- ggplot2::ggplot(top.df, ggplot2::aes(x = feature, y = loading, fill = block)) +
                  ggplot2::geom_col(width = 0.7) +
                  ggplot2::coord_flip() +
                  ggplot2::labs(
                    title = paste("DIABLO Loadings - Component", nc),
                    x = NULL,
                    y = "Loading"
                  ) +
                  ggplot2::theme_bw(base_size = 10) +
                  ggplot2::theme(
                    legend.position = "bottom",
                    plot.title = ggplot2::element_text(hjust = 0.5)
                  )
                print(p)
              }, error = function(e) {
                message("[custom loading comp=", nc, "] ", e$message)
                plot.new()
                text(0.5, 0.5, paste0("Component ", nc, " failed:\n", e$message), cex = 1)
              })
            }
          })

          # 3. perf — skipped here, runs in separate PlotDiagnostic call (clean environment)
          diablo.comp <- NA

          # 4. circos JSON for interactive chord viewer
          # Can be disabled by setting diablo.skip.circos <<- TRUE from Java
          skip_circos <- if (exists("diablo.skip.circos")) diablo.skip.circos else FALSE
          if (!skip_circos) {
          message("[DIABLO] Computing correlations for circos plot...")
          circos_start <- Sys.time()
          tryCatch({
            # Compute cross-block correlations from projected data
            block_names <- names(model$X)
            message("[DIABLO circos] Filtering features with non-zero variance...")
            X_proj <- lapply(model$X, function(x) {
              x[, which(apply(x, 2, var) > 0), drop = FALSE]
            })

            # PRE-FILTER: Limit features to top N by loading weights to speed up correlation
            max_features <- if (exists("diablo.circos.max.features")) diablo.circos.max.features else 1000
            X_proj_filtered <- lapply(names(X_proj), function(block_name) {
              x <- X_proj[[block_name]]
              if (ncol(x) > max_features) {
                # Select top features by sum of absolute loadings across components
                loading_mat <- as.matrix(model$loadings[[block_name]])
                loading_scores <- rowSums(abs(loading_mat))
                top_idx <- order(loading_scores, decreasing = TRUE)[1:max_features]
                message("[DIABLO circos] Pre-filtering ", block_name, ": ", ncol(x), " -> ", max_features, " features")
                x[, top_idx, drop = FALSE]
              } else {
                x
              }
            })
            names(X_proj_filtered) <- names(X_proj)
            X_proj <- X_proj_filtered

            message("[DIABLO circos] Computing correlation matrix (", ncol(X_proj[[1]]), " x ", ncol(X_proj[[2]]), " = ", ncol(X_proj[[1]]) * ncol(X_proj[[2]]), " correlations)...")
            cor_cross <- cor(X_proj[[1]], X_proj[[2]])
            message("[DIABLO circos] Correlation matrix computed in ", round(difftime(Sys.time(), circos_start, units="secs"), 1), " seconds")
            # Use cutoff from Java (default 0.5)
            cutoff <- if (exists("diablo.circos.cutoff")) diablo.circos.cutoff else 0.5
            max_edges <- if (exists("diablo.circos.max.edges")) diablo.circos.max.edges else 100

            message("[DIABLO circos] Filtering edges with cutoff ", cutoff, "...")
            sig_idx <- which(abs(cor_cross) > cutoff, arr.ind = TRUE)
            message("[DIABLO circos] Found ", nrow(sig_idx), " edges above cutoff")

            # ALWAYS enforce maximum edges limit to prevent freeze
            if (nrow(sig_idx) > max_edges) {
              message("[DIABLO circos] Too many edges (", nrow(sig_idx), "), limiting to top ", max_edges, " by correlation...")
              # Get correlation values for edges above cutoff
              edge_cors <- abs(cor_cross[sig_idx])
              # Sort and keep only top max_edges
              top_indices <- order(edge_cors, decreasing = TRUE)[1:max_edges]
              sig_idx <- sig_idx[top_indices, , drop = FALSE]
              message("[DIABLO circos] Reduced to ", nrow(sig_idx), " edges")
            }

            if (nrow(sig_idx) > 0) {
              edges <- lapply(1:nrow(sig_idx), function(i) {
                list(source = rownames(cor_cross)[sig_idx[i,1]],
                     target = colnames(cor_cross)[sig_idx[i,2]],
                     corr = round(cor_cross[sig_idx[i,1], sig_idx[i,2]], 4),
                     type1 = block_names[1], type2 = block_names[2],
                     label1 = rownames(cor_cross)[sig_idx[i,1]],
                     label2 = colnames(cor_cross)[sig_idx[i,2]])
              })
            } else {
              # If no edges above cutoff, use top N by absolute correlation
              top_n <- min(max_edges, length(cor_cross))
              message("[DIABLO circos] No edges above cutoff, selecting top ", top_n, " edges by correlation...")
              top_idx <- order(abs(cor_cross), decreasing = TRUE)[1:top_n]
              message("[DIABLO circos] Top edges selected")
              sig_idx <- arrayInd(top_idx, dim(cor_cross))
              edges <- lapply(1:nrow(sig_idx), function(i) {
                list(source = rownames(cor_cross)[sig_idx[i,1]],
                     target = colnames(cor_cross)[sig_idx[i,2]],
                     corr = round(cor_cross[sig_idx[i,1], sig_idx[i,2]], 4),
                     type1 = block_names[1], type2 = block_names[2],
                     label1 = rownames(cor_cross)[sig_idx[i,1]],
                     label2 = colnames(cor_cross)[sig_idx[i,2]])
              })
            }
            message("[DIABLO circos] Creating edge list with ", length(edges), " edges...")
            circos_json <- list(DIABLO = edges)
            message("[DIABLO circos] Writing JSON to file...")
            jsonlite::write_json(circos_json, "diablo_circos.json", auto_unbox = TRUE, pretty = FALSE)
            message("[DIABLO circos] JSON completed in ", round(difftime(Sys.time(), circos_start, units="secs"), 1), " seconds")
          }, error = function(e) message("[DIABLO circosJSON] ", e$message))

          # Also generate static circos PNG as fallback (skip if too many features to avoid freeze)
          n_features <- sum(sapply(model$X, ncol))
          skip_png <- n_features > 1000  # Skip PNG if more than 1000 total features

          if (skip_png) {
            message("[DIABLO] Skipping circos PNG rendering (", n_features, " features is too many, would freeze)")
            message("[DIABLO] Creating placeholder PNG instead...")
            # Create a simple placeholder PNG
            .safe_plot(input$circos_img, 8, 8, function() {
              plot.new()
              text(0.5, 0.5, paste0("Circos plot not rendered\n(", n_features, " features)\n\nUse interactive viewer instead"),
                   cex = 1.5, col = "gray40")
            })
          } else {
            message("[DIABLO] Rendering circos PNG (", n_features, " features, max ", max_edges, " edges)...")
            png_start <- Sys.time()
            .safe_plot(input$circos_img, 8, 8, function() {  # Reduced from 10×10 to 8×8 for faster rendering
              cutoff_val <- if (exists("diablo.circos.cutoff")) diablo.circos.cutoff else 0.5
              size_val <- if (exists("diablo.circos.feature.size")) diablo.circos.feature.size else 0.25
              # Use higher cutoff for PNG if many features to prevent freeze
              png_cutoff <- if (n_features > 500) max(0.7, cutoff_val) else cutoff_val
              if (png_cutoff > cutoff_val) {
                message("[DIABLO circos PNG] Using higher cutoff (", png_cutoff, ") for static plot")
              }
              circosPlot(model, cutoff = png_cutoff, size.variables = size_val, showIntraLinks = FALSE)
            })
            message("[DIABLO] Circos PNG completed in ", round(difftime(Sys.time(), png_start, units="secs"), 1), " seconds")
          }
          } else {
            message("[DIABLO] Circos plot generation skipped (diablo.skip.circos = TRUE)")
            # Create empty JSON so viewer doesn't break
            jsonlite::write_json(list(DIABLO = list()), "diablo_circos.json", auto_unbox = TRUE, pretty = FALSE)
          }

          gc(verbose = FALSE, full = TRUE)
          message("[DIABLO] All plots completed. Total time: ", round(difftime(Sys.time(), start_time, units="secs"), 1), " seconds")
          list(pos.xyz = pos.xyz, pos.xyz2 = pos.xyz2,
               loading.pos.xyz = loading.pos.xyz, loadings = loadings,
               var.vec = var.vec, diablo.comp = diablo.comp)
        },
        args = list(work_dir = work_dir),
        timeout_sec = 600
      )
      unlink("diablo_input.qs")
      if (!is.null(diablo_result$diablo.comp) && !is.na(diablo_result$diablo.comp)) {
        diablo.comp <<- diablo_result$diablo.comp
      }

      # Model stays in diablo.res.qs (single qs, never read into master)
      # Read back for dim.res field needed by diablo.res structure
      res[[l]] <- qs::qread("diablo.res.qs")$dim.res[[1]]
      pos.xyz[[l]] <- diablo_result$pos.xyz
      pos.xyz2[[l]] <- diablo_result$pos.xyz2
      loading.pos.xyz[[l]] <- diablo_result$loading.pos.xyz
      var.vec <- diablo_result$var.vec

      loadingNames[[l]]=rownames(loading.pos.xyz[[l]])
      names = lapply(pos.xyz, function(x) rownames(x))
      newmeta = current.proc$meta_para$sample_data

      misc[[l]]$pct2[["microbiome"]] = unname(signif(as.numeric(var.vec[['mic']],4)))*100
      misc[[l]]$pct2[["metabolomics"]] = unname(signif(as.numeric(var.vec[['met']],4)))*100

    }

  }else if(reductionOpt == "procrustes"){
    # vegan quarantined — matching original callr pattern
    mic_data_list <- d.list$mic$data.proc
    met_data <- d.list$met$data.proc
    work_dir <- getwd()
    qs::qsave(list(mic = mic_data_list, met = met_data), "procrustes_input.qs")

    proc_result <- run_func_via_rsclient(
      func = function(work_dir) {
        setwd(work_dir)
        require(vegan)
        input <- qs::qread("procrustes_input.qs")
        ndat1 <- lapply(input$mic, function(x) vegan::decostand(t(x), method = "standardize"))
        pca.dat1 <- lapply(ndat1, function(x) vegan::rda(x))
        ndat2 <- vegan::decostand(t(input$met), method = "standardize")
        pca.dat2 <- vegan::rda(ndat2)
        proc_res <- lapply(pca.dat1, function(x) vegan::procrustes(x, pca.dat2, choices = c(1,2,3), symmetric = TRUE, scale = TRUE))
        prot_res <- lapply(pca.dat1, function(x) vegan::protest(X = x, Y = pca.dat2, scores = "sites", permutations = 999))
        misc <- lapply(prot_res, function(x) list(`Sum of Squares` = x$ss, Significance = x$signif, Correlation = x$scale))
        pos.xyz <- lapply(proc_res, function(x) rbind(x$X, x$Yrot))
        gc(verbose = FALSE, full = TRUE)
        list(proc_res = proc_res, prot_res = prot_res, misc = misc, pos.xyz = pos.xyz)
      },
      args = list(work_dir = work_dir),
      timeout_sec = 300
    )
    unlink("procrustes_input.qs")

    misc <- proc_result$misc
    pos.xyz <- proc_result$pos.xyz

    names = lapply(pos.xyz, function(x) make.unique(as.character(rownames(x))))
    newmeta$omics[c(1:(length(names[[1]])/2))] = "microbiome"
    newmeta$omics[c(((length(names[[1]])/2)+1):(length(names[[1]])))] = "metabolomics"
    dim.res <- vector("list", length=length(proc_result$proc_res))
    names(dim.res) <- names(misc) <- names(proc_result$proc_res)
    for(i in 1:length(dim.res)){
      dim.res[[i]] <- list(proc_result$proc_res[[i]], proc_result$prot_res[[i]])
    }

    procrustes.res <- list(misc=misc, dim.res=dim.res)
    procrustes.res$pos.xyz = lapply(pos.xyz,function(x) unitAutoScale(x))
    procrustes.res$newmeta = newmeta
    combined.res$meta = newmeta
    shadow_save(combined.res,"combined.res.qs")
    shadow_save(procrustes.res,"procrustes.res.qs")

  } else if(reductionOpt == "mofa"){
    # Save input for external MOFA script (avoids HDF5Array conflicts in Rserve)
    mic_mat <- as.matrix(d.list$mic$data.proc[[1]])
    met_mat <- as.matrix(d.list$met$data.proc)
    tax_name <- names(d.list$mic$data.proc)[1]
    met_features <- rownames(d.list$met$data.proc)

    mofa.input <- list(
      data.list = list(mic = mic_mat, met = met_mat),
      tax_name = tax_name,
      met_features = met_features
    )
    shadow_save(mofa.input, "mofa_input.qs")
    shadow_save(combined.res, "combined.res.qs")

    # Return 2 to signal Java to run _perform_mofa.R externally
    return(2)
  }

    pos.xyz <- lapply(pos.xyz,function(x) as.data.frame(x)[,c(1:3)]);
    pos.xyz <- lapply(pos.xyz,function(x) unitAutoScale(x));
    pos.xyz <- lapply(pos.xyz, "rownames<-", names[[1]]);

  if(reductionOpt %in% c("diablo")){
    names2 <- lapply(pos.xyz2, function(x) rownames(x))
    pos.xyz <- lapply(pos.xyz,function(x) as.data.frame(x)[,c(1:3)]);
    pos.xyz2 <- lapply(pos.xyz2, function(x) unitAutoScale(x));
    pos.xyz2 <- lapply(pos.xyz2, "rownames<-", names2[[1]]);
  }

  
  
  
  if(reductionOpt != "procrustes"){
    hit.inx <- mapply(function(x,y){
      match(x,y)
    },loadingNames, combined.res$enrich_ids);
    loadingSymbols <- mapply(function(x,y){
      y[x]
    },hit.inx, combined.res$enrich_ids);
    if(micDataType=="ko"){
      loadingSymbols=list(OTU=loadingSymbols)
    } else if(!is.list(loadingSymbols)){
      loadingSymbols=setNames(list(loadingSymbols), names(loadingNames))
    }
  }
  if(reductionOpt == "mofa"){
    loading.pos.xyz = lapply(loading.pos.xyz, as.data.frame)
    loading.pos.xyz <- lapply(loading.pos.xyz, unitAutoScale)
    mofa.res$pos.xyz <- pos.xyz
    mofa.res$loading.pos.xyz <- loading.pos.xyz
    mofa.res$loadingNames <- loadingNames
    mofa.res$loading.enrich <- loadingSymbols
    shadow_save(mofa.res, "mofa.res.qs")
  }
  if(reductionOpt == "diablo"){
    loading.pos.xyz = lapply(loading.pos.xyz,as.data.frame)
    loading.pos.xyz <- lapply(loading.pos.xyz,unitAutoScale);

    diablo.res$dim.res <- res
    diablo.res$pos.xyz <- pos.xyz
    diablo.res$pos.xyz2 <- pos.xyz2
    diablo.res$misc <- misc
    diablo.res$loading.pos.xyz <- loading.pos.xyz
    diablo.res$loadingNames <- loadingNames
    diablo.res$loading.enrich = loadingSymbols
    diablo.res$newmeta=newmeta
    shadow_save(combined.res,"combined.res.qs")
    shadow_save(diablo.res,"diablo.res.qs")
    
  }
  reductionOptGlobal <<- reductionOpt;
  
  mbSetObj$analSet$dimrdc<-reductionOpt;
  
  .set.mbSetObj(mbSetObj);

  return(1)
}
