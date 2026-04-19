
#'Plot PCA plot for meta-analysis samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'

PlotMetaPCA <- function(imgNm, dpi, format, factor="NA"){
  require('ggplot2');
  require('ggsci');

  microbiome.meta <- ov_qs_read("microbiome_meta.qs");
  x <- microbiome.meta[["data"]];
  dpi <- as.numeric(dpi);
  imgNm <- paste(imgNm, ".", format, sep="");

  Conditions <- factor(microbiome.meta$cls.lbl);
  Datasets <- factor(microbiome.meta$data.lbl);
  dataset_names <- levels(Datasets);

  .make_pca <- function(dat, conditions, title){
    dat <- na.omit(dat)
    vars <- apply(dat, 1, var, na.rm=TRUE)
    dat <- dat[vars > 0, , drop=FALSE]
    if(nrow(dat) < 3 || ncol(dat) < 3) return(ggplot() + ggtitle(title) + theme_void())
    pca <- prcomp(t(dat), center=TRUE, scale.=TRUE)
    imp <- summary(pca)$importance
    xl <- paste0("PC1 (", 100*round(imp[2,1], 3), "%)")
    yl <- paste0("PC2 (", 100*round(imp[2,2], 3), "%)")
    df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Conditions=conditions)
    ggplot(df, aes(x=PC1, y=PC2, color=Conditions)) +
      geom_point(size=3, alpha=0.7) +
      scale_color_npg() +
      xlab(xl) + ylab(yl) + ggtitle(title) +
      theme_bw() +
      theme(plot.title=element_text(size=11, face="bold"))
  }

  # Two columns: Before (individual) | After (individual from merged)
  # One row per dataset + one row for combined
  fig.list <- list()
  idx <- 1

  if(file.exists("microbiome_meta_prenorm.qs")){
    meta.before <- ov_qs_read("microbiome_meta_prenorm.qs")
  } else {
    meta.before <- microbiome.meta  # fallback: before = after
  }

  valid_datasets <- c()
  for(ds in dataset_names){
    ds_inx <- which(Datasets == ds)
    if(length(ds_inx) == 0) next  # skip if no samples

    valid_datasets <- c(valid_datasets, ds)
    # Before: from pre-normalization/pre-batch snapshot
    ds_before <- meta.before$data[, ds_inx, drop=FALSE]
    ds_cond <- Conditions[ds_inx]
    fig.list[[idx]] <- .make_pca(ds_before, ds_cond, paste0(ds, " (Before)"))
    idx <- idx + 1

    # After: from current (post-normalization/batch-corrected)
    ds_after <- x[, ds_inx, drop=FALSE]
    fig.list[[idx]] <- .make_pca(ds_after, ds_cond, paste0(ds, " (After)"))
    idx <- idx + 1
  }

  # Combined row
  fig.list[[idx]] <- .make_pca(meta.before$data, Conditions, "Combined (Before)")
  idx <- idx + 1
  fig.list[[idx]] <- .make_pca(x, Conditions, "Combined (After)")

  nrows <- length(valid_datasets) + 1
  h <- 4 * nrows

  # ggpubr quarantined — ggarrange in subprocess
  rsclient_isolated_exec(
    func_body = function(d) {
      require(ggpubr)
      Cairo::Cairo(file = d$imgNm, width = 10, height = d$h, type = d$format, bg = "white", unit = "in", dpi = 96)
      p <- ggpubr::ggarrange(plotlist = d$fig_list, ncol = 2, nrow = d$nrows, common.legend = TRUE, legend = "bottom")
      print(p)
      dev.off()
    },
    input_data = list(imgNm = imgNm, h = h, format = format, fig_list = fig.list, nrows = nrows),
    packages = c("ggpubr", "Cairo", "qs"), timeout = 120, output_type = "qs"
  )
}

#'Plot density plot of datasets in meta-analysis module
#'@description 
#'@param imgNm name of the image to output
#'@param format File format of output image (default to png)
#'@param dpi dots per inch resolution of image (default 72)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotMetaDensity <- function(imgNm, dpi=default.dpi, format="png", factor=""){
  require("ggplot2")

  microbiome.meta <- ov_qs_read("microbiome_meta.qs");
  x <- microbiome.meta$data;
  imgNm <- paste(imgNm, ".", format, sep="");
  dpi <- as.numeric(dpi);

  Conditions <- factor(microbiome.meta$cls.lbl);
  Datasets <- factor(microbiome.meta$data.lbl);
  dataset_names <- levels(Datasets);

  .make_density <- function(dat, conditions, title){
    if(ncol(dat) < 2) return(ggplot() + ggtitle(title) + theme_void())
    df <- data.frame(dat, stringsAsFactors=FALSE, check.names=FALSE)
    df <- stack(df)
    if(nrow(df) == 0 || all(is.na(df$values))) return(ggplot() + ggtitle(title) + theme_void())
    df$ind <- as.character(df$ind)
    conv <- data.frame(ind=colnames(dat), class=conditions, stringsAsFactors=FALSE)
    # Normalize both sides for matching
    df$ind <- gsub("-", ".", df$ind)
    conv$ind <- gsub("-", ".", conv$ind)
    df1 <- merge(df, conv, by="ind")
    if(nrow(df1) == 0) return(ggplot() + ggtitle(title) + theme_void())
    ggplot(df1, aes(x=values)) +
      geom_line(aes(color=class, group=ind), stat="density", alpha=0.3) +
      geom_line(aes(color=class), stat="density", alpha=0.6, linewidth=1) +
      ggtitle(title) +
      theme_bw() +
      theme(plot.title=element_text(size=11, face="bold"))
  }

  if(file.exists("microbiome_meta_prenorm.qs")){
    meta.before <- ov_qs_read("microbiome_meta_prenorm.qs")
  } else {
    meta.before <- microbiome.meta  # fallback: before = after
  }

  fig.list <- list()
  idx <- 1
  valid_datasets <- c()
  for(ds in dataset_names){
    ds_inx <- which(Datasets == ds)
    if(length(ds_inx) == 0) next  # skip if no samples

    valid_datasets <- c(valid_datasets, ds)
    ds_cond <- Conditions[ds_inx]

    # Before
    ds_before <- meta.before$data[, ds_inx, drop=FALSE]
    fig.list[[idx]] <- .make_density(ds_before, ds_cond, paste0(ds, " (Before)"))
    idx <- idx + 1

    # After
    ds_after <- x[, ds_inx, drop=FALSE]
    fig.list[[idx]] <- .make_density(ds_after, ds_cond, paste0(ds, " (After)"))
    idx <- idx + 1
  }

  # Combined row
  fig.list[[idx]] <- .make_density(meta.before$data, Conditions, "Combined (Before)")
  idx <- idx + 1
  fig.list[[idx]] <- .make_density(x, Conditions, "Combined (After)")

  nrows <- length(valid_datasets) + 1
  h <- 3.5 * nrows

  rsclient_isolated_exec(
    func_body = function(d) {
      require(ggpubr)
      Cairo::Cairo(file = d$imgNm, width = 10, height = d$h, type = d$format, bg = "white", dpi = 96, unit = "in")
      p <- ggpubr::ggarrange(plotlist = d$fig_list, ncol = 2, nrow = d$nrows, common.legend = TRUE, legend = "bottom")
      print(p)
      dev.off()
    },
    input_data = list(imgNm = imgNm, h = h, format = format, fig_list = fig.list, nrows = nrows),
    packages = c("ggpubr", "Cairo", "qs"), timeout = 120, output_type = "qs"
  )
}


#'Perform batch correction of meta-analysis datasets
#'@description 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
ApplyMetaAutoScale <- function(apply="true"){
  microbiome.meta <- ov_qs_read("microbiome_meta.qs");
  if(apply == "true"){
    microbiome.meta$data <- t(scale(t(microbiome.meta$data), center=TRUE, scale=TRUE))
    # Replace NaN from zero-variance features with 0
    microbiome.meta$data[is.nan(microbiome.meta$data)] <- 0
  }
  shadow_save(microbiome.meta, "microbiome_meta.qs")
  return(1)
}

PerformBatchCorrection <- function(){
    # MMUPHin batch correction — runs in isolated subprocess for memory safety
    bridge_in <- paste0(tempdir(), "/bridge_", paste0(sample(letters,6,replace=TRUE), collapse=""), "_in.qs")
    bridge_out <- sub("_in.qs", "_out.qs", bridge_in)
    ov_qs_save(list(placeholder = TRUE), bridge_in, preset = "fast")
    on.exit(unlink(c(bridge_in, bridge_out)), add = TRUE)

    run_func_via_rsclient(
      func = function(wd, bridge_in, bridge_out) {
        setwd(wd)
        require(MMUPHin)
        require(phyloseq)
        phyobj <- ov_qs_read("merged.data.qs")
        sam.data <- as.data.frame(as.matrix(sample_data(phyobj)))

        cov <- colnames(sam.data)[1]
        otu.tbl <- otu_table(phyobj)

        microbiome.meta <- ov_qs_read("microbiome_meta.qs")

        fit_adjust_batch <- adjust_batch(feature_abd = otu.tbl,
                                         batch = "dataset",
                                         covariates = cov,
                                         data = sam.data,
                                         control = list(verbose = FALSE))

        adj.otu.tbl <- fit_adjust_batch$feature_abd_adj
        phyobj@otu_table <- otu_table(adj.otu.tbl, taxa_are_rows = TRUE)

        ov_qs_save(phyobj, "merged.data.norm.qs")
        microbiome.meta$data <- adj.otu.tbl
        ov_qs_save(microbiome.meta, "microbiome_meta.qs")
        merged.data <- transform_sample_counts(phyobj, function(x) x / sum(x))
        ov_qs_save(merged.data, "merged.data.qs")
        ov_qs_save(TRUE, bridge_out, preset = "fast")
      },
      args = list(wd = getwd(), bridge_in = bridge_in, bridge_out = bridge_out),
      timeout_sec = 600
    )
    return(1);
}

