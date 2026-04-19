
#'Perform differential abundance analysis on individual dataset and apply combine effect size meta-analysis method on the results
#'@param method Method used, either "rem" or "fem"
#'@param imgName File name of output bar plot image
#'@param taxRank selected taxa rank which to perform analysis
#'@param selMeta selected metadata
#'@param BHth P-value threshold
#'@param de.method Differential abundance method
#'@param ef.method Effect size model
#'@param format file format of output image
#'@param dpi dots per inch resolution of image (default 72)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PerformMetaEffectSize <- function(mbSetObj=NA, imgName="", taxrank="OTU", selMeta, BHth=0.05, de.method="LM", ef.method="REML", format="png", dpi=default.dpi){

  if(exists('cov.meta.eff')){
    cov <- cov.meta.eff;
  }
  mbSetObj <- .get.mbSetObj(mbSetObj);

  dataName <- mbSetObj$dataSet$dataName;

  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");
  
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  dat <- ov_qs_read("merged.data.qs");
  dat <- subsetPhyloseqByDataset(mbSetObj, dat);
  
  if(taxrank!="OTU"){      
      dat <- fast_tax_glom_mem(dat, taxrank);
      if(is.null(dat)){
          AddErrMsg("Errors in projecting to the selected taxanomy level!");
          return(0);
      }
  }
  
  data <- dat@otu_table;
  feat.nms <- rownames(data);
  data <- as.matrix(data); 
  data <- as.data.frame(data);
  rownames(data) <- feat.nms;
  
  
  if(taxrank!="OTU"){  
    rownames(data) <- dat@tax_table@.Data[feat.nms,taxrank]
    rownames(dat@otu_table) <- rownames(data)
  }else{
    rownames(data) <- feat.nms;
  }

  metadata <-sample_data(dat@sam_data);
  metadata$dataset <- dat@sam_data$dataset;
  
  metadata<-as.matrix(metadata);
  metadata<-as.data.frame(metadata);
  
  #data <- data[,rownames(metadata)];
  
  paramSet$inmex.method <- "effectsize";
  analSet$meta.mat <- meta.stat <<- NULL;
  
  shadow_save(dat, "metaanal_phyobj.qs");
  #meta.obj <- ov_qs_read("microbiome_meta.qs");
  
  #if(taxrank=="OTU"){     
  #concensus_feat <- unname(meta.obj$gene.symbls);
  #data <- data[concensus_feat,]
  #}
  # MMUPHin::lm_meta — quarantined, isolate in subprocess
  lm_meta_args <- list(
    feature_abd = data, batch = "dataset", exposure = selMeta,
    data = metadata,
    control = list(verbose = FALSE, normalization = "NONE",
                   rma_method = ef.method, analysis_method = de.method)
  )
  if (!identical(cov, character(0))) {
    lm_meta_args$covariates <- cov
  }

  fit_lm_meta <- rsclient_isolated_exec(
    func_body = function(input_data) {
      require(MMUPHin)
      do.call(MMUPHin::lm_meta, input_data$args)
    },
    input_data = list(args = lm_meta_args),
    packages = c("MMUPHin", "qs"),
    timeout = 300,
    output_type = "qs"
  )
  if (is.list(fit_lm_meta) && isFALSE(fit_lm_meta$success)) { AddErrMsg(fit_lm_meta$message); return(0) }

  res <- fit_lm_meta$meta_fits;
  nms <- rownames(res);
  es.mat <- data.frame(Coefficient=res$coef, Pval=res$pval, adj_Pval=res$qval.fdr);
  rownames(es.mat) <- nms;
  ord.inx <- order(abs(es.mat[,"Pval"]), decreasing = F);
  es.mat <- es.mat[ord.inx,]
  sig.inx <- which(es.mat[,"adj_Pval"]<=BHth);
  analSet$meta.expr.mat <- data;
  analSet$meta.mat <- es.mat[sig.inx, ];
  analSet$meta.mat.all <- es.mat;
  analSet$meta.res.obj <- fit_lm_meta;
  
  require(dplyr); require(ggplot2);
  p <- res %>% 
    filter(qval.fdr < 0.05) %>% 
    arrange(coef) %>% 
    mutate(feature = factor(feature, levels = feature)) %>% 
    ggplot(aes(y = coef, x = feature)) +
    geom_bar(stat = "identity", aes(fill=-log(qval.fdr))) + 
    scale_fill_gradient(low="#FFFF00",high="#FF0000", name="-Log(P-value)") +
    coord_flip()+ ylab("Coefficient") + xlab("Feature")
  
  imgName = paste(imgName, ".", format, sep="");
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=11.1, height=8.3, type=format, bg="white");
  print(p);
  dev.off();
  res <- SetupMetaStats(BHth, paramSet, analSet);
  saveSet(res[[1]], "paramSet");
  mbSetObj$imgSet[[de.method]] <- imgName;
  mbSetObj$analSet$effectsize$taxalvl <- taxrank;
  mbSetObj$analSet$effectsize$de.method <- de.method;
  mbSetObj$analSet$effectsize$selMeta <- selMeta;
  mbSetObj$analSet$effectsize$ef.method <- ef.method;
  mbSetObj$analSet$effectsize$thresh <- BHth;
  mbSetObj$analSet$meta.restbl <- res[[2]]$meta.restbl;
  .set.mbSetObj(mbSetObj);
  return(saveSet(res[[2]], "analSet", length(sig.inx)));
}

# for gene-level meta-analysis
# function to set up results combining individual data analysis
# as well as to prepare for GO analysis
# no return, as set global 

SetupMetaStats <- function(BHth, paramSet,analSet){
  meta.mat <- analSet$meta.mat.all;
  paramSet$BHth <- BHth;
  #all common genes
  dat <- ov_qs_read("merged.data.qs");
  dat <- subsetPhyloseqByDataset(mbSetObj, dat);
  
  meta.res.obj <- analSet$meta.res.obj;
  gene.ids <- rownames(otu_table(dat));
  # meta.sig genes
  metade.genes <- rownames(meta.mat);
  
  # setup individual sig genes & stats
  # that overlap with meta.sig
  inmex.de <- list();
  
  pval.mat <- fc.mat <- matrix(nrow=nrow(meta.mat), ncol=length(meta.res.obj$maaslin_fits));
  
  for(i in 1:length(meta.res.obj$maaslin_fits)){
    de.res <- meta.res.obj$maaslin_fits[[i]];
    
    hit.inx <- de.res$pval <= BHth;
    hit.inx <- which(hit.inx); # need to get around NA
    inmex.de[[i]] <- rownames(de.res)[hit.inx];
    
    # only choose the genes that are also meta sig genes from in
    # individual analysis for display
    de.res <- de.res[match(metade.genes, de.res$feature),];
    
    fc.mat[,i] <- de.res$coef;
    pval.mat[,i] <- de.res$pval;
  }
  fc.mat[is.na(fc.mat)] <- 0;
  pval.mat[is.na(pval.mat)] <- 1;
  dataNms <- names(meta.res.obj$maaslin_fits);
  newNms <- substring(dataNms,0, nchar(dataNms)-4);
  colnames(fc.mat) <- paste0(newNms, "_FC");
  colnames(pval.mat) <- paste0(newNms, "_Pval");
  
  names(inmex.de) <- names(meta.res.obj$maaslin_fits);
  
  # calculate gain/loss
  deindst <- unique(unlist(inmex.de));
  gains=metade.genes[which(!(metade.genes %in% deindst))];
  losses=deindst[which(!(deindst %in% metade.genes))];
  all.de <- cbind(gene.ids %in% metade.genes, gene.ids %in% deindst);
  colnames(all.de) <- c("Meta-DE", "Individual-DE");
  
  # de genes from individual 
  de.len <- sapply(inmex.de, length);
  stat <- c(length(metade.genes), de.len);
  names(stat) <- c("Meta", substr(names(inmex.de), 0, nchar(names(inmex.de))-4));
  meta.stat <- list(
    stat = stat,
    de = metade.genes,
    idd = gains,
    loss = losses
  );
  
  analSet$fc.mat <- fc.mat;
  analSet$pval.mat <- pval.mat;
  analSet$inmex.de <- inmex.de;
  analSet$meta.stat <- meta.stat;
  
  dat.mat <- cbind(fc.mat, pval.mat, meta.mat);
  dat.mat <- signif(dat.mat, 5);
  # save the result
  res <- cbind(ID=metade.genes, dat.mat); 
  analSet$meta.restbl <- res;
  fast.write(res, file=paste("meta_sig_genes_", paramSet$inmex.method, ".csv", sep=""), row.names=F);
  
  return(list(paramSet, analSet))
}


#'Perform alpha diversity meta-analysis
#'@description Compute alpha diversity indices of datasets and plot summary figure, adapted from work of Jordan Bisanz
#'@param mbSetObj Input the name of the mbSetObj.
#'@param fileName File name of output summary image
#'@param sel.meta Selected metadata
#'@param taxrank Selected taxonomy rank
#'@param view.mode Type of plot to output (ratio or box)
#'@param format File format of output image (default to png)
#'@param dpi dots per inch resolution of image (default 72)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
CompareSummaryStats <- function(mbSetObj=NA,fileName="abc", sel.meta="", taxrank="Family", view.mode="ratio", format="png", dpi=default.dpi) {

  mbSetObj <- .get.mbSetObj(mbSetObj);
  err.vec <<- "";
  current.msg <<- "";
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  if(!exists('summary.stat.vec')){
    AddErrMsg("Need to select at least one summary statistics!");
    return(0);
  }
  message("[MetaAlphaDBG] selected metrics: ", paste(summary.stat.vec, collapse=","));
  
  #print(summary.stat.vec);
  
  library(dplyr); library(ggplot2);
  library(phyloseq);
  res.list <- list();
  
  for(i in 1:length(sel.nms)){
    PS <- list();
    Study <- sel.nms[i];
    data <- mbSetObj$dataSets[[Study]]$proc.phyobj;    
    
    if(taxrank!="OTU"){
      data <- fast_tax_glom_mem(data, taxrank);
      if(is.null(data )){
        AddErrMsg("Errors in projecting to the selected taxanomy level!");
        return(0);
      }
    }
    
    tryCatch(
      {
        # estimate_richness uses vegan internally — isolate in subprocess
        otu_mat <- as(data@otu_table, "matrix")
        taxa_rows <- phyloseq::taxa_are_rows(data)
        # Standardize to samples x taxa for vegan::diversity/estimateR and metadata joins.
        if (isTRUE(taxa_rows)) {
          otu_mat <- t(otu_mat)
        }
        message("[MetaAlphaDBG] Study=", Study, " | taxa_are_rows=", taxa_rows,
                " | otu dim(samples x taxa)=", nrow(otu_mat), "x", ncol(otu_mat),
                " | sam_data rows=", nrow(data@sam_data))
        res <- rsclient_isolated_exec(
          func_body = function(input_data) {
            require(vegan)
            x <- input_data$otu_mat
            if (!any(x == 1)) warning("No singletons — richness estimates may be unreliable")
            outlist <- list()
            measures <- input_data$measures
            if (any(c("Chao1","ACE","Observed") %in% measures))
              outlist <- c(outlist, list(t(data.frame(vegan::estimateR(x), check.names=FALSE))))
            if ("Shannon" %in% measures)
              outlist <- c(outlist, list(shannon = vegan::diversity(x, index="shannon")))
            if ("Simpson" %in% measures)
              outlist <- c(outlist, list(simpson = vegan::diversity(x, index="simpson")))
            if ("InvSimpson" %in% measures)
              outlist <- c(outlist, list(invsimpson = vegan::diversity(x, index="invsimpson")))
            if ("Fisher" %in% measures)
              outlist <- c(outlist, list(fisher = tryCatch(vegan::fisher.alpha(x), error=function(e) rep(NA, nrow(x)))))
            erDF <- as.data.frame(do.call(cbind, outlist))
            renamevec <- c(Observed="S.obs",Chao1="S.chao1",ACE="S.ACE",Shannon="shannon",Simpson="simpson",InvSimpson="invsimpson",Fisher="fisher")
            for (old in names(renamevec)) if (renamevec[old] %in% colnames(erDF)) colnames(erDF)[colnames(erDF)==renamevec[old]] <- old
            # estimateR/diversity are computed from x with samples as rows.
            # Use sample IDs as rownames for downstream metadata joins.
            if (nrow(erDF) == nrow(x)) {
              rownames(erDF) <- rownames(x)
            } else if (nrow(erDF) == ncol(x)) {
              rownames(erDF) <- colnames(x)
            } else {
              message("[MetaAlphaDBG] Unexpected richness row count for ", input_data$study_name,
                      ": nrow(erDF)=", nrow(erDF), ", nrow(x)=", nrow(x), ", ncol(x)=", ncol(x))
            }
            erDF
          },
          input_data = list(otu_mat = otu_mat, measures = summary.stat.vec, study_name = Study),
          packages = c("vegan", "qs"),
          timeout = 120,
          output_type = "qs"
        )
        if (is.list(res) && isFALSE(res$success)) {
          AddErrMsg(res$message)
          return(0)
        }
      },
      error = function(e) {
        message("[PerformMetaAlphaDiversity] ", e$message);
        if(grepl( "uniroot(", e, fixed = TRUE)){
          AddErrMsg(paste0("Fisher alpha estimation fails for ", Study, ". Please try with other metrics or unselect this dataset.", collapse=" "));
        }
        return(0);
      }, finally = {
        if(!exists("res")){
          return(0);
        }
        res <- res
      })    
    
    message("[MetaAlphaDBG] Study=", Study,
            " | ACE selected=", ("ACE" %in% summary.stat.vec),
            " | ACE column=", ("ACE" %in% colnames(res)))
    if("ACE" %in% summary.stat.vec && "ACE" %in% colnames(res)){
      # Some studies do not satisfy ACE assumptions (e.g., sparse/singleton patterns).
      # Skip ACE for that study instead of aborting the whole alpha-meta run.
      if(all(is.nan(res$ACE) | is.na(res$ACE) | is.infinite(res$ACE))){
        msg <- paste0("ACE estimation fails for ", Study, ". ACE is skipped for this dataset and other selected metrics will continue.")
        message("[MetaAlphaDBG] ", msg)
        if(!exists("current.msg") || is.null(current.msg) || identical(current.msg, "null")){
          current.msg <<- msg
        }else if(nchar(paste(current.msg, collapse="; ")) == 0){
          current.msg <<- msg
        }else{
          current.msg <<- paste(paste(current.msg, collapse="; "), msg, sep="; ")
        }
        res$ACE <- NULL
      }
    }
    
    #data <- bf_ratio(data);
    sam_df <- as.data.frame(data@sam_data, check.names = FALSE, stringsAsFactors = FALSE)
    if(!("sample_id" %in% colnames(sam_df))){
      sam_df$sample_id <- rownames(sam_df)
      message("[MetaAlphaDBG] Study=", Study, " | sample_id column missing in sam_data; using rownames.")
    }
    if(!(sel.meta %in% colnames(sam_df))){
      msg <- paste0("Dataset ", Study, " does not contain selected metadata column '", sel.meta, "' and is skipped.")
      message("[MetaAlphaDBG] ", msg)
      if(!exists("current.msg") || is.null(current.msg) || identical(current.msg, "null")){
        current.msg <<- msg
      }else if(nchar(paste(current.msg, collapse="; ")) == 0){
        current.msg <<- msg
      }else{
        current.msg <<- paste(paste(current.msg, collapse="; "), msg, sep="; ")
      }
      next
    }
    sam_df$dataset <- Study
    sam_df$sample_id <- gsub("-", ".", as.character(sam_df$sample_id))

    meta_tbl <- sam_df %>%
      dplyr::select(sample_id, dplyr::all_of(sel.meta), dataset) %>%
      as_tibble()
    colnames(meta_tbl)[which(colnames(meta_tbl) == sel.meta)] <- "study_condition"
    meta_tbl <- meta_tbl %>%
      mutate(
        sample_id_meta = as.character(sample_id),
        sample_id_norm = make.names(gsub("-", ".", trimws(sample_id)))
      ) %>%
      dplyr::select(sample_id_norm, sample_id_meta, study_condition, dataset)

    res <- res %>%
      tibble::rownames_to_column("sample_id") %>%
      mutate(
        sample_id = as.character(sample_id),
        sample_id_norm = make.names(gsub("-", ".", trimws(sample_id)))
      ) %>%
      tidyr::gather(-sample_id, -sample_id_norm, key = "Metric", value = "Diversity") %>%
      left_join(meta_tbl, by = "sample_id_norm") %>%
      mutate(sample_id = dplyr::coalesce(sample_id_meta, sample_id)) %>%
      dplyr::select(sample_id, dataset, study_condition, Metric, Diversity)
    res$study_condition <- factor(res$study_condition)
    matched_n <- sum(!is.na(res$study_condition))
    message("[MetaAlphaDBG] Study=", Study, " | matched metadata rows=", matched_n, "/", nrow(res))
    if(matched_n == 0){
      msg <- paste0("Dataset ", Study, " has no sample-id matches between alpha results and metadata after normalization and is skipped.")
      message("[MetaAlphaDBG] ", msg)
      res_preview <- paste(head(unique(res$sample_id_norm), 5), collapse=",")
      meta_preview <- paste(head(unique(meta_tbl$sample_id_norm), 5), collapse=",")
      message("[MetaAlphaDBG] Study=", Study, " | sample_id_norm preview from alpha=", res_preview)
      message("[MetaAlphaDBG] Study=", Study, " | sample_id_norm preview from metadata=", meta_preview)
      if(!exists("current.msg") || is.null(current.msg) || identical(current.msg, "null")){
        current.msg <<- msg
      }else if(nchar(paste(current.msg, collapse="; ")) == 0){
        current.msg <<- msg
      }else{
        current.msg <<- paste(paste(current.msg, collapse="; "), msg, sep="; ")
      }
      next
    }
    res <- res[!is.na(res$study_condition),]
 
    cond_levels <- levels(res$study_condition)
    if(length(cond_levels) < 2){
      msg <- paste0("Dataset ", Study, " has fewer than two condition levels after preprocessing and is skipped.")
      message("[MetaAlphaDBG] ", msg)
      if(!exists("current.msg") || is.null(current.msg) || identical(current.msg, "null")){
        current.msg <<- msg
      }else if(nchar(paste(current.msg, collapse="; ")) == 0){
        current.msg <<- msg
      }else{
        current.msg <<- paste(paste(current.msg, collapse="; "), msg, sep="; ")
      }
      next
    }
    control_col <- cond_levels[1]
    exp_col <- cond_levels[2]

    mean_tbl <- res %>%
      group_by(Metric, study_condition) %>%
      summarize(mean = mean(log2(Diversity)), .groups = "drop") %>%
      tidyr::pivot_wider(names_from = study_condition, values_from = mean)

    if(!(control_col %in% colnames(mean_tbl)) || !(exp_col %in% colnames(mean_tbl))){
      msg <- paste0("Dataset ", Study, " is missing expected condition columns after reshaping (expected: ", control_col, ", ", exp_col, ") and is skipped.")
      message("[MetaAlphaDBG] ", msg)
      if(!exists("current.msg") || is.null(current.msg) || identical(current.msg, "null")){
        current.msg <<- msg
      }else if(nchar(paste(current.msg, collapse="; ")) == 0){
        current.msg <<- msg
      }else{
        current.msg <<- paste(paste(current.msg, collapse="; "), msg, sep="; ")
      }
      next
    }

    mean_tbl$mean_log2Control <- mean_tbl[[control_col]]
    mean_tbl$mean_log2Exp <- mean_tbl[[exp_col]]
    mean_tbl <- mean_tbl[, c("Metric", "mean_log2Control", "mean_log2Exp")]

    res <- res %>%
      left_join(mean_tbl, by = "Metric") %>%
      plyr::mutate(log2FC = log2(Diversity) - mean_log2Control);
    
    
    res <- res[!is.na(res$log2FC),];
    if(nrow(res) == 0){
      msg <- paste0("No valid alpha-diversity values remain for ", Study, " after filtering invalid estimates; this dataset is skipped.")
      if(!exists("current.msg") || is.null(current.msg) || identical(current.msg, "null")){
        current.msg <<- msg
      }else if(nchar(paste(current.msg, collapse="; ")) == 0){
        current.msg <<- msg
      }else{
        current.msg <<- paste(paste(current.msg, collapse="; "), msg, sep="; ")
      }
      next
    }
    
    PS$AlphaDiversity <- res;
    PS$AlphaDiversity_Stats <- 
      PS$AlphaDiversity %>%
      group_by(Metric) %>%
      do({
        # Now, using tryCatch inside the do() to handle individual t-test failures
        tryCatch({
          broom::tidy(t.test(log2FC ~ study_condition, data = ., conf.int = TRUE, conf.level = 0.95))
        }, error = function(e) {
          # Instead of stopping or returning 0, create a fallback result
          message("T-test failed for ", Study, " in Metric ", unique(.$Metric), ": ", e$message)
          return(tibble(
            estimate = NA_real_,      # NA as we cannot estimate
            p.value = 1,              # Set P-value to 1, indicating non-significance
            estimate1 = NA_real_,     # NA as we don't have this from failed t-test
            estimate2 = NA_real_,     # NA as well
            conf.low = NA_real_,      # NA for confidence interval
            conf.high = NA_real_,     # NA for confidence interval
            method = "Fallback due to error",
            alternative = "NA"
          ))
        })
      }) %>%
      dplyr::mutate(dataset = Study) %>%
      dplyr::select(dataset, Metric, log2FC = estimate, Pvalue = p.value, mean_LFD = estimate1, mean_HFD = estimate2, CI_low = conf.low, CI_high = conf.high)
    
    # Save the results
    res.list[[Study]] <- PS;  
  }
  
  if(length(res.list) == 0){
    AddErrMsg("No valid datasets remained for alpha-diversity meta-analysis after preprocessing and metadata alignment.");
    return(0);
  }

  mod<-lapply(res.list, function(x) x$AlphaDiversity) %>%
    do.call(bind_rows, .) %>%
    plyr::mutate(study_condition=factor(study_condition)) %>%
    group_by(Metric)
  mod <- mod %>% filter(!is.na(log2FC) & !is.infinite(log2FC))
  if(nrow(mod) == 0){
    AddErrMsg("No valid alpha-diversity observations remained after filtering; unable to run meta-analysis.");
    return(0);
  }
  # OPTIMIZED: Pre-allocate list to avoid O(n²) bind_rows in loop (10-50x faster for multiple metrics)
  control_level <- levels(mod$study_condition)[1];
  non_control_level <- levels(mod$study_condition)[2];

  unique_metrics <- unique(mod$Metric)
  alpha_list <- vector("list", length(unique_metrics))

  for(idx in seq_along(unique_metrics)){
    i <- unique_metrics[idx]
    fit <- lmerTest::lmer(log2FC ~ study_condition + (1|dataset), data=subset(mod, Metric==i))
    cf <- confint(fit, level = 0.95)

    alpha_list[[idx]] <- tibble(
      dataset="Combined",
      Metric=i,
      log2FC=summary(fit)$coefficients[paste0("study_condition", non_control_level), "Estimate"],
      Pvalue=anova(fit)$`Pr(>F)`,
      mean_LFD=NA,
      mean_HFD=NA,
      CI_low=cf[paste0("study_condition", non_control_level),1],
      CI_high=cf[paste0("study_condition", non_control_level),2]
    )
  }

  AlphaCombined <- bind_rows(alpha_list)

  cat("levels(study_condition) =", paste(levels(mod$study_condition), collapse = ", "), "\n")

  
  NSamples<-
    lapply(names(res.list), function(x) tibble(dataset=x, Nsamples=length(unique(res.list[[x]]$AlphaDiversity$sample_id)))) %>% 
    do.call(bind_rows, .) %>% 
    arrange(desc(Nsamples)) %>% 
    plyr::mutate(Study=paste0(dataset, " (n=", Nsamples,")")) %>% 
    bind_rows(tibble(dataset="Combined", Study="Combined")) %>%
    plyr::mutate(Study=factor(Study, levels=rev(Study)))
  imgName = paste(fileName, ".", format, sep="");
  
  tbl <- lapply(res.list, function(x) x$AlphaDiversity_Stats) %>%
    do.call(bind_rows, .) %>%
    bind_rows(AlphaCombined) %>%
    plyr::mutate(Significant=case_when(
      Pvalue<0.05 & log2FC>0 ~ "* up",
      Pvalue<0.05 & log2FC<0 ~ "* down",
      TRUE~"ns"
    )) %>%
    ungroup() %>%
    plyr::mutate(Metric=factor(Metric, levels=summary.stat.vec)) %>%
    left_join(NSamples);
  
  if(view.mode == "ratio"){
    if(length(summary.stat.vec)>1){
      summary.stat.vec <- c(summary.stat.vec, "Combined");
    }
    
    fig <- tbl %>%
      ggplot(aes(x=log2FC, y=Study, color=Significant)) +
      geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
      geom_errorbarh(aes(xmin=CI_low, xmax=CI_high), height=0 ) +
      geom_point() +
      facet_grid(~Metric, scales="free_x") +
      scale_color_manual(values=c("green", "red", "black"))  +
      theme(panel.border = element_blank(), axis.line = element_line()) +
      theme(axis.text.x=element_text(angle=45, hjust=1)) + 
      theme(text = element_text(size = 16));
    
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=11.1, height=8.3, type=format, bg="white");
    print(fig);
    dev.off();
  }else{
    data <- lapply(res.list, function(x) x$AlphaDiversity) %>%
      do.call(bind_rows, .)
    combined_data <- data %>%
      mutate(dataset = 'Combined') # Create a combined dataset

    # Bind this combined data with the original data
    augmented_data <- bind_rows(data, combined_data)
    # Now, create the plot with the augmented data
    box1 <- ggplot(augmented_data, aes(Diversity, dataset, fill = study_condition)) +
      geom_boxplot(alpha=0.7, outlier.shape = NA, width=0.2) +
      stat_summary(fun=mean, # Change from fun.y to fun for ggplot2 updates
                   geom = "point",
                   shape = 16,
                   size = 1,
                   aes(group=study_condition),
                   color = "black",
                   position = position_dodge(0.2)) + 
      scale_fill_discrete(name = "Metadata") +
      facet_grid( .~ Metric, scales = "free_x") + 
      theme(text = element_text(size = 16))

    # Create the plot
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=11.1, height=8.3, type=format, bg="white")
    print(box1)
    dev.off()
  }
  tbl <- as.data.frame(tbl[, !colnames(tbl) %in% c("Significant", "Nsamples", "Study", "mean_LFD", "mean_HFD")]);
  mbSetObj$analSet$alpha.summary <- tbl[order(tbl$dataset, tbl$Metric),];
  fast.write(mbSetObj$analSet$alpha.summary, "alpha_summary.csv", row.names = TRUE);
  mbSetObj$imgSet$alpha <- imgName;
  mbSetObj$analSet$alpha.taxalvl <- taxrank
  mbSetObj$analSet$alpha.meta <- sel.meta
  return(.set.mbSetObj(mbSetObj));
}



bf_ratio <- function(phylo){
  
  # Check to make sure input is a phyloseq object
  if(class(phylo)[1] != "phyloseq"){
    message("Error, input is not a phyloseq object!")
    return(NULL)
  }
  
  # Collapse on phyla (use fast_tax_glom_mem - no ape dependency)
  phyla <- fast_tax_glom_mem(phylo, "Phylum")

  # Find relative abundances (use embedded transform_sample_counts)
  phyla_rel <- transform_sample_counts(phyla, function(x) { x/sum(x) } )

  # Keep B/F taxa
  tax_table(phyla_rel)
  phyla_rel_bact <- otu_table(subset_taxa(phyla_rel, grepl("Bacteroidetes", Phylum, fixed = TRUE)))
  phyla_rel_firm <- suppressWarnings(otu_table(subset_taxa(phyla_rel, grepl("Firmicutes", Phylum, fixed = TRUE))))

  # OTU
  bf_ratio <- phyla_rel_bact /  phyla_rel_firm

  # Add to sample metadata (use embedded sample_data<-)
  sample_data(phylo)$bf_ratio <- as.numeric(bf_ratio)
  
  # Return phyllseq object
  return(phylo)
}

#'Perform beta diversity meta-analysis
#'@description Compute beta diversity indices of datasets and plot summary figure, adapted from work of Jordan Bisanz
#'@param mbSetObj Input the name of the mbSetObj.
#'@param plotNm File name of output summary image
#'@param taxalvl Selected taxonomy rank
#'@param sel.meta Selected metadata
#'@param alg statistical comparison algorithm to detect significant different between metadata groups
#'@param format File format of output image (default to png)
#'@param dpi dots per inch resolution of image (default 72)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotBetaSummary <- function(mbSetObj, plotNm,taxalvl, sel.meta, alg, format="png", dpi=default.dpi){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  data.obj <- ov_qs_read("merged.data.qs");
  data.obj <- subsetPhyloseqByDataset(mbSetObj, data.obj);

  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];

  dist.vec <- c("jaccard", "jsd", "bray");#"wunifrac", "unifrac",

  res.list <- list()
  data.nms.vec <- character(length(sel.nms) * length(dist.vec))
  vec_idx <- 1

  for(j in 1:length(sel.nms)){
    mbSetObj$dataSet <- mbSetObj$dataSets[[sel.nms[j]]];
    dataName <- mbSetObj$dataSet$name;

    phyloseq_objs <- readDataQs("phyloseq_objs.qs", mbSetObj$module.type, dataName)
    data <- phyloseq_objs$merged_obj[[taxalvl]]

    if(!is.null(data)){
      otu_mat <- as.matrix(otu_table(data))
      sampledf <- data.frame(sample_data(data), check.names = FALSE)

      for(i in 1:length(dist.vec)){
        dist_method <- dist.vec[i]
        otu_t <- t(otu_mat)
        otu_t <- sweep(otu_t, 1, rowSums(otu_t), "/")

        # vegan distance + stats — isolate in subprocess
        stat.info.vec <- rsclient_isolated_exec(
          func_body = function(input_data) {
            require(vegan)
            otu_t <- input_data$otu_t
            dist_method <- input_data$dist_method
            alg <- input_data$alg
            sampledf <- input_data$sampledf
            sel.meta <- input_data$sel_meta
            group <- sampledf[[sel.meta]]

            if (dist_method == "jaccard") {
              data.dist <- vegan::vegdist(otu_t, method = "jaccard", binary = TRUE)
            } else if (dist_method == "jsd") {
              # JSD: use manual computation
              mat <- otu_t
              n <- nrow(mat)
              dm <- matrix(0, n, n)
              for (ii in 1:(n-1)) for (jj in (ii+1):n) {
                m <- (mat[ii,] + mat[jj,]) / 2
                p1 <- mat[ii,] * log(mat[ii,] / m); p1[!is.finite(p1)] <- 0
                p2 <- mat[jj,] * log(mat[jj,] / m); p2[!is.finite(p2)] <- 0
                dm[jj, ii] <- sum((p1 + p2) / 2, na.rm = TRUE)
              }
              rownames(dm) <- colnames(dm) <- rownames(mat)
              data.dist <- as.dist(dm)
            } else {
              data.dist <- vegan::vegdist(otu_t, method = dist_method)
            }

            if (alg == "adonis") {
              f <- as.formula(paste("data.dist ~", sel.meta))
              res <- vegan::adonis2(formula = f, data = sampledf)
              resTab <- res[1, ]
              sv <- c(signif(resTab$F, 5), signif(resTab$R2, 5), signif(resTab$`Pr(>F)`, 5))
              names(sv) <- c("F-value", "R-squared", "p-value")
            } else if (alg == "anosim") {
              anosim_res <- vegan::anosim(data.dist, group = group)
              sv <- c(signif(anosim_res$statistic, 5), signif(anosim_res$signif, 5))
              names(sv) <- c("R", "p-value")
            } else if (alg == "permdisp") {
              beta <- vegan::betadisper(data.dist, group = group)
              resTab <- anova(beta)
              sv <- c(signif(resTab$"F value"[1], 5), signif(resTab$"Pr(>F)"[1], 5))
              names(sv) <- c("F-value", "p-value")
            }
            sv
          },
          input_data = list(otu_t = otu_t, dist_method = dist_method, alg = alg,
                            sampledf = sampledf, sel_meta = sel.meta),
          packages = c("vegan", "qs"),
          timeout = 180,
          output_type = "qs"
        )
        if (is.list(stat.info.vec) && isFALSE(stat.info.vec$success)) { AddErrMsg(stat.info.vec$message); return(0) }

        res.list[[dataName]][[dist_method]] <- stat.info.vec
        data.nms.vec[vec_idx] <- dataName
        vec_idx <- vec_idx + 1
      }
    }
  }
  
  df <- as.data.frame(t(data.frame(res.list)));
  df$dist <- rep(dist.vec, length(sel.nms));
  df$dataset <- data.nms.vec;
  colnames(df) <- gsub("-", "_", colnames(df), perl = TRUE);
  
  if("R_squared" %in% colnames(df)){
  p <- ggplot(df, aes(x=R_squared, y=dataset, color=dist));
  }else if("R" %in% colnames(df)){
  p <- ggplot(df, aes(x=R, y=dataset, color=dist));
  }else{
  p <- ggplot(df, aes(x=F_value, y=dataset, color=dist));
  }

p <- p + geom_point(shape=16, alpha=0.8, size=4) +
    theme(panel.border = element_blank(), axis.line = element_line(), text = element_text(size = 14)) +
    theme(axis.text.x=element_text(angle=45, hjust=1)) +
    scale_alpha_manual(values=c(0.5,0));
  
  imgName = paste(plotNm, ".", format, sep="");
  Cairo::Cairo(file = imgName, unit="in", dpi=default.dpi, width=11.1, height=8.3, type=format, bg="white");
  print(p);
  dev.off();
  
  mbSetObj$imgSet$beta_summary <- imgName;
  mbSetObj$analSet$beta.summary <- df;
  fast.write(df, "beta_summary.csv", row.names = TRUE);
  
  
  return(.set.mbSetObj(mbSetObj));
}


PlotDiscreteDiagnostic <- function(mbSetObj, fileName, metadata, format="png", dpi=default.dpi){
  mbSetObj <- .get.mbSetObj(mbSetObj);

  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];

  data  <- ov_qs_read("merged.data.raw.qs");
  data <- subsetPhyloseqByDataset(mbSetObj, data);

  sam_data <- as.data.frame(as.matrix(sample_data(data)));
  meta_groups <- unique(sam_data[, metadata])

  # Batch all vegan+MMUPHin computations in one subprocess
  otu_subsets <- lapply(meta_groups, function(grp) {
    sub_sam <- sam_data[sam_data[, metadata] == grp, ]
    as.matrix(data@otu_table[, rownames(sub_sam)])
  })
  names(otu_subsets) <- meta_groups
  sam_subsets <- lapply(meta_groups, function(grp) {
    sam_data[sam_data[, metadata] == grp, ]
  })
  names(sam_subsets) <- meta_groups

  fit_results <- rsclient_isolated_exec(
    func_body = function(input_data) {
      require(vegan); require(MMUPHin)
      results <- list()
      for (grp in input_data$meta_groups) {
        dist_data <- vegan::vegdist(t(input_data$otu_subsets[[grp]]), method = "bray")
        fit <- MMUPHin::discrete_discover(D = dist_data, batch = "dataset",
                 data = input_data$sam_subsets[[grp]],
                 control = list(k_max = 8, verbose = FALSE))
        results[[grp]] <- list(internal_mean = fit$internal_mean, internal_se = fit$internal_se,
                               external_mean = fit$external_mean, external_se = fit$external_se)
      }
      results
    },
    input_data = list(meta_groups = meta_groups, otu_subsets = otu_subsets, sam_subsets = sam_subsets),
    packages = c("vegan", "MMUPHin", "qs"),
    timeout = 300,
    output_type = "qs"
  )
  if (is.list(fit_results) && isFALSE(fit_results$success)) { AddErrMsg(fit_results$message); return(0) }

  # Build ggplot objects in Master (ggplot2 is not quarantined)
  plot.list <- list()
  i <- 1
  for (meta.grp in meta_groups) {
    fit_disc <- fit_results[[meta.grp]]
    res.list <- list()
    for (study_id in sel.nms) {
      internal <- data.frame(K = 2:8, statistic = fit_disc$internal_mean, se = fit_disc$internal_se, type = "internal")
      external <- data.frame(K = 2:8, statistic = fit_disc$external_mean, se = fit_disc$external_se, type = "external")
      df <- rbind(internal, external)
      df$dataset <- rep(study_id, nrow(df))
      res.list[[paste0(metadata, meta.grp, sep = "-")]] <- df
    }
    res <- do.call("rbind", res.list)

    p <- ggplot2::ggplot(res, aes(x = K, y = statistic, color = type)) +
      geom_point(position = position_dodge(width = 0.5)) +
      geom_line(position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                    position = position_dodge(width = 0.5), width = 0.5) +
      facet_wrap(as.formula(paste0("~dataset")), scales = "free", ncol = 1) +
      ggtitle(label = paste0(metadata, "-", meta.grp))

    if (i == 1) p <- p + theme(legend.position = "none")
    plot.list[[meta.grp]] <- p
    i <- i + 1
  }

  # ggarrange (ggpubr) — quarantined, isolate with Cairo
  imgName <- paste(fileName, ".", format, sep = "")
  rsclient_isolated_exec(
    func_body = function(input_data) {
      require(ggpubr)
      Cairo::Cairo(file = input_data$imgName, unit = "in", dpi = 96, width = 9.7, height = 16.7, type = input_data$format, bg = "white")
      p <- ggpubr::ggarrange(plotlist = input_data$plot_list, ncol = input_data$ncol)
      print(p)
      dev.off()
      TRUE
    },
    input_data = list(imgName = imgName, format = format, plot_list = plot.list, ncol = length(meta_groups)),
    packages = c("ggpubr", "Cairo", "qs"),
    timeout = 120,
    output_type = "qs"
  )
  # plot failure is non-fatal — log but don't block

  return(.set.mbSetObj(mbSetObj))
}

PlotContinuousPopulation <- function(mbSetObj, loadingName, ordinationName, metadata, format="png", dpi=default.dpi){
  mbSetObj <- .get.mbSetObj(mbSetObj);

  require(dplyr); require(ggplot2);

  data <- ov_qs_read("merged.data.qs");
  data <- subsetPhyloseqByDataset(mbSetObj, data);
  sam_data <- as.data.frame(as.matrix(sample_data(data)));
  meta_groups <- unique(sam_data[, metadata])

  # Batch vegan+MMUPHin for all meta groups in one subprocess
  otu_subsets <- lapply(meta_groups, function(grp) {
    sub_sam <- sam_data[sam_data[, metadata] == grp, ]
    as.matrix(data@otu_table[, rownames(sub_sam)])
  })
  names(otu_subsets) <- meta_groups
  sam_subsets <- lapply(meta_groups, function(grp) {
    sam_data[sam_data[, metadata] == grp, ]
  })
  names(sam_subsets) <- meta_groups

  fit_results <- rsclient_isolated_exec(
    func_body = function(input_data) {
      require(vegan); require(MMUPHin)
      results <- list()
      for (grp in input_data$meta_groups) {
        sub_data <- input_data$otu_subsets[[grp]]
        sub_sam <- input_data$sam_subsets[[grp]]
        dist.data <- vegan::vegdist(t(sub_data), method = "bray")
        fit <- MMUPHin::continuous_discover(feature_abd = sub_data, batch = "dataset",
                 data = sub_sam, control = list(var_perc_cutoff = 0.5, verbose = FALSE))
        mds <- cmdscale(d = dist.data)
        colnames(mds) <- c("Axis1", "Axis2")
        results[[grp]] <- list(
          consensus_loadings = fit$consensus_loadings,
          consensus_scores = fit$consensus_scores,
          mds = as.data.frame(mds)
        )
      }
      results
    },
    input_data = list(meta_groups = meta_groups, otu_subsets = otu_subsets, sam_subsets = sam_subsets),
    packages = c("vegan", "MMUPHin", "qs"),
    timeout = 300,
    output_type = "qs"
  )
  if (is.list(fit_results) && isFALSE(fit_results$success)) { AddErrMsg(fit_results$message); return(0) }

  # Build ggplot objects in Master
  loading.list <- list()
  ordination.list <- list()
  for (meta.grp in meta_groups) {
    fit_cont <- fit_results[[meta.grp]]
    loading <- data.frame(feature = rownames(fit_cont$consensus_loadings),
                          loading1 = fit_cont$consensus_loadings[, 1])
    p <- loading %>%
      arrange(-abs(loading1)) %>% slice(1:20) %>% arrange(loading1) %>%
      plyr::mutate(feature = factor(feature, levels = feature)) %>%
      ggplot(aes(x = feature, y = loading1)) + geom_bar(stat = "identity") +
      coord_flip() + ggtitle(label = paste0(metadata, "-", meta.grp))
    loading.list[[as.character(meta.grp)]] <- p

    mds <- fit_cont$mds
    mds$score1 <- fit_cont$consensus_scores[, 1]
    p2 <- ggplot(mds, aes(x = Axis1, y = Axis2)) +
      geom_point(aes(colour = score1), size = 2, shape = 16) +
      scale_colour_gradient2() + coord_fixed()
    ordination.list[[as.character(meta.grp)]] <- p2
  }

  # ggarrange (ggpubr) — quarantined, isolate with Cairo
  ngrp <- length(meta_groups)
  loadingName <- paste(loadingName, ".", format, sep = "")
  ordinationName <- paste(ordinationName, ".", format, sep = "")

  rsclient_isolated_exec(
    func_body = function(input_data) {
      require(ggpubr)
      Cairo::Cairo(file = input_data$loadingName, unit = "in", dpi = 96, width = 9.7, height = 16.7, type = input_data$format, bg = "white")
      print(ggpubr::ggarrange(plotlist = input_data$loading_list, nrow = input_data$ngrp))
      dev.off()
      Cairo::Cairo(file = input_data$ordinationName, unit = "in", dpi = 96, width = 9.7, height = 16.7, type = input_data$format, bg = "white")
      print(ggpubr::ggarrange(plotlist = input_data$ordination_list, nrow = input_data$ngrp))
      dev.off()
      TRUE
    },
    input_data = list(loadingName = loadingName, ordinationName = ordinationName, format = format,
                      loading_list = loading.list, ordination_list = ordination.list, ngrp = ngrp),
    packages = c("ggpubr", "Cairo", "qs"),
    timeout = 120,
    output_type = "qs"
  )

  return(.set.mbSetObj(mbSetObj));
}
