
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
PerformMetaEffectSize <- function(mbSetObj=NA, imgName="", taxrank="OTU", selMeta, BHth=0.05, de.method="LM", ef.method="REML", format="png", dpi=100){

  if(exists('cov.meta.eff')){
    cov <- cov.meta.eff;
  }
  mbSetObj <- .get.mbSetObj(mbSetObj);

  dataName <- mbSetObj$dataSet$dataName;

  paramSet <- readSet(paramSet, "paramSet");
  analSet <- readSet(analSet, "analSet");
  
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  dat <- qs::qread("merged.data.qs");
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
  
  qs::qsave(dat, "metaanal_phyobj.qs");
  #meta.obj <- qs::qread("microbiome_meta.qs");
  
  #if(taxrank=="OTU"){     
  #concensus_feat <- unname(meta.obj$gene.symbls);
  #data <- data[concensus_feat,]
  #}
  library(MMUPHin);
  
    if(identical(cov, character(0))){
      fit_lm_meta <- lm_meta(feature_abd = data,
                             batch = "dataset",
                             exposure = selMeta,
                             data = metadata,
                             control = list(verbose = FALSE, normalization="NONE", 
                            #transform="NONE", 
                            rma_method=ef.method, analysis_method=de.method));
    }else{
      fit_lm_meta <- lm_meta(feature_abd = data,
                             batch = "dataset",
                             exposure = selMeta,
                             data = metadata,
                             covariates = cov,
                             control = list(verbose = FALSE, normalization="NONE", 
                          #  transform="NONE", 
                            rma_method=ef.method, analysis_method=de.method));
    }
  
  
  
  
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
  
  require(tidyverse);
  p <- res %>% 
    filter(qval.fdr < 0.05) %>% 
    arrange(coef) %>% 
    mutate(feature = factor(feature, levels = feature)) %>% 
    ggplot(aes(y = coef, x = feature)) +
    geom_bar(stat = "identity", aes(fill=-log(qval.fdr))) + 
    scale_fill_gradient(low="#FFFF00",high="#FF0000", name="-Log(P-value)") +
    coord_flip()+ ylab("Coefficient") + xlab("Feature")
  
  imgName = paste(imgName, ".", format, sep="");
  Cairo::Cairo(file = imgName, unit="px", dpi=dpi, width=800, height=600, type=format, bg="white");
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
  dat <- qs::qread("merged.data.qs");
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
CompareSummaryStats <- function(mbSetObj=NA,fileName="abc", sel.meta="", taxrank="Family", view.mode="ratio", format="png", dpi=100) {
  #save.image("alpha.RData");

cat("R version  :", R.version.string,                "\n")
cat("tidyverse  :", as.character(packageVersion("tidyverse")), "\n")
cat("broom      :", as.character(packageVersion("broom")),     "\n")
cat("Cairo OK?  :", capabilities("cairo"),           "\n")

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  if(!exists('summary.stat.vec')){
    AddErrMsg("Need to select at least one summary statistics!");
    return(0);
  }
  
  #print(summary.stat.vec);
  
  library(tidyverse);
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
        res <- estimate_richness(data@otu_table, measures=summary.stat.vec);
      },
      error = function(e) {
        print(e);
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
    
    if("ACE" %in% colnames(res)){
      if(all(res$ACE %in% NaN)){
        AddErrMsg(paste0("ACE estimation fails for ", Study, ". Please try with other metrics or unselect this dataset."));
        return(0);
      }
    }
    
    #data <- bf_ratio(data);
    data@sam_data$dataset <- rep(Study, nrow(data@sam_data));
    data@sam_data$sample_id <- gsub("-", ".", data@sam_data$sample_id);
    
    
    res <- res %>% rownames_to_column("sample_id") %>%
      select(sample_id, everything()) %>% 
      left_join(data@sam_data[, c("sample_id")]) %>%
      as_tibble() %>%
      gather(-sample_id, key="Metric", value="Diversity") %>%
      left_join(data@sam_data[,c("sample_id",sel.meta, "dataset")]) %>%
      dplyr::select(sample_id, dataset, matches(sel.meta), Metric, Diversity);
    
    colnames(res)[which(colnames(res) == sel.meta)] <- "study_condition";
 
    #R package conflict with rename
    res <-
      res %>%
      left_join(
        res %>%
          group_by(Metric, study_condition) %>%
          summarize(mean=mean(log2(Diversity))) %>%
          spread(key=study_condition, value=mean) %>%
          dplyr::rename(mean_log2Control := !!quo_name(levels(res$study_condition)[1]), mean_log2Exp := !!quo_name(levels(res$study_condition)[2]))
      )  %>%
      plyr::mutate(log2FC=log2(Diversity)-mean_log2Control);
    
    
    res <- res[!is.na(res$log2FC),];
    
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
  
  mod<-lapply(res.list, function(x) x$AlphaDiversity) %>%
    do.call(bind_rows, .) %>%
    plyr::mutate(study_condition=factor(study_condition)) %>%
    group_by(Metric)
  mod <- mod %>% filter(!is.na(log2FC) & !is.infinite(log2FC))
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
    
    Cairo::Cairo(file = imgName, unit="px", dpi=dpi, width=800, height=600, type=format, bg="white");
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
    Cairo::Cairo(file = imgName, unit="px", dpi=dpi, width=800, height=600, type=format, bg="white")
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
  
  # Collapse on phyla
  phyla <- phyloseq::tax_glom(physeq = phylo, taxrank = "Phylum")
  
  # Find relative abundances
  phyla_rel <- phyloseq::transform_sample_counts(phyla, function(x) { x/sum(x) } )
  
  # Keep B/F taxa
  tax_table(phyla_rel)
  phyla_rel_bact <- otu_table(subset_taxa(phyla_rel, grepl("Bacteroidetes", Phylum, fixed = TRUE)))
  phyla_rel_firm <- suppressWarnings(otu_table(subset_taxa(phyla_rel, grepl("Firmicutes", Phylum, fixed = TRUE))))
  
  # OTU
  bf_ratio <- phyla_rel_bact /  phyla_rel_firm
  
  # Add to sample metadata
  phyloseq::sample_data(phylo)$bf_ratio <- as.numeric(bf_ratio)
  
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
PlotBetaSummary <- function(mbSetObj, plotNm,taxalvl, sel.meta, alg, format="png", dpi=100){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  data.obj <- qs::qread("merged.data.qs");
  data.obj <- subsetPhyloseqByDataset(mbSetObj, data.obj);
  
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  dist.vec <- c("jaccard", "jsd", "bray");#"wunifrac", "unifrac",
  res.list <- list();
  # OPTIMIZED: Pre-allocate vector to avoid O(n²) c() growing in nested loop
  data.nms.vec <- character(length(sel.nms) * length(dist.vec))
  vec_idx <- 1
  for(j in 1:length(sel.nms)){
    mbSetObj$dataSet <- mbSetObj$dataSets[[sel.nms[j]]];
    dataName <- mbSetObj$dataSet$name;
    .set.mbSetObj(mbSetObj);
    res.list[[dataName]] <- list();

    for(i in 1:length(dist.vec)){
      distName <- dist.vec[i];
      PerformCategoryComp(mbSetObj, taxalvl, alg,distName, sel.meta);
      mbSetObj <- .get.mbSetObj(mbSetObj);
      res <- mbSetObj$analSet$stat.info.vec;
      res.list[[dataName]][[dist.vec[i]]] <- res;
      data.nms.vec[vec_idx] <- dataName
      vec_idx <- vec_idx + 1
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
  Cairo::Cairo(file = imgName, unit="px", dpi=72, width=800, height=600, type=format, bg="white");
  print(p);
  dev.off();
  
  mbSetObj$imgSet$beta_summary <- imgName;
  mbSetObj$analSet$beta.summary <- df;
  fast.write(df, "beta_summary.csv", row.names = TRUE);
  
  
  return(.set.mbSetObj(mbSetObj));
}


PlotDiscreteDiagnostic <- function(mbSetObj, fileName, metadata, format="png", dpi=72){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  
  require(MMUPHin);
  require(ggpubr);
  require(vegan);
  
  data  <- qs::qread("merged.data.raw.qs");
  data <- subsetPhyloseqByDataset(mbSetObj, data);
  
  sam_data <- as.data.frame(as.matrix(sample_data(data)));
  
  plot.list <- list();
  i <- 1;
  for(meta.grp in unique(sam_data[,metadata])){
    sub_sam_data <- sam_data[which(sam_data[,metadata] == meta.grp), ]
    sub_data <- data@otu_table[, rownames(sub_sam_data)]
    
    dist_data <- vegdist(t(sub_data))
    
    fit_discrete <- discrete_discover(D = dist_data,
                                      batch = "dataset",
                                      data = sub_sam_data,
                                      control = list(k_max = 8, verbose = F));
    nms <- sel.nms;
    res.list <- list();
    for(study_id in nms){
    
    internal <- data.frame(K = 2:8,
                           statistic = fit_discrete$internal_mean,
                           se = fit_discrete$internal_se,
                           type = "internal")
    
    external <- data.frame(K = 2:8,
                           statistic = fit_discrete$external_mean,
                           se = fit_discrete$external_se,
                           type = "external")
    
    df <- rbind(internal, external);
    df$dataset <- rep(study_id, nrow(df));
    res.list[[paste0(metadata, meta.grp, sep="-")]] <- df;
    }  
    res <- do.call("rbind", res.list);
    
    p <- ggplot2::ggplot(res, aes(x = K, y = statistic, color = type)) +
      geom_point(position = position_dodge(width = 0.5)) + 
      geom_line(position = position_dodge(width = 0.5)) +
      geom_errorbar(aes(ymin = statistic - se, ymax = statistic + se),
                    position = position_dodge(width = 0.5), width = 0.5) + 
      facet_wrap(as.formula(paste0("~dataset")), scales="free", ncol= 1) +
      ggtitle(label=paste0(metadata, "-", meta.grp));
    
    if(i == 1){
      p <- p + theme(legend.position = "none")
    }
    
    plot.list[[meta.grp]] <- p;
  }

  imgName = paste(fileName, ".", format, sep="");
  
  Cairo::Cairo(file=imgName, width=700, height=1200, type=format, bg="white",dpi=dpi);
  p <- ggarrange(plotlist=plot.list, ncol = length(unique(sam_data[,metadata])));
  print(p);
  dev.off();
  return(.set.mbSetObj(mbSetObj))
}

PlotContinuousPopulation <- function(mbSetObj, loadingName, ordinationName, metadata, format="png", dpi=72){
  mbSetObj <- .get.mbSetObj(mbSetObj);

  require(vegan);
  require(MMUPHin);
  require(tidyverse);
  
  loading.list <- list();
  ordination.list <- list();
  
  data <- qs::qread("merged.data.qs");
  data <- subsetPhyloseqByDataset(mbSetObj, data);
  
  sam_data <- as.data.frame(as.matrix(sample_data(data)));
  for(meta.grp in unique(sam_data[,metadata])){
    sub_sam_data <- sam_data[which(sam_data[,metadata] == meta.grp), ]
    sub_data <- data@otu_table[, rownames(sub_sam_data)];
    dist.data <- vegdist(t(sub_data));
    
    fit_continuous <- continuous_discover(feature_abd = sub_data,
                                          batch = "dataset",
                                          data = sub_sam_data,
                                          control = list(var_perc_cutoff = 0.5,
                                                         verbose = FALSE));
    
    loading <- data.frame(feature = rownames(fit_continuous$consensus_loadings),
                          loading1 = fit_continuous$consensus_loadings[, 1])
    
    p <- loading %>%
      arrange(-abs(loading1)) %>%
      slice(1:20) %>%
      arrange(loading1) %>%
      plyr::mutate(feature = factor(feature, levels = feature)) %>%
      ggplot(aes(x = feature, y = loading1)) +
      geom_bar(stat = "identity") +
      coord_flip() +
      ggtitle(label=paste0(metadata, "-", meta.grp));
    
    loading.list[[meta.grp]] <- p;
    
    mds <- cmdscale(d = dist.data)
    colnames(mds) <- c("Axis1", "Axis2");
    mds <- as.data.frame(mds);
    mds$score1 <- fit_continuous$consensus_scores[, 1];
    
    p2 <- ggplot(mds, aes(x = Axis1, y = Axis2)) +
      geom_point(aes(colour = score1), size=2, shape=16) +
      scale_colour_gradient2() +
      coord_fixed()
    
    ordination.list[[meta.grp]] <- p2;
  }
  
  
  
  loadingName <-paste(loadingName, ".", format, sep="");
  Cairo::Cairo(file=loadingName, width=700, height=1200, type=format, bg="white",dpi=dpi);
  pl1 <- ggarrange(plotlist=loading.list, nrow = length(unique(sam_data[,metadata])));
  print(pl1);
  dev.off();
  
  ordinationName <-paste(ordinationName, ".", format, sep="");
  Cairo::Cairo(file=ordinationName, width=700, height=1200, type=format, bg="white",dpi=dpi);
  pl2 <- ggarrange(plotlist=ordination.list, nrow = length(unique(sam_data[,metadata])));
  print(pl2);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj));
}

handle_t_test <- function(data) {
  result <- tryCatch(
    {
      broom::tidy(t.test(log2FC ~ study_condition, data = data, conf.int = TRUE, conf.level = 0.95))
    },
    error = function(e) {
      message("Error in t.test: ", e$message)  # Print the error message
      
      # Return a tibble/data.frame with the same structure as broom::tidy(t.test) but with Pvalue set to 1
      tibble::tibble(
        estimate = NA_real_,  # NA values for estimates since the test did not properly execute
        estimate1 = NA_real_,
        estimate2 = NA_real_,
        statistic = NA_real_,
        p.value = 1,  # Set P-value to 1 indicating no significant difference
        parameter = NA_real_,
        conf.low = NA_real_,
        conf.high = NA_real_,
        conf.level = 0.95  # The confidence level you've used
      )
    }
  )
  return(result)
}