#' Perform Partial Correlation Analysis
#' @description Function to perform and plot
#' partial correlations between all taxonomic features,
#' the outcome, and selected confounders.
#' NOTE: All metadata must be numeric
#' @param mbSetObj Input the name of the mbSetObj.
#' @param taxa.lvl Character, input the taxonomic level
#' to perform partial correlation analysis.
#' @param variable Character, input the selected variable.
#' @param alg Use "kendall" or "spearman" for non-parametric and 
#' "pearson" for parametric.
#' @export
#' @import ppcor
PerformPartialCorr <- function(mbSetObj, taxa.lvl="Phylum", variable=NA, alg = "pearson", pval.cutoff = 0.05){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  suppressMessages(library(ppcor));
  suppressMessages(library(viridis));
  
  # retrieve sample info
  metadata <- data.frame(sample_data(mbSet$dataSet$proc.phyobj), check.names=F, stringsAsFactors = FALSE,check.names=FALSE);
  confounders <- mbSetObj$dataSet$confs
  
  # get list of confounders to consider
  if(length(confounders)==0){
    current.msg <<- "No confounders inputted!"
    return(0);
  }
  
  #check that confounders do not include variable
  check <- variable %in% confounders
  
  if(check){
    current.msg <<- "Invalid confounders! Variable included."
    return(0)
  }
  
  check2 <- "NA" %in% confounders
  
  if(check2){
    current.msg <<- "NA included as a confounder!"
    return(0)
  }
  
  # create list of confounders
  meta.subset <- metadata[, confounders]
  
  if("data.frame" %in% class(meta.subset)){
    meta.subset <- data.frame(apply(meta.subset, 2, function(x) as.numeric(as.character(x))),check.names=FALSE);
    meta.subset <- as.list(meta.subset)
  }else{
    meta.subset <- as.numeric(meta.subset)
    meta.subset <- as.list(as.data.frame(meta.subset,check.names=FALSE));
  }
  
  # check variable is numeric
  var <- metadata[,variable]
  
  if(class(levels(var))=="character"){
    var <- as.numeric(var)
  }
  
  # now get otu data
  if(taxa.lvl=="OTU"){
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
    data1 <- as.matrix(otu_table(data));
  }else{
    #get otu table
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data,taxa.lvl);
    nm <- as.character(tax_table(data)[,taxa.lvl]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    #all NA club together
    data1 <- as.matrix(t(sapply(by(data1, rownames(data1), colSums), identity)));
  }
  
  otu.table <- t(data1);
  mbSetObj$analSet$abund_data <- otu.table;
  
  #replace 0s and NAs with small number
  otu.table[otu.table==0|is.na(otu.table)] <- .00001
  
  #more than 1 imp.feat, convert to named list of vectors then perform partial correlation
  if(class(otu.table)=="numeric"){
    otu.subset <- otu.table
    # first calculate corr
    cor.result <- cor.test(otu.subset, var, method=alg)
    cor.results <- data.frame(cor_pval=cor.result$p.value, cor_est=cor.result$estimate,check.names=FALSE)
    # calculate pcorr
    pcor.results <- ppcor::pcor.test(otu.subset, var, meta.subset, method = alg)
    pcor.results <- cbind(cor.results, pcor.results)
  }else{
    otu.subset <- lapply(seq_len(ncol(otu.table)), function(i) otu.table[,i])
    
    # first calculate corr
    cor.results <- do.call(rbind, lapply(otu.subset, function(x){
      cor.result <- cor.test(x, var, method = alg);
      data.frame(cor_pval=cor.result$p.value, cor_est=cor.result$estimate,check.names=FALSE)}))  
    row.names(cor.results) <- colnames(otu.table)
    
    # calculate pcorr
    pcor.fun <- function(feats){ ppcor::pcor.test(feats, var, meta.subset, method = alg) }
    pcor.results <- do.call("rbind", lapply(otu.subset, pcor.fun))
    pcor.results <- cbind(cor.results, pcor.results)
  }
  
  #order results by p.value
  row.names(pcor.results) <- colnames(otu.table)
  resTable <- as.data.frame(pcor.results,check.names=FALSE)
  ord.inx <- order(resTable$p.value);
  resTable <- resTable[ord.inx, , drop=FALSE];
  fast.write(resTable, "partial_corr.csv", row.names = TRUE);
  resTable$taxarank = row.names(pcor.results)
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}

#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxrank Character, input the taxonomic level to use. Default is "OTU".
#'@param metadata Either a character vector containing the metadata variables to be extracted
#'from the uploaded metadata table or a metadata matrix containing different information such
#'as Fatty Acid Profile or Transcriptomics data for the same samples. If a matrix, ensure that
#'the metaDistMatrix is set to TRUE. 
#'@param metaDistMatrix Boolean. If TRUE, the object metadata is a matrix, such as metabolomics data
#'of the same matched samples to perform the Mantel test. If FALSE, the object metadata is a vector
#'of selected variables from the metadata to extract from the uploaded metadata table to use for the 
#'Mantel test.
#'@param microDist Character, input the dissimilarity index used by vegan::vegdist to create
#'a dissimilarity indice for the microbial data matrix. Default is "bray".
#'@param metaDist Character, input the dissimilarity index used by vegan::vegdist to create
#'a dissimilarity indice for the meta-data matrix. Default is "jaccard".
#'@param nperm  Numeric, input the number of permutations to be performed for the Mantel test.
#'Default is set to 100.
#'@param method Correlation method used, as returned by cor.test.
#'@export
PerformMantelTest <- function(mbSet, taxrank = "OTU", microDist = "bray", 
                              metadata, metaDistMatrix = FALSE, metaDist = "jaccard",  
                              nperm = 100, method = "spearman"){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  suppressMessages(library(data.table));
  suppressMessages(library(viridis));
  
  metadata <- metadata;
  
  # using normalized data
  if(taxrank=="OTU"){
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
  }else{
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
    # merging at taxonomy levels
    data <- fast_tax_glom_mem(data,taxrank);
    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
  }
  
  micro.dist <- phyloseq::distance(data, microDist)
  
  micro_meta <- sample_data(data);
  valid.meta <- all(metadata %in% colnames(micro_meta))
  
  # now make meta-data distance matrix
  if(valid.meta & !metaDistMatrix){
    meta_data <- data.frame(micro_meta,check.names=FALSE)[, metadata, drop = FALSE];
    meta_data[] <- lapply(meta_data, as.numeric)
    meta_dist <- vegan::vegdist(meta_data, metaDist)
  }else if(metaDistMatrix){
    meta_data[] <- lapply(metadata, as.numeric)
    meta_dist <- vegan::vegdist(meta_data, metaDist)
  }else{
    AddErrMsg("Invalid metadata! Ensure all metadata match the uploaded metadata file!");
    return(0);
  }
  
  mantel_results <- vegan::mantel(micro.dist, meta_dist, permutations = nperm, method = method)
  
  if(!.on.public.web){
    print(mantel_results)
  }
  
  mantel_df <- data.frame(Mantel_Statistic = mantel_results$statistic, Mantel_Signif = mantel_results$signif, 
                          nPerm = mantel_results$permutations, Method = mantel_results$method,check.names=FALSE)
  
  fast.write(mantel_df, "mantel_results.csv")
  mbSetObj$analSet$mantel <- mantel_df;
  
  return(.set.mbSetObj(mbSetObj))
}

#'@description This function runs a general linear model
#'to investigate if the continuous variable of interest
#'impacts the calculated alpha diversity index.
#'@param variable Input the name of the experimental factor to group the samples.
#'@param var_distribution Input the type of model to be used for the glm model.
#'By default the model will assume that the variable of interest follows a normal
#'distribution, hence "gaussian". Other distributions (non-normal) include "poisson",
#'"quasipoisson" - see ?family for more details.
#'"poisson" - usually count data. If you have overdispersion (if residual deviance is much larger than degrees of freedom), 
#'you may want to use quasipoisson() instead of poisson().
#'@param covariates Boolean. TRUE to consider covariates, FALSE to only consider the
#'main effect. Only valid for PERMANOVA.
#'@param cov.vec Character vector. Input the names of the covariates to consider in the 
#'PERMANOVA model.
#'@param model.additive Boolean. If TRUE, the model will be additive (i.e. data ~ var1 + var2),
#'making the assumption that the two factor variables are independent. 
#'However, if FALSE, the model will consider the synergistic effects of the variables - interaction (i.e. data ~ var1*var2).
#'"Different explanatory variables the effects on the outcome of a change in one variable may either not depend on the
#'level of the other variable (additive model) or it may depend on the level of the other variable (interaction model)." 
#'@usage mbSet <- PerformAlphaGLM(mbSet, "SampleType", covariates = T, cov.vec = c("sex", "age_group"))

PerformAlphaGLM <- function(mbSet, variable, var_distribution="gaussian",
                            covariates = FALSE, cov.vec = NA, model.additive = TRUE){

  mbSetObj <- .get.mbSetObj(mbSet);
  
  # Set up alpha-diversity
  alpha.div <- mbSetObj$analSet$alpha$value
  sampledf <- mbSetObj$analSet$alpha[, 1:(length(mbSetObj$analSet$alpha)-3)]
  
  if(covariates){
    
    valid.meta <- all(cov.vec %in% colnames(sampledf))
    
    if(valid.meta){
      outcome <- "alpha.div"
      variables <- c(variable, cov.vec)
      
      if(model.additive){
        f <- as.formula(
          paste(outcome,
                paste(variables, collapse = " + "),
                sep = " ~ ")
        )
      }else{
        f <- as.formula(
          paste(outcome,
                paste(variables, collapse = " * "),
                sep = " ~ ")
        )
      }
      res <- glm(formula = f, data = sampledf, family = var_distribution)
    }else{
      AddErrMsg("Invalid metadata! Ensure all metadata match the uploaded metadata file!");
      return(0);
    }
  }else{
    res <- glm(alpha.div ~ sampledf[ ,variable], data = sampledf, family = var_distribution);
  }
  
  mbSetObj$analSet$alpha_glm <- res
  results_df <- summary.glm(res)$coefficients
  fast.write(results_df, file="glm_alpha_diversity.csv");
  
  # Now check goodness of fit
  # Code adapted from jtools R package
  
  # https://pj.freefaculty.org/guides/stat/Regression-GLM/GLM2-SigTests/GLM-2-guide.pdf
  # analyze the variance between the null model + full model
  # how much better your model is (residual deviance) than just the intercept (null deviance)
  anal_var <- anova(res, test = "Chisq")
  
  # Model only explains X% of the deviance.
  pseudo_rsquared <- 1-(res$deviance/res$null.deviance)
  
  # Cragg-Uhler r2
  llh <- logLik(res)
  llhNull <- logLik(update(res, ~ 1, data = model.frame(res)))
  n <- nobs(res)
  G2 <- -2*(llhNull-llh)
  r2ML <- 1 - exp(-G2/n)
  r2ML.max <- 1 - exp(llhNull*2/n)
  r2CU <- r2ML/r2ML.max
  
  # McFadden r2
  # Values from 0.2-0.4 indicate (in McFadden's words) excellent model fit
  McFadden <- 1 - llh/llhNull
  
  glm_res <- list(n = n,
                  dep_var = names(res$model[1]),
                  lmFamily = res$family,
                  table = results_df,
                  chisq_dev = na.omit(anal_var$Deviance),
                  chisq_p = na.omit(anal_var$`Pr(>Chi)`),
                  pseudo_r2 = pseudo_rsquared,
                  mcf_r2 = McFadden,
                  cu_r2 = r2CU,
                  aic = res$aic)
  
  mbSetObj$analSet$alpha_glm_res <- glm_res
  
  if(!.on.public.web){
    print_glm_stats(glm_res)
  }
  
  return(.set.mbSetObj(mbSetObj));
}

print_glm_stats <- function(glm_res){
  
  if("crayon" %in% rownames(installed.packages())){
    
    cat(underline("MODEL INFO:"), "\n",
        italic("Number of Observations:"), " ",  glm_res$n, "\n",
        italic("Dependent Variable:"), " ", glm_res$dep_var, "\n", 
        italic("Type:"), " ",  glm_res$lmFamily$family, "\n\n", sep = "")
    
    cat(underline("MODEL FIT:"), "\n",
        italic("ANOVA X\u00B2 Deviance:"), " ",  glm_res$chisq_dev, italic(", p:"), " ", glm_res$chisq_p, "\n",
        italic("Pseudo-R\u00B2 (Deviance):"), " ", glm_res$pseudo_r2, "\n", 
        italic("Pseudo-R\u00B2 (Cragg-Uhler):"), " ", glm_res$mcf_r2, "\n", 
        italic("Pseudo-R\u00B2 (McFadden):"), " ", glm_res$cu_r2, "\n", 
        italic("AIC:"), " ",  glm_res$aic, "\n\n", sep = "")
    
    print(glm_res$table)
    
  }else{
    cat(paste("MODEL INFO:",
              paste("Number of Observations:", " ",  glm_res$n),
              paste("Dependent Variable:", " ", glm_res$dep_var), 
              paste("Type:", " ",  glm_res$lmFamily$family), sep = "\n"))
    
    cat(paste("MODEL FIT:", 
              paste("ANOVA X\u00B2 Deviance:", " ",  glm_res$chisq_dev,", p:", " ", glm_res$chisq_p),
              paste("Pseudo-R\u00B2 (Deviance):", " ", glm_res$pseudo_r2), 
              paste("Pseudo-R\u00B2 (Cragg-Uhler):", " ", glm_res$mcf_r2), 
              paste("Pseudo-R\u00B2 (McFadden):", " ", glm_res$cu_r2), 
              paste("AIC:", " ",  glm_res$aic), sep = "\n"))
    
    print(glm_res$table)
  }
}

# Plots based on this thread
# https://stats.stackexchange.com/questions/121490/interpretation-of-plot-glm-model

# Residual 𝑟𝑖 is defined as the difference between observed and predicted values

# https://stats.stackexchange.com/questions/76226/interpreting-the-residuals-vs-fitted-values-plot-for-verifying-the-assumptions
# https://stats.stackexchange.com/questions/101274/how-to-interpret-a-qq-plot

# https://gist.github.com/atyre2/ff4e1ec24e42adda8dbd43cda99d6282
# https://github.com/yeukyul/lindia/tree/master/R
# https://rpubs.com/therimalaya/43190
# https://stats.stackexchange.com/questions/70558/diagnostic-plots-for-count-regression

PlotAlphaGLM <- function(mbSet){
  
  mbSetObj <- .get.mbSetObj(mbSet)
  
  plot(mbSet$analSet$alpha_glm)
  
  residuals_plot(mbSetObj)
  glm_qqplot(mbSetObj)
  
}

residuals_plot <- function(mbSet, dot_size = 2){
  
  model <- mbSetObj$analSet$alpha_glm
  
  library(ggplot2)
  
  p <- ggplot(model, aes(x = .fitted, y = .resid)) + geom_point(size = dot_size) +
    stat_smooth(method="loess", color = "black", size = 0.5) + 
    geom_hline(yintercept=0, col="red", linetype="dashed", size = 0.5) +
    labs(x = "Predicted Values", y = "Residuals") +
    ggtitle("Residuals vs Fitted Plot") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill="NA"))
  p
}

glm_qqplot <- function(mbSet, dot_size = 2, line_size = 1) {
  
  model <- mbSetObj$analSet$alpha_glm
  # extract residuals from lm object
  res = residuals(model)
  
  # calculate slope and interncept for qqline
  slope = (quantile(res, .75) - quantile(res, .25)) / (qnorm(.75) - qnorm(.25))
  intercept = quantile(res,.25) - slope*qnorm(.25)
  qq_line = data.frame(intercept = intercept, slope = slope,check.names=FALSE);
  
  library(ggplot2)
  # generate ggplot for qqplot
  qq_plot <- ggplot(data = model) +
    stat_qq(aes(sample = res), size = dot_size) +
    labs(x = "Theoretical Quantiles", y = "Standardized Residuals") +
    geom_abline(data = qq_line, aes(intercept = intercept, slope = slope), color = "indianred3", size = line_size) +
    ggtitle("Normal-QQ Plot") + 
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          panel.border = element_rect(fill="NA"))
  
  qq_plot
}

## Utility functions ##

#'Function to update confounders used for partial correlation
#'@description This function updates which confounders will be
#'used to calculate partial correlation.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
UpdateConfItems <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(!exists("conf.vec")){
    current.msg <<- "Cannot find the current list of available metadata!";
    return (0);
  }
  
  #double check validity of metadata
  metadata <- colnames(mbSetObj$dataSet$sample_data)
  check <- metadata[(which(metadata %in% conf.vec))]
  mbSetObj$dataSet$confs <- check
  
  current.msg <<- "Successfully updated selected confounders!";
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}
