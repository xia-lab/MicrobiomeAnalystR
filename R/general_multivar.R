#'Function to create box plots of important features
#'@description This functions plots box plots of a selected feature.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param boxplotName Character, input the name of the 
#'box plot.
#'@param feat Character, input the name of the selected 
#'feature.
#'@param format Character, by default the plot format
#'is "png".
#'@param dpi Dots per inch. Numeric, by default
#'it is set to 72.
#'@parm colorPal Character, input the name of the preferred color palette.
#'Use "default" for the RColor brewer Set1 palette, "virdis" for the viridis color palette, and
#'"dark" for the RColor brewer Dark2 palette.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import grid
#'@import gridExtra
PlotBoxDataCorr<-function(mbSetObj, boxplotName, feat, format="png", dpi=72, colorPal = "dark"){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_ggplot();
  load_grid();
  load_gridExtra();
  
  variable <-  mbSetObj$analSet$var.typecor 
  
  data <- mbSetObj$analSet$boxdatacor;
  a <- data[,feat];
  ind <- which(a=="0");
  a[ind] <- 0.1;
  data$log_feat <- log(a);
  boxplotName = paste(boxplotName,".",format, sep="");
  
  numGrps <- length(levels(data$class))
  
  if(numGrps == 2){
    width <- 325
  }else if(numGrps < 4){
    width <- 350
  }else if(numGrps < 6){
    width <- 375
  }else{
    width <- 400
  }
  
  Cairo::Cairo(file=boxplotName, width=width, height=300, type=format, bg="white", dpi=dpi);
  
  box <- ggplot(data, aes(x=data$class, y = data$log_feat, fill=as.factor(class))) + stat_boxplot(geom ='errorbar') + 
    geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + labs(y="Log-transformed Counts\n", x=paste0("\n",variable), fill=variable) +
    ggtitle(feat) + theme(plot.title = element_text(hjust=0.5, size=13, face="bold"), axis.title=element_text(size=11), legend.title=element_text(size=11), axis.text=element_text(size=10));
  #remove grid
  box <- box + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "#787878", fill=NA, size=0.5))
  
  if(colorPal == "viridis"){
    box <- box + scale_fill_viridis_d()
  }else if(colorPal == "set1"){
    box <- box + scale_fill_brewer(palette="Set1")
  }else if(colorPal == "dark"){
    box <- box + scale_fill_brewer(palette="Dark2")
  }
  
  print(box)
  dev.off();
  return(.set.mbSetObj(mbSetObj))
}

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
  
  load_ppcor()
  load_viridis()
  
  # retrieve sample info
  metadata <- data.frame(sample_data(mbSet$dataSet$proc.phyobj), check.names=F, stringsAsFactors = FALSE);
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
  
  if(class(meta.subset) == "data.frame"){
    meta.subset <- data.frame(apply(meta.subset, 2, function(x) as.numeric(as.character(x))))
    meta.subset <- as.list(meta.subset)
  }else{
    meta.subset <- as.numeric(meta.subset)
    meta.subset <- as.list(as.data.frame(meta.subset))
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
    cor.results <- data.frame(cor_pval=cor.result$p.value, cor_est=cor.result$estimate)
    # calculate pcorr
    pcor.results <- ppcor::pcor.test(otu.subset, var, meta.subset, method = alg)
    pcor.results <- cbind(cor.results, pcor.results)
  }else{
    otu.subset <- lapply(seq_len(ncol(otu.table)), function(i) otu.table[,i])
    
    # first calculate corr
    cor.results <- do.call(rbind, lapply(otu.subset, function(x){
      cor.result <- cor.test(x, var, method = alg);
      data.frame(cor_pval=cor.result$p.value, cor_est=cor.result$estimate)}))  
    row.names(cor.results) <- colnames(otu.table)
    
    # calculate pcorr
    pcor.fun <- function(feats){ ppcor::pcor.test(feats, var, meta.subset, method = alg) }
    pcor.results <- do.call("rbind", lapply(otu.subset, pcor.fun))
    pcor.results <- cbind(cor.results, pcor.results)
  }
  
  #order results by p.value
  row.names(pcor.results) <- colnames(otu.table)
  resTable <- as.data.frame(pcor.results)
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
