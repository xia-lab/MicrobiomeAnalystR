##################################################
## R script for MicrobiomeAnalyst
## Description: comparative analysis methods for microbiome data
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#######################################
########### Random Forest #############
#######################################

#'Performs Random Forest Analysis
#'@description This functions performs the Random Forest (RF) analysis.
#'@param mbSetObj Input the name of the mbSetObj. By default,
#'the name should be mbSet.
#'@param treeNum Numeric, input the number of trees to grow.
#'@param tryNum Numeric, input the number of predictors to try.
#'@param randomOn Randomness setting. 1 is on, 0 is off.
#'@param variable Character, input the experimental factor for classification.
#'@param taxrank Character, input the taxonomic level to perform
#'classification. For instance, "OTU-level" to use OTUs.
#'@param datatype Character, "16S" if performing RF on 16S rRNA
#'marker gene data and "metageno" if performing RF on shotgun metagenomic data.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import randomForest
RF.Anal <- function(mbSetObj, treeNum, tryNum, randomOn, variable, taxrank){

  load_randomforest();
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # set up random numbers
  if(is.null(mbSetObj$analSet$random.seeds)){
    mbSetObj$analSet$random.seeds <- GetRandomNumbers();
    mbSetObj$analSet$cur.inx <- 0;
    mbSetObj$analSet$rn.seed <- mbSetObj$analSet$random.seeds[1];
  }
  
  if(randomOn == -1){
    rn.sd <- 123456;
  }else if(randomOn == 0){ # keep current
    rn.sd <- mbSetObj$analSet$rn.seed;
  }else{ # random on
    cur.inx <- mbSetObj$analSet$cur.inx + 1;
    rn.sd <- mbSetObj$analSet$random.seeds[cur.inx];
    mbSetObj$analSet$cur.inx <- cur.inx;
  }
  
  set.seed(rn.sd);
  # save the seed
  mbSetObj$analSet$rn.seed <- rn.sd;
  
  if(mbSetObj$module.type=="sdp"){
    taxrank<-"OTU";
    data <- mbSetObj$dataSet$norm.phyobj;
    data1 <- as.data.frame(t(otu_table(data)),check.names=FALSE);
  }else{
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    }
    
    if(taxrank=="OTU"){
      data1 <- t(phyloseq_objs$count_tables$OTU)
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- t(phyloseq_objs$count_tables[[taxrank.inx]])
    } 
  }

  meta.info <- sample_data(mbSetObj$dataSet$norm.phyobj)
  if(length(meta.vec.rf) == 0){
    sel.meta.vecs <- "NA"
  }else{
    sel.meta.vecs <- meta.info[, meta.vec.rf]
  }

  data1 <- data1[match(rownames(meta.info), rownames(data1)), ] # make sure the metadata and metabolites have the same row order
  if(length(meta.vec.rf) >0) {
    data1 <- cbind(sel.meta.vecs, data1);
    colnames(data1)[1:length(meta.vec.rf)] <- meta.vec.rf;
  }
  
  data.impfeat <<- data1;
  cls <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  variable <<- variable;
  rf_out <- randomForest(data1,cls, ntree = treeNum, mtry = tryNum, importance = TRUE, proximity = TRUE);
  # set up named sig table for display
  impmat <- rf_out$importance;
  impmat <- impmat[rev(order(impmat[,"MeanDecreaseAccuracy"])),]
  sigmat <- impmat[,"MeanDecreaseAccuracy", drop=F];
  sigmat <- signif(sigmat, 5);
  fast.write(sigmat,file="randomforests_sigfeatures.csv");
  
  mbSetObj$analSet$cls <- cls;
  mbSetObj$analSet$rf <- rf_out;
  mbSetObj$analSet$rf.sigmat <- sigmat;
  return(.set.mbSetObj(mbSetObj))
}

#'Plot Random Forest Classification
#'@description This functions plots the classification of samples
#'from the Random Forest analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param feature Numeric, input the number of important features. 
#'This is more for RF.VIP to make the sizes of the plot consistent.
#'@param imgName Character, input the name of the plot.
#'@param format Character, by default the plot is .png format.
#'@param dpi The dots per inch. Numeric, by default it is set to 72.
#'@param width Width of the plot. Numeric, by default it is set to NA.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotRF.Classify<-function(mbSetObj, feature, imgName, format="png", dpi=72, width=NA){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  imgName = paste(imgName,".", format, sep="");
  
  mbSetObj$imgSet$rf.cls <- imgName;

  if (feature < 5) {
    h <- feature * 1.2
  } else if (feature < 10) {
    h <- feature * 1.4
  } else if (feature < 15) {
    h <- feature / 1.6
  } else if (feature < 20) {
    h <- feature / 1.8
  } else if (feature < 25) {
    h <- feature / 2
  } else if (feature < 30) {
    h <- feature / 2.2
  } else if (feature < 40) {
    h <- feature / 2.5
  } else {
    h <- feature / 10
  }
  
  if(is.na(width)){
    w <- 9;
  }else if(width == 0){
    w <- 8;
  }else{
    w <- width;
  }

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(mar=c(4,4,3,2));
  cols <- grDevices::rainbow(length(levels(mbSetObj$analSet$cls))+1);
  plot(mbSetObj$analSet$rf, main="Random Forest Classification", col=cols);
  legend("topright", legend = c("Overall", levels(mbSetObj$analSet$cls)), lty=2, lwd=1, col=cols);
  dev.off();
  return(.set.mbSetObj(mbSetObj))
}

#'Plot variable importance ranked by MeanDecreaseAccuracy
#'@description This functions plot variable importance ranked by MeanDecreaseAccuracy
#'from the random forest analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param feature Numeric, input the number of important features. 
#'This is more for RF.VIP to make the sizes of the plot consistent.
#'@param imgName Character, input the name of the plot.
#'@param format Character, by default the plot is .png format.
#'@param dpi The dots per inch. Numeric, by default it is set to 72.
#'@param width Width of the plot. Numeric, by default it is set to NA.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotRF.VIP<-function(mbSetObj, feature, imgName, format="png", dpi=72, width=NA){

  load_phyloseq();
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  imgName = paste(imgName, ".", format, sep="");
  mbSetObj$imgSet$rf.imp <- imgName;

  if(is.na(width)){
    w <- 9;
  }else if(width == 0){
    w <- 8;
  }else{
    w <- width;
  }
  h <- w-1; # margin is big

  vip.score <- rev(sort(mbSetObj$analSet$rf$importance[,"MeanDecreaseAccuracy"]));
  cls <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  cls.len <- length(levels(cls));

  # re-adjust width based on group size
  rt.mrg <- cls.len + 3;
  w <- w + rt.mrg/72; # convert to inch

  Cairo::Cairo(file = imgName,  unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  op <- par(mar=c(5,7,3,rt.mrg)); # update right side margin with the number of class
  PlotImpVar(mbSetObj, vip.score,"MeanDecreaseAccuracy", feature);
  par(op);
  dev.off();
  return(.set.mbSetObj(mbSetObj))
}

#'Helper function to plot variable importance (for RF and LEfSe)
#'@description This functions plots the variable importance 
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
PlotImpVar <- function(mbSetObj, imp.vec, xlbl, feature, color.BW=FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  cls <- as.factor(mbSetObj$analSet$cls);
  cls.lbl <- levels(cls);
  cls.len <- length(levels(cls));
  
  feat.num = feature;
  if(feat.num <= 0){
    feat.num = 15;
  }
  
  if(feat.num > length(imp.vec)){
    feat.num <- length(imp.vec);
  }
  
  # first get the top subset
  imp.vec <- rev(sort(imp.vec))[1:feat.num];
  
  # reverser the order for display
  imp.vec <- sort(imp.vec);
  
  # as data should already be normalized, use mean/median should be the same
  # mns is a list contains means of all vars at each level
  # conver the list into a matrix with each row contains var averages across different lvls
  data1 <- data.impfeat;
  mns <- by(data1[, names(imp.vec)], cls,
            function(x){ # inner function note, by send a subset of dataframe
              apply(x, 2, mean, trim=0.1)
            });
  mns <- t(matrix(unlist(mns), ncol=feat.num, byrow=TRUE));
  
  vip.nms <- names(imp.vec);
  names(imp.vec) <- NULL;
  
  # modified for B/W color
  dotcolor <- ifelse(color.BW, "darkgrey", "#585855");
  dotchart(imp.vec, bg=dotcolor, xlab= xlbl, cex=1.35);
  
  mtext(side=2, at=1:feat.num, vip.nms, las=2, line=1, cex=1.1)
  
  axis.lims <- par("usr"); # x1, x2, y1 ,y2
  
  # get character width
  shift <- 2*par("cxy")[1];
  lgd.x <- axis.lims[2] + shift;
  
  x <- rep(lgd.x, feat.num);
  y <- 1:feat.num;
  par(xpd=T);
  
  load_rcolorbrewer();
  
  nc <- ncol(mns);
  
  # modified for B/W color
  colorpalette <- ifelse(color.BW, "Greys", "RdYlBu");
  col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, colorpalette))(nc); # set colors for each class
  if(color.BW) col <- rev(col);
  
  # calculate background
  bg <- matrix("", nrow(mns), nc);
  for (m in 1:nrow(mns)){
    bg[m,] <- (col[nc:1])[rank(mns[m,])];
  }
 
  
  for (n in 1:ncol(mns)){
    points(x,y, bty="n", pch=22, bg=bg[,n], cex=3);
    # now add label
    text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5), cex=1.1);
    # shift x, note, this is good for current size
    x <- x + shift/1.25;
  }
  
  # now add color key, padding with more intermediate colors for contiuous band
  col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(25, colorpalette))(50)
  
  if(color.BW) col <- rev(col);
  
  nc <- length(col);
  x <- rep(x[1] + shift, nc);
  
  shifty <- (axis.lims[4]-axis.lims[3])/3;
  starty <- axis.lims[3] + shifty;
  endy <- axis.lims[3] + 2*shifty;
  y <- seq(from = starty, to = endy, length = nc);
  
  points(x,y, bty="n", pch=15, col=rev(col), cex=2);
  
  text(x[1], endy+shifty/8, "High", cex=1.1);
  text(x[1], starty-shifty/8, "Low", cex=1.1);

}

#######################################
###########Univariate Analysis ########
#######################################

#'Main function to perform classical univariate analysis.
#'@description This functions performs classical univariate analysis
#' on the normalized microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param variable Character, input the name of the experimental factor.
#'@param p.lvl Numeric, input the adjusted p-value cutoff.
#'@param datatype Character, input whether the data is marker gene
#'data ("16S") or metagenomic data ("metageno").
#'@param shotgunid If 16S, it is set to "NA".
#'@param taxrank Character, input the taxonomic level to perform
#'univariate analysis.
#'@param statOpt Character, input "nonpar" to use non-paramentric tests including
#'Mann-Whitney for two-group comparisons or Kruskall-Wallis for multiple group comparisons. 
#'Input "tt" to use parametric tests including T-tests for two-group comparisons and ANOVA
#'for multiple group comparisons.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import MASS

PerformUnivarTest <- function(mbSetObj=NA, variable, p.lvl=0.05, shotgunid=NA, taxrank, statOpt, fc.thresh=0){
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  #rather than whole name from taxonomy just last name.
  cls <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  lvl <- length(levels(cls));
  
  if(mbSetObj$module.type=="sdp"){
    taxrank<-"OTU";
    data <- mbSetObj$dataSet$norm.phyobj;
    data1 <- as.data.frame(t(otu_table(data)),check.names=FALSE);
  }else{
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    }
    
    if(taxrank=="OTU"){
      data1 <- t(phyloseq_objs$count_tables$OTU)
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- t(phyloseq_objs$count_tables[[taxrank.inx]])
    } 
  }
  
  isNonPar <- statOpt=="nonpar"
  if(ncol(data1) > 1000 & !isNonPar){
    resTable <- PerformFastUnivTests(data1, cls);
  }else{
    resTable <- PerformUnivTests(cls, data1, nonpar=isNonPar);
  }
  
  rownames(resTable) <- colnames(data1);
  colnames(resTable) <- c("Statistics","Pvalues");
  resTable <- data.frame(resTable,check.names=FALSE);
  resTable$FDR <- p.adjust(resTable$Pvalues,method ="fdr");
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  resTable <- signif(resTable, digits = 5);
  resTable <- resTable[complete.cases(resTable), ];
  resTable <- resTable[,c(2,3,1)];
  fast.write(resTable, file="univar_test_output.csv");
  
  #getting log2fc from edgeR
  if(!is.null(mbSetObj$analSet$rnaseq$resTable.edger.all)){
  edger_df <- mbSetObj$analSet$rnaseq$resTable.edger.all
  resTable_logFC <- edger_df[match(rownames(resTable), rownames(edger_df)), 'log2FC']
  resTable_logCPM <- edger_df[match(rownames(resTable), rownames(edger_df)), 'logCPM']
  resTable <- data.frame(log2FC = resTable_logFC, resTable,logCPM=resTable_logCPM)
  }

  sigHits <- (resTable$FDR < p.lvl & abs(resTable$log2FC) > fc.thresh);
  resTable <- resTable[order(-sigHits, resTable$Pvalues), ]
  de.Num <- sum(sigHits);
  
  if(de.Num == 0){
    current.msg <<- "No significant features were identified using the given p value cutoff.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features were identified!");
  }
  
  if(nrow(resTable) > 500){
    resTable <- resTable[1:500, ];
  }
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  sigfeat <- rownames(resTable);
  
  #for individual box plot, filtered data instead of normalizated data will be used.
  #following process is just for box plot.
  ##############################################
  taxrank_boxplot <- taxrank;
  claslbl_boxplot <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  # build phyloslim obj in fly
  filt.dataphy <- mbSetObj$dataSet$filt.data;
  filt.dataphy <- apply(filt.dataphy, 2, as.integer);
  filt.dataphy <- otu_table(filt.dataphy, taxa_are_rows =TRUE);
  sample_table_boxplot <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  filt.dataphy <- merge_phyloseq(filt.dataphy, sample_table_boxplot);
  taxa_names(filt.dataphy) <- rownames(mbSetObj$dataSet$filt.data);
  data_boxplot <- filt.dataphy;
  
  if(mbSetObj$module.type=="mdp"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data_boxplot <- merge_phyloseq(data_boxplot, mbSetObj$dataSet$taxa_table);
  }
  
  #using by default names for shotgun data
  if(mbSetObj$module.type=="sdp"){
    taxrank_boxplot<-"OTU";
  }
  
  if(taxrank_boxplot!="OTU"){
    #merging at taxonomy levels
    data_boxplot <- fast_tax_glom_mem(data_boxplot, taxrank_boxplot);
    if(is.null(data_boxplot)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
    
    nm_boxplot <- tax_table(data_boxplot)[,taxrank_boxplot];
    nm_boxplot <- as.character(nm_boxplot);
    #converting NA values to unassigned
    nm_boxplot[is.na(nm_boxplot)] <- "Not_Assigned";
    data1_boxplot <- as.matrix(otu_table(data_boxplot));
    rownames(data1_boxplot) <- nm_boxplot;
    
    #all NA club together
    data1_boxplot <- as.matrix(t(sapply(by(data1_boxplot, rownames(data1_boxplot), colSums), identity)));
    data1_boxplot <- otu_table(data1_boxplot,taxa_are_rows=T);
    data_boxplot <- merge_phyloseq(data1_boxplot, sample_data(data_boxplot));
  }
  
  nm_boxplot <- taxa_names(data_boxplot);
  dat3t_boxplot <- as.data.frame(t(otu_table(data_boxplot)),check.names=FALSE);
  colnames(dat3t_boxplot) <- nm_boxplot; 


  #individual boxplot for features
  box_data <- data.frame(dat3t_boxplot[,sigfeat %in% nm_boxplot],check.names=FALSE);
  box_data$class <- claslbl_boxplot;
  mbSetObj$analSet$boxdata <- box_data;
  fast.write(t(box_data), "uni_abund_data.csv")

  mbSetObj$analSet$anal.type <- "tt";
  mbSetObj$analSet$var.type <- variable;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$univar$resTable <- mbSetObj$analSet$resTable <- resTable;
  mbSetObj$analSet$univar.taxalvl <- taxrank;

  # record parameters
  mbSetObj$paramSet$univar <- list(
        exp.factor = variable,
        anal.type = "t-tests",
        method = statOpt,
        taxalvl = taxrank,
        p.lvl = p.lvl,
        fc.thresh= fc.thresh
    );


  return(.set.mbSetObj(mbSetObj));
}

#######################################
###########MetagenomeSeq ##############
#######################################

#'Main function to perform metagenome seq analysis
#'@description This functions performs metagenome seq analysis on the microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param variable Character, input the name of the experimental factor.
#'@param p.lvl Numeric, input the adjusted p-value cutoff.
#'@param datatype Character, input whether the data is marker gene
#'data ("16S") or metagenomic data ("metageno").
#'@param shotgunid Only valid for SDP module, set to "NA".
#'@param taxrank Character, input the taxonomic level to perform
#'univariate analysis.
#'@param model Character, input the name of the statistical
#'model to fit the data. Use "zigfit" for zero-inflated Gaussian fit
#'and "ffm" for the fitFeatureModel.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import metagenomeSeq

PerformMetagenomeSeqAnal<-function(mbSetObj, variable, p.lvl, shotgunid, taxrank, model, fc.thresh=0){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_metagenomeseq();
  current.msg<<-"null"
  filt.dataphy <- mbSetObj$dataSet$filt.data;
  filt.dataphy <- apply(filt.dataphy,2,as.integer);
  filt.dataphy <- otu_table(filt.dataphy,taxa_are_rows =TRUE);
  sample_table <- sample_data(mbSetObj$dataSet$norm.phyobj, errorIfNULL=TRUE);
  filt.dataphy <- merge_phyloseq(filt.dataphy, sample_table);
  taxa_names(filt.dataphy) <- rownames(mbSetObj$dataSet$filt.data);
  data <- filt.dataphy;
  
  #data<-dataSet$norm.phyobj;
  cls <- as.factor(sample_data(data)[[variable]]);
  lvl <- length(levels(cls));

  if(mbSetObj$module.type=="mdp"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
  }else{ #using by default names for shotgun data
    taxrank <- "OTU";
  }
  if(taxrank!="OTU"){
    #merging at taxonomy levels
    data <- fast_tax_glom_mem(data, taxrank);
    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
    nm <- as.character(tax_table(data)[,taxrank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    #all NA club together
    data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    data1 <- otu_table(data1, taxa_are_rows=T);
    data <- merge_phyloseq(data1, sample_data(data));
    nm <- taxa_names(data);
  }
  tree_data <<- data;
  data <- phyloseq_to_metagenomeSeq(data);
  data <- cumNorm(data, p=cumNormStat(data));
  mod <- model.matrix(~phenoData(data)@data[,variable]);

  if(model=="zigfit"){
    tryCatch(
      {
        fit <- fitZig(data, mod);  
      }, warning = function(w){ print() },
      error = function(e) {
        current.msg <<- paste( "fitZig model failed to fit to your data! Consider a different model or further filtering your dataset!");
      }, finally = {
        if(!exists("fit")){
          return(0);
        }
        fit <- fit
      })
  }else{
    if(length(levels(cls)) > 2){
      current.msg <<- paste( "More than two groups present in your experimental factor. This model can only be used with two groups.");
      return(0);
    }else{
      tryCatch(
        {
          fit <-fitFeatureModel(data, mod);
        }, warning = function(w){ print() },
        error = function(e) {
          current.msg <<- paste( "fitFeatureModel model failed to fit to your data! Consider a different model or further filtering your dataset!");
        }, finally = {
          if(!exists("fit")){
            return(0);
          }
          fit <- fit
        })
    }
  }

  x <- MRfulltable(fit, number = nrow(assayData(data)$counts));
  x <- x[!is.na(rownames(x)), ];
  
  rownames(x) <- gsub(":1", "", x = rownames(x), fixed = TRUE);
  x$OTUnames <- as.character(rownames(x))
  
  if (!is.null(tax_table(data, errorIfNULL = FALSE))) {
    #Attach the bacterial taxonomy to the table, if available
    TAX = data.frame(tax_table(data),check.names=FALSE);
    TAX$OTUnames <- as.character(rownames(TAX));
    res = merge(x, TAX, by = "OTUnames")
  } else {
    res = x;
  }
  
  if(model=="ffm"){
    resTable <- res[,c("pvalues","adjPvalues","logFC")];
    resTable <- signif(resTable[,c(3,1,2)], digits = 5);
    colnames(resTable) <- c("logFC","Pvalues","FDR");
  }else{
    resTable <-res[,c("pvalues","adjPvalues")];
    resTable <- signif(resTable[,1:2], digits = 5);
    colnames(resTable) <- c("Pvalues","FDR");
  }

  #getting log2fc from edgeR
  if(!is.null(mbSetObj$analSet$rnaseq$resTable.edger.all)){
  edger_df <- mbSetObj$analSet$rnaseq$resTable.edger.all
  resTable_logFC <- edger_df[match(rownames(resTable), rownames(edger_df)), 'log2FC']
  resTable_logCPM <- edger_df[match(rownames(resTable), rownames(edger_df)), 'logCPM']
  resTable <- data.frame(log2FC = resTable_logFC, resTable,logCPM=resTable_logCPM)
  }

  sigHits <- (resTable$FDR < p.lvl & abs(resTable$log2FC) > fc.thresh);
  resTable <- resTable[order(-sigHits, resTable$Pvalues), , drop=FALSE]
  de.Num <- length(which(sigHits));
  print(de.Num)
  if(de.Num == 0){
    current.msg <<- paste( "No significant features were identified using the given p value cutoff. Please change the cutoff limit.");
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features were identified!")
  }

  fast.write(resTable, file="metageno_de_output.csv");
  
  if(nrow(resTable) > 500){
    resTable <- resTable[1:500, ];
  }
  
  mbSetObj$analSet$metagenoseq$resTable <- mbSetObj$analSet$resTable <- data.frame(resTable,check.names=FALSE);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  sigfeat <- rownames(resTable);
  
  #prepare individual boxplot
  box_data <- MRcounts(data);
  #subset only diff. Abundant features
  box_data <- box_data[sigfeat, ];
  
  #samples in rows
  box_data <- t(box_data);
  box_data <- data.frame(box_data,check.names=FALSE);
  colnames(box_data) <- sigfeat;
  
  claslbl <- pData(data)[ ,variable];
  box_data$class <- unlist(claslbl);

  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$anal.type <- "metagseq";
  mbSetObj$analSet$var.type <- variable;
  mbSetObj$analSet$metageno.taxalvl <- taxrank;
  mbSetObj$analSet$id.type <- shotgunid;
  # record parameters
  mbSetObj$paramSet$metagenoseq <- list(
        exp.factor = variable,
        anal.type = "metagseq",
        method = model,
        taxalvl = taxrank,
        p.lvl = p.lvl,
        fc.thresh = fc.thresh
    );

  
  return(.set.mbSetObj(mbSetObj));
}

####################################
##############LEfSe ################
####################################

#'Main function to perform LEfSe analysis
#'@description This functions performs LEfSe analysis on the microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param p.lvl Numeric, input the adjusted p-value cutoff.
#'@param lda.lvl Numeric, input the Log LDA score cutoff.
#'@param variable Character, input the name of the experimental factor.
#'@param isfunc Logical, default set to "F".
#'@param datatype Character, input whether the data is marker gene
#'data ("16S") or metagenomic data ("metageno").
#'@param shotgunid Only valid for SDP module, set to "NA".
#'@param taxrank Character, input the taxonomic level to perform
#'univariate analysis.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import MASS

PerformLefseAnal <- function(mbSetObj, p.lvl, pvalOpt="fdr", lda.lvl, variable, isfunc, shotgunid, taxrank){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_MASS();
 
  claslbl <<- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);

  if(mbSetObj$module.type=="sdp"){
    taxrank<-"OTU";
    data <- mbSetObj$dataSet$norm.phyobj;
    dat3t <- as.data.frame(t(otu_table(data)),check.names=FALSE);
  }else{
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    }
    if(taxrank=="OTU"){
      dat3t <- t(phyloseq_objs$count_tables$OTU)
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      dat3t <- t(phyloseq_objs$count_tables[[taxrank.inx]])
    } 
  }
  data.impfeat_lefse <<- dat3t;
  
  set.seed(56290);
  #KW rank sum test
  rawpvalues <- apply(dat3t, 2, function(x) kruskal.test(x, claslbl)$p.value);
  ord.inx <- order(rawpvalues);
  rawpvalues <- rawpvalues[ord.inx];
  clapvalues <- p.adjust(rawpvalues, method ="fdr");
  dat3t <- dat3t[,ord.inx];
  
  if(length(rawpvalues) > 500){
    rawpvalues <- rawpvalues[1:500];
    clapvalues <- clapvalues[1:500];
    dat3t <- dat3t[,1:500];
  };
  
  wil_datadf <- as.data.frame(dat3t,check.names=FALSE);
  
  #if no subclass within classes then no wilcoxon rank sum test  
  #Linear Discriminant analysis (LDA)
  ldares <- lda(claslbl~ .,data = wil_datadf);
  ldamean <- as.data.frame(t(ldares$means),check.names=FALSE);
  class_no <<- length(unique(claslbl));
  if(class_no==2){
  ldamean$LDAscore <- sign(ldamean[,1]-ldamean[,2])*signif(log10(1+abs(ldamean[,1]-ldamean[,2])/2),digits=3);
  }else{
   ldamean$max <- apply(ldamean[,1:class_no],1,max);
  ldamean$min <- apply(ldamean[,1:class_no],1,min);
  ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);
  }
  ldamean$Pvalues <- signif(rawpvalues,digits=5);
  ldamean$FDR <- signif(clapvalues,digits=5);
  resTable <- ldamean;

  # it seems lda add ` around names containing dash "-", need to strip this off
  rawNms <- rownames(resTable);
  rownames(resTable) <- gsub("`", '', rawNms);
  if(pvalOpt == "raw"){
    de.Num <- sum(rawpvalues<=p.lvl & abs(ldamean$LDAscore)>=lda.lvl)
  }else{
    de.Num <- sum(clapvalues<=p.lvl & abs(ldamean$LDAscore)>=lda.lvl)
  }
  
  if(de.Num == 0){
    current.msg <<- "No significant features were identified with given criteria.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
  }
  
  # sort by p value
  ord.inx <- order(resTable$Pvalues, rev(abs(resTable$LDAscore)));
  resTable <- resTable[ord.inx, ,drop=FALSE];
  #p-values column to appear first; then FDR and then others
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];
  
  # need to control digits
  resTable <- signif(resTable, 5);
  
  #only getting the names of DE features
  if(pvalOpt == "raw"){
    diff_ft <<- rownames(resTable)[(abs(resTable$LDAscore) > lda.lvl & resTable$Pvalues < p.lvl)]
  }else{
    diff_ft <<- rownames(resTable)[(abs(resTable$LDAscore) > lda.lvl & resTable$FDR < p.lvl)]
  }

  resTable$max <- resTable$min <- NULL;
  #if only two groups are present in sample variable
  cls.lbls <- levels(claslbl);
  grp.dat <- resTable[, cls.lbls];
  res.cls <- as.character(apply(grp.dat, 1, function(x){cls.lbls[which.max(x)]}));

  fast.write(resTable, file="lefse_de_output.csv");
  mbSetObj$analSet$lefse$resTable <- mbSetObj$analSet$resTable <- resTable;
  
  #subset dataset for bar plot visualization (LDA Score)
  ldabar <- as.data.frame(rownames(resTable),check.names=FALSE);
  ldabar[,2] <- resTable$LDAscore;
  ldabar[,3] <- res.cls;
  if(pvalOpt == "raw"){
    ldabar[,4] <- resTable$Pvalues;
    ldabar.sub <- subset(ldabar, (abs(ldabar[,2]) > lda.lvl & ldabar[,4] < p.lvl));
  }else{
    ldabar[,4] <- resTable$FDR;
    ldabar.sub <- subset(ldabar, (abs(ldabar[,2]) > lda.lvl & ldabar[,4] < p.lvl));
  }
  
  #visualizing top features based on LDA score
  ldabar <<- ldabar.sub[,-4];
  
  #preparing data for indvidual box plot
  sigfeat <<- rownames(resTable);
  taxrank <<- taxrank;
  box_data <- as.data.frame(wil_datadf[, sigfeat],check.names=FALSE);
  colnames(box_data) <- sigfeat;
  box_data$class <- claslbl;
  
  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$anal.type <- "lefse";
  mbSetObj$analSet$lefse.taxalvl <- taxrank;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$meta <- variable;
  
  # record parameters
  
  mbSetObj$paramSet$lefse <- list(
        taxalvl = taxrank,
        p.val = p.lvl,
        p.type = pvalOpt,
        lda.lvl = lda.lvl,
        factor = variable
   );

  return(.set.mbSetObj(mbSetObj));
}

#'Plot LEfSe summary
#'@description This functions graphically summarizes the LEfSe results
#'in a VIP plot.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param ldaFeature Numeric, input the number of top
#'features to include in the plot.
#'@param layoutOptlf Character, input "dot" to create the VIP plot 
#'and "bar" to create the bar graph.
#'@param imgName Character, input the name
#'of the plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param width Numeric, input the width of the plot. By
#'default it is set to NA.
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@param colOpt Character, "default", "viridis",
#'"cividis", or "plasma".
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotLEfSeSummary <- function(mbSetObj, ldaFeature, layoutOptlf, imgName, format="png", width = NA, dpi=72, colOpt="default") {

  load_phyloseq();
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  set.seed(280561493);
  
  imgName = paste(imgName, ".", format, sep="");
  mbSetObj$analSet$lefse_plot <- imgName;

  ldabar <- ldabar;
  
  if(nrow(ldabar) == 0){
    Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=7, height=7, type=format, bg="white");
    empty <- ggplot() + theme_void()
    print(empty)
    dev.off();
    return(.set.mbSetObj(mbSetObj))
  }
  
  ldabar <-  ldabar[order(-abs(ldabar[[2]])), ];
  
  if(ldaFeature < nrow(ldabar)) {
    ldabar <- ldabar[1:ldaFeature,];
  };
  
  vip.score <- ldabar[[2]];
  names(vip.score) <- ldabar[[1]];
  
    sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
    meta <- mbSetObj$analSet$meta
    cls.len <- length(levels(as.factor(sample_table[[meta]])));

  if(is.na(width)){
    
    if(cls.len <= 4){
      w <- 9.25;
    }else if(cls.len < 6){
      w <- 10.25;
    }else if(cls.len < 8){
      w <- 11.25;
    }else if(cls.len < 10){
      w <- 12.25;
    }else if(cls.len < 12){
      w <- 13.25;
    }else{
      w <- 12;
    }
    
    if(length(vip.score) < 5 ){
      h <- 6;
    } else if (length(vip.score) < 10){
      h <- length(vip.score)/1.3;
    } else if (length(vip.score) < 15){
      h <- length(vip.score)/1.6;
    } else if (length(vip.score) < 20){
      h <- length(vip.score)/1.8;
    } else if (length(vip.score) < 25){
      h <- length(vip.score)/2;
    } else if (length(vip.score) < 30){
      h <- length(vip.score)/2.2;
    } else if (length(vip.score) < 40){
      h <- length(vip.score)/2.5;
    } else {
      h <- length(vip.score)/4.5;
    };
  }else if(width == 0){
    w <- 8;
  }else{
    w <- width;
  }
  
  # re-adjust width based on group size to accomodate miniheatmap
  if(layoutOptlf == "dot") {
       rt.mrg <- cls.len + 3;
       w <- w + rt.mrg/72; # convert to inch
       Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
       op <- par(mar=c(5,7,3,rt.mrg)); # update right side margin with the number of class
   }else{
       Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
   }
   PlotImpVarLEfSe(mbSetObj, vip.score, layoutOptlf, mbSetObj$analSet$meta, colOpt);
   dev.off();
   return(.set.mbSetObj(mbSetObj))
}

#'Plot LEfSe summary
#'@description This functions graphically summarizes the LEfSe results
#'using a bargraph.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotImpVarLEfSe <- function(mbSetObj, imp.vec, layoutOptlf, meta, colOpt="default", color.BW=FALSE){

  load_phyloseq();
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(layoutOptlf == "dot") {
    
    sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
    
    feat.num <- length(imp.vec);
    
    # first get the top subset
    imp.vec <- rev(sort(imp.vec))[1:feat.num];
    
    # reverser the order for display
    imp.vec <- sort(imp.vec);
    
    # as data should already be normalized, use mean/median should be the same
    # mns is a list contains means of all vars at each level
    # conver the list into a matrix with each row contains var averages across different lvls
    data1 <- data.impfeat_lefse;
    mns <- by(data1[, names(imp.vec)], 
              sample_table[[meta]],
              function(x){ # inner function note, by send a subset of dataframe
                apply(x, 2, mean, trim=0.1)
              });
    mns <- t(matrix(unlist(mns), ncol=feat.num, byrow=TRUE));
    
    vip.nms <- names(imp.vec);
    vip.clean.nms <- CleanTaxaNames(mbSetObj, vip.nms)
    names(imp.vec) <- NULL;
    
    # modified for B/W color
    dotcolor <- ifelse(color.BW, "darkgrey", "#585855");
    dotchart(imp.vec, bg=dotcolor, xlab= "LDA score", cex=1.35);
    
    mtext(side=2, at=1:feat.num, vip.clean.nms, las=2, line=1, cex=1.1)
    
    axis.lims <- par("usr"); # x1, x2, y1 ,y2
    
    # get character width
    shift <- 2*par("cxy")[1];
    lgd.x <- axis.lims[2] + shift;
    
    x <- rep(lgd.x, feat.num);
    y <- 1:feat.num;
    par(xpd=T);
    
    load_rcolorbrewer();
    
    nc <- ncol(mns);
    
    # modified for B/W color
    colorpalette <- ifelse(color.BW, "Greys", "RdYlBu");
    col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, colorpalette))(nc); # set colors for each class
    if(color.BW) col <- rev(col);
    
    # calculate background
    bg <- matrix("", nrow(mns), nc);
    for (m in 1:nrow(mns)){
      bg[m,] <- (col[nc:1])[rank(mns[m,])];
    }
    
    cls.lbl <- levels(as.factor(sample_table[[meta]]));
    
    for (n in 1:ncol(mns)){
      points(x,y, bty="n", pch=22, bg=bg[,n], cex=3);
      # now add label
      text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5), cex=1.1);
      # shift x, note, this is good for current size
      x <- x + shift/1.25;
    }
    
    # now add color key, padding with more intermediate colors for contiuous band
    col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(25, colorpalette))(50)
    
    if(color.BW) col <- rev(col);
    
    nc <- length(col);
    x <- rep(x[1] + shift, nc);
    
    shifty <- (axis.lims[4]-axis.lims[3])/3;
    starty <- axis.lims[3] + shifty;
    endy <- axis.lims[3] + 2*shifty;
    y <- seq(from = starty, to = endy, length = nc);
    
    points(x,y, bty="n", pch=15, col=rev(col), cex=2);
    
    text(x[1], endy+shifty/8, "High", cex=1.1);
    text(x[1], starty-shifty/8, "Low", cex=1.1);
    
    #par(op);
    
  } else {
    
    ldabar <-  ldabar[order(-abs(ldabar[[2]])), ];
    feat.num <- length(imp.vec);
    if(feat.num > nrow(ldabar)){
      ldabar <- ldabar;
    }else{
      ldabar <- ldabar[1:feat.num, ];
    }
    
    # trim levels for ASV
    levels(ldabar[,1]) <- strtrim(levels(ldabar[,1]), 15)
    box<- ggplot(ldabar, aes(x=reorder(ldabar[,1],ldabar[,2]), y=ldabar[,2], fill=ldabar[,3]))+ 
      geom_bar(stat="identity",width=0.8) + coord_flip() + labs(y="LDA score", x="Features", fill="Class") + 
      theme_bw() + scale_fill_brewer(palette="Set1") + 
      theme(axis.title=element_text(size=14), axis.text=element_text(size=12), legend.text = element_text(size=12));
    print(box);
  }
}



#######################################
###########EdgeR/DESeq2################
#######################################

#'Main function to perform RNAseq analysis
#'@description This functions performs RNAseq analysis on the microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param opts Character, input "EdgeR" to use the edgeR algorithm and
#'"DESeq2" to use the DESeq2 algorithm.
#'@param p.lvl Numeric, input the adjusted p-value cutoff.
#'@param variable Character, input the experimental factor.
#'@param datatype Character, input "16S" if the data is marker gene
#'data and "metageno" if it is metagenomic data.
#'@param shotgunid Only valid for SDP module, set to "NA".
#'@param taxrank Character, input the taxonomic level
#'to use for RNAseq analysis.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import DESeq2

PerformRNAseqDE<-function(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank, fc.thresh=0, comp1=1, comp2=2){
  if(opts=="DESeq2"){
    .prepare.deseq(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank, fc.thresh, comp1, comp2);
    .perform.computing();
    res = .save.deseq.res();
  }else{
    data <- .prepare_rnaseq(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank, fc.thresh, comp1, comp2)
    res <- .perform_edger(variable, data, p.lvl, fc.thresh, comp1, comp2)
  }
  return(1)
}

.prepare.deseq<-function(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank,fc.thresh=0, comp1=1, comp2=2){
  data <- .prepare_rnaseq(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank,fc.thresh, comp1, comp2)
  mbSetObj = .get.mbSetObj(mbSetObj);
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  if(length(claslbl) > 120){
    AddErrMsg("Only EdgeR is supported for sample size over 120.");
    return(0);
  }else{
    dat.in <- list(data=data, variable=variable);
    my.fun <- function(){
      # create formula based on user selection
      variable <- dat.in$variable;
      my.formula <- as.formula(paste("~", variable));
      #converting from phyloslim object to deseq
      library(phyloseq);
      library(DESeq2);
      diagdds <- phyloseq_to_deseq2(dat.in$data, my.formula);
      geoMeans <- apply(counts(diagdds), 1, 
                        function(x, na.rm=TRUE){exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))}
      );
      diagdds <- estimateSizeFactors(diagdds, geoMeans = geoMeans);
      diagdds <- DESeq(diagdds, test="Wald", fitType="parametric");

      if(comp1 %in% claslbl && comp2 %in% claslbl ){
      res <- results(diagdds, independentFiltering = FALSE, cooksCutoff = Inf, contrast=c(variable, comp1,comp2));
      }else{
      res <- results(diagdds, independentFiltering = FALSE, cooksCutoff = Inf);
      }
      # make sure it is basic R, not DESeq2 obj
      resTable <- data.frame(res[,c("log2FoldChange" ,"lfcSE","pvalue","padj")],check.names=FALSE); 
      return(resTable);
    }
    dat.in <- list(data=data, variable=variable, my.fun=my.fun);
    qs::qsave(dat.in, file="dat.in.qs");
    return(1)
  }
return(1)
}

.save.deseq.res <- function(){
  mbSetObj <-.get.mbSetObj("NA");
  dat3t<- mbSetObj$analSet$rnaseq$data.rnaseq
  
  dat.in <- qs::qread("dat.in.qs"); 
  variable = dat.in$variable
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]])
  p.lvl <- mbSetObj$analSet$rnaseq.plvl
  fc.thresh <- mbSetObj$analSet$rnaseq.fcthresh
  
  resTable <- dat.in$my.res;
  resTable <- signif(data.matrix(resTable), digits=5);
  colnames(resTable) <- c("log2FC","lfcSE","Pvalues","FDR");
  mbSetObj$analSet$anal.type <- "deseq";
  
  resTable <- as.data.frame(resTable,check.names=FALSE);
  sigHits <- resTable$FDR < p.lvl & abs(resTable$log2FC) > fc.thresh;
  de.Num <- length(which(sigHits));
  if(de.Num == 0){
    current.msg <<- "No significant features were identified using the given p value cutoff.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features were identified!");
  }
  print(current.msg);
  
  ord.inx <- order(-sigHits, resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  fast.write(resTable, file="rnaseq_de.csv");
  
  if(nrow(resTable) > 500){
    resTable<-resTable[1:500, ];
  }
  
  mbSetObj$analSet$rnaseq$resTable <- mbSetObj$analSet$resTable <- as.data.frame(resTable,check.names=FALSE);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[sigHits]
  
  #individual boxplot for features
  sigfeat <- rownames(resTable);
  box_data <- as.data.frame(dat3t[ ,sigfeat],check.names=FALSE);
  colnames(box_data) <- sigfeat;
  box_data$class <- claslbl;
  
  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;      
  tree_data <<- data;
  .set.mbSetObj(mbSetObj)
  return(1)
}

.prepare_rnaseq<-function(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank, fc.thresh=0, comp1=1, comp2=2){

  load_phyloseq();

  mbSetObj <-.get.mbSetObj(mbSetObj);

  taxrank <<- taxrank;
  # build phyloslim obj in fly

  filt.dataphy <- mbSetObj$dataSet$filt.data;
  filt.dataphy <- apply(filt.dataphy, 2, as.integer);
  filt.dataphy <- otu_table(filt.dataphy, taxa_are_rows =TRUE);
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  filt.dataphy <- merge_phyloseq(filt.dataphy, sample_table);
  taxa_names(filt.dataphy) <- rownames(mbSetObj$dataSet$filt.data);
  data <- filt.dataphy;
  
  if(mbSetObj$module.type=="mdp"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
  }else{
    data <- data;
  }
  
  #using by default names for shotgun data
  if(mbSetObj$module.type=="sdp"){
    taxrank<-"OTU";
  }
  
  if(taxrank=="OTU"){
    data <- data;
    nm <- taxa_names(data);
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_mem(data, taxrank);
    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
    nm <- as.character(tax_table(data)[,taxrank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    
    #all NA club together
    data1 <- as.matrix(t(sapply(by(data1, rownames(data1), colSums), identity)));
    data1 <- otu_table(data1,taxa_are_rows=T);
    data <- merge_phyloseq(data1, sample_data(data));
    nm <- taxa_names(data);
  }

  dat3t <- as.data.frame(t(otu_table(data)),check.names=FALSE);
  colnames(dat3t) <- nm;
  
  mbSetObj$analSet$rnaseq$data.rnaseq <- dat3t
  mbSetObj$analSet$var.type <- variable;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$rnaseq.taxalvl <- taxrank;
  mbSetObj$analSet$rnaseq.meth <- opts;
  mbSetObj$analSet$rnaseq.plvl <- p.lvl
  mbSetObj$analSet$rnaseq.fcthresh <- fc.thresh

  # record parameters
  mbSetObj$paramSet$rnaseq <- list(
        exp.factor = variable,
        anal.type = "rnaseq",
        method = opts,
        taxalvl = taxrank,
        p.lvl = p.lvl,
        fc.thresh = fc.thresh,
        comp1 = comp1,
        comp2 = comp2
    );
  .set.mbSetObj(mbSetObj)
  return(data)
}

.perform_edger <- function(variable, data, p.lvl=0.05, fc.thresh=0, comp1="", comp2="") {
  #save.image("edger.RData");
  mbSetObj <- .get.mbSetObj(mbSetObj)
  dat3t <- mbSetObj$analSet$rnaseq$data.rnaseq
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]])

  # Using filtered data, RLE normalization within it
  dge <- phyloseq_to_edgeR(data, variable)
  if(comp1 %in% claslbl && comp2 %in% claslbl ){
  et <- edgeR::exactTest(dge, c(comp1,comp2))
  }else{
  et <- edgeR::exactTest(dge)
  }
  tt <- edgeR::topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue")
  res <- tt@.Data[[1]]
  sigHits <-res$FDR < p.lvl & abs(res$logFC) > fc.thresh;
  de.Num <- sum(sigHits);
 
  if (de.Num == 0) {
    current.msg <<- "No significant features were identified using the given p value cutoff."
  } else {
    current.msg <<- paste("A total of", de.Num, "significant features were identified!")
  }

  resTable <- res[, c("logFC", "logCPM", "PValue", "FDR")]
  resTable <- signif(data.matrix(resTable), digits = 5)
  colnames(resTable) <- c("log2FC", "logCPM", "Pvalues", "FDR")
  mbSetObj$analSet$anal.type <- "edgr"

  resTable <- as.data.frame(resTable, check.names=FALSE)
  ord.inx <- order(-sigHits, resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE]
  mbSetObj$analSet$rnaseq$resTable.edger.all <- resTable

  fast.write(resTable, file="rnaseq_de.csv")

  if (nrow(resTable) > 500) {
    resTable <- resTable[1:500, ]
  }

  mbSetObj$analSet$rnaseq$resTable <- mbSetObj$analSet$resTable <- as.data.frame(resTable, check.names=FALSE)

  # Only getting the names of DE features
  diff_ft <<- rownames(resTable)[which(resTable$FDR < p.lvl & abs(resTable$log2FC) > fc.thresh)]
  taxrank <<- taxrank  

  # Individual boxplot for features
  sigfeat <- diff_ft  # Use the DE features identified
  box_data <- as.data.frame(dat3t[, sigfeat], check.names=FALSE)
  colnames(box_data) <- sigfeat
  box_data$class <- claslbl

  mbSetObj$analSet$boxdata <- box_data
  mbSetObj$analSet$sig.count <- de.Num

  tree_data <<- data
  .set.mbSetObj(mbSetObj)

  return(1)
}



########################################################
###########)Permanova_Pairwise##########################
########################################################
###adopted from ecole package https://rdrr.io/github/phytomosaic/ecole/
.permanova_pairwise <- function(x,
                                 grp,
                                 permutations = 999,
                                 method = 'bray',
                                 padj = 'fdr', ...) {
  f     <- grp
  if (!all(table(f) > 1)) warning('factor has singletons! perhaps lump them?')
  co    <- combn(unique(as.character(f)),2)
  nco   <- NCOL(co)
  out   <- data.frame(matrix(NA, nrow=nco, ncol=5))
  dimnames(out)[[2]] <- c('pairs', 'SumOfSqs', 'F.Model', 'R2', 'pval')
  if (!inherits(x, 'dist')) {
    D <- vegan::vegdist(x, method=method)
  } else {
    D <- x
  }
  #cat('Now performing', nco, 'pairwise comparisons. Percent progress:\n')
  for(j in 1:nco) {
    cat(round(j/nco*100,0),'...  ')
    ij  <- f %in% c(co[1,j],co[2,j])
    Dij <- as.dist(as.matrix(D)[ij,ij])
    fij <- data.frame(fij = f[ij])
    a   <- vegan::adonis2(Dij ~ fij, data=fij, permutations = permutations, ...);
    out[j,1] <- paste(co[1,j], 'vs', co[2,j])
    out[j,2] <- a$SumOfSqs[1]
    out[j,3] <- a$F[1]
    out[j,4] <- a$R2[1]
    out[j,5] <- a$`Pr(>F)`[1]
  }
  #cat('\n')
  out$p.adj <- p.adjust(out$pval, method=padj)
  out$SumOfSqs <-NULL
  #attr(out, 'p.adjust.method') <- padj
  #cat('\np-adjust method:', padj, '\n\n');
  return(out)
}


########################################################
###########Correlation Analysis (Pattern Hunter)########
########################################################

#'Main function to calculate correlation of all other feature to a given feature name
#'@description This functions calculate correlation of all other feature to a given feature name.
#'This is used in the Pattern Search analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dist.name Input the distance measure. Input "pearson" for 
#'Pearson R, "spearman" for Spearman rank correlation, and "kendall"
#'for Kendall rank correlation.
#'@param taxrank Character, input the taxonomic level to use.
#'@param taxa If the pattern is defined by using a specific taxon, input the name
#'of that taxa here.
#'@param variable Input the name of the experimental factor.
#'@param datatype If the data is marker gene, input "16S". If the
#'data is metagenomic, use "metageno".
#'@param shotfeat Only valid for SDP module, set to "null".
#'@param shotgunid Only valid for SDP module, set to "NA".
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

FeatureCorrelation <- function(mbSetObj, dist.name, taxrank, feat){

  load_phyloseq();
   
  mbSetObj <- .get.mbSetObj(mbSetObj);
 
  mbSetObj$dataSet$proc.phyobj@tax_table <- mbSetObj$dataSet$proc.phyobj@tax_table[,!is.na(colnames(mbSetObj$dataSet$proc.phyobj@tax_table))]
  mbSetObj$dataSet$norm.phyobj@tax_table <- mbSetObj$dataSet$norm.phyobj@tax_table[,!is.na(colnames(mbSetObj$dataSet$norm.phyobj@tax_table))]
 
  if(mbSetObj$module.type=="mdp"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
    if(taxrank=="OTU"){
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      data1 <- as.matrix(otu_table(data));
    }else{
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      #merging at taxonomy levels
      data <- fast_tax_glom_mem(data, taxrank);
      if(is.null(data)){
        AddErrMsg("Errors in projecting to the selected taxanomy level!");
        return(0);
      }
      
      nm <- as.character(tax_table(data)[,taxrank]);
      #converting NA values to unassigned
      nm[is.na(nm)] <- "Not_Assigned";
      data1 <- as.matrix(otu_table(data));
      rownames(data1) <- nm;
      #all NA club together
      data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
    }
  }else{
    data <- mbSetObj$dataSet$norm.phyobj;
    taxrank <- "OTU";
    data1 <- as.matrix(otu_table(data));
  }
  
  if(feat %in% rownames(data1)){
    feat_data <- as.numeric(data1[feat,]);
  }else{
    AddErrMsg("Cannot find the given feature name in the data!");
    return(0);
  }
  
  data1 <- t(data1);
  qs::qsave(data1, "match_data.qs")
  #making boxplot data
  
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  variable <- colnames(sample_table)[1];
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  boxdata <- as.data.frame(data1,check.names=FALSE);
  boxdata$class <- clslbl;
  mbSetObj$analSet$boxdata <- boxdata;
 
  if(dist.name == "sparcc"){
    
    permNum <- 100
    pvalCutoff <- 1 
    corrCutoff <- 0
    
    if(.on.public.web){
      sparcc_results <- RunFastSpar_mem(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network", "feat")
    }else{
      sparcc_results <- RunFastSpar(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network", "feat")
    }
    cbtempl.results <- sparcc_results[which(sparcc_results$Taxon1==feat), 2:4]
    fdr.col <- p.adjust(cbtempl.results[,3], "fdr");
    rownames(cbtempl.results) <- cbtempl.results[,1]
    
    cor.res <- signif(cbind(cbtempl.results[,-1], fdr.col),5)
    ord.inx <- order(cor.res[,2])
    sig.mat <- cor.res[ord.inx,]
    colnames(sig.mat) <- c("Correlation", "P-Value", "FDR");
  }else{
    cbtempl.results <- apply(data1, 2, template.match, feat_data, dist.name);
    cor.res <- t(cbtempl.results);
    fdr.col <- p.adjust(cor.res[,3], "fdr");
    cor.res <- cbind(cor.res, fdr.col);
    colnames(cor.res) <- c("Correlation", "T-Stat", "P-Value", "FDR");
    ord.inx <- order(cor.res[,3])
    sig.mat <-signif(cor.res[ord.inx,],5);
  } 
  fileName <- "correlation_feature.csv";
  fast.write(sig.mat,file=fileName);
  mbSetObj$analSet$resTable <- as.data.frame(sig.mat,check.names=FALSE);
  #removing Inf values from table
  is.na(mbSetObj$analSet$resTable) <- sapply(mbSetObj$analSet$resTable, is.infinite);
  mbSetObj$analSet$resTable[is.na(mbSetObj$analSet$resTable)]<-0;
  
  mbSetObj$analSet$pattern <- feat;
  mbSetObj$analSet$taxrank <- taxrank;
  mbSetObj$analSet$feat.corr <- TRUE;
  
  sig.nm <<- fileName;
  mbSetObj$analSet$cor.mat <- sig.mat;
  mbSetObj$analSet$corph.taxalvl <- taxrank;
  mbSetObj$analSet$corph.meth <- dist.name;
  mbSetObj$analSet$sig.count <- 0; # note, not a DE analysis here
  
  return(.set.mbSetObj(mbSetObj));
  
}

#'Plot Pattern Search
#'@description This functions plots a bargraph of the
#' correlation between a specific taxon to other taxa in the Pattern
#' Search analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param imgName Character, input the name
#'of the plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param width Numeric, input the width of the plot. By
#'default it is set to NA.
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotCorr <- function(mbSetObj, imgName, format="png", dpi=72,appendnm, width=NA){

  load_phyloseq();
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  variable <- mbSetObj$analSet$pattern.var;
  data <- qs::qread("match_data.qs")
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  
  if(is.null(variable)){
    variable <- colnames(sample_table)[1];
  }
  
  cor.res <- mbSetObj$analSet$cor.mat;
  pattern <- mbSetObj$analSet$pattern;
  title <- paste(mbSetObj$analSet$taxrank, "correlated with", pattern);
  
  if(nrow(cor.res) > 25){
    cor.res <- cor.res[1:26, ];
  }
  
  # first get most signficant ones (p value)
  ord.inx <- order(cor.res[,3]);
  cor.res <- cor.res[ord.inx, ];
  # then order by their direction (correlation)
  ord.inx <- order(cor.res[,1]);
  
  if(sum(cor.res[,1] > 0) == 0){ # all negative correlation
    ord.inx <- rev(ord.inx);
  }
  
  cor.res <- cor.res[ord.inx, ];
  # need to remove pattern feature from corr matrix (always 1)
  if(mbSetObj$analSet$feat.corr & mbSetObj$analSet$corph.meth != "sparcc"){
    match.inx <- match(pattern, rownames(cor.res))
    cor.res <- cor.res[-match.inx,]
  }else{
    if(nrow(cor.res) > 25){
      cor.res <- cor.res[-(nrow(cor.res)),]
    }
  }


  
  title <- paste("Top",nrow(cor.res), tolower(mbSetObj$analSet$taxrank), "correlated with", pattern);
  imgName = paste(imgName, ".", format, sep="");
  mbSetObj$imgSet$cor.ph <- imgName;
  feat.num <- nrow(cor.res);
  
  if(is.na(width)){
    h <- w <- 10
  }else if(width == 0){
    h <- w <- 10
  }else{
    w <- width;
    h <- w/0.85
  }
  
  cls.len <- length(levels(as.factor(sample_table[[variable]])));
  
  if(cls.len == 2){
    rt.mrg <- 7.75;
  }else if(cls.len == 3){
    rt.mrg <- 8.75;
    w <- w + 1
  }else if(cls.len == 4){
    rt.mrg <- 9.75;
    w <- w + 1.5
  }else if(cls.len == 5){
    rt.mrg <- 10.75;
    w <- w + 2
  }else if(cls.len == 6){
    rt.mrg <- 11.75;
    w <- w + 2.5
  }else{
    rt.mrg <- 12;
    w <- w + 3
  }
  
  # for heatmap
  mns <- by(data[, rownames(cor.res)], 
            sample_table[[variable]],
            # sample_table[["SampleType"]],
            function(x){ # inner function note, by send a subset of dataframe
              apply(x, 2, mean, trim=0.1)
            });
  mns <- t(matrix(unlist(mns), ncol=nrow(cor.res), byrow=TRUE));
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  
  op <- par(mar=c(5, 10, 3, rt.mrg));
  
  feat.nms <- substr(rownames(cor.res), 1, 20);
  feat.num <- length(feat.nms)
  rownames(cor.res) <- NULL;
  
  cols <- case_when(
    cor.res > 0.75 ~ "#d47273",
    cor.res > 0.5 & cor.res < 0.75 ~ "#dd9091",
    cor.res > 0.25 & cor.res < 0.5 ~ "#e6aeaf",
    cor.res > 0 & cor.res < 0.25 ~ "#f0cccd",
    cor.res < 0 & cor.res > -0.25 ~ "#ccdef0",
    cor.res > -0.25 ~ "#ccdef0",
    cor.res > -0.5 ~ "#aecbe6",
    cor.res > -0.75 ~ "#90b7dd",
    TRUE ~ "#72a4d4"
  )
  
  dotchart(cor.res[,1], pch="", xlim=c(-1,1), xlab="Correlation Coefficients", main=title, cex = 1.05);
  mtext(side=2, at=1:feat.num, feat.nms, las=2, line=1, cex=1.05)
  barplot(cor.res[,1], space=c(0.5, rep(0, nrow(cor.res)-1)), xlim=c(-1,1), xaxt="n", col = cols, add=T,horiz=T);
  
  axis.lims <- par("usr"); # x1, x2, y1 ,y2
  
  # get character width
  shift <- 2*par("cxy")[1];
  lgd.x <- axis.lims[2] + shift;
  
  x <- rep(lgd.x, feat.num);
  y <- 1:feat.num;
  par(xpd=T);
  
  load_rcolorbrewer();
  nc <- ncol(mns);
  
  col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(nc); # set colors for each class
  
  # calculate background
  bg <- matrix("", nrow(mns), nc);
  for (m in 1:nrow(mns)){
    bg[m,] <- (col[nc:1])[rank(mns[m,])];
  }
  
  cls.lbl <- levels(as.factor(sample_table[[variable]]));
  
  for (n in 1:ncol(mns)){
    points(x,y, bty="n", pch=22, bg=bg[,n], cex=3);
    # now add label
    text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2, 0.2), cex=1.05);
    # shift x, note, this is good for current size
    x <- x + shift/1.05;
  }
  
  # now add color key, padding with more intermediate colors for contiuous band
  col <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(25, "RdYlBu"))(50)
  
  nc <- length(col);
  x <- rep(x[1] + shift, nc);
  
  shifty <- (axis.lims[4]-axis.lims[3])/3;
  starty <- axis.lims[3] + shifty;
  endy <- axis.lims[3] + 2*shifty;
  y <- seq(from = starty, to = endy, length = nc);
  
  points(x,y, bty="n", pch=15, col=rev(col), cex=2);
  
  text(x[1], endy+shifty/8, "High", cex=1.05);
  text(x[1], starty-shifty/8, "Low", cex=1.05);
  
  par(op)
  dev.off();
  
  return(.set.mbSetObj(mbSetObj))
}

#'Match Patterns Function
#'@description This functions matches patterns (when
#'the defined pattern is set to a predefined profile or custom
#'profile) in the Pattern Search analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dist.name Input the distance measure. Input "pearson" for 
#'Pearson R, "spearman" for Spearman rank correlation, "kendall"
#'for Kendall rank correlation and "sparcc" for SparCC.
#'@param pattern Character, input the specified pattern.
#'@param taxrank Character, input the taxonomic level to use.
#'@param taxa If the pattern is defined by using a specific taxon, input the name
#'of that taxa here.
#'@param variable Input the name of the experimental factor.
#'@param datatype If the data is marker gene, input "16S". If the
#'data is metagenomic, use "metageno".
#'@param shotfeat Only valid for SDP module, set to "null".
#'@param shotgunid Only valid for SDP module, set to "NA".
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export


Match.Pattern <- function(mbSetObj, dist.name="pearson", pattern=NULL, taxrank, variable, appendname){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(pattern %in% names(mbSetObj$dataSet$sample_data)){
    new.template = as.numeric(sample_data(mbSetObj$dataSet$norm.phyobj)[[pattern]])
  }else{
    if(!.on.public.web){
      clslbl <- sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]];
    }
    
    clslbl <- as.factor(clslbl);
    if(is.null(pattern)){
      pattern <- paste(1:length(levels(clslbl)), collapse="-");
    }
    
    templ <- as.numeric(ClearStrings(strsplit(pattern, "-")[[1]]));

    if (all(na.omit(templ) == na.omit(templ)[1])) {
      AddErrMsg("Cannot calculate correlation on constant values!")
      return(0)
    }

    
    new.template <- vector(mode="numeric", length=length(clslbl))
    # expand to match each levels in the analSet$cls
    all.lvls <- levels(clslbl);

    if(length(templ)!=length(all.lvls)){
      AddErrMsg("Wrong template - must the same length as the group number!");
      return(0);
    }
    
    for(i in 1:length(templ)){
      hit.inx <- clslbl == all.lvls[i]
      new.template[hit.inx] = templ[i];
    }
 
  }
 
  if(mbSetObj$module.type=="mdp"){
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    }
    
    if(taxrank=="OTU"){
      data <- phyloseq_objs$count_tables$OTU
    }else{
      
      data <- phyloseq_objs$merged_obj[[taxrank]]
      if(is.null(data)){
        AddErrMsg("Errors in projecting to the selected taxanomy level!");
        return(0);
      }
      
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      nm <- as.character(tax_table(data)[,taxrank]);
      nm.idx = match(rownames(phyloseq_objs$count_tables[[taxrank.inx]]),nm);
      y <- which(is.na(nm)==TRUE);
      
      #converting NA values to unassigned
      
      nm[y] <- "Not_Assigned";
      data1 <- as.matrix(otu_table(data));
      
      if(appendname=="T"){
        
        all_nm <- colnames(tax_table(data));
        hg_nmindx <- which(all_nm==taxrank)-1;
        
        if(hg_nmindx!=0){
          nma <- as.character(tax_table(data)[,hg_nmindx]);
          y1 <- which(is.na(nma)==TRUE);
          nma[y1] <- "Not_Assigned";
          nm <- paste0(nma,"_",nm);
          ind <- which(nm=="Not_Assigned_Not_Assigned");
          nm[ind] <- "Not_Assigned";
          nm <- gsub("_Not_Assigned", "",nm, perl = TRUE);
        }

        idx.table =  data.frame(original = rownames(phyloseq_objs$count_tables[[taxrank.inx]]),
                                taxPrepend = nm[nm.idx],stringsAsFactors = F)
        idx.table$taxPrepend[is.na(idx.table$taxPrepend)] <-  'Not_Assigned'
        data <- phyloseq_objs$count_tables[[taxrank.inx]]
        rownames(data) <- idx.table$taxPrepend
      }else{
        data <- phyloseq_objs$count_tables[[taxrank.inx]]
      }
      
    }
  }else{
    taxrank <- "OTU";
    data <- as.matrix(otu_table(mbSetObj$dataSet$norm.phyobj));
  }

  data <- t(data);
  qs::qsave(data, "match_data.qs")
 
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  boxdata <- as.data.frame(data,check.names=FALSE);
  boxdata$class <- clslbl;
  mbSetObj$analSet$boxdata <- boxdata;
  
  if(dist.name == "sparcc"){
    
    pattern.data <- cbind(data, new.template)
    qs::qsave(t(pattern.data), "pattern_data.qs")
    permNum <- 100
    pvalCutoff <- 1 
    corrCutoff <- 0
    
    if(.on.public.web){
      sparcc_results <- RunFastSpar_mem(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network", "pattern")
    }else{
      sparcc_results <- RunFastSpar(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network", "pattern")
    }
    
    cbtempl.results <- sparcc_results[which(sparcc_results$Taxon1=="new.template"), 2:4]
    fdr.col <- p.adjust(cbtempl.results[,3], "fdr");
    rownames(cbtempl.results) <- cbtempl.results[,1]
    cor.res <- signif(cbind(cbtempl.results[,-1], fdr.col),5)
    ord.inx <- order(cor.res[,2])
    sig.mat <- cor.res[ord.inx,]
    colnames(sig.mat) <- c("Correlation", "P-Value", "FDR");
    
  }else if(dist.name %in% c("spearman", "kendall", "pearson")){
    cbtempl.results <- apply(data, 2, template.match, new.template, dist.name);
    cor.res <- t(cbtempl.results);
    fdr.col <- p.adjust(cor.res[,3], "fdr");
    cor.res <- cbind(cor.res, fdr.col);
    colnames(cor.res) <- c("Correlation", "T-Stat", "P-Value", "FDR");
    ord.inx <- order(cor.res[,3]);
    sig.mat <- signif(cor.res[ord.inx,],5);
  }else{
    AddErrMsg("Distance measure invalid!")
  }

  fileName <- "correlation_feature.csv";
  fast.write(sig.mat, file=fileName);
  mbSetObj$analSet$resTable <- as.data.frame(sig.mat,check.names=FALSE);
  #removing Inf values from table
  is.na(mbSetObj$analSet$resTable) <- sapply(mbSetObj$analSet$resTable, is.infinite);
  mbSetObj$analSet$resTable[is.na(mbSetObj$analSet$resTable)]<-0;
  
  if(mbSetObj$module.type=="sdp"){
    mbSetObj$analSet$pattern<-pattern;
  }else{
    mbSetObj$analSet$pattern<-pattern;
    mbSetObj$analSet$taxrank<-taxrank;
  }
  
  sig.nm <<- fileName;
  mbSetObj$analSet$corph.meth <- dist.name;
  mbSetObj$analSet$feat.corr <- FALSE;
  mbSetObj$analSet$cor.mat<-sig.mat;
  mbSetObj$analSet$pattern.var <- variable
  
  return(.set.mbSetObj(mbSetObj));
}

## plots for multi-factor linear model
# IMPORTANT NOTE: data must be log transformed, and ideally TSS scaled as well (follow Maaslin defaults).
# There is no way of ensuring this from current tool design. User must choose either TSS or LOG, or they could choose none for both. 
# We could take original data, but then it is not filtered and Maaslin takes a long time to run.
# Best case: have filtered but raw data, and enforce both TSS and LOG in this function
# filt.data.orig <- qs::qread("filt.data.orig")@.Data %>% as.data.frame() <- this is what we want for OTU. Trying to find something similar for other taxa.
# See note on general_io.R -> "JE note" etc. Filtered, unnormalized tables do not exist other than OTU. Not sure if this is correct.
ProcessMaaslin <- function(
    mbSetObj,
    analysis.var,
    is.norm = "false",
    comp = NULL,
    ref = NULL,
    block = "NA",
    taxrank = "NA",
    model="LM",
    imgNm = "NA",
    thresh = 0.05){

  if(.on.public.web){
    # make this lazy load
    if(!exists(".perform.my.maaslin")){ # public web on same user dir
      .load.scripts.on.demand("utils_maaslin.Rc");    
    }
  }
  .perform.my.maaslin(mbSetObj,analysis.var,is.norm,comp, ref, block, taxrank,model,imgNm, thresh);
}
setMaaslinModel <- function(model){
print("model")
model<<-model
}
###################
## Transformation #
###################
transformFeatures = function(features, transformation) {
  if (transformation == 'LOG')     {
    features <- apply(features, 2, LOG)
  }
  return(features)
}


##################
## Normalization #
##################

normalizeFeatures = function(features, normalization) { 
  if (normalization == 'TSS')
  {
    features <- TSSnorm(features)
  }
  
  if (normalization == 'NONE')
  {
    features <- features
  }
 
  return(features)
}

######################
## TSS Normalization #
######################

# Apply TSS Normalization To A Dataset

TSSnorm = function(features) {
  # Convert to Matrix from Data Frame
  features_norm = as.matrix(features)
  dd <- colnames(features_norm)
  
  # TSS Normalizing the Data
  features_TSS <-
    vegan::decostand(
      features_norm,
      method = "total",
      MARGIN = 1,
      na.rm = TRUE)
  
  # Convert back to data frame
  features_TSS <- as.data.frame(features_TSS)
  
  # Rename the True Positive Features - Same Format as Before
  colnames(features_TSS) <- dd
  
  
  # Return
  return(features_TSS)
}


######################
# Log Transformation #
######################

# Log Transformation
LOG <- function(x) {
  y <- replace(x, x == 0, min(x[x>0]) / 2)
  return(log2(y))
}

###################################
########## Util Functions #########
###################################

# use memoise tech to avoid repeating
RunFastSpar_mem <- function(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output="network", opt="corr"){
  if(!exists("fastsparr_mem")){
    require("memoise");
    fastsparr_mem <<- memoise(RunFastSpar); 
  }
  
  res <- fastsparr_mem(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output=output, opt=opt);
  return(res);
}

#'Function to call FastSpar
#'@description This function runs the fastspar 
#'@param mbSetObj Input the name of the mbSetObj.
#'@import igraph 
#'@export
RunFastSpar <- function(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output="network", opt="corr"){
  
    if(!exists("my.fast.spar")){ # public web on same user dir
        .load.scripts.on.demand("utils_fastspar.Rc");    
    }
    return(my.fast.spar(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output, opt));
}

# Set of functions to perform cor.test
PearsonCorrFunc <- function(var1, var2, data){
  result <- cor.test(data[,var1], data[,var2])
  data.frame(var1, var2, result[c("estimate", "p.value", "statistic", "method")], stringsAsFactors = FALSE,check.names=FALSE)
}

SpearmanCorrFunc <- function(var1, var2, data){
  result <- cor.test(data[,var1], data[,var2], method = "spearman", exact = FALSE)
  data.frame(var1, var2, result[c("estimate", "p.value", "statistic", "method")], stringsAsFactors = FALSE,check.names=FALSE)
}

KendallCorrFunc <- function(var1, var2, data){
  result <- cor.test(data[,var1], data[,var2], method = "kendall", exact = FALSE)
  data.frame(var1, var2, result[c("estimate", "p.value", "statistic", "method")], stringsAsFactors = FALSE,check.names=FALSE)
}

#'Function to generate templates
#'@description This functions generate templates
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

GenerateTemplates <- function(mbSetObj, variable){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  print(variable)
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  level.len <- length(levels(clslbl));
  # only specify 4: increasing, decreasing, mid high, mid low, constant
  incs <- 1:level.len;
  desc <- level.len:1;
  
  if(level.len > 2){
    # use ceiling, so that the peak will be right for even length
    mid.pos <- ceiling((level.len+1)/2);
    mid.high <- c(1:mid.pos, seq(mid.pos-1,by=-1,length.out=level.len-mid.pos));
    mid.low <- c(mid.pos:1, seq(2, length.out=level.len-mid.pos));
    res <- rbind(incs, desc, mid.high, mid.low); # add the constant one
  }else{
    res <- rbind(incs, desc);
  }
  
  # turn into string
  res <- apply(res, 1, paste, collapse="-");
  # add the legends
  res <- c(paste(levels(clslbl), collapse="-"), res);
  clslbl<<-clslbl;
  return(res);
}

# helper function
# calculate geometric means prior to estimate size factors
gm_mean = function(x, na.rm=TRUE){
  exp(sum(log(x[x > 0]), na.rm=na.rm) / length(x))
}

# helper function
phyloseq_to_edgeR = function(physeq, group, method="RLE", ...){
  
  load_edgeR();
  load_phyloseq();
  
  # Enforce orientation.
  if( !taxa_are_rows(physeq) ){ physeq <- t(physeq) }
  x = as(otu_table(physeq), "matrix")
  # Add one to protect against overflow, log(0) issues.
  x = x + 1
  # Check `group` argument
  if( identical(all.equal(length(group), 1), TRUE) & nsamples(physeq) > 1 ){
    # Assume that group was a sample variable name (must be categorical)
    group = phyloseq::get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"),check.names=FALSE);
  }
  
  # Now turn into a DGEList
  y = edgeR::DGEList(counts=x, group=group, genes=taxonomy, remove.zeros = TRUE, ...)
  z = edgeR::calcNormFactors(y, method="RLE");
  # Estimate dispersions
  return(edgeR::estimateTagwiseDisp(edgeR::estimateCommonDisp(z)))
}

phyloseq_to_metagenomeSeq = function (physeq, ...) 
{
  if (!taxa_are_rows(physeq)) {
    physeq <- t(physeq)
  }
  countData = round(as(otu_table(physeq), "matrix"), digits = 0)
  if (!is.null(sample_data(physeq, FALSE))) {
    ADF = Biobase::AnnotatedDataFrame(data.frame(sample_data(physeq),check.names=FALSE));
  }else {
    ADF = NULL;
  }
  taxonomy = tax_table(physeq, errorIfNULL=FALSE);
  if (!is.null(taxonomy)) {
    TDF = Biobase::AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), data.frame(as(taxonomy, "matrix")), row.names = taxa_names(physeq),check.names=FALSE))
    #  data.frame(tax_table(physeq)@.Data), row.names = taxa_names(physeq)))
    #  data.frame(tax_table(physeq)), row.names = taxa_names(physeq))) ## bugs with tax_table
  }else {
    TDF = Biobase::AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), row.names = taxa_names(physeq),check.names=FALSE));
  }
  if (requireNamespace("metagenomeSeq")) {
    mrobj = metagenomeSeq::newMRexperiment(counts = countData, 
                                           phenoData = ADF, featureData = TDF, ...)
    if (sum(colSums(countData > 0) > 1) < ncol(countData)) {
      p = suppressMessages(metagenomeSeq::cumNormStat(mrobj))
    }
    else {
      p = suppressMessages(metagenomeSeq::cumNormStatFast(mrobj))
    }
    mrobj = metagenomeSeq::cumNorm(mrobj, p = p)
    return(mrobj)
  }
}

# Helper function
# Run template on all the high region effect genes
template.match <- function(x, template, dist.name) {
  k<-cor.test(x,template, method=dist.name);
  c(k$estimate, k$stat, k$p.value)
}

# Helper function for volcano plot
GetTopInx <- function (vec, n, dec = T) {
  inx <- order(vec, decreasing = dec)[1:n];
  vec <- rep(F, length = length(vec));
  vec[inx] <- T;
  return(vec);
}

####################################
############ Get Funs ##############
####################################

GetMDI <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$mdi)
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetSigTable.UNIVAR<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  GetSigTable(mbSetObj$analSet$univar$resTable, "Univariate analysis");
}


GetSigTable.Corr<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  GetSigTable(mbSetObj$analSet$cor.mat, "Pattern search using correlation analysis");
}


GetSigTable.METAGENOSEQ<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  GetSigTable(mbSetObj$analSet$metagenoseq$resTable, "metagenomeSeq");
}


GetSigTable.RNASeq<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  GetSigTable(mbSetObj$analSet$rnaseq$resTable, "RNAseq method");
}


GetSigTable.LEFSE<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  GetSigTable(mbSetObj$analSet$lefse$resTable, "LEfSe");
}


GetRFConf.Table<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  load_xtable();
  print(xtable::xtable(mbSetObj$analSet$rf$confusion, caption="Random Forest Classification Performance"), size="\\scriptsize");
}


GetMapTable<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  load_xtable();
  print(xtable::xtable(mbSetObj$analSet$mapTable, caption="Result from Taxa Name Mapping"),
        tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}

####
#### WEB UTILS
####

# get the OOB error for the last signif
GetRFOOB<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  errors = mbSetObj$analSet$rf$err.rate;
  nrow = dim(errors)[1];
  signif(errors[nrow, 1],3);
}

# significance measure, double[][]
GetRFSigMat<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(CleanNumber(mbSetObj$analSet$rf.sigmat))
}

GetRFSigRowNames<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rownames(mbSetObj$analSet$rf.sigmat));
}

GetRFSigColNames<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(colnames(mbSetObj$analSet$rf.sigmat));
}

# return double[][] confusion matrix
GetRFConfMat<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(signif(mbSetObj$analSet$rf$confusion,3));
}

GetRFConfRowNames<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rownames(mbSetObj$analSet$rf$confusion));
}

GetRFConfColNames<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(colnames(mbSetObj$analSet$rf$confusion));
}

GetSigTable.MMPMic<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  GetSigTable(mbSetObj$analSet$resTable, "MaasLin");
}

GetSigTable.MMPMet<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  GetSigTable(mbSetObj$dataSet$metabolomics$resTable, "Limma");
}

GetMMPMicTable<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  load_xtable();
  print(xtable::xtable(mbSetObj$analSet$mic.map, caption="Result from Taxa Name Mapping"),
        tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}


GetMMPMetTable<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  load_xtable();
  print(xtable::xtable(mbSetObj$analSet$met.map, caption="Result from Metabolite Name Mapping"),
        tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}

#Generate json file for plotly for comparison tests
GenerateCompJson <- function(mbSetObj = NA, fileName, type, mode = 1, taxlvl, parent = "Phylum", sigLevel = 0.05, fcLevel = 0) {
  library(RColorBrewer)
  library("dplyr")
  print(c(taxlvl, parent, type,sigLevel,fcLevel))
  sigLevel <- as.numeric(sigLevel)
    fcLevel <- as.numeric(fcLevel)
  mbSetObj <- .get.mbSetObj(mbSetObj)
  resList <- ""
  if (type %in% c("tt", "nonpar")) {
    resTable <- mbSetObj$analSet$univar$resTable
    resTable$id <- rownames(resTable)
    resList <- list(data = resTable, param = mbSetObj$paramSet$univar)
  } else if (type %in% c("zigfit", "ffm")) {
    resTable <- mbSetObj$analSet$metagenoseq$resTable
    resTable$id <- rownames(resTable)
    resList <- list(data = resTable, param = mbSetObj$paramSet$metagenoseq)
  } else if (type %in% c("EdgeR", "DESeq2")) {
    resTable <- mbSetObj$analSet$rnaseq$resTable
    resTable$id <- rownames(resTable)
    resList <- list(data = resTable, param = mbSetObj$paramSet$rnaseq)
  }

 if(mbSetObj[["module.type"]]=="mdp"){
tax_table <- data.frame(mbSetObj[["dataSet"]][["proc.phyobj"]]@tax_table)
  
  if (taxlvl == "OTU") {
    resTable[["parent"]] <- tax_table[[parent]][match(resTable$id, rownames(tax_table))]
  } else {
    resTable[["parent"]] <- tax_table[[parent]][match(resTable$id, tax_table[[taxlvl]])]
  }
  resTable[["parent"]][is.na(resTable[["parent"]])] <- "Not_Assigned"


  resTable <- resTable %>%
    group_by(parent) %>%
    group_modify(~ .x %>%
      mutate(len = if (n() > 3) {
        sample(1:(n() %/% 3), size = n(), replace = TRUE)
      } else {
        seq_len(n())
      }))


  don <- data.frame(resTable %>%
    group_by(parent) %>%
    summarise(chr_len = max(len)) %>%
    mutate(tot = cumsum(chr_len) - chr_len) %>%
    select(-chr_len) %>%
    left_join(resTable, ., by = c("parent" = "parent")) %>%
    arrange(parent, len) %>%
    mutate(BPcum = len + tot))

  if (type %in% c("EdgeR", "DESeq2")) {
    don$shape <- ifelse(don$log2FC > 0, "triangle-up", "triangle-down")
    don$shape[don$FDR > sigLevel | abs(don$log2FC) < fcLevel] <- "circle"
    resList$param$multigroup <- FALSE
  } else {
    resList$param$multigroup <- TRUE
    don$shape <- "circle"
  }

  colors <- setNames(colorRampPalette(brewer.pal(8, "Set1"))(length(unique(resTable$parent))), unique(resTable$parent))
  don$color <- unname(colors[match(don$parent, names(colors))])
  don$color[don$FDR > sigLevel | abs(don$log2FC ) <fcLevel] <- "#808080"
  don$size <- 7 + (don$logCPM - min(don$logCPM)) * 7 / (max(don$logCPM) - min(don$logCPM))

  axisdf <- don %>%
    group_by(parent) %>%
    summarize(
      center = (max(BPcum) + min(BPcum)) / 2,
      start = min(BPcum) - 0.5, end = max(BPcum) + 0.5
    )
  resList$data <- don
  resList$axisdf <- axisdf
  resList$param$parent <- parent

  if (mode == 2) {
    axisdf$color <- scales::alpha(unname(colors[match(axisdf$parent, names(colors))]), 0.1)

    areadf <- lapply(1:nrow(axisdf), function(x) {
      info <- list(
        x = c(axisdf$start[x], axisdf$start[x], axisdf$end[x], axisdf$end[x], axisdf$start[x]),
        y = c(0, ceiling(max(-log10(resTable$Pvalues))), ceiling(max(-log10(resTable$Pvalues))), 0, 0),
        mode = "line",
        fill = "toself",
        fillcolor = axisdf$color[x],
        line = list(
          color = "rgba(255,255,255,0.1)",
          width = 0,
          dash = "solid"
        ),
        showlegend = F,
        name = axisdf$parent[x],
        legendgroup = axisdf$parent[x],
        hoverinfo = ""
      )
      return(info)
    })

    resList$areadf <- areadf
  }

 }

  json.obj <- rjson::toJSON(resList)
  sink(fileName)
  cat(json.obj)
  sink()
  return(1)
}


PlotlyCompRes <- function(mbSetObj = NA, type="", fileName="") {
  #save.image("comp.res");
  library(htmlwidgets)
  library(plotly)

  mbSetObj <- .get.mbSetObj(mbSetObj)
  if (type %in% c("univ")) {
    resTable <- mbSetObj$analSet$univar$resTable
    params <- mbSetObj$paramSet$univar
  } else if (type %in% c("metagenome")) {
    resTable <- mbSetObj$analSet$metagenoseq$resTable
    params <- mbSetObj$paramSet$metagenoseq
  } else { 
    resTable <- mbSetObj$analSet$rnaseq$resTable
    params <- mbSetObj$paramSet$rnaseq
  }
  p.lvl <- params$p.lvl
  fc.thresh <- params$fc.thresh

  resTable$id <- rownames(resTable)
  raw_data <- resTable
  
  if ("log2FC" %in% names(raw_data)) {

    p <- plot_ly(
      data = raw_data, 
      x = ~log2FC, 
      y = -log10(raw_data$Pvalues), 
      type = 'scatter',
      mode = 'markers',
      marker = list(
        color = mapply(getColor, raw_data$FDR, raw_data$log2FC, p.lvl, fc.thresh), # getColor function should be defined
        size=mapply(getSizeForPValue, raw_data$FDR, p.lvl),
        line = list(color = 'white', width = 0.8)
      ),
      text = ~paste("Feature ID: ", id, 
                    "<br>P-value: ", format(Pvalues, scientific = TRUE),
                    "<br>FDR: ", format(FDR, scientific = TRUE)),
      hoverinfo = 'text',
      label= ~id
    )
  layout <- list(
    xaxis = list(title = "log2FC"),
    yaxis = list(title = '-log10(FDR)')
  )

  } else {
    p <- plot_ly(
      data = raw_data, 
      x = seq_along(raw_data$id),
      y = ~-log10(FDR),
      type = 'scatter',
      mode = 'markers',
      marker = list(
        color = mapply(getColorForPValue, raw_data$FDR, p.lvl), # getColorForPValue function should be defined
        size = mapply(getSizeForPValue, raw_data$FDR, p.lvl),
        line = list(color = 'white', width = 0.8)
      ),
      text = ~paste("Feature ID: ", id, 
                    "<br>P-value: ", format(Pvalues, scientific = TRUE),
                    "<br>FDR: ", format(FDR, scientific = TRUE)),
      hoverinfo = 'text',
      label= ~id
    )
  layout <- list(
    xaxis = list(title = "Features"),
    yaxis = list(title = '-log10(FDR)')
  )
  
  }
  
  p <- p %>% layout(layout);
  p <- p %>% onRender("
  function(el, x) {
    el.on('plotly_click', function(data) {
        var pointIndex = data.points[0].pointIndex;
        console.log(data.points[0].data.label)
        parent.window.document.getElementById('form1:selectedVarInput').value = data.points[0].data.label[pointIndex];
        parent.window.plotFeature();
    });
  }
");
  return(p)
}

getColor <- function(pValue, log2fc, pval.thresh, log2fcThreshold) {
  pValueThreshold <- pval.thresh# Replace with the actual way to access this value in R
  
  if (pValue < pValueThreshold && abs(log2fc) > log2fcThreshold) {
    if (log2fc > 0) {
      return('red') # Interpolate red for positive log2fc
    } else {
      return('blue') # Interpolate blue for negative log2fc
    }
  } else {
    return('grey') # Non-significant
  }
}

getColorForPValue <- function(val, pval.thresh) {
  if (abs(val) >= pval.thresh) {
    return('grey') # Non-significant
  } else {
    return('blue') # Significant (simplified example)
  }
}

getSizeForPValue <- function(val, pval.thresh) {
  if (abs(val) >= pval.thresh) {
    return(5) # Non-significant
  } else {
    return(10) # Significant (simplified example)
  }
}

PlotCompRes <- function(mbSetObj = NA, type = "", imgName = "") {
  library(ggplot2)
  library(Cairo)  # For high-quality image output
  
  # Process the data similar to the original function
  mbSetObj <- .get.mbSetObj(mbSetObj)  # Assuming this is a function defined elsewhere
  fileName <- paste0(imgName, ".png");
  if (type %in% c("univ")) {
    resTable <- mbSetObj$analSet$univar$resTable
    params <- mbSetObj$paramSet$univar
    mbSetObj$imgSet$univar.plot <- fileName;
  } else if (type %in% c("metagenome")) {
    resTable <- mbSetObj$analSet$metagenoseq$resTable
    params <- mbSetObj$paramSet$metagenoseq
    mbSetObj$imgSet$metagenoseq.plot <- fileName;
  } else { 
    resTable <- mbSetObj$analSet$rnaseq$resTable
    params <- mbSetObj$paramSet$rnaseq
    mbSetObj$imgSet$rnaseq.plot <- fileName;
  }
  
  p.lvl <- params$p.lvl
  fc.thresh <- params$fc.thresh
  
  resTable$id <- rownames(resTable)
  raw_data <- resTable
  
  # Define color and size mapping based on your getColor and getSizeForPValue functions
  raw_data$color <- mapply(getColor, raw_data$FDR, raw_data$log2FC, MoreArgs = list(p.lvl, fc.thresh))
  raw_data$size <- mapply(getSizeForPValue, raw_data$FDR, MoreArgs = list(p.lvl))
  
  # Create the ggplot object based on the type of data
  if ("log2FC" %in% names(raw_data)) {
    p <- ggplot(raw_data, aes(x = log2FC, y = -log10(Pvalues), color = color, size = size)) +
      geom_point(alpha = 0.6) +
      scale_size_continuous(range = c(1, 10)) +
      scale_color_identity() +
      labs(x = "log2FC", y = "-log10(P-value)", title = "Volcano Plot") +
      theme_minimal()
  } else {
    p <- ggplot(raw_data, aes(x = seq_along(id), y = -log10(FDR), color = color, size = size)) +
      geom_point(alpha = 0.6) +
      scale_size_continuous(range = c(1, 10)) +
      scale_color_identity() +
      labs(x = "Feature Index", y = "-log10(FDR)", title = "Feature Plot") +
      theme_minimal()
  }
  
  # Save the plot using Cairo for high quality output
  if (fileName != "") {
    Cairo::Cairo(file = fileName, width = 1000, height = 800, dpi=72)
    print(p)
    dev.off()
  } else {
    print(p)
  }
  .set.mbSetObj(mbSetObj)
  return(mbSetObj)
}


#univ/metagenome/rnaseq/lefse/maaslin/list/global
SetKeggProjectionType <- function(mbSetObj=NA, type){
    mbSetObj <- .get.mbSetObj(mbSetObj);
    mbSetObj$paramSet$koProj.type <- type;
    return(.set.mbSetObj(mbSetObj));
}

SetKEGGNetVisOpt <- function(mbSetObj, nm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$analSet$keggnet$background <- nm;
  return(.set.mbSetObj(mbSetObj));
}