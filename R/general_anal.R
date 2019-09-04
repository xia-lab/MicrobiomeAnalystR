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
    data1 <- as.data.frame(t(otu_table(data)));
  }else{
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- readRDS("phyloseq_objs.RDS")
    }
    
    if(taxrank=="OTU"){
      data1 <- t(phyloseq_objs$count_tables$OTU)
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- t(phyloseq_objs$count_tables[[taxrank.inx]])
    } 
  }
  
  data.impfeat <<- data1;
  cls <- sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]];
  variable <<- variable;
  rf_out <- randomForest(data1,cls, ntree = treeNum, mtry = tryNum, importance = TRUE, proximity = TRUE);
  # set up named sig table for display
  impmat <- rf_out$importance;
  impmat <- impmat[rev(order(impmat[,"MeanDecreaseAccuracy"])),]
  sigmat <- impmat[,"MeanDecreaseAccuracy", drop=F];
  sigmat <- signif(sigmat, 5);
  write.csv(sigmat,file="randomforests_sigfeatures.csv");
  
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
  
  if(is.na(width)){
    if(feature < 5 ){
      h <- feature*1.2;
      w <- 9;
    } else if (feature < 10){
      h <- feature*1.4;
      w <- 9;
    } else if (feature < 15){
      h <- feature/1.6;
      w <- 9;
    } else if (feature < 20){
      h <- feature/1.8;
      w <- 9;
    } else if (feature < 25){
      h <- feature/2;
      w <- 9;
    } else if (feature < 30){
      h <- feature/2.2;
      w <- 9;
    } else if (feature < 40){
      h <- feature/2.5;
      w <- 9;
    } else {
      h <- feature/10;
      w <- 9;
    }
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
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  imgName = paste(imgName, ".", format, sep="");
  mbSetObj$imgSet$rf.imp <- imgName;
  vip.score <- rev(sort(mbSetObj$analSet$rf$importance[,"MeanDecreaseAccuracy"]));
  cls <- sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]];
  cls.length <- length(levels(cls))
  cls.width <- cls.length/2
  
  if(is.na(width)){
    if(feature < 5 ){
      h <- feature*1.2;
    } else if (feature < 10){
      h <- feature*1.4;
    } else if (feature < 15){
      h <- feature/1.6;
    } else if (feature < 20){
      h <- feature/1.8;
    } else if (feature < 25){
      h <- feature/2;
    } else if (feature < 30){
      h <- feature/2.2;
    } else if (feature < 40){
      h <- feature/2.5;
    } else {
      h <- feature/10;
    }
    
    if(cls.length < 5){
      w <- 9.25;
    } else if(cls.length < 10){
      w <- 11.5;
    } else if(cls.length < 15){
      w <- 12.5;
    }else{
      w <- 15;
    }
  }else if(width == 0){
    w <- 8;
  }else{
    w <- width;
  }
  
  Cairo::Cairo(file = imgName,  unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  PlotImpVar(mbSetObj, vip.score,"MeanDecreaseAccuracy", feature);
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
  
  cls.len <- length(levels(mbSetObj$analSet$cls));
  
  if(cls.len == 2){
    rt.mrg <- 5;
  }else if(cls.len == 3){
    rt.mrg <- 6;
  }else if(cls.len == 4){
    rt.mrg <- 7;
  }else if(cls.len == 5){
    rt.mrg <- 8;
  }else if(cls.len == 6){
    rt.mrg <- 9;
  }else{
    rt.mrg <- 11;
  }
  op <- par(mar=c(5,9,2,rt.mrg)); # set right side margin with the number of class
  
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
  mns <- by(data1[, names(imp.vec)], mbSetObj$analSet$cls,
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
  
  cls.lbl <- levels(mbSetObj$analSet$cls);
  
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
  par(op);
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

PerformUnivarTest <- function(mbSetObj, variable, p.lvl, shotgunid, taxrank, statOpt){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  #rather than whole name from taxonomy just last name.
  cls <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  lvl <- length(levels(cls));
  
  if(mbSetObj$module.type=="sdp"){
    taxrank<-"OTU";
    data <- mbSetObj$dataSet$norm.phyobj;
    data1 <- as.data.frame(t(otu_table(data)));
  }else{
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- readRDS("phyloseq_objs.RDS")
    }
    
    if(taxrank=="OTU"){
      data1 <- t(phyloseq_objs$count_tables$OTU)
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- t(phyloseq_objs$count_tables[[taxrank.inx]])
    } 
  }
  
  isNonPar <- statOpt=="nonpar"
  
  if(length(levels(cls)) > 2){
    resTable <- data.frame(GetAovRes(cls, data1, nonpar=isNonPar));
  }else{
    resTable <- data.frame(GetTtestRes(cls, data1, nonpar=isNonPar));
  }
  
  rownames(resTable) <- colnames(data1);
  colnames(resTable) <- c("Statistics","Pvalues");
  resTable$FDR <- p.adjust(resTable$Pvalues,method ="fdr");
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  resTable <- signif(resTable, digits = 5);
  resTable <- resTable[complete.cases(resTable), ];
  resTable <- resTable[,c(2,3,1)];
  write.csv(resTable, file="univar_test_output.csv");
  
  sigHits <- resTable$FDR < p.lvl;
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
  dat3t_boxplot <- as.data.frame(t(otu_table(data_boxplot)));
  colnames(dat3t_boxplot) <- nm_boxplot; 
  
  #individual boxplot for features
  box_data <- data.frame(dat3t_boxplot[,sigfeat %in% nm_boxplot]);
  box_data$class <- claslbl_boxplot;
  mbSetObj$analSet$boxdata <- box_data;
  write.csv(t(box_data), "uni_abund_data.csv")
  
  mbSetObj$analSet$anal.type <- "tt";
  mbSetObj$analSet$var.type <- variable;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$Univar$resTable <- mbSetObj$analSet$resTable <- resTable;
  mbSetObj$analSet$univar.taxalvl <- taxrank;
  
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

PerformMetagenomeSeqAnal<-function(mbSetObj, variable, p.lvl, shotgunid, taxrank, model){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_metagenomeseq();
  
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
    fit <- fitZig(data, mod);  
  }else{
    if(length(levels(cls)) > 2){
      current.msg <<- paste( "More than two groups present in your experimental factor. This model can only be used with two groups.");
      return(0);
    }else{
      fit <-fitFeatureModel(data, mod);
    }
  }
  
  x <- MRfulltable(fit, number = nrow(assayData(data)$counts));
  x <- x[!is.na(rownames(x)), ];
  
  rownames(x) <- gsub(":1", "", x = rownames(x), fixed = TRUE);
  x$OTUnames <- as.character(rownames(x))
  
  if (!is.null(tax_table(data, errorIfNULL = FALSE))) {
    #Attach the bacterial taxonomy to the table, if available
    TAX = data.frame(tax_table(data));
    TAX$OTUnames <- as.character(rownames(TAX));
    res = merge(x, TAX, by = "OTUnames")
  } else {
    res = x;
  }
  
  # Sort and return #sighits return TRUE or FALSE
  sigHits <-res$adjPvalues<=p.lvl;
  de.Num <- length(which(sigHits));
  
  if(de.Num == 0){
    current.msg <<- paste( "No significant features were identified using the given p value cutoff. Please change the cutoff limit.");
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features were identified!")
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
  
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  write.csv(resTable, file="metageno_de_output.csv");
  
  if(nrow(resTable) > 500){
    resTable <- resTable[1:500, ];
  }
  
  mbSetObj$analSet$metagenoseq$resTable <- mbSetObj$analSet$resTable <- data.frame(resTable);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  sigfeat <- rownames(resTable);
  
  #prepare individual boxplot
  box_data <- MRcounts(data);
  #subset only diff. Abundant features
  box_data <- box_data[sigfeat, ];
  
  #samples in rows
  box_data <- t(box_data);
  box_data <- data.frame(box_data);
  colnames(box_data) <- sigfeat;
  
  claslbl <- pData(data)[ ,variable];
  box_data$class <- unlist(claslbl);
  
  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$anal.type <- "metagseq";
  mbSetObj$analSet$var.type <- variable;
  mbSetObj$analSet$metageno.taxalvl <- taxrank;
  mbSetObj$analSet$id.type <- shotgunid;
  
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
    dat3t <- as.data.frame(t(otu_table(data)));
  }else{
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- readRDS("phyloseq_objs.RDS")
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
  
  wil_datadf <- as.data.frame(dat3t);
  
  #if no subclass within classes then no wilcoxon rank sum test  
  #Linear Discriminant analysis (LDA)
  ldares <- lda(claslbl~ .,data = wil_datadf);
  ldamean <- as.data.frame(t(ldares$means));
  class_no <<- length(unique(claslbl));
  ldamean$max <- apply(ldamean[,1:class_no],1,max);
  ldamean$min <- apply(ldamean[,1:class_no],1,min);
  ldamean$LDAscore <- signif(log10(1+abs(ldamean$max-ldamean$min)/2),digits=3);
  ldamean$Pvalues <- signif(rawpvalues,digits=5);
  ldamean$FDR <- signif(clapvalues,digits=5);
  resTable <- ldamean;
  
  # it seems lda add ` around names containing dash "-", need to strip this off
  rawNms <- rownames(resTable);
  rownames(resTable) <- gsub("`", '', rawNms);
  
  if(pvalOpt == "raw"){
    de.Num <- sum(rawpvalues<=p.lvl & ldamean$LDAscore>=lda.lvl)
  }else{
    de.Num <- sum(clapvalues<=p.lvl & ldamean$LDAscore>=lda.lvl)
  }
  
  if(de.Num == 0){
    current.msg <<- "No significant features were identified with given criteria.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features with given criteria.")
  }
  
  # sort by p value
  ord.inx <- order(resTable$Pvalues, resTable$LDAscore);
  resTable <- resTable[ord.inx, ,drop=FALSE];
  #p-values column to appear first; then FDR and then others
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];
  resTable <- resTable[,c(ncol(resTable),1:(ncol(resTable)-1))];
  
  # need to control digits
  resTable <- signif(resTable, 5);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  resTable$max <- resTable$min <- NULL;
  #if only two groups are present in sample variable
  cls.lbls <- levels(claslbl);
  grp.dat <- resTable[, cls.lbls];
  res.cls <- as.character(apply(grp.dat, 1, function(x){cls.lbls[which.max(x)]}));
  
  if(class_no==2){
    indx <- which(res.cls==unique(claslbl)[1]);
    resTable$LDAscore[indx]<--abs(resTable$LDAscore[indx]);
  }
  
  write.csv(resTable, file="lefse_de_output.csv");
  mbSetObj$analSet$lefse$resTable <- mbSetObj$analSet$resTable <- resTable;
  
  #subset dataset for bar plot visualization (LDA Score)
  ldabar <- as.data.frame(rownames(resTable));
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
  box_data <- as.data.frame(wil_datadf[, sigfeat]);
  colnames(box_data) <- sigfeat;
  box_data$class <- claslbl;
  
  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$anal.type <- "lefse";
  mbSetObj$analSet$lefse.taxalvl <- taxrank;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$meta <- variable;
  
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
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  set.seed(280561493);
  mbSetObj$analSet$lefse_plot <- imgName;
  imgName = paste(imgName, ".", format, sep="");
  ldabar <- ldabar;
  
  ldabar <-  ldabar[order(-abs(ldabar[[2]])), ];
  
  if(ldaFeature < nrow(ldabar)) {
    ldabar <- ldabar[1:ldaFeature,];
  };
  
  vip.score <- ldabar[[2]];
  names(vip.score) <- ldabar[[1]];
  
  if(is.na(width)){
    if(length(vip.score) < 5 ){
      h <- 6;
      w <- 9.25;
    } else if (length(vip.score) < 10){
      h <- length(vip.score)/1.3;
      w <- 9.25;
    } else if (length(vip.score) < 15){
      h <- length(vip.score)/1.6;
      w <- 9.25;
    } else if (length(vip.score) < 20){
      h <- length(vip.score)/1.8;
      w <- 9.25;
    } else if (length(vip.score) < 25){
      h <- length(vip.score)/2;
      w <- 9.25;
    } else if (length(vip.score) < 30){
      h <- length(vip.score)/2.2;
      w <- 9.25;
    } else if (length(vip.score) < 40){
      h <- length(vip.score)/2.5;
      w <- 9.25;
    } else {
      h <- length(vip.score)/4.5;
      w <- 9.25;
    };
  }else if(width == 0){
    w <- 8;
  }else{
    w <- width;
  }
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
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
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(layoutOptlf == "dot") {
    
    sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
    cls.len <- length(levels(sample_table[[meta]]));
    
    if(cls.len == 2){
      rt.mrg <- 5;
    }else if(cls.len == 3){
      rt.mrg <- 6;
    }else if(cls.len == 4){
      rt.mrg <- 7;
    }else if(cls.len == 5){
      rt.mrg <- 8;
    }else if(cls.len == 6){
      rt.mrg <- 9;
    }else{
      rt.mrg <- 11;
    }
    op <- par(mar=c(5,9,2,rt.mrg)); # set right side margin with the number of class
    
    feat.num <- length(imp.vec);
    
    # first get the top subset
    imp.vec <- rev(sort(imp.vec))[1:feat.num];
    
    # reverser the order for display
    imp.vec <- sort(imp.vec);
    
    # as data should already be normalized, use mean/median should be the same
    # mns is a list contains means of all vars at each level
    # conver the list into a matrix with each row contains var averages across different lvls
    
    #data.impfeat <- as.matrix(t(otu_table(mbSetObj$dataSet$filt.data)))
    data1 <- data.impfeat_lefse;
    mns <- by(data1[, names(imp.vec)], 
              sample_table[[meta]],
              # sample_table[["SampleType"]],
              function(x){ # inner function note, by send a subset of dataframe
                apply(x, 2, mean, trim=0.1)
              });
    mns <- t(matrix(unlist(mns), ncol=feat.num, byrow=TRUE));
    
    vip.nms <- names(imp.vec);
    names(imp.vec) <- NULL;
    
    # modified for B/W color
    dotcolor <- ifelse(color.BW, "darkgrey", "#585855");
    dotchart(imp.vec, bg=dotcolor, xlab= "LDA score", cex=1.35);
    
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
    
    cls.lbl <- levels(sample_table[[meta]]);
    
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
    
    par(op);
    
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

PerformRNAseqDE<-function(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  
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
  
  dat3t <- as.data.frame(t(otu_table(data)));
  colnames(dat3t) <- nm;
  
  mbSetObj$analSet$rnaseq$data.rnaseq <- dat3t
  
  if(opts=="DESeq2"){
    # only for small data set (< 100)
    if(length(claslbl) > 100){
      current.msg <<- "Only EdgeR is supported for sample size over 100.";
      return(0);
    }else{
      
      load_deseq();
      
      # create formula based on user selection
      my.formula <- as.formula(paste("~", variable));
      
      #converting from phyloslim object to deseq
      diagdds = phyloseq_to_deseq2(data, my.formula);
      geoMeans = apply(counts(diagdds), 1, gm_mean);
      diagdds = DESeq2::estimateSizeFactors(diagdds, geoMeans = geoMeans);
      diagdds = DESeq2::DESeq(diagdds, test="Wald", fitType="parametric");
      res = DESeq2::results(diagdds, independentFiltering = FALSE, cooksCutoff =  Inf);
      sigHits <- which(res$padj < p.lvl);
      de.Num <- length(sigHits);
      
      if(de.Num == 0){
        current.msg <<- "No significant features were identified using the given p value cutoff.";
      }else{
        current.msg <<- paste("A total of", de.Num, "significant features were identified!");
      }
      
      resTable <- res[,c("log2FoldChange" ,"lfcSE","pvalue","padj")];
      resTable <- signif(data.matrix(resTable), digits=5);
      colnames(resTable) <- c("log2FC","lfcSE","Pvalues","FDR");
      mbSetObj$analSet$anal.type <- "deseq";
    }
  }else{
    #using by filtered data ,RLE normalization within it.
    dge <- phyloseq_to_edgeR(data, variable);
    et = edgeR::exactTest(dge);
    tt = edgeR::topTags(et, n=nrow(dge$table), adjust.method="BH", sort.by="PValue");
    res = tt@.Data[[1]];
    de.Num <- sum(res$FDR < p.lvl);
    
    if(de.Num == 0){
      current.msg <<- "No significant features were identified using the given p value cutoff.";
    }else{
      current.msg <<- paste("A total of", de.Num, "significant features were identified!");
    }
    
    resTable <- res[,c("logFC","logCPM","PValue","FDR")];
    resTable <- signif(data.matrix(resTable), digits = 5);
    colnames(resTable) <- c("log2FC","logCPM","Pvalues","FDR");
    mbSetObj$analSet$anal.type <- "edgr";
  }
  
  resTable <- as.data.frame(resTable);
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  write.csv(resTable, file="rnaseq_de.csv");
  
  if(nrow(resTable) > 500){
    resTable<-resTable[1:500, ];
  }
  
  mbSetObj$analSet$rnaseq$resTable <- mbSetObj$analSet$resTable <- as.data.frame(resTable);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  taxrank <<- taxrank;
  
  #individual boxplot for features
  sigfeat <- rownames(resTable);
  box_data <- as.data.frame(dat3t[ ,sigfeat]);
  colnames(box_data) <- sigfeat;
  box_data$class <- claslbl;
  
  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$var.type <- variable;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$rnaseq.taxalvl <- taxrank;
  mbSetObj$analSet$rnaseq.meth <- opts;
  
  tree_data <<- data;
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to plot volacano analysis results from RNAseq analysis
#'@description This functions plots the RNAseq analysis 
#'results from the microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param fc.cutoff Numeric, input the fold-change cutoff.
#'@param p.cutoff Numeric, input the p-value cutoff.
#'@param colpal Character, "default" for a blue and red color
#'palette, "pog" for the purple and orange color palette, "bog"
#'for the blue and orange color palette, and "rgg" for the red and 
#'green color palette.
#'@param imgName Character, input the name of the plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import ggplot2
#'@import ggrepel
#'@import viridis
PlotRNAseqVolcano <- function(mbSetObj, fc.cutoff = 2, p.cutoff = 0.05, colpal = "default", imgName="volcano",
                              format = "png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_ggplot();
  load_ggrepel();
  
  rna_results <- mbSetObj$analSet$rnaseq$resTable
  
  pvals <- rna_results$Pvalues
  inx.p <- pvals <= p.cutoff
  pvals[pvals == 0] <- .Machine$double.xmin
  logpvals <- -log10(pvals)
  
  df <- data.frame(id = rownames(rna_results), logfc = rna_results$log2FC, logp = logpvals, inx.p = inx.p)
  df$sig <- ifelse(df$logfc >= fc.cutoff,"A", ifelse(df$logfc<=-fc.cutoff , "B", "C"))
  
  # label significant features
  imp.inx <- (df$sig == "A" | df$sig == "B") & df$inx.p;
  
  sig.inx <- imp.inx;
  p.topInx <- GetTopInx(df$logp, 5, T) & (df$sig == "A");
  fc.leftInx <- GetTopInx(df$logfc, 5, F);
  lblInx.up <-  sig.inx & (p.topInx | fc.leftInx);
  text.lbls.up <- rownames(rna_results)[lblInx.up]
  
  p.topInx <- GetTopInx(df$logp, 5, T) & (df$sig == "B");
  fc.rtInx <- GetTopInx(df$logfc, 5, T);
  lblInx.down <- sig.inx & (p.topInx | fc.rtInx);
  text.lbls.dwn <- rownames(rna_results)[lblInx.down]
  
  text.lbls.all <- c(text.lbls.dwn, text.lbls.up)
  indices <- which(df$id %in% text.lbls.all)
  
  df$siglabels <- ifelse(df$id %in% text.lbls.all, "TRUE", "FALSE")
  
  if(colpal=="default"){
    cols <- c("red", "blue", "lightgrey")
  }else if(colpal=="pog"){ # orange, purple, grey
    cols <- c("#FA9E3BFF", "#47039FFF", "lightgrey")
  }else if(colpal=="bog"){ # orange, blue, grey
    cols <- c("#ED7953FF", "#1F9E89FF", "lightgrey")
  }else if(colpal=="rgg"){ # red, green, grey
    cols <- c("#F04C3BFF", "#1D9A6CFF", "lightgrey")
  }
  
  plotname <- paste(imgName, ".", format, sep="")
  Cairo::Cairo(file=plotname, width=550, height=550, type=format, bg="white", dpi=dpi);
  mbSetObj$imgSet$volcano <- plotname;
  
  # plot volcano
  p <- ggplot(df, aes(x = logfc, y = logp)) + geom_point(aes(color=sig), size = 3.5, alpha=0.75) + 
    scale_color_manual(values = c("A" = cols[1], "B" = cols[2], "C" = cols[3])) +
    theme_bw() + labs(x = "\nLog2 Fold Change", y = "-Log10 P-Value") + theme(legend.position = "none") + 
    theme(axis.title = element_text(size = 13), axis.text = element_text(size = 11)) 
  
  if(.on.public.web){
    p <- p + geom_text(data = subset(df, siglabels == TRUE), aes(label = subset(df, siglabels == TRUE)[,'id']), angle = 45, hjust = "inward", vjust="inward", check_overlap = TRUE)
  }else{
    p <- p + geom_text_repel(data = subset(df, siglabels == TRUE), aes(label = subset(df, siglabels == TRUE)[,'id']))
  }
  
  print(p);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj))
  
}

#'Plot fold-change dot plot
#'@description This functions creates a dot plot of the RNAseq analysis 
#'results from the microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param top.inx Numeric, input the number of features to include in the dot plot.
#'@param colpal Character, "default", "viridis" for the 
#'default viridis color palette, "plasma"
#'for the plasma color palette, and "cividis" for the  
#'cividis color palette.
#'@param imgName Character, input the name of the plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@import ggplot2
#'@import viridis
PlotRNASeqDotPlot <- function(mbSetObj, top.inx = 15, colpal = "plasma", imgName="rnaseq_dot",
                              format = "png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  if(mbSetObj$module.type == "sdp"){
    return(0)
  }
  
  taxrank <- taxrank
  rna_results <- mbSetObj$analSet$rnaseq$resTable
  
  # results to plot
  feats <- rna_results[1:top.inx,]
  cols <- colnames(feats)
  feats <- cbind(rownames(feats), feats)
  colnames(feats) <- c(taxrank, cols) 
  
  # create taxonomy table
  if(taxrank != "OTU"){
    # dirty fix to deal w. multiple matching taxa showing up
    glom <- fast_tax_glom_mem(mbSetObj$dataSet$proc.phyobj, taxrank);
    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
    taxa_table <- as(tax_table(glom), "matrix")
    t.first <- taxa_table[match(unique(taxa_table[,taxrank]), taxa_table[,taxrank]),]
    
    matches <- which(t.first[,taxrank] %in% rownames(feats))
    taxa_table_sub <- t.first[matches,]
    merged_table <- merge(feats, taxa_table_sub, taxrank)
  }else{
    glom <- mbSetObj$dataSet$proc.phyobj
    t.first <- as(tax_table(glom), "matrix")
    taxa <- rownames(t.first)
    oldcolnames <- colnames(t.first)
    t.first.merged <- cbind(taxa, t.first)
    colnames(t.first.merged) <- c("OTU", oldcolnames)
    merged_table <- merge(feats, t.first.merged, "OTU")
  }
  
  merged_table[,1] <- factor(merged_table[,1], levels = merged_table[,1])
  merged_table[,1] <- strtrim(merged_table[,1], 20)
  merged_table <- merged_table[order(merged_table$FDR),]
  
  size <- length(rownames(feats))
  
  if(size <= 20){
    h <- 600
    w <- 600
  }else if(size <= 35){
    h <- 700
    w <- 600
  }else{
    h <- 800
    w <- 650
  }
  
  plotname <- paste(imgName, ".", format, sep="")
  Cairo::Cairo(file=plotname, width=w, height=h, type=format, bg="white", dpi=dpi);
  mbSetObj$imgSet$dotplot <- plotname;
  
  p <- ggplot(merged_table, aes(x = merged_table[,1], y = log2FC, color = Phylum)) + geom_point(size = 6, alpha = 0.75) +
    theme_bw() + theme(axis.text.x = element_text(angle = -90, hjust = 0, vjust = 0.5, size = 11), axis.text.y = element_text(size = 11)) +
    labs(y = "\nLog2 Fold Change", x = paste(taxrank)) + coord_flip() + 
    theme(legend.text = element_text(size = 11), axis.title = element_text(size = 13), legend.title = element_text(size = 13));
  
  if(colpal == "viridis"){
    p <- p + viridis::scale_color_viridis(discrete = TRUE);
  }else if(colpal == "plasma"){
    p <- p + viridis::scale_color_viridis(discrete = TRUE, option="C");
  }else if(colpal == "cividis"){
    p <- p + viridis::scale_color_viridis(discrete = TRUE, option="E");
  }
  
  print(p)
  dev.off()
  
  return(.set.mbSetObj(mbSetObj))
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
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
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
  #making boxplot data
  boxdata <- as.data.frame(data1);
  mbSetObj$analSet$boxdata <- boxdata;
  cbtempl.results <- apply(data1, 2, template.match, feat_data, dist.name);
  cor.res <- t(cbtempl.results);
  fdr.col <- p.adjust(cor.res[,3], "fdr");
  cor.res <- cbind(cor.res, fdr.col);
  colnames(cor.res) <- c("correlation", "t-stat", "p-value", "FDR");
  ord.inx <- order(cor.res[,3])
  sig.mat <-signif(cor.res[ord.inx,],5);
  fileName <- "correlation_feature.csv";
  write.csv(sig.mat,file=fileName);
  mbSetObj$analSet$resTable <- as.data.frame(sig.mat);
  #removing Inf values from table
  is.na(mbSetObj$analSet$resTable) <- sapply(mbSetObj$analSet$resTable, is.infinite);
  mbSetObj$analSet$resTable[is.na(mbSetObj$analSet$resTable)]<-0;
  
  mbSetObj$analSet$pattern <- feat;
  mbSetObj$analSet$taxrank <- taxrank;
  
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

PlotCorr <- function(mbSetObj, imgName, format="png", dpi=72, width=NA){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  cor.res <- mbSetObj$analSet$cor.mat;
  pattern <- mbSetObj$analSet$pattern;
  title <- paste(mbSetObj$analSet$taxrank, "correlated with the", pattern);
  
  if(nrow(cor.res) > 25){
    cor.res <- cor.res[1:25, ];
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
  title <- paste("Top",nrow(cor.res), tolower(mbSetObj$analSet$taxrank), "correlated with", pattern);
  imgName = paste(imgName, ".", format, sep="");
  mbSetObj$imgSet$cor.ph <- imgName;
  
  if(is.na(width)){
    w <- h <- 7.2;
  }else if(width == 0){
    w <- 7.2;
  }else{
    w <- h <- width;
  }
  
  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(mar=c(5,6,4,3))
  rownames(cor.res)<-substr(rownames(cor.res), 1, 18);
  cols <- ifelse(cor.res[,1] >0, "mistyrose","lightblue");
  dotchart(cor.res[,1], pch="", xlim=c(-1,1), xlab="Correlation Coefficients", main=title);
  rownames(cor.res) <- NULL;
  barplot(cor.res[,1], space=c(0.5, rep(0, nrow(cor.res)-1)), xlim=c(-1,1), xaxt="n", col = cols, add=T,horiz=T);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj))
}

#'Match Patterns Function
#'@description This functions matches patterns (when
#'the defined pattern is set to a predefined profile or custom
#'profile) in the Pattern Search analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dist.name Input the distance measure. Input "pearson" for 
#'Pearson R, "spearman" for Spearman rank correlation, and "kendall"
#'for Kendall rank correlation.
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

Match.Pattern <- function(mbSetObj, dist.name="pearson", pattern=NULL, taxrank, variable){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(!.on.public.web){
    clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  }
  
  if(is.null(pattern)){
    pattern <- paste(1:length(levels(clslbl)), collapse="-");
  }
  
  templ <- as.numeric(ClearStrings(strsplit(pattern, "-")[[1]]));
  
  if(all(templ==templ[1])){
    AddErrMsg("Cannot calculate correlation on constant values!");
    return(0);
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
  
  if(mbSetObj$module.type=="mdp"){
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- readRDS("phyloseq_objs.RDS")
    }
    
    if(taxrank=="OTU"){
      data <- phyloseq_objs$count_tables$OTU
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data <- phyloseq_objs$count_tables[[taxrank.inx]]
    }
  }else{
    taxrank <- "OTU";
    data <- as.matrix(otu_table(mbSetObj$dataSet$norm.phyobj));
  }
  
  data <- t(data);
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  boxdata <- as.data.frame(data);
  boxdata$class <- clslbl;
  mbSetObj$analSet$boxdata <- boxdata;
  cbtempl.results <- apply(data, 2, template.match, new.template, dist.name);
  cor.res <- t(cbtempl.results);
  fdr.col <- p.adjust(cor.res[,3], "fdr");
  cor.res <- cbind(cor.res, fdr.col);
  colnames(cor.res) <- c("correlation", "t-stat", "p-value", "FDR");
  ord.inx <- order(cor.res[,3]);
  sig.mat <- signif(cor.res[ord.inx,],5);
  fileName <- "correlation_feature.csv";
  write.csv(sig.mat, file=fileName);
  mbSetObj$analSet$resTable <- as.data.frame(sig.mat);
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
  mbSetObj$analSet$cor.mat<-sig.mat;
  
  return(.set.mbSetObj(mbSetObj));
}

#'Function to create correlation heat map
#'@description This function creates the correlation heat map
#'@param mbSetObj Input the name of the mbSetObj.
#'@param imgName Character, input the name
#'of the plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param width Numeric, input the width of the plot. By
#'default it is set to NA.
#'@param cor.method Input the distance measure. Input "pearson" for 
#'Pearson R, "spearman" for Spearman rank correlation, and "kendall"
#'for Kendall rank correlation.
#'@param colors_cntrst Set the colors of the heatmap. By default it 
#'is set to "bwm", blue, white, to red. Use "gbr" for green, black, red, use
#'"heat" for red to yellow, "topo" for blue to yellow, "gray" for 
#'white to black, and "byr" for blue, yellow, red. Use "viridis"
#'for the viridis color palette and "plasma" for the plasma color palette from
#'the viridis R package.
#'@param viewOpt Character, "overview" to view an overview
#'of the heatmap, and "detail" to iew a detailed view of the
#'heatmap (< 1500 features).
#'@param taxrank Character, input the taxonomic level to perform
#'classification. For instance, "OTU-level" to use OTUs.
#'@param fix.col Logical, default set to FALSE.
#'@param no.clst Logical, default set to FALSE.
#'@param top Logical, default set to FALSE.
#'@param topNum Numeric, set the number of top features
#'to include.
#'@param datatype If the data is marker gene, input "16S". If the
#'data is metagenomic, use "metageno".
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import RColorBrewer
#'@import gplots
#'@import pheatmap

PlotCorrHeatMap <- function(mbSetObj, imgName, format="png", cor.method,
                            colors_cntrst, viewOpt, taxrank, 
                            topNum, permNum=100, pvalCutoff=0.05, corrCutoff=0.3){
  
  
  fix.col="F"; 
  no.clst="F"; 
  top="F"
  width ="NA"
  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_viridis();
  
  main <- xlab <- ylab <- NULL;
  
  if(mbSetObj$module.type=="sdp"){
    data <- mbSetObj$dataSet$norm.phyobj;
    data1 <- as.matrix(otu_table(data));
    taxrank <- "OTU";
  }else if(mbSetObj$module.type=="mdp"){
    if(taxrank=="OTU"){
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      data1 <- as.matrix(otu_table(data));
    }else{
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      #merging at taxonomy levels
      data <- fast_tax_glom_mem(data,taxrank);
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
    }
  }
  
  data <- t(data1);
  mbSetObj$analSet$abund_data<-data;
  
  if(ncol(data) > 1000){
    filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
    rk <- rank(-filter.val, ties.method='random');
    data <- as.data.frame(data[,rk <=1000]);
  }
  
  colnames(data)<-substr(colnames(data), 1, 18);
  
  if(cor.method=="sparcc"){
    #corr.mat <- RunFastSpar(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "heatmap")
    if(.on.public.web){
        corr.mat <- RunFastSpar_mem(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "heatmap")
    }else{
        corr.mat <- RunFastSpar(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "heatmap")
    }
  }else{
    corr.mat<-cor(data, method=cor.method);
  }
  
  # use total abs(correlation) to select
  if(top){
    cor.sum <- apply(abs(corr.mat), 1, sum);
    cor.rk <- rank(-cor.sum);
    var.sel <- cor.rk <= topNum;
    corr.mat <- corr.mat[var.sel, var.sel];
  }
  
  # set up parameter for heatmap
  load_rcolorbrewer();
  load_gplots();
  
  if(colors_cntrst=="gbr"){
    colors <- grDevices::colorRampPalette(c("green", "black", "red"), space="rgb")(256);
  }else if(colors_cntrst == "heat"){
    colors <- grDevices::heat.colors(256);
  }else if(colors_cntrst == "topo"){
    colors <- grDevices::topo.colors(256);
  }else if(colors_cntrst == "gray"){
    colors <- grDevices::colorRampPalette(c("grey90", "grey10"))(256);
  }else if(colors_cntrst == "byr"){
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256));
  }else if(colors_cntrst == "viridis") {
    colors <- rev(viridis::viridis(10))
  }else if(colors_cntrst == "plasma") {
    colors <- rev(viridis::plasma(10))
  }else{
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256));
  }
  
  #used for network edges
  first_color <- colors[1];
  last_color <- colors[length(colors)];
  imgName = paste(imgName, ".", format, sep="");
  
  if(viewOpt == "overview"){
    if(is.na(width)){
      w <- 9;
    }else if(width == 0){
      w <- 7.2;
      mbSetObj$imgSet$heatmap <-imgName;
    }else{
      w <- 7.2;
    }
    h <- w;
  }else{
    if(ncol(corr.mat) > 50){
      myH <- ncol(corr.mat)*12 + 40;
    }else if(ncol(corr.mat) > 20){
      myH <- ncol(corr.mat)*12 + 60;
    }else{
      myH <- ncol(corr.mat)*12 + 120;
    }
    h <- round(myH/72,2);
    
    if(is.na(width)){
      w <- h;
    }else if(width == 0){
      w <- h <- 7.2;
      mbSetObj$imgSet$corr.heatmap <-imgName;
    }else{
      w <- h <- 7.2;
    }
    
    # to prevent too small
    min.w <- 4.8;
    
    if(w < min.w){
      w <- h <- min.w;
    }
  }
  
  if(format=="pdf"){
    grDevices::pdf(file = imgName, width=w, height=h, bg="white", onefile=FALSE);
  }else{
    #Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
    Cairo::Cairo(file = imgName, unit="in", dpi=72, width=w, height=h, type=format, bg="white");
  }
  
  mbSetObj$imgSet$cor.heat<-imgName;
  
  if(no.clst){
    rowv = FALSE;
    colv = FALSE;
    dendro = "none";
  }else{
    rowv = TRUE;
    colv = TRUE;
    dendro = "both";
  }
  
  load_pheatmap();
  
  if(fix.col){
    breaks <- seq(from = -1, to = 1, length = 257);
    pheatmap::pheatmap(corr.mat, fontsize=8, fontsize_row=8, cluster_rows = colv,
                       cluster_cols = rowv, color = colors, breaks = breaks);
  }else{
    pheatmap::pheatmap(corr.mat, fontsize=8, fontsize_row=8, cluster_rows = colv,
                       cluster_cols = rowv, color = colors);
  }
  
  dev.off();
  
  mbSetObj$analSet$cor.heatmat<-corr.mat;
  mbSetObj$analSet$colors_cntrst<-colors_cntrst;
  mbSetObj$analSet$edge_col_low<-first_color;
  mbSetObj$analSet$edge_col_high<-last_color;
  mbSetObj$analSet$colors<-colors;
  mbSetObj$analSet$corheat.taxalvl<-taxrank;
  mbSetObj$analSet$corheat.meth<-cor.method;
  write.csv(signif(corr.mat,5), file="correlation_table.csv");
  
  if(cor.method=="sparcc"){
    method = "SparCC"
  }else{
    method = "Correlation"
  }
  
  current.msg <<- paste(method, "performed successfully!") 
  return(.set.mbSetObj(mbSetObj));
}

#'Function to call for correlation network
#'@description This function runs the fastspar or cor.test 
#'@param mbSetObj Input the name of the mbSetObj.
#'@import ppcor
#'@import igraph
#'@export
PerformNetworkCorrelation <- function(mbSetObj, taxrank, cor.method="pearson", colorOpt="expr", 
                                      permNum=100, pvalCutoff=0.05, corrCutoff=0.3){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(.on.public.web){
    load_ppcor();
    load_igraph();
  }

  mbSetObj$analSet$corr_color_opt <- colorOpt;
  
  if(cor.method == "sparcc"){
    if(.on.public.web){
        sparcc_results <- RunFastSpar_mem(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network")
    }else{
        sparcc_results <- RunFastSpar(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network")
    }
  
    if(nrow(sparcc_results)==0){
      AddErrMsg("No correlations meet the p-value and correlation thresholds!")
      return(0)
    }else{
      mbSetObj$analSet$network_cor <- sparcc_results
    }
  }else{
    
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- readRDS("phyloseq_objs.RDS")
    }
    
    if(taxrank=="OTU"){
      data1 <- phyloseq_objs$count_tables$OTU
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- phyloseq_objs$count_tables[[taxrank.inx]]
    }
    
    data <- t(data1);
    data <- data[which(rownames(data) %in% mbSetObj$dataSet$selected.grps),]
    data[data==0|is.na(data)] <- .00001
    
    if(ncol(data) > 1000){
      filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
      rk <- rank(-filter.val, ties.method='random');
      data <- as.data.frame(data[,rk <=1000]);
    }
    
      mbSetObj$analSet$netcorr_data<-data;
      vars = data.frame(t(combn(colnames(data), 2)), stringsAsFactors = FALSE)
      
      if(cor.method == "spearman"){
        cor.results <- do.call(rbind, mapply(SpearmanCorrFunc, vars[,1], vars[,2], MoreArgs = list(data=data), SIMPLIFY = FALSE))
      }else if(cor.method == "kendall"){
        cor.results <- do.call(rbind, mapply(KendallCorrFunc, vars[,1], vars[,2], MoreArgs = list(data=data), SIMPLIFY = FALSE))
      }else if(cor.method == "pearson"){
        cor.results <- do.call(rbind, mapply(PearsonCorrFunc, vars[,1], vars[,2], MoreArgs = list(data=data), SIMPLIFY = FALSE))
      }else{
        AddErrMsg("Invalid correlation method!")
        return(0)
      }
      
      cor.results.filt <- cor.results[(abs(cor.results[,3]) > corrCutoff & cor.results[,4] < pvalCutoff),]
      colnames(cor.results.filt) <- c("Taxon1", "Taxon2", "Correlation", "P.value", "Statistic", "Method")
      cor.results.filt[,3] <- round(cor.results.filt[,3], digits=4)
      cor.results.filt[,4] <- round(cor.results.filt[,4], digits=4)
      write.csv(cor.results.filt, "correlation_table.csv", row.names = FALSE)
      mbSetObj$analSet$network_cor <- cor.results.filt;
    
  }
  
  #network building only needed for web
  if(.on.public.web){
    all.taxons = unique(c(mbSetObj$analSet$network_cor[,1] , mbSetObj$analSet$network_cor[,2]))
    taxColNms = vector();
    
    if(taxrank == "OTU"){
      tax.tbl = as.data.frame(matrix(mbSetObj$dataSet$taxa_table[,1:7],ncol=7))
      colnames(tax.tbl) = colnames(mbSetObj$dataSet$taxa_table[,1:7])
      colnames(tax.tbl)[7] = "OTU"
      tax.tbl[,7]=rownames(mbSetObj$dataSet$taxa_table)
      taxColNms = colnames(tax.tbl)
      tax.tbl = as.data.frame(tax.tbl[which(tax.tbl[,simpleCap(taxrank)] %in% all.taxons),])
      colnames(tax.tbl) = taxColNms
    }else{
      tax.tbl = data.frame(mbSetObj$dataSet$taxa_table[,1:which( colnames(mbSetObj$dataSet$taxa_table)== simpleCap(taxrank))])
      taxColNms = colnames(mbSetObj$dataSet$taxa_table[,1:which( colnames(mbSetObj$dataSet$taxa_table)== simpleCap(taxrank))])
      colnames(tax.tbl) = taxColNms
      tax.tbl = as.data.frame(tax.tbl[which(tax.tbl[,simpleCap(taxrank)] %in% all.taxons),])
      colnames(tax.tbl) = taxColNms
    }
    
    inx = !duplicated(tax.tbl[,simpleCap(taxrank)])
    filt.taxa.table = data.frame(tax.tbl[inx,])
    colnames(filt.taxa.table) = taxColNms
    mbSetObj$analSet$filt.taxa.table = filt.taxa.table 
    .set.mbSetObj(mbSetObj);
    res <- SparccToNet(mbSetObj)
    if(res == 0){
      AddErrMsg("Errors during creating correlation network!")
      return(0);
    }
  }
  
  if(cor.method=="sparcc"){
    method = "SparCC"
  }else{
    method = paste(cor.method, "correlation"); 
  }
  
  # calculate microbial dysbiosis index
  feats_used <- unique(mbSetObj$analSet$network_cor[,1], mbSetObj$analSet$network_cor[,2]) 
  fc <- mbSetObj$analSet$diff_table
  taxa_lvl <- mbSetObj$analSet$network_taxalvl
  clean_feats_used <- sub("^[^_]*__", "", unique(feats_used))
  match.inx <- which(fc$tax_name %in% clean_feats_used)
  
  if(length(match.inx)!=0){
    fc.filt <- fc[match.inx,]
    
    if(length(unique(fc.filt$tax_name)) < length(fc.filt$tax_name)){
      # need to delete doubles in fc.filt
      keep <- c(TRUE, head(fc.filt$tax_name, -1) != tail(fc.filt$tax_name, -1))
      fc.filt <- fc.filt[keep,]
    }
    
    # need to filter again
    if(cor.method=="sparcc"){
      abund_data <- readRDS("sparcc_data.RDS")
    }else{
      abund_data <- t(as.matrix(mbSetObj$analSet$netcorr_data))
    }
    
    abund.inx <- which(rownames(abund_data) %in% feats_used)
    abund_data.filt <- abund_data[abund.inx,]
    
    log2fc <- fc.filt$log2_median_ratio
    names(log2fc) <- fc.filt$tax_name
    
    inc <- which(log2fc >= 0) # total abundance increased in X
    dec <- which(log2fc < 0) # total abundance decreased in X
    
    if(length(inc)==0 | length(dec)==0){
      mbSetObj$analSet$mdi <- "MD-index could not be calculated."
    }else{
      abund <- c(sum(abund_data.filt[inc,]), sum(abund_data.filt[dec,]))
      abund[abund==0] <- 0.1
      
      MDI1 <- round(log(abund[1]/abund[2]), digits = 4) 
      first <- unique(fc.filt$treatment_1)
      second <- unique(fc.filt$treatment_2)
      groups <- paste(first, "/", second, sep="")
      
      mbSetObj$analSet$mdi <- paste(groups, "MD-index:", MDI1)
    }
  }
  
  if(current.msg == "Only the top 500 features are kept, ranked by their variance!"){
    current.msg <<- paste(current.msg, method, "performed successfully!")
  }else{
    current.msg <<- paste(method, "performed successfully!")
  }
  
  return(.set.mbSetObj(mbSetObj));
}

###################################
########## Util Functions #########
###################################

# use memoise tech to avoid repeating
RunFastSpar_mem <- function(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output="network", palette="viridis"){
  if(!exists("fastsparr_mem")){
    require("memoise");
    fastsparr_mem<<- memoise(RunFastSpar); 
  }
  
  res <- fastsparr_mem(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output="network", palette="viridis");
  return(res);
}

#'Function to call FastSpar
#'@description This function runs the fastspar 
#'@param mbSetObj Input the name of the mbSetObj.
#'@import igraph 
#'@export
RunFastSpar <- function(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output="network", palette="viridis"){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_igraph();
  }
  
  path_fastspar <- "";
  path_fastspar_bs <- "";
  path_fastspar_pvals <- "";
  
  # for development purposes
  if(.on.public.web){ 
    if(file.exists("/home/glassfish/libraries/fastspar")){ #server
      path_fastspar <- "/home/glassfish/libraries/fastspar/fastspar"
      path_fastspar_bs <- "/home/glassfish/libraries/fastspar/fastspar_bootstrap"
      path_fastspar_pvals <- "/home/glassfish/libraries/fastspar/fastspar_pvalues"
    }else if(file.exists("~/Downloads/fastspar-0.0.10_linux/fastspar")){ #jasmine workstation
      path_fastspar <- "~/Downloads/fastspar-0.0.10_linux/fastspar"
      path_fastspar_bs <- "~/Downloads/fastspar-0.0.10_linux/fastspar_bootstrap"
      path_fastspar_pvals <- "~/Downloads/fastspar-0.0.10_linux/fastspar_pvalues"
    }else{
      AddErrMsg("Cannot find path to fastspar! Please intall fastspar and add path here!");
      return(0);
    }
  }else{
    if(file.exists("~/fastspar")){ # need to adapt per user
      path_fastspar <- "~/Downloads/fastspar-0.0.10_linux/fastspar"
      path_fastspar_bs <- "~/Downloads/fastspar-0.0.10_linux/fastspar_bootstrap"
      path_fastspar_pvals <- "~/Downloads/fastspar-0.0.10_linux/fastspar_pvalues"
    }else{
      AddErrMsg("Cannot find path to fastspar! If fastspar is not downloaded, please install!");
      return(0);
    }
  }
  
  # first get OTU table to BIOM TSV format
  if(mbSetObj$module.type=="sdp"){
    data <- mbSetObj$dataSet$norm.phyobj;
    data1 <- as.matrix(otu_table(data));
    taxrank <- "OTU";
  }else if(mbSetObj$module.type=="mdp"){
    
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- readRDS("phyloseq_objs.RDS")
    }
    
    if(taxrank=="OTU"){
      data1 <- phyloseq_objs$count_tables$OTU
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- phyloseq_objs$count_tables[[taxrank.inx]]
    }
  }
  
  if(!is.null(mbSetObj$dataSet$selected.grps)){
    data1 <- data1[,which(colnames(data1) %in% mbSetObj$dataSet$selected.grps)]
  }
  
  # filter data if necessary
  abund_filtered <- mbSetObj$dataSet$ab.filtered
  var_filtered <- mbSetObj$dataSet$var.filtered
  
  feat_num <- nrow(data1)
  
  # if number of feats is greater than 500, do filtering of poorly represented OTUs
  # first, do low count filtering 
  if(feat_num > 500 & !abund_filtered){
    minLen <- 0.2*ncol(data1); # filter out feats that do not have a minimum count of 4 in at least 20 percent of samples
    kept.inx <- apply(data1, MARGIN = 1,function(x) {sum(x >= 4) >= minLen});
    data1 <- data1[kept.inx, ];
  }
  
  # second, do low variance filtering
  if(feat_num > 500 & !var_filtered){
    filter.val <- apply(data1, 1, IQR, na.rm=T);
    rk <- rank(-filter.val, ties.method='random');
    remain <- rk < feat_num*(1-0.1); # filter out 10% of low variance feats based on IQR
    data1 <- data1[remain,];
  }
  
  # third, if still over 500 feats, rank and keep top 500
  if(feat_num > 500){
    filter.val <- apply(data1, 1, IQR, na.rm=T);
    rk <- rank(-filter.val, ties.method='random');
    remain <- rk < 500;
    data1 <- data1[remain,];
    current.msg <<- "Only the top 500 features are kept, ranked by their variance!"
  }
  
  # replace zeros with random numbers based on lowest non-zero count
  
  set.seed(12345)
  lowest.count <- min(apply(data1, 1, FUN = function(x) {min(x[x > 0])}))
  
  # replacements
  replacements <- apply(data1, 1, function(x) { (sample(1:lowest.count, size=sum(x==0), replace = F))/lowest.count} )
  
  # this works for rows
  replaced <- vector(length = nrow(data1), "list")
  
  for(i in 1:nrow(data1)){
    replaced[[i]] <- replace(as.vector(data1[i,]), as.vector(data1[i,]==0), replacements[[i]])
  }
  
  zero.output <- matrix(unlist(replaced), ncol = ncol(data1), byrow = TRUE)
  colnames(zero.output) <- colnames(data1)
  rownames(zero.output) <- rownames(data1)
  
  saveRDS(zero.output, "sparcc_data.RDS")
  
  # add header
  new_col <- c("#OTU ID", colnames(zero.output))
  data <- cbind(rownames(zero.output), zero.output)
  colnames(data) <- new_col
  
  write.table(data, file="otu_table_corr.tsv", quote=FALSE, row.names = FALSE, sep="\t")
  
  otu_table <- "otu_table_corr.tsv"
  corr_output <- "sparcc_median_correlation.tsv"
  cov_output <- "sparcc_median_covariance.tsv"
  bootstrap_counts <- "bootstrap_counts"
  counts_prefix <- "bootstrap_counts/boot_data"
  bootstrap_correlation <- "bootstrap_correlation"
  corr_bs_output <- "${jnew}_median_correlation.tsv"
  cov_bs_output <- "${k}_median_covariance.tsv"
  corr_prefix <- "cor_boot_data_"
  test <- "text.txt"
  pval_output <- "sparcc_pvals.tsv"
  
  system(paste(path_fastspar, "--otu_table", otu_table, "--correlation", corr_output, "--covariance", cov_output, ";",
               "mkdir", bootstrap_counts, ";",
               path_fastspar_bs,  "--otu_table", otu_table, "--number", permNum, "--prefix", counts_prefix, ";",
               "find", bootstrap_counts, ">>", test, ";",
               "cat", test, "| while read line; do echo $line; j=$(basename ${line}); jnew=$(echo cor_${j}); k=$(echo cov_${j}); jnew=$(echo ${jnew} | sed 's/\\.tsv//'); k=$(echo ${k} | sed 's/\\.tsv//');", path_fastspar, "--otu_table ${line} --correlation",
               corr_bs_output, "--covariance", cov_bs_output, "-i 5 ; done;",
               path_fastspar_pvals,  "--otu_table", otu_table, "--correlation", corr_output, "--prefix", corr_prefix, "--permutations", permNum, "--outfile", pval_output, ";",
               "rm cor_boot_data*.tsv; rm cov_boot_data*.tsv;", "rm -rf", bootstrap_counts, "rm -rf text.txt"))
  
  sparcc_results <- read.table(file = corr_output, sep="\t", stringsAsFactors = FALSE)
  names <- sparcc_results[,1]
  sparcc_results_new <- as.matrix(sparcc_results[,-1])
  colnames(sparcc_results_new) <- rownames(sparcc_results_new) <- names
  
  sparcc_pvals <- read.table(file = pval_output, sep="\t", stringsAsFactors = FALSE)
  sparcc_pvals_new <- as.matrix(sparcc_pvals[,-1])
  colnames(sparcc_pvals_new) <- rownames(sparcc_pvals_new) <- names
  
  if(output == "network"){
    sparcc_results_new[sparcc_results_new==0] <- .00001
    sparcc_corr <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(sparcc_results_new, weighted = TRUE))
    sparcc_pvals_new[sparcc_pvals_new==0] <- .00001
    sparcc_pvals <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(sparcc_pvals_new, weighted = TRUE))
    
    sparcc_combo <- cbind(sparcc_corr, sparcc_pvals[,3])
    sparcc_combo <- sparcc_combo[(abs(sparcc_combo[,3]) > corrCutoff & sparcc_combo[,4] < pvalCutoff),]
    colnames(sparcc_combo) <- c("Taxon1", "Taxon2", "Correlation", "P.Value")
    write.csv(sparcc_combo, "correlation_table.csv", row.names = FALSE)
    .set.mbSetObj(mbSetObj)
    return(sparcc_combo)
  }else if(output=="heatmap"){
    .set.mbSetObj(mbSetObj)
    return(sparcc_results_new)
  }
}

# Set of functions to perform cor.test
PearsonCorrFunc <- function(var1, var2, data){
  result <- cor.test(data[,var1], data[,var2])
  data.frame(var1, var2, result[c("estimate", "p.value", "statistic", "method")], stringsAsFactors = FALSE)
}

SpearmanCorrFunc <- function(var1, var2, data){
  result <- cor.test(data[,var1], data[,var2], method = "spearman", exact = FALSE)
  data.frame(var1, var2, result[c("estimate", "p.value", "statistic", "method")], stringsAsFactors = FALSE)
}

KendallCorrFunc <- function(var1, var2, data){
  result <- cor.test(data[,var1], data[,var2], method = "kendall", exact = FALSE)
  data.frame(var1, var2, result[c("estimate", "p.value", "statistic", "method")], stringsAsFactors = FALSE)
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
    group = get_variable(physeq, group)
  }
  # Define gene annotations (`genes`) as tax_table
  taxonomy = tax_table(physeq, errorIfNULL=FALSE)
  if( !is.null(taxonomy) ){
    taxonomy = data.frame(as(taxonomy, "matrix"))
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
    ADF = Biobase::AnnotatedDataFrame(data.frame(sample_data(physeq)))
  }else {
    ADF = NULL
  }
  taxonomy = tax_table(physeq, errorIfNULL=FALSE);
  if (!is.null(taxonomy)) {
    TDF = Biobase::AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                                 data.frame(as(taxonomy, "matrix")), row.names = taxa_names(physeq)))
    #  data.frame(tax_table(physeq)@.Data), row.names = taxa_names(physeq)))
    #  data.frame(tax_table(physeq)), row.names = taxa_names(physeq))) ## bugs with tax_table
  }else {
    TDF = Biobase::AnnotatedDataFrame(data.frame(OTUname = taxa_names(physeq), 
                                                 row.names = taxa_names(physeq)))
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
  inx <- order(vec, decreasing = dec)[1:n]
  vec <- rep(F, length = length(vec))
  vec[inx] <- T
  return(vec)
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
  GetSigTable(mbSetObj$analSet$Univar$resTable, "Univariate analysis");
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetSigTable.Corr<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  GetSigTable(mbSetObj$analSet$cor.mat, "Pattern search using correlation analysis");
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetSigTable.METAGENOSEQ<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  GetSigTable(mbSetObj$analSet$metagenoseq$resTable, "metagenomeSeq");
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetSigTable.RNASeq<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  GetSigTable(mbSetObj$analSet$rnaseq$resTable, "RNAseq method");
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetSigTable.LEFSE<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  GetSigTable(mbSetObj$analSet$lefse$resTable, "LEfSe");
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import xtable
GetRFConf.Table<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  load_xtable();
  print(xtable::xtable(mbSetObj$analSet$rf$confusion, caption="Random Forest Classification Performance"), size="\\scriptsize");
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import xtable
GetMapTable<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  load_xtable();
  print(xtable::xtable(mbSetObj$analSet$mapTable, caption="Result from Taxa Name Mapping"),
        tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import xtable
GetORATable<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  load_xtable();   
  res <- mbSetObj$analSet$ora.mat;
  print(xtable::xtable(res, caption="Result from Over Representation Analysis"),
        tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}

# utility method to get anova results
GetAovRes <- function(cls, data, nonpar=F){
  
  if(nonpar){
    anova.res <- apply(as.matrix(data), 2, function(x){kruskal.test(x ~ cls)});
    res <- unlist(lapply(anova.res, function(x) {c(x$statistic, x$p.value)}));
    return(matrix(res, nrow=length(anova.res), byrow=T));
  }else{
    aov.res <- apply(as.matrix(data), 2, function(x){aov(x ~ cls)});
    anova.res <- lapply(aov.res, anova);
    res <- unlist(lapply(anova.res, function(x) { c(x["F value"][1,], x["Pr(>F)"][1,])}));
    return(matrix(res, nrow=length(aov.res), byrow=T));
  }
}

# utility method to get p values
#'@import genefilter
GetTtestRes<- function(cls, data, paired=FALSE, equal.var=TRUE, nonpar=F){
  
  if(nonpar){
    inx1 <- which(cls==levels(cls)[1]);
    inx2 <- which(cls==levels(cls)[2]);
    res <- apply(data, 2, function(x) {
      tmp <- try(wilcox.test(x[inx1], x[inx2], paired = paired));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    })
  }else{
    if(ncol(data) < 1000){
      data <- na.omit(data);
      inx1 <- which(cls==levels(cls)[1]);
      inx2 <- which(cls==levels(cls)[2]);
      
      res <- apply(data, 2, function(x) {
        tmp <- try(t.test(x[inx1], x[inx2], paired = paired, var.equal = equal.var, silent=TRUE));
        if(class(tmp) == "try-error") {
          return(c(NA, NA));
        }else{
          return(c(tmp$statistic, tmp$p.value)); }})
    }else{ # use fast version
      load_genefilter();
      res <- try(rowttests((t(data)), cls));
      if(class(res) == "try-error") {
        res <- c(NA, NA);
      }else{
        res <- t(cbind(res$statistic, res$p.value));
      }
    }
  }
  return(t(res));
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

SparccToNet <- function(mbSetObj=NULL){
  library(igraph)
  mbSetObj <- .get.mbSetObj(mbSetObj);
  edge.list= mbSetObj$analSet$network_cor
  edge.list = edge.list[,c("Taxon1", "Taxon2"), drop=FALSE]
  g = graph_from_data_frame(edge.list, directed = FALSE, vertices = NULL)
  E(g)$weight = abs(mbSetObj$analSet$network_cor[, "Correlation"])
  E(g)$correlation = mbSetObj$analSet$network_cor[, "Correlation"]
  corr_matrix = get.adjacency(g,sparse=FALSE)
  mbSetObj$analSet$corr_matrix = corr_matrix
  
  colnames(edge.list) = c("Source", "Target")
  nodes = unique(c(edge.list[,1], edge.list[,2]))
  node.list = data.frame(Id=nodes, Name=nodes)
  overall.graph <- g
  
  ppi.net <- list(
    db.type="abc",
    order=1, 
    seeds=" ", 
    table.nm=" ", 
    node.data = node.list, 
    edge.data = edge.list,
    require.exp = FALSE,
    min.score = 400
  );
  
  ppi.net <<-ppi.net
  
  ppi.comps <- list();
  ppi.comps[["subnetwork"]] <- overall.graph;
  ppi.comps <<- ppi.comps
  mbSetObj <- .set.mbSetObj(mbSetObj);
  corr.net.name <- paste0("correlation_net_", corr.net.count, ".json")
  corr.net.name <<- corr.net.name
  return(PrepareNetworkExp(mbSetObj, "subnetwork",corr.net.name));
}


PrepareNetworkExp <- function(mbSetObj=NULL, net.nm, json.nm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$analSet$corr_net_nm <- json.nm
  if(is.null(net.nm)){
    net.nm <- names(ppi.comps)[1];
  }
  
  taxa = mbSetObj$analSet$filt.taxa.table
  
  my.ppi <- ppi.comps[[net.nm]];
  nd.nms <- V(my.ppi)$name;
  
  entrezIDs <- nd.nms
  names(entrezIDs) <- nd.nms;
  current.anot <<- entrezIDs;
  current.net.nm <<- net.nm; 
  
  return(convertIgraph2JSONExp(net.nm, json.nm, FALSE));
}


convertIgraph2JSONExp <- function(net.nm, filenm, thera=FALSE){
  g <- ppi.comps[[net.nm]];
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  abundance_table = mbSet$analSet$boxdatacor[,-length(colnames(mbSet$analSet$boxdatacor))]
  ab = apply(abundance_table,2,function(x){sum(x)/length(x)})  
  taxa = mbSetObj$analSet$filt.taxa.table
  
  colorOpt <- mbSetObj$analSet$corr_color_opt;
  if(!colorOpt %in% colnames(taxa) && colorOpt != "expr"){
    current.msg <<- "Invalid taxa is selected for coloring (must be same or higher taxonomy level of selected taxa used for correlation calculation)"
    return(0);
  }
  
  # annotation
  nms <- V(g)$name;
  inx <- !nms %in% taxa[,length(colnames(taxa))]
  toDelete = nms[inx]
  g <- delete_vertices(g, toDelete);
  nms <- V(g)$name;
  taxa=taxa[match(nms, taxa[,length(colnames(taxa))]),]
  hit.inx <- match(nms, ppi.net$node.data[,1]);
  lbls <- ppi.net$node.data[hit.inx, 2];
  
  # setup shape (gene circle, other squares)
  shapes <- rep("circle", length(nms));
  itypes <- rep("circle", length(nms));
  seeds <- rep("circle", length(nms));
  
  # get edge data
  edge.mat <- get.edgelist(g);

  edge.mat1 = data.frame(edge.mat)
  edge.mat1$color = ComputeColorGradientCorr(E(g)$correlation);
  edge.mat1 = as.matrix(edge.mat1)
  
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2], color = edge.mat1[,3], weight=E(g)$weight, correlation = E(g)$correlation);
  
  # now get coords
  pos.xy <- PerformLayOut(net.nm, "Default");
  # get the note data
  node.btw <- as.numeric(betweenness(g));
  node.eig <- eigen_centrality(g);
  node.eig = as.numeric(node.eig$vector);
  node.tra=transitivity(g,type=c("local"))
  node.dgr <- as.numeric(degree(g));
  node.exp <- as.numeric(get.vertex.attribute(g, name="abundance", index = V(g)));
  
  # node size to abundance values
  if(vcount(g)>500){
    min.size = 1;
  }else if(vcount(g)>200){
    min.size = 2;
  }else{
    min.size = 2;
  }
  
  #modules = FindCommunities("walktrap", FALSE);
  abundance = unname(ab[nms])
  node.sizes <- as.numeric(rescale2NewRange((log(abundance+1))^2, min.size, 15));
  
  centered = T;
  notcentered = F;
  # update node color based on betweenness
  require("RColorBrewer");
  topo.val <- log(node.btw+1);
  exp_table = mbSetObj$analSet$diff_table
  #log(exp_table$median_diff + min(exp_table$median_diff/2))
  colVecNms = rep("NA", length(topo.val))
  nms2=nms
  if(colorOpt == "expr"){
    nms1 = strsplit(nms, "__")
    if(length(nms1[[2]])>1){
    l = unlist(lapply(nms1, function(x) unlist(x[2])));
    inx = is.na(l)
if(length(which(inx == T))>0){
    l[inx]= paste0(nms1[[which(inx == T)]][1],"__")
}
    nms2=l
    }else{
    nms2 = nms
    }
    exp_table = exp_table[match(nms2, exp_table$tax_name), ]
    topo.colsb <- topo.colsb1 <-ComputeColorGradientCorr(exp_table$median_diff);
    topo.colsw <- topo.colsw1 <-ComputeColorGradientCorr(exp_table$median_diff);
  }else{
    color_var <- as.character(unique(taxa[,colorOpt]));
    x <- length(color_var);
    if(x<10){
      cols = brewer.pal(n = 9, name = "Set3")
      colVec <- rep(cols,length.out=x);
    }else if(x<21){
      colVec <- rep(custom_col21,length.out=x);
    }else{
      colVec <- rep(custom_col42,length.out=x);
    }
    names(colVec) = unique(taxa[,colorOpt])
    colVecNms<-names(colVec[as.character(taxa[,colorOpt])])
    topo.colsb  <- unname(colVec[as.character(taxa[,colorOpt])])
    topo.colsw  <- topo.colsb
  }
  #topo.colsb <- topo.colsb1 <- topo.val
  #topo.colsw <- topo.colsw1 <- topo.val
  
  # color based on expression
  bad.inx <- is.na(node.exp) | node.exp==0;
  if(!all(bad.inx)){
    exp.val <- node.exp;
    node.colsb.exp <- ComputeColorGradientCorr(exp.val); 
    node.colsw.exp <- ComputeColorGradientCorr(exp.val);
    node.colsb.exp[bad.inx] <- "#d3d3d3"; 
    node.colsw.exp[bad.inx] <- "#c6c6c6"; 
  }else{
    node.colsb.exp <- rep("#D3D3D3",length(node.dgr)); 
    node.colsw.exp <- rep("#C6C6C6",length(node.dgr)); 
  }
  
  # now update for bipartite network
  
  #topo.colsw <- rep("#FF8484", length(nms))
  #topo.colsb <- rep("#FF8484", length(nms))
  
  types_arr <<- c("gene", "interactor")
  
  shapes<- "gene";
  
  mat_source <<- edge.mat[,"source"]
  mat_target <<- edge.mat[,"target"]
  
  network_prop = list();
  for(i in 1:length(node.sizes)){
    network_prop[[i]]  <- list(
      eigen = node.eig[i],
      transitivity = node.tra[i]
    )
  }
  
  type <- rep(FALSE,length(node.dgr))
  
  lblsu <<- lbls;
  node_attr = list.vertex.attributes(g);
  node_attr = node_attr[which(node_attr!="names")] 
  
  # now create the json object
  nodes <- vector(mode="list");
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(
      id=nms[i], 
      idx=i,
      label=nms2[i],
      size=node.sizes[i], 
      tsize=node.sizes[i],
      type=shapes[i],
      itype=itypes[i],
      taxon=colVecNms[i],
      color=topo.colsb[i],
      colorb=topo.colsb[i],
      colorw=topo.colsw[i],
      topocolb=topo.colsb[i],
      topocolw=topo.colsw[i],
      expcolb=node.colsb.exp[i],
      expcolw=node.colsw.exp[i],
      #x = pos.xy[i,1],
      #y = pos.xy[i,2],
      x=pos.xy[i,1],
      y= pos.xy[i,2],
      user =network_prop[[i]],
      attributes=list(
        degree=node.dgr[i], 
        between=node.btw[i],
        expr = exp_table$mean_diff[i],
        eigen = node.eig[i],
        transitivity = node.tra[i]
      )
    );
  }
  summary = vector(mode="list");
  summary[[1]] = list(nodes=vcount(g));
  summary[[2]] = list(edges=ecount(g));
  #if(nrow(edge.mat) > 2*length(nms)){
  #    edge.mat <- edge.mat[-order(edge.mat[,5),]
  #    edge.mat <- edge.mat[c(1:length(nms)*2),]
  #}
  # save node table
  nd.tbl <- data.frame(Id=nms, Label=lbls, Degree=node.dgr, Betweenness=round(node.btw,2));
  # order 
  ord.inx <- order(nd.tbl[,3], nd.tbl[,4], decreasing = TRUE)
  nd.tbl <- nd.tbl[ord.inx, ];
  write.csv(nd.tbl, file="node_table.csv", row.names=FALSE);
  edge.matall <- edge.mat
  # covert to json
  require(RJSONIO);
  taxNodes = list(tax=colVecNms,colors=topo.colsb)
  mbSetObj <- .set.mbSetObj(mbSetObj);
  netData <- list(nodes=nodes, edges=edge.mat, taxa=taxa, taxaNodes = taxNodes)#, modules = modules, taxaNodes = list(tax=colVecNms,colors=topo.colsb));
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(1)
}

GetColorGradient <- function(background, center){
  if(background == "black"){
    if(center){
      return(c(colorRampPalette(c("#31A231", "#5BC85B", "#90EE90", "#C1FFC1"))(50), colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50)));
    }else{
      return(c(colorRampPalette(c("#edf8fb","#b2e2e2","#66c2a4","#2ca25f","#006d2c"))(100)));
    }
  }else{ # white background
    if(center){
      return(c(colorRampPalette(c("#137B13", "#31A231", "#5BC85B", "#90EE90"))(50), colorRampPalette(c("#FF7783", "#E32636", "#BD0313", "#96000D"))(50)));
    }else{
      return(colorRampPalette(hsv(h = seq(0.72, 1, 0.035), s = 0.72, v = 1))(100));
    }
  }
}

ComputeColorGradient <- function(nd.vec, background="black", centered){
  require("RColorBrewer");
  minval = min(nd.vec, na.rm=TRUE);
  maxval = max(nd.vec, na.rm=TRUE);
  res = maxval-minval;
  
  if(res == 0){
    return(rep("#0b6623", length(nd.vec)));
  }
  
  if(sum(nd.vec<0, na.rm=TRUE) > 0){ 
    centered <- T;
  }else{
    centered <- F;
  }
  color <- GetColorGradient(background, centered);
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}


PerformLayOut <- function(net.nm, algo){
  g <- ppi.comps[[net.nm]];
  vc <- vcount(g);
  if(algo == "Default"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else if(vc > 1000) {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }else if(vc < 150){
      pos.xy <- layout.kamada.kawai(g);
    }else{
      pos.xy <- layout.fruchterman.reingold(g);
    }
  }else if(algo == "FrR"){
    pos.xy <- layout.fruchterman.reingold(g);
  }else if(algo == "random"){
    pos.xy <- layout.random(g);
  }else if(algo == "lgl"){
    if(vc > 3000) {
      pos.xy <- layout.lgl(g, maxiter = 100);
    }else if(vc > 2000) {
      pos.xy <- layout.lgl(g, maxiter = 150);
    }else {
      pos.xy <- layout.lgl(g, maxiter = 200);
    }
  }else if(algo == "gopt"){
    # this is a slow one
    if(vc > 3000) {
      maxiter = 50;
    }else if(vc > 2000) {
      maxiter = 100;
    }else if(vc > 1000) {
      maxiter = 200;
    }else{
      maxiter = 500;
    }
    pos.xy <- layout.graphopt(g, niter=maxiter);
  }
  pos.xy;
}

GetCorrNetName <- function(){
  
  return(corr.net.name)
}

rescale2NewRange <- function(qvec, a, b){
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}

PrepareCorrExpValues <- function(mbSetObj, meta, taxalvl, color, layoutOpt, comparison, 
                                 wilcox.cutoff){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_metacoder();
  tax_o <- taxalvl;
  dm <- mbSetObj$dataSet$proc.phyobj;
  dims <- ncol(tax_table(dm))
  tax_table_new = data.frame("Kingdom" = "Root", as(tax_table(dm), "matrix")[, 1:dims]) # add root to tax table
  tax_table(dm) <- as.matrix(tax_table_new);
  
  dm_samples = as(sample_data(dm), "data.frame");
  dm_samples <- cbind.data.frame("sample_id" = row.names(dm_samples), dm_samples); # add sample_id to sample table
  row.names(dm_samples) <- c();
  dm_samples[] <- lapply(dm_samples, as.character);
  
  grp.nms <- strsplit(comparison, "_vs_")[[1]];
  dm_samples_cmf <- dm_samples[dm_samples[[meta]] %in% grp.nms, ]; #subset sample data by meta variable
  
  otu_dm <- as.data.frame(as(otu_table(dm), "matrix"));
  tax_dm <- as.data.frame(as(tax_table(dm), "matrix"));
  tax_dm[] <- lapply(tax_dm, as.character);#make sure characters in tax_dm;
  
  depth <- ncol(tax_dm)
  rank_dm <- c("r", "p", "c", "o", "f", "g", "s");
  names(tax_dm) <- rank_dm[1:depth];
  
  for (i in 1:ncol(tax_dm)){
    for (j in 1:nrow(tax_dm)){
      if (is.na(tax_dm[j, i])){
        tax_dm[j, i] <- "";
      } else {
        tax_dm[j, i] <- paste(names(tax_dm)[i],
                              tax_dm[j, i],
                              sep = "__");
      }
    }
  } #add __ to tax table
  
  if(taxalvl == "Phylum"){
    tax <- "p";
  } else if(taxalvl == "Class"){
    tax <- "c";
  } else if(taxalvl == "Order"){
    tax <- "o";
  } else if(taxalvl == "Family"){
    tax <- "f";
  } else if(taxalvl == "Genus"){
    tax <- "g"
  } else {
    tax <- "s";
  }; # get tax rank for heat tree
  
  tax_dm <- tax_dm[, 1:which(names(tax_dm) == tax)]; #subset tax table
  rank_dm_new <- rank_dm[1:which(rank_dm == tax)];
  tax_dm$lineage <- apply(tax_dm[, rank_dm_new], 1, paste, collapse = ";"); #collapse all tax ranks
  dm_otu <- cbind.data.frame("otu_id" = row.names(tax_dm),
                             "lineage" = tax_dm$lineage,
                             otu_dm); #make new otu table
  row.names(dm_otu) <- c();
  dm_otu$lineage <- gsub(";+$", "", dm_otu$lineage); #remove empty tax names
  
  dm_otu_cmf <- dm_otu[, c("otu_id", "lineage", dm_samples_cmf$sample_id)]; #make otu table ready for heat tree  
  
  PrepareHeatTreePlotDataParse_cmf_res <- PrepareHeatTreePlotDataParse_cmf(dm_otu_cmf, dm_samples_cmf, meta);
  dm_obj_cmf = PrepareHeatTreePlotDataParse_cmf_res;
  mbSetObj$dataSet$selected.grps = dm_samples_cmf$sample_id
  mbSetObj$dataSet$comparison = grp.nms 
  mbSetObj$dataSet$meta = meta
  cls = data.frame(dm_obj_cmf$data$class_data)
  cls = cls[!duplicated(cls$taxon_id),]
  cls = cls[,c(1,4)]
  diff_table = data.frame(dm_obj_cmf$data$diff_table)
  df = merge(diff_table, cls, by="taxon_id")
  df = df[which(df$tax_name != "Root"),]
  df = df[order(df$tax_name),]
  mbSetObj$analSet$diff_table <- df;
  #PerformUnivarTestForCorrBoxPlot(mbSetObj, meta,0.05,"NA",taxalvl,"par");
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj);
    PrepareBoxPlot(taxalvl, mbSetObj$dataSet$meta);
    return(1)
  }else{
    return(.set.mbSetObj(mbSetObj));
  }
}

PrepareBoxPlot <- function(taxrank, variable){
  mbSetObj <- .get.mbSetObj();
  
  selSamples = mbSetObj$dataSet$selected.grps
  grps = mbSetObj$dataSet$comparison
  taxrank_boxplot <- taxrank;
  claslbl_boxplot <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  claslbl_boxplot <- claslbl_boxplot[which(claslbl_boxplot %in% grps)];
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
  dat3t_boxplot <- as.data.frame(t(otu_table(data_boxplot)));
  colnames(dat3t_boxplot) <- nm_boxplot; 
  
  #individual boxplot for features
  box_data <- as.data.frame(dat3t_boxplot[which(rownames(dat3t_boxplot) %in% selSamples),]);
  box_data$class <- claslbl_boxplot;
  mbSetObj$analSet$boxdatacor <- box_data;
  mbSetObj$analSet$var.typecor <- variable
  mbSet<<- mbSetObj
  return(.set.mbSetObj(mbSetObj));
}
