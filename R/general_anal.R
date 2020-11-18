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
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    }
    
    if(taxrank=="OTU"){
      data1 <- t(phyloseq_objs$count_tables$OTU)
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- t(phyloseq_objs$count_tables[[taxrank.inx]])
    } 
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
  cls <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
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
  cls <- as.factor(mbSetObj$analSet$cls);
  cls.lbl <- levels(cls);
  cls.len <- length(levels(cls));
  
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
  resTable <- data.frame(resTable);
  resTable$FDR <- p.adjust(resTable$Pvalues,method ="fdr");
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  resTable <- signif(resTable, digits = 5);
  resTable <- resTable[complete.cases(resTable), ];
  resTable <- resTable[,c(2,3,1)];
  fast.write(resTable, file="univar_test_output.csv");
  
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
  fast.write(t(box_data), "uni_abund_data.csv")
  
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
  fast.write(resTable, file="metageno_de_output.csv");
  
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
  
  if(class_no==2){
    indx <- which(res.cls==unique(claslbl)[1]);
    resTable$LDAscore[indx]<--abs(resTable$LDAscore[indx]);
  }
  
  fast.write(resTable, file="lefse_de_output.csv");
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
  
  if(is.na(width)){
    
    sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
    meta <- mbSetObj$analSet$meta
    cls.len <- length(levels(as.factor(sample_table[[meta]])));
    
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
    cls.len <- length(levels(as.factor(sample_table[[meta]])));
    
    if(cls.len <= 3){
      rt.mrg <- 5;
    }else if(cls.len <= 5){
      rt.mrg <- 6.5;
    }else if(cls.len <= 7){
      rt.mrg <- 8;
    }else if(cls.len <= 9){
      rt.mrg <- 9.5;
    }else if(cls.len <= 11){
      rt.mrg <- 11;
    }else if(cls.len <= 12){
      rt.mrg <- 12.5;
    }else{
      rt.mrg <- 15;
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
    data1 <- data.impfeat_lefse;
    mns <- by(data1[, names(imp.vec)], 
              sample_table[[meta]],
              # sample_table[["SampleType"]],
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
  if(opts=="DESeq2"){
    .prepare.deseq(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank);
    .perform.computing();
    res = .save.deseq.res();
  }else{
    data <- .prepare_rnaseq(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank)
    res <- .perform_edger(variable, data, p.lvl)
    
  }
  return(1)
}

.prepare.deseq<-function(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank){
  data <- .prepare_rnaseq(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank)
  mbSetObj = .get.mbSetObj(mbSetObj);
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  if(length(claslbl) > 100){
    AddErrMsg("Only EdgeR is supported for sample size over 100.");
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
      res <- results(diagdds, independentFiltering = FALSE, cooksCutoff = Inf);
      # make sure it is basic R, not DESeq2 obj
      resTable <- data.frame(res[,c("log2FoldChange" ,"lfcSE","pvalue","padj")]); 
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
  resTable <- dat.in$my.res;
  sigHits <- which(resTable$padj < p.lvl);
  de.Num <- length(sigHits);
  if(de.Num == 0){
    current.msg <<- "No significant features were identified using the given p value cutoff.";
  }else{
    current.msg <<- paste("A total of", de.Num, "significant features were identified!");
  }
  print(current.msg);
  
  resTable <- signif(data.matrix(resTable), digits=5);
  colnames(resTable) <- c("log2FC","lfcSE","Pvalues","FDR");
  mbSetObj$analSet$anal.type <- "deseq";
  
  resTable <- as.data.frame(resTable);
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  fast.write(resTable, file="rnaseq_de.csv");
  
  if(nrow(resTable) > 500){
    resTable<-resTable[1:500, ];
  }
  
  mbSetObj$analSet$rnaseq$resTable <- mbSetObj$analSet$resTable <- as.data.frame(resTable);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  
  #individual boxplot for features
  sigfeat <- rownames(resTable);
  box_data <- as.data.frame(dat3t[ ,sigfeat]);
  colnames(box_data) <- sigfeat;
  box_data$class <- claslbl;
  
  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;      
  tree_data <<- data;
  .set.mbSetObj(mbSetObj)
  return(1)
}


.prepare_rnaseq<-function(mbSetObj, opts, p.lvl, variable, shotgunid, taxrank){
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

  dat3t <- as.data.frame(t(otu_table(data)));
  colnames(dat3t) <- nm;
  
  mbSetObj$analSet$rnaseq$data.rnaseq <- dat3t
  mbSetObj$analSet$var.type <- variable;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$rnaseq.taxalvl <- taxrank;
  mbSetObj$analSet$rnaseq.meth <- opts;
  mbSetObj$analSet$rnaseq.plvl <- p.lvl
  .set.mbSetObj(mbSetObj)
  return(data)
}

.perform_edger<-function(variable,data, p.lvl){
  mbSetObj <- .get.mbSetObj(mbSetObj)
  dat3t<- mbSetObj$analSet$rnaseq$data.rnaseq
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  
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
  
  resTable <- as.data.frame(resTable);
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  fast.write(resTable, file="rnaseq_de.csv");
  
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
  
  tree_data <<- data;
  .set.mbSetObj(mbSetObj)
  return(1);
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
  qs::qsave(data1, "match_data.qs")
  #making boxplot data
  
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  variable <- colnames(sample_table)[1];
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  boxdata <- as.data.frame(data1);
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
  mbSetObj$analSet$resTable <- as.data.frame(sig.mat);
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
PlotCorr <- function(mbSetObj, imgName, format="png", dpi=72, width=NA){
  
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
    h <- w <- 8 
  }else if(width == 0){
    h <- w <- 8 
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
  
  feat.nms <- substr(rownames(cor.res), 1, 18);
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

Match.Pattern <- function(mbSetObj, dist.name="pearson", pattern=NULL, taxrank, variable){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(!.on.public.web){
    clslbl <- sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]];
  }
  
  clslbl <- as.factor(clslbl);
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
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
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
  qs::qsave(data, "match_data.qs")
  
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  boxdata <- as.data.frame(data);
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
  mbSetObj$analSet$corph.meth <- dist.name;
  mbSetObj$analSet$feat.corr <- FALSE;
  mbSetObj$analSet$cor.mat<-sig.mat;
  mbSetObj$analSet$pattern.var <- variable
  
  return(.set.mbSetObj(mbSetObj));
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
  
  # only works for unix
  my.os <- tolower(GetOS());
  
  if(my.os == "mac" | my.os == "windows"){
    AddErrMsg("The built-in fastspar only works for unix system. Please visit https://github.com/scwatts/fastspar for more details.");
    return(0);
  }
  
  if(.on.public.web){
    spar_home <- "../../lib/fastspar/";
  }else if (file.exists("/home/jasmine/Downloads/fastspar/fastspar")){ #jas local
    spar_home <- "/home/jasmine/Downloads/fastspar/"
  }else{
    spar_home <- "https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/fast_spar/"
  }
  
  path_fastspar <- paste(spar_home, "fastspar", sep="");
  path_fastspar_bs <- paste(spar_home, "fastspar_bootstrap", sep="");
  path_fastspar_pvals <- paste(spar_home, "fastspar_pvalues", sep="");
  
  # make sure they are executable
  my.cmd <- paste("chmod a+x", path_fastspar, path_fastspar_bs, path_fastspar_pvals);
  system(my.cmd);
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(.on.public.web){
    load_igraph();
  }
  
  if(opt == "pattern"){
    data1 <- qs::qread("pattern_data.qs")
  }else{
    # first get OTU table to BIOM TSV format
    if(mbSetObj$module.type=="sdp"){
      data1 <- as.matrix(otu_table(mbSetObj$dataSet$norm.phyobj));
      taxrank <- "OTU";
    }else if(mbSetObj$module.type=="mdp"){
      if(!exists("phyloseq_objs")){
        phyloseq_objs <- qs::qread("phyloseq_objs.qs")
      }
      
      if(taxrank=="OTU"){
        data1 <- phyloseq_objs$count_tables$OTU
      }else{
        taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
        data1 <- phyloseq_objs$count_tables[[taxrank.inx]]
      }
    }
  }
  
  if(opt=="corr"){
    qs::qsave(data1, "sparcc_full_data.qs")  
  }
  
  if(!is.null(mbSetObj$dataSet$selected.grps) & opt == "corr"){
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
  
  if(opt %in% c("corr", "feat")){
    lowest.count <- min(apply(data1, 1, FUN = function(x) {min(x[x > 0])}))
  }else{ #only pattern searching
    lowest.count <- min(apply(data1[-nrow(data1),], 1, FUN = function(x) {min(x[x > 0])}))
  }
  
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
  
  if(opt=="corr"){
    qs::qsave(zero.output, "sparcc_data.qs")
  }
  
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
  
  invisible(system(paste(path_fastspar, "--otu_table", otu_table, "--correlation", corr_output, "--covariance", cov_output, ";",
                         "mkdir", bootstrap_counts, ";",
                         path_fastspar_bs,  "--otu_table", otu_table, "--number", permNum, "--prefix", counts_prefix, ";",
                         "find", bootstrap_counts, ">>", test, ";",
                         "cat", test, "| while read line; do echo $line; j=$(basename ${line}); jnew=$(echo cor_${j}); k=$(echo cov_${j}); jnew=$(echo ${jnew} | sed 's/\\.tsv//'); k=$(echo ${k} | sed 's/\\.tsv//');", path_fastspar, "--otu_table ${line} --correlation",
                         corr_bs_output, "--covariance", cov_bs_output, "-i 5 ; done;",
                         path_fastspar_pvals, "--otu_table", otu_table, "--correlation", corr_output, "--prefix", corr_prefix, "--permutations", permNum, "--outfile", pval_output, ";",
                         "rm cor_boot_data*.tsv; rm cov_boot_data*.tsv;", "rm -rf", bootstrap_counts, "rm -rf text.txt"), intern=TRUE, ignore.stdout=TRUE))
  
  sparcc_results <- read.table(file = "sparcc_median_correlation.tsv", sep="\t", stringsAsFactors = FALSE)
  names <- sparcc_results[,1]
  sparcc_results_new <- as.matrix(sparcc_results[,-1])
  colnames(sparcc_results_new) <- rownames(sparcc_results_new) <- names
  
  sparcc_pvals <- read.table(file = "sparcc_pvals.tsv", sep="\t", stringsAsFactors = FALSE)
  sparcc_pvals_new <- as.matrix(sparcc_pvals[,-1])
  colnames(sparcc_pvals_new) <- rownames(sparcc_pvals_new) <- names
  
  if(output == "network"){
    sparcc_results_new[sparcc_results_new==0] <- .00001
    sparcc_corr <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(sparcc_results_new, weighted = TRUE))
    sparcc_pvals_new[sparcc_pvals_new==0] <- .00001
    sparcc_pvals <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(sparcc_pvals_new, weighted = TRUE))
    
    sparcc_combo <- cbind(sparcc_corr, sparcc_pvals[,3])
    colnames(sparcc_combo) <- c("Taxon1", "Taxon2", "Correlation", "P.Value")
    qs::qsave(sparcc_combo, "network_correlation.qs")
    
    sparcc_combo <- sparcc_combo[(abs(sparcc_combo[,3]) > corrCutoff & sparcc_combo[,4] < pvalCutoff),]
    fast.write(sparcc_combo, "correlation_table.csv", row.names = FALSE)
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

#'Function to create Network Summary plot
#'@description This function outputs the network
#'summary plot for a specific taxa.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param colOpt Character, "default", "viridis",
#'"cividis", or "plasma".
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotNetworkSummary <- function(mbSetObj, imgName, taxa, format, width=NA, dpi=72,
                               colOpt="default", num_taxa=10){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  group <- mbSetObj$dataSet$meta
  pvalCutoff <- mbSetObj$dataSet$corr.pval.cutoff
  
  mbSetObj$analSet$network_sum <- imgName;
  imgName <- paste(imgName, ".", format, sep="");
  
  cor.method <- mbSetObj$dataSet$cor.method
  
  network_results <- qs::qread("network_correlation.qs")
  
  if(cor.method != "sparcc"){
    network_results <- network_results[,1:4]
    data <- qs::qread("network_cor_data.qs") 
    taxa_inx1 <- which(network_results$Taxon1 %in% taxa)
    taxa_inx2 <- which(network_results$Taxon2 %in% taxa)
    taxa_results1 <- network_results[taxa_inx1,]
    taxa_results2 <- network_results[taxa_inx2,]
    taxa_results2 <- taxa_results2[, c(2,1,3,4)]
    colnames(taxa_results2) <- colnames(taxa_results1)
    taxa_results <- rbind(taxa_results1, taxa_results2) 
  }else{
    data <- qs::qread("sparcc_full_data.qs") 
    data <- t(data)
    taxa_inx <- which(network_results$Taxon1 %in% taxa)
    taxa_results <- network_results[taxa_inx,]
  }
  
  taxa_results <- taxa_results[!taxa_results$Taxon1 == taxa_results$Taxon2,] 
  taxa_results <- taxa_results[order(-abs(taxa_results[[3]])), ];
  
  if(num_taxa < nrow(taxa_results)) {
    taxa_results <- taxa_results[1:num_taxa,]
  }
  
  # reverse order
  taxa_results <- taxa_results[nrow(taxa_results):1,]
  
  # next create dot-size by p-value
  pval <- taxa_results[[4]]
  corr <- taxa_results[[3]]
  sig.names <- names(corr) <- taxa_results[[2]]
  
  sig.inx <- which(pval < pvalCutoff)
  sig.names[sig.inx] <- paste0(sig.names[sig.inx], "*")
  
  if(sum(corr > 0) == 0){ # all negative correlation
    corr <- rev(corr);
  }
  
  feat.num <- length(corr);
  
  if(is.na(width)){
    if(length(corr) < 5 ){
      h <- 6;
      w <- 8;
    } else if (length(corr) < 10){
      h <- length(corr)/1.3;
      w <- 8;
    } else if (length(corr) < 15){
      h <- length(corr)/1.6;
      w <- 8;
    } else if (length(corr) < 20){
      h <- length(corr)/1.8;
      w <- 8;
    } else if (length(corr) < 25){
      h <- length(corr)/2;
      w <- 8;
    } else if (length(corr) < 30){
      h <- length(corr)/2.2;
      w <- 8;
    } else if (length(corr) < 40){
      h <- length(corr)/2.5;
      w <- 8;
    } else {
      h <- length(corr)/4.5;
      w <- 8;
    };
  }else if(width == 0){
    w <- 8.5;
  }else{
    w <- width;
  }
  
  Cairo::Cairo(file = imgName, unit="in", dpi=72, width=w, height=h, type=format, bg="white");
  
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  keep <- mbSetObj$dataSet$selected.grps
  keep.inx <- rownames(sample_table) %in% keep
  sample_table$Keep <- keep.inx
  sample_table <- eval(parse(text = paste("phyloseq:::subset_samples(sample_table,", "Keep == TRUE)", sep="")))
  
  cls.len <- length(mbSetObj$dataSet$comparison)
  
  if(cls.len == 2){
    rt.mrg <- 7;
  }else if(cls.len == 3){
    rt.mrg <- 8;
  }else if(cls.len == 4){
    rt.mrg <- 9;
  }else if(cls.len == 5){
    rt.mrg <- 10;
  }else if(cls.len == 6){
    rt.mrg <- 11;
  }else{
    rt.mrg <- 12;
  }
  
  op <- par(mar=c(5,10,3,rt.mrg)); # set right side margin with the number of class (always 2)
  
  # as data should already be normalized, use mean/median should be the same
  # mns is a list contains means of all vars at each level
  # convert the list into a matrix with each row contains var averages across different lvls
  data1 <- data[which(rownames(data) %in% mbSetObj$dataSet$selected.grps),]
  mns <- by(data1[, names(corr)], 
            sample_table[[group]],
            function(x){ # inner function note, by send a subset of dataframe
              apply(x, 2, mean, trim=0.1)
            });
  mns <- t(matrix(unlist(mns), ncol=feat.num, byrow=TRUE));
  
  vip.nms <- sig.names;
  names(corr) <- NULL;
  
  cols <- case_when(
    corr > 0.75 ~ "#d47273",
    corr > 0.5 & corr < 0.75 ~ "#dd9091",
    corr > 0.25 & corr < 0.5 ~ "#e6aeaf",
    corr > 0 & corr < 0.25 ~ "#f0cccd",
    corr < 0 & corr > -0.25 ~ "#ccdef0",
    corr > -0.25 ~ "#ccdef0",
    corr > -0.5 ~ "#aecbe6",
    corr > -0.75 ~ "#90b7dd",
    TRUE ~ "#72a4d4"
  )
  
  title <- taxa
  
  dotchart(corr, xlab="Correlation", xlim=c(-1, 1), pch="", main=title, cex.lab = 1.15);
  barplot(corr, space=c(0.5, rep(0, length(corr)-1)), xlim=c(-1,1), xaxt="n", col = cols, add=T, horiz=T);
  
  mtext(side=2, at=1:feat.num, vip.nms, las=2, line=1) # set y-axis (species names)
  
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
  
  cls.lbl <- levels(sample_table[[group]]);
  
  for (n in 1:ncol(mns)){
    points(x,y, bty="n", pch=22, bg=bg[,n], cex=3);
    # now add label
    text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5));
    # shift x, note, this is good for current size
    x <- x + shift/1.25;
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
  
  text(x[1], endy+shifty/8, "High");
  text(x[1], starty-shifty/8, "Low");
  
  par(op);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj));
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

