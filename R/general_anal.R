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
      res <- results(diagdds, independentFiltering = FALSE, cooksCutoff = Inf);
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
  
  resTable <- as.data.frame(resTable,check.names=FALSE);
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  fast.write(resTable, file="rnaseq_de.csv");
  
  if(nrow(resTable) > 500){
    resTable<-resTable[1:500, ];
  }
  
  mbSetObj$analSet$rnaseq$resTable <- mbSetObj$analSet$resTable <- as.data.frame(resTable,check.names=FALSE);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  
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

  dat3t <- as.data.frame(t(otu_table(data)),check.names=FALSE);
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
  
  resTable <- as.data.frame(resTable,check.names=FALSE);
  ord.inx <- order(resTable$Pvalues);
  resTable <- resTable[ord.inx, , drop=FALSE];
  fast.write(resTable, file="rnaseq_de.csv");
  
  if(nrow(resTable) > 500){
    resTable<-resTable[1:500, ];
  }
  
  mbSetObj$analSet$rnaseq$resTable <- mbSetObj$analSet$resTable <- as.data.frame(resTable,check.names=FALSE);
  
  #only getting the names of DE features
  diff_ft <<- rownames(resTable)[1:de.Num];
  taxrank <<- taxrank;
  
  #individual boxplot for features
  sigfeat <- rownames(resTable);
  box_data <- as.data.frame(dat3t[ ,sigfeat],check.names=FALSE);
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


Match.Pattern <- function(mbSetObj, dist.name="pearson", pattern=NULL, taxrank, variable,appendname){
  
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
PerformMaaslin <- function(
    mbSetObj,
    analysis.var,
    is.norm = "false",
    comp = NULL,
    ref = NULL,
    block = "NA",
    taxrank = "NA",
    imgNm = "NA",
    thresh = 0.05){

  require(dplyr);
  require(R.utils); 
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if (!exists('adj.vec')) {
    adj.bool <- F;
  } else {
    if (length(adj.vec) > 0) {
      adj.bool <- T;
    } else {
      adj.bool <- F;
    }
  }


  thresh <- as.numeric(thresh);
  adj.vars <- adj.vec;
  
  if (mbSetObj$module.type=="sdp"){
    taxrank<-"OTU";
    if (is.norm == "false"){
      phyloseq_objs <- qs::qread("phyloseq_prenorm_objs.qs")
      input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
    } else {
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
      input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
    }
  } else {
    if (is.norm == "false"){
      phyloseq_objs <- qs::qread("phyloseq_prenorm_objs.qs")
      
      if (taxrank=="OTU"){
        input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
      } else {
        taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
        input.data <- phyloseq_objs$count_tables[[taxrank.inx]] %>% as.data.frame()
      } 
    } else {
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
      if (taxrank=="OTU") {
        input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
      } else {
        taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
        input.data <- phyloseq_objs$count_tables[[taxrank.inx]] %>% as.data.frame()
      } 
    }
  }
  
  meta.nms <- colnames(mbSetObj$dataSet$sample_data)
  input.meta <- mbSetObj$dataSet$sample_data@.Data %>% as.data.frame()
  colnames(input.meta) <- meta.nms
  rownames(input.meta) <- input.meta$sample_id
  if(adj.bool){
    fixed.effects <- c(analysis.var, adj.vars)
    fixed.types <- mbSetObj$dataSet$meta.types[names(mbSetObj$dataSet$meta.types) %in% fixed.effects]
    fixed.types <- fixed.types[match(fixed.effects, names(fixed.types))]
  } else { # to do still
    fixed.effects <- analysis.var
    fixed.types <- mbSetObj$dataSet$meta.types[names(mbSetObj$dataSet$meta.types) == analysis.var]
  }
  
  analysis.type <- fixed.types[fixed.effects == analysis.var]
  disc.effects <- fixed.effects[fixed.types == "disc"]
  mbSetObj$analSet$adj.bool = adj.bool
  mbSetObj$analSet$block = block
  mbSetObj$analSet$analysis.type = analysis.type
  mbSetObj$analSet$disc.effects = disc.effects
  # build refs vector (may need to add for blocking too)
  if(length(disc.effects) > 0){
    if(analysis.type == "disc"){
      refs <- paste0(analysis.var, ",", ref)
      if(length(disc.effects) > 1){
        for(i in c(2:length(disc.effects))){
          ref.temp <- paste0(disc.effects[i], ",", levels(unlist(c(input.meta[,disc.effects[i]])))[1])
          refs <- c(refs, ref.temp)
        }
      }
    } else {
      refs <- c()
      if(length(disc.effects) > 1){
        for(i in c(1:length(disc.effects))){
          ref.temp <- paste0(disc.effects[i], ",", levels(unlist(c(input.meta[,disc.effects[i]])))[1])
          refs <- c(refs, ref.temp)
        }
      }
    }
  }

  # MaAslin does not require samples or orders to exactly match - it takes care of this
  # set normalized/transformation parameters
  if(is.norm == "false"){
    norm.method = "TSS"
    trans.method = "LOG"
  } else {
    norm.method = "NONE"
    trans.method = "NONE"
  }
 
###set adjust paprameter
   if((!adj.bool) & (block == "NA")){
    maaslin.para.adj<<-0
  } else {
    refs <- refs[grep(paste0(analysis.var, ","), refs)];
    
   maaslin.para.noadj<<- list(
       input_data = input.data, 
       input_metadata = input.meta, 
      fixed_effects = c(analysis.var),
      reference = c(refs),
       max_significance = 0.05,
      min_abundance = 0.0,
     min_prevalence = 0.0,
      min_variance = 0.0,
      normalization = norm.method,
     transform = trans.method);
  }

.set.mbSetObj(mbSetObj)
  if(block == "NA"){
    if(length(disc.effects) > 0){ # case: discrete variables, no blocking factor
     maaslin.para<<- list(input_data = input.data, 
        input_metadata = input.meta, 
        fixed_effects = c(fixed.effects),
        reference = c(refs),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method)
      return(1)
      maaslin <- Maaslin2.MicrobiomeAnalyst(
        input_data = input.data, 
        input_metadata = input.meta, 
        fixed_effects = c(fixed.effects),
        reference = c(refs),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method);
    } else { # case: no discrete variables, no blocking factor
    maaslin.para<<- list(input_data = input.data, 
        fixed_effects = c(fixed.effects),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method)
      return(1)
      maaslin <- Maaslin2.MicrobiomeAnalyst(
        input_data = input.data, 
        fixed_effects = c(fixed.effects),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method); 
    }
  } else { # case: discrete variables, blocking factor (blocking factor must be discrete)
  maaslin.para <<-list(check= list(input_data = input.data[1,], 
      input_metadata = input.meta, 
      fixed_effects = c(fixed.effects),
      random_effects = c(block),
      reference = c(refs),
      max_significance = 0.05,
      min_abundance = 0.0,
      min_prevalence = 0.0,
      min_variance = 0.0,
      normalization = norm.method,
      transform = trans.method),
     test=list(input_data = input.data, 
        input_metadata = input.meta, 
        fixed_effects = c(fixed.effects),
        random_effects = c(block),
        reference = c(refs),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method)
     )
    return(2)
    check.rank <- capture.output(Maaslin2.MicrobiomeAnalyst(
      input_data = input.data[1,], 
      input_metadata = input.meta, 
      fixed_effects = c(fixed.effects),
      random_effects = c(block),
      reference = c(refs),
      max_significance = 0.05,
      min_abundance = 0.0,
      min_prevalence = 0.0,
      min_variance = 0.0,
      normalization = norm.method,
      transform = trans.method), type=c("message"));
    
    if((length(grep("rank deficient", check.rank)) + length(grep("singular", check.rank))) > 0){
      # often random effects model matrix are rank deficient - check this way and return 
      # feedback that the experimental design does not support using a blocking factor.
      return(-2)
    } else {

   maaslin <- Maaslin2.MicrobiomeAnalyst(
        input_data = input.data, 
        input_metadata = input.meta, 
        fixed_effects = c(fixed.effects),
        random_effects = c(block),
        reference = c(refs),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method);
    }
  }

  res <- maaslin$results
  inds <- !(res$feature %in% rownames(input.data)); # rownames that are all integers have "X" appended to front
  res$feature[inds] <- substring(res$feature[inds], 2);
  
  # get unadjusted results
  if((!adj.bool) & (block == "NA")){
    res.noadj <- res;
  } else {
    refs <- refs[grep(paste0(analysis.var, ","), refs)];
    
    maaslin.noadj <- Maaslin2.MicrobiomeAnalyst(
      input_data = input.data, 
      input_metadata = input.meta, 
      fixed_effects = c(analysis.var),
      reference = c(refs),
      max_significance = 0.05,
      min_abundance = 0.0,
      min_prevalence = 0.0,
      min_variance = 0.0,
      normalization = norm.method,
      transform = trans.method);
    
    res.noadj <- maaslin.noadj$results;
  }
  
  inds <- !(res.noadj$feature %in% rownames(input.data)) # rownames that are all integers have "X" appended to front
  res.noadj$feature[inds] <- substring(res.noadj$feature[inds], 2)
  
  # filter results to get only ones related to analysis var
  res <- res[res$metadata == analysis.var, ];
  
  # make res pretty
  res$coef <- signif(res$coef, digits = 3);
  res$stderr <- signif(res$stderr, digits = 3);
  res$pval <- signif(res$pval, digits = 3);
  res$qval <- signif(res$qval, digits = 3);
  
  if(analysis.type == "disc"){
    res <- res[res$value == comp, ];
    res.noadj <- res.noadj[res.noadj$value == comp, ];
    rownames(res) <- res$feature;
    rownames(res.noadj) <- res.noadj$feature;
    res <- res[ ,c("coef", "stderr", "pval", "qval")];
    res.noadj <- res.noadj[ ,c("coef", "stderr", "pval", "qval")];
    colnames(res) <- c("Log2FC", "St.Error", "P-value", "FDR");
    colnames(res.noadj) <- c("Log2FC", "St. Error", "P-value", "FDR");
  } else {
    rownames(res) <- res$feature;
    rownames(res.noadj) <- res.noadj$feature;
    res <- res[ ,c("coef", "stderr", "pval", "qval")];
    res.noadj <- res.noadj[ ,c("coef", "stderr", "pval", "qval")];
    colnames(res) <- c("Coefficient", "St.Error", "P-value", "FDR");
    colnames(res.noadj) <- c("Coefficient", "St.Error", "P-value", "FDR");
  }
  
  # write out/save results
  fileName <- "multifac_output.csv";
  fast.write(res, file = fileName);
  
  # put results in mbSetObj, learn pattern of analysis set
  sigfeat <- rownames(res)[res$FDR < thresh];
  sig.count <- length(sigfeat);
  if(sig.count == 0){
    current.msg <<- "No significant features were identified using the given p value cutoff.";
  }else{
    current.msg <<- paste("A total of", sig.count, "significant features were identified!");
  }

  # process data for individual feature boxplot
  taxrank_boxplot <- taxrank;
  claslbl_boxplot <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[analysis.var]]);
  nm_boxplot <- rownames(input.data);
  dat3t_boxplot <- as.data.frame(t(input.data),check.names=FALSE);
  colnames(dat3t_boxplot) <- nm_boxplot; 
  box_data <- dat3t_boxplot;
  box_data$class <- claslbl_boxplot;
  box_data$norm <- is.norm;
  
  
  # for graphical summary
  adj.mat <- res[, c("P-value", "FDR")]
  noadj.mat <- res.noadj[, c("P-value", "FDR")]
  
  colnames(adj.mat) <- c("pval.adj", "fdr.adj")
  colnames(noadj.mat) <- c("pval.no", "fdr.no")
  
  both.mat <- merge(adj.mat, noadj.mat, by = "row.names")
  both.mat$pval.adj <- -log10(both.mat$pval.adj)
  both.mat$fdr.adj <- -log10(both.mat$fdr.adj)
  both.mat$pval.no <- -log10(both.mat$pval.no)
  both.mat$fdr.no <- -log10(both.mat$fdr.no)
  
  rownames(both.mat) = both.mat[,1]
  
  # for plotting adjp vs p
  jsonNm <- gsub(".png", ".json", imgNm);
  jsonObj <- RJSONIO::toJSON(both.mat);
  sink(jsonNm);
  cat(jsonObj);
  sink();
  
  mbSetObj$analSet$cov.mat <- both.mat; 
  mbSetObj$analSet$multiboxdata <- box_data;
  mbSetObj$analSet$sig.count <- sig.count;
  mbSetObj$analSet$resTable <- res;
  mbSetObj$analSet$maas.resnoadj <- res.noadj;
  
  return(.set.mbSetObj(mbSetObj))
}


Maaslin2.MicrobiomeAnalyst <-
    function(
        input_data,
        input_metadata,
        min_abundance = 0.0,
        min_prevalence = 0.1,
        min_variance = 0.0,
        normalization = "TSS",
        transform = "LOG",
        analysis_method = "LM",
        max_significance = 0.25,
        random_effects = NULL,
        fixed_effects = NULL,
        correction = "BH",
        standardize = TRUE,
        cores = 1,
        reference = NULL)
    {
      
      require('data.table')
      require('dplyr') 
        # Allow for lower case variables
        normalization <- toupper(normalization);
        transform <- toupper(transform);
        analysis_method <- toupper(analysis_method);
        correction <- toupper(correction);

        #################################################################
        # Read in the data and metadata, create output folder, init log #
        #################################################################
        # if a character string then this is a file name, else it 
        # is a data frame
        if (is.character(input_data)) {
            data <- data.frame(data.table::fread(input_data, header = TRUE, sep = "\t"), row.names = 1);
            if (nrow(data) == 1) {
                # read again to get row name
                data <- read.table(input_data, header = TRUE, row.names = 1);
            }
        } else {
            data <- input_data;
        }
        if (is.character(input_metadata)) {
            metadata <- data.frame(data.table::fread(input_metadata, header = TRUE, sep = "\t"), row.names = 1);
            if (nrow(metadata) == 1) {
                metadata <- read.table(input_metadata, header = TRUE, row.names = 1);
            }
        } else {
            metadata <- input_metadata;
        }

        ###############################################################
        # Determine orientation of data in input and reorder to match #
        ###############################################################
        
        samples_row_row <- intersect(rownames(data), rownames(metadata))
        if (length(samples_row_row) > 0) {
            # this is the expected formatting so do not modify data frames
        } else {
            samples_column_row <- intersect(colnames(data), rownames(metadata));

            if (length(samples_column_row) == 0) {
                # modify possibly included special chars in sample names in metadata
                rownames(metadata) <- make.names(rownames(metadata));
                samples_column_row <- intersect(colnames(data), rownames(metadata));
            }

            if (length(samples_column_row) > 0) {
                # transpose data frame so samples are rows
                data <- as.data.frame(t(data));
            } else {
                samples_column_column <- intersect(colnames(data), colnames(metadata));
                if (length(samples_column_column) > 0) {
                    data <- as.data.frame(t(data));
                    metadata <- type.convert(as.data.frame(t(metadata)));
                } else {
                    samples_row_column <- intersect(rownames(data), colnames(metadata));

                    if (length(samples_row_column) == 0) {
                        # modify possibly included special chars in sample names in data
                        rownames(data) <- make.names(rownames(data));
                        samples_row_column <- intersect(rownames(data), colnames(metadata));
                    }

                    if (length(samples_row_column) > 0) {
                        metadata <- type.convert(as.data.frame(t(metadata)));
                    } else {
                        stop()
                    }
                }
            }
        }
       
        # replace unexpected characters in feature names
        colnames(data) <- make.names(colnames(data))
 
        # check for samples without metadata
        extra_feature_samples <- setdiff(rownames(data), rownames(metadata))
        # check for metadata samples without features
        extra_metadata_samples <- setdiff(rownames(metadata), rownames(data))
        
        # get a set of the samples with both metadata and features
        intersect_samples <- intersect(rownames(data), rownames(metadata))
        
        # now order both data and metadata with the same sample ordering
        data <- data[intersect_samples, , drop = FALSE]
        metadata <- metadata[intersect_samples, , drop = FALSE]
        
        ###########################################
        # Compute the formula based on user input #
        ###########################################
        
        random_effects_formula <- NULL
        # use all metadata if no fixed effects are provided
        if (is.null(fixed_effects)) {
            fixed_effects <- colnames(metadata)
        } else {
            fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
            # remove any fixed effects not found in metadata names
            to_remove <- setdiff(fixed_effects, colnames(metadata))
            if (length(to_remove) > 0)
            fixed_effects <- setdiff(fixed_effects, to_remove)
            if (length(fixed_effects) == 0) {
                stop()
            }
        }
        
        if (!is.null(random_effects)) {
            random_effects <- unlist(strsplit(random_effects, ",", fixed = TRUE))
            
            # subtract random effects from fixed effects
            fixed_effects <- setdiff(fixed_effects, random_effects)
            
            # remove any random effects not found in metadata
            to_remove <- setdiff(random_effects, colnames(metadata))
            
            if (length(to_remove) > 0)
            random_effects <- setdiff(random_effects, to_remove)
            
            # create formula
            if (length(random_effects) > 0) {
                random_effects_formula_text <- paste("expr ~ (1 | ",
                        paste(random_effects, ")", sep = '', collapse = " + (1 | "), sep = '');
                
                random_effects_formula <- tryCatch(as.formula(random_effects_formula_text),
                                                   error = function(e)
                                                     stop(paste("Invalid formula for random effects: ",
                                                                random_effects_formula_text)));
            }
        }
        
        # reduce metadata to only include fixed/random effects in formula
        effects_names <- union(fixed_effects, random_effects)
        metadata <- metadata[, effects_names, drop = FALSE]
        
        # create the fixed effects formula text
        formula_text <- paste("expr ~ ", paste(fixed_effects, collapse = " + "));
        formula <- tryCatch(as.formula(formula_text),
                error = function(e)
                    stop(paste("Invalid formula.",
                               "Please provide a different formula: ",
                               formula_text)));
        #########################################################
        # Filter data based on min abundance and min prevalence #
        #########################################################

        # use ordered factor for variables with more than two levels
        # find variables with more than two levels
        if (is.null(reference)) {reference <- ","}

        for ( i in colnames(metadata) ) {
            mlevels <- unique(na.omit(metadata[,i]));
            numeric_levels <- grep('^-?[0-9.]+[eE+-]?', mlevels, value = T);
            if ( ( length(mlevels[! (mlevels %in% c("UNK"))]) > 1 ) &&  # modification to allow setting reference when only two classes in metadata
                 ( i %in% fixed_effects ) &&
                 ( length(numeric_levels) == 0)) {
                    split_reference <- unlist(strsplit(reference, "[,;]"));
                if (! i %in% split_reference ) {
                    stop(paste("Please provide the reference for the variable '",
                        i, "' which includes more than 2 levels: ",
                        paste(as.character(mlevels), collapse=", "), ".", sep=""));
                } else {
                    ref <- split_reference[match(i,split_reference)+1];
                    other_levels <- as.character(mlevels)[! as.character(mlevels) == ref];
                    metadata[,i] <- factor(metadata[,i], levels=c(ref, other_levels));
                }
            }
        }       
 
        unfiltered_data <- data;
        unfiltered_metadata <- metadata;
        
        # require at least total samples * min prevalence values 
        # for each feature to be greater than min abundance
        total_samples <- nrow(unfiltered_data);
        min_samples <- total_samples * min_prevalence;
        
        # Filter by abundance using zero as value for NAs
        data_zeros <- unfiltered_data;
        data_zeros[is.na(data_zeros)] <- 0;
        filtered_data <- unfiltered_data[, colSums(data_zeros > min_abundance) > min_samples, drop = FALSE];
        total_filtered_features <- ncol(unfiltered_data) - ncol(filtered_data);
        filtered_feature_names <- setdiff(names(unfiltered_data), names(filtered_data));
        
        #################################
        # Filter data based on variance #
        #################################
        
        sds <- apply(filtered_data, 2, sd)
        variance_filtered_data <- filtered_data[, which(sds > min_variance), drop = FALSE]
        variance_filtered_features <- ncol(filtered_data) - ncol(variance_filtered_data)
        variance_filtered_feature_names <- setdiff(names(filtered_data), names(variance_filtered_data))
        filtered_data <- variance_filtered_data
       
        ######################
        # Normalize features #
        ######################
        
        filtered_data_norm <- normalizeFeatures(filtered_data, normalization = normalization)
        
        ################################
        # Standardize metadata, if set #
        ################################
        
        if (standardize) {
            metadata <- metadata %>% dplyr::mutate_if(is.numeric, scale)
        }
        
        ############################
        # Transform and run method #
        ############################
       
        # transform features
        filtered_data_norm_transformed <- transformFeatures(filtered_data_norm, transformation = transform)
        
        # apply the method to the data with the correction

    print( sum(memory.profile()) )   
#Rprof(line.profiling = TRUE, memory.profiling = TRUE)

 fit_data <- fit.data(
                filtered_data_norm_transformed,
                metadata,
                analysis_method,
                formula = formula,
                random_effects_formula = random_effects_formula,
                correction = correction,
                cores = cores);
 
    print( sum(memory.profile()) )  

#Rprof(NULL)
#print(summaryRprof(lines = "both", memory = "both"))

    print( sum(memory.profile()) )  
        
        fit_data$results$N <- apply(fit_data$results, 1, FUN = function(x)
                    length(filtered_data_norm[, x[1]]));
        
        fit_data$results$N.not.zero <- apply(fit_data$results, 1, FUN = function(x)
                    length(which(filtered_data_norm[, x[1]] > 0)));
        
        fit_data$input$data <- filtered_data_norm_transformed
        fit_data$input$metadata <- metadata
        fit_data$input$analysis_method <- analysis_method
        fit_data$input$formula <- formula
        fit_data$input$random_effects_formula <- random_effects_formula
        fit_data$input$correction <- correction
        
        return(fit_data)
   
    }

############ use RserveMicro to perform MaasLin2
.prepare.maaslin2<-function(case,input_data,
                            input_metadata,
                            min_abundance = 0.0,
                            min_prevalence = 0.1,
                            min_variance = 0.0,
                            normalization = "TSS",
                            transform = "LOG",
                            analysis_method = "LM",
                            max_significance = 0.25,
                            random_effects = NULL,
                            fixed_effects = NULL,
                            correction = "BH",
                            standardize = TRUE,
                            cores = 1,
                            reference = NULL){
  require('data.table')
  require('dplyr')
 
if(case==1){
  input_data = maaslin.para$input_data
  if(exists("input_metadata",where = maaslin.para)){
    input_metadata = maaslin.para$input_metadata
    fixed_effects = maaslin.para$fixed_effects
  }
  reference = maaslin.para$reference
  max_significance = maaslin.para$max_significance
  min_abundance = maaslin.para$min_abundance
  min_prevalence = maaslin.para$min_prevalence
  min_variance = maaslin.para$min_variance
  normalization = maaslin.para$normalization
  transform = maaslin.para$transform
  
}else if(case==2){
  
  input_data = input_data
  input_metadata = maaslin.para$check$input_metadata
  fixed_effects = maaslin.para$check$fixed_effects
  random_effects = maaslin.para$check$random_effects
  reference = maaslin.para$check$reference
  max_significance = maaslin.para$check$max_significance
  min_abundance = maaslin.para$check$min_abundance
  min_prevalence = maaslin.para$check$min_prevalence
  min_variance = maaslin.para$check$min_variance
  normalization = maaslin.para$check$normalization
  transform = maaslin.para$check$transform
  
}else if(case==3){
  
  input_data = input_data
  input_metadata = maaslin.para$test$input_metadata
  fixed_effects = maaslin.para$test$fixed_effects
  random_effects = maaslin.para$test$random_effects
  reference = maaslin.para$test$reference
  max_significance = maaslin.para$test$max_significance
  min_abundance = maaslin.para$test$min_abundance
  min_prevalence = maaslin.para$test$min_prevalence
  min_variance = maaslin.para$test$min_variance
  normalization = maaslin.para$test$normalization
  transform = maaslin.para$test$transform
  
}else if(case==4){
  
  input_data = maaslin.para.noadj$input_data 
  input_metadata = maaslin.para.noadj$input_metadata
  fixed_effects = maaslin.para.noadj$fixed_effects
  reference = maaslin.para.noadj$reference
  max_significance = maaslin.para.noadj$max_significance
  min_abundance = maaslin.para.noadj$min_abundance
  min_prevalence = maaslin.para.noadj$min_prevalence
  min_variance = maaslin.para.noadj$min_variance
  normalization = maaslin.para.noadj$normalization
  transform = maaslin.para.noadj$transform
}

  # Allow for lower case variables
  normalization <- toupper(normalization);
  transform <- toupper(transform);
  analysis_method <- toupper(analysis_method);
  correction <- toupper(correction);
  
  #################################################################
  # Read in the data and metadata, create output folder, init log #
  #################################################################
  # if a character string then this is a file name, else it 
  # is a data frame
  if (is.character(input_data)) {
    data <- data.frame(data.table::fread(input_data, header = TRUE, sep = "\t"), row.names = 1);
    if (nrow(data) == 1) {
      # read again to get row name
      data <- read.table(input_data, header = TRUE, row.names = 1);
    }
  } else {
    data <- input_data;
  }
  if (is.character(input_metadata)) {
    metadata <- data.frame(data.table::fread(input_metadata, header = TRUE, sep = "\t"), row.names = 1);
    if (nrow(metadata) == 1) {
      metadata <- read.table(input_metadata, header = TRUE, row.names = 1);
    }
  } else {
    metadata <- input_metadata;
  }
  
  ###############################################################
  # Determine orientation of data in input and reorder to match #
  ###############################################################
  
  samples_row_row <- intersect(rownames(data), rownames(metadata))
  if (length(samples_row_row) > 0) {
    # this is the expected formatting so do not modify data frames
  } else {
    samples_column_row <- intersect(colnames(data), rownames(metadata));
    
    if (length(samples_column_row) == 0) {
      # modify possibly included special chars in sample names in metadata
      rownames(metadata) <- make.names(rownames(metadata));
      samples_column_row <- intersect(colnames(data), rownames(metadata));
    }
    
    if (length(samples_column_row) > 0) {
      # transpose data frame so samples are rows
      data <- as.data.frame(t(data));
    } else {
      samples_column_column <- intersect(colnames(data), colnames(metadata));
      if (length(samples_column_column) > 0) {
        data <- as.data.frame(t(data));
        metadata <- type.convert(as.data.frame(t(metadata)));
      } else {
        samples_row_column <- intersect(rownames(data), colnames(metadata));
        
        if (length(samples_row_column) == 0) {
          # modify possibly included special chars in sample names in data
          rownames(data) <- make.names(rownames(data));
          samples_row_column <- intersect(rownames(data), colnames(metadata));
        }
        
        if (length(samples_row_column) > 0) {
          metadata <- type.convert(as.data.frame(t(metadata)));
        } else {
          stop()
        }
      }
    }
  }
  
  # replace unexpected characters in feature names
  colnames(data) <- make.names(colnames(data))
  
  # check for samples without metadata
  extra_feature_samples <- setdiff(rownames(data), rownames(metadata))
  # check for metadata samples without features
  extra_metadata_samples <- setdiff(rownames(metadata), rownames(data))
  
  # get a set of the samples with both metadata and features
  intersect_samples <- intersect(rownames(data), rownames(metadata))
  
  # now order both data and metadata with the same sample ordering
  data <- data[intersect_samples, , drop = FALSE]
  metadata <- metadata[intersect_samples, , drop = FALSE]
  
  ###########################################
  # Compute the formula based on user input #
  ###########################################
  
  random_effects_formula <- NULL
  # use all metadata if no fixed effects are provided
  if (is.null(fixed_effects)) {
    fixed_effects <- colnames(metadata)
  } else {
    fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
    # remove any fixed effects not found in metadata names
    to_remove <- setdiff(fixed_effects, colnames(metadata))
    if (length(to_remove) > 0)
      fixed_effects <- setdiff(fixed_effects, to_remove)
    if (length(fixed_effects) == 0) {
      stop()
    }
  }
  
  if (!is.null(random_effects)) {
    random_effects <- unlist(strsplit(random_effects, ",", fixed = TRUE))
    
    # subtract random effects from fixed effects
    fixed_effects <- setdiff(fixed_effects, random_effects)
    
    # remove any random effects not found in metadata
    to_remove <- setdiff(random_effects, colnames(metadata))
    
    if (length(to_remove) > 0)
      random_effects <- setdiff(random_effects, to_remove)
    
    # create formula
    if (length(random_effects) > 0) {
      random_effects_formula_text <- paste("expr ~ (1 | ",
                                           paste(random_effects, ")", sep = '', collapse = " + (1 | "), sep = '');
      
      random_effects_formula <- tryCatch(as.formula(random_effects_formula_text),
                                         error = function(e)
                                           stop(paste("Invalid formula for random effects: ",
                                                      random_effects_formula_text)));
    }
  }
  
  # reduce metadata to only include fixed/random effects in formula
  effects_names <- union(fixed_effects, random_effects)
  metadata <- metadata[, effects_names, drop = FALSE]
  
  # create the fixed effects formula text
  formula_text <- paste("expr ~ ", paste(fixed_effects, collapse = " + "));
  formula <- tryCatch(as.formula(formula_text),
                      error = function(e)
                        stop(paste("Invalid formula.",
                                   "Please provide a different formula: ",
                                   formula_text)));
  #########################################################
  # Filter data based on min abundance and min prevalence #
  #########################################################
  
  # use ordered factor for variables with more than two levels
  # find variables with more than two levels
  if (is.null(reference)) {reference <- ","}
  
  for ( i in colnames(metadata) ) {
    mlevels <- unique(na.omit(metadata[,i]));
    numeric_levels <- grep('^-?[0-9.]+[eE+-]?', mlevels, value = T);
    if ( ( length(mlevels[! (mlevels %in% c("UNK"))]) > 1 ) &&  # modification to allow setting reference when only two classes in metadata
         ( i %in% fixed_effects ) &&
         ( length(numeric_levels) == 0)) {
      split_reference <- unlist(strsplit(reference, "[,;]"));
      if (! i %in% split_reference ) {
        stop(paste("Please provide the reference for the variable '",
                   i, "' which includes more than 2 levels: ",
                   paste(as.character(mlevels), collapse=", "), ".", sep=""));
      } else {
        ref <- split_reference[match(i,split_reference)+1];
        other_levels <- as.character(mlevels)[! as.character(mlevels) == ref];
        metadata[,i] <- factor(metadata[,i], levels=c(ref, other_levels));
      }
    }
  }       
  
  unfiltered_data <- data;
  unfiltered_metadata <- metadata;
  
  # require at least total samples * min prevalence values 
  # for each feature to be greater than min abundance
  total_samples <- nrow(unfiltered_data);
  min_samples <- total_samples * min_prevalence;
  
  # Filter by abundance using zero as value for NAs
  data_zeros <- unfiltered_data;
  data_zeros[is.na(data_zeros)] <- 0;
  filtered_data <- unfiltered_data[, colSums(data_zeros > min_abundance) > min_samples, drop = FALSE];
  total_filtered_features <- ncol(unfiltered_data) - ncol(filtered_data);
  filtered_feature_names <- setdiff(names(unfiltered_data), names(filtered_data));
  
  #################################
  # Filter data based on variance #
  #################################
  
  sds <- apply(filtered_data, 2, sd)
  variance_filtered_data <- filtered_data[, which(sds > min_variance), drop = FALSE]
  variance_filtered_features <- ncol(filtered_data) - ncol(variance_filtered_data)
  variance_filtered_feature_names <- setdiff(names(filtered_data), names(variance_filtered_data))
  filtered_data <- variance_filtered_data
  
  ######################
  # Normalize features #
  ######################
  
  filtered_data_norm <- normalizeFeatures(filtered_data, normalization = normalization)
  
  ################################
  # Standardize metadata, if set #
  ################################
  
  if (standardize) {
    metadata <- metadata %>% dplyr::mutate_if(is.numeric, scale)
  }
  
  ############################
  # Transform and run method #
  ############################
  
  # transform features
  filtered_data_norm_transformed <- transformFeatures(filtered_data_norm, transformation = transform)
  
  dat.in <- list(features=filtered_data_norm_transformed,
                 metadata=metadata, 
                 model=analysis_method, 
                 formula=formula, 
                 random_effects_formula=random_effects_formula,
                 correction=correction,
                 cores=cores);

  my.fun <- function(){
    features=dat.in$features
    metadata=dat.in$metadata
    model=dat.in$model
    formula=dat.in$formula
    random_effects_formula=dat.in$random_effects_formula
    correction=dat.in$correction
    cores=dat.in$cores
    
    if (is.null(formula))
      formula <-
        as.formula(paste(
          "expr ~ ", 
          paste(colnames(metadata), 
                collapse = "+")))
    
    if (!(is.null(random_effects_formula))) {
      formula <-
        paste(
          '. ~', 
          paste(all.vars(formula)[-1], collapse = ' + '), 
          '.', 
          sep = ' + ')
      formula <- update(random_effects_formula, formula)
    }
    
    #############################################################
    # Determine the function and summary for the model selected #
    #############################################################
    
    ################
    # Linear Model #
    ################
    
    if (model == "LM") {
      if (is.null(random_effects_formula)) {
        model_function <-
          function(formula, data, na.action) {
            return(glm(
              formula,
              data = data,
              family = 'gaussian',
              na.action = na.action
            ))
          }
        summary_function <- function(fit) {
          lm_summary <- summary(fit)$coefficients
          para <- as.data.frame(lm_summary)[-1, -3]
          para$name <- rownames(lm_summary)[-1]
          return(para)
        }
      } else {
        ranef_function <- lme4::ranef
        model_function <-
          function(formula, data, na.action) {
            return(lmerTest::lmer(
              formula, 
              data = data, 
              na.action = na.action))
          }
        summary_function <- function(fit) {
          lm_summary <- coef(summary(fit))
          para <- as.data.frame(lm_summary)[-1, -c(3:4)]
          para$name <- rownames(lm_summary)[-1]
          return(para)
        }
      }
    }
    
    #######################################
    # Init cluster for parallel computing #
    #######################################
    
    cluster <- NULL
    if (cores > 1)
    {
      cluster <- parallel::makeCluster(cores)
    }
    
    ##############################
    # Apply per-feature modeling #
    ##############################

    outputs <-
      pbapply::pblapply(seq_len(ncol(features)), cl = cluster, function(x) {
        # Extract Features One by One
        featuresVector <- features[, x]
        
        dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
        fit <- tryCatch({
          fit1 <-
            model_function(
              formula, 
              data = dat_sub, 
              na.action = na.exclude)
        }, error = function(err) {
          fit1 <-
            try({
              model_function(
                formula, 
                data = dat_sub, 
                na.action = na.exclude)
            })
          return(fit1)
        })
        
        # Gather Output
        output <- list()
        if (all(!inherits(fit, "try-error"))) {
          output$para <- summary_function(fit)
          if (!(is.null(random_effects_formula))) {
            l <- ranef_function(fit)
            d<-as.vector(unlist(l))
            names(d)<-unlist(lapply(l, row.names))
            output$ranef<-d
          }
        }
        else
        {
          output$para <-
            as.data.frame(matrix(NA, 
                                 nrow = ncol(metadata), ncol = 3))
          output$para$name <- colnames(metadata)
          if (!(is.null(random_effects_formula))) output$ranef <- NA
        }
        colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
        output$para$feature <- colnames(features)[x]
        return(output)
      })
    
    # stop the cluster
    if (!is.null(cluster))
      parallel::stopCluster(cluster)
    
    # bind the results for each feature
    paras <-
      do.call(rbind, lapply(outputs, function(x) {
        return(x$para)
      }))
    
    if (!(is.null(random_effects_formula))) {
      ranef <-
        do.call(rbind, lapply(outputs, function(x) {
          return(x$ranef)
        }))
      row.names(ranef) <- colnames(features) 
    }
    
    ################################
    # Apply correction to p-values #
    ################################
    
    paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
    
    #####################################################
    # Determine the metadata names from the model names #
    #####################################################
    
    metadata_names <- colnames(metadata)
    # order the metadata names by decreasing length
    metadata_names_ordered <-
      metadata_names[order(
        nchar(metadata_names), decreasing = TRUE)]
    # find the metadata name based on the match 
    # to the beginning of the string
    extract_metadata_name <- function(name) {
      return(metadata_names_ordered[mapply(
        startsWith, 
        name, 
        metadata_names_ordered)][1])
    }
    paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))
    # compute the value as the model contrast minus metadata
    paras$value <-
      mapply(function(x, y) {
        if (x == y)
          x
        else
          gsub(x, "", y)
      }, paras$metadata, paras$name)
    
    ##############################
    # Sort by decreasing q-value #
    ##############################
    
    paras <- paras[order(paras$qval, decreasing = FALSE), ]
    paras <-
      dplyr::select(
        paras,
        c('feature', 'metadata', 'value'),
        dplyr::everything())
    
    rownames(paras)<-NULL
    
    if (!(is.null(random_effects_formula))) {
      return(list("results" = paras, "ranef" = ranef))
    } else {
      return(list("results" = paras))
    }
  }
  dat.in <- list(features=filtered_data_norm_transformed,
                 metadata=metadata, 
                 model=analysis_method, 
                 formula=formula, 
                 random_effects_formula=random_effects_formula,
                 correction=correction,
                 cores=cores,
                 my.fun=my.fun,
                filtered_data_norm=filtered_data_norm);
  qs::qsave(dat.in, file="dat.in.qs");
  return(1);
}

.save.maaslin.res <- function(mbSetObj,case){
  mbSetObj <- .get.mbSetObj(mbSetObj);
 dat.in <- qs::qread("dat.in.qs"); 
 filtered_data_norm <- dat.in$filtered_data_norm
 fit_data <- dat.in$my.res
  fit_data$results$N <- apply(fit_data$results, 1, FUN = function(x)
    length(filtered_data_norm[, x[1]]));
  
  fit_data$results$N.not.zero <- apply(fit_data$results, 1, FUN = function(x)
    length(which(filtered_data_norm[, x[1]] > 0)));
  
  fit_data$input$data <- dat.in$filtered_data_norm_transformed
  fit_data$input$metadata <- dat.in$metadata
  fit_data$input$analysis_method <- dat.in$model
  fit_data$input$formula <- dat.in$formula
  fit_data$input$random_effects_formula <- dat.in$random_effects_formula
  fit_data$input$correction <- dat.in$correction
  if(case==4){
    mbSetObj$analSet$maaslin.noadj<-fit_data
  }else if(case !=2){
    mbSetObj$analSet$maaslin<-fit_data
  }
  .set.mbSetObj(mbSetObj)
  if(case==1|case==4){
    print("here")
      return(1)
  }else if(case==2){
    check.rank <- capture.output(fit_data, type=c("message"));
    
    if((length(grep("rank deficient", check.rank)) + length(grep("singular", check.rank))) > 0){
      # often random effects model matrix are rank deficient - check this way and return 
      # feedback that the experimental design does not support using a blocking factor.
      return(-2)
    }else{
      return(3)
    }
   
  }

}
PostProcessMaaslin <- function(mbSetObj,analysis.var,comp=NULL, thresh = 0.05,taxrank,is.norm,imgNm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
input.data<-maaslin.para$input_data
  res <- mbSetObj$analSet$maaslin$results
inds <- !(res$feature %in% rownames(input.data)); # rownames that are all integers have "X" appended to front
res$feature[inds] <- substring(res$feature[inds], 2);

# get unadjusted results
if((!mbSetObj$analSet$adj.bool) & (mbSetObj$analSet$block == "NA")){
  res.noadj <- res;
} else {
  
  res.noadj <- mbSetObj$analSet$maaslin.noadj$results;
}

inds <- !(res.noadj$feature %in% rownames(input.data)) # rownames that are all integers have "X" appended to front
res.noadj$feature[inds] <- substring(res.noadj$feature[inds], 2)

# filter results to get only ones related to analysis var
res <- res[res$metadata == analysis.var, ];

# make res pretty
res$coef <- signif(res$coef, digits = 3);
res$stderr <- signif(res$stderr, digits = 3);
res$pval <- signif(res$pval, digits = 3);
res$qval <- signif(res$qval, digits = 3);

if(  mbSetObj$analSet$analysis.type == "disc"){
  res <- res[res$value == comp, ];
  res.noadj <- res.noadj[res.noadj$value == comp, ];
  rownames(res) <- res$feature;
  rownames(res.noadj) <- res.noadj$feature;
  res <- res[ ,c("coef", "stderr", "pval", "qval")];
  res.noadj <- res.noadj[ ,c("coef", "stderr", "pval", "qval")];
  colnames(res) <- c("Log2FC", "St.Error", "P-value", "FDR");
  colnames(res.noadj) <- c("Log2FC", "St. Error", "P-value", "FDR");
} else {
  rownames(res) <- res$feature;
  rownames(res.noadj) <- res.noadj$feature;
  res <- res[ ,c("coef", "stderr", "pval", "qval")];
  res.noadj <- res.noadj[ ,c("coef", "stderr", "pval", "qval")];
  colnames(res) <- c("Coefficient", "St.Error", "P-value", "FDR");
  colnames(res.noadj) <- c("Coefficient", "St.Error", "P-value", "FDR");
}

# write out/save results
fileName <- "multifac_output.csv";
fast.write(res, file = fileName);
thresh <- as.numeric(thresh);
# put results in mbSetObj, learn pattern of analysis set
sigfeat <- rownames(res)[res$FDR < thresh];
sig.count <- length(sigfeat);
if(sig.count == 0){
  current.msg <<- "No significant features were identified using the given p value cutoff.";
}else{
  current.msg <<- paste("A total of", sig.count, "significant features were identified!");
}

# process data for individual feature boxplot
taxrank_boxplot <- taxrank;
claslbl_boxplot <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[analysis.var]]);
nm_boxplot <- rownames(input.data);
dat3t_boxplot <- as.data.frame(t(input.data),check.names=FALSE);
colnames(dat3t_boxplot) <- nm_boxplot; 
box_data <- dat3t_boxplot;
box_data$class <- claslbl_boxplot;
box_data$norm <- is.norm;


# for graphical summary
adj.mat <- res[, c("P-value", "FDR")]
noadj.mat <- res.noadj[, c("P-value", "FDR")]

colnames(adj.mat) <- c("pval.adj", "fdr.adj")
colnames(noadj.mat) <- c("pval.no", "fdr.no")

both.mat <- merge(adj.mat, noadj.mat, by = "row.names")
both.mat$pval.adj <- -log10(both.mat$pval.adj)
both.mat$fdr.adj <- -log10(both.mat$fdr.adj)
both.mat$pval.no <- -log10(both.mat$pval.no)
both.mat$fdr.no <- -log10(both.mat$fdr.no)

rownames(both.mat) = both.mat[,1]

# for plotting adjp vs p
jsonNm <- gsub(".png", ".json", imgNm);
jsonObj <- RJSONIO::toJSON(both.mat);
sink(jsonNm);
cat(jsonObj);
sink();

mbSetObj$analSet$cov.mat <- both.mat; 
mbSetObj$analSet$multiboxdata <- box_data;
mbSetObj$analSet$sig.count <- sig.count;
mbSetObj$analSet$resTable <- res;
mbSetObj$analSet$maas.resnoadj <- res.noadj;

return(.set.mbSetObj(mbSetObj))
}
############## Fit Function #############
fit.data <- function(
    features,
    metadata,
    model,
    formula = NULL,
    random_effects_formula = NULL,
    correction = "BH",
    cores = 1) {
  # Load Required Packages
 # for (lib in 
       # c('dplyr','pbapply','lmerTest','car','parallel','glmmTMB','MASS','cplm','pscl')) {
   # suppressPackageStartupMessages(require(lib, character.only = TRUE))
  # }
  
  # set the formula default to all fixed effects if not provided
  if (is.null(formula))
    formula <-
      as.formula(paste(
        "expr ~ ", 
        paste(colnames(metadata), 
              collapse = "+")))
  
  if (!(is.null(random_effects_formula))) {
    formula <-
      paste(
        '. ~', 
        paste(all.vars(formula)[-1], collapse = ' + '), 
        '.', 
        sep = ' + ')
    formula <- update(random_effects_formula, formula)
  }
  
  #############################################################
  # Determine the function and summary for the model selected #
  #############################################################
  
  ################
  # Linear Model #
  ################
  
  if (model == "LM") {
    if (is.null(random_effects_formula)) {
      model_function <-
        function(formula, data, na.action) {
          return(glm(
            formula,
            data = data,
            family = 'gaussian',
            na.action = na.action
          ))
        }
      summary_function <- function(fit) {
        lm_summary <- summary(fit)$coefficients
        para <- as.data.frame(lm_summary)[-1, -3]
        para$name <- rownames(lm_summary)[-1]
        return(para)
      }
    } else {
      ranef_function <- lme4::ranef
      model_function <-
        function(formula, data, na.action) {
          return(lmerTest::lmer(
            formula, 
            data = data, 
            na.action = na.action))
        }
      summary_function <- function(fit) {
        lm_summary <- coef(summary(fit))
        para <- as.data.frame(lm_summary)[-1, -c(3:4)]
        para$name <- rownames(lm_summary)[-1]
        return(para)
      }
    }
  }
  
  #######################################
  # Init cluster for parallel computing #
  #######################################
  
  cluster <- NULL
  if (cores > 1)
  {
    cluster <- parallel::makeCluster(cores)
  }
  
  ##############################
  # Apply per-feature modeling #
  ##############################
    outputs <-
   lapply(seq_len(ncol(features)), function(x) {
      # Extract Features One by One
      featuresVector <- features[, x]      
      dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata)
      fit <- tryCatch({
        fit1 <-
          model_function(
            formula, 
            data = dat_sub, 
            na.action = na.exclude)
      }, error = function(err) {
        fit1 <-
          try({
            model_function(
              formula, 
              data = dat_sub, 
              na.action = na.exclude)
          })
        return(fit1)
      })
      
      # Gather Output
      output <- list()
      if (all(!inherits(fit, "try-error"))) {
        output$para <- summary_function(fit)
        if (!(is.null(random_effects_formula))) {
          l <- ranef_function(fit)
          d<-as.vector(unlist(l))
          names(d)<-unlist(lapply(l, row.names))
          output$ranef<-d
        }
      }
      else
      {
        output$para <-
          as.data.frame(matrix(NA, 
                               nrow = ncol(metadata), ncol = 3))
        output$para$name <- colnames(metadata)
        if (!(is.null(random_effects_formula))) output$ranef <- NA
      }
      colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
      output$para$feature <- colnames(features)[x]
      return(output)
    })
  
  # stop the cluster
  if (!is.null(cluster))
    parallel::stopCluster(cluster)
  
  # bind the results for each feature
  paras <-
    do.call(rbind, lapply(outputs, function(x) {
      return(x$para)
    }))
  
  if (!(is.null(random_effects_formula))) {
    ranef <-
      do.call(rbind, lapply(outputs, function(x) {
        return(x$ranef)
      }))
    row.names(ranef) <- colnames(features) 
  }
  
  ################################
  # Apply correction to p-values #
  ################################
  
  paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
  
  #####################################################
  # Determine the metadata names from the model names #
  #####################################################
  
  metadata_names <- colnames(metadata)
  # order the metadata names by decreasing length
  metadata_names_ordered <-
    metadata_names[order(
      nchar(metadata_names), decreasing = TRUE)]
  # find the metadata name based on the match 
  # to the beginning of the string
  extract_metadata_name <- function(name) {
    return(metadata_names_ordered[mapply(
      startsWith, 
      name, 
      metadata_names_ordered)][1])
  }
  paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))
  # compute the value as the model contrast minus metadata
  paras$value <-
    mapply(function(x, y) {
      if (x == y)
        x
      else
        gsub(x, "", y)
    }, paras$metadata, paras$name)
  
  ##############################
  # Sort by decreasing q-value #
  ##############################
  
  paras <- paras[order(paras$qval, decreasing = FALSE), ]
  paras <-
    dplyr::select(
      paras,
      c('feature', 'metadata', 'value'),
      dplyr::everything())
  
  rownames(paras)<-NULL

  if (!(is.null(random_effects_formula))) {
    return(list("results" = paras, "ranef" = ranef))
  } else {
    return(list("results" = paras))
  }
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
  GetSigTable(mbSetObj$analSet$univar$resTable, "Univariate analysis");
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

