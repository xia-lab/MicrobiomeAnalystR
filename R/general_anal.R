# current.selected.tax
SetCurrentSelectedTaxLevel<-function(taxLvl){
   current.selected.tax <<- taxLvl;
}

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
RF.Anal <- function(mbSetObj, treeNum, tryNum, randomOn, variable, taxrank, datatype){
    
  if(.on.public.web){
    load_randomforest();
  }

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

  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj,tax_table(mbSetObj$dataSet$proc.phyobj));
  }else{
    data <- mbSetObj$dataSet$norm.phyobj;
  }
  
  #using by default names for shotgun data
  if(datatype=="metageno"){
    taxrank <- "OTU";
  }

  if(taxrank=="OTU"){
    data1 <- as.matrix(t(otu_table(data)));
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data,taxrank);
    nm <- as.character(tax_table(data)[,taxrank]);
    y <- which(is.na(nm)==TRUE);
    #converting NA values to unassigned
    nm[y] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    #all NA club together
    data1 <- sapply(by(data1,rownames(data1),colSums),identity);
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
  plot(mbSetObj$analSet$rf, main="Random Forest classification", col=cols);
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
  dotcolor <- ifelse(color.BW, "darkgrey", "blue");
  dotchart(imp.vec, bg=dotcolor, xlab= xlbl, cex=1.3);

  mtext(side=2, at=1:feat.num, vip.nms, las=2, line=1)

  axis.lims <- par("usr"); # x1, x2, y1 ,y2

  # get character width
  shift <- 2*par("cxy")[1];
  lgd.x <- axis.lims[2] + shift;

  x <- rep(lgd.x, feat.num);
  y <- 1:feat.num;
  par(xpd=T);
  
  if(.on.public.web){
    load_rcolorbrewer();
  }

  nc <- ncol(mns);

  # modified for B/W color
  colorpalette <- ifelse(color.BW, "Greys", "RdYlGn");
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
    text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5));
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

  text(x[1], endy+shifty/8, "High");
  text(x[1], starty-shifty/8, "Low");
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

PerformUnivarTest <- function(mbSetObj, variable, p.lvl, datatype, shotgunid, taxrank, statOpt){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  #rather than whole name from taxonomy just last name.
  cls <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  lvl <- length(levels(cls));
  data <- mbSetObj$dataSet$norm.phyobj;

  #using just normalized abundant data
  if(datatype=="16S"){
    # dynamically add taxa table to norm.phyobj
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
  }else{
    data<-data;
  }

  #using by default names for shotgun data
  if(datatype=="metageno"){
    taxrank<-"OTU";
  }

  if(taxrank=="OTU"){
    data <- data;
    nm <- taxa_names(data);
    data1 <- as.data.frame(t(otu_table(data)));
  }else{
    #merging at taxonomy levels
    data <- tax_glom(data, taxrank);
    nm <- as.character(tax_table(data)[,taxrank]);

    #converting NA values to unassigned
    nm[is.na(nm)] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    #all NA club together
    data1 <- t(t(sapply(by(data1,rownames(data1),colSums),identity)));
    nm <- colnames(data1);
  }
  
  data.classical.univar <<- data1 

  colnames(data1)<-nm;
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
  
  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data_boxplot <- merge_phyloseq(data_boxplot, mbSetObj$dataSet$taxa_table);
  }else{
    data_boxplot <- data_boxplot;
  }
  
  #using by default names for shotgun data
  if(datatype=="metageno"){
    taxrank_boxplot<-"OTU";
  }
  
  if(taxrank_boxplot=="OTU"){
    data_boxplot <- data_boxplot;
    nm_boxplot <- taxa_names(data_boxplot);
  }else{
    #merging at taxonomy levels
    data_boxplot <- fast_tax_glom_first(data_boxplot, taxrank_boxplot);
    nm_boxplot <- as.character(tax_table(data_boxplot)[,taxrank_boxplot]);
    #converting NA values to unassigned
    nm_boxplot[is.na(nm_boxplot)] <- "Not_Assigned";
    data1_boxplot <- as.matrix(otu_table(data_boxplot));
    rownames(data1_boxplot) <- nm_boxplot;
    
    #all NA club together
    data1_boxplot <- as.matrix(t(sapply(by(data1_boxplot, rownames(data1_boxplot), colSums), identity)));
    data1_boxplot <- otu_table(data1_boxplot,taxa_are_rows=T);
    data_boxplot <- merge_phyloseq(data1_boxplot, sample_data(data_boxplot));
    nm_boxplot <- taxa_names(data_boxplot);
  }
  
  dat3t_boxplot <- as.data.frame(t(otu_table(data_boxplot)));
  colnames(dat3t_boxplot) <- nm_boxplot;

  #individual boxplot for features
  box_data <- as.data.frame(dat3t_boxplot[ ,sigfeat]);
  colnames(box_data) <- sigfeat;
  box_data$class <- claslbl_boxplot;
  mbSetObj$analSet$boxdata <- box_data;
  write.csv(t(box_data), "uni_abund_data.csv")
  
  mbSetObj$analSet$datatype <- datatype;
  mbSetObj$analSet$anal.type <- "tt";
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$Univar$resTable <- mbSetObj$analSet$resTable <- resTable;
  mbSetObj$analSet$univar.taxalvl <- taxrank;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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

PerformMetagenomeSeqAnal<-function(mbSetObj, variable, p.lvl, datatype, shotgunid, taxrank, model){

  variable <<-  variable;
  p.lvl <<-  p.lvl; 
  datatype <<-  datatype; 
  shotgunid <<-  shotgunid; 
  taxrank <<-  taxrank; 
  model <<- model; 
 
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_metagenomeseq();
  }
  
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
    
  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
  }else{
    data <- data;
  }

  #using by default names for shotgun data
  if(datatype=="metageno"){
    taxrank <- "OTU";
  }

  if(taxrank=="OTU"){
    data <- data;
    nm <- taxa_names(data);
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data, taxrank);
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
      current.msg <<- paste( "More than two group present in sample variable. This model can only be used with two groups.");
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
  mbSetObj$analSet$datatype <- datatype;
  mbSetObj$analSet$anal.type <- "metagseq";
  mbSetObj$analSet$metageno.taxalvl <- taxrank;
  mbSetObj$analSet$id.type <- shotgunid;
    
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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

PerformLefseAnal <- function(mbSetObj, p.lvl, lda.lvl, variable, isfunc, datatype, shotgunid, taxrank){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  datatype <- datatype;
  
  if(.on.public.web){
    load_mass();
  }
  claslbl <<- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
    
  #normalized data
  data <- mbSetObj$dataSet$norm.phyobj;
    
  if(datatype=="16S"){
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, taxa_table);
  }else{
    data <- data;
  }

  #using by default names for shotgun data
    
  if(datatype=="metageno"){
    taxrank <- "OTU";
  }
    
  if(taxrank=="OTU"){
    data <- data;
    tax_orig <<- nm <- taxa_names(data);
    dat3t <- as.data.frame(t(otu_table(data)));
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data, taxrank);
    tax_orig <<- taxa_names(data);
    nm <- as.character(tax_table(data)[,taxrank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    dat3t <- as.data.frame(t(t(sapply(by(data1, rownames(data1), colSums), identity))));
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
  
  wil_data <- as.data.frame(dat3t);
  #if no subclass within classes then no wilcoxon rank sum test
  wil_datadf <- wil_data;

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
  de.Num <- sum(clapvalues<=p.lvl & ldamean$LDAscore>=lda.lvl)
    
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
    
  #visualizing top features based on LDA score
  ldabar <<- ldabar;
 
  #preparing data for indvidual box plot
  sigfeat <<- rownames(resTable);
  taxrank <<- taxrank;
  box_data <- as.data.frame(wil_datadf[, sigfeat]);
  colnames(box_data) <- sigfeat;
  box_data$class <- claslbl;
    
  mbSetObj$analSet$boxdata <- box_data;
  mbSetObj$analSet$sig.count <- de.Num;
  mbSetObj$analSet$datatype <- datatype;
  mbSetObj$analSet$anal.type <- "lefse";
  mbSetObj$analSet$lefse.taxalvl <- taxrank;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$meta <- variable;

  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}

#'Plot LEfSe summary
#'@description This functions graphically summarizes the LEfSe results
#'in a VIP plot.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param ldaFeature Numeric, input the number of top
#'features to include in the plot.
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
PlotLEfSeSummary <- function(mbSetObj, ldaFeature, imgName, format="png", width = NA, dpi=72) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  set.seed(280561493);
  imgName = paste(imgName, ".", format, sep="");
  ldabar <- ldabar;
  ldabar <-  ldabar[order(-ldabar[[2]]), ];
  if(ldaFeature < nrow(ldabar)) {
     ldabar <- ldabar[1:ldaFeature,];
  };

  vip.score <- ldabar[[2]];
  names(vip.score) <- ldabar[[1]];

  if(is.na(width)){
      if(length(vip.score) < 5 ){
            h <- length(vip.score)/1.2;
            w <- 9;
        } else if (length(vip.score) < 10){
            h <- length(vip.score)/1.4;
            w <- 9;
        } else if (length(vip.score) < 15){
            h <- length(vip.score)/1.6;
            w <- 9;
        } else if (length(vip.score) < 20){
            h <- length(vip.score)/1.8;
            w <- 9;
        } else if (length(vip.score) < 25){
            h <- length(vip.score)/2;
            w <- 9;
        } else if (length(vip.score) < 30){
            h <- length(vip.score)/2.2;
            w <- 9;
        } else if (length(vip.score) < 40){
            h <- length(vip.score)/2.5;
            w <- 9;
        } else {
        h <- length(vip.score)/5;
        w <- 9;
        };
    }else if(width == 0){
      w <- 8;
    }else{
      w <- width;
  }
  

  Cairo::Cairo(file = imgName, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  PlotImpVarLEfSe(mbSetObj, vip.score, mbSetObj$analSet$meta);
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

PlotImpVarLEfSe <- function(mbSetObj, imp.vec, meta, color.BW=FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
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
  
  #if(feat.num > length(imp.vec)){
    #feat.num <- length(imp.vec);
  #}
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
  dotcolor <- ifelse(color.BW, "darkgrey", "blue");
  dotchart(imp.vec, bg=dotcolor, xlab= "LDA score", cex=1.3);
  
  mtext(side=2, at=1:feat.num, vip.nms, las=2, line=1)
  
  axis.lims <- par("usr"); # x1, x2, y1 ,y2
  
  # get character width
  shift <- 2*par("cxy")[1];
  lgd.x <- axis.lims[2] + shift;
  
  x <- rep(lgd.x, feat.num);
  y <- 1:feat.num;
  par(xpd=T);
  
  if(.on.public.web){
    load_rcolorbrewer();
  }
  
  nc <- ncol(mns);
  
  # modified for B/W color
  colorpalette <- ifelse(color.BW, "Greys", "RdYlGn");
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
    text(x[1], axis.lims[4], cls.lbl[n], srt=45, adj=c(0.2,0.5));
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
  
  text(x[1], endy+shifty/8, "High");
  text(x[1], starty-shifty/8, "Low");
  
  par(op);
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

PerformRNAseqDE<-function(mbSetObj, opts, p.lvl, variable, datatype, shotgunid, taxrank){
  
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

  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
  }else{
    data <- data;
  }

  #using by default names for shotgun data
  if(datatype=="metageno"){
    taxrank<-"OTU";
  }

  if(taxrank=="OTU"){
    data <- data;
    nm <- taxa_names(data);
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data, taxrank);
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
  
  data.rnaseq <<- dat3t

  if(opts=="DESeq2"){
    # only for small data set (< 100)
    if(length(claslbl) > 100){
      current.msg <<- "Only EdgeR is supported for sample size over 100.";
      return(0);
    }else{
      
      if(.on.public.web){
        load_deseq();
      }
      
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
  mbSetObj$analSet$datatype <- datatype;
  mbSetObj$analSet$id.type <- shotgunid;
  mbSetObj$analSet$rnaseq.taxalvl <- taxrank;
  mbSetObj$analSet$rnaseq.meth <- opts;
  
  tree_data <<- data;

  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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
  
  if(.on.public.web){
    load_ggplot();
    load_ggrepel();
  }
  
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
  
  anal.type <- anal.type
  print(anal.type)
    
  if(anal.type == "shotgun"){
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
    glom <- MicrobiomeAnalystR:::fast_tax_glom_first(mbSetObj$dataSet$proc.phyobj, taxrank);
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
  merged_table <- merged_table[order(merged_table$FDR),]
  
  size <- length(rownames(feats))
  
  if(size <= 20){
    h <- 550
    w <- 500
  }else if(size <= 35){
    h <- 650
    w <- 500
  }else{
    h <- 800
    w <- 500
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

##################################
###########3D PCoA/PCA############
##################################

#'Main function to perform PCoA analysis
#'@description This functions creates a 3D PCoA plot from the microbiome data.
#'This is used by the Beta-Diversity analysis.
#'The 3D interactive visualization is on the web.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param ordMeth Character, input the name
#'of the ordination method. "PCoA" for principal coordinate analysis and "NMDS" for 
#'non-metric multidimensional scaling.
#'@param distName Character, input the name of the distance method.
#'@param datatype Character, input "16S" if the data is marker
#'gene data and "metageno" if it is metagenomic data.
#'@param taxrank Character, input the taxonomic
#'level for beta-diversity analysis.
#'@param colopt Character, color the data points by the experimental factor,
#'the taxon abundance of a selected taxa, or alpha diversity.
#'@param variable Character, input the name of the experimental factor.
#'@param taxa Character, if the data points are colored by taxon abundance, 
#'input the name of the selected taxa.
#'@param alphaopt Character, if the data points are colored by alpha-diversity, 
#'input the preferred alpha-diversity measure.
#'@param jsonNm Character, input the name of the json file to output.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import vegan
#'@import RJSONIO

PCoA3D.Anal <- function(mbSetObj, ordMeth, distName, datatype, taxrank, colopt, variable, taxa, alphaopt, jsonNm){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_vegan();
  }
  
  variable <<- variable;

  if(taxrank=="OTU"){
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
  }else{
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data, taxrank)
  }

  if(colopt=="taxa"){
    if(taxrank=="OTU"){
      data1 <- as.matrix(otu_table(data));
      feat_data <- as.numeric(data1[taxa,]);
    }else{
      nm <- as.character(tax_table(data)[,taxrank]);
      #converting NA values to unassigned
      nm[is.na(nm)] <- "Not_Assigned";
      data1 <- as.matrix(otu_table(data));
      rownames(data1) <- nm;
      #all NA club together
      data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
      feat_data <- data1[taxa,];
    }
    sample_data(data)$taxa <- feat_data;
    indx <- which(colnames(sample_data(data))=="taxa");
    colnames(sample_data(data))[indx] <- taxa;
  }else if(colopt=="alphadiv"){
    data1 <- mbSetObj$dataSet$proc.phyobj;
    box <- plot_richness(data1, measures = alphaopt);
    alphaboxdata <- box$data;
    sam_nm <- sample_names(data);
    alphaboxdata <- alphaboxdata[alphaboxdata$samples %in% sam_nm,];
    alphaval <- alphaboxdata$value;
    sample_data(data)$alphaopt <- alphaval;
    indx <- which(colnames(sample_data(data))=="alphaopt");
    colnames(sample_data(data))[indx]<-alphaopt;
  }else{
    data<-data;
  }

  datacolby <<- data;

  if(distName=="wunifrac"){
    pg_tree <- readRDS("tree.RDS");
    pg_tb <- tax_table(data);
    pg_ot <- otu_table(data);
    pg_sd <- sample_data(data);
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
    data <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree);

    if(!is.rooted(phy_tree(data))){
        pick_new_outgroup <- function(tree.unrooted){
            treeDT <- cbind(cbind(data.table(tree.unrooted$edge),data.table(length = tree.unrooted$edge.length))[1:Ntip(tree.unrooted)],
                    data.table(id = tree.unrooted$tip.label));
                    new.outgroup <- treeDT[which.max(treeDT$length), ]$id
                    return(new.outgroup);
            }
            new.outgroup <- pick_new_outgroup(phy_tree(data));
            phy_tree(data) <- ape::root(phy_tree(data),
                                outgroup = new.outgroup,
                                resolve.root=TRUE)
    }
    GP.ord <-ordinate(data,ordMeth,"unifrac",weighted=TRUE);
  } else if (distName=="wunifrac"){
    pg_tree <- readRDS("tree.RDS");
    pg_tb <- tax_table(data);
    pg_ot <- otu_table(data);
    pg_sd <- sample_data(data);
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
    data <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree);

    if(!is.rooted(phy_tree(data))){
        pick_new_outgroup <- function(tree.unrooted){
            treeDT <- cbind(cbind(data.table(tree.unrooted$edge),data.table(length = tree.unrooted$edge.length))[1:Ntip(tree.unrooted)],
            data.table(id = tree.unrooted$tip.label));
            new.outgroup <- treeDT[which.max(treeDT$length), ]$id
            return(new.outgroup);
        }
        new.outgroup <- pick_new_outgroup(phy_tree(data));
        phy_tree(data) <- ape::root(phy_tree(data),
                                outgroup = new.outgroup,
                                resolve.root=TRUE)
        }
    GP.ord <-ordinate(data,ordMeth,"unifrac",weighted=FALSE);
  }else{
    GP.ord <- ordinate(data,ordMeth,distName);
  }
    
  # obtain variance explained
  sum.pca <- GP.ord;
  imp.pca <- sum.pca$values;
  std.pca <- imp.pca[1,]; # eigen values
  var.pca <- imp.pca[,2]; # variance explained by each PC
  cum.pca <- imp.pca[5,]; # cummulated variance explained
  sum.pca <- append(sum.pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca));

  pca3d <- list();
  
  if(ordMeth=="NMDS"){
    pca3d$score$axis <- paste("NMDS", 1:3 , sep="");
    coord<-sum.pca$points;
    write.csv(signif(coord,5), file="pcoa_score.csv");
    list2 <- rep(as.numeric(0),nrow(coord));
    coord <- cbind(coord, list2);
    coords <- data.frame(t(signif(coord[,1:3], 5)));
  }else{
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(sum.pca$variance[1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(sum.pca$vectors[,1:3], 5)));
    write.csv(signif(sum.pca$vectors,5), file="pcoa_score.csv");
  }
    
  colnames(coords) <- NULL;
  pca3d$score$xyz <- coords;
  pca3d$score$name <- sample_names(mbSetObj$dataSet$norm.phyobj);
  col.type <- "factor";
    
  if(colopt=="taxa"){
    cls <- sample_data(data)[[taxa]];
    col.type <- "gradient"
    cols <- ComputeColorGradient(cls);
  }else if(colopt=="alphadiv") {
    cls <- sample_data(data)[[alphaopt]];
    col.type <- "gradient";
    cols <- ComputeColorGradient(cls);
  }else{
    cls <- factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
    # now set color for each group
    cols <- unique(as.numeric(cls)) + 1;
  }

  pca3d$score$type <- col.type;
  pca3d$score$facA <- cls;
  rgbcols <- col2rgb(cols);
  cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
  pca3d$score$colors <- cols;
    
  if(.on.public.web){
    load_rjsonio();
  }
    
  json.obj <- RJSONIO::toJSON(pca3d);
  sink(jsonNm);
  cat(json.obj);
  sink();

  if(.on.public.web){
    .set.mbSetObj(mbSetObj) #may need to delete?
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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

FeatureCorrelation <- function(mbSetObj, dist.name, taxrank, taxa, variable, 
                               datatype, shotfeat, shotgunid){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
        
    if(taxrank=="OTU"){
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      data1 <- as.matrix(otu_table(data));
      feat_data <- as.numeric(data1[taxa,]);
    }else{
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      #merging at taxonomy levels
      data <- fast_tax_glom_first(data, taxrank);
      nm <- as.character(tax_table(data)[,taxrank]);
      #converting NA values to unassigned
      nm[is.na(nm)] <- "Not_Assigned";
      data1 <- as.matrix(otu_table(data));
      rownames(data1) <- nm;
      #all NA club together
      data1 <- as.matrix(t(sapply(by(data1,rownames(data1),colSums),identity)));
      feat_data <- data1[taxa,];
    }
  }else{
    data <- mbSetObj$dataSet$norm.phyobj;
  }
    
  if(datatype=="metageno"){
    taxrank <- "OTU";
    data1 <- as.matrix(otu_table(data));
    feat_data <- as.numeric(data1[shotfeat,]);
  }
    
  data1 <- t(data1);
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  
  #making boxplot data
  boxdata <- as.data.frame(data1);
  boxdata$class <- clslbl;
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
  
  if(datatype=="metageno"){
    mbSetObj$analSet$pattern <- shotfeat;
  }else{
    mbSetObj$analSet$pattern <- taxa;
    mbSetObj$analSet$taxrank <- taxrank;
  }
  
  sig.nm <<- fileName;
  mbSetObj$analSet$cor.mat <- sig.mat;
  mbSetObj$analSet$corph.taxalvl <- taxrank;
  mbSetObj$analSet$corph.meth <- dist.name;
  mbSetObj$analSet$sig.count <- 0; # note, not a DE analysis here

  if(.on.public.web){
    .set.mbSetObj(mbSetObj) 
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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
  title <- paste("Top",nrow(cor.res), tolower(mbSetObj$analSet$taxrank), "correlated with the", pattern);
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
  dotchart(cor.res[,1], pch="", xlim=c(-1,1), xlab="Correlation coefficients", main=title);
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

Match.Pattern <- function(mbSetObj, dist.name="pearson", pattern=NULL, taxrank, 
                          taxa, variable, datatype, shotfeat, shotgunid){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
    
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

  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
        
    if(taxrank=="OTU"){
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      data1 <- as.matrix(otu_table(data));
      feat_data <- as.numeric(data1[taxa,]);
        
    }else{
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      #merging at taxonomy levels
      data <- fast_tax_glom_first(data, taxrank);
      nm <- as.character(tax_table(data)[,taxrank]);
      #converting NA values to unassigned
      nm[is.na(nm)] <- "Not_Assigned";
      data1 <- as.matrix(otu_table(data));
      rownames(data1)<-nm;
      #all NA club together
      data1<-as.matrix(t(sapply(by(data1, rownames(data1), colSums), identity)));
      feat_data<-data1[taxa,];
    }
  }else{
    data <- mbSetObj$dataSet$norm.phyobj;
  }
    
  if(datatype=="metageno"){
    taxrank <- "OTU";
    data1 <- as.matrix(otu_table(data));
    feat_data <- as.numeric(data1[shotfeat,]);
  }
    
  data1 <- t(data1);
  clslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  boxdata <- as.data.frame(data1);
  boxdata$class <- clslbl;
  mbSetObj$analSet$boxdata <- boxdata;
  cbtempl.results <- apply(data1, 2, template.match, new.template, dist.name);
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
    
  if(datatype=="metageno"){
    mbSetObj$analSet$pattern<-shotfeat;
  }else{
    mbSetObj$analSet$pattern<-taxa;
    mbSetObj$analSet$taxrank<-taxrank;
  }
  
  sig.nm <<- fileName;
  mbSetObj$analSet$cor.mat<-sig.mat;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj) 
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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

PlotCorrHeatMap <- function(mbSetObj, imgName, format="png", width=NA, cor.method,
                            colors_cntrst, viewOpt, taxrank, fix.col, no.clst, top, 
                            topNum, datatype){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  main <- xlab <- ylab <- NULL;
    
  if(datatype=="metageno"){
    data <- mbSetObj$dataSet$norm.phyobj;
    data1 <- as.matrix(otu_table(data));
    taxrank <- "OTU";
  }
    
  if(datatype=="16S"){
    if(taxrank=="OTU"){
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      data1 <- as.matrix(otu_table(data));
    }else{
      taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
      data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
      #merging at taxonomy levels
      data <- fast_tax_glom_first(data,taxrank);
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
  corr.mat<-cor(data, method=cor.method);

  # use total abs(correlation) to select
  if(top){
    cor.sum <- apply(abs(corr.mat), 1, sum);
    cor.rk <- rank(-cor.sum);
    var.sel <- cor.rk <= topNum;
    corr.mat <- corr.mat[var.sel, var.sel];
  }

  # set up parameter for heatmap
  
  if(.on.public.web){
    load_rcolorbrewer();
    load_gplots();
  }

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
  
  if(.on.public.web){
    load_pheatmap();
  }
  
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
  return(.set.mbSetObj(mbSetObj))
}

###################################
########## Util Functions #########
###################################

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
  
  if(.on.public.web){
    load_edgeR();
    load_phyloseq();
  }

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
  if(.on.public.web){
    load_xtable();
  }
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
  if(.on.public.web){
    load_xtable();
  }
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
  if(.on.public.web){
    load_xtable();
  }    
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
      if(.on.public.web){
        load_genefilter();
      }
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
