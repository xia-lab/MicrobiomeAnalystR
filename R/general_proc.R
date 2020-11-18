##################################################
## R script for MicrobiomeAnalyst
## Description: filtering and normalization functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#######################################
# filter zeros and singletons this is not good enough ==> single occurance
# reads occur in only one sample should be treated as artifacts and removed

#'Main function to sanity check on uploaded data
#'@description This function performs a sanity check on the uploaded
#'user data. It checks the grouping of samples, if a phylogenetic tree
#'was uploaded.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param filetype Character, "biom" if the uploaded data
#'was .biom format, "mothur" if mothur format, or "txt" if
#'as .txt or .csv format.
#'@param disableFilter Boolean. Set to TRUE to bypass the hard-coded filter
#'to remove features occuring in <2 samples.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SanityCheckData <- function(mbSetObj, filetype, disableFilter = FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  feat.sums <- apply(mbSetObj$dataSet$data.orig, 1, function(x){sum(x>0, na.rm=T)});

  if(disableFilter){
    data.proc <- mbSetObj$dataSet$data.orig
  }else{
    gd.inx <- feat.sums > 1; # occur in at least 2 samples
    if(length(which(gd.inx=="TRUE"))==0){
      AddErrMsg("Reads occur in only one sample.  All these are considered as artifacts and have been removed from data. No data left after such processing.");
      return(0);
    }
    data.proc <- mbSetObj$dataSet$data.orig[gd.inx, ]; 
  }

  # filtering the constant features here
  # check for columns with all constant (var=0)
  varCol <- apply(data.proc, 1, var, na.rm=T);
  constCol <- varCol == 0 | is.na(varCol);
  
  # making copy of data.proc and proc.phyobj(phyloseq)
  data.proc <- data.proc[!constCol, ];
  
  if(length(data.proc)==0){
    AddErrMsg("All features are found to be constant and have been removed from data. No data left after such processing.");
    return(0);
  }
  
  qs::qsave(data.proc, file="data.proc.orig"); # save an copy
  qs::qsave(data.proc, file="data.prefilt"); # save an copy
  
  # now get stats
  taxa_no <- nrow(mbSetObj$dataSet$data.orig);
  
  # from abundance table
  sample_no <- ncol(mbSetObj$dataSet$data.orig);
  
  if(filetype=="biom"||filetype=="mothur"){
    samplemeta_no <- nrow(mbSetObj$dataSet$sample_data);
  }else{
    samplemeta_no <- sample_no;
  }
  
  if(mbSetObj$module.type!="ppd"){
    mbSetObj$dataSet$sample_data <- mbSetObj$dataSet$sample_data[sapply(mbSetObj$dataSet$sample_data, function(col) length(unique(col))) > 1];
  }
  
  if(ncol(mbSetObj$dataSet$sample_data)==0){
    AddErrMsg("No sample variable have more than one group. Please provide variables with at least two groups.");
    return(0);
  }
  
  #converting all sample variables to factor type
  character_vars <- sapply(mbSetObj$dataSet$sample_data, is.character);
  
  if(any(character_vars)=="TRUE"){
    mbSetObj$dataSet$sample_data[, character_vars] <- sapply(mbSetObj$dataSet$sample_data[, character_vars], as.factor);
  }
  
  num_vars <- sapply(mbSetObj$dataSet$sample_data, is.numeric);
  
  if(any(num_vars)=="TRUE"){
    mbSetObj$dataSet$sample_data[, num_vars] <- sapply(mbSetObj$dataSet$sample_data[, num_vars],as.factor);
  }
  
  if(file.exists("tree.qs")){
    tree_exist <- 1
  } else {
    tree_exist <- 0
  }
  
  if(identical(sort(row.names(mbSetObj$dataSet$sample_data)),
               sort(colnames(data.proc)))){
    samname_same <- 1;
  } else {
    samname_same <- 0;
  }
  
  sample_no_in_outfile <- ncol(data.proc);
  samname_same_number <- sum(row.names(mbSetObj$dataSet$sample_data) %in% colnames(data.proc));
  
  # now store data.orig to RDS
  qs::qsave(mbSetObj$dataSet$data.orig, file="data.orig");
  qs::qsave(data.proc, file="data.proc");
  qs::qsave(mbSetObj$dataSet$sample_data, file = "data.sample_data")
  mbSetObj$dataSet$data.orig <- NULL;
  mbSetObj$dataSet$tree <- tree_exist
  
  vari_no <- ncol(mbSetObj$dataSet$sample_data);
  disc_no <- sum(mbSetObj$dataSet$meta_info$disc.inx);
  cont_no <- sum(mbSetObj$dataSet$meta_info$cont.inx);
  smpl.sums <- apply(data.proc, 2, sum);
  tot_size <- sum(smpl.sums);
  smin <- min(smpl.sums)
  smean <- mean(smpl.sums)
  smax <- max(smpl.sums);
  gd_feat <- nrow(data.proc);
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(c(1,taxa_no,sample_no,vari_no, smin,smean,smax,gd_feat,samplemeta_no,tot_size, tree_exist, samname_same, samname_same_number, sample_no_in_outfile, disc_no, cont_no));
  }else{
    print("Sanity check passed!")
    return(.set.mbSetObj(mbSetObj))
  }
}

#'Function to filter uploaded data
#'@description This function filters data based on low counts in high percentage samples.
#'Note, first is abundance, followed by variance.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param filt.opt Character, input the low count filter option. "prevalence" to
#'filter based on prevalence in samples, "mean" to filter by the mean abundance value, and
#'"median" to filter by median abundance value.
#'@param count Numeric, input the minimum count. Set to 0 to disable the low count filter.
#'@param smpl.perc Numeric, input the percentage of samples for which to filter low counts.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

ApplyAbundanceFilter <- function(mbSetObj, filt.opt, count, smpl.perc){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data <- qs::qread("data.prefilt");
  msg <- NULL;
  
  #this data is used for sample categorial comparision further
  rmn_feat <- nrow(data);
  
  mbSetObj$dataSet$ab.filtered <- FALSE
  
  if(count==0){# no low-count filtering
    kept.inx <- rep(TRUE, rmn_feat);
  }else{
    mbSetObj$dataSet$ab.filtered <- TRUE
    if(filt.opt == "prevalence"){
      rmn_feat <- nrow(data);
      minLen <- smpl.perc*ncol(data);
      kept.inx <- apply(data, MARGIN = 1,function(x) {sum(x >= count) >= minLen});
    }else if (filt.opt == "mean"){
      filter.val <- apply(data, 1, mean, na.rm=T);
      kept.inx <- filter.val >= count;
    }else if (filt.opt == "median"){
      filter.val <- apply(data, 1, median, na.rm=T);
      kept.inx <- filter.val >= count;
    }
  }
  
  mbSetObj$dataSet$filt.data <- data[kept.inx, ];
  qs::qsave(mbSetObj$dataSet$filt.data, file="filt.data.orig"); # save an copy
  current.msg <<- paste("A total of ", sum(!kept.inx), " low abundance features were removed based on ", filt.opt, ".", sep="");
  
  return(.set.mbSetObj(mbSetObj));
}

#'Function to filter uploaded data
#'@description This function filters data based on low abundace or variance.
#'Note, this is applied after abundance filter.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param filtopt Character, input the low variance filter option. "iqr" for 
#'inter-quantile range, "sd" for standard deviation, and "cov" for coefficient of variation.
#'@param filtPerct Numeric, input the percentage cutoff for low variance. 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

ApplyVarianceFilter <- function(mbSetObj, filtopt, filtPerct){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data <- mbSetObj$dataSet$filt.data;
  msg <- NULL;
  
  rmn_feat <- nrow(data);
  
  filter.val <- nm <- NULL;
  
  mbSetObj$dataSet$var.filtered <- FALSE
  
  if(filtPerct==0){# no low-count filtering
    remain <- rep(TRUE, rmn_feat);
  }else{
    mbSetObj$dataSet$var.filtered <- TRUE
    if (filtopt == "iqr"){
      filter.val <- apply(data, 1, IQR, na.rm=T);
      nm <- "IQR";
    }else if (filtopt == "sd"){
      filter.val <- apply(data, 1, sd, na.rm=T);
      nm <- "standard deviation";
    }else if (filtopt == "cov"){
      sds <- apply(data, 1, sd, na.rm=T);
      mns <- apply(data, 1, mean, na.rm=T);
      filter.val <- abs(sds/mns);
      nm <- "Coeffecient of variation";
    }
    # get the rank
    rk <- rank(-filter.val, ties.method='random');
    var.num <- nrow(data);
    remain <- rk < var.num*(1-filtPerct);
  }
  
  data <- data[remain,];
  mbSetObj$dataSet$filt.data <- data;
  qs::qsave(mbSetObj$dataSet$filt.data, file="filt.data.orig"); # save an copy
  rm.msg1 <- paste("A total of ", sum(!remain), " low variance features were removed based on ", filtopt, ".", sep="");
  rm.msg2 <- paste("The number of features remains after the data filtering step:", nrow(data));
  current.msg <<- paste(c(current.msg, rm.msg1, rm.msg2), collapse=" ");
  mbSetObj$dataSet$filt.msg <- current.msg;
  
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to update samples
#'@description This function prunes samples to be included 
#'for downstream analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
UpdateSampleItems <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(!exists("smpl.nm.vec")){
    AddErrMsg("Cannot find the current sample names!");
    return(0);
  }
  
  # read from saved original copy
  data.proc.orig  <- qs::qread("data.proc.orig");
  proc.phyobj.orig <- qs::qread("proc.phyobj.orig");
  hit.inx <- colnames(data.proc.orig) %in% smpl.nm.vec;
  
  #preserving the rownames
  # read from saved prefilt copy
  prefilt.data <- qs::qread("data.prefilt");
  tax_nm <- rownames(prefilt.data);
  data.proc <- qs::qread("data.proc");
  modi_data <- data.proc[,hit.inx];
  prefilt.data <- modi_data;
  
  #doesn't need to change the names object
  rownames(prefilt.data) <- tax_nm;
  #getting unmatched sample names
  unmhit.indx <- which(hit.inx==FALSE);
  allnm <- colnames(data.proc.orig);
  mbSetObj$dataSet$remsam <- allnm[unmhit.indx];
  
  mbSetObj$dataSet$proc.phyobj <- prune_samples(colnames(prefilt.data),proc.phyobj.orig);
  qs::qsave(prefilt.data, file="data.prefilt")
  current.msg <<- "Successfully updated the sample items!";
  
  # need to update metadata info after removing samples
  my.meta <- sample_data(mbSetObj$dataSet$proc.phyobj)
  disc.inx <- GetDiscreteInx(my.meta);
  
  if(sum(disc.inx) == 0){
    mbSetObj$poor.replicate <- TRUE;
  }else{
    mbSetObj$poor.replicate <- FALSE;
  }
  
  mbSetObj$dataSet$sample_data <- my.meta
  
  return(.set.mbSetObj(mbSetObj));
  
}


#'Function to perform normalization
#'@description This function performs normalization on the uploaded
#'data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param rare.opt Character, "rarewi" to rarefy the data and 
#'"none" to not.
#'@param scale.opt Character, input the name of the data scaling option.
#'"colsum" for total sum scaling, "CSS" for cumulative sum scaling,
#' "upperquartile" for upper-quartile normalization, and "none" to none.
#'@param transform.opt Character, input the name of the data transformation
#'to be applied. "rle" for relative log expression, "TMM" for trimmed mean of 
#'M-values, "clr" for centered log ratio, and "none" for none.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import edgeR
#'@import metagenomeSeq
PerformNormalization <- function(mbSetObj, rare.opt, scale.opt, transform.opt){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  data <- qs::qread("filt.data.orig");
  tax_nm <- rownames(data);
  msg <- NULL;
  
  if(rare.opt != "none"){
    data <- PerformRarefaction(mbSetObj, data, rare.opt);
    tax_nm <- rownames(data);
    msg <- c(msg, paste("Performed data rarefaction."));
    
    ###### note, rarefying (should?) affect filt.data in addition to norm
    mbSetObj$dataSet$filt.data <- data;
  }else{
    msg <- c(msg, paste("No data rarefaction was performed."));
  }
  
  qs::qsave(mbSetObj$dataSet$proc.phyobj, file="orig.phyobj"); # save original phylo.obj
  # now proc.phyobj is now filtered data
  mbSetObj$dataSet$proc.phyobj <- merge_phyloseq(otu_table(data,taxa_are_rows =TRUE), mbSetObj$dataSet$sample_data, mbSetObj$dataSet$taxa_table);
  
  if(scale.opt != "none"){
    if(scale.opt=="colsum"){
      data <- sweep(data, 2, colSums(data), FUN="/")
      data <- data*10000000;
      msg <- c(msg, paste("Performed total sum normalization."));
    }else if(scale.opt=="upperquartile"){
      load_edgeR();
      otuUQ <- edgeRnorm(data,method="upperquartile");
      data <- as.matrix(otuUQ$counts);
      msg <- c(msg, paste("Performed upper quartile normalization"));
    }else if(scale.opt=="CSS"){
      load_metagenomeseq();
      #biom and mothur data also has to be in class(matrix only not in phyloseq:otu_table)
      data1 <- as(data,"matrix");
      dataMR <- newMRexperiment(data1);
      data <- cumNorm(dataMR,p=cumNormStat(dataMR));
      data <- MRcounts(data,norm = T);
      msg <- c(msg, paste("Performed cumulative sum scaling normalization"));
    }else{
      print(paste("Unknown scaling parameter:", scale.opt));
    }
  }else{
    msg <- c(msg, paste("No data scaling was performed."));
  }
  
  if(transform.opt != "none"){
    if(transform.opt=="rle"){
      load_edgeR();            
      otuRLE <- edgeRnorm(data,method="RLE");
      data <- as.matrix(otuRLE$counts);
      msg <- c(msg, paste("Performed RLE Normalization"));
    }else if(transform.opt=="TMM"){
      load_edgeR();   
      otuTMM <- edgeRnorm(data,method="TMM");
      data <- as.matrix(otuTMM$counts);
      msg <- c(msg, paste("Performed TMM Normalization"));
    }else if(transform.opt=="clr"){
      data <- apply(data, 2, clr_transform);
      msg <- "Performed centered-log-ratio normalization.";
    }else{
      print(paste("Unknown scaling parameter:", transform.opt));
    }
  }else{
    msg <- c(msg, paste("No data transformation was performed."));
  }
  #this step has to be done after all the normalization step
  otu.tab <- otu_table(data,taxa_are_rows =TRUE);
  taxa_names(otu.tab) <- tax_nm;
  
  # create phyloseq obj
  mbSetObj$dataSet$sample_data$sample_id <- rownames(mbSetObj$dataSet$sample_data);
  sample_table <- sample_data(mbSetObj$dataSet$sample_data, errorIfNULL=TRUE);
  phy.obj <- merge_phyloseq(otu.tab, sample_table);
  
  #using this object for plotting
  mbSetObj$dataSet$norm.phyobj <- phy.obj;
  current.msg <<- paste(msg, collapse=" ");
  mbSetObj$dataSet$norm.msg <- current.msg;
  
  #make hierarchies
  if(mbSetObj$module.type=="sdp"){
    ranks <- "OTU"
  }else{
    ranks <- c(GetMetaTaxaInfo(mbSetObj), "OTU")
    ranks <- unique(ranks)
  } 
  
  # start with lowest
  data.list <- list()
  data.list$merged_obj <- vector(length = length(ranks), "list")
  data.list$count_tables <- vector(length = length(ranks), "list")
  names(data.list$count_tables) <- names(data.list$merged_obj) <- ranks
  
  for(i in 1:length(ranks)){
    phyloseq.obj <- UtilMakePhyloseqObjs(mbSetObj, ranks[i])
    data.list$merged_obj[[i]] <- phyloseq.obj
    count.table <- UtilMakeCountTables(phyloseq.obj, ranks[i])
    data.list$count_tables[[i]] <- count.table
  }
  qs::qsave(data.list, "phyloseq_objs.qs")
  
  return(.set.mbSetObj(mbSetObj));
}

#'Utility function to perform rarefraction (used by PerformNormalization)
#'@description This function performs rarefraction on the uploaded
#'data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param data Input the data.
#'@param rare.opt Input the option for rarefying the microbiome data.
#'"rarewi" to rarefy with replacement to the minimum library depth and 
#'"rarewo" to rarefy without replacemet to the minimum library depth.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PerformRarefaction <- function(mbSetObj, data, rare.opt){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data <- data.matrix(data);
  tax_nm<-rownames(data);
  
  # data must be count data (not contain fractions)
  data <- round(data);
  
  # create phyloseq obj
  otu.tab<-otu_table(data,taxa_are_rows =TRUE);
  taxa_names(otu.tab)<-tax_nm;
  
  mbSetObj$dataSet$sample_data$sample_id<-rownames(mbSetObj$dataSet$sample_data);
  sample_table<-sample_data(mbSetObj$dataSet$sample_data, errorIfNULL=TRUE);
  phy.obj<-merge_phyloseq(otu.tab,sample_table);
  msg<-NULL;
  
  # first doing rarefaction, this is on integer or count data
  if(rare.opt=="rarewi"){
    phy.obj <- rarefy_even_depth(phy.obj, replace=TRUE,rngseed = T)
    msg <- c(msg, paste("Rarefy with replacement to minimum library depth."));
  }else if(rare.opt=="rarewo"){ # this is not shown on web due to computational issue
    phy.obj <- rarefy_even_depth(phy.obj, replace=FALSE,rngseed = T);
    msg <- c(msg, paste("Rarefaction without replacement to minimum library depth."));
  }
  
  cleanMem();
  otu.tab <- otu_table(phy.obj);
  return(otu.tab);
}

#'Function to plot rarefraction curves
#'@description This function plots rarefraction curves from the uploaded
#'data for both sample-wise or group-wise
#'@param mbSetObj Input the name of the mbSetObj.
#'@param graphName Input the name of the plot.
#'@param variable Input the experimental factor.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import vegan

PlotRareCurve <- function(mbSetObj, graphName, variable){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  set.seed(13789);
  load_vegan();
  
  data <- data.matrix(mbSetObj$dataSet$filt.data);
  rarefaction_curve_data<-as.matrix(otu_table(data));
  
  Cairo::Cairo(file=graphName, width=900, height=480, type="png", bg="white");
  
  #sample-wise
  subsmpl=50;
  
  if (ncol(rarefaction_curve_data)>subsmpl) {
    ss  = sample(ncol(rarefaction_curve_data), subsmpl);
    rarefaction_curve_data = rarefaction_curve_data[,ss,drop=FALSE];
    #retaining the column names
    data<-prune_samples(colnames(rarefaction_curve_data),data);
  } else {
    rarefaction_curve_data = rarefaction_curve_data;
  }
  
  sam_data<-sample_data(data);
  raremax <- min(rowSums(t(rarefaction_curve_data)));
  
  #getting colors
  grp.num <- length(levels(as.factor(sam_data[[variable]])));
  
  if(grp.num<9){
    dist.cols <- 1:grp.num + 1;
    lvs <- levels(as.factor(sam_data[[variable]]));
    colors <- vector(mode="character", length=length(as.factor(sam_data[[variable]])));
    
    for(i in 1:length(lvs)){
      colors[as.factor(sam_data[[variable]]) == lvs[i]] <- dist.cols[i];
    }
  }else{
    colors<-"blue";
  }
  
  rarecurve(t(rarefaction_curve_data), step = 20, sample = raremax, col = colors, cex = 0.6, xlab = "Sequencing depth", ylab = "Observed Species");
  legend("bottomright",lvs,lty = rep(1,grp.num),col = dist.cols);
  
  plot.data <- rarefaction_curve_data;
  sam_data <- sample_data(data);
  cls.lbls <- as.factor(sam_data[[variable]]);
  grp.lvls <- levels(cls.lbls);
  grp.num <- length(grp.lvls);
  grp.data <- list();
  
  for(lvl in grp.lvls){
    grp.data[[lvl]] <- rowSums(plot.data[, cls.lbls == lvl])
  }
  
  raremax <- min(unlist(lapply(grp.data, sum)));
  grp.data <- data.frame(grp.data);
  dist.cols <- 1:grp.num + 1;
  # now plot
  rarecurve(t(grp.data), step = 20, sample = raremax, col = dist.cols, lwd = 2, cex=1.5, xlab = "Sequencing depth", ylab = "Observed Species");
  legend("bottomright",grp.lvls,lty = rep(1,grp.num),col = dist.cols);
  
  dev.off();
  return(.set.mbSetObj(mbSetObj))
}

#'Function to plot library size 
#'@description This function creates a plot summarizing the library
#'size of microbiome dataset.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param imgName Character, input the name
#'of the plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotLibSizeView <- function(mbSetObj, imgName, format="png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data.proc <- qs::qread("data.proc");
  data_bef <- data.matrix(data.proc);
  smpl.sums <- colSums(data_bef);
  names(smpl.sums) <- colnames(data_bef);
  smpl.sums <- sort(smpl.sums);
  smpl.sums <- rev(smpl.sums);
  vip.nms <- names(smpl.sums);
  names(smpl.sums) <- NULL;
  vip.nms <- substr(vip.nms, 1, 16);
  
  myH <- ncol(data_bef)*25 + 50;
  imgName = paste(imgName,".", format, sep="");
  mbSetObj$imgSet$lib.size<-imgName;
  Cairo::Cairo(file=imgName, width=580, height=myH, type=format, bg="white",dpi=dpi);
  xlim.ext <- GetExtendRange(smpl.sums, 10);
  par(mar=c(4,7,2,2));
  dotchart(smpl.sums, col="forestgreen", xlim=xlim.ext, pch=19, xlab="Read Counts", main="Library Size Overview");
  mtext(side=2, at=1:length(vip.nms), vip.nms, las=2, line=1)
  text(x=smpl.sums,y=1:length(smpl.sums),labels= round(smpl.sums), col="blue", pos=4, xpd=T);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj))
}

################################################
###########Phyloseq object creation#############
################################################

#'Function to recreate phyloseq object
#'@description This function recreates the phyloseq object.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param type Input the format of the uploaded files. "text" for 
#'.txt or .csv files, "biom" for .biom, and "mothur" for mothur files.
#'@param taxa_type Input the taxonomy labels. "SILVA" for SILVA taxonomy,
#'"Greengenes" for Greengenes taxonomy, "QIIME" for QIIME taxonomy,
#'"GreengenesID" for Greengenes OTU IDs, and "Others/Not_specific" for others.
#'@param taxalabel Logical, T if taxonomy labels were already included in
#'the OTU table and F if not.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
CreatePhyloseqObj<-function(mbSetObj, type, taxa_type, taxalabel){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_phyloseq();
  load_splitstackshape();
  
  data.proc <- qs::qread("data.proc");
  
  # do some sanity check here on sample and feature names
  smpl.nms <- colnames(data.proc);
  taxa.nms <- row.names(data.proc);
  
  # check for uniqueness of dimension name
  if(length(unique(smpl.nms))!=length(smpl.nms)){
    dup.nms <- paste(smpl.nms[duplicated(smpl.nms)], collapse="; ");
    AddErrMsg(paste(c("Duplicate sample names are not allowed:"), dup.nms, collapse=" "));
    return(0);
  }
  
  # check for uniqueness of taxon names for metagenomic data only
  if(mbSetObj$module.type == "sdp"){
    if(length(unique(taxa.nms))!=length(taxa.nms)){
      dup.tax.nms <- paste(taxa.nms[duplicated(taxa.nms)], collapse="; ");
      AddErrMsg(paste(c("Duplicate taxon names are not allowed:"), dup.tax.nms, collapse=" "));
      return(0);
    }
  }else{
    # for ASV data, replace sequences with ASV1, ASV2 ....
    # check is avegage name length is > 75
    nm.ln <- sum(nchar(taxa.nms)/length(taxa.nms));
    if(nm.ln > 75){
      mbSetObj$is.ASV <-TRUE;
    }
  }
  
  # now check for special characters in the data labels
  if(sum(is.na(iconv(smpl.nms)))>0){
    na.inx <- is.na(iconv(smpl.nms));
    nms <- paste(smpl.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in sample names:", nms, collapse=" "));
    return(0);
  }
  
  var.nms <- rownames(data.proc);
  
  if(sum(is.na(iconv(var.nms)))>0){
    na.inx <- is.na(iconv(var.nms));
    nms <- paste(var.nms[na.inx], collapse="; ");
    AddErrMsg(paste("No special letters (i.e. Latin, Greek) are allowed in feature or taxa names:", nms, collapse=" "));
    return(0);
  }
  
  #standard name to be used
  classi.lvl<- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species");
  
  if(mbSetObj$module.type == "mdp" | mbSetObj$module.type == "ppd"){
    if(type=="text"){
      # prepare data for phyloseq visualization.
      # if features names are present in specific taxonomy format (greengenes or silva).
      if(taxalabel=="T"){ #taxa labels present and need to be parsed in OTU table
        feat_nm <- rownames(data.proc);
        mbSetObj$dataSet$feat_nm <- feat_nm;
        
        if(taxa_type=="SILVA"){
          
          load_splitstackshape();
          
          feat_nm<-data.frame(mbSetObj$dataSet$feat_nm);
          names(feat_nm)<-"Rank";
          taxonomy<-splitstackshape::cSplit(feat_nm,"Rank",";");
          taxmat= data.frame(matrix(NA, ncol = 7, nrow = nrow(taxonomy)));
          colnames(taxmat) <- classi.lvl;
          taxmat[,1:ncol(taxonomy)]<-taxonomy;
          taxmat<-taxmat[colSums(!is.na(taxmat))>0];
          taxmat<-as.matrix(taxmat);
          rownames(taxmat)<-c(1:nrow(taxmat));
          #phyloseq taxonomy object
          taxa_table <- tax_table(taxmat);
          taxa_names(taxa_table)<-rownames(taxmat);
        }else if(taxa_type=="Greengenes"||taxa_type=="QIIME"){
          my.parser<-lapply(feat_nm, parse_taxonomy_qiime); # note, this functions will remove empty strings, replace with NA!?
          taxa_table<-build_tax_table(my.parser);
          
          if(taxa_type=="Greengenes"){
            #additional columns are added(starting with Rank1 etc; has to be removed)
            indx<-grep("^Rank",colnames(taxa_table));
            # need to test if still columns left (could be that users specified wrong tax format)
            if(ncol(taxa_table) == length(indx)){
              AddErrMsg("Error with parsing Greengenes Taxonomy Labels. Please make sure the taxonomy label was specified correctly during data upload.");
              return(0);
            }
            if(length(indx) > 0){
              taxa_table<-taxa_table[,-indx];
            }
            #if all taxa ranks of an taxa are NA with pase_taxonomy_qiime function #covert to Not_assigned
            ind<-which(apply(taxa_table,1, function(x)all(is.na(x)))==TRUE);
            taxa_table[ind,]<-rep("Not_Assigned",ncol(taxa_table));
          }
          taxa_names(taxa_table)<-c(1:nrow(taxa_table));
        }else if(taxa_type=="GreengenesID"||taxa_type=="Others/Not_specific"){
          
          # need to parse Taxonomy still!
          load_splitstackshape();               
          
          feat_nm <- data.frame(mbSetObj$dataSet$feat_nm);
          feat_nm <- data.frame(apply(feat_nm, 1, function(x) gsub(";\\s;", ";", x))) # remove empty taxa
          names(feat_nm) <- "Rank";
          taxonomy <- splitstackshape::cSplit(feat_nm,"Rank",";");
          taxmat <- data.frame(matrix(NA, ncol = 7, nrow = nrow(taxonomy)));
          colnames(taxmat) <- classi.lvl;
          taxmat[,1:ncol(taxonomy)] <- taxonomy;
          taxmat <- taxmat[colSums(!is.na(taxmat))>0];
          taxmat <- as.matrix(taxmat);
          rownames(taxmat) <- c(1:nrow(taxmat));
          #phyloseq taxonomy object
          taxa_table <- tax_table(taxmat);
          taxa_names(taxa_table) <- rownames(taxmat);
        }
        
        # making unique id for each OTU consist of lowest taxonomy level present followed by row number
        tmat <- as(tax_table(taxa_table), "matrix")
        rowLn <- nrow(tmat);
        colLn <- ncol(tmat);
        num.mat <- matrix(0, nrow=rowLn, ncol=colLn);
        
        for(i in 1:colLn){
          num.mat[,i] <- rep(i, rowLn);
        }
        
        na.inx <- is.na(tmat);
        num.mat[na.inx] <- 0;
        maxInx <- apply(num.mat, 1, max);
        gd.coords <- cbind(1:rowLn, maxInx);
        gd.nms <- tmat[gd.coords];
        mynames <- make.unique(gd.nms, sep="_");
        
        if(length(mynames) != nrow(tmat)){
          AddErrMsg(paste("Errors with the taxonomy table! Please check if there are any names with all NA."));
          return(0);
        }
        
        taxa_names(taxa_table)<-mynames;
        feat_nm<-NULL;
        mbSetObj$dataSet$taxa_table<-taxa_table;
      }else{
        #sanity check: features name of taxonomy table should match feature names in OTU abundance table,
        # Only features in filtered OTU tables are selected additional are removed.
        taxa_table <- mbSetObj$dataSet$taxa_table;
        indx<-match(rownames(data.proc), rownames(taxa_table));
        
        if(sum(is.na(indx)) > 0){
          na.nms <- rownames(data.proc)[is.na(indx)];
          if(length(na.nms)>10){
            na.nms <- na.nms[1:10];
          }
          AddErrMsg(paste("The following names cannot be found in your taxanomy table (showing 10 max):", paste(na.nms, collapse="; ")));
          return(0);
        }
        
        # let's address ASV sequence ID here
        if(mbSetObj$is.ASV){
          # do the conversion and save 
          new.nms <- paste("ASV", 1:length(taxa.nms), sep="_");
          master.nms <- cbind("NewID"=new.nms, "OriginalID"=taxa.nms);
          fast.write(master.nms, file="ASV_ID_mapping.csv");
          
          # update the taxa and sample
          names(new.nms) <- taxa.nms;
          rownames(data.proc) <- new.nms[rownames(data.proc)];
          rownames(taxa_table) <- new.nms[rownames(taxa_table)];
          
          # update tree file if uploaded
          if(mbSetObj$tree.uploaded){
            pg_tree <- qs::qread("tree.qs");
            pg_tree$tip.label <- new.nms[pg_tree$tip.label];
            qs::qsave(pg_tree, "tree.qs");
          }
        }
        
        nm<-colnames(taxa_table);
        taxa_table<-as.matrix(taxa_table[indx,]);
        colnames(taxa_table)<-nm;
        #converting taxonomy file(data frame) to phyloseq taxonomy object;keeping the taxonomy names as it is.
        mbSetObj$dataSet$taxa_table<-tax_table(taxa_table);
        taxa_nms <- rownames(taxa_table)
        
        if(length(unique(taxa_nms))!=length(taxa_nms)){
          dup.tax.nms <- paste(taxa_nms[duplicated(taxa_nms)], collapse="; ");
          AddErrMsg(paste(c("Duplicate taxon names are not allowed:"), dup.tax.nms, collapse=" "));
          return(0);
        }
        taxa_names(mbSetObj$dataSet$taxa_table)<-rownames(taxa_table);
      }
      
      # creating phyloseq object
      data.proc<-apply(data.proc,2,as.integer);
      data.proc<-otu_table(data.proc,taxa_are_rows =TRUE);
      taxa_names(data.proc)<-taxa_names(mbSetObj$dataSet$taxa_table);
      # removing constant column and making standard names
      ind<-apply(mbSetObj$dataSet$taxa_table[,1], 2, function(x) which(x %in% c("Root", "root")));
      
      if(length(ind)>0){
        mbSetObj$dataSet$taxa_table<-mbSetObj$dataSet$taxa_table[,-1];
        #dataSet$taxa_table<-dataSet$taxa_table #newnew
      }
      
      # removing first KINGDOM OR DOMAIN column
      indx<-apply(mbSetObj$dataSet$taxa_table[,1], 2, function(x) which(x %in% c("k__Fungi", "Fungi", "k__Viruses", "Viruses",
                                                                                 "Bacteria", "Archaea", "k__Bacteria", "k__Archea", 
                                                                                 "D_0__Bacteria", "D_0__Archea", "d__Bacteria")));
      
      if(length(indx)>0){
        mbSetObj$dataSet$taxa_table<-mbSetObj$dataSet$taxa_table[,-1];
      }
      
      classi.lvl<- c("Phylum", "Class", "Order", "Family", "Genus", "Species");
      colnames(mbSetObj$dataSet$taxa_table)<-classi.lvl[1:ncol(mbSetObj$dataSet$taxa_table)]; 
      #sanity check: no of sample of both abundance and metadata should match.
      indx<-match(colnames(data.proc), rownames(mbSetObj$dataSet$sample_data));
      
      if(all(is.na(indx))){
        AddErrMsg("Please make sure that sample names and their number are same in metadata and OTU abundance files.");
        return(0);
      }else{
        mbSetObj$dataSet$sample_data<-mbSetObj$dataSet$sample_data[indx, ,drop=FALSE];
      }
      
      mbSetObj$dataSet$sample_data<-sample_data(mbSetObj$dataSet$sample_data, errorIfNULL=TRUE);
      #cleaning up the names#deleting [] if present; and substituting (space,.,/ with underscore(_))
      mbSetObj$dataSet$taxa_table<- gsub("[[:space:]./_-]", "_",mbSetObj$dataSet$taxa_table);
      mbSetObj$dataSet$taxa_table<- gsub("\\[|\\]","",mbSetObj$dataSet$taxa_table);
      
      #sometimes after removal of such special characters rownames beacame non unique; so make it unique
      mynames1<- gsub("[[:space:]./_-]", "_",taxa_names(mbSetObj$dataSet$taxa_table));
      mynames1<-gsub("\\[|\\]","",mynames1);
      
      if(length(unique(mynames1))< length(taxa_names(mbSetObj$dataSet$taxa_table))){
        mynames1 <- make.unique(mynames1, sep="_");
      }
      
      #after cleaning, names might now be NA
      missing.inx <- which(is.na(mynames1))
      missing_vector <- paste0("New_OTU", seq_along(missing.inx))
      mynames1[missing.inx] <- missing_vector
      
    } else if(type == "biom"||type == "mothur"){
      # creating phyloseq object
      feat_nm<-rownames(data.proc);
      data.proc<-apply(data.proc,2,as.integer);
      data.proc<-otu_table(data.proc,taxa_are_rows =TRUE);
      taxa_names(data.proc)<-feat_nm;
      #sanity check: features name of taxonomy table should match feature names in OTU abundance table,Only features in filtered OTU tables are selected additional are removed.
      indx<-match(taxa_names(data.proc),taxa_names(mbSetObj$dataSet$taxa_table));
      
      if(all(is.na(indx))){
        AddErrMsg("Please make sure that features name of the taxonomy table match the feature names in the OTU abundance table.");
        return(0);
      }
      
      mbSetObj$dataSet$taxa_table<-mbSetObj$dataSet$taxa_table[indx,];
      
      # removing constant column (root) if present otherwise just kingdom column from taxonomy table
      ind<-apply(mbSetObj$dataSet$taxa_table[,1], 2, function(x) which(x %in% c("Root", "root")));
      
      if(length(ind)>0){
        mbSetObj$dataSet$taxa_table<-mbSetObj$dataSet$taxa_table[,-1];
        #dataSet$taxa_table<-dataSet$taxa_table #newnew
      }
      
      # removing first KINGDOM OR DOMAIN column
      indx<-apply(mbSetObj$dataSet$taxa_table[,1], 2, function(x) which(x %in% c("k__Fungi", "Fungi", "k__Viruses", "Viruses", 
                                                                                 "Bacteria","Archaea","k__Bacteria","k__Archea", 
                                                                                 "D_0__Bacteria","D_0__Archea", "d__Bacteria")));
      
      if(length(indx)>0){
        mbSetObj$dataSet$taxa_table<-mbSetObj$dataSet$taxa_table[,-1];
      }
      
      classi.lvl<- c("Phylum", "Class", "Order", "Family", "Genus", "Species"); 
      colnames(mbSetObj$dataSet$taxa_table)<-classi.lvl[1:ncol(mbSetObj$dataSet$taxa_table)];
      #sanity check: no of sample of both abundance and metadata should match.
      indx<-match(colnames(data.proc), rownames(mbSetObj$dataSet$sample_data));
      
      if(all(is.na(indx))){
        AddErrMsg("Please make sure that sample names and their number are the same in metadata and OTU abundance files.");
        return(0);
      }else{
        mbSetObj$dataSet$sample_data<-mbSetObj$dataSet$sample_data[indx, ,drop=FALSE];
      }
      
      #cleaning up the names#deleting [] if present; and substituting (space,.,/ with underscore(_))
      mbSetObj$dataSet$taxa_table <- gsub("[[:space:]./_-]", "_",mbSetObj$dataSet$taxa_table);
      mbSetObj$dataSet$taxa_table<- gsub("\\[|\\]","",mbSetObj$dataSet$taxa_table);
      
      #sometimes after removal of such special characters rownames beacame non unique; so make it unique
      mynames1 <- gsub("[[:space:]./_-]", "_",taxa_names(mbSetObj$dataSet$taxa_table));
      mynames1<-gsub("\\[|\\]","",mynames1);
      
      if(length(unique(mynames1))< length(taxa_names(mbSetObj$dataSet$taxa_table))){
        mynames1 <- make.unique(mynames1, sep="_");
      }
      
      #after cleaning, names might now be NA
      missing.inx <- which(is.na(mynames1))
      missing_vector <- paste0("New_OTU", seq_along(missing.inx))
      mynames1[missing.inx] <- missing_vector
    }
    
    taxa_names(mbSetObj$dataSet$taxa_table)<-taxa_names(data.proc)<-mynames1;
    mbSetObj$dataSet$sample_data<-sample_data(mbSetObj$dataSet$sample_data, errorIfNULL = TRUE);
    
    sd_names <- sample_names(mbSetObj$dataSet$sample_data)
    name.match.inx <- sd_names %in% smpl.nms
    
    if(sum(as.numeric(name.match.inx)) == 0){
      AddErrMsg("Issues with sample names in your files! Make sure names are not purely numeric (i.e. 1, 2, 3).");
      return(0);
    }
    
    mbSetObj$dataSet$proc.phyobj <- merge_phyloseq(data.proc, mbSetObj$dataSet$sample_data, mbSetObj$dataSet$taxa_table);
    
    if(length(rank_names(mbSetObj$dataSet$proc.phyobj)) > 7){
      tax_table(mbSetObj$dataSet$proc.phyobj) <- tax_table(mbSetObj$dataSet$proc.phyobj)[, 1:7]
    }
    
    #also using unique names for our further data
    prefilt.data <- qs::qread("data.prefilt");
    rownames(prefilt.data)<-taxa_names(mbSetObj$dataSet$taxa_table);
    qs::qsave(prefilt.data, file="data.prefilt")
    
  }else if(mbSetObj$module.type == "sdp"){
    #constructing phyloseq object for aplha diversity.
    data.proc<-otu_table(data.proc, taxa_are_rows = TRUE);
    taxa_names(data.proc)<-rownames(data.proc);
    
    #sanity check: sample names of both abundance and metadata should match.
    smpl_var<-colnames(mbSetObj$dataSet$sample_data);
    indx<-match(colnames(data.proc), rownames(mbSetObj$dataSet$sample_data));
    
    if(all(is.na(indx))){
      AddErrMsg("Please make sure that sample names are matched between both files.");
      return(0);
    }
    mbSetObj$dataSet$sample_data<-sample_data(mbSetObj$dataSet$sample_data, errorIfNULL=TRUE);
    mbSetObj$dataSet$proc.phyobj <- merge_phyloseq(data.proc, mbSetObj$dataSet$sample_data);
  }
  
  qs::qsave(mbSetObj$dataSet$proc.phyobj, file="proc.phyobj.orig");
  #do some data cleaning
  mbSetObj$dataSet$data.prefilt <- NULL;
  mbSetObj$dataSet$taxa_type <- taxa_type;
  
  return(.set.mbSetObj(mbSetObj));
}

###########################################
## Utilty function to perform edgeR norm ##
###########################################
edgeRnorm = function(x,method){
  # Enforce orientation.
  # See if adding a single observation, 1,
  # everywhere (so not zeros) prevents errors
  # without needing to borrow and modify
  # calcNormFactors (and its dependent functions)
  # It did. This fixed all problems.
  # Can the 1 be reduced to something smaller and still work?
  x = x + 1;
  # Now turn into a DGEList
  y = DGEList(counts=x, remove.zeros=TRUE);
  # Perform edgeR-encoded normalization, using the specified method (...)
  z = edgeR::calcNormFactors(y, method=method);
  # A check that we didn't divide by zero inside `calcNormFactors`
  if( !all(is.finite(z$samples$norm.factors)) ){
    AddErrMsg(paste("Something wrong with edgeR::calcNormFactors on this data, non-finite $norm.factors."));
    return(0);
  }
  return(z)
}

###########################################
##################Getters##################
###########################################

GetSamplNames<-function(){
  prefilt.data <- qs::qread("data.prefilt");
  return(colnames(prefilt.data));
}

GetGroupNames<-function(){
  prefilt.data <- qs::qread("data.prefilt");
  return(colnames(prefilt.data));
}

GetSampleNamesaftNorm<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(sample_names(mbSetObj$dataSet$norm.phyobj));
}

