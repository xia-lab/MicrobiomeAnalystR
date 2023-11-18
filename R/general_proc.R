##################################################
## R script for MicrobiomeAnalyst
## Description: filtering and normalization functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#######################################
# filter zeros and singletons this is not good enough ==> single occurance
# reads occur in only one sample should be treated as artifacts and removed

#'Main function to sanity check on uploaded microbiome data
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
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type;
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
  
  saveDataQs(data.proc, "data.proc.orig", module.type, dataName);
  saveDataQs(data.proc, "data.prefilt", module.type, dataName);

  if(module.type == "meta"){
    sanCheck <- .checkSampleConsistencyMeta(mbSetObj,dataName);
    if(sanCheck == 0){
        return(0);
    }
  }

  # now get stats
  taxa_no <- nrow(mbSetObj$dataSet$data.orig);
  
  # from abundance table
  sample_no <- ncol(mbSetObj$dataSet$data.orig);
  
  if(any(c(filetype=="biom", filetype=="mothur"))){
    samplemeta_no <- nrow(mbSetObj$dataSet$sample_data);
  }else{
    samplemeta_no <- sample_no;
  }
  
  if(!mbSetObj$module.type %in% c("ppd") ){
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
    tree<-qs::qread("tree.qs")  
    
    if(length(intersect(rownames(data.proc),tree$tip.label))==nrow(data.proc)){
      tree_tip <- 1
    }else{
      tree_tip <- 0
      AddErrMsg("The tip labels of the tree are not matched with your feature names in the abundance table!");
      
    }
  } else {
    tree_exist <- 0
    tree_tip <- 0
  }
  if(identical(sort(row.names(mbSetObj$dataSet$sample_data)),
               sort(colnames(data.proc)))){
    samname_same <- 1;
  } else {
    samname_same <- 0;
  }
  
  sample_no_in_outfile <- ncol(data.proc);
  samname_same_number <- sum(row.names(mbSetObj$dataSet$sample_data) %in% colnames(data.proc));
   if(samname_same_number){
      mis.num <- length(which(!row.names(mbSetObj$dataSet$sample_data) %in% colnames(data.proc)))
     if(mis.num==1){
     current.msg <<- paste0("One sample name was not included in the OTU/ASV abundance table.");

   }else{
   current.msg <<- paste0("A total of ",mis.num," sample names were not included in the OTU/ASV abundance table.");

   }

  }
  
  # now store data.orig to RDS
  saveDataQs(mbSetObj$dataSet$data.orig, "data.orig", module.type, dataName);
  saveDataQs(data.proc,"data.proc", module.type, dataName);
  saveDataQs(mbSetObj$dataSet$sample_data, "data.sample_data", module.type, dataName);
  
  #mbSetObj$dataSet$data.orig <- NULL;
  mbSetObj$dataSet$tree <- tree_exist
  
  vari_no <- ncol(mbSetObj$dataSet$sample_data);
  disc_no <- sum(mbSetObj$dataSet$meta_info$disc.inx);
  cont_no <- sum(mbSetObj$dataSet$meta_info$cont.inx);
  vari_dup_no <- disc_no + cont_no
  smpl.sums <- apply(data.proc, 2, sum);
  tot_size <- sum(smpl.sums);
  smin <- min(smpl.sums)
  smean <- mean(smpl.sums)
  smax <- max(smpl.sums);
  gd_feat <- nrow(data.proc);
  
  if(exists("current.proc")){
  current.proc$mic$data.proc<<-data.proc
}

  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(c(1,taxa_no,sample_no,vari_no, smin,smean,smax,gd_feat,samplemeta_no,tot_size, tree_exist, samname_same, samname_same_number, sample_no_in_outfile, disc_no, cont_no,vari_dup_no,tree_tip));
  }else{
    print("Sanity check passed!")
    
    return(.set.mbSetObj(mbSetObj))
  }
}

#'Function to sanity check on uploaded metabolomics data
#'@description This function performs a sanity check on the uploaded
#'user data. It checks the grouping of samples, if a phylogenetic tree
#'was uploaded.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param disableFilter Booleanmetabo. Set to TRUE to bypass the hard-coded filter
#'to remove features occuring in <2 samples.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export


SanityCheckMetData <- function(mbSetObj,isNormMetInput, disableFilter = FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  feat.sums <- apply(mbSetObj$dataSet$metabolomics$data.orig, 1, function(x){sum(x>0, na.rm=T)});
  
  if(disableFilter){
    data.proc <-mbSetObj$dataSet$metabolomics$data.orig
  }else{
    gd.inx <- feat.sums > 1; # occur in at least 2 samples
    if(length(which(gd.inx=="TRUE"))==0){
      AddErrMsg("Reads occur in only one sample.  All these are considered as artifacts and have been removed from data. No data left after such processing.");
      return(0);
    }
    data.proc <- mbSetObj$dataSet$metabolomics$data.orig[gd.inx, ]; 
  }
  
  # filtering the constant features here
  # check for columns with all constant (var=0)
  varCol <- apply(data.proc, 1, var, na.rm=T);
  constCol <- varCol == 0 | is.na(varCol);
  
  # making copy of data.proc
  data.proc <- data.proc[!constCol, ];
  
  if(length(data.proc)==0){
    AddErrMsg("All features are found to be constant and have been removed from data. No data left after such processing.");
    return(0);
  }
  
  # check zero, NA values
  totalCount <- nrow(data.proc)*ncol(data.proc);
  naCount <- sum(is.na(data.proc));
  naPercent <- round(100*naCount/totalCount,1)
  #  print(naCount)
  mbSetObj$dataSet$metabolomics$missingCount <- naCount;
  
  msg<-paste("A total of ", naCount, " (", naPercent, "%) missing values were detected.", sep="");
  
  
  #Reset to default
  mbSetObj$dataSet$metabolomics$data.filt<-  mbSetObj$dataSet$metabolomics$edit <- NULL;
  
  # replace zero and missing values using Detection Limit for each variable 
  int.mat <- t(ReplaceMissingByLoD(t(data.proc)));   ##notice that samples in column and features in row
  
  mbSetObj$dataSet$metabolomics$proc.feat.num <- ncol(int.mat);
  
  msg<-c(msg, paste("Zero or missing values were replaced by 1/5 of the min positive value for each variable."));
  
  
  mbSetObj$dataSet$metabolomics$check.msg <- msg;
  mbSetObj$dataSet$metabolomics$isNormInput <- isNormMetInput;
  samname_same_number <- sum(row.names(mbSetObj$dataSet$sample_data) %in% colnames(data.proc));
  # now store data.orig to RDS
  qs::qsave(mbSetObj$dataSet$metabolomics$data.orig, file="metabo.data.orig.qs");
  qs::qsave(int.mat, file="metabo.data.init.qs");
  current.proc$met$data.proc<<-int.mat

  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(c(1,length(feat.sums),nrow(data.proc),naPercent,samname_same_number));
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
  
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type;
  data <- readDataQs("data.prefilt", module.type, dataName);
  #data <- qs::qread("data.prefilt");
  
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
  
  if(sum(kept.inx)==0){
    AddErrMsg(paste0("No samples have counts less than ", count, "! If your data has already been normalized, please turn off count filter (0)."))
    return(0)
  }
  
  mbSetObj$dataSet$filt.data <- data[kept.inx, ];
  
  if(exists("current.proc")){
    
    current.proc$mic$data.proc<<-mbSetObj$dataSet$filt.data

  }
  
  
  saveDataQs(mbSetObj$dataSet$filt.data, "filt.data.orig", module.type, dataName);
  #qs::qsave(mbSetObj$dataSet$filt.data, file="filt.data.orig"); # save an copy
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
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type;
  data <- mbSetObj$dataSet$filt.data;
  
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
  
  if(nrow(data) == 0){
    AddErrMsg("No samples remaining after variance filtering!")
    return(0)
  }
  
  mbSetObj$dataSet$filt.data <- data;
  if(exists("current.proc")){ 
    current.proc$mic$data.proc<<-data
  }
  
  #qs::qsave(mbSetObj$dataSet$filt.data, file="filt.data.orig"); # save an copy
  saveDataQs(mbSetObj$dataSet$filt.data, "filt.data.orig", module.type, dataName);
  
  rm.msg1 <- paste("A total of ```", sum(!remain), "``` low variance features were removed based on ```", filtopt, "```.", sep="");
  rm.msg2 <- paste("The number of features remains after the data filtering step: ```", nrow(data), "```.");
  current.msg <<- c(current.msg, rm.msg1, rm.msg2);
  mbSetObj$dataSet$filt.msg <- current.msg;
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to filter uploaded metabolomic data
#'@description .
#'Note, this is applied after abundance filter.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param filtopt Character, input the low variance filter option. "iqr" for 
#'inter-quantile range, "sd" for standard deviation, and "cov" for coefficient of variation.
#'@param filtPerct Numeric, input the percentage cutoff for low variance. 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

ApplyMetaboFilter <- function(mbSetObj=NA, filter,  rsd){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data <- mbSetObj$dataSet$metabolomics$data.orig;
  int.mat <- t(data) ;
  
  
  feat.num <- ncol(int.mat);
  feat.nms <- colnames(int.mat);
  nm <- NULL;
  msg <- "";
  if(all(c(filter == "none", feat.num < 5000))) { # only allow for less than 4000
    remain <- rep(TRUE, feat.num);
    msg <- paste(msg, "No filtering was applied");
  }else {
    if (filter == "rsd"){
      sds <- apply(int.mat, 2, sd, na.rm=T);
      mns <- apply(int.mat, 2, mean, na.rm=T);
      filter.val <- abs(sds/mns);
      nm <- "Relative standard deviation";
    }else if (filter == "nrsd" ){
      mads <- apply(int.mat, 2, mad, na.rm=T);
      meds <- apply(int.mat, 2, median, na.rm=T);
      filter.val <- abs(mads/meds);
      nm <- "Non-paramatric relative standard deviation";
    }else if (filter == "mean"){
      filter.val <- apply(int.mat, 2, mean, na.rm=T);
      nm <- "mean";
    }else if (filter == "sd"){
      filter.val <- apply(int.mat, 2, sd, na.rm=T);
      nm <- "standard deviation";
    }else if (filter == "mad"){
      filter.val <- apply(int.mat, 2, mad, na.rm=T);
      nm <- "Median absolute deviation";
    }else if (filter == "median"){
      filter.val <- apply(int.mat, 2, median, na.rm=T);
      nm <- "median";
    }else{ # iqr
      filter.val <- apply(int.mat, 2, IQR, na.rm=T);
      nm <- "Interquantile Range";
    }
    
    # get the rank of the filtered variables
    rk <- rank(-filter.val, ties.method='random');
    
    
    if(feat.num < 250){ # reduce 5%
      remain <- rk < feat.num*0.95;
      msg <- paste(msg, "Further feature filtering based on", nm,". A total of ",sum(!remain), "features were removed.");
    }else if(feat.num < 500){ # reduce 10%
      remain <- rk < feat.num*0.9;
      msg <- paste(msg, "Further feature filtering based on", nm,". A total of ",sum(!remain), "features were removed.");
    }else if(feat.num < 1000){ # reduce 25%
      remain <- rk < feat.num*0.75;
      msg <- paste(msg, "Further feature filtering based on", nm,". A total of ",sum(!remain), "features were removed.");
    }else{ # reduce 40%, if still over 5000, then only use top 5000
      remain <- rk < feat.num*0.6;
      msg <- paste(msg, "Further feature filtering based on", nm);
      
      if(sum(remain) > 5000){
        remain <- rk < 5000;
        msg <- paste(msg, paste("Reduced to", 5000, "features based on", nm,""));
      }
    }
    
  }
  
  # save a copy for user 
  #fast.write.csv(cbind(filter=filter.val, t(int.mat)), file=paste0("data_prefilter_", filter, ".csv"));
  
  
  filt.res=int.mat[, remain]
  mbSetObj$dataSet$metabolomics$filt.data <- t(filt.res); 
  current.proc$met$data.proc<<-t(filt.res)
  qs::qsave(mbSetObj$dataSet$metabolomics$filt.data, file="metabo.filt.data"); # save an copy
  current.msg <<- msg
 mbSetObj$dataSet$metabolomics$filt.msg <- current.msg;
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
  
  load_phyloseq();
  mbSetObj <- .get.mbSetObj(mbSetObj);
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type;
  
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
  prefilt.data <- readDataQs("data.prefilt", module.type, dataName);
  tax_nm <- rownames(prefilt.data);
  data.proc <- readDataQs("data.proc", module.type, dataName);
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
PerformNormalization <- function(mbSetObj, rare.opt, scale.opt, transform.opt,isAutoScale=F){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type;
  data <- readDataQs("filt.data.orig", module.type, dataName);
  saveDataQs(mbSetObj$dataSet$proc.phyobj, "orig.phyobj", module.type, dataName);  
  # now proc.phyobj is now filtered data
  tax_nm <- rownames(data);
  msg <- NULL;
  
  #rarefying (should?) affect filt.data in addition to norm 
  
  if(rare.opt != "none"){
    data <- PerformRarefaction(mbSetObj, data, rare.opt);
    tax_nm <- rownames(data);
    msg <- c(msg, paste("Performed data rarefaction."));
    
    ###### note, rarefying (should?) affect filt.data in addition to norm
    mbSetObj$dataSet$filt.data <- data;
  }else{
    msg <- c(msg, paste("No data rarefaction was performed."));
  }

  # create phyloseq obj
  mbSetObj$dataSet$sample_data$sample_id <- rownames(mbSetObj$dataSet$sample_data);
  sample_table <- sample_data(mbSetObj$dataSet$sample_data, errorIfNULL=TRUE);
  mbSetObj$dataSet$proc.phyobj<- merge_phyloseq(otu_table(data,taxa_are_rows =TRUE), sample_table, mbSetObj$dataSet$taxa_table);
 
  
  #make hierarchies
  if(mbSetObj$module.type=="sdp" | mbSetObj$micDataType=="ko"){
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
  
  ##### The prenorm object is the filter result of the raw counts and also affected by the rarefy option 
### Note the counts table have different row numbers with the merged_obj as the otus without taxonomy assignment at a certain level will be merged in the count table.
  saveDataQs(data.list, "phyloseq_prenorm_objs.qs", module.type, dataName);  
 
  
  for(i in 1:length(ranks)){
    
    data <- data.list$count_tables[[i]]
     if(scale.opt != "none"){
      if(scale.opt=="colsum"){
        data <- sweep(data, 2, colSums(data), FUN="/")
        data <- data*10000000;
        if(i==1){
          msg <- c(msg, paste("Performed ```total sum scaling``` normalization."));
        }
      }else if(scale.opt=="upperquartile"){
        load_edgeR();
        otuUQ <- edgeRnorm(data,method="upperquartile");
        data <- as.matrix(otuUQ$counts);
        if(i==1){
          msg <- c(msg, paste("Performed ```upper quartile``` normalization"));
        }
      }else if(scale.opt=="CSS"){
        load_metagenomeseq();
        #biom and mothur data also has to be in class(matrix only not in phyloseq:otu_table)
        data1 <- as(data,"matrix");
        dataMR <- newMRexperiment(data1);
        data <- cumNorm(dataMR,p=cumNormStat(dataMR));
        data <- MRcounts(data,norm = T);
        if(i==1){
          msg <- c(msg, paste("Performed ```cumulative sum scaling``` normalization"));
        }
      }else{
        if(scale.opt!="auto" & i==1){
          print(paste("Unknown scaling parameter:", scale.opt));
        }
      }
    }else{
      if(i==1){
        msg <- c(msg, paste("No data scaling was performed."));
      }
    }
    
    if(transform.opt != "none"){
      if(transform.opt=="rle"){
        load_edgeR();            
        otuRLE <- edgeRnorm(data,method="RLE");
        data <- as.matrix(otuRLE$counts);
        if(i==1){        
          msg <- c(msg, paste("Performed ```RLE``` Normalization"));
        }
      }else if(transform.opt=="TMM"){
        load_edgeR();   
        otuTMM <- edgeRnorm(data,method="TMM");
        data <- as.matrix(otuTMM$counts);
        if(i==1){        
          msg <- c(msg, paste("Performed ```TMM``` Normalization"));
        }
      }else if(transform.opt=="clr"){
        data <- apply(data, 2, clr_transform);
        if(i==1){        
          msg <- "Performed ```centered-log-ratio (CLR)``` normalization."
        }
      }else{
        if(scale.opt!="auto" & i==1){
          print(paste("Unknown scaling parameter:", transform.opt));
        }
      }
    }else{
      if(i==1){        
        msg <- c(msg, paste("No data transformation was performed."));
      }
    }
   data.list$count_tables[[i]] <- data

    if(ranks[i]=="OTU"){
      otu.tab <- otu_table(data,taxa_are_rows =TRUE);
      data.list$merged_obj[[i]] <- merge_phyloseq(otu.tab, sample_data(data.list$merged_obj[[i]]), mbSetObj$dataSet$taxa_table);
    }else{
      otu.tab <- otu_table(data,taxa_are_rows =TRUE);
      tax.tab <- tax_table(data.list$merged_obj[[i]])
      nm <- as.character(tax.tab[,ranks[i]])
      nm[is.na(nm)|nm==""|nm==" "] <- "Not_Assigned";
      rownames(tax.tab) <- nm
      tax.tab <- tax.tab[!(duplicated(nm)),] 
      tax.tab[which(rownames(tax.tab)=="Not_Assigned"),] <- NA
      data.list$merged_obj[[i]] <- merge_phyloseq(otu.tab, sample_data(data.list$merged_obj[[i]]), tax.tab);
    }
  }

  #using the OTU level for plotting when first initialize
  mbSetObj$dataSet$norm.phyobj <- data.list$merged_obj[["OTU"]];
  mbSetObj$dataSet$norm.msg <- msg;
  
  ### prepare for mmp module
  if(exists("current.proc")){
    qs::qsave(data.list,"prescale.phyobj.qs")
    if(isAutoScale=="true" | scale.opt == "auto"){
      data.list$count_tables <- lapply(data.list$count_tables,function(df) t(AutoScale(t(df))))
    }  
    current.proc$mic$data.proc<<- data.list$count_tables[["OTU"]]
    
  }
  
  ### Normalized object for each taxonomy level
  saveDataQs(data.list, "phyloseq_objs.qs", module.type, dataName); 

  return(.set.mbSetObj(mbSetObj));
}


#'Function to perform normalization on metabolomics data
#'@description This function performs normalization on the uploaded metabolomics
#'data.
#@param rowNorm Select the option for row-wise normalization, "QuantileNorm" for Quantile Normalization, 
#'"SumNorm" for Normalization to constant sum, 
#'"MedianNorm" for Normalization to sample median, and 
#'@param transNorm Select option to transform the data, "LogNorm" for Log Normalization,
#'and "CrNorm" for Cubic Root Transformation. 
#'@param scaleNorm Select option for scaling the data, "MeanCenter" for Mean Centering,
#'"AutoNorm" for Autoscaling, "ParetoNorm" for Pareto Scaling, amd "RangeNorm" for Range Scaling.
#'@param ref Input the name of the reference sample or the reference feature, use " " around the name.  
#'@param ratio This option is only for biomarker analysis.
#'@param ratioNum Relevant only for biomarker analysis.  
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import edgeR
#'@import metagenomeSeq
PerformMetaboNormalization <- function(mbSetObj, rowNorm, transNorm, scaleNorm,isAutoScale){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  data <-t(mbSetObj$dataSet$metabolomics$filt.data);
  # now proc.phyobj is now filtered data

  msg <- "";
  
  colNames <- colnames(data);
  rowNames <- rownames(data);
  
  # row-wise normalization
  if(rowNorm=="QuantileNorm"){
    
    data<-t(preprocessCore::normalize.quantiles(t(data), copy=FALSE));
    # this can introduce constant variables if a variable is 
    # at the same rank across all samples (replaced by its average across all)
    varCol <- apply(data, 2, var, na.rm=T);
    constCol <- (varCol == 0 | is.na(varCol));
    constNum <- sum(constCol, na.rm=T);
    if(constNum > 0){
      print(paste("After quantile normalization", constNum, "features with a constant value were found and deleted."));
      data <- data[,!constCol, drop=FALSE];
      colNames <- colnames(data);
      rowNames <- rownames(data);
    }
    rownm<-"Quantile Normalization";
  }else if(rowNorm=="SumNorm"){
    data<-t(apply(data, 1, function(x){
      1000*x/sum(x, na.rm=T);
    }));
    rownm<-"Normalization to constant sum";
  }else if(rowNorm=="MedianNorm"){
    data<-t(apply(data, 1, function(x){
      x/median(x, na.rm=T);
    }));
    rownm<-"Normalization to sample median";
  }else{
    # nothing to do
    rownm<-"N/A";
  }
  
  if(rownm!="N/A"){
    msg <- paste(rownm,"performed.")
    
  }
  # use apply will lose dimension info (i.e. row names and colnames)
  rownames(data)<-rowNames;
  colnames(data)<-colNames;
  
  # record row-normed data for fold change analysis (b/c not applicable for mean-centered data)
  row.norm <- as.data.frame(CleanData(data, T, T)); #moved below ratio 
  qs::qsave(row.norm, file="metabo.row.norm.qs");
  # transformation
  # may not be able to deal with 0 or negative values
  if(transNorm=='LogNorm'){
    min.val <- min(abs(data[data!=0]))/10;
    data<-apply(data, 2, LogNorm, min.val);
    transnm<-"Log10 Normalization";
  }else if(transNorm=='SrNorm'){
    min.val <- min(abs(data[data!=0]))/10;
    data<-apply(data, 2, SquareRootNorm, min.val);
    transnm<-"Square Root Transformation";
  }else if(transNorm=='CrNorm'){
    norm.data <- abs(data)^(1/3);
    norm.data[data<0] <- - norm.data[data<0];
    data <- norm.data;
    transnm<-"Cubic Root Transformation";
  }else{
    transnm<-"N/A";
  }
  
  if(transnm!="N/A"){
    msg <- paste(msg,paste(transnm,"performed."))
    
  }
  # scaling
  if(scaleNorm=='MeanCenter'){
    data<-apply(data, 2, function(x){
      x - mean(x);
    });
    scalenm<-"Mean Centering";
  }else if(scaleNorm=='AutoNorm'| isAutoScale == "true"){
    data<-apply(data, 2, function(x){
      (x - mean(x))/sd(x, na.rm=T);
    });
    scalenm<-"Autoscaling";
  }else if(scaleNorm=='ParetoNorm'){
    data<-apply(data, 2, function(x){
      (x - mean(x))/sqrt(sd(x, na.rm=T));
    });
    scalenm<-"Pareto Scaling";
  }else if(scaleNorm=='RangeNorm'){
    data<-apply(data, 2, function(x){
      if(max(x) == min(x)){
        x;
      }else{
        (x - mean(x))/(max(x)-min(x));
      }
    });
    scalenm<-"Range Scaling";
  }else{
    scalenm<-"N/A";
  }
  
  if(scalenm!="N/A"){
    msg <- paste(msg,paste(scalenm,"performed."))
    
  }
  # note after using "apply" function, all the attribute lost, need to add back
  rownames(data)<-rowNames;
  colnames(data)<-colNames;

  # need to do some sanity check, for log there may be Inf values introduced
  data <- CleanData(data, T, F);
  
  
  mbSetObj$dataSet$metabolomics$norm.data <- t(as.data.frame(data)); ## feature in row and sample in column
  
  qs::qsave( mbSetObj$dataSet$metabolomics$norm.data, file="metabo.complete.norm.qs");
  
  current.proc$met$data.proc<<-t(as.data.frame(data))
  
  mbSetObj$dataSet$metabolomics$rownorm.method <- rownm;
  mbSetObj$dataSet$metabolomics$trans.method <- transnm;
  mbSetObj$dataSet$metabolomics$scale.method <- scalenm;
  
  #using this object for plotting
  if(rownm=="N/A" & transnm=="N/A" & scalenm=="N/A"){
    current.msg <<- "No normalization was performed";
  }else{
    current.msg <<-msg
  }
  
  mbSetObj$dataSet$metabolomics$norm.msg <- current.msg;
  
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
  msg <- paste(msg, collapse=" ");
  #print(msg);
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
  grp.data <- data.frame(grp.data,check.names=FALSE);
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
PlotLibSizeView <- function(mbSetObj, origImgName="",format="png", dpi=72, dataName="", interactive=F){
  mbSetObj <- .get.mbSetObj(mbSetObj)
  library(ggplot2)
  library(dplyr)
  library(Cairo)
  
  if(dataName != ""){
    ind <- TRUE
  } else {
    dataName <- mbSetObj$dataSet$name
    ind <- FALSE
  }
  
  module.type <- mbSetObj$module.type
  
  if(all(c(mbSetObj$module.type == "meta", !ind))){
    sums.list <- list()
    for(i in 1:length(mbSetObj$dataSets)){
      dataName <- mbSetObj$dataSets[[i]]$name
      data.proc <- readDataQs("data.proc", module.type, dataName)
      data_bef <- data.matrix(data.proc)
      smpl.sums <- colSums(data_bef)
      names(smpl.sums) <- colnames(data_bef)
      smpl.sums <- sort(smpl.sums)
      sums.list[[i]] <- smpl.sums
    }
    smpl.sums <- unlist(sums.list)
  } else {
    data.proc <- readDataQs("data.proc", module.type, dataName)
    data_bef <- data.matrix(data.proc)
    smpl.sums <- colSums(data_bef)
    names(smpl.sums) <- colnames(data_bef)
    smpl.sums <- sort(smpl.sums)
  }
  
  # Create a data frame for ggplot
  library_size_data <- data.frame(Sample = names(smpl.sums), LibrarySize = smpl.sums)
  library_size_data <- library_size_data %>% arrange(desc(LibrarySize))
  library_size_data$Sample = factor(library_size_data$Sample, levels=library_size_data[order(library_size_data$LibrarySize), "Sample"])
  
  # Plotting with ggplot2
  imgName <- paste0(origImgName, ".", format)
  g <- ggplot(library_size_data, aes(x = Sample, y = LibrarySize)) +
    geom_point(color = "forestgreen") +
    geom_hline(yintercept = mean(library_size_data$LibrarySize), color = "blue", linetype = "dashed") +
    labs(x = "Sample", y = "Read Counts") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  Cairo::Cairo(file=imgName, width=800, height=600, type=format, bg="white",dpi=dpi);
  print(g);
  dev.off();
  
  mbSetObj$imgSet$lib.size <- imgName;

  mean_line <- mean(library_size_data$LibrarySize)
  annotation_offset <- 0.05 * (max(library_size_data$LibrarySize) - min(library_size_data$LibrarySize))

  tooltips <- paste("Sample:", library_size_data$Sample, "<br>Library Size:", library_size_data$LibrarySize)

  # Convert data to JSON for Plotly
  plot_data <- list(
  data = list(
    list(
      x = seq_along(library_size_data$Sample),
      y = library_size_data$LibrarySize,
      type = 'scatter',
      mode = 'markers',
      marker = list(
        size = 10,
        line = list(color = 'white', width = 0.8)
      ),
      text = tooltips,
      hoverinfo = 'text'
    )
  ),
  layout = list(
    title = 'Library Size Overview',
    xaxis = list(title = 'Sample'),
    yaxis = list(title = 'Read Counts'),
    shapes = list(
      list(
        type = 'line',
        x0 = 0,  # Start at the left edge of the plot area
        y0 = mean_line,  # The y-value at which the line is placed
        x1 = 1,  # End at the right edge of the plot area
        y1 = mean_line,  # Same y-value to keep the line horizontal
        line = list(
          color = 'blue',  # Line color
          width = 3  # Line width
        ),
        xref = 'paper',  # Use the 'paper' reference for the x-axis
        yref = 'y'  # Use the y-axis data reference for the y-axis
      )
    ),
    annotations = list(
       list(
        x = 1,
        y = mean_line + annotation_offset,
        xref = 'paper',
        yref = 'y',
        text = paste('Mean:', signif(mean_line, 4)),
        showarrow = FALSE,
        font = list(family = 'Arial', size = 12, color = 'black')
      )
    )
  )
)
  
  jsonFileName <- paste0(origImgName, ".json")
  json.obj <- rjson::toJSON(plot_data );
  sink(jsonFileName);
  cat(json.obj);
  sink();

  if(interactive){
    library(plotly);
    tooltips <- paste("Sample: ", library_size_data$Sample, "\nSize:", library_size_data$LibrarySize)
    annotation_offset <- 50 # Adjust as needed

    fig <- plot_ly() %>%
      add_trace(data = library_size_data, x = seq_along(library_size_data$Sample), y = ~LibrarySize, type = 'scatter', mode = 'markers', text = ~tooltips, hoverinfo = 'text', marker = list(size = 10, line = list(color = 'white', width = 0.8)))

    fig <- fig %>% layout(
     title = plot_data$data$layout$title,
     xaxis = plot_data$layout$xaxis,
     yaxis = plot_data$layout$yaxis,
     shapes = plot_data$layout$shapes,
     annotations = plot_data$annotations
   )

   fig <- fig %>% config(displayModeBar = FALSE)

    return(fig);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}

#'Plot PCA plot for multi-omics samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'


PlotPCAView <- function(imgName, format="png", dpi=72,init){
  require("ggsci")
  
  #dpi<-as.numeric(dpi)
  imgName = paste(imgName,".", format, sep="");
  
  fig.list <- list()
  pca.list<- list()
  pct <- list();
  
  for(i in 1:2){
    if(i==1){
      data.proc <- current.proc$mic$data.proc;
      x <- data.matrix(data.proc);
      type="Microbiome"
    }else{
      data.proc <- current.proc$met$data.proc;
      x <- data.matrix(data.proc);
      type="Metabolomics"
    }
    pca <- prcomp(t(na.omit(x)), center=T, scale=T);
    imp.pca<-summary(pca)$importance;
    xlabel <- paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
    ylabel <- paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
    names <- colnames(x);
    pca.res <- as.data.frame(pca$x);
    pca.res <- pca.res[,c(1,2)]
    
    xlim <- GetExtendRange(pca.res$PC1);
    ylim <- GetExtendRange(pca.res$PC2);
    
    metaFile <- qs::qread("data.sample_data")
    Factor <- metaFile[,1]
    pca.rest <- pca.res
    pca.rest$Conditions <- Factor
    pca.rest$names <- rownames(pca.res)
    
    pcafig <- ggplot(pca.rest, aes(x=PC1, y=PC2,  color=Conditions)) +
      geom_point(size=3, alpha=0.7) + 
      scale_color_npg()+
      xlim(xlim) + 
      ylim(ylim) + 
      xlab(xlabel) + 
      ylab(ylabel) +
      theme_bw()
    fig.list[[i]] <- pcafig
    
    pos.xyz = data.frame(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3]);
    loading.pos.xyz = data.frame(pca$rotation[,c(1:3)]);
    loadingNames = rownames(pca$rotation);
    
    pos.xyz <- as.data.frame(pos.xyz);
    pos.xyz <- unitAutoScale(pos.xyz);
    
    loadingNames <- rownames(loading.pos.xyz);
    loading.pos.xyz <- as.data.frame(loading.pos.xyz);
    loading <- unitAutoScale(loading.pos.xyz);
    rownames(loading) <- loadingNames;
    nm <- paste0("pca_", type);
    pca.list[[nm]][["score"]] <- pos.xyz * 1000;
    pca.list[[nm]][["loading"]] <- loading* 1000;    
    
    pct[[nm]] <- unname(round(imp.pca[2,],3))[c(1:3)]*100;
  }
  pca.list$pct2 <- pct;
  qs::qsave(pca.list, file="pca.scatter.qs");
  
  h<-6*round(length(fig.list)/2)
  Cairo::Cairo(file=imgName, width=14, height=h, type=format, bg="white", unit="in", dpi=dpi);
  library("ggpubr")
  p1 <- ggarrange(plotlist=fig.list, ncol = 2, nrow = round(length(fig.list)/2), labels=c("Microbiome","Metabolomics"))
  print(p1)
  dev.off();
  return(1)
  
}


PlotDensityView <- function(imgName, format="png",  init=0, autoScale=T,dpi=72){

 # dpi <- as.numeric(dpi)
  imgName = paste(imgName,".", format, sep="");
  fig.list <- list()
  merged.df <- data.frame
   
  df.list <- list()
  for(i in 1:2){
    
    if(i==1){
      data.proc <- current.proc$mic$data.proc;
      dat <- data.matrix(data.proc);
      if(init==0){dat= t(AutoScale(t(data.proc))) }
      st <- stack(as.data.frame(dat))
      st$type <- "Microbiome"
      merged.df <- st
    }else{
      data.proc <- current.proc$met$data.proc;
      dat <- data.matrix(data.proc);
      if(init==0){dat = t(AutoScale(t(data.proc)))}
      st <- stack(as.data.frame(dat))
      st$type <- "Metabolomics"
      merged.df <-rbind(merged.df, st)
    } 
    
  }
  

 
  type<-merged.df$type
   merged.df$ind <- paste0(merged.df$ind, "_", merged.df$type)
  #if(init==0){
  # merged.df$values[which(merged.df$values==0)] <- 10^-3
  # merged.df$values <- log((merged.df$values),base=2)

# }

  Cairo::Cairo(file=imgName, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  g <- ggplot(merged.df, aes(x=values)) + 
    geom_line(aes(color=type, group=ind), stat="density", alpha=0.1) + 
    geom_line(aes(color=type), stat="density", alpha=0.7, linewidth=3) +
    scale_color_manual(values=c("#E64B35B2","#00A087B2"))+
    theme_bw()
  print(g)
  dev.off();
  

  return(1)
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
CreatePhyloseqObj<-function(mbSetObj, type, taxa_type, taxalabel,isNormInput){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type;
  load_phyloseq();
  load_splitstackshape();
  
  data.proc <- readDataQs("data.proc", module.type, dataName);
  
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
  if(mbSetObj$module.type == "sdp" | mbSetObj$micDataType =="ko"){
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
  if(any(c(mbSetObj$module.type == "mdp", mbSetObj$micDataType =="otu", mbSetObj$module.type == "ppd", mbSetObj$module.type == "meta"))){
    if(type=="text"){
      # prepare data for phyloseq visualization.
      # if features names are present in specific taxonomy format (greengenes or silva).
      if(taxalabel=="T"){ #taxa labels present and need to be parsed in OTU table
        feat_nm <- rownames(data.proc);
        mbSetObj$dataSet$feat_nm <- feat_nm;
        
        if(!(is.null(mbSetObj$dataSet$sglTax))){
          idxsgl <- match(toupper(mbSetObj$dataSet$sglTax),toupper(classi.lvl))
          feat_nm<-data.frame(mbSetObj$dataSet$feat_nm,check.names=FALSE);
          names(feat_nm)<-classi.lvl[idxsgl];
          #sanity check: no of sample of both abundance and metadata should match.
          taxmat<-as.matrix(feat_nm);
          rownames(taxmat)<-c(1:nrow(taxmat));
          
          taxa_table <- tax_table(taxmat);
          taxa_names(taxa_table)<-rownames(taxmat);
          
        }else{
          
       
          if(taxa_type=="SILVA"){
            
            load_splitstackshape();
            feat_nm<-data.frame(mbSetObj$dataSet$feat_nm,check.names=FALSE);
            names(feat_nm)<-"Rank";
            taxonomy<-splitstackshape::cSplit(feat_nm,"Rank",";");
            taxmat= data.frame(matrix(NA, ncol = 7, nrow = nrow(taxonomy)),check.names=FALSE);
            colnames(taxmat) <- classi.lvl;
            taxmat[,1:ncol(taxonomy)]<-taxonomy;
            taxmat<-taxmat[colSums(!is.na(taxmat))>0];
            taxmat<-as.matrix(taxmat);
            rownames(taxmat)<-c(1:nrow(taxmat));
            #phyloseq taxonomy object
            taxa_table <- tax_table(taxmat);
            taxa_names(taxa_table)<-rownames(taxmat);
          }else if(any(c(taxa_type=="Greengenes", taxa_type=="QIIME"))){
            my.parser<-lapply(feat_nm, parse_taxonomy_qiime); # note, this functions will remove empty strings, replace with NA!?
            taxa_table<-build_tax_table(my.parser);
            
            if(taxa_type=="Greengenes"){
              #additional columns are added(starting with Rank1 etc; has to be removed)
              indx<-grep("^Rank",colnames(taxa_table));
              # need to test if still columns left (could be that users specified wrong tax format)
              if(ncol(taxa_table) == length(indx)){
                AddErrMsg("Error parsing Greengenes Taxonomy Labels - you may try <b>Greengenes OTU ID</b> or <b>Not Specific /Other</b> in data upload.");
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
          }else if(any(c(taxa_type=="GreengenesID", taxa_type=="Others/Not_specific"))){
            
            # need to parse Taxonomy still!
            load_splitstackshape();               
            
            feat_nm <- data.frame(mbSetObj$dataSet$feat_nm,check.names=FALSE);
            feat_nm <- data.frame(apply(feat_nm, 1, function(x) gsub(";\\s;", ";", x)),check.names=FALSE) # remove empty taxa
            names(feat_nm) <- "Rank";
            taxonomy <- splitstackshape::cSplit(feat_nm,"Rank",";");
            taxmat <- data.frame(matrix(NA, ncol = 7, nrow = nrow(taxonomy)),check.names=FALSE);
            colnames(taxmat) <- classi.lvl;
            taxmat[,1:ncol(taxonomy)] <- taxonomy;
            taxmat <- taxmat[colSums(!is.na(taxmat))>0];
            taxmat <- as.matrix(taxmat);
            rownames(taxmat) <- c(1:nrow(taxmat));
            #phyloseq taxonomy object
            taxa_table <- tax_table(taxmat);
            taxa_names(taxa_table) <- rownames(taxmat);
          }
          
          
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
        blank.inx <- which(tmat=="")
        num.mat[na.inx] <- 0;
        num.mat[blank.inx] <- 0;
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
          mbSetObj$dataSet$master.nms = master.nms;
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
      if(isNormInput=="false"){
       data.proc<-apply(data.proc,2,as.integer);
      }   
      data.proc<-otu_table(data.proc,taxa_are_rows =TRUE);
      taxa_names(data.proc)<-taxa_names(mbSetObj$dataSet$taxa_table);
      # removing constant column and making standard names
      ind<-apply(mbSetObj$dataSet$taxa_table[,1], 2, function(x) which(x %in% c("Root", "root")));
      
      if(length(ind)>0){
        mbSetObj$dataSet$taxa_table<-mbSetObj$dataSet$taxa_table[,-1];
        #dataSet$taxa_table<-dataSet$taxa_table #newnew
      }
      
      # removing first KINGDOM OR DOMAIN column
      indx<-apply(mbSetObj$dataSet$taxa_table[,1], 2, function(x) which(x %in% c("k__Fungi", "Fungi", "k__Viruses", "Viruses","k_Bacteria",
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
      mbSetObj$dataSet$taxa_table[which(mbSetObj$dataSet$taxa_table==''|grepl("_$",mbSetObj$dataSet$taxa_table))]=NA;
       
if(taxalabel=="T"){# for mapping to the database in mmp module
  if("Species" %in% colnames( mbSetObj$dataSet$taxa_table)){
    sps = data.frame(mbSetObj$dataSet$taxa_table@.Data[,c('Genus','Species')])
    if(grepl("^g__",sps[1,"Genus"])){
      sps$Species = gsub("^s__","",sps$Species)
      if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
        sps$Genus =gsub("^g__","",sps$Genus)
        idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
        sps$Species[idxsp] = paste0("s__",sps$Genus[idxsp],"_",sps$Species[idxsp])
      }
    }else if(grepl("^g_",sps[1,"Genus"])){
      sps$Species = gsub("^s_","",sps$Species)
      if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
        sps$Genus =gsub("^g_","",sps$Genus)
        idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
        sps$Species[idxsp] = paste0("s__",sps$Genus[idxsp],"_",sps$Species[idxsp])
      }
    }else if(grepl("^g_0__",sps[1,"Genus"])){
      sps$Species = gsub("^s_0__","",sps$Species)
      if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
        sps$Genus =gsub("^g_0__","",sps$Genus)
        idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
        sps$Species[idxsp] = paste0("s_0__",sps$Genus[idxsp],"_",sps$Species[idxsp])
      }
  }else{
    if(all(!grepl(" ",sps$Species)) & all(!grepl("_",sps$Species))){
      idxsp = which(!is.na(sps$Genus) & !is.na(sps$Species) & sps$Genus!="" & sps$Species!="")
      sps$Species[idxsp] = paste0(sps$Genus[idxsp],"_",sps$Species[idxsp])
    }
  }
  mbSetObj$dataSet$taxa_table@.Data[,'Species']  <- sps$Species
}

}
#sometimes after removal of such special characters rownames beacame non unique; so make it unique
      mynames1<- gsub("[[:space:]./_-]", "_",taxa_names(mbSetObj$dataSet$taxa_table));
      mynames1<- gsub("\\[|\\]","",mynames1);
      
      if(length(unique(mynames1))< length(taxa_names(mbSetObj$dataSet$taxa_table))){
        mynames1 <- make.unique(mynames1, sep="_");
      }
      
      #after cleaning, names might now be NA
      missing.inx <- which(is.na(mynames1))
      missing_vector <- paste0("New_OTU", seq_along(missing.inx))
      mynames1[missing.inx] <- missing_vector
      
    } else if(any(c(type == "biom", type == "mothur"))){
      # creating phyloseq object
      feat_nm<-rownames(data.proc);
       if(isNormInput=="false"){
       data.proc<-apply(data.proc,2,as.integer);
      }
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
      indx<-apply(mbSetObj$dataSet$taxa_table[,1], 2, function(x) which(x %in% c("k__Fungi", "Fungi", "k__Viruses", "Viruses", "k_Bacteria",
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
      mbSetObj$dataSet$taxa_table[which(mbSetObj$dataSet$taxa_table==''|grepl("_$",mbSetObj$dataSet$taxa_table))]=NA;
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
    prefilt.data <- readDataQs("data.prefilt", module.type, dataName);
    rownames(prefilt.data)<-taxa_names(mbSetObj$dataSet$taxa_table);
   
    saveDataQs(prefilt.data, "data.prefilt", module.type, dataName);
    
  }else if(mbSetObj$module.type == "sdp"| mbSetObj$micDataType =="ko"){
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
  mbSetObj$dataSet$proc.phyobj@tax_table <- mbSetObj$dataSet$proc.phyobj@tax_table[,!is.na(colnames(mbSetObj$dataSet$proc.phyobj@tax_table))]
  saveDataQs(mbSetObj$dataSet$proc.phyobj, "proc.phyobj.orig", module.type, dataName);
  #do some data cleaning
  mbSetObj$dataSet$data.prefilt <- NULL;
  mbSetObj$dataSet$taxa_type <- taxa_type;
  
  if(module.type == "meta"){
    mbSetObj$dataSet$proc.phyobj <- mbSetObj$dataSet$proc.phyobj;
    mbSetObj$dataSet$norm.phyobj <-mbSetObj$dataSet$proc.phyobj;
  }

  return(.set.mbSetObj(mbSetObj));
}


###########################################################################
## Function to make the fake objects for normalized input ##
###########################################################################
### This function is used to create fake objects for
### normalized input to make sure all the function can work properly

CreateFakeFile <- function(mbSetObj,isNormalized="true",isNormalizedMet,module.type){

mbSetObj <- .get.mbSetObj(mbSetObj);

if(isNormalized=="false" & isNormalizedMet=="false"){

 AddErrMsg("Please make sure your data has been normalized properly!");
return(0)

}
mbSetObj$dataSet$filt.data <- mbSetObj$dataSet$data.orig
mbSetObj$dataSet$filt.msg <- "No data filtering has been performed since the input data has been transformed."

mbSetObj$dataSet$norm.phyobj <- mbSetObj$dataSet$proc.phyobj
mbSetObj$dataSet$norm.msg <- "No normalization has been performed since the input data has been transformed."

    #make hierarchies
  if(mbSetObj$module.type=="sdp" | mbSetObj$micDataType=="ko"){
    ranks <- "OTU"
  }else if(mbSetObj$module.type=="meta"){ #temporary
    ranks <- "OTU"
  }else{
    ranks <- c(GetMetaTaxaInfo(mbSetObj), "OTU")
    ranks <- unique(ranks)
  } 


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

if(exists("current.proc")){
  qs::qsave(data.list,"prescale.phyobj.qs")

current.proc$mic$data.proc<<- data.list$count_tables[["OTU"]]

}
 mbSetObj$dataSet$sample_data$sample_id <- rownames(mbSetObj$dataSet$sample_data);

saveDataQs(data.list, "phyloseq_prenorm_objs.qs",module.type, mbSetObj$dataSet$name);  
  saveDataQs(data.list, "phyloseq_objs.qs",module.type, mbSetObj$dataSet$name);  

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
##################Scale##################
###########################################
scale_rows = function(x){
  m = apply(x, 1, mean, na.rm = T)
  s = apply(x, 1, sd, na.rm = T)
  return((x - m) / s)
}
scale_mat = function(mat, scale){
  if(!(scale %in% c("none", "row", "column"))){
    stop("scale argument shoud take values: 'none', 'row' or 'column'")
  }
  mat = switch(scale, none = mat, row = scale_rows(mat), column = t(scale_rows(t(mat))))
  return(mat)
}

GetExtendRange<-function(vec, unit=10){
  var.max <- max(vec, na.rm=T);
  var.min <- min(vec, na.rm=T);
  exts <- (var.max - var.min)/unit;
  c(var.min-exts, var.max+exts);
}

unitAutoScale <- function(df){
  df <- as.data.frame(df)
  row.nms <- rownames(df);
  col.nms <- colnames(df);
  df<-apply(df, 2, AutoNorm);
  rownames(df) <- row.nms;
  colnames(df) <- col.nms;
  maxVal <- max(abs(df))
  df<- df/maxVal
  return(df)
}
AutoNorm<-function(x){
  (x - mean(x))/sd(x, na.rm=T);
}

#######################################################################
###Help function for processing metabolomics data##################
#######################################################################

# Limit of detection (1/5 of min for each var)
.replace.by.lod <- function(x){
  lod <- min(x[x>0], na.rm=T)/5;
  x[x==0|is.na(x)] <- lod;
  return(x);
}

ReplaceMissingByLoD <- function(int.mat=data.proc){
  int.mat <- as.matrix(int.mat);
  
  rowNms <- rownames(int.mat);
  colNms <- colnames(int.mat);
  int.mat <- apply(int.mat, 2, .replace.by.lod);
  rownames(int.mat) <- rowNms;
  colnames(int.mat) <- colNms;
  return (int.mat);
}


#'Perform data cleaning
#'@description Cleans data and removes -Inf, Inf, NA, negative and 0s.
#'@param bdata Input data to clean
#'@param removeNA Logical, T to remove NAs, F to not. 
#'@param removeNeg Logical, T to remove negative numbers, F to not. 
#'@param removeConst Logical, T to remove samples/features with 0s, F to not. 
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: MIT License
#'

CleanData <-function(bdata, removeNA=T, removeNeg=T, removeConst=T){
  
  if(sum(bdata==Inf, na.rm=TRUE)>0){
    inx <- bdata == Inf;
    bdata[inx] <- NA;
    bdata[inx] <- max(bdata, na.rm=T)*2
  }
  if(sum(bdata==-Inf, na.rm=TRUE)>0){
    inx <- bdata == -Inf;
    bdata[inx] <- NA;
    bdata[inx] <- min(bdata, na.rm=T)/2
  }
  if(removeNA){
    if(sum(is.na(bdata))>0){
      bdata[is.na(bdata)] <- min(bdata, na.rm=T)/2
    }
  }
  if(removeNeg){
    if(sum(as.numeric(bdata<=0)) > 0){
      inx <- bdata <= 0;
      bdata[inx] <- NA;
      bdata[inx] <- min(bdata, na.rm=T)/2
    }
  }
  if(removeConst){
    varCol <- apply(data.frame(bdata), 2, var, na.rm=T); # getting an error of dim(X) must have a positive length, fixed by data.frame
    constCol <- (varCol == 0 | is.na(varCol));
    constNum <- sum(constCol, na.rm=T);
    if(constNum > 0){
      bdata <- data.frame(bdata[,!constCol, drop=FALSE], check.names = F); # got an error of incorrect number of dimensions, added drop=FALSE to avoid vector conversion
    }
  }
  bdata;
}

LogNorm<-function(x, min.val){
  log10((x + sqrt(x^2 + min.val^2))/2)
}


# square root, tolerant to negative values
SquareRootNorm<-function(x, min.val){
  ((x + sqrt(x^2 + min.val^2))/2)^(1/2);
}

###########################################
##################Getters##################
###########################################

GetSamplNames<-function(dataName, module.type){
  prefilt.data <- readDataQs("data.prefilt", module.type, dataName);
  return(colnames(prefilt.data));
}

GetSampleNamesaftNorm<-function(mbSetObj, dataName){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$module.type == "meta"){
    if(dataName != ""){
      data <- readDataset(dataName);
    }else{
      data <- qs::qread("merged.data.raw.qs");
      data <- subsetPhyloseqByDataset(mbSetObj, data);
    }
    return(sample_names(data));
  }else{
    return(sample_names(mbSetObj$dataSet$norm.phyobj));
  }
}
