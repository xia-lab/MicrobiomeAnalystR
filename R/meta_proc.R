##################################################
## R scripts for MicrobiomeAnalyst
## Description: Data IO functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#' Init resources & global settings & variables.
#'
#' @return mbSetObj object, Input the name of the mbSetObj. By default,
# 'the name should be mbSet.
#' @export
#' @examples
Init.DataMeta <- function(mbSetObj, onWeb=T){ 
  isKo <<- F
  reductionOptGlobal <<- "pca"
  dataSet <- list(annotated=FALSE);
  dataSet <<- dataSet;
  
  paramSet <- list(annotated=FALSE);
  paramSet$on.public.web <- onWeb;
  .on.public.web <<- onWeb;

  dataSet <- list(annotated=FALSE);
  dataSet <<- dataSet;
  analSet <<- list(annotated=FALSE);
  paramSet <<- list(annotated=FALSE);
  msgSet <<- list(annotated=FALSE);
  cmdSet <<- list(annotated=FALSE);
  mdata.all <<- list();
  msg.vec <<- vector(mode="character");
  lib.path <<- "../../data/";
  data.org <<- NULL;
  if(onWeb){
  Set.Config("web");
  }else{
  Set.Config("api");
  }

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets <- list();
  mbSetObj$mdata.all <- list();

  saveSet(paramSet, "paramSet");
  saveSet(msgSet, "msgSet");
  saveSet(analSet, "analSet");
  saveSet(cmdSet, "cmdSet");

  return(.set.mbSetObj(mbSetObj));
}

# only for switching single expression data results
SetCurrentData <- function(nm){
  if(dataSet$name != nm){
    dataSet <<- readRDS(nm);
  }
  return(1);
}

# remove data object, the current dataSet will be the last one by default 
RemoveData <- function(dataName){
  mbSetObj <- .get.mbSetObj(NA);
  if(!is.null(mbSetObj$dataSets[[dataName]])){
    mbSetObj$dataSets[[dataName]] <- NULL;
    unlink(paste0(dataName, "_data"), recursive = TRUE)
  }
 if(mbSetObj$module.type=="meta"){
  if(!is.null(mdata.all[[dataName]])){
    mdata.all[[dataName]] <<- NULL;
  }
 }
  
  return(.set.mbSetObj(mbSetObj));
}

# users can select one or more data for analysis
# note, we use 1 to indicate this is selected
# and by default is all selected. 
SelectData <- function(){
  if(!exists('nm.vec')){
    current.msg <<-"No dataset is selected for analysis!";
    print(current.msg);
    return(0);
  }

  mbSetObj <- .get.mbSetObj(NA);
  mdata.all <- mbSetObj$mdata.all;

  all.nms <- names(mdata.all);
  for(nm in all.nms){
    if(nm %in% nm.vec){
      mdata.all[[nm]] <<- 1;
    }else{
      mdata.all[[nm]] <<- 0;
    }
  }
  
  rm('nm.vec', envir = .GlobalEnv);

  mbSetObj$mdata.all <- mdata.all;
  mdata.all <<- mdata.all;
  
  return(1);
}


GetAllDataNames <- function(){
  names(mdata.all);
}


PlotDataProfile<-function(mbSetObj=NA, dataName, type, libSizeName, pcaName){
  mbSetObj <- .get.mbSetObj(NA);
  dataSet <- readDataset(dataName);
  PlotLibSizeView(mbSetObj, libSizeName, "png", 72, dataName);
  qc.pcaplot2(as.matrix(otu_table(dataSet$proc.phyobj)), pcaName);
  return(.set.mbSetObj(mbSetObj));
}


qc.boxplot2 <- function(dat, imgNm){
  require('lattice');
  require('Cairo');
  imgNm = paste(imgNm, ".png", sep="");
  subgene=10000;
  if (nrow(dat)>subgene) {
    set.seed(28051968);
    sg  = sample(nrow(dat), subgene)
    Mss = dat[sg,,drop=FALSE]
  } else {
    Mss = dat
  }
  
  subsmpl=100;
  if (ncol(Mss)>subsmpl) {
    set.seed(28051968);
    ss  = sample(ncol(Mss), subsmpl)
    Mss = Mss[,ss,drop=FALSE]
  } else {
    Mss = Mss
  }
  
  sample_id = rep(seq_len(ncol(Mss)), each = nrow(Mss));
  values  = as.numeric(Mss)
  formula = sample_id ~ values
  
  box = bwplot(formula, groups = sample_id, layout = c(1,1), as.table = TRUE,
               strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
               horizontal = TRUE,
               pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
               xlab = "", ylab = "Samples",
               fill = "#1c61b6AA",
               panel = panel.superpose,
               scales = list(x=list(relation="free"), y=list(axs="i")),
               ylim = c(ncol(Mss)+0.7,0.3),
               prepanel = function(x, y) {
                 list(xlim = quantile(x, probs = c(0.01, 0.99), na.rm=TRUE))
               },
               panel.groups = function(x, y, ...) {
                 panel.bwplot(x, y, ...)
               })
  
  Cairo(file=imgNm, width=460, height=420, type="png", bg="white");
  print(box);
  dev.off();
}

qc.pcaplot2 <- function(x, imgNm){
  imgNm = paste(imgNm, ".png", sep="");
  require('lattice');
  pca <- prcomp(t(na.omit(x)));
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  
  # increase xlim ylim for text label
  GetExtendRange<-function(vec, unit=10){
    var.max <- max(vec);
    var.min <- min(vec);
    exts <- (var.max - var.min)/unit;
    return(c(var.min-exts, var.max+exts));
  }
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);
  pcafig = xyplot(PC2~PC1, data=pca.res, pch=19, cex=1,aspect = "iso", xlim = xlim, ylim=ylim,
                  panel=function(x, y, ...) {
                    panel.xyplot(x, y, ...);
                    ltext(x=x, y=y, labels=names, pos=1, offset=1, cex=0.8, col="magenta")
                  })
  
  Cairo(file=imgNm, width=480, height=480, type="png", bg="white");
  print(pcafig);
  dev.off();
}

###
### Data utilities
###

# read the uploaded data into memory
# return the meta-data information (multiple groups)
ReadDataForMetaInfo<-function(dataName){
    dataSet <- readDataset(dataName);
    return(colnames(dataSet$sample_data));
}


GetFeatureNum <-function(dataName){
  dataSet <- readDataset(dataName);
  if("norm.phyobj" %in% names(dataSet)){
    res <- nrow(dataSet$norm.phyobj@otu_table);
  }else{
    res <- nrow(dataSet$data.orig);

  }
  return(res);
}


PlotDataOverview<-function(imgNm){
    dat <- dataSet$data.anot;
    library('lattice');
    subgene=10000;
    if (nrow(dat)>subgene) {
        set.seed(28051968);
        sg  = sample(nrow(dat), subgene)
        Mss = dat[sg,,drop=FALSE]
    } else {
        Mss = dat
    }

    subsmpl=100;
    if (ncol(Mss)>subsmpl) {
        set.seed(28051968);
        ss  = sample(ncol(Mss), subsmpl)
        Mss = Mss[,ss,drop=FALSE]
    } else {
        Mss = Mss
    }

    sample_id = rep(seq_len(ncol(Mss)), each = nrow(Mss));
    values  = as.numeric(Mss)
    formula = sample_id ~ values

  Cairo(file=imgNm, width=460, height=420, type="png", bg="white");
    box = bwplot(formula, groups = sample_id, layout = c(1,1), as.table = TRUE,
        strip = function(..., bg) strip.default(..., bg ="#cce6ff"),
        horizontal = TRUE,
        pch = "|",  col = "black", do.out = FALSE, box.ratio = 2,
        xlab = "", ylab = "Samples",
        fill = "#1c61b6AA",
        panel = panel.superpose,
        scales = list(x=list(relation="free"), y=list(axs="i")),
        ylim = c(ncol(Mss)+0.7,0.3),
        prepanel = function(x, y) {
          list(xlim = quantile(x, probs = c(0.01, 0.99), na.rm=TRUE))
        },
        panel.groups = function(x, y, ...) {
          panel.bwplot(x, y, ...)
        })
  print(box);
  dev.off();
}

#'Function to read 16S Abundance data in meta-analysis module
#'@description This is the main function read in the 16S data
#'into the mbSetObj.
# ismetafile: whether meta-data is given (for BIOM format only)
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName File name of selected dataset
#'@param format File type, only supports ".txt"
#'@param taxa_type Taxonomy type (i.e.
#'@param ismetafile Whether meta-data is given (for BIOM format only)
#'@param is.normalized Specify whether uploaded data is normalized or raw counts
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
Read16SAbundDataMeta <- function(mbSetObj, dataName, format, taxa_type, ismetafile, is.normalized=F) {
    Read16SAbundData(mbSetObj, dataName, format, taxa_type, ismetafile, is.normalized=F);
    mbSetObj <- .get.mbSetObj(NA);
    .set.mbSetObj(mbSetObj);
    mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
    dir.create(paste0(dataName, "_data"), showWarnings = FALSE);
    return(.set.mbSetObj(mbSetObj));
}

#' SanityCheckSampleData
#'
#' @param mbSetObj microbiomeanalyst object, initialized by InitDataObjects("pktable", "mf", FALSE)
#' @param dataName File name of selected dataset
#' @export
SanityCheckSampleDataMeta <- function(mbSetObj=NA, dataName=""){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  SanityCheckSampleData(mbSetObj);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet; 

  return(.set.mbSetObj(mbSetObj)); 
}

#'Function to read in sample data in meta-analysis module
#'@description This functions reads in sample data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the file name of selected dataset
#'@param metaName Input the file name of meta-data file name.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
ReadSampleTableMetaInd <- function(mbSetObj, dataName, metaName) {
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  ReadSampleTable(mbSetObj, metaName);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;  
  return(.set.mbSetObj(mbSetObj));
}

#'Function to recreate phyloseq object in meta-analysis module
#'@description This function recreates the phyloseq object.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the file name of selected dataset
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
CreatePhyloseqObjMeta<-function(mbSetObj, dataName, type, taxa_type, taxalabel,isNormInput){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  CreatePhyloseqObj(mbSetObj, type, taxa_type, taxalabel,isNormInput);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
  return(.set.mbSetObj(mbSetObj));
}

#'Function to filter uploaded data in meta-analysis module
#'@description This function filters data based on low counts in high percentage samples.
#'Note, first is abundance, followed by variance.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the name of selected dataset data file
#'@param filt.opt Character, input the low count filter option. "prevalence" to
#'filter based on prevalence in samples, "mean" to filter by the mean abundance value, and
#'"median" to filter by median abundance value.
#'@param count Numeric, input the minimum count. Set to 0 to disable the low count filter.
#'@param smpl.perc Numeric, input the percentage of samples for which to filter low counts.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

ApplyAbundanceFilterMeta <- function(mbSetObj, dataName, filt.opt, count, smpl.perc){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  ApplyAbundanceFilter(mbSetObj, filt.opt, count, smpl.perc);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;  
  return(.set.mbSetObj(mbSetObj));
}

#'Function to filter uploaded data
#'@description This function filters data based on low abundace or variance.
#'Note, this is applied after abundance filter.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the name of selected dataset data file
#'@param filtopt Character, input the low variance filter option. "iqr" for 
#'inter-quantile range, "sd" for standard deviation, and "cov" for coefficient of variation.
#'@param filtPerct Numeric, input the percentage cutoff for low variance. 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
ApplyVarianceFilterMeta <- function(mbSetObj, dataName, filtopt, filtPerct){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  ApplyVarianceFilter(mbSetObj, filtopt, filtPerct);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;  
  return(.set.mbSetObj(mbSetObj));
}



#'Function to perform normalization
#'@description This function performs normalization on the uploaded
#'data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the name of selected dataset data file
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
PerformNormalizationMeta <- function(mbSetObj, dataName, rare.opt, scale.opt, transform.opt){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  PerformNormalization(mbSetObj, rare.opt, scale.opt, transform.opt);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;  
  return(.set.mbSetObj(mbSetObj));

}

#'Individual sanity check on selected datase in meta-analysis module
#'@description This function performs a sanity check on the uploaded
#'user data. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the name of selected dataset data file
#'@param filetype Character, "biom" if the uploaded data
#'was .biom format, "mothur" if mothur format, or "txt" if
#'as .txt or .csv format.
#'@param disableFilter Boolean. Set to TRUE to bypass the hard-coded filter
#'to remove features occuring in <2 samples.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SanityCheckDataMeta <- function(mbSetObj, dataName, filetype, disableFilter = FALSE){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  res <- SanityCheckData(mbSetObj, filetype, disableFilter);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;  
  .set.mbSetObj(mbSetObj);
  return(res);
}


#'Function to read 16S taxonomy table from txt format in meta-analysis module
#'@description This function reads in the 16S taxonomy table from txt format
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the name of selected dataset data file
#'@param fileName Input the name of selected dataset's taxa file
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
Read16STaxaTableMeta <- function(mbSetObj, dataName, fileName) {
    cat(fileName, dataName, "\n");
    mbSetObj <- .get.mbSetObj(mbSetObj);  
    mbSetObj$dataSet <- readDataset(dataName);
    .set.mbSetObj(mbSetObj);  
    Read16STaxaTable(mbSetObj, fileName);
    mbSetObj <- .get.mbSetObj(mbSetObj);
    .set.mbSetObj(mbSetObj);
    mbSetObj$dataSet$name <- dataName;
    mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
    return(.set.mbSetObj(mbSetObj));
}


#'Function to read in sample data in meta-analysis module (merged metadata file)
#'@description This functions reads in sample data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param fileName Input the name of merged metadata file
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
ReadSampleTableMeta <- function(mbSetObj, fileName) {
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_phyloseq();
  mydata <- .readDataTable(fileName);

  mydata[is.na(mydata)] <- "NA";
  if(any(c(any(is.na(mydata)), class(mydata) == "try-error"))){
    AddErrMsg("Failed to read in the metadata! Please make sure that the metadata file is in the right format and does not have empty cells or NA.");
    return(0);
  }
  
  # look for #NAME, store in a list
  #sam.nm <- substr(colnames(mydata[1]),1,5);
  #sam.nm <- tolower(sam.nm);
  #sam.inx <- grep("#name",sam.nm);
  #if(length(sam.inx) > 0){

   # first col is sample
    smpl_nm<-mydata[,1];
    smpl_var<-colnames(mydata[-1]);
  #}else{
  #  AddErrMsg("Please make sure you have the label #NAME in your sample data file!");
  #  return(0);
  #}
  
  # converting to character matrix as duplicate row names not allowed in data frame.
  mydata <- as.matrix(mydata[,-1]);
  
  if(nrow(mydata)==1){
    AddErrMsg("Only one sample in the dataset or the metadata file must be transposed!")
    return(0);
  }
  
  rownames(mydata) <- smpl_nm;
  colnames(mydata) <- smpl_var;
  
  # empty cell or NA cannot be tolerated in metadata
  na.inx  <- is.na(mydata);
  na.msg <-na.msg1 <- NULL;
  
  if(sum(na.inx) > 0){
    mydata[na.inx] <- "Unknown";
    na.msg1 <- paste("A total of", sum(na.inx), "empty or NA values were replaced by 'Unknown'.");
  }
  
  #Check group label names for spaces and replace with underscore
  my.meta <- data.frame(mydata,check.names=FALSE);
  my.meta.blank <- grep("[[:blank:]]", my.meta);
  if(length(my.meta.blank) >=  1){
    names(my.meta) <- gsub("\\s+","_", names(my.meta));
    na.msg1 <- c(na.msg1, "Blank spaces in group names are replaced with underscore '_'");
  }
  
  for(i in 1:length(mbSetObj$dataSets)){
  dataSet <- mbSetObj$dataSets[[i]];
  dataSet$group_names <- colnames(my.meta)
  #na.msg <- paste(na.msg, "The sample data contains a total of ", nrow(mydata), "samples and  ", ncol(mydata), " sample variables.", collapse=" ");
  
  # as most functions are for discrete groups (not continuous values)
  # require at least one column contains discrete factors with at least two replicates 
  disc.inx <- GetDiscreteInx(my.meta);
  if(sum(disc.inx) == 0){ # all class labels are unique! 
    na.msg <- c("<font style=\"color:red\"><b> No.</b></font>", "Cannot find a suitable variable as the primary metadata! MicrobiomeAnalyst requires some biological replicates for robust analysis!");
    mbSetObj$poor.replicate <- TRUE;
    dataSet$sample_data <- my.meta[which(rownames(my.meta) %in% rownames(dataSet$data.orig)),]
    dataSet$meta_info$disc.inx <- 0;
    dataSet$meta_info$cont.inx <- 0;
  }else{
    
    na.msg <- c(na.msg, "<font style=\"color:green\"><b> Yes. </b></font>");
    if(sum(disc.inx) == length(disc.inx)){
      na.msg <- c(na.msg,"All metadata columns are discrete.")
    }else{
      cont.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
      na.msg <- c(na.msg, paste0("Continuous variable(s): <b>", cont.meta, "</b>"));
    }
    
    dataSet$meta_info$disc.inx <- disc.inx;
    dataSet$sample_data <- my.meta[,disc.inx, drop=FALSE];
    
    cont.inx <- GetNumbericalInx(my.meta);
    cont.inx <- !disc.inx & cont.inx; # discrete is first
    dataSet$meta_info$cont.inx <- cont.inx;
    
    if(sum(cont.inx)>0){
      # make sure the discrete data is on the left side
      dataSet$sample_data <- cbind(dataSet$sample_data, my.meta[,cont.inx, drop=FALSE]);
    }
    
  }

  dataSet$smpl.msg <- c(na.msg,na.msg1);  

  mbSetObj$dataSets[[i]] <- dataSet; 
  .set.mbSetObj(mbSetObj);
  qs::qsave(dataSet$sample_data, paste0("sampleData",i, ".qs"));
  }

  mdata.all <- list();
  for(i in 1:length(mbSetObj$dataSets)){
    mdata.all[[names(mbSetObj$dataSets)[i]]] <- 1;
  }

  mbSetObj$mdata.all <- mdata.all;
  return(.set.mbSetObj(mbSetObj));
}

GetMetaInfoMeta <- function(mbSetObj, dataName, type="disc"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);
  return(GetMetaInfo(mbSetObj, type));
}


#'Check if data are ready for meta-analysis
#'@description This function determines if all annotated data are ready for meta-analysis
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

CheckMetaDataIntegrity <- function(mbSetObj, taxo_type="OTU", sample_var="NA"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  paramSet <- readSet(paramSet, "paramSet");
  paramSet$performedDE <- FALSE;
  mdata.all <- mbSetObj$mdata.all;
  for(i in 1:length(mdata.all)){
    mdata.all[i] <- 1;
  }
  
  sel.nms <- names(mbSetObj$dataSets);
  
  # first check that all class labels are consistent
  dataSet <- readDataset(sel.nms[1]);
  lvls <- levels(dataSet$proc.phyobj@sam_data[,1][[1]]);
  id.type <- ""; #dataSet$id.type;
  nms <- rownames(dataSet$proc.phyobj@otu_table);
  shared.nms <- nms;
  shared.meta.nms <- colnames(dataSet$proc.phyobj@sam_data);
  if(length(sel.nms)>1){
    for(i in 2:length(sel.nms)){
      dataSet <- readDataset(sel.nms[i]);
      data <- dataSet$proc.phyobj@otu_table;
      # check if class label is consistent
      if(!all(levels(dataSet$sample_data) == lvls)){
        msgSet$current.msg <- paste(sel.nms[i], "has different group labels", 
                                    paste(levels(dataSet$sample_data), collapse=":"), 
                                    "from", sel.nms[1], paste(lvls, collapse=":"));
        saveSet(msgSet, "msgSet");    
        print(msgSet$current.msg);                  
        return(0);
      }
      
      # check and record if there is common genes            
      shared.nms <- intersect(shared.nms, rownames(data));
      if(length(shared.nms) < 10){
        msgSet$current.msg <- paste(sel.nms[i], "has less than 10 common genes/probes from previous data sets");
        saveSet(msgSet, "msgSet");                      
        print(msgSet$current.msg);
        return(0);
      }
      
      shared.meta.nms <- intersect(shared.meta.nms, colnames(dataSet$proc.phyobj@sam_data));
      if(length(shared.meta.nms) < 1){
        msgSet$current.msg <- paste(sel.nms[i], "has no shared metadata with previous dataset.");
        saveSet(msgSet, "msgSet");                      
        print(msgSet$current.msg);
        return(0);
      }
    }   
  }
  
  print("Passed exp condition check!");
  
  # now construct a common matrix to faciliated plotting across all studies
  dataName <- sel.nms[1];
  dataSet <- readDataset(dataName);
  data <- dataSet$norm.phyobj@otu_table;
  common.matrix <- data[rownames(data) %in% shared.nms, ];
  nms.vec <- rownames(data);
  smps.vec <- colnames(data);
  data.lbl <- rep(dataName, ncol(common.matrix));
  sam_data <- dataSet$norm.phyobj@sam_data;
  dataSet$norm.phyobj@sam_data <- sam_data[,shared.meta.nms]
  dataSet$proc.phyobj@sam_data <- sam_data[,shared.meta.nms]
  mbSetObj$dataSets[[dataSet$name]] <- dataSet;
  if(length(sel.nms)>1){
    for(i in 2:length(sel.nms)){
      dataName <- sel.nms[i];
      dataSet <- readDataset(dataName);
      data <- dataSet$norm.phyobj@otu_table;
      
      ndat <- data[rownames(data) %in% shared.nms, ];
      nms.vec <- c(nms.vec, rownames(data));
      smps.vec <- c(smps.vec, colnames(data));
      plot.ndat <- t(scale(t(ndat)));
      common.matrix <- cbind(common.matrix, ndat);
      data.lbl <- c(data.lbl, rep(dataName, ncol(data)));
      dataSet$norm.phyobj@sam_data <- dataSet$norm.phyobj@sam_data[,shared.meta.nms];
      dataSet$proc.phyobj@sam_data <- dataSet$proc.phyobj@sam_data[,shared.meta.nms];
      sam_data <- rbind(sam_data[,shared.meta.nms], dataSet$norm.phyobj@sam_data[,shared.meta.nms]);
      mbSetObj$dataSets[[dataSet$name]] <- dataSet;
    }
  }
  .set.mbSetObj(mbSetObj);
  #levels(cls.lbl) <- lvls;
  
  smps.nms <- colnames(common.matrix)
  if(length(unique(smps.nms)) != length(smps.nms)){
    data.nb <- length(unique(data.lbl));
    data.vec <- vector()
    for(i in 1:data.nb){
      data.vec[i] <- paste0("d", i);
    }
    levels(data.lbl) <- data.vec;
    colnames(common.matrix) <- make.unique(paste(data.vec, smps.nms, sep="_"));
    
    dataSet <- readDataset(sel.nms[1]);
    colnames(data) <- paste("d1", colnames(data), sep="_");
    dataSet$norm.phyobj@otu_table <- data;
    mbSetObj$dataSets[[dataSet$name]] <- dataSet;
    
    for(i in 2:length(sel.nms)){
      dataSet <- readDataset(sel.nms[i]);
      colnames(data) <- paste0("d",i,"_",colnames(data));
      # check if class label is consistent
      dataSet$norm.phyobj@otu_table <- data;
      mbSetObj$dataSets[[dataSet$name]] <- dataSet;
    }
    smps.vec <- smps.nms;
    msgSet$current.msg <- paste("Duplicated sample names detected, samples have been renamed to make them unique.");
  }
  
  # note: index by entrez, gene symbol DO have duplicates
  rownames(common.matrix) <- shared.nms;
  
  
  if(sample_var == "null"){
    sample_var = colnames(sam_data)[1];
  }
  
  if(length(sel.nms) > 1){
    data <- mbSetObj$dataSets[[1]];
    data2 <- mbSetObj$dataSets[[2]];
    
    PerformDataMerging(mbSetObj, data, data2, taxo_type, sample_var, T);
    if(length(sel.nms)>2){
      for(i in 3:length(sel.nms)){
        dataX <- mbSetObj$dataSets[[i]];
        merged.data <- qs::qread("merged.data.raw.qs");
        orig.proc <-  mbSetObj$dataSet$proc.phyobj ;
        mbSetObj$dataSet$proc.phyobj <- merged.data;
        PerformDataMerging(mbSetObj, mbSetObj$dataSet, dataX, taxo_type, sample_var, F);
        mbSetObj$dataSet$proc.phyobj <- orig.proc;
      }
    }
    merged.data <- qs::qread("merged.data.raw.qs");
  }else{
    samdata <- as.data.frame(as.matrix(sample_data(mbSetObj$dataSet$proc.phyobj)));
    samdata$dataset <- rep(mbSetObj$dataSet$name,nrow(sample_data(mbSetObj$dataSet$proc.phyobj)));
    sample_data(mbSetObj$dataSet$proc.phyobj) <- samdata;
    merged.data <- mbSetObj$dataSet$proc.phyobj;
    qs::qsave(merged.data, "merged.data.raw.qs");

  }
  .set.mbSetObj(mbSetObj);
  
  ref.nms <- mbSetObj$ref.nms;
  if(length(ref.nms) > 0){
    for(i in 1:length(ref.nms)){
      orig.proc <-  mbSetObj$dataSet$proc.phyobj ;
      res <- PerformRefDataMapping(mbSetObj, ref.nms[1], taxo_type, sample_var, "");
      if(res == 0){
        return(0);
      }
      mbSetObj$dataSet$proc.phyobj <- qs::qread("merged.data.raw.qs");
      .set.mbSetObj(mbSetObj)
    }
  }
  #merged.data <- qs::qread("merged.data.raw.qs");
  
  meta.keep.inx <- sapply(sample_data(merged.data), function(col) length(unique(col))) > 1;
  merged.data@sam_data<- sample_data(merged.data)[, meta.keep.inx];
  sam_data <- as.data.frame(as.matrix(merged.data@sam_data));
  rownames(sam_data) <- rownames(sample_data(merged.data));
  cls.lbl <- sam_data[,1];
  
  
  #qs::qsave(merged.data, "merged.data.raw.qs");
  # if entrez, get symbols for display
  shared.nms <- rownames(common.matrix);
  symbols <- shared.nms;
  msg <- vector();
  names(symbols) <- shared.nms;
  meta_info <- list();
  disc.inx <- GetDiscreteInx(sam_data);
  if(sum(disc.inx) == 0){ # all class labels are unique! 
    msg <- c(msg, "It seems that all your meta data values are unique! MicrobiomeAnalyst requires some biological replicates for robust analysis");
    mbSetObj$poor.replicate <- TRUE;
  }else{
    if(sum(disc.inx) == length(disc.inx)){
      msg <- c(msg, "All metadata columns are OK.");
    }else{
      bad.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
      msg <- c(msg, paste0("The following metadata columns are excluded: <font style=\"color:red\"><b>", bad.meta, "</b></font>"));
    }
    meta_info$disc.inx <- disc.inx;
    sam_data <- sam_data[,disc.inx, drop=FALSE];
    
    cont.inx <- GetNumbericalInx(sam_data);
    cont.inx <- !disc.inx & cont.inx; # discrete is first
    meta_info$cont.inx <- cont.inx;
    
    if(sum(cont.inx)>0){
      # make sure the discrete data is on the left side
      sam_data <- cbind(sam_data, sam_data[,cont.inx, drop=FALSE]);
    }
  }
  
  # resort data, first on data.lbl, then on class lbl
  if("dataset" %in% colnames(sample_data(merged.data))){
    data.lbl <- sample_data(merged.data)$dataset;
  }else{
    data.lbl <- rep(mbSetObj$dataSet$name, nrow(sample_data(merged.data)));
  }
    ord.inx <- order(data.lbl, sam_data[,1]);
    common.matrix <- data.matrix(common.matrix[,ord.inx]);
    sam_data <- sam_data[ord.inx, ];
    data.lbl <- data.lbl[ord.inx];
    smps.vec <- smps.vec[ord.inx];
    nms.vec <- unique(nms.vec)

  if(ncol(common.matrix) > 1000){  # max sample number allow 1000
    msgSet$current.msg <- paste("Total combined sample #:", ncol(common.matrix), "(exceed the limit: 1000!)");
    saveSet(msgSet, "msgSet");
    return(0);
  }
  
  # save the meta-dataset
  res <- data.frame(colnames(common.matrix), cls.lbl, data.lbl, t(common.matrix));
  colnames(res) <- c('#NAME', '#CLASS.condition', paste('#CLASS.dataset',id.type, paramSet$data.org, sep="."), rownames(common.matrix));
  write.table(t(res), sep="\t", file="MicrobiomeAnalyst_merged_data.txt", col.names=F, quote=FALSE);
  
  # need to set up the data for plotting (boxplot, heatmpa) so 
  # we need to scale row for each dataset in order to elimiate the maganitude difference 
  plot.matrix <- matrix(NA, nrow=nrow(common.matrix), ncol=ncol(common.matrix));
  rownames(plot.matrix) <- rownames(common.matrix);
  colnames(plot.matrix) <- colnames(common.matrix);
  #for(i in 1:length(sel.nms)){
  #  data.inx <- data.lbl == sel.nms[i];
  #  plot.matrix[,data.inx] <- t(scale(t(common.matrix[,data.inx])));
  #}

  if(class(sam_data) != "data.frame"){
    sam_data <- as.data.frame(sam_data);
    colnames(sam_data) <- "CLASS";
  }
  meta.types <- rep("disc", ncol(sam_data));
  meta.types[meta_info$cont.inx] <- "cont";
  names(meta.types) <- colnames(sam_data);
  meta.types<-meta.types[!is.na(meta.types)];  

  microbiome.meta <- list(data=common.matrix,
                          plot.data=plot.matrix,
                          id.type = id.type,
                          gene.symbls = symbols,
                          cls.lbl=factor(cls.lbl),
                          meta.data=sam_data,
                          meta.info=meta_info,
                          meta.types=meta.types,
                          data.lbl=data.lbl);
  
  qs::qsave(microbiome.meta, "microbiome_meta.qs");
  paramSet$smps.vec <- colnames(common.matrix);
  
  # setup common stats gene number, smpl number, grp info
  msgSet$current.msg <- paste("Sample #:", ncol(microbiome.meta$data),
                              "Common ID #:", nrow(microbiome.meta$data), 
                              "Condition:", paste(levels(microbiome.meta$cls.lbl), collapse=" vs. "));
  
  saveSet(msgSet, "msgSet");
  saveSet(paramSet, "paramSet");
  mbSetObj$mdata.all <- rep(1, length(mbSetObj$dataSets));
  names(mbSetObj$mdata.all ) <- names(mbSetObj$dataSets);
  return(.set.mbSetObj(mbSetObj));
}

#'Merge two or more datasets into merged phyloseq object
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
MergeDatasets <- function(mbSetObj, taxo_type, sample_var){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(sample_var == "null"){
    #sample_var = colnames(sam_data)[1];
  }
  sel.nms <- names(mbSetObj$dataSets);
  if(length(sel.nms) > 1){
    data <- mbSetObj$dataSets[[1]];
    data2 <- mbSetObj$dataSets[[2]];
    PerformDataMerging(mbSetObj, data, data2, taxo_type, sample_var, T, "norm");
    PerformDataMerging(mbSetObj, data, data2, taxo_type, sample_var, T, "raw");

    if(length(sel.nms)>2){
      for(i in 3:length(sel.nms)){
        dataX <- mbSetObj$dataSets[[i]];
        merged.data <- qs::qread("merged.data.norm.qs");
        merged.data.raw <- qs::qread("merged.data.raw.qs");

        orig.norm <-  mbSetObj$dataSet$norm.phyobj;
        orig.proc <-  mbSetObj$dataSet$proc.phyobj;

        mbSetObj$dataSet$norm.phyobj <- merged.data;
        mbSetObj$dataSet$proc.phyobj <- merged.data.raw;

        PerformDataMerging(mbSetObj, mbSetObj$dataSet, dataX, taxo_type, sample_var, F, "norm");
        PerformDataMerging(mbSetObj, mbSetObj$dataSet, dataX, taxo_type, sample_var, F, "raw");

        mbSetObj$dataSet$norm.phyobj <- orig.norm;
        mbSetObj$dataSet$proc.phyobj <- orig.proc;
      }
    }

    phyobj <- qs::qread("merged.data.norm.qs");
    microbiome.meta <- qs::qread("microbiome_meta.qs");
    microbiome.meta$data <- as.data.frame(otu_table(phyobj));
    qs::qsave(microbiome.meta, "microbiome_meta.qs");

    .set.mbSetObj(mbSetObj);
  }else{
    mbSetObj$dataSet <- mbSetObj$dataSets[[1]];
    samdata <- as.data.frame(as.matrix(sample_data(mbSetObj$dataSet$proc.phyobj)));
    samdata$dataset <- rep(mbSetObj$dataSet$name,nrow(sample_data(mbSetObj$dataSet$proc.phyobj)));
    sample_data(mbSetObj$dataSet$proc.phyobj) <- samdata;
    merged.data <- mbSetObj$dataSet$proc.phyobj;
    qs::qsave(merged.data, "merged.data.raw.qs");
    .set.mbSetObj(mbSetObj);
  }
  
  ref.nms <- mbSetObj$ref.nms;

  if(length(ref.nms) > 0){
    for(i in 1:length(ref.nms)){
      res <- PerformRefDataMapping(mbSetObj, ref.nms[1], taxo_type, sample_var, "");
      if(res == 0){
        return(0);
      }
      mbSetObj$dataSet$proc.phyobj <- qs::qread("merged.data.raw.qs");
      .set.mbSetObj(mbSetObj)
    }
  }
  merged.data <- qs::qread("merged.data.raw.qs");

  sam_data <- as.data.frame(as.matrix(merged.data@sam_data));
  rownames(sam_data) <- rownames(sample_data(merged.data));
  cls.lbl <- sam_data[,1];
  qs::qsave(merged.data, "merged.data.raw.qs");
  norm.data <- transform_sample_counts(merged.data, function(x) x / sum(x) );
  qs::qsave(norm.data, "merged.data.qs");


  return(1);
}

PerformDataMerging <- function(mbSetObj, data1, data2, taxo_type, sample_var, init=T, type="proc"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  sel.nms <- names(mbSetObj$dataSets);
  load_phyloseq();
  #reading user data
  if(type == "proc"){
  data <- data1$proc.phyobj;
  current.refset <- data2$proc.phyobj;
current.sample <- current.refset@sam_data
  }else if(type == "filt"){
  data <- data1$filt.data;
  current.refset <- data2$filt.data;
current.sample <- data2$proc.phyobj@sam_data
  }else{
  data <- data1$norm.phyobj;
  current.refset <- data2$norm.phyobj;
current.sample <- current.refset@sam_data
  }

  #sample data
  otu_no <- ntaxa(data);
  
  #since we use modified name for each OTU(to make it unique and convenient);in order to do mapping we need complete and original mapping label.
  #since our reference data has Greengenes OTU Ids as onomy identifier, so if data is already in this format, we will directly merge both of them(reference and users data;kepping all OTU's unique to user data)
  #if taxo_type is SILVA; first we have to map it with Greengenes OTU Ids from mapping file;#then we have to compare Greengenes OTU Ids OTUs between user and reference data
  if(taxo_type=="SILVA"){
    taxa_names(data) <- mbSetObj$dataSet$comp_taxnm;
    a <- taxa_names(data);
    #taxonomy mapping file
    
    if(.on.public.web){
      otu.dic <<- qs::qread("../../lib/picrust/greengenes_taxmap.qs");
    }else{
      otu.dic <<- readRDS("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/picrust/greengenes_taxmap.rds");
    }
    
    #returns a vector of the positions of (first) matches of its first argument in user data(second argument).
    match_ind <- match(unique(otu.dic[ ,taxo_type]),a);
    #getting all the indices of user data that match(removing NA)
    match_ind1 <- match_ind[!is.na(match_ind)];
    
    if(length(match_ind1)==0){
      AddErrMsg("No SILVA labels match between user and taxonomy mapping file. Please check taxonomy labels ");
      return(0);
    }
    
    #getting only matched taxa in user data and then have to prune other taxa;
    match_taxa <- a[match_ind1];
    #now replacing the taxonomy names with Greengenes OTU Idss
    match_taxagg <- otu.dic[match(match_taxa,otu.dic[ ,taxo_type]),1];
    match_taxagg <- as.character(match_taxagg);
    
    #getting only matched taxa in user data and then have to prune other taxa;followed by replacing it with Greengenes OTU Ids;
    data = prune_taxa(match_taxa,data);
    taxa_names(data) <- match_taxagg;
  }
  
  # have to check whether 20% of  OTUs common in user and reference data or not;
  taxa_ind <- match(taxa_names(data),taxa_names(current.refset));
  
  if(length(which(taxa_ind != "NA")) < 0.20*otu_no){
    AddErrMsg(paste(c("Less than 20 percent OTU match between user and reference data.", length(which(taxa_ind != "NA")), "% match!"), collapse=" "));
    return(0);
  } 
  
  #pruning reference dataset
  current.refset = prune_taxa(taxa_names(data), current.refset);
  
  if(taxo_type=="SILVA"){
    #for SILVA id we can compare only between common id that matches between user and reference data;pruning others OTUs in reference data.
    data <- prune_taxa(taxa_names(current.refset),data);
  }
  
  #merging both the user and reference phyloslim objects(sample data will be merged seperately)
  #storing taxonomic rank for users data
  #all reference data have kingdom column which can be removed;
  
  userdatarank <<- rank_names(current.refset);
  current.ref_userdata <- merge_phyloseq(otu_table(data),otu_table(current.refset),tax_table(data),tax_table(current.refset));
  #dummy variable for showing different shape for user and reference data
  
  current.sample$dataset <- rep(data2$name,nrow(current.sample));
  
  sam_data <- as.data.frame(sample_data(data),check.names=FALSE);
  if(init){
    sam_data$dataset <- rep(data1$name,nrow(sample_data(data)));
  }
  shared.nms <- intersect(colnames(current.sample), colnames(sam_data));
  current.sample <- current.sample[, shared.nms];
  sam_data <- sam_data[, shared.nms];
  #selecting primary variable data
  #user_sam<-sam_data[sample_var];
  #making the name of same variable same for reference sample data and then merging them
  current.ref_usersamdata <- rbind(current.sample,sam_data);
  current.ref_usersamdata$dataset <- as.factor(current.ref_usersamdata$dataset);    
  current.ref_usersamdata <- sample_data(current.ref_usersamdata);
  #storing the taxonomy rank from users data
  #final phyloslim object;(taxonomy table after merging will get distorted if taxa_ranks are not same in both user and reference data;but it's of no use)
  merged.data <- merge_phyloseq(current.ref_userdata, current.ref_usersamdata);
  
  #change column name to dataset due to issue in plotting later

  if(type == "proc"){
  qs::qsave(merged.data, "merged.data.raw.qs");
  }else if(type == "filt"){
  qs::qsave(merged.data, "merged.data.filt.qs");
  }else{
  qs::qsave(merged.data, "merged.data.norm.qs");
  }
  
  #data filteration and transformation
  merged.data <- transform_sample_counts(merged.data, function(x) x / sum(x) );
  
  qs::qsave(merged.data, "merged.data.qs");
  #current.msg <<- paste(msg, collapse=".");
  mbSetObj$dataSet$lib.msg <- current.msg;
  
  return(.set.mbSetObj(mbSetObj));
  
}

SelectRefData <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$ref.nms <- ref.nms.vec;
  return(.set.mbSetObj(mbSetObj));
}


.checkSampleConsistencyMeta <- function(mbSetObj, dataName){
    # first check that all class labels are consistent
    sel.nms <- names(mbSetObj$dataSets);
    if(length(sel.nms) > 1 ){
      dataSet <- readDataset(dataName);
      lvls <- levels(dataSet$sample_data[,1]);
      id.type <- ""; #dataSet$id.type;
      nms <- rownames(dataSet$data.orig);
      shared.nms <- nms;
      shared.meta.nms <- colnames(dataSet$sample_data); 
      subset.sel.nms <- sel.nms[sel.nms != dataName]
      is.norm <- dataSet$is.normalized;
        for(i in 1:length(subset.sel.nms)){
          dataSet <- readDataset(subset.sel.nms[i]);
          data <- dataSet$data.orig;
          if(dataSet$is.normalized != is.norm){
            AddErrMsg("Please make sure uploaded datasets are consistently raw counts or normalized.");
            return(0);
          }
          # check if class label is consistent
          if(!all(levels(dataSet$sample_data) == lvls)){
            AddErrMsg(paste(subset.sel.nms[i], "has different group labels", 
                                        paste(levels(dataSet$sample_data), collapse=":"), 
                                        "from", subset.sel.nms[1], paste(lvls, collapse=":")));
            return(0);
          }
          
          # check and record if there is common genes      

          shared.nms <- intersect(shared.nms, rownames(data));
          if(length(shared.nms) < 10){
            AddErrMsg(paste(subset.sel.nms[i], "has less than 10 common genes/probes from previous data sets"));
            return(0);
          }
          
          shared.meta.nms <- intersect(shared.meta.nms, colnames(dataSet$sample_data));
          if(length(shared.meta.nms) < 1){
            AddErrMsg(paste(subset.sel.nms[i], "has no shared metadata with previous dataset."));
            return(0);
          }

      }
    }
    return(1)
}

PlotLibSizeHistogram <- function(mbSetObj, imgName,format="png", dpi=72, dataName=""){
  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(dataName != ""){
    # plot single dataset or not, in metaanal
    ind <- T; 
  }else{
    dataName <- mbSetObj$dataSet$name;
    ind <- F;
  }
  module.type <- mbSetObj$module.type;
  
  if(all(c(mbSetObj$module.type == "meta", !ind))){
    sums.list <- list();
    for(i in 1:length(mbSetObj$dataSets)){
      dataName <- mbSetObj$dataSets[[i]]$name;
      data.proc <- readDataQs("data.proc", module.type, dataName);
      data_bef <- data.matrix(data.proc);
      smpl.sums <- colSums(data_bef);
      names(smpl.sums) <- colnames(data_bef);
      smpl.sums <- sort(smpl.sums);
      sums.list[[i]] <- smpl.sums
    }
    smpl.sums  <- unlist(sums.list)
  }else{
    data.proc <- readDataQs("data.proc", module.type, dataName);
    data_bef <- data.matrix(data.proc);
    smpl.sums <- colSums(data_bef);
    names(smpl.sums) <- colnames(data_bef);
    smpl.sums <- sort(smpl.sums);
  }
  
  # save the full lib size 
  fast.write(cbind(Size=smpl.sums), file="norm_libsizes.csv");
  
  smpl.sums <- sort(smpl.sums);
  vip.nms <- names(smpl.sums);
  vip.nms <- substr(vip.nms, 1, 16);
  
  myH <- ncol(data_bef)*25 + 50;
  imgName = paste(imgName,".", format, sep="");
  mbSetObj$imgSet$lib.size<-imgName;
  
  df <- data.frame(sampleName=vip.nms, counts=smpl.sums, index=seq(1,length(smpl.sums)))
  p1 <- ggplot(df) +
    geom_histogram(aes(x = counts), color = "black", fill = "gray", bins = 30) +
    labs(x = "Library size", y = "Frequency (n)") + 
    # scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), 
    # labels = scales::trans_format("log10", scales::math_format(10^.x))) +
    theme_bw() +
    theme(panel.grid.major = element_blank(), # Removes the grid
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "black")) # Adds y-axis
  
  Cairo::Cairo(file=imgName, width=580, height=450, type=format, bg="white",dpi=dpi);
  print(p1)
  dev.off();

  return(.set.mbSetObj(mbSetObj))
}

subsetPhyloseqByDataset <- function(mbSetObj, phyobj){
   mbSetObj <- .get.mbSetObj(mbSetObj);
   mdata.all <- mbSetObj$mdata.all;
   sel.nms <- names(mdata.all)[mdata.all==1];

  if (is.null(sample_data(phyobj))) {
    cat("Nothing subset. No sample_data in phyobj.\n")
    return(phyobj)
  }
  else {
    oldDF <- as(sample_data(phyobj), "data.frame")
    newDF <- oldDF[oldDF$dataset %in% sel.nms, ]
    if (class(phyobj) == "sample_data") {
      return(sample_data(newDF))
    }
    else {
      sample_data(phyobj) <- sample_data(newDF)
      return(phyobj)
    }
  }
}


##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################

GetOrigSmplNmsInd <-function(mbSetObj=NA, dataName=""){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);  

  return(rownames(mbSetObj$dataSet$sample_data));
}

GetOrigSmplGroupNamesPerMetaInd <- function(mbSetObj=NA,dataName="", metaname="NA"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);  

  return(GetOrigSmplGroupNamesPerMeta(mbSetObj, metaname));
}

UpdateSampleGroupsInd<-function(mbSetObj=NA,dataName, metadata="NA"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);  
  UpdateSampleGroups(mbSetObj, metadata);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  .set.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
  .set.mbSetObj(mbSetObj);
  return(SanityCheckDataMeta(mbSetObj, dataName, "txt", FALSE));
}

RemoveSelectedMetaInd <- function(mbSetObj=NA, dataName, meta){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);  
  RemoveSelectedMeta(mbSetObj, meta);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  .set.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
  return(.set.mbSetObj(mbSetObj));
}

UpdateMetaColNameInd <- function(mbSetObj=NA, dataName, oldName, newName){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);  
  UpdateMetaColName(mbSetObj, oldName, newName);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  .set.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
  return(.set.mbSetObj(mbSetObj));

}

UpdateMetaLevelsInd <- function(mbSetObj=NA,dataName, metaNm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);  
  UpdateMetaLevels(mbSetObj, metaNm);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  .set.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
  return(.set.mbSetObj(mbSetObj));
}

UpdateMetaOrderInd <- function(mbSetObj=NA,dataName, metaName){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet <- readDataset(dataName);
  .set.mbSetObj(mbSetObj);  
  UpdateMetaOrder(mbSetObj, metaName);
  mbSetObj <- .get.mbSetObj(mbSetObj);
  .set.mbSetObj(mbSetObj);
  mbSetObj$dataSets[[dataName]] <- mbSetObj$dataSet;
  return(.set.mbSetObj(mbSetObj));
}
