##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# This is only for web version
#.on.public.web <- TRUE; # only TRUE when on MicrobiomeAnalyst web server

# note, this is usually used at the end of a function
# for local, return itself; for web, push to global environment
.set.mbSetObj <- function(mbSetObj=NA){
  if(.on.public.web){
    mbSet <<- mbSetObj;
    return (1);
  }
  return(mbSetObj);
}

.get.mbSetObj <- function(mbSetObj=NA){
  if(.on.public.web){
    return(mbSet);
  }else{
    return(mbSetObj);
  }
}

#'Constructs a mbSet object for storing data 
#'@description This functions handles the construction of a mbSetObj object for storing data 
#'produced by MicrobiomeAnalystR for further processing and analysis.
#'@usage Init.mbSetObj()
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
Init.mbSetObj <- function(){

  rm(list = ls(all.names = TRUE))
  dataSet <- list();
  analSet <- list();
  imgSet <- list();
   if(!exists(".on.public.web")){
    .on.public.web <<- FALSE
  }  

  mbSetObj <- list();
  mbSetObj$dataSet <- dataSet;
  mbSetObj$analSet <- analSet;
  mbSetObj$imgSet <- imgSet;
  mbSetObj$is.ASV <- FALSE;
  mbSetObj$poor.replicate <- FALSE; # flag to show if all unique values
  mbSetObj$tree.uploaded <- FALSE;
  mbSetObj$cmdSet <- vector(mode="character"); # store R command
  # set global variables
  current.msg <<- "";
  err.vec <<- "";
  current.selected.tax <<- "NA";
  enrich.type <<- "hyper";
  
  # to control parallel computing for some packages
  BiocParallel::register(BiocParallel::SerialParam());
  Sys.setenv("OMP_NUM_THREADS" = 2); 
  Sys.setenv("OPENBLAS_NUM_THREADS" = 2);
  
  # preload some general package
  load_cairo();
  load_ggplot();
  Cairo::CairoFonts("Arial:style=Regular","Arial:style=Bold","Arial:style=Italic","Helvetica","Symbol")
  print("Init MicrobiomeAnalyst!");
  return(.set.mbSetObj(mbSetObj))
}

#' Read RDS files from the internet
#' @description Function downloads the required file and reads it only if not already in working directory.
#' Need to specify the file URL and the destfile. 
#' @param filenm Input the name of the file to download.
#' @param opt Default set to "none".
#' @param ref Default set to "NA".

# read binary RDS files
.read.microbiomeanalyst.lib.rds <- function(filenm, sub.dir = NULL, ref = NA){
  
  if(.on.public.web){
    if(is.null(sub.dir)){
      lib.path <- paste("../../lib/", filenm, sep="");
    }else{
      lib.path <- paste("../../lib/", sub.dir, "/", filenm,  sep="");
    }
    return(readRDS(lib.path));
  }else{
    lib.download <- FALSE;
    file_name <- basename(filenm)
    if(!file.exists(filenm)){
      lib.download <- TRUE;
    }else{
      time <- file.info(filenm)
      diff_time <- difftime(Sys.time(), time[,"mtime"], unit="days") 
      if(diff_time>30){
        lib.download <- TRUE;
      }
    }
    # Deal with curl issues
    if(lib.download){
      if(sub.dir == "tsea"){
        lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/tsea/", filenm, sep="");
      }else if(sub.dir == "picrust"){
        lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/picrust/", filenm, sep="");
      }else if(sub.dir == "ppd"){
        lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/ppd", filenm, sep = "");
      }else{
        lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/", filenm, sep="");
      }
      tryCatch(
        {
          download.file(lib.url, destfile=file_name, method="curl")
        }, warning = function(w){ print() },
        error = function(e) {
          print("Download unsucceful. Ensure that curl is downloaded on your computer.")
          print("Attempting to re-try download using libcurl...")
          download.file(lib.url, destfile=file_name, method="libcurl")
        }
      )
    }
    lib.path <- filenm;
  }
  
  # Deal w. corrupt downloaded files
  tryCatch({
    my.lib <- readRDS(file_name); # this is a returned value, my.lib never called outside this function, should not be in global env.
    print("Loaded files from MetaboAnalyst web-server.")
  },
  warning = function(w) { print("Warning, files not successfully downloaded from web.") },
  error = function(err) {
    print("Reading data unsuccessful, attempting to re-download file...")
    tryCatch(
      {
        if(sub.dir == "tsea"){
          lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/tsea/", filenm, sep="");
        }else if(sub.dir == "picrust"){
          lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/picrust/", filenm, sep="");
        }else if(sub.dir == "ppd"){
          lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/ppd", filenm, sep = "");
        }else{
          lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/", filenm, sep="");
        }
        
        download.file(lib.url, destfile=file_name, method="curl")
        my.lib <- readRDS(file_name);
        print("Loaded necessary files.")
      },
      warning = function(w) { print() },
      error = function(err) {
        print("Loading files from server unsuccessful. Ensure curl is downloaded on your computer.")
      }
    )
  })
  return(my.lib)
}

#'Function to set analysis type
#'@description This functions sets the module name.
#'@param analType Input the analysis type. If the data is marker gene data, 
#'use "mdp", if the data is shotgun metagenomics or transcriptomics data, 
#'use ""sdp". If performing the Projection with Public Data module, use "ppd". 
#'If performing Taxon Set Enrichment Analysis, use "tsea". 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SetModuleType <- function(mbSetObj, nm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$module.type <- nm;
  if(mbSetObj$module.type!="mmp"){
  mbSetObj$micDataType = "na"
  }
  return(.set.mbSetObj(mbSetObj));
}

#'Function to read in sample data
#'@description This functions reads in sample data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the sample data file name.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
ReadSampleTable <- function(mbSetObj, fileName) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  module.type <- mbSetObj$module.type;
  load_phyloseq();
  mydata <- .readDataTable(fileName);
  mydata[is.na(mydata)] <- "Not Available";
  if(any(is.na(mydata)) || class(mydata) == "try-error"){
    AddErrMsg("Failed to read in the metadata! Please make sure that the metadata file is in the right format and does not have empty cells or NA.");
    return(0);
  }
  
  # look for #NAME, store in a list
  #sam.nm <- substr(colnames(mydata[1]),1,5);
  #sam.nm <- tolower(sam.nm);
  #sam.inx <- grep("#name",sam.nm);
  #if(length(sam.inx) > 0){
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
  my.meta.blank <- any(grepl("[[:blank:]]", my.meta)| grepl("[[:blank:]]", names(my.meta)));
  if(my.meta.blank){
    names(my.meta) <- gsub("\\s+","_", names(my.meta));
    rownms <- rownames(my.meta)
    my.meta <- data.frame(sapply(my.meta, function(x) gsub(" ","_",x)));
   rownames(my.meta) <- rownms
    na.msg1 <- c(na.msg1, "Blank spaces in group names are replaced with underscore '_'");
  }

  mbSetObj$dataSet$group_names <- colnames(my.meta)
  #na.msg <- paste(na.msg, "The sample data contains a total of ", nrow(mydata), "samples and  ", ncol(mydata), " sample variables.", collapse=" ");

  # as most functions are for discrete groups (not continuous values)
  # require at least one column contains discrete factors with at least two replicates 
  disc.inx <- GetDiscreteInx(my.meta);

  if(sum(disc.inx) == 0){ # No discrete column detected
    AddErrMsg("Metadata Table: make sure there is at least one column contains experimental design for group comparisons (i.e., the primary metadata), with each group contains at least 3 replicate. No unique values are allowed in the primary metadata column.");
    na.msg <- c("<font style=\"color:red\"><b> No.</b></font>", "It seems that some of your metadata values are unique! Please make sure at least 3 replicate per group!");
    mbSetObj$poor.replicate <- TRUE;
    mbSetObj$dataSet$sample_data <- my.meta
    mbSetObj$dataSet$meta_info$disc.inx <- 0;
    mbSetObj$dataSet$meta_info$cont.inx <- 0;
    return(0);
  }else{
    
    na.msg <- c(na.msg, "<font style=\"color:green\"><b> Yes. </b></font>");
    if(sum(disc.inx) == length(disc.inx)){
      na.msg <- c(na.msg,"All metadata columns are OK!")
    }else{
      bad.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
      na.msg <- c(na.msg, paste0("<font style=\"color:red\">Detected presence of unique values, the following metadata columns are excluded: <b>", bad.meta, "</b></font>"));
    }
    
    mbSetObj$dataSet$meta_info$disc.inx <- disc.inx;
    mbSetObj$dataSet$sample_data <- my.meta[,disc.inx, drop=FALSE];
    
    cont.inx <- GetNumbericalInx(my.meta);
    cont.inx <- !disc.inx & cont.inx; # discrete is first
    mbSetObj$dataSet$meta_info$cont.inx <- cont.inx;
    
    if(sum(cont.inx)>0){
      # make sure the discrete data is on the left side
      mbSetObj$dataSet$sample_data <- cbind(mbSetObj$dataSet$sample_data, my.meta[,cont.inx, drop=FALSE]);
    }
    
  }
  mbSetObj$dataSet$smpl.msg <- c(na.msg,na.msg1);
  return(.set.mbSetObj(mbSetObj));
}

#'Function to read in tree files.
#'@description This functions reads in tree files.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param dataName Input the tree file name.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
ReadTreeFile <- function(mbSetObj, fileName, dataName="",module.type) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_phyloseq();
  tree <- tryCatch({
    read_tree(fileName);
  })  
  if(!is.null(tree)){
    saveDataQs(tree ,"tree.qs", module.type, dataName);
    mbSetObj$tree.uploaded <- TRUE;
    return(.set.mbSetObj(mbSetObj));
  }else{
    AddErrMsg("Failed to parse tree file data!")
    return(0)
  }
}

IsTreeUploaded <- function(mbSetObj) {
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$tree.uploaded){
    return(1);
  }else{
    return(0)
  }
}

RecordRCommand <- function(mbSetObj=NA, cmd){
  write(cmd, file = "Rhistory.R", append = TRUE);
  mbSetObj <- .get.mbSetObj(mbSetObj); 
  mbSetObj$cmdSet <- c(mbSetObj$cmdSet, cmd);
  return(.set.mbSetObj(mbSetObj));
}

GetRCommandHistory <- function(mbSetObj=NA){
  mbSetObj <- .get.mbSetObj(mbSetObj); 
  if(is.null(mbSetObj$cmdSet)){
    return("NA");
  }
  return(mbSetObj$cmdSet);
}

ClearRCommandHistory <- function(mbSetObj=NA){
  mbSetObj <- .get.mbSetObj(mbSetObj); 
  mbSetObj$cmdSet <- c();
  return(.set.mbSetObj(mbSetObj));
}

####################################
############ Get Funs ##############
####################################
GetResRowNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rownames(mbSetObj$analSet$resTable));
}

GetResColNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(colnames(mbSetObj$analSet$resTable));
}

GetMetaboResRowNames <- function(){
  return(rownames(current.proc$met$res_deAnal));
}

GetMetaboResColNames <- function(){
  return(colnames(current.proc$met$res_deAnal));
}


GetMetadataMsg <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(paste0(mbSetObj$dataSet$smpl.msg, collapse="; "));
}


IsPoorReplicate <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$poor.replicate){
    return(1);
  }
  return(0);
}

GetResMat <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(as.matrix(mbSetObj$analSet$resTable));
}

GetResMetabo <- function(){
  return(as.matrix(current.proc$met$res_deAnal));
}


# type can be all, discrete or continuous
GetMetaInfo <- function(mbSetObj, type="disc"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$module.type == "meta"){
    microbiome.meta <- qs:::qread("microbiome_meta.qs");
    meta.types <- microbiome.meta$meta.types;
  }else{
    meta.types <- mbSetObj$dataSet$meta.types
  }
  meta.nms <- names(meta.types)
  if(type=="all"){

  }else if(type=="disc"){
    meta.nms <- meta.nms[meta.types == "disc"]
  }else if(type=="cont"){
    meta.nms <- meta.nms[meta.types == "cont"]
    if(length(meta.nms) == 0){
      return("NA");
    }
  }
    return(meta.nms);

}

GetMetaTaxaInfo <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  proc.phyobj <- mbSetObj$dataSet$proc.phyobj;

  if(mbSetObj$module.type == "meta"){
    if(file.exists("merged.data.qs")){
     proc.phyobj <- qs::qread("merged.data.qs");
     proc.phyobj <- subsetPhyloseqByDataset(mbSetObj, proc.phyobj);
     }
  }

  #check that each rank has >2 groups
  taxa.tbl <- as(tax_table(proc.phyobj), "matrix")

  if(ncol(taxa.tbl)==1){
    taxa.nms <- "Phylum"
    return(taxa.nms)
  }
  
  #drop taxa with only 1 level (i.e. Viruses at Phylum)
  gd.inx <- apply(taxa.tbl, 2, function(x) length(unique(x))!=1);
  taxa.tbl.update <- taxa.tbl[,gd.inx, drop=FALSE];
  
  if(ncol(taxa.tbl.update) == 0){
    current.msg <<- c("All taxa info for the remaining features are the same!")
    return("OTU")
  }
  
  taxa.nms <- rank_names(taxa.tbl.update);
  return(taxa.nms[!is.na(taxa.nms)]);
}


GetSampleGrpInfo <- function(mbSetObj, clsLbl){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(levels(factor(get_variable(mbSetObj$dataSet$norm.phyobj, clsLbl))));
}

GetSampleGrpNo <- function(mbSetObj, clsLbl){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  #Issue with phyloslim (after merging into phyloslim object the sample variable are converted to numeric again rather than factor)
  return(length(levels(factor(get_variable(mbSetObj$dataSet$norm.phyobj, clsLbl)))));
}

GetTaxaNames<- function(mbSetObj, taxlvl){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
   if(!exists("phyloseq_objs")){
 phyloseq_objs <- readDataQs("phyloseq_objs.qs",mbSetObj$module.type,mbSetObj$dataSet$name)
  }
  
  nm = taxa_names(phyloseq_objs$merged_obj[[taxlvl]])
  
  if(sum(is.na(nm))/length(nm) > 0.7){
      AddErrMsg("More than 70% values are missing at this level!");
      return(0);
    }
    indx<-which(is.na(nm)==TRUE);
    nm[indx]<-"Not_Assigned";
    return(nm);

}

GetTaxaFeatSize<- function(mbSetObj, taxlvl){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(!exists("phyloseq_objs")){
 phyloseq_objs <- readDataQs("phyloseq_objs.qs",mbSetObj$module.type,mbSetObj$dataSet$name)
  }

   feat.size <- length(unique(as.character(rownames(phyloseq_objs$count_tables[[taxlvl]]))));

  return(feat.size);
}

GetLowerTaxaLvlNm<- function(mbSetObj, taxrank){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  all.rnks <- colnames(tax_table(mbSetObj$dataSet$proc.phyobj));
  bad.inx <- is.na(all.rnks) | nchar(all.rnks)==0;
  gd.rnks <- all.rnks[!bad.inx];
  indx <- which(gd.rnks==taxrank);
  low.rnks <- gd.rnks[(indx+1):length(gd.rnks)];
  return(low.rnks);
}

GetHighTaxaLvlNm<- function(mbSetObj, taxrank){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(taxrank=="OTU"){
    return(colnames(tax_table(mbSetObj$dataSet$proc.phyobj))[1:length(colnames(tax_table(mbSetObj$dataSet$proc.phyobj)))]);
  }else{
    indx <- which(colnames(tax_table(mbSetObj$dataSet$proc.phyobj))==taxrank);
    rem <- (indx):1;
    return(colnames(tax_table(mbSetObj$dataSet$proc.phyobj))[rev(rem)]);
  }
}

##########################
######## Checks ##########
##########################

ValidateFeatureName <- function(mbSetObj, taxlvl, nm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(taxlvl=="OTU"){
    tax.nms <- taxa_names(mbSetObj$dataSet$norm.phyobj);
  }else{
    taxa_table<-tax_table(mbSetObj$dataSet$proc.phyobj);
    data<-merge_phyloseq(mbSetObj$dataSet$norm.phyobj,taxa_table);
    tax.nms<-unique(as.character(tax_table(data)[,taxlvl]));
    if(sum(is.na(tax.nms))/length(tax.nms) > 0.7){
      AddErrMsg("More than 70% values are missing at this level!");
      return(0);
    }
  }
  if(nm %in% tax.nms){
    return(1);
  }else{
    return(0);
  }
}

# save the processed data with class names
PrepareDownloadData <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$module.type == "mdp" || mbSetObj$module.type == "sdp"){
    if(!is.null(mbSetObj$dataSet$data.orig)){
      fast.write(mbSetObj$dataSet$data.orig, file="data_original.csv");
    }
    if(!is.null(mbSetObj$dataSet$filt.data)){
      fast.write(mbSetObj$dataSet$filt.data, file="data_filtered.csv");
    }
    if(!is.null(otu_table(mbSetObj$dataSet$norm.phyobj,as.matrix()))){
      fast.write(otu_table(mbSetObj$dataSet$norm.phyobj,as.matrix()), file="data_normalized.csv");
    }
  }else if(mbSetObj$module.type == "ppd"){
    fast.write(otu_table(userrefdata,as.matrix()), file="merge_otutable.csv");
    fast.write(as.matrix(tax_table(userrefdata)), file="merge_taxtable.csv");
    fast.write(as.data.frame(sample_data(userrefdata),check.names=FALSE), file="merge_sampletable.csv");
  }
  return(.set.mbSetObj(mbSetObj));
};

## utility functions to create phyloseq obs + count tables

UtilMakePhyloseqObjs <- function(mbSetObj, taxrank){

# if(mbSetObj$module.type=="mdp" |mbSetObj$micDataType=="otu" ){
#   
#   taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
#   data <- merge_phyloseq(mbSetObj$dataSet$proc.phyobj, taxa_table);
# }else{
#   data <- mbSetObj$dataSet$proc.phyobj; #for shotgun
# }

  data <- mbSetObj$dataSet$proc.phyobj;
  if(taxrank!="OTU"){
    #merging at taxonomy levels
    # JE note Dec 14, 2022: We are condensing normalized, filtered data. Is this correct? How are the data condensed? 
    # If a sum/mean operation, condensing before/after normalization could have a big impact.
    data <- fast_tax_glom_mem(data, taxrank);

    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
  }
  return(data)
}

UtilMakeCountTables <- function(phyloseq.obj, taxrank){

  if(taxrank=="OTU"){
    data1 <- as(otu_table(phyloseq.obj), "matrix");

  }else{
    nm <- as.character(tax_table(phyloseq.obj)[,taxrank]);
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Not_Assigned";
    data1 <- as(otu_table(phyloseq.obj), "matrix");
    rownames(data1) <- nm;
    #all NA club together
    data1 <- as.matrix(t(sapply(by(data1, rownames(data1), colSums), identity)));
  }

  return(data1)
}

saveSet <- function(obj=NA, set="", output=1){
  
  if(globalConfig$anal.mode == "api"){ 
    qs:::qsave(obj, paste0(set, ".qs"));
    return(output)
  }else{
    if(set == "dataSet"){
      dataSet <<- obj;
    }else if(set == "analSet"){
      analSet <<- obj;
    }else if(set == "imgSet"){
      imgSet <<- obj;
    }else if(set == "paramSet"){
      head(paramSet);
      paramSet <<- obj;
    }else if(set == "msgSet"){
      msgSet <<- obj;
    }else if(set == "cmdSet"){
      cmdSet <<- obj;
    }
    
    if(globalConfig$anal.mode == "web"){
      return(output);
    }else{
      return(obj);
    }
  }
  
}

readSet <- function(obj=NA, set=""){
  if(globalConfig$anal.mode == "api"){
    path <- "";
    if(exists('user.path')){
      path <- user.path;
    }
    
    if(path != ""){
      obj <- load_qs(paste0(path, set, ".qs"));
    }else{
      obj <- qs:::qread(paste0(set, ".qs"));
    }
  }
  return(obj);
}

load_qs <- function(url) qs::qdeserialize(curl::curl_fetch_memory(url)$content)

readDataset <- function(fileName="", obj=NA){

  if(globalConfig$anal.mode == "api"){
    if(exists('user.path')){
      path <- user.path;
      dat <- load_qs(paste0(path, fileName));
    }else{
      dat <- qs:::qread(fileName);
    }
  }else{
    if(!is.na(obj)){
      dat <- obj$dataSets[[fileName]];
      
    }else{
      mbSetObj <- .get.mbSetObj(NA);
      dat <- mbSetObj$dataSets[[fileName]];
    }
  }
  
  return(dat);
}

Set.Config <-function(anal.mode="web"){
  
  globalConfig <- list();
  globalConfig$anal.mode <- anal.mode
  globalConfig <<- globalConfig;
}

saveDataQs <-function(data, name, module.nm, dataName){
  if(module.nm == "meta"){
    qs::qsave(data, file=paste0(dataName, "_data/", name));
  }else{
    qs::qsave(data, file=name);
  }
}

readDataQs <-function(name, module.nm, dataName){
  if(module.nm == "meta"){
    dat <- qs::qread(file=paste0(dataName, "_data/", name));
  }else{
    dat <- qs::qread(file=name);
  }
  return(dat);
}
