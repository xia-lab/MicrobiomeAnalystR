##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# This is only for web version
.on.public.web <- FALSE; # only TRUE when on MicrobiomeAnalyst web server

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
  
  dataSet <- list();
  analSet <- list();
  imgSet <- list();
  
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
  current.selected.tax <<- "NA";
  enrich.type <<- "hyper";
  
  load_cairo();
  load_ggplot();
  BiocParallel::register(BiocParallel::SerialParam());
  
  # preload some general package
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
ReadSampleTable<- function(mbSetObj, dataName) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_phyloseq();
  
  msg <- NULL;
  mydata <- .readDataTable(dataName);
  
  if(any(is.na(mydata)) || class(mydata) == "try-error"){
    AddErrMsg("Failed to read in the metadata! Please make sure that the metadata file is in the right format and does not have empty cells or NA.");
    return(0);
  }
  
  # look for #NAME, store in a list
  sam.nm <- substr(colnames(mydata[1]),1,5);
  sam.nm <- tolower(sam.nm);
  sam.inx <- grep("#name",sam.nm);
  
  if(length(sam.inx) > 0){
    smpl_nm<-mydata[,1];
    smpl_var<-colnames(mydata[-1]);
  }else{
    AddErrMsg("Please make sure you have the label #NAME in your sample data file!");
    return(0);
  }
  
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
  na.msg <- NULL;
  
  if(sum(na.inx) > 0){
    mydata[na.inx] <- "Unknown";
    na.msg <- paste("A total of", sum(na.inx), "empty or NA values were replaced by 'Unknown'.");
  }
  
  # as most functions are for discrete groups (not continuouse values
  # require at least one column contains discrete factors with at least two replicates 
  my.meta <- data.frame(mydata);
  mbSetObj$dataSet$group_names <- colnames(my.meta)

  disc.inx <- GetDiscreteInx(my.meta);
  if(sum(disc.inx) == 0){ # all class labels are unique! 
    na.msg <- c(na.msg, "It seems that all your meta data values are unique! MicrobiomeAnalyst requires some biological replicates for robust analysis");
    mbSetObj$poor.replicate <- TRUE;
    mbSetObj$dataSet$sample_data <- my.meta
  }else{
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

  current.msg <<- paste(na.msg, "The sample data contains a total of ", nrow(mydata), "samples and  ", ncol(mydata), " sample variables.", collapse=" ");
  mbSetObj$dataSet$smpl.msg <- current.msg;
  
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
ReadTreeFile <- function(mbSetObj, dataName) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_phyloseq();
  msg <- NULL;
  
  tree <- tryCatch({
    read_tree(dataName);
  })
  
  if(!is.null(tree)){
    qs::qsave(tree, "tree.qs");
    mbSetObj$tree.uploaded <- TRUE;
    return(.set.mbSetObj(mbSetObj));
  }else{
    print("Failed to parse tree file data!")
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
  mbSetObj <- .get.mbSetObj(mSetObj); 
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

# type can be all, discrete or continuous
GetMetaInfo <- function(mbSetObj, type="disc"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  all.nms <- mbSetObj$dataSet$group_names
  if(type=="all"){
    return(all.nms);
  }else if(type=="disc"){
    return(all.nms[mbSetObj$dataSet$meta_info$disc.inx]);
  }else if(type=="cont"){
    cont.inx <- mbSetObj$dataSet$meta_info$cont.inx;
    if(sum(cont.inx) == 0){
      return("NA");
    }
    return(all.nms[mbSetObj$dataSet$meta_info$cont.inx]);
  }else{
    print(paste("Unknown type:", type));
    return("NA");
  }
}

GetMetaTaxaInfo <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  #check that each rank has >2 groups
  taxa.tbl <- as(tax_table(mbSetObj$dataSet$proc.phyobj), "matrix")
  
  if(ncol(taxa.tbl)==1){
    taxa.nms <- "Phylum"
    return(taxa.nms)
  }
  
  #drop taxa with only 1 level (i.e. Viruses at Phylum)
  gd.inx <- apply(taxa.tbl, 2, function(x) length(unique(x))!=1);
  taxa.tbl.update <- taxa.tbl[,gd.inx, drop=FALSE];

  if(ncol(taxa.tbl.update) == 0){
    msg <- c("All taxa info for the remaining features are the same!")
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
  
  if(taxlvl=="OTU"){
    return(taxa_names(mbSetObj$dataSet$norm.phyobj));
  }else{
    taxa_table<-tax_table(mbSetObj$dataSet$proc.phyobj);
    data<-merge_phyloseq(mbSetObj$dataSet$norm.phyobj,taxa_table);
    nm<-unique(as.character(tax_table(data)[,taxlvl]));
    if(sum(is.na(nm))/length(nm) > 0.7){
      AddErrMsg("More than 70% values are missing at this level!");
      return(0);
    }
    indx<-which(is.na(nm)==TRUE);
    nm[indx]<-"Not_Assigned";
    return(nm);
  }
}

GetTaxaFeatSize<- function(mbSetObj, taxlvl){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(taxlvl=="OTU"){
    feat.size <- length(taxa_names(mbSetObj$dataSet$norm.phyobj));
  }else{
    taxa_table<-tax_table(mbSetObj$dataSet$proc.phyobj);
    data<-merge_phyloseq(mbSetObj$dataSet$norm.phyobj,taxa_table);
    feat.size <- length(unique(as.character(tax_table(data)[,taxlvl])));
  }
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
    fast.write(as.data.frame(sample_data(userrefdata)), file="merge_sampletable.csv");
  }
  return(.set.mbSetObj(mbSetObj));
};

## utility functions to create phyloseq obs + count tables

UtilMakePhyloseqObjs <- function(mbSetObj, taxrank){
  
  if(mbSetObj$module.type=="mdp"){
    taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(mbSetObj$dataSet$norm.phyobj, taxa_table);
  }else{
    data <- mbSetObj$dataSet$norm.phyobj; #for shotgun
  }
  
  if(taxrank!="OTU"){
    #merging at taxonomy levels
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

# only used when data is norm.phyobj
MakeRankedCountTables <- function(mbSetObj){
  
  # make hierarchies
  ranks <- c(GetMetaTaxaInfo(mbSetObj), "OTU")
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
  return(.set.mbSetObj(mbSetObj))
}
