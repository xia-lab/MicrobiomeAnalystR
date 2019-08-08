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
    return(mbSet)
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
  mbSetObj$cmdSet <- vector(mode="character"); # store R command

  # set global variables
  current.msg <<- "";
  current.selected.tax <<- "NA";
  enrich.type <<- "hyper";
  msg.vec <<- vector(mode="character");

  if(.on.public.web){
    load_cairo();
    load_ggplot();
    load_biocparallel();
    BiocParallel::register(BiocParallel::SerialParam());
  }

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
.read.microbiomeanalyst.lib <- function(filenm, opt = "none", ref = NA){
  
  if(.on.public.web){
    if(opt=="tsea"){
      lib.path <- paste("../../lib/tsea/", filenm, sep="");
    }else if(opt=="ppd"){
      lib.path <- paste(paste("../../lib/ref_data/", ref, sep=""), filenm, sep="/");
    }else{
      lib.path <- paste("../../lib/", filenm, sep="");
    }
    return(readRDS(lib.path));
  }else{
    lib.download <- FALSE;
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
      if(opt == "tsea"){
        lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/tsea/", filenm, sep="");
      }else if(opt == "picrust"){
        lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/picrust/", filenm, sep="");
      }else if(opt == "ppd"){
        lib.url <- paste(paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/ref_data/", ref, sep=""), filenm, sep = "/");
      }else{
        lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/", filenm, sep="");
      }
      tryCatch(
        {
          download.file(lib.url, destfile=filenm, method="curl")
        }, warning = function(w){ print() },
        error = function(e) {
          print("Download unsucceful. Ensure that curl is downloaded on your computer.")
          print("Attempting to re-try download using libcurl...")
          download.file(lib.url, destfile=filenm, method="libcurl")
        }
      )
    }
    lib.path <- filenm;
  }
  
  # Deal w. corrupt downloaded files
  tryCatch({
    my.lib <- readRDS(lib.path); # this is a returned value, my.lib never called outside this function, should not be in global env.
    print("Loaded files from MetaboAnalyst web-server.")
  },
  warning = function(w) { print("Warning, files not successfully downloaded from web.") },
  error = function(err) {
    print("Reading data unsuccessful, attempting to re-download file...")
    tryCatch(
      {
        if(opt == "tsea"){
          lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/tsea/", filenm, sep="");
        }else if(opt == "picrust"){
          lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/picrust/", filenm, sep="");
        }else if(opt == "ppd"){
          lib.url <- paste(paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/ref_data/", ref, sep=""), filenm, sep = "/");
        }else{
          lib.url <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/", filenm, sep="");
        }

        download.file(lib.url, destfile=filenm, method="curl")
        my.lib <- readRDS(filenm);
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

# read binary RDA files (old style should be all RDS)
# type should mset or kegg
.load.microbiomeanalyst.lib <- function(libname){
  
  destfile <- libname;
  if(.on.public.web){
    destfile <- paste("../../lib/", libname, sep="");
  }else{
    lib.download <- FALSE;
    if(!file.exists(destfile)){
      lib.download <- TRUE;
    }else{
      time <- file.info(destfile)
      diff_time <- difftime(Sys.time(), time[,"mtime"], unit="days") 
      if(diff_time>30){
        lib.download <- TRUE;
      }
    }
    if(lib.download){
      libPath <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/", libname, sep="");
      download.file(libPath, destfile);
    }
  }
  load(destfile, .GlobalEnv);  
}

#'Function to set analysis type
#'@description This functions sets the analysis mode.
#'@param analType Input the analysis type. If the data is marker gene data, 
#'use "markergene", if the data is shotgun metagenomics or transcriptomics data, 
#'use ""shotgun". If performing the Projection with Public Data module, use "dataprojection". 
#'If performing Taxon Set Enrichment Analysis, use "species". 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SetAnalType <- function(analType){
  anal.type <<- analType;
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
  
  if(.on.public.web){
    load_phyloseq();
  }
  
  msg <- NULL;
  mydata <- .readDataTable(dataName);
  
  if(any(is.na(mydata)) || class(mydata) == "try-error"){
    current.msg <<- "Failed to read in the metadata! Please make sure that the metadata file is in the right format and does not have empty cells or NA.";
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
    current.msg <<- "Please make sure you have the label #NAME in your sample data file!";
    return(0);
  }
  
  # converting to character matrix as duplicate row names not allowed in data frame.
  mydata <- as.matrix(mydata[,-1]);
  rownames(mydata) <- smpl_nm;
  colnames(mydata) <- smpl_var;
  
  # empty cell or NA cannot be tolerated in metadata
  na.inx  <- is.na(mydata);
  na.msg <- NULL;
  
  if(sum(na.inx) > 0){
    mydata[na.inx] <- "Unknown";
    na.msg <- paste("A total of", sum(na.inx), "empty or NA values were replaced by 'Unknown'.");
  }
  
  mbSetObj$dataSet$sample_data <- data.frame(mydata);
  current.msg <<- paste(na.msg, "The sample data contains a total of ", nrow(mydata), "samples and  ", ncol(mydata), " sample variables.", collapse=" ");
  mbSetObj$dataSet$smpl.msg <- current.msg;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)  
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj));
  }
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
ReadTreeFile<- function(mbSetObj, dataName) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(.on.public.web){
    load_phyloseq();
  }
  
  msg <- NULL;
  
  tree <- tryCatch(
    {read_tree(dataName)}
  )
  
  if(!is.null(tree)){
    saveRDS(tree, "tree.RDS");
  }else{
    return(0)
  }
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)  
    return(1);
  }else{
    print("Tree file successfully uploaded!")
    return(.set.mbSetObj(mbSetObj));
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
  return(mbSetObj$cmdSet);
}

####################################
############ Get Funs ##############
####################################

GetNameMapCol <-function(mbSetObj, colInx){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(.on.public.web){
    return(mbSetObj$analSet$resTable[,colInx]);
  }else{
    print(mbSetObj$analSet$resTable[,colInx]);
    return(.set.mbSetObj(mbSetObj))
  }
}

GetResRowNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rownames(mbSetObj$analSet$resTable));
}

GetTseaRowNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rownames(mbSetObj$analSet$tseaInfo));
}

GetTseaCol <-function(mbSetObj, colInx){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$resTable[,colInx]);
}

GetMirResCol <-function(mbSetObj, colInx){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$resTable[,colInx]);
}

GetResColNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(colnames(mbSetObj$analSet$resTable));
}

GetResMat <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(as.matrix(mbSetObj$analSet$resTable));
}

GetMetaInfo <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(colnames(mbSetObj$dataSet$sample_data));
}

GetConfounderOpts <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  metadata <- colnames(mbSetObj$dataSet$sample_data)
  if(length(metadata)==1){
    return("NA")
  }else{
    return(metadata)
  }
}

GetMetaTaxaInfo <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rank_names(mbSetObj$dataSet$proc.phyobj));
}

GetSampleGrpInfo <- function(mbSetObj, class){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(levels(factor(get_variable(mbSetObj$dataSet$norm.phyobj, class))));
}

GetSampleGrpNo <- function(mbSetObj, class){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  #Issue with phyloslim (after merging into phyloslim object the sample variable are converted to numeric again rather than factor)
  return(length(levels(factor(get_variable(mbSetObj$dataSet$norm.phyobj, class)))));
}

GetAllSampleGrpInfo <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  #sample variable having more than one group will be selected as default
  sam_var<-which(sapply(sample_data(mbSetObj$dataSet$norm.phyobj)[,sapply(sample_data(mbSetObj$dataSet$norm.phyobj), is.factor)], nlevels)>1);
  if(length(sam_var)>0){
    return(names(sam_var[1]));
  } else {
      return(NULL);
  }
}

GetTaxaFeatName<- function(mbSetObj, taxlvl){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
    
  if(taxlvl=="OTU"){
    return(taxa_names(mbSetObj$dataSet$norm.phyobj));
  }else{
    taxa_table<-tax_table(mbSetObj$dataSet$proc.phyobj);
    data<-merge_phyloseq(mbSetObj$dataSet$norm.phyobj,taxa_table);
    nm<-unique(as.character(tax_table(data)[,taxlvl]));
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
  indx <- which(colnames(tax_table(mbSetObj$dataSet$proc.phyobj))==taxrank);
  rem <- ncol(tax_table(mbSetObj$dataSet$proc.phyobj))-indx;
  return(colnames(tax_table(mbSetObj$dataSet$proc.phyobj))[indx+1:rem]);
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

GetSampleGrpUser <- function(mbSetObj, class){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(levels(get_variable(mbSetObj$dataSet$proc.phyobj, class)));
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
  }
  
  if(nm %in% tax.nms){
    if(.on.public.web){
      .set.mbSetObj(mbSetObj)  
      return(1);
    }else{
      return(.set.mbSetObj(mbSetObj));
    }
  }
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)  
    return(0);
  }else{
    return(.set.mbSetObj(mbSetObj));
  }
}

#check whether sample variable is continuous and give warning abt it in metagenomeSeq and LEfSe
CheckContSampleVar <- function(mbSetObj, variable){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  cls<-as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  
  if(length(cls)/length(levels(cls)) < 2){
    current.msg <<-"Sample variable having continuous values. This method is only applicable for variables containing discrete values.";
    if(.on.public.web){
      .set.mbSetObj(mbSetObj)  
      return(0);
    }else{
      return(.set.mbSetObj(mbSetObj));
    }
  } else {
    if(.on.public.web){
      .set.mbSetObj(mbSetObj)  
      return(1);
    }else{
      return(.set.mbSetObj(mbSetObj));
    }
  }
}

# save the processed data with class names
SaveData <- function(mbSetObj, anal.type){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(anal.type == "16S" || anal.type == "metageno"){
    if(!is.null(mbSetObj$dataSet$data.orig)){
      write.csv(mbSetObj$dataSet$data.orig, file="data_original.csv");
    }
    if(!is.null(mbSetObj$dataSet$filt.data)){
      write.csv(mbSetObj$dataSet$filt.data, file="data_filtered.csv");
    }
    if(!is.null(otu_table(mbSetObj$dataSet$norm.phyobj,as.matrix()))){
      write.csv(otu_table(mbSetObj$dataSet$norm.phyobj,as.matrix()), file="data_normalized.csv");
    }
  }else if(anal.type == "16S_ref"){
    write.csv(otu_table(userrefdata,as.matrix()), file="merge_otutable.csv");
    write.csv(as.matrix(tax_table(userrefdata)), file="merge_taxtable.csv");
    write.csv(as.data.frame(sample_data(userrefdata)), file="merge_sampletable.csv");
  }
  return(.set.mbSetObj(mbSetObj));
}