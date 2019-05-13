##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

# This is only for web version
.on.public.web <- TRUE; # only TRUE when on metaboanalyst web server

# note, this is usually used at the end of a function
# for local, return itself; for web, push to global environment
.set.microSet <- function(microSetObj=NA){
  if(.on.public.web){
    microSet <<- microSetObj;
    return (1);
  }
  return(microSetObj);
}

.get.microSet <- function(microSetObj=NA){
  if(.on.public.web){
    return(microSet)
  }else{
    return(microSetObj);
  }
}

#'Constructs a microSet object for storing data 
#'@description This functions handles the construction of a microSetObj object for storing data 
#'produced by MicrobiomeAnalystR for further processing and analysis.
#'@usage Init.Data()
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
Init.Data <- function(){
  
  dataSet <- list();
  analSet <- list();
  imgSet <- list();
  
  microSetObj <- list();
  microSetObj$dataSet <- dataSet;
  microSetObj$analSet <- analSet;
  microSetObj$imgSet <- imgSet;
  
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
  return(.set.microSet(microSetObj))
}

SetAnalType <- function(analType){
  anal.type <<- analType;
}

GetNameMapCol <-function(microSetObj, colInx){
  microSetObj <- .get.microSetObj(microSetObj);
  return(microSetObj$analSet$resTable[,colInx]);
}

GetResRowNames <- function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(rownames(microSetObj$analSet$resTable));
}

GetTseaRowNames <- function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(rownames(microSetObj$analSet$tseaInfo));
}

GetTseaCol <-function(microSetObj, colInx){
  microSetObj <- .get.microSetObj(microSetObj);
  return(microSetObj$analSet$resTable[,colInx]);
}

GetMirResCol <-function(microSetObj, colInx){
  microSetObj <- .get.microSetObj(microSetObj);
  return(microSetObj$analSet$resTable[,colInx]);
}

GetResColNames <- function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(colnames(microSetObj$analSet$resTable));
}

GetResMat <- function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(as.matrix(microSetObj$analSet$resTable));
}

GetMetaInfo <- function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(colnames(microSetObj$dataSet$sample_data));
}

GetMetaTaxaInfo <- function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(rank_names(microSetObj$dataSet$proc.phyobj));
}

GetSampleGrpInfo <- function(microSetObj, class){
  microSetObj <- .get.microSetObj(microSetObj);
  return(levels(factor(get_variable(microSetObj$dataSet$norm.phyobj, class))));
}

GetSampleGrpNo <- function(microSetObj, class){
  microSetObj <- .get.microSetObj(microSetObj);
  #Issue with phyloslim (after merging into phyloslim object the sample variable are converted to numeric again rather than factor)
  return(length(levels(factor(get_variable(microSetObj$dataSet$norm.phyobj, class)))));
}

GetAllSampleGrpInfo <- function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  #sample variable having more than one group will be selected as default
  sam_var<-which(sapply(sample_data(microSetObj$dataSet$norm.phyobj)[,sapply(sample_data(microSetObj$dataSet$norm.phyobj), is.factor)], nlevels)>1);
  if(length(sam_var)>0){
    return(names(sam_var[1]));
  } else {
      return(NULL);
  }
}

GetTaxaFeatName<- function(microSetObj, taxlvl){
  
  microSetObj <- .get.microSetObj(microSetObj);
    
  if(taxlvl=="OTU"){
    return(taxa_names(microSetObj$dataSet$norm.phyobj));
  }else{
    taxa_table<-tax_table(microSetObj$dataSet$proc.phyobj);
    data<-merge_phyloslim(microSetObj$dataSet$norm.phyobj,taxa_table);
    nm<-unique(as.character(tax_table(data)[,taxlvl]));
    indx<-which(is.na(nm)==TRUE);
    nm[indx]<-"Not_Assigned";
    return(nm);
  }
}

GetTaxaFeatSize<- function(microSetObj, taxlvl){

  microSetObj <- .get.microSetObj(microSetObj);
    
  if(taxlvl=="OTU"){
    feat.size <- length(taxa_names(microSetObj$dataSet$norm.phyobj));
  }else{
    taxa_table<-tax_table(microSetObj$dataSet$proc.phyobj);
    data<-merge_phyloslim(microSetObj$dataSet$norm.phyobj,taxa_table);
    feat.size <- length(unique(as.character(tax_table(data)[,taxlvl])));
  }
  return(feat.size);
}

ValidateFeatureName<- function(microSetObj, taxlvl, nm){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  if(taxlvl=="OTU"){
    tax.nms <- taxa_names(microSetObj$dataSet$norm.phyobj);
  }else{
    taxa_table<-tax_table(microSetObj$dataSet$proc.phyobj);
    data<-merge_phyloslim(microSetObj$dataSet$norm.phyobj,taxa_table);
    tax.nms<-unique(as.character(tax_table(data)[,taxlvl]));
  }
  
  if(nm %in% tax.nms){
    if(.on.public.web){
      .set.microSet(microSetObj)  
      return(1);
    }else{
      return(.set.microSet(microSetObj));
    }
  }
  
  if(.on.public.web){
    .set.microSet(microSetObj)  
    return(0);
  }else{
    return(.set.microSet(microSetObj));
  }
}

GetLowerTaxaLvlNm<- function(microSetObj, taxrank){
  microSetObj <- .get.microSetObj(microSetObj);
  indx <- which(colnames(tax_table(microSetObj$dataSet$proc.phyobj))==taxrank);
  rem <- ncol(tax_table(microSetObj$dataSet$proc.phyobj))-indx;
  return(colnames(tax_table(microSetObj$dataSet$proc.phyobj))[indx+1:rem]);
}

GetHighTaxaLvlNm<- function(microSetObj, taxrank){
  microSetObj <- .get.microSetObj(microSetObj);
  if(taxrank=="OTU"){
    return(colnames(tax_table(microSetObj$dataSet$proc.phyobj))[1:length(colnames(tax_table(microSetObj$dataSet$proc.phyobj)))]);
  }else{
    indx <- which(colnames(tax_table(microSetObj$dataSet$proc.phyobj))==taxrank);
    rem <- (indx):1;
    return(colnames(tax_table(microSetObj$dataSet$proc.phyobj))[rev(rem)]);
  }
}

GetSampleGrpUser <- function(microSetObj, class){
  microSetObj <- .get.microSetObj(microSetObj);
  return(levels(get_variable(microSetObj$dataSet$proc.phyobj, class)));
}

#check whether sample variable is continuous and give warning abt it in metagenomeSeq and LEfSe
CheckContSampleVar <- function(microSetObj, variable){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  cls<-as.factor(sample_data(microSetObj$dataSet$norm.phyobj)[[variable]]);
  
  if(length(cls)/length(levels(cls)) < 2){
    current.msg <<-"Sample variable having continuous values. This method is only applicable for variables containing discrete values.";
    if(.on.public.web){
      .set.microSet(microSetObj)  
      return(0);
    }else{
      return(.set.microSet(microSetObj));
    }
  } else {
    if(.on.public.web){
      .set.microSet(microSetObj)  
      return(1);
    }else{
      return(.set.microSet(microSetObj));
    }
  }
}

# save the processed data with class names
SaveData <- function(microSetObj, anal.type){
  microSetObj <- .get.microSetObj(microSetObj);
  if(anal.type == "16S" || anal.type == "metageno"){
    if(!is.null(dataSet$data.orig)){
      write.csv(microSetObj$dataSet$data.orig, file="data_original.csv");
    }
    if(!is.null(dataSet$filt.data)){
      write.csv(microSetObj$dataSet$filt.data, file="data_filtered.csv");
    }
    if(!is.null(otu_table(microSetObj$dataSet$norm.phyobj,as.matrix()))){
      write.csv(otu_table(microSetObj$dataSet$norm.phyobj,as.matrix()), file="data_normalized.csv");
    }
  }else if(anal.type == "16S_ref"){
    write.csv(otu_table(userrefdata,as.matrix()), file="merge_otutable.csv");
    write.csv(as.matrix(tax_table(userrefdata)), file="merge_taxtable.csv");
    write.csv(as.data.frame(sample_data(userrefdata)), file="merge_sampletable.csv");
  }
  return(.set.microSet(microSetObj));
}

ReadSampleTable<- function(microSetObj, dataName) {
  
  microSetObj <- .get.microSetObj(microSetObj);
  
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
  
  microSetObj$dataSet$sample_data <- data.frame(mydata);
  current.msg <<- paste(na.msg, "The sample data contains a total of ", nrow(mydata), "samples and  ", ncol(mydata), " sample variables.", collapse=" ");
  microSetObj$dataSet$smpl.msg <- current.msg;
  
  if(.on.public.web){
    .set.microSet(microSetObj)  
    return(1);
  }else{
    return(.set.microSet(microSetObj));
  }
}

ReadTreeFile<- function(dataName) {
  suppressMessages(library("phyloslimR"));
  msg <- NULL;
  tree <- phyloslimR::read_tree(dataName);
  saveRDS(tree, "tree.RDS");
  return(1)
}
