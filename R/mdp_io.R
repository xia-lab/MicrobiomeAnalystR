
# from clicking a interactive PCoA 3D to generate pie chart
SetCurrentSelectedTaxLevel<-function(taxLvl){
  current.selected.tax <<- taxLvl;
}

#######################################
########### I/O for 16S data ##########
########### Used by MDP & PPD #########
#######################################

#'Main function to read 16S data
#'@description This is the main function read in the 16S data
#'into the mbSetObj.
# ismetafile: whether meta-data is given (for BIOM format only)
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
Read16SAbundData <- function(mbSetObj, dataName, format, taxa_type, ismetafile, is.normalized) {

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(format=="text"){
    if(.on.public.web){
      if(Read16STabData(mbSetObj, dataName)==0){
        return(0);
      }
    }else{
      # offline
      mbSetObj <- Read16STabData(mbSetObj, dataName)
      if(!mbSetObj$dataSet$read){
        return(0);
      }
    }
    
  }else if(format=="biom"){
    if(Read16SBiomData(mbSetObj, dataName, taxa_type, ismetafile)==0){
      return(0);
    } 
  }
  
  #need to refetch mbSetObj
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data.orig <- mbSetObj$dataSet$data.orig;
  
  if(mode(data.orig)=="character"){
    AddErrMsg(paste("Errors in parsing your data as numerical - possible reason: comma as decimal separator?"));
    return(0);
  }
  if(mbSetObj$module.type=="mmp"){
   mbSetObj$micDataType = "otu"
  }else{
  mbSetObj$micDataType = "na"
}

  #storing whole taxonomy label(used for PICRUST & Tax4Fun; reference data mapping)
  if(is.null(mbSetObj$dataSet$comp_taxnm)){
  mbSetObj$dataSet$comp_taxnm <- rownames(data.orig)
  }
  current.msg <<- paste("A total of",ncol(data.orig) ,"samples and ", nrow(data.orig), "features or taxa are present.");
  mbSetObj$dataSet$read.msg <- current.msg;
  mbSetObj$dataSet$data.type <- format;
  mbSetObj$dataSet$taxa.type <- taxa_type;
  mbSetObj$dataSet$is.normalized <- is.normalized;

  return(.set.mbSetObj(mbSetObj));

}

#'Function to read 16S data
#'@description This is the main function read in the 16S data in txt format
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

Read16STabData <- function(mbSetObj, dataName) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  msg <- NULL;
  mydata <- .readDataTable(dataName);

  if(any(c(any(is.na(mydata)), class(mydata) == "try-error"))){
    AddErrMsg("Failed to read in the OTU abundance data! Please make sure the data is in the right format and do not have empty cells or NA.");
    return(0);
  }
 
  # note dateSet$sample_data already set
  # getting NAME label
  if(mbSetObj$module.type=="mmp" |mbSetObj$module.type=="meta" ){
     sam.nm <- substr(colnames(mydata[1]),1,5);
    sam.nm <- tolower(sam.nm);
    sam.inx <- grep("^#name",sam.nm);
    
    if(length(sam.inx) == 0){
      
      sam.nm <- tolower(colnames(mydata[1]));
      if(sam.nm %in% c("#phylum","#class","#order","#family","#genus","#species") ){
        mbSetObj$dataSet$sglTax <- gsub("#","",sam.nm)
      }else{
  
        AddErrMsg("No labels ID found in your otu/asv table!");
        return(0);
      }
    
    } 
  
  }else{
    
    sam.nm <- substr(colnames(mydata[1]),1,5);
    sam.nm <- tolower(sam.nm);
    sam.inx <- grep("^#name",sam.nm);
    
    if(length(sam.inx) == 0){
      AddErrMsg("No labels #NAME found in your data!");
      return(0);
    }  
  }
  
  smpl_nm <- colnames(mydata[-1]);

  mydata <- .to.numeric.mat(mydata);

  # empty cell or NA cannot be tolerated in metadata
  na.inx  <- is.na(mydata);
  
  if(sum(na.inx) > 0){
    mydata[na.inx] <- "Unknown";
    msg <- c(msg, paste("A total of", sum(na.inx), "empty or NA values were replaced by zeros."));
  }
  
  current.msg <<- paste(msg, collapse="; ");
  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$smpl_nm <- smpl_nm;
  mbSetObj$dataSet$data.orig <- data.matrix(mydata);
  
  mbSetObj$dataSet$read <- TRUE
  return(.set.mbSetObj(mbSetObj))
}

#'Function to read 16S data in biom format
#'@description This function reads in the 16S data from biom format
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
Read16SBiomData <- function(mbSetObj, dataName, taxa_type, ismetadata){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  suppressMessages(library(biomformat));

  msg <- NULL;

  if(any(c(taxa_type=="Greengenes", taxa_type=="GreengenesID"))){
    
    mydata <- tryCatch( 
      
      { # try
        import_biom(dataName,parseFunction = parse_taxonomy_greengenes)
      },
      
      error = function(error_cond){
        
        test <- read_biom(dataName)
        test <- is.null(observation_metadata(test))
        
        if(test){
          AddErrMsg("Biom file does not contain taxonomic information!");
          return(0);
        }
      })
    
  }else if(taxa_type == "SILVA"){
    
    mydata <- tryCatch( 
      
      { # try
        import_biom(dataName,parseFunction = parse_taxonomy_silva_128);
      },
      
      error = function(error_cond){
        
        test <- read_biom(dataName)
        test <- is.null(observation_metadata(test))
        
        if(test){
          AddErrMsg("Biom file does not contain taxonomic information!");
          return(0);
        }
      })

  }else{
    
    mydata = tryCatch( 
      
      { # try
        import_biom(dataName,parseFunction = parse_taxonomy_default);
      },
      
      error = function(error_cond){
        
        test <- read_biom(dataName)
        test <- is.null(observation_metadata(test))
        
        if(test){
          AddErrMsg("Biom file does not contain taxonomic information!");
          return(0);
        }
      }
    )
  }

  if(class(mydata) == "numeric"){
    AddErrMsg("Verify that the uploaded biom file has both abundance and taxonomic information!");
    return(0)
  }
    
  #reading data separately(abundance;taxonomy is must and sample data can be uploaded separately)
  otu.dat <- otu_table(mydata,taxa_are_rows = TRUE,errorIfNULL = FALSE);
    
  if(length(otu.dat)==0){
    AddErrMsg("File does not contain 16S abundance data");
    return(0);
  }
    
  otu.dat <- as.matrix(otu.dat);
  msg <- c(msg, "Abundance data present.");

  #taxonomy table
  taxa_table <- tax_table(mydata,errorIfNULL = FALSE);
    
  if(length(taxa_table)==0){
    AddErrMsg("File does not contain taxonomy data.");
    return(0);
  }
    
  taxa_table<-as.matrix(taxa_table);

  if(ncol(taxa_table)==1){
    AddErrMsg("Taxonomy table in .biom file is invalid!")
    return(0)
  }

  msg <- c(msg, "Taxonomy file is detected.");

  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$data.orig <- otu.dat;
  mbSetObj$dataSet$taxa_table <- taxa_table;
    
  #sample data
  if(ismetadata=="T"){
    sample_data <- sample_data(mydata,errorIfNULL = FALSE);
    sample_data <- as.data.frame(sample_data,check.names=FALSE);
    if(length(sample_data)==0){
      AddErrMsg("Metadata file not detected in your biom file. Please upload metadata file separately.");
      return(0);
    }
        
    mbSetObj$dataSet$sample_data <- sample_data;
    msg <- c(msg, "Metadata file is detected in your biom file.");
    
    mbSetObj$dataSet$group_names <- colnames(sample_data)
    
    disc.inx <- GetDiscreteInx(sample_data);
    if(sum(disc.inx) == 0){ # all class labels are unique! 
       
        AddErrMsg("Metadata Table: make sure there is at least one column contains experimental design for group comparisons (i.e., the primary metadata), with each group contains at least 3 replicate. No unique values are allowed in the primary metadata column.");
        msg <- c(msg, "It seems that some of your metadata values are unique! MicrobiomeAnalyst requires some biological replicates for robust analysis");
        mbSetObj$poor.replicate <- TRUE;
        mbSetObj$dataSet$sample_data <- sample_data
       return(0)
    }else{
        if(sum(disc.inx) == length(disc.inx)){
            msg <- c(msg, "All metadata columns are OK.");
        }else{
            bad.meta<- paste(names(disc.inx)[!disc.inx], collapse="; ");
            msg <- c(msg, paste0("The following metadata columns are excluded: <font style=\"color:red\"><b>", bad.meta, "</b></font>"));
        }
        mbSetObj$dataSet$meta_info$disc.inx <- disc.inx;
        mbSetObj$dataSet$sample_data <- sample_data[,disc.inx, drop=FALSE];
      
        cont.inx <- GetNumbericalInx(sample_data);
        cont.inx <- !disc.inx & cont.inx; # discrete is first
        mbSetObj$dataSet$meta_info$cont.inx <- cont.inx;
      
        if(sum(cont.inx)>0){
            # make sure the discrete data is on the left side
            mbSetObj$dataSet$sample_data <- cbind(mbSetObj$dataSet$sample_data, sample_data[,cont.inx, drop=FALSE]);
        }
        }
        msg <- paste(msg, "The sample data contains a total of ", nrow(sample_data), "samples and  ", ncol(sample_data), " sample variables.", collapse=" ");
    }
    current.msg <<- paste(msg, collapse="; ");
    mbSetObj$dataSet$smpl.msg <- current.msg;
    mbSetObj$dataSet$read <- TRUE;
    return(.set.mbSetObj(mbSetObj));
  
}

#'Function to read 16S data in mothur format
#'@description This function reads in the 16S data from mothur format
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
#'@import data.table

ReadMothurData<-function(mbSetObj, dataName, taxdataNm, taxa_type){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  suppressMessages(library(data.table));
  
  msg <- NULL;

  indx <- grep(".shared$", dataName);
    
  if(length(indx)!=1){
    AddErrMsg("Data format error. Make sure that abundance file ends with (.shared) extension");
    return(0);
  }
    
  indx2 <- grep(".taxonomy$", taxdataNm);
    
  if(length(indx2)!=1){
    AddErrMsg("Data format error. Make sure that taxonomy file ends with (.taxonomy) extension");
    return(0);
  }

  #reading mothur data using phyloseq seperately.
  mydata <- import_mothur(mothur_shared_file = dataName);
    
  if(length(mydata)==0){
    AddErrMsg("Data format error. Make sure the file is not empty and is in mothur shared file format(.shared).");
    return(0);
  }

  msg <- c(msg, "OTU abundance table is present.");

  #taxonomy table
    
  if(any(c(taxa_type=="Greengenes", taxa_type=="GreengenesID"))){
    tax_data <- import_mothur(mothur_constaxonomy_file = taxdataNm,parseFunction = parse_taxonomy_greengenes);
  }else{
    tax_data <- import_mothur(mothur_constaxonomy_file = taxdataNm,parseFunction = parse_taxonomy_default);
  }

  #getting orignal label from .taxonomy file (Tax4Fun)
  if(nrow(tax_data)!=nrow(mydata)){
    AddErrMsg("Please make sure that number of features and their row names are same in OTU abundance table and taxonomy table.");
    return(0);
  }
    
  taxa_names(tax_data) <- taxa_names(mydata);
  tax_tabdata <- .readDataTable(taxdataNm);
  comp_taxnm <- tax_tabdata$Taxonomy;
  taxa_table <- tax_table(tax_data,errorIfNULL = FALSE);
    
  if(length(taxa_table)==0){
    AddErrMsg("Make sure the file is not empty and is in Mothur taxonomy file format (should have header).");
    return(0);
  }
    
  msg <- c(msg,paste("Taxonomy file is also detected."));

  #storing object for other purpose
  taxa_table <- as.matrix(taxa_table);

  #sanity check: features name of taxonomy table should match feature names in OTU abundance table, only features in filtered OTU tables are selected additional are removed.
  indx<-match(rownames(mydata), rownames(taxa_table));
    
  if(all(is.na(indx))){
    AddErrMsg("Please make sure that features name of taxonomy table should match feature names in OTU abundance table.");
    return(0);
  }

  msg <- c(msg,paste("A total of ",ncol(mydata) , " samples and ", nrow(mydata), "features were found."));
  current.msg <<- paste(msg, collapse="; ");
    
  taxa_table <- taxa_table[indx,];
  rownames(taxa_table) <- rownames(mydata);
  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$data.orig <- mydata;
  mbSetObj$dataSet$taxa_table <- taxa_table;
  mbSetObj$dataSet$comp_taxnm <- comp_taxnm;
  mbSetObj$dataSet$read.msg <- current.msg;
  mbSetObj$dataSet$data.type <- "mothur";
  mbSetObj$dataSet$taxa.type <- taxa_type;
  
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to read 16S taxonomy table from txt format
#'@description This function reads in the 16S taxonomy table from txt format
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
Read16STaxaTable <- function(mbSetObj, dataName) {
  mydata <- .readDataTable(dataName);
  if(any(c(is.null(mydata), class(mydata) == "try-error"))){
    AddErrMsg("Failed to read in the taxonomic data! Please make sure the data is in the right format.");
    return(0);
  }
  # look for #TAXONOMY, store in a list
  sam.nm <- substr(colnames(mydata[1]),1,9);
  sam.nm <- tolower(sam.nm);
  sam.inx <- grep("^#taxonomy",sam.nm);
    
  if(length(sam.inx) > 0){
    tax_nm <- mydata[,1];
    tax_rank <- colnames(mydata[-1]);
  }else{
    AddErrMsg("Failed to read in the taxonomic data! Please make sure the data is in the right format.");
    return(0);
  }

  # converting to character matrix as duplicate row names not allowed in data frame.
  mydata <- as.matrix(mydata[,-1]);
  rownames(mydata) <- tax_nm;
  colnames(mydata) <- tax_rank;
  current.msg <<- paste("Taxonomy file has total of ",nrow(mydata),"features and",ncol(mydata), "taxonomic rank. Additional features which are not present in abundance table has been removed if present.");

  mbSetObj <- .get.mbSetObj(mbSetObj);   
 

if("Species" %in% colnames( mydata)){# for mapping to the database in mmp module and save the taxanm for function pediction
  sps = data.frame(mydata[,c('Genus','Species')])
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
  mydata[,'Species'] <- sps$Species
  
}

 mbSetObj$dataSet$comp_taxnm <- apply(mydata,1,function(x) paste(x[!is.na(x)],collapse = ";"))
  mbSetObj$dataSet$taxa_table <- mydata;

  return(.set.mbSetObj(mbSetObj));
}

#'Main function to read metabolomics table
#'@description This is the main function read in the metabolomics data
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'@param metType whether a peak or metabolite intensity table
#'@param idType One of compound name,KEGG or HMDB ID. Only applicable when the metType is metabolite. 
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

ReadMetabolicTable <- function(mbSetObj, dataName, metType, idType) {

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet$metabolomics <- list()

  current.msg <<- NULL;
  mydata <- .readDataTable(dataName); 

  if( class(mydata) == "try-error"){
    AddErrMsg("Failed to read in the metabolomics abundance table! Please make sure the data is in the right format.");
    return(0);
  }
  
 # getting NAME label
  sam.nm <- substr(colnames(mydata[1]),1,5);
  sam.nm <- tolower(sam.nm);
  sam.inx <- grep("^#name",sam.nm);
  
  if(length(sam.inx) == 0){
    
    AddErrMsg("No labels #NAME found in your data!");
    return(0);
  }   
  
  #smpl_nm <- colnames(mydata[-1]);
  
  mydata <- .to.numeric.mat(mydata); 
  row.nas <- apply(is.na(mydata)|is.null(mydata), 1, sum);
  good.inx<- row.nas/ncol(mydata) < 0.5;
  if(sum(!good.inx) > 0){
    mydata <- mydata[good.inx,];
    current.msg <<- c(current.msg, paste("removed ", sum(!good.inx), " features with over 50% missing values"));
  }
 
  minVal <- min(mydata, na.rm=T);
  na.inx <- is.na(mydata);
  if(sum(na.inx) > 0){
    mydata[na.inx] <- minVal/2;
    current.msg <<- c(current.msg, "the remaining", sum(na.inx), "missing variables were replaced with data min");
  }

  mbSetObj$dataSet$metabolomics$name <- basename(dataName);
  #mbSetObj$dataSet$metabolomics$smpl_nm <- smpl_nm;
  mbSetObj$dataSet$metabolomics$data.orig <- data.matrix(mydata);
  mbSetObj$dataSet$metabolomics$comp_metnm <- rownames(mydata);
  mbSetObj$dataSet$metabolomics$feature.type <- metType;
  mbSetObj$dataSet$metabolomics$id.type <- idType;
  mbSetObj$dataSet$metabolomics$read <- TRUE;
  
  current.msg <<- c(current.msg,paste("A total of",ncol(mydata) ,"metabolomics samples and ", nrow(mydata), metType," are present."));
  mbSetObj$dataSet$metabolomics$read.msg <- current.msg
 

  return(.set.mbSetObj(mbSetObj));
  
}


ReadPeakList <- function(mbSetObj=NA, dataName = NA, rankOpt="pval",rtOpt,mode,instrumentOpt) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  mydata <- data.frame(.readDataTable(dataName));

  if(any(c(is.null(mydata), class(mydata) == "try-error"))){
    AddErrMsg("Failed to read in the peak list! Please make sure the data is in the right format.");
    return(0);
  }
  
  user_cols <- gsub("[^[:alnum:]]", "", colnames(mydata));
  mummi.cols <- c("m.z", "p.value", "t.score", "r.t","logfc")

  # next check what column names are there
  hit <- "mz" %in% user_cols;
  
  if(sum(hit) < 1){
    AddErrMsg("Missing information, data must contain a 'm.z' column!");
    return(0);
  }
  
  
  if(current.proc$mumRT  & ! "rt" %in% user_cols){
     AddErrMsg("Missing information, r.t was not found in your data!");
    return(0);
  }

  if(current.proc$mumRT  & ! "rt" %in% user_cols){
    AddErrMsg("Missing information, r.t was not found in your data!");
    return(0);
  }
  
  if(current.proc$mode =="mixed"){
    if(! mode %in% user_cols){
      AddErrMsg("Missing information, mode was not found in your data which is mandatory for mixed mode!");
      return(0);
    }else{
      current.proc$pos_inx = mydata$mode == "positive"
      mydata <- subset(mydata, select=-mode)
      user_cols <- gsub("[^[:alnum:]]", "", colnames(mydata))
    }
   
  }

  
  if(rankOpt == "pvalue"){
    if( ! "pvalue" %in% user_cols){
      AddErrMsg("Missing information, p.values were not found in your data!");
      return(0);
    }else{
      current.proc$rankopt = "p.value"
    }
  }
  
  if(rankOpt == "tscore"){
    if( ! "tscore" %in% user_cols){
      AddErrMsg("Missing information, t-scores were not found in your data!");
      return(0);
    }else{
      current.proc$rankopt = "tscore"
    }
  }
  if(rankOpt == "fc"){
    if( ! "logfc" %in% user_cols){
      AddErrMsg("Missing information, fold changes were not found in your data!");
      return(0);
    }else{
      current.proc$rankopt = "tscore"
    }
  }
 
  mbSetObj$dataSet$metabolomics$name <- basename(dataName);
  mbSetObj$dataSet$metabolomics$data.orig <- data.matrix(mydata);
  mbSetObj$dataSet$metabolomics$feature.type <- "peak";
  mbSetObj$dataSet$metabolomics$read <- TRUE

return(.set.mbSetObj(mbSetObj));
}

#'Function to create a summary of a sample using a piechart at different tax level.
#'@description This function creates a piechart summary of a sample at a specific
#'taxonomic level. Used by PPD and MDP modules.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotSelectedSample <-function(mbSetObj, imgNm, smplID, OtuIdType, rel_perct, format="png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(current.selected.tax == "NA"){
    txlvl <- "Phylum";
  }else{
    txlvl <- current.selected.tax;
  }
 
  if(mbSetObj$module.type=="mdp"){
    data <- mbSetObj$dataSet$norm.phyobj;
   }else{
    data <- userrefdata;
  }
    
  #prune data sample wise
  data_s <- prune_samples(smplID,data);
    
  #creating a new pruned phyloseq object
  data_newph <- data_s;
  data_tax <- tax_table(data_newph);
  piedata_new <- GetDataForPie(data_newph,data_tax,txlvl,OtuIdType,feat_cnt);
  h <- 460;
  a <- length(unique(data_tax[,txlvl]));
  x.colors <- rep(col_vector,length.out=a);
  rowNum <- ceiling(a/3);
  myH <- rowNum*18 + h;
  imgNm = paste(imgNm,".",format, sep="");
  Cairo::Cairo(file=imgNm,width=420, height=myH, type=format, bg="white",dpi=dpi);
  piedata_rel <- transform(transform(piedata_new, value=value/sum(value)));
  ind <- which(piedata_rel[,"value"]>rel_perct);
  ind1 <- which(piedata_rel[,"value"]<rel_perct);
    
  if(length(ind)==0){
    AddErrMsg("All features have lower relative abundance than given minimum abundance. Please lower the cut off for relative abundance.");
    return(0);
  }
    
  if(length(ind1)>0){
    levels(piedata_rel$variable)[ind1] <- "Others";
    piedata_rel <- aggregate(value~variable, piedata_rel, FUN=sum);
  }
    
  ind_zero <- which(piedata_rel[,"value"]==0);
    
  if(length(ind_zero)>0){
    piedata_rel <- piedata_rel[-ind_zero,];
  }
    
  box = ggplot(piedata_rel,
        aes(x="", y = value, fill=variable)) +geom_bar(width = 1, stat = "identity") + theme_bw() +coord_polar(theta = "y") +
        geom_text(aes(x=1.7, label = scales::percent(value)), check_overlap = T,size=3, position = position_stack(vjust = 0.5)) +
        theme(axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank())+
        guides(fill=guide_legend(ncol=3)) + labs(x="", y="",fill ="") +scale_fill_manual(values=c(x.colors))+
        ggtitle(paste(smplID, " [", txlvl, "]", sep="")) + theme(legend.position="bottom", plot.title = element_text(hjust=0.5, face="bold"));
  print(box);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj));
}

#'Utility function to exract data for piechart.
#'@description This function extracts necessary data to create a piechart.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import stringr
GetDataForPie<-function(data_n, datataxa, txlvl, OtuIdType, feat_cnt){
  
  #using reduce names
  data_new <- otu_table(data_n);
  data_new <- data.frame(data_new,check.names=FALSE);
  data_new <- t(data_new);

  if(txlvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    taxa_nm <- as.data.frame(datataxa[,txlvl],check.names=FALSE);
  }
    
  if(OtuIdType=="GreengenesID"){
    
    suppressMessages(library(stringr));
    
    nm<-c("k__","p__","c__","o__","f__","g__","s__");
    first.10 <- substr(taxa_nm[,1], start=1, stop=3);
    mat<-match(first.10,nm);
    ind<-which(!is.na(mat));
    
    if(length(ind)!=0){
      y<-str_sub(taxa_nm[ind,1], 4, -1);
      indx <- which(y=="");
      y[indx]<-"NA";
      taxa_nm[,1]<-y;
    }
  }
    
  taxa_nm<-as.matrix(taxa_nm);
  y <- which(is.na(taxa_nm)==TRUE);
  #converting NA values to unassigned
  taxa_nm[y] <- "Not_Assigned";
  colnames(data_new)<-taxa_nm[,1];
  nms <-colnames(data_new);
  data_new<-data_new %*% sapply(unique(nms),"==",nms);
  data_new<-data.frame(data_new,check.names=FALSE);
  data_new$step<-factor(rownames(data_new));
  data_new<-reshape2::melt(data_new,id='step');
  data_new$step<-as.numeric(data_new$step);
  piedata_new<-data_new[-1];
  piedata_new<-aggregate(. ~variable , data=piedata_new, FUN=sum);
  piedata_new<-piedata_new;
  return(piedata_new);
}

################################# 
######## Utility Function #######
################################# 

# utilty function to remove from, within, leading and trailing spaces
# also removes /
ClearFactorStrings<-function(cls.nm, query){
  # remove leading and trailing space
  query<- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);

  # kill multiple white space
  query <- gsub(" +","_",query);
  # remove non alphabets and non numbers
  query <- gsub("[^[:alnum:] ]", "_", query);

  # test all numbers (i.e. Time points)
  chars <- substr(query, 0, 1);
  num.inx<- chars >= '0' & chars <= '9';
    
  if(all(num.inx)){
    query = as.numeric(query);
    nquery <- paste(cls.nm, query, sep="_");
    query <- factor(nquery, levels=paste(cls.nm, sort(unique(query)), sep="_"));
  }else{
    query[num.inx] <- paste(cls.nm, query[num.inx], sep="_");
    query <- factor(query);
  }
  return (query);
}

# deal with issues with SILVA taxonomy parsing for import_biom
# https://gist.github.com/grabear/018e86413b19b62a6bb8e72a9adba349
parse_taxonomy_silva_128 <- function(char.vec){
  # Use default to assign names to elements in case problem with greengenes prefix
  char.vec = parse_taxonomy_default(char.vec)
  # Check for unassigned taxa
  if (char.vec["Rank1"] == "Unassigned") {
    char.vec <- c(Rank1="D_0__Unassigned", Rank2="D_1__Unassigned", Rank3="D_2__Unassigned", Rank4="D_3__Unassigned",
                  Rank5="D_4__Unassigned", Rank6="D_5__Unassigned", Rank7="D_6__Unassigned")
  }
  # Define the meaning of each prefix according to SILVA taxonomy
  Tranks = c(D_0="Kingdom", D_1="Phylum", D_2="Class", D_3="Order", D_4="Family", D_5="Genus", D_6="Species")
  # Check for prefix using regexp, warn if there were none. trim indices, ti
  ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
  if( length(ti) == 0L ){
    warning(
      "No silva prefixes were found. \n",
      "Consider using parse_taxonomy_delfault() instead if true for all OTUs. \n",
      "Dummy ranks may be included among taxonomic ranks now."
    )
    # Will want to return without further modifying char.vec
    taxvec = char.vec
    # Replace names of taxvec according to prefix, if any present...
  } else {
    # Format character vectors for Ambiguous taxa
    if( length(ti) < 7 ){
      for (key in names(char.vec)) {
        if ( char.vec[key] == "Ambiguous_taxa" ) {
          tax_no <- (as.numeric(substr(key, 5, 5)) - 1)
          char.vec[key] = sprintf("D_%s__Ambiguous_taxa", tax_no)
        }
      }
      # Reset the trimmed indicies if Ambiguous taxa
      ti = grep("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", char.vec)
    }
    # Remove prefix using sub-"" regexp, call result taxvec
    taxvec = gsub("[[:alpha:]]\\_[[:digit:]]{1}\\_\\_", "", char.vec)
    # Define the ranks that will be replaced
    repranks = Tranks[substr(char.vec[ti], 1, 3)]
    # Replace, being sure to avoid prefixes notK present in Tranks
    names(taxvec)[ti[!is.na(repranks)]] = repranks[!is.na(repranks)]
  }
  return(taxvec)
}
