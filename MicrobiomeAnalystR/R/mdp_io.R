#######################################
########### I/O for 16S data ##########
########### Used by MDP & PPD #########
#######################################

#'Main function to read 16S data
#'@description This is the main function read in the 16S data
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
Read16SAbundData <- function(mbSetObj, dataName, type, taxalabel, taxa_type,
                             ismetafile, module.type) {
  
  if(.on.public.web){
    load_phyloseq();
  }
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(type=="text"){
    
    if(.on.public.web){
      if(!Read16STabData(mbSetObj, dataName, type, ismetafile)){
        return(0);
      }
    }else{
      # offline
      mbSetObj <- Read16STabData(mbSetObj, dataName, type, ismetafile)
      if(!mbSetObj$dataSet$read){
        return(0);
      }
    }
    
  }else if(type=="biom"){
    
    if(.on.public.web){
      if(!Read16SBiomData(mbSetObj, dataName, type, taxa_type, ismetafile)){
        return(0);
      } 
    }else{
      # offline
      mbSetObj <- Read16SBiomData(mbSetObj, dataName, type, taxa_type, ismetafile)
      if(!mbSetObj$dataSet$read){
        return(0);
      }
    }
  }
  
  #need to refetch mbSetObj
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data.orig <- mbSetObj$dataSet$data.orig;
  
  if(mode(data.orig)=="character"){
    current.msg <<- paste("Errors in parsing your data as numerical - possible reason: comma as decimal separator?");
    return(0);
  }
  
  #storing whole taxonomy label(used for PICRUST & Tax4Fun; reference data mapping)
  mbSetObj$dataSet$comp_taxnm <- rownames(data.orig);
  current.msg <<- paste("A total of",ncol(data.orig) ,"samples and ", nrow(data.orig), "features or taxa are present.");
  mbSetObj$dataSet$read.msg <- current.msg;
  mbSetObj$dataSet$data.type <- type;
  mbSetObj$dataSet$taxa.type <- taxa_type;
  mbSetObj$dataSet$module.type <- module.type;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}

#'Function to read 16S data
#'@description This is the main function read in the 16S data in txt format
#'into the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

Read16STabData <- function(mbSetObj, dataName, type, ismetafile) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  msg <- NULL;
  mydata <- .readDataTable(dataName);

  if(any(is.na(mydata)) || class(mydata) == "try-error"){
    current.msg <<- "Failed to read in the OTU abundance data! Please make sure the data is in the right format and do not have empty cells or NA.";
    return(FALSE);
  }

  # if metadata file present seperately
  if(ismetafile=="T"){
    # note dateSet$sample_data already set
    # getting NAME label
    sam.nm <- substr(colnames(mydata[1]),1,5);
    sam.nm <- tolower(sam.nm);
    sam.inx <- grep("^#name",sam.nm);
        
    if(length(sam.inx) == 0){
      current.msg <<- "No labels #NAME found in your data!";
      return(FALSE);
    }
        
    smpl_nm <- colnames(mydata[-1]);
    
  }else{
    # look for #CLASS,have class labels, store in a list
    meta.info <- list();
    #getting Group label
    sam.nm <- substr(mydata[1,1],1,6);
    sam.nm <- tolower(sam.nm);
    cls.inx <- grep("^#class",sam.nm);
        
    if(length(cls.inx) > 0){
      for(i in 1:length(cls.inx)){
        inx <- cls.inx[i];
        cls.nm <- substring(mydata[inx, 1],2); # discard the first char #
        if(nchar(cls.nm) > 6){
          cls.nm <- substring(cls.nm, 7); # remove class
        }
        cls.lbls <- mydata[inx, -1];
        # test NA
        na.inx <- is.na(cls.lbls);
        cls.lbls[na.inx] <- "NA";
        cls.lbls <- ClearFactorStrings(cls.nm, cls.lbls);
        meta.info[[cls.nm]] <- cls.lbls;
      }
    }else{
      current.msg <<- "No metadata labels #CLASS found in your data!";
      return(FALSE);
    }
        
    #creating phyloseq(sample_data)object
    sample_data <- as.data.frame(cls.lbls);
    names(sample_data) <- "CLASS";
    rownames(sample_data) <- smpl_nm <- colnames(mydata)[-1];
    mbSetObj$dataSet$sample_data <- sample_data; # set up sample_data
    mydata<-mydata[-1, ];
  }

  # remove extra comments if any
  dat.nms <- mydata[,1];
  comment.inx <- grep("^#",dat.nms);
    
  if(sum(comment.inx) > 0){
    mydata < -mydata[-comment.inx, ];
    msg <- c(msg, paste("A total of", sum(comment.inx), "comment rows were removed from OTU table"));
    dat.nms <- mydata[,1];
  }

  mydata <- data.matrix(mydata[,-1]);
  rownames(mydata) <- dat.nms;

  # empty cell or NA cannot be tolerated in metadata
  na.inx  <- is.na(mydata);
    
  if(sum(na.inx) > 0){
    mydata[na.inx] <- "Unknown";
    msg <- c(msg, paste("A total of", sum(na.inx), "empty or NA values were replaced by zeros."));
  }

  msg <<- paste(msg, collapse="; ");

  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$smpl_nm <- smpl_nm;
  mbSetObj$dataSet$data.orig <- data.matrix(mydata);
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(TRUE);
  }else{
    mbSetObj$dataSet$read <- TRUE
    return(.set.mbSetObj(mbSetObj))
  }
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
Read16SBiomData <- function(mbSetObj, dataName, type, taxa_type, ismetadata){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_phyloseq();
    load_biomformat();
  }
  
  msg <- NULL;

  if(taxa_type=="Greengenes"||taxa_type=="GreengenesID"){
    mydata <- import_biom(dataName,parseFunction = parse_taxonomy_greengenes);
  }else{
    mydata <- import_biom(dataName,parseFunction = parse_taxonomy_default);
  }
    
  #reading data separately(abundance;taxonomy is must and sample data can be uploaded separately)
  otu.dat <- otu_table(mydata,taxa_are_rows = TRUE,errorIfNULL = FALSE);
    
  if(length(otu.dat)==0){
    current.msg <<-"File does not contain 16S abundance data";
    return(FALSE);
  }
    
  otu.dat <- as.matrix(otu.dat);
  msg <- c(msg, "Abundance data present.");

  #taxonomy table
  taxa_table <- tax_table(mydata,errorIfNULL = FALSE);
    
  if(length(taxa_table)==0){
    current.msg <<-"File does not contain taxonomy data.";
    return(FALSE);
  }
    
  taxa_table<-as.matrix(taxa_table);
  msg <- c(msg, "Taxonomy file is detected.");

  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$data.orig <- otu.dat;
  mbSetObj$dataSet$taxa_table <- taxa_table;
    
  #sample data
  if(ismetadata=="T"){
    sample_data <- sample_data(mydata,errorIfNULL = FALSE);
    sample_data <- as.data.frame(sample_data);
        
    if(length(sample_data)==0){
      current.msg <<- "Metadata file not detected in your biom file. Please upload metadata file seperately.";
      ismetafile<-"F";
      return(FALSE);
    }
        
    mbSetObj$dataSet$sample_data <- as.data.frame(sample_data);
    msg <- c(msg, "Metadata file is detected in your biom file.");
  }
    
  current.msg <<- paste(msg, collapse=". ");
  mbSetObj$dataSet$taxa.type <- taxa_type;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(TRUE);
  }else{
    mbSetObj$dataSet$read <- TRUE
    return(.set.mbSetObj(mbSetObj))
  }
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

ReadMothurData<-function(mbSetObj, dataName, taxdataNm, taxa_type, module.type){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_phyloseq();
    load_datatable();
  }
  
  msg <- NULL;

  indx <- grep(".shared$", dataName);
    
  if(length(indx)!=1){
    current.msg <<-"Data format error. Make sure that abundance file ends with (.shared) extension";
    return(0);
  }
    
  indx2 <- grep(".taxonomy$", taxdataNm);
    
  if(length(indx2)!=1){
    current.msg <<-"Data format error. Make sure that taxonomy file ends with (.taxonomy) extension";
    return(0);
  }

  #reading mothur data using phyloseq seperately.
  mydata <- import_mothur(mothur_shared_file = dataName);
    
  if(length(mydata)==0){
    current.msg <<-"Data format error. Make sure the file is not empty and is in mothur shared file format(.shared).";
    return(0);
  }

  msg <- c(msg, "OTU abundance table is present.");

  #taxonomy table
    
  if(taxa_type=="Greengenes"||taxa_type=="GreengenesID"){
    tax_data <- import_mothur(mothur_constaxonomy_file = taxdataNm,parseFunction = parse_taxonomy_greengenes);
  }else{
    tax_data <- import_mothur(mothur_constaxonomy_file = taxdataNm,parseFunction = parse_taxonomy_default);
  }

  #getting orignal label from .taxonomy file (Tax4Fun)
  if(nrow(tax_data)!=nrow(mydata)){
    current.msg <<-"Please make sure that number of features and their row names are same in OTU abundance table and taxonomy table.";
    return(0);
  }
    
  taxa_names(tax_data) <- taxa_names(mydata);
  tax_tabdata <- fread(taxdataNm, header=TRUE, check.names=FALSE, data.table=FALSE);
  comp_taxnm <- tax_tabdata$Taxonomy;
  taxa_table <- tax_table(tax_data,errorIfNULL = FALSE);
    
  if(length(taxa_table)==0){
    current.msg <- "Make sure the file is not empty and is in Mothur taxonomy file format (should have header).";
    return(0);
  }
    
  msg <- c(msg,paste("Taxonomy file is also detected."));

  #storing object for other purpose
  taxa_table <- as.matrix(taxa_table);

  #sanity check: features name of taxonomy table should match feature names in OTU abundance table, only features in filtered OTU tables are selected additional are removed.
  indx<-match(rownames(mydata), rownames(taxa_table));
    
  if(all(is.na(indx))){
    current.msg <<-"Please make sure that features name of taxonomy table should match feature names in OTU abundance table.";
    return(0);
  }
    
  taxa_table <- taxa_table[indx,];
  rownames(taxa_table) <- rownames(mydata);

  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$data.orig <- mydata;
  mbSetObj$dataSet$taxa_table <- taxa_table;
  mbSetObj$dataSet$comp_taxnm <- comp_taxnm;
  msg <- c(msg,paste("A total of ",ncol(mydata) , " samples and ", nrow(mydata), "features were found."));
  current.msg <<- paste(msg, collapse="; ");
  mbSetObj$dataSet$read.msg <- current.msg;
  mbSetObj$dataSet$data.type <- "mothur";
  mbSetObj$dataSet$taxa.type <- taxa_type;
  mbSetObj$dataSet$module.type <- module.type;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  msg <- NULL;
  mydata <- .readDataTable(dataName);
    
  if(is.null(mydata) || class(mydata) == "try-error"){
    current.msg <<- "Failed to read in the taxonomic data! Please make sure the data is in the right format.";
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
    current.msg <<- "No labels #TAXONOMY found in your data!";
    return(0);
  }
    
  # converting to character matrix as duplicate row names not allowed in data frame.
  mydata <- as.matrix(mydata[,-1]);
  rownames(mydata) <- tax_nm;
  colnames(mydata) <- tax_rank;
  current.msg <<- paste("Taxonomy file has total of ",nrow(mydata),"features and",ncol(mydata), "taxonomic rank. Additional features which are not present in abundance table has been removed if present.");
    
  mbSetObj$dataSet$taxa_table <- mydata;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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
PlotSelectedSample <-function(mbSetObj, imgNm, smplID, idtype, OtuIdType, rel_perct,
                              format="png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_reshape();
  }

  if(current.selected.tax == "NA"){
    txlvl <- "Phylum";
  }else{
    txlvl <- current.selected.tax;
  }

  if(idtype=="16S"){
    data <- mbSetObj$dataSet$norm.phyobj;
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
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
    current.msg<<-"All features have lower relative abundance than given minimum abundance. Please lower the cut off for relative abundance.";
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
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(1);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
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
  data_new <- data.frame(data_new);
  data_new <- t(data_new);
  taxa_nm <- as.data.frame(datataxa[,txlvl]);
    
  if(OtuIdType=="GreengenesID"){
    
    if(.on.public.web){
      load_stringr();
    }
    
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
  data_new<-data.frame(data_new);
  data_new$step<-factor(rownames(data_new));
  data_new<-melt(data_new,id='step');
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
