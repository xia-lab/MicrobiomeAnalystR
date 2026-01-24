#' SanityCheckSampleData
#'
#' @param mbSetObj microbiomeanalyst object, initialized by InitDataObjects("pktable", "mf", FALSE)
#' @export
SanityCheckSampleData <- function(mbSetObj=NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type; 
  metadata <- mbSetObj$dataSet$sample_data
  smpl.nms <- rownames(metadata)
  
  # are sample names identical to data$orig
  dat.table <- readDataQs("data.proc", module.type, dataName);
  dat.table <- dat.table@.Data
  data.smpl.nms <- colnames(dat.table)
  
  nm.hits <- data.smpl.nms %in% smpl.nms
  if(!all(nm.hits)){
   perct <- round(sum(!nm.hits)/length(data.smpl.nms)*100, 3);
    AddErrMsg(paste0("A total of ", sum(!nm.hits), " or (", perct, "%) sample names are not present in the metadata file!" ));
    mis.nms <- data.smpl.nms[!nm.hits];
    AddErrMsg(paste0("Sample names missing in metadata file:", paste(mis.nms, collapse="; ")));
   return(0)
}

  
  nm.hits2 <- which(smpl.nms %in% data.smpl.nms);
  if(dim(metadata)[2] == 1){
    meta.col <- colnames(metadata)[1]
    meta.row <- rownames(metadata)[nm.hits2]
    metadata1 <- data.frame(V1 = metadata[nm.hits2, 1])
    colnames(metadata1) <- meta.col
    rownames(metadata1) <- meta.row
  } else {
    metadata1 <- metadata[nm.hits2,];
  }
  
  smpl.nms <- smpl.nms[nm.hits2];
  metadata1[] <- lapply( metadata1, factor);
  metadata <- metadata1
  
  # ensure order is the same
  if(dim(metadata)[2] == 1){
    meta.col <- colnames(metadata)[1]
    meta.row <- rownames(metadata)[match(data.smpl.nms, rownames(metadata))]
    metadata <- data.frame(V1 = metadata[match(data.smpl.nms, rownames(metadata)),])
    colnames(metadata) <- meta.col
    rownames(metadata) <- meta.row
  } else {
    metadata <- metadata[match(data.smpl.nms, rownames(metadata)),]
  }
  
  # write info to mbSetObj
  mbSetObj$dataSet$smpl_nm <- data.smpl.nms;
  mbSetObj$dataSet$meta.types <- rep("disc", ncol(metadata));
  mbSetObj$dataSet$meta.status <- rep("OK", ncol(metadata));
  names(mbSetObj$dataSet$meta.status) <- colnames(metadata);
  mbSetObj$dataSet$sample_data <- mbSetObj$dataSet$orig.sample_data <- metadata;
  mbSetObj$dataSet$types.cls.lbl <- sapply(metadata, function(x) class(x) ) 
  mbSetObj$dataSet$orig.cls <- mbSetObj$dataSet$cls <- metadata[,1]
  names(mbSetObj$dataSet$meta.types) <- colnames(mbSetObj$dataSet$sample_data)
  names(mbSetObj$dataSet$types.cls.lbl) <- colnames(mbSetObj$dataSet$sample_data)
  mbSetObj$dataSet$type.cls.lbl <- class(mbSetObj$dataSet$sample_data[,1]);
  return(.set.mbSetObj(mbSetObj));
}

GetMetaByCol <- function(mbSetObj=NA, metaname){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  metaData <- mbSetObj$dataSet$sample_data;
  return(metaData[,metaname]);
}

UpdateMetaData <- function(mbSetObj=NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(!exists('meta.nm.vec')){
    sel.meta.vecs <- mbSetObj$dataSet$orig.sample_data
  }else{
    sel.meta.vecs <- mbSetObj$dataSet$orig.sample_data[, meta.nm.vec]
  }
  mbSetObj$dataSet$sample_data <- sel.meta.vecs;
  return(.set.mbSetObj(mbSetObj))
}

RemoveSelectedMeta <- function(mbSetObj=NA, meta){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  inx <- which(colnames(mbSetObj$dataSet$sample_data) == meta);
  mbSetObj$dataSet$sample_data <- mbSetObj$dataSet$sample_data[, -inx];
  if(mbSetObj$module.type == "meta"){
    mbSetObj$dataSet$norm.phyobj@sam_data <- mbSetObj$dataSet$sample_data;
    mbSetObj$dataSet$proc.phyobj@sam_data <- mbSetObj$dataSet$sample_data;
  }

  inx1 <- which(names(mbSetObj$dataSet$meta.types) == meta);
  mbSetObj$dataSet$meta.types <- mbSetObj$dataSet$meta.types[-inx1];
  inx2 <- which(names(mbSetObj$dataSet$meta.status) == meta);
  mbSetObj$dataSet$meta.status <- mbSetObj$dataSet$meta.status[-inx2];
  inx3 <- which(names(mbSetObj$dataSet$types.cls.lbl) == meta);
  mbSetObj$dataSet$types.cls.lbl <- mbSetObj$dataSet$types.cls.lbl[-inx3];
  return(.set.mbSetObj(mbSetObj))
}

# Detect status and variable type of each metadata
#' SetMetaAttributes'
#' @param mbSetObj microbiomeanalyst object
#' @param init can be 0 or 1
#' @export
SetMetaAttributes <- function(mbSetObj=NA, init = 1){
  msg <- NULL;
  mbSetObj <- .get.mbSetObj(mbSetObj);
  sample_data  <- mbSetObj$dataSet$sample_data
  
  init <- as.numeric(init);
  cls.lbl <- sample_data[,1];
  cls.num <- length(levels(cls.lbl));
  msg <- c(msg, paste0("A total of ", length(colnames(mbSetObj$dataSet$sample_data)), " metadata factors were detected: ", paste0(colnames(mbSetObj$dataSet$sample_data), collapse=", "), "."));
  msg <- c(msg, paste0("The primary metadata factor is: ", colnames(mbSetObj$dataSet$sample_data)[1], ", which contains ", cls.num, " groups."));

  check.inx <-apply(sample_data , 2, function(x){(sum(is.na(x))/length(x) + sum(x=="NA", na.rm=TRUE)/length(x) + sum(x=="", na.rm=TRUE)/length(x)) >0})

  if(sum(check.inx)>0){
    if(init == 0){
      mbSetObj$dataSet$sample_data <- sample_data[,!check.inx]
    }else{
      mbSetObj$dataSet$meta.status[check.inx] = "<font color='red'>Missing values</font>";
    }
    msg <- c(msg, paste0( "<b>",paste0(colnames(sample_data)[check.inx], collapse=", "),"</b>", " meta-data factors have missing values."));
  }
  #print(msg)
  cls.vec <- vector()
  lowrep.vec <- vector()
  toolow.vec <- vector();
  for(i in 1:ncol(sample_data)){
      cls.Clean <- sample_data[,i];
      meta.name <- colnames(sample_data)[i]
      min.grp.size <- min(table(cls.Clean));
      cls.num <- length(levels(cls.Clean));


    # checking if too many groups but a few samples in each group
      if(cls.num/min.grp.size > 3){
        mbSetObj$dataSet$small.smpl.size <- 1;
        if(init == 1){
           isNum <- grepl("^-?[0-9.]+$", cls.Clean);
           if(all(isNum)){
             mbSetObj$dataSet$meta.types[i] = "cont";
             cls.vec <- c(cls.vec, meta.name)
           }else{
             if(!check.inx[i]){
             toolow.vec <- c(toolow.vec, meta.name)
             }
           }
        }
    # checking if some groups have low replicates
      } else if(min.grp.size < 3 | cls.num < 2){
        if(init == 1){
           isNum <- grepl("^-?[0-9.]+$", cls.Clean);
           if(all(isNum)){
             mbSetObj$dataSet$meta.types[i] = "cont";
             cls.vec <- c(cls.vec, meta.name)
           }else{
             if(all(c(!check.inx[i], !meta.name %in% toolow.vec))){
             lowrep.vec <- c(lowrep.vec, meta.name)
             }
           }
        }
      }
  
  }
 
  if(length(cls.vec) > 0){
    if(init == 0){
        mbSetObj$dataSet$sample_data <- sample_data[,which(!colnames(sample_data) %in% cls.vec)]
    }else if(init == 1){
        msg <- c(msg, paste0( "<b>",paste0(cls.vec, collapse=", "),"</b>", " meta-data factors are assigned to be continuous and remaining are categorical."));
        msg <- c(msg, "For categorical metadata, at least <b>two</b> groups with <b>three replicates</b> per group are required.");
        msg <- c(msg, "Please double check if these auto-assigned metadata types are correct.");
        msg <- c(msg, "You can manually update the metadata using the table below.");
    }
  }


  if(all(c(length(toolow.vec)>0, init == 1))){
    msg <- c(msg, paste0( "<b>",paste0(toolow.vec, collapse=", "),"</b>", " meta-data factors have too many groups with low replicates (less than 3) per group."));
  }

  if(all(c(length(lowrep.vec)>0, init == 1))){
    msg <- c(msg, paste0( "<b>",paste0(lowrep.vec, collapse=", "),"</b>", " meta-data factors have some groups with low replicates (less than 3) per group."));
  }

 if(init == 1){
  mbSetObj$msgSet$metacheck.msg <- msg;
  if(length(cls.vec)==1){
    sample_data[,cls.vec] = as.numeric(as.character(sample_data[,cls.vec]))
  }else{
    sample_data[,cls.vec] <- apply(sample_data[,cls.vec], 2, function(x){as.numeric(as.character(x))})
  }
  mbSetObj$dataSet$sample_data <- sample_data
}
  return(.set.mbSetObj(mbSetObj));
}


#' SetDataTypeOfMeta
#' @param mbSetObj microbiomeanalyst object
#' @export
SetDataTypeOfMeta <- function(mbSetObj=NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  sample_data  <- mbSetObj$dataSet$sample_data
  for(i in 1:ncol(sample_data)){
    metaNm <- colnames(sample_data)[i];
    if(mbSetObj$dataSet$meta.types[metaNm] == "cont"){
        sample_data[,i] = as.numeric(as.character(sample_data[,i] ));
    }
  }
  mbSetObj$dataSet$sample_data <- sample_data
  return(.set.mbSetObj(mbSetObj));
}

UpdateMetaColName <- function(mbSetObj=NA, oldName, newName){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  sample_data  <- mbSetObj$dataSet$sample_data
  inx <- which(colnames(sample_data) == oldName);
  mbSetObj$dataSet$meta.types <- mbSetObj$dataSet$meta.types[colnames(sample_data)]
  mbSetObj$dataSet$meta.status <- mbSetObj$dataSet$meta.status[colnames(sample_data)]
  mbSetObj$dataSet$types.cls.lbl <- mbSetObj$dataSet$types.cls.lbl[colnames(sample_data)]
  colnames(sample_data)[inx] <- newName
  names(mbSetObj$dataSet$meta.types)[inx] <- newName
  names(mbSetObj$dataSet$meta.status)[inx] <- newName
  names(mbSetObj$dataSet$types.cls.lbl)[inx] <- newName
  mbSetObj$dataSet$sample_data <- sample_data
  return(.set.mbSetObj(mbSetObj));

}

UpdateMetaLevels <- function(mbSetObj=NA,metaNm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(exists("meta.lvls")){ 
      sample_data  <- mbSetObj$dataSet$sample_data
      sel.meta <- sample_data[, metaNm]
      levels(sel.meta) <- meta.lvls
      sample_data <- cbind(sample_data, sel.meta);
      nms.vec <- colnames(sample_data)
      nms.vec[length(nms.vec)] <- metaNm
      nms.vec <- make.unique(nms.vec)
      colnames(sample_data) <- nms.vec
      new.nm <- colnames(sample_data)[length(nms.vec)]

     inx1 <- which(names(mbSetObj$dataSet$meta.types) == metaNm);
     mbSetObj$dataSet$meta.types <- c(mbSetObj$dataSet$meta.types, mbSetObj$dataSet$meta.types[inx1])
     names(mbSetObj$dataSet$meta.types)[length(nms.vec)] <- new.nm

     
     inx2 <- which(names(mbSetObj$dataSet$meta.status) == metaNm);
     mbSetObj$dataSet$meta.status <- c(mbSetObj$dataSet$meta.status, mbSetObj$dataSet$meta.status[inx2])
     names(mbSetObj$dataSet$meta.status)[length(nms.vec)] <- new.nm


     inx3 <- which(names(mbSetObj$dataSet$types.cls.lbl) == metaNm);
     mbSetObj$dataSet$types.cls.lbl <- c(mbSetObj$dataSet$types.cls.lbl, mbSetObj$dataSet$types.cls.lbl[inx3]);
     names(mbSetObj$dataSet$types.cls.lbl)[length(nms.vec)] <- new.nm

     mbSetObj$dataSet$sample_data <- sample_data;
     #print(head(sample_data));
  }
  return(.set.mbSetObj(mbSetObj));
}

GetMetaDataGroups <- function(){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$dataSet$sample_data)
}

UpdateMetaType <- function(mbSetObj=NA, metadata="NA", type="disc"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  meta.dat <- mbSetObj$dataSet$sample_data[[metadata]]
  if(type == "disc"){
    # check for at least two replicates in each level, then convert to factor
    meta.dat <- as.factor(as.character(meta.dat))
    reps <- as.data.frame(table(meta.dat))
    if(any(reps$Freq < 2)){
      return("You must have at least two replicates of each group to convert to a categorical variable!")
    } else {
      mbSetObj$dataSet$sample_data[[metadata]] <- meta.dat
      mbSetObj$dataSet$meta.types[metadata] <- type
      mbSetObj$dataSet$types.cls.lbl[metadata] <- "factor"
      return("Meta data type updated from continuous to categorical!")
    }
  } else {
    # check for alphabetic chars, then convert to number
    if(any(is.na(as.numeric(levels(meta.dat))))){
      return("You cannot convert a variable with non-numeric characters to a continuous variable!")
    } else {
      meta.dat <- as.numeric(levels(meta.dat))[meta.dat]
      mbSetObj$dataSet$sample_data[[metadata]] <- meta.dat
      mbSetObj$dataSet$meta.types[metadata] <- type
      mbSetObj$dataSet$types.cls.lbl[metadata] <- "numeric"
      return("Meta data type updated from categorical to continuous!")
    }
  }
  return(.set.mbSetObj(mbSetObj));
}

GetMetaTypes <- function(mbSetObj=NA, dataName){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$module.type == "meta"){
     mbSetObj$dataSet <- readDataset(dataName);
  }
  return(unname(mbSetObj$dataSet$meta.types));
}

GetMetaDataGroups <- function(mbSetObj=NA, dataName){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$module.type == "meta"){
     if(dataName == ""){
        mbSetObj$dataSet <- mbSetObj$dataSets[[1]]
     }else{
     mbSetObj$dataSet <- readDataset(dataName);
    }
  }
  colnms = colnames(mbSetObj$dataSet$sample_data)

  #print(head(mbSetObj$dataSet$sample_data))
  return(colnms[colnms!="sample_id"]);
  
}

GetMetaDataStatus <- function(mbSetObj=NA, dataName){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$module.type == "meta"){
     mbSetObj$dataSet <- readDataset(dataName);
  }
  return(unname(mbSetObj$dataSet$meta.status));
}


SetMetaTypes <- function(mbSetObj=NA, dataName, metaTypes.vec){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(mbSetObj$module.type == "meta"){
     mbSetObj$dataSet <- readDataset(dataName);
  }
  names(metaTypes.vec) <- colnames(mbSetObj$dataSet$sample_data)
  mbSetObj$dataSet$meta.types <- metaTypes.vec;
  return(.set.mbSetObj(mbSetObj));
}

UpdateMetaOrder <- function(mbSetObj=NA, metaName){
  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(exists('meta.ord.vec')){
   
if( isS4(mbSetObj$dataSet$sample_data)){
  metadata <- sample_data(mbSetObj$dataSet$sample_data)
  metadata[[metaName]] <- factor(metadata[[metaName]], levels = meta.ord.vec)
  mbSetObj$dataSet$sample_data <- sample_data(metadata)
}else{
  metadata <- mbSetObj$dataSet$sample_data[,metaName];
  mbSetObj$dataSet$sample_data[,metaName] <- factor(as.character(metadata), levels=meta.ord.vec)
}
 
  }

    return(.set.mbSetObj(mbSetObj));
}

GetUniqueMetaNames <-function(mbSetObj=NA, metadata){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(metadata == "NA"){
    metadata <- names(mbSetObj[["dataSet"]][["meta.types"]])[1];
  }

  data.type <- mbSetObj[["dataSet"]][["meta.types"]][metadata];
  
  if(data.type == "cont"){
    return("--- NA ---");
  } else {
    meta.dat <- mbSetObj$dataSet$sample_data
    if(class(meta.dat) == "data.frame"){
      return(levels(meta.dat[,metadata]))
    } else {
      return(levels(meta.dat[[metadata]]))
    }
  }
}

##############################################
##############################################
########## Utilities for web-server ##########
##############################################
##############################################

GetOrigSmplNms <-function(mbSetObj=NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$dataSet$smpl_nm);
}

GetOrigSmplGroupNamesPerMeta <- function(mbSetObj=NA, metaname="NA"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
if(class(mbSetObj$dataSet$sample_data) == "sample_data"){
  if(metaname %in% colnames(mbSetObj$dataSet$sample_data)){
    as.vector(mbSetObj$dataSet$sample_data[[metaname]]);
  }else{
    as.vector(mbSetObj$dataSet$sample_data[[1]]);
  }
}else{
  if(metaname %in% colnames(mbSetObj$dataSet$sample_data)){
    as.vector(mbSetObj$dataSet$sample_data[,metaname]);
  }else if(!is.null(mbSetObj$dataSet$sample_data[,1])) {
    as.vector(mbSetObj$dataSet$sample_data[,1]);
  }else {
    return(GetOrigSmplGroupNames(NA))
  }
}
}

UpdateSampleGroups<-function(mbSetObj=NA, metadata="NA"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  cls.lbl <- ClearStrings(as.vector(grp.vec));
  if(is.null(mbSetObj$dataSet$sample_data)) {
    mbSetObj$dataSet$sample_data <- matrix(nrow = length(cls.lbl))
  }
  sample_data <- mbSetObj$dataSet$sample_data;
  inx <- 1;
  if(metadata %in% colnames(sample_data)){
    inx <- which(colnames(sample_data) == metadata);
    type <- mbSetObj$dataSet$meta.types[inx];
    x <- cls.lbl

    if(type == "cont"){
        is.num <- T
        if(type == "cont"){
           isNum <- grepl("^-?[0-9.]+$", x);
           if(!all(isNum)){
                is.num <- F;
           }
        }

      if(!is.num){
          mbSetObj$dataSet$meta.status[inx] <- "<font color='red'>Not all numeric</font>"
      }else{
          mbSetObj$dataSet$meta.status[inx] <- "OK"
      }
    }else{

        containsMissing <-  sum(is.na(x))/length(x) + sum(x=="NA")/length(x) + sum(x=="")/length(x) + sum(x=="-")/length(x)  >0
        qb.inx <- tolower(cls.lbl) %in% c("qc", "blank");
        if(sum(qb.inx) > 0){
            cls.Clean <- as.factor(as.character(cls.lbl[!qb.inx])); # make sure drop level
        } else {
            cls.Clean <- as.factor(cls.lbl);
        }
        meta.name <- colnames(sample_data)[inx];
        min.grp.size <- min(table(cls.Clean));
        cls.num <- length(levels(cls.Clean));
        lowReplicate <- min.grp.size < 3 | cls.num < 2
        tooManyLow <- cls.num/min.grp.size > 4
        if(containsMissing){
            mbSetObj$dataSet$meta.status[inx] <- "<font color='red'>Missing values</font>"
        }else if (tooManyLow){
            mbSetObj$dataSet$meta.status[inx] <- "<font color='red'>Too many low replicates</font>"
        }else if (lowReplicate){
            mbSetObj$dataSet$meta.status[inx] <- "<font color='darkorange'>Low replicates</font>"
        }else{
            mbSetObj$dataSet$meta.status[inx] <- "OK"
        }
    }
  }else{
    mbSetObj$dataSet$orig.cls <- mbSetObj$dataSet$proc.cls <- mbSetObj$dataSet$prenorm.cls <- mbSetObj$dataSet$cls <- as.factor(cls.lbl);
  }
  if(mbSetObj$module.type == "meta"){
  mbSetObj$dataSet$sample_data[[inx]] = as.factor(cls.lbl);
  mbSetObj$dataSet$norm.phyobj@sam_data = mbSetObj$dataSet$sample_data;
  mbSetObj$dataSet$proc.phyobj@sam_data = mbSetObj$dataSet$sample_data;
  }else{
  mbSetObj$dataSet$sample_data[,inx] = as.factor(cls.lbl);
  }

  return(.set.mbSetObj(mbSetObj));
} 


##############################################
##############################################
################  GETTERS ####################
##############################################
##############################################

GetSampleGroups <- function(mbSetObj){
    mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$dataSet$group_names);
}


GetMetaDataCol <- function(mbSetObj,colnm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(colnm == "NA"){
    colnm <- 1;
  }

  df <- sample_data(mbSetObj$dataSet$sample_data)
  cls <- levels(df[[colnm]])
  cls <- cls[cls!="NA"];
  return(cls);
}
  
