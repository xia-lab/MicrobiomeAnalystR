##################################################
## R script for ExpressAnalyst
## Description: Functions related to web interface
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
##################################################

IncludeDataset <- function(mbSetObj=NA, dataName, include=1){
    mbSetObj <- .get.mbSetObj(mbSetObj);
    mbSetObj$mdata.all[dataName] <- include;
    #print(mbSetObj$mdata.all)
    return(.set.mbSetObj(mbSetObj));
}

GetSigGeneCount <- function(dataName){
  analSet <- readSet(analSet, "analSet");
  return(analSet$sig.gene.count);
}

CheckRawDataAlreadyNormalized <- function(dataName=""){
  dataSet <- readDataset(dataName);
  data <- dataSet$data.anot
  if(sum(data > 100) > 100){ # now we think it is raw counts
    return(0);
  }else{
    return(1);
  }
}

GetMetaCol<- function(dataName=""){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");;
  anal.type <- paramSet$anal.type;
  if(anal.type == "onedata"){
    colNms <- colnames(dataSet$comp.res);
    if (dataSet$de.method=="limma"){
      inx <- match("AveExpr", colNms)
    } else if (dataSet$de.method=="deseq2"){
      inx <- match("baseMean", colNms)
      return(colnames(dataSet$contrast.matrix));
    } else {
      inx <- match("logCPM", colNms)
    }
    resT <- dataSet$comp.res;
    if(inx > 2){
      resT <- resT[,1:inx-1];
      nms <- gsub("logFC.", "logFC_", colnames(resT));
      nms <- gsub("\\.", " vs ", nms);
      return(as.vector(nms));
    }else{
      return(dataSet$par1);
    }
  }else{
    nms <- paste(unique(dataSet$cls), collapse=" vs ");
    return(nms);
  }
}

GetSummaryData <- function(){
  msgSet <- readSet(msgSet, "msgSet");
  return(msgSet$summaryVec);
}

GetMetaColLength<- function(dataName=""){
  dataSet <- readDataset(dataName);
  paramSet <- readSet(paramSet, "paramSet");;

  if (dataSet$de.method=="limma"){
    inx <- match("AveExpr", colnames(dataSet$comp.res))
  } else if (dataSet$de.method=="deseq2"){
    inx <- match("baseMean", colnames(dataSet$comp.res))
    if(dataSet$contrast.type == "default"){
        return(dim(dataSet$contrast.matrix)[2]);
    }
  } else {
    inx <- match("logCPM", colnames(dataSet$comp.res))
  }
  resT <- dataSet$comp.res;
  resT <- resT[,1:inx-1]
  return(length(colnames(resT)));
}

GetInitLib <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  init.lib <- paramSet$init.lib;
  return(init.lib)
}

GetMetaDatasets<- function(){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];
  return(sel.nms);
}

SetSelMetaData<- function(selNm){
    paramSet <- readSet(paramSet, "paramSet");;
    paramSet$selDataNm <- selNm;
    paramSet$jsonNms$dataName <- selNm;
    saveSet(paramSet, "paramSet");
}

# only for switching single expression data results
SetCurrentData <- function(nm){
  dataSet <- readDataset(nm);
  
  return(1);
}

GetOmicsDataDims <- function(dataName){

dataSet <- readDataset(dataName);
  if("norm.phyobj" %in% names(dataSet)){
    dm <- dim(dataSet$norm.phyobj@otu_table);
  }else if("proc.phyobj" %in% names(dataSet)){
    dm <- dim(dataSet$proc.phyobj@otu_table);
  }else{
    dm <- dim(dataSet$data.orig);
  }
  naNum <- sum(is.na(dataSet$data.orig));

  return(c(dm, naNum));
} 
 

# given dataSet Name, sample name, and class name, do update
# note, for multiple #class, this set which one to use in the subsequent steps
# last one wins

# read in the data and perform
# gene ID mapping using built in libraries
# matchMin is minimal matched probe (%)
# return the total matched gene number

# obtain sample names and their class labels
GetSampleInfo <- function(dataName, clsLbl){
    dataSet <- readDataset(dataName);
    grpInfo <- dataSet$sample_data[[clsLbl]];
    grpLbls <- paste(levels(grpInfo), collapse="\n");
    smplInfo <- paste(Sample = colnames(dataSet$data.orig), "\t", Class=grpInfo, collapse="\n");
    return(c(grpLbls, smplInfo));
}

GetMetaSummary<- function(mbSetObj=NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
    feat.nums.by.dat <- c();
    total.col.num <- 0;
    for(i in 1:length(mbSetObj$dataSets)){
        dataSet <- readDataset(names(mbSetObj$dataSets)[i]);
        feat.nums.by.dat <- c(feat.nums.by.dat, nrow(otu_table(dataSet$proc.phyobj)));
        total.col.num  <- total.col.num + ncol(otu_table(dataSet$proc.phyobj));
    }
    feat.nums.by.dat <- paste(feat.nums.by.dat ,collapse="; ")
    mdata.all <- mbSetObj$mdata.all;
    sel.nms <- names(mdata.all)[mdata.all==1];
    sel.nms <- paste(sel.nms, collapse="; ")
    merged.data <- qs::qread("merged.data.raw.qs");
    cls.lbls <- sort(unname(unique(as.matrix(as.data.frame(sample_data(merged.data)))[,1])))
    cls.lbls <- paste(cls.lbls, collapse="; ")
    return(c(total.col.num, feat.nums.by.dat, sel.nms, cls.lbls));
}

GetDatasetNamesString <- function(){
    microbiome.meta <- qs::qread("microbiome_meta.qs");
    paste(unique(microbiome.meta$data.lbl), collapse="||");
}

##Single matrix
GetSampleNumber <-function(){
  data.proc <- qs::qread("data.proc.qs");
  return(ncol(data.proc));
}


GetFilesToBeSaved <-function(naviString){
  paramSet <- readSet(paramSet, "paramSet");
  return(unique(paramSet$partialToBeSaved));
}

###Gene list
GetNumOfLists <- function(){
  paramSet <- readSet(paramSet, "paramSet");
  return(paramSet$numOfLists)
}

GetMetaSigGeneCount <- function(){
  analSet <- readSet(analSet, "analSet");
  return(nrow(analSet$meta.mat));
}


##metastat



GetMetaGeneIDType<-function(){
  microbiome.meta <- qs::qread("microbiome_meta.qs");
  return(microbiome.meta$id.type);
}

GetMetaResultGeneIDs<-function(){
  analSet <- readSet(analSet, "analSet");

  rnms <- rownames(as.matrix(analSet$meta.mat.all));# already sorted based on meta-p values
  if(length(rnms) > 1000){
    rnms <- rnms[1:1000];
  }
  return(rnms);
}

# note, due to limitation of get/post
# maximum gene symb for list is top 500
GetMetaResultGeneSymbols<-function(){
  analSet <- readSet(analSet, "analSet");

  ids <- rownames(as.matrix(analSet$meta.mat.all));
  if(length(ids) > 1000){
    ids <- ids[1:1000];
  }
  microbiome.meta <- qs::qread("microbiome_meta.qs");
  if(microbiome.meta$id.type == "entrez"){ # row name gene symbols
    ids <- microbiome.meta$gene.symbls[ids];
  }
  return(ids);
}

GetMetaResultPathSymbols<-function(){
  analSet <- readSet(analSet, "analSet");
  meta.mat.all <- analSet$meta.mat.all;
  return(rownames(meta.mat.all));
}

GetMetaResultGeneIDLinks <- function(){
  paramSet <- readSet(paramSet, "paramSet");;
  analSet <- readSet(analSet, "analSet");

  ids <- rownames(as.matrix(analSet$meta.mat.all));
  if(length(ids) > 1000){
    ids <- ids[1:1000];
  }

  # set up links to genbank
  if(paramSet$data.org == "ko"){
    annots <- paste("<a href='https://www.genome.jp/dbget-bin/www_bget?", ids, "' target='_blank'>KEGG</a>", sep="");  
  }else{
    annots <- paste("<a href='http://www.ncbi.nlm.nih.gov/gene?term=", ids, "' target='_blank'>NCBI</a>", sep="");
  }
  return(annots);
}

GetMetaResultColNames<-function(){
  paramSet <- readSet(paramSet, "paramSet");
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mdata.all <- mbSetObj$mdata.all;
  sel.nms <- names(mdata.all)[mdata.all==1];

  # note, max 9 data columns can be displayed
  if(length(sel.nms) + ncol(analSet$meta.mat.all) > 9){
    max.col <- 9 - ncol(analSet$meta.mat.all);
    sel.nms <- sel.nms[1:max.col];
  }

  #no individual DE stat for now
  c(substring(sel.nms, 0, nchar(sel.nms)-4), colnames(analSet$meta.mat));
  #c(colnames(analSet$meta.mat))
}

# single.type return logFC or p value for individual data analysis
GetMetaResultMatrix<-function(single.type="fc"){

  analSet <- readSet(analSet, "analSet");

  if(single.type == "fc"){
    dat.mat <- analSet$fc.mat;
  }else{
    dat.mat <- analSet$pval.mat;
  }

  # note, max 9 data columns can be displayed
  if(ncol(dat.mat) + ncol(analSet$meta.mat.all) > 9){
    max.col <- 9 - ncol(analSet$meta.mat.all);
    dat.mat <- dat.mat[,1:max.col];
  }
  meta.mat2 <- cbind(dat.mat, analSet$meta.mat.all);

  # display at most 1000 genes
  if(nrow(meta.mat2) > 1000){
    meta.mat2 <- meta.mat2[1:1000,]; # already sorted based on meta-p values
  }

  meta.mat2 <-signif(as.matrix(meta.mat2), 5);
  meta.mat2;
  
}

GetMetaStat<-function(){
  analSet <- readSet(analSet, "analSet");
  return (analSet$meta.stat$stat);
}

GetMetaStatNames<-function(){
  analSet <- readSet(analSet, "analSet");
  return (names(analSet$meta.stat$stat));
}


# here should first try to load the original data
# the data in the memory could be changed
GetGroupNamesMeta <- function(dataName, meta="NA"){
    dataSet <- readDataset(dataName);
    if(meta == "NA"){
        return(levels(factor(dataSet$sample_data[,1])));
    }else{
        return(levels(factor(dataSet$sample_data[,meta])));
    }

}

setIncludeMeta <- function(metaBool){
    paramSet <- readSet(paramSet, "paramSet");
    paramSet$meta.selected <- metaBool;
    saveSet(paramSet, "paramSet");
}
PlotSelectedFeature<-function(mbSetObj, imgName, feat.id, format="png", sel.meta="", dpi=72){
  load_ggplot();
  load_grid();
  load_gridExtra();

  imgName <- paste(imgName,".", format, sep="");
  singleCol <- F;
  mbSetObj <- .get.mbSetObj(mbSetObj);
  phyobj <- qs::qread("metaanal_phyobj.qs");
  
  plot.data <- as.data.frame(as.matrix(otu_table(phyobj)))[feat.id,];
  a <- as.numeric(plot.data);
  min.val <- min(abs(a[a!=0]))/5;
  plot.data.norm <- log2((a + sqrt(a^2 + min.val))/2);
  
  plot.data.list <- list(plot.data, plot.data.norm);
  names(plot.data.list) <- c("Fitlered counts", "Log-transformed Count");
  
  sam_data <- as.data.frame(as.matrix(sample_data(phyobj)));
  cls.lbl <- sam_data[, sel.meta];
  num <- length(mbSetObj$dataSets);
  cmpdNm <- feat.id;
  # calculate width based on the dateset number
  if(num == 1){
    col <- unique(GetColorSchemaFromFactor(as.character(cls.lbl)));   
    df.norm <- data.frame(value=plot.data[feat.id,], name = as.character(cls.lbl))
    p.norm <- ggplot2::ggplot(df.norm, aes(x=name, y=value, fill=name))  
    p.norm <- p.norm + + geom_violin(trim = FALSE, aes(color = name), show.legend = FALSE) + geom_jitter(height = 0, width = 0.05, show.legend = FALSE)  + theme_bw()
    p.norm <- p.norm + theme(axis.title.x = element_blank(), axis.title.y = element_blank(), legend.position = "none")
    p.norm <- p.norm + stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)
    p.norm <- p.norm + scale_fill_manual(values=col) + 
      scale_color_manual(values=col) +
      ggtitle(cmpdNm) + theme(axis.text.x = element_text(angle=90, hjust=1))
    p.norm <- p.norm + theme(plot.title = element_text(size = 11, hjust=0.5)) 
    p.norm <- p.norm + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
    myplot <- p.norm + theme(plot.margin = margin(t=0.35, r=0.25, b=0.15, l=0.5, "cm"), axis.text = element_text(size=10))
    Cairo(file = imgName, width=280, height=320, type="png", bg="white");
    print(myplot);
    dev.off();
  }else{
    
    if(singleCol){
      layout <- c(1, num);
      height <- 200*num;
      width <- 280;
      row.num <- num;
    }else{
      layout <- c(1, num);
      height <- 200*num;
      width <- 400*length(plot.data.list);
    } 
    plot.list <- list();
    for(i in 1:length(plot.data.list)){
      data.lbl <- as.character(sample_data(phyobj)$dataset);
      data.lbl <- substr(data.lbl, 0, nchar(data.lbl)-4);
      
      # get counts in each data, same order as a levels
      counts <- table(data.lbl);
      # back to factor 
      data.lbl <- factor(data.lbl);
      # get new lbls to cut potential long names, and add sample numbers
      nlbls <- data.lbl;
      levels(nlbls) <- abbreviate(levels(nlbls),9);
      nlbls <- paste(levels(nlbls), "( n=", as.vector(counts), ")");
      # update labels
      data.lbl <- factor(data.lbl, labels=nlbls);
      # some time the transformed plot.data can switch class label, use the original data, need to be similar scale
      p_all <- list()
      out.fac <- data.lbl;
      in.fac <- cls.lbl;
      xlab <- sel.meta;
      
      col <- unique(GetColorSchemaFromFactor(as.character(cls.lbl)));   
      for(lv in levels(out.fac)){
        inx <- out.fac == lv;
        df.orig <- data.frame(facA = lv, value = unname(unlist(plot.data.list[[i]][inx])), name = in.fac[inx])
        p_all[[lv]] <- df.orig
      }
      
      alldata <- do.call(rbind, p_all);
      alldata$facA <- factor(as.character(alldata$facA), levels=levels(out.fac));
      
      p.time <- ggplot2::ggplot(alldata, aes(x=name, y=value, fill=name)) + stat_boxplot(geom ='errorbar') + 
        geom_boxplot(outlier.shape = NA) +  geom_jitter(height = 0, width = 0.05, show.legend = FALSE) + theme_bw()
      p.time <- p.time + facet_wrap(~facA, nrow = num) + theme(axis.title.x = element_blank(), legend.position = "none")
      p.time <- p.time + scale_fill_manual(values=col) + 
        scale_color_manual(values=col) +
        theme(axis.text.x = element_text(angle=90, hjust=1))
      p.time <- p.time + ggtitle(names(plot.data.list)[i]) + theme(plot.title = element_text(size = 11, hjust=0.5, face = "bold")) + ylab("Expression")
      p.time <- p.time + theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank()) # remove gridlines
      myplot <- p.time + theme(plot.margin = margin(t=0.15, r=0.25, b=0.15, l=0.25, "cm"), axis.text = element_text(size=10)) 
      myplot <- myplot + coord_flip()
      plot.list[[i]] <- myplot;
    }

    Cairo(file = imgName, width=width, height=height, type="png", bg="white");
    grid.arrange(plot.list[[1]], plot.list[[2]],ncol=2,nrow=1,top=feat.id);
    dev.off();
  }
  return(.set.mbSetObj(mbSetObj));
}

GetResColNamesMeta <- function(mbSetObj, type="alpha"){
    mbSetObj <- .get.mbSetObj(mbSetObj);

    if(type == "alpha"){
        vec <- colnames(mbSetObj$analSet$alpha.summary)[-c(1,2)];
    }else if (type == "beta"){
        vec <- colnames(mbSetObj$analSet$beta.summary)[-c(4,5)];
    }
    return(vec)
}
GetResRowNames <- function(mbSetObj, type="alpha"){
    mbSetObj <- .get.mbSetObj(mbSetObj);

    if(type == "alpha"){
        vec <- mbSetObj$analSet$alpha.summary[,1];
    }else if (type == "beta"){
        vec <- mbSetObj$analSet$beta.summary$dataset;
    }

    return(vec)
}

GetResMatrix<- function(mbSetObj, type="alpha"){
    mbSetObj <- .get.mbSetObj(mbSetObj);

    if(type == "alpha"){
        df <- mbSetObj$analSet$alpha.summary[,-c(1:2)];
    }else if (type == "beta"){
        df <- mbSetObj$analSet$beta.summary[,-c(4,5)];
    }

    return(signif(as.matrix(df), 5))
}

GetResMetric<- function(mbSetObj, type="alpha"){
    mbSetObj <- .get.mbSetObj(mbSetObj);

    if(type == "alpha"){
        vec <- mbSetObj$analSet$alpha.summary$Metric;
    }else if (type == "beta"){
        vec <- mbSetObj$analSet$beta.summary$dist;
    }

    return(vec)
}
