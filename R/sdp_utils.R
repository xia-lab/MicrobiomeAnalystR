##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

###########################
########### I/O ###########
###########################

#'Function to get gene list statistics
#'@description This function gets the gene list stats.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

GetGeneListStat <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  range <- c(0, 0);
    
  if(mbSetObj$analSet$gene.only){
    val.range <- c(0, 0);
  }else{
    val.range <- range(mbSetObj$analSet$ko.mapped[,1])
  }
  
  return(c(nrow(mbSetObj$analSet$ko.orig), nrow(mbSetObj$analSet$ko.mapped), mbSetObj$analSet$gene.only, val.range));
}

#'Function to filter list data based on a minimum count
#'@description This function filters the inputted list data
#'based on a minimum count.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

UpdateListInput <- function(mbSetObj, minL, maxL=Inf){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  hit.inx <- mbSetObj$analSet$ko.mapped >= minL & mbSetObj$analSet$ko.mapped <= maxL;
    
  if(sum(hit.inx) > 0){
    current.msg <<- paste("A total of unqiue", sum(hit.inx), "KO genes were selected!");
    mbSetObj$analSet$data <- mbSetObj$analSet$ko.mapped[hit.inx, , drop=FALSE];
    return(.set.mbSetObj(mbSetObj));
  }else{
    AddErrMsg("No genes were selected in this range!");
    return(0);
  }
}

#'Main function to read in shotgun data
#'@description This function reads in shotgun data
#'in a tabular format.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
ReadShotgunTabData <- function(mbSetObj, dataName, geneidtype, datatype) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  mydata <- .readDataTable(dataName);
    
  if(any(is.na(mydata)) || class(mydata) == "try-error"){
    AddErrMsg("Failed to read in the abundance data! Please make sure the gene abundance table is in the right format and do not have empty cells or NA.");
    return(0);
  }

  # look for #NAME, store in a list
  #just to check the name of labels.
  sam.nm <- substr(colnames(mydata[1]),1,5);
  sam.nm <- tolower(sam.nm);
  sam.inx <- grep("^#name",sam.nm);
    
  if(length(sam.inx) > 0){
    smpl_nm<-colnames(mydata[-1]);
  }else{
    AddErrMsg("No labels #NAME found in your data!");
    return(0);
  }

  dat.nms <- mydata[,1];
  mydata <- as.matrix(mydata[,-1]);
    
  if(mode(mydata)=="character"){
    AddErrMsg(paste("Errors in parsing your data as numerics - possible reason: comma as decimal separator?"));
    return(0);
  }
  
  rownames(mydata) <- dat.nms;

  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$data.orig <- mydata;
  mbSetObj$dataSet$smpl_nm <- smpl_nm;
  current.msg <<- paste("A total of ",ncol(mydata) , " samples and ", nrow(mydata), " metagenomic features are present.");
  mbSetObj$dataSet$read.msg<-current.msg;
  mbSetObj$dataSet$data.type<-"text";
  mbSetObj$dataSet$gene.id<-geneidtype;
  
  return(.set.mbSetObj(mbSetObj));

}

#'Main function to read in shotgun data
#'@description This function reads in shotgun data
#'in biom format.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import biomformat
ReadShotgunBiomData <- function(mbSetObj, dataName, geneidtype, module.type, ismetadata) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  #to support biom file produced by picrust, (read_biom:phyloseq) function doesn't work
  msg <- NULL;
    
  if(on.public.web){
    load_biomformat();
  }
  
  #reading .biom file using phyloseq
  mydata <- biomformat::read_biom(dataName);

  #converting to matrix for further use
  mydata = as(biom_data(mydata), "matrix")
  otu.dat = otu_table(mydata, taxa_are_rows=TRUE)
    
  if(length(otu.dat)==0){
    AddErrMsg("Biom file does not contain abundance information");
    return(0);
  }
    
  otu.dat<-as.matrix(otu.dat);
  msg <- c(msg, "Abundance data present.");

  mbSetObj$dataSet$name <- basename(dataName);
  mbSetObj$dataSet$data <- otu.dat;

  #sample file if present within
  if(ismetadata=="T"){
    sample_data<-sample_data(mydata,errorIfNULL = FALSE);
    if(length(sample_data)==0){
      AddErrMsg("Metadata file not detected in your biom file. Please upload metadeta file seperately.");
      ismetafile<-"F";
      return(0);
    }
    mbSetObj$dataSet$sample_data<-as.data.frame(sample_data);
    msg <- c(msg, "Metadata file is also detected in your biom file.");
  }

  msg <- c(msg, paste("A total of ",ncol(mydata) , " samples and ", nrow(mydata), " metagenomic features were found."));
  current.msg <<- paste(msg, collapse="; ");
  mbSetObj$dataSet$read.msg<-current.msg;
  mbSetObj$dataSet$data.type<-"biom";
  mbSetObj$module.type<-module.type;
  mbSetObj$dataSet$data.orig <- otu.dat;
  mbSetObj$dataSet$gene.id<-geneidtype;
  
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to prepare shotgun data for PCA.
#'@description This function formats shotgun data for PCA.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import RJSONIO
#'@import ggfortify
PreparePCA4Shotgun <- function(mbSetObj, imgName,imgName2, format="json", inx1, inx2, inx3,
                              variable, showlabel, format2d="png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_rjsonio();
  load_ggfortify();
  
  imgName2 = paste(imgName2, ".", format2d, sep="");
  mbSetObj$imgSet$pca<-imgName2;

  dat <- as.matrix(otu_table(mbSetObj$dataSet$norm.phyobj));
  pca3d <- list();
  pca <- prcomp(t(dat), center=T, scale=T);
  imp.pca <- summary(pca)$importance;
  fast.write(signif(pca$x,5), file="pca_score.csv");
  pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(imp.pca[2,][1:3], 3), "%)", sep="");
  
  #3D
  coords <- data.frame(t(signif(pca$x[,1:3], 5)));
  colnames(coords) <- NULL;
  pca3d$score$type <- "factor";
  pca3d$score$xyz <- coords;
  pca3d$score$name <- sample_names(mbSetObj$dataSet$norm.phyobj);
  sam_data <- data.frame(sample_data(mbSetObj$dataSet$norm.phyobj));
  cls <- as.character(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  clsLbl <- sam_data[[variable]];
  pca3d$score$facA <- cls;
  variable <<- variable;
  
  # now set color for each group
  cols <- unique(as.numeric(factor(cls)))+1;
  rgbcols <- col2rgb(cols);
  cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
  pca3d$score$colors <- cols;
  
  json.obj <- RJSONIO::toJSON(pca3d);
  sink(imgName);
  cat(json.obj);
  sink();

  #2D
  Cairo::Cairo(file=imgName2, width=720, height=500, type=format2d, bg="white",dpi=dpi);
  label = FALSE;
    
  if(showlabel=="samnm"){
    label = TRUE;
    box <- autoplot(pca,data=sam_data,colour=variable,size=4,alpha =0.8,label = label) + theme_bw()
    box <- box + stat_ellipse(type="norm", linetype=2, geom = "polygon",alpha = 0.2, aes_string(fill = quo(clsLbl)), show.legend=FALSE)
    box <- box + labs(x = pca3d$score$axis[1], y = pca3d$score$axis[2]);
  } else if(showlabel=="none") {
    box <- autoplot(pca,data=sam_data,colour=variable,size=4,alpha =0.8,label = label) + theme_bw()
    box <- box + stat_ellipse(type="norm", linetype=2, geom = "polygon",alpha = 0.2, aes_string(fill = quo(clsLbl)), show.legend=FALSE)
    box <- box + labs(x = pca3d$score$axis[1], y = pca3d$score$axis[2]);
  } else {
    grplbl <<- sam_data[ ,showlabel];
    clsLbl <<- clsLbl;
    box <- autoplot(pca,data=sam_data,colour=variable,size=4,alpha =0.8,label = label) + geom_text(aes(label=grplbl,colour=clsLbl))
    box <- box + theme_bw()+ stat_ellipse(type="norm", linetype=2, geom = "polygon",alpha = 0.2, aes_string(fill = quo(clsLbl)), show.legend=FALSE)
    box <- box + labs(x = pca3d$score$axis[1], y = pca3d$score$axis[2]);
  }
  print(box);
  dev.off();
  mbSetObj$analSet$pca<-pca;
  return(.set.mbSetObj(mbSetObj))
}

#############################################
###########Functional Profiling #############
#############################################

#'Function to plot stacked bar chart of functional data.
#'@description This function plots stacked bar charts of functional data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotFunctionStack<-function(mbSetObj, summaryplot, functionlvl, abundcal, geneidtype, metadata,
                            colpalopt, format="png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  summaryplot <- paste(summaryplot, ".", format, sep="");
  
  load_reshape();

  data <- mbSetObj$dataSet$proc.phyobj;
  smpl_nm <- sample_names(data);
  clsLbl <- factor(sample_data(data)[[metadata]]);
  if(length(levels(clsLbl)) > 9 && min(table(clsLbl)) < 3){
    AddErrMsg("Too many facets to be displayed - please select a more meaningful facet option with at least 3 samples per group.");
    return(0);
  }
  #reorder data based on groups
  ord.inx <- order(clsLbl);
  clsLbl <- clsLbl[ord.inx];
  colvec <- as.numeric(clsLbl) + 1;
  smpl_nm <- smpl_nm[ord.inx];
  query <- otu_table(data)[,ord.inx];

  if(geneidtype=="cog"){
    ko_higher_path <- .read.microbiomeanalyst.lib.rds("cog_functioncount.rds", "ko");
  }else {
    if(functionlvl=="KEGG metabolism"){
      ko_higher_path<-.read.microbiomeanalyst.lib.rds("ko_higherpathway.rds", "ko");
    }else if(functionlvl=="KEGG pathways"){
      ko_higher_path<-.read.microbiomeanalyst.lib.rds("ko_pathwaycount.rds", "ko");
    }else if(functionlvl=="KEGG modules"){
      ko_higher_path<-.read.microbiomeanalyst.lib.rds("ko_modulecount.rds", "ko");
    }else if(functionlvl=="COG Functional category"){
      ko_higher_path<-.read.microbiomeanalyst.lib.rds("ko_cogfunction.rds", "ko");
    }
  }

  clsLbl_new <- as.character(clsLbl);

  #sample,categories name
  samplenm <- colnames(query);
  categnm <- colnames(ko_higher_path)<-gsub("\\.", " ",colnames(ko_higher_path));

  #merging user KO query with the require KO file
  result2 <- merge(query,ko_higher_path, by ="row.names");
  indx2 <- match(samplenm,colnames(result2), nomatch = NA_integer_, incomparables = NULL);
  indx1 <- match(categnm,colnames(result2), nomatch = NA_integer_, incomparables = NULL);

  #actual weight of node/no of pathways in which it is present
  if(abundcal=="weig_hit"){
    result2[indx1] <- result2[indx1]/rowSums(result2[,indx1]);
    #removing NA introduced
    result2 <- replace(result2, is.na(result2), 0);
  }
    
  myList <- vector('list', length(indx2));
  for (i in 1:length(indx2)) {
    myList[[i]] <- as.data.frame(colSums(result2[indx1]*result2[,indx2[i]]));
  }

  MyMerge<- function(x, y){
    df <- merge(x, y, by= "row.names", all.x= F, all.y= F);
    rownames(df) <- df$Row.names
    df$Row.names <- NULL
    return(df)
  }
    
  result <- Reduce(MyMerge, myList);

  #orignal class label
  colnames(result) <- samplenm;
  result <- result[rowSums(result)!=0, ];

  if(abundcal=="norm_hit"){
    #getting the size of categories
    categ_size <- as.data.frame(colSums(ko_higher_path));
    colnames(categ_size) <- "size";
    result <- merge(result,categ_size, by ="row.names");
    result[indx2] <- result[indx2]/result[['size']];

    #removing the extra (size) column
    result <- result[ ,c(1,indx2)];
    rownames(result) <- result[,1];
    result <- result[,-1];

    #removing zero abundance KEGG pathways, metabolism and modules
    result <- result[ rowSums(result)!=0, ];
  }
    
  fast.write(result, file="funcprof_abund.csv");

  # now plotting
  w <- 1000;
    
  if(nrow(result)>100){
    w<-w+200;
  }else if(nrow(result)>150){
    w<-w+300;
  }
  ##selecting 250 samples
  subsmpl <- 250;
    
  if (ncol(result)>subsmpl) {
    ss  <- sample(ncol(result), subsmpl);
    result <- result[,ss,drop=FALSE];
  } else {
    result <- result;
  }

  data<-t(result);

  nms <- colnames(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data);

  data$facetOpt <- as.character(clsLbl_new);
  data$step <- factor(rownames(data), levels = rownames(data));
  data <- melt(data,id=c('step', 'facetOpt'));


  data$step <- as.numeric(data$step);
  data$variable <- gsub("\\.", " ",data$variable); # remove the dot introduced for readability
  #color schema
  color_var <- levels(factor(data$variable));
  x <- length(color_var);
    
  if(colpalopt=="grad"){
    indx <- which(color_var=="NA");
    #color schema for ggplot
    x.colors <- hcl(h=seq(15,375,length=(x+1)),l=65,c=100)[1:x];
    x.colors[indx] <- "#666666";
  }else if (colpalopt=="cont21"){
    x.colors <- rep(custom_col21,length.out=x);
  }else if (colpalopt=="cont28"){
    x.colors <- rep(custom_final28,length.out=x);
  }else {
    x.colors <- rep(custom_col42,length.out=x);
  }
    
  Cairo::Cairo(file=summaryplot,width=w, height=600, type=format, bg="white",dpi=dpi);
  mbSetObj$imgSet$func.prof<-summaryplot;

  box <- ggplot(data,aes(x=step,y=value)) + 
    facet_grid(~ facetOpt, space = "free", scales = "free") +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust =1, vjust=0.5)) +
    geom_area(aes(fill=variable),position='fill') +
    scale_x_continuous(breaks=seq(1,length(unique(data$step)),1),labels=smpl_nm) +
    #scale_fill_manual(values=c(x.colors))+
    labs(y=" Relative Abundance",fill=functionlvl) +
    theme(axis.text.y = element_text(size = 10)) +
    theme(legend.text=element_text(size=10), strip.text = element_text(size = 12)) +
    theme(axis.text.x = element_text(colour="black", size = 10), axis.title.x=element_blank()); 
   
  if(colpalopt=="set3"){
    cols.needed <- length(unique(data$variable))
    if(cols.needed > 12){
      col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      box <- box + scale_fill_manual(values=col.func(cols.needed),
                                     #guide = guide_legend(direction = "horizontal", ncol = 5)
                                    )
    } else {
      box <- box + scale_fill_brewer(palette = "Set3",
                                     #guide = guide_legend(direction = "horizontal",ncol = 5)
                                    )                           
    }
  } else {
    box <- box + scale_fill_manual(values=c(x.colors))
  }

  print(box);
  dev.off();
  mbSetObj$analSet$func.prof<-data;
  mbSetObj$analSet$func.lvl<-functionlvl;
  return(.set.mbSetObj(mbSetObj))
}

#'Function to perform KO mapping.
#'@description This function performs KO mapping.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PerformKOmapping <- function(mbSetObj, geneIDs, type){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  mbSetObj$analSet <- list();
  mbSetObj$analSet$orig <- geneIDs;
  current.msg <<- NULL;

  lines <- unlist(strsplit(geneIDs, "\r|\n|\r\n")[1]);
    
  if(substring(lines[1],1,1)=="#"){
    lines <- lines[-1];
  }
    
  gene.lists <- strsplit(lines, "\\s+");
  gene.mat <- do.call(rbind, gene.lists);

  if(dim(gene.mat)[2] == 1){ # add 1
    gene.only <- 1;
    gene.mat <- cbind(gene.mat, rep(1, nrow(gene.mat)));
  }else {
    gene.only <- 0;
    gene.mat <- gene.mat[,1:2];
  }

  rownames(gene.mat) <- gene.mat[,1];
  gene.mat <- gene.mat[,-1, drop=F];
  mbSetObj$analSet$ko.orig <- gene.mat;
  mbSetObj$analSet$gene.only <- gene.only;
  gene.mat <- RemoveDuplicates(gene.mat, "sum", quiet=F);
  mbSetObj$analSet$ko.uniq <- gene.mat;

  # now get input that are in the lib
  kos <-  doKOFiltering(rownames(gene.mat), type);

  if(sum(!is.na(kos)) < 2){
    AddErrMsg("Less than two hits found in the database. ");
    return(0);
  }else{
    rownames(gene.mat) <- kos;
    gd.inx <- (!is.na(kos)) & gene.mat[,1] > 0;
    gene.mat <- gene.mat[gd.inx, ,drop=F];
    mbSetObj$analSet$ko.mapped <- mbSetObj$analSet$data <- gene.mat; # data will be updated, ko.map will keep intact
    current.msg <<- paste("A total of unique", nrow(gene.mat), "KO genes were mapped to KEGG network!");

    return(.set.mbSetObj(mbSetObj));
  }
}

#'Function to prepare query for JSON.
#'@description This function prepares the data for JSON.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import RJSONIO

PrepareQueryJson <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_rjsonio();

  if(enrich.type == "hyper"){
    exp.vec <- mbSetObj$analSet$data[,1]; # drop dim for json
  }else{
    # for global test, all KO measured should be highlighted
    genemat <- as.data.frame(t(otu_table(mbSetObj$dataSet$norm.phyobj)));
    exp.vec <- rep(2, ncol(genemat));
    names(exp.vec) <- colnames(genemat);
  }
    
  edge.mat <- MapKO2KEGGEdges(exp.vec);
  row.names(edge.mat) <- eids <- rownames(edge.mat);
  query.ko <- edge.mat[,1];
  net.orig <- edge.mat[,2];
  query.res <- edge.mat[,3];# abundance
  names(query.res) <- eids; # named by edge

  json.mat <- RJSONIO::toJSON(query.res, .na='null');
  sink("network_query.json");
  cat(json.mat);
  sink();

  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to prepare KO enrichment analysis.
#'@description This function performs KO enrichment analysis
#'using the KO01100 map.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PerformKOEnrichAnalysis_KO01100 <- function(mbSetObj, category, file.nm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
    
  if(enrich.type == "hyper"){
  LoadKEGGKO_lib(category);
    PerformKOEnrichAnalysis_List(mbSetObj, file.nm);
  }else{
    .prepare.global(mbSetObj, category, file.nm);
    .perform.computing();
    .save.global.res();
  }
  return(.set.mbSetObj(mbSetObj))
}

.prepare.global<-function(mbSetObj, category, file.nm){
  LoadKEGGKO_lib(category);
  mbSetObj <- .get.mbSetObj(mbSetObj);

  phenotype <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[selected.meta.data]]);
  genemat <- as.data.frame(t(otu_table(mbSetObj$dataSet$norm.phyobj)));
  # first, get the matched entries from current.geneset
  hits <- lapply(current.geneset, function(x){x[x %in% colnames(genemat)]});
  set.num <- unlist(lapply(current.geneset, length), use.names = FALSE);
  dat.in <- list(cls=phenotype, data=genemat, subsets=hits, set.num=set.num, filenm=file.nm);

  my.fun <- function(){
        gt.obj <- globaltest::gt(dat.in$cls, dat.in$data, subsets=dat.in$subsets);
        gt.res <- globaltest::result(gt.obj);

        match.num <- gt.res[,5];
        if(sum(match.num>0)==0){
            return(NA);
        }
        raw.p <- gt.res[,1];

        # add adjust p values
        bonf.p <- p.adjust(raw.p, "holm");
        fdr.p <- p.adjust(raw.p, "fdr");

        res.mat <- cbind(set.num, match.num, gt.res[,2], gt.res[,3], raw.p, bonf.p, fdr.p);
        rownames(res.mat) <- names(hits);
        colnames(res.mat) <- c("Size", "Hits", "Statistic Q", "Expected Q", "Pval", "Holm p", "FDR");
        hit.inx <- res.mat[,2]>0;
        res.mat <- res.mat[hit.inx, ];
        ord.inx <- order(res.mat[,5]);
        res.mat <- res.mat[ord.inx,];
        return(res.mat);
    }
    dat.in <- list(cls=phenotype, data=genemat, subsets=hits, set.num=set.num, filenm=file.nm , my.fun=my.fun);
    qs::qsave(dat.in, file="dat.in.qs");
    return(1);
}

.save.global.res <- function(){
  mbSetObj = .get.mbSetObj("NA");
  dat.in <- qs::qread("dat.in.qs"); 
  hits = dat.in$subsets
  file.nm = dat.in$filenm;
  my.res <- dat.in$my.res;
  if(length(my.res)==1 && is.na(my.res)){
        AddErrMsg("No match was found to the selected metabolite set library!");
        return(0);
  }

    # in R, sort list is by its name!, using pos order has issues!
    nms <- rownames(my.res);
    hits <- hits[nms];

    Save2KEGGJSON(hits, my.res, file.nm);
    .set.mbSetObj(mbSetObj);
    return(1);
}

#'Perform KO Enrichment Analysis 
#'@description This functions performs KO enrichment analysis
#'on a tabled input.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PerformKOEnrichAnalysis_Table <- function(mbSetObj, file.nm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  phenotype <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[selected.meta.data]]);
  genemat <- as.data.frame(t(otu_table(mbSetObj$dataSet$norm.phyobj)));
  # first, get the matched entries from current.geneset
  hits <- lapply(current.geneset, function(x){x[x %in% colnames(genemat)]});
  set.num <- unlist(lapply(current.geneset, length), use.names = FALSE);

   # now, perform the enrichment analysis
   library(RSclient);
    rsc <- RS.connect();
    RS.assign(rsc, "my.dir", getwd()); 
    RS.eval(rsc, setwd(my.dir));

    gt.out <- list(cls=phenotype, data=genemat, subsets=hits, set.num=set.num);
    RS.assign(rsc, "gt.in", gt.out); 

    # there are more steps, better drop a function to compute in the remote env.
    my.fun <- function(){
        gt.obj <- globaltest::gt(gt.in$cls, gt.in$data, subsets=gt.in$subsets);
        gt.res <- globaltest::result(gt.obj);

        match.num <- gt.res[,5];
        if(sum(match.num>0)==0){
            return(NA);
        }
        raw.p <- gt.res[,1];

        # add adjust p values
        bonf.p <- p.adjust(raw.p, "holm");
        fdr.p <- p.adjust(raw.p, "fdr");

        res.mat <- cbind(set.num, match.num, gt.res[,2], gt.res[,3], raw.p, bonf.p, fdr.p);
        rownames(res.mat) <- names(hits);
        colnames(res.mat) <- c("Size", "Hits", "Statistic Q", "Expected Q", "Pval", "Holm p", "FDR");
        hit.inx <- res.mat[,2]>0;
        res.mat <- res.mat[hit.inx, ];
        ord.inx <- order(res.mat[,5]);
        res.mat <- res.mat[ord.inx,];
        return(res.mat);
    }
    RS.assign(rsc, my.fun);
    my.res <- RS.eval(rsc, my.fun());
    RS.close(rsc);

    if(length(my.res)==1 && is.na(my.res)){
        AddErrMsg("No match was found to the selected metabolite set library!");
        return(0);
    }

    # in R, sort list is by its name!, using pos order has issues!
    nms <- rownames(my.res);
    hits <- hits[nms];

    Save2KEGGJSON(hits, my.res, file.nm);
    return(.set.mbSetObj(mbSetObj));
}

# Utility function
LoadKEGGKO_lib<-function(category){
    
  if(category == "module"){
    kegg.anot <- .read.microbiomeanalyst.lib.rds("ko_modules.rds", "ko")
    current.setlink <- kegg.anot$link;
    current.mset <- kegg.anot$sets$"Pathway module";
  }else{
    kegg.anot <- .read.microbiomeanalyst.lib.rds("ko_pathways.rds", "ko")
    current.setlink <- kegg.anot$link;
    current.mset <- kegg.anot$sets$Metabolism;
  }
    
  # now need to update the msets to contain only those in ko01100 map
  if(!exists("ko.edge.map")){
    
    if(.on.public.web){
      ko.edge.path <- paste("../../lib/ko/ko_edge.csv", sep="");
    }else{
      ko.edge.path <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/ko/ko_edge.csv", sep="");
    }
    ko.edge.map <<- .readDataTable(ko.edge.path);
  }

  kos.01100 <- ko.edge.map$gene[ko.edge.map$net == "ko01100"];
  current.mset <- lapply(current.mset, function(x) {
            as.character(unique(x[x %in% kos.01100]))});
  # remove those empty ones
  mset.ln <- lapply(current.mset, length);
  current.mset <- current.mset[mset.ln > 0];
  set.ids <- names(current.mset);
  names(set.ids) <- names(current.mset) <- kegg.anot$term[set.ids];

  current.setlink <<- current.setlink;
  current.setids <<- set.ids;
  current.geneset <<- current.mset;
  current.universe <<- unique(unlist(current.mset));
}

###################################################
########### Metabolic Pathway projection ##########
###################################################

#'Function to perform KO projection
#'@description This function projects user-uploaded
#'KOs onto the KEGG Global Metabolic Network.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PerformKOProjection <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  datatype <- mbSetObj$analSet$datatype;
  de.Num <- mbSetObj$analSet$sig.count;
  resTable <- mbSetObj$analSet$resTable[1:de.Num,];
    
  if(nrow(resTable) > 5000){
    resTable <- resTable[1:5000,];
    current.msg <<- paste(" Due to computational constraints, only top 5000 features will be used. ", collapse="\n");
  }

  mbSetObj$analSet$sig.mat <- resTable;
  gene.mat <- data.matrix(-log10(resTable$Pvalues));
  rownames(gene.mat) <- row.names(resTable);
  gene.mat <- RemoveDuplicates(gene.mat, "sum", quiet=F);
  mbSetObj$analSet$ko.uniq <- gene.mat;
  id.type <- mbSetObj$analSet$id.type;
    
  if(is.na(id.type) || id.type == "NA"){
    id.type <- "ko";
  }
    
  kos <- doKOFiltering(rownames(gene.mat),id.type);
    
  if(sum(!is.na(kos)) < 2){
    AddErrMsg("Less than two hits found in the database.");
    return(0);
  }else{
    rownames(gene.mat) <- kos;
    gd.inx <- (!is.na(kos)) & gene.mat[,1] > 0;
    gene.mat <- gene.mat[gd.inx, ,drop=F];
    mbSetObj$analSet$ko.mapped <- mbSetObj$analSet$data <- gene.mat; # data will be updated, ko.map will keep intact
    current.msg <<- paste("A total of unqiue", nrow(gene.mat), "KO genes were mapped to KEGG network!");
  }

  return(.set.mbSetObj(mbSetObj));
  
}

##############################################
############## Utility Functions #############
##############################################

# Utility function
MapKO2KEGGEdges<- function(kos, net="ko01100"){
  
  if(!exists("ko.edge.map")){
    
    if(.on.public.web){
      ko.edge.path <- paste("../../lib/ko/ko_edge.csv", sep="");
    }else{
      ko.edge.path <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/ko/ko_edge.csv", sep="");
    }
    ko.edge.map <<- .readDataTable(ko.edge.path);
  }
  
  all.hits <- ko.edge.map$gene %in% names(kos) & ko.edge.map$net == net;
  my.map <- ko.edge.map[all.hits, ];
  q.map <- data.frame(gene=names(kos), expr=as.numeric(kos));
  
  # first merge to get ko abundance to each edge
  dat <- merge(my.map, q.map, by="gene");
  
  # now merge duplicated edge to sum
  dup.inx <- duplicated(dat[,2]);
  dat <- dat[!dup.inx,];
  rownames(dat) <- dat[,2];
  
  return(dat[,-2]);
}

# Utility function
# note: only return hits in this map KO01100
PerformKOEnrichAnalysis_List <- function(mbSetObj, file.nm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # prepare for the result table
  set.size <- length(current.geneset);
  res.mat <- matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat) <- names(current.geneset);
  colnames(res.mat) <- c("Total", "Expected", "Hits", "Pval", "FDR");
  
  # prepare query
  ora.vec <- NULL;
  exp.vec <- mbSetObj$analSet$data[,1]; # drop dim for json
  ora.vec <- names(exp.vec);
  
  # need to cut to the universe covered by the pathways, not all genes
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  #ora.nms <- ora.nms[hits.inx];
  q.size <- length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.geneset, function(x) {
    ora.vec[ora.vec%in%unlist(x)];});
  names(hits.query) <- names(current.geneset);
  
  hit.num <- unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
  
  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.geneset, length));
  
  res.mat[,1] <- set.size;
  res.mat[,2] <- q.size*(set.size/uniq.count);
  res.mat[,3] <- hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4] <- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx <- order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    imp.inx <- res.mat[,4] <= 0.05;
    
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }
  Save2KEGGJSON(hits.query, res.mat, file.nm);
  return(.set.mbSetObj(mbSetObj));
}

# Utility function
# for KO01100
#'@import RJSONIO
Save2KEGGJSON <- function(hits.query, res.mat, file.nm){
  
  load_rjsonio();
  
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat);
  current.msg <<- "Functional enrichment analysis was completed";
  
  if(!exists("ko.edge.map")){
    
    if(.on.public.web){
      ko.edge.path <- paste("../../lib/ko/ko_edge.csv", sep="");
    }else{
      ko.edge.path <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/ko/ko_edge.csv", sep="");
    }
    ko.edge.map <- .readDataTable(ko.edge.path);
    ko.edge.map <- ko.edge.map[ko.edge.map$net=="ko01100",];  #only one map
    ko.edge.map <<- ko.edge.map;
  }
  
  hits.edge <- lapply(hits.query, function(x) {
    as.character(unique(ko.edge.map$edge[ko.edge.map$gene%in%unlist(x)]));});
  
  # only keep hits with edges in the map
  hits.inx <- unlist(lapply(hits.edge, length))>0;
  hits.query <- hits.query[hits.inx];
  resTable <- resTable[hits.inx, ];
  
  # write json
  fun.pval = resTable$Pval; if(length(fun.pval) ==1) { fun.pval <- matrix(fun.pval) };
  hit.num = resTable$Hits; if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  fun.ids <- as.vector(current.setids[names(hits.query)]); if(length(fun.ids) ==1) { fun.ids <- matrix(fun.ids) };
  
  json.res <- list(hits.query = hits.query,
                   hits.edge = hits.edge,
                   path.id = fun.ids,
                   fun.pval = fun.pval,
                   hit.num = hit.num);
  
  json.mat <- RJSONIO::toJSON(json.res, .na='null');
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();
  
  # write csv
  fast.write(resTable, file=paste(file.nm, ".csv", sep=""), row.names=F);
}

# Utility function
SetCurrentMetaData <- function(meta.data, type){
  selected.meta.data <<- meta.data
  enrich.type <<- type;
}

# Given a data with duplicates, dups is the one with duplicates
RemoveDuplicates <- function(data, lvlOpt, quiet=T){

  all.nms <- rownames(data);
  colnms <- colnames(data);
  dup.inx <- duplicated(all.nms);
  dim.orig  <- dim(data);
  data <- apply(data, 2, as.numeric); # force to be all numeric
  dim(data) <- dim.orig; # keep dimension (will lost when only one item)
  rownames(data) <- all.nms;
  colnames(data) <- colnms;
    
  if(sum(dup.inx) > 0){
    uniq.nms <- all.nms[!dup.inx];
    uniq.data <- data[!dup.inx,,drop=F];

    dup.nms <- all.nms[dup.inx];
    uniq.dupnms <- unique(dup.nms);
    uniq.duplen <- length(uniq.dupnms);

    for(i in 1:uniq.duplen){
      nm <- uniq.dupnms[i];
      hit.inx.all <- which(all.nms == nm);
      hit.inx.uniq <- which(uniq.nms == nm);

      # average the whole sub matrix
      if(lvlOpt == "mean"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, mean, na.rm=T);
      }else if(lvlOpt == "median"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, median, na.rm=T);
      }else if(lvlOpt == "max"){
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, max, na.rm=T);
      }else{ # sum
        uniq.data[hit.inx.uniq, ]<- apply(data[hit.inx.all,,drop=F], 2, sum, na.rm=T);
      }
    }
        
    if(!quiet){
      current.msg <<- paste(current.msg, paste("A total of ", sum(dup.inx), " of duplicates were replaced by their ", lvlOpt, ".", sep=""), collapse="\n");
    }
    return(uniq.data);
  }else{
    if(!quiet){
      current.msg <<- paste(current.msg, "All IDs are unique.", collapse="\n");
    }
    return(data);
  }
}

#############################
##### Utility Functions #####
#############################

# return matched KO in the same order (NA if no match)
doKOFiltering <- function(ko.vec, type){
  
  ko.dic <- .read.microbiomeanalyst.lib.rds("ko_dic.rds", "ko");
  
  if(type == "ko"){
    hit.inx <- match(ko.vec, ko.dic$KO);
    return(ko.dic$KO[hit.inx]);
  }
  
  if(type == "ec"){
    hit.inx = vector(mode='numeric', length=length(ko.vec)); # record hit index, initial 0
    match.values = rep(NA, length=length(ko.vec)); # the best matched values (hit names), initial ""
    
    if(substring(ko.vec[1], 0, 2) == "EC"){
      ko.vec <- gsub("^EC", "", ko.vec)
    }
    
    # first find exact match to the common compound names
    hit.inx <- match(ko.vec, ko.dic[, "EC"]);
    match.values <- ko.dic$KO[hit.inx];
    todo.inx <-which(is.na(hit.inx));
    
    if(length(todo.inx) > 0){
      multi.ec.inx <- grep("; ", ko.dic$EC);
      m.dic <- ko.dic[multi.ec.inx,];
      ec.list <- strsplit(m.dic$EC, "; ");
      ecs2ko <- m.dic$KO;
      
      for(i in 1:length(ec.list)){
        hitInx <- match(ko.vec[todo.inx], ec.list[[i]]);
        hitPos <- which(!is.na(hitInx));
        orig.inx<-todo.inx[hitPos];
        match.values[orig.inx] <- ecs2ko[i];
      }
    }
    return(match.values);
  }
  
  if(type == "cog"){
    hit.inx = vector(mode='numeric', length=length(ko.vec)); # record hit index, initial 0
    match.values = rep(NA, length=length(ko.vec)); # the best matched values (hit names), initial ""
    # first find exact match to the common compound names
    hit.inx <- match(ko.vec, ko.dic[, "COG"]);
    match.values <- ko.dic$KO[hit.inx];
    todo.inx <-which(is.na(hit.inx));
    
    if(length(todo.inx) > 0){
      multi.cog.inx <- grep(";", ko.dic$COG);
      m.dic <- ko.dic[multi.cog.inx,];
      cog.list <- strsplit(m.dic$COG, ";");
      cog2ko <- m.dic$KO;
      
      for(i in 1:length(cog.list)){
        hitInx <- match(ko.vec[todo.inx], cog.list[[i]]);
        hitPos <- which(!is.na(hitInx));
        orig.inx<-todo.inx[hitPos];
        match.values[orig.inx] <- cog2ko[i];
      }
    }
    return(match.values);
  }
}