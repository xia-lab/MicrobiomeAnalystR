##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

GetBetaDiversityStats<-function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(microSetObj$analSet$stat.info);
}

GetStressNMDS<-function(microSetObj){
  microSetObj <- .get.microSetObj(microSetObj);
  return(microSetObj$analSet$beta.stress);
}

#'Perform alpha diversity
#'@description This functions performs alpha diversity.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PerformAlphaDiversityComp<-function(microSetObj, opt, metadata){
  
  microSetObj <- .get.microSetObj(microSetObj); 
    
  data <- microSetObj$analSet$alpha;
  cls <- data[,metadata];
  x <- data$value;
  stat.info <- NULL;
    
  if(length(levels(cls)) > 2){
    if(opt=="tt"){
      res <- anova(aov(x ~ cls));
      stat.info <- paste("p-value: ", signif(res$"Pr(>F)"[1], 5), "; [ANOVA] F-value: ", signif(res$"F value"[1], 5), sep="");
    }else{
      res <- kruskal.test(x ~ cls);
      stat.info <- paste("p-value: ", signif(res$p.value, 5), "; [Kruskal-Wallis] statistic: ", signif(res$statistic, 5) , sep="");
    }
  }else{
    inx1 <- which(cls==levels(cls)[1]);
    inx2 <- which(cls==levels(cls)[2]);
        
    if(opt=="tt"){
      res <- t.test(x[inx1], x[inx2]);
      stat.info <- paste("p-value: ", signif(res$p.value, 5), "; [T-test] statistic: ", signif(res$statistic, 5), sep="");
    }else{
      res <- wilcox.test(x[inx1], x[inx2]);
      stat.info <- paste("p-value: ", signif(res$p.value, 5), "; [Mann-Whitney] statistic: ", signif(res$statistic, 5), sep="");
    }
  }
  
  microSetObj$analSet$alpha.stat.info <- stat.info;
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(stat.info);
  }else{
    print(stat.info);
    return(.set.microSet(microSetObj))
  }
}

#######################################
###########Core Microbiome#############
#######################################

#'Perform core microbiome analysis
#'@description This functions performs core microbiome analysis.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import RColorBrewer
CoreMicrobeAnalysis<-function(microSetObj, imgName, preval, detection, taxrank,
                              palette, viewOpt, format="png", dpi=72, width=NA){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  #library("tidyr");
  data <- microSetObj$dataSet$proc.phyobj;
    
  if(taxrank=="OTU"){
    data <- otu_table(microSetObj$dataSet$proc.phyobj,taxa_are_rows=T);
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_first(microSetObj$dataSet$proc.phyobj, taxrank);
    tax_table <- tax_table(data);
    nm <- as.character(tax_table(data)[,taxrank]);
    y <- which(is.na(nm)==TRUE);
    #converting NA values to unassigned
    nm[y] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    #all NA club together
    data1 <- as.matrix(t(sapply(by(data1, rownames(data1), colSums), identity)));
    data <- otu_table(data1, taxa_are_rows=T);
  }
    
  data.compositional <- transform_sample_counts(data,function(x) x / sum(x));
  data.core <- core(data.compositional, detection = detection, prevalence = preval);
  core.nm <- data.frame(prevalence(data.compositional, detection = detection, sort = TRUE));
  colnames(core.nm)[1] <- "Prevelance";
  fileName <- "core_microbiome.csv";
  write.csv(core.nm,file=fileName);
  
  imgName = paste(imgName, ".", format, sep="");
  microSetObj$imgSet$core <- imgName;

  #if more than 1500 features will be present;subset to most abundant=>1500 features.
  #OTUs already in unique names;
  if(ntaxa(data.core)>1500){
    data.core = prune_taxa(names(sort(taxa_sums(data.core), TRUE))[1:1500], data.core);
    viewOpt == "overview";
  }
    
  #setting the size of plot
  if(is.na(width)){
    minW <- 800;
    myW <- 10*18 + 200;
    if(myW < minW){
      myW <- minW;
    }
    w <- round(myW/72,2);
  }
    
  myH <- nrow(data.core)*18 + 150;
  h <- round(myH/72,2);
    
  if(viewOpt == "overview"){
    if(is.na(width)){
      if(w >9.3){
        w <- 9.3;
      }
    }
      if(h > w){
        h <- w;
      }
  }
    
  Cairo::Cairo(file=imgName, unit="in",width=w, height=h, type=format, bg="white",dpi=dpi);
    
  # set up colors for heatmap
  if(palette=="gbr"){
    colors <- grDevices::colorRampPalette(c("green", "black", "red"), space="rgb")(10);
  }else if(palette == "heat"){
    colors <- heat.colors(10);
  }else if(palette == "topo"){
    colors <- topo.colors(10);
  }else if(palette == "gray"){
    colors <- grDevices::colorRampPalette(c("grey90", "grey10"), space="rgb")(10);
  }else{
    
    if(.on.public.web){
      load_rcolorbrewer();
    }
    
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(10));
  }
    
  p <- plot_core(data.core,plot.type = "heatmap",colours = colors,prevalences = seq(.05, 1, .05),detections = 10^seq(log10(1e-3), log10(0.2), length = 10)) + xlab("Detection Threshold (Relative Abundance (%))")+guides(fill = guide_legend(keywidth = 1.5, keyheight = 1));
  print(p);
  microSetObj$analSet$core<-as.matrix(core.nm);
  microSetObj$analSet$core.taxalvl<-taxrank;
  dev.off();
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

##############Adapted from the Microbiome R package#########################
core<-function(x, detection, prevalence, include.lowest=FALSE) {
    xorig <- x
    # TODO: add optional renormalization such that the core member
    # abundances would
    # sum up to 1 ?
    taxa <- core_members(x, detection, prevalence,
                       include.lowest=include.lowest)
    prune_taxa(taxa, xorig)
}

core_members<-function(x, detection=detection, prevalence=prevalence,
                         include.lowest=FALSE) {
  # Pick taxa x samples matrix
    x <- abundances(x)
    if (include.lowest) {
        taxa <- names(which(prevalence(x, detection,
                                   include.lowest=include.lowest) >= prevalence))
    } else {
    taxa <- names(which(prevalence(x, detection,
                                   include.lowest=include.lowest) > prevalence))
    }
    taxa
}

abundances<-function(x, transform="identity") {
  # Pick the OTU data
    if (class(x) == "phyloslim") {
    # Pick OTU matrix
        otu <- get_taxa(x)
    # Ensure that taxa are on the rows
        if (!taxa_are_rows(x) && ntaxa(x) > 1 && nsamples(x) > 1) {
            otu <- t(otu)
        }
        if (ntaxa(x) == 1) {
            otu <- matrix(otu, nrow=1)
            rownames(otu) <- taxa(x)
            colnames(otu) <- sample_names(x)
        }
        if (nsamples(x) == 1) {
            otu <- matrix(otu, ncol=1)
            rownames(otu) <- taxa(x)
            colnames(otu) <- sample_names(x)
        }
    } else if (is.vector(x)) {
        otu <- as.matrix(x, ncol=1)
    } else {
        # If x is not vector or phyloslim object then let us assume it is a
        # taxa x samples
        # count matrix
        otu <- x
    }
  # Apply the indicated transformation
  if (!transform == "identity") {
        otu <- transform(otu, transform)
  }
  otu
}

plot_core<-function(x, prevalences=seq(.1, 1, 0.1), detections=20,
                      plot.type="lineplot", colours=NULL, # gray(seq(0, 1, length=5)),
                      min.prevalence=NULL, taxa.order=NULL, horizontal=FALSE) {
    if (length(detections) == 1) {
     detections <- 10^seq(log10(0.001), log10(max(abundances(x),
                                                 na.rm=TRUE)), length=detections)
    }
    if (plot.type == "lineplot") {
        # Calculate the core matrix (prevalences x abundance thresholds)
        coremat <- core_matrix(x, prevalences, detections)
        res <- core_lineplot(coremat)
    } else if (plot.type == "heatmap") {
        # Here we use taxon x abundance thresholds table indicating prevalences
        res <- core_heatmap(abundances(x),
                        dets=detections, cols=colours,
                        min.prev=min.prevalence, taxa.order=taxa.order)
    }
    p <- res$plot;
    if (horizontal) {
        p <- p + coord_flip() + theme(axis.text.x=element_text(angle=90))
    }
    p
}

#'@import data.table
core_heatmap<-function(x, dets, cols, min.prev, taxa.order){
    
  if(.on.public.web){
    load_datatable();
  }
  
    data <- x
    #colours <- gray(seq(0, 1, length=5)),
    DetectionThreshold <- Taxa <- Prevalence <- NULL
    # Prevalences with varying dets
    prev <- lapply(dets, function(th) {
        prevalence(data, detection=th)
    })
    prev <- do.call("cbind", prev)
    colnames(prev) <- as.character(dets)

    # # Exclude rows and cols that never exceed the given prevalence
    if (!is.null(min.prev)) {
        prev <- prev[rowMeans(prev > min.prev) > 0,
                   colMeans(prev > min.prev) > 0]
    }

    df <- as.data.frame(prev)
    df$ID <- rownames(prev)
    df <- melt(df, "ID")
    names(df) <- c("Taxa", "DetectionThreshold", "Prevalence")
    df$DetectionThreshold <- as.numeric(as.character(df$DetectionThreshold))
    df$Prevalence <- as.numeric(as.character(df$Prevalence))
    if (is.null(taxa.order)) {
        o <- names(sort(rowSums(prev)))
    } else {
        o <- taxa.order
    }
    df$Taxa <- factor(df$Taxa, levels=o)
    theme_set(theme_bw(10))
    p <- ggplot(df, aes(x=DetectionThreshold, y=Taxa, fill=Prevalence))
    p <- p + geom_tile()
    p <- p + xlab("Detection Threshold")
    p <- p + scale_x_log10()
    if (!is.null(cols)) {
        p <- p + scale_fill_gradientn("Prevalence",
                breaks=seq(from=0, to=1, by=0.1),
                colours=cols,
                limits=c(0, 1))
    }
    return(list(plot=p, data=df))
}

prevalence <- function(x, detection=0, sort=FALSE, count=FALSE,
    include.lowest=FALSE) {
    if (is.null(detection)) {
        detection <- (-Inf)
    }
    if (is.null(x)) {
        warning("x is NULL - returning NULL")
        return(NULL)
    }
    # Convert phyloslim to matrix
    if (class(x) == "phyloslim") {
        x <- abundances(x)
    }
    if (is.vector(x)) {
        if (include.lowest) {
            prev <- sum(x >= detection)
        } else {
            prev <- sum(x > detection)
        }
    } else if (is.matrix(x) || is.data.frame(x)) {

        if (include.lowest) {
            prev <- rowSums(x >= detection)
        } else {
            prev <- rowSums(x > detection)
        }
    }
    if (!count) {
        prev <- prev/prevalence_nsamples(x)
    }
    if (sort) {
        prev <- rev(sort(prev))
    }
    prev
}

# Internal auxiliary function
prevalence_nsamples <- function(x) {
    if (is.vector(x)) {
        n <- length(x)
    } else if (is.matrix(x) || is.data.frame(x)) {
        n <- ncol(x)
    }
    n
}

#'Main function to plot pie graphs of microbiome data.
#'@description This functions plots pie graphs of microbiome data. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape

PlotOverallPieGraph<-function(microSetObj, taxalvl, feat_cnt, calcmeth){

  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }
   
  set.seed(28053448);

  #using filtered data
  data <- microSetObj$dataSet$filt.data;
   
  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data <- otu_table(data, taxa_are_rows =TRUE);
  }
   
  #data <- prune_taxa(taxa_sums(dataSet$filt.data)>feat_cnt,dataSet$filt.data);
  sample_table <- sample_data(microSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  datapie <<- merge_phyloslim(data, tax_table(microSetObj$dataSet$proc.phyobj), sample_table);
   
  #using reduce names
  data <- otu_table(datapie);
  data <- data.frame(data);
  data_tax <- tax_table(datapie);

  #reshaping data
  data <- t(data);
  taxa_nm <- data.matrix(data_tax[,taxalvl]);
  y <- which(is.na(taxa_nm)==TRUE);
   
  #converting NA values to unassigned
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data.frame(data %*% sapply(unique(nms),"==",nms));
  colnames(data)<- gsub("\\."," ",colnames(data));
   
  if(calcmeth=="sum"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
  }else{
    dt <- apply(data,2,median);
    ind <- which(dt>feat_cnt);
    ind1 <- which(dt<feat_cnt);
  }
   
  if(length(ind)==0){
    current.msg<<-"All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
    return(0);
  }
   
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)));
  }
   
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  piedata <- aggregate(. ~variable , data=data[-1], FUN=sum);

  # order by abundance
  ord.inx <- order(piedata$value, decreasing = TRUE);
  piedata <<- piedata[ord.inx,];
  piedata$percentage <- round((piedata$value / sum(piedata$value) * 100), digits = 2)
  colnames(piedata) <- c("Taxa", "Abundance", "Percentage")
  write.csv(piedata, "piechart_abundances.csv")
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

# getter
GetPieTaxaNames<- function(){
    res<-as.character(piedata$variable);
    return(res);
}

# getter
GetPieTaxaAbund<- function(){
    return(piedata$value);
}

#'Main function to plot pie graphs of microbiome data.
#'@description This functions plots pie graphs of microbiome data. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotGroupPieGraph <- function(microSetObj, taxalvl, metadata, clslevel,
                              feat_cnt, calcmeth){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }

  set.seed(28053443);

  #using filtered data
  data <- microSetObj$dataSet$filt.data;
    
  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data<-otu_table(data,taxa_are_rows =TRUE);
  }
    
  smpl <- data.frame(sample_data(microSetObj$dataSet$proc.phyobj));
  smpl1 <- sample_data(subset(smpl,smpl[metadata]==clslevel));
  data <- prune_samples(sample_names(smpl1), data);
  data <- merge_phyloslim(data, smpl1);
  microSetObj$dataSet$taxa_table <- tax_table(microSetObj$dataSet$proc.phyobj);
  datapie <<- merge_phyloslim(data, microSetObj$dataSet$taxa_table);

  #using reduce names
  data <- t(data.frame(otu_table(datapie)));
  data_tax <- tax_table(datapie);
  taxa_nm <- data.matrix(data_tax[,taxalvl]);

  #converting NA values to unassigned
  y <- which(is.na(taxa_nm)==TRUE);
  taxa_nm[y] <- "Not_Assigned";
  nms <- colnames(data) <- taxa_nm[,1];
  data <- data.frame(data %*% sapply(unique(nms),"==",nms));
  colnames(data) <- gsub("\\."," ",colnames(data));
    
  if(calcmeth=="sum"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
  }else{
    dt <- apply(data,2,median);
    ind <- which(dt>feat_cnt);
    ind1 <- which(dt<feat_cnt);
  }
    
  if(length(ind)==0){
    current.msg<<-"All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
    return(0);
  }
    
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    smpl_nm <- rownames(data);
    data <- by(t(data),INDICES=names(data),FUN=colSums);
    data <- as.data.frame(do.call(cbind,as.list(data)));
    rownames(data) <- smpl_nm;
  }
    
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  piedata <- aggregate(. ~variable , data=data[-1], FUN=sum);

  # order by abundance
  ord.inx <- order(piedata$value, decreasing = TRUE);
  piedata <<- piedata[ord.inx,];
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

#'Main function to plot sample-wise pie graphs of microbiome data.
#'@description This functions plots sample-wise pie graphs of microbiome data. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape

PlotSamplePieGraph<-function(microSetObj, taxalvl, smplnm, feat_cnt){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }
  
  set.seed(28053443);
  
  #using filtered data
  data <- microSetObj$dataSet$filt.data;
   
  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data <- otu_table(data,taxa_are_rows =TRUE);
  }
   
  sample_table <- sample_data(microSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data <- merge_phyloslim(data, tax_table(microSetObj$dataSet$proc.phyobj), sample_table);
  data <- prune_samples(smplnm, data);
  microSetObj$dataSet$taxa_table <- tax_table(microSetObj$dataSet$proc.phyobj);
  datapie <<- merge_phyloslim(data, microSetObj$dataSet$taxa_table);

  #using reduce names
  data <- otu_table(datapie);
  data <- data.frame(data);
  data_tax <- tax_table(datapie);

  #reshaping data
  data <- t(data);
  taxa_nm <- data.matrix(data_tax[,taxalvl]);
  y <- which(is.na(taxa_nm)==TRUE);
  #converting NA values to unassigned
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data.frame(data %*% sapply(unique(nms),"==",nms));
  colnames(data) <- gsub("\\."," ",colnames(data));
  ind <- which(colSums(data)>feat_cnt);
  ind1 <- which(colSums(data)<feat_cnt);
    
  if(length(ind)==0){
    current.msg<<-"All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
    return(0);
  }
    
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data<-as.data.frame(t(rowsum(t(data),group = colnames(data))));
  }
   
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  piedata <- aggregate(. ~variable , data=data[-1], FUN=sum);

  # order by abundance
  ord.inx <- order(piedata$value, decreasing = TRUE);
  piedata <<- piedata[ord.inx,];
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

#'Function to plot pie-chart data.
#'@description This functions plots pie charts of microbiome data. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape

PlotDataPieFromPie<-function(microSetObj, taxalvl, metadata, clslevel,
                             taxaposn, lowtaxa){

  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }
  
  set.seed(280534431);
  high_taxa <<- as.character(piedata$variable[taxaposn]);
  
  if(high_taxa=="Not_Assigned"){
    return(0);
  }
   
  lowlvl_nm <<- lowtaxa;

  datataxa <- as.matrix(tax_table(datapie));
  subsettax_table <- tax_table(subset(datataxa,datataxa[,taxalvl]==high_taxa));
  data1 <- prune_taxa(taxa_names(subsettax_table),datapie);
  datapietaxatab <<- data_tax <- tax_table(data1);
  datapietaxa <<- data <- t(data.frame(otu_table(data1)));
  taxa_nm <- data.matrix(data_tax[,lowlvl_nm]);

  #converting NA values to unassigned
  y <- which(is.na(taxa_nm)==TRUE);
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data.frame(data %*% sapply(unique(nms),"==",nms));
   
  if(length(nms)==1){
    colnames(data)<-nms;
  }

  colnames(data) <- gsub("\\."," ",colnames(data));
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  color_var <- levels(factor(data$variable));
  x <- length(color_var);
  x.colors <<- rep(col_vector,length.out=x);

  piedata2 <- aggregate(. ~variable , data=data[-1], FUN=sum);
  # order by abundance
  ord.inx <- order(piedata2$value, decreasing = TRUE);
  piedata2 <<- piedata2[ord.inx,];
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

#'Function to update pie-chart data.
#'@description This functions updates pie chart data. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
UpdatePieData<-function(microSetObj, lowtaxa){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }
  
  set.seed(280534432);
   
  data <- datapietaxa;
  data_tax <- datapietaxatab;
  high_taxa <- high_taxa;
  lowlvl_nm <<- lowtaxa;
  taxa_nm <- as.data.frame(data_tax[,lowlvl_nm]);
  taxa_nm <- as.matrix(taxa_nm);
  y <- which(is.na(taxa_nm)==TRUE);
   
  #converting NA values to unassigned
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data);
   
  if(length(nms)==1){
    colnames(data)<-nms;
  }
   
  colnames(data) <- gsub("\\."," ",colnames(data));
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  color_var <- levels(factor(data$variable));
  x <- length(color_var);
  x.colors <<- rep(col_vector,length.out=x);
  # piedatas is two column stats (variable and value)
  piedata2 <<- aggregate(. ~variable , data=data[-1], FUN=sum);

  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

#'Function to save pie-chart
#'@description This functions saves created pie chart plot. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SavePiechartImg <- function(microSetObj, taxalvl, pieName, format="png", dpi=72) {
    
  microSetObj <- .get.microSetObj(microSetObj);
  set.seed(280);
    
  pieName = paste(pieName,".", format, sep="");
  piedata <- transform(transform(piedata, value=value/sum(value)));
    
  #rownames are still arranged by decending order
  piedataimg <- piedata;
  microSetObj$imgSet$pie <- pieName;

  row.names(piedataimg) <- NULL;
  x.cols <- pie.cols;
    
  # java color code to R color code
  x.cols <- paste("#",x.cols, sep="");
  Cairo(file=pieName, width=630, height=500, type=format, bg="white", dpi=dpi);
    
  box <- ggplot(piedataimg, aes(x="", y = value, fill=reorder(variable,-value))) +
         geom_bar(width = 1, stat = "identity") + theme_bw() +
         coord_polar(theta = "y",direction=-1,start = 4.71239) + scale_fill_manual(values=c(x.cols))+
         geom_text(aes(x=1.6,label = scales::percent(round(value,2))), check_overlap = T,size=3,position = position_stack(vjust = 0.5),color="grey48") +
         theme(legend.position="right",axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), plot.title = element_text(hjust=0.5, face="bold"),legend.text=element_text(color="grey48")) +
         labs(x="", y="",fill="");
  print(box);
  microSetObj$analSet$pie<-piedataimg;
  microSetObj$analSet$pie.taxalvl<-taxalvl;
  dev.off();
  return(.set.microSet(microSetObj))
}

#'Function to create pie-chart plot.
#'@description This functions creates the pie chart plot. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotPiechart <- function(microSetObj, rel_perct, pieName, format="png", dpi=72) {
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  set.seed(28056188); 
  pieName = paste(pieName,".", format, sep="");
  piedata2 <- transform(transform(piedata2, value=value/sum(value)));
  ind <- which(piedata2[,"value"]>rel_perct);
  ind1 <- which(piedata2[,"value"]<rel_perct);
    
  if(length(ind)==0){
    current.msg <<- "All features have lower relative abundance than given minimum abundance. Please lower the cut off for relative abundance.";
    return(0);
  }
    
  if(length(ind1)>0){
    levels(piedata2$variable)[ind1] <- "Others";
    piedata2 <- aggregate(value~variable, piedata2, FUN=sum);
  }
    
  ind_zero <- which(piedata2[,"value"]==0);
    
  if(length(ind_zero)>0){
    piedata2 <- piedata2[-ind_zero,];
  }
    
  Cairo::Cairo(file=pieName,width=630, height=400, type=format, bg="white",dpi=dpi);
    
  box=ggplot(piedata2, aes(x="", y = value, fill=variable)) +
        geom_bar(width = 1, stat = "identity") + theme_bw() +
        coord_polar(theta = "y") + scale_fill_manual(values=c(x.colors))+
        geom_text(aes(x=1.7, label = scales::percent(value)), check_overlap = T,size=3, position = position_stack(vjust = 0.5)) +
        theme(legend.position="bottom",axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), plot.title = element_text(hjust=0.5, face="bold")) +
        labs(x="", y="",fill =lowlvl_nm) + ggtitle(high_taxa);
  print(box);
    
  dev.off();
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

#######################################
###########Alpha-diversity#############
#######################################

#'Function to plot alpha-diversity analysis.
#'@description This functions creates a plot for alpha-diversity. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotAlphaData<-function(microSetObj, data.src, bargraphName, distName,
                        metadata, taxrank, format="png", dpi=72){
  
  microSetObj <- .get.microSetObj(microSetObj);  
  
  set.seed(13133);
    
  if(data.src == "orig"){
    data <- readRDS("orig.phyobj");
  }else{
    data <- microSetObj$dataSet$proc.phyobj;
  }

  if(taxrank!="OTU"){
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data, taxrank);
  }

  bargraphName = paste(bargraphName, ".", format, sep="");
  microSetObj$imgSet$alpha <- bargraphName;

  Cairo::Cairo(file=bargraphName, width=900, height=450, type=format, bg="white",dpi=dpi);
    
  #reordering the sample (getting reordered rownames)
  sam <- sample_data(data);
  sam <- sam[order(sam[[metadata]])];
  smplord <- rownames(sam);

  box = plot_richness(data,color=metadata,measures = distName)+scale_x_discrete(limits=c(smplord));
  microSetObj$analSet$alpha<-box$data;
  write.csv(microSetObj$analSet$alpha, file="alphadiversity.csv");
  box = box + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust= 1));
  box$layers <- box$layers[-1];
  box <- box + geom_point(size=2);
  #getting scale for plot (using same for boxplot also)
  ylimits <<- layer_scales(box)$y$range$range;
  print(box);
    
  dev.off();
  microSetObj$analSet$alpha.meth <- distName;
  microSetObj$analSet$alpha.metadata <- metadata;
  microSetObj$analSet$alpha.taxalvl <- taxrank;
  return(.set.microSet(microSetObj))
}

#'Function to create bar plots of selected taxa level.
#'@description This functions creates bar plots of a selected taxa level. 
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotTaxaAlphaBarSam<-function(microSetObj, barplotName, taxalvl, samplnm,
                              imgOpt, feat_cnt, format="png", dpi=72){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }

  #using filtered data
  data <- microSetObj$dataSet$filt.data;
    
  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data<-otu_table(data, taxa_are_rows =TRUE);
  }
    
  sample_table <- sample_data(microSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data1 <- merge_phyloslim(data, tax_table(microSetObj$dataSet$proc.phyobj), sample_table);

  yLbl <- "Actual Abundance";

  data <- as.matrix(otu_table(data1));
  data <- data[,samplnm,drop=FALSE];
  data <- t(data);
  data_tax <- tax_table(data1);
    
  #reshaping data
  taxa_nm <- as.data.frame(data_tax[,taxalvl]);
  taxa_nm <- as.matrix(taxa_nm);
  y <- which(is.na(taxa_nm)==TRUE);

  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data);
  data <- data[ , order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  colnames(data)[indx] <- "NA";
  ind <- which(colSums(data)>feat_cnt);
  ind1 <- which(colSums(data)<feat_cnt);
    
  if(length(ind)==0){
    current.msg <<- "All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
    return(0);
  }
    
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data <- as.data.frame(t(rowsum(t(data),group = colnames(data))));
  }
    
  if(imgOpt=="barnorm"){
    data <- as.data.frame(apply(data,1, function(x) x/sum(x)));
    data <- as.data.frame(t(data));
    yLbl <- "Relative Abundance";
  }
    
  feat_no<-ncol(data);
  #adjust height according to number of legends
    
  w<-600;
    
  if(feat_no < 10){
    w<-w;
  } else if (feat_no < 20){
    w<-w+100;
  } else if (feat_no < 50){
    w<-w+200;
  } else if (feat_no < 100){
    w<-w+400;
  } else if (feat_no > 100){
    w<-w+500;
  }
  
  write.csv(t(data), file="taxa_abund.csv");
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  data <- data[order(data[,2]),];
  data <- data[,-1];
  a <- feat_no;
    
  if(length(a)<50){
    h<-feat_no*50;
  }else{
    h<-feat_no*25;
  }
    
  #sorting by descending
  barplotName = paste(barplotName, ".",format, sep="");
  microSetObj$imgSet$stack<-barplotName;

  Cairo::Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
  box<-ggplot(data, aes(x=reorder(variable,value),y=value))+geom_bar(stat="identity",width=0.6,fill="steelblue")+theme_bw()+
        theme(axis.text.x = element_text(angle = 0,vjust=0.5))+
        labs(y=yLbl,x="",fill=taxalvl)+coord_flip()+
        theme(panel.background = element_blank(),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none");
  print(box);
  dev.off();
  microSetObj$analSet$stack<-data;
  microSetObj$analSet$stack.taxalvl<-taxalvl;
  microSetObj$analSet$plot<-"Stacked Bar";

  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

###########Beta-diversity#############
#######################################

#'Function to plot beta diversity.
#'@description This functions creates beta diversity plots.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import data.table
#'@import ape
PlotBetaDiversity<-function(microSetObj, plotNm, ordmeth, distName, colopt, metadata, 
                            showlabel, taxrank, taxa, alphaopt, ellopt, format="png", dpi=72){

  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_datatable();
  }
  
  set.seed(13134);
  metadata <- metadata;
    
  #using normalized data
  if(taxrank=="OTU"){
    taxa_table <- tax_table(microSetObj$dataSet$proc.phyobj);
    data <- merge_phyloslim(microSetObj$dataSet$norm.phyobj, taxa_table);
  }else{
    taxa_table <- tax_table(microSetObj$dataSet$proc.phyobj);
    data <- merge_phyloslim(microSetObj$dataSet$norm.phyobj, taxa_table);
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data,taxrank);
  }
    
  if(colopt=="taxa"){
    if(taxrank=="OTU"){
      data1 <- as.matrix(otu_table(data));
      feat_data <- as.numeric(data1[taxa,]);
    }else{
      nm <- as.character(tax_table(data)[,taxrank]);
      #converting NA values to unassigned
      nm[is.na(nm)] <- "Not_Assigned";
      data1 <- as.matrix(otu_table(data));
      rownames(data1) <- nm;
      #all NA club together
      data1 <- as.matrix(t(sapply(by(data1, rownames(data1), colSums), identity)));
      feat_data <- data1[taxa,];
    }
    sample_data(data)$taxa <- feat_data;
    indx <- which(colnames(sample_data(data))=="taxa");
    colnames(sample_data(data))[indx] <- taxa;
    taxa1 <- colnames(sample_data(data))[indx];
    taxaorig <- taxa;
    #if the taxa names are numeric then X is appending to column name
        
    if(!is.na(as.numeric(taxa))=="TRUE"){
      taxa <- paste("X",taxa, sep = "", collapse = NULL);
    }
  }else if(colopt=="alphadiv") {
    data1 <- microSetObj$dataSet$proc.phyobj;
    box <- plot_richness(data1,measures = alphaopt);
    alphaboxdata <- box$data;
    sam_nm <- sample_names(data);
    alphaboxdata <- alphaboxdata[alphaboxdata$samples %in% sam_nm,];
    alphaval <- alphaboxdata$value;
    sample_data(data)$alphaopt <- alphaval;
    indx <- which(colnames(sample_data(data))=="alphaopt");
    colnames(sample_data(data))[indx] <- alphaopt;
  }else{
    data<-data;
  }

  if(distName=="wunifrac"){
    
    if(.on.public.web){
      load_ape();
    }
    
    pg_tree <- readRDS("tree.RDS");
    pg_tb <- tax_table(data);
    pg_ot <- otu_table(data);
    pg_sd <- sample_data(data);
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
    data <- merge_phyloslim(pg_tb, pg_ot, pg_sd, pg_tree);

    if(!ape::is.rooted(phy_tree(data))){
      pick_new_outgroup <- function(tree.unrooted){
        treeDT <- cbind(cbind(data.table(tree.unrooted$edge),
                  data.table(length = tree.unrooted$edge.length))[1:Ntip(tree.unrooted)],
                  data.table(id = tree.unrooted$tip.label));
                  new.outgroup <- treeDT[which.max(treeDT$length), ]$id
                  return(new.outgroup);
      }
      new.outgroup <- pick_new_outgroup(phy_tree(data));
      phy_tree(data) <- ape::root(phy_tree(data),
                        outgroup = new.outgroup,
                        resolve.root=TRUE)
      }
        
    saveRDS(data, "data_unifra.RDS");
    ord <- phyloslimR::ordinate(data,method = ordmeth,"unifrac",weighted=TRUE);

  } else if (distName=="unifrac") {
      
    if(.on.public.web){
      load_ape();
    }        
      
    pg_tree <- readRDS("tree.RDS");
    pg_tb <- tax_table(data);
    pg_ot <- otu_table(data);
    pg_sd <- sample_data(data);
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
    data <- merge_phyloslim(pg_tb, pg_ot, pg_sd, pg_tree);

    if(!ape::is.rooted(phy_tree(data))){
        pick_new_outgroup <- function(tree.unrooted){
        treeDT <- cbind(cbind(data.table(tree.unrooted$edge),
        data.table(length = tree.unrooted$edge.length))[1:Ntip(tree.unrooted)],
        data.table(id = tree.unrooted$tip.label));
        new.outgroup <- treeDT[which.max(treeDT$length), ]$id
        return(new.outgroup);}

        new.outgroup <- pick_new_outgroup(phy_tree(data));
        phy_tree(data) <- ape::root(phy_tree(data),
                          outgroup = new.outgroup,
                          resolve.root=TRUE);
    }
    
    saveRDS(data, "data_unifra.RDS");
    ord<-phyloslimR::ordinate(data,method = ordmeth,"unifrac");

  }else{
    ord<-phyloslimR::ordinate(data,method = ordmeth,distName);
  }
  
  plotNm = paste(plotNm, ".", format, sep="");
  microSetObj$imgSet$beta2d<-plotNm;

  Cairo::Cairo(file=plotNm, width=720, height=500, type=format, bg="white",dpi=dpi);
    
  if(colopt=="taxa"){
    box = phyloslimR::plot_ordination(data,ord,color=taxa)+ labs(aesthetic=taxaorig)+scale_colour_gradient(low="green", high="red");
  }else if(colopt=="alphadiv") {
    box = phyloslimR::plot_ordination(data,ord,color=alphaopt)+ scale_colour_gradient(low="green", high="red");
  }else{
    box = phyloslimR::plot_ordination(data,ord,color=metadata);
  }
    
  box$layers <- box$layers[-1];
    
  if(showlabel=="samnm"){
    box = box + geom_text(aes(label=sample_id), hjust=0.5, vjust=2, size=3, fontface="bold");
    box = box + geom_point(size=4, alpha=0.6) + theme_bw();
  }else if(showlabel=="none"){
    box=box+geom_point(size=4, alpha=0.8) + theme_bw();
  }else{
    showlabel <<- showlabel;
    bx_data <<- data.frame(box$data);
    box = box + geom_text(aes(label=bx_data[ ,showlabel]), hjust=0.5, vjust=2, size=3, fontface="bold");
    box = box + geom_point(size=4, alpha=0.6) + theme_bw();
  }
  
  #used for area color for ellipse
  if(colopt=="expfac"){
    sam_data <- as.data.frame(sample_data(data));
    clsLbl <- sam_data[[metadata]];

    if (ellopt=="yes"){
      box = box + stat_ellipse(type="norm", linetype=2, geom = "polygon", alpha = 0.2, aes_string(fill = clsLbl), show.legend=FALSE);
    }
  }

  print(box);
  dev.off();

  #saving info for report generation
  microSetObj$analSet$beta <- data;
  microSetObj$analSet$beta.meth <- ordmeth;
  microSetObj$analSet$beta.dist <- distName;
  microSetObj$analSet$beta.taxalvl <- taxrank;

  # if it is NMDS, show stress value in the plot
  if(ordmeth == "NMDS"){
    microSetObj$analSet$beta.stress <- paste("[NMDS] Stress =", signif(ord$stress, 5));
  }
  return(.set.microSet(microSetObj))
}

#'Plot functional annotation summary
#'@description This functions plots the functional annotation summary.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotFunAnotSummary<-function(microSetObj, imgName, format="png", dpi=72){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  set.seed(280561499);
    
  if(is.null(microSetObj$analSet$func.pred)){
    result <- readRDS("func.pred");
  }else{
    result <- microSetObj$analSet$func.pred;
    saveRDS(microSetObj$analSet$func.pred, file="func.pred");
    microSetObj$analSet$func.pred <- NULL;
  }
    
  imgName = paste(imgName, ".",format, sep="");
  microSetObj$imgSet$func.pred <- imgName;
  Cairo::Cairo(file=imgName, width=900, height=480, type=format, bg="white",dpi=dpi);
  box <- ggplot(stack(log(result)), aes(x = factor(ind, levels = names(result)), y=values)) + labs(x=NULL, y="log (KO Counts)") + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(box);
  dev.off();
  
  return(.set.microSet(microSetObj))
}

#'Plot functional annotation summary
#'@description This functions plots the functional annotation summary.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotTaxaAlphaArea<-function(microSetObj, barplotName, viewOpt, taxalvl, metadata,
                            imgOpt, feat_cnt, colpalopt, calcmeth, format="png", dpi=72){

  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }

  #using filtered data
  data <- microSetObj$dataSet$filt.data;

  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data <- otu_table(data,taxa_are_rows =TRUE);
  }

  sample_table <- sample_data(microSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data1 <- merge_phyloslim(data, tax_table(microSetObj$dataSet$proc.phyobj), sample_table);

  if(viewOpt=="smpl_grp"){
    data <- as.data.frame(t(otu_table(data1)));
    data <- cbind(data,variable=data1@sam_data[[metadata]]);
    data <- aggregate(. ~variable,data,sum);
    gp_nm <- rownames(data)<-data[,1];
    data <- data[,-1];
    data <- data[order(rownames(data)),];
    clsLbl <- sort(unique(factor(data1@sam_data[[metadata]])));
    colvec <- NULL;
  }else {
    smpl_nm <- sample_names(data1);
    data <- data.frame(otu_table(data1));

    # reorder data based on groups
    sam <- sample_data(data1);
    clsLbl <- factor(sam[[metadata]]);
    ord.inx <- order(clsLbl);
    smpl_nm <- smpl_nm[ord.inx];
    clsLbl <- clsLbl[ord.inx];
    colvec <- as.numeric(clsLbl)+1;
    data <- t(data[,ord.inx]);
  }

  data_tax <- tax_table(data1);
  #reshaping data
  taxa_nm <- as.data.frame(data_tax[,taxalvl]);
  taxa_nm <- as.matrix(taxa_nm);
  y <- which(is.na(taxa_nm)==TRUE);

  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- as.matrix(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data);
  data <- data[ ,order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  colnames(data)[indx] <- "NA";

  if(calcmeth=="sum"){
    ind<-which(colSums(data)>feat_cnt);
    ind1<-which(colSums(data)<feat_cnt);
  } else {
    dt<-apply(data,2,median);
    ind<-which(dt>feat_cnt);
    ind1<-which(dt<feat_cnt);
  }
  
  if(length(ind)==0){
    current.msg<<-"All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
    return(0);
  }
  
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data<- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)));
  }
  
  feat_no<-ncol(data);
  #adjust height according to number of legends
  h<-540;
  if(feat_no < 10){
    h<-h;
  } else if (feat_no < 20){
    h<-h+100;
  } else if (feat_no < 50){
    h<-h+200;
  } else if (feat_no < 100){
    h<-h+400;
  } else if (feat_no > 100){
    h<-h+500;
  }

  # width calculation
  a <- nsamples(data1);
  min.w <- 540;
  if(length(a)<50){
    w <- a*35;
  }else{
    w <- a*25;
  }
  if(w < min.w){
    w <- min.w;
  }

  write.csv(t(data), file="taxa_abund.csv");
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  data <- data[order(data[,2]),];
  data <- data[,-1];
  if(viewOpt=="smpl_grp"){
    data$step <- rep(1:length(gp_nm),feat_no);
    lbl <- gp_nm;
  }else {
    data$step <- rep(1:a,feat_no);
    lbl <- smpl_nm;
  }

  tmp_df <- aggregate(data$value, by=list(data$variable), FUN=mean)
  var_level <- tmp_df[order(tmp_df$x, decreasing = TRUE), ][[1]]
  data$variable <- factor(data$variable,
                          levels = var_level) # change the factor level of taxa

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

  barplotName = paste(barplotName, ".",format, sep="");
  microSetObj$imgSet$stack<-barplotName;
  Cairo::Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
  box <- ggplot(data,aes(x=step,y=value)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust =1,vjust=0.5)) +
    geom_area(aes(fill=variable), position='fill') +
    scale_x_continuous(breaks=seq(1,length(unique(data$step)),1),labels=lbl) +
    labs(y="Relative abundance",x="",fill=taxalvl) +
    theme(legend.position="bottom",legend.box = "vertical") +
    theme(axis.text.x = element_text(colour=colvec),axis.title.x=element_blank());

  if(colpalopt=="set3"){
    cols.needed <- length(unique(data$variable))
    if(cols.needed > 12){
      col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      box <- box + scale_fill_manual(values=col.func(cols.needed),
                                     guide = guide_legend(direction = "horizontal",
                                                          ncol = 5))
    } else {
      box <- box + scale_fill_brewer(palette = "Set3",
                                     guide = guide_legend(direction = "horizontal",
                                                          ncol = 5))
    }
  } else {
    box <- box + scale_fill_manual(values=c(x.colors))
  }

  print(box);
  dev.off();
  microSetObj$analSet$stack<-data;
  microSetObj$analSet$stack.taxalvl<-taxalvl;
  microSetObj$analSet$plot<-"Stacked Area";
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

#'Function to plot bar charts for alpha diversity
#'@description This functions plots bar charts of different taxonomic levels
#'for alpha diversity.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotTaxaAlphaBar<-function(microSetObj, barplotName, taxalvl, metadata, facet, imgOpt, 
                           feat_cnt, colpalopt, calcmeth, format="png", dpi=72){

  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
    load_ggplot();
  }

  data <- microSetObj$dataSet$filt.data;
  
  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data <- otu_table(data, taxa_are_rows =TRUE);
  }
  
  sample_table <- sample_data(microSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data1 <- merge_phyloslim(data, tax_table(microSetObj$dataSet$proc.phyobj), sample_table);

  data <- as.data.frame(otu_table(data1));
  data_tax <- tax_table(data1);

  # reorder data based on groups
  sam <- sample_data(data1);
  
  if(facet == "none"){
    sam$newnewnew <- rep("one", nrow(sam));
    facet <- "newnewnew";
  }
  
  sample_data(data1) <- sam

  metalp <- as(sample_data(data1), "data.frame") #extract metadata table

  smpl_nm <- sample_names(data1);
  clsLbl <- factor(sam[[facet]]);
  
  if(min(table(clsLbl)) < 2){
    current.msg<<-"Too many facets to be displayed - please select a more meaninful facet option with at least 3 samples per group.";
    return(0);
  }

  ord.inx <- order(clsLbl);
  smpl_nm <- smpl_nm[ord.inx];

  clsLbl <- clsLbl[ord.inx];
  colvec <- as.numeric(clsLbl)+1;
  data <- t(data[,ord.inx]);
  #reshaping data
  taxa_nm <- as.data.frame(data_tax[,taxalvl]);
  taxa_nm <- as.matrix(taxa_nm);
  y <- which(is.na(taxa_nm)==TRUE);

  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- as.matrix(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data);
  data <- data[ , order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  colnames(data)[indx] <- "NA";

  if(calcmeth=="sum"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
  } else {
    dt <- apply(data,2,median);
    ind <- which(dt>feat_cnt);
    ind1 <- which(dt<feat_cnt);
  }
  
  if(length(ind)==0){
    current.msg<<-"All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
    return(0);
  }
  
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)));
  }
  
  yLbl <- "Actual Abundance";
  
  if(imgOpt=="barnorm"){
    data <- as.data.frame(apply(data,1, function(x) x/sum(x)));
    data<-as.data.frame(t(data));
    yLbl <- "Relative Abundance";
  }

  # height according to number of legends
  feat_no <- ncol(data);
  h <- 540;
  if(feat_no < 10){
    h <- h;
  } else if (feat_no < 20){
    h <- h+100;
  } else if (feat_no < 50){
    h <- h+200;
  } else if (feat_no < 100){
    h <- h+400;
  } else if (feat_no > 100){
    h <- h+500;
  }

  # width calculation
  a <- nsamples(data1);
  min.w <- 540;
  if(length(a)<50){
    w <- a*35;
  }else{
    w <- a*25;
  }
  if(w < min.w){
    w <- min.w;
  }

  data <- data[row.names(metalp), ]
  data[[get("facet")]] <- metalp[[get("facet")]]
  data$sample <- row.names(data);
  write.csv(t(data), file="taxa_abund.csv");
  data <- melt(data, id = c("sample", get("facet")))
  tmp_df <- aggregate(data$value, by=list(data$variable), FUN=mean)
  var_level <- tmp_df[order(tmp_df$x, decreasing = TRUE), ][[1]]
    as.character()# get the order of taxa relative abundance

  data$variable <- factor(data$variable,
                          levels = var_level) # change the factor level of taxa

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

  barplotName = paste(barplotName, ".",format, sep="");
  microSetObj$imgSet$stack <- barplotName;

  Cairo::Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);

  box <- ggplot(data = data,
                aes(x = sample,
                    y = value,
                    fill = variable)) +
    geom_bar(stat="identity", position="stack") +
    facet_grid(as.formula(paste("~", facet)), scales = "free", space = "free") +
    labs(y=yLbl,x="", fill = taxalvl)+
    theme(legend.position="bottom",legend.box = "vertical",
          axis.text.x = element_text(angle = 45,
                                     hjust =1,vjust=1), #remove colour=colvec
          axis.title.x=element_blank(),
          panel.background = element_rect(fill = "grey90",
                                colour = "grey90",
                                size = 0.5, linetype = "solid"),
          axis.line = element_line(colour = "lightgrey"));
  
  if(colpalopt=="set3"){
    cols.needed <- length(unique(data$variable))
      if(cols.needed > 12){
        col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
        box <- box + scale_fill_manual(values=col.func(cols.needed), guide = guide_legend(direction = "horizontal", ncol = 5))
      }else{
        box <- box + scale_fill_brewer(palette = "Set3", guide = guide_legend(direction = "horizontal", ncol = 5))
      }
  }else{
    box <- box + scale_fill_manual(values=c(x.colors))
  }

  if(facet == "newnewnew"){
    box <- box + theme(strip.text.x = element_blank())
  }

  print(box);
  dev.off();
  microSetObj$analSet$stack <- data;
  microSetObj$analSet$stack.taxalvl <- taxalvl;
  microSetObj$analSet$plot <- "Stacked Bar";
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

#'Function to perform categorical comparison.
#'@description This functions performs categorical comparisons.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PerformCategoryComp <- function(microSetObj, method, distnm, variable){

  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_vegan();
  }

  if(distnm %in% c("wunifrac", "unifrac")) {
    data <- readRDS("data_unifra.RDS");
  } else {
    data <- microSetObj$dataSet$proc.phyobj;
    data <- merge_phyloslim(data);
  }

  data <- transform_sample_counts(data, function(x) x/sum(x));
  data.dist <- phyloslimR::distance(data, method=distnm);
  group <- get_variable(data,variable);
  stat.info <- "";
  resTab <- list();
    
  if(method=="adonis"){
    sampledf <- data.frame(sample_data(data),check.names=F);
    res <- adonis(data.dist~sampledf[ ,variable],data = sampledf);
    resTab <- res$aov.tab[1,];
    stat.info <- paste("[PERMANOVA] F-value: ", signif(resTab$F.Model, 5),  "; R-squared: ", signif(resTab$R2, 5), "; p-value < ", signif(resTab$Pr, 5), sep="");
  }else if(method=="anosim"){
    anosim <- anosim(data.dist,group=group);
    resTab$Rval <- anosim$statistic;
    resTab$pval <- anosim$signif;
    stat.info <- paste("[ANOSIM] R: ", signif(resTab$Rval, 5), "; p-value < ", signif(resTab$pval, 5), sep="");
  }else if (method=="permdisp"){
    beta <- betadisper(data.dist,group=group);
    resTab <- anova(beta);
    stat.info <- paste("[PERMDISP] F-value: ", signif(resTab$"F value"[1], 5), "; p-value: ", signif(resTab$"Pr(>F)"[1], 5), sep="");
  }

  microSetObj$analSet$stat.info <- stat.info;
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}

###################################################################
# generate figure barplot with group#############################
###################################################################

#'Function to plot group-wise bar charts.
#'@description This functions plots group-wise bar charts of a specified 
#'taxa for alpha diversity analysis.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PlotTaxaAlphaBarSamGrp<-function(microSetObj, barplotName, taxalvl, metadata, imgOpt,
                                 feat_cnt, colpalopt, calcmeth, format="png", dpi=72){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  if(.on.public.web){
    load_reshape();
  }
  
  #using filtered data
  data <- microSetObj$dataSet$filt.data;
    
  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data<-otu_table(data,taxa_are_rows =TRUE);
  }
    
  sample_table <- sample_data(microSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data1 <- merge_phyloslim(data, tax_table(microSetObj$dataSet$proc.phyobj), sample_table);

  yLbl <- "Actual Abundance";

  data <- as.data.frame(t(otu_table(data1)));
  data <- cbind(data, variable=data1@sam_data[[metadata]]);
  data <- aggregate(. ~variable,data,sum);
  gp_nm <- rownames(data) <- data[,1];
  data <- data[,-1];
  data <- data[order(rownames(data)),];
  clsLbl <- sort(unique(factor(data1@sam_data[[metadata]])));
  data_tax <- tax_table(data1);
    
  #reshaping data
  taxa_nm <- as.data.frame(data_tax[,taxalvl]);
  taxa_nm <- as.matrix(taxa_nm);
  y <- which(is.na(taxa_nm)==TRUE);

  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- as.matrix(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data);
  data <- data[ , order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  colnames(data)[indx] <- "NA";
    
  if(calcmeth=="sum"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
  } else {
    dt <- apply(data,2,median);
    ind <- which(dt>feat_cnt);
    ind1 <- which(dt<feat_cnt);
  }
    
  if(length(ind)==0){
    current.msg <<- "All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.";
    return(0);
  }
    
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)));
  }

  if(imgOpt=="barnorm"){
    data <- as.data.frame(apply(data,1, function(x) x/sum(x)));
    data <- as.data.frame(t(data));
    yLbl <- "Relative Abundance";
  }
    
  feat_no <- ncol(data);

  #adjust height according to number of legends
  w <- 600;
    
  if(feat_no < 10){
    w<-w;
  } else if (feat_no < 20){
    w<-w+100;
  } else if (feat_no < 50){
    w<-w+200;
  } else if (feat_no < 100){
    w<-w+400;
  } else if (feat_no > 100){
    w<-w+500;
  }
    
  write.csv(t(data), file="taxa_abund.csv");
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  data <- data[order(data[,2]),];
  data <- data[,-1];
  data$step <- rep(gp_nm,feat_no);

  tmp_df <- aggregate(data$value, by=list(data$variable), FUN=mean)
  var_level <- tmp_df[order(tmp_df$x, decreasing = TRUE), ][[1]]
  data$variable <- factor(data$variable, levels = var_level) # change the factor level of taxa
   
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
    x.colors<-rep(custom_col42,length.out=x);
  }

  h<-length(clsLbl)*100;

  barplotName = paste(barplotName, ".",format, sep="");
  microSetObj$imgSet$stack<-barplotName;

  Cairo::Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
  box<-ggplot(data,aes(x = step, y = value, fill = variable))+
            geom_bar(stat="identity", position="stack", width = 0.4)+
            #scale_y_continuous(expand = c(0, 0, 0.3, 0)) +
            theme_bw() + theme(legend.position="bottom",legend.box = "vertical")+
            labs(y=yLbl,x="",fill=taxalvl)+
            theme(axis.text.x = element_text(angle = 0,vjust=0.5))+
            coord_flip()+
            theme(panel.background = element_rect(fill = "grey90",
                                          colour = "grey90",
                                          size = 0.5, linetype = "solid"),
                  panel.border = element_blank())

  if(colpalopt=="set3"){
    cols.needed <- length(unique(data$variable))
    if(cols.needed > 12){
      col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(12, "Set3"))
      box <- box + scale_fill_manual(values=col.func(cols.needed),
                                     guide = guide_legend(direction = "horizontal",
                                                          ncol = 5))
    } else {
      box <- box + scale_fill_brewer(palette = "Set3",
                                     guide = guide_legend(direction = "horizontal",
                                                          ncol = 5))
        }
  } else {
    box <- box + scale_fill_manual(values=c(x.colors))
  }
    
  print(box);
  dev.off();
  microSetObj$analSet$stack<-data;
  microSetObj$analSet$stack.taxalvl<-taxalvl;
  microSetObj$analSet$plot<-"Stacked Bar";
  
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}


###############################################
###########alpha diversity boxplot#############
###############################################

#'Function to create box plots for alpha diversity analysis
#'@description This functions performs metagenome seq analysis on the microbiome data.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotAlphaBoxData<-function(microSetObj, boxplotName, distName, metadata, format="png", dpi=72){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  set.seed(1313397);
  data <- microSetObj$analSet$alpha;
  CLASS <- data[,metadata];
  boxplotName = paste(boxplotName, ".", format, sep="");
  microSetObj$imgSet$alpha.box<-boxplotName;

  Cairo::Cairo(file=boxplotName,width=500, height=400, type=format, bg="white",dpi=dpi);
  
  box1 = ggplot(data,
              aes(data[,metadata],data$value,
                  color = CLASS))+#change fill to color
    stat_boxplot(geom ='errorbar',width=0.2)+
    geom_boxplot(aes(outlier.shape=1),
                 position = position_dodge(width = 0),width=0.3)+
    geom_jitter(width = 0.1) + #reduce the width of jitter
    stat_summary(fun.y=mean, #add mean point
                 geom = "point",
                 shape = 18,
                 size = 4,
                 color = "black")+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 45,vjust = 1, hjust = 1),
          plot.title = element_text(color="black", face="bold", hjust = 0.5)) + #ajust title
    labs(title = "",
         y= paste("Alpha-diversity index", as.character(data$variable[[1]]), sep = " "),
         x="") + #remove x = CLASS, add title name, change y name
    coord_cartesian(ylim = c(ylimits[1],ylimits[2]));
  print(box1);
  dev.off();
  return(.set.microSet(microSetObj))
}

#'Function to create rarefraction curves of microbiome data
#'@description This functions plots rarefraction curves.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotRarefactionCurve<-function(microSetObj, data.src, linecolor, linetype, facet, step, imgName, format="png",dpi=72){

  microSetObj <- .get.microSetObj(microSetObj);
    
  # should use unfiltered data
  if(data.src == "orig"){
    data_rare <- readRDS("orig.phyobj");
  }else{
    data_rare <- microSetObj$dataSet$proc.phyobj;
  }

  #get good's coverage index
  goods_coverage <- ComputeGoods(data_rare)
  write.csv(goods_coverage, "goods_coverage.csv", row.names = FALSE, quote = FALSE);

  feat_no <- nsamples(data_rare);
    
  #adjust height according to number of legends
  w <- 600;
    
  if(feat_no < 10){
    w<-w;
  } else if (feat_no < 20){
    w<-w+100;
  } else if (feat_no < 50){
    w<-w+200;
  } else if (feat_no < 100){
    w<-w+400;
  } else if (feat_no > 100){
    w<-w+500;
  }
    
  imgName = paste(imgName,".", format, sep="");
  Cairo::Cairo(file=imgName, width = 1.5 * w, height = 540, type=format, bg="white", dpi=dpi);
  linecolor <- ifelse(linecolor == "none", "NULL", linecolor);
  linetype <- ifelse(linetype == "none", "NULL", linetype);

  box <- ggrare2(data_rare,
                data.src = data.src,
                color = linecolor,
                label = "Sample",
                linetype = linetype,
                se = FALSE,  # this is not to meaningful
                step = step);

  if(!is.null(facet) & facet != "none"){
    box <- box + facet_wrap(as.formula(paste("~", facet)))
  } else {
    box <- box;
  }
  print(box);
  dev.off();
  return(.set.microSet(microSetObj))
}

# Utility function
# get goods's coverage
ComputeGoods <-function(physeq_object){
  com <- t(get_sample(physeq_object))
  no.seqs <- rowSums(com)
  sing <- com==1
  no.singleton <- apply(sing, 1, sum)
  goods <- 100*(1-no.singleton/no.seqs)
  sample <- row.names(com)
  goods.sum <- cbind(sample, no.singleton, no.seqs, goods)
  goods.sum <- as.data.frame(goods.sum)
  row.names(goods.sum) <- c()
  return(goods.sum)
}

# Utility function that performs rarefaction
ggrare2 <- function(physeq_object, data.src, label = NULL, color = NULL, plot = TRUE, linetype = NULL, se = FALSE, step=5) {

  x <- methods::as(phyloslimR::otu_table(physeq_object), "matrix")
  
  if (phyloslimR::taxa_are_rows(physeq_object)) { x <- t(x) }

    ## This script is adapted from vegan `rarecurve` function
    tot <- rowSums(x)
    S <- rowSums(x > 0)
    nr <- nrow(x)
    step_new = floor(max(tot)/as.integer(step))
    rarefun <- function(i) {
      #cat(paste("rarefying sample", rownames(x)[i]), sep = "\n")
      n <- seq(1, tot[i], by = step_new)
      if (n[length(n)] != tot[i]) {
        n <- c(n, tot[i])
      }
      y <- vegan::rarefy(x[i, ,drop = FALSE], n, se = se)
      if (nrow(y) != 1) {
        rownames(y) <- c(".S", ".se")
        return(data.frame(t(y), Size = n, Sample = rownames(x)[i]))
      } else {
        return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i]))
      }
    }

    f_n <- paste(data.src, step, "rds", sep = ".");

    if(file.exists(f_n)){
        df <- readRDS(file = f_n);
    } else {
        out <- lapply(seq_len(nr), rarefun)
        df <- do.call(rbind, out);
        saveRDS(df, file = f_n);
    }

    # Get sample data
    if (!is.null(phyloslimR::sample_data(physeq_object, FALSE))) {
      sdf <- methods::as(phyloslimR::sample_data(physeq_object), "data.frame")
      sdf$Sample <- rownames(sdf)
      data <- merge(df, sdf, by = "Sample")
      labels <- data.frame(x = tot, y = S, Sample = rownames(x))
      labels <- merge(labels, sdf, by = "Sample")
    }


    # Add, any custom-supplied plot-mapped variables
    if (length(color) > 1) {
      data$color <- color
      names(data)[names(data) == "color"] <- deparse(substitute(color))
      color <- deparse(substitute(color))
    }

    if (length(label) > 1) {
      labels$label <- label
      names(labels)[names(labels) == "label"] <- deparse(substitute(label))
      label <- deparse(substitute(label))
    }

    if (length(linetype) > 1) {
      data$linetype <- linetype
      names(data)[names(data) == "linetype"] <- deparse(substitute(linetype))
      linetype <- deparse(substitute(linetype))
    }

    p <- ggplot2::ggplot(data = data,
                         ggplot2::aes_string(x = "Size",
                                             y = ".S",
                                             group = "Sample",
                                             color = color,
                                             linetype = linetype))

    p <- p + ggplot2::labs(x = "Sequence Sample Size", y = "Species Richness")

    if (!is.null(label)) {
       p <- p + ggplot2::geom_text(data = labels,
                              ggplot2::aes_string(x = "x",
                                                  y = "y",
                                                  label = label,
                                                  color = color),
                              size = 4, hjust = 0) +
        scale_x_continuous(expand = c(0, 0, 0.3, 0))
    }

    p <- p + ggplot2::geom_line()
    if (se) { ## add standard error if available
      p <- p +
        ggplot2::geom_ribbon(ggplot2::aes_string(ymin = ".S - .se",
                                                 ymax = ".S + .se",
                                                 color = NULL,
                                                 fill = color),
                             alpha = 0.2)
    }
    invisible(p)
}

#'Function to prepare data for phylogenetic tree.
#'@description This functions prepares the data to plot the phylogenetic tree.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PreparePhylogeneticTreePlot <-function(microSetObj, color, shape, taxa, treeshape, imgName, format="png", dpi=72){
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  data <- microSetObj$dataSet$filt.data;
  
  if(class(microSetObj$dataSet$filt.data)=="matrix"){
    data <- otu_table(data,taxa_are_rows =TRUE);
  }
  
  sample_table <- sample_data(microSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data1 <- merge_phyloslim(data,tax_table(microSetObj$dataSet$proc.phyobj), sample_table);
  pg_tree <- readRDS("tree.RDS");
  pg_tb <- tax_table(data1);
  pg_ot <- otu_table(data1);
  pg_sd <- sample_data(data1);
  pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
  data_tree <- merge_phyloslim(pg_tb, pg_ot, pg_sd, pg_tree);

  for (i in sample_variables(data_tree)){
    if(is.integer(sample_data(data_tree)[[i]]) | is.character(sample_data(data_tree)[[i]])){
      sample_data(data_tree)[[i]] <- factor(sample_data(data_tree)[[i]])
    }
  }

  feat_no <- length(phy_tree(tax_glom(data_tree, taxa))$tip.label);

  if(feat_no > 50){
    current.msg <<-"There are too many tree tips to be display, please select a higher taxonomy level";
    return(0);
  }
  
  #adjust height according to number of legends
  w<-600;
  if(feat_no < 10){
    w<-w;
  } else if (feat_no < 20){
    w<-w+100;
  } else if (feat_no < 50){
    w<-w+200;
  } else if (feat_no < 100){
    w<-w+400;
  } else if (feat_no > 100){
    w<-w+500;
  }
  
  imgName = paste(imgName,".", format, sep="");
  box <- plot_tree(tax_glom(data_tree, taxa),
                   ladderize = "left",
                   color = color,
                   shape = shape,
                   label.tips = taxa,
                   size = "abundance");

  if(treeshape == "rectangular"){
    box <- box;
    Cairo(file=imgName, width = 1.5 * w, height = w, type=format, bg="white", dpi=dpi);
    print(box);
    dev.off();
  }else{
    box <- box + coord_polar(theta="y");
    Cairo::Cairo(file=imgName, width = 1.5*w, height = 1.5*w, type=format, bg="white", dpi=dpi);
    print(box);
    dev.off();
  }
  return(.set.microSet(microSetObj))
}

#######################################
###########Tax4Fun/PICRUSt ############
#######################################

#'Main function to perform 16S functional annotation
#'@description This is the main function to perform either PICRUSt or
#'SILVA for functional annotation on the microSetObj.
#'@param microSetObj Input the name of the microSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import Tax4Fun
Perform16FunAnot<-function(microSetObj, type, pipeline) {
  
  microSetObj <- .get.microSetObj(microSetObj);
  
  microSetObj$dataSet$type <- type;
  merge.otu <- readRDS("data.orig");
  merge.otu <- apply(merge.otu, 2, as.numeric);
  rownames(merge.otu) <- microSetObj$dataSet$comp_taxnm;

  #getting whole taxa labels back
  if(type=="SILVA"){
    
    func.meth<-"Tax4Fun";

    if(.on.public.web){
      load_tax4fun();
    }
    
    if(pipeline=="qi_silva"){
            ModSilvaIds <- gsub("uncultured archaeon","",rownames(merge.otu));
            ModSilvaIds <- gsub("uncultured organism","",ModSilvaIds);
            ModSilvaIds <- gsub("uncultured bacterium","",ModSilvaIds);
            ModSilvaIds <- gsub("uncultured crenarchaeote","",ModSilvaIds);
            ModSilvaIds <- gsub("uncultured euryarchaeote","",ModSilvaIds);
            ModSilvaIds <- gsub("; ",";",ModSilvaIds);
            rownames(merge.otu)<-ModSilvaIds;
            merge.otu <- rowsum(as.data.frame(merge.otu),ModSilvaIds);
        }
        data<-list(sampleNames=colnames(merge.otu),otuTable=merge.otu);
        
        if(.on.public.web){
          folderReferenceData <- "../../lib/SILVA123";
        }else{
          folderReferenceData <- "https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/SILVA123";
        }
        Tax4FunOutput <- Tax4Fun(data, folderReferenceData, fctProfiling = TRUE, refProfile = "UProC", shortReadMode = TRUE, normCopyNo = TRUE);
        result2<-as.data.frame(Tax4FunOutput[1]);
        #removing unnecessary data from column name.
        colnames(result2)<-sub("Tax4FunProfile.","",colnames(result2));
        colnames(result2)<-substr(colnames(result2),1,6);
        result<-t(result2);
        result<-as.data.frame(round(1000000*result));  # get to integers
  } else {
    func.meth<-"PICRUSt";
        
    if(type == "Greengenes"){
      a<-rownames(merge.otu);
      rownames(merge.otu)<-NULL;
      merge.otu<-data.frame(merge.otu);
      merge.otu <- cbind(a,merge.otu)
      colnames(merge.otu)[1]<-"X.OTU_IDs";
      
      if(.on.public.web){
        otu.dic <<- readRDS("../../lib/greengenes_taxmap.rds");
      }else{
        otu.dic <<- .read.microbiomeanalyst.lib("greengenes_taxmap.rds");
      }
      
      merge.otu[,1]<-as.character(merge.otu[,1]);
      merge.otu[,1]<-otu.dic[match(merge.otu$X.OTU_IDs,otu.dic$Greengenes),1];
      merge.otu<-merge.otu[!is.na(merge.otu$X.OTU_IDs),];
      nm<-rownames(merge.otu)<-merge.otu[,1];
      merge.otu<-merge.otu[,-1];
      merge.otu<-apply(merge.otu, 2, as.numeric);
      rownames(merge.otu)<-nm;
    }

    query<-as.data.frame(merge.otu);
    samplenm<-colnames(query);

    #normalizing by 16S copy number
    
    if(.on.public.web){
      copyno <- readRDS("../../lib/16S_copyno.rds");
    }else{
      copyno <- .read.microbiomeanalyst.lib("16S_copyno.rds");
    }
  
    result2<-merge(query,copyno, by ="row.names");
    index1<-match(samplenm,colnames(result2), nomatch = NA_integer_, incomparables = NULL);
    result2[index1]<-result2[index1]/result2[['X16S_rRNA_Count']];
    result2[index1]<-round(result2[index1],2);
    rownames(result2)<-result2[,1];
    result2<-result2[,-1];

    # need to fetch and merge results from 5 parts of picrust to get around memory issue (from 2G => 400M)
    res <- data.frame(matrix(0, nrow=nrow(result2), ncol= 6885));
    row.names(res) <- row.names(result2);
    
    if(.on.public.web){
      pc.lib <- readRDS("../../lib/picrust/picrust_part1.rds");
    }else{
      pc.lib <- .read.microbiomeanalyst.lib("picrust_part1.rds", "picrust");
    }
    
    row.hits <- match(row.names(result2), rownames(pc.lib));
    res[, 1:1377] <- pc.lib[row.hits,];
    colnames(res)[1:1377] <- colnames(pc.lib);
    
    if(.on.public.web){
      pc.lib <- readRDS("../../lib/picrust/picrust_part2.rds");
    }else{
      pc.lib <- .read.microbiomeanalyst.lib("picrust_part2.rds", "picrust");
    }
    
    res[, 1378:2754] <- pc.lib[row.hits,];
    colnames(res)[1378:2754] <- colnames(pc.lib);

    if(.on.public.web){
      pc.lib <- readRDS("../../lib/picrust/picrust_part3.rds");
    }else{
      pc.lib <- .read.microbiomeanalyst.lib("picrust_part3.rds", "picrust");
    }
    
    res[, 2755:4131] <- pc.lib[row.hits,];
    colnames(res)[2755:4131] <- colnames(pc.lib);

    if(.on.public.web){
      pc.lib <- readRDS("../../lib/picrust/picrust_part4.rds");
    }else{
      pc.lib <- .read.microbiomeanalyst.lib("picrust_part4.rds", "picrust");
    }
    
    res[, 4132:5508] <- pc.lib[row.hits,];
    colnames(res)[4132:5508] <- colnames(pc.lib);

    if(.on.public.web){
      pc.lib <- readRDS("../../lib/picrust/picrust_part5.rds");
    }else{
      pc.lib <- .read.microbiomeanalyst.lib("picrust_part5.rds", "picrust");
    }
    
    res[, 5509:6885] <- pc.lib[row.hits,];
    colnames(res)[5509:6885] <- colnames(pc.lib);

    ko.nms <- colnames(res);
    res<-merge(result2,res, by ="row.names");

    index2<-match(samplenm,colnames(res), nomatch = NA_integer_, incomparables = NULL);
    index3<-match(ko.nms,colnames(res), nomatch = NA_integer_, incomparables = NULL);

    myList <- vector('list', length(index2));
        
    for (i in 1:length(index2)) {
      myList[[i]]<-data.frame(colSums(res[index3]*res[,index2[i]]));
    }

    res <- NULL;
    gc();

    MyMerge <- function(x, y){
      df<- merge(x, y, by= "row.names", all.x= F, all.y= F);
      rownames(df) <- df$Row.names
      df$Row.names <- NULL
      return(df)
    }
        
    result<- Reduce(MyMerge, myList);
    #orignal class label
    colnames(result)<-samplenm;
  }
    
  result <- round(result, digits =0);
  write.csv(result, file="functionalprof_output.csv");
  result <- result[!apply(result==0,1,all), ]; #filtering zero counts across all

  # save as RDS for memory saving
  microSetObj$analSet$func.pred<-result;
  microSetObj$analSet$func.meth<-func.meth;
    
  if(.on.public.web){
    .set.microSet(microSetObj)
    return(1);
  }else{
    return(.set.microSet(microSetObj))
  }
}
