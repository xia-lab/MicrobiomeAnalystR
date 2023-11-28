##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#'Perform alpha diversity
#'@description This functions performs alpha diversity.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param opt Character, input the name of the statistical method
#'used to calculate the significance of the alpha diversity measure.
#'"tt" for T-test or ANOVA, and "nonpar" for Mann-Whitney or
#'Kruskall-Wallis.
#'@param metadata Character, input the name of the experimental factor
#'to group the samples.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PerformAlphaDiversityComp<-function(mbSetObj, opt, metadata, pair.wise = "false"){
  
  mbSetObj <- .get.mbSetObj(mbSetObj); 
  data <- mbSetObj$analSet$alpha;
  cls <- as.factor(data[,metadata]);
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

    if(pair.wise != "false"){
        # get pairs
          co    <- combn(unique(as.character(cls)),2);
          nco   <- NCOL(co);
          out   <- data.frame(matrix(NA, nrow=nco, ncol=3));
          dimnames(out)[[2]] <- c('pairs', 'statistic', 'pval');

          for(j in 1:nco) {

            inx1 <- cls %in% co[1,j];
            inx2 <- cls %in% co[2,j];

            if(opt=="tt"){
                res <- t.test(x[inx1], x[inx2]);
            }else{
                res <- wilcox.test(x[inx1], x[inx2]);
            }

            out[j,1] <- paste(co[1,j], 'vs', co[2,j]);
            out[j,2] <- res$statistic;
            out[j,3] <- res$p.value;
         }
         
         # add adjusted p value
         out$p.adj <- p.adjust(out$pval, method="fdr");
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

    if(pair.wise != "false"){ # same for two group case
        out <- list();
        out$pairs <- paste(levels(cls)[1], 'vs', levels(cls)[2]);
        out$statistic <- signif(res$statistic, 5);
        out$p.adj <- out$p.value <- signif(res$p.value, 5);
        out<- data.frame(out);
    }
  }
  
  mbSetObj$analSet$alpha.stat.info <- stat.info;
  if(pair.wise != "false"){ 
    rownames(out) <- out$pairs;
    out$pairs <- NULL;
    mbSetObj$analSet$alpha.stat.pair <- mbSetObj$analSet$resTable <- signif(out,5);
    fast.write(mbSetObj$analSet$resTable, file="pairwise_alpha.csv");
  }

  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(stat.info);
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}
####################################
###########Correlation Networks###########
#####################################


#'Function to call for correlation network
#'@description This function runs the fastspar or cor.test 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxrank Character, input the taxonomic level
#' to perform partial correlation analysis.
#'@param cor.method Character, input the correlation method.
#'Supported methods are pearson, spearman, kendall and sparcc.
#'@param colorOpt Character, input what to color the nodes by. Default
#'the nodes will be colored by their expression levels.
#'@param permNum Numeric, input the number of permutations to perform. Default
#'is set to 100.
#'@param pvalCutoff Numeric, input the p-value cutoff. Default is set to 0.05.
#'@param corrCutoff Numeric, input the correlation cutoff. Default is set to 0.3.
#'@param abundOpt Character, default is set to "mean".
#'@param corr.net.name Character, input the name of the plot to save.
#'@param plotNet Boolean. Set to TRUE if you would like a network visualization
#'outputted to your current working directory. 
#'@param networkType Character, "static" to create a static image or
#'"interactive" to create an interactive network saved as an html 
#'in your working directory.
#'@param netLayout Character, layout from ggraph. "kk" for the spring-based algorithm by Kamada and Kawai
#'as default. "drl" for force directed algorithm from the DrL toolbox. "lgl" for Large Graph Layout. "fr" for
#'force-directed of Fruchterman and Reingold.
#'@param netTextSize Numeric, input the preferred font size to be used in the network
#'plot.
#'@import ppcor
#'@import igraph
#'@export

PerformNetworkCorrelation <- function(mbSetObj, taxrank, cor.method="pearson", colorOpt="expr", 
                                      permNum=100, pvalCutoff=0.05, corrCutoff=0.3, abundOpt="mean", 
                                      corr.net.name, plotNet = FALSE, netType="static", netLayout="kk",
                                      netTextSize = 2.5){
  if(!exists("my.corr.net")){ # public web on same user dir
    .load.scripts.on.demand("utils_corrnet.Rc");    
  }
  return(my.corr.net(mbSetObj, taxrank, cor.method, colorOpt, permNum, pvalCutoff, corrCutoff, abundOpt, corr.net.name, plotNet, netType, netLayout,netTextSize));
}

#'Function to create box plots of important features
#'@description This functions plots box plots of a selected feature.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param boxplotName Character, input the name of the 
#'box plot.
#'@param feat Character, input the name of the selected 
#'feature.
#'@param format Character, by default the plot format
#'is "png".
#'@param dpi Dots per inch. Numeric, by default
#'it is set to 72.
#'@parm colorPal Character, input the name of the preferred color palette.
#'Use "default" for the RColor brewer Set1 palette, "virdis" for the viridis color palette, and
#'"dark" for the RColor brewer Dark2 palette.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import grid
#'@import gridExtra
PlotBoxDataCorr<-function(mbSetObj, boxplotName, feat, format="png", sel.meta="", dpi=72){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  colorPal <- "dark";
  
  load_ggplot();
  load_grid();
  load_gridExtra();
  if(mbSetObj$module.type == "meta"){
    merged <- qs::qread("merged.data.qs");
    merged <- subsetPhyloseqByDataset(mbSetObj, merged);
    data <- t(as.data.frame(as(otu_table(merged), "matrix"),check.names=FALSE));
    data <- as.data.frame(data);
    variable <- sel.meta;
    data$class <- as.vector(as.matrix(sample_data(merged))[,sel.meta])
    
    a <- as.numeric(data[,feat]);
  }else{
    variable <- mbSetObj$analSet$var.typecor;
    data <- mbSetObj$analSet$boxdatacor;
    a <- as.numeric(data[,feat]);
  }
  min.val <- min(abs(a[a!=0]))/5;
  data$log_feat <- log2((a + sqrt(a^2 + min.val))/2);
  boxplotName = paste(boxplotName,".",format, sep="");
  
  numGrps <- length(levels(data$class))
  
  if(numGrps == 2){
    width <- 325
  }else if(numGrps < 4){
    width <- 350
  }else if(numGrps < 6){
    width <- 375
  }else{
    width <- 400
  }
  
  Cairo::Cairo(file=boxplotName, width=width, height=300, type=format, bg="white", dpi=dpi);
  
  box <- ggplot(data, aes(x=data$class, y = data$log_feat, fill=as.factor(class))) + stat_boxplot(geom ='errorbar') + 
    geom_boxplot(outlier.shape = NA) + geom_jitter() + theme_bw() + labs(y="Log-transformed Counts\n", x=paste0("\n",variable), fill=variable) +
    ggtitle(feat) + theme(plot.title = element_text(hjust=0.5, size=13, face="bold"), axis.title=element_text(size=11), legend.title=element_text(size=11), axis.text=element_text(size=10));
  #remove grid
  box <- box + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank(), panel.border = element_rect(colour = "#787878", fill=NA, size=0.5))
  
  if(colorPal == "viridis"){
    box <- box + scale_fill_viridis_d()
  }else if(colorPal == "set1"){
    box <- box + scale_fill_brewer(palette="Set1")
  }else if(colorPal == "dark"){
    box <- box + scale_fill_brewer(palette="Dark2")
  }
  
  print(box)
  dev.off();
  return(.set.mbSetObj(mbSetObj))
}

SparccToNet <- function(mbSetObj=NULL, corr.net.name, networkType="static", netLayout="kk", netTextSize){
  if(!exists("my.sparcc.net")){ # public web on same user dir
    .load.scripts.on.demand("utils_sparcc.Rc");    
  }
  return(my.sparcc.net(mbSetObj, corr.net.name, networkType, netLayout, netTextSize));
}

simpleCap <- function(x) {
  s <- strsplit(x, " ")[[1]]
  paste(toupper(substring(s, 1,1)), substring(s, 2),
        sep="", collapse=" ")
}

PerformLayOut <- function(g){
  vc <- vcount(g);
  if(vc > 3000) {
    pos.xy <- layout.lgl(g, maxiter = 100);
  }else if(vc > 2000) {
    pos.xy <- layout.lgl(g, maxiter = 150);
  }else if(vc > 1000) {
    pos.xy <- layout.lgl(g, maxiter = 200);
  }else if(vc < 150){
    pos.xy <- layout.kamada.kawai(g);
  }else{
    pos.xy <- layout.fruchterman.reingold(g);
  }
  pos.xy;
}

rescale2NewRange <- function(qvec, a, b){
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}

PrepareCorrExpValues <- function(mbSetObj, meta, taxalvl, color, layoutOpt, comparison, wilcox.cutoff){
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_metacoder();

  tax_o <- taxalvl;
  dm <- mbSetObj$dataSet$proc.phyobj;
  dims <- ncol(tax_table(dm))
  tax_table_new = data.frame("Kingdom" = "Root", as(tax_table(dm), "matrix")[, 1:dims],check.names=FALSE) # add root to tax table
  tax_table(dm) <- as.matrix(tax_table_new);
  
  dm_samples = as(sample_data(dm), "data.frame");
  dm_samples <- cbind.data.frame("sample_id" = row.names(dm_samples), dm_samples); # add sample_id to sample table
  row.names(dm_samples) <- c();
  dm_samples[] <- lapply(dm_samples, as.character);
  
  if(comparison != "all"){
    grp.nms <- mbSetObj$corr.net$comparison;
    if(is.null(grp.nms)){
      AddErrMsg("Please specify groups first!");
      return(0);
    }
    dm_samples_cmf <- dm_samples[dm_samples[[meta]] %in% grp.nms, ]; #subset sample data by meta variable
  }else{
    grp.nms <- unique(dm_samples[[meta]]) # all 
    dm_samples_cmf <- dm_samples
  }
  
  otu_dm <- as.data.frame(as(otu_table(dm), "matrix"),check.names=FALSE);
  tax_dm <- as.data.frame(as(tax_table(dm), "matrix"),check.names=FALSE);
  tax_dm[] <- lapply(tax_dm, as.character);#make sure characters in tax_dm;
  
  depth <- ncol(tax_dm)
  rank_dm <- c("r", "p", "c", "o", "f", "g", "s");
  names(tax_dm) <- rank_dm[1:depth];
  
  for (i in 1:ncol(tax_dm)){
    for (j in 1:nrow(tax_dm)){
      if (is.na(tax_dm[j, i])){
        tax_dm[j, i] <- "";
      } else {
        tax_dm[j, i] <- paste(names(tax_dm)[i],
                              tax_dm[j, i],
                              sep = "__");
      }
    }
  } #add __ to tax table
  
  if(taxalvl == "Phylum"){
    tax <- "p";
  } else if(taxalvl == "Class"){
    tax <- "c";
  } else if(taxalvl == "Order"){
    tax <- "o";
  } else if(taxalvl == "Family"){
    tax <- "f";
  } else if(taxalvl == "Genus"){
    tax <- "g"
  } else {
    tax <- "s";
  }; # get tax rank for heat tree
  
  tax_dm <- tax_dm[, 1:which(names(tax_dm) == tax)]; #subset tax table
  rank_dm_new <- rank_dm[1:which(rank_dm == tax)];
  tax_dm$lineage <- apply(tax_dm[, rank_dm_new], 1, paste, collapse = ";"); #collapse all tax ranks
  dm_otu <- cbind.data.frame("otu_id" = row.names(tax_dm),
                             "lineage" = tax_dm$lineage,
                             otu_dm); #make new otu table
  row.names(dm_otu) <- c();
  dm_otu$lineage <- gsub(";+$", "", dm_otu$lineage); #remove empty tax names
  
  dm_otu_cmf <- dm_otu[, c("otu_id", "lineage", dm_samples_cmf$sample_id)]; #make otu table ready for heat tree  
  
  PrepareHeatTreePlotDataParse_cmf_res <- PrepareHeatTreePlotDataParse_cmf(dm_otu_cmf, dm_samples_cmf, meta);
  dm_obj_cmf = PrepareHeatTreePlotDataParse_cmf_res;
  mbSetObj$dataSet$selected.grps = dm_samples_cmf$sample_id
  mbSetObj$dataSet$comparison = grp.nms 
  mbSetObj$dataSet$meta = meta
  cls = data.frame(dm_obj_cmf$data$class_data,check.names=FALSE)
  cls = cls[!duplicated(cls$taxon_id),]
  cls = cls[,c(1,4)]
  diff_table = data.frame(dm_obj_cmf$data$diff_table,check.names=FALSE)
  df = merge(diff_table, cls, by="taxon_id")
  df = df[which(df$tax_name != "Root"),]
  df = df[order(df$tax_name),]
  mbSetObj$analSet$diff_table <- df;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj);
    PrepareBoxPlot(mbSetObj, taxalvl, mbSetObj$dataSet$meta);
    return(1)
  }else{
    mbSetObj <- PrepareBoxPlot(mbSetObj, taxalvl, mbSetObj$dataSet$meta);
    return(.set.mbSetObj(mbSetObj));
  }
}

PrepareBoxPlot <- function(mbSetObj, taxrank, variable){

  load_phyloseq();
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  selSamples <- mbSetObj$dataSet$selected.grps
  
  grps = mbSetObj$dataSet$comparison
  taxrank_boxplot <- taxrank;
  claslbl_boxplot <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  claslbl_boxplot <- claslbl_boxplot[which(claslbl_boxplot %in% grps)];
  # build phyloslim obj in fly
  filt.dataphy <- mbSetObj$dataSet$filt.data;
  filt.dataphy <- apply(filt.dataphy, 2, as.integer);
  filt.dataphy <- otu_table(filt.dataphy, taxa_are_rows =TRUE);
  sample_table_boxplot <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  filt.dataphy <- merge_phyloseq(filt.dataphy, sample_table_boxplot);
  taxa_names(filt.dataphy) <- rownames(mbSetObj$dataSet$filt.data);
  data_boxplot <- filt.dataphy;
  
  if(mbSetObj$module.type=="mdp"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data_boxplot <- merge_phyloseq(data_boxplot, mbSetObj$dataSet$taxa_table);
  }
  
  #using by default names for shotgun data
  if(mbSetObj$module.type=="sdp"){
    taxrank_boxplot<-"OTU";
  }
  
  if(taxrank_boxplot!="OTU"){
    #merging at taxonomy levels
    data_boxplot <- fast_tax_glom_mem(data_boxplot, taxrank_boxplot);
    if(is.null(data_boxplot)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
    
    nm_boxplot <- tax_table(data_boxplot)[,taxrank_boxplot];
    nm_boxplot <- as.character(nm_boxplot);
    #converting NA values to unassigned
    nm_boxplot[is.na(nm_boxplot)] <- "Not_Assigned";
    data1_boxplot <- as.matrix(otu_table(data_boxplot));
    rownames(data1_boxplot) <- nm_boxplot;
    
    #all NA club together
    data1_boxplot <- as.matrix(t(sapply(by(data1_boxplot, rownames(data1_boxplot), colSums), identity)));
    data1_boxplot <- otu_table(data1_boxplot,taxa_are_rows=T);
    data_boxplot <- merge_phyloseq(data1_boxplot, sample_data(data_boxplot));
  }
  
  nm_boxplot <- taxa_names(data_boxplot);
  dat3t_boxplot <- as.data.frame(t(otu_table(data_boxplot)),check.names=FALSE);
  colnames(dat3t_boxplot) <- nm_boxplot; 
  
  #individual boxplot for features
  box_data <- as.data.frame(dat3t_boxplot[which(rownames(dat3t_boxplot) %in% selSamples),],check.names=FALSE);
  box_data$class <- claslbl_boxplot;
  mbSetObj$analSet$boxdatacor <- box_data;
  mbSetObj$analSet$var.typecor <- variable;
  
  return(.set.mbSetObj(mbSetObj));
}

#####################################
###########Core Microbiome###########
#####################################

#'Perform core microbiome analysis
#'@description This functions performs core microbiome analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param imgName Input the name of the core microbiome analysis plot.
#'@param preval Numeric, input the sample prevalence cutoff. This is
#'the minimum number of samples that share this microbe.
#'@param detection Numeric, input the relative abundance cutoff. This
#'is the minimum detection threshold of X percent relative abundance.
#'@param palette Character, input the color palette.
#'"bwm" for the default color palette /(blue, white, red/), "gbr"
#'for the red and green color palette, "heat" for the yellow to red color
#'palette, "topo" for the blue, green, and yellow color palette, 
#'"gray" for the gray color palette, "byr" for the blue, yellow and red color
#'palette, "viridis" for the viridis default color palette and "plasma"
#'for the plasma viridis color palette.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@param interactive Boolean, if set to TRUE, saves the plot
#'as an interactive html plot.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import RColorBrewer
#'@import viridis
CoreMicrobeAnalysis<-function(mbSetObj, imgName, preval, detection, taxrank,
                              palette, viewOpt, analOpt, expFact, group, 
                              format="png", dpi=72, width=NA, interactive = FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data <- mbSetObj$dataSet$proc.phyobj;
  #print(expFact)
  expFact <- expFact
  group <- group
  
  if(!analOpt == "all_samples"){
    data <- eval(parse(text = paste("phyloseq:::subset_samples(data,", expFact, "==", "\"", group, "\"", ")", sep="")))
    
    # check min 2 reps 
    samples_left <- nsamples(data)
    
    if(samples_left<2){
      AddErrMsg("More than 2 replicates are required in your group!")
      return(0)
    }
  }
  
  if(taxrank=="OTU"){
    data <- otu_table(data,taxa_are_rows=T);
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_mem(data, taxrank);
    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
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
  
  #transform data to relative abundances, then obtain full phyloseq obj of just core microbiota
  data.compositional <- transform_sample_counts(data,function(x) x / sum(x));
  data.core <- core(data.compositional, detection = detection, prevalence = preval);
  core.nm <- data.frame(prevalence(data.compositional, detection = detection, sort = TRUE),check.names=FALSE);
  colnames(core.nm)[1] <- "Prevelance";
  fileName <- "core_microbiome.csv";
  fast.write(core.nm, file=fileName);
  
  imgName = paste(imgName, ".", format, sep="");
  mbSetObj$imgSet$core <- imgName;
  
  #if more than 1500 features will be present; subset to most abundant=>1500 features.
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
  h <- round(myH/65);
  
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
  }else if(palette == "viridis") {
    colors <- rev(viridis::viridis(10))
  }else if(palette == "plasma") {
    colors <- rev(viridis::plasma(10))
  }else {
    load_rcolorbrewer();
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(10));
  }
  
  p <- plot_core(data.core, plot.type = "heatmap", colours = colors, prevalences = seq(.05, 1, .05), 
                 detections = 10^seq(log10(detection), log10(max(abundances(data.core))), length = 10)) + 
    ylab(paste0("\n", taxrank)) + xlab("\nDetection Threshold (Relative Abundance (%))") + 
    guides(fill = guide_legend(keywidth = 1.5, keyheight = 1)) + 
    theme(axis.text=element_text(size=10), axis.title=element_text(size=11.5), legend.title=element_text(size=11));
  
  print(p);
  dev.off();
  
  if(interactive){
    library(plotly)
    library(htmlwidgets)
    ax <- list(
      zeroline = FALSE,
      showline = FALSE,
      showgrid = FALSE
    )
    p <- plotly::ggplotly(p)
    p <- p %>% layout(xaxis = ax, yaxis = ax)
    htmlwidgets::saveWidget(p, "core_interactive.html")
  }
  mbSetObj$analSet$core<-signif(as.matrix(core.nm),5);
  mbSetObj$analSet$core.taxalvl<-taxrank;
  mbSetObj$paramSet$core <- list(
        taxalvl=taxrank,
        preval=preval, 
        detection=detection,
        analOpt = analOpt,
        expFact = expFact, 
        group = group
    );
  return(.set.mbSetObj(mbSetObj));
}

######################Adapted from the Microbiome R package#########################

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
  if (class(x) == "phyloseq") {
    # Pick OTU matrix
    otu <- get_taxa(x)
    # Ensure that taxa are on the rows
    if (all(c(!taxa_are_rows(x), ntaxa(x) > 1, nsamples(x) > 1))) {
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
  
  load_datatable();
  
  data <- x;
  
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
  
  df <- as.data.frame(prev,check.names=FALSE)
  df$ID <- rownames(prev)
  df <- data.table::melt(df, "ID");
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
  p <- p + scale_x_log10(breaks=round(as.numeric(names(split(df$DetectionThreshold, sort(df$DetectionThreshold%%3)))),3))
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
  if (class(x) == "phyloseq") {
    x <- abundances(x)
  }
  if (is.vector(x)) {
    if (include.lowest) {
      prev <- sum(x >= detection)
    } else {
      prev <- sum(x > detection)
    }
  } else if (any(c(is.matrix(x), is.data.frame(x)))) {
    
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
  } else if (any(c(is.matrix(x), is.data.frame(x)))) {
    n <- ncol(x)
  }
  n
}

###########################
### Pie Chart Functions ###
###########################

#'Main function to plot pie graphs of microbiome data.
#'@description This functions plots pie graphs of microbiome data. 
#'In particular, it plots the overall pie chart (all samples). The
#'data used to calculate the pie-chart is saved in your working 
#'directory as ""piechart_abundances.csv".
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param feat_cnt Set the minimum feature count that is used to "bin"
#'together small taxa. 
#'@param calcmeth Merge small taxa based on the sum /("sum"/)
#'or median /("med"/) counts across all samples or groups.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotOverallPieGraph<-function(mbSetObj, taxalvl, feat_cnt, calcmeth, 
                              toptaxapie, pietoptaxaopt){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_reshape();
  load_phyloseq();

  set.seed(28053448);
  
  #using filtered data
  data <- mbSetObj$dataSet$filt.data;
  
  if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
    data <- otu_table(data, taxa_are_rows =TRUE);
  }
  
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  datapie <<- merge_phyloseq(data, tax_table(mbSetObj$dataSet$proc.phyobj), sample_table);
  
  #using reduce names
  data <- otu_table(datapie);
  data <- data.frame(data,check.names=FALSE);
  data_tax <- tax_table(datapie);
  
  #reshaping data
  data <- t(data);
  
  if(taxalvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    taxa_nm <- data.matrix(data_tax[,taxalvl]);
  }
  
  y <- which(is.na(taxa_nm)==TRUE);
  
  #converting NA values to unassigned
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data.frame(data %*% sapply(unique(nms),"==",nms),check.names=FALSE);
  colnames(data)<- gsub("\\."," ",colnames(data));
  
  if(pietoptaxaopt == "bottom"){
    if(calcmeth=="sum"){
      ind <- which(colSums(data)>feat_cnt);
      ind1 <- which(colSums(data)<feat_cnt);
    }else{
      dt <- apply(data,2,median);
      ind <- which(dt>feat_cnt);
      ind1 <- which(dt<feat_cnt);
    }
    if(length(ind)==0){
      AddErrMsg("All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.");
      return(0);
    }
    
    if(length(ind1)>0){
      colnames(data)[ind1] <- "Others";
      data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)),check.names=FALSE);
    }
    
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step);
    piedata <- aggregate(. ~variable , data=data[-1], FUN=sum);
    
    # order by abundance
    ord.inx <- order(piedata$value, decreasing = TRUE);
    piedata <- piedata[ord.inx,];
    piedata_write <- piedata;
  } else {
    
    if(calcmeth=="sum"){
      order_taxa <- order(sapply(data, sum), decreasing = TRUE)
    } else {
      order_taxa <- order(sapply(data, median), decreasing = TRUE)
    }
    
    data <- data[, order_taxa];
    data <- data.frame("variable" = names(data),
                       "value" = apply(data, 2, sum),check.names=FALSE);
    row.names(data) <- NULL;
    data$variable <- as.character(data$variable);
    
    piedata_write <- data[order(data$value, decreasing = TRUE), ]
    piedata <- data;
    contain_notassign <- any(grepl("^Not_Assigned$", piedata$variable));
    if(nrow(piedata) < toptaxapie){
      piedata <- piedata;
    } else {
      if (! contain_notassign){
        top_piedata <- head(piedata, n = toptaxapie);
        bottom_piedata <- tail(piedata, n = nrow(piedata) - toptaxapie);
        other_piedata <- data.frame("variable" = "Others", 
                                    "value" = sum(bottom_piedata$value),check.names=FALSE);
        piedata <- rbind.data.frame(top_piedata, other_piedata);
      } else {
        if(nrow(piedata) == toptaxapie + 1){
          piedata <- piedata;
        } else {
          not_assigned_piedata <- piedata[which(piedata$variable == "Not_Assigned"), ];
          piedata <- piedata[which(piedata$variable != "Not_Assigned"), ];
          top_piedata <- head(piedata, n = toptaxapie);
          bottom_piedata <- tail(piedata, n = nrow(piedata) - toptaxapie);
          other_piedata <- data.frame("variable" = "Others", 
                                      "value" = sum(bottom_piedata$value),check.names=FALSE);
          piedata <- rbind.data.frame(top_piedata, other_piedata, not_assigned_piedata);
        }
      }
    }
    piedata <- piedata[order(piedata$value, decreasing = TRUE), ];
  }
  
  piedata_write$percentage <- round((piedata_write$value / sum(piedata_write$value) * 100), digits = 2)
  colnames(piedata_write) <- c("Taxa", "Abundance", "Percentage")
  fast.write(piedata_write, "piechart_abundances.csv");
  
  piedata <<- piedata;
  
  return(.set.mbSetObj(mbSetObj));
}

#'Main function to plot pie graphs of microbiome data.
#'@description This functions plots pie graphs of microbiome data. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param metadata Character, select which grouping to use.
#'@param clslevel Character, input which group to plot.
#'@param feat_cnt Set the minimum feature count that is used to "bin"
#'together small taxa. 
#'@param calcmeth Merge small taxa based on the sum /("sum"/)
#'or median /("med"/) counts across all samples or groups.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotGroupPieGraph <- function(mbSetObj, taxalvl, metadata, clslevel,
                              feat_cnt, calcmeth, toptaxapie, pietoptaxaopt){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_reshape();
  load_phyloseq();
  
  set.seed(28053443);
  
  #using filtered data
  data <- mbSetObj$dataSet$filt.data;
  
  if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
    data<-otu_table(data,taxa_are_rows =TRUE);
  }
  
  smpl <- data.frame(sample_data(mbSetObj$dataSet$proc.phyobj),check.names=FALSE);
  smpl1 <- sample_data(subset(smpl,smpl[metadata]==clslevel));
  data <- prune_samples(sample_names(smpl1), data);
  data <- merge_phyloseq(data, smpl1);
  mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
  datapie <<- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
  
  #using reduce names
  data <- t(data.frame(otu_table(datapie),check.names=FALSE));
  data_tax <- tax_table(datapie);
  
  if(taxalvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    taxa_nm <- data.matrix(data_tax[,taxalvl]);
  }
  
  #converting NA values to unassigned
  y <- which(is.na(taxa_nm)==TRUE);
  taxa_nm[y] <- "Not_Assigned";
  nms <- colnames(data) <- taxa_nm[,1];
  data <- data.frame(data %*% sapply(unique(nms),"==",nms),check.names=FALSE);
  colnames(data) <- gsub("\\."," ",colnames(data));
  
  if(pietoptaxaopt == "bottom"){
    if(calcmeth=="sum"){
      ind <- which(colSums(data)>feat_cnt);
      ind1 <- which(colSums(data)<feat_cnt);
    }else{
      dt <- apply(data,2,median);
      ind <- which(dt>feat_cnt);
      ind1 <- which(dt<feat_cnt);
    }
    if(length(ind)==0){
      AddErrMsg("All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.");
      return(0);
    }
    
    if(length(ind1)>0){
      colnames(data)[ind1] <- "Others";
      data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)),check.names=FALSE);
    }
    
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step);
    piedata <- aggregate(. ~variable , data=data[-1], FUN=sum);
    
    # order by abundance
    ord.inx <- order(piedata$value, decreasing = TRUE);
    piedata <- piedata[ord.inx,];
    piedata_write <- piedata;
  } else {
    
    if(calcmeth=="sum"){
      order_taxa <- order(sapply(data, sum), decreasing = TRUE)
    } else {
      order_taxa <- order(sapply(data, median), decreasing = TRUE)
    }
    
    data <- data[, order_taxa];
    data <- data.frame("variable" = names(data),
                       "value" = apply(data, 2, sum),check.names=FALSE);
    row.names(data) <- NULL;
    data$variable <- as.character(data$variable);
    
    piedata_write <- data[order(data$value, decreasing = TRUE), ]
    piedata <- data;
    contain_notassign <- any(grepl("^Not_Assigned$", piedata$variable));
    if(nrow(piedata) < toptaxapie){
      piedata <- piedata;
    } else {
      if (! contain_notassign){
        top_piedata <- head(piedata, n = toptaxapie);
        bottom_piedata <- tail(piedata, n = nrow(piedata) - toptaxapie);
        other_piedata <- data.frame("variable" = "Others", 
                                    "value" = sum(bottom_piedata$value),check.names=FALSE);
        piedata <- rbind.data.frame(top_piedata, other_piedata);
      } else {
        if(nrow(piedata) == toptaxapie + 1){
          piedata <- piedata;
        } else {
          not_assigned_piedata <- piedata[which(piedata$variable == "Not_Assigned"), ];
          piedata <- piedata[which(piedata$variable != "Not_Assigned"), ];
          top_piedata <- head(piedata, n = toptaxapie);
          bottom_piedata <- tail(piedata, n = nrow(piedata) - toptaxapie);
          other_piedata <- data.frame("variable" = "Others", 
                                      "value" = sum(bottom_piedata$value),check.names=FALSE);
          piedata <- rbind.data.frame(top_piedata, other_piedata, not_assigned_piedata);
        }
      }
    }
    piedata <- piedata[order(piedata$value, decreasing = TRUE), ];
  }
  piedata_write$percentage <- round((piedata_write$value / sum(piedata_write$value) * 100), digits = 2)
  colnames(piedata_write) <- c("Taxa", "Abundance", "Percentage")
  fast.write(piedata_write, "piechart_abundances.csv");
  piedata <<- piedata;
  return(.set.mbSetObj(mbSetObj));
}

#'Main function to plot sample-wise pie graphs of microbiome data.
#'@description This functions plots sample-wise pie graphs of microbiome data. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param feat_cnt Set the minimum feature count that is used to "bin"
#'together small taxa. 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@exportmbSetObj, taxalvl, smplmbSetObj, taxalvl, smplnm, feat_cnt, calcmethnm, feat_cnt, calcmeth
#'@import reshape
PlotSamplePieGraph<-function(mbSetObj, taxalvl, smplnm, feat_cnt, toptaxapie, pietoptaxaopt){
  mbSetObj <- .get.mbSetObj(mbSetObj);

  load_reshape();
  load_phyloseq();
  set.seed(28053443);

  #using filtered data
  data <- mbSetObj$dataSet$filt.data;
  
  if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
    data <- otu_table(data,taxa_are_rows =TRUE);
  }
  
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data <- merge_phyloseq(data, tax_table(mbSetObj$dataSet$proc.phyobj), sample_table);
  data <- prune_samples(smplnm, data);
  mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
  #datapie <<- merge_phyloseq(data, mbSetObj$dataSet$taxa_table)
  pie_try <- try(merge_phyloseq(data, mbSetObj$dataSet$taxa_table));
  if(class(pie_try) == "try-error"){
    datapie <<- data;
  } else {
    datapie <<- merge_phyloseq(data, mbSetObj$dataSet$taxa_table);
  }
  
  #using reduce names
  data <- otu_table(datapie);
  data <- data.frame(data,check.names=FALSE);
  zero_row <- which(data[[1]] != 0);
  data_tmp <- as.data.frame(data[zero_row, ],check.names=FALSE)
  row.names(data_tmp) <- row.names(data)[zero_row];
  names(data_tmp) <- smplnm;
  data <- data_tmp;
  
  data_tax <- tax_table(datapie);
  data_tax <- data_tax[zero_row, ]
  
  #reshaping data
  data <- t(data);
  
  if(taxalvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    taxa_nm <- data.matrix(data_tax[,taxalvl]);
  }
  
  y <- which(is.na(taxa_nm)==TRUE);
  #converting NA values to unassigned
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data.frame(data %*% sapply(unique(nms),"==",nms),check.names=FALSE);
  colnames(data) <- gsub("\\."," ",colnames(data));
  
  if(pietoptaxaopt == "bottom"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
    
    if(length(ind)==0){
      AddErr("All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.");
      return(0);
    }
    
    if(length(ind1)>0){
      colnames(data)[ind1] <- "Others";
      data<-as.data.frame(t(rowsum(t(data),group = colnames(data))),check.names=FALSE);
    }
    
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step);
    piedata <- aggregate(. ~variable , data=data[-1], FUN=sum);
    
    # order by abundance
    ord.inx <- order(piedata$value, decreasing = TRUE);
    piedata <- piedata[ord.inx,];
  } else {
    piedata <- data.frame("variable" = colnames(data),
                          "value" = t(data)[, 1],check.names=FALSE)
    # order by abundance
    ord.inx <- order(piedata$value, decreasing = TRUE);
    piedata <- piedata[ord.inx,];
    contain_notassign <- any(grepl("^Not_Assigned$", piedata$variable));
    if(nrow(piedata) <= toptaxapie){
      piedata <- piedata;
    } else {
      if (! contain_notassign){
        top_piedata <- head(piedata, n = toptaxapie);
        bottom_piedata <- tail(piedata, n = nrow(piedata) - toptaxapie);
        other_piedata <- data.frame("variable" = "Others", 
                                    "value" = sum(bottom_piedata$value),check.names=FALSE);
        piedata <- rbind.data.frame(top_piedata, other_piedata);
      } else {
        if(nrow(piedata) == toptaxapie + 1){
          piedata <- piedata;
        } else {
          not_assigned_piedata <- piedata[which(piedata$variable == "Not_Assigned"), ];
          piedata <- piedata[which(piedata$variable != "Not_Assigned"), ];
          top_piedata <- head(piedata, n = toptaxapie);
          bottom_piedata <- tail(piedata, n = nrow(piedata) - toptaxapie);
          other_piedata <- data.frame("variable" = "Others", 
                                      "value" = sum(bottom_piedata$value),check.names=FALSE);
          piedata <- rbind.data.frame(top_piedata, other_piedata, not_assigned_piedata);
        }
      }
    }
  }
  piedata <<- piedata;
  return(.set.mbSetObj(mbSetObj));
}

#'Function to plot pie-chart data.
#'@description This functions plots pie charts of microbiome data. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape

PlotDataPieFromPie<-function(mbSetObj, taxalvl, metadata, clslevel,
                             taxaposn, lowtaxa){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_reshape();
  set.seed(280534431);
  high_taxa <<- as.character(piedata$variable[taxaposn]);
  
  if(high_taxa=="Not_Assigned"){
    AddErrMsg("Cannot map to Not_Assigned level!");
    return(0);
  }
  
  lowlvl_nm <<- lowtaxa;
  
  datataxa <- as.matrix(tax_table(datapie));
  subsettax_table <- tax_table(subset(datataxa,datataxa[,taxalvl]==high_taxa));
  data1 <- prune_taxa(taxa_names(subsettax_table),datapie);
  datapietaxatab <<- data_tax <- tax_table(data1);
  datapietaxa <<- data <- t(data.frame(otu_table(data1),check.names=FALSE));
  taxa_nm <- data.matrix(data_tax[,lowlvl_nm]);
  
  #converting NA values to unassigned
  y <- which(is.na(taxa_nm)==TRUE);
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data.frame(data %*% sapply(unique(nms),"==",nms),check.names=FALSE);
  
  if(length(nms)==1){
    colnames(data)<-nms;
  }
  
  colnames(data) <- gsub("\\."," ",colnames(data));
  data$step <- factor(rownames(data));
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  
  fact <- factor(data$variable)
  levels(fact) <- sub("^X", "", levels(fact))
  
  color_var <- levels(fact);
  x <- length(color_var);
  x.colors <<- rep(col_vector,length.out=x);
  
  piedata2 <- aggregate(. ~variable , data=data[-1], FUN=sum);
  # order by abundance
  ord.inx <- order(piedata2$value, decreasing = TRUE);
  piedata2 <<- piedata2[ord.inx,];
  
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to update pie-chart data.
#'@description This functions updates pie chart data. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
UpdatePieData<-function(mbSetObj, lowtaxa){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_reshape();
  
  set.seed(280534432);
  
  data <- datapietaxa;
  data_tax <- datapietaxatab;
  high_taxa <- high_taxa;
  lowlvl_nm <<- lowtaxa;
  taxa_nm <- as.data.frame(data_tax[,lowlvl_nm],check.names=FALSE);
  taxa_nm <- as.matrix(taxa_nm);
  y <- which(is.na(taxa_nm)==TRUE);
  
  #converting NA values to unassigned
  taxa_nm[y] <- "Not_Assigned";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data,check.names=FALSE);
  
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
  
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to save pie-chart
#'@description This functions saves created pie chart plot. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param pieName Character, input the name of the pie chart plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SavePiechartImg <- function(mbSetObj, taxalvl, pieName, format="png", dpi=72) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  set.seed(280);
  
  pieName = paste(pieName,".", format, sep="");
  piedata <- transform(transform(piedata, value=value/sum(value)));
  
  #rownames are still arranged by decending order
  piedataimg <- piedata;
  mbSetObj$imgSet$pie <- pieName;
  
  row.names(piedataimg) <- NULL;
  x.cols <- pie.cols;
  
  # java color code to R color code
  x.cols <- paste("#",x.cols, sep="");
  Cairo::Cairo(file=pieName, width=630, height=500, type=format, bg="white", dpi=dpi);
  
  box <- ggplot(piedataimg, aes(x="", y = value, fill=reorder(variable,-value))) +
    geom_bar(width = 1, stat = "identity") + theme_bw() +
    coord_polar(theta = "y",direction=-1,start = 4.71239) + scale_fill_manual(values=c(x.cols))+
    geom_text(aes(x=1.6,label = scales::percent(round(value,2))), check_overlap = T,size=3,position = position_stack(vjust = 0.5),color="grey48") +
    theme(legend.position="right",axis.text = element_blank(),axis.ticks = element_blank(),panel.grid  = element_blank(), plot.title = element_text(hjust=0.5, face="bold"),legend.text=element_text(color="grey48")) +
    labs(x="", y="",fill="");
  
  print(box);
  mbSetObj$analSet$pie<-piedataimg;
  mbSetObj$analSet$pie.taxalvl<-taxalvl;
  dev.off();
  
  return(.set.mbSetObj(mbSetObj))
}

#'Function to create pie-chart plot.
#'@description This functions creates the pie chart plot. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotPiechart <- function(mbSetObj, rel_perct, pieName, format="png", dpi=72) {
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  set.seed(28056188); 
  pieName = paste(pieName,".", format, sep="");
  piedata2 <- transform(transform(piedata2, value=value/sum(value)));
  ind <- which(piedata2[,"value"]>rel_perct);
  ind1 <- which(piedata2[,"value"]<rel_perct);
  
  if(length(ind)==0){
    AddErrMsg("All features have lower relative abundance than given minimum abundance. Please lower the cut off for relative abundance.");
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
  
  return(.set.mbSetObj(mbSetObj));
}

#######################################
###########Alpha-diversity#############
#######################################

#'Function to plot alpha-diversity analysis.
#'@description This functions creates a plot for alpha-diversity. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param data.src Character, input whether alpha diversity 
#'is calculated using the filtered /("filt"/) or raw data /("orig"/).
#'@param bargraphName Character, input the name of the plot.
#'@param distName Character, input the diversity measure
#'to calculate alpha-diversity. "Chao1", "Observed", "ACE", "Shannon", 
#'"Simpson", or "Fisher".
#'@param metadata Character, input the name of the experimental
#'factor to group the samples.
#'@param taxrank Character, input the taxonomic level to calculate
#'alpha-diversity. "Phylum", "Class", "Order",
#'"Family", "Genus", "Species" or "OTU".   
#'@param group Boolean, input whether or not to group the samples in the 
#'dot plot. 0 or 1.
#'@param colors Character, set to "default" to use the default colors, 
#'"viridis" to use the viridis color palette, "plasma" to use the plasma
#'color palette, "magma" to use the magma color palette, or
#'"magma" to use the magma color palette.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import viridis
PlotAlphaData<-function(mbSetObj, data.src, bargraphName, distName, metadata, 
                        taxrank, colors = "default", format="png", dpi=72, interactive = FALSE){

  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);  
  set.seed(13133);
  
  dataName <- mbSetObj$dataSet$name;
  module.type <- mbSetObj$module.type;
  
  if(module.type == "meta"){
    data <- qs::qread("merged.data.raw.qs");
    data <- subsetPhyloseqByDataset(mbSetObj, data);

  }else{
    if(data.src == "orig"){
      data <- readDataQs("orig.phyobj", module.type, dataName);
    }else{
      data <- mbSetObj$dataSet$proc.phyobj;
    }
  }
  
  if(taxrank!="OTU"){
    #merging at taxonomy levels
    data <- fast_tax_glom_mem(data, taxrank);
    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
  } 
  
  #reordering the sample (getting reordered rownames)
  sam <- sample_data(data);
  sam <- sam[order(sam[[metadata]])];
  smplord <- rownames(sam);
  
  smpl.num <- length(smplord);
  if(smpl.num < 25){
    width <- 600
  }else if(smpl.num >= 25 & smpl.num <= 50){
    width <- 750
  }else{
    width <- 900
  }
  
  bargraphName = paste(bargraphName, ".", format, sep="");
  mbSetObj$imgSet$alpha <- bargraphName;  
  
  box = plot_richness(data, color = metadata, measures = distName) + scale_x_discrete(limits=c(smplord));
  
  if(colors == "viridis"){
    box = box + viridis::scale_color_viridis(discrete=TRUE)
  }else if(colors %in% c("magma","plasma","inferno")){
    box <- box + viridis::scale_color_viridis(option=colors, discrete=TRUE)
  }
  
  mbSetObj$analSet$alpha <- box$data;
  fast.write(mbSetObj$analSet$alpha, file="alphadiversity.csv");
  box = box + theme_bw() + theme(axis.text.x = element_text(angle = 45, hjust = 1,vjust= 1));
  box$layers <- box$layers[-1];
  box <- box + geom_point(size=3, alpha=0.7);
  #getting scale for plot (using same for boxplot also)
  ylimits <<- layer_scales(box)$y$range$range;
  
  Cairo::Cairo(file=bargraphName, width, height=450, type=format, bg="white", dpi=dpi);

  if(mbSetObj$module.type == "meta"){
     box <- box + facet_grid(. ~ dataset, scales = "free", space = "free");
  }
  print(box);
  dev.off();
  
  mbSetObj$analSet$alpha.meth <- distName;
  mbSetObj$analSet$alpha.metadata <- metadata;
  mbSetObj$analSet$alpha.taxalvl <- taxrank;
  
  if(interactive){
    library(plotly)
    library(htmlwidgets)
    p <- plotly::ggplotly(box)
    p[["x"]][["layout"]][["annotations"]][[1]][["y"]] <- -0.075 # Samples
    p[["x"]][["layout"]][["annotations"]][[2]][["x"]] <- -0.06 # Alpha
    p <- p %>% layout(margin = list(l = 20, r = 20, b = 20, t = 40))
    htmlwidgets::saveWidget(p, "alpha_richness_interactive.html")
  }
  
  return(.set.mbSetObj(mbSetObj))
}

#'Function to create bar plots of selected taxa level.
#'@description This functions creates bar plots of a selected taxa level. 
#'@param mbSetObj Input the name of the mbSetObj.
#'@param barplotName Character, input the name of the bar plot.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotSampleTaxaAundanceBar<-function(mbSetObj, barplotName, taxalvl, samplnm,
                                    imgOpt, feat_cnt, toptaxa, abunTopTaxaOpt, 
                                    appendnm, format="png", dpi=72){
  save.image("taxsample.RData")
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_reshape();
    load_phyloseq();
  }
  
  #using filtered data
  data <- mbSetObj$dataSet$filt.data;
  
  if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
    data<-otu_table(data, taxa_are_rows =TRUE);
  }
  
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  data1 <- merge_phyloseq(data, tax_table(mbSetObj$dataSet$proc.phyobj), sample_table);
  
  yLbl <- "Actual Abundance";
  
  data <- as.matrix(otu_table(data1));
  data <- data[,samplnm,drop=FALSE];
  data <- t(data);
  data_tax <- tax_table(data1);
  
  #reshaping data
  if(taxalvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    
    taxa_nm <- as.character(data_tax[,taxalvl]);
    
    y <- which(is.na(taxa_nm)==TRUE);
    #converting NA values to unassigned
    taxa_nm[y] <- "Not_Assigned";
    
    if(appendnm=="T"){
      all_nm <- colnames(tax_table(data1));
      hg_nmindx <- which(all_nm==taxalvl)-1;
      
      if(hg_nmindx!=0){
        nma <- as.character(tax_table(data1)[,hg_nmindx]);
        y1 <- which(is.na(nma)==TRUE);
        nma[y1] <- "Not_Assigned";
        nm <- paste0(nma,"_",taxa_nm);
        ind <- which(nm=="Not_Assigned_Not_Assigned");
        nm[ind] <- "Not_Assigned";
        nm <- gsub("_Not_Assigned", "",nm, perl = TRUE);
        taxa_nm <- nm;
      }
    }
  }
  
  taxa_nm <- as.matrix(taxa_nm);
  
  if(appendnm=="T"){
    y <- which(grepl("Not_Assigned", taxa_nm))
  }else{
    y <- which(is.na(taxa_nm)==TRUE);
  }
  
  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data,check.names=FALSE);
  data <- data[ , order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  colnames(data)[indx] <- "NA";
  
  if(abunTopTaxaOpt == "bottom"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
    
    if(length(ind)==0){
      AddErrMsg("All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.");
      return(0);
    }
    
    if(length(ind1)>0){
      colnames(data)[ind1] <- "Others";
      data <- as.data.frame(t(rowsum(t(data),group = colnames(data))),check.names=FALSE);
    }
    
    if(imgOpt=="barnorm"){
      data <- as.data.frame(apply(data,1, function(x) x/sum(x)),check.names=FALSE);
      data <- as.data.frame(t(data),check.names=FALSE);
      yLbl <- "Relative Abundance";
    }
    
    feat_no<-ncol(data);
    fast.write(t(data), file="taxa_abund.csv");
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step);
    data <- data[order(data[,2]),];
    data <- data[,-1];
  } else {
    if(imgOpt=="barnorm"){
      data <- as.data.frame(apply(data,1, function(x) x/sum(x)),check.names=FALSE);
      data <- as.data.frame(t(data),check.names=FALSE);
      yLbl <- "Relative Abundance";
    }
    
    fast.write(t(data), file="taxa_abund.csv");
    data$step <- factor(rownames(data));
    data <- melt(data,id='step');
    data$step <- as.numeric(data$step);
    data <- data[order(data[,2]),];
    data <- data[,-1];
    
    data <- data[which(data$value != 0), ]
    data$variable <- as.character(data$variable);
    data <- data[order(data$value, decreasing = TRUE), ];
    
    contain_notassign <- any(grepl('^NA$', data$variable));
    if(nrow(data) <= toptaxa){
      data <- data;
    } else {
      if (! contain_notassign){
        data <- head(data, n = toptaxa);
        bottom_data <- tail(data, n = nrow(data) - toptaxa);
        other_data <- data.frame("variable" = "Others", 
                                 "value" = sum(bottom_data$value),check.names=FALSE);
        data <- rbind.data.frame(top_data, other_data);
      } else {
        if(nrow(data) == toptaxa + 1){
          data <- data;
        } else {
          not_assigned_data <- data[which(data$variable == "NA"), ];
          data <- data[which(data$variable != "NA"), ];
          top_data <- head(data, n = toptaxa);
          bottom_data <- tail(data, n = nrow(data) - toptaxa);
          other_data <- data.frame("variable" = "Others", 
                                   "value" = sum(bottom_data$value),check.names=FALSE);
          data <- rbind.data.frame(top_data, other_data, not_assigned_data);
        }
      }
    }
    feat_no<-nrow(data);
  }
  
  w<-600;
  
  if(feat_no < 3){
    w<-w;
  } else if (feat_no <5){
    w<-w+100;
  } else if (feat_no < 10){
    w<-w+200;
  } else if (feat_no < 20){
    w<-w+400;
  } else if (feat_no > 20){
    w<-w+500;
  }
  
  a <- feat_no;
  
  if(length(a)<50){
    h<-feat_no*50;
  }else{
    h<-feat_no*25;
  }
  
  #sorting by descending
  rdaName = paste(barplotName, ".rda", sep="");
  jsonName = paste(barplotName, ".json", sep="")
  barplotName = paste(barplotName, ".",format, sep="");
  mbSetObj$imgSet$stack<-barplotName;
  mbSetObj$imgSet$stackRda <-rdaName;
  mbSetObj$imgSet$stackType <- "sample";
  
  Cairo::Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
  box <- ggplot(data, aes(x=reorder(variable,value),y=value))+geom_bar(stat="identity",width=0.6,fill="steelblue")+theme_bw()+
    theme(axis.text.x = element_text(angle = 0,vjust=0.5))+
    labs(y=yLbl,x="",fill=taxalvl)+coord_flip()+
    theme(panel.background = element_blank(),
          panel.grid.major = element_blank(), panel.grid.minor = element_blank(),legend.position="none");
  print(box);
  dev.off();

  save(box,file=rdaName);

  mbSetObj$analSet$stack<-data;
  mbSetObj$analSet$stack.taxalvl<-taxalvl;
  mbSetObj$analSet$plot<-"Stacked Bar";
  return(.set.mbSetObj(mbSetObj));
}

#'Function to create box plots for alpha diversity analysis
#'@description This functions performs metagenome seq analysis on the microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param boxplotName Character, input the name of the boxplot.
#'@param distName Character, input the diversity measure
#'to calculate alpha-diversity. "Chao1", "Observed", "ACE", "Shannon", 
#'"Simpson", or "Fisher".
#'@param metadata Input the name of the experimental factor to group the samples.
#'@param colors Use "default", "viridis" to 
#'use the viridis color palette, "plasma" to use the plasma
#'color palette, "magma" to use the magma color palette, or
#'"magma" to use the magma color palette.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import ggplot2
#'@import viridis
  
  PlotAlphaBoxData<-function(mbSetObj, boxplotName, distName, metadata, 
                             colors="default", format="png", dpi=72){
    mbSetObj <- .get.mbSetObj(mbSetObj);
    load_viridis();
    
    set.seed(1313397);
    data <- mbSetObj$analSet$alpha;
    CLASS <- data[,metadata];
    vals <- data[,"value"]
    boxplotName = paste(boxplotName, ".", format, sep="");
    mbSetObj$imgSet$alpha.box<-boxplotName;
    
    grp_size <- length(levels(CLASS))
    
    if(grp_size <= 2){
      width <- 500;
    }else if(grp_size >= 3 & grp_size < 6){
      width <- 600;
    }else{
      width <- 700;
    }
    
    if(colors %in% c("magma","plasma","inferno","viridis")){
      box1 = ggplot(data, aes(CLASS, vals, fill = CLASS))
    }else{
      box1 = ggplot(data, aes(CLASS, vals, color = CLASS))
    }
    
    box1 = box1 + stat_boxplot(geom ='errorbar', width=0.2) +
      geom_boxplot(alpha=0.7, outlier.shape = NA,
                   position = position_dodge(width = 0), width=0.3) +
      geom_jitter(width = 0.1) + #reduce the width of jitter
      stat_summary(fun.y=mean, #add mean point
                   geom = "point",
                   shape = 18,
                   size = 4,
                   color = "black") +
      theme_bw() +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 11),
            axis.text.y = element_text(size = 11),
            legend.text = element_text(size = 11)) + #adjust titles
      labs(title = "",
           y= paste("Alpha-diversity Index:", as.character(data$variable[[1]]), sep = " "),
           x="") + #remove x = CLASS, add title name, change y name
      theme(axis.title.y = element_text(size=14)) + 
      coord_cartesian(ylim = c(ylimits[1], ylimits[2]));
    
    if(colors == "viridis"){
      box1 <- box1 + viridis::scale_fill_viridis(discrete = TRUE);
    }else if(colors %in% c("magma","plasma","inferno")){
      box1 <- box1 + viridis::scale_fill_viridis(option=colors, discrete=TRUE)
    }
    
    Cairo::Cairo(file=boxplotName, width=width, height=400, type=format, bg="white", dpi=dpi);
    if(mbSetObj$module.type == "meta"){
      box1 <- box1 + facet_grid(. ~ dataset, scales = "free", space = "free");
    }
    print(box1);
    dev.off();
    
    return(.set.mbSetObj(mbSetObj))
  }
  

######################################
###########Beta-diversity#############
######################################

#'Function to plot beta diversity.
#'@description This functions creates beta diversity plots.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param plotNm Character, input the name of the beta-diversity score plot.
#'@param ordmeth Ordination method. Character, input "PCoA" to create a PCoA plot and
#'"NMDS" to create a NMDS plot.
#'@param distName Character, input the name of the distance method. "bray" for
#'Bray-Curtis Index, "jsd" for Jensen-Shannon Divergence, "jaccard" for Jaccard Index,
#'"unifrac" for Unweighted Unifrac Distance and "wunifrac" for Weighted Unifrac Distance.
#'@param colopt Character, input whether the data points should be colored by
#'experimental factor with "expfac", taxon abundance with "taxa" or alpha diversity 
#'with "alphadiv".
#'@param metadata Input the name of the preferred experimental factor, only used if 
#'colopt is "expfac".
#'@param showlabel Character, input whether or not to label samples in the plot.
#'"none" to label no samples, "samnm" to label samples by their name, and "Class" for
#'their group classification.
#'@param taxrank Character, input the taxonomic level to calculate
#'beta-diversity.
#'@param taxa Character, input the specific taxon used to color the data points. Only
#'used if colopt is set to "taxa".
#'@param alphaopt Character, input the name of the alpha-diversity metric. "Chao1",
#' "Observed", "ACE", "Shannon", "Simpson", or "Fisher".
#'@param ellopt Character, input "yes" to show ellipses and "no" to not.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@param custom_col Set to "none" to use the default color palette, "viridis" to use the 
#'default viridis color palette, "plasma" to use the plasma viridis color palette and
#'"magma" to use the magma color palette.
#'@param interactive Boolean, if set to TRUE, saves the plot
#'as an interactive html plot.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import data.table
#'@import ape

PerformBetaDiversity <- function(mbSetObj, plotNm, ordmeth, distName, colopt, metadata, 
                                 showlabel, taxrank, taxa, alphaopt, ellopt, comp.method, format="png", dpi=72,
                                 custom_col = "none",pairwise, interactive = FALSE){
 
  combined <- F;
  mbSetObj <- .get.mbSetObj(mbSetObj);
  module.type <- mbSetObj$module.type;

  load_datatable();
  load_viridis();
  load_phyloseq();

  set.seed(13134);
  if(all(c(module.type == "meta", !combined))){
    mdata.all <- mbSetObj$mdata.all;
    sel.nms <- names(mdata.all)[mdata.all==1];
    dataNames <- sel.nms;
  }else{
    dataNames <- mbSetObj$dataSet$name
  }
  ord.list <- list();
  for(i in 1:length(dataNames)){
    dataName <- dataNames[i];
    if(module.type == "meta"){
      mbSetObj$dataSet <- readDataset(dataName);
    }
    if(all(c(module.type == "meta", combined))){
      proc.phyobj <- qs::qread("merged.data.raw.qs");
      norm.phyobj <- qs::qread("merged.data.qs");
      proc.phyobj <- subsetPhyloseqByDataset(mbSetObj, proc.phyobj);
      norm.phyobj <- subsetPhyloseqByDataset(mbSetObj, norm.phyobj);
      
    }else{
      proc.phyobj <- mbSetObj$dataSet$proc.phyobj;
      norm.phyobj <- mbSetObj$dataSet$norm.phyobj;
    }    
    #using normalized data
    
    phyloseq_objs <- readDataQs("phyloseq_objs.qs",mbSetObj$module.type,dataName)

    data <- phyloseq_objs$merged_obj[[taxrank]]
    if(is.null(data)){
      AddErrMsg("Errors in projecting to the selected taxanomy level!");
      return(0);
    }
    
    if(colopt=="taxa"){
      if(taxrank=="OTU"){
        data1 <- as.matrix(otu_table(data));
        feat_data <- as.numeric(data1[taxa,]);
      }else{
        nm <- as.character(tax_table(data)[,taxrank]);
        if(sum(is.na(nm))/length(nm) > 0.7){
          AddErrMsg("More than 70% values are missing at this taxa!");
          return(0);
        }
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
      data1 <- proc.phyobj;
      box <- plot_richness(data1, measures = alphaopt);
      alphaboxdata <- box$data;
      sam_nm <- sample_names(data);
      alphaboxdata <- alphaboxdata[alphaboxdata$samples %in% sam_nm,];
      alphaval <- alphaboxdata$value;
      sample_data(data)$alphaopt <- alphaval;
      indx <- which(colnames(sample_data(data))=="alphaopt");
      colnames(sample_data(data))[indx] <- alphaopt;
    }else if(colopt=="continuous") {
      require("MMUPHin");
      require("vegan");
      
      data <- proc.phyobj;
      #sub_sam_data <- sam_data[which(sam_data[,metadata] == meta.grp), ]
      #sub_data <- data@otu_table[, rownames(sub_sam_data)];
      sub_sam_data <- sample_data(data); 
      
      sub_data <- as.matrix(otu_table(data)); #data@otu_table[, rownames(sub_sam_data)];
      sub_data <- apply(sub_data, 2, function(x) x / sum(x))
      dist.data <- vegdist(t(sub_data));
      
      fit_continuous <- continuous_discover(feature_abd = sub_data,
                                            batch = "dataset",
                                            data = sub_sam_data,
                                            control = list(var_perc_cutoff = 0.5,
                                                           verbose = FALSE));
      
      loading <- data.frame(feature = rownames(fit_continuous$consensus_loadings),
                            loading1 = fit_continuous$consensus_loadings[, 1]);
      sample_data(data)$loading <- fit_continuous$consensus_scores[, 1];
      oldDF <- as(sample_data(data), "data.frame");
      newDF <- oldDF[oldDF$dataset %in% dataName, ];
      sample_data(data) <- sample_data(newDF);
      
    }else{
      data<-data;
    }
    if(distName=="wunifrac"){
      
      load_ape();
      pg_tree <- qs::qread("tree.qs");
      pg_tb <- tax_table(data);
      pg_ot <- otu_table(data);
      pg_sd <- sample_data(data);
      pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
      data <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree);
      # check if pg_tree is null
      if(is.null(data@phy_tree)){
        AddErrMsg("Tree tip labels do not match feature names in the OTU/taxonomic tables!");
        return(0)
      }
      
      qs::qsave(data, "data_unifra.qs");
      ord <- ordinate(data,method = ordmeth,"unifrac",weighted=TRUE);
      
    } else if (distName=="unifrac") {
      
      load_ape();      
      pg_tree <- qs::qread("tree.qs");
      pg_tb <- tax_table(data);
      pg_ot <- otu_table(data);
      pg_sd <- sample_data(data);
      pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
      data <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree);
      
      # check if pg_tree is null
      if(is.null(data@phy_tree)){
        AddErrMsg("Tree tip labels do not match feature names in the OTU/taxonomic tables!");
        return(0)
      }
      
      qs::qsave(data, "data_unifra.qs");
      ord <- ordinate(data, method = ordmeth,"unifrac");
    }else{
      ord <- ordinate(data, method = ordmeth,distName);
    }

      if(ordmeth == "NMDS"){
        ord$vectors <- ord$points;
        colnames(ord$vectors) <- c("Axis.1", "Axis.2")
      }

    if(colopt == "continuous"){
      ord$vectors <- as.data.frame(ord$vectors);
      ord$vectors <- cbind(sample_data(data)$loading, ord$vectors);
      colnames(ord$vectors)[1] <- "loading"
    }else if(all(c(colopt == "alphadiv", module.type == "meta"))){
      ord$vectors <- as.data.frame(ord$vectors);
      
      ord$vectors <- cbind(sample_data(data)[,alphaopt], ord$vectors);
      colnames(ord$vectors)[1] <- alphaopt
      
    }
    .set.mbSetObj(mbSetObj);
    PerformCategoryComp(mbSetObj, taxrank, comp.method,distName, metadata,pairwise);
    mbSetObj <- .get.mbSetObj(mbSetObj);
    ord$stat.info <- mbSetObj$analSet$stat.info;
    
    ord.list[[dataName]] <- ord;
    
  }
  plotNm = paste(plotNm, ".", format, sep="");
  mbSetObj$imgSet$beta2d<-plotNm;
  
  if(all(c(module.type == "meta", !combined))){
    dat.num <- length(ord.list);
    if(dat.num > 3){
      colNum=2;
    }else{
      colNum=1;
    }
    require(plyr);
    merged.data <- qs::qread("merged.data.qs");
    merged.data <- subsetPhyloseqByDataset(mbSetObj, merged.data);
    sam_data <- sample_data(merged.data);
    pdataframe = ldply(ord.list, function(x){
      df = sam_data;
      return(x$vectors)
    })
    pdataframe <- cbind(sam_data, pdataframe);
    res <- lapply(ord.list, function(x){ return(x$stat.info)})
    stats <- unlist(res);
    for(i in 1:length(stats)){
      stats[i] <- paste0(dataNames[i], ": ", stats[i]);
    }
    names(stats) <- dataNames;
    
    if(colopt =="continuous") {
      box <- ggplot(pdataframe, aes(Axis.1, Axis.2, color=loading)) + scale_colour_gradient2(low="green",mid="black", high="red");
    }else if(colopt =="alphadiv"){
      box <- ggplot(pdataframe, aes(Axis.1, Axis.2, color=pdataframe[,alphaopt])) + scale_colour_gradient2(low="green", high="red") + labs(color = alphaopt);
    }else{
      box <- ggplot(pdataframe, aes(Axis.1, Axis.2, color=study_condition, fill=study_condition));
    }
    box <- box + geom_point(size=4) +
      facet_wrap(~dataset, scales="free", ncol= colNum, labeller = labeller(dataset = stats));
    if(colNum == 1){
      height.multiplier <- dat.num;
      width <- 720;
    }else{
      height.multiplier <- ceiling(dat.num/2)
      width <- 1220;
    }
    height <- 100 + 400 * height.multiplier;
    
  }else{
    sam_data <- as.data.frame(sample_data(data),check.names=FALSE);
    if(colopt=="taxa"){
      box = plot_ordination(data, ord, color=taxa) + labs(aesthetic=taxaorig) + scale_colour_gradient(low="green", high="red");
    }else if(colopt=="alphadiv") {
      box = plot_ordination(data, ord, color=alphaopt)+ scale_colour_gradient(low="green", high="red");
    }else if(colopt=="continuous") {
      box = plot_ordination(data, ord, color="loading")+ scale_colour_gradient2(low="green", high="red");
    }else{
      box = plot_ordination(data, ord, color=metadata);
    }  
    height <- 500;
    width <- 720;
  }
  
  box$layers <- box$layers[-1];
  
  if(showlabel=="samnm"){
    box = box + geom_text(aes(label=sample_id), hjust=0.5, vjust=2, size=3, fontface="bold");
    box = box + geom_point(size=4, alpha=0.6) + theme_bw();
  }else if(showlabel=="none"){
    box=box+geom_point(size=4, alpha=0.8) + theme_bw();
  }else{
    showlabel <<- showlabel;
    bx_data <<- data.frame(box$data,check.names=FALSE);
    box = box + geom_text(aes(label=bx_data[ ,showlabel]), hjust=0.5, vjust=2, size=3, fontface="bold");
    box = box + geom_point(size=4, alpha=0.6) + theme_bw();
  }
  
  #used for area color for ellipse
  if(colopt=="expfac"){
    clsLbl <- quo(sam_data[[metadata]]);
    
    if (ellopt=="yes"){
      box = box + stat_ellipse(type="norm", linetype=2, geom = "polygon", alpha = 0.2, aes(fill = !!clsLbl), show.legend=FALSE);
    }
  }
  
  if(custom_col == "viridis"){
    box = box + viridis::scale_color_viridis(discrete = TRUE) + viridis::scale_fill_viridis(discrete = TRUE)
  }else if(custom_col == "plasma"){
    box = box + viridis::scale_color_viridis(option="plasma", discrete = TRUE) + viridis::scale_fill_viridis(option="plasma", discrete = TRUE)
  }else if(custom_col == "magma"){
    box = box + viridis::scale_color_viridis(option="magma", discrete = TRUE) + viridis::scale_fill_viridis(option="magma", discrete = TRUE)
  }
  
  box = box + theme(strip.text.x = element_text(size = 12));

  Cairo::Cairo(file=plotNm, width=width, height=height, type=format, bg="white",dpi=dpi);
  print(box);
  dev.off();
  
  #saving info for report generation
  mbSetObj$analSet$beta <- data;
  mbSetObj$analSet$beta.meth <- ordmeth;
  mbSetObj$analSet$beta.dist <- distName;
  mbSetObj$analSet$beta.taxalvl <- taxrank;
  
  # if it is NMDS, show stress value in the plot
  if(ordmeth == "NMDS"){
    mbSetObj$analSet$beta.stress <- paste("[NMDS] Stress =", signif(ord$stress, 5));
  }
  return(.set.mbSetObj(mbSetObj))
}

#'Plot functional annotation summary
#'@description This functions plots the functional annotation summary.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotFunAnotSummary<-function(mbSetObj, imgName, format="png",funanno, dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  set.seed(280561499);
   
  if(funanno=="gg12"){
    func.file = "func.pred.gg12"
    mbSetObj$analSet$funcode <- "picrust_ko";
  }else if(funanno=="gg13"){
    func.file = "func.pred.gg13"
    mbSetObj$analSet$funcode <- "picrust_ko";
  }else if(funanno =="tax4fun2"){
    func.file = "func.pred.tax4fun2";
    mbSetObj$analSet$funcode <- "tax4fun2_ko";
  }else{
    func.file = "func.pred.tax4fun"
    mbSetObj$analSet$funcode <- "tax4fun_ko";
  }
  
  if(is.null(mbSetObj$analSet[[func.file]])){
    result <- qs::qread(func.file);
  }else{
    result <- mbSetObj$analSet[[func.file]];
    qs::qsave(mbSetObj$analSet[[func.file]], file=func.file);
    mbSetObj$analSet[[func.file]] <- NULL;
  }
 
  imgName = paste(imgName, ".",format, sep="");
  mbSetObj$imgSet[[func.file]] <- imgName;
  Cairo::Cairo(file=imgName, width=900, height=480, type=format, bg="white",dpi=dpi);
  box <- ggplot(stack(log(result)), aes(x = factor(ind, levels = names(result)), y=values)) + labs(x=NULL, y="log (KO Counts)") + geom_boxplot() + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  print(box);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj))
}

#'Plot functional annotation summary
#'@description This functions plots the functional annotation summary.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param barplotName Character, input the name of the barplot.
#'@param viewOpt Character.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param metadata If users wish to merge samples to groups in the stacked bar plot,
#'set this to the preferred grouping.  
#'@param feat_cnt Set the minimum feature count that is used to "bin"
#'together small taxa. 
#'@param colpalopt Select the color palette options. "set3", 
#'which is the Set3 from the R Color Brewer, "cont21" which is
#'a set of 21 colors, "cont28" which is a set of 28 colors, and 
#'"cont42" which is a set of 42 colors. For users who wish to use
#'a color palette robust to colorblindness 
#'/(https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html/), 
#'use "viridis", "magma", "plasma", or "inferno".
#'@param calcmeth Merge small taxa based on the sum /("sum"/)
#'or median /("med"/) counts across all samples or groups.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
#'@import viridis
PlotTaxaAbundanceArea<-function(mbSetObj, barplotName, viewOpt, taxalvl, metadata,
                                feat_cnt, colpalopt, calcmeth, toptaxa, abunTopTaxaOpt, 
                                appendnm=FALSE, format="png", dpi=72){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_reshape();
  load_viridis();
  load_phyloseq();
  
  if(mbSetObj$module.type == "meta"){
    if(viewOpt=="smpl_grp"){
      AddErrMsg("Taxa area can not be plotted for merged samples.");
      return(0);
    }
    data1 <- qs::qread("merged.data.raw.qs");
    data1 <- subsetPhyloseqByDataset(mbSetObj, data1);

    sample_table <- sample_data(data1);
    data <- as.data.frame(t(otu_table(data1)),check.names=FALSE);
  }else{
    #using filtered data
    data <- mbSetObj$dataSet$filt.data;
    if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
      data <- otu_table(data,taxa_are_rows =TRUE);
    }
    sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
    data1 <- merge_phyloseq(data, tax_table(mbSetObj$dataSet$proc.phyobj), sample_table);
  }
  
  if(viewOpt=="smpl_grp"){
    data <- as.data.frame(t(otu_table(data1)),check.names=FALSE);
    data <- cbind(data,variable=data1@sam_data[[metadata]]);
    data <- aggregate(. ~variable,data,sum);
    gp_nm <- rownames(data)<-data[,1];
    data <- data[,-1];
    data <- data[order(rownames(data)),];
    clsLbl <- sort(unique(factor(data1@sam_data[[metadata]])));
    colvec <- NULL;
  }else {
    data <- data.frame(otu_table(data1),check.names=FALSE);
    
    # reorder data based on groups
    sam <- sample_data(data1);
    smpl_nm <- row.names(sam);
    
    if(metadata == "none"){
      sam$newnewnew <- rep("one", nrow(sam));
      metadata <- "newnewnew";
    }
    clsLbl <- factor(sam[[metadata]]);
    if(all(c(length(levels(clsLbl)) > 9, min(table(clsLbl)) < 3))){
      AddErrMsg("Too many facets to be displayed - please select a more meaningful facet option with at least 3 samples per group.");
      return(0);
    }
    ord.inx <- order(clsLbl);
    smpl_nm <- smpl_nm[ord.inx];
    clsLbl <- clsLbl[ord.inx];
    colvec <- as.numeric(clsLbl)+1;
    data <- t(data[,ord.inx]);
  }
  
  data_tax <- tax_table(data1);
  
  #reshaping data
  if(taxalvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    
    taxa_nm <- as.character(data_tax[,taxalvl]);
    
    y <- which(is.na(taxa_nm)==TRUE);
    #converting NA values to unassigned
    taxa_nm[y] <- "Not_Assigned";
    
    if(appendnm=="T"){
      all_nm <- colnames(tax_table(data1));
      hg_nmindx <- which(all_nm==taxalvl)-1;
      
      if(hg_nmindx!=0){
        nma <- as.character(tax_table(data1)[,hg_nmindx]);
        y1 <- which(is.na(nma)==TRUE);
        nma[y1] <- "Not_Assigned";
        nm <- paste0(nma,"_",taxa_nm);
        ind <- which(nm=="Not_Assigned_Not_Assigned");
        nm[ind] <- "Not_Assigned";
        nm <- gsub("_Not_Assigned", "",nm, perl = TRUE);
        taxa_nm <- nm;
      }
    }
  }
  
  taxa_nm <- as.matrix(taxa_nm);
  
  if(appendnm=="T"){
    y <- which(grepl("Not_Assigned", taxa_nm))
  }else{
    y <- which(is.na(taxa_nm)==TRUE);
  }
  
  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- as.matrix(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data,check.names=FALSE);
  data <- data[ ,order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  colnames(data)[indx] <- "NA";
  
  if(abunTopTaxaOpt == "top"){
    if(calcmeth=="sum"){
      feat_cnt <- sort(sapply(data, sum), decreasing = TRUE)[toptaxa];
    }else {
      feat_cnt <- sort(sapply(data, median), decreasing = TRUE)[toptaxa];
    }
  }  
  
  if(calcmeth=="sum"){
    ind<-which(colSums(data)>feat_cnt);
    ind1<-which(colSums(data)<feat_cnt);
  } else {
    dt<-apply(data,2,median);
    ind<-which(dt>feat_cnt);
    ind1<-which(dt<feat_cnt);
  }
  
  if(length(ind)==0){
    AddErrMsg("All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.");
    return(0);
  }
  
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data<- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)),check.names=FALSE);
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
  #a <- nsamples(data1);
  a <- nrow(sample_table);
  min.w <- 540;
  if(length(a)<50){
    w <- a*35;
  }else{
    w <- a*25;
  }
  if(w < min.w){
    w <- min.w;
  }
  
  fast.write(t(data), file="taxa_abund.csv");
  data$facetOpt <- as.character(clsLbl);
  data$step <- factor(rownames(data), levels = rownames(data));
  data <- melt(data,id=c('step', 'facetOpt'));
  data$step <- as.numeric(data$step);
  # data$sample <- data$step;
  data <- data[order(data[,3]),];
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
  levels(data$variable) <- sub("^X", "", levels(data$variable)) 
  
  
  data$facetOpt <- factor(data$facetOpt,levels=(unique(data$facetOpt)) )
  
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
  
  jsonName = paste(barplotName, ".json", sep="");
  rdaName = paste(barplotName, ".rda", sep="");
  barplotName = paste(barplotName, ".",format, sep="");
  mbSetObj$imgSet$stack <- barplotName;
  mbSetObj$imgSet$stackRda <-rdaName;
  mbSetObj$imgSet$stackType <- "area";

  Cairo::Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
  
  box <- ggplot(data,aes(x=step,y=value)) + theme_bw() +
    theme(axis.text.x = element_text(angle = 90, hjust =1,vjust=0.5)) +
    geom_area(aes(fill=variable), position='fill') +
    facet_grid(~ facetOpt, space = "free", scales = "free") +
    scale_x_continuous(breaks=seq(1,length(unique(data$step)),1),labels=lbl) +
    labs(y="Relative abundance",x="",fill=taxalvl) +
    theme(legend.position="bottom",legend.box = "vertical") +
    theme(axis.text.x = element_text(colour="black"),axis.title.x=element_blank());
  
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
  }else if(colpalopt=="viridis"){
    box <- box + viridis::scale_fill_viridis(discrete=TRUE)
  }else if(colpalopt %in% c("magma","plasma","inferno")){
    box <- box + viridis::scale_fill_viridis(option=colpalopt, discrete=TRUE)
  }else{
    box <- box + scale_fill_manual(values=c(x.colors))
  }
  
  if(metadata == "newnewnew"){
    box <- box + theme(strip.text.x = element_blank())
  }
  if(mbSetObj$module.type == "meta"){
    #box <- box + facet_grid(. ~ dataset  , scales = "free", space = "free");
  }
  
  print(box);
  dev.off();
  # for plotly
  save(box,file=rdaName);
  
  mbSetObj$analSet$stack<-data;
  mbSetObj$analSet$stack.taxalvl<-taxalvl;
  mbSetObj$analSet$plot<-"Stacked Area";
  
  return(.set.mbSetObj(mbSetObj));
  
}

#'Function to plot bar charts for alpha diversity
#'@description This functions plots bar charts of different taxonomic levels
#'for alpha diversity.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param barplotName Character, input the name of the barplot.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param metadata If users wish to merge samples to groups in the stacked bar plot,
#'set this to the preferred grouping.  
#'@param facet Character, set the group to separate the bar plots.
#'@param imgOpt Character, set the graph type. Stacked bar
#'using the actual abundance, use "barraw". Stacked bar using
#'the percentage abundance, use "barnorm".
#'@param feat_cnt Set the minimum feature count that is used to "bin"
#'together small taxa. 
#'@param colpalopt Select the color palette options. "set3", 
#'which is the Set3 from the R Color Brewer, "cont21" which is
#'a set of 21 colors, "cont28" which is a set of 28 colors, and 
#'"cont42" which is a set of 42 colors. For users who wish to use
#'a color palette robust to colorblindness 
#'/(https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html/), 
#'use "viridis", "magma", "plasma", or "inferno".
#'@param calcmeth Merge small taxa based on the sum /("sum"/)
#'or median /("med"/) counts across all samples or groups.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@param interactive Boolean, if TRUE, will save an interactive
#'version of the Stacked Bar plot using plotly.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
#'@import ggplot2
#'@import viridis

PlotTaxaAundanceBar<-function(mbSetObj, barplotName, taxalvl, facet, facet2, imgOpt, 
                              feat_cnt, colpalopt, calcmeth, toptaxa, abunTopTaxaOpt, 
                              appendnm, format="png", dpi=72, interactive = FALSE){
  save.image("tax.RData");
  load_reshape();
  load_ggplot();
  load_viridis();
  load_phyloseq(); 

  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(mbSetObj$module.type == "meta"){
    data1 <- qs::qread("merged.data.raw.qs");
    #data1 <- subsetPhyloseqByDataset(mbSetObj, data1);
    sample_table <- sample_data(data1);
    data <- as.data.frame(otu_table(data1),check.names=FALSE);
    #facet2="dataset";
  }else{
    data <- mbSetObj$dataSet$filt.data;
    if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
      data <- otu_table(data, taxa_are_rows =TRUE);
    }
    sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
    data1 <- merge_phyloseq(data, tax_table(mbSetObj$dataSet$proc.phyobj), sample_table);
    data <- as.data.frame(otu_table(data1),check.names=FALSE);
  }
  # reorder data based on groups
  sam <- sample_data(data1);
  
  if(facet == "none"){
    sam$newnewnew <- rep("one", nrow(sam));
    facet <- "newnewnew";
  }
  
  clsLbl <- factor(sam[[facet]]);
  if(length(clsLbl)==0){
    AddErrMsg("Invalid class label! Ensure there are no special characters (i.e. punctuation marks) in your group names!")
    return(0)
  }
  
  if(all(c(length(levels(clsLbl)) > 9, min(table(clsLbl)) < 3))){
    AddErr("Too many facets to be displayed - please select a more meaningful facet option with at least 3 samples per group.");
    return(0);
  }
  
  sample_data(data1) <- sam;
  
  metalp <- as(sample_data(data1), "data.frame") #extract metadata table
  smpl_nm <- sample_names(data1);
  
  ord.inx <- order(clsLbl);
  smpl_nm <- smpl_nm[ord.inx];
  clsLbl <- clsLbl[ord.inx];
  colvec <- as.numeric(clsLbl)+1;
  data <- t(data[,ord.inx]);
  
  data_tax <- tax_table(data1);

  if(taxalvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    
    taxa_nm <- as.character(data_tax[,taxalvl]);
    y <- which(is.na(taxa_nm)==TRUE | taxa_nm == "");
    
    #converting NA and empty values to unassigned
    taxa_nm[y] <- "Not_Assigned";
    
    if(appendnm=="T"){
      all_nm <- colnames(tax_table(data1));
      hg_nmindx <- which(all_nm==taxalvl)-1;
      
      if(hg_nmindx!=0){
        nma <- as.character(tax_table(data1)[,hg_nmindx]);    
        y1 <- which(is.na(nma)==TRUE | nma == "");
        nma[y1] <- "Not_Assigned";
        nm <- paste0(nma,"_",taxa_nm);
        ind <- which(nm=="Not_Assigned_Not_Assigned");
        nm[ind] <- "Not_Assigned";
        nm <- gsub("_Not_Assigned", "",nm, perl = TRUE);
        taxa_nm <- nm;
      }
    }
  }
  
  #reshaping data
  taxa_nm <- as.matrix(taxa_nm);
  
  if(appendnm=="T"){
    y <- which(grepl("Not_Assigned", taxa_nm))
  }else{
    y <- which(is.na(taxa_nm)==TRUE | taxa_nm== "");
  }
  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- as.matrix(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data,check.names=FALSE);
  
  data <- data[ , order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  
  colnames(data)[indx] <- "NA";
  
  
  if(abunTopTaxaOpt == "top"){
    if(calcmeth=="sum"){
      feat_cnt <- sort(sapply(data, sum), decreasing = TRUE)[toptaxa];
    }else {
      feat_cnt <- sort(sapply(data, median), decreasing = TRUE)[toptaxa];
    }
  }
  
  if(calcmeth=="sum"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
  } else {
    dt <- apply(data,2,median);
    ind <- which(dt>feat_cnt);
    ind1 <- which(dt<feat_cnt);
  }
  
  if(length(ind)==0){
    AddErrMsg("All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.");
    return(0);
  }
  
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)),check.names=FALSE);
  }
  
  yLbl <- "Actual Abundance";
  
  if(imgOpt=="barnorm"){
    data <- as.data.frame(apply(data,1, function(x) x/sum(x)),check.names=FALSE);
    data<-as.data.frame(t(data),check.names=FALSE);
    yLbl <- "Relative Abundance";
  }
  
  feat_no <- ncol(data);
  h <- 630;
  if(feat_no < 10){
    h <- h;
  } else if (feat_no < 20){
    h <- h+100;
  } else if (feat_no < 50){
    h <- h+200;
  } else if (feat_no < 100){
    h <- h+300;
  } else {
    h <- 1100;
  }
  w <- h + 100;

  data <- data[row.names(metalp),]
  data[[get("facet")]] <- metalp[[get("facet")]]
  data$sample <- row.names(data);
  fast.write(t(data), file="taxa_abund.csv");
  data <- melt(data, id = c("sample", get("facet")))
  tmp_df <- aggregate(data$value, by=list(data$variable), FUN=mean)
  var_level <- tmp_df[order(tmp_df$x, decreasing = TRUE), ][[1]]
  
  data$variable <- factor(data$variable, levels = var_level) # change the factor level of taxa
  levels(data$variable) <- sub("^X", "", levels(data$variable));
  
  if(all(c(facet2 != "null", facet2 != "none", facet2 != facet))){
    data$variable2 <- sam[[facet2]][match(data$sample, row.names(sam))]; 
  }
  
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
  
  jsonName = paste(barplotName,".json", sep="");
  rdaName = paste(barplotName,".rda", sep="");
  barplotName = paste(barplotName, ".",format, sep="");
  
  mbSetObj$imgSet$stack <- barplotName;
  mbSetObj$imgSet$stackRda <-rdaName;
  mbSetObj$imgSet$stackType <- "default";

  
  if(length(unique(data$sample)) <= 10){
    guide_num = 3
  }else if (length(unique(data$sample)) <= 20){
    guide_num = 4
  }else{
    guide_num = 5
  }
  
  box <- ggplot(data = data,
                aes(x = sample,
                    y = value,
                    fill = variable)) +
    geom_bar(stat="identity", position="stack") +
    labs(y=yLbl,x="", fill = taxalvl)+
    theme(legend.position="bottom",legend.box = "vertical",
          axis.text.x = element_text(angle = 45,
                                     hjust =1,vjust=1), #remove colour=colvec
          axis.title.x=element_blank(),
          panel.background = element_rect(fill = "grey90",
                                          colour = "grey90",
                                          size = 0.5, linetype = "solid"),
          axis.line = element_line(colour = "lightgrey"));
  
  if(any(c(facet2 == "null", facet2 == "none", facet2 == facet))){
    box <- box +
      facet_grid(as.formula(paste("~", facet)), scales = "free", space = "free") 
  } else {
    box <- box +
      facet_grid(as.formula(paste("variable2", "~", facet)), scales = "free", space = "free") 
  }
  
  if(colpalopt %in% c("set1", "set2", "set3")){
    
    cols.needed <- length(unique(data$variable))
    colpalopt <- tools::toTitleCase(tolower(colpalopt))
    
    if(colpalopt == "Set1"){
      n.cols <- 9
    }else if(colpalopt == "Set2"){
      n.cols <- 8
    }else{
      n.cols <- 12
    }
    
    if(cols.needed > n.cols){
      col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n.cols, colpalopt))
      box <- box + scale_fill_manual(values=col.func(cols.needed), guide = guide_legend(direction = "horizontal", ncol=guide_num))
    }else{
      box <- box + scale_fill_brewer(palette = colpalopt, guide = guide_legend(direction = "horizontal", ncol=guide_num))
    }
  } else if(colpalopt=="viridis"){
    box <- box + viridis::scale_fill_viridis(discrete=TRUE) + guides(color=guide_legend(ncol=guide_num)) 
  }else if(colpalopt %in% c("magma","plasma","inferno")){
    box <- box + viridis::scale_fill_viridis(option=colpalopt, discrete=TRUE) + guides(color=guide_legend(ncol=guide_num))
  }else{
    box <- box + scale_fill_manual(values=c(x.colors)) + guides(color=guide_legend(ncol=guide_num))
  }
  
  if(facet == "newnewnew"){
    box <- box + theme(strip.text.x = element_blank())
  }
  
  save(box,file=rdaName);

  mbSetObj$analSet$stack <- data;
  mbSetObj$analSet$stack.taxalvl <- taxalvl;
  mbSetObj$analSet$plot <- "Stacked Bar";
  
    return(.set.mbSetObj(mbSetObj));
  }


#'Function to perform categorical comparison.
#'@description This functions performs categorical comparisons.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxaLvl Input the name of the taxonomic level to calculate Beta-diversity comparison.
#'@param method Statistical method to calculate
#'beta-diversity significance. Use "adonis" for Permutational MANOVA, "anosim" for
#'Analysis of Group Similarities and "permdisp" for
#'Homogeneity of Group Dispersions.
#'@param distnm Character, input the name of the distance method. "bray" for
#'Bray-Curtis Index, "jsd" for Jensen-Shannon Divergence, "jaccard" for Jaccard Index,
#'"unifrac" for Unweighted Unifrac Distance and "wunifrac" for Weighted Unifrac Distance.
#'@param variable Input the name of the experimental factor to group the samples.
#'@param covariates Boolean. TRUE to consider covariates, FALSE to only consider the
#'main effect. Only valid for PERMANOVA.
#'@param cov.vec Character vector. Input the names of the covariates to consider in the 
#'PERMANOVA model.
#'@param model.additive Boolean. If TRUE, the model will be additive (i.e. data ~ var1 + var2),
#'making the assumption that the two factor variables are independent. 
#'However, if FALSE, the model will consider the synergistic effects of the variables - interaction (i.e. data ~ var1*var2).
#'"Different explanatory variables the effects on the outcome of a change in one variable may either not depend on the
#'level of the other variable (additive model) or it may depend on the level of the other variable (interaction model)." 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import vegan
PerformCategoryComp <- function(mbSetObj, taxaLvl, method, distnm, variable, pairwise=FALSE,
                                covariates = FALSE, cov.vec = NA, model.additive = TRUE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_vegan();

  if(distnm %in% c("wunifrac", "unifrac")) {
    data <- qs::qread("data_unifra.qs");
  } else {
     if(!(exists("phyloseq_objs"))){
      phyloseq_objs <- readDataQs("phyloseq_objs.qs",mbSetObj$module.type,mbSetObj$dataSet$name)
   }
     data <- phyloseq_objs$merged_obj[[taxaLvl]]
     if(is.null(data)){
        AddErrMsg("Errors in projecting to the selected taxanomy level!");
        return(0);
      }

  }
  
  data <- transform_sample_counts(data, function(x) x/sum(x));

  data.dist <- phyloseq::distance(data, method=distnm);
  group <- phyloseq::get_variable(data,variable);
  stat.info <- "";
  resTab <- list();
 
  if(method=="adonis"){
    sampledf <- data.frame(sample_data(data), check.names=FALSE); 
    outcome <- "data.dist";
    variables <- variable;
    f <- as.formula(paste(outcome,paste(variables, collapse = " * "), sep = " ~ "));
    
    # advanced feature to include co-variates and different model options (NOT implemented in web?)
    if(covariates){
      variables <- c(variable, cov.vec)
      
      if(model.additive){
        f <- as.formula(
          paste(outcome,
                paste(variables, collapse = " + "),
                sep = " ~ "))
      }else{
        f <- as.formula(
          paste(outcome,
                paste(variables, collapse = " * "),
                sep = " ~ "))
      }
    }
    res <- adonis2(formula = f, data = sampledf);
    resTab <- res[1,];
    stat.info <- paste("[PERMANOVA] F-value: ", signif(resTab$F, 5),  "; R-squared: ", signif(resTab$R2, 5), "; p-value: ", signif(resTab$Pr, 5), sep="");   
    stat.info.vec <- c(signif(resTab$F, 5), signif(resTab$R2, 5), signif(resTab$Pr, 5));
    names(stat.info.vec) <- c("F-value", "R-squared", "p-value");
 
    if(pairwise != "false"){
        grp <- sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]
       if(length(levels(grp)) == 2){ # same as before  
            res <- data.frame(t(c(stat.info.vec, p.adj=signif(resTab$Pr, 5))));
            rownames(res) <- paste(levels(grp)[1], ".vs.", levels(grp)[2]);
       }else{
            res <- .permanova_pairwise(x = data.dist,  grp);
            rownames(res) <- res$pairs;
            res$pairs <- NULL;
      }

      mbSetObj$analSet$beta.stat.pair <- mbSetObj$analSet$resTable <- signif(res,5);
      fast.write(mbSetObj$analSet$resTable, file="pairwise_permanova.csv");
   }
 
}else if(method=="anosim"){ # just one group
    anosim <- anosim(data.dist, group=group);
    resTab$Rval <- anosim$statistic;
    resTab$pval <- anosim$signif;
    stat.info <- paste("[ANOSIM] R: ", signif(resTab$Rval, 5), "; p-value < ", signif(resTab$pval, 5), sep="");
    stat.info.vec <- c(signif(resTab$Rval, 5), signif(resTab$pval, 5));
    names(stat.info.vec) <- c("R", "p-value");

  }else if (method=="permdisp"){ # just one group
    beta <- betadisper(data.dist, group=group);
    resTab <- anova(beta);
    stat.info <- paste("[PERMDISP] F-value: ", signif(resTab$"F value"[1], 5), "; p-value: ", signif(resTab$"Pr(>F)"[1], 5), sep="");
    stat.info.vec <- c(signif(resTab$"F value"[1], 5), signif(resTab$"Pr(>F)"[1], 5));
    names(stat.info.vec) <- c("F-value", "p-value");
  }else if(method=="mirkat"){
    
    if(!exists("MiRKAT")){ # public web on same user dir
      .load.scripts.on.demand("utils_mirkat.Rc"); 
    }
    exclud <- which(colnames(data@sam_data) %in% variable)
    X = data@sam_data[,-exclud]
    X = apply(X, 2, function(x) as.numeric(factor(x)))
    K.weighted <- D2K(as.matrix(data.dist))
    Gp <- as.numeric(factor(group))
    res<- MiRKAT(y =Gp, X = X, Ks = K.weighted, out_type = "C", method = "davies", returnKRV = TRUE, returnR2 = TRUE)
    resTab <- do.call(cbind,res);
    stat.info <- paste("[MiRKAT] R-squared: ", signif(resTab[1,"R2"], 5), "; p-value: ", signif(resTab[1,"p_values"], 5),"; KRV: ",signif(resTab[1,"KRV"], 5), sep="");
    stat.info.vec <- signif(resTab[1, c("R2", "p_values", "KRV")], 5);
    names(stat.info.vec) <- c("R-squared", "p-value", "KRV");
  }

  mbSetObj$analSet$stat.info <- stat.info;
  mbSetObj$analSet$stat.info.vec <- stat.info.vec;
  return(.set.mbSetObj(mbSetObj));
}

################################################
###### generate figure barplot with group#######
################################################
#'Function to plot group-wise bar charts.
#'@description This functions plots group-wise bar charts of a specified 
#'taxa for alpha diversity analysis.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import reshape
PlotTaxaAbundanceBarSamGrp<-function(mbSetObj, barplotName, taxalvl, metadata, facet2, imgOpt,
                                     feat_cnt, colpalopt, calcmeth, toptaxa,abunTopTaxaOpt, 
                                     appendnm, format="png", dpi=80, interactive = FALSE){
  save.image("taxgrp.RData")
  load_phyloseq();
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_reshape();
    load_viridis();
  }
  
  if(mbSetObj$module.type == "meta"){
    data1 <- qs::qread("merged.data.raw.qs");
    sample_table <- sample_data(data1);
    data <- as.data.frame(t(otu_table(data1)),check.names=FALSE);
    facet2="dataset";
  }else{
    #using filtered data
    data <- mbSetObj$dataSet$filt.data;
    if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
      data<-otu_table(data,taxa_are_rows =TRUE);
    }
  
    sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
    data1 <- merge_phyloseq(data, tax_table(mbSetObj$dataSet$proc.phyobj), sample_table);
  }

  yLbl <- "Actual Abundance";
  data <- as.data.frame(t(otu_table(data1)),check.names=FALSE);
  flag <- (facet2 == "null" || facet2 == "none" || facet2 == metadata);
  
  if(flag){
    data <- cbind(data, variable=data1@sam_data[[metadata]]);
    clsLbl <- unique(factor(data1@sam_data[[metadata]]));
  } else {
    data$variable <- paste(data1@sam_data[[metadata]], data1@sam_data[[facet2]], sep = ";");
    clsLbl <- unique(factor(paste(data1@sam_data[[metadata]],data1@sam_data[[facet2]], sep = ";")));
  }

  data <- aggregate(. ~variable,data,sum);
  rownames(data) <- data[,1];
  data <- data[,-1];
  data <- data[rownames(data),];
  
  data_tax <- tax_table(data1);
  #reshaping data
  if(taxalvl=="OTU"){
    taxa_nm <- as.matrix(colnames(data));
    rownames(taxa_nm) <- colnames(data);
    rownames(taxa_nm) <- sub("^X", "", rownames(taxa_nm))
  }else{
    taxa_nm <- as.character(data_tax[,taxalvl]);
    y <- which(is.na(taxa_nm)==TRUE);
    #converting NA values to unassigned
    taxa_nm[y] <- "Not_Assigned";
    
    if(appendnm=="T"){
      all_nm <- colnames(tax_table(data1));
      hg_nmindx <- which(all_nm==taxalvl)-1;
      
      if(hg_nmindx!=0){
        nma <- as.character(tax_table(data1)[,hg_nmindx]);
        y1 <- which(is.na(nma)==TRUE);
        nma[y1] <- "Not_Assigned";
        nm <- paste0(nma,"_",taxa_nm);
        ind <- which(nm=="Not_Assigned_Not_Assigned");
        nm[ind] <- "Not_Assigned";
        nm <- gsub("_Not_Assigned", "",nm, perl = TRUE);
        taxa_nm <- nm;
      }
    }
  }
  
  taxa_nm <- as.matrix(taxa_nm);
  
  if(appendnm=="T"){
    y <- which(grepl("Not_Assigned", taxa_nm))
  }else{
    y <- which(is.na(taxa_nm)==TRUE);
  }
  
  #converting NA values to unassigned; before order it to last position using ZZZ as its name
  taxa_nm[y] <- "ZZZ";
  colnames(data) <- taxa_nm[,1];
  nms <- colnames(data);
  data <- as.matrix(data);
  data <- data %*% sapply(unique(nms),"==",nms);
  data <- data.frame(data,check.names=FALSE);
  data <- data[ , order(names(data))];
  indx <- which(colnames(data)=="ZZZ");
  colnames(data)[indx] <- "NA";
  if(abunTopTaxaOpt == "top"){
    if(calcmeth=="sum"){
      feat_cnt <- sort(sapply(data, sum), decreasing = TRUE)[toptaxa];
    }else {
      feat_cnt <- sort(sapply(data, median), decreasing = TRUE)[toptaxa];
    }
  }
  if(calcmeth=="sum"){
    ind <- which(colSums(data)>feat_cnt);
    ind1 <- which(colSums(data)<feat_cnt);
  } else {
    dt <- apply(data,2,median);
    ind <- which(dt>feat_cnt);
    ind1 <- which(dt<feat_cnt);
  }
  
  if(length(ind)==0){
    AddErrMsg("All features have lower read count than given minimum count filter. Please lower the cut off for minimum count.");
    return(0);
  }
  
  if(length(ind1)>0){
    colnames(data)[ind1] <- "Others";
    data <- as.data.frame(do.call(cbind, by(t(data),INDICES=names(data),FUN=colSums)),check.names=FALSE);
  }
  
  if(imgOpt=="barnorm"){
    data <- as.data.frame(apply(data,1, function(x) x/sum(x)),check.names=FALSE);
    data <- as.data.frame(t(data),check.names=FALSE);
    yLbl <- "Relative Abundance";
  }
  
  w <- 650;
  #adjust height according to number of taxa 
  feat_no <- ncol(data); 
  if(feat_no < 10){
    h<-w;
  } else if (feat_no < 20){
    w<-w+100;
    h<-w/2+100;
  } else if (feat_no < 50){
    w<-w+200;
    h<-w/2+150;
  } else if (feat_no < 100){
    w<-w+350
    h<-w/2+200;
  } else if (feat_no > 100){
    w<-w+500
    h<-w/2+300;
  }
  
  if(flag){
  } else {
    feat_no2 <- nlevels(as.factor(data1@sam_data[[facet2]]));
    if(feat_no2 < 3){
      w <- w;
    } else if (feat_no2 < 4){
      w <- w * 1.1
    } else if (feat_no2 < 6){
      w <- w * 1.2;
    } else if (feat_no2 < 8){
      w <- w * 1.3;
    } else  if (feat_no2 < 10) {
      w <- w * 1.4;
    } else  if (feat_no2 < 12) {
      w <- w * 1.5;
    } else  if (feat_no2 < 14) {
      w <- w * 1.6;
    } else {
      w <- w * 2;
    }
  }
  # final check
  if(w > 1800){
    w <- 1800;
  }
  if(h < w/3){
    h <- w/3
  }
  if(h < 480){
    h <- 480;
  }
  
  fast.write(t(data), file="taxa_abund.csv");
  data$step <- factor(rownames(data), levels = rownames(data));
  gp_nm <- data$step
  data <- melt(data,id='step');
  data$step <- as.numeric(data$step);
  data <- data[order(data[,2]),];
  data <- data[,-1];
  data$step <- rep(gp_nm,feat_no);
  
  tmp_df <- aggregate(data$value, by=list(data$variable), FUN=mean)
  var_level <- tmp_df[order(tmp_df$x, decreasing = TRUE), ][[1]]
  #as.character()# get the order of taxa relative abundance
  data$variable <- factor(data$variable, levels = var_level) # change the factor level of taxa
  levels(data$variable) <- sub("^X", "", levels(data$variable)) 
  
  if(flag){
  }else {
    data$variable2 <- sapply(strsplit(as.character(data$step), ";", fixed = TRUE), "[", 2);
    data$step <- sapply(strsplit(as.character(data$step), ";", fixed = TRUE), "[", 1);
  }
  
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
  
  if(x <= 10){
    guide_num = 3
  }else if (x <= 20){
    guide_num = 4
  }else{
    guide_num = 5
  }

  rdaName = paste(barplotName, ".rda", sep="");
  jsonName = paste(barplotName, ".json", sep="");
  barplotName = paste(barplotName, ".",format, sep="");
  mbSetObj$imgSet$stack<-barplotName;
  mbSetObj$imgSet$stackRda <-rdaName;
  mbSetObj$imgSet$stackType <- "group";

  stackdata <<- data
  Cairo::Cairo(file=barplotName,width=w, height=h, type=format, bg="white",dpi=dpi);
  box <- ggplot(data,aes(x = step, y = value, fill = variable))+
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
  
  if(flag){
    box <- box;
  } else {
    box <- box + facet_grid(variable2 ~ ., scales = "free");
  }
  
  if(colpalopt %in% c("set1", "set2", "set3")){
    
    cols.needed <- length(unique(data$variable))
    colpalopt <- tools::toTitleCase(tolower(colpalopt))
    
    if(colpalopt == "Set1"){
      n.cols <- 9
    }else if(colpalopt == "Set2"){
      n.cols <- 8
    }else{
      n.cols <- 12
    }
    
    if(cols.needed > n.cols){
      col.func <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(n.cols, colpalopt))
      box <- box + scale_fill_manual(values=col.func(cols.needed),
                                     guide = guide_legend(direction = "horizontal",
                                                          ncol = guide_num))
    } else {
      box <- box + scale_fill_brewer(palette = colpalopt,
                                     guide = guide_legend(direction = "horizontal",
                                                          ncol = guide_num))
    }
  }else if(colpalopt %in% c("magma","plasma","inferno","viridis")){
    box <- box + viridis::scale_fill_viridis(option=colpalopt, discrete=TRUE) + guides(fill=guide_legend(ncol=guide_num))
  }else{
    box <- box + scale_fill_manual(values=c(x.colors)) + guides(fill=guide_legend(ncol=guide_num))
  }

  if(mbSetObj$module.type == "meta"){
      box <- box + facet_grid(variable2 ~ . , scales = "free", space = "free");
  }

  save(box,file=rdaName);
  print(box+guides(fill=guide_legend(ncol=3)));
  dev.off();
  
 

  mbSetObj$analSet$stack<-data;
  mbSetObj$analSet$stack.taxalvl<-taxalvl;
  mbSetObj$analSet$plot<-"Stacked Bar";
  return(.set.mbSetObj(mbSetObj));
}

#'Function to create rarefraction curves of microbiome data
#'@description This functions plots rarefraction curves.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param data.src Character, use "orig" to use the original data or
#'"proc" to use the processed data.
#'@param linecolor Character, input the metadata you wish to 
#'group the sample line colors.
#'@param linetype Character, input the metadata you wish to 
#'group the sample line types.
#'@param facet Character, input the metadata you wish to 
#'group the samples.
#'@param step Numeric, input the step number. Default is set to 5.
#'@param imgName Character, input the name of the plot to be saved.
#'@param format Character, input the type of plot to be saved. Default is 
#'set to png. 
#'@param dpi Numeric, input the dpi for the plot to be saved. Default
#'is set to 72.
#'@param interactive Boolean, if set to TRUE, saves the plot
#'as an interactive html plot.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotRarefactionCurve <- function(mbSetObj, data.src, linecolor, linetype, facet, step=5, 
                                 imgName, format="png", dpi=72, interactive = FALSE){

  load_phyloseq();
 current.msg<<-"";
  mbSetObj <- .get.mbSetObj(mbSetObj);
  # should use unfiltered data
  if(data.src == "orig"){
    data_rare <- qs::qread("orig.phyobj");
  }else{
    data_rare <- mbSetObj$dataSet$proc.phyobj;
  }
  
  sample_table_msg <- sample_data(data_rare, errorIfNULL=TRUE);
  
  if(min(table(factor(sample_table_msg [[linecolor]]))) < 3 | min(table(factor(sample_table_msg [[linetype]]))) < 3 | min(table(factor(sample_table_msg [[facet]]))) < 3){
    AddErrMsg("Too many groups to be displayed - please select a more meaningful group option with at least 3 samples per group.");
    return(0);
  }
  
  #get good's coverage index
  goods_coverage <- ComputeGoods(data_rare)
 
  fast.write(goods_coverage, "goods_coverage.csv", row.names = FALSE, quote = FALSE);
  
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
  
  mbSetObj$analSet$rarefaction_curve <- imgName;
  mbSetObj$analSet$rarefaction_curve_data.src = ifelse(data.src == "orig", "original", "filtered");
  
  if(interactive){
    library(plotly)
    library(htmlwidgets)
    ax <- list(
      zeroline = FALSE,
      showline = FALSE,
      showgrid = FALSE
    )
    p <- plotly::ggplotly(box)
    p <- p %>% layout(xaxis = ax, yaxis = ax)
    htmlwidgets::saveWidget(p, "rarefaction_interactive.html")
  }
  
  return(.set.mbSetObj(mbSetObj))
}

#'Function to prepare data for phylogenetic tree.
#'@description This functions prepares the data to plot the phylogenetic tree.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotPhylogeneticTree <-function(mbSetObj, color, shape, taxa, treeshape, imgName, format="png", dpi=72){

  load_phyloseq();
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  data <- mbSetObj$dataSet$filt.data;
  
  if("matrix" %in% class(mbSetObj$dataSet$filt.data)){
    data <- otu_table(data,taxa_are_rows =TRUE);
  }
  
  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL=TRUE);
  
  if(min(table(factor(sample_table [[color]]))) < 3 | min(table(factor(sample_table [[shape]]))) < 3){
    AddErrMsg("Too many groups to be displayed - please select a more meaningful group option with at least 3 samples per group.");
    return(0);
  }
  
  data1 <- merge_phyloseq(data,tax_table(mbSetObj$dataSet$proc.phyobj), sample_table);
  pg_tree <- qs::qread("tree.qs");
  pg_tb <- tax_table(data1);
  pg_ot <- otu_table(data1);
  pg_sd <- sample_data(data1);
  pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
  data_tree <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree);
  for (i in sample_variables(data_tree)){
    if(is.integer(sample_data(data_tree)[[i]]) | is.character(sample_data(data_tree)[[i]])){
      sample_data(data_tree)[[i]] <- factor(sample_data(data_tree)[[i]])
    }
  }
  
  data_tree <- fast_tax_glom_mem(data_tree, taxa);
  if(is.null(data_tree)){
    AddErrMsg("Errors in projecting to the selected taxanomy level!");
    return(0);
  } 
  
  feat_no <- length(phy_tree(data_tree)$tip.label);
  if(feat_no > 50){
    AddErrMsg("There are too many tree tips to be display, please select a higher taxonomy level");
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
  box <- plot_tree(data_tree,
                   ladderize = "left",
                   color = color,
                   shape = shape,
                   label.tips = taxa,
                   size = "abundance");
  
  if(treeshape == "rectangular"){
    Cairo(file=imgName, width = 1.5 * w, height = w, type=format, bg="white", dpi=dpi);
  }else{
    box <- box + coord_polar(theta="y");
    Cairo(file=imgName, width = 1.5*w, height = 1.5*w, type=format, bg="white", dpi=dpi);
  }
  print(box);
  dev.off();
  
  mbSetObj$analSet$phylogenetic_tree_curve <- imgName;
  mbSetObj$analSet$phylogenetic_tree_curve_tax_level <- taxa;
  
  return(.set.mbSetObj(mbSetObj))
}

#########################
#####plot heat tree######
#########################

GetHtGroupItems <- function(mbSetObj, meta){
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  dm <- mbSetObj$dataSet$proc.phyobj;  
  dm_samples = as(sample_data(dm), "data.frame");
  grp.nms <- as.character(unique(dm_samples[[meta]]));
  return(grp.nms);
}

SetGroupItems <- function(mbSetObj, groups){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  num.groups <- length(groups)
  
  if(num.groups < 2){
    AddErrMsg("Minimum 2 groups for comparison!")
    return(0)
  }
  
  mbSetObj$corr.net$comparison <- groups
  return(.set.mbSetObj(mbSetObj))
}

GetHtMetaCpInfo <- function(mbSetObj, meta){
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  dm <- mbSetObj$dataSet$proc.phyobj;  
  dm_samples = as(sample_data(dm), "data.frame");
  
  grp.nms <- as.character(unique(dm_samples[[meta]]));
  myargs <- list();
  inx = 0
  for (m in 1:(length(grp.nms) - 1)) {
    for (n in (m + 1):length(grp.nms)) {
      inx <- inx + 1
      myargs[[inx]] <- paste(grp.nms[m], "_vs_", grp.nms[n], sep = "");
    }
  }
  return(unlist(myargs));
}

#'Function to plot heat tree
#'@description This functions plots the heat tree.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param meta Input the name of the group.
#'@param taxalvl Character, input the taxonomic level to perform
#'classification. For instance, "Genus" to use the Genus level.
#'@param color Character, input the color palette code. "dbgr"
#'for dark blue, grey and red. "bgy" for greenblue, grey and yellow.
#'"ggr" for green, grey and red. "pgy" for purple, gray and yellow. 
#'"tgr" for teal, grey and red. "ggg" for green, grey and gold. Additionally,
#'the viridis R package can be used to generate color schemes. "plasma" for 
#'the plasma color scheme, "viridis" for the default viridis color scheme and
#'"cividis" for the cividis color scheme.
#'@param layoutOpt Character, input the layout of the heat tree. "dft"
#'for the default layout and "reda" for reingold-tilford.
#'@param comparison Character, input the group comparisons.
#'@param wilcox.cutoff Numeric, input the Wilcoxan p-value cutoff
#'for significant node labels.
#'@param imgName Character, input the name of the heat tree plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import metacoder
PrepareHeatTreePlot <- function(mbSetObj, meta, taxalvl, color, layoutOpt, comparison, 
                                wilcox.cutoff, switchCmpDirection, colorMode, showLabels, imgName, format="png", dpi=72){
  load_metacoder();
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);  
  tax_o <- taxalvl;
  dm <- mbSetObj$dataSet$proc.phyobj;
  tax_table_new = data.frame("Kingdom" = "Root", as(tax_table(dm), "matrix")[, 1:ncol(tax_table(dm))],check.names=FALSE) # add root to tax table
  
  tax_table(dm) <- as.matrix(tax_table_new);
  
  dm_samples = as(sample_data(dm), "data.frame");
  dm_samples <- cbind.data.frame("sample_id" = row.names(dm_samples), dm_samples); # add sample_id to sample table
  row.names(dm_samples) <- c();
  dm_samples[] <- lapply(dm_samples, as.character);
  
  otu_dm <- as.data.frame(as(otu_table(dm), "matrix"),check.names=FALSE);
  tax_dm <- as.data.frame(as(tax_table(dm), "matrix"),check.names=FALSE);
  tax_dm[] <- lapply(tax_dm, as.character);#make sure characters in tax_dm;
  depth <- ncol(tax_dm)
  rank_dm <- c("r", "p", "c", "o", "f", "g", "s");
  names(tax_dm) <- rank_dm[1:depth];
  
  for (i in 1:ncol(tax_dm)){
    for (j in 1:nrow(tax_dm)){
      if (is.na(tax_dm[j, i])){
        tax_dm[j, i] <- "";
      } else {
        tax_dm[j, i] <- paste(names(tax_dm)[i],
                              tax_dm[j, i],
                              sep = "__");
      }
    }
  } #add __ to tax table
  
  if(taxalvl == "Phylum"){
    tax <- "p";
  } else if(taxalvl == "Class"){
    tax <- "c";
  } else if(taxalvl == "Order"){
    tax <- "o";
  } else if(taxalvl == "Family"){
    tax <- "f";
  } else if(taxalvl == "Genus"){
    tax <- "g"
  } else {
    tax <- "s";
  }; # get tax rank for heat tree
  
  tax_dm <- tax_dm[, 1:which(names(tax_dm) == tax)]; #subset tax table
  rank_dm_new <- rank_dm[1:which(rank_dm == tax)];
  tax_dm$lineage <- apply(tax_dm[, rank_dm_new], 1, paste, collapse = ";"); #collapse all tax ranks
  dm_otu <- cbind.data.frame("otu_id" = row.names(tax_dm),
                             "lineage" = tax_dm$lineage,
                             otu_dm); #make new otu table
  row.names(dm_otu) <- c();
  dm_otu$lineage <- gsub(";+$", "", dm_otu$lineage); #remove empty tax names
  
  grp.nms <- strsplit(comparison, "_vs_")[[1]];
  #dm_samples_cmf <- dm_samples[dm_samples[[meta]] %in% grp.nms, ]; #subset sample data by meta variable  
  grpDF1 <- dm_samples[dm_samples[[meta]] == grp.nms[1], ];
  grpDF2 <- dm_samples[dm_samples[[meta]] == grp.nms[2], ];
  if(switchCmpDirection == "false"){
    dm_samples_cmf <- rbind.data.frame(grpDF1, grpDF2);
  } else {
    dm_samples_cmf <- rbind.data.frame(grpDF2, grpDF1);
    comparison <- paste0(grp.nms[2], "_vs_", grp.nms[1]);
  }
  dm_otu_cmf <- dm_otu[, c("otu_id", "lineage", dm_samples_cmf$sample_id)]; #make otu table ready for heat tree  
  
  PrepareHeatTreePlotDataParse_cmf_res <- PrepareHeatTreePlotDataParse_cmf(dm_otu_cmf, dm_samples_cmf, meta); #parse otu table to heat tree object
  PrepareHeatTreePlotDataParse_cmf_diff_table_res <- PrepareHeatTreePlotDataParse_cmf_diff_table(PrepareHeatTreePlotDataParse_cmf_res); #generate diff table
  PrepareHeatTreePlotDataParse_cmf_res <<- PrepareHeatTreePlotDataParse_cmf_res; #generate heat tree
  
  PrepareHeatTreePlotDataParse_cmf_plot(mbSetObj, color, layoutOpt, comparison, wilcox.cutoff,colorMode, showLabels, imgName, format, dpi=72);
  
  #below is for PDF reporter
  mbSetObj$analSet$heat_tree_plot <- paste0(imgName,".",format); 
  mbSetObj$analSet$heat_tree_meta <- meta;
  mbSetObj$analSet$heat_tree_tax <- tax_o;
  mbSetObj$analSet$heat_tree_comparison <- comparison;
  current.msg<<-"Heat tree analysis successful!";
  return(.set.mbSetObj(mbSetObj));
};

#'Function to prepare heat tree data
#'@description This functions plots the heat tree.
#'@param meta Input the name of the group.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import metacoder
#'@import taxa
PrepareHeatTreePlotDataParse_cmf <- function(dm_otu_cmf,
                                             dm_samples_cmf,
                                             meta){
  dm_obj_cmf <- metacoder::parse_tax_data(dm_otu_cmf,
                                          class_cols = "lineage",
                                          class_sep = ";",
                                          class_regex = "^(.+)__(.+)$",
                                          class_key = c(tax_rank = "info",
                                                        tax_name = "taxon_name"));
  
  dm_obj_cmf$data$tax_data <- calc_obs_props(dm_obj_cmf, "tax_data"); # normalization
  
  dm_obj_cmf$data$tax_abund <- calc_taxon_abund(dm_obj_cmf, "tax_data",
                                                cols = dm_samples_cmf$sample_id); #calculate taxon abundance
  dm_obj_cmf$data$diff_table <- compare_groups(dm_obj_cmf, dataset = "tax_abund", #make diff table
                                               cols = dm_samples_cmf$sample_id, # What columns of sample data to use
                                               groups = dm_samples_cmf[[meta]]); # What category each sample is assigned to;
  return(dm_obj_cmf);
};

#'Function to prepare heat tree data
#'@description This functions plots the heat tree.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import metacoder
PrepareHeatTreePlotDataParse_cmf_diff_table <- function(PrepareHeatTreePlotDataParse_cmf_res){#generate diff table for downloading
  
  dm_obj_cmf = PrepareHeatTreePlotDataParse_cmf_res;
  table_tax_dm <- dm_obj_cmf$data$class_data[-2];
  table_tax_dm <- table_tax_dm[!duplicated(table_tax_dm$taxon_id), ]; #remove redundancy
  
  tax_diff_dm <- merge(dm_obj_cmf$data$diff_table, 
                       table_tax_dm,
                       by = "taxon_id"); #add tax names
  tax_diff_dm <- cbind.data.frame("id" = tax_diff_dm[, 1],
                                  tax_diff_dm[, 8:9],
                                  tax_diff_dm[, 2:7]);
  
  fast.write(tax_diff_dm, 
             "tax_diff_dm.csv", 
             row.names = FALSE, quote = FALSE);
  
  return(tax_diff_dm);
};

#'Function to prepare heat tree data
#'@description This functions plots the heat tree.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param color Character, input the color palette code. "dbgr"
#'for dark blue, grey and red. "bgy" for greenblue, grey and yellow.
#'"ggr" for green, grey and red. "pgy" for purple, gray and yellow. 
#'"tgr" for teal, grey and red. "ggg" for green, grey and gold. Additionally,
#'the viridis R package can be used to generate color schemes. "plasma" for 
#'the plasma color scheme, "viridis" for the default viridis color scheme and
#'"cividis" for the cividis color scheme.
#'@param layoutOpt Character, input the layout of the heat tree. "dft"
#'for the default layout and "reda" for reingold-tilford.
#'@param comparison Character, input the group comparisons.
#'@param wilcox.cutoff Numeric, input the Wilcoxan p-value cutoff
#'for significant node labels.
#'@param imgName Character, input the name of the heat tree plot.
#'@param format Character, input the preferred
#'format of the plot. By default it is set to "png".
#'@param dpi Numeric, input the dots per inch. By default
#'it is set to 72.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import metacoder
#'@import viridis
PrepareHeatTreePlotDataParse_cmf_plot <- function(mbSetObj, color, layoutOpt, comparison, wilcox.cutoff, colorMode, showLabels, imgName, format, dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_viridis();
  
  dm_obj_cmf = PrepareHeatTreePlotDataParse_cmf_res;
  if(color == "ggr"){
    color_new <- c("#006B30", "gray", "#E21818");
  } else if(color == "dbgr") {
    color_new <- c("#2A4196", "gray", "#F71735");
  } else if(color == "tgr"){
    color_new <- c("#007777", "gray", "#F71735");
  }else if(color == "bgy"){
    color_new <- c("#2EAA9C", "gray", "#FCB932");
  }else if(color == "ggg"){
    color_new <- c("#007777", "gray", "#DAA520");
  }else if(color == "plasma") {
    color_new <- rev(c(viridis::plasma(10)[4:9], "lightgrey"))
  }else if(color == "viridis") {
    color_new <- rev(c(viridis::viridis(10)[4:9], "lightgrey"))
  }else if(color == "cividis") {
    color_new <- rev(c(viridis::cividis(10)[4:9], "lightgrey"))
  }else {
    color_new <- c("#764b93", "gray", "#F0C808");
  };
  
  mbSetObj$analSet$heat_tree_plot <- imgName; # for PDF reporter
  
  set.seed(56784);
  
  Cairo::Cairo(file=paste0(imgName, ".", format), height = 875, width = 1000, type=format, bg="white", dpi=96);
  wilcox.cutoff <- as.numeric(wilcox.cutoff)

   if(colorMode=="sig"){
    dm_obj_cmf$data$diff_table$log2_median_ratio[dm_obj_cmf$data$diff_table$wilcox_p_value>wilcox.cutoff]  <- 0
   }

if(showLabels=="true"){
  if(layoutOpt == "reda"){# two layouts are provided
    box <- heat_tree(dm_obj_cmf,
                     node_label =  ifelse(wilcox_p_value < wilcox.cutoff, taxon_names, NA),  #taxon names
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                     node_color_interval = c(-8, 8), # The range of `log2_median_ratio` to display
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "Abundance level",
                     node_color_axis_label = "Median ratio (log2)",
                     layout = "davidson-harel", # The primary layout algorithm
                     initial_layout = "reingold-tilford",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);
  } else {
    box <- heat_tree(dm_obj_cmf,
                     node_label =  ifelse(wilcox_p_value < wilcox.cutoff, taxon_names, NA),
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                     node_color_interval = c(-8, 8), # The range of `log2_median_ratio` to display
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "Abundance level",
                     node_color_axis_label = "Median ratio (log2)",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);

  }
}else{
if(layoutOpt == "reda"){# two layouts are provided
    box <- heat_tree(dm_obj_cmf,
                     node_label =  NA,  #taxon names
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                     node_color_interval = c(-8, 8), # The range of `log2_median_ratio` to display
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "Abundance level",
                     node_color_axis_label = "Median ratio (log2)",
                     layout = "davidson-harel", # The primary layout algorithm
                     initial_layout = "reingold-tilford",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);
  } else {
    box <- heat_tree(dm_obj_cmf,
                     node_label =  NA,
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = log2_median_ratio, # A column from `obj$data$diff_table`
                     node_color_interval = c(-8, 8), # The range of `log2_median_ratio` to display
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "Abundance level",
                     node_color_axis_label = "Median ratio (log2)",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);

  }

};
  print(box);
  dev.off();
  
  return(.set.mbSetObj(mbSetObj))
}

PlotGroupDataHeattree <- function(mbSetObj, meta, comparison, taxalvl, color, layoutOpt, 
                                  showLabels,imgName, format="png", dpi=72){
  load_metacoder();
  load_phyloseq();
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  tax_o <- taxalvl;
  
  dm <- mbSetObj$dataSet$proc.phyobj;
  tax_table_new = data.frame("Kingdom" = "Root", as(tax_table(dm), "matrix")[, 1:6],check.names=FALSE) # add root to tax table
  tax_table(dm) <- as.matrix(tax_table_new);
  
  dm_samples = as(sample_data(dm), "data.frame");
  dm_samples <- cbind.data.frame("sample_id" = row.names(dm_samples), dm_samples); # add sample_id to sample table
  row.names(dm_samples) <- c();
  dm_samples[] = lapply(dm_samples, FUN = function(x){gsub(" ", "", x)});
  comparison <- gsub(" ", "", comparison);
  if(sum(grepl("^[0-9]{1,}$", dm_samples[[meta]])) == nrow(dm_samples)){
    dm_samples[[meta]] <- paste(meta, dm_samples[[meta]], sep = "_");
    comparison <- paste(meta, comparison, sep = "_");
  };
  dm_samples_cmf <- dm_samples[dm_samples[[meta]] %in% comparison, ]; #subset sample data by meta variable
  dm_samples_cmf$meta_com <- comparison;
  
  flag <- FALSE;
  
  PrepareHeatTreePlotAbR(dm, tax_dm, taxalvl, dm_samples_cmf, color, showLabels,imgName, format, layoutOpt, comparison, flag);
  
  #below is for PDF reporter
  mbSetObj$analSet$heat_tree_plot <- imgName; 
  #mbSetObj$analSet$heat_tree_plot <- paste(imgName,".", format, sep="")
  mbSetObj$analSet$heat_tree_meta <- meta;
  mbSetObj$analSet$heat_tree_tax <- tax_o;
  mbSetObj$analSet$heat_tree_comparison <- comparison;
  current.msg<<-"Heat tree analysis successful!";
  return(.set.mbSetObj(mbSetObj));
};

PlotSampleDataHeattree <- function(mbSetObj,comparison, taxalvl, color, layoutOpt, 
                                  showLabels, imgName, format="png", dpi=72){
  load_metacoder();
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  tax_o <- taxalvl;
  
  dm <- mbSetObj$dataSet$proc.phyobj;
  tax_table_new = data.frame("Kingdom" = "Root", as(tax_table(dm), "matrix")[, 1:6],check.names=FALSE) # add root to tax table
  tax_table(dm) <- as.matrix(tax_table_new);
  
  dm_samples = as(sample_data(dm), "data.frame");
  dm_samples <- cbind.data.frame("sample_id" = row.names(dm_samples), dm_samples); # add sample_id to sample table
  row.names(dm_samples) <- c();
  dm_samples[] <- lapply(dm_samples, as.character);
  
  if(sum(grepl("^[0-9]{1,}$", dm_samples$sample_id)) == nrow(dm_samples)){
    dm_samples$sample_id <- paste("sample", dm_samples$sample_id, sep = "_");
    comparison <- paste("sample", comparison, sep = "_"); 
    flag <- TRUE;
  } else {
    flag <- FALSE; 
  };  
  
  dm_samples_cmf <- dm_samples[dm_samples$sample_id %in% comparison, ];
  dm_samples_cmf$meta_com <- comparison;
  
  PrepareHeatTreePlotAbR(dm, tax_dm, taxalvl, dm_samples_cmf, color,showLabels, imgName, format, layoutOpt, comparison, flag);
  
  #below is for PDF reporter
  mbSetObj$analSet$heat_tree_plot <- imgName; 
  mbSetObj$analSet$heat_tree_meta <- comparison;
  mbSetObj$analSet$heat_tree_tax <- tax_o;
  mbSetObj$analSet$heat_tree_comparison <- comparison;
  current.msg<<-"Heat tree analysis successful!";
  return(.set.mbSetObj(mbSetObj));
};

PlotOverviewDataHeattree <- function(mbSetObj, taxalvl, color, layoutOpt, 
                                    showLabels, imgName, format="png", dpi=72){
  load_metacoder();
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  tax_o <- taxalvl;
  
  dm <- mbSetObj$dataSet$proc.phyobj;
  tax_table_new = data.frame("Kingdom" = "Root", as(tax_table(dm), "matrix")[, 1:6],check.names=FALSE) # add root to tax table
  tax_table(dm) <- as.matrix(tax_table_new);
  
  dm_samples = as(sample_data(dm), "data.frame");
  dm_samples <- cbind.data.frame("sample_id" = row.names(dm_samples), dm_samples); # add sample_id to sample table
  row.names(dm_samples) <- c();
  dm_samples[] <- lapply(dm_samples, as.character);
  
  dm_samples_cmf <- dm_samples;
  dm_samples_cmf$meta_com <- "All_samples";
  comparison <- "All_samples";
  
  flag <- FALSE;
  
  PrepareHeatTreePlotAbR(dm, tax_dm, taxalvl, dm_samples_cmf, color,showLabels, imgName, format, layoutOpt, comparison, flag);  
  
  #below is for PDF reporter
  mbSetObj$analSet$heat_tree_plot <- imgName; 
  mbSetObj$analSet$heat_tree_meta <- "All samples";
  mbSetObj$analSet$heat_tree_tax <- tax_o;
  mbSetObj$analSet$heat_tree_comparison <- "All samples";
  current.msg<<-"Heat tree analysis successful!";
  return(.set.mbSetObj(mbSetObj));
};

PrepareHeatTreePlotAbR <- function(dm = dm, tax_dm = tax_dm, taxalvl = taxalvl, dm_samples_cmf = dm_samples_cmf, 
                                   color = color,showLabels, imgName = imgName, format = format, layoutOpt= layoutOpt, comparison = comparison, flag = flag){
  
  otu_dm <- as.data.frame(as(otu_table(dm), "matrix"),check.names=FALSE);
  tax_dm <- as.data.frame(as(tax_table(dm), "matrix"),check.names=FALSE);
  tax_dm[] <- lapply(tax_dm, as.character);#make sure characters in tax_dm;
  depth <- ncol(tax_dm)
  rank_dm <- c("r", "p", "c", "o", "f", "g", "s");
  names(tax_dm) <- rank_dm[1:depth];
  
  for (i in 1:ncol(tax_dm)){
    for (j in 1:nrow(tax_dm)){
      if (is.na(tax_dm[j, i])){
        tax_dm[j, i] <- "";
      } else {
        tax_dm[j, i] <- paste(names(tax_dm)[i],
                              tax_dm[j, i],
                              sep = "__");
      }
    }
  } #add __ to tax table
  
  if(taxalvl == "Phylum"){
    tax <- "p";
  } else if(taxalvl == "Class"){
    tax <- "c";
  } else if(taxalvl == "Order"){
    tax <- "o";
  } else if(taxalvl == "Family"){
    tax <- "f";
  } else if(taxalvl == "Genus"){
    tax <- "g"
  } else {
    tax <- "s";
  }; # get tax rank for heat tree
  
  tax_dm <- tax_dm[, 1:which(names(tax_dm) == tax)]; #subset tax table
  rank_dm_new <- rank_dm[1:which(rank_dm == tax)];
  tax_dm$lineage <- apply(tax_dm[, rank_dm_new], 1, paste, collapse = ";"); #collapse all tax ranks
  dm_otu <- cbind.data.frame("otu_id" = row.names(tax_dm),
                             "lineage" = tax_dm$lineage,
                             otu_dm); #make new otu table
  row.names(dm_otu) <- c();
  if(flag){
    names(dm_otu)[-c(1:2)] <- paste("sample", names(dm_otu)[-c(1:2)], sep = "_");
  };
  dm_otu$lineage <- gsub(";+$", "", dm_otu$lineage); #remove empty tax names
  
  dm_otu_cmf <- dm_otu[, c("otu_id", "lineage", dm_samples_cmf$sample_id)]; #make otu table ready for heat tree  
  dm_obj_cmf <- metacoder::parse_tax_data(dm_otu_cmf,
                                          class_cols = "lineage",
                                          class_sep = ";",
                                          class_regex = "^(.+)__(.+)$",
                                          class_key = c(tax_rank = "info",
                                                        tax_name = "taxon_name"));
  
  dm_obj_cmf$data$tax_data <- calc_obs_props(dm_obj_cmf, "tax_data"); # normalization
  
  dm_obj_cmf$data$tax_abund <- calc_taxon_abund(dm_obj_cmf, "tax_data",
                                                cols = dm_samples_cmf$sample_id); #calculate taxon abundance
  
  dm_obj_cmf$data$tax_occ <- calc_n_samples(dm_obj_cmf, "tax_abund", groups = dm_samples_cmf$meta_com)
  
  if(color == "ggr"){
    color_new <- c("#006B30", "gray", "#E21818");
  } else if(color == "dbgr") {
    color_new <- c("#2A4196", "gray", "#F71735");
  } else if(color == "tgr"){
    color_new <- c("#007777", "gray", "#F71735");
  }else if(color == "bgy"){
    color_new <- c("#2EAA9C", "gray", "#FCB932");
  }else if(color == "ggg"){
    color_new <- c("#007777", "gray", "#DAA520");
  }else if(color == "plasma") {
    color_new <- rev(c(viridis::plasma(10)[4:9], "lightgrey"))
  }else if(color == "viridis") {
    color_new <- rev(c(viridis::viridis(10)[4:9], "lightgrey"))
  }else if(color == "cividis") {
    color_new <- rev(c(viridis::cividis(10)[4:9], "lightgrey"))
  }else {
    color_new <- c("#764b93", "gray", "#F0C808");
  };
  
  
  set.seed(56784);
  
  Cairo::Cairo(file=paste0(imgName, ".", format), height = 875, width = 1000, type=format, bg="white", dpi=96);
  
  if(showLabels=="true"){
  if(layoutOpt == "reda"){# two layouts are provided
    box <- heat_tree(dm_obj_cmf,
                     node_label = taxon_names,
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = dm_obj_cmf$data$tax_occ[[2]],# A column from `obj$data$tax_occ`
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "OTU count",
                     node_color_axis_label = "Samples with reads",
                     layout = "davidson-harel", # The primary layout algorithm
                     initial_layout = "reingold-tilford",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);
  } else {
    box <- heat_tree(dm_obj_cmf,
                     node_label = taxon_names,
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = dm_obj_cmf$data$tax_occ[[2]],# A column from `obj$data$tax_occ`
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "OTU count",
                     node_color_axis_label = "Samples with reads",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);
  }
}else{
if(layoutOpt == "reda"){# two layouts are provided
    box <- heat_tree(dm_obj_cmf,
                     node_label = NA,
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = dm_obj_cmf$data$tax_occ[[2]],# A column from `obj$data$tax_occ`
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "OTU count",
                     node_color_axis_label = "Samples with reads",
                     layout = "davidson-harel", # The primary layout algorithm
                     initial_layout = "reingold-tilford",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);
  } else {
    box <- heat_tree(dm_obj_cmf,
                     node_label = NA,
                     node_size = n_obs, # n_obs is a function that calculates, in this case, the number of OTUs per taxon
                     node_color = dm_obj_cmf$data$tax_occ[[2]],# A column from `obj$data$tax_occ`
                     node_color_range = color_new, # The color palette used
                     node_size_axis_label = "OTU count",
                     node_color_axis_label = "Samples with reads",
                     title = comparison,
                     title_size = 0.05,
                     node_label_size_range = c(0.02, 0.05),
                     output_file = NULL);
  }

}
   ;
  print(box);
  dev.off();
}

#######################################
###########Tax4Fun/PICRUSt ############
#######################################

#'Main function to perform 16S functional annotation
#'@description This is the main function to perform either PICRUSt or
#'SILVA for functional annotation on the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import Tax4Fun

# use memoise tech to avoid repeating
Perform16FunAnot_mem <- function(mbSetObj, type, pipeline,ggversion){
  if(!exists("perform16funanot_mem")){
    require("memoise");
    perform16funanot_mem <<- memoise(Perform16FunAnot); 
  }
  
  res <- perform16funanot_mem(mbSetObj, type, pipeline,ggversion);
  return(res);
}


Perform16FunAnot<-function(mbSetObj, type, pipeline,ggversion) {
  if(!exists("my.16sfun.anot")){ # public web on same user dir
    .load.scripts.on.demand("utils_16s2fun.Rc");    
  }
  return(my.16sfun.anot(mbSetObj, type, pipeline,ggversion));
}

######################################
########## Getter Functions ##########
######################################

GetBetaDiversityStats<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$stat.info);
}

GetStressNMDS<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$beta.stress);
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


######################################
######### Utility Functions ##########
######################################

# Get goods's coverage
ComputeGoods <-function(physeq_object){
  com <- t(get_sample(physeq_object))
  no.seqs <- rowSums(com)
  sing <- com==1
  no.singleton <- apply(sing, 1, sum)
  goods <- 100*(1-no.singleton/no.seqs)
  sample <- row.names(com)
  goods.sum <- cbind(sample, no.singleton, no.seqs, goods)
  goods.sum <- as.data.frame(goods.sum,check.names=FALSE)
  row.names(goods.sum) <- c()
  return(goods.sum)
}

# Utility function that performs rarefaction
ggrare2 <- function(physeq_object, data.src, label = NULL, color = NULL, plot = TRUE, linetype = NULL, se = FALSE, step=5) {
  
  x <- methods::as(otu_table(physeq_object), "matrix")
  
  if (taxa_are_rows(physeq_object)) { x <- t(x) }
  
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
    
if(length(y)==1){
 current.msg <<- paste0("All the feature counts in sample ",rownames(x)[i]," is less than 5 which is necessary for rarefy. Please check the data or change the method for filteration.")
 return(0)
}
    if (nrow(y) != 1) {
      rownames(y) <- c(".S", ".se")
      return(data.frame(t(y), Size = n, Sample = rownames(x)[i],check.names=FALSE));
    } else {
      return(data.frame(.S = y[1, ], Size = n, Sample = rownames(x)[i],check.names=FALSE));
    }
  }
  
  f_n <- paste(data.src, step, "rds", sep = ".");
  print(nr)
  out <- lapply(seq_len(nr), rarefun)
  df <- do.call(rbind, out);
  
  # Get sample data
  if (!is.null(sample_data(physeq_object, FALSE))) {
    sdf <- methods::as(sample_data(physeq_object), "data.frame")
    sdf$Sample <- rownames(sdf);
    data <- merge(df, sdf, by = "Sample")
    labels <- data.frame(x = tot, y = S, Sample = rownames(x),check.names=FALSE)
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

generateColorArr <- function(grp.num, filenm=NULL) {
  grp.num <- as.numeric(grp.num);
  library(RColorBrewer);
  colArr<-brewer.pal(grp.num,"Dark2")
  if(is.null(filenm)){
    return(colArr);
  }else{
    sink(filenm);
    cat(rjson::toJSON(colArr));
    sink();
    return(filenm);
  }
}

PlotlyTaxaAbundance <- function(rdaName, type){
load(rdaName);
p <- plotly::ggplotly(box, width=1000, height=800);

if(type=="area"){
   narm <- p[["x"]][["data"]]
   for(i in 1:length(narm)){
      narm[[i]]$y[is.na(narm[[i]]$y)]=0;
   }
   p[["x"]][["data"]] <- narm
}

return(p);
}