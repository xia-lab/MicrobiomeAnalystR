##################################################
## R script for MicrobiomeAnalyst
## Description: generating different graphics output
## Authors: Achal Dhariwal(achal.dhariwal@mail.mcgill.ca), Jeff Xia (jeff.xia@mcgillca)
## Last updated: 09/01/2017
###################################################

#'Main function to plot tree graphics
#'@description This functions creates tree plots from the microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotTreeGraph<-function(mbSetObj, plotNm, distnm, clstDist, metadata, datatype,
                        taxrank, format="png", dpi=72, width=NA){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  set.seed(2805619);
  plotNm <- paste(plotNm,".", format, sep="");
  variable<<-metadata;
    
  data <- mbSetObj$dataSet$norm.phyobj;
    
  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloslim(data,mbSetObj$dataSet$taxa_table);
  }else{
    data <- data;
  }
    
  #using by default names for shotgun data
  if(datatype=="metageno"){
    taxrank<-"OTU";
  }

  if(taxrank=="OTU"){
    data<-data;
  }else{
    #merging at taxonomy levels
    data<-fast_tax_glom_first(data,taxrank);
  }

  hc.cls <-as.factor(sample_data(data)[[variable]]);

  # must call distance within the phyloslim package
  #dist.mat<-phyloslim::distance(data,distnm,type = "samples");
  if(distnm == "unifrac" | distnm == "wunifrac"){
    pg_tree <- readRDS("tree.RDS");
    pg_tb <- tax_table(data);
    pg_ot <- otu_table(data);
    pg_sd <- sample_data(data);
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
    data <- merge_phyloslim(pg_tb, pg_ot, pg_sd, pg_tree);

    if(!is.rooted(phy_tree(data))){
      pick_new_outgroup <- function(tree.unrooted){
        treeDT <-
          cbind(cbind(
            data.table(tree.unrooted$edge),
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
    dist.mat<-distance(data,distnm,type = "samples")
  } else {
    dist.mat<-distance(data,distnm,type = "samples")
  }

  # build the tree
  hc_tree<-hclust(dist.mat, method=clstDist);
  mbSetObj$imgSet$tree<-plotNm;

  if(is.na(width)){
    w <- minH <- 650;
    myH <- nsamples(data)*16 + 150;
        
    if(myH < minH){
      myH <- minH;
    }
    w <- round(w/72,2);
    h <- round(myH/72,2);
  }
    
  Cairo::Cairo(file=plotNm, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  par(mar=c(4,2,2,10));
  clusDendro<-as.dendrogram(hc_tree);
  cols <- GetColorSchema();
  names(cols) <- sample_names(data);
  labelColors <- cols[hc_tree$order];
    
  colLab <- function(n){
        if(is.leaf(n)) {
            a <- attributes(n);
            labCol <- labelColors[a$label];
            attr(n, "nodePar") <-
                        if(is.list(a$nodePar)){
                            c(a$nodePar,lab.col = labCol,pch=NA)
                        }else{
                            list(lab.col = labCol,pch=NA)
                        }
        }
        n
  }
    
  clusDendro<-dendrapply(clusDendro, colLab);
  plot(clusDendro,horiz=T,axes=T);
  par(cex=1);
  legend.nm <- unique(as.character(hc.cls));
  legend.nm <-gsub("\\.", " ",legend.nm)
  legend("topleft", legend = legend.nm, pch=15, col=unique(cols), bty = "n");
  dev.off();
  mbSetObj$analSet$tree<-hc_tree;
  mbSetObj$analSet$tree.dist<-distnm;
  mbSetObj$analSet$tree.clust<-clstDist;
  mbSetObj$analSet$tree.taxalvl<-taxrank;
  return(.set.mbSetObj(mbSetObj))
}

#######################################
###########DE(feature boxplot)#########
#######################################

#'Function to create box plots of important features
#'@description This functions plots box plots of a selected feature.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import grid
#'@import gridExtra
PlotBoxData<-function(mbSetObj, boxplotName, feat, format="png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_ggplot();
    load_grid();
    load_gridExtra();
  }
  
  data <- mbSetObj$analSet$boxdata;
  a <- data[,feat];
  ind <- which(a=="0");
  a[ind] <- 0.1;
  data$log_feat <- log(a);
  boxplotName = paste(boxplotName,".",format, sep="");
  Cairo::Cairo(file=boxplotName,width=720, height=360, type=format, bg="white",dpi=dpi);
    
  box=ggplot(data,aes(x=data$class, y = data[,feat])) + stat_boxplot(geom ='errorbar') + geom_boxplot(aes(color=class), outlier.shape = NA) +
      geom_boxplot(aes(fill=class), outlier.shape = NA) + geom_jitter() + theme_bw() + labs(y="Abundance",x="class") +
      ggtitle("Orignal Count") + theme(plot.title = element_text(hjust=0.5), legend.position="none");
    
  box1=ggplot(data,aes(x=data$class, y = data$log_feat)) + stat_boxplot(geom ='errorbar') + geom_boxplot(aes(color=class), outlier.shape = NA) +
      geom_boxplot(aes(fill=class), outlier.shape = NA) + geom_jitter() + theme_bw() + labs(y="",x="class") +
      ggtitle("Log-transformed Count") + theme(plot.title = element_text(hjust=0.5));
  
  grid.arrange(ggplotGrob(box), ggplotGrob(box1),ncol=2,nrow=1,top=feat);
  dev.off();
  return(.set.mbSetObj(mbSetObj))
}

###############################
###########Heatmap#############
###############################

#'Main function to plot heatmap.
#'@description This functions plots a heatmap from the mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import pheatmap

PlotHeatmap<-function(mbSetObj, plotNm, smplDist, clstDist, palette, metadata,
                      taxrank, datatype, viewOpt, doclust, format="png", showfeatname,
                      appendnm, rowV=F, colV=T, var.inx=NA, border=T, width=NA, dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(.on.public.web){
    load_pheatmap();
  }
  
  set.seed(2805614);
  #used for color pallete
  variable <<- metadata;
  data <- mbSetObj$dataSet$norm.phyobj;
    
  if(datatype=="16S"){
    mbSetObj$dataSet$taxa_table <- tax_table(mbSetObj$dataSet$proc.phyobj);
    data <- merge_phyloslim(data, mbSetObj$dataSet$taxa_table);
  }else{
    data <- data;
  }

  #using by default names for shotgun data
  if(datatype=="metageno"){
    taxrank <- "OTU";
  }

  #if more than 1500 features will be present;subset to most abundant=>1500 features.
  #OTUs already in unique names;
  if(ntaxa(data)>1500){
    data = prune_taxa(names(sort(taxa_sums(data), TRUE))[1:1500], data);
    viewOpt == "overview";
  }

  if(taxrank=="OTU"){
    data <- data;
    nm <- taxa_names(data);
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
  }else{
    #merging at taxonomy levels
    data <- fast_tax_glom_first(data,taxrank);
    nm <- as.character(tax_table(data)[,taxrank]);
    y <- which(is.na(nm)==TRUE);
    #converting NA values to unassigned
    nm[y] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
        
    if(appendnm=="T"){
      all_nm <- colnames(tax_table(data));
      hg_nmindx <- which(all_nm==taxrank)-1;
            
      if(hg_nmindx!=0){
        nma <- as.character(tax_table(data)[,hg_nmindx]);
        y1 <- which(is.na(nma)==TRUE);
        nma[y1] <- "Not_Assigned";
        nm <- paste0(nma,"_",nm);
        ind <- which(nm=="Not_Assigned_Not_Assigned");
        nm[ind] <- "Not_Assigned";
        nm <- gsub("_Not_Assigned", "",nm, perl = TRUE);
      }
    }
      
    rownames(data1) <- nm;
    #all NA club together
    data1 <- (t(sapply(by(data1,rownames(data1),colSums),identity)));
    nm <- rownames(data1);
  }

  rownames(data1) <- nm;
  # arrange samples on the basis of slected experimental factor and using the same for annotation also
  annotation <- data.frame(sample_data(data));

  ind <- which(colnames(annotation)!=metadata && colnames(annotation)!="sample_id");
    
  if(length(ind)>0){
    ind1 <- ind[1];
    annotation <- annotation[order(annotation[,metadata],annotation[,ind1]),];
  }else{
    annotation <- annotation[order(annotation[,metadata]),];
  }

  #there is an additional column sample_id which need to be removed first
  annotation <- subset(annotation,select = -sample_id);
  sam.ord <- rownames(annotation);
  data1 <- data1[,sam.ord];

  # set up colors for heatmap
  if(palette=="gbr"){
    colors <- grDevices::colorRampPalette(c("green", "black", "red"), space="rgb")(256);
  }else if(palette == "heat"){
    colors <- grDevices::heat.colors(256);
  }else if(palette == "topo"){
    colors <- grDevices::topo.colors(256);
  }else if(palette == "gray"){
    colors <- grDevices::colorRampPalette(c("grey90", "grey10"), space="rgb")(256);
  }else{
    
    if(.on.public.web){
      load_rcolorbrewer();
    }
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(256));
  }

  #setting the size of plot
  if(is.na(width)){
    minW <- 800;
    myW <- ncol(data1)*18 + 200;
    if(myW < minW){
      myW <- minW;
    }
    w <- round(myW/72,2);
  }

  myH <- nrow(data1)*18 + 150;
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

  if(border){
    border.col<-"grey60";
  }else{
    border.col <- NA;
  }
  
  plotNm = paste(plotNm, ".", format, sep="");
  mbSetObj$imgSet$heatmap<-plotNm;

  if(format=="pdf"){
    grDevices::pdf(file = plotNm, width=w, height=h, bg="white", onefile=FALSE);
  }else{
    Cairo::Cairo(file = plotNm, unit="in", dpi=dpi, width=w, height=h, type=format, bg="white");
  }
    
  # set up color schema for samples
  if(palette== "gray"){
    cols <- GetColorSchema(T);
    uniq.cols <- unique(cols);
  }else{
    cols <- GetColorSchema();
    uniq.cols <- unique(cols);
  }

  if(doclust=="T"){
    rowV<-T;
  }
    
  if(showfeatname=="T"){
    showfeatname<-T;
  } else {
    showfeatname<-F;
  }

  pheatmap::pheatmap(data1,
        annotation=annotation,
        fontsize=8, fontsize_row=8,
        clustering_distance_rows = smplDist,
        clustering_distance_cols = smplDist,
        clustering_method = clstDist,
        show_rownames = showfeatname,
        border_color = border.col,
        cluster_rows = colV,
        cluster_cols = rowV,
        scale= "row",
        color = colors
        );
  
  dev.off();
    
  # storing for Report Generation
  mbSetObj$analSet$heatmap<-data1;
  mbSetObj$analSet$heatmap.dist<-smplDist;
  mbSetObj$analSet$heatmap.clust<-clstDist;
  mbSetObj$analSet$heat.taxalvl<-taxrank;
  return(.set.mbSetObj(mbSetObj))
}

#'Function to get color palette for graphics.
#'@description This function is called to create a color palette
#'based on the number of groups. It returns a vector of color
#'hex codes based on the number of groups.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetColorSchema <- function(mbSetObj, grayscale=F){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
    
  # test if total group number is over 9
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]]);
  grp.num <- length(levels(claslbl));

  if(grayscale){
    dist.cols <- grDevices::colorRampPalette(c("grey90", "grey30"))(grp.num);
    lvs <- levels(claslbl);
    colors <- vector(mode="character", length=length(claslbl));
        
    for(i in 1:length(lvs)){
      colors[mbSetObj$analSet$cls == lvs[i]] <- dist.cols[i];
    }
  }else if(grp.num > 9){
    pal12 = c("#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
              "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
              "#FFFF99", "#B15928");
    dist.cols <- grDevices::colorRampPalette(pal12)(grp.num);
    lvs <- levels(claslbl);
    colors <- vector(mode="character", length=length(mbSetObj$analSet$cls));

    for(i in 1:length(lvs)){
      colors[claslbl == lvs[i]] <- dist.cols[i];
    }
  }else{
    if(exists("colVec") && !any(colVec =="#NA") ){
      cols <- vector(mode="character", length=length(claslbl));
      clsVec <- as.character(claslbl);
      grpnms <- names(colVec);
            
      for(i in 1:length(grpnms)){
        cols[clsVec == grpnms[i]] <- colVec[i];
      }
      colors <- cols;
    }else{
      colors <- as.numeric(claslbl)+1;
    }
  }
  return (colors);
}
