
#'Plot PCA plot for meta-analysis samples
#'@description 
#'@param imgNm name of the image to output
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'

PlotMetaPCA <- function(imgNm, dpi, format,factor="NA"){
  microbiome.meta <- qs::qread("microbiome_meta.qs");
  x <- microbiome.meta[["data"]];
  dpi <- as.numeric(dpi);
  #imgNm <- paste(imgNm, "_dpi", dpi, ".", format, sep="");
  imgNm <- paste(imgNm, ".", format, sep="");
  require('lattice');
  require('ggplot2');
  pca <- prcomp(t(na.omit(x)));
  imp.pca<-summary(pca)$importance;
  xlabel <- paste0("PC1"," (", 100*round(imp.pca[2,][1], 3), "%)")
  ylabel <- paste0("PC2"," (", 100*round(imp.pca[2,][2], 3), "%)")
  names <- colnames(x);
  pca.res <- as.data.frame(pca$x);
  # increase xlim ylim for text label
  xlim <- GetExtendRange(pca.res$PC1);
  ylim <- GetExtendRange(pca.res$PC2);

  #print(dim(pca.res));

  if(factor == "NA"){
    Conditions <- factor(microbiome.meta$cls.lbl);
  }else{
    Conditions <- factor(microbiome.meta$cls.lbl); #temporarily
  }
  Datasets <- factor(microbiome.meta$data.lbl)
  pca.res$Conditions <- Conditions;
  pca.res$Datasets <- Datasets;
  pcafig <- ggplot(pca.res, aes(x=PC1, y=PC2,  color=Conditions ,shape=Datasets)) +
    geom_point(size=4, alpha=0.5) + 
    xlim(xlim)+ ylim(ylim) + 
    xlab(xlabel) + ylab(ylabel) + 
    theme_bw()
  
  Cairo(file=imgNm, width=9, height=8, type=format, bg="white", unit="in", dpi=dpi);
  print(pcafig);
  dev.off();
  
}

#'Plot density plot of datasets in meta-analysis module
#'@description 
#'@param imgNm name of the image to output
#'@param format File format of output image (default to png)
#'@param dpi dots per inch resolution of image (default 72)
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PlotMetaDensity<- function(imgNm, dpi=72, format="png", factor=""){
  require("ggplot2")
  microbiome.meta <- qs::qread("microbiome_meta.qs");
  dat <- microbiome.meta$data;
  #imgNm <- paste(imgNm, "_dpi", dpi, ".", format, sep="");
  imgNm <- paste(imgNm, ".", format, sep="");
  dpi <- as.numeric(dpi);
  
  df <- data.frame(microbiome.meta$data, stringsAsFactors = FALSE);
  df <- stack(df);

  Factor <- microbiome.meta$data.lbl;
  
  conv <- data.frame(ind=colnames(microbiome.meta$data), class=Factor);
  conv$ind <- gsub("-", ".", conv$ind);
  df1 <- merge(df, conv, by="ind");
  Cairo(file=imgNm, width=10, height=6, type=format, bg="white", dpi=dpi, unit="in");
  g =ggplot(df1, aes(x=values)) + 
        geom_line(aes(color=class, group=ind), stat="density", alpha=0.3) + 
        geom_line(aes(color=class), stat="density", alpha=0.6, size=1.5) +
        theme_bw()
  print(g);
  dev.off();
}


#'Perform batch correction of meta-analysis datasets
#'@description 
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'
PerformBatchCorrection <- function(){
    .prepare.batch();
    .perform.computing();
    return(dataSets);
    # no need to , already done
}


.prepare.batch<-function(){
  my.fun <- function(){
    require('MMUPHin');
    require('phyloseq');
    phyobj <- qs::qread("merged.data.qs");
    #phyobj <- subsetPhyloseqByDataset(NA, phyobj);
    sam.data <- as.data.frame(as.matrix(sample_data(phyobj)));
    
    cov <- colnames(sam.data)[1];
    otu.tbl <- otu_table(phyobj);
    #print(cov);

    microbiome.meta <- qs::qread("microbiome_meta.qs");
    data.lbl <- microbiome.meta$data.lbl;
    
    fit_adjust_batch <- adjust_batch(feature_abd = otu.tbl,
                                     batch = "dataset",
                                     covariates = cov,
                                     data = sam.data,
                                     control = list(verbose = FALSE))
    
    adj.otu.tbl <- fit_adjust_batch$feature_abd_adj;
    phyobj@otu_table <- otu_table(adj.otu.tbl, taxa_are_rows = TRUE);
    
    shadow_save(phyobj, "merged.data.norm.qs");
    microbiome.meta$data <- adj.otu.tbl;
    shadow_save(microbiome.meta, "microbiome_meta.qs");
    merged.data <- transform_sample_counts(phyobj, function(x) x / sum(x) );
    shadow_save(merged.data, "merged.data.qs");
  }
  dat.in <- list(my.fun=my.fun);
  qs:::qsave(dat.in, file="dat.in.qs");
  return(1);
}

