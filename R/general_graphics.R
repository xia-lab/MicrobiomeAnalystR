##################################
########### 3D PCoA/PCA############
##################################

#' Main function to perform PCoA analysis
#' @description This functions creates a 3D PCoA plot from the microbiome data.
#' This is used by the Beta-Diversity analysis.
#' The 3D interactive visualization is on the web.
#' @param mbSetObj Input the name of the mbSetObj.
#' @param ordMeth Character, input the name
#' of the ordination method. "PCoA" for principal coordinate analysis and "NMDS" for
#' non-metric multidimensional scaling.
#' @param distName Character, input the name of the distance method.
#' @param datatype Character, input "16S" if the data is marker
#' gene data and "metageno" if it is metagenomic data.
#' @param taxrank Character, input the taxonomic
#' level for beta-diversity analysis.
#' @param colopt Character, color the data points by the experimental factor,
#' the taxon abundance of a selected taxa, or alpha diversity.
#' @param variable Character, input the name of the experimental factor.
#' @param taxa Character, if the data points are colored by taxon abundance,
#' input the name of the selected taxa.
#' @param alphaopt Character, if the data points are colored by alpha-diversity,
#' input the preferred alpha-diversity measure.
#' @param jsonNm Character, input the name of the json file to output.
#' @author Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import vegan

PCoA3D.Anal <- function(mbSetObj, ordMeth, distName, taxrank, colopt, variable, taxa, alphaopt, jsonNm) {
  if (!exists("my.pcoa.3d")) { # public web on same user dir
    .load.scripts.on.demand("utils_pcoa3d.Rc")
  }

  return(my.pcoa.3d(mbSetObj, ordMeth, distName, taxrank, colopt, variable, taxa, alphaopt, jsonNm))
}

#' Function to plot tree graphics for dendogram.
#' @description This functions creates dendogram tree plots.
#' @param mbSetObj Input the name of the mbSetObj.
#' @param plotNm Character, input the name of the plot.
#' @param distnm Character, input the name of the selected
#' distance measure. "bray" for Bray-Curtis Index, "jsd" for
#' Jensen-Shannon Divergence, "jaccard" for Jaccard Index,
#' "unifrac" for Unweighted Unifrac Distance, and "wunifrac" for weighted
#' unifrac distance.
#' @param clstDist Character, input the name of the
#' selected clustering algorithm. "ward" for Ward, "average" for Average,
#' "complete" for Complete, and "single" for Single.
#' @param metadata Character, input the name of the experimental factor.
#' @param datatype Character, "16S" if marker gene data and
#' "metageno" if shotgun metagenomic data.
#' @param colorOpts Character, "default" or "viridis".
#' @param taxrank Character, input the taxonomic level to perform
#' classification. For instance, "OTU-level" to use OTUs.
#' @param format Character, by default the plot is .png format.
#' @param dpi The dots per inch. Numeric, by default it is set to 72.
#' @param width Width of the plot. Numeric, by default it is set to NA.
#' @param plotType Character, default is set to "rectangle". Alternative is
#' "triangle".
#' @author Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import ape
#' @import viridis
PlotTreeGraph <- function(mbSetObj, plotNm, distnm, clstDist, metadata,
                          taxrank, colorOpts, format = "png", dpi = 72, width = NA,
                          plotType = "rectangle") {
  load_ape()
  load_viridis()
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj)

  set.seed(2805619)
  plotNm <- paste(plotNm, ".", format, sep = "")
  variable <<- metadata

  data <- mbSetObj$dataSet$norm.phyobj

  if (mbSetObj$module.type == "mdp") {
    if (!(exists("phyloseq_objs"))) {
      phyloseq_objs <- readDataQs("phyloseq_objs.qs", mbSetObj$module.type, dataName)
    }
    data <- phyloseq_objs$merged_obj[[taxrank]]
    if (is.null(data)) {
      AddErrMsg("Errors in projecting to the selected taxanomy level!")
      return(0)
    }
  } else {
    data <- data
  }

  # using by default names for shotgun data
  if (mbSetObj$module.type == "sdp") {
    taxrank <- "OTU"
  }



  hc.cls <- as.factor(sample_data(data)[[variable]])

  # must call distance within the phyloslim package
  if (distnm == "unifrac" | distnm == "wunifrac") {
    pg_tree <- qs::qread("tree.qs")
    pg_tb <- tax_table(data)
    pg_ot <- otu_table(data)
    pg_sd <- sample_data(data)
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree)
    data <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree)

    if (!is.rooted(phy_tree(data))) {
      pick_new_outgroup <- function(tree.unrooted) {
        treeDT <-
          cbind(
            cbind(
              data.table(tree.unrooted$edge),
              data.table(length = tree.unrooted$edge.length)
            )[1:Ntip(tree.unrooted)],
            data.table(id = tree.unrooted$tip.label)
          )
        new.outgroup <- treeDT[which.max(treeDT$length), ]$id
        return(new.outgroup)
      }
      new.outgroup <- pick_new_outgroup(phy_tree(data))
      phy_tree(data) <- ape::root(phy_tree(data),
        outgroup = new.outgroup,
        resolve.root = TRUE
      )
    }
    dist.mat <- distance(data, distnm, type = "samples")
  } else {
    dist.mat <- distance(data, distnm, type = "samples")
  }

  # build the tree
  hc_tree <- hclust(dist.mat, method = clstDist)
  mbSetObj$imgSet$tree <- plotNm

  if (is.na(width)) {
    w <- minH <- 650
    myH <- nsamples(data) * 16 + 150

    if (myH < minH) {
      myH <- minH
    }
    w <- round(w / 72, 2)
    h <- round(myH / 72, 2)
  }

  Cairo::Cairo(file = plotNm, unit = "in", dpi = dpi, width = w, height = h, type = format, bg = "white")
  par(mar = c(4, 2, 2, 10))
  clusDendro <- as.dendrogram(hc_tree)

  if (colorOpts == "default") {
    cols <- GetColorSchema(mbSetObj)
  } else {
    claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]])
    grp.num <- length(levels(claslbl))

    if (colorOpts == "viridis") {
      cols <- viridis::viridis(grp.num)
    } else if (colorOpts == "plasma") {
      cols <- viridis::plasma(grp.num)
    } else if (colorOpts == "cividis") {
      cols <- viridis::cividis(grp.num)
    }

    lvs <- levels(claslbl)
    colors <- vector(mode = "character", length = length(mbSetObj$analSet$cls))

    for (i in 1:length(lvs)) {
      colors[claslbl == lvs[i]] <- cols[i]
    }

    cols <- colors
  }

  names(cols) <- sample_names(data)
  labelColors <- cols[hc_tree$order]

  colLab <- function(n) {
    if (is.leaf(n)) {
      a <- attributes(n)
      labCol <- labelColors[a$label]
      attr(n, "nodePar") <-
        if (is.list(a$nodePar)) {
          c(a$nodePar, lab.col = labCol, pch = NA)
        } else {
          list(lab.col = labCol, pch = NA)
        }
    }
    n
  }

  clusDendro <- dendrapply(clusDendro, colLab)
  plot(clusDendro, horiz = T, axes = T, type = plotType)
  par(cex = 1)
  legend.nm <- unique(as.character(hc.cls))
  legend.nm <- gsub("\\.", " ", legend.nm)
  legend("topleft", legend = legend.nm, pch = 15, col = unique(cols), bty = "n")
  dev.off()

  mbSetObj$analSet$tree <- hc_tree
  mbSetObj$analSet$tree.dist <- distnm
  mbSetObj$analSet$tree.clust <- clstDist
  mbSetObj$analSet$tree.taxalvl <- taxrank
  return(.set.mbSetObj(mbSetObj))
}

#######################################
########### DE(feature boxplot)#########
#######################################

#' Function to create box plots of important features
#' @description This functions plots box plots of a selected feature.
#' @param mbSetObj Input the name of the mbSetObj.
#' @param boxplotName Character, input the name of the
#' box plot.
#' @param feat Character, input the name of the selected
#' feature.
#' @param format Character, by default the plot format
#' is "png".
#' @param dpi Dots per inch. Numeric, by default
#' it is set to 72.
#' @author Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import grid
#' @import gridExtra
PlotBoxData <- function(mbSetObj, boxplotName, feat, plotType, format = "png", dpi = 72) {
  mbSetObj <- .get.mbSetObj(mbSetObj)

  load_ggplot()
  load_grid()
  load_gridExtra()
  load_phyloseq()

  print(paste("plottype==", plotType))
  print(feat)
  variable <- mbSetObj$analSet$var.type

  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL = TRUE)

  if (is.null(variable)) {
    variable <- colnames(sample_table)[1]
  }
  data <- mbSetObj$analSet$boxdata
  a <- as.numeric(data[, feat])
  min.val <- min(abs(a[a != 0])) / 5
  data$log_feat <- log2((a + sqrt(a^2 + min.val)) / 2)
  boxplotName <- paste(boxplotName, ".", format, sep = "")
  if (is.null(mbSetObj$imgSet$boxplots)) {
    mbSetObj$imgSet$boxplots <- boxplotName
    mbSetObj$analSet$boxplot_feats <- feat
  } else {
    idxft <- which(mbSetObj$analSet$boxplot_feats == feat)
    if (length(idxft) == 0) {
      mbSetObj$analSet$boxplot_feats <- c(mbSetObj$analSet$boxplot_feats, feat)
      mbSetObj$imgSet$boxplots <- c(mbSetObj$imgSet$boxplots, boxplotName)
    } else {
      mbSetObj$imgSet$boxplots[idxft] <- boxplotName
    }
  }

  Cairo::Cairo(file = boxplotName, width = 720, height = 360, type = format, bg = "white", dpi = dpi)

    box <- ggplot2::ggplot(data, aes(x = data$class, y = data[, feat], fill = data$class))  
    box1 <- ggplot2::ggplot(data, aes(x = data$class, data$log_feat, fill = data$class)) 
    
    if(plotType == "violin"){
          box <- box + geom_violin(trim=FALSE) 
         box1 <- box1 + geom_violin(trim=FALSE) 
        } else {
          box <- box + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA)
          box1 <- box1 + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA)
        }

     box <- box+theme_bw() + geom_jitter(size=1) + 
            stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)+
             labs(y = "", x = variable, fill = variable) + theme(axis.text.x = element_text(angle=90, hjust=1))+
              ggtitle("Filtered Count") +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
             theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

    
     box1 <- box1 +theme_bw() + geom_jitter(size=1) + 
            stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)+
             labs(y = "", x = variable, fill = variable) + theme(axis.text.x = element_text(angle=90, hjust=1))+
             ggtitle("Log-transformed Count") +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
             theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

 
  grid.arrange(ggplotGrob(box), ggplotGrob(box1), ncol = 2, nrow = 1, top = feat)
  dev.off()
  return(.set.mbSetObj(mbSetObj))
}

PlotBoxMultiData <- function(mbSetObj, boxplotName, analysis.var, feat, plotType, format = "png", dpi = 72) {
  mbSetObj <- .get.mbSetObj(mbSetObj)

  load_ggplot()
  load_grid()
  load_gridExtra()
  load_phyloseq()

  variable <- analysis.var
  var.type <- mbSetObj$dataSet$meta.types[names(mbSetObj$dataSet$meta.types) == variable]

  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL = TRUE)

  data <- mbSetObj$analSet$multiboxdata
  is.norm <- unique(data[, "norm"])
  a <- as.numeric(data[, feat])
  min.val <- min(abs(a[a != 0])) / 5
  data$log_feat <- log2((a + sqrt(a^2 + min.val)) / 2)
  boxplotName <- paste(boxplotName, ".", format, sep = "")

  if (is.norm == "false") {
    feat2 <- "log_feat"
    plot.title <- "Log-transformed Counts"
  } else {
    feat2 <- feat
    plot.title <- paste0(feat, ": Normalized Counts")
  }

  


  if (var.type == "disc") {
     plot1 <- ggplot2::ggplot(data, aes(x = data$class, y = data[, feat], fill = data$class))  
     plot2 <- ggplot2::ggplot(data, aes(x = data$class, data$log_feat, fill = data$class)) 
    
    if(plotType == "violin"){
        plot1 <- plot1 + geom_violin(trim=FALSE) 
         plot2 <- plot2 + geom_violin(trim=FALSE) 
        } else {
          plot1 <- plot1 + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA)
          plot2 <- plot2 + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA)
        }

     plot1 <- plot1+theme_bw() + geom_jitter(size=1) + 
            stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)+
             labs(y = "", x = variable, fill = variable) + theme(axis.text.x = element_text(angle=90, hjust=1))+
              ggtitle("Filtered Count") +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
             theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    
     plot2 <- plot2 +theme_bw() + geom_jitter(size=1) + 
            stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)+
             labs(y = "", x = variable, fill = variable) + theme(axis.text.x = element_text(angle=90, hjust=1))+
             ggtitle("Log-transformed Count") +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
             theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())

  } else {
    plot1 <- ggplot(data, aes(x = as.numeric(data$class), y = .data[[feat]])) +
      geom_point() +
      geom_smooth(method = lm) +
      theme_bw() +
      labs(y = "Abundance", x = variable) +
      ggtitle("Raw Counts") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

    plot2 <- ggplot(data, aes(x = as.numeric(data$class), y = .data[[feat2]])) +
      geom_point() +
      geom_smooth(method = lm) +
      theme_bw() +
      labs(y = "", x = variable, fill = variable) +
      ggtitle(plot.title) +
      theme(plot.title = element_text(hjust = 0.5))
  }

  if (is.norm == "false") {
    Cairo::Cairo(file = boxplotName, width = 720, height = 360, type = format, bg = "white", dpi = dpi)
    grid.arrange(ggplotGrob(plot1), ggplotGrob(plot2), ncol = 2, nrow = 1, top = feat)
    dev.off()
  } else {
    Cairo::Cairo(file = boxplotName, width = 720, height = 360, type = format, bg = "white", dpi = dpi)
    grid.arrange(ggplotGrob(plot2), ncol = 1, nrow = 1)
    dev.off()
  }

  return(.set.mbSetObj(mbSetObj))
}



PlotBoxMultiMetabo <- function(mbSetObj, boxplotName, analysis.var, feat,plotType, format = "png", dpi = 72) {
  mbSetObj <- .get.mbSetObj(mbSetObj)

  load_ggplot()
  load_grid()
  load_gridExtra()
  load_phyloseq()

  variable <- analysis.var
  var.type <- mbSetObj$dataSet$meta.types[names(mbSetObj$dataSet$meta.types) == variable]

  sample_table <- sample_data(mbSetObj$dataSet$proc.phyobj, errorIfNULL = TRUE)

  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[analysis.var]])

  dt.orig <- t(mbSetObj$dataSet$metabolomics$filt.data)
  dt.norm <- t(mbSetObj$dataSet$metabolomics$norm.data)

  boxplotName <- paste(boxplotName, ".", format, sep = "")


  if (var.type == "disc") {
    df.orig <- data.frame(value = as.numeric(dt.orig[, feat]), name = claslbl)
   df.norm <- data.frame(value = as.numeric(dt.norm[, feat]), name = claslbl)

    plot1 <- ggplot2::ggplot(data = df.orig, aes(x =name, y = value, fill = name))  
     plot2 <- ggplot2::ggplot(data = df.norm, aes(x =name, y = value, fill = name)) 
    
    if(plotType == "violin"){
        plot1 <- plot1 + geom_violin(trim=FALSE) 
         plot2 <- plot2 + geom_violin(trim=FALSE) 
        } else {
          plot1 <- plot1 + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA)
          plot2 <- plot2 + geom_boxplot(notch=FALSE, outlier.shape = NA, outlier.colour=NA)
        }

     plot1 <- plot1+theme_bw() + geom_jitter(size=1) + 
            stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)+
             labs(y = "", x = variable, fill = variable) + theme(axis.text.x = element_text(angle=90, hjust=1))+
                  ggtitle("Original Data") +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
             theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())
    
     plot2 <- plot2 +theme_bw() + geom_jitter(size=1) + 
            stat_summary(fun=mean, colour="yellow", geom="point", shape=18, size=3, show.legend = FALSE)+
             labs(y = "", x = variable, fill = variable) + theme(axis.text.x = element_text(angle=90, hjust=1))+
                  ggtitle("Normalized Data") +
            theme(plot.title = element_text(hjust = 0.5), legend.position = "none")+
             theme(panel.grid.minor = element_blank(), panel.grid.major = element_blank())


  } else {
    df.orig <- data.frame(value = as.numeric(dt.orig[, feat]), as.numeric(as.character(claslbl)))

    plot1 <- ggplot(df.orig, aes(x = name, y = value)) +
      geom_point() +
      geom_smooth(method = lm) +
      theme_bw() +
      labs(y = "Abundance", x = variable) +
      ggtitle("Raw Counts") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")

    df.norm <- data.frame(value = as.numeric(dt.norm[, feat]), as.numeric(as.character(claslbl)))

    plot1 <- ggplot(df.norm, aes(x = name, y = value)) +
      geom_point() +
      geom_smooth(method = lm) +
      theme_bw() +
      labs(y = "Abundance", x = variable) +
      ggtitle("Raw Counts") +
      theme(plot.title = element_text(hjust = 0.5), legend.position = "none")
  }

  if (mbSetObj$dataSet$metabolomics$isNormInput == "false") {
    Cairo::Cairo(file = boxplotName, width = 720, height = 360, type = format, bg = "white", dpi = dpi)
    grid.arrange(ggplotGrob(plot1), ggplotGrob(plot2), ncol = 2, nrow = 1, top = feat)
    dev.off()
  } else {
    Cairo::Cairo(file = boxplotName, width = 720, height = 360, type = format, bg = "white", dpi = dpi)
    grid.arrange(ggplotGrob(plot2), ncol = 1, nrow = 1)
    dev.off()
  }

  message("metabo box plot done")
  return(.set.mbSetObj(mbSetObj))
}

###############################
########### Heatmap#############
###############################

#' Main function to plot heatmap.
#' @description This functions plots a heatmap from the mbSetObj.
#' @param mbSetObj Input the name of the mbSetObj.
#' @param plotNm Character, input the name
#' of the plot.
#' @param smplDist Input the distance measure. "euclidean" for
#' Euclidean distance, "correlation" for Pearson, and "minkowski"
#' for Minkowski.
#' @param clstDist Character, input the name of the
#' selected clustering algorithm. "ward" for Ward, "average" for Average,
#' "complete" for Complete, and "single" for Single.
#' @param palette Set the colors of the heatmap. By default it
#' is set to "bwm", blue, white, to red. Use "gbr" for green, black, red, use
#' "heat" for red to yellow, "topo" for blue to yellow, "gray" for
#' white to black, and "byr" for blue, yellow, red.
#' @param metadata Character, input the name of the experimental factor
#' to cluster samples by.
#' @param taxrank Character, input the taxonomic level to perform
#' classification. For instance, "OTU-level" to use OTUs.
#' @param viewOpt Character, "overview" to view an overview
#' of the heatmap, and "detail" to iew a detailed view of the
#' heatmap (< 1500 features).
#' @param doclust Logicial, default set to "F".
#' @param format Character, input the preferred
#' format of the plot. By default it is set to "png".
#' @param colname Logical, specify the if column name is shown, default set to "T",
#' @param rowname Logical, specify the if row name is shown, default set to "T",
#' @param fontsize_row Numeric, fontsize for rownames
#' @param fontsize_col Numeric, fontsize for colnames
#' @param appendnm Logical, "T" to prepend higher taxon names.
#' @param rowV Logical, default set to "F".
#' @param colV Logical, default set to "T".
#' @param var.inx Default set to NA.
#' @paraboxdatam border Logical, show cell borders, default set to "T".
#' @param width Numeric, input the width of the plot. By
#' default it is set to NA.
#' @param dpi Numeric, input the dots per inch. By default
#' it is set to 72.
#' @author Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
#' @import pheatmap
#' @import viridis

PlotHeatmap <- function(mbSetObj, plotNm, dataOpt = "norm",
                        scaleOpt = "row", smplDist, clstDist, palette, metadata,
                        taxrank, viewOpt, doclust, format = "png", showColnm, showRownm,
                        unitCol, unitRow, fzCol, fzRow, annoPer, fzAnno,
                        appendnm = "F", rowV = F, colV = T, var.inx = NA, border = T, width = NA, dpi = 72) {
  mbSetObj <- .get.mbSetObj(mbSetObj)
  load_iheatmapr()
  load_rcolorbrewer()
  load_viridis()
  load_phyloseq()

  print(c(unitCol, unitRow, fzCol, fzRow, annoPer, fzAnno))
  print(c(showColnm, showRownm))
  set.seed(2805614)
  # used for color pallete
  variable <<- metadata
  if (dataOpt == "norm") {
    if (mbSetObj$module.type != "mdp") {
      taxrank <- "OTU"
    }

    if (!(exists("phyloseq_objs"))) {
      phyloseq_objs <- readDataQs("phyloseq_objs.qs", mbSetObj$module.type, dataName)
    }

    data <- phyloseq_objs$merged_obj[[taxrank]]
    if (is.null(data)) {
      AddErrMsg("Errors in projecting to the selected taxanomy level!")
      return(0)
    }
  } else {
    data <- mbSetObj$dataSet$proc.phyobj
    if (mbSetObj$module.type != "mdp") {
      taxrank <- "OTU"
    }
  }

  # if more than 500 features will be present;subset to most abundant=>500 features.
  # OTUs already in unique names;
  if (ntaxa(data) > 500) {
    data <- prune_taxa(names(sort(taxa_sums(data), TRUE))[1:500], data)
    viewOpt == "overview"
  }

  if (taxrank == "OTU") {
    data1 <- as.matrix(otu_table(data))
    rownames(data1) <- taxa_names(data)
  } else {
    # merging at taxonomy levels
    data <- fast_tax_glom_mem(data, taxrank)
    if (is.null(data)) {
      AddErrMsg("Errors in projecting to the selected taxanomy level!")
      return(0)
    }
    nm <- as.character(tax_table(data)[, taxrank])
    y <- which(is.na(nm) == TRUE)
    # converting NA values to unassigned
    nm[y] <- "Not_Assigned"
    data1 <- as.matrix(otu_table(data))

    if (appendnm == "T") {
      all_nm <- colnames(tax_table(data))
      hg_nmindx <- which(all_nm == taxrank) - 1

      if (hg_nmindx != 0) {
        nma <- as.character(tax_table(data)[, hg_nmindx])
        y1 <- which(is.na(nma) == TRUE)
        nma[y1] <- "Not_Assigned"
        nm <- paste0(nma, "_", nm)
        ind <- which(nm == "Not_Assigned_Not_Assigned")
        nm[ind] <- "Not_Assigned"
        nm <- gsub("_Not_Assigned", "", nm, perl = TRUE)
      }
    }

    rownames(data1) <- nm
    # all NA club together
    data1 <- (t(sapply(by(data1, rownames(data1), colSums), identity)))
    nm <- rownames(data1)
  }

  # arrange samples on the basis of selected experimental factor and using the same for annotation also
  annotation <- data.frame(sample_data(data))

  ind <- which(!colnames(annotation) %in% c(metadata, "sample_id"))

  if (length(ind) > 0) {
    ind1 <- ind[1]
    annotation <- annotation[order(annotation[, metadata], annotation[, ind1]), ]
  } else {
    annotation <- annotation[order(annotation[, metadata]), ]
  }

  # remove those columns that all values are unique (continuous or non-factors)
  # uniq.inx <- apply(annotation, 2, function(x){length(unique(x)) == length(x)});
  # there is an additional column sample_id which need to be removed first
  # get only good meta-data
  good.inx <- GetDiscreteInx(annotation)
  if (sum(good.inx) > 0) {
    annotation <- annotation[, good.inx, drop = FALSE]
    sam.ord <- rownames(annotation)
    data1 <- data1[, sam.ord]
  } else {
    annotation <- NA
  }
  # set up colors for heatmap

  if (palette == "gbr") {
    colors <- grDevices::colorRampPalette(c("green", "black", "red"), space = "rgb")(256)
  } else if (palette == "heat") {
    colors <- grDevices::heat.colors(256)
  } else if (palette == "topo") {
    colors <- grDevices::topo.colors(256)
  } else if (palette == "gray") {
    colors <- grDevices::colorRampPalette(c("grey90", "grey10"), space = "rgb")(256)
  } else if (palette == "byr") {
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256))
  } else if (palette == "viridis") {
    colors <- rev(viridis::viridis(10))
  } else if (palette == "plasma") {
    colors <- rev(viridis::plasma(10))
  } else if (palette == "npj") {
    colors <- c("#00A087FF", "white", "#E64B35FF")
  } else if (palette == "aaas") {
    colors <- c("#4DBBD5FF", "white", "#E64B35FF")
  } else if (palette == "d3") {
    colors <- c("#2CA02CFF", "white", "#FF7F0EFF")
  } else {
    colors <- c("#0571b0", "#92c5de", "white", "#f4a582", "#ca0020")
  }

  plotjs <- paste0(plotNm, ".json")
  plotwidget <- paste0(plotNm, ".rda")
  plotNmPdf <- paste(plotNm, ".", "pdf", sep = "")
  # print(plotNm)
  mbSetObj$imgSet$heatmap <- plotNmPdf

  if (doclust == "T") {
    rowV <- TRUE
  }

  if (!(is.null(mbSetObj$analSet$heat.taxalvl))) {
    if (mbSetObj$analSet$heat.taxalvl != taxrank) {
      # annoPer <- fontsize_col <- fontsize_row <- fzAnno <- "0.0"
      mbSetObj$analSet$heatVar <- 1
    } else {
      mbSetObj$analSet$heatVar <- 0
    }
  } else {
    mbSetObj$analSet$heatVar <- 0
  }

  data1sc <- as.matrix(apply(data1, 2, as.numeric))
  rownames(data1sc) <- rownames(data1)
  data1sc <- scale_mat(data1sc, scaleOpt)


  data1sc <- round(data1sc, 5)
  w <- min(1500, ncol(data1sc) * unitCol + 50)
  h <- min(1800, nrow(data1sc) * unitRow + 50)

  if (ncol(data1sc) < 100) {
    w <- w + (100 - ncol(data1sc)) * 6
  }

  if (nrow(data1sc) < 100) {
    h <- h + (100 - nrow(data1sc)) * 5
  }

  print(c(w, h))
  idx <- which(mbSetObj$dataSet$meta.types[names(annotation)] == "disc")
  if (any(apply(annotation[, names(idx), drop = F], 2, function(x) length(unique(x))) > 10)) {
    if (h < 750) {
      nr <- 2
      ys <- 0.85
    } else if (h < 1500) {
      nr <- 4.5
      ys <- 0.9
    } else {
      nr <- 9
      ys <- 0.95
    }
  } else {
    if (h < 750) {
      nr <- 3
      ys <- 0.85
    } else if (h < 1500) {
      nr <- 6
      ys <- 0.95
    } else {
      nr <- 11
      ys <- 0.95
    }
  }
  sz <- max(as.numeric(annoPer) / 100, 0.015)
  bf <- min(0.01, (sz / 3))


  dend_row <- hclust(dist(data1sc, method = smplDist), method = clstDist)

  p <- iheatmap(data1sc,
    name = "Abundance", x_categorical = TRUE,
    layout = list(font = list(size = fzAnno)),
    colors = colors,
    colorbar_grid = setup_colorbar_grid(nrows = nr, x_start = 1.1, y_start = ys, x_spacing = 0.15)
  ) %>%
    add_col_annotation(annotation,
      side = "top", size = sz, buffer = bf, inner_buffer = bf / 3
    ) %>%
    add_row_dendro(dend_row, side = "right")


  if (showColnm) {
    p <- p %>%
      add_col_labels(size = 0.2, font = list(size = fzCol))
  }

  if (showRownm) {
    p <- p %>%
      add_row_labels(size = 0.2, font = list(size = fzRow), side = "left")
  }

  if (doclust == "T") {
    dend_col <- hclust(dist(t(data1), method = smplDist), method = clstDist)
    p <- p %>% add_col_dendro(dend_col)
  }

  as_list <- to_plotly_list(p)



  as_list[["layout"]][["width"]] <- w
  as_list[["layout"]][["height"]] <- h


  as_json <- attr(as_list, "TOJSON_FUNC")(as_list)
  as_json <- paste0("{ \"x\":", as_json, ",\"evals\": [],\"jsHooks\": []}")

  print(plotjs)
  write(as_json, plotjs)


  write(as_json, plotjs)

  # storing for Report Generation
  mbSetObj$analSet$heatmap <- data1
  mbSetObj$analSet$heatmap.dist <- smplDist
  mbSetObj$analSet$heatmap.clust <- clstDist
  mbSetObj$analSet$heat.taxalvl <- taxrank
  mbSetObj$imgSet$heatmap_int <- plotwidget
  save(p, file=plotwidget);

  #pstatic <- CreateStaticHeatmap(data1sc, fzAnno, colors, nrows, x_start, y_start, x_spacing, annotation, sz, bf, showColnm, showRownm, doclust, smplDist, clstDist, fzCol, fzRow)

  #Cairo::Cairo(file = paste0(plotNm, ".png"), unit="px", dpi=72, width=w, height=h, type="png");    
  #print(pstatic)
  #dev.off()

  return(.set.mbSetObj(mbSetObj))
}



#########################
### Utility Functions ###
#########################

#' Function to get color palette for graphics.
#' @description This function is called to create a color palette
#' based on the number of groups. It returns a vector of color
#' hex codes based on the number of groups.
#' @param mbSetObj Input the name of the mbSetObj.
#' @param grayscale Logical, default set to F.
#' @author Jeff Xia \email{jeff.xia@mcgill.ca}
#' McGill University, Canada
#' License: GNU GPL (>= 2)
#' @export
GetColorSchema <- function(mbSetObj, grayscale = F) {
  mbSetObj <- .get.mbSetObj(mbSetObj)
  load_phyloseq();

  # test if total group number is over 9
  claslbl <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]])
  grp.num <- length(levels(claslbl))

  if (grayscale) {
    dist.cols <- grDevices::colorRampPalette(c("grey90", "grey30"))(grp.num)
    lvs <- levels(claslbl)
    colors <- vector(mode = "character", length = length(claslbl))

    for (i in 1:length(lvs)) {
      colors[mbSetObj$analSet$cls == lvs[i]] <- dist.cols[i]
    }
  } else if (grp.num > 9) {
    pal12 <- c(
      "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
      "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
      "#FFFF99", "#B15928"
    )
    dist.cols <- grDevices::colorRampPalette(pal12)(grp.num)
    lvs <- levels(claslbl)
    colors <- vector(mode = "character", length = length(mbSetObj$analSet$cls))

    for (i in 1:length(lvs)) {
      colors[claslbl == lvs[i]] <- dist.cols[i]
    }
  } else {
    colors <- as.numeric(claslbl) + 1
    if (exists("colVec")) {
      if (any(colVec == "#NA")) {
        cols <- vector(mode = "character", length = length(claslbl))
        clsVec <- as.character(claslbl)
        grpnms <- names(colVec)

        for (i in 1:length(grpnms)) {
          cols[clsVec == grpnms[i]] <- colVec[i]
        }
        colors <- cols
      }
    }
  }
  return(colors)
}

CleanTaxaNames <- function(mbSetObj, names) {
  mbSetObj <- .get.mbSetObj(mbSetObj)
  # first get taxa type
  type <- mbSetObj$dataSet$taxa_type

  if (type == "QIIME") {
    new <- gsub("D.*__", "", names)
  } else if (type == "Greengenes") {
    new <- gsub(".*__", "", names)
  } else {
    new <- names
  }
  return(new)
}

GetColorSchemaFromFactor <- function(my.grps) {
  # test if total group number is over 9
  my.grps <- as.factor(my.grps)
  grp.num <- length(levels(my.grps))

  # if(grp.num > 9){
  pal12 <- c(
    "#A6CEE3", "#1F78B4", "#B2DF8A", "#33A02C", "#FB9A99",
    "#E31A1C", "#FDBF6F", "#FF7F00", "#CAB2D6", "#6A3D9A",
    "#FFFF99", "#B15928"
  )
  dist.cols <- colorRampPalette(pal12)(grp.num)
  lvs <- levels(my.grps)
  colors <- vector(mode = "character", length = length(my.grps))
  for (i in 1:length(lvs)) {
    colors[my.grps == lvs[i]] <- dist.cols[i]
  }
  # }else{
  #  colors <- as.numeric(my.grps)+1;
  # }
  return(colors)
}


PlotCovariateMap<- function(mbSetObj, thresh = "0.05", theme="default", imgName="NA", format="png", dpi=72, interactive=T){

  thresh <- as.numeric(thresh);
  mbSetObj <- .get.mbSetObj(mbSetObj); 
  both.mat <- mbSetObj$analSet$cov.mat;
  both.mat <- both.mat[order(-both.mat[,"pval.adj"]),];
  logp_val <- -log10(thresh);
  load_ggplot();
  library(ggrepel);
  topFeature <- 5;

  if (nrow(both.mat) < topFeature) {
    topFeature <- nrow(both.mat)
  }
  if (theme == "default") {
    p <- ggplot(both.mat, mapping = aes(x = fdr.no, y = fdr.adj, label = Row.names)) +
      geom_rect(
        mapping = aes(
          xmin = logp_val, xmax = Inf,
          ymin = logp_val, ymax = Inf
        ),
        fill = "#6699CC"
      ) +
      geom_rect(
        mapping = aes(
          xmin = -Inf, xmax = logp_val,
          ymin = -Inf, ymax = logp_val
        ),
        fill = "grey"
      ) +
      geom_rect(
        mapping = aes(
          xmin = logp_val, xmax = Inf,
          ymin = -Inf, ymax = logp_val
        ),
        fill = "#E2808A"
      ) +
      geom_rect(
        mapping = aes(
          xmin = -Inf, xmax = logp_val,
          ymin = logp_val, ymax = Inf
        ),
        fill = "#94C973"
      ) +
      guides(size = "none") +
      geom_point(aes(size = pval.adj), alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
      xlab("-log10(Adj. P-value): no covariate adjustment") +
      ylab("-log10(Adj. P-value): covariate adjustment") +
      geom_text_repel(
        data = both.mat[c(1:topFeature), ],
        aes(x = pval.no, y = pval.adj, label = Row.names)
      ) +
      theme_bw()
  } else {
    p <- ggplot(both.mat, mapping = aes(x = fdr.no, y = fdr.adj, label = Row.names)) +
      guides(size = "none") +
      geom_point(aes(size = pval.adj), alpha = 0.5) +
      geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "red", size = 1) +
      geom_vline(xintercept = logp_val) +
      geom_hline(yintercept = logp_val) +
      xlab("-log10(Adj. P-value): no covariate adjustment") +
      ylab("-log10(Adj. P-value): covariate adjustment") +
      geom_text_repel(
        data = both.mat[c(1:topFeature), ],
        aes(x = pval.no, y = pval.adj, label = Row.names)
      )
  }
  fileName <- paste0(imgName, ".", format);
  mbSetObj$imgSet$covAdj <- fileName;

  width <- 8;
  height <- 8.18;

    Cairo::Cairo(file = fileName, unit="in", dpi=dpi, width=width, height=height, type=format);    
    print(p)
    dev.off()

  if(interactive){
    library(plotly);
    ggp_build <- layout(ggplotly(p,width = 800, height = 600, tooltip = c("text")), autosize = FALSE, margin = mbSetObj$imgSet$margin.config)
    .set.mbSetObj(mbSetObj)
    return(ggp_build);
  }else{
    return(.set.mbSetObj(mbSetObj));
  }
}

CreateStaticHeatmap <- function(data1sc, fzAnno, colors, nrows, x_start, y_start, x_spacing, annotation, sz, bf, showColnm, showRownm, doclust, smplDist, clstDist, fzCol, fzRow) {
  library(pheatmap);
  # Note: In pheatmap, annotations are usually a data frame where each column is a different annotation
  # You might need to adjust this part based on your actual data structure for annotations
  
  # Prepare the color mapping
  color_breaks <- seq(min(data1sc), max(data1sc), length.out = length(colors) + 1)
  color_mapping <- colorRampPalette(colors)(length(colors))

  # Prepare clustering if needed
  clustering_distance_rows <- if(doclust == "T") smplDist else "none"
  clustering_distance_cols <- if(doclust == "T") clstDist else "none"
  
  # Create the heatmap
  p <- pheatmap(data1sc,
           color = color_mapping,
           breaks = color_breaks,
           cluster_rows = doclust == "T",
           cluster_cols = doclust == "T",
           clustering_distance_rows = clustering_distance_rows,
           clustering_distance_cols = clustering_distance_cols,
           fontsize = fzAnno, # Set the font size for annotations
           fontsize_row = fzRow,
           fontsize_col = fzCol,
           show_rownames = F,
           show_colnames = T,
           annotation = annotation
  )

  return(p)
}
