##################################################
## R script for MicrobiomeAnalyst
## Description: Functions to create PCoA3D plot
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
##################################################

my.pcoa.3d <- function(mbSetObj, ordMeth, distName, taxrank, colopt, variable, taxa, alphaopt, jsonNm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  suppressMessages(library(vegan));
  
  variable <<- variable;

  if(!exists("phyloseq_objs")){

 phyloseq_objs <- readDataQs("phyloseq_objs.qs",mbSetObj$module.type,mbSetObj$dataSet$name)
  }
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
      #converting NA values to unassigned
      nm[is.na(nm)] <- "Not_Assigned";
      data1 <- as.matrix(otu_table(data));
      rownames(data1) <- nm;
      #all NA club together
      data1 <- rowsum(as.matrix(data1), rownames(data1));
      feat_data <- data1[taxa,];
    }
    sample_data(data)$taxa <- feat_data;
    indx <- which(colnames(sample_data(data))=="taxa");
    colnames(sample_data(data))[indx] <- taxa;
  }else if(colopt=="alphadiv"){
    data1 <- mbSetObj$dataSet$proc.phyobj;
    box <- plot_richness(data1, measures = alphaopt);
    alphaboxdata <- box$data;
    sam_nm <- sample_names(data);
    alphaboxdata <- alphaboxdata[alphaboxdata$samples %in% sam_nm,];
    alphaval <- alphaboxdata$value;
    sample_data(data)$alphaopt <- alphaval;
    indx <- which(colnames(sample_data(data))=="alphaopt");
    colnames(sample_data(data))[indx]<-alphaopt;
  }else{
    data<-data;
  }

  if(distName=="wunifrac"){
    pg_tree <- qs::qread("tree.qs");
    pg_tb <- tax_table(data);
    pg_ot <- otu_table(data);
    pg_sd <- sample_data(data);
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
    data <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree);
    
    if(!is.rooted(phy_tree(data))){
      pick_new_outgroup <- function(tree.unrooted){
        treeDT <- cbind(cbind(data.table(tree.unrooted$edge),data.table(length = tree.unrooted$edge.length))[1:Ntip(tree.unrooted)],
                        data.table(id = tree.unrooted$tip.label));
        new.outgroup <- treeDT[which.max(treeDT$length), ]$id
        return(new.outgroup);
      }
      new.outgroup <- pick_new_outgroup(phy_tree(data));
      phy_tree(data) <- ape::root(phy_tree(data),
                                  outgroup = new.outgroup,
                                  resolve.root=TRUE)
    }
    if(ordMeth=="PCA"){
             GP.ord  <- prcomp(t(data@otu_table@.Data), center=TRUE, scale=F)
      }else{
       
    GP.ord <-ordinate(data,ordMeth,"unifrac",weighted=TRUE);
      }
  } else if (distName=="unifrac"){
    pg_tree <- qs::qread("tree.qs");
    pg_tb <- tax_table(data);
    pg_ot <- otu_table(data);
    pg_sd <- sample_data(data);
    pg_tree <- prune_taxa(taxa_names(pg_ot), pg_tree);
    data <- merge_phyloseq(pg_tb, pg_ot, pg_sd, pg_tree);
    
    if(!is.rooted(phy_tree(data))){
      pick_new_outgroup <- function(tree.unrooted){
        treeDT <- cbind(cbind(data.table(tree.unrooted$edge),data.table(length = tree.unrooted$edge.length))[1:Ntip(tree.unrooted)],
                        data.table(id = tree.unrooted$tip.label));
        new.outgroup <- treeDT[which.max(treeDT$length), ]$id
        return(new.outgroup);
      }
      new.outgroup <- pick_new_outgroup(phy_tree(data));
      phy_tree(data) <- ape::root(phy_tree(data),
                                  outgroup = new.outgroup,
                                  resolve.root=TRUE)
    }
   
 if(ordMeth=="PCA"){
       GP.ord  <- prcomp(t(data@otu_table@.Data), center=TRUE, scale=F)
      }else{
       
       GP.ord <-ordinate(data,ordMeth,"unifrac",weighted=FALSE);
      }
  }else{
 
    
      if(ordMeth=="PCA"){
             GP.ord  <- prcomp(t(data@otu_table@.Data), center=TRUE, scale=F)
      }else{
       
       GP.ord <- ordinate(data,ordMeth,distName);
      }
  }
 
  # obtain variance explained

  if(ordMeth=="PCA"){
     sum.pca <-  GP.ord;
     sum.pca$vectors <- sum.pca$x;
  imp.pca <- summary(GP.ord)$importance;
  std.pca <- imp.pca[1,]; # standard devietation
  var.pca <- imp.pca[2,]; # variance explained by each PC
  cum.pca <- imp.pca[3,]; # cummulated variance explained
  

  }else{
    sum.pca <- GP.ord;
 imp.pca <- sum.pca$values;
  std.pca <- imp.pca[1,]; # eigen values
  var.pca <- imp.pca[,2]; # variance explained by each PC
  cum.pca <- imp.pca[5,]; # cummulated variance explained
   }
 sum.pca <- append(sum.pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca));
  
  pca3d <- list();
  
  if(ordMeth=="NMDS"){
    pca3d$score$axis <- paste("NMDS", 1:3 , sep="");
    coord<-sum.pca$points;
    fast.write(signif(coord,5), file="pcoa_score.csv");
    list2 <- rep(0,nrow(coord));
    coord <- cbind(coord, list2);
    coords <- data.frame(t(signif(coord[,1:3], 5)),check.names=FALSE);
  }else{
    pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(sum.pca$variance[1:3], 3), "%)", sep="");
    coords <- data.frame(t(signif(sum.pca$vectors[,1:3], 5)),check.names=FALSE);
    fast.write(signif(sum.pca$vectors,5), file="pcoa_score.csv");
  }
  
  colnames(coords) <- NULL;
  pca3d$score$xyz <- coords;
  pca3d$score$name <- sample_names(mbSetObj$dataSet$norm.phyobj);
  col.type <- "factor";
  
  if(colopt=="taxa"){
    cls <- sample_data(data)[[taxa]];
    col.type <- "gradient"
    cols <- ComputeColorGradient(cls);
  }else if(colopt=="alphadiv") {
    cls <- sample_data(data)[[alphaopt]];
    col.type <- "gradient";
    cols <- ComputeColorGradient(cls);
  }else{
    cls <- sample_data(mbSetObj$dataSet$norm.phyobj)[[variable]];
    # now set color for each group
    cols <- unique(as.numeric(as.factor(cls))) + 1;
  }
  cls <- as.factor(cls);
  pca3d$score$type <- col.type;
  pca3d$score$facA <- cls;
  rgbcols <- col2rgb(cols);
  cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
  pca3d$score$colors <- cols;

  #json.obj <- rjson::toJSON(pca3d);
  #sink(jsonNm);
  #cat(json.obj);
  #sink();

  mbSetObj$analSet$pca<-pca3d;

  if(!exists("my.json.scatter")){
    .load.scripts.on.demand("utils_scatter3d.Rc");    
  }

  .set.mbSetObj(mbSetObj)
  my.json.scatter(mbSetObj, jsonNm, F);

  return(.set.mbSetObj(mbSetObj));
}
