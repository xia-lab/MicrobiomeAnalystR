##################################################
## R script for MicrobiomeAnalyst
## Description: Compute json file for 3D scatter plot
## Author: Guangyan Zhou, guangyan.zhou@mail.mcgill.ca
##################################################

#dataName is used for onedata dataset
my.json.scatter.meta <- function(mbSetObj=NA, filenm="abc",selMeta=""){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  paramSet <- readSet(paramSet, "paramSet");
  mdata.all <- paramSet$mdata.all;
  anal.type <- mbSetObj$module.type;
  
  data <- qs::qread("merged.data.qs");
  #data <- subsetPhyloseqByDataset(mbSetObj, data);
  
  meta <- sample_data(data);
  
  Sys.setenv(RGL_USE_NULL = TRUE);
  require(rgl);
  require(igraph);
  pos.xyz <-qs::qread("score_pos_xyz.qs");
  nodes <- vector(mode="list");
  names <- c(rownames(pos.xyz));
  metadf <- meta;
  
  
  a<-list();
  a$objects <- "NA";
  meshes<-"NA";
  
  col <- vector();
  
  
  # can be selected meta as well if = dataSet$sel.meta
  meta.vec <- as.vector(metadf[,selMeta][[1]]);
  meta.vec.num <- as.integer(as.factor(metadf[,selMeta][[1]]));
  dataset.vec <- metadf$dataset;
  
  col.s <- gg_color_hue_scatter(length(unique(meta.vec)), "green");
  for(i in 1:length(meta.vec.num)){
    col[i] <- col.s[meta.vec.num[i]];
  }
  color = col;
  nodeSize = 18;
  if(length(names)>200){
    nodeSize = 16;
  }
  seriesList <- list();
  for(j in 1:length(unique(dataset.vec))){
    nodes <- list();
    k = 1;
    for(i in 1:length(names)){
      if(unique(dataset.vec)[j] == dataset.vec[i]){
        nodes[[k]] <- c(
          #id=
          names[i],
          #size=
          nodeSize,
          #meta=
          meta.vec[i],
          #cluster=
          meta.vec.num[i],
          #fx = 
          unname(pos.xyz[i,1]),
          #fy = 
          unname(pos.xyz[i,2]),
          #fz = 
          unname(pos.xyz[i,3]),
          #color=
          color[i],
          dataset.vec[i]
        );
        k = k+1;
      }
      
    }
    
    seriesList[[j]] <- nodes;
    
  }
  
  edge.mat = "NA";
  modules = "NA";
  # save node table
  ellipse ="NA";  
  require(RJSONIO);
  
  netData <- list( seriesList=seriesList,colArr=col.s,metaVecUniq=unique(meta.vec.num), metaVec = meta.vec, edges=edge.mat, modules=modules, objects=a$objects, ellipse=meshes, meta=metadf, loading="NA", reductionOpt="pca" ,omicstype=c("rna.b"));
  
  netData[["misc"]] <- "";
  paramSet$partialToBeSaved <- c(paramSet$partialToBeSaved, c(filenm));
  paramSet$jsonNms["scatter3d"] <- filenm;
  
  saveSet(paramSet, "paramSet");
  
  sink(filenm);
  cat(toJSON(netData));
  sink();
  
  return(1);
}
