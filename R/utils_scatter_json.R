##################################################
## R script for MicrobiomeAnalyst
## Description: Compute 3D scatter plot from dimension reduction results
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## J. Xia, jeff.xia@mcgill.ca
###################################################

my.json.scatter <- function(filenm,analysisVar){
  omicstype.vec <- c("microbiome","metabolomics");

  library(RJSONIO)
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
  }
print(micDataType)

   metdat <- current.proc$met$data.proc
  if(micDataType=="ko"){
   sig.mic <- current.proc$mic$res_deAnal[which(current.proc$mic$res_deAnal$P_value<current.proc$mic$plvl),]
   sig.mic <- list(OTU=sig.mic)
  }else{
   sig.mic <- lapply(phyloseq_objs$res_deAnal,function(x) x[which(x$P_value<current.proc$mic$plvl),])
  } 
  sig.met <- current.proc$met$res_deAnal
  sig.met <- sig.met[which(sig.met$P_value<current.proc$met$plvl),]
  sig.met$ids <- rownames(sig.met)
  sig.mats <-vector("list",length=length(sig.mic))
  
  names(sig.mats)<-names(sig.mic)
  for(s in 1:length(sig.mats)){
      sig.mic[[s]]$ids<-rownames(sig.mic[[s]])
     sig.mats[[s]][["microbiome"]] <- sig.mic[[s]]      
     sig.mats[[s]][["metabolomics"]] <- sig.met
  }

  seeds <- lapply(sig.mic,function(x) c(rownames(x),rownames(sig.met)))
  #meta <- meta ### for other methods not import here
  
  
  Sys.setenv(RGL_USE_NULL = TRUE)
  library(rgl)
  library(igraph)
  
  if(reductionOptGlobal == "procrustes"){
    if(!exists("procrustes.res")){
      procrustes.res <- qs::qread("procrustes.res.qs")
    }
    pos.xyz.all =  procrustes.res$pos.xyz
    metadf = procrustes.res$newmeta
  }else if(reductionOptGlobal == "diablo"){
    if(!exists("diablo.res")){
      diablo.res <- qs::qread("diablo.res.qs")
    }
    if(!exists("combined.res")){
      combined.res <- qs::qread("combined.res.qs")
    }
    pos.xyz.all =  diablo.res$pos.xyz
    metadf = diablo.res$newmeta
  }
 
  if("omics" %in% colnames(metadf)){
    sel.meta = "omics"
  }
  netData <- vector("list",length=length(pos.xyz.all))
  names(netData) <- names(pos.xyz.all)
  for(tax in names(pos.xyz.all)){
    pos.xyz <-   pos.xyz.all[[tax]]
    nodes <- vector(mode="list");
    names <-  make.unique(as.character(rownames(pos.xyz)))
    
    a=list();
    a$objects = "NA";
    meshes="NA"
    col = vector();
    
    # can be selected meta as well if = dataSet$sel.meta
    if(!is.null(analysisVar)){
      meta.vec = as.vector(metadf[,analysisVar])
      meta.vec.num = as.integer(as.factor(metadf[,analysisVar]))
    }else{
      meta.vec = as.vector(metadf[,1])
      meta.vec.num = as.integer(as.factor(metadf[,1]))
    }
    col.s <- gg_color_hue(length(unique(meta.vec)), "green")
    for(i in 1:length(meta.vec.num)){
      col[i] = col.s[meta.vec.num[i]];
    }
    color = col
    nodeSize = 18;
    if(length(names)>200){
      nodeSize = 16;
    }
    
    
    for(i in 1:length(names)){
      nodes[[i]] <- list(
        id=names[i],
        label=names[i],
        size=nodeSize,
        meta=meta.vec[i],
        cluster=meta.vec.num[i],
        fx = unname(pos.xyz[i,1])*1000,
        fy = unname(pos.xyz[i,2])*1000,
        fz = unname(pos.xyz[i,3])*1000,
        colorb=color[i],
        colorw=color[i],
        topocolb=color[i],
        topocolw=color[i],
        expcolb=color[i],
        expcolw=color[i],
        attributes=list(
          expr = 1,
          degree=1,
          between=1
        )
      );
    }
    
    if(reductionOptGlobal == "procrustes"){
      pos.xyz.length = nrow(pos.xyz)
      edge.mat <- cbind(id=c(1:(pos.xyz.length)), source=names[c(1:(pos.xyz.length/2))], target=names[c(((pos.xyz.length/2)+1):pos.xyz.length) ], opacity = 0);
      procrustes.res$pos.xyz[[tax]] = pos.xyz
      modules = "NA"
      # save node table
      ellipse ="NA"
      netData[[tax]] <- list(nodes=nodes, edges=edge.mat, modules=modules, objects=a$objects,
                             ellipse=meshes, meta=metadf, loading="NA", reductionOpt=reductionOptGlobal, 
                             misc=procrustes.res$misc[[tax]], omicstype=c("microbiome","metabolomics"), sigMat=sig.mats[[tax]]);
      
      pca.scatter <- generatePCAscatter(phyloseq_objs$count_tables[[tax]],
                                        metdat,tax)
      
      procrustes.res$misc[[tax]]$pct2 <- pca.scatter$pct2
      
      netData[[tax]][["misc"]] <- procrustes.res$misc[[tax]]

     qs::qsave(procrustes.res,"procrustes.res.qs")
      
    }else if(reductionOptGlobal == "diablo"){
    
      pos.xyz2 <-   diablo.res$pos.xyz2[[tax]]
      names <- rownames(pos.xyz2)
      nodes_samples2 <- vector(mode="list");
      for(i in 1:length(names)){
        nodes_samples2[[i]] <- list(
          id=names[i],
          label=names[i],
          size=nodeSize,
          meta=meta.vec[i],
          cluster=1,
          fx = unname(pos.xyz2[i,1])*1000,
          fy = unname(pos.xyz2[i,2])*1000,
          fz = unname(pos.xyz2[i,3])*1000,
          colorb=color[i],
          colorw=color[i],
          topocolb=color[i],
          topocolw=color[i],
          expcolb=color[i],
          expcolw=color[i],
          attributes=list(
            expr = 1,
            degree=1,
            between=1
          )
        );
      }
      
      edge.mat = "";
      modules = "NA"
      ellipse ="NA"
    
      diablo.res$pos.xyz[[tax]] = pos.xyz;
          
      loading.data = diablo.res$loading.pos.xyz[[tax]];
      cluster = diablo.res$loadingCluster[[tax]];
      aLoading=list();
      aLoading$objects = "NA";
      
      names = diablo.res$loading.enrich[[tax]]
      ids = diablo.res$loadingNames[[tax]]
      rownames(loading.data) = names
      de = combined.res$comp.res[[tax]]
      de = de[which(rownames(de) %in% ids),]
      de[de == "NaN"] = 1
      pv = as.numeric(de[,"P_value"])
      pv_no_zero = pv[pv != 0]
      minval = min(pv_no_zero)
      pv[pv == 0] = minval/2
      pvals <<- -log10(pv);
      type.vec <- pvals;
      if(exists("comp.res.inx",combined.res)){
        for(i in 1:length(unique(combined.res$comp.res.inx[[tax]]))){
          inx = combined.res$comp.res.inx[[tax]] == i
          type.vec[inx] <- omicstype.vec[i]
        }
      }
      colors<- ComputeColorGradient(pvals,  F);
      colorb <- colors;
      sizes <- as.numeric(rescale2NewRange(-log10(pv), 15, 25));
      nodes2 <- vector(mode="list");
      
      #loading.data = loading.data[which(rownames(loading.data) %in% as.character(ids)),];
      
      seed.inx <- names %in% unique(seeds[[tax]]);
      seed_arr <- rep("notSeed",length(names));
      seed_arr[seed.inx] <- "seed";
      
      for(i in 1:length(pvals)){
        nodes2[[i]] <- list(
          id=ids[i],
          label=names[i],
          size=sizes[i],
          cluster=1,
          omicstype=type.vec[i],
          fx = unname(loading.data[i,1])*1000,
          fy = unname(loading.data[i,2])*1000,
          fz = unname(loading.data[i,3])*1000,
          seedArr = seed_arr[i],
          colorb=colorb[i],
          colorw=colorb[i],
          topocolb="#ffa500",
          topocolw="#ffa500",
          expcolb="#ffa500",
          expcolw="#ffa500",
          attributes=list(
            expr = 1,
            degree=1,
            between=1
          )
        );
      }
      netData[[tax]] <- list(omicstype=omicstype.vec, nodes=nodes, edges=edge.mat, modules=modules, objects=a$objects, ellipse=meshes, meta=metadf, loading=nodes2, reductionOpt=reductionOptGlobal , objectsLoading=aLoading$objects, sigMat=sig.mats[[tax]]);
   
      type <- omicstype.vec[2]
      netData[[tax]][[ type]] <- nodes_samples2;
      
      
      pca.scatter <- generatePCAscatter(phyloseq_objs$count_tables[[tax]],
                                        metdat,tax)
    
      diablo.res$pca.scatter[[tax]] <- pca.scatter
      for(i in 1:length(omicstype.vec)){
        pos<-pca.scatter[[paste0("pca_", omicstype.vec[i])]]$score
        
        pca_nodes <- nodes[c(1:nrow(pos))]
        for(j in 1:nrow(pos)){
          pca_nodes[[j]][["id"]] <-rownames(pos)[j]
          pca_nodes[[j]][["label"]] <-rownames(pos)[j]
          pca_nodes[[j]][["fx"]] <-pos[j,1]
          pca_nodes[[j]][["fy"]] <-pos[j,2]
          pca_nodes[[j]][["fz"]] <-pos[j,3]
        }
        nm <- paste0("pca_", omicstype.vec[[i]])
        netData[[tax]][[nm]] <- pca_nodes;
        
        loading.data<-pca.scatter[[paste0("pca_", omicstype.vec[i])]]$loading
        #print(head(loading.data))
        loadingNames <- rownames(loading.data);
        enrich_ids <- intersect(combined.res$enrich_ids[[tax]],loadingNames)
        loading.enrich = enrich_ids[order(match(enrich_ids, loadingNames))]
        ids = loadingNames
        names = loading.enrich
        
        
        pca_loading <- nodes2[c(1:nrow(loading.data))];
        count = 1
        for(k in 1:length(nodes2)){
          if(nodes2[[k]][["id"]] %in% rownames(loading.data) || nodes2[[k]][["label"]] %in% rownames(loading.data)){
            pca_loading[[count]] <- nodes2[[k]]
            pca_loading[[count]][["id"]]=ids[k]
            pca_loading[[count]][["label"]]=names[k]
            pca_loading[[count]][["omicstype"]]=type.vec[k]
            inx = which(rownames(loading.data) == nodes2[[k]][["id"]]);
            if(length(inx) == 0){
              inx = which(rownames(loading.data) == nodes2[[k]][["label"]]);
            }
            pca_loading[[count]][["fx"]] = loading.data[inx,1]
            pca_loading[[count]][["fx"]] = loading.data[inx,2]
            pca_loading[[count]][["fx"]] = loading.data[inx,3]
            count = count +1;
            
          }
        }
        nm <- paste0("pca_", omicstype.vec[[i]], "_loading")
        netData[[tax]][[nm]] <- pca_loading;
      }
      
      diablo.res$misc[[tax]]$pct2 <- c( diablo.res$misc[[tax]]$pct2, pca.scatter$pct2)
      
      netData[[tax]][["misc"]] <- diablo.res$misc[[tax]]

      
       qs::qsave(diablo.res,"diablo.res.qs")
    }

  }
 
  jsonNms_scatter <<- filenm;
  
  sink(filenm);
  cat(toJSON(netData));
  sink();
  return(1)

}

generatePCAscatter <- function(micdat,metdat,tax){
  
  fig.list <- list()
  pca.list<- list()
  pct <- list();
  
  for(i in 1:2){
    if(i==1){

      x <- data.matrix(micdat);
      type= "microbiome"
    }else{
     
      x <- data.matrix(metdat);
      type="metabolomics"
    }
    pca <- prcomp(t(na.omit(x)), center=T, scale=T);
    imp.pca<-summary(pca)$importance;
    nm <- paste0("pca_", type);
 
  pos.xyz = data.frame(x=pca$x[,1], y=pca$x[,2], z=pca$x[,3]);
  pos.xyz <- unitAutoScale(pos.xyz);

 loading.pos.xyz = data.frame(pca$rotation[,c(1:3)]);
  loadingNames <- rownames(loading.pos.xyz);
  loading <- unitAutoScale(loading.pos.xyz);
  rownames(loading) <- loadingNames;
  

  pca.list[[nm]][["score"]] <- pos.xyz * 1000;
  pca.list[[nm]][["loading"]] <- loading* 1000;    
  
  pct[[nm]] <- unname(round(imp.pca[2,],3))[c(1:3)]*100;

  }
pca.list$pct2 <- pct;
  return(pca.list)
}

UpdateSampleBasedOnLoading<-function(filenm, gene.id, omicstype){
  if(omicstype != "NA"){
  sel.nms <- names(mdata.all)[mdata.all==1];
    for(i in 1:length(sel.nms)){
      dat = readRDS(sel.nms[i])
      if(dat$type == omicstype){
        dataSet <- dat;
      }
    }
  }else{
    dataSet <- .get.rdt.set();
  }
  
  inx <- which(dataSet$enrich_ids == gene.id)
  id <- unname(dataSet$enrich_ids[inx])
  vec = as.vector(dataSet$data.proc[rownames(dataSet$data.proc) == gene.id,])
  colors<- ComputeColorGradient(as.numeric(vec), "black", F, F);
  sink(filenm);
  cat(toJSON(colors));
  sink();
}

ComputeEncasing <- function(filenm, type, names.vec, level=0.95, omics="NA",taxalvl="OTU"){
  Sys.setenv(RGL_USE_NULL = TRUE)
  level <- as.numeric(level)
  names = strsplit(names.vec, "; ")[[1]]
  
  if(reductionOptGlobal %in% c("diablo", "spls") || omics != "NA"){
    if(!exists("diablo.res")){
      diablo.res <- qs::qread("diablo.res.qs")
    }
    if(grepl("pca_", omics, fixed=TRUE)){
      pos.xyz<-diablo.res$pca.scatter[[taxalvl]]$omics$score/1000
    }else{
      
      if(omics == "microbiome"){
        pos.xyz = diablo.res$pos.xyz[[taxalvl]]
      }else{
        pos.xyz = diablo.res$pos.xyz2[[taxalvl]]
      }
    }
    
  }else{
    procrustes.res <- qs::qread("procrustes.res.qs")
    pos.xyz = procrustes.res$pos.xyz[[taxalvl]]
  }
  
  inx = rownames(pos.xyz) %in% names;
  
  coords = as.matrix(pos.xyz[inx,c(1:3)])
  mesh = list()
  if(type == "alpha"){
    library(alphashape3d)
    library(rgl)
    sh=ashape3d(coords, 1.0, pert = FALSE, eps = 1e-09);
    mesh[[1]] = as.mesh3d(sh, triangles=T);
  }else if(type == "ellipse"){
    library(rgl);
    pos=cov(coords, y = NULL, use = "everything");
    mesh[[1]] = ellipse3d(x=as.matrix(pos), level=level);
  }else{
    library(ks);
    res=kde(coords);
    r = plot(res, cont=level*100, display="rgl");
    sc = scene3d();
    mesh = sc$objects;
  }
  library(RJSONIO);
  sink(filenm);
  cat(toJSON(mesh));
  sink();
  return(filenm);
}
