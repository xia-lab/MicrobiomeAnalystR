##################################################
## R script for MicrobiomeAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## Y. Lu   
## J. Xia, jeff.xia@mcgill.ca
###################################################

my.reduce.dimension <- function(mbSetObj, reductionOpt= "procrustes", method="globalscore", dimn=10,analysisVar){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(method == ""){
    method="globalscore"
  }

  dimn = as.numeric(dimn);
  if(dimn == 2){
    dimn = 3;
  }
 
  d.list = vector("list",length = 2)
  names(d.list) <- c("mic","met")
  omics.type = vector();
  featureNms <- vector();
  combined.res <- list()
  
  omics.type <- c("microbiome","metabolomics")
  
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
  }
  
  d.list[["mic"]] = list()
  d.list[["mic"]][["data.proc"]] = phyloseq_objs$count_tables
  d.list[["mic"]][["comp.res"]] = lapply(phyloseq_objs$res_deAnal, function(x){names(x)[1] ="T.Stats"; return(x[,c(3,4,1)])})
  d.list[["mic"]][["enrich.nms"]] = lapply(phyloseq_objs$res_deAnal ,function(x) rownames(x))
  d.list[["mic"]][["meta"]] = data.frame(mbSetObj$dataSet$sample_data)
  
  d.list[["met"]] = list()
  d.list[["met"]][["data.proc"]] =   current.proc$met$data.proc
  d.list[["met"]][["comp.res"]] =   current.proc$met$res_deAnal[,c(1:3)] #comp.res
  d.list[["met"]][["enrich.nms"]] = rownames(current.proc$met$res_deAnal)
  d.list[["met"]][["meta"]] = data.frame(mbSetObj$dataSet$sample_data)
  
  
  newmeta <- rbind(  d.list[["mic"]][["meta"]],d.list[["met"]][["meta"]])
  comp.res1 = lapply( d.list[["mic"]][["comp.res"]],function(x) rbind(x, d.list[["met"]][["comp.res"]]) )
  
  enrich.nms1 = lapply(d.list[["mic"]][["enrich.nms"]],function(x) c(x, d.list[["met"]][["enrich.nms"]]) )
  comp.res.inx1 =lapply(d.list[["mic"]][["comp.res"]],function(x) c(rep(1,nrow(x)), rep(2,nrow( d.list[["met"]][["comp.res"]])) ));
  featureNms <- lapply( d.list[["mic"]][["data.proc"]] ,function(x) c(rownames(x), rownames(d.list[["met"]][["data.proc"]]) ));
  
  
  
  combined.res$comp.res = comp.res1
  combined.res$enrich_ids = enrich.nms1
  combined.res$comp.res.inx = comp.res.inx1
  combined.res$meta = newmeta
  
  if(reductionOpt == "diablo"){
    library(mixOmics)
    if(is.null(analysisVar)){
    Y <- current.proc$meta_para$sample_data[,analysisVar]
    }else{
     Y <- current.proc$meta_para$sample_data[,1]
    }
    diablo.res <- list()
    dats <- vector("list",length=length(d.list$mic$data.proc))
    names(dats) <- names(d.list$mic$data.proc)
    
    res <- pos.xyz <- pos.xyz2 <- misc <- loading.pos.xyz <- loadingNames <- vector("list",length=length(d.list$mic$data.proc))
    names(res) <- names(pos.xyz) <- names(pos.xyz2) <- names(misc) <- 
      names(loading.pos.xyz) <- names(loadingNames) <- names(d.list$mic$data.proc)

    for(l in 1:length(dats)){
     dats[[l]][["mic"]] <- d.list$mic$data.proc[[l]]
     dats[[l]][["met"]] <- d.list$met$data.proc
     dats[[l]]= lapply(dats[[l]], function(x){
       x <- data.frame(x, stringsAsFactors = FALSE);
       x <- t(x);
     })
     ncomp = min(ncol(dats[[l]]$mic),dimn)
     design = matrix(1, ncol = length(dats[[l]]), nrow = length(dats[[l]]), 
                     dimnames = list(names(dats[[l]]), names(dats[[l]])))
     res[[l]] = mixOmics:::block.splsda(X = dats[[l]], Y = Y, ncomp = ncomp, design = design)
     pos.xyz[[l]] <- res[[l]]$variates[[1]]
     pos.xyz2[[l]] <- res[[l]]$variates[[2]]
     
     for(i in 1:length(res[[l]]$loadings)){
       pos = as.data.frame(res[[l]]$loadings[[i]])
       rn <- rownames(res[[l]]$loadings[[i]])
       pos <- unitAutoScale(pos);
       res[[l]]$loadings[[i]] <- pos
       rownames(res[[l]]$loadings[[i]]) <- rn
     }
     
     loading.pos.xyz[[l]] <- rbind(res[[l]]$loadings[[1]], res[[l]]$loadings[[2]])
     rownames(loading.pos.xyz[[l]]) = c(rownames(d.list$mic$data.proc[[l]]), rownames(d.list$met$data.proc))
     loadingNames[[l]]=rownames(loading.pos.xyz[[l]])
     
     names = lapply(pos.xyz, function(x) rownames(x))
     newmeta = current.proc$meta_para$sample_data
      
      if("prop_expl_var" %in% names(res[[l]])){
         var.vec <- res[[l]]$prop_expl_var
       }else if("explained_variance" %in% names(res[[l]])){
         var.vec <- res[[l]]$explained_variance
       }else{
         var.vec <- 0;
       }
       misc[[l]]$pct2[["microbiome"]] = unname(signif(as.numeric(var.vec[['mic']],4)))*100
       misc[[l]]$pct2[["metabolomics"]] = unname(signif(as.numeric(var.vec[['met']],4)))*100

    }
    
  }else if(reductionOpt == "procrustes"){
    library(vegan)
    ndat1 <- lapply(d.list$mic$data.proc,function(x) decostand(t(x), method = "standardize"))
    pca.dat1 <- lapply(ndat1,function(x) rda(x))
    ndat2 <- decostand(t(d.list$met$data.proc), method = "standardize")
    pca.dat2 <- rda(ndat2)
    
    if(dimn == 2){
      choicesVec = c(1,2)
    }else{
      choicesVec = c(1,2,3)
    }
    res= lapply(pca.dat1,function(x) procrustes(x, pca.dat2, choices=choicesVec, symmetric = T, scale = TRUE))
    res2= lapply(pca.dat1,function(x) protest(X = x, Y = pca.dat2, scores = "sites", permutations = 999) )
    
    misc <- list(`Sum of Squares`=lapply(res2,function(x) x$ss),
                 Significance = lapply(res2,function(x) x$signif),
                 Correlation =lapply(res2,function(x) x$scale))
    
    
    pos.xyz = lapply(res,function(x) rbind(x$X, x$Yrot))
    names = lapply(pos.xyz,function(x) make.unique(as.character(rownames(x))))
    newmeta$omics[c(1:(length(names[[1]])/2))] = "microbiome"
    newmeta$omics[c( ((length(names[[1]])/2)+1) : (length(names[[1]])))] = "metabolomics"
    dim.res <-  misc <- vector("list",length=length(res))
    names(dim.res) <- names(misc) <- names(res)
    for(i in 1:length(dim.res)){
      dim.res[[i]] <- list(res[[i]],res2[[i]])
      misc[[i]] <- list(`Sum of Squares`=res2[[i]]$ss,
                        Significance =  res2[[i]]$signif,
                        Correlation =res2[[i]]$scale)
      
    }
    
    procrustes.res <- list(misc=misc,dim.res=dim.res)
    procrustes.res$pos.xyz = lapply(pos.xyz,function(x) unitAutoScale(x))
    procrustes.res$newmeta = newmeta
    combined.res$meta = newmeta
    qs::qsave(combined.res,"combined.res.qs")
    qs::qsave(procrustes.res,"procrustes.res.qs")
  }
  if(reductionOpt != "mcia"){ # had to rescale with segment points
    pos.xyz <- lapply(pos.xyz,function(x) as.data.frame(x));
    pos.xyz <- lapply(pos.xyz,function(x) unitAutoScale(x));
    pos.xyz <- lapply(pos.xyz, "rownames<-", names[[1]]);
  }
  if(reductionOpt %in% c("diablo", "spls")){
    names2 <- lapply(pos.xyz2, function(x) rownames(x))
    pos.xyz2 <- lapply(pos.xyz2,function(x) as.data.frame(x))
    pos.xyz2 <- lapply(pos.xyz2, function(x) unitAutoScale(x));
    pos.xyz2 <- lapply(pos.xyz2, "rownames<-", names2[[1]]);
  }
  
  if(reductionOpt != "procrustes"){
    hit.inx <- mapply(function(x,y){
      match(x,y)
    },loadingNames, combined.res$enrich_ids);
    loadingSymbols <- mapply(function(x,y){
      y[x]
    },hit.inx, combined.res$enrich_ids);
  }
  if(reductionOpt %in% c("o2pls", "mcia","rcca", "mbpca", "diablo", "rcca", "spls")){
    loading.pos.xyz = lapply(loading.pos.xyz,as.data.frame)
    loading.pos.xyz <- lapply(loading.pos.xyz,unitAutoScale);
    
    #rownames(loading.pos.xyz) = loadingNames
   
    # for(i in 1:length(names(data.list))){
    #   loadingNames <- gsub(paste0("_", names(data.list)[i]), "", loadingNames)
    # }

    diablo.res$dim.res <- res
    diablo.res$pos.xyz <- pos.xyz
    diablo.res$pos.xyz2 <- pos.xyz2
    diablo.res$misc <- misc
    diablo.res$loading.pos.xyz <- loading.pos.xyz
    diablo.res$loadingNames <- loadingNames
    diablo.res$loading.enrich = loadingSymbols
    diablo.res$newmeta=newmeta
    qs::qsave(combined.res,"combined.res.qs")
    qs::qsave(diablo.res,"diablo.res.qs")
    
    
    if(reductionOpt == "o2pls"){ #second loading plot
      loading.pos.xyz2 <- as.data.frame(loading.pos.xyz2)
      loading.pos.xyz2 <- unitAutoScale(loading.pos.xyz2);
      rownames(loading.pos.xyz2) = loadingNames2
      dataSet$loading.names2 = loadingNames2
      dataSet$loading.pos.xyz2 = loading.pos.xyz2
    }
    
  }
  reductionOptGlobal <<- reductionOpt
  
  .set.mbSetObj(mbSetObj);
  
  if(dimn == 2){
    netData <-list()
    if(reductionOpt == "mcia"){
      pos.xyz = rbind(pos.xyz, dataSet$mcia.seg.points);
    }
    netData$score=pos.xyz
    if(reductionOpt %in% loadingOpts){
      netData$loading = list(pos=loading.pos.xyz, name=loadingSymbols);
    }else{
      netData$loading = "NA";
    }
    sink(TwoDJsonName);
    cat(toJSON(netData));
    sink();
  }
  return(1)
}


