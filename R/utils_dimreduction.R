##################################################
## R script for MicrobiomeAnalyst
## Description: Compute dimension reduction
## Authors: 
## G. Zhou (guangyan.zhou@mail.mcgill.ca) 
## Y. Lu   
## J. Xia, jeff.xia@mcgill.ca
###################################################
#procrustes or diablo
my.reduce.dimension <- function(mbSetObj, reductionOpt= "procrustes", method="globalscore", dimn=10,analysisVar, diabloPar=0.2){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(method == ""){
    method="globalscore"
  }
  dimn = as.numeric(dimn);
  
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
  if(micDataType=='ko'){
    d.list[["mic"]][["comp.res"]]  = current.proc$mic$res_deAnal[,c(3,4,1)]
    names(d.list[["mic"]][["comp.res"]])[3] = "T.Stats"
    d.list[["mic"]][["comp.res"]] = list(OTU=d.list[["mic"]][["comp.res"]])
    d.list[["mic"]][["enrich.nms"]] = list(OTU=rownames(current.proc$mic$res_deAnal))
    
  }else{
    d.list[["mic"]][["comp.res"]] = lapply(phyloseq_objs$res_deAnal, function(x){names(x)[1] ="T.Stats"; return(x[,c(3,4,1)])})
    d.list[["mic"]][["enrich.nms"]] = lapply(phyloseq_objs$res_deAnal ,function(x) rownames(x))
  }
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
    # NOTE: Run DIABLO analysis in callr subprocess (mixOmics isolated)
    diablo.meta.type <- mbSetObj$dataSet$meta.types[analysisVar]
    diabloPar <- as.numeric(diabloPar); #default diabloPar was 0.2
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

      # Prepare Y and design
      if(diablo.meta.type == "disc"){
        Y <-  current.proc$meta_para$sample_data[,analysisVar];
      } else {
        meta.var <- current.proc$meta_para$sample_data[,analysisVar];
        Y <- matrix(as.numeric(as.character(meta.var)));
        rownames(Y) <- rownames(current.proc$meta_para$sample_data);
      }
      design = matrix(diabloPar, ncol = length(dats[[l]]), nrow = length(dats[[l]]),
                      dimnames = list(names(dats[[l]]), names(dats[[l]])));
      diag(design) = 0;

      # Run DIABLO in callr subprocess (mixOmics isolated)
      if (exists("diablo_analysis_isolated", mode = "function")) {
        diablo_result <- diablo_analysis_isolated(
          dats = dats[[l]],
          Y = Y,
          ncomp = ncomp,
          design = design,
          diablo.meta.type = diablo.meta.type
        )
        res[[l]] <- diablo_result$model
        pos.xyz[[l]] <- diablo_result$pos.xyz
        pos.xyz2[[l]] <- diablo_result$pos.xyz2
        loading.pos.xyz[[l]] <- diablo_result$loading.pos.xyz
        var.vec <- diablo_result$var.vec
      } else {
        # Fallback to direct call (for non-Pro)
        library(mixOmics)
        if(diablo.meta.type == "disc"){
          res[[l]] = block.splsda(X = dats[[l]], Y = Y, ncomp = ncomp, design = design, near.zero.var = T)
        } else {
          res[[l]] = block.spls(X = dats[[l]], Y = Y, ncomp = ncomp, design = design, mode = "regression", near.zero.var = T)
        }
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

        if("prop_expl_var" %in% names(res[[l]])){
          var.vec <- res[[l]]$prop_expl_var
        }else if("explained_variance" %in% names(res[[l]])){
          var.vec <- res[[l]]$explained_variance
        }else{
          var.vec <- list(mic = 0, met = 0);
        }
      }

      loadingNames[[l]]=rownames(loading.pos.xyz[[l]])
      names = lapply(pos.xyz, function(x) rownames(x))
      newmeta = current.proc$meta_para$sample_data

      misc[[l]]$pct2[["microbiome"]] = unname(signif(as.numeric(var.vec[['mic']],4)))*100
      misc[[l]]$pct2[["metabolomics"]] = unname(signif(as.numeric(var.vec[['met']],4)))*100

    }

  }else if(reductionOpt == "procrustes"){
    # NOTE: Run entire procrustes workflow in SINGLE callr subprocess
    # This replaces 4N+2 individual callr calls with just 1
    proc_result <- procrustes_analysis_isolated(
      mic_data_list = d.list$mic$data.proc,
      met_data = d.list$met$data.proc,
      choices = c(1,2,3),
      permutations = 999
    )

    # Unpack results from isolated execution
    res <- proc_result$proc_res
    res2 <- proc_result$prot_res
    misc <- proc_result$misc
    pos.xyz <- proc_result$pos.xyz

    names = lapply(pos.xyz,function(x) make.unique(as.character(rownames(x))))
    newmeta$omics[c(1:(length(names[[1]])/2))] = "microbiome"
    newmeta$omics[c( ((length(names[[1]])/2)+1) : (length(names[[1]])))] = "metabolomics"
    dim.res <- vector("list",length=length(res))
    names(dim.res) <- names(misc) <- names(res)
    for(i in 1:length(dim.res)){
      dim.res[[i]] <- list(res[[i]],res2[[i]])
    }

    procrustes.res <- list(misc=misc,dim.res=dim.res)
    procrustes.res$pos.xyz = lapply(pos.xyz,function(x) unitAutoScale(x))
    procrustes.res$newmeta = newmeta
    combined.res$meta = newmeta
    shadow_save(combined.res,"combined.res.qs")
    shadow_save(procrustes.res,"procrustes.res.qs")
  }

    pos.xyz <- lapply(pos.xyz,function(x) as.data.frame(x)[,c(1:3)]);
    pos.xyz <- lapply(pos.xyz,function(x) unitAutoScale(x));
    pos.xyz <- lapply(pos.xyz, "rownames<-", names[[1]]);

  if(reductionOpt %in% c("diablo")){
    names2 <- lapply(pos.xyz2, function(x) rownames(x))
    pos.xyz <- lapply(pos.xyz,function(x) as.data.frame(x)[,c(1:3)]);
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
    if(micDataType=="ko"){
      loadingSymbols=list(OTU=loadingSymbols)
    }
  }
  if(reductionOpt == "diablo"){
    loading.pos.xyz = lapply(loading.pos.xyz,as.data.frame)
    loading.pos.xyz <- lapply(loading.pos.xyz,unitAutoScale);
    
    diablo.res$dim.res <- res
    diablo.res$pos.xyz <- pos.xyz
    diablo.res$pos.xyz2 <- pos.xyz2
    diablo.res$misc <- misc
    diablo.res$loading.pos.xyz <- loading.pos.xyz
    diablo.res$loadingNames <- loadingNames
    diablo.res$loading.enrich = loadingSymbols
    diablo.res$newmeta=newmeta
    shadow_save(combined.res,"combined.res.qs")
    shadow_save(diablo.res,"diablo.res.qs")
    
  }
  reductionOptGlobal <<- reductionOpt;
  
  mbSetObj$analSet$dimrdc<-reductionOpt;
  
  .set.mbSetObj(mbSetObj);

  return(1)
}



