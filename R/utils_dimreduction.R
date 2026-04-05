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
  # Use normalized/auto-scaled data for integration
  if(micDataType=='ko'){
    d.list[["mic"]][["data.proc"]] = list(OTU = current.proc$mic$data.proc)
  } else {
    # Wrap as named list matching the analyzed taxonomy level
    analyzed_lvl <- names(which(!sapply(phyloseq_objs$res_deAnal, is.null)))
    if(length(analyzed_lvl) == 0) analyzed_lvl <- names(phyloseq_objs$count_tables)
    d.list[["mic"]][["data.proc"]] = setNames(list(current.proc$mic$data.proc), analyzed_lvl[1])
  }
  if(micDataType=='ko'){
    d.list[["mic"]][["comp.res"]]  = current.proc$mic$res_deAnal[,c(3,4,1)]
    names(d.list[["mic"]][["comp.res"]])[3] = "T.Stats"
    d.list[["mic"]][["comp.res"]] = list(OTU=d.list[["mic"]][["comp.res"]])
    d.list[["mic"]][["enrich.nms"]] = list(OTU=rownames(current.proc$mic$res_deAnal))
    
  }else{
    d.list[["mic"]][["comp.res"]] = lapply(phyloseq_objs$res_deAnal, function(x){names(x)[1] ="T.Stats"; return(x[,c(3,4,1), drop=FALSE])})
    d.list[["mic"]][["enrich.nms"]] = lapply(phyloseq_objs$res_deAnal ,function(x) rownames(x))
  }
  d.list[["mic"]][["meta"]] = data.frame(mbSetObj$dataSet$sample_data)
  
  d.list[["met"]] = list()
  d.list[["met"]][["data.proc"]] = if(!is.null(current.proc$met$data.norm)) current.proc$met$data.norm else current.proc$met$data.proc
  d.list[["met"]][["comp.res"]] =   current.proc$met$res_deAnal[,c(1:3), drop=FALSE] #comp.res
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
    meta.df <- current.proc$meta_para$sample_data
    if(!(analysisVar %in% colnames(meta.df))){
      stop(paste0("[DIABLO] Analysis variable '", analysisVar, "' is not available in sample metadata."))
    }

    diablo.meta.type <- mbSetObj$dataSet$meta.types[analysisVar]
    diablo.meta.type <- if(length(diablo.meta.type) > 0) as.character(diablo.meta.type[[1]]) else NA_character_
    if(is.na(diablo.meta.type) || !(diablo.meta.type %in% c("disc", "cont"))){
      meta.var.raw <- meta.df[, analysisVar, drop = TRUE]
      if(is.numeric(meta.var.raw)){
        diablo.meta.type <- "cont"
      } else {
        meta.num <- suppressWarnings(as.numeric(as.character(meta.var.raw)))
        if(sum(!is.na(meta.num)) == length(meta.var.raw) && length(unique(meta.num)) > 5){
          diablo.meta.type <- "cont"
        } else {
          diablo.meta.type <- "disc"
        }
      }
      message("[DIABLO] Metadata type for '", analysisVar, "' was missing/invalid; inferred as '", diablo.meta.type, "'.")
    }

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

      # Prepare Y and remove invalid samples for this analysis variable
      if(isTRUE(diablo.meta.type == "disc")){
        Y_factor <- as.factor(meta.df[, analysisVar, drop = TRUE])
        keep.inx <- !is.na(Y_factor) & as.character(Y_factor) != ""
        Y_factor <- droplevels(Y_factor[keep.inx])
      } else {
        meta.var <- meta.df[, analysisVar, drop = TRUE]
        Y_numeric <- suppressWarnings(as.numeric(as.character(meta.var)))
        keep.inx <- is.finite(Y_numeric)
        Y_numeric <- Y_numeric[keep.inx]
      }

      dats[[l]] <- lapply(dats[[l]], function(x) x[keep.inx, , drop = FALSE])

      n_samples <- nrow(dats[[l]]$mic)
      if(n_samples < 2){
        stop(paste0("[DIABLO] Not enough valid samples for '", analysisVar, "' after filtering missing/non-numeric values."))
      }

      max_ncomp_by_samples <- n_samples - 1  # Maximum components = n_samples - 1

      # Calculate appropriate ncomp considering sample size and class size
      if(isTRUE(diablo.meta.type == "disc")){
        if(length(unique(Y_factor)) < 2){
          stop(paste0("[DIABLO] Analysis variable '", analysisVar, "' needs at least 2 groups after filtering."))
        }
        min_class_size <- min(table(Y_factor))
        max_ncomp_by_classes <- min_class_size - 1  # Each class needs enough samples
        max_ncomp <- max(1, min(max_ncomp_by_samples, max_ncomp_by_classes, ncol(dats[[l]]$mic), dimn))
        Y <- Y_factor
      } else {
        max_ncomp <- max(1, min(max_ncomp_by_samples, ncol(dats[[l]]$mic), dimn))
        Y <- matrix(Y_numeric, ncol = 1)
        rownames(Y) <- rownames(dats[[l]]$mic)
      }

      ncomp <- max_ncomp

      design = matrix(diabloPar, ncol = length(dats[[l]]), nrow = length(dats[[l]]),
                      dimnames = list(names(dats[[l]]), names(dats[[l]])));
      diag(design) = 0;

      library(mixOmics)

      # Determine if we can use near.zero.var (requires sufficient samples for internal CV)
      # near.zero.var uses internal 10-fold CV by default
      use_nzv <- TRUE
      if(isTRUE(diablo.meta.type == "disc")){
        min_class_size <- min(table(Y))
        # Need at least 10 samples per class for default 10-fold CV
        if(n_samples < 10 || min_class_size < 10) {
          use_nzv <- FALSE
          message("[DIABLO] Disabling near-zero variance filtering due to small sample size (n=", n_samples, ", min class=", min_class_size, ")")
        }
      } else {
        # For regression, need at least 10 total samples
        if(n_samples < 10) {
          use_nzv <- FALSE
          message("[DIABLO] Disabling near-zero variance filtering due to small sample size (n=", n_samples, ")")
        }
      }

      if(isTRUE(diablo.meta.type == "disc")){
        res[[l]] = block.splsda(X = dats[[l]], Y = Y, ncomp = ncomp, design = design, near.zero.var = use_nzv)
      } else {
        res[[l]] = block.spls(X = dats[[l]], Y = Y, ncomp = ncomp, design = design, mode = "regression", near.zero.var = use_nzv)
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

      loadingNames[[l]]=rownames(loading.pos.xyz[[l]])
      names = lapply(pos.xyz, function(x) rownames(x))
      newmeta = current.proc$meta_para$sample_data

      misc[[l]]$pct2[["microbiome"]] = unname(signif(as.numeric(var.vec[['mic']],4)))*100
      misc[[l]]$pct2[["metabolomics"]] = unname(signif(as.numeric(var.vec[['met']],4)))*100

    }

  }else if(reductionOpt == "procrustes"){
    require(vegan)

    mic_data_list <- d.list$mic$data.proc
    met_data <- d.list$met$data.proc

    ndat1 <- lapply(mic_data_list, function(x) decostand(t(x), method = "standardize"))
    pca.dat1 <- lapply(ndat1, function(x) rda(x))
    ndat2 <- decostand(t(met_data), method = "standardize")
    pca.dat2 <- rda(ndat2)

    res <- lapply(pca.dat1, function(x) procrustes(x, pca.dat2, choices = c(1,2,3), symmetric = TRUE, scale = TRUE))
    res2 <- lapply(pca.dat1, function(x) protest(X = x, Y = pca.dat2, scores = "sites", permutations = 999))

    misc <- lapply(res2, function(x) list(`Sum of Squares` = x$ss, Significance = x$signif, Correlation = x$scale))
    pos.xyz <- lapply(res, function(x) rbind(x$X, x$Yrot))

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

  } else if(reductionOpt == "mofa"){
    # Save input for external MOFA script (avoids HDF5Array conflicts in Rserve)
    mic_mat <- as.matrix(d.list$mic$data.proc[[1]])
    met_mat <- as.matrix(d.list$met$data.proc)
    tax_name <- names(d.list$mic$data.proc)[1]
    met_features <- rownames(d.list$met$data.proc)

    mofa.input <- list(
      data.list = list(mic = mic_mat, met = met_mat),
      tax_name = tax_name,
      met_features = met_features
    )
    shadow_save(mofa.input, "mofa_input.qs")
    shadow_save(combined.res, "combined.res.qs")

    # Return 2 to signal Java to run _perform_mofa.R externally
    return(2)
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
    } else if(!is.list(loadingSymbols)){
      loadingSymbols=setNames(list(loadingSymbols), names(loadingNames))
    }
  }
  if(reductionOpt == "mofa"){
    loading.pos.xyz = lapply(loading.pos.xyz, as.data.frame)
    loading.pos.xyz <- lapply(loading.pos.xyz, unitAutoScale)
    mofa.res$pos.xyz <- pos.xyz
    mofa.res$loading.pos.xyz <- loading.pos.xyz
    mofa.res$loadingNames <- loadingNames
    mofa.res$loading.enrich <- loadingSymbols
    shadow_save(mofa.res, "mofa.res.qs")
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


