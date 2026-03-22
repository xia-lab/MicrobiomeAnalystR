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
    if(!exists("run_mofa")){
      .load.scripts.on.demand("mofa_core.Rc");
      .load.scripts.on.demand("util_mofa.Rc");
    }
    library(MOFA2)
    library(reshape2)

    mic_mat <- as.matrix(d.list$mic$data.proc[[1]])
    met_mat <- as.matrix(d.list$met$data.proc)
    data.list.mofa <- list(mic = mic_mat, met = met_mat)

    MOFAobject <- create_mofa_from_matrix(data.list.mofa)
    data_opts <- get_default_data_options(MOFAobject)
    model_opts <- get_default_model_options(MOFAobject)
    model_opts$num_factors <- 5
    train_opts <- get_default_training_options(MOFAobject)
    MOFAobject <- prepare_mofa(object = MOFAobject, data_options = data_opts,
                                model_options = model_opts, training_options = train_opts)
    model <- run_mofa(MOFAobject, save_data = FALSE)

    factors <- get_factors(model, as.data.frame = TRUE)
    mofa_pos <- reshape2::dcast(factors, sample ~ factor, value.var = "value")
    rownames(mofa_pos) <- mofa_pos$sample
    mofa_pos <- mofa_pos[, -1]

    weights <- get_weights(model, as.data.frame = TRUE)
    mofa_loading <- reshape2::dcast(weights, feature ~ factor, value.var = "value")
    mofa_loading$ids <- as.character(mofa_loading$feature)
    mofa_loading <- mofa_loading[, -1]

    var.exp <- model@cache[["variance_explained"]][["r2_per_factor"]][[1]] / 100
    var.exp <- round(var.exp, digits = 3)

    tax_name <- names(d.list$mic$data.proc)[1]
    mofa_pos_xyz <- mofa_pos[, 1:min(3, ncol(mofa_pos)), drop=FALSE]
    pos.xyz <- setNames(list(mofa_pos_xyz), tax_name)
    names <- list(rownames(mofa_pos_xyz))
    names(names) <- tax_name

    loading_ids <- mofa_loading$ids
    loading_xyz <- mofa_loading[, 1:min(3, ncol(mofa_loading) - 1), drop=FALSE]
    loading.pos.xyz <- setNames(list(as.data.frame(loading_xyz)), tax_name)
    loadingNames <- setNames(list(loading_ids), tax_name)

    met_features <- rownames(d.list$met$data.proc)
    omics_vec <- ifelse(loading_ids %in% met_features, "metabolomics", "microbiome")

    misc <- setNames(list(list(
      pct2 = list(microbiome = var.exp[, 1] * 100, metabolomics = var.exp[, 2] * 100),
      var.exp = var.exp
    )), tax_name)

    mofa.res <- list()
    mofa.res$pos.xyz <- pos.xyz
    mofa.res$loading.pos.xyz <- loading.pos.xyz
    mofa.res$loadingNames <- loadingNames
    mofa.res$misc <- misc
    mofa.res$newmeta <- newmeta
    mofa.res$omics_vec <- omics_vec
    shadow_save(combined.res, "combined.res.qs")
    shadow_save(mofa.res, "mofa.res.qs")
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



