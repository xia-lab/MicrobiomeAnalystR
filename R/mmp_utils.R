##############################################################
## R script for MicrobiomeAnalyst
## Description: Functions for microbiome metabolomics analysis
## Author: Jeff Xia, jeff.xia@mcgill.ca
################################################################



#####################################################
#############processing functions####################
#####################################################

CreateMMPFakeFile <- function(mbSetObj,isNormalized="true",isNormalizedMet="true",module.type){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(isNormalized=="false" & isNormalizedMet=="false"){
    
    AddErrMsg("Please make sure your data has been normalized properly!");
    return(0)
    
  }
  current.msg <<- ""
  
  if(isNormalized=="true"){
    
    mbSetObj$dataSet$filt.data <- mbSetObj$dataSet$data.orig
    mbSetObj$dataSet$filt.msg <- "No data filtering has been performed for microbiome data since it has been transformed."
    
    mbSetObj$dataSet$norm.phyobj <- mbSetObj$dataSet$proc.phyobj
    mbSetObj$dataSet$norm.msg <- "No normalization has been performed for microbiome data since it has been transformed."
    
    #make hierarchies
    
    ranks <- c(GetMetaTaxaInfo(mbSetObj), "OTU")
    ranks <- unique(ranks)
    data.list <- list()
    data.list$merged_obj <- vector(length = length(ranks), "list")
    data.list$count_tables <- vector(length = length(ranks), "list")
    names(data.list$count_tables) <- names(data.list$merged_obj) <- ranks
    
    
    for(i in 1:length(ranks)){
      phyloseq.obj <- UtilMakePhyloseqObjs(mbSetObj, ranks[i])
      data.list$merged_obj[[i]] <- phyloseq.obj
      count.table <- UtilMakeCountTables(phyloseq.obj, ranks[i])
      data.list$count_tables[[i]] <- count.table
    }
    
    qs::qsave(data.list,"prescale.phyobj.qs")
    current.proc$mic$data.proc<<- data.list$count_tables[["OTU"]]
    
    saveDataQs(data.list, "phyloseq_prenorm_objs.qs",module.type, mbSetObj$dataSet$name);  
    saveDataQs(data.list, "phyloseq_objs.qs",module.type, mbSetObj$dataSet$name);  
    
  }
  
  if(isNormalizedMet=="true"){
    
    mbSetObj$dataSet$metabolomics$filt.data <- mbSetObj$dataSet$metabolomics$norm.data <- mbSetObj$dataSet$metabolomics$data.orig; ## feature in row and sample in column
    
    qs::qsave( mbSetObj$dataSet$metabolomics$norm.data, file="metabo.complete.norm.qs");
    
    current.proc$met$data.proc<<-mbSetObj$dataSet$metabolomics$data.orig
    
    mbSetObj$dataSet$metabolomics$norm.msg <- "No normalization has been performed for metabolomics data since it has been transformed."
    
  }
  
  mbSetObj$dataSet$sample_data$sample_id <- rownames(mbSetObj$dataSet$sample_data);
  return(.set.mbSetObj(mbSetObj));
}

#####################################################
#############differential analysis###################
#####################################################
#'Perform differential analysis
#'@description This functions performs alpha diversity.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param taxalvl Character, input taxonomy level
#'@param metadata Character, input the name of the experimental factor
#'to group the samples.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

PerformDEAnalyse<- function(mbSetObj, taxalvl="Genus",netType="gem",overlay,initDE=1,
                            analysisVar="CLASS",adjustedVar,alg="limma",plvl=0.05, fc.lvl=1, selected="NA",nonpar=FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  require(dplyr)
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
  }
  
  plvl<<-plvl
  analysisVar<<-analysisVar
  alg<<-alg
  metdat <- mbSetObj$dataSet$metabolomics$norm.data
  sample_data <-  mbSetObj$dataSet$sample_data
  sample_type <- mbSetObj$dataSet$meta_info
  
  metdat.de <- performLimma(metdat,sample_data,sample_type,analysisVar)
  
  if(initDE=="1"){
    phyloseq_objs[["res_deAnal"]] <- vector("list",length=length(phyloseq_objs$count_tables))
    names( phyloseq_objs[["res_deAnal"]]) <- names( phyloseq_objs$count_tables)
    
    micdat <- phyloseq_objs$count_tables
    micdat.de <- lapply(micdat,function(x) performLimma(x,sample_data,sample_type,analysisVar))
    predres.met <- qs::qread(paste0("m2m_pred_",predDB,".qs"))
    predres.met <- lapply(predres.met,function(x) return(x$fun_prediction_met))
    predDE<- vector("list",length=length(predres.met))
    pred.dat<- vector("list",length=length(predres.met))
    if(netType=="gem" & overlay =="true"){
      if(!(file.exists(paste0("m2m_pred_",predDB,".qs")))){
        
        AddErrMsg("Cannot import the prediction result!")
      }else{
        
        for(i in 1:length( predres.met)){
          taxalvl2<- names( predres.met)[i]
          m2m_pred <-  predres.met[[taxalvl2]]
          m2m_pred <- lapply(m2m_pred,function(x) reshape2::melt(x) )
          m2m_for_de <- mapply(`[<-`, m2m_pred, 'sample', value = names(m2m_pred), SIMPLIFY = FALSE)
          m2m_for_de <- do.call(rbind,m2m_for_de)
          rownames(m2m_for_de) <- NULL
          m2m_for_de$pair <- apply(m2m_for_de[,1:2],1,function(x) paste(x,collapse = ";;"))
          tokeep <- aggregate(m2m_for_de$value,list(m2m_for_de$pair),sum) %>% filter(x!=0)
          m2m_for_de <- m2m_for_de %>% filter(pair %in% tokeep$Group.1)
          m2m_pair_dat <- reshape2::dcast(m2m_for_de,pair~sample,value.var = "value")
          pred.dat[[taxalvl2]] <- m2m_pair_dat
          #qs::qsave(list(m2m_for_de=m2m_for_de,m2m_pair_dat=m2m_pair_dat),"m2m_pair_pred.qs")
          rownames(m2m_pair_dat) <- m2m_pair_dat$pair
          m2m_pair_dat$pair <- NULL
          m2m_pair_de <- performLimma(m2m_pair_dat,sample_data,sample_type,analysisVar)
          m2m_pair_de$mic <- m2m_for_de$Var1[match(rownames(m2m_pair_de),m2m_for_de$pair)]
          m2m_pair_de$met <- m2m_for_de$Var2[match(rownames(m2m_pair_de),m2m_for_de$pair)]
          predDE[[taxalvl2]]<-m2m_pair_de
          
        }
        
        qs::qsave(list(pred.dat=pred.dat,predDE=predDE),"m2m_pair_de.qs")
      }
      
    }else if(netType=="kegg"){
      
    }
    
    
  }else{
    
    micdat <- phyloseq_objs$count_tables[[taxalvl]]
    micdat.de <- performLimma(micdat,sample_data,sample_type,analysisVar)
    
    if(netType=="gem" & overlay =="true" &  taxalvl != "OTU"){
      if(!(file.exists(paste0(tolower(taxalvl),"_metabolite_pred_pair.qs")))){
        
        AddErrMsg("Cannot import the prediction result!")
      }else{
        m2m_pred <- qs::qread(paste0(tolower(taxalvl),"_metabolite_pred_pair.qs"))
        m2m_pred <- lapply(m2m_pred,function(x) reshape2::melt(x) )
        m2m_for_de <- mapply(`[<-`, m2m_pred, 'sample', value = names(m2m_pred), SIMPLIFY = FALSE)
        m2m_for_de <- do.call(rbind,m2m_for_de)
        rownames(m2m_for_de) <- NULL
        m2m_for_de$pair <- apply(m2m_for_de[,1:2],1,function(x) paste(x,collapse = ";;"))
        tokeep <- aggregate(m2m_for_de$value,list(m2m_for_de$pair),sum) %>% filter(x!=0)
        m2m_for_de <- m2m_for_de %>% filter(pair %in% tokeep$Group.1)
        m2m_pair_dat <- reshape2::dcast(m2m_for_de,pair~sample,value.var = "value")
        qs::qsave(list(m2m_for_de=m2m_for_de,m2m_pair_dat=m2m_pair_dat),"m2m_pair_pred.qs")
        rownames(m2m_pair_dat) <- m2m_pair_dat$pair
        m2m_pair_dat$pair <- NULL
        m2m_pair_de <- performLimma(m2m_pair_dat,sample_data,sample_type,analysisVar)
        m2m_pair_de$mic <- m2m_for_de$Var1[match(rownames(m2m_pair_de),m2m_for_de$pair)]
        m2m_pair_de$met <- m2m_for_de$Var2[match(rownames(m2m_pair_de),m2m_for_de$pair)]
        
        mbSetObj$analSet$m2m_pair_de <- m2m_pair_de
        qs::qsave(m2m_pair_de,"m2m_pair_de.qs")
      }
      
    }
  }
  #fast.write(micdat.de, file=paste0(taxalvl,adjustedVar,"_",alg,"_Res.csv"));
  fast.write(metdat.de, file=paste0("metabolite",'_',analysisVar,"_",alg,"_Res.csv"));
  
  phyloseq_objs$res_deAnal <- micdat.de
  mbSetObj$dataSet$metabolomics$res_deAnal <- metdat.de
  qs::qsave(phyloseq_objs,"phyloseq_objs.qs")
  return(.set.mbSetObj(mbSetObj))
}


PerformPairDEAnalyse <- function(mbSetObj, taxalvl, analysisVar,alg="limma",plvl=1, selected="NA",nonpar=FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  require(dplyr)
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
  }
  
  if(taxalvl=="null" | is.null(taxalvl)){
    if(exist(current.proc$taxalvl)){
      
      taxalvl = current.proc$taxalvl
    }else{
      taxalvl=names(phyloseq_objs$count_tables)[length(phyloseq_objs$count_tables)-1]
    }
    
  }
  
  analysisVar <- current.proc$meta_para$analysis.var
  
  #if(analysisVar=="null" | is.null(analysisVar)){
  # analysisVar = names(current.proc$sample)[1]
  # }
  sample_data <-  mbSetObj$dataSet$sample_data
  sample_type <- mbSetObj$dataSet$meta_info
  
  #tempnm <- paste0(analysisVar,"_",alg)
  predDB<- current.proc$predDB
  if(micDataType=="ko"){
    AddErrMsg("Prediction is supportive for KO abundance table! Please check your data type")
  }else{
    
    if(!(file.exists(paste0("m2m_pred_",predDB,".qs")))){
      
      AddErrMsg("Cannot import the prediction result!")
    }else{
      if(taxalvl=="all"){
        predres.met <- qs::qread(paste0("m2m_pred_",predDB,".qs"))
        predres.met <- lapply(predres.met,function(x) return(x$fun_prediction_met))
        
        predDE<- vector("list",length=length(predres.met))
        pred.dat<- vector("list",length=length(predres.met))
        names(predDE) <- names(pred.dat) <- names(predres.met)
        for(tax in names( predres.met)){
          m2m_pred <-  predres.met[[tax]]
          m2m_pred <- lapply(m2m_pred,function(x) reshape2::melt(x) )
          m2m_for_de <- mapply(`[<-`, m2m_pred, 'sample', value = names(m2m_pred), SIMPLIFY = FALSE)
          m2m_for_de <- do.call(rbind,m2m_for_de)
          rownames(m2m_for_de) <- NULL
          m2m_for_de$pair <- apply(m2m_for_de[,1:2],1,function(x) paste(x,collapse = ";;"))
          tokeep <- aggregate(m2m_for_de$value,list(m2m_for_de$pair),sum) %>% filter(x!=0)
          m2m_for_de <- m2m_for_de %>% filter(pair %in% tokeep$Group.1)
          m2m_pair_dat <- reshape2::dcast(m2m_for_de,pair~sample,value.var = "value")
          # m2m_pair_dat[,-1] <- t(apply(m2m_pair_dat[,-1],1,function(x) ReScale(x,0,1)))
          pred.dat[[tax]] <- m2m_pair_dat
          rownames(m2m_pair_dat) <- m2m_pair_dat$pair
          m2m_pair_dat$pair <- NULL
          m2m_pair_de <- performLimma(m2m_pair_dat,sample_data,sample_type,analysisVar)
          m2m_pair_de$mic <- m2m_for_de$Var1[match(rownames(m2m_pair_de),m2m_for_de$pair)]
          m2m_pair_de$met <- m2m_for_de$Var2[match(rownames(m2m_pair_de),m2m_for_de$pair)]
          predDE[[tax]]<-m2m_pair_de
          
        }
        m2m_pair_de <- list()
        m2m_pair_de$pred.dat <- pred.dat
        m2m_pair_de$predDE <- predDE
        qs::qsave(m2m_pair_de,"m2m_pair_de.qs")
        
        
      }else{
        if(!exists("predres",current.proc)){
          predres.met <- qs::qread(paste0("m2m_pred_",predDB,".qs"))
          m2m_pred <- predres.met[[taxalvl]]$fun_prediction_met
        }else{
          m2m_pred <- current.proc$predres$fun_prediction_met
        }
        m2m_pred <- lapply(m2m_pred,function(x) reshape2::melt(x) )
        m2m_for_de <- mapply(`[<-`, m2m_pred, 'sample', value = names(m2m_pred), SIMPLIFY = FALSE)
        m2m_for_de <- do.call(rbind,m2m_for_de)
        rownames(m2m_for_de) <- NULL
        m2m_for_de$pair <- apply(m2m_for_de[,1:2],1,function(x) paste(x,collapse = ";;"))
        tokeep <- aggregate(m2m_for_de$value,list(m2m_for_de$pair),sum) %>% filter(x!=0)
        m2m_for_de <- m2m_for_de %>% filter(pair %in% tokeep$Group.1)
        m2m_pair_dat <- reshape2::dcast(m2m_for_de,pair~sample,value.var = "value")
        # m2m_pair_dat[,-1] <- t(apply(m2m_pair_dat[,-1],1,function(x) ReScale(x,0,1)))
        pred.dat <- m2m_pair_dat
        rownames(m2m_pair_dat) <- m2m_pair_dat$pair
        m2m_pair_dat$pair <- NULL
        m2m_pair_de <- performLimma(m2m_pair_dat,sample_data,sample_type,analysisVar)
        m2m_pair_de$mic <- m2m_for_de$Var1[match(rownames(m2m_pair_de),m2m_for_de$pair)]
        m2m_pair_de$met <- m2m_for_de$Var2[match(rownames(m2m_pair_de),m2m_for_de$pair)]
        current.proc$pred.dat <<- pred.dat
        current.proc$predDE <<- m2m_pair_de 
        fast.write(m2m_pair_de, file=paste0("prediction_differential.csv"));
      }
      
    }
    
  }
  
  return(.set.mbSetObj(mbSetObj))
  
}

resMsg <<- "";

CompareMic <- function(mbSetObj, taxalvl,initDE=1,
                       analysis.var,comp = NULL,ref = NULL,block = "NA",
                       alg="maaslin",plvl=0.05, 
                       is.norm=FALSE,selected="NA",nonpar=FALSE){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  current.proc$mic$alg<<-alg
  current.proc$mic$plvl<<-plvl
  sample_data <-  data.frame(mbSetObj$dataSet$sample_data)
  sample_type <- mbSetObj$dataSet$meta_info
  meta_type <- mbSetObj$dataSet$meta.types

  if (!exists('adj.vec')) {
    adj.bool = F;
  } else {
    if (length(adj.vec) > 0) {
        if(length(adj.vec)==1){
            if(adj.vec == ""){
                adj.bool = F;
            } else {
                adj.bool = T;
            }
        } else {
           adj.bool = T;
        }      
    } else {
      adj.bool = F;
    }
  }
  adj.vars <- adj.vec;
  
  
  meta.nms <- colnames(sample_data)
    library(dplyr)
  input.meta <-sample_data@.Data %>% as.data.frame()
  colnames(input.meta) <- meta.nms
  rownames(input.meta) <- input.meta$sample_id
  
  if(adj.bool){
    fixed.effects <- c(analysis.var, adj.vars)
    fixed.types <- meta_type[names(meta_type) %in% fixed.effects]
    fixed.types <- fixed.types[match(fixed.effects, names(fixed.types))]
  } else { # to do still
    fixed.effects <- analysis.var
    fixed.types <- meta_type[names(meta_type) == analysis.var]
  }
  
  analysis.type <- fixed.types[fixed.effects == analysis.var]
  disc.effects <- fixed.effects[fixed.types == "disc"]
  
  # build refs vector (may need to add for blocking too)
  if(length(disc.effects) > 0){
    if(analysis.type == "disc"){
      refs <- paste0(analysis.var, ",", ref)
      if(length(disc.effects) > 1){
        for(i in c(2:length(disc.effects))){
          ref.temp <- paste0(disc.effects[i], ",", levels(unlist(c(input.meta[,disc.effects[i]])))[1])
          refs <- c(refs, ref.temp)
        }
      }
    } else {
      refs <- c()
      if(length(disc.effects) > 1){
        for(i in c(1:length(disc.effects))){
          ref.temp <- paste0(disc.effects[i], ",", levels(unlist(c(input.meta[,disc.effects[i]])))[1])
          refs <- c(refs, ref.temp)
        }
      }
    }
  }
  
  # MaAslin does not require samples or orders to exactly match - it takes care of this
  # set normalized/transformation parameters
  if(is.norm == "false"){
    phyloseq_objs <- qs::qread("phyloseq_prenorm_objs.qs")
    norm.method = "TSS"
    trans.method = "LOG"
  } else {
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    norm.method = "NONE"
    trans.method = "NONE"
  }
  
  current.proc$meta_para<<-list(analysis.var=analysis.var,sample_data=sample_data,
                                sample_type=sample_type, input.meta=input.meta,
                                fixed.effects=fixed.effects,analysis.type=analysis.type,
                                disc.effects=disc.effects,comp=comp,ref=ref,refs=refs,block=block,
                                norm.method=norm.method,trans.method =trans.method)
  
  if(micDataType=="ko"){
    micdat <- phyloseq_objs$count_tables[["OTU"]]
    micdat.de <- doMaAslin(micdat,plvl)
  }else{
    if(initDE=="1"|taxalvl=="all"){
      micdat <- phyloseq_objs$count_tables
      micdat.de <- lapply(micdat,function(x) doMaAslin(x,plvl))
    }else{
      micdat <- phyloseq_objs$count_tables[[taxalvl]]
      micdat.de <- doMaAslin(micdat,plvl)
    }
    
  }
  return(micdat.de) 
}



CompareMet <- function(mbSetObj, analysisVar,
                       alg="limma",plvl=0.05,ref,compr, selected="NA",nonpar=FALSE){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  # current.proc$sample<<-data.frame(mbSetObj$dataSet$sample_data)
  require(dplyr)

  if(analysisVar=="null" | is.null(analysisVar)){
    analysisVar = names(current.proc$meta_para$sample_data)[1]
  }
  current.proc$met$plvl<<-plvl
  current.proc$met$alg<<-alg
 
  sample_data <-  mbSetObj$dataSet$sample_data
  sample_data <- sample_data[sample_data[[analysisVar]] %in% c(ref,compr),]
  sample_type <- mbSetObj$dataSet$meta_info
  metdat <-current.proc$met$data.proc %>% 
         .[,colnames(.)%in%sample_data$sample_id]
  metdat.de <- performLimma(metdat,sample_data,sample_type,analysisVar)
  fast.write(metdat.de, file="limma_output.csv");
  current.proc$met$res_deAnal <<- metdat.de
  mbSetObj$dataSet$metabolomics$resTable <- metdat.de
  
  sigfeat <- rownames(metdat.de)[metdat.de$FDR < plvl];
  sig.count <- length(sigfeat);
  if(sig.count == 0){
    current.msg <<- "No significant metabolomic features were identified using the given p value cutoff.";
  }else{
    if(metDataType=="peak"){
      current.msg <<-  paste(sig.count, "significant peaks were identified!");
    }else{
      current.msg <<- paste(sig.count, "significant metabolites were identified!");
    }
    
  }

  mbSetObj$dataSet$metabolomics$sigfeat <- sigfeat
  mbSetObj$dataSet$metabolomics$sig.count <- sig.count
  current.proc$met$sigfeat <<- sigfeat 
  resMsg<<- paste(resMsg,current.msg)
  return(.set.mbSetObj(mbSetObj))
  
}


performLimma <-function(data,sample_data,sample_type,analysisVar){
  require(limma);
  covariates <- data.frame(sample_data)
  if(is.null(covariates$sample_id)){
    covariates$sample_id <- rownames(covariates)
  }
  
  sample_type <- lapply(sample_type, function(x) return(x[x]))
  
  for(i in 1:(ncol(covariates)-1)){ # ensure all columns are the right type
    if(names(covariates)[i] %in% names(sample_type[["disc.inx"]])){
      covariates[,i] <- covariates[,i] %>% make.names() %>% factor()
    } else {
      covariates[,i] <- covariates[,i] %>% as.character() %>% as.numeric()
    }
  }
  
  covariates <- data.frame(covariates[,-ncol(covariates),drop=F])
  
  if(!exists('adj.vec')){
    adj.bool = F;
  }else{
    if(length(adj.vec) > 0){
      adj.bool = T;
      adj.vars <- adj.vec;
      if(length(adj.vec) == 1){
        if(adj.vec == ""){
            adj.bool = F;
        }
      }
    }else{
      adj.bool = F;
    }
  }
  
  feature_table = as.matrix(data[,which(colnames(data) %in% rownames(covariates))])
  covariates <- covariates[match(colnames(feature_table), rownames(covariates)),,drop=F]
  analysis.var <- analysisVar
  analysis.type <- ifelse(analysis.var %in% names(sample_type[["disc.inx"]]),"disc","count")
  if(adj.bool){
    vars <- c(analysis.var, adj.vars)
  }else{
    vars <- analysis.var
  }
  if(analysis.type == "disc"){
    covariates[, analysis.var] <- covariates[, analysis.var] %>% make.names() %>% factor();
    grp.nms <- unique(c(current.proc$meta_para$comp,current.proc$meta_para$ref,levels(covariates[, analysis.var])))
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data =covariates );
    if(adj.bool){
      
      nms=sapply(seq(adj.vars), function(x) nms= levels(covariates[,adj.vars[x]])[-1])
      nms=sapply(seq(adj.vars), function(x) {
        if(!(is.null(levels(covariates[,adj.vars[x]])))){
          return(levels(covariates[,adj.vars[x]])[-1])
        }else{
          return(adj.vars[x])
        }
      })
      colnames(design) = c(grp.nms[order(grp.nms)],unlist(nms))
    }else{
      colnames(design) =  grp.nms[order(grp.nms)]
      
    }
    inx = 0;
    myargs <- list();
    for(m in 1:(length(grp.nms)-1)){
      for(n in (m+1):length(grp.nms)){
        inx <- inx + 1;
        myargs[[inx]] <- paste(grp.nms[m], "-", grp.nms[n], sep="")
      }
    }
    
    myargs[["levels"]] <- design;
    
    contrast.matrix <- do.call(makeContrasts, myargs);
    fit <- lmFit(feature_table, design)
    fit <- contrasts.fit(fit, contrast.matrix);
    fit <- eBayes(fit);
    if(length(levels( covariates[, analysis.var]))==2){
      topFeatures <- topTable(fit, number = Inf);
      res = data.frame( P_value=signif(topFeatures[,"P.Value"] , digits = 3), 
                        FDR=signif(topFeatures[,"adj.P.Val"], digits = 3),
                        T.Stats=signif(topFeatures[,"t"], digits = 3),
                        Log2FC=signif(topFeatures[,"logFC"], digits = 3))
      rownames(res) <- rownames(topFeatures);

    }else{
      
      res <- data.frame(P_value=signif(fit$p.value[,1],digits = 3),
                        FDR=signif(p.adjust(fit$p.value[,1],"fdr"),digits = 3),
                        T.Stats=signif(fit$t[,1],digits = 3),
                        F.Stats=signif(fit$F,digits = 3),
                        F.Pval=signif(fit$F.p.value,digits = 3))
      rownames(res) <- rownames(fit$p.value);
      
    }
    
  } else { 
    
    covariates[, analysis.var] <- covariates[, analysis.var] %>% as.numeric();
    design <- model.matrix(formula(paste0("~ 0", paste0(" + ", vars, collapse = ""))), data = covariates);
    fit <- eBayes(lmFit(feature_table, design));
    rest <- topTable(fit, number = Inf, coef = analysis.var);
    colnames(rest)[1] <- analysis.var;
    ### get results with no adjustment
    # design <- model.matrix(formula(paste0("~ 0", paste0(" + ", analysis.var, collapse = ""))), data = covariates);
    # fit <- eBayes(lmFit(feature_table, design));
    # topFeatures <- topTable(efit, number = Inf, adjust.method = "fdr");
  }
  if(length(which(duplicated(rownames(fit$p.value))))>0){
    current.msg<<-"Duplicate features names are not allowed! Please double check your input!"
    return()
  }else{
    rownames(fit$p.value[order(fit$p.value),,drop=F])
  } 

  res <- na.omit(res)
  res <- res[order(res[,2], decreasing=FALSE),]
  res[res == "NaN"] = 1
  return(res)
  
}


doMaAslin <- function(input.data,thresh = 0.05,adj.bool=F){
  
  require(dplyr);
  require(R.utils);
  if(.on.public.web){
    # make this lazy load
    if(!exists(".prepare.maaslin2")){ # public web on same user dir
      .load.scripts.on.demand("utils_maaslin.Rc");    
    }
  }
  thresh <- as.numeric(thresh);
  input.data <- as.data.frame(input.data)
  block <- current.proc$meta_para$block
  disc.effects <- current.proc$meta_para$disc.effects
  analysis.var <- current.proc$meta_para$analysis.var
  if(block == "NA"){
    if(length(disc.effects) > 0){ # case: discrete variables, no blocking factor
      maaslin.para<<- list(input_data = input.data, 
                           input_metadata = current.proc$meta_para$input.meta, 
                           fixed_effects =  current.proc$meta_para$fixed.effects,
                           reference = current.proc$meta_para$refs,
                           max_significance = 0.05,
                           min_abundance = 0.0,
                           min_prevalence = 0.0,
                           min_variance = 0.0,
                           normalization = current.proc$meta_para$norm.method,
                           transform = current.proc$meta_para$trans.method)
      return(1)
    } else { # case: no discrete variables, no blocking factor
      maaslin.para<<- list(input_data = input.data, 
                           fixed_effects = current.proc$meta_para$fixed.effects,
                           max_significance = 0.05,
                           min_abundance = 0.0,
                           min_prevalence = 0.0,
                           min_variance = 0.0,
                           normalization = current.proc$meta_para$norm.method,
                           transform = trans.method)
      return(1)
    }
  } else { # case: discrete variables, blocking factor (blocking factor must be discrete)
    maaslin.para <<-list(check= list(input_data = input.data[1,], 
                                     input_metadata = current.proc$meta_para$input.meta, 
                                     fixed_effects = current.proc$meta_para$fixed.effects,
                                     random_effects = block,
                                     reference =current.proc$meta_para$refs,
                                     max_significance = 0.05,
                                     min_abundance = 0.0,
                                     min_prevalence = 0.0,
                                     min_variance = 0.0,
                                     normalization = current.proc$meta_para$norm.method,
                                     transform = current.proc$meta_para$trans.method),
                         test=list(input_data = input.data, 
                                   input_metadata = current.proc$meta_para$input.meta, 
                                   fixed_effects = current.proc$meta_para$fixed.effects,
                                   random_effects = block,
                                   reference = current.proc$meta_para$refs,
                                   max_significance = 0.05,
                                   min_abundance = 0.0,
                                   min_prevalence = 0.0,
                                   min_variance = 0.0,
                                   normalization = current.proc$meta_para$norm.method,
                                   transform = current.proc$meta_para$trans.method)
    )
    return(2)
    
  }
  
}



PrepareResTable <- function(mbSetObj,micDataType,taxalvl,is.norm=F){
  load_phyloseq();

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$analSet$maaslin$taxalvl <- taxalvl;
  if(micDataType=="otu"){
    if(is.null(taxalvl)|taxalvl=="null"){
      taxalvl = colnames(mbSetObj$dataSet$taxa_table)[length(colnames(mbSetObj$dataSet$taxa_table))]
    }
    
    resTab = qs::qread("phyloseq_objs.qs")$res_deAnal[[taxalvl]]
    sigfeat <- qs::qread("phyloseq_objs.qs")$sigfeat[[taxalvl]]
    fileName <- paste0(taxalvl,"_maaslin_output.csv");
  }else{
    taxalvl =="OTU"
    resTab = current.proc$mic$res_deAnal
    sigfeat <- current.proc$mic$sigfeat
    fileName <- paste0("maaslin_output.csv");
  }
  
  sig.count <- length(sigfeat);
  if(sig.count == 0){
    resMsg <<- "No significant microbiome features were identified using the given p value cutoff.";
  }else{
    resMsg <<- paste("A total of", sig.count,"significant", ifelse(taxalvl=="OTU",taxalvl,tolower(taxalvl)) , "were identified, ");
  }
  
  if(is.norm){
    phylonm <- "phyloseq_objs.qs"
  }else{
    phylonm <- "phyloseq_prenorm_objs.qs"
  }
  
  input.data = qs::qread(phylonm)$count_tables[[taxalvl]]
  analysis.var = current.proc$meta_para$analysis.var
  # put results in mbSetObj, learn pattern of analysis set
  
  fast.write(resTab, file = fileName);
  compMicFile<<-fileName
  
  
  # process data for individual feature boxplot
  taxrank_boxplot <- taxalvl;
  claslbl_boxplot <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[analysis.var]]);
  nm_boxplot <- rownames(input.data);
  dat3t_boxplot <- as.data.frame(t(input.data),check.names=FALSE);
  colnames(dat3t_boxplot) <- nm_boxplot; 
  box_data <- dat3t_boxplot;
  box_data$class <- claslbl_boxplot;
  box_data$norm <- is.norm;
  
  
  message("Result table done")
  
  mbSetObj$analSet$multiboxdata <- box_data;
  mbSetObj$analSet$sig.count <- sig.count;
  mbSetObj$analSet$resTable <- resTab;
  resMsg<<- paste0(resMsg,current.msg)
  return(.set.mbSetObj(mbSetObj))
}




ProcessMaaslinRes <- function(mbSetObj,taxalvl,analysis.var,thresh){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  input.data<-maaslin.para$input_data
  res <- mbSetObj$analSet$maaslin$results
  inds <- !(res$feature %in% rownames(input.data)); 
  # filter results to get only ones related to analysis var
  res <- res[res$metadata == analysis.var, ];
  
  # make res pretty
  res$coef <- signif(res$coef, digits = 3);
  res$stderr <- signif(res$stderr, digits = 3);
  res$pval <- signif(res$pval, digits = 3);
  res$qval <- signif(res$qval, digits = 3);
  if(current.proc$meta_para$analysis.type == "disc"){
    res <- res[res$value == current.proc$meta_para$comp, ];
    rownames(res) <- res$feature;
    res <- res[ ,c("coef", "stderr", "pval", "qval")];
    colnames(res) <- c("Log2FC", "St.Error", "P_value", "FDR");
  } else {
    rownames(res) <- res$feature;
    res <- res[ ,c("coef", "stderr", "pval", "qval")];
    colnames(res) <- c("Coefficient", "St.Error", "P_value", "FDR");
  }
  
  res = res[order(res$P_value),] 
  # write out/save results
  fileName <-paste0(taxalvl,"_maaslin_output.csv");
  fast.write(res, file = fileName);
  
  plvl<- current.proc$mic$plvl
  if(micDataType=="ko"){
    current.proc$mic$res_deAnal <<- res
    current.proc$mic$sigfeat <<-  rownames(current.proc$mic$res_deAnal)[current.proc$mic$res_deAnal$FDR< plvl]
  }else{ 
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    phyloseq_objs$res_deAnal[[taxalvl]] <- res
    phyloseq_objs$sigfeat[[taxalvl]] <- rownames(phyloseq_objs$res_deAnal[[taxalvl]])[phyloseq_objs$res_deAnal[[taxalvl]]$FDR< plvl]
    qs::qsave(phyloseq_objs,"phyloseq_objs.qs")
    
  }
   mbSetObj$analSet$maaslin$taxalvl <- "OTU"
  return(.set.mbSetObj(mbSetObj))
  
}

#####################################################
##################Prediction#########################
#####################################################
#####################################################


#####lib path
lib.path.mmp <<- "../../lib/mmp/"
  

MetaboIDmap <- function(netModel,predDB,IDtype,met.vec=NA){
  
  # met.vec <- rownames(qs::qread("metabo.complete.norm.qs"))
  if(inputType=="table"){
    met.vec <- rownames(current.proc$met$data.proc)
  }else{
    met.vec <- met.vec
  }
  
  if(netModel=="gem"){
    if(predDB=="agora"){
      metdb <- qs::qread(paste0(lib.path.mmp,"agora.met.qs"))
    }else if(predDB=="embl"){
      metdb <- qs::qread(paste0(lib.path.mmp,"embl.met.qs"))
    }
    
    if(IDtype=="name"){
      metInfo <- qs::qread(paste0(lib.path.mmp,"synonymGem.qs"));
      met.map <- data.frame(Query=met.vec,Match=met.vec,stringsAsFactors = F)
      met.map$Match <-  metInfo$metID[match(tolower(met.map$Query),tolower(metInfo$Name))]
      met.map <- met.map[which(met.map$Match %in% metdb),]
      map.l <- length(unique(met.map$Match))
    }else if(IDtype=="kegg"){
      metInfo <- qs::qread(paste0(lib.path.mmp,"gem2kegg.qs"));
      met.map <- data.frame(Query=met.vec,Match=met.vec,stringsAsFactors = F)
      met.map$Match <-  metInfo$metID[match(met.map$Query,metInfo$KEGG)]
      met.map <- met.map[which(met.map$Match %in% metdb),]
      map.l <- length(unique(met.map$KEGG))
    }else if(IDtype=="hmdb"){
      metInfo <- qs::qread(paste0(lib.path.mmp,"gem2hmdb.qs"));
      met.map <- data.frame(Query=met.vec,Match=met.vec,stringsAsFactors = F)
      met.map$Match <-  metInfo$metID[match(met.map$Query,metInfo$HMDB)]
      met.map <- met.map[which(met.map$Match %in% metdb),]
      map.l <- length(unique(met.map$HMDB))
    }
  }else if(netModel=="keggNet"){
    
    if(IDtype=="name"){
      metInfo <- qs::qread(paste0(lib.path.mmp,"general_kegg2name.qs"));
      met.map <- data.frame(Query=met.vec,Match=met.vec,stringsAsFactors = F)
      met.map$Match <-  metInfo$ID[match(tolower(met.map$Query),tolower(metInfo$Name))]
      met.map$Name <- met.map$Query
      met.map$Node <-  metInfo$id[match(met.map$Query,metInfo$Name)]
      met.map <- met.map[!(is.na(met.map$Match)),]
      map.l <- length(unique(met.map$Match))
    }else if(IDtype=="kegg"){
      
      metInfo <- qs::qread(paste0(lib.path.mmp,"general_kegg2name.qs"));
      met.map <- data.frame(Query=met.vec,Match=met.vec,Name=met.vec,stringsAsFactors = F)
      met.map$Name <-  metInfo$Name[match(met.map$Query,metInfo$ID)]
      met.map$Node <-  metInfo$node[match(met.map$Query,metInfo$ID)]
      met.map <- met.map[!(is.na(met.map$Name)),]
      map.l <- length(unique(met.map$Match))
    }
    
  }
  
  
  if(inputType=="table"){
    qs::qsave(met.map,paste0(netModel,".met.map.qs"))
    
    current.proc[[netModel]]<<-met.map
    fast.write(met.map, file=paste0(netModel,"_metabo_match_result.csv"))
    return(1)
  }else{
    return(met.map)
  }
  
}



MicIDmap <- function(netModel,predDB,taxalvl="all"){
  load_stringr();
  
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
  }
  
  mic.vec <- lapply(phyloseq_objs$count_tables, function(x) return(list(rownames(x))))
  
  mic.vec[["OTU"]] <- NULL
  lvlnm <- c("phylum","class","order","family","genus","species")
  lvlidx <- match(tolower(names(mic.vec)),lvlnm)
  lvlppl <- c("p__","c__","o__","f__","g__","s__")
  lvlppl2 <- c("p_","c_","o_","f_","g_","s_")
  for(i in 1:length(mic.vec)){
    mic.vec[[i]][[2]] <- gsub(paste0("^",lvlppl[lvlidx[i]]),"",    mic.vec[[i]][[1]])
    mic.vec[[i]][[2]] <- gsub(paste0("^",lvlppl2[lvlidx[i]]),"",    mic.vec[[i]][[2]])
    mic.vec[[i]][[2]] <- str_trim(mic.vec[[i]][[2]],side="both")
    # mic.vec[[i]][[2]] <- gsub("_"," ",mic.vec[[i]][[2]])
    # mic.vec[[i]][[2]] <- gsub("\\.","",mic.vec[[i]][[2]])
  }
  
  
  if(netModel=="gem"){
    if(taxalvl=="all"){
      mic.map <- list()
      for(i in 1:length(mic.vec)){
        mic.map[[i]] <- doGemNameMatch(mic.vec[[i]],lvlidx[i],predDB)
      }
      names(mic.map)<-lvlnm[lvlidx]
      # map_num <- orig.num <-setNames(rep(0,6),lvlnm)
      # orig.num[lvlidx] <- unlist(lapply(mic.map,function(x) length(unique(x$Query))))
      # map_num[lvlidx] <- unlist(lapply(mic.map,function(x) length(which(!(is.na(x$Match))))))
      # 
    }
    
  }else if(netModel=="keggNet"){
    if(taxalvl=="default"){
      taxalvl = names(mic.vec)[length(mic.vec)]
    }
    
    if(file.exists("kegg.mic.map.qs")){
      mic.map = qs::qread("kegg.mic.map.qs")
    }else{
      mic.map <- list()
    }
    if(is.null(mic.map[[taxalvl]]) | length(mic.map[[taxalvl]])==0){
      mic.map <- doKeggNameMatch(mic.vec[[taxalvl]],taxalvl)
      
    }
    if(any(c(is.null(mic.map), length(mic.map)==0))){
      current.msg<<-paste0("No ",taxalvl, " was found in kegg database!")
      return(0)
    }
    sig.mic <- phyloseq_objs$sigfeat[[taxalvl]]
    sig.mic <- mic.map$Match[match(sig.mic,mic.map$Query)][!is.na( mic.map$Match)]
    sig.mic<<-unlist(strsplit(sig.mic,split=";"))
    
    
    if(any(c(is.null(sig.mic), length(sig.mic)==0))){
      current.msg<<-paste0("No significant ",taxalvl, " was found in kegg database! Taxonomy level can be change on comparison analysis page!")
    }
  }
  
  
  qs::qsave(mic.map,paste0(netModel,".mic.map.qs"))
  return(1)
}


doGemNameMatch <- function(qvec,l,predDB){
  
  taxMapLong <- qs::qread(paste0(lib.path.mmp,predDB,"_tax.qs"))[[l]]
  names(taxMapLong)[1] <- "taxa"
  res <- data.frame(Query=qvec[[1]],Qtrans=qvec[[2]],stringsAsFactors = F)
  taxMapLong$taxa2<-  gsub("[[:space:]./_-]", "_",taxMapLong$taxa)
  taxMapLong$taxa2<-  gsub("\\[|\\]","",taxMapLong$taxa2)
  res$Match <- taxMapLong[match(tolower(res$Qtrans),tolower(taxMapLong$taxa2)),1]
  fast.write(res, paste("gem_taxa_match_result.csv"));
  return(res)
}

doKeggNameMatch <- function(qvec,taxalvl){
  taxalvl = tolower(taxalvl)
  taxalvl<<- taxalvl
  taxMapKEGG <- qs::qread(paste0(lib.path.mmp,"taxMapKEGG.qs"))[[taxalvl]]
  taxnms <- gsub("[[:space:]./_-]", "_",names(taxMapKEGG)[-1])
  taxnms<-  gsub("\\[|\\]","",taxnms)
  names(taxnms) <- names(taxMapKEGG)[-1]
  res <- data.frame(Query=qvec[[1]],Qtrans=qvec[[2]],stringsAsFactors = F)
  nmsidx= which(taxnms %in% res$Qtrans)
  mtchidx <-  taxMapKEGG[names(taxnms)[nmsidx]]
  mtcls <<- unique(unlist(mtchidx))
  mtchidx <- unlist(lapply(mtchidx, function(x) paste(unique(x),collapse = ";")))
  res$Match <- mtchidx[match(res$Qtrans,as.character(taxnms[names(mtchidx)]))]
  fast.write(res, paste("kegg_taxa_match_result.csv"));
  message("kegg taxonomy mapping done!")
  return(res)
}


CreatPathwayLib <- function(contain){
  #for LTS
  mbSetObj <- .get.mbSetObj(mbSet);
  if(!is.null(mbSetObj$paramSet$includeInfoFileNm)){
    includeInfoNm <- mbSetObj$paramSet$includeInfoFileNm;
  }else{
    includeInfoNm <- "includeInfo";
  }

  if(contain=="usrbac"){
    mtcls = mtcls
  }else if(contain=="sigbac"){
    mtcls = sig.mic[sig.mic!="NA"]
  }
  bacpath <- qs::qread(paste0(lib.path.mmp,"bacpathway.met.qs"))[mtcls]
  bacpath <- bacpath[!is.na(names(bacpath))]
  paths <- unique(unlist(lapply(bacpath,function(x) names(x))))
  
  current.lib = vector("list",length=length(paths))
  names(current.lib) = paths
  for(p in paths){
    pth = lapply(bacpath, function(x) x[[p]])
    current.lib[[p]] =unique(unlist(pth))
  }
  
  includeInfo = list(nodes=unique(unlist(current.lib)))
  edges.bc = qs::qread(paste0(lib.path.mmp,"edge.bac.qs"))
  edges  = data.frame(edge=edges.bc$id_rxn,cpd = edges.bc$met)
  edges = edges[which(edges$cpd %in% includeInfo$nodes),]
  edges = edges.bc[which(edges.bc$id_rxn %in% edges$edge),]
  
  includeInfo$edges = edges
  
  qs::qsave(current.lib,paste0(taxalvl,".current.lib.qs"))
  
  json.mat <- rjson::toJSON(includeInfo);
  sink(paste0(includeInfoNm, ".json"));
  cat(json.mat);
  sink();
  
}



M2Mprediction<- function(model,predDB,taxalvl,psc=0.5,metType="metabolite"){
  
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("prescale.phyobj.qs")
  }
  
  if(predDB=="null"| is.null(predDB) | predDB==""){
    predDB <- "agora"
  }
  current.proc$predDB <<-predDB
  
  if(is.null(taxalvl)|taxalvl=="null"){
    taxalvl=names(phyloseq_objs$count_tables)[length(phyloseq_objs$count_tables)-1]
  }
  
  if(taxalvl=="all"){
    lvlnm <-  names(phyloseq_objs$count_tables)
    lvlnm <- lvlnm[lvlnm!="OTU"]
    #taxalvls <- lvlnm[length(lvlnm)]
    predres <- vector('list',length=length(lvlnm))
    names(predres)<- lvlnm
    for(taxalvl in lvlnm){
      OTUtab <<- phyloseq_objs$count_tables[[taxalvl]]
      predres[[taxalvl]] <- doGemPrediction(predDB,taxalvl,psc,metType)
    }
  }else{
    
    OTUtab <<- phyloseq_objs$count_tables[[taxalvl]]
    predres <- doGemPrediction(predDB,taxalvl,psc,metType)
    current.proc$taxalvl <<-taxalvl
    current.proc$predres<<-predres
    
  }
  # met.map <- qs::qread("met.map.qs")
  # predres <-predres[rownames(predres) %in% met.map$Match,]
  # mbSetObj$analSet$m2m.pred <- predres
  #mbSetObj$imgSet$m2m.pred <- imgName;
  
  qs::qsave(predres,paste0("m2m_pred_",predDB,".qs"))
  
  message("Prediction completed!")
  return(1)
}

doGemPrediction <- function(predDB,taxalvl,psc=0.5,metType,matchonly=T,sigonly=T){
  #print(c(predDB,taxalvl,metType))
  require(reshape2)
  message('Loading the model database..')
  psc <- as.numeric(psc)
  taxalvl<-tolower(taxalvl) 
  tax_map <- qs::qread("gem.mic.map.qs")[[taxalvl]]
  tax_map <- tax_map[which(!is.na(tax_map$Match)),]
  m2m_ls <- qs::qread(paste0(lib.path.mmp,predDB,".qs"))[[taxalvl]]
  names(m2m_ls)[1] <- "taxa"
  m2m_ls <- m2m_ls[which(m2m_ls$potential>=psc),]
  m2m_ls <- m2m_ls[which(m2m_ls$taxa %in% tax_map$Match),]
  m2m_ls$taxa <- tax_map$Query[match(m2m_ls$taxa,tax_map$Match)]
  
  if(metType=="metabolite"){
    met.map<- current.proc$gem
    m2m_ls <- m2m_ls[which(m2m_ls$metID %in% met.map$Match),]
    m2m_ls$metabolite <- met.map$Query[match(m2m_ls$metID,met.map$Match)]
  }
  
  m2m_db <- dcast(m2m_ls,taxa~metabolite,value.var="potential")
  
  m2m_db[is.na(m2m_db)] <- 0
  dbnorm <- as.matrix(m2m_db[,-1])
  ##filter otu table
  OTUtab <- OTUtab[which(rownames(OTUtab) %in% tax_map$Query),]
  if(!(all( rownames(OTUtab) ==m2m_db$taxa))){
    AddErrMsg("Names not match!");
    return(0);
  }
  OTUtab <- apply(OTUtab, 2, function(x) ReScale(rank(x),0,1) )
  fun_prediction = NULL
  fun_m2m_pair <- list()
  message('Generating metabolic profile..')
  
  rownames(dbnorm) <- m2m_db$taxa
  for(sample in 1:ncol(OTUtab)){
    fun_prediction_sample = dbnorm * as.numeric(OTUtab[,sample])
    #fun_prediction_sample <- t(preprocessCore::normalize.quantiles(t(fun_prediction_sample), copy=FALSE))
    #??zero should be back transfer??
    fun_m2m_pair[[sample]] <- fun_prediction_sample
    fun_prediction_sample = colMeans(fun_prediction_sample)
    fun_prediction_sample = fun_prediction_sample / sum(fun_prediction_sample)
    if(is.na(sum(fun_prediction_sample))) fun_prediction_sample[1:ncol(dbnorm)] = 0
    fun_prediction = cbind(fun_prediction, fun_prediction_sample)
  }
  
  message("Prediction done")
  
  names(fun_m2m_pair) <-colnames(fun_prediction) <-  colnames(OTUtab)
  fun_m2m_pair <- lapply(fun_m2m_pair, function(p){
    rownames(p)=rownames(OTUtab)
    return(p) })
  
  keep = which(rowSums(fun_prediction) > 0)
  
  if (length(keep) == 0) stop("No functional prediction possible!\nEither no nearest neighbor found or your table is empty!")
  fun_prediction_final = fun_prediction[unname(keep),]
  
  fast.write(fun_prediction_final, paste0(taxalvl,"_prediction.csv"))
  return(list(fun_prediction_sample=fun_prediction_final,fun_prediction_met=fun_m2m_pair))
  
}



###########################################################
####################Prediction && correlation heatmap######
###########################################################
###########################################################

DoM2Mcorr <- function(mic.sig,met.sig,cor.method="univariate",cor.stat="pearson",taxalvl){
  labels <- c(rownames(mic.sig),rownames(met.sig))
  nan.msg<<-"null"
  if(cor.method == "univariate"){
    require(psych)
    res <- corr.test(cbind(t(mic.sig), t(met.sig)),method=cor.stat);
    rowidx <- which(rownames(res$r) %in% rownames(mic.sig))
    colidx <- which(colnames(res$r) %in% rownames(met.sig))
    corr.mat <- res$r[rowidx,colidx];
    corr.pval <- res$p[rowidx,colidx]
    
  }else if(cor.method == "MI"){
    library(parmigene)
    res = knnmi.all(rbind(mic.sig, met.sig), k=5)
    scale = 1/max(res)
    corr.mat = res * scale
    
  }else if(cor.method=="discor"){
    require(energy)
    corr.mat <- matrix(data=NA,nrow=nrow(mic.sig),ncol = nrow(met.sig))
    corr.pval <- matrix(data=NA,nrow=nrow(mic.sig),ncol = nrow(met.sig))
    for(row in 1:nrow(mic.sig)){
      res<-lapply(1:nrow(met.sig), function(y) {
        corr=dcor.test(mic.sig[row,],met.sig[y,],R=100)
        return(corr)
      })
      corr.mat[row,] <- unlist(lapply(res,function(z) return(z[["statistic"]])))
      corr.pval[row,] <- unlist(lapply(res,function(z) return(z[["p.value"]])))
    }
    colnames(corr.mat) <- colnames(corr.pval) <- rownames(met.sig)
    rownames(corr.mat) <- rownames(corr.pval) <- rownames(mic.sig)   
  }else{
    library(ppcor);
    sel.res <-  cbind(t(mic.sig), t(met.sig))
    
    res  <- tryCatch( 
      { pcor(sel.res, method=cor.stat);
      },
      error = function(error_cond){
        current.msg <<-"Fail to perform patial correlation";
      })
    
    if(!is.null(dim(res$estimate))){
      corr.mat <- res$estimate;
      corr.pval <- res$p.value;
      if(any(is.nan(corr.pval))){
        corr.pval=0
        nan.msg <<-"NaNs produced in p_value calculation using current correlation parameters! ";
        rownames(corr.mat) <-  colnames(sel.res)
        colnames(corr.mat) <- colnames(sel.res)
        rowidx <- which(rownames(corr.mat) %in% rownames(mic.sig))
        colidx <- which(colnames(corr.mat) %in% rownames(met.sig))
        corr.mat <- corr.mat[rowidx,colidx];
        
      }else{
        
        rownames(corr.mat) <- rownames(corr.pval) <- colnames(sel.res)
        colnames(corr.mat) <-colnames(corr.pval) <- colnames(sel.res)
        rowidx <- which(rownames(corr.mat) %in% rownames(mic.sig))
        colidx <- which(colnames(corr.mat) %in% rownames(met.sig))
        corr.mat <- corr.mat[rowidx,colidx];
        corr.pval <- corr.pval[rowidx,colidx]
        
      }
    }else{
      corr.mat=0
      corr.pval=0
    }
  }
  return(list(corr.mat=corr.mat,corr.pval=corr.pval));
}

performeCorrelation <- function(mbSetObj,taxalvl,initDE,cor.method="univariate",cor.stat="pearson",sign, cor.thresh=0.5,
                                corp.thresh=0.05){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(!exists("phyloseq_objs")){
    phyloseq_objs <- qs::qread("phyloseq_objs.qs")
  }
  
  micdat <- phyloseq_objs$count_tables[[taxalvl]]
  metdat <- current.proc$met$data.proc
  if(micDataType=="ko"){
    lbl.mic <-current.proc$mic$sigfeat
  }else{
    lbl.mic <- phyloseq_objs$sigfeat[[taxalvl]] 
  }
  lbl.met <- current.proc$met$sigfeat
  if(length(lbl.mic) >100){lbl.mic= lbl.mic[1:100]}
  if(length(lbl.met) >100){lbl.met= lbl.met[1:100]}
  
  mic.sig <- micdat[which(rownames(micdat) %in% lbl.mic),]
  met.sig <- metdat[which(rownames(metdat) %in% lbl.met),match(colnames(mic.sig),colnames(metdat)), drop = FALSE]
  res.corr <- DoM2Mcorr(mic.sig,met.sig,cor.method,cor.stat,taxalvl)
  corr.mat <- res.corr$corr.mat
  if(is.null(dim(res.corr$corr.mat))){
    corr.mat <- 0
    corr.pval <- 0
    return(0)
  }
  res.corr.filt <- doCorrelationFilt(res.corr,cor.thresh,corp.thresh,sign)
  if(is.null(dim(res.corr.filt$corr.mat))){
    corr.mat <- 0
    corr.pval <- 0
    return(0)
  }
  output.dat <- reshape2::melt(res.corr$corr.mat,value.name = "correlation")
  output.p <- reshape2::melt(res.corr$corr.pval,value.name = "pval")
  if(!is.null(output.p)&nrow(output.p)>0){
    output.dat <- merge(output.dat,output.p)
    output.dat <- output.dat[order(output.dat$pval),]
  }else{
    output.dat <- output.dat[order(abs(output.dat$correlation),decreasing = T),]
  }
  fast.write(output.dat, file=paste("correlation", "_",cor.method,"_",cor.stat,".csv", sep=""),row.names=F);
  corrNm<<-paste("correlation", "_",cor.method,"_",cor.stat,".csv", sep="")
  
  mbSetObj$analSet$corr.method <- paste0(cor.method,"_",cor.stat)
  current.proc$corr.mat <<- res.corr.filt$corr.mat
  current.proc$corr.pval <<- res.corr.filt$corr.pval
  message("correlation completed")
  return(.set.mbSetObj(mbSetObj));
}


doCorrelationFilt <- function( res.corr,cor.thresh,corp.thresh,sign){
  
  corr.mat <- res.corr$corr.mat
  corr.pval <- res.corr$corr.pval
  if(any(is.nan(corr.pval))|is.null(dim(corr.pval))){
    corr.pval=0
    nan.msg <<-"NaNs produced in p_value calculation using current correlation parameters! ";
  }else{
    keepidx1p <- apply(corr.pval,2,function(x) sum(x<corp.thresh))>0
    keepidx2p <- apply(corr.pval[,keepidx1p],1,function(x) sum(x<corp.thresh))>0
    corr.pval <- corr.pval[keepidx2p,keepidx1p]; if(length(corr.pval) ==1) { corr.pval <- matrix(corr.pval) }
    corr.mat <- corr.mat[keepidx2p,keepidx1p]
  }
  
  if(sign=="positive"){
    keepidx1 <- apply(corr.mat,2,function(x) sum(x>cor.thresh))>0
    keepidx2 <- apply(corr.mat[,keepidx1],1,function(x) sum(x>cor.thresh))>0
  }else if(sign=="negative"){
    keepidx1 <- apply(corr.mat,2,function(x) sum(x<(-cor.thresh)))>0
    keepidx2 <- apply(corr.mat[,keepidx1],1,function(x) sum(x<(-cor.thresh)))>0
  }else{
    keepidx1 <- apply(corr.mat,2,function(x) sum(abs(x)>cor.thresh))>0
    keepidx2 <- apply(corr.mat[,keepidx1],1,function(x) sum(abs(x)>cor.thresh))>0
  }
  
  corr.mat <- corr.mat[keepidx2,keepidx1]
  if(!is.null(dim(corr.pval))){
    
    corr.pval <- corr.pval[keepidx2,keepidx1]
  }
  
  return(list(corr.mat=corr.mat,corr.pval=corr.pval))
  
}


CreatM2MHeatmap<-function(mbSetObj,htMode,overlay, taxalvl, plotNm,  format="png", 
                          smplDist="euclidean", clstDist="ward.D", palette="npj",viewOpt="barraw", 
                          clustRow="T", clustCol="T", 
                          colname="T",rowname="T", fontsize_col=10, fontsize_row=10,
                          sign, cor.thresh=0.5,corp.thresh=0.05,
                          potential.thresh=0.5,predpval.thresh=0.05,
                          var.inx=NA, border=T, width=NA, dpi=72){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  load_iheatmapr();
  load_rcolorbrewer();
  load_viridis();
  current.msg<<- NULL
  set.seed(2805614);

  #used for color pallete
  ######set up plot
  #colors for heatmap
  if(palette=="gbr"){
    colors <- grDevices::colorRampPalette(c("green", "black", "red"), space="rgb")(256);
  }else if(palette == "heat"){
    colors <- grDevices::heat.colors(256);
  }else if(palette == "topo"){
    colors <- grDevices::topo.colors(256);
  }else if(palette == "gray"){
    colors <- grDevices::colorRampPalette(c("grey90", "grey10"), space="rgb")(256);
  }else if(palette == "byr"){
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256));
  }else if(palette == "viridis") {
    colors <- rev(viridis::viridis(10))
  }else if(palette == "plasma") {
    colors <- rev(viridis::plasma(10))
  }else  if(palette == "npj"){
    colors <- c("#00A087FF","white","#E64B35FF")
  }else  if(palette == "aaas"){
    colors <- c("#4DBBD5FF","white","#E64B35FF");
  }else  if(palette == "d3"){
    colors <- c("#2CA02CFF","white","#FF7F0EFF");
  }else {
    colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")), alpha=0.8)(100)
    #c("#0571b0","#92c5de","white","#f4a582","#ca0020");
  }
  
  plotjs = paste0(plotNm, ".json");
  plotwidget <- paste0(plotNm, ".rda")
  plotNm = paste(plotNm, ".", format, sep="");
  mbSetObj$imgSet$IntegrationHeatmap<-plotNm;
  if(htMode=="predht"){  ####using  prediction pair pval
    pred.dat <- current.proc$pred.dat
    predDE <- current.proc$predDE
    
    
    data.abd <- data.frame(mic=as.character(predDE$mic[match(pred.dat$pair,rownames(predDE))]),
                           met=as.character(predDE$met[match(pred.dat$pair,rownames(predDE))]),
                           var = rowMeans(pred.dat[,-1]),                        
                           value = predDE$P_value)
    
    data.abd <- data.abd[order(data.abd$value,-(data.abd$var)),]
    
    if(length(unique(data.abd$mic))>100){
      micnms <- unique(data.abd$mic)[1:100]
    }else{
      micnms <- unique(data.abd$mic)
    }
    
    if(length(unique(data.abd$met))>100){
      metnms <- unique(data.abd$met)[1:100]
    }else{
      metnms <- unique(data.abd$met)
    }
    
    data <- data.abd[which(data.abd$mic %in% micnms & data.abd$met %in% metnms),-4]
    
    data <- reshape2::dcast(data,mic~met)
    data[is.na(data)] <-0
    data.mtr <- data[,-1]
    micnms <- data$mic
    metnms <- colnames(data.mtr)
    nameHt <- "Ave.Potential"
    
    anno.mat0 <- data.abd[which(data.abd$mic %in% micnms & data.abd$met %in% metnms),-3]
    
    anno.mat0 <- anno.mat0[which(anno.mat0$value<predpval.thresh),]
    if(nrow(anno.mat0)==0){   
      current.msg <<- paste("No significant prediction was detected using current parameters!");
      
    }else{
      #anno.mat$value <- as.character(round(anno.mat$value,2))
      names(anno.mat0) <- c("Var1","Var2","value")
    }
    
    if(overlay=="true"){
      corr.mat <- current.proc$corr.mat
      corr.pval <- current.proc$corr.pval
      if(nrow(corr.mat)==0){
        current.msg <<- paste("No statistical correlation was detected using current parameters!");
      }else{
        corr.mat <- corr.mat[which(rownames(corr.mat) %in% as.character(micnms)),
                             which(colnames(corr.mat) %in% metnms)]
        anno.mat <- reshape2::melt(corr.mat,value.name = "correlation")
        anno.mat <- anno.mat[which(anno.mat$correlation>cor.thresh),]
        if(is.null(anno.mat) | nrow(anno.mat)==0){
          current.msg <<- paste("No significant statistical correlation was detected using current parameters!");          
        }else{
          if(is.null(dim(corr.pval))){            
            current.msg <<- paste("No significant statistical correlation was detected! The triangle only show the ones pass the correlation thresh hold!");            
          }else{
            corr.pval <- corr.pval[which(rownames(corr.pval) %in% as.character(micnms)),
                                   which(colnames(corr.pval) %in% metnms)]
            anno.pval <- reshape2::melt(corr.pval,value.name="pval")
            anno.pval <- anno.pval[which(anno.pval$pval<corp.thresh),]
            if(nrow(anno.pval)==0| is.null(anno.pval)){
              current.msg <<- paste("No statistical correlation pass the significance thresh hold using current parameters! The triangle show the ones pass the correlation thresh hold!");
              
            }else{
              anno.mat <- unique(left_join(anno.mat,anno.pval))
              if(all(is.na(anno.mat$pval))){
                anno.mat$pval<- NULL
                current.msg <<- paste("No statistical correlation pass the significance thresh hold using current parameters! The triangle show the ones pass the correlation thresh hold!");
              }
            }          
          }
          anno.mat <- unique(left_join(anno.mat0,anno.mat))
          if(all(is.na(anno.mat$correlation))){
            anno.mat$correlation<- NULL
            current.msg <<- paste("No statistical correlation was detected using current parameters");
          }
          anno.mat$size <- as.numeric(ReScale(-log(anno.mat$value),8,12))
          annols <- vector("list",length=nrow( anno.mat))
        }   
      }  
      mbSetObj$analSet$integration$corr<- cor.thresh
      mbSetObj$analSet$integration$corrPval<- corp.thresh
    }else{
      anno.mat <- anno.mat0
      anno.mat$size <- as.numeric(ReScale(-log(anno.mat$value),8,12))
      annols <- vector("list",length=nrow( anno.mat))
    }    
    mbSetObj$analSet$integration$potential<- potential.thresh
    mbSetObj$analSet$integration$predPval<- predpval.thresh
  }else if(htMode=="corrht"){ 
    data.mtr <- current.proc$corr.mat
    corr.pval <- current.proc$corr.pval
    
    if(nrow(data.mtr)==0){      
      current.msg <<- paste("No statistical correlation was detected using current parameters!");
      return(0)      
    }else{
      micnms <- rownames(data.mtr)
      metnms <- colnames(data.mtr)
      nameHt <- "Correlation"
      anno.mat0 <- reshape2::melt(data.mtr)
      if(is.null(dim(corr.pval))){
        
        current.msg <<- paste("No significant correlation was detected using current parameters!");
      }else{
        
        #### fro annotation using pval
        corr.pval <- corr.pval[which(rownames(corr.pval) %in% as.character(micnms)),
                               which(colnames(corr.pval) %in% metnms)]
        
        
        if(sign=="positive"){
          anno.mat0 <- anno.mat0[which(anno.mat0$value>cor.thresh),] 
        }else if(sign=="negative"){
          anno.mat0 <- anno.mat0[which(anno.mat0$value< (-cor.thresh)),] 
        }else{
          anno.mat0 <- anno.mat0[which(abs(anno.mat0$value)>cor.thresh),] 
        }
        
        anno.pval <- reshape2::melt(corr.pval,value.name = "pval")
        anno.pval <- anno.pval[which(anno.pval$pval<corp.thresh),]
        anno.mat <- unique(left_join(anno.mat0,anno.pval))
        anno.mat <- anno.mat[!(is.na(anno.mat$pval)),]
        if(nrow(anno.mat)==0){
          current.msg <<- paste("No significant correlation was detected using current parameters!");          
        }else{
          anno.mat$size <- as.numeric(ReScale(-log(anno.mat$pval),8,12))
          annols <- vector("list",length=nrow( anno.mat))
        }
        
      }
      
      # anno.mat$value <- as.character(round(anno.mat$value,2))      
    }
    
    if(overlay=="true"){
      pred.de <- current.proc$predDE[,c(5,6,1)] %>% filter(P_value <predpval.thresh)
      rownames(pred.de) <- NULL
      names(pred.de)[1:2] <- names(anno.mat0)[1:2]
      if(exists("anno.mat")){
        if(nrow(anno.mat)>0){
          anno.mat <- unique(left_join(anno.mat,pred.de))
        }else{
          anno.mat <- anno.mat0
          anno.mat <- unique(left_join(anno.mat,pred.de)) %>% filter(!(is.na(P_value)))
        }        
      }else{
        anno.mat <- anno.mat0
        anno.mat <- unique(left_join(anno.mat,pred.de)) %>% filter(!(is.na(P_value)))        
      }
      if(nrow(anno.mat)==0){
        if(!is.null(current.msg)){
          current.msg <<- c(current.msg," No overlay prediction result was detected using current parameters!")
        }else{
          current.msg <<- paste("No overlay prediction result was detected using current parameters!");          
        }        
      }else{
        anno.mat$size <- as.numeric(ReScale(-log(anno.mat$pval),8,12))
        annols <- vector("list",length=nrow( anno.mat))
      }
      mbSetObj$analSet$integration$potential<- potential.thresh
      mbSetObj$analSet$integration$predPval<- predpval.thresh
    }
    mbSetObj$analSet$integration$corr<- cor.thresh
    mbSetObj$analSet$integration$corrPval<- corp.thresh
  }
  
  data1 <- data.mtr;
  data1sc <- as.matrix(apply(data1, 2, as.numeric))
  rownames(data1sc) <- micnms
  #data1sc <- scale_mat(data1sc, scaleOpt)
  
  fzCol <- round(as.numeric(fontsize_col), 1)
  fzRow <- round(as.numeric(fontsize_row), 1)
  map.height=nrow(data1)*30
  map.width=ncol(data1)*30
  
  #cb_grid <- setup_colorbar_grid(nrow = 100, x_start = 1.1, y_start = 0.95, x_spacing = 0.15)
  dend_row <- hclust(dist(data1sc, method = smplDist), method = clstDist)
  p <- iheatmap(data1sc,
                # colorbar_grid = cb_grid, 
                name = nameHt, x_categorical = TRUE,
                layout = list(font = list(size = 10)),
                colors = colors
  )
  
  if (clustRow == "true") {
    p <- p %>% add_row_dendro(dend_row, side = "right")
  }
  
  if (colname == "true" ){    
    p <- p %>%  add_col_labels(size = 0.1, font = list(size = fzCol))
  }
  
  if (colname == "true" ){    
    p <- p %>% add_row_labels(size = 0.1, font = list(size = fzRow), side = "left")
  }
    
  if (clustCol == "true") {
    dend_col <- hclust(dist(t(data1), method = smplDist), method = clstDist)
    p <- p %>% add_col_dendro(dend_col)
  }
  
  pwidget <- to_widget(p)
  save(pwidget, file = plotwidget)
  
  as_list <- to_plotly_list(p)
  ### add the layer for  annotation
  if(exists("annols")){    
    annht <- as_list$layout$annotations
    annht <- data.frame(label=unlist(lapply(annht,function(x) x[["text"]])),
                        X= unlist(lapply(annht,function(x)  x[["x"]])),
                        Y=unlist(lapply(annht,function(x)  x[["y"]])))
    for(i in 1:nrow(anno.mat)){
      annols[[i]]$text <- "*"
      annols[[i]]$x <- annht$X[match(anno.mat$Var2[i],annht$label)]
      annols[[i]]$y <- annht$Y[match(anno.mat$Var1[i],annht$label)]
      annols[[i]][["font"]][["size"]] <- anno.mat$size[i]
      annols[[i]]$showarrow <- FALSE
    }
    
    if(htMode=="predht"&overlay=="true"){
      if(!(is.null(anno.mat$pval))){
        anno.mat$pval[is.na(anno.mat$pval)]=1
        anno.mat$size2 <- as.numeric(ReScale(-log(anno.mat$pval),10,14))
        for(i in 1:nrow(anno.mat)){
          if(anno.mat$pval[i]!=1){
            annols[[i]]$text <- "۞"
            annols[[i]][["font"]][["size"]] <- anno.mat$size2[i]
          }       
        } 
      }else{
        if(!(is.null(anno.mat$correlation))){
          anno.mat$correlation[is.na(anno.mat$correlation)]=0
          anno.mat$size2 <- as.numeric(ReScale(anno.mat$correlation,10,14))           
          for(i in 1:nrow(anno.mat)){
            if(anno.mat$correlation[i]!=0){
              annols[[i]]$text <- "۞"
              annols[[i]][["font"]][["size"]] <- anno.mat$size2[i]
            }       
          } 
        }
      }
    }
    
    if(!(is.null(anno.mat$P_value))&htMode=="corrht"&overlay=="true"){
      anno.mat$P_value[is.na(anno.mat$P_value)]=1
      anno.mat$size2 <- as.numeric(ReScale(-log(anno.mat$P_value),10,14))  
      for(i in 1:nrow(anno.mat)){
        if(anno.mat$P_value[i]!=1){
          annols[[i]]$text <- "۞"
          annols[[i]][["font"]][["size"]] <- 12
        }       
      } 
    }

    as_list$layout$annotations <- c(as_list$layout$annotations,annols)
  }
  
  overlyNum = length(which(unlist(lapply(annols,function(x) x[["text"]]=="۞"))))
  if (viewOpt != "overview") {
    as_list[["layout"]][["width"]] <- max(map.width,1000)
    as_list[["layout"]][["height"]] <- max(map.height,800)
  } else {
    as_list[["layout"]][["width"]] <- 1200
    as_list[["layout"]][["height"]] <- map.height
  }
  
  if(exists("id2nm",where=current.proc)){
    for(i in 1:ncol(data1sc)){
      
      as_list$layout$annotations[[i]]$text = unname(current.proc$id2nm[as_list$layout$annotations[[i]]$text])
    }
  }
  
  as_json <- attr(as_list, "TOJSON_FUNC")(as_list)
  as_json <- paste0("{ \"x\":", as_json, ",\"evals\": [],\"jsHooks\": []}")
  
  write(as_json, plotjs)
  
  if(is.null(current.msg)){
    current.msg<<-"null"
  }
  # storing for Report Generation
  mbSetObj$analSet$integration$heatmap <- data1sc
  mbSetObj$analSet$integration$heatmap.dist <- smplDist
  mbSetObj$analSet$integration$heatmap.clust <- clstDist
  mbSetObj$analSet$integration$taxalvl <- taxalvl
  mbSetObj$analSet$integration$overlay <- overlay
  mbSetObj$analSet$integration$htMode <- htMode
  mbSetObj$analSet$integration$sign <- sign
  mbSetObj$imgSet$heatmap_cormmp <- plotwidget
  message("heatmap done")
  .set.mbSetObj(mbSetObj)
  return(overlyNum)
}

###########################################################
####################KEGG Metabolism Network################
###########################################################
###########################################################


PrepareOTUQueryJson <- function(mbSetObj,taxalvl,contain="bac"){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(contain=="bac"| contain=="hsabac"|contain=="all"|contain=="hsa"){
    met.map <- qs::qread("keggNet.met.map.qs")
    query.res <- rep(2,length(unique(met.map$id[!(is.na(met.map$id))])))
    names(query.res) <- unique(met.map$id[!(is.na(met.map$id))])
    
  }else{
    doNetpathFilt(contain)
  }
  json.mat <- rjson::toJSON(query.res);
  sink("network_query_met.json");
  cat(json.mat);
  sink();
  
  return(.set.mbSetObj(mbSetObj));
  
}



PerformTuneEnrichAnalysis <- function(mbSetObj, dataType,category, file.nm,contain="hsabac",enrich.type){
  mbSetObj <- .get.mbSetObj(mbSetObj);

  if(enrich.type == "hyper"){
    if(dataType=="metabolite"){
      mbSetObj <- PerformMetListEnrichment(mbSetObj, contain, file.nm);
    }else{
      MicrobiomeAnalystR:::LoadKEGGKO_lib(category);
      PerformKOEnrichAnalysis_List(mbSetObj, file.nm);
      mbSetObj <- .get.mbSetObj(mbSet);
    }
    
  }else if(enrich.type =="global"){
    if(contain=="usrbac" & micDataType=="ko"){
      tuneKOmap()
      contain = "bac"
    }
    .prepare.global.tune(mbSetObj, dataType, category, file.nm,contain);
    .perform.computing();
    
    if(dataType=="ko"){
      mbSetObj = .save.global.res();
      taxalvl = "ko"
    }else if(dataType=="metabolite"){
      
      mbSetObj=enrich2json()
    }
    
  }else if(enrich.type =="mummichog"){
    
    if(!exists("performPeakEnrich")){ # public web on same user dir
      .load.scripts.on.demand("utils_peak2fun.Rc");    
    }
    mbSetObj=performPeakEnrich(lib=contain);
    mbSetObj <- .get.mbSetObj(mbSet);
  }
  if(!exists("taxalvl")){taxalvl = "ko"}
  mbSetObj$analSet$keggnet$background <- contain
  mbSetObj$analSet$keggnet$taxalvl <- taxalvl
  return(.set.mbSetObj(mbSetObj))
}


.prepare.global.tune<-function(mbSetObj, dataType,category, file.nm,contain){
  
  load_phyloseq();
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  phenotype <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[selected.meta.data]]);
  
  if(dataType=="metabolite"){
    
    if(contain=="bac"){
      current.set <- qs::qread(paste0(lib.path.mmp,"kegg_bac_mummichog.qs"))$pathways$cpds
      
    }else if(contain=="hsabac"){
      current.set <- qs::qread(paste0(lib.path.mmp,"kegg_hsa_bac_mummichog.qs"))$pathways$cpds
      
    }else if(contain=="hsa"){
      current.set <- qs::qread(paste0(lib.path.mmp,"kegg_hsa_mummichog.qs"))$pathways$cpds
      
    }else if(contain=="all"){
      current.set <- qs::qread(paste0(lib.path.mmp,"kegg_all_mummichog.qs"))$pathways$cpds
      
    }else{
      current.set <- qs::qread(paste0(taxalvl,".current.lib.qs"))
      
    }
    

    metmat <-  t(current.proc$met$data.proc)          
    met.map <-  qs::qread("keggNet.met.map.qs")
    met.map <- met.map[!(is.na(met.map$Node)),]
    met.map$include = ifelse(met.map$Match %in% unique(unlist(current.set)),T,F)
    qs::qsave(met.map,"keggNet.met.map.qs")
    colnames(metmat) <- met.map$Match[match(colnames(metmat),met.map$Query)]
    datmat <- metmat[,which(colnames(metmat)!='')]
    
    hits <- lapply(current.set, function(x){x[x %in% colnames(datmat)]});
    set.num <- unlist(lapply(current.set, length), use.names = FALSE);
    dat.in <- list(cls=phenotype, data=datmat, subsets=hits, set.num=set.num, filenm=file.nm);
    
  }else if(dataType=="ko"){
    if(contain=="bac"){
      current.set <- qs::qread(paste0(lib.path.mmp,"ko_set_bac.qs"))
      
    }else if(contain=="hsabac"){
      current.set <- qs::qread(paste0(lib.path.mmp,"ko_set_hsa_bac.qs"))
      
    }else if(contain=="hsa"){
      current.set <- qs::qread(paste0(lib.path.mmp,"ko_set_hsa.qs"))
      
    }else if(contain=="all"){
      kegg.anot <- .read.microbiomeanalyst.lib.rds("ko_pathways.rds", "ko")
      current.setlink <- kegg.anot$link;
      current.set <- kegg.anot$sets$Metabolism;
    }else{
      current.set <- qs::qread(paste0(lib.path.mmp,"ko_set_bac.qs"))
    }
    
    set2nm <-  qs::qread("../../lib/mmp/set2nm.qs")[["pathway"]];
    set.ids <- names(current.set);
    names(set.ids) <- names(current.set)<-  set2nm[set.ids];
    
    current.setids <<-  set.ids;
    
    datmat <- as.data.frame(t(otu_table(mbSetObj$dataSet$norm.phyobj)),check.names=FALSE);
    # first, get the matched entries from current.set
    
    hits <- lapply(current.set, function(x){x[x %in% colnames(datmat)]});
    set.num <- unlist(lapply(current.set, length), use.names = FALSE);
    dat.in <- list(cls=phenotype, data=datmat, subsets=hits, set.num=set.num, filenm=file.nm);
  }
  
  
  my.fun <- function(){
    gt.obj <- globaltest::gt(dat.in$cls, dat.in$data, subsets=dat.in$subsets);
    gt.res <- globaltest::result(gt.obj);
    
    match.num <- gt.res[,5];
    
    if(sum(match.num>0)==0){
      return(NA);
    }
    raw.p <- gt.res[,1];
    
    # add adjust p values
    bonf.p <- p.adjust(raw.p, "holm");
    fdr.p <- p.adjust(raw.p, "fdr");
    
    res.mat <- cbind(set.num, match.num, gt.res[,2], gt.res[,3], raw.p, bonf.p, fdr.p);
    rownames(res.mat) <- names(hits);
    colnames(res.mat) <- c("Size", "Hits", "Statistic Q", "Expected Q", "Pval", "Holm p", "FDR");
    hit.inx <- res.mat[,2]>0;
    res.mat <- res.mat[hit.inx, ];
    ord.inx <- order(res.mat[,5]);
    res.mat <- res.mat[ord.inx,];
    
    return(res.mat);
  }
  
  
  dat.in <- list(cls=phenotype, data=datmat, subsets=hits, set.num=set.num, filenm=file.nm , my.fun=my.fun);
  
  qs::qsave(dat.in, file="dat.in.qs");
  return(1);
}

enrich2json <- function(){
  dat.in <- qs::qread("dat.in.qs"); 
  hits = dat.in$subsets
  file.nm = dat.in$filenm;
  my.res <- dat.in$my.res;
  if(all(c(length(my.res)==1, is.na(my.res)))){
    AddErrMsg("No match was found to the selected metabolite set library!");
    return(0);
  }
  nms <- rownames(my.res);
  hits <- hits[nms];
  resTable <- data.frame(Pathway=rownames(my.res), my.res,check.names=FALSE);
  current.msg <<- "Functional enrichment analysis was completed";
  met.map <-  qs::qread("keggNet.met.map.qs")
  hits.met <- lapply(hits, function(x){
    x=met.map$Name[match(x,met.map$Match)]  
    return(x)
  })
  hits.node <- lapply(hits, function(x){
    x=met.map$Node[match(x,met.map$Match)]  
    return(x)
  })
  # write json 
  path.pval = resTable$Pval; if(length(path.pval) ==1) { path.pval <- matrix(path.pval) };
  hit.num = paste0(resTable$Hits,"/",resTable$Size); if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  path.nms <- resTable$Pathway; if(length(path.nms) ==1) { path.nms <- matrix(path.nms) };
  hits.query <- convert2JsonList(hits)
  path.fdr <- resTable$FDR; if(length(path.fdr) ==1) { path.fdr <- matrix(path.fdr) };
  sig.path <- resTable$Pathway[which(resTable$FDR<0.05)];
  sig.path <- hits.node[names(hits.node) %in% sig.path]; if(length(sig.path) ==0) { sig.path <-0};
  expr.mat = current.proc$met$res_deAnal
  expr.mat$kegg = met.map$Match[match(rownames(expr.mat),met.map$Query)]
  expr.mat <- expr.mat[expr.mat$kegg %in% met.map$Match[met.map$include],]
  expr.mat <- split(expr.mat$T.Stats,expr.mat$kegg)
  json.res <- list(hits.query = hits.query,
                   path.nms = path.nms,
                   path.pval = path.pval,
                   hit.num = hit.num,
                   path.fdr = path.fdr,
                   hits.query.nm = convert2JsonList(hits.met),
                   hits.node = convert2JsonList(hits.node),
                   expr.mat=convert2JsonList(expr.mat),
                   sig.path=sig.path );
  
  json.mat <- RJSONIO::toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();

  mbSetObj <- .get.mbSetObj(NA);
  
  mbSetObj <- recordEnrTable(mbSetObj, "mmp", resTable, "KEGG", "Global Test");

  # write csv
  fast.write(resTable, file=paste(file.nm, ".csv", sep=""), row.names=F);
  return(mbSetObj)
}


GetAssociationPlot <- function(type,keggid,koid,micDataType,metIDType,taxalvl,imgNm,corrMethod,corrSign,corrThresh,corrPval,topNum=10,subset="false"){
  current.msg<<-"null"
  topNum = as.numeric(topNum)
  current.proc$keggmap$corrMethod <<-corrMethod
  current.proc$keggmap$corrSign <<-corrSign
  current.proc$keggmap$corrThresh <<-corrThresh
  current.proc$keggmap$corrPval <<-corrPval
  
  if(type=="corr"){
    corrThresh <- as.numeric(corrThresh)
    corrPval <- as.numeric(corrPval)
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    }
    if(!(keggid %in% current.proc$keggNet$Match)){
      current.msg <<- "The selected compound is not provided in your input data! Please choose the related compounds!"
      return(0)
    }
    if(metDataType == "peak"){
      qvecidx =  which(current.proc$keggNet$Match==keggid)
      adduct = current.proc$keggNet$adduct[qvecidx]
      if(length(interaction(adduct,primary_ions))>0){
        qvec = current.proc$keggNet$Query[which(current.proc$keggNet$Match==keggid)][which(adduct %in% primary_ions)]
      }else{
        qvec =  current.proc$keggNet$Query[which(current.proc$keggNet$Match==keggid)] 
      }
      print(paste0(length(qvec)," peaks related! The plots take longer time."))
    }else{
      qvec =  current.proc$keggNet$Query[which(current.proc$keggNet$Match==keggid)] 
      
    }
    
    
    current.proc$keggmap$current_qvec<<-qvec
    
    micdat <- phyloseq_objs$count_tables[[taxalvl]]
    
    metdat <- current.proc$met$data.proc[qvec,,drop=FALSE]
    if(grepl("u-",corrMethod)){
      cor.method = "univariate"
      cor.stat = gsub("u-","",corrMethod)
      
    }else if(grepl("p-",corrMethod)){
      cor.method = "partial"
      cor.stat = gsub("p-","",corrMethod)
      
    }else if(corrMethod=="discor"){
      cor.method = "discor"
      cor.stat = "discor"
    }
    
    if(nrow(micdat)>2000){
      if(micDataType=="ko"){
        keepft = rownames(current.proc$mic$res_deAnal)[1:2000]
      }else{
        keepft =phyloseq_objs$res_deAnal[[taxalvl]] [1:2000]
      }
      micdat <-  micdat[rownames( micdat) %in%keepft, ]
      
    }
    
    res.corr <- DoM2Mcorr(micdat,metdat,cor.method,cor.stat,taxalvl)
    
    if(is.null(res.corr$corr.mat)| length(res.corr$corr.mat)==1){
      corr.mat <- 0
      corr.pval <- 0
      current.msg<<-"No correlation is detected using the selected parameters! Please adjust the parameters!"
      return(0)
    }else{
      
      if(!is.matrix(res.corr$corr.mat)){
        corr.mat <- as.matrix(res.corr$corr.mat)
      }else{
        corr.mat <- res.corr$corr.mat
      }
      colnames(corr.mat) <- qvec
      corr.mat <- reshape2::melt(corr.mat,value.name = "correlation")
      
      if(!is.null(res.corr$corr.pval) & length(res.corr$corr.pval)>1){
        if(!is.matrix(res.corr$corr.pval)){
          corr.pval <- as.matrix(res.corr$corr.pval)
        }else{
          corr.pval <- res.corr$corr.pval
        }
        corr.pval <-  reshape2::melt(corr.pval,value.name = "pval")
        corr.mat$pval <- corr.pval$pval
        #  corr.mat <- corr.mat[which(corr.mat$pval < corrPval),]
      }
      fast.write(corr.mat, file=paste(imgNm, ".csv", sep=""), row.names=F);
      current.proc$keggmap$corrplot <<- corr.mat
      if(corrSign=="positive"){
        corr.mat <- corr.mat[which(corr.mat$correlation>0),]
      }else if(corrSign=="negative"){
        corr.mat <- corr.mat[which(corr.mat$correlation< 0),]  
      }
      
      if(nrow(corr.mat)>0 & length(unique(corr.mat$Var2))==1){
        corr.mat <- corr.mat[1:min(topNum,nrow(corr.mat)),]
        names(corr.mat)[1:2] <- c("mic","met")
        colnm = 1
        wb=6
        hb = max(0.25*nrow(corr.mat),2.5)
        wc=7
        hc = 8
      }else if(nrow(corr.mat)>0 & length(unique(corr.mat$Var2))>1){
        library(dplyr)
        names(corr.mat)[1:2] <- c("mic","met")
        corr.mat <-  corr.mat[order(corr.mat$correlation,decreasing = T),] %>% 
          group_by(met) %>%
          mutate(idx=1:length(mic)) %>% 
          filter(idx<(min(topNum,nrow(corr.mat))+1))
        colnm = 2
        wb= 8
        hb = min(0.25*nrow(corr.mat)/2,60)
        wc= 12
        hc = 3*length(qvec)
      }else{
        current.msg<<-"No correlation is detected using the selected parameters! Please adjust the parameters!"
        return(0)
      }
      ylim0 = min(min(corr.mat$correlation)-0.1,-0.1)
      ylim1 = max(corr.mat$correlation)+0.1
      
      require("Cairo");
      library(ggplot2);
      library(viridis);
      library(geomtextpath)
      
      barplot <- vector("list",length=length(qvec))
      circleplot <- vector("list",length=length(qvec))
      for(pl in 1:length(qvec)){
        plot.df <- corr.mat %>% filter(met==qvec[pl]) 
        plot.df <-   plot.df[order(abs(plot.df$correlation)),]
        plot.df$mic <- factor(plot.df$mic,levels = unique(plot.df$mic))
        plot.df$met <- factor(plot.df$met,levels = unique(plot.df$met))
        
        if("pval" %in% colnames(plot.df)){
          
          barplot[[pl]] <-  ggplot(plot.df, aes(x=mic, y=correlation,fill= pval)) + 
            scale_fill_viridis_c(option = "plasma",alpha = 0.8)+
            geom_bar(stat = "identity") + 
            ggtitle(unname(current.proc$id2nm[qvec[pl]]))+
            xlab("")+theme_minimal()+coord_flip() 
          if(length(qvec)>1){
            barplot[[pl]] <-  barplot[[pl]] +theme(legend.key.size = unit(0.45, 'cm'))
          }
          if(nrow(plot.df)>2 ){
            angle <-  90 - 360 * (1:nrow(plot.df)-0.5) /nrow(plot.df)
            plot.df$hjust<-ifelse( angle < -90, 1, 0)
            plot.df$angle<-ifelse(angle < -90, angle+180, angle) 
            plot.df$yh <- 0.05
            if(length(which(plot.df$correlation>0))>0){
              pidx <- which(plot.df$correlation>0)
              lidx <- which(plot.df$correlation>0 & nchar(as.character(plot.df$mic))>20)
              plot.df$yh[pidx] <-  max(plot.df$correlation)-nchar(as.character(plot.df$mic[pidx]))/100
              sidx <- which((plot.df$yh[pidx]-(plot.df$correlation[pidx]+0.05))>0)
              plot.df$yh[pidx[sidx]] <- plot.df$correlation[pidx[sidx]]+0.05
              plot.df$yh[lidx] <- max(plot.df$correlation)-0.2-nchar(as.character(plot.df$mic[lidx]))/100
            }
            if(pl%%2==0){
              yl=""
            }else{
              yl="correlation"
            }
            circleplot[[pl]] <-  ggplot(plot.df,aes(x=mic, y=correlation,fill= pval))+
              geom_bar(stat="identity", color="black")+
              ylim(ylim0,ylim1) + xlab("")+ylab(yl)+
              theme_minimal() + ggtitle(unname(current.proc$id2nm[qvec[pl]]))+
              geom_text(data=plot.df, aes(x=mic, y=yh, label=mic, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= plot.df$angle, inherit.aes = FALSE )+ 
              coord_polar(start = 0) +scale_fill_viridis_c(option = "plasma",alpha = 0.7)+
              theme(
                axis.text.x = element_blank(),
                # axis.text.y = element_text(hjust = -15),
                # plot.margin = margin(1, 1, 1, 1, "cm")
              )  
            if(length(qvec)>1){
              circleplot[[pl]] <-  circleplot[[pl]] + theme(legend.key.size = unit(0.4, 'cm'))
            }
          }else{
            current.msg<<-"Circle plot is not supported when associated taxa are less than 3!" 
          }
          
          
        }else{
          
          barplot[[pl]] <-  ggplot(plot.df, aes(x=mic, y=correlation,fill= correlation)) + 
            scale_fill_viridis_c(option = "plasma",alpha = 0.8)+
            geom_bar(stat = "identity") + xlab("")+ theme_minimal()+
            coord_flip() + ggtitle(unname(current.proc$id2nm[qvec[pl]]))
          
          if(nrow(plot.df)>2){
            angle <-  90 - 360 * (1:nrow(plot.df)-0.5) /nrow(plot.df)
            plot.df$hjust<-ifelse( angle < -90, 1, 0)
            plot.df$angle<-ifelse(angle < -90, angle+180, angle)
            plot.df$yh <- 0.05
            if(length(which(plot.df$correlation>0))>0){
              pidx <- which(plot.df$correlation>0)
              lidx <- which(plot.df$correlation>0 & nchar(as.character(plot.df$mic))>20)
              plot.df$yh[pidx] <-  max(plot.df$correlation)-nchar(as.character(plot.df$mic[pidx]))/100
              sidx <- which((plot.df$yh[pidx]-(plot.df$correlation[pidx]+0.05))>0)
              plot.df$yh[pidx[sidx]] <- plot.df$correlation[pidx[sidx]]+0.05
              plot.df$yh[lidx] <- max(plot.df$correlation)-0.2-nchar(as.character(plot.df$mic[lidx]))/100
            }
            circleplot[[pl]] <-  ggplot(plot.df,aes(x=mic, y=correlation,fill= correlation))+
              geom_bar(stat="identity", color="black")+
              ylim(ylim0,ylim1) + xlab("")+
              theme_minimal() +
              geom_text(data=plot.df, aes(x=mic, y=yh, label=mic, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= plot.df$angle, inherit.aes = FALSE )+ 
              coord_polar(start = 0) +scale_fill_viridis_c(option = "plasma",alpha = 0.7)+
              theme(
                axis.text.x = element_blank()
              ) + ggtitle(unname(current.proc$id2nm[qvec[pl]]))
            
          }else{
            current.msg<<-"Circle plot is not supported when associated taxa are less than 3!"
            
          }
          
        }
      }
      library(grid)
      library(gridExtra)
      library(gridGraphics);
      library(cowplot)
      imgNm.bar <- paste("barplot_",imgNm,  ".png",sep="");
      imgNm.circle <- paste("circleplot_",imgNm,  ".png",sep="");
      Cairo(file=imgNm.bar, width=wb, height=hb, type="png", bg="white", unit="in", dpi=100);
      grid.arrange(grobs =barplot, ncol=colnm)
      dev.off();
      if(exists("circleplot")){
        Cairo(file=imgNm.circle, width=wc, height=hc, type="png", bg="white", unit="in", dpi=100);
        grid.arrange(grobs =circleplot, ncol=colnm)
        dev.off();
      }
      #print("accociation plot done")
      
      # write json
      mic = corr.mat$mic; if(length(mic) ==1) { mic <- matrix(mic) };
      met = corr.mat$met; if(length(met) ==1) { met <- matrix(met) };
      correlation = corr.mat$correlation; if(length(correlation) ==1) { correlation <- matrix(correlation) };
      pval = corr.mat$pval; if(length(pval) ==1) { pval <- matrix(pval) };
      fdr <- p.adjust(corr.mat$pval, "fdr"); if(length(fdr) ==1) { fdr <- matrix(fdr) };
      
      json.res <- list(mic = mic,
                       met = met,
                       correlation = correlation,
                       pval = pval,
                       fdr = fdr);
      
      json.mat <- RJSONIO::toJSON(json.res);
      json.nm <- paste(imgNm, ".json", sep="");
      sink(json.nm)
      cat(json.mat);
      sink();
      
    }
    
  }else{
    
    if(micDataType=="otu" & !is.null(keggid)){
      if(!(exists("gem",current.proc))){
        met.map <- MetaboIDmap(netModel="gem", predDB="agora", IDtype=metIDType);
      }
      if(!exists("gem.mic.map.qs")){
        mic.map <- MicIDmap(predModel="gem",predDB="agora",taxalvl="")
      }
      
      M2Mprediction(model="gem", predDB="agora",taxalvl="all",metType=metIDType)
      current.proc$meta_para$analysis.var<<-selected.meta.data
      PerformFeatDEAnalyse(NA,taxalvl=taxalvl, 0, analysisVar=selected.meta.data, adjustedVar=NULL);
      
      PerformPairDEAnalyse(NA,taxalvl=taxalvl, overlay="T", initDE="0", analysisVar=selected.meta.data, adjustedVar=NULL);
      
      predDE <- current.proc$predDE
      qvec =  current.proc$keggNet$Query[match(keggid,current.proc$keggNet$Match)] 
      
      data.plot = predDE[which(predDE$met==qvec),c(4,5,2)]
    }else{
      
      
      PerformFeatDEAnalyse(NA,taxalvl="OTU", 0,  analysisVar=selected.meta.data, adjustedVar=NULL);
      performeCorrelation(NA,taxalvl="OTU",initDE=0)
      
    }
    
    
    if(nrow(data.plot)==0){
      return(0)
      current.msg<<-"no mtached metabolite was found"
    }
    if(nrow(data.plot)>20){
      data.plot = data.plot[1:20,]
    }
    
    data.plot$`-log(p)` = -log(data.plot$p_value)
    data.plot <- data.plot[order(data.plot$p_value,decreasing = T),]
    data.plot$mic <- factor(data.plot$mic,levels = data.plot$mic)
    
    require("Cairo");
    library(ggplot2);
    library(viridis);
    imgNm <- paste(imgNm,  ".png",sep="");
    Cairo(file=imgNm, width=5, height=5, type="png", bg="white", unit="in", dpi=100);
    p1 <- ggplot(data.plot, aes(x=mic, y=`-log(p)`,fill=`-log(p)`)) + 
      scale_fill_viridis_c(option = "plasma",alpha = 0.8)+
      geom_bar(stat = "identity") + xlab("")+
      coord_flip()
    print(p1)
    dev.off();
    
  }
  
  
  print("plot done")
  return(1)
}



UpdateAssociationPlot <- function(imgNm,topNum=10){
  current.msg<<-"null"
  topNum = as.numeric(topNum)
  corr.mat <-current.proc$keggmap$corrplot
  corrSign<-current.proc$keggmap$corrSign
  qvec <- current.proc$keggmap$current_qvec
  library(dplyr)
  if(corrSign=="positive"){
    corr.mat <- corr.mat[which(corr.mat$correlation>0),]
  }else if(corrSign=="negative"){
    corr.mat <- corr.mat[which(corr.mat$correlation< 0),]
  }
  
  if(is.null(corr.mat) | nrow(corr.mat)==0){
    current.msg<<-"Correlation failed! Please check the parameters"
    return(0)
  }
  names(corr.mat)[1:2] <- c("mic","met")
  if(nrow(corr.mat)>0 & length(unique(corr.mat$met))==1){
    corr.mat <- corr.mat[1:min(topNum,nrow(corr.mat)),]
    colnm = 1
    wb=6
    hb = max(0.25*nrow(corr.mat),2.5)
    wc=6
    hc = 7
  }else if(nrow(corr.mat)>0 & length(unique(corr.mat$met))>1){
    
    corr.mat <-  corr.mat[order(corr.mat$correlation,decreasing = T),] %>% 
      group_by(met) %>%
      mutate(idx=1:length(mic)) %>% 
      filter(idx<(min(topNum,nrow(corr.mat))+1))
    colnm = 2
    wb= 8
    hb = min(0.25*nrow(corr.mat)/2,60)
    wc= 12
    hc = 3*length(qvec)
  }else{
    current.msg<<-"No correlation is detected using the selected parameters! Please adjust the parameters!"
    return(0)
  }
  ylim0 = min(min(corr.mat$correlation)-0.1,-0.1)
  ylim1 = max(corr.mat$correlation)+0.1
  require("Cairo");
  library(ggplot2);
  library(viridis);
  library(geomtextpath)
  
  barplot <- vector("list",length=length(qvec))
  circleplot <- vector("list",length=length(qvec))
  for(pl in 1:length(qvec)){
    plot.df <- corr.mat %>% filter(met==qvec[pl]) 
    plot.df <-   plot.df[order(abs(plot.df$correlation)),]
    plot.df$mic <- factor(plot.df$mic,levels = unique(plot.df$mic))
    plot.df$met <- factor(plot.df$met,levels = unique(plot.df$met))
    
    if("pval" %in% colnames(plot.df)){
      
      barplot[[pl]] <-  ggplot(plot.df, aes(x=mic, y=correlation,fill= pval)) + 
        scale_fill_viridis_c(option = "plasma",alpha = 0.8)+
        geom_bar(stat = "identity") + 
        ggtitle(qvec[pl])+
        xlab("")+theme_minimal()+coord_flip() 
      if(length(qvec)>1){
        barplot[[pl]] <-  barplot[[pl]] +theme(legend.key.size = unit(0.45, 'cm'))
      }
      if(nrow(plot.df)>2 ){
        angle <-  90 - 360 * (1:nrow(plot.df)-0.5) /nrow(plot.df)
        plot.df$hjust<-ifelse( angle < -90, 1, 0)
        plot.df$angle<-ifelse(angle < -90, angle+180, angle) 
        plot.df$yh <- 0.05
        if(length(which(plot.df$correlation>0))>0){
          pidx <- which(plot.df$correlation>0)
          lidx <- which(plot.df$correlation>0 & nchar(as.character(plot.df$mic))>20)
          plot.df$yh[pidx] <-  max(plot.df$correlation)-nchar(as.character(plot.df$mic[pidx]))/100
          sidx <- which((plot.df$yh[pidx]-(plot.df$correlation[pidx]+0.05))>0)
          plot.df$yh[pidx[sidx]] <- plot.df$correlation[pidx[sidx]]+0.05
          plot.df$yh[lidx] <- max(plot.df$correlation)-0.2-nchar(as.character(plot.df$mic[lidx]))/100
        }
        if(pl%%2==0){
          yl=""
        }else{
          yl="correlation"
        }
        circleplot[[pl]] <-  ggplot(plot.df,aes(x=mic, y=correlation,fill= pval))+
          geom_bar(stat="identity", color="black")+
          ylim(ylim0,ylim1) + xlab("")+ylab(yl)+
          theme_minimal() +  ggtitle(qvec[pl])+
          geom_text(data=plot.df, aes(x=mic, y=yh, label=mic, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= plot.df$angle, inherit.aes = FALSE )+ 
          coord_polar(start = 0) +scale_fill_viridis_c(option = "plasma",alpha = 0.7)+
          theme(
            axis.text.x = element_blank(),
            # axis.text.y = element_text(hjust = -15),
            # plot.margin = margin(1, 1, 1, 1, "cm")
          )  
        if(length(qvec)>1){
          circleplot[[pl]] <-  circleplot[[pl]] + theme(legend.key.size = unit(0.4, 'cm'))
        }
      }else{
        current.msg<<-"Circle plot is not supported when associated taxa are less than 3!" 
      }
      
      
    }else{
      
      barplot[[pl]] <-  ggplot(plot.df, aes(x=mic, y=correlation,fill= correlation)) + 
        scale_fill_viridis_c(option = "plasma",alpha = 0.8)+
        geom_bar(stat = "identity") + xlab("")+ theme_minimal()+
        coord_flip() + ggtitle(qvec[pl])
      
      if(nrow(plot.df)>2){
        angle <-  90 - 360 * (1:nrow(plot.df)-0.5) /nrow(plot.df)
        plot.df$hjust<-ifelse( angle < -90, 1, 0)
        plot.df$angle<-ifelse(angle < -90, angle+180, angle)
        plot.df$yh <- 0.05
        if(length(which(plot.df$correlation>0))>0){
          pidx <- which(plot.df$correlation>0)
          lidx <- which(plot.df$correlation>0 & nchar(as.character(plot.df$mic))>20)
          plot.df$yh[pidx] <-  max(plot.df$correlation)-nchar(as.character(plot.df$mic[pidx]))/100
          sidx <- which((plot.df$yh[pidx]-(plot.df$correlation[pidx]+0.05))>0)
          plot.df$yh[pidx[sidx]] <- plot.df$correlation[pidx[sidx]]+0.05
          plot.df$yh[lidx] <- max(plot.df$correlation)-0.2-nchar(as.character(plot.df$mic[lidx]))/100
        }
        circleplot[[pl]] <-  ggplot(plot.df,aes(x=mic, y=correlation,fill= correlation))+
          geom_bar(stat="identity", color="black")+
          ylim(ylim0,ylim1) + xlab("")+
          theme_minimal() +
          geom_text(data=plot.df, aes(x=mic, y=yh, label=mic, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= plot.df$angle, inherit.aes = FALSE )+ 
          coord_polar(start = 0) +scale_fill_viridis_c(option = "plasma",alpha = 0.7)+
          theme(
            axis.text.x = element_blank()
          ) + ggtitle(qvec[pl])
        
      }else{
        current.msg<<-"Circle plot is not supported when associated taxa are less than 3!"
        
      }
      
    }
  }
  
  library(grid)
  library(gridExtra)
  library(gridGraphics);
  library(cowplot)
  imgNm.bar <- paste("barplot_",imgNm,  ".png",sep="");
  imgNm.circle <- paste("circleplot_",imgNm,  ".png",sep="");
  Cairo(file=imgNm.bar, width=wb, height=hb, type="png", bg="white", unit="in", dpi=100);
  grid.arrange(grobs =barplot, ncol=colnm)
  dev.off();
  if(exists("circleplot")){
    Cairo(file=imgNm.circle, width=wc, height=hc, type="png", bg="white", unit="in", dpi=100);
    grid.arrange(grobs =circleplot, ncol=colnm)
    dev.off();
  }
  print("Update accociation plot done")
  
  # write json
  mic = corr.mat$mic; if(length(mic) ==1) { mic <- matrix(mic) };
  met = corr.mat$met; if(length(met) ==1) { met <- matrix(met) };
  correlation = corr.mat$correlation; if(length(correlation) ==1) { correlation <- matrix(correlation) };
  pval = corr.mat$pval; if(length(pval) ==1) { pval <- matrix(pval) };
  fdr <- p.adjust(corr.mat$pval, "fdr"); if(length(fdr) ==1) { fdr <- matrix(fdr) };
  
  json.res <- list(mic = mic,
                   met = met,
                   correlation = correlation,
                   pval = pval,
                   fdr = fdr);
  
  json.mat <- RJSONIO::toJSON(json.res);
  json.nm <- paste(imgNm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();
  
}

tuneKOmap <- function(){
  edges.ko = qs::qread(paste0(lib.path.mmp,"ko.info.qs"))
  include = rownames(current.proc$mic$data.proc)
  edges.ko = edges.ko[which(edges.ko$ko %in% include),]
  includeInfo = list(edges=edges.ko)
  includeInfo$nodes = unique(edges.ko$met)
  
  json.mat <- rjson::toJSON(includeInfo);
  sink("includeInfo.json");
  cat(json.mat);
  sink();
  
}

###########################################################
####################3DScatter Plot#########################
###########################################################
###########################################################


DoDimensionReductionIntegrative <- function(mbSetObj, reductionOpt, method="globalscore", dimn,analysisVar,diabloPar=0.2){
  if(!exists("my.reduce.dimension")){ # public web on same user dir
    .load.scripts.on.demand("utils_dimreduction.Rc");    
  }
  if(analysisVar=="null"){
    analysisVar = current.proc$meta_para$analysis.var
  }
  return(my.reduce.dimension(mbSetObj, reductionOpt, method,dimn, analysisVar,diabloPar));
}

doScatterJsonPair <- function(filenm,analysisVar,taxrank="genus"){
   
  if(!exists("my.json.scatter.pair")){ # public web on same user dir
    .load.scripts.on.demand("utils_scatter_json.Rc");    
  }
  
  return(my.json.scatter.pair(filenm,current.proc$meta_para$analysis.var, taxrank));
}

DoStatComparisonVis <- function(filenm, alg, meta, selected, meta.vec, omicstype, taxalvl,nonpar=FALSE){
  
  mbSetObj <- .get.mbSetObj(NA);
  
  if(taxalvl=="null"|taxalvl=="NULL"|is.null(taxalvl)){
    taxalvl="OTU"
  }
  
  if(meta == "null"){
    meta = 1;
  }
  combined.res <- qs::qread("combined.res.qs")
  if(meta.vec == "NA"){ # process page
    #if(dataSet$name != filenm){
    dataSet <- readRDS(filenm);
    #}
  }else{
    
    if(omicstype != "NA"){
      
      if(omicstype %in% c("microbiome","mic")){
        data <-  qs::qread("phyloseq_objs.qs") 
        data<- data$count_tables[[taxalvl]]
        
      }else{
        data<- current.proc$met$data.proc
        
      }
      
      
    }else{
      data.mic <- current.proc$mic$data.proc
      data.met <- current.proc$met$data.proc
      
    }
  }
  
  if(meta.vec == "NA"){ # process page
    if(meta == ""){
      meta=1;
    }
    metavec = current.proc$meta_para$sample_data
    sel = unique(metavec)
  }else{
    metavec <- strsplit(meta.vec, "; ")[[1]];
    sel <- strsplit(selected, "; ")[[1]];
  }
  
  combined.res$meta$newcolumn = metavec
  metadf = combined.res$meta
  metadf = metadf[which(metadf$omics==omicstype),]
  
  
  sel_meta1 = metadf[which(metadf[,"newcolumn"] %in% sel[1]),]
  sel_meta2 = metadf[which(metadf[,"newcolumn"] %in% sel[2]),] 
  
  nms1 <- rownames(sel_meta1)
  nms2 <- rownames(sel_meta2)
  sel_meta_more_than_2 = metadf[which(metadf[,"newcolumn"] %in% sel),]
  nms <- rownames(sel_meta_more_than_2)
  sample_type <- mbSetObj$dataSet$meta_info
  
  newcol_disc <- c(T);
  names(newcol_disc) <- "newcolumn";
  newcol_cont <- c(F);
  names(newcol_cont) <- "newcolumn";
  sample_type$disc_inx <- c(sample_type$disc.inx, newcol_disc);
  sample_type$cont_inx <- c(sample_type$cont.inx, newcol_cont);
  
  if(alg=="ttest"){
    res <- PerformFastUnivTests(data,factor(metadf[,"newcolumn"]),F,F)
  }else if(alg =="kruskal"){
    res <- PerformFastUnivTests(data,factor(metadf[,"newcolumn"]),F,T)
  }else if(alg =="limma"){
    metadf[,meta] <- metadf[,"newcolumn"]
    res <- performLimma(data,metadf,sample_type,meta)
  }
  res = res[,c(1,2)]
  rownames(res) = rownames(data)
  colnames(res) = c("stat", "p_value")
  
  res = na.omit(res)
  res = res[order(res[,2], decreasing=FALSE),]
  pvals <- p.adjust(res[,"p_value"],method="BH");
  res = cbind(res, pvals)
  res = cbind(res, rownames(res))
  de = res
  de[de == "NaN"] = 1
  pv = as.numeric(de[,"p_value"])
  pv_no_zero = pv[pv != 0]
  minval = min(pv_no_zero)
  pv[pv == 0] = minval/2
  pvals <- -log10(pv);
  colorb<- ComputeColorGradient(pvals, "black");
  sizes <- as.numeric(rescale2NewRange(-log10(pv), 15, 35));
  res = cbind(res, colorb);
  res = cbind(res, sizes);
  
  #ids <- names(dataSet$enrich_ids[order(match(combined.res$enrich_ids,rownames(res)))])
  res = cbind(res, rownames(res));
  colnames(res) = c("stat", "p_value", "p_adj", "ids", "color", "size","name");
  res= as.matrix(res)
  library(RJSONIO)
  sink(filenm);
  cat(toJSON(res));
  sink();
  
  if(meta.vec == "NA"){
    filenm = "OK";
    combined.res[[omicstype]]$comp.res = de;
    combined.res$sel.meta = meta
  }
  
  return(filenm)
}


.calNMI<-function (x, y){
  x <- as.vector(x)
  y <- as.vector(y)
  mutual.info <- (.mutualInformation(x, y)/sqrt(.entropy(x) * 
                                                  .entropy(y)))
  return(max(0, mutual.info, na.rm = TRUE))
}


###########################################################
####################list input################################
###########################################################
###########################################################


PrepareListInput<-function(mbSetObj, qvec, omic){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$inputType <- "list";
  mbSetObj$micDataType <- "otu"
  lines <- unlist(strsplit(qvec, "\r|\n|\r\n")[1]);
  if(substring(lines[1],1,1)=="#"){
    lines <- lines[-1];
  }
  mbSetObj$dataSet[[omic]][["original"]] <- lines;
  inputType <<- "list"
  return(.set.mbSetObj(mbSetObj))
}

SetListInfo <-function(mbSetObj,taxalvl){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$paraSet<-list()
  mbSetObj$paraSet$taxalvl <-taxalvl
  mbSetObj$paraSet$metDataType <-metDataType
  mbSetObj$paraSet$metIDType <-metIDType
  taxalvl<<-taxalvl
  return(.set.mbSetObj(mbSetObj))
}

PerformMicNameMap <- function(mbSetObj,taxalvl){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mic.map = data.frame(Query=mbSetObj$dataSet$mic$original,agora=NA,embl=NA,kegg=NA,ncbi=NA,stringsAsFactors = F)
  taxMapLong <- qs::qread(paste0(lib.path.mmp,"agora_tax.qs"))[[taxalvl]]
  names(taxMapLong)[1] <- "taxa"
  
  mic.map$agora <- taxMapLong[match(mic.map$Query,taxMapLong$taxa),1]
  
  mic.map$ncbi <- taxMapLong$ncbi_id[match(mic.map$agora,taxMapLong$taxa)]
  
  taxMapLong <- qs::qread(paste0(lib.path.mmp,"embl_tax.qs"))[[taxalvl]]
  names(taxMapLong)[1] <- "taxa"
  mic.map$embl <- taxMapLong[match(mic.map$Query,taxMapLong$taxa),1]
  
  mic.map$ncbi[is.na(mic.map$ncbi)] <- taxMapLong$ncbi_id[match(mic.map$embl[is.na(mic.map$ncbi)],taxMapLong$taxa)]
  
  taxMapKEGG <- qs::qread(paste0(lib.path.mmp,"taxMapKEGG.qs"))[[taxalvl]]
  taxMapLong <- taxMapKEGG[["info"]]
  names(taxMapLong)[1] <- "taxa"
  mic.map$kegg <- taxMapLong[match(mic.map$Query,taxMapLong$taxa),1]
  
  mic.map$ncbi[is.na(mic.map$ncbi)] <- taxMapLong$id[match(mic.map$kegg[is.na(mic.map$ncbi)],taxMapLong$taxa)]
  
  mtchidx <-  taxMapKEGG[which(names(taxMapKEGG) %in% mic.map$Query)]
  mtcls <<- unique(unlist(mtchidx))
  
  fast.write(mic.map, paste("taxa_match_result.csv"));
  
  mbSetObj$analSet$mic.map = mic.map
  mbSetObj$analSet$mtcls = mtcls
  
  return(.set.mbSetObj(mbSetObj))
}


PerformMetNameMap <- function(mbSetObj,metIDType="kegg"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  met.map = data.frame(Query=mbSetObj$dataSet$met$original,agora=NA,embl=NA,kegg=NA,stringsAsFactors = F)
  
  res = MetaboIDmap("gem","agora",metIDType,met.map$Query)
  met.map$agora_id = res$Match[match( met.map$Query,res$Query)]
  res = MetaboIDmap("gem","embl",metIDType,met.map$Query)
  met.map$embl_id = res$Match[match( met.map$Query,res$Query)]
  
  if(metIDType !='name'){
    metInfo <- qs::qread(paste0(lib.path.mmp,"synonymGem.qs"));
    met.map$agora = metInfo$Name[match(met.map$agora_id,metInfo$metID)]
    met.map$embl = metInfo$Name[match(met.map$embl_id,metInfo$metID)]
  }
  
  if(metIDType=="kegg"){
    met.map$kegg = met.map$Query
    metInfo <- qs::qread(paste0(lib.path.mmp,"general_kegg2name.qs"));
    
    met.map$name <- metInfo$Name[match(met.map$kegg,metInfo$ID)]
    met.map$node <- metInfo$node[match(met.map$kegg,metInfo$ID)]
  }else{
    res = MetaboIDmap("keggNet","kegg",metIDType,met.map$Query)
    met.map$kegg = res$Match[match( met.map$Query,res$Query)]
    met.map$name = res$Name
    met.map$node = res$Node
  }
  
  fast.write(met.map, paste("metabolite_match_result.csv"));
  
  mbSetObj$analSet$met.map = met.map
  
  return(.set.mbSetObj(mbSetObj))
}

GetMicMapCol <-function(mbSetObj, colInx){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(colInx==0){
    return(rownames(mbSetObj$analSet$mic.map));
  }else{
    return(mbSetObj$analSet$mic.map[,colInx]);
  }
  
}

GetMetMapCol <-function(mbSetObj, colInx){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(colInx==0){
    return(rownames(mbSetObj$analSet$met.map));
  }else{
    return(mbSetObj$analSet$met.map[,colInx]);
  }
  
}

PerformMetListEnrichment <- function(mbSetObj, contain,file.nm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(contain=="bac"){
    current.set <- qs::qread(paste0(lib.path.mmp,"kegg_bac_mummichog.qs"))$pathways$cpds
    
  }else if(contain=="hsabac"){
    current.set <- qs::qread(paste0(lib.path.mmp,"kegg_hsa_bac_mummichog.qs"))$pathways$cpds
    
  }else if(contain=="hsa"){
    current.set <- qs::qread(paste0(lib.path.mmp,"kegg_hsa_mummichog.qs"))$pathways$cpds
    
  }else if(contain=="all"){
    current.set <- qs::qread(paste0(lib.path.mmp,"kegg_all_mummichog.qs"))$pathways$cpds
    
  }else{
    current.set <- qs::qread(paste0(taxalvl,".current.lib.qs"))
    
  }
  
  current.universe <- unique(unlist(current.set));
  met.map <- mbSetObj$analSet$met.map
  
  # prepare for the result table
  set.size <- length(current.set);
  res.mat <- matrix(0, nrow=set.size, ncol=5);
  rownames(res.mat) <- names(current.set);
  colnames(res.mat) <- c("Total", "Expected", "Hits", "Pval", "FDR");
  
  # prepare query
  ora.vec <- NULL;
  ora.vec <- met.map$kegg[!is.na(met.map$kegg)];
  
  # need to cut to the universe covered by the pathways, not all genes
  hits.inx <- ora.vec %in% current.universe;
  ora.vec <- ora.vec[hits.inx];
  #ora.nms <- ora.nms[hits.inx];
  q.size <- length(ora.vec);
  
  # get the matched query for each pathway
  hits.query <- lapply(current.set, function(x) {
    ora.vec[ora.vec%in%unlist(x)];});
  names(hits.query) <- names(current.set);
  
  hit.num <- unlist(lapply(hits.query, function(x){length(x)}), use.names=FALSE);
  
  # total unique gene number
  uniq.count <- length(current.universe);
  
  # unique gene count in each pathway
  set.size <- unlist(lapply(current.set, length));
  
  res.mat[,1] <- set.size;
  res.mat[,2] <- q.size*(set.size/uniq.count);
  res.mat[,3] <- hit.num;
  
  # use lower.tail = F for P(X>x)
  raw.pvals <- phyper(hit.num-1, set.size, uniq.count-set.size, q.size, lower.tail=F);
  res.mat[,4] <- raw.pvals;
  res.mat[,5] <- p.adjust(raw.pvals, "fdr");
  
  # now, clean up result, synchronize with hit.query
  res.mat <- res.mat[hit.num>0,,drop = F];
  hits.query <- hits.query[hit.num>0];
  
  if(nrow(res.mat)> 1){
    # order by p value
    ord.inx <- order(res.mat[,4]);
    res.mat <- signif(res.mat[ord.inx,],3);
    hits.query <- hits.query[ord.inx];
    imp.inx <- res.mat[,4] <= 0.05;
    
    if(sum(imp.inx) < 10){ # too little left, give the top ones
      topn <- ifelse(nrow(res.mat) > 10, 10, nrow(res.mat));
      res.mat <- res.mat[1:topn,];
      hits.query <- hits.query[1:topn];
    }else{
      res.mat <- res.mat[imp.inx,];
      hits.query <- hits.query[imp.inx];
      
      if(sum(imp.inx) > 120){
        # now, clean up result, synchronize with hit.query
        res.mat <- res.mat[1:120,];
        hits.query <- hits.query[1:120];
      }
    }
  }
  fast.write(res.mat, file=paste(file.nm, ".csv", sep=""), row.names=F);
  
  resTable <- data.frame(Pathway=rownames(res.mat), res.mat,check.names=FALSE);
  
  path.pval = resTable$Pval; if(length(path.pval) ==1) { path.pval <- matrix(path.pval) };
  hit.num = paste0(resTable$Hits,"/",resTable$Total); if(length(hit.num) ==1) { hit.num <- matrix(hit.num) };
  path.nms <- resTable$Pathway; if(length(path.nms) ==1) { path.nms <- matrix(path.nms) };
  path.fdr <- resTable$FDR; if(length(path.fdr) ==1) { path.fdr <- matrix(path.fdr) };
  expr.mat = setNames(rep(2,q.size),ora.vec)
  hits.met <- lapply(hits.query, function(x){
    x=met.map$name[match(x,met.map$kegg)]  
    return(x)
  })
  hits.node <- lapply(hits.query, function(x){
    x=met.map$node[match(x,met.map$kegg)]  
    return(x)
  })
  
  mbSetObj <- recordEnrTable(mbSetObj, "mmp", resTable, "KEGG", "Overrepresentation Analysis");

  json.res <- list(hits.query =convert2JsonList(hits.query),
                   path.nms = path.nms,
                   path.pval = path.pval,
                   hit.num = hit.num,
                   path.fdr = path.fdr,
                   hits.query.nm = convert2JsonList(hits.met),
                   hits.node = convert2JsonList(hits.node),
                   expr.mat=expr.mat);
  
  json.mat <-  RJSONIO::toJSON(json.res);
  json.nm <- paste(file.nm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();

  .set.mbSetObj(mbSetObj)
  return(mbSetObj);
}



M2MPredictionList<- function(mbSetObj,model,predDB,psc=0.5,metType="metabolite",taxalvl){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if(predDB=="null"| is.null(predDB) | predDB==""){
    predDB <- "agora"
  }
  mbSetObj$paraSet$gemdb<-predDB
  
  require(reshape2)
  message('Loading the model database..')
  psc <- as.numeric(psc)
  taxalvl<-tolower(taxalvl) 
  
  tax_map <- mbSetObj$analSet$mic.map
  tax_map <- tax_map[which(!is.na(tax_map[,predDB])),]
  m2m_ls <- qs::qread(paste0(lib.path.mmp,predDB,".qs"))[[taxalvl]]
  names(m2m_ls)[1] <- "taxa"
  m2m_ls <- m2m_ls[which(m2m_ls$potential>=psc),]
  m2m_ls <- m2m_ls[which(m2m_ls$taxa %in% tax_map[,predDB]),]
  m2m_ls$taxa <- tax_map$Query[match(m2m_ls$taxa,tax_map[,predDB])]
  
  if(metType=="metabolite"){
    met.map<-  mbSetObj$analSet$met.map
    m2m_ls <- m2m_ls[which(m2m_ls$metID %in% met.map[,paste0(predDB,"_id")]),]
    if(mbSetObj$paraSet$metIDType=="name"){
      m2m_ls$metabolite <- met.map$Query[match(m2m_ls$metID,met.map[,paste0(predDB,"_id")])]
    }
  }
  
  
  mbSetObj$analSet$predres<-m2m_ls
  qs::qsave(m2m_ls,paste0("m2m_pred_",predDB,".qs"))
  mbSetObj$analSet$integration$db <- predDB;
  message("Prediction completed!")
  return(.set.mbSetObj(mbSetObj));
}




CreatM2MHeatmapList<-function(mbSetObj, plotNm,  format="png", 
                              smplDist="euclidean", clstDist="ward.D", palette="npj",viewOpt="barraw", 
                              clustRow="T", clustCol="T", 
                              colname="T",rowname="T", fontsize_col=10, fontsize_row=10,
                              potential.thresh=0.5,
                              var.inx=NA, border=T, width=NA, dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  load_iheatmapr();
  load_rcolorbrewer();
  load_viridis();
  current.msg<<- NULL
  set.seed(2805614);
  #used for color pallete
  ######set up plot
  #colors for heatmap
  
  if(palette=="gbr"){
    colors <- grDevices::colorRampPalette(c("green", "black", "red"), space="rgb")(256);
  }else if(palette == "heat"){
    colors <- grDevices::heat.colors(256);
  }else if(palette == "topo"){
    colors <- grDevices::topo.colors(256);
  }else if(palette == "gray"){
    colors <- grDevices::colorRampPalette(c("grey90", "grey10"), space="rgb")(256);
  }else if(palette == "byr"){
    colors <- rev(grDevices::colorRampPalette(RColorBrewer::brewer.pal(10, "RdYlBu"))(256));
  }else if(palette == "viridis") {
    colors <- rev(viridis::viridis(10))
  }else if(palette == "plasma") {
    colors <- rev(viridis::plasma(10))
  }else  if(palette == "npj"){
    colors <- c("#00A087FF","white","#E64B35FF")
  }else  if(palette == "aaas"){
    colors <- c("#4DBBD5FF","white","#E64B35FF");
  }else  if(palette == "d3"){
    colors <- c("#2CA02CFF","white","#FF7F0EFF");
  }else {
    colors <- colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")), alpha=0.8)(100)
    #c("#0571b0","#92c5de","white","#f4a582","#ca0020");
  }
  
  plotjs = paste0(plotNm, ".json");
  plotNm = paste(plotNm, ".", format, sep="");
  mbSetObj$imgSet$IntegrationHeatmap<-plotNm;
  ####using  prediction pair pval
  pred.dat <-mbSetObj$analSet$predres
  data.abd <- data.frame(mic=pred.dat$taxa,
                         met=pred.dat$metabolite,
                         var = pred.dat$potential)
  
  data <- data.abd[order(data.abd$var,decreasing = T),]
  
  if(nrow(data)>1000){
    thresh <- data$var[1000]
  }else{
    thresh<- potential.thresh
  }
  
  data <- reshape2::dcast(data,mic~met)
  data[is.na(data)] <-0
  data.mtr <- data[,-1]
  
  micnms <- data$mic
  metnms <- colnames(data.mtr)
  nameHt <- "Potential Score"
  
  #### fro annotation using pval
  
  data1 <- data.mtr;
  data1sc <- as.matrix(apply(data1, 2, as.numeric))
  rownames(data1sc) <- micnms
  #data1sc <- scale_mat(data1sc, scaleOpt)
  
  fzCol <- round(as.numeric(fontsize_col), 1)
  fzRow <- round(as.numeric(fontsize_row), 1)
  map.height=nrow(data1)*30
  map.width=ncol(data1)*30
  
  #cb_grid <- setup_colorbar_grid(nrow = 100, x_start = 1.1, y_start = 0.95, x_spacing = 0.15)
  dend_row <- hclust(dist(data1sc, method = smplDist), method = clstDist)
  p <- iheatmap(data1sc,
                # colorbar_grid = cb_grid, 
                name = nameHt, x_categorical = TRUE,
                layout = list(font = list(size = 10)),
                colors = colors
  )
  
  if (clustRow == "true") {
    p <- p %>% add_row_dendro(dend_row, side = "right")
  }
  
  if (colname == "true" ){
    
    p <- p %>%  add_col_labels(size = 0.1, font = list(size = fzCol))
  }
  
  if (colname == "true" ){
    
    p <- p %>% add_row_labels(size = 0.1, font = list(size = fzRow), side = "left")
  }
  
  
  if (clustCol == "true") {
    dend_col <- hclust(dist(t(data1), method = smplDist), method = clstDist)
    p <- p %>% add_col_dendro(dend_col)
  }
  
  
  as_list <- to_plotly_list(p)

  plotwidget <- paste0(plotNm, ".rda")
  pwidget <- to_widget(p)
  save(pwidget, file = plotwidget)
  
  if (viewOpt != "overview") {
    as_list[["layout"]][["width"]] <- max(map.width,1000)
    as_list[["layout"]][["height"]] <- max(map.height,800)
  } else {
    as_list[["layout"]][["width"]] <- 1200
    as_list[["layout"]][["height"]] <- map.height
  }
  
  
  
  as_json <- attr(as_list, "TOJSON_FUNC")(as_list)
  as_json <- paste0("{ \"x\":", as_json, ",\"evals\": [],\"jsHooks\": []}")
  
  write(as_json, plotjs)
  
  if(is.null(current.msg)){
    current.msg<<-"null"
  }
  # storing for Report Generation
  mbSetObj$analSet$integration$heatmap <- data1sc
  mbSetObj$analSet$integration$heatmap.dist <- smplDist
  mbSetObj$analSet$integration$heatmap.clust <- clstDist
  mbSetObj$analSet$integration$taxalvl <- taxalvl
  mbSetObj$analSet$integration$overlay <- "false"
  mbSetObj$analSet$integration$htMode <- "prediction"
  mbSetObj$analSet$integration$potential <- potential.thresh;
  mbSetObj$imgSet$heatmap_predmmp <- plotwidget;
  message("heatmap done")
  return(.set.mbSetObj(mbSetObj))
}


GetPredictionPlot <- function(mbSetObj, keggid,imgNm,predDB="agora",potentialThresh=0.5,topNum=10,subset="false"){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  current.msg<<-"null"
  topNum = as.numeric(topNum)
  potentialThresh = as.numeric(potentialThresh)
  mic.map = mbSetObj$analSet$mic.map
  met.map = mbSetObj$analSet$met.map
  macthed_cmpd = met.map$kegg[which( !is.na(met.map[,predDB]))]
  metType = mbSetObj$paraSet$metDataType
  
  if(!keggid %in% met.map$kegg){
    current.msg <<- "The selected compound is not provided in your input data! Please choose the related compounds!"
    return(0)
  }else if(!keggid %in% macthed_cmpd){
    current.msg <<- "The selected compound is not supported by the selected GEM database!"
    return(0)
  }
  
  if(metType=="metabolite"){
    
    qvec =  met.map$Query[which(met.map$kegg==keggid)] 
    
  }else{
    qvecidx =  which(current.proc$keggNet$Match==keggid)
    adduct = current.proc$keggNet$adduct[qvecidx]
    if(length(interaction(adduct,primary_ions))>0){
      qvec = current.proc$keggNet$Query[which(current.proc$keggNet$Match==keggid)][which(adduct %in% primary_ions)]
    }else{
      qvec =  current.proc$keggNet$Query[which(current.proc$keggNet$Match==keggid)] 
    }
    print(paste0(length(qvec)," peaks related! The plots take longer time."))
    
  }
  
  mbSetObj$analSet$current_qvec <- qvec
  message('Loading the model database..')
  taxalvl<-tolower(taxalvl) 
  
  tax_map <- mic.map[which(!is.na(mic.map[,predDB])),]
  m2m_ls <- qs::qread(paste0(lib.path.mmp,predDB,".qs"))[[taxalvl]]
  names(m2m_ls)[1] <- "taxa"
  m2m_ls <- m2m_ls[which(m2m_ls$potential>=potentialThresh),]
  m2m_ls <- m2m_ls[which(m2m_ls$taxa %in% tax_map[,predDB]),]
  m2m_ls$taxa <- tax_map$Query[match(m2m_ls$taxa,tax_map[,predDB])]
  
  if(metType=="metabolite"){
    m2m_ls <- m2m_ls[which(m2m_ls$metID %in% met.map[which(met.map$kegg==keggid),paste0(predDB,"_id")]),]
    m2m_ls$met <- met.map$Query[match(m2m_ls$metID,met.map[,paste0(predDB,"_id")])]
    predres <- data.frame(mic = m2m_ls$taxa,met=m2m_ls$met,potentialSore=m2m_ls$potential,stringsAsFactors = F)
  }
  library(dplyr)
  if(nrow(predres)>0 & length(unique(predres$met))==1){
    predres <-  predres[order(predres$potentialSore,decreasing = T),]
    predres <- predres[1:min(topNum,nrow(predres)),]
    colnm = 1
    wb=6
    hb = max(0.25*nrow(predres),2.5)
    wc=7
    hc = 7
  }else if(nrow(predres)>0 & length(unique(predres$met))>1){
    predres <-  predres[order(predres$potentialSore,decreasing = T),] %>% 
      group_by(met) %>%
      mutate(idx=1:length(mic)) %>% 
      filter(idx<(min(topNum,nrow(corr.mat))+1))
    colnm = 2
    wb= 8
    hb = min(0.25*nrow(predres)/2,60)
    wc= 12
    hc = 3*length(qvec)
  }else{
    current.msg<<-"No prediction is detected using the selected parameters! Please adjust the parameters!"
    return(0)
  }
  ylim0 = min(min(predres$potentialSore)-0.1,-0.1)
  ylim1 = max(predres$potentialSore)+0.1
  
  
  
  require("Cairo");
  library(ggplot2);
  library(viridis);
  library(geomtextpath)
  
  barplot <- vector("list",length=length(qvec))
  circleplot <- vector("list",length=length(qvec))
  for(pl in 1:length(qvec)){
    plot.df <- predres %>% filter(met==qvec[pl]) 
    plot.df <-   plot.df[order(abs(plot.df$potentialSore)),]
    plot.df$mic <- factor(plot.df$mic,levels = unique(plot.df$mic))
    plot.df$met <- factor(plot.df$met,levels = unique(plot.df$met))
    
    barplot[[pl]] <-   ggplot(plot.df, aes(x=mic, y=potentialSore,fill= potentialSore)) + 
      scale_fill_viridis_c(option = "plasma",alpha = 0.8)+
      geom_bar(stat = "identity") + xlab("")+ ylab("Potential Sore")+ theme_minimal()+
      coord_flip() + ggtitle(qvec[pl])+labs(fill =  "Potential Sore")
    
    if(nrow(plot.df)>2){
      angle <-  90 - 360 * (1:nrow(plot.df)-0.5) /nrow(plot.df)
      plot.df$hjust<-ifelse( angle < -90, 1, 0)
      plot.df$angle<-ifelse(angle < -90, angle+180, angle)
      plot.df$yh <- 0.05
      lidx <- which(nchar(as.character(plot.df$mic))>15)
      plot.df$yh <-  max(plot.df$potentialSore)-nchar(as.character(plot.df$mic))/90
      sidx <- which((plot.df$yh-(plot.df$potentialSore+0.05))>0)
      plot.df$yh[sidx] <- plot.df$potentialSore[sidx]+0.05
      plot.df$yh[lidx] <- max(plot.df$potentialSore)-0.2-nchar(as.character(plot.df$mic[lidx]))/100 
      
      circleplot[[pl]] <-  ggplot(plot.df,aes(x=mic, y=potentialSore,fill= potentialSore))+
        geom_bar(stat="identity", color="black")+
        ylim(ylim0,ylim1) + xlab("")+ ylab("Potential Sore")+
        theme_minimal() +
        geom_text(data=plot.df, aes(x=mic, y=yh, label=mic, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= plot.df$angle, inherit.aes = FALSE )+ 
        coord_polar(start = 0) +scale_fill_viridis_c(option = "plasma",alpha = 0.7)+labs(fill =  "Potential Sore")+
        theme(
          axis.text.x = element_blank()
        ) + ggtitle(qvec[pl])
      
    }else{
      current.msg<<-"Circle plot is not supported when associated taxa are less than 3!"
      
    }
    
  }
  
  library(grid)
  library(gridExtra)
  library(gridGraphics);
  library(cowplot)
  imgNm.bar <- paste("barplot_",imgNm,  ".png",sep="");
  imgNm.circle <- paste("circleplot_",imgNm,  ".png",sep="");
  Cairo(file=imgNm.bar, width=wb, height=hb, type="png", bg="white", unit="in", dpi=100);
  grid.arrange(grobs =barplot, ncol=colnm)
  dev.off();
  if(exists("circleplot")){
    Cairo(file=imgNm.circle, width=wc, height=hc, type="png", bg="white", unit="in", dpi=100);
    grid.arrange(grobs =circleplot, ncol=colnm)
    dev.off();
  }
  print("prediction plot done")
  
  mbSetObj$analSet$current_predres <- predres
  # write json
  mic = predres$mic; if(length(mic) ==1) { mic <- matrix(mic) };
  met = predres$met; if(length(met) ==1) { met <- matrix(met) };
  potentialSore = predres$potentialSore; if(length(potentialSore) ==1) { potentialSore <- matrix(potentialSore) };
  
  json.res <- list(mic = mic,
                   met = met,
                   potentialSore = potentialSore);
  
  json.mat <- RJSONIO::toJSON(json.res);
  json.nm <- paste(imgNm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();
  return(.set.mbSetObj(mbSetObj))
}



UpdatePredictionPlot <- function(mbSetObj,imgNm,topNum=10){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  current.msg<<-"null"
  topNum = as.numeric(topNum)
  
  qvec <-  mbSetObj$analSet$current_qvec 
  predres <- mbSetObj$analSet$current_predres
  
  library(dplyr)
  if(nrow(predres)>0 & length(unique(predres$met))==1){
    predres <-  predres[order(predres$potentialSore,decreasing = T),]
    predres <- predres[1:min(topNum,nrow(predres)),]
    colnm = 1
    wb=6
    hb = max(0.25*nrow(predres),2.5)
    wc=7
    hc = 7
  }else if(nrow(predres)>0 & length(unique(predres$met))>1){
    predres <-  predres[order(predres$potentialSore,decreasing = T),] %>% 
      group_by(met) %>%
      mutate(idx=1:length(mic)) %>% 
      filter(idx<(min(topNum,nrow(corr.mat))+1))
    colnm = 2
    wb= 8
    hb = min(0.25*nrow(predres)/2,60)
    wc= 12
    hc = 3*length(qvec)
  }else{
    current.msg<<-"No prediction is detected using the selected parameters! Please adjust the parameters!"
    return(0)
  }
  ylim0 = min(min(predres$potentialSore)-0.1,-0.1)
  ylim1 = max(predres$potentialSore)+0.1
  
  
  
  require("Cairo");
  library(ggplot2);
  library(viridis);
  library(geomtextpath)
  
  barplot <- vector("list",length=length(qvec))
  circleplot <- vector("list",length=length(qvec))
  for(pl in 1:length(qvec)){
    plot.df <- predres %>% filter(met==qvec[pl]) 
    plot.df <-   plot.df[order(abs(plot.df$potentialSore)),]
    plot.df$mic <- factor(plot.df$mic,levels = unique(plot.df$mic))
    plot.df$met <- factor(plot.df$met,levels = unique(plot.df$met))
    
    barplot[[pl]] <-   ggplot(plot.df, aes(x=mic, y=potentialSore,fill= potentialSore)) + 
      scale_fill_viridis_c(option = "plasma",alpha = 0.8)+
      geom_bar(stat = "identity") + xlab("")+ ylab("Potential Sore")+ theme_minimal()+
      coord_flip() + ggtitle(qvec[pl])+labs(fill =  "Potential Sore")
    
    if(nrow(plot.df)>2){
      angle <-  90 - 360 * (1:nrow(plot.df)-0.5) /nrow(plot.df)
      plot.df$hjust<-ifelse( angle < -90, 1, 0)
      plot.df$angle<-ifelse(angle < -90, angle+180, angle)
      plot.df$yh <- 0.05
      lidx <- which(nchar(as.character(plot.df$mic))>15)
      plot.df$yh <-  max(plot.df$potentialSore)-nchar(as.character(plot.df$mic))/90
      sidx <- which((plot.df$yh-(plot.df$potentialSore+0.05))>0)
      plot.df$yh[sidx] <- plot.df$potentialSore[sidx]+0.05
      plot.df$yh[lidx] <- max(plot.df$potentialSore)-0.2-nchar(as.character(plot.df$mic[lidx]))/100 
      
      circleplot[[pl]] <-  ggplot(plot.df,aes(x=mic, y=potentialSore,fill= potentialSore))+
        geom_bar(stat="identity", color="black")+
        ylim(ylim0,ylim1) + xlab("")+ ylab("Potential Sore")+
        theme_minimal() +
        geom_text(data=plot.df, aes(x=mic, y=yh, label=mic, hjust=hjust), color="black", fontface="bold",alpha=0.8, size=3, angle= plot.df$angle, inherit.aes = FALSE )+ 
        coord_polar(start = 0) +scale_fill_viridis_c(option = "plasma",alpha = 0.7)+labs(fill =  "Potential Sore")+
        theme(
          axis.text.x = element_blank()
        ) + ggtitle(qvec[pl])
      
    }else{
      current.msg<<-"Circle plot is not supported when associated taxa are less than 3!"
      
    }
    
  }
  
  library(grid)
  library(gridExtra)
  library(gridGraphics);
  library(cowplot)
  imgNm.bar <- paste("barplot_",imgNm,  ".png",sep="");
  imgNm.circle <- paste("circleplot_",imgNm,  ".png",sep="");
  Cairo(file=imgNm.bar, width=wb, height=hb, type="png", bg="white", unit="in", dpi=100);
  grid.arrange(grobs =barplot, ncol=colnm)
  dev.off();
  if(exists("circleplot")){
    Cairo(file=imgNm.circle, width=wc, height=hc, type="png", bg="white", unit="in", dpi=100);
    grid.arrange(grobs =circleplot, ncol=colnm)
    dev.off();
  }
  print("prediction plot done")
  
  mbSetObj$analSet$current_predres <- predres
  # write json
  mic = predres$mic; if(length(mic) ==1) { mic <- matrix(mic) };
  met = predres$met; if(length(met) ==1) { met <- matrix(met) };
  potentialSore = predres$potentialSore; if(length(potentialSore) ==1) { potentialSore <- matrix(potentialSore) };
  
  json.res <- list(mic = mic,
                   met = met,
                   potentialSore = potentialSore);
  
  json.mat <- RJSONIO::toJSON(json.res);
  json.nm <- paste(imgNm, ".json", sep="");
  sink(json.nm)
  cat(json.mat);
  sink();
  
  return(.set.mbSetObj(mbSetObj))
}
###########################################################
####################Plots################################
###########################################################
###########################################################

PlotCorrHistogram <- function(imgNm, dpi=72, format="png"){
  dpi<-as.numeric(dpi)
  imgNm <- paste(imgNm, "dpi", dpi, ".", format, sep="");
  
  library(ggplot2)
  reductionSet <- .get.rdt.set();
  
  fig.list <- list();
  for( i in 1:2){
    if(i == 1){
      cors <- reductionSet$corr.mat.inter
      titleText <- "Between-omics correlation"
    }else{
      cors <- reductionSet$corr.mat.intra
      titleText <- "Intra-omics correlation"
    }
    cor.data <- cors[upper.tri(cors, diag = FALSE)]  # we're only interested in one of the off-diagonals, otherwise there'd be duplicates
    cor.data <- as.data.frame(cor.data)  # that's how ggplot likes it
    summary(cor.data)
    colnames(cor.data) <- "coefficient"
    coefficient.p <- function(r, n) {
      pofr <- ((1-r^2)^((n-4)/2))/beta(a = 1/2, b = (n-2)/2)
      return(pofr)
    }
    
    fig.list[[i]] <- ggplot(cors, aes(x=correlation)) + geom_histogram() +
      xlim(-1,1) +
      theme_bw()
  }
  library(Cairo)
  library(ggpubr)
  Cairo(file=imgNm, width=10, height=8, unit="in", type="png", bg="white", dpi=dpi);
  p1 <- ggarrange(plotlist=fig.list, ncol = 1, labels=c("Between-omics correlation", "Intra-omics correlation"))
  print(p1)
  dev.off();
  
}


PlotDiagnostic <- function(imgName, dpi=72, format="png",alg, taxrank="OTU"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  dpi <- as.numeric(dpi);
  imgNm <- paste(imgName,  ".", format, sep="");
  require("Cairo");
  if(alg %in% c("snf", "spectrum") ){
    h=8
    fig.list <- list()
    library(ggpubr);
  }else{
    h=8
  }
  
  Cairo(file=imgNm, width=10, height=h, type=format,unit="in", bg="white", dpi=dpi);
  if(alg == "procrustes"){
    procrustes.res <- qs::qread("procrustes.res.qs")
    if(length(procrustes.res$dim.res) == 1){
    res <- procrustes.res$dim.res[[1]]
    }else{
    res <- procrustes.res$dim.res[[taxrank]]
    }
    error = residuals(res[[1]])
    require("ggrepel")
    
    error.df = data.frame(Samples=names(error), Procrustes_residual=unname(error))
    
    rankres <- rank(-abs(error), ties.method="random");
    inx.x <- which(rankres < 6);
    inx.y <- error[inx.x];
    nms <- names(error)[inx.x];
    subsetdf <- error.df[which(error.df$Samples %in% nms),]
    
    p = ggplot(error.df, aes(x = Samples, y = Procrustes_residual)) + geom_col() + geom_hline(yintercept = summary(error)[c(2,3,5)])+
      theme(axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            text = element_text(size = 16)) +
      geom_text_repel(
        data = subsetdf,
        aes(label = Samples),
        size = 5,
        box.padding = unit(0.35, "lines"),
        point.padding = unit(0.3, "lines")
      ) +
      theme_bw()
    print(p)
    mbSetObj$imgSet$procrustes$diagnostic <- imgNm
  }else if(alg == "diablo"){
    require(mixOmics)
    diablo.res <- qs::qread("diablo.res.qs")
    res <- diablo.res$dim.res[[length(diablo.res$dim.res)]]
    set.seed(123) # for reproducibility, only when the `cpus' argument is not used
    # this code takes a couple of min to run
    perf.res <- mixOmics:::perf(res, validation = 'Mfold', folds = 10, nrepeat = 1, dist="max.dist",near.zero.var=T)
    diablo.comp <<- median(perf.res$choice.ncomp$WeightedVote)
    plot(perf.res) 
    mbSetObj$imgSet$diablo$diagnostic <- imgNm
  }
  dev.off();
  .set.mbSetObj(mbSetObj)
  return(1);
}


PlotDiagnosticPca <- function(imgNm, dpi=72, format="png",type="diablo", taxrank="OTU"){
  #save.image("diag.RData");
  mbSetObj <- .get.mbSetObj(mbSetObj);
  require("Cairo");
  library(ggplot2);
  dpi<-as.numeric(dpi)
  imgNm<- paste(imgNm, ".", format, sep="");
  fig.list <- list()
  if(type == "diablo"){ 
    library(grid)
    library(gridExtra)
    library(gridGraphics);
    library(cowplot)
    diablo.res <- qs::qread("diablo.res.qs")
    if(length(diablo.res$dim.res) == 1){
        dim.res <- diablo.res$dim.res[[length(diablo.res$dim.res)]]
    }else{
        dim.res <- diablo.res$dim.res[[taxrank]]
    }
    fig.list[[1]] <- as_grob(function(){
      plotDiablo(dim.res, ncomp = 1)
    })
    
    fig.list[[2]] <- as_grob(function(){
      plotDiablo(dim.res, ncomp = 2)
    })
    
    fig.list[[3]] <- as_grob(function(){
      plotDiablo(dim.res, ncomp = 3)
    })
    h<-8*round(length(fig.list))
    
    Cairo(file=imgNm, width=10, height=h, type=format, bg="white", unit="in", dpi=dpi);
    
    grid.arrange(grobs =fig.list, nrow=length(fig.list))
    dev.off(); 
    mbSetObj$imgSet$diablo$pca <- imgNm;
  }else if(type == "procrustes"){
    library(ggplot2)
    library(grid)
    procrustes.res <- qs::qread("procrustes.res.qs")
    if(length(procrustes.res$dim.res) == 1){
    pro.test <- procrustes.res$dim.res[[1]][[1]]
    }else{
    pro.test <- procrustes.res$dim.res[[taxrank]][[1]]
    }
    pct <- pro.test$svd$d
    ctest <- data.frame(rda1=pro.test$Yrot[,1], rda2=pro.test$Yrot[,2], xrda1=pro.test$X[,1],
                        xrda2=pro.test$X[,2],Type=procrustes.res$newmeta[,"omics"], Conditions = procrustes.res$newmeta[,1])
    xlabel <- paste0("Component 1 ", "(" , signif(pct[1],4), ")")
    ylabel <- paste0("Component 2 ", "(" , signif(pct[2],4), ")")
    
    p <- ggplot(ctest) +
      geom_point(aes(x=rda1, y=rda2, colour=Conditions, shape=Type)) +
      geom_point(aes(x=xrda1, y=xrda2, colour=Conditions, shape=Type)) +
      geom_segment(aes(x=rda1,y=rda2,xend=xrda1,yend=xrda2,colour=Conditions), alpha=0.4,arrow=arrow(length=unit(0.1,"cm"))) + 
      xlab(xlabel) + ylab(ylabel) +
      theme_bw()
    Cairo(file=imgNm, width=10, height=10, type=format, bg="white", unit="in", dpi=dpi);
    print(p)
    dev.off();
    mbSetObj$imgSet$procrustes$pca <- imgNm
  }
  .set.mbSetObj(mbSetObj)
}


PlotDiagnosticLoading <- function(imgNm, dpi=72, format="png",type="diablo",taxrank="OTU"){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  require("Cairo");
  library(ggplot2)
  dpi <- as.numeric(dpi);
  imgNm <- paste(imgNm,  ".", format, sep="");
  mbSetObj$imgSet$diablo$loading <- imgNm
  if(type == "diablo"){
    library(grid)
    library(gridExtra)
    library(cowplot)
    fig.list <- list()
    diablo.res <- qs::qread("diablo.res.qs")
    if(length(diablo.res$dim.res) == 1){
    dim.res <- diablo.res$dim.res[[1]]
    }else{
    dim.res <- diablo.res$dim.res[[taxrank]]
    }
    fig.list[[1]] <- as_grob(function(){
      plotLoadings(dim.res, ndisplay=10, comp = 1, contrib="max", method="median", size.name=1.1, legend=T)
    })
    fig.list[[2]] <- as_grob(function(){
      plotLoadings(dim.res, ndisplay=10, comp = 2, contrib="max", method="median", size.name=1.1, legend=T)
    })
    fig.list[[3]] <-as_grob(function(){
      plotLoadings(dim.res, ndisplay=10, comp = 3, contrib="max", method="median", size.name=1.1, legend=T)
    })
    h <- 8*round(length(fig.list))
    Cairo(file=imgNm, width=13, height=h, type=format, bg="white", unit="in", dpi=dpi);
    grid.arrange(grobs =fig.list, nrow=length(fig.list))
    dev.off();
  }
  .set.mbSetObj(mbSetObj)
}

GetDiagnosticSummary<- function(type){
  if(type %in% c("perturbation", "spectrum", "snf")){
    reductionSet <- .get.rdt.set();
    clustNum <- length(unique(reductionSet$clustVec))
    return(c(clustNum, signif(reductionSet$clustNmi)))
  }else if(type == "procrustes"){
    procrustes.res <- qs::qread("procrustes.res.qs")
    res <-procrustes.res$dim.res[[length(procrustes.res$dim.res)]][[2]];
    return(c(signif(res$ss,4), signif(res$scale,4)));
  }else{
    return(c("","") )
  }
}


gg_color_hue <- function(grp.num, type="green", filenm=NULL) {
  
  grp.num <- as.numeric(grp.num)
  if(type == "green"){
    pal18 <- c("#e6194B", "#3cb44b", "#4363d8", "#ffff00", "#f032e6", "#ffe119", "#911eb4", "#f58231", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#000075");
  }else{
    pal18 <- c( "#4363d8","#e6194B" , "#3cb44b", "#f032e6", "#ffe119", "#e6194B", "#f58231", "#bfef45", "#fabebe", "#469990", "#e6beff", "#9A6324", "#800000", "#aaffc3", "#808000", "#ffd8b1", "#42d4f4","#000075", "#ff4500");
  }
  if(grp.num <= 18){ # update color and respect default
    colArr <- pal18[1:grp.num];
  }else{
    colArr <- colorRampPalette(pal18)(grp.num);
  }
  if(is.null(filenm)){
    return(colArr);
  }else{
    library(RJSONIO)
    sink(filenm);
    cat(toJSON(colArr));
    sink();
    return(filenm);
  }
}


###########################################################
####################Helpper functions######################
###########################################################
###########################################################

SetMMPDataType <- function(inputType,micDataType,metDataType,metIDType){
  inputType<<-inputType
  micDataType<<-micDataType
  metDataType<<-metDataType
  metIDType<<-metIDType
  
  return(1)
}

InitCurrentProc <-function(){
  current.proc<- vector("list",length=2)
  names(current.proc)<-c("mic","met")
  moduleType <<-"mmp"
  if(metIDType=="kegg"){
    metInfo <- qs::qread(paste0(lib.path.mmp,"general_kegg2name.qs"));
    mbSetObj <- .get.mbSetObj(mbSetObj);
    keggids <- rownames(mbSetObj$dataSet$metabolomics$data.orig)
    nms<- metInfo$Name[match(keggids, metInfo$ID)]
    current.proc$id2nm <- setNames(nms,keggids)
  }else{
    mbSetObj <- .get.mbSetObj(mbSetObj);
    nms <- rownames(mbSetObj$dataSet$metabolomics$data.orig)
    current.proc$id2nm <- setNames(nms,nms)
  }
  
  mbSetObj$inputType <- inputType;
  mbSetObj$micDataType <-micDataType;
  mbSetObj$metDataType <-metDataType;
  mbSetObj$metIDType <-metIDType;
  
  current.proc<<-current.proc
  .set.mbSetObj(mbSetObj)
  return(1)
  
}

SetPeakParameter<-function(rtOpt,mode,instrumentOpt){
  if(rtOpt=="no"){
    current.proc$mumRT <- FALSE
    current.proc$mumRT.type <- NA;
  }else{
    current.proc$mumRT <- TRUE
    current.proc$mumRT.type <- rtOpt;
  }
  current.proc$mode <- mode;
  current.proc$adducts <- NA;
  current.proc$peakFormat <- "mpt";
  current.proc$instrument <- as.numeric(instrumentOpt);
  current.proc$rt_frac <- 0.02
  current.proc$primary_ion <- "yes"
  current.proc <<- current.proc
  return(1);
}


RemoveData <- function(dataName){
  mbSetObj <- .get.mbSetObj(NA);
  if(!is.null(mbSetObj$dataSets[[dataName]])){
    mbSetObj$dataSets[[dataName]] <- NULL;
    unlink(paste0(dataName, "_data"), recursive = TRUE)
  }
  if(mbSetObj$module.type=="meta"){
    if(!is.null(mdata.all[[dataName]])){
      mdata.all[[dataName]] <<- NULL;
    }
  }
  
  return(.set.mbSetObj(mbSetObj));
}

GetMicMetDataDims <- function(dataType,dataName){
  if(is.null(current.proc$mic)){
    
    data<-data.table::fread(dataName)
    dm <- dim(data);
    dm[2] <- dm[2]
    naNum <- sum(is.na(data));
    
    
  }else{
    
    if(dataType=="mic"){
      dm <- current.proc$mic$data.proc
      naNum <- sum(is.na(current.proc$mic$data.proc));
      
    }else{
      dm <- current.proc$met$data.proc
      naNum <- sum(is.na(current.proc$met$data.proc));
      
    }
  }
  
  return(c(dm, naNum));
}


GetMetaTaxaInfoMMP <- function(mbSetObj,istaxalbl){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  proc.phyobj <- mbSetObj$dataSet$proc.phyobj;
  
  #check that each rank has >2 groups
  taxa.tbl <- as(tax_table(proc.phyobj), "matrix")
  
  if(ncol(taxa.tbl)==1){
    taxa.nms <- "Phylum"
    return(taxa.nms)
  }
  
  #drop taxa with only 1 level (i.e. Viruses at Phylum)
  gd.inx <- apply(taxa.tbl, 2, function(x) length(unique(x))!=1);
  taxa.tbl.update <- taxa.tbl[,gd.inx, drop=FALSE];
  
  if(ncol(taxa.tbl.update) == 0){
    current.msg <<- c("All taxa info for the remaining features are the same!")
    return("OTU")
  }
  
  taxa.nms <- rank_names(taxa.tbl.update);
  
  return(c(taxa.nms[!is.na(taxa.nms)],"OTU"));
  
}


CleanMMP<- function(){
  rm(list = ls())
}


ReScale <- function(x,first,last){(last-first)/(max(x)-min(x))*(x-min(x))+first}

rescale2NewRange <- function(qvec, a, b){ReScale
  qvec = replace(qvec, qvec == 0, 1)
  q.min <- min(qvec);
  q.max <- max(qvec);
  if(length(qvec) < 50){
    a <- a*2;
  }
  if(q.max == q.min){
    new.vec <- rep(8, length(qvec));
  }else{
    coef.a <- (b-a)/(q.max-q.min);
    const.b <- b - coef.a*q.max;
    new.vec <- coef.a*qvec + const.b;
  }
  return(new.vec);
}


GetNetsName <- function(){
  rownames(net.stats);
}

GetNetsEdgeNum <- function(){
  as.numeric(net.stats$Edge);
}

GetNetsNodeNum <- function(){
  as.character(net.stats$Node);
}

GetNetsQueryNum <- function(){
  as.numeric(net.stats$Query);
}

GetMicNm <- function(){
  as.character(corrRes$source);
}

GetMetNm <- function(){
  as.character(corrRes$target);
}

GetCorrIdx <- function(){
  as.numeric(corrRes$corrindex);
}

GetCorrPval <- function(){
  as.numeric(corrRes$corrpval);
}

GetPredPval <- function(){
  as.numeric(corrRes$predpval);
}


GetNetsNameString <- function(){
  nms <- paste(rownames(net.stats), collapse="||");
  return(nms);
}

GetCorrRes <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  corr.pval <- mbSetObj$analSet$corr.pval
  idxkeep <- rownames(corr.pval) %in% mbSetObj$analSet$m2mNet_graphlist$nodes$name
  corr.pval <- reshape2::melt(corr.pval[idxkeep,idxkeep])
  corrRes <- mbSetObj$analSet$m2mNet_graphlist$links
  corrRes <- dplyr::left_join(corrRes,corr.pval,by=c("source"="Var1","target"="Var2"))
  corrRes <- corrRes[,c(1,2,3,7,4)]
  names(corrRes)[3:5] <- c("corrindex","corrpval","predpval")
  corrRes<<-corrRes
  return(1)
}


GetCurrentScatter  <- function(){
  
  return(jsonNms_scatter)
}


convert2JsonList <- function(my.list){
  lapply(my.list, function(x){
    if(length(x) == 1){
      matrix(x);
    }else{
      x;
    }
  });
}


ComputeEncasingDiablo <- function(filenm, type, names.vec, level=0.95, omics="NA",taxalvl="OTU"){
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


ComputeEncasing <- function(filenm, type, names.vec, level=0.95, omics="NA"){
  
  
  level <- as.numeric(level)
  names = strsplit(names.vec, "; ")[[1]]
  
  pos.xyz <- qs::qread("pos.xyz.qs");
  #print(head(pos.xyz));
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
