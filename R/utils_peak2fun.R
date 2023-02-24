performPeakEnrich <- function(lib){

  mz.de <- current.proc$met$res_deAnal
  current.proc$pos_inx <- rep(as.numeric(as.logical(current.proc$mode =="positive")) == 1,nrow(mz.de));

  if( current.proc$mumRT){
    
    spl <- ifelse(grepl("__",rownames(mz.de)[1]),"__","@")
    mzls <- strsplit(rownames(mz.de),split = spl)
    current.proc$ref_mz  <- lapply(mzls,function(x) x[1])
    current.proc$ret_time <- unlist(lapply(mzls,function(x) x[2]))
    #names(current.proc$ret_time ) <-  current.proc$ref_mz 
    rt_tol <- max(as.numeric(current.proc$ret_time)) *current.proc$rt_frac 
    current.proc$rt_tol <-  rt_tol
    }else{
    
    current.proc$ref_mz  <- rownames(mz.de)
    
  }
  current.proc$expr_dic <-  current.proc$met$res_deAnal[,1]
  names( current.proc$expr_dic)<- current.proc$ref_mz 
  
   sig.idx <- min(length(current.proc$met$sigfeat),0.1*nrow(mz.de))
  current.proc$sig.feat <- unlist(current.proc$ref_mz[1:round(sig.idx)])
  current.proc<<-current.proc

  .setup.library(lib);
  
  if(current.proc$mumRT){
     .init.RT.Permutations(permNum)
  } else {
     .init.Permutations(permNum)
  }
  
  return(1);
}



.setup.library <- function( lib, minLib=3, adduct.custom=NULL){
  
if(lib=="usrbac"|lib=="sigbac"){

pathways = qs::qread(paste0(taxalvl,".current.lib.qs"))
allmet = qs::qread(paste0(lib.path.mmp,"keggmet_allpathway.qs"))

cpd.lib = allmet[!is.na(allmet$mw) & allmet$id %in% unique(unlist(pathways)),]

cpd.lib <- list(
  id = cpd.lib$id,
  name = cpd.lib$name,
  mw = as.numeric(cpd.lib$mw)
)

create.adducts(cpd.lib)

lib <- qs::qread("cmn_compound_lib.qs")
mummichog.lib<-list()
mummichog.lib$pathways$cpds <- pathways
mummichog.lib$pathways$name <- names(pathways)
mummichog.lib$cpd.tree <- lib$cpd.tree
mummichog.lib$cpd.lib <- lib$cpd.lib
qs::qsave(mummichog.lib,"current_mummichog_lib.qs")

}else{

  if(lib=="bac"){
     filenm <- paste0(lib.path.mmp,"kegg_bac_mummichog.qs")

    }else if(lib=="hsa"){
     filenm <- paste0(lib.path.mmp,"kegg_hsa_mummichog.qs")
    }else if(lib=="all"){
      filenm <- paste0(lib.path.mmp,"kegg_all_mummichog.qs")
    }else if(lib=="hsabac"){
    filenm <- paste0(lib.path.mmp,"kegg_hsa_bac_mummichog.qs")


}
  mummichog.lib <- qs::qread(filenm);

}


   if(!is.null(adduct.custom)){
    mw <- mummichog.lib$cpd.lib$mw;
    new_adducts <- new_adduct_mzlist(current.proc, mw);
    
    cpd.lib <- list(
      mz.matp = new_adducts$pos,
      mz.matn = new_adducts$neg,
      mw = mummichog.lib$cpd.lib$mw,
      id = mummichog.lib$cpd.lib$id,
      name = mummichog.lib$cpd.lib$name
    );
    
  }else{
    cpd.lib <- list(
      mz.matp = mummichog.lib$cpd.lib$adducts[["positive"]],
      mz.matn = mummichog.lib$cpd.lib$adducts[["negative"]],
      mw = mummichog.lib$cpd.lib$mw,
      id = mummichog.lib$cpd.lib$id,
      name = mummichog.lib$cpd.lib$name
    );
  }
  
  cpd.treep <- mummichog.lib$cpd.tree[["positive"]];
  cpd.treen <- mummichog.lib$cpd.tree[["negative"]];
  
  # build empirical compound library after
 # path.length <- sapply(mummichog.lib$pathways$cpds, length)
  
  
  # min.inx <- which(path.length >= minLib)
  # 
  # cleaned.pathways <- vector("list")
  # cleaned.pathways$cpds <- mummichog.lib$pathways$cpds[min.inx]
  # cleaned.pathways$id <- mummichog.lib$pathways$id[min.inx]
  # cleaned.pathways$name <- mummichog.lib$pathways$name[min.inx]
  # 
  current.proc$pathways <<- mummichog.lib$pathways;
  # 
   .search.compoundLib(cpd.lib, cpd.treep, cpd.treen);

  
  if(current.proc$mumRT ){
    # only for empirical compounds
    cleaned.pathways<-current.proc$pathways
      # map cpds to empirical cpds
      cleaned.pathways$emp_cpds <- lapply(cleaned.pathways$cpds, 
                                          function(x) {
                                            unique(unlist(current.proc$cpd_ecpd_dict[na.omit(match(x, names(current.proc$cpd_ecpd_dict)))]))
                                          })
      
      # delete emp_cpds, cpds and names with no emp_cpds
      null.inx <- sapply(cleaned.pathways$emp_cpds, is.null)
      
      new_pathways <- vector("list");
      
      new_pathways$cpds <- cleaned.pathways$cpds[!null.inx]
      new_pathways$name <- cleaned.pathways$name[!null.inx]
      new_pathways$emp_cpds <- cleaned.pathways$emp_cpds[!null.inx]
      
      current.proc$pathways <<- new_pathways;
    
  }
  
  return(1);
}



.search.compoundLib <- function( cpd.lib, 
                                cpd.treep, 
                                cpd.treen){
  
  
  ref_mzlist <- as.numeric(current.proc$ref_mz);
  print(paste0("Got ", length(ref_mzlist), " mass features."))
  pos_inx <- current.proc$pos_inx;
  ref_mzlistp <- ref_mzlist[pos_inx];
  ref_mzlistn <- ref_mzlist[!pos_inx];
  
  # for empirical compounds
  if( current.proc$mumRT ){
    ord_rt <- rank( current.proc$ret_time, ties.method = "random")
    ret_time_pos <- current.proc$ret_time[pos_inx];
    ret_time_rank_pos <- ord_rt[pos_inx];
    ret_time_neg <- current.proc$ret_time[!pos_inx];
    ret_time_rank_neg <- ord_rt[!pos_inx];
    rt_tol <- current.proc$rt_tol;
    rt_tol_rank <- length(ref_mzlist)*current.proc$rt_frac;
  } else {
    # add fake RT
    ret_time_pos <- rep(1, length(ref_mzlistp))
    ret_time_rank_pos <- rep(1, length(ref_mzlistp))
    ret_time_neg <- rep(1, length(ref_mzlistn))
    ret_time_rank_neg <- rep(1, length(ref_mzlistn))
  }
  
  modified.statesp <- colnames(cpd.lib$mz.matp);
  modified.statesn <- colnames(cpd.lib$mz.matn);
  my.tolsp <- mz_tolerance(ref_mzlistp, current.proc$instrument);
  my.tolsn <- mz_tolerance(ref_mzlistn, current.proc$instrument);
  
  # get mz ladder (pos index)
  self.mzsp <- floor(ref_mzlistp);
  all.mzsp <- cbind(self.mzsp-1, self.mzsp, self.mzsp+1);
  
  self.mzsn <- floor(ref_mzlistn);
  all.mzsn <- cbind(self.mzsn-1, self.mzsn, self.mzsn+1);
  
  # matched_res will contain detailed result (cmpd.id. query.mass, mass.diff) for all mz;
  # use a high-performance variant of list
  matched_resp <- myFastList();
  matched_resn <- myFastList();
  
  if(current.proc$mode != "negative"){
    for(i in 1:length(ref_mzlistp)){
      mz <- ref_mzlistp[i];
      rt <- ret_time_pos[i];
      rt_rank <- ret_time_rank_pos[i];
      my.tol <- my.tolsp[i];
      all.mz <- all.mzsp[i,];
      pos.all <- as.numeric(unique(unlist(cpd.treep[all.mz])));
      
      for(pos in pos.all){
        id <- cpd.lib$id[pos];
        mw.all <- cpd.lib$mz.matp[pos,]; #get modified mzs
        diffs <- abs(mw.all - mz); #modified mzs - mz original
        hit.inx <- which(diffs < my.tol);
        if(length(hit.inx)>0){
          for(spot in 1:length(hit.inx)){
            hit.pos <- hit.inx[spot];# need to match all
            index <- paste(mz, id, rt, hit.pos, sep = "___");
            matched_resp$add(index, c(i, id, mz, rt, rt_rank, mw.all[hit.pos], modified.statesp[hit.pos], diffs[hit.pos])); #replaces previous when hit.inx>1
          }
        }
      }
    }
  }
  
  all.mzsn <<- all.mzsn
  
  if (current.proc$mode != "positive") {
    for(i in 1:length(ref_mzlistn)){
      mz <- ref_mzlistn[i];
      rt <- ret_time_neg[i];
      rt_rank <- ret_time_rank_neg[i];
      my.tol <- my.tolsn[i];
      all.mz <- all.mzsn[i,];
      pos.all <- as.numeric(unique(unlist(cpd.treen[all.mz])));
      
      for(pos in pos.all){
        id <- cpd.lib$id[pos]; # position of compound in cpd.tree
        mw.all <- cpd.lib$mz.matn[pos,]; #get modified mzs
        diffs <- abs(mw.all - mz); #modified mzs - mz original
        hit.inx <- which(diffs < my.tol);
        if(length(hit.inx)>0){
          for(spot in 1:length(hit.inx)){
            hit.pos <- hit.inx[spot];# need to match all
            index <- paste(mz, id, rt, hit.pos, sep = "___"); #name in fast_list
            matched_resn$add(index, c(i, id, mz, rt, rt_rank, mw.all[hit.pos], modified.statesn[hit.pos], diffs[hit.pos])); #replaces previous when hit.inx>1
          }
        }
      }
    }
  }
  
  # convert to regular list
  if (current.proc$mode == "mixed") {
    
    matched_resn <- matched_resn$as.list();
    matched_resp <- matched_resp$as.list();
    
    neg_matches <- length(matched_resn) > 0
    pos_matches <- length(matched_resp) > 0
    
    if(!neg_matches & !pos_matches){
      msg.vec <<- "No compound matches from upload peak list!"
      return(0)
    }
    
    if(neg_matches){
      matched_resn <- data.frame(matrix(unlist(matched_resn), nrow=length(matched_resn), byrow=T), stringsAsFactors = FALSE);
    }
    
    if(pos_matches){
      matched_resp <- data.frame(matrix(unlist(matched_resp), nrow=length(matched_resp), byrow=T), stringsAsFactors = FALSE);
    }
    
    if(neg_matches & pos_matches){ # both w. matches
      matched_res <- rbind(matched_resp, matched_resn)
    }else if(neg_matches & !pos_matches){ # only neg w. matches
      matched_res <- matched_resn
    }else{ # only pos w. matches
      matched_res <- matched_resp
    }
    
  } else if(current.proc$mode == "positive") {
    matched_resp <- matched_resp$as.list();
    
    if(is.null(unlist(matched_resp))){
      msg.vec <<- "No compound matches from upload peak list!"
      return(0)
    }
    
    matched_resp <- data.frame(matrix(unlist(matched_resp), nrow=length(matched_resp), byrow=T), stringsAsFactors = FALSE);
    matched_res <- matched_resp;
    
  } else {
    matched_resn <- matched_resn$as.list();
    if(is.null(unlist(matched_resn))){
      msg.vec <<- "No compound matches from upload peak list!"
      return(0)
    }
    
    matched_resn <- data.frame(matrix(unlist(matched_resn), nrow=length(matched_resn), byrow=T), stringsAsFactors = FALSE);
    matched_res <- matched_resn
  }
  
  # re-order columns for output
  matched_res <- matched_res[, c(3,2,7,8,4,5)];
  colnames(matched_res) <- c("Query.Mass", "Matched.Compound", "Matched.Form", "Mass.Diff", "Retention.Time", "RT.Rank");
  
  if(!current.proc$mumRT ){
    matched_res <- matched_res[,-(5:6)]
  }
  
  #print(paste0(length(unique(matched_res[,2])), " matched compounds! cpd2mz"))
  
  # now create empirical compounds if necessary!
  # 1 compound matches to multiple m/z, filter by RT 
  if(current.proc$mumRT){
    start <- Sys.time()
    # mz, ion
    empirical.cpd.list <- split(matched_res[,c(1,3,5,6)], matched_res[,2]); # split mz, ion and rt by compound
    empirical.cpds2cpds <- vector(length=(length(empirical.cpd.list)), "list")
    names(empirical.cpds2cpds) <- names(empirical.cpd.list)
    
    # for each compound, if multiple matches, split into ECpds if > RT tolerance - rt_tol
    for(i in 1:length(empirical.cpd.list)){
      
      mzs <- empirical.cpd.list[[i]]$Query.Mass
      ions <- empirical.cpd.list[[i]]$Matched.Form
      rts <- empirical.cpd.list[[i]]$Retention.Time
      rt.rank <- empirical.cpd.list[[i]]$RT.Rank
      cpds <- names(empirical.cpd.list)[i]
      
      # first, for each compound, determine ECs among matched ions
      if(length(mzs)>1){ # if multiple ECs per compound
        
        # first group together to create empirical cpds by rt
        rts <- as.numeric(rts)
        names(rts) <- paste0(mzs, ";", ions, ";", rts, ";", cpds)
        rts <- sort(rts)
        
        # second, group together to create empirical cpds by rt rank
        rt.ranks <- as.numeric(rt.rank)
        names(rt.ranks) <- paste0(mzs, ";", ions, ";", rts, ";", cpds)
        rt.ranks <- sort(rt.ranks)
        
        split.inx <- c(0, cumsum(Reduce("&", list(abs(diff(rts)) > rt_tol, abs(diff(rt.ranks)) > rt_tol_rank) )))
        
        # need to deal w. multiple rts but only 1 EC
        if(length(unique(split.inx)) > 1){
          e.cpds <- split(rts, split.inx)
          empirical.cpds2cpds[[i]] <- lapply(e.cpds, names)
        }else{
          empirical.cpds2cpds[[i]] <- paste0(names(rts), collapse="__")
        }
        
      }else{ # if only 1 EC per compound
        empirical.cpds2cpds[[i]] <- paste0(mzs, ";", ions, ";", rts, ";", cpds)
      }
    }
    
    initial_ecs <- unlist(empirical.cpds2cpds, recursive=FALSE)
    names(initial_ecs) <- paste0("EC", 1:length(initial_ecs))
    print(paste0(length(initial_ecs), " initial ECs created!"))
    
    # second, merge ECs if same m/z and form - append compounds
    try <- reshape2::melt(initial_ecs)
    try2 <- strsplit(as.character(try[,1]), split="__", fixed=TRUE) # deals with multiple rts belonging to 1 EC
    try2.df <- data.frame(value=unlist(try2), L1 = rep(try$L1, sapply(try2, length)))
    
    info <- strsplit(as.character(try2.df[,1]), split=";")
    df_ecs <- data.frame(ec=as.character(try2.df[,2]), mz=sapply(info, `[[`, 1), form=sapply(info, `[[`, 2), rt = sapply(info, `[[`, 3), cpd = sapply(info, `[[`, 4), stringsAsFactors = F)
    df_ecs$str_row_inx <- paste(df_ecs$mz, df_ecs$form, df_ecs$rt, sep = "___")
    qs::qsave(df_ecs, "initial_ecs.qs")
    merged_ecs <- aggregate(. ~ str_row_inx, df_ecs, paste, collapse=";")
    
    # cleaning the df
    # merged_ecs$ec <- sapply(strsplit(merged_ecs$ec, ";", fixed=TRUE), function(x) unlist(x)[1]) - keep as long name
    merged_ecs$mz <- sapply(strsplit(merged_ecs$mz, ";", fixed=TRUE), function(x) unique(unlist(x)))
    merged_ecs$form <- sapply(strsplit(merged_ecs$form, ";", fixed=TRUE), function(x) unique(unlist(x)))
    merged_ecs$rt <- sapply(strsplit(merged_ecs$rt, ";", fixed=TRUE), function(x) unique(unlist(x)))
    print(paste0(length(unique(merged_ecs$ec)), " merged ECs identified!"))
    
    # third, check if primary ion is present
    # needs to be per EC!
    if(current.proc$primary_ion=="yes"){
   
      ecs <- unique(merged_ecs$ec);
      
      # function to group ECs and verify if contains primary ion
      new_info <- lapply(ecs, function(x) { 
        new_info <- merged_ecs[which(merged_ecs$ec == x),] # subset merged_ecs to rows containing ECx
        primary.inx <- length(intersect(new_info$form, primary_ions))
        
        if(primary.inx>0){
          new_info <- new_info
        }else{
          new_info <- NULL
        }
        new_info
      })  
      
      final_ecs <- do.call(args=new_info, what=rbind)[,-1]
      
    }else{
      final_ecs <- merged_ecs[,-1]
    }
    
    colnames(final_ecs) <- c("Empirical.Compound", "Query.Mass", "Matched.Form", "Retention.Time", "Matched.Compound")
    
    # transform to long format
    cpd_split <- strsplit(as.character(final_ecs$Matched.Compound), ";", fixed=TRUE)
    reps <- pmax(lengths(cpd_split))
    df2 <- final_ecs[rep(1:nrow(final_ecs), reps), 1:4]
    df2$Matched.Compound <- unlist(mapply(function(x,y) c(x, rep(NA, y)), cpd_split, reps-lengths(cpd_split)))
    
    matched_res <- merge(matched_res, df2)
    matched_res <- matched_res[,-6] #rm rt rank
    matched_res[,6] <- as.character(matched_res[,6])
    
    # now deal with the fact that if at least one EC overlap, need to count as same EC per compound...
    my_final_cpds <- aggregate(. ~ Matched.Compound, matched_res, paste, collapse="___")
    my_final_cpds_list <- lapply(split(my_final_cpds$Empirical.Compound, my_final_cpds$Matched.Compound), unlist)
    
    cpd2ec1 <- lapply(seq_along(my_final_cpds_list), function(x) { # function used to make grouping of ecs per cpd
      
      ecs <- unlist(strsplit(my_final_cpds_list[[x]], "___", fixed=TRUE))
      
      if(length(ecs) > 1){
        ecs.list <- as.list(strsplit(ecs, ";", fixed=TRUE))
        library(igraph)
        m = sapply(ecs.list, function(x) sapply(ecs.list, function(y) length(intersect(x,y))>0))
        g = igraph::groups(components(graph_from_adjacency_matrix(m)))
        ecs <- paste0(sapply(g, function(z) paste0(ecs[z], collapse = "|") ), collapse = "___")
      }
      ecs
    })
    
    names(cpd2ec1) <- names(my_final_cpds_list)
    
    update_ecs <- lapply(seq_along(cpd2ec1), function(z) {
      
      ecs.old <- unlist(strsplit(my_final_cpds_list[[z]], "___", fixed=TRUE))
      ecs.new <- unlist(strsplit(cpd2ec1[[z]], "___", fixed=TRUE))
      
      for(i in seq_along(ecs.new)){
        pattern <- ecs.new[i]
        pattern_vec <- unlist(strsplit(pattern, "\\|"))
        up.pattern <- paste0(unique(pattern_vec), collapse = "|")
        ecs.old[ ecs.old %in% pattern_vec  ] <- up.pattern
      }
      
      ecs.old <- paste0(ecs.old, collapse = "___")
      ecs.old
    })
    
    updated_ecs <- do.call(rbind, update_ecs)
    my_final_cpds$Empirical.Compound <- updated_ecs
    
    new_dt <- data.table::data.table(my_final_cpds)
    new_dt <- new_dt[, list(Query.Mass = unlist(strsplit(as.character(Query.Mass), "___", fixed=TRUE)), 
                            Matched.Form = unlist(strsplit(as.character(Matched.Form), "___", fixed=TRUE)),
                            Retention.Time = unlist(strsplit(as.character(Retention.Time), "___", fixed=TRUE)),
                            Mass.Diff = unlist(strsplit(as.character(Mass.Diff), "___", fixed=TRUE)),
                            Empirical.Compound = unlist(strsplit(as.character(Empirical.Compound), "___", fixed=TRUE))),
                     by = Matched.Compound]
    
    matched_res <- data.frame(Query.Mass = new_dt$Query.Mass, Matched.Compound = new_dt$Matched.Compound, Matched.Form = new_dt$Matched.Form,
                              Retention.Time = new_dt$Retention.Time, Mass.Diff = new_dt$Mass.Diff, Empirical.Compound = new_dt$Empirical.Compound, stringsAsFactors = FALSE)
    
    # make EC names
    ec <- matched_res$Empirical.Compound
    ec.unique <- unique(matched_res$Empirical.Compound)
    
    for(i in seq_along(ec.unique)){
      ec <- replace(ec, grep(paste0("\\b", ec.unique[i], "\\b"), ec, perl=TRUE), paste0("EC000", i))
    }
    
    matched_res$Empirical.Compound <- gsub("\\|.*", "", ec)
    end <- Sys.time()
    totaltime <- end-start
    print(paste0(length(unique(matched_res$Empirical.Compound)), " empirical compounds identified in ", totaltime, " seconds."))
  }
  
  #matched_res <- matched_res[which(matched_res$Matched.Form %in% primary_ions),]
  fast.write(matched_res, file="mummichog_matched_compound_all.csv", row.names=FALSE);
  qs::qsave(matched_res, "mum_res.qs");
  
  # now update expr. profile
  matched_mz <- matched_res[,1];
  matched_ts <- current.proc$expr_dic[matched_mz];
  
  if(current.proc$mumRT) { # RT need to be in EC space
    # first create ecpd to expression dict
    ec.exp.mat <- data.frame(key=matched_res[,6], value=as.numeric(matched_ts), stringsAsFactors = F)
    ec_exp_dict <- Convert2Dictionary(ec.exp.mat);
    ec.exp.vec <- unlist(lapply(ec_exp_dict, max));
    
    # also need to make cpd_exp_dict for KEGG network view
    exp.mat <- data.frame(key=matched_res[,2], value=as.numeric(matched_ts));
    cpd_exp_dict <- Convert2Dictionary(exp.mat);
    
    # ecpd to cpd dict
    cpd_ecpd_dict <- Convert2Dictionary(matched_res[,c(2,6)])
    ecpd_cpd_dict <- Convert2Dictionary(matched_res[,c(6,2)])
    
    # now mz 2 ecpd dict
    mz2cpd_dict <- Convert2Dictionary(matched_res[,c(1,2)]); #indexed/named by mz
    mz2ec_dict <- Convert2Dictionary(matched_res[,c(1,6)])
    ec2mz_dict <- Convert2Dictionary(matched_res[,c(6,1)])
    
   # label.mat <- data.frame(key=matched_res[,2], label= apply(matched_res[,c(1,4)], 1, function(x) paste(x,collapse = "__")));
   # cpd_label <- Convert2Dictionary(label.mat)
    # save to current.proc
    current.proc$ec_exp_dict <- ec_exp_dict
    current.proc$cpd_exp_dict <- cpd_exp_dict;
    current.proc$ec_exp <- ec.exp.vec
    current.proc$mz2cpd_dict <- mz2cpd_dict;
    current.proc$mz2ec_dict <- mz2ec_dict
    current.proc$ec2mz_dict <- ec2mz_dict
    current.proc$ecpd_cpd_dict <- ecpd_cpd_dict
    current.proc$cpd_ecpd_dict <- cpd_ecpd_dict
    current.proc$cpd_ecpd_counts <- cpd2ec1
    current.proc$keggNet$Query <- apply(matched_res[,c(1,4)], 1, function(x) paste(x,collapse = "__"))
    current.proc$keggNet$Match <- matched_res[,2]
    current.proc$keggNet$adduct <- matched_res[,3]
    # now do matching to identify significant input_ecpdlist
    refmz <- names(mz2ec_dict)
    hits.index <- which(refmz %in% as.character(current.proc$sig.feat));
    ec1 <- unique(unlist(mz2ec_dict[hits.index]));
    current.proc$input_ecpdlist <- ec1;
    current.proc$total_matched_ecpds <- unique(as.vector(matched_res$Empirical.Compound));
    
  } else {
    # get the expression profile for each 
    exp.mat <- data.frame(key=matched_res[,2], value=as.numeric(matched_ts));
    cpd_exp_dict <- Convert2Dictionary(exp.mat);
    # create average exp
    exp.vec <- unlist(lapply(cpd_exp_dict, mean));
    
    # now need to get the mapping from mz to compound id (one mz can have 0, 1, or more id hits)
    mz2cpd_dict <- Convert2Dictionary(matched_res[,c(1,2)]); #indexed/named by mz
    cpd2mz_dict <- Convert2Dictionary(matched_res[,c(2,1)]); # indexed/named by id
    
    # now do matching to identify significant input_cpdlist
    refmz <- names(mz2cpd_dict)
    hits.index <- which(refmz %in% as.character(current.proc$ref_mz));
    cpd1 <- unique(unlist(mz2cpd_dict[hits.index]));
    
    current.proc$mz2cpd_dict <- mz2cpd_dict;
    current.proc$cpd_exp_dict <- cpd_exp_dict;
    current.proc$cpd_exp <- exp.vec;
    current.proc$cpd2mz_dict <- cpd2mz_dict;
    current.proc$input_cpdlist <- cpd1;
    current.proc$total_matched_cpds <- unique(as.vector(matched_res$Matched.Compound));
    current.proc$keggNet$Query <- matched_res[,1]
    current.proc$keggNet$Match <- matched_res[,2]
    current.proc$keggNet$adduct <- matched_res[,3]
  }
  
  form.mat <- cbind(matched_res[,2], matched_res[,3]);
  cpd_form_dict <- Convert2Dictionary(form.mat);
  current.proc$cpd_form_dict <- cpd_form_dict;
  
  current.proc<<-current.proc
  return(1);
}

#' Internal function to perform PSEA, with RT
#' @noRd
.init.RT.Permutations <- function(permNum){
  current.proc <- .perform.mummichogRTPermutations(current.proc,permNum=100);
  current.proc <<- .compute.mummichogRTSigPvals( current.proc);
  
}

# Internal function for permutation
.perform.mummichogRTPermutations <- function(current.proc, permNum){
  
  print(paste('Resampling, ', permNum, 'permutations to estimate background ...'));
  permutation_hits <- permutation_record <- vector("list", permNum);
  matched_res <- qs::qread("mum_res.qs");
  set.seed(123)
  for(i in 1:permNum){ # for each permutation, create list of input emp compounds and calculate pvalues for each pathway
    input_mzlist <- unlist(sample(current.proc$ref_mz, length(current.proc$sig.feat)))
    t <- make_ecpdlist(input_mzlist);
    perm <- ComputeMummichogRTPermPvals(t, current.proc$total_matched_ecpds, current.proc$pathways, matched_res, input_mzlist);
    permutation_record[[i]] <- perm[1]
    permutation_hits[[i]] <- perm[2]
  }
  
  # append new info
  current.proc$perm_record <- permutation_record;
  current.proc$perm_hits <- permutation_hits;
  print("perm done")
  return(current.proc);
}


#' Internal function to perform PSEA, no retention time
#' @importFrom data.table setDF
#' @noRd
.init.Permutations <- function(permNum){

  .perform.mummichogPermutations(permNum);
  .compute.mummichogSigPvals();

}
 

# Calculate p-values for each Lperm
# Used in higher mummichogR functions w. RT
ComputeMummichogRTPermPvals <- function(input_ecpdlist, total_matched_ecpds, pathways, matches.res, input_mzlist){

  ora.vec <- input_ecpdlist; #Lperm
  query_set_size <- length(ora.vec) # query set size
  current.mset <- pathways$emp_cpds; #all
  total_ecpds <- unique(total_matched_ecpds) # matched number of empirical compounds
  total_feature_num <- length(total_ecpds)
  
  size <- negneg <- vector(mode="list", length=length(current.mset));
  
  ecpds <- lapply(current.mset, function(x) intersect(x, total_ecpds)); # pathways & all ref ecpds
  feats <- lapply(current.mset, function(x) intersect(x, ora.vec)); #pathways & query ecpds (perm lsig)
  feat_len <- unlist(lapply(feats, length)); # length of overlap features
  set.num <- unlist(lapply(ecpds, length)); #cpdnum
  
  negneg <- sizes <- vector(mode="list", length=length(current.mset));
  
  for(i in seq_along(current.mset)){ # for each pathway
    sizes[[i]] <- feat_len[i] # for ecs, just use length of overlap feats - overlap_size
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size;
  }
  
  unsize <- as.integer(unlist(sizes))
  res.mat <- matrix(0, nrow=length(current.mset), ncol=1)
  fishermatrix <- cbind(unsize-1, set.num, (query_set_size + unlist(negneg)), query_set_size)
  res.mat[,1] <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail=FALSE));
  perm_records <- list(res.mat, as.matrix(unsize));
  return(perm_records);
}


.compute.mummichogRTSigPvals <- function(current.proc){

  qset <- unique(unlist(current.proc$input_ecpdlist)); #Lsig ora.vec
  query_set_size <- length(qset); #q.size
  input_cpd <- unique(unlist(current.proc$ecpd_cpd_dict[qset]));
  
  total_ecpds <- unique(current.proc$total_matched_ecpds) #all matched compounds
  total_feature_num <- length(total_ecpds)
  total_cpds <- unique(unlist(current.proc$ecpd_cpd_dict));
  
  current.mset <- current.proc$pathways$emp_cpds; #all compounds per pathway
  path.num <- unlist(lapply(current.mset, length));
  
  ecpds <- lapply(current.mset, function(x) intersect(x, total_ecpds)); #pathways & all ref ecpds
  set.num <- unlist(lapply(ecpds, length)); # total ecpd num in pathway
  cpd.num<- unlist(lapply(lapply(current.proc$pathways$cpds, function(x) intersect(x, total_cpds)), length)); 
  
  feats <- lapply(current.mset, function(x) intersect(x, qset)); #pathways & lsig
  feat_len <- unlist(lapply(feats, length)); # length of overlap features
  cpd_len<- lapply(current.proc$pathways$cpds, function(x) intersect(x, input_cpd))
  feat_vec <- sapply(cpd_len, function(x) paste(x, collapse=";"))
  
  negneg <- sizes <- vector(mode="list", length=length(current.mset)); #empty lists
  
  for(i in seq_along(current.mset)){ # for each pathway
    sizes[[i]] <- feat_len[i] # overlap size
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size; # failure in left part
  }
  
  #error fixing for negatives, problem occurs when total_feat_num and query_set_size too close (lsig too close to lall)
  negneg <- rapply(negneg, function(x) ifelse(x<0,0,x), how = "replace") 
  
  unsize <- as.integer(unlist(sizes));
  
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));
  
  # prepare for the result table
  res.mat <- matrix(0, nrow=length(current.mset), ncol=8);
  
  #fishermatrix for phyper
  fishermatrix <- cbind(unsize-1, set.num, (query_set_size + unlist(negneg) - unsize), query_set_size);
  first <- unlist(lapply(sizes, function(x) max(0, x-1)));
  easematrix <- cbind(first, (set.num - unsize + 1), (query_set_size - unsize), unlist(negneg)); 
  
  
  res.mat[,1] <- unlist(lapply(current.proc$pathways$cpds, length));  
  res.mat[,2] <- cpd.num;
  res.mat[,3] <- as.integer(unlist(lapply(cpd_len, length)));
  res.mat[,4] <- query_set_size*(path.num/uniq.count); #expected
  res.mat[,5] <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail=FALSE));
  res.mat[,6] <- apply(easematrix, 1, function(x) fisher.test(matrix(x, nrow=2), alternative = "greater")$p.value);
  res.mat[,7] <- set.num
  res.mat[,8] <- unsize
  
  colnames(res.mat) <- c("Pathway total", "Hits.total", "Hits.sig", "Expected", "FET", "EASE","Total.EC","Sig.EC");
  rownames(res.mat) <- names(current.proc$pathways$emp_cpds)
  
  current.proc$pvals <- res.mat;
  permutations_hits <- matrix(unlist(current.proc$perm_hits), nrow=length(current.proc$perm_hits), byrow=TRUE);
  sig_hits <- unsize; # sighits
  sigpvalue <- res.mat[,5]; # EASE scores
  
  perm_record <- unlist(current.proc$perm_record);
  perm_minus <- abs(0.9999999999 - perm_record);
  
  if(length(sig_hits[sig_hits!=0]) < round(length(sig_hits)*0.05)){ # too few hits that can't calculate gamma dist!
    if(!exists("adjustedp")){
      adjustedp <- rep(NA, length = length(res.mat[,1]))
    }
    res.mat <- cbind(res.mat, Gamma=adjustedp);
  }else{
    tryCatch({
      fit.gamma <- fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle", lower = c(0, 0), start = list(scale = 1, shape = 1));
      rawpval <- as.numeric(sigpvalue);
      adjustedp <- 1 - (pgamma(1-rawpval, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["scale"]));
    }, error = function(e){
  
      print(e)   
    }, finally = {
      if(!exists("adjustedp")){
        adjustedp <- rep(NA, length = length(res.mat[,1]))
      }
      res.mat <- cbind(res.mat, Gamma=adjustedp);
    })
  }
  
  
  #calculate empirical p-values
  record <- current.proc$perm_record
  
  
  fisher.p <- as.numeric(res.mat[,5])
  
  #pathway in rows, perms in columns
  record_matrix <- do.call(cbind, do.call(cbind, record))
  num_perm <- ncol(record_matrix)
  
  #number of better hits for web
  better.hits <- sapply(seq_along(record_matrix[,1]), function(i) sum(record_matrix[i,] <= fisher.p[i])  )
  
  #account for a bias due to finite sampling - Davison and Hinkley (1997)
  emp.p <- sapply(seq_along(record_matrix[,1]), function(i) (sum(record_matrix[i,] <= fisher.p[i])/num_perm) )
  
  res.mat <- cbind(res.mat, Emp.Hits = better.hits, Empirical = emp.p, Sig.Cmpd.Hits = feat_vec)
  # remove pathways with no hits
  hit.inx <- as.numeric(as.character(res.mat[,8])) > 0;
  res.mat <- res.mat[hit.inx, , drop=FALSE];
  
  if(nrow(res.mat) <= 1){
    AddErrMsg("Not enough m/z to compound hits for pathway analysis! Try Version 1 (no RT considerations)!")
    return(0)
  }
  
  # prepare json element for network
  # need to convert ecpds to cpds
  # and get average expression based on ec
  cpds <- lapply(ecpds, function(x) unique(unlist(current.proc$ecpd_cpd_dict[match(x, names(current.proc$ecpd_cpd_dict))])) )  
  cpd.exp.vec <- sapply(ecpds, function(x) mean(current.proc$ec_exp[match(x, names(current.proc$ec_exp))]) )
  cpds_feats <- lapply(feats, function(x) unique(unlist(current.proc$ecpd_cpd_dict[match(x, names(current.proc$ecpd_cpd_dict))])) )  
  
  # now make exp vec for all compounds
  cpds2ec <- current.proc$cpd_ecpd_dict
  cpds.all <- unique(unlist(current.proc$ecpd_cpd_dict[match(total_ecpds, names(current.proc$ecpd_cpd_dict))]))
  cpd.exp.vec <- sapply(cpds.all, function(x) sapply(seq_along(x), function(i) mean(current.proc$ec_exp[match(unique(unlist(cpds2ec[match(x[[i]], names(cpds2ec))])), names(current.proc$ec_exp))]) ) )
  
  hits.all <- cpds[hit.inx];
  hits.sig <- cpds_feats[hit.inx];  
  path.nms <- names(current.proc$pathways$emp_cpds)[hit.inx];
  
  metInfo <- qs::qread(paste0(lib.path.mmp,"general_kegg2name.qs"));
  hits.nms <- lapply(hits.all, function(x){
  x=metInfo$Name[match(x,metInfo$ID)]  
  return(x)
  })
  hits.node <- lapply(hits.all, function(x){
  x=metInfo$id[match(x,metInfo$ID)]  
  return(x)
  })

 # order by p-values
  if(length(na.omit(res.mat[,9])) == 0){
    ord.inx <- order(res.mat[,5]); # order by FET if gamma not able to be calc
  }else{
    ord.inx <- order(res.mat[,9]); # order by gamma
  }
  
  Sig.Cmpd.Hits = res.mat[ord.inx, 12]
  res.mat <- signif(apply(as.matrix(res.mat[ord.inx, 1:11]), 2, as.numeric), 5); # loop through columns and keep rownames
  rownames(res.mat) <- path.nms[ord.inx]
  hit.num = paste0(res.mat[,2],"/",res.mat[,1]);if(length(hit.num) ==1) { hit.num <- matrix(hit.num) }
  
  current.proc$mummi.resmat <- res.mat[,-11];
  current.proc$path.nms <- path.nms[ord.inx]
  current.proc$path.hits <- convert2JsonList(hits.all[ord.inx])
  current.proc$path.pval <- as.numeric(res.mat[,5])
  matched_res <- qs::qread("mum_res.qs");
  expr.mat <- lapply(current.proc$cpd_exp_dict,function(x) mean(as.numeric(x)))
  
  json.res <- list(
    hits.query = hits.all[ord.inx],
    path.nms = path.nms[ord.inx],
    path.pval = as.numeric(res.mat[,5]),
    hit.num = hit.num,
    path.fdr =  as.numeric(res.mat[,9]),
    hits.query.nm = hits.nms[ord.inx],
    hits.node = hits.node[ord.inx],
    hits.sig = convert2JsonList(hits.sig[ord.inx]),
    peakToMet = current.proc$cpd_label,
     expr.mat=convert2JsonList(expr.mat)
  );
  
  matri = cbind(res.mat[,-c(7:8)], paste0("P", seq.int(1, nrow(res.mat))))
  colnames(matri)[ncol(matri)] = "Pathway Number"
  matri <- cbind(matri, Sig.Cmpd.Hits)
  fast.write(matri, file="peak_enrichment_res.csv", row.names=TRUE);
  json.mat <- RJSONIO::toJSON(json.res);
  sink("network_enrichment_peak_0.json");
  cat(json.mat);
  sink();
  return(current.proc);
}
# Internal function for significant p value 
.compute.mummichogSigPvals <- function(current.proc){
  
  qset <- unique(unlist(current.proc$input_cpdlist)); #Lsig ora.vec
  query_set_size <- length(qset); #q.size
  
  total_cpds <- unique(current.proc$total_matched_cpds) #all matched compounds
  total_feature_num <- length(total_cpds)
  
  current.mset <- current.proc$pathways$cpds; #all compounds per pathway
  path.num <- unlist(lapply(current.mset, length));
  
  cpds <- lapply(current.mset, function(x) intersect(x, total_cpds)); #pathways & all ref cpds
  set.num <- unlist(lapply(cpds, length)); #cpdnum
  
  feats <- lapply(current.mset, function(x) intersect(x, qset)); #pathways & lsig
  feat_len <- unlist(lapply(feats, length)); # length of overlap features
  feat_vec <- sapply(feats, function(x) paste(x, collapse=";"))
  
  negneg <- sizes <- vector(mode="list", length=length(current.mset)); #empty lists
  
  for(i in seq_along(current.mset)){ # for each pathway
    sizes[[i]] <- min(feat_len[i], count_cpd2mz(current.proc$cpd2mz_dict, unlist(feats[i]), current.proc$dataSet$input_mzlist)) #min overlap or mz hits
    negneg[[i]] <- total_feature_num + sizes[[i]] - set.num[i] - query_set_size; # failure in left part
  }
  
  #error fixing for negatives, problem occurs when total_feat_num and query_set_size too close (lsig too close to lall)
  negneg <- rapply(negneg, function(x) ifelse(x<0,0,x), how = "replace") 
  
  unsize <- as.integer(unlist(sizes));
  
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));
  
  # prepare for the result table
  res.mat <- matrix(0, nrow=length(current.mset), ncol=8);
  
  #fishermatrix for phyper
  fishermatrix <- cbind(unsize-1, set.num, (query_set_size + unlist(negneg) - unsize), query_set_size); 
  first <- unlist(lapply(sizes, function(x) max(0, x-1)));
  easematrix <- cbind(first, (set.num - unsize + 1), (query_set_size - unsize), unlist(negneg)); 
  
  res.mat[,1] <- path.num;  
  res.mat[,2] <- set.num;
  res.mat[,3] <- unsize;
  res.mat[,4] <- query_set_size*(path.num/uniq.count); #expected
  res.mat[,6] <- apply(easematrix, 1, function(x) fisher.test(matrix(x, nrow=2), alternative = "greater")$p.value);
  res.mat[,5] <- apply(fishermatrix, 1, function(x) phyper(x[1], x[2], x[3], x[4], lower.tail=FALSE));
  res.mat[,7] <- set.num
  res.mat[,8] <- unsize
  colnames(res.mat) <- c("Pathway total", "Hits.total", "Hits.sig", "Expected", "FET", "EASE","Total","Sig");
  rownames(res.mat) <- current.proc$pathways$name
  
  current.proc$pvals <- res.mat;
  permutations_hits <- matrix(unlist(current.proc$perm_hits), nrow=length(current.proc$perm_hits), byrow=TRUE);
  sig_hits <- res.mat[,3]; # sighits
  sigpvalue <- res.mat[,6]; # EASE scores
  
  perm_record <- unlist(current.proc$perm_record);
  perm_minus <- abs(0.9999999999 - perm_record);
  
  if(length(sig_hits[sig_hits!=0]) < round(length(sig_hits)*0.05)){ # too few hits that can't calculate gamma dist!
    if(!exists("adjustedp")){
      adjustedp <- rep(NA, length = length(res.mat[,1]))
    }
    res.mat <- cbind(res.mat, Gamma=adjustedp);
  }else{
    tryCatch({
      fit.gamma <- fitdistrplus::fitdist(perm_minus, distr = "gamma", method = "mle", lower = c(0, 0), start = list(scale = 1, shape = 1));
      rawpval <- as.numeric(sigpvalue);
      adjustedp <- 1 - (pgamma(1-rawpval, shape = fit.gamma$estimate["shape"], rate = fit.gamma$estimate["scale"]));
    }, error = function(e){
      print(e)   
    }, finally = {
      if(!exists("adjustedp")){
        adjustedp <- rep(NA, length = length(res.mat[,1]))
      }
      res.mat <- cbind(res.mat, Gamma=adjustedp);
    })
  }
  
  #calculate empirical p-values
  record <- current.proc$perm_record
  fisher.p <- as.numeric(res.mat[,5])
  
  #pathway in rows, perms in columns
  record_matrix <- do.call(cbind, do.call(cbind, record))
  num_perm <- ncol(record_matrix)
  
  #number of better hits for web
  better.hits <- sapply(seq_along(record_matrix[,1]), function(i) sum(record_matrix[i,] <= fisher.p[i])  )
  
  #account for a bias due to finite sampling - Davison and Hinkley (1997)
  emp.p <- sapply(seq_along(record_matrix[,1]), function(i) (sum(record_matrix[i,] <= fisher.p[i])/num_perm) )
  
  res.mat <- cbind(res.mat, Emp.Hits=better.hits, Empirical=emp.p, Cpd.Hits = feat_vec)
  
  # remove those no hits
  hit.inx <- as.numeric(as.character(res.mat[,3])) > 0;
  res.mat <- res.mat[hit.inx, , drop=FALSE];
  
  if(nrow(res.mat) <= 1){
    AddErrMsg("Not enough m/z to compound hits for pathway analysis!")
    return(0)
  }

  # prepare json element for network
  hits.all <- cpds[hit.inx];
  hits.sig <- feats[hit.inx];  
  path.nms <- names(current.proc$pathways$cpds)[hit.inx];
  
  # order by p-values
  ord.inx <- order(res.mat[,9]);
  
  Cpd.Hits <- res.mat[ord.inx, 12]
  res.mat <- signif(apply(as.matrix(res.mat[ord.inx, 1:11, drop=FALSE]), 2, as.numeric), 5);
  rownames(res.mat) <- path.nms[ord.inx]
  current.proc$mummi.resmat <- res.mat[,-11];
  
  current.proc$path.nms <- path.nms[ord.inx]
  current.proc$path.hits <- convert2JsonList(hits.all[ord.inx])
  current.proc$path.pval <- as.numeric(res.mat[,9])
  matched_res <- qs::qread("mum_res.qs"); 
    json.res <- list(
    hits.query = hits.all[ord.inx],
    path.nms = path.nms[ord.inx],
    path.pval = as.numeric(res.mat[,5]),
    hit.num = hit.num,
    path.fdr =  as.numeric(res.mat[,9]),
    hits.sig = convert2JsonList(hits.sig[ord.inx]),
    peakToMet = current.proc$cpd_label
  )
    
  matri = cbind(res.mat[,-c(7:8)], paste0("P", seq.int(1, nrow(res.mat))))
  colnames(matri)[ncol(matri)] = "Pathway Number"
  matri <- cbind(matri, Cpd.Hits)
  fast.write(matri, file="peak_enrichment_res.csv", row.names=TRUE);
  json.mat <- RJSONIO::toJSON(json.res);
  sink("network_enrichment_peak_0.json");
  cat(json.mat);
  sink();
  return(current.proc);
}




mz_tolerance <- function(mz, ms.type){
  return(ms.type*1e-06*mz)
}

myFastList <- function(capacity = 50) {
  buffer <- vector('list', capacity)
  names <- character(capacity)
  length <- 0
  methods <- list()
  
  methods$double.size <- function() {
    buffer <<- c(buffer, vector('list', capacity))
    names <<- c(names, character(capacity))
    capacity <<- capacity * 2
  }
  
  methods$add <- function(name, val) {
    if(length == capacity) {
      methods$double.size()
    }
    
    length <<- length + 1
    buffer[[length]] <<- val
    names[length] <<- name
  }
  
  methods$as.list <- function() {
    b <- buffer[0:length]
    names(b) <- names[0:length]
    return(b)
  }
  
  methods
}

Convert2Dictionary <- function(data, quiet=T){
  
  all.ids <- data[,1];
  dup.inx <- duplicated(all.ids);
  if(sum(dup.inx) > 0){
    uniq.ids <- all.ids[!dup.inx];
    uniq.vals <- data[!dup.inx,2];
    
    # convert two-col data it to list (vals as list values, ids as list names)
    uniq.list <- split(uniq.vals, uniq.ids)
    
    # the list element orde will be sorted by the names alphabetically, need to get updated ones
    uniq.id.list <- names(uniq.list)
    
    dup.ids <- all.ids[dup.inx];
    uniq.dupids <- unique(dup.ids);
    uniq.duplen <- length(uniq.dupids);
    
    for(id in uniq.dupids){ # only update those with more than one hits
      hit.inx.all <- which(all.ids == id);
      hit.inx.uniq <- which(uniq.id.list == id);
      uniq.list[[hit.inx.uniq]]<- data[hit.inx.all,2];
    }
    
    current.msg <- paste("A total of ", sum(dup.inx), " of duplicates were merged.", sep="");
    return(uniq.list);
  }else{
    current.msg <- "All IDs are unique.";
    uniq.list <- split(data[,2], data[,1]);
    return(uniq.list);
  }
}

make_ecpdlist <- function(input_mzs){
  ecpd <- unique(unlist(current.proc$mz2ec_dict[input_mzs]));
  ecpd <- ecpd[!is.null(ecpd)];
  return(ecpd);
}



 primary_ions <- c('M+H[1+]', 'M+Na[1+]', 'M-H2O+H[1+]', 'M-H[-]', 'M-2H[2-]', 'M-H2O-H[-]',
                      'M+H [1+]', 'M+Na [1+]', 'M-H2O+H [1+]', 'M-H [1-]', 'M-2H [2-]', 'M-H2O-H [1-]')
  
create.adducts <- function(cpd.lib){
  
  cpd.lib <- cpd.lib;
  ms_modes <- c('dpj_positive', 'positive', 'negative');
  adducts <- list();
  
  for (ms_mode in ms_modes){
    adducts[[ms_mode]] <- Compound_function_mzlist(ms_mode, cpd.lib$mw);
  }
  
  cpd.lib$adducts <- adducts;
  
  # create a dictionary for look up in the range of 50-2000
  # now need to create ladder (tree) for each new mz
  # key is the mass 50 to 2000, values are the compounds (if any of their modified mw gives the value)
  # now create cpd tree for each mass pos
  # note, this can be slow, but this can be created before hand
  # for each species and for each mode
  # note l2 only stores the index of the cpd.lib
  
  cpd.tree <- list();
  for (ms_mode in ms_modes){
    l2 <- list();
    l2[[49]] <- "";
    l2[[2001]] <- "";
    mz.mat <- cpd.lib$adducts[[ms_mode]];
    floor.mzs <- floor(mz.mat);
    for(i in 1:nrow(floor.mzs)){
      neighbourhood <- floor.mzs[i,];
      for(n in neighbourhood){
        if((n>50) & (n<2000)){
          l2[[n]] <- append(l2[[n]], i);
        }
      }
    }
    cpd.tree[[ms_mode]] <- lapply(l2, unique);
  }
  
  # set up the variables
  mummichog.lib <- list(
    cpd.tree = cpd.tree,
    cpd.lib = cpd.lib
  )
  
  message("Library of CMN compound universe created!")
  file_name <- "cmn_compound_lib.qs"
  
  qs::qsave(mummichog.lib, file=file_name);
}

Compound_function_mzlist <- function(ms_mode, mw){
  
 library(stringr)
  
  PROTON <- 1.00727646677;
  mw_modified <- NULL;
  
  if (ms_mode == "dpj_positive"){
    mw_modified <- cbind(mw, mw + PROTON, mw/2 + PROTON, mw +1.0034 + PROTON, mw/2 + 0.5017 + PROTON, mw +1.9958 + PROTON, mw +1.9972 + PROTON, mw + 21.9820 + PROTON, mw/2 + 10.991 + PROTON, mw + 37.9555 + PROTON, mw + 67.9874 + PROTON, mw + 83.9613 + PROTON);
    colnames(mw_modified) <- c('M[1+]', 'M+H[1+]', 'M(C13)+H[1+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 'M(S34)+H[1+]', 'M(Cl37)+H[1+]', 'M+Na[1+]', 'M+H+Na[2+]', 'M+K[1+]', 'M+HCOONa[1+]', 'M+HCOOK[1+]');
    
  }else if (ms_mode == "positive" | ms_mode == 'generic'){
    mw_modified <- cbind(mw, mw + PROTON, mw/2 + PROTON, mw/3 + PROTON, mw +1.0034 + PROTON, mw/2 + 0.5017 + PROTON, mw/3 + 0.3344 + PROTON, mw +1.9958 + PROTON, mw +1.9972 + PROTON, mw + 21.9820 + PROTON, mw/2 + 10.991 + PROTON, mw + 37.9555 + PROTON, mw + 18.0106 + PROTON, mw - 18.0106 + PROTON, mw - 36.0212 + PROTON, mw - 17.0265 + PROTON, mw - 27.9950 + PROTON, mw - 43.9898 + PROTON, mw - 46.0054 + PROTON, mw + 67.9874 + PROTON, mw - 67.9874 + PROTON, mw + 57.9586 + PROTON, mw - 72.0211 + PROTON, mw + 83.9613 + PROTON, mw - 83.9613 + PROTON);
    colnames(mw_modified) <- c('M[1+]', 'M+H[1+]', 'M+2H[2+]', 'M+3H[3+]', 'M(C13)+H[1+]', 'M(C13)+2H[2+]', 'M(C13)+3H[3+]', 'M(S34)+H[1+]', 'M(Cl37)+H[1+]', 'M+Na[1+]', 'M+H+Na[2+]', 'M+K[1+]', 'M+H2O+H[1+]', 'M-H2O+H[1+]', 'M-H4O2+H[1+]', 'M-NH3+H[1+]', 'M-CO+H[1+]', 'M-CO2+H[1+]', 'M-HCOOH+H[1+]', 'M+HCOONa[1+]', 'M-HCOONa+H[1+]', 'M+NaCl[1+]', 'M-C3H4O2+H[1+]', 'M+HCOOK[1+]', 'M-HCOOK+H[1+]');
    
  }else if (ms_mode == "negative"){
    mw_modified <- cbind(mw - PROTON, mw/2 - PROTON, mw + 1.0034 - PROTON, mw + 1.9958 - PROTON, mw + 1.9972 - PROTON, mw + 21.9820 - 2*PROTON, mw + 37.9555 - 2*PROTON, mw - 18.0106 - PROTON, mw + 34.9689, mw + 36.9659, mw + 78.9183, mw + 80.9163, mw + 2*12 + 3*1.007825 + 14.00307 - PROTON, mw + 1.007825 + 12 + 2*15.99491, mw + 3*1.007825 + 2*12 + 2*15.99491, mw - PROTON + 15.99491);
    colnames(mw_modified) <- c('M-H[-]', 'M-2H[2-]', 'M(C13)-H[-]', 'M(S34)-H[-]', 'M(Cl37)-H[-]', 'M+Na-2H[-]', 'M+K-2H[-]', 'M-H2O-H[-]', 'M+Cl[-]', 'M+Cl37[-]', 'M+Br[-]', 'M+Br81[-]', 'M+ACN-H[-]', 'M+HCOO[-]', 'M+CH3COO[-]', 'M-H+O[-]');
    
  }else{
    print("Unrecognized mode of instrumentation.")
  }
  
  return(mw_modified);
}
