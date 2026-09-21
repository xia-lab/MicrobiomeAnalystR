############################
###########TSEA#############
############################

GetNameMapCol <-function(mbSetObj, colInx){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$resTable[,colInx]);
}

GetTseaRowNames <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(rownames(mbSetObj$analSet$tseaInfo));
}

GetTseaCol <-function(mbSetObj, colInx){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$resTable[,colInx]);
}

#'Function to set up data for TSEA
#'@description This function sets up data for TSEA.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
Setup.MapData<-function(mbSetObj, qvec, method="ora"){

  lines <- unlist(strsplit(qvec, "\r|\n|\r\n")[1]);
  if(substring(lines[1],1,1)=="#"){
    lines <- lines[-1];
  }
  lines <- lines[lines != ""];
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet$tsea.method <- method;
  if(method == "gsea"){
    parts <- strsplit(lines, "\t");
    mbSetObj$dataSet$species <- trimws(sapply(parts, `[`, 1));
    mbSetObj$dataSet$scores <- as.numeric(sapply(parts, `[`, 2));
  } else {
    mbSetObj$dataSet$species <- lines;
  }
  return(.set.mbSetObj(mbSetObj))
}

#'Getter function
#'@description This function retrieves table from mbSetObj.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import xtable
GetORATable<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);  
  suppressMessages(library(xtable));   
  res <- mbSetObj$analSet$ora.mat;
  print(xtable::xtable(res, caption="Result from Over Representation Analysis"),
        tabular.environment = "longtable", caption.placement="top", size="\\scriptsize");
}

#'Perform cross referencing.
#'@description This function performs cross referencing of user's data
#'with the MicrobiomeAnalyst database. Given a list of species names or ids, 
#'it finds matched names or ids from selected internal databases.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export

CrossReferencing <- function(mbSetObj, q.type){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  # record all the data
  name.map <<- list();
  
  # distribute job
  mbSetObj$dataSet$q.type <- q.type;  
  .set.mbSetObj(mbSetObj)
  
  qvec <- mbSetObj$dataSet$species;
  resTable <- SpeciesMappingExact(qvec, q.type);
  mbSetObj$analSet$resTable <- resTable;
  mbSetObj$analSet$mapTable <- cbind(Query=qvec, resTable[,2:ncol(resTable)]);

  # do some sanity check. note name.map is on global env.
  if(length(which(is.na(name.map$hit.inx)))/length(name.map$hit.inx) > 0.75){
    nmcheck.msg <<- c(1, "Over 3/4 of the IDs could not be matched to our database. Please make 
                        sure that correct taxonomy IDs or common taxa names are used.");        
  }else{
    nmcheck.msg <<- c(1, "Name matching OK, please inspect (and manual correct) the results then proceed.");   
  }  
  return(.set.mbSetObj(mbSetObj))
}

# Utility function
# Mapping from different metabolite IDs
# For compound names to other id, can do exact or approximate match
# For other IDs, except HMDB ID, all other may return multiple /non-unique hits
# multiple hits or non-unique hits will all users to manually select
SpeciesMappingExact<-function(qvec, q.type){
       
  # local variable to save memory
  species.db <- .read.microbiomeanalyst.lib.rds("microbe_db_new.rds", "tsea")
       
  # variables to record results
  hit.inx = vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
  match.values = vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
  match.state = vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0 
       
  if(q.type == "gold"){
    hit.inx <- match(tolower(qvec), tolower(species.db$GOLD_ID));
    match.values <- species.db$GOLD_ID[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type %in% c("mixed","species","strain")){
    hit.inx <- match(tolower(qvec), tolower(species.db$taxa));
    # not spelled exactly as in the table: try the tidied / current NCBI name
    miss <- is.na(hit.inx);
    if(any(miss)){
      can <- .canonical.taxon.name(qvec[miss]);
      alt <- match(tolower(can), tolower(species.db$canonical));
      alt[is.na(alt)] <- match(tolower(can), tolower(species.db$taxa))[is.na(alt)];
      hit.inx[miss] <- alt;
    }
    match.values <- species.db$taxa[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else if(q.type == "ncbitax"){
    hit.inx <- match(tolower(qvec), tolower(species.db$NCBITAX));
    match.values <- species.db$NCBITAX[hit.inx];
    match.state[!is.na(hit.inx)] <- 1;
  }else{
    print(paste("Unknown species ID type:", q.type));
  }
      
  # empty memory
  name.map$hit.inx <- hit.inx;
  name.map$hit.values <- match.values;
  name.map$match.state <- match.state;
  name.map <<- name.map;

  # style for highlighted background for unmatched names
  pre.style <- NULL;
  post.style <- NULL;

  # style for no matches
  no.prestyle <- "<strong style=\"background-color:yellow; font-size=125%; color=\"black\">";
  no.poststyle <- "</strong>";
    
  hit.inx <- name.map$hit.inx;
  hit.values <- name.map$hit.values;
  match.state <- name.map$match.state;

  # construct the result table with cells wrapped in html tags
  # the unmatched will be highlighted in different background
  html.res <- matrix("", nrow=length(qvec), ncol=6);
  colnames(html.res) <- c("Query", "Match", "Species", "Genus", "NCBI_Taxonomy_ID", "GOLDSTAMP_ID");

  for (i in 1:length(qvec)){
    if(match.state[i]==1){
      pre.style <- "";
      post.style = "";
    }else{ # no matches
      pre.style <- no.prestyle;
      post.style <- no.poststyle;
    }
           
    hit <- species.db[hit.inx[i], ,drop=F];

    html.res[i, ] <- c(paste(pre.style, qvec[i], post.style, sep=""),
                            paste(ifelse(match.state[i]==0, "", hit.values[i]), sep=""),
                            paste(ifelse(match.state[i]==0 || is.na(hit$species) ||is.null(hit$species) || hit$species=="" || hit$species=="NA","-",hit$species),sep=""),
                            paste(ifelse(match.state[i]==0 || is.na(hit$genus) ||is.null(hit$genus) || hit$genus=="" || hit$genus=="NA","-",hit$genus),  sep=""), 
                            paste(ifelse(match.state[i]==0 || is.na(hit$NCBITAX) ||is.null(hit$NCBITAX) || hit$NCBITAX=="" || hit$NCBITAX=="NA","-", paste("<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=", hit$NCBITAX," target='_blank'>", hit$NCBITAX,"</a>", sep="")),  sep=""),
                            paste(ifelse(match.state[i]==0 || is.na(hit$GOLDMAPID) ||is.null(hit$GOLDMAPID) || hit$GOLDMAPID=="" || hit$GOLDMAPID=="NA", "-", paste("<a href=https://gold.jgi.doe.gov/project?id=", hit$GOLDMAPID," target='_blank'>", hit$GOLDMAPID,"</a>", sep="")), sep=""))
  }

  return(data.frame(html.res,check.names=FALSE));
}

#'Calculate enrichment score.
#'@description This function calculates the enrichment score for TSEA.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
CalculateHyperScore <- function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  nm.map <- GetFinalNameMap(mbSetObj);
  valid.inx <- !(is.na(nm.map$Strain)| duplicated(nm.map$Strain));
  ora.vec <- nm.map$Strain[valid.inx];

  q.size <- length(ora.vec);
  if(q.size==0) {
    print("Query taxa is missing!");
    return(0);
  }

  if(all(is.na(ora.vec))) {
    print("Query taxa are all NA!");
    return(0);
  }
    
  # total uniq cmpds in the current mset lib
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));
  set.size<-length(current.mset);
    
  if(set.size ==1){
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!");
    return(0);
  }

  hits <- lapply(current.mset, function(x){x[x %in% ora.vec]});
  hit.num <- unlist(lapply(hits, function(x) length(x)), use.names = FALSE);
    
  if(sum(hit.num>0)==0){

    AddErrMsg("No matches were found in the selected taxon set library!");

    if(grepl("_species$", mbSetObj$dataSet$tset.type)){
      AddErrMsg("Species-level taxa set was selected: verify that your list contains species names!");
    }else if(grepl("_strain$", mbSetObj$dataSet$tset.type)){
      AddErrMsg("Strain-level taxa set was selected: verify that your list contains strain names!");
    }else{
      AddErrMsg("Mixed-level taxa set was selected!");
    }
    return(0);
  }

  set.num<-unlist(lapply(current.mset, length), use.names = FALSE);

  # prepare for the result table
  res.mat<-matrix(NA, nrow=set.size, ncol=6);        
  rownames(res.mat)<-names(current.mset);
  colnames(res.mat)<-c("total", "expected", "hits", "Raw p", "Holm p", "FDR");

  res.mat[,1]<-set.num;
  res.mat[,2]<-q.size*(set.num/uniq.count);
  res.mat[,3]<-hit.num;
  res.mat[,4]<-phyper(hit.num-1, set.num, uniq.count-set.num, q.size, lower.tail=F);

  # Multiple-testing correction is applied only to sets with >= 2 hits. A single
  # matched taxon is not meaningful evidence of enrichment, and counting every
  # 1-hit set as a test inflates the penalty on the sets that matter (the
  # per-study libraries make this noticeable). 1-hit sets keep their raw p but
  # are reported with NA for Holm/FDR.
  test.inx <- hit.num >= 2;
  if(any(test.inx)){
    res.mat[test.inx, 5] <- p.adjust(res.mat[test.inx, 4], "holm");
    res.mat[test.inx, 6] <- p.adjust(res.mat[test.inx, 4], "fdr");
  }

  res.mat <- res.mat[hit.num>0,];

  # fix error when only 1 hit (not sig), no longer a matrix
  if(!("matrix" %in% class(res.mat))){
    AddErrMsg("No significant hits found using enrichment analysis!");
    return(0);
  }

  ord.inx<-order(res.mat[,4]);

  # download result
  mbSetObj$analSet$ora.mat = signif(res.mat[ord.inx,],3);
  mbSetObj$analSet$ora.hits = hits;
  fast.write(mbSetObj$analSet$ora.mat, file="tsea_ora_result.csv");

  # Safe-Handshake: Arrow save with verification
  tryCatch({
    ExportResultMatArrow(mbSetObj$analSet$ora.mat, "ora_mat");
  }, error = function(e) {
    warning(paste("Arrow save failed for ora_mat:", e$message));
  });

  return(.set.mbSetObj(mbSetObj));

}

#'Calculate GSEA-based enrichment score.
#'@description This function calculates the enrichment score using fgsea for TSEA.
#'@param mbSetObj Input the name of the mbSetObj.
#'@export
CalculateGseaScore <- function(mbSetObj){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  nm.map <- GetFinalNameMap(mbSetObj);
  valid.inx <- !(is.na(nm.map$Strain) | duplicated(nm.map$Strain));
  mapped.names <- nm.map$Strain[valid.inx];
  scores <- mbSetObj$dataSet$scores[valid.inx];

  # remove entries with NA scores
  keep <- !is.na(scores);
  mapped.names <- mapped.names[keep];
  scores <- scores[keep];

  if(length(mapped.names) == 0){
    AddErrMsg("No valid scored taxa after mapping!");
    return(0);
  }

  ranked.vec <- scores;
  names(ranked.vec) <- mapped.names;
  # resolve duplicates by keeping highest absolute score
  ranked.vec <- ranked.vec[order(abs(ranked.vec), decreasing = TRUE)];
  ranked.vec <- ranked.vec[!duplicated(names(ranked.vec))];
  ranked.vec <- sort(ranked.vec, decreasing = TRUE);

  set.size <- length(current.mset);
  if(set.size <= 1){
    AddErrMsg("Cannot perform enrichment analysis on a single taxon set!");
    return(0);
  }

  gsea.res <- run_func_via_microservice(
    func = function(pathways, stats){
      res <- fgsea::fgsea(pathways = pathways, stats = stats,
                          minSize = 3, maxSize = 500, scoreType = "std");
      as.data.frame(res);
    },
    args = list(pathways = current.mset, stats = ranked.vec),
    timeout_sec = 180
  );

  # filter to sets with leading edge
  le.counts <- sapply(gsea.res$leadingEdge, length);
  keep.inx <- le.counts > 0;
  if(sum(keep.inx) == 0){
    AddErrMsg("No enriched taxon sets found!");
    return(0);
  }
  gsea.res <- gsea.res[keep.inx, ];
  le.counts <- le.counts[keep.inx];

  # build result matrix matching ORA column layout
  res.mat <- matrix(NA, nrow = nrow(gsea.res), ncol = 6);
  rownames(res.mat) <- gsea.res$pathway;
  colnames(res.mat) <- c("total", "expected", "hits", "Raw p", "Holm p", "FDR");

  res.mat[, 1] <- gsea.res$size;        # total = set_size
  res.mat[, 2] <- gsea.res$NES;         # expected = NES
  res.mat[, 3] <- le.counts;            # hits = leading_edge_count
  res.mat[, 4] <- gsea.res$pval;        # Raw p
  res.mat[, 5] <- p.adjust(gsea.res$pval, "holm");
  res.mat[, 6] <- gsea.res$padj;        # FDR (fgsea's BH-adjusted)

  ord.inx <- order(res.mat[, 4]);
  res.mat <- signif(res.mat[ord.inx, , drop = FALSE], 3);

  # build hits list (leading edge taxa per set)
  hits <- lapply(seq_len(nrow(gsea.res)), function(i){
    gsea.res$leadingEdge[[i]];
  });
  names(hits) <- gsea.res$pathway;
  hits <- hits[rownames(res.mat)];

  mbSetObj$analSet$ora.mat <- res.mat;
  mbSetObj$analSet$ora.hits <- hits;
  fast.write(res.mat, file = "tsea_gsea_result.csv");

  tryCatch({
    ExportResultMatArrow(res.mat, "ora_mat");
  }, error = function(e){
    warning(paste("Arrow save failed for ora_mat:", e$message));
  });

  return(.set.mbSetObj(mbSetObj));
}

#'Getter to return final map
#'@description This function returns the final (after user selection) map as a dataframe.
#'Consists of two columns, original name and strain.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
GetFinalNameMap<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  enrtype <- mbSetObj$dataSet$q.type;
  qvec <- mbSetObj$dataSet$species;
  nm.mat <- matrix(nrow=length(qvec), ncol=2);
  colnames(nm.mat) <- c("query", "Strain");
    
  if(enrtype=="taxa"){
    for (i in 1:length(qvec)){
      nm.mat[i, ]<-c(qvec[i],qvec[i]);
    }
  }else{
    hit.inx <- name.map$hit.inx;
    hit.values <- name.map$hit.values;
    match.state <- name.map$match.state;
    species.db <- .read.microbiomeanalyst.lib.rds("microbe_db_new.rds", "tsea");
        
    # enrichment is done on canonical names: the mapped table entry's current NCBI name when the
    # query matched, otherwise the tidied query itself (it may still equal a set member verbatim)
    can <- .canonical.taxon.name(qvec);
    for (i in 1:length(qvec)){
      if(match.state[i]==1 && !is.na(hit.inx[i])){
        hit <- species.db[hit.inx[i], , drop=FALSE];
        nm <- if(!is.null(hit$canonical) && !is.na(hit$canonical) && nchar(hit$canonical) > 0) hit$canonical else hit$taxa;
      }else{
        nm <- can[i];
      }
      nm.mat[i, ] <- c(qvec[i], nm);
    }
  }
  return(as.data.frame(nm.mat,check.names=FALSE));
}

#'Function to prepare data for enrichment network.
#'@description This function prepares data for enrichment network.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PrepareEnrichNet<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
    
  #calculate the enrichment fold change
  folds <- mbSetObj$analSet$ora.mat[,3]/mbSetObj$analSet$ora.mat[,2];
  names(folds) <- GetShortNames(rownames(mbSetObj$analSet$ora.mat));
  hits <- mbSetObj$analSet$ora.mat[,3];
  pvals <- mbSetObj$analSet$ora.mat[,4];
  PlotEnrichNet.Overview(hits, pvals);
}

#'Set the microbe set library
#'@description This function sets the microbe
# ---- taxon name canonicalisation -------------------------------------------------------------
# The libraries mix nomenclature eras (Firmicutes / Bacillota, Ruminococcus gnavus /
# Mediterraneibacter gnavus ...). Both query taxa and set members are mapped to the current
# NCBI scientific name through tsea_name_synonyms.rds (built from NCBI names.dmp for every
# taxid in microbe_db_new.rds) so that either spelling reaches the same sets. Names NCBI does
# not know (most strain designations) are left as they are.
.tsea.cache <- new.env(parent = emptyenv())    # an environment, so it also works inside a locked package namespace
.get.tsea.synonyms <- function(){
  if(is.null(.tsea.cache$synonyms)){
    syn <- .read.microbiomeanalyst.lib.rds("tsea_name_synonyms.rds", "tsea");
    .tsea.cache$synonyms <- setNames(syn$canonical, syn$name);
  }
  .tsea.cache$synonyms
}

# tidy a user-supplied or library taxon name: rank prefixes (s__Genus_species), underscores,
# surrounding quotes/whitespace, repeated spaces
.tidy.taxon.name <- function(x){
  x <- sub("^[a-z]__", "", trimws(as.character(x)));
  x <- gsub("_", " ", x, fixed=TRUE);
  x <- gsub("^[\"']+|[\"']+$", "", x);
  gsub("[[:space:]]+", " ", x)
}

.canonical.taxon.name <- function(x){
  x <- .tidy.taxon.name(x);
  syn <- .get.tsea.synonyms();
  can <- syn[tolower(x)];
  ifelse(is.na(can), x, can)
}

# Taxon-set library file for each library key used by the web interface.
.tsea.lib.files <- c(
  # mixed level
  host_int          = "tsea_host_int.csv",
  host_ext          = "tsea_host_ext.csv",
  host_diet         = "tsea_host_diet_lifestyle.csv",
  host_drug         = "tsea_host_medication.csv",
  env               = "tsea_environment.csv",
  mic_met           = "taxon_metabolite_tsea.csv",
  mic_int           = "tsea_microbiome_int.csv",
  gene              = "tsea_host_snps_new.csv",
  food_matrix       = "tsea_food_matrix.csv",
  # species level
  host_int_species  = "tsea_host_int_species.csv",
  host_ext_species  = "tsea_host_ext_species.csv",
  host_diet_species = "tsea_host_diet_lifestyle_species.csv",
  host_drug_species = "tsea_host_medication_species.csv",
  env_species       = "tsea_environment_species.csv",
  food_matrix_species = "tsea_food_matrix_species.csv",
  # strain level
  host_int_strain   = "tsea_host_int_strain.csv",
  env_strain        = "tsea_environment_strain.csv",
  mic_int_strain    = "tsea_microbiome_int_strain.csv",
  food_matrix_strain = "tsea_food_matrix_strain.csv"
);

# Resolve a taxon-set library file to a local path. On the web server the file is
# read from the resources tree. In package mode it is downloaded once into the
# working directory and re-used for 30 days (same policy as
# .read.microbiomeanalyst.lib.rds), instead of being fetched on every call.
.get.tsea.lib.path <- function(filenm){
  if(.on.public.web){
    return(paste0(rpath, "libs/tsea/", filenm));
  }
  stale <- !file.exists(filenm) ||
           difftime(Sys.time(), file.info(filenm)$mtime, units="days") > 30;
  if(stale){
    lib.url <- paste0("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/libs/tsea/", filenm);
    # download to a temp file first so an interrupted transfer never masquerades as a cached library
    tmp <- paste0(filenm, ".part");
    ok <- tryCatch({ download.file(lib.url, destfile=tmp, method="curl", quiet=TRUE); TRUE },
                   error=function(e){
                     tryCatch({ download.file(lib.url, destfile=tmp, method="libcurl", quiet=TRUE); TRUE },
                              error=function(e2) FALSE)
                   });
    if(ok && file.exists(tmp) && file.info(tmp)$size > 0){
      file.rename(tmp, filenm);
    }else{
      if(file.exists(tmp)) unlink(tmp);
      if(!file.exists(filenm)){
        AddErrMsg(paste("Could not download taxon set library", filenm, "- check your internet connection."));
        return(NULL);
      }
      # else: keep using the stale copy
    }
  }
  filenm
}

#'set library for TSEA.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
SetTaxonSetLib <- function(mbSetObj, tset.type){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  mbSetObj$dataSet$tset.type <- tset.type
  
  filenm <- .tsea.lib.files[tset.type];
  if(is.na(filenm)){
    AddErrMsg(paste("Unknown taxon set library:", tset.type));
    return(0);
  }
  libPath <- .get.tsea.lib.path(filenm);
  if(is.null(libPath)){
    return(0);
  }

  current.msetlib <<- .readDataTable(libPath);
  ms.list <- strsplit(current.msetlib[,2],"; ");
  names(ms.list) <- current.msetlib[,1];
  # Some taxon-set libraries repeat a set name across rows. Consolidate duplicates
  # here, at load, by unioning their members into one entry. Leaving duplicates in
  # place (a) double-counts the same set in the hypergeometric test and (b) makes a
  # later as.data.frame() of the ORA result matrix call make.names(unique = TRUE),
  # which mangles every set name (e.g. "Crohn Disease (decrease)" ->
  # "Crohn.Disease..decrease."). Keeping set names unique at the source avoids both.
  if (anyDuplicated(names(ms.list))) {
    nms    <- names(ms.list);
    keep   <- !duplicated(nms);
    merged <- lapply(nms[keep], function(nm) unique(unlist(ms.list[nms == nm], use.names = FALSE)));
    names(merged) <- nms[keep];
    lib2 <- current.msetlib[keep, , drop = FALSE];
    lib2[, 2] <- vapply(merged, paste, character(1), collapse = "; ");
    current.msetlib <<- lib2;
    ms.list <- merged;
  }
  # members under old nomenclature join their current-name equivalents (see .canonical.taxon.name)
  set.nms <- names(ms.list);
  all.can <- .canonical.taxon.name(unlist(ms.list, use.names = FALSE));
  grp <- factor(rep(seq_along(ms.list), lengths(ms.list)), levels = seq_along(ms.list));   # keeps empty sets in place
  ms.list <- lapply(split(all.can, grp), unique);
  names(ms.list) <- set.nms;
  current.mset <<- ms.list;
  # total uniq cmpds in the mset lib
  uniq.count <<- length(unique(unlist(current.mset, use.names = FALSE)));
  return(.set.mbSetObj(mbSetObj))
}

#'Create network for enrichmnet overview
#'@description This function creates the plot
#'for the enrichent network overview.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import igraph
#'@import reshape
PlotEnrichNet.Overview<-function(hits, pvals){
  
  suppressMessages(library(igraph));

  # due to space limitation, plot top 50 if more than 50 were given
  title <- "Taxon Set Enrichment Network Overview";
  if(length(hits) > 50){
    hits <- hits[1:50];
    pvals <- pvals[1:50];
    title <- "Enrichment Overview (top 50)";
  }
  
  pvalue <- pvals;
  id <- names(pvalue);
  geneSets <- current.mset;
  n <- length(pvalue);
  w <- matrix(NA, nrow=n, ncol=n);
  colnames(w) <- rownames(w) <- id;

  # OPTIMIZED: Vectorized computation of overlap matrix
  # Instead of nested loops, compute upper triangle indices and vectorize
  upper_tri_indices <- which(upper.tri(w, diag = TRUE), arr.ind = TRUE);

  # Vectorized computation using mapply
  overlap_values <- mapply(
    function(i_idx, j_idx) {
      overlap_ratio(geneSets[id[i_idx]], geneSets[id[j_idx]])
    },
    upper_tri_indices[, 1],
    upper_tri_indices[, 2],
    SIMPLIFY = TRUE
  );

  # Fill upper triangle with computed values
  w[upper_tri_indices] <- overlap_values;

  # Mirror to lower triangle (overlap is symmetric)
  w[lower.tri(w)] <- t(w)[lower.tri(w)];

  wd <- reshape2::melt(w);
  wd <- wd[wd[,1] != wd[,2],];
  wd <- wd[!is.na(wd[,3]),];
  g <- graph_from_data_frame(wd[,-3], directed=F);
  E(g)$width <- sqrt(wd[,3]*20);
  g <- delete_edges(g, E(g)[wd[,3] < 0.2]);
  idx <- unlist(sapply(V(g)$name, function(x) match(x,id)));
  pvalue <- pvalue[idx]
  cols <- color_scale("red", "#E5C494");
  V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))];

  cnt <- hits + 2;
  names(cnt) <- id;
  cnt2 <- cnt[V(g)$name];
  V(g)$size <- log(cnt2, base=10) * 10; ## cnt2/sum(cnt2) * 100;
  #V(g)$size <- cnt2/sum(cnt2) * 10;
    
  # layout
  pos.xy <- layout_with_fr(g);

  # now create the json object
  nodes <- vector(mode="list");
  node.nms <- V(g)$name;
  node.sizes <- V(g)$size;
  node.cols <- V(g)$color;
    
  for(i in 1:length(node.sizes)){
    nodes[[i]] <- list(id = node.nms[i],
                  label=node.nms[i], 
                  size=node.sizes[i], 
                  color=node.cols[i],
                  x = pos.xy[i,1],
                  y = pos.xy[i,2]);
  }
    
  edge.mat <- as_edgelist(g);
  edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);

  # covert to json
  netData <- list(nodes=nodes, edges=edge.mat);
  sink("tsea_network.json");
  cat(RJSONIO::toJSON(netData));
  sink();
}

#################################
########## Utility Fx ###########
#################################

# Getter for ORA matrix
GetORA.rowNames<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  nms <- rownames(mbSetObj$analSet$ora.mat);
  
  if(is.null(nms)){
    return("NA");
  }
  return(nms);
}

# Getter
GetORA.mat<-function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  ora_mat <- mbSetObj$analSet$ora.mat;

  # Safe-Handshake: Arrow save with verification
  tryCatch({
    ExportResultMatArrow(ora_mat, "ora_mat");
  }, error = function(e) {
    warning(paste("Arrow save failed for ora_mat:", e$message));
  });

  return(ora_mat);
}

# Getter
GetORA.colorBar<-function(mbSetObj){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  len <- nrow(mbSetObj$analSet$ora.mat);
  
  if(len > 60){
    ht.col <- c(substr(heat.colors(50), 0, 7), rep("#FFFFFF", len-50));
  }else{
    # reduce to hex by remove the last character so HTML understand
    ht.col <- substr(heat.colors(len), 0, 7);
  }

  return (ht.col);
}

# methods to return the selected metset require to java for display
GetMsetNames<-function(){
  return(current.msetlib$name);
}

GetMsetMembers<-function(){
  return(current.msetlib$member);
}

GetMsetReferences<-function(){
  return(current.msetlib$reference);
}

GetShortNames<-function(nm.vec, max.len= 45){
  new.nms <- vector(mode="character", length=length(nm.vec));
  for(i in 1:length(nm.vec)){
    nm <- nm.vec[i];
    # NA-safe: nchar(NA) is NA_integer_ which makes `if (NA <= max.len)`
    # throw "missing value where TRUE/FALSE needed" — happens when the
    # taxon-set library has NA-rownamed entries that propagate into
    # mbSet$analSet$ora.mat. Substitute a placeholder and continue.
    if(is.na(nm) || is.null(nm)){
      new.nms[i] <- "(unnamed)";
      next;
    }
    if(nchar(nm) <= max.len){
      new.nms[i] <- nm;
    }else{
      wrds <- strsplit(nm, "[[:space:]]+")[[1]];
      new.nm <- "";
      if(length(wrds)>1){
        for(m in 1:length(wrds)){
          wrd <- wrds[m];
          if(nchar(new.nm)+4+nchar(wrd) <= max.len){
            new.nm <- paste(new.nm, wrd);
          }else{
            new.nms[i] <- paste (new.nm, "...", sep="");
            break;
          }
        }
      }else{
        new.nms[i] <- paste (substr(nm, 0, 21), "...", sep="");
      }
    }
  }
  return (new.nms);
}

overlap_ratio <- function(x, y) {
  x <- unlist(x)
  y <- unlist(y)
  length(intersect(x, y))/length(unique(c(x,y)))
}

#' @importFrom grDevices colorRampPalette
color_scale <- function(c1="grey", c2="red") {
  pal <- grDevices::colorRampPalette(c(c1, c2))
  colors <- pal(100)
  return(colors)
}

getIdx <- function(v, MIN, MAX) {
  if ( MIN == MAX ) {
    return(100)
  }
  intervals <- seq(MIN, MAX, length.out=100)
  max(which(intervals <= v))
}

GetCurrentImg <- function(){
  return (current.img);
}

# given a metset inx, return hmtl highlighted metset cmpds and references
GetHTMLMetSet<-function(mbSetObj, msetNm){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  hits <- mbSetObj$analSet$ora.hits;
  # highlighting with different colors
  mset <- current.mset[[msetNm]];
  red.inx <- which(mset %in% hits[[msetNm]]);
  mset[red.inx] <- paste("<font color=\"red\">", "<b>", mset[red.inx], "</b>", "</font>",sep="");

  grey.inx <- which(!(mset %in% current.mset[[msetNm]]));
  mset[grey.inx] <- paste("<font color=\"grey\">", "<b>", mset[grey.inx], "</b>", "</font>",sep="");

  # get references
  matched.inx <- match(tolower(msetNm), tolower(current.msetlib$name))[1];

  return(cbind(msetNm, paste(mset, collapse="; "), current.msetlib$reference[matched.inx]));
}

GetMsetPval<-function(mbSetObj, msetNm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  return(mbSetObj$analSet$ora.mat[msetNm, "Raw p"]);
}

GetMsetEvidence <- function(mbSetObj, msetNm){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  matched.inx <- match(tolower(msetNm), tolower(current.msetlib$name))[1];
  #print(current.msetlib$evidence[matched.inx])
  return(current.msetlib$evidence[matched.inx])
}

