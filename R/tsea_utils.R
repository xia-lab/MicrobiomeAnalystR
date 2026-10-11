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
  mbSetObj$dataSet$tsea.input <- NULL;   # a new list starts without a chosen input kind (SetTseaListInput)
  if(method == "gsea"){
    parts <- strsplit(lines, "\t");
    mbSetObj$dataSet$species <- trimws(sapply(parts, `[`, 1));
    mbSetObj$dataSet$scores <- as.numeric(sapply(parts, `[`, 2));
    # optional third and fourth columns: the list's p-value and adjusted p-value (FDR); blank where it has none
    p.col <- suppressWarnings(as.numeric(sapply(parts, `[`, 3)));
    f.col <- suppressWarnings(as.numeric(sapply(parts, `[`, 4)));
    mbSetObj$dataSet$p.values <- if(any(is.finite(p.col))) p.col else NULL;
    mbSetObj$dataSet$fdr.values <- if(any(is.finite(f.col))) f.col else NULL;
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
  }else if(q.type %in% c("genus","species","strain","mixed")){   # one table holds the names of every rank, so a genus name is mapped as any name
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

#'Say what kind of list was entered, and where a complete list is cut.
#'@description A taxon list is either a set of significant taxa (the default: tested against the taxa
#'measured in the study, see \code{CalculateHyperScore}) or a COMPLETE list, every taxon measured with a
#'score. A complete list is its own background: \code{CalculateHyperScore} then tests the part of it chosen
#'by the cutoff against the whole list (set sizes are the list's members, sets with at least two of them
#'are tested, the adjustment runs over all of those), and \code{CalculateGseaScore} scores the whole ranking.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param input "significant" or "complete".
#'@param cutoff.by "top" (the top fraction of the list by the size of its score, or by its order when it has
#'no score), "score" (the taxa whose |score| passes the cutoff in the direction \code{cutoff.dir}; refused when the
#'list has no score), "fdr" or "p" (the taxa below a threshold of the list's own adjusted p / p column; refused when
#'the list has no such column).
#'@param cutoff the fraction for "top" (strictly between 0 and 1, default 0.10), the threshold for "score" / "fdr" / "p".
#'@param cutoff.dir for "score" only: ">=" (|score| at least the cutoff: a fold change or a statistic) or "<="
#'(|score| at most the cutoff: a p-value-like score). When not given it follows the kind of score the list carries
#'(\code{.tsea.score.kind}). The direction used is recorded in \code{analSet$ora.background$cutoff.dir}.
#'@export
SetTseaListInput <- function(mbSetObj = NA, input = "significant", cutoff.by = "top", cutoff = 0.10, cutoff.dir = NULL){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(!(length(input) == 1 && input %in% c("complete", "significant"))){
    AddErrMsg("The list is either 'complete' or 'significant'.");
    return(0);
  }
  if(!(length(cutoff.by) == 1 && cutoff.by %in% c("top", "score", "fdr", "p"))){
    AddErrMsg("The cutoff is a top fraction ('top'), a score ('score'), an adjusted p-value ('fdr') or a p-value ('p').");
    return(0);
  }
  if(!(is.numeric(cutoff) && length(cutoff) == 1 && is.finite(cutoff) && cutoff > 0 && (cutoff.by != "top" || cutoff < 1))){
    AddErrMsg(if(cutoff.by == "top") "The top fraction must be more than 0 and less than 1 (0.10 is the top 10%)." else "The threshold must be more than 0.");
    return(0);
  }
  if(!is.null(cutoff.dir) && !(length(cutoff.dir) == 1 && cutoff.dir %in% c(">=", "<="))){
    AddErrMsg("The direction of a score cutoff is '>=' (|score| at least the cutoff) or '<=' (|score| at most the cutoff).");
    return(0);
  }
  mbSetObj$dataSet$tsea.input <- input;
  mbSetObj$dataSet$tsea.cutoff.by <- cutoff.by;
  mbSetObj$dataSet$tsea.cutoff <- cutoff;
  mbSetObj$dataSet$tsea.cutoff.dir <- if(cutoff.by == "score") cutoff.dir else NULL;
  return(.set.mbSetObj(mbSetObj));
}

# What kind of score a complete list carries, by its values: "p" when every value lies in [0, 1] (a p-value-like score: small is
# strong), else "effect" (a fold change or a statistic: any negative value, or a value beyond [0, 1]; large in size is strong).
.tsea.score.kind <- function(sc){
  sc <- sc[is.finite(sc)];
  if(length(sc) == 0) return(NA_character_);
  if(all(sc >= 0 & sc <= 1)) "p" else "effect"
}

# The query and the background of a COMPLETE list: the whole mapped list is the background, the part chosen by the cutoff is the query.
# Returns list(query, background, kind, cutoff), or NULL with the reason recorded.
.tsea.complete.query <- function(mbSetObj, nm.map, valid.inx){
  all.vec <- nm.map$Strain[valid.inx];
  pick <- function(v) if(length(v) == length(nm.map$Strain)) v[valid.inx] else NULL;
  by <- mbSetObj$dataSet$tsea.cutoff.by; if(is.null(by)) by <- "top";
  cutoff <- mbSetObj$dataSet$tsea.cutoff; if(is.null(cutoff)) cutoff <- 0.10;
  if(by == "top"){
    sc <- pick(mbSetObj$dataSet$scores);
    n.q <- max(1, min(length(all.vec), floor(length(all.vec) * cutoff)));
    # the largest change in either direction first (as MetaboAnalyst reads a complete list); the list order when it has no score
    ord <- if(is.null(sc)) seq_along(all.vec) else order(-abs(ifelse(is.finite(sc), sc, NA)), na.last = TRUE);
    query <- all.vec[ord][seq_len(n.q)];
  } else if(by == "score"){
    sc <- pick(mbSetObj$dataSet$scores);
    if(is.null(sc) || !any(is.finite(sc))){
      AddErrMsg("This list has no score, so a score cutoff cannot be used; use the top fraction (for example the top 10%) instead.");
      return(NULL);
    }
    kind <- .tsea.score.kind(sc);
    # the direction is the user's when given; otherwise the kind of score decides (large is strong in a fold change or a statistic, small in a p-value)
    dir <- mbSetObj$dataSet$tsea.cutoff.dir; if(is.null(dir)) dir <- if(kind == "p") "<=" else ">=";
    pass <- is.finite(sc) & (if(dir == "<=") abs(sc) <= cutoff else abs(sc) >= cutoff);
    query <- all.vec[pass];
    if(length(query) == 0){
      AddErrMsg(paste0("No taxon of the list has |score| ", dir, " ", cutoff, "."));
      return(NULL);
    }
    return(list(query = unique(query), background = unique(all.vec), kind = "score", cutoff = cutoff, score.kind = kind, cutoff.dir = dir));
  } else {
    v <- pick(if(by == "fdr") mbSetObj$dataSet$fdr.values else mbSetObj$dataSet$p.values);
    col <- if(by == "fdr") "FDR" else "p-value";
    if(is.null(v) || all(is.na(v))){
      AddErrMsg(paste0("This list has no ", col, " column, so a ", col, " cutoff cannot be used; use the top fraction (for example the top 10%) instead."));
      return(NULL);
    }
    query <- all.vec[!is.na(v) & v < cutoff];
    if(length(query) == 0){
      AddErrMsg(paste0("No taxon of the list has a ", col, " below ", cutoff, "."));
      return(NULL);
    }
  }
  list(query = unique(query), background = unique(all.vec), kind = if(by == "top") "top_fraction" else by, cutoff = cutoff);
}

# ---- rolling a list up to the level of the library --------------------------------------------------------------------------------------
# A library holds the names of one level, so a list that names finer taxa matches almost none of it: a genus-level library has no place for a
# species name. A list can be tested against a library of its own level or a coarser one, never a finer one, so for such a library the list
# is rolled up: a strain name to its species for a species-level library, a species or a strain name to its genus for a genus-level one. One
# row stays per name, a name keeping the score of its strongest member (the largest size of score, its sign kept; the smallest value for a
# list of p values and for the p and FDR columns), and the taxa measured are rolled up the same way.

# The name a species or a strain name rolls up to at a level ("genus" or "species"): its first word, or its first two, without a bracket or a
# "Candidatus" in front ("[Eubacterium] siraeum", "Escherichia coli K-12"). A name of a rank that is not finer than the level is returned as it is.
.tsea.to.level <- function(x, level){
  x <- as.character(x);
  if(!(level %in% c("genus", "species"))) return(x);
  rk <- .tsea.taxon.rank(x);
  down <- !is.na(rk) & rk %in% (if(level == "genus") c("species", "strain") else "strain");
  if(any(down)){
    y <- .tidy.taxon.name(x[down]);
    y <- trimws(gsub("[\\[\\]]", "", sub("^Candidatus +", "", y, ignore.case = TRUE), perl = TRUE));
    keep <- if(level == "genus") 1 else 2;
    x[down] <- .canonical.taxon.name(vapply(strsplit(y, " ", fixed = TRUE), function(v) paste(utils::head(v, keep), collapse = " "), ""));
  }
  x
}

# The level of the library a list is tested against ("genus", "species" or "strain"; NA for a library of the mixed level): the library the user
# chose, steered by the rank of the list as CalculateHyperScore and CalculateGseaScore steer it (ResolveTaxonSetLib).
.tsea.library.level <- function(mbSetObj, nms){
  req <- mbSetObj$dataSet$tset.type;
  if(length(req) != 1 || is.na(req) || !(req %in% names(.tsea.lib.files))) return(NA_character_);
  nms <- nms[!is.na(nms) & nzchar(nms)];
  if(length(nms) == 0) return(NA_character_);
  used <- ResolveTaxonSetLib(req, TaxonRankSummary(nms)$rank);
  if(grepl("_(genus|species|strain)$", used)) sub("^.*_", "", used) else NA_character_
}

# The name map of a list, and the list's scores, rolled up to the level of the library; NULL when there is nothing to roll up.
# Returns list(nm.map, set, from, to, scored, level): the map with the rolled name as its Strain, a COPY of the data set holding the rolled
# scores (the list as typed is not changed, so the same list can be run against another library), the names before and after, whether it had
# scores, and the level. The value of a name sits on its first row; the other rows of that name are the duplicates the callers already drop.
.tsea.roll.up <- function(mbSetObj, nm.map){
  nms <- as.character(nm.map$Strain);
  level <- .tsea.library.level(mbSetObj, nms);
  if(is.na(level) || level == "strain") return(NULL);
  g <- .tsea.to.level(nms, level);
  if(!any(!is.na(nms) & !is.na(g) & g != nms)) return(NULL);
  strongest <- function(v, small = FALSE){
    v <- v[!is.na(v)];
    if(length(v) == 0) NA_real_ else if(small) min(v) else v[which.max(abs(v))]
  };
  roll <- function(v, small = FALSE){
    if(is.null(v) || length(v) != length(g)) return(v);
    for(grp in split(seq_along(g)[!is.na(g)], g[!is.na(g)])){
      if(length(grp) > 1) v[grp[1]] <- strongest(v[grp], small);
    }
    v
  };
  qset <- mbSetObj;
  sc <- mbSetObj$dataSet$scores;
  qset$dataSet$scores <- roll(sc, small = identical(.tsea.score.kind(sc), "p"));
  qset$dataSet$p.values <- roll(mbSetObj$dataSet$p.values, small = TRUE);
  qset$dataSet$fdr.values <- roll(mbSetObj$dataSet$fdr.values, small = TRUE);
  nm2 <- nm.map;
  nm2$Strain <- g;
  list(nm.map = nm2, set = qset, from = length(unique(nms[!is.na(nms)])), to = length(unique(g[!is.na(g)])), scored = !is.null(sc), level = level);
}

#'The level a list of taxon names is read at, and the libraries it can be tested against.
#'@description A list can be tested against the libraries of the finest level it names or of a coarser one (a strain name rolls up to its
#'species and to its genus, a genus name cannot roll down). The finest level is the rank of the list when its names are of one rank
#'(\code{TaxonRankSummary}), and the finest rank present when they are of several.
#'@param mbSetObj Input the name of the mbSetObj.
#'@return Three strings: the finest level ("genus", "species", "strain"; "higher" for a list of families and above; "mixed" when it cannot be
#'read), the rank of the list ("mixed" when its names are of several ranks), and a phrase with the number of names of each rank.
#'@export
GetTseaListLevel <- function(mbSetObj){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  rs <- TaxonRankSummary(mbSetObj$dataSet$species);
  nb <- rs$n.by.rank;
  if(is.na(rs$rank) || length(nb) == 0) return(c("mixed", "mixed", ""));
  levels.fine <- c("genus", "species", "strain");
  present <- intersect(levels.fine, names(nb));
  finest <- if(rs$rank %in% levels.fine) rs$rank else if(length(present) == 0) "higher" else present[length(present)];
  parts <- paste(nb, ifelse(names(nb) == "higher", "higher-rank", names(nb)));
  txt <- paste0(if(length(parts) > 1) paste(paste(parts[-length(parts)], collapse = ", "), "and", parts[length(parts)]) else parts, " names");
  c(finest, rs$rank, txt)
}

#'Calculate enrichment score.
#'@description This function calculates the enrichment score for TSEA.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param background Optional character vector of the taxa measured in the study. When given, chance is
#'defined by what could have been seen: the background is the measured taxa that are also in the library,
#'a set's size is its measured members, the query is the query taxa that are in that background, and
#'only sets with at least two measured members are tested (the multiple-testing correction runs over all
#'of them). When NULL (the default) or empty, the whole library is the background. A "source" attribute on
#'the vector (a short phrase naming where the list came from) is kept with the result, in
#'\code{analSet$ora.background}.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
CalculateHyperScore <- function(mbSetObj, background = NULL){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  nm.map <- GetFinalNameMap(mbSetObj);
  # a library holds names of one level: the list is rolled up to it, never down (qset: what the list is read from)
  roll <- .tsea.roll.up(mbSetObj, nm.map);
  qset <- mbSetObj;
  if(!is.null(roll)){ nm.map <- roll$nm.map; qset <- roll$set; }
  valid.inx <- !(is.na(nm.map$Strain)| duplicated(nm.map$Strain));
  ora.vec <- nm.map$Strain[valid.inx];

  q.size <- length(ora.vec);
  if(q.size==0) {
    print("Query taxa is missing!");
    return(0);
  }
  list.names <- ora.vec;      # the list as given: the query of a significant list, the whole of a complete one (its rank is stated, and decides the library variant)

  if(all(is.na(ora.vec))) {
    print("Query taxa are all NA!");
    return(0);
  }

  # A COMPLETE list (SetTseaListInput) is its own background: the query is the part of it the cutoff picks, the whole mapped list is
  # what was measured. Whatever background was passed in is not used then.
  complete.sel <- NULL;
  if(identical(mbSetObj$dataSet$tsea.input, "complete")){
    complete.sel <- .tsea.complete.query(qset, nm.map, valid.inx);
    if(is.null(complete.sel)) return(0);
    ora.vec <- complete.sel$query;
    q.size <- length(ora.vec);
    background <- complete.sel$background;
    attr(background, "source") <- paste0("all ", length(background), " taxa of your list");
  }
    
  # the list's rank decides which variant of the chosen library it is tested against (ResolveTaxonSetLib); the rank and the variant are recorded
  st <- .tsea.steer(mbSetObj, list.names, if(is.null(complete.sel)) background else NULL);
  mbSetObj <- st$set;

  # total uniq cmpds in the current mset lib
  uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));

  # With a list of the taxa measured in the study, chance is defined by what could have been seen: the
  # background is the measured taxa that are in the library, a set counts only its measured members, and the
  # query is the query taxa inside that background. Without one, the whole library is the background.
  q.input <- q.size;
  bg.source <- attr(background, "source");
  use.bg <- length(background) > 0;
  if(use.bg){
    background <- unique(as.character(background));
    if(!identical(mbSetObj$dataSet$q.type, "taxa")){
      background <- unique(.canonical.taxon.name(background));  # the same names the query and the sets are read in
    }
    bg.level <- .tsea.library.level(mbSetObj, list.names);
    if(!is.na(bg.level) && bg.level != "strain") background <- unique(.tsea.to.level(background, bg.level));   # measured species count as their genus
    current.mset <- lapply(current.mset, function(x){x[x %in% background]});
    uniq.count <- length(unique(unlist(current.mset, use.names = FALSE)));
    ora.vec <- ora.vec[ora.vec %in% unlist(current.mset, use.names = FALSE)];
    q.size <- length(ora.vec);
    if(q.size==0) {
      AddErrMsg("None of the query taxa was both measured in the study and present in the selected taxon set library!");
      return(0);
    }
  }

  # A set of fewer than two members is not tested (one taxon cannot be "enriched"). Filtered by size alone, before
  # any test, so the family the p-values are adjusted over never depends on the hits; the background is counted above.
  current.mset <- current.mset[lengths(current.mset) >= 2];
  set.size<-length(current.mset);

  if(use.bg && set.size==0){
    AddErrMsg("No taxon set has two or more of the measured taxa, so there is nothing to test!");
    return(0);
  }

  if(set.size ==1){
    AddErrMsg("Cannot perform enrichment analysis on a single taxon set!");
    return(0);
  }

  hits <- lapply(current.mset, function(x){x[x %in% ora.vec]});
  hit.num <- unlist(lapply(hits, function(x) length(x)), use.names = FALSE);
    
  if(sum(hit.num>0)==0){

    AddErrMsg("No matches were found in the selected taxon set library!");

    if(grepl("_genus$", mbSetObj$dataSet$tset.type)){
      AddErrMsg("Genus-level taxa set was selected: verify that your list contains genus names!");
    }else if(grepl("_species$", mbSetObj$dataSet$tset.type)){
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

  # Multiple-testing correction over EVERY set tested (a set with no hit has p = 1). Adjusting only over the sets
  # with >= 2 hits chose the family by the outcome, which understates the adjusted p-values; the sets too small to
  # be evidence are left out by their SIZE above, before any test.
  res.mat[,5] <- p.adjust(res.mat[,4], "holm");
  res.mat[,6] <- p.adjust(res.mat[,4], "fdr");

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
  # the members each reported set was tested with (its measured members): the set's view marks the others as not tested
  mbSetObj$analSet$ora.tested = current.mset[rownames(mbSetObj$analSet$ora.mat)];
  # what the test was run against, so every result surface can say it
  mbSetObj$analSet$ora.background = list(mode = ifelse(use.bg, "measured", "library"),
                                         background = uniq.count, query = q.size, query.input = q.input, sets.tested = set.size,
                                         source = if(use.bg && !is.null(bg.source)) as.character(bg.source)[1] else "");
  mbSetObj$analSet$ora.background[names(st$rec)] <- st$rec;     # the rank tested and the library variant used
  if(!is.null(roll)) mbSetObj$analSet$ora.background[c("rollup.from", "rollup.to", "rollup.scored", "rollup.level")] <- list(roll$from, roll$to, roll$scored, roll$level);
  # how big the tested sets are: a test over sets of two or three members has little power, and the result can say so
  mbSetObj$analSet$ora.background$set.size.median <- stats::median(set.num);
  mbSetObj$analSet$ora.background$sets.small <- sum(set.num < 5);
  if(!is.null(complete.sel)){
    # what was tested, for pages and manifests: "the top 8% (5) of your 62 taxa against all 62"
    mbSetObj$analSet$ora.background$kind <- complete.sel$kind;
    mbSetObj$analSet$ora.background$cutoff <- complete.sel$cutoff;
    mbSetObj$analSet$ora.background$list.size <- length(complete.sel$background);
    if(!is.null(complete.sel$score.kind)){
      mbSetObj$analSet$ora.background$score.kind <- complete.sel$score.kind;
      mbSetObj$analSet$ora.background$cutoff.dir <- complete.sel$cutoff.dir;
    }
  }
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
#'@description This function calculates the enrichment score using fgsea for TSEA: the complete ranked list
#'(every taxon measured, each with a score) is the input, no significance cutoff is applied and no background has
#'to be chosen, because the ranking itself is the universe. A set's size is the number of its members that are in
#'the ranking; sets with fewer than three are not tested. A ranking whose scores are all non-negative (a list
#'ordered by importance, or -log10 p-values) is scored one-sided, so that only the top of the list can be enriched.
#'The multiple-testing correction runs over every set tested. The permutations are seeded, so the same ranking
#'always gives the same table. What was tested is recorded in \code{analSet$ora.background}.
#'@param mbSetObj Input the name of the mbSetObj.
#'@export
CalculateGseaScore <- function(mbSetObj){

  mbSetObj <- .get.mbSetObj(mbSetObj);

  nm.map <- GetFinalNameMap(mbSetObj);
  # a library holds names of one level: the ranking is rolled up to it, never down (qset: what the ranking is read from)
  roll <- .tsea.roll.up(mbSetObj, nm.map);
  qset <- mbSetObj;
  if(!is.null(roll)){ nm.map <- roll$nm.map; qset <- roll$set; }
  valid.inx <- !(is.na(nm.map$Strain) | duplicated(nm.map$Strain));
  mapped.names <- nm.map$Strain[valid.inx];
  scores <- qset$dataSet$scores[valid.inx];

  # remove entries with NA scores
  keep <- !is.na(scores);
  mapped.names <- mapped.names[keep];
  scores <- scores[keep];
  n.scored <- length(scores);

  if(length(mapped.names) == 0){
    AddErrMsg("No valid scored taxa after mapping!");
    return(0);
  }

  # the ranking's rank decides which variant of the chosen library it is scored against (ResolveTaxonSetLib); the rank and the variant are recorded
  st <- .tsea.steer(mbSetObj, mapped.names);
  mbSetObj <- st$set;

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
      set.seed(1);  # the same ranking gives the same table, so a recorded run replays
      # all scores non-negative: the list is ordered from most to least important, only its top can be enriched
      res <- fgsea::fgsea(pathways = pathways, stats = stats,
                          minSize = 3, maxSize = 500, scoreType = ifelse(all(stats >= 0), "pos", "std"));
      as.data.frame(res);
    },
    args = list(pathways = current.mset, stats = ranked.vec),
    timeout_sec = 180
  );

  # adjusted over EVERY set tested, before the sets with no leading edge are left out of the table
  n.tested <- nrow(gsea.res);
  tested.sizes <- gsea.res$size;
  gsea.res$holm <- p.adjust(gsea.res$pval, "holm");

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
  res.mat[, 5] <- gsea.res$holm;
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
  # the members each reported set was tested with (those in the ranking): the set's view marks the others as not tested
  mbSetObj$analSet$ora.tested <- setNames(lapply(rownames(res.mat), function(nm) current.mset[[nm]][current.mset[[nm]] %in% names(ranked.vec)]), rownames(res.mat));
  # what the test was run over, so every result surface can say it: the ranking itself, no background
  mbSetObj$analSet$ora.background <- list(mode = "ranked", background = length(ranked.vec), query = length(ranked.vec),
                                          query.input = n.scored, sets.tested = n.tested, source = "the ranked list itself",
                                          set.size.median = stats::median(tested.sizes), sets.small = sum(tested.sizes < 5));
  mbSetObj$analSet$ora.background[names(st$rec)] <- st$rec;     # the rank tested and the library variant used
  if(!is.null(roll)) mbSetObj$analSet$ora.background[c("rollup.from", "rollup.to", "rollup.scored", "rollup.level")] <- list(roll$from, roll$to, roll$scored, roll$level);
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
  # named by set: on a one-row matrix m[,3] is a bare number, and the overview then had no set ids
  m <- mbSetObj$analSet$ora.mat;
  hits <- setNames(m[,3], rownames(m));
  pvals <- setNames(m[,4], rownames(m));
  PlotEnrichNet.Overview(hits, pvals);
}

# What belongs to one result of a taxon set test: the table, the hits, the members each set was tested with, and what it was tested against.
.tsea.kept.fields <- c("ora.mat", "ora.hits", "ora.tested", "ora.background");

#'Keep the result in state under a name.
#'@description A complete list is tested two ways and the page shows one result at a time, so the result in state is kept under a name
#'before the second test overwrites it, and put back with \code{ShowTseaResult} when its tab is chosen.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param key The name of the result, for example "top" or "rank".
#'@export
KeepTseaResult <- function(mbSetObj, key){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  keep <- mbSetObj$analSet$tsea.keep;
  if(is.null(keep)) keep <- list();
  keep[[key]] <- mbSetObj$analSet[intersect(names(mbSetObj$analSet), .tsea.kept.fields)];
  mbSetObj$analSet$tsea.keep <- keep;
  .set.mbSetObj(mbSetObj)
}

#'Make a kept result the one in state.
#'@description The details of a taxon set (its hits, its p value, the members it was tested with) and the enrichment network are read from the
#'result in state. A part the kept result does not have is removed from state, not left from another result.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param key The name the result was kept under.
#'@export
ShowTseaResult <- function(mbSetObj, key){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  k <- mbSetObj$analSet$tsea.keep[[key]];
  if(is.null(k)) stop("no such result");
  for(n in .tsea.kept.fields) mbSetObj$analSet[[n]] <- k[[n]];
  .set.mbSetObj(mbSetObj)
}

#'Prepare the data of the enrichment network viewer from a taxon set result.
#'@description The enrichment network viewer (the one MetaboAnalyst opens for its enrichment results) draws two networks from one JSON
#'and switches between them with its own View option: the sets only (nodes, edges: two sets are joined when they share at least 20% of
#'their taxa) and the sets with the taxa in them (bnodes, bedges: a line from a set to each of its members, the library's members the set
#'network joins its sets by, so two sets that share taxa are joined there too; a taxon node says whether it hit and whether it was tested). It is written from
#'the result in state (\code{analSet$ora.mat} and \code{ora.hits}: the over-representation result and the rank-based one alike, the
#'leading-edge taxa being the hits of the latter), for the top 50 sets as the page's own network. The JSON has the shape
#'MetaboAnalyst's builder writes. The same call writes the result table the viewer offers to download.
#'@param mbSetObj Input the name of the mbSetObj.
#'@param netNm The name of the JSON file without its extension; the result table is named from it (enrichNet_ -> enrichment_).
#'@return 1 when the network was written, 0 when it could not be (the reason is in the error messages).
#'@export
PrepareEnrichNetViewer <- function(mbSetObj, netNm = "enrichNet_tsea"){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  enr.mat <- mbSetObj$analSet$ora.mat;
  set.hits <- mbSetObj$analSet$ora.hits;
  set.tested <- mbSetObj$analSet$ora.tested;      # the members each set was tested with; NULL for a result saved before they were recorded
  if(is.null(enr.mat) || nrow(enr.mat) == 0){
    AddErrMsg("No taxon set enrichment result to draw a network from!");
    return(0);
  }
  if(nrow(enr.mat) > 50){
    enr.mat <- enr.mat[1:50, , drop=FALSE];      # the same top 50 as the page's own network
  }
  if(nrow(enr.mat) < 2){
    AddErrMsg("Fewer than two taxon sets to display — aborting network build.");
    return(0);
  }

  set.nms <- rownames(enr.mat);
  n <- length(set.nms);
  pvals <- setNames(enr.mat[, "Raw p"], set.nms);
  hit.num <- setNames(enr.mat[, "hits"], set.nms);
  hit.lists <- lapply(set.nms, function(s){
    h <- as.character(set.hits[[s]]);
    unique(h[!is.na(h) & nzchar(h)]);
  });
  names(hit.lists) <- set.nms;

  # the set network: two sets are joined when they share at least 20% of the taxa of the library sets
  members <- lapply(set.nms, function(s){
    m <- current.mset[[s]];
    if(is.null(m)) hit.lists[[s]] else m
  });
  pairs <- utils::combn(n, 2);
  ovl <- mapply(function(i, j) overlap_ratio(members[[i]], members[[j]]), pairs[1, ], pairs[2, ]);
  keep <- !is.na(ovl) & ovl >= 0.2;
  e.df <- data.frame(from = set.nms[pairs[1, keep]], to = set.nms[pairs[2, keep]], weight = ovl[keep], stringsAsFactors = FALSE);
  g <- igraph::graph_from_data_frame(e.df, directed = FALSE, vertices = data.frame(name = set.nms, stringsAsFactors = FALSE));

  # colour by the p value, size by the number of hits
  set.cols <- .enrnet.gradient(-log10(pmax(pvals, 1e-300) + min(pvals) / 2));
  set.sizes <- .enrnet.rescale(log10(pmax(hit.num, 0) + 1), 8, 20, 12);
  pos.xy <- .enrnet.layout(g, "Default");

  nodes <- lapply(seq_len(n), function(i){
    list(id = set.nms[i], pwId = set.nms[i], label = set.nms[i], size = set.sizes[i], true_size = set.sizes[i],
         molType = "pathway", colorb = set.cols[i], colorw = set.cols[i], pval = unname(pvals[i]),
         hits = unname(hit.num[i]), posx = pos.xy[i, 1], posy = pos.xy[i, 2]);
  });
  edges <- lapply(seq_len(nrow(e.df)), function(k){
    list(id = k, source = e.df$from[k], target = e.df$to[k], weight = e.df$weight[k], width = 0.5);
  });

  # the bipartite network: the sets and the taxa in them (a taxon named like a set is left out: node ids are unique). A set is drawn with
  # every member the library holds for it, as its own view lists them; the taxa the result calls hits are marked, and so are the members
  # that were not in the list tested (their sets were not tested on them)
  set.members <- lapply(seq_len(n), function(i){
    m <- c(members[[i]], hit.lists[[i]]);
    unique(m[!is.na(m) & nzchar(m)])
  });
  taxa <- setdiff(unique(unlist(set.members, use.names = FALSE)), set.nms);
  b.from <- unlist(lapply(seq_len(n), function(i) rep(set.nms[i], length(intersect(set.members[[i]], taxa)))), use.names = FALSE);
  b.to <- unlist(lapply(seq_len(n), function(i) intersect(set.members[[i]], taxa)), use.names = FALSE);
  if(is.null(b.from)){ b.from <- character(0); b.to <- character(0); }
  bg <- igraph::graph_from_data_frame(data.frame(from = b.from, to = b.to, stringsAsFactors = FALSE), directed = FALSE,
                                      vertices = data.frame(name = c(set.nms, taxa), stringsAsFactors = FALSE));
  deg <- as.numeric(igraph::degree(bg));
  is.set <- seq_along(deg) <= n;
  b.nms <- igraph::V(bg)$name;
  hit.taxa <- unique(unlist(hit.lists, use.names = FALSE));
  tested.taxa <- if(is.null(set.tested)) NULL else unique(unlist(set.tested[set.nms], use.names = FALSE));
  b.hit <- !is.set & b.nms %in% hit.taxa;
  b.tested <- is.set | is.null(tested.taxa) | b.nms %in% tested.taxa;
  # the sets keep the colour they have in the set network; the taxa: a hit strong, a tested member light, a member not tested grey
  b.cols <- c(set.cols, ifelse(b.hit[!is.set], "#1B9E9E", ifelse(b.tested[!is.set], "#6CD0D0", "#C8C8C8")));
  b.sizes <- .enrnet.rescale(log10(deg + 1), 6, 20, 13);
  b.sizes[is.set] <- set.sizes;
  b.pos <- .enrnet.layout(bg, "Default");
  bnodes <- lapply(seq_along(b.nms), function(i){
    list(id = b.nms[i], label = b.nms[i], size = b.sizes[i], colorb = b.cols[i], colorw = b.cols[i], true_size = b.sizes[i],
         molType = if(is.set[i]) "pathway" else "compound", isHit = as.integer(b.hit[i]), tested = as.integer(b.tested[i]),
         posx = b.pos[i, 1], posy = b.pos[i, 2]);
  });
  bedges <- lapply(seq_along(b.from), function(k){
    list(id = paste0("b", k), source = b.to[k], target = b.from[k]);
  });

  # the table of the viewer and the file it offers to download
  csv.mat <- enr.mat;
  colnames(enr.mat) <- gsub("[ .]", "_", colnames(enr.mat));
  colnames(enr.mat)[colnames(enr.mat) == "hits"] <- "Hits";
  fast.write(enr.mat, file = paste0(gsub("enrichNet", "enrichment", netNm), ".csv"), row.names = TRUE);
  enr <- apply(enr.mat, 1, as.list);

  # RJSONIO writes the names of a list as keys without escaping them: a set named with a quote or a backslash would break the JSON
  hits.json <- hit.lists;
  names(hits.json) <- gsub("([\"\\\\])", "\\\\\\1", names(hits.json));
  netData <- list(nodes = nodes, edges = edges, bnodes = bnodes, bedges = bedges, enr = unname(enr), id = names(enr),
                  sigCmpds = "", sizes = NULL, hits = hits.json, genelist = NULL, analType = "tsea", pwType = "taxon",
                  org = "NA", backgroundColor = list("#f5f5f5", "#0066CC"), naviString = "EnrichNetwork");
  writeLines(enc2utf8(RJSONIO::toJSON(netData)), paste0(netNm, ".json"), useBytes = TRUE);

  # kept for the viewer's Layout option (UpdateNetworkLayout)
  pheno.comps <- list();
  pheno.comps[["enrichNet"]] <- g;
  pheno.comps[["enrichNet_bipartite"]] <- bg;
  pheno.comps <<- pheno.comps;

  return(1);
}

#'Positions for the enrichment network on screen, by the layout the viewer's Layout option names.
#'@description The viewer asks for a new layout of the network it shows (the sets only or the sets with their taxa) and draws the
#'positions it gets back. The networks are those PrepareEnrichNetViewer kept.
#'@param algo The layout: Default, FrR, circle, random, lgl or gopt.
#'@param filenm The JSON file to write: the id and the position of every node.
#'@param curr.nm The network: enrichNet or enrichNet_bipartite.
#'@return The name of the file written.
#'@export
UpdateNetworkLayout <- function(algo, filenm, curr.nm = "enrichNet"){

  g <- pheno.comps[[curr.nm]];
  pos.xy <- .enrnet.layout(g, algo);
  nms <- igraph::V(g)$name;
  nodes <- lapply(seq_along(nms), function(i) list(id = nms[i], x = pos.xy[i, 1], y = pos.xy[i, 2]));
  writeLines(enc2utf8(RJSONIO::toJSON(list(nodes = nodes))), filenm, useBytes = TRUE);
  return(filenm);
}

# orange to fire-brick, the viewer's colours for a p value; values that are all the same get the middle colour (a gradient needs a range)
.enrnet.gradient <- function(x){
  col <- grDevices::colorRampPalette(c("#FFA500", "#FF4500", "#B22222"))(100);
  if(length(unique(x)) < 2) return(rep(col[50], length(x)));
  col[as.integer(cut(x, breaks = seq(min(x), max(x), length.out = 101), include.lowest = TRUE))]
}

# x rescaled to [from, to]; all the same: `same`
.enrnet.rescale <- function(x, from, to, same){
  x <- unname(x);
  if(length(x) == 0 || max(x) == min(x)) return(rep(same, length(x)));
  (x - min(x)) / (max(x) - min(x)) * (to - from) + from
}

# the same picture for the same network: a layout that draws on the random numbers is run under a fixed seed, and the caller's random state is kept
.enrnet.layout <- function(g, algo = "Default"){
  vc <- igraph::vcount(g);
  if(igraph::ecount(g) == 0 || vc < 3 || algo == "circle") return(igraph::layout_in_circle(g));
  rs <- if(exists(".Random.seed", envir = globalenv(), inherits = FALSE)) get(".Random.seed", envir = globalenv()) else NULL;
  on.exit(if(is.null(rs)) suppressWarnings(rm(".Random.seed", envir = globalenv())) else assign(".Random.seed", rs, envir = globalenv()));
  set.seed(8574);
  if(algo == "FrR"){
    igraph::layout_with_fr(g)
  }else if(algo == "random"){
    igraph::layout_randomly(g)
  }else if(algo == "lgl"){
    igraph::layout_with_lgl(g)
  }else if(algo == "gopt"){
    igraph::layout_with_graphopt(g)
  }else if(vc > 5000){
    igraph::layout_with_lgl(g)
  }else if(vc < 100){
    igraph::layout_with_kk(g)
  }else{
    igraph::layout_with_fr(g)
  }
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
.tsea.cache$missing <- list()                   # library files asked for quietly and not there
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

# ---- taxon ranks ---------------------------------------------------------------------------------
# Most taxon-set libraries MIX ranks in one set (genus, species and strain names side by side); their _species / _strain variants are
# single-rank. The rank a list is written at decides which variant it is tested against, and the result states the rank it tested.
# The rank is read from the NAME (the name-mapping database has no rank column to trust: its 'strain' column repeats species names):
# one word is a genus; two words a species; more a strain ("Escherichia coli K-12", "Streptococcus phocae subsp. salmonis"). A
# placeholder species is a species however long its tail is ("Eubacterium sp CAG 180", "Clostridium sp. CAG:127", "Firmicutes bacterium
# CAG:110"), a genus label of the marker-gene databases is a genus however many words it has ("Clostridium sensu stricto 1",
# "Ruminococcaceae UCG-014", "Lachnospiraceae NK4A136 group", "Prevotella 9", "Pseudoramibacter Eubacterium": a second word that is not a
# lower-case epithet), and a single word that ends as a family, order, class or phylum name does
# is "higher" (above genus: "Lachnospiraceae", "Bacteroidota", "Firmicutes"). "Candidatus" and square brackets are not words of the name.
.tsea.rank.levels <- c("genus", "species", "strain", "higher")     # the order everywhere: Genus, Species, Strain (then above genus; mixed is last)
.tsea.variant.ranks <- c("genus", "species", "strain")                # the ranks a library can have a single-rank variant of (suffix _genus / _species / _strain)
.tsea.rank.pure <- 0.9;     # a list is "at" a rank when at least this share of its names are written at it; otherwise it is mixed

.tsea.taxon.rank <- function(x){
  x <- .tidy.taxon.name(x);
  x <- trimws(gsub("[\\[\\]]", "", sub("^Candidatus +", "", x, ignore.case = TRUE), perl = TRUE));
  w <- lengths(strsplit(x, " ", fixed = TRUE));
  rk <- ifelse(w <= 1, "genus", ifelse(w == 2, "species", "strain"));
  rk[w >= 3 & grepl("^[^ ]+ (sp\\.?|bacterium|cf\\.?|aff\\.?)( |$)", x)] <- "species";
  rk[w >= 2 & grepl("( group| clade| cluster)$|sensu stricto|incertae sedis| UCG[- ]?[0-9]+$", x, ignore.case = TRUE)] <- "genus";
  rk[w == 2 & grepl("^[^ ]+ ([A-Z]|[0-9])", x)] <- "genus";      # a second word that is not a lower-case epithet: "Pseudoramibacter Eubacterium", "Prevotella 9"
  rk[w <= 1 & grepl("(aceae|ales|ineae|ota|etes|bacteria)$", x, ignore.case = TRUE)] <- "higher";
  rk[is.na(x) | !nzchar(x)] <- NA_character_;
  rk
}

#'The rank a list of taxon names is written at.
#'@description Reads each name's rank from the name itself (one word a genus, two a species, more a strain; see the notes above
#'\code{.tsea.taxon.rank}). A list is "at" a rank when at least 90 percent of its names are; otherwise it is "mixed" and the counts say its mix.
#'@param x character vector of taxon names (duplicates are counted once; blanks and NA are left out).
#'@return list(rank = "genus" | "species" | "strain" | "higher" | "mixed" (NA for an empty list), n.by.rank = named integer counts of the
#'ranks present, n = number of distinct names).
#'@export
TaxonRankSummary <- function(x){
  x <- unique(as.character(x));
  rk <- .tsea.taxon.rank(x);
  rk <- rk[!is.na(rk)];
  n.by <- vapply(.tsea.rank.levels, function(l) sum(rk == l), integer(1));
  n.by <- n.by[n.by > 0];
  n <- sum(n.by);
  rank <- if(n == 0) NA_character_ else if(max(n.by) / n >= .tsea.rank.pure) names(n.by)[which.max(n.by)] else "mixed";
  list(rank = rank, n.by.rank = n.by, n = n)
}

#'The variant of a taxon-set library a list of this rank is tested against.
#'@description The ONE place that decides it, called by \code{CalculateHyperScore} and \code{CalculateGseaScore} (so the web app and
#'the AI workflows steer alike): a genus-, species- or strain-level list goes to the "_genus", "_species" or "_strain" variant of the
#'library chosen, when the library has one (one rule on the rank suffix; a variant whose file is not shipped is tried and, when it is
#'not there, the mixed library is used and the result says so). A mixed or higher-rank list, a library with no such variant, an unknown
#'library, and a level the user chose themselves (a library that already ends in _genus, _species or _strain) stay as they are.
#'@param tset.type the library key the user chose (see \code{SetTaxonSetLib}).
#'@param rank the list's rank, from \code{TaxonRankSummary}.
#'@return the library key to test against.
#'@export
ResolveTaxonSetLib <- function(tset.type, rank){
  if(length(tset.type) != 1 || is.na(tset.type) || !(tset.type %in% names(.tsea.lib.files))) return(tset.type);
  if(length(rank) != 1 || is.na(rank)) return(tset.type);
  if(grepl("_(genus|species|strain)$", tset.type)) return(tset.type);
  if(!(rank %in% .tsea.variant.ranks)) return(tset.type);
  cand <- paste0(tset.type, "_", rank);
  if(cand %in% names(.tsea.lib.files)) cand else tset.type
}

# Steer the library to the list's rank and say what was tested: returns list(set = mbSetObj, rec = the fields to add to analSet$ora.background).
# list.names = the list tested (the query of a significant list, the whole of a complete or ranked list); background.names = a separate
# measured list, whose own mix of ranks is recorded apart. Reloads the library globals only when the variant changes.
.tsea.steer <- function(mbSetObj, list.names, background.names = NULL){
  rs <- TaxonRankSummary(list.names);
  rec <- list(rank = rs$rank, n.by.rank = rs$n.by.rank);
  req <- mbSetObj$dataSet$tset.type;
  if(length(req) == 1 && !is.na(req) && req %in% names(.tsea.lib.files)){
    used <- ResolveTaxonSetLib(req, rs$rank);
    have <- mbSetObj$dataSet$tset.used; if(is.null(have)) have <- req;
    if(!identical(used, have)){
      if(isTRUE(tryCatch(.tsea.load.lib(used, quiet = TRUE), error = function(e) FALSE))) have <- used     # a variant that cannot be read (not shipped yet) leaves the library as it was, and the record says which was wanted
      else rec$library.variant.missing <- used;
    } else if(rs$rank %in% .tsea.variant.ranks && identical(used, req) && !grepl("_(genus|species|strain)$", req)){
      rec$library.variant.missing <- paste0(req, "_", rs$rank);      # this library has no variant at the list's rank at all
    }
    mbSetObj$dataSet$tset.used <- have;
    rec$library <- have;
    rec$library.requested <- req;
    rec$library.steered <- !identical(have, req);
    rec$library.sets <- length(current.mset);
    rec$library.rank <- TaxonRankSummary(unlist(current.mset, use.names = FALSE))$rank;
  }
  if(length(background.names) > 0){
    bs <- TaxonRankSummary(background.names);
    rec$background.rank <- bs$rank;
    rec$background.n.by.rank <- bs$n.by.rank;
  }
  list(set = mbSetObj, rec = rec)
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
  mic_int_species   = "tsea_microbiome_int_species.csv",
  food_matrix_species = "tsea_food_matrix_species.csv",
  # strain level
  host_int_strain   = "tsea_host_int_strain.csv",
  host_ext_strain   = "tsea_host_ext_strain.csv",
  host_diet_strain  = "tsea_host_diet_lifestyle_strain.csv",
  host_drug_strain  = "tsea_host_medication_strain.csv",
  env_strain        = "tsea_environment_strain.csv",
  mic_int_strain    = "tsea_microbiome_int_strain.csv",
  food_matrix_strain = "tsea_food_matrix_strain.csv",
  # genus level: the library curation ships a _genus file for every category (mic_met and gene included). A genus list is tested against
  # the variant when the file is there and against the mixed library, said so, when it is not.
  host_int_genus    = "tsea_host_int_genus.csv",
  host_ext_genus    = "tsea_host_ext_genus.csv",
  host_diet_genus   = "tsea_host_diet_lifestyle_genus.csv",
  host_drug_genus   = "tsea_host_medication_genus.csv",
  env_genus         = "tsea_environment_genus.csv",
  mic_met_genus     = "taxon_metabolite_tsea_genus.csv",
  mic_int_genus     = "tsea_microbiome_int_genus.csv",
  gene_genus        = "tsea_host_snps_genus.csv",
  food_matrix_genus = "tsea_food_matrix_genus.csv"
);

# Resolve a taxon-set library file to a local path. On the web server the file is
# read from the resources tree. In package mode it is downloaded once into the
# working directory and re-used for 30 days (same policy as
# .read.microbiomeanalyst.lib.rds), instead of being fetched on every call.
.get.tsea.lib.path <- function(filenm, quiet = FALSE){
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
        if(!quiet) AddErrMsg(paste("Could not download taxon set library", filenm, "- check your internet connection."));
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
  if(!.tsea.load.lib(tset.type)){
    return(0);
  }
  # the variant the list is tested against: what the user chose, until CalculateHyperScore / CalculateGseaScore steer it (ResolveTaxonSetLib)
  mbSetObj$dataSet$tset.used <- tset.type
  return(.set.mbSetObj(mbSetObj))
}

# Read one taxon-set library into the globals the tests use (current.msetlib, current.mset, uniq.count). TRUE when it was read.
# quiet = TRUE: a library that is not there is not an error (a variant the curation has not shipped yet); the caller falls back.
.tsea.load.lib <- function(tset.type, quiet = FALSE){

  filenm <- .tsea.lib.files[tset.type];
  if(is.na(filenm)){
    if(!quiet) AddErrMsg(paste("Unknown taxon set library:", tset.type));
    return(FALSE);
  }
  # package mode: a variant already found missing in this session is not downloaded for again (on the web server it is a cheap file check)
  if(quiet && !.on.public.web && isTRUE(.tsea.cache$missing[[filenm]])) return(FALSE);
  libPath <- .get.tsea.lib.path(filenm, quiet);
  if(is.null(libPath) || (.on.public.web && !file.exists(libPath))){
    if(quiet) .tsea.cache$missing[[filenm]] <- TRUE;
    if(!quiet && !is.null(libPath)) AddErrMsg(paste("Taxon set library", filenm, "is not available."));
    return(FALSE);
  }

  lib <- .readDataTable(libPath);
  # the reader makes the text "NA" a missing value, and a set named NA would answer to no name: the set keeps the text
  if(anyNA(lib[, 1])) lib[is.na(lib[, 1]), 1] <- "NA";
  current.msetlib <<- lib;
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
  return(TRUE)
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
  # one row per pair of sets: the mirrored matrix lists A-B and B-A, which made two edges for every link
  wd <- wd[match(as.character(wd[,1]), id) < match(as.character(wd[,2]), id),];
  wd <- wd[!is.na(wd[,3]),];
  g <- graph_from_data_frame(wd[,-3], directed=F, vertices=id);
  E(g)$width <- sqrt(wd[,3]*20);
  g <- delete_edges(g, E(g)[wd[,3] < 0.2]);
  idx <- unlist(sapply(V(g)$name, function(x) match(x,id)));
  pvalue <- pvalue[idx]
  cols <- color_scale("red", "#E5C494");
  V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))];

  # a set's size is its number of hits, spread over the range the enrichment network viewer draws (8 to 20), so the sets with the most and the
  # fewest hits differ whatever the counts are
  sz <- .enrnet.rescale(log10(hits + 1), 8, 20, 12);
  names(sz) <- id;
  V(g)$size <- sz[V(g)$name];
    
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
  # seq_len, not 1:nrow: when no two sets overlap enough to keep an edge, 1:0 is c(1, 0) and the file carried two
  # edges with no ends, so the network drew nothing. No edge left is an empty list; the sets still show as nodes.
  edges <- if (nrow(edge.mat) == 0) list() else
    cbind(id=seq_len(nrow(edge.mat)), source=edge.mat[,1], target=edge.mat[,2]);

  # covert to json
  netData <- list(nodes=nodes, edges=edges);
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
  # the members the set was tested with; NULL for a result saved before they were recorded (then none can be called untested)
  tested <- mbSetObj$analSet$ora.tested[[msetNm]];
  # highlighting with different colors: hits red, members that were not tested grey
  members <- current.mset[[msetNm]];
  mset <- members;
  red.inx <- which(members %in% hits[[msetNm]]);
  grey.inx <- if(is.null(tested)) integer(0) else setdiff(which(!(members %in% tested)), red.inx);
  mset[red.inx] <- paste("<font color=\"red\">", "<b>", members[red.inx], "</b>", "</font>",sep="");
  mset[grey.inx] <- paste("<font color=\"grey\" title=\"not tested\">", "<b>", members[grey.inx], "</b>", "</font>",sep="");

  # the library keeps every known member of a set; the table counts only those in the list tested, so the view says both
  n.mem <- length(members);
  said <- paste0("The library holds ", n.mem, if(n.mem == 1) " member" else " members", " of this set");
  if(is.null(tested)){
    said <- paste0(said, ".");
  }else{
    n.tested <- sum(members %in% tested);
    said <- paste0(said, "; ", n.tested, " of them ", if(n.tested == 1) "is" else "are", " in the list tested here, and only those count in the table (Total tested and Hits). ",
                   "Hits are red; grey members are not in the tested list, so they were not tested.");
  }

  # get references
  matched.inx <- match(tolower(msetNm), tolower(current.msetlib$name))[1];

  return(cbind(msetNm, paste0("<p style=\"font-size:12px\"><i>", said, "</i></p>", paste(mset, collapse="; ")), current.msetlib$reference[matched.inx]));
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

