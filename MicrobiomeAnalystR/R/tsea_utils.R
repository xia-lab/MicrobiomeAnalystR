######################################
###########TSEA#############
#######################################

Setup.MapData<-function(qvec){
    dataSet$species <- qvec;
    dataSet <<- dataSet;
}

# given a list of species names or ids, find matched name or ids from selected databases
CrossReferencing <- function(q.type, module.type){
   # record all the data
   name.map <<- list();
   # distribute job
   dataSet$q.type <- q.type;
   dataSet$module.type<-module.type;
   dataSet <<- dataSet;
   SpeciesMappingExact(q.type);

    # do some sanity check
   if(length(which(is.na(name.map$hit.inx)))/length(name.map$hit.inx) > 0.75){
        nmcheck.msg <<- c(1, "Over 3/4 of the IDs could not be matched to our database. Please make 
                        sure that correct taxonomy IDs or common taxa names are used.");        
    }else{
        nmcheck.msg <<- c(1, "Name matching OK, please inspect (and manual correct) the results then proceed.");   
    }  
}

# Mapping from different metabolite IDs
# For compound names to other id, can do exact or approximate match
# For other IDs, except HMDB ID, all other may return multiple /non-unique hits
# multiple hits or non-unique hits will all users to manually select
SpeciesMappingExact<-function(q.type){
       qvec <- dataSet$species;
       # local variable to save memory
       species.db <- readRDS("../../lib/tsea/microbe_db.rds")
       # variables to record results
       hit.inx = vector(mode='numeric', length=length(qvec)); # record hit index, initial 0
       match.values = vector(mode='character', length=length(qvec)); # the best matched values (hit names), initial ""
       match.state = vector(mode='numeric', length=length(qvec));  # match status - 0, no match; 1, exact match; initial 0 
       if(q.type == "gold"){
            hit.inx <- match(tolower(qvec), tolower(species.db$GOLD_ID));
            match.values <- species.db$GOLD_ID[hit.inx];
            match.state[!is.na(hit.inx)] <- 1;
       }else if(q.type == "taxa"){
            hit.inx <- match(tolower(qvec), tolower(species.db$taxa));
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

       qvec <- dataSet$species;
       # style for highlighted background for unmatched names
       pre.style<-NULL;
       post.style<-NULL;

       # style for no matches
       no.prestyle<-"<strong style=\"background-color:yellow; font-size=125%; color=\"black\">";
       no.poststyle<-"</strong>";
    
       hit.inx<-name.map$hit.inx;
       hit.values<-name.map$hit.values;
       match.state<-name.map$match.state;

       # construct the result table with cells wrapped in html tags
       # the unmatched will be highlighted in different background
       html.res<-matrix("", nrow=length(qvec), ncol=6);
       colnames(html.res)<-c("Query", "Match","Species","Genus","NCBI_Taxonomy_ID","GOLDSTAMP_ID");

       for (i in 1:length(qvec)){
         if(match.state[i]==1){
           pre.style<-"";
           post.style="";
           
           }else{ # no matches
               pre.style<-no.prestyle;
               post.style<-no.poststyle;
           }
           hit <-species.db[hit.inx[i], ,drop=F];

           html.res[i, ]<-c(paste(pre.style, qvec[i], post.style, sep=""),
                            paste(ifelse(match.state[i]==0, "", hit.values[i]), sep=""),
                            paste(ifelse(match.state[i]==0 || is.na(hit$species) ||is.null(hit$species) || hit$species=="" || hit$species=="NA","-",hit$species),sep=""),
                            paste(ifelse(match.state[i]==0 || is.na(hit$genus) ||is.null(hit$genus) || hit$genus=="" || hit$genus=="NA","-",hit$genus),  sep=""), 
                            paste(ifelse(match.state[i]==0 || is.na(hit$NCBITAX) ||is.null(hit$NCBITAX) || hit$NCBITAX=="" || hit$NCBITAX=="NA","-", paste("<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=", hit$NCBITAX," target='_blank'>", hit$NCBITAX,"</a>", sep="")),  sep=""),
                            paste(ifelse(match.state[i]==0 || is.na(hit$GOLDMAPID) ||is.null(hit$GOLDMAPID) || hit$GOLDMAPID=="" || hit$GOLDMAPID=="NA", "-", paste("<a href=https://gold.jgi.doe.gov/project?id=", hit$GOLDMAPID," target='_blank'>", hit$GOLDMAPID,"</a>", sep="")), sep=""))
        }

     analSet$resTable <- data.frame(html.res);
     analSet$mapTable<-cbind(Query=qvec,analSet$resTable[,2:ncol(analSet$resTable)]);
     analSet <<- analSet;
}

CalculateHyperScore <- function(){

    nm.map <- GetFinalNameMap();
    valid.inx <- !(is.na(nm.map$Strain)| duplicated(nm.map$Strain));
    ora.vec <- nm.map$Strain[valid.inx];
    q.size <- length(ora.vec);

    if(is.na(ora.vec) || q.size==0) {
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

        if(dataSet$tset.type=="host_int_species"|dataSet$tset.type=="env_species"|dataSet$tset.type=="host_ext_species"){
           AddErrMsg("Species-level taxa set was selected: verify that your list contains species names!");
        }else if(dataSet$tset.type=="host_int_strain"|dataSet$tset.type=="env_strain"|dataSet$tset.type=="mic_int_strain"){
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

    # adjust for multiple testing problems
    res.mat[,5] <- p.adjust(res.mat[,4], "holm");
    res.mat[,6] <- p.adjust(res.mat[,4], "fdr");

    res.mat <- res.mat[hit.num>0,];

    # fix error when only 1 hit (not sig), no longer a matrix
    if(class(res.mat)!= "matrix"){
      AddErrMsg("No significant hits found using enrichment analysis!");
      return(0);
    }

    ord.inx<-order(res.mat[,4]);

    # download result
    
    analSet$ora.mat = signif(res.mat[ord.inx,],3);
    analSet$ora.hits = hits;
    write.csv(analSet$ora.mat, file="tsea_ora_result.csv");
    analSet <<- analSet;

    return(1);
}

# return the final (after user selection) map as dataframe
# two col, original name and strain ,
# for Enrichment and pathway analysis, respectively
GetFinalNameMap<-function(){
    enrtype <- dataSet$q.type;
    qvec <- dataSet$species;
    nm.mat<-matrix(nrow=length(qvec), ncol=2);
    colnames(nm.mat)<-c("query", "Strain");
    if(enrtype=="taxa"){
        for (i in 1:length(qvec)){
            nm.mat[i, ]<-c(qvec[i],qvec[i]);
        }
        return(as.data.frame(nm.mat));
    }else{
        hit.inx<-name.map$hit.inx;
        hit.values<-name.map$hit.values;
        match.state<-name.map$match.state;
        species.db <- readRDS("../../lib/tsea/microbe_db.rds");
        for (i in 1:length(qvec)){
            hit <-species.db[hit.inx[i], ,drop=F];
            if(match.state[i]==0){
                kegg.hit <- NA;
            }else{
                kegg.hit <- ifelse(nchar(hit$organism.name)==0, NA, hit$organism.name);
            }
            nm.mat[i, ]<-c(qvec[i], kegg.hit);
        }
        return(as.data.frame(nm.mat));
    }
}

PrepareEnrichNet<-function(){
    #calculate the enrichment fold change
    folds <- analSet$ora.mat[,3]/analSet$ora.mat[,2];
    names(folds)<-GetShortNames(rownames(analSet$ora.mat));
    hits <-  analSet$ora.mat[,3];
    pvals <- analSet$ora.mat[,4];
    PlotEnrichNet.Overview(hits, pvals);
}

GetORA.rowNames<-function(){
    nms <- rownames(analSet$ora.mat);
    if(is.null(nms)){
        return("NA");
    }
    return (nms);
}

GetORA.mat<-function(){
    return(analSet$ora.mat);
}

GetORA.colorBar<-function(){
    len <- nrow(analSet$ora.mat);
    if(len > 50){
        ht.col <- c(substr(heat.colors(50), 0, 7), rep("#FFFFFF", len-50));
    }else{
        # reduce to hex by remove the last character so HTML understand
        ht.col <- substr(heat.colors(len), 0, 7);
    }
    return (ht.col);
}

SetMsetLib <- function(tset.type){

    dataSet$tset.type <- tset.type
    dataSet <<- dataSet;

    if(tset.type=="host_int"){
      libPath <- "../../lib/tsea/tsea_host_int.csv";
    }else if(tset.type=="host_ext"){
      libPath <- "../../lib/tsea/tsea_host_ext.csv";
    }else if(tset.type=="env"){
      libPath <- "../../lib/tsea/tsea_environment.csv";
    }else if(tset.type=="mic_int"){
      libPath <- "../../lib/tsea/tsea_microbiome_int.csv";
    }else if(tset.type=="gene"){
      libPath <- "../../lib/tsea/tsea_host_snps_new.csv";
    }else if(tset.type=="host_int_species") {
      libPath <- "../../lib/tsea/tsea_host_int_species.csv";
    }else if(tset.type=="env_species"){
      libPath <- "../../lib/tsea/tsea_environment_species.csv";
    }else if(tset.type=="host_ext_species"){
      libPath <- "../../lib/tsea/tsea_host_ext_species.csv";
    }else if(tset.type=="host_int_strain"){
      libPath <- "../../lib/tsea/tsea_host_int_strain.csv";
    }else if(tset.type=="env_strain"){
      libPath <- "../../lib/tsea/tsea_environment_strain.csv";
    }else{
      libPath <- "../../lib/tsea/tsea_microbiome_int_strain.csv";
    }
      
    current.msetlib <<- .readDataTable(libPath);
    
    ms.list<-strsplit(current.msetlib[,2],"; ");
    names(ms.list)<-current.msetlib[,1];
    current.mset <<- ms.list;
    # total uniq cmpds in the mset lib
    uniq.count <<- length(unique(unlist(current.mset, use.names = FALSE)));
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

# barplot height is enrichment fold change
# color is based on p values
PlotEnrichNet.Overview<-function(hits, pvals){

    # due to space limitation, plot top 50 if more than 50 were given
    title <- "Taxon Set Enrichment Network Overview";
    if(length(hits) > 50){
        hits <- hits[1:50];
        pvals <- pvals[1:50];
        title <- "Enrichment Overview (top 50)";
    }

    library(igraph);
    library(reshape);
    pvalue <- pvals;
    id <- names(pvalue);
    geneSets <- current.mset;
    n <- length(pvalue);
    w <- matrix(NA, nrow=n, ncol=n);
    colnames(w) <- rownames(w) <- id;

    for (i in 1:n) {
        for (j in i:n) {
            w[i,j] = overlap_ratio(geneSets[id[i]], geneSets[id[j]])
        }
    }

    wd <- melt(w);
    wd <- wd[wd[,1] != wd[,2],];
    wd <- wd[!is.na(wd[,3]),];
    g <- graph.data.frame(wd[,-3], directed=F);
    E(g)$width <- sqrt(wd[,3]*20);
    g <- delete.edges(g, E(g)[wd[,3] < 0.2]);
    idx <- unlist(sapply(V(g)$name, function(x) which(x == id)));
    cols <- color_scale("red", "#E5C494");
    V(g)$color <- cols[sapply(pvalue, getIdx, min(pvalue), max(pvalue))];

    cnt <- hits + 2;
    names(cnt) <- id;
    cnt2 <- cnt[V(g)$name];
    V(g)$size <- log(cnt2, base=10) * 10; ## cnt2/sum(cnt2) * 100;
    #V(g)$size <- cnt2/sum(cnt2) * 10;
    
    # layout
    pos.xy <- layout.fruchterman.reingold(g);

    # now create the json object
    nodes <- vector(mode="list");
    node.nms <- V(g)$name;
    node.sizes <- V(g)$size;
    node.cols <- V(g)$color;
    for(i in 1:length(node.sizes)){
        nodes[[i]] <- list(
                  id = node.nms[i],
                  label=node.nms[i], 
                  size=node.sizes[i], 
                  color=node.cols[i],
                  x = pos.xy[i,1],
                  y = pos.xy[i,2]
                );
    }
    edge.mat <- get.edgelist(g);
    edge.mat <- cbind(id=1:nrow(edge.mat), source=edge.mat[,1], target=edge.mat[,2]);
    # covert to json
    library(RJSONIO);
    netData <- list(nodes=nodes, edges=edge.mat);
    sink("tsea_network.json");
    cat(RJSONIO::toJSON(netData));
    sink();
}

overlap_ratio <- function(x, y) {
    x <- unlist(x)
    y <- unlist(y)
    length(intersect(x, y))/length(unique(c(x,y)))
}

##' @importFrom grDevices colorRampPalette
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

GetCurrentImg<-function(){
    return (current.img);
}

# given a metset inx, return hmtl highlighted metset cmpds and references
GetHTMLMetSet<-function(msetNm){
    hits <- analSet$ora.hits;
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

GetMsetPval<-function(msetNm){
    return(analSet$ora.mat[msetNm, "Raw p"]);
}

# given a metset inx, return all info
GetTaxaSet<-function(msetNm){
  mset <- subset(current.msetlib, name == msetNm)
  analSet$tseaInfo <- t(rbind(colnames(current.msetlib), mset))
  print(analSet$tseaInfo)
  analSet <<- analSet
  return(1);
}
