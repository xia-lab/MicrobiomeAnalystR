##################################################
## R scripts for MicrobiomeAnalyst
## Various utility methods
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

GetRandomNumbers <- function(){
  rm(.Random.seed);
  runif(1);
  return(.Random.seed[3:626]);
}

# need to obtain the full path to convert (from imagemagik) for cropping images
GetBashFullPath<-function(){
  path <- system("which bash", intern=TRUE);
    
  if((length(path) == 0) && (typeof(path) == "character")){
    print("Could not find bash in the PATH!");
    return("NA");
  }
  return(path);
}

# tell if window/mac/unix
GetOS <- function(){
  sysinf <- Sys.info()
  if (!is.null(sysinf)){
    os <- sysinf['sysname']
    if (os == 'Darwin')
      os <- "mac"
  } else { ## mystery machine
    os <- .Platform$OS.type
    if (grepl("^darwin", R.version$os))
      os <- "mac"
    if (grepl("linux-gnu", R.version$os))
      os <- "linux"
  }
  tolower(os)
}

UnzipUploadedFile<-function(zip_file){
    dataName <- unzip(zip_file, list = TRUE)$Name;
    if(length(dataName) > 1){
      # test if "__MACOSX" or ".DS_Store"
      osInx <- grep('MACOSX',dataName,perl=TRUE);
      if(length(osInx) > 0){
        dataName <- dataName[-osInx];
      }
      dsInx <- grep('DS_Store',dataName,perl=TRUE);
      if(length(dsInx) > 0){
        dataName <- dataName[-dsInx];
      }
      if(length(dataName) != 1){
        current.msg <<- "More than one data files found in the zip file.";
        print(dataName);
        return("NA");
      }
    }
    a <- try(unzip(zip_file));
    if(class(a) == "try-error" | length(a)==0){
      current.msg <<- "Failed to unzip the uploaded files!";
      return ("NA");
    }
    return(dataName);
}

# get qualified inx with at least number of replicates
GetDiscreteInx <- function(my.dat, min.rep=2){
  good.inx <- apply(my.dat, 2, function(x){
                    good1.inx <- length(x) > length(unique(x));
                    good2.inx <- min(table(x)) >= min.rep;
                    return (good1.inx & good2.inx);
            });
   return(good.inx);
}

# get columns that are "most likely" continuous values
GetNumbericalInx <- function(my.dat){
  good.inx <- apply(my.dat, 2, function(x){
                return(all(!is.na(as.numeric(as.character(x))))); 
            });
   return(good.inx);
}

###########
# improved list of objects
.ls.objects <- function (pos = 1, pattern, order.by,
                        decreasing=FALSE, head=FALSE, n=5) {
  
  napply <- function(names, fn) sapply(names, function(x)
                                         fn(get(x, pos = pos)))
  names <- ls(pos = pos, pattern = pattern)
  obj.class <- napply(names, function(x) as.character(class(x))[1])
  obj.mode <- napply(names, mode)
  obj.type <- ifelse(is.na(obj.class), obj.mode, obj.class)
  obj.prettysize <- napply(names, function(x) {
                           capture.output(format(utils::object.size(x), units = "auto")) })
  obj.size <- napply(names, object.size)
  obj.dim <- t(napply(names, function(x)
                        as.numeric(dim(x))[1:2]))
  vec <- is.na(obj.dim)[, 1] & (obj.type != "function")
  obj.dim[vec, 1] <- napply(names, length)[vec]
  out <- data.frame(obj.type, obj.size, obj.prettysize, obj.dim)
  names(out) <- c("Type", "Size", "PrettySize", "Rows", "Columns")
    
  if (!missing(order.by))
      out <- out[order(out[[order.by]], decreasing=decreasing), ]
  if (head)
      out <- head(out, n)
  out
}

# shorthand
ShowMemoryUse <- function(..., n=20) {
 
  sink(); # make sure print to screen
  print(pryr::mem_used());
  print(sessionInfo());
  print(.ls.objects(..., order.by="Size", decreasing=TRUE, head=TRUE, n=n));
  print(warnings());
}

GetExtendRange<-function(vec, unit=10){
  var.max <- max(vec, na.rm=T);
  var.min <- min(vec, na.rm=T);
  exts <- (var.max - var.min)/unit;
  c(var.min-exts, var.max+exts);
}

ComputeColorGradient <- function(nd.vec, centered=TRUE){
  
  load_rcolorbrewer();
  if(sum(nd.vec<0, na.rm=TRUE) > 0){
    centered <- T;
  }else{
    centered <- F;
  }
    
  color <- grDevices::colorRampPalette(c("green", "yellow", "red"))(100);
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}

generate_breaks = function(x, n, center = F){
    
  if(center){
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }else{
        res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1)
  }
  return(res)
}

scale_vec_colours = function(x, col = rainbow(10), breaks = NA){
  breaks <- sort(unique(breaks));
  return(col[as.numeric(cut(x, breaks = breaks, include.lowest = T))])
}

scale_colours = function(mat, col = rainbow(10), breaks = NA){
  mat = as.matrix(mat)
  return(matrix(scale_vec_colours(as.vector(mat), col = col, breaks = breaks), nrow(mat), ncol(mat), dimnames = list(rownames(mat), colnames(mat))))
}

#'Read data table
#'@description Function to read in a data table. First, it will try to use fread, however, it has issues with 
#'some windows 10 files. In such case, use the slower read.table method.
#'@param fileName Input filename
#'@author Jeff Xia\email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
.readDataTable <- function(fileName){
  dat <- try(data.table::fread(fileName, header=TRUE, check.names=FALSE, blank.lines.skip=TRUE, data.table=FALSE));
  if(class(dat) == "try-error" || any(dim(dat) == 0)){
      print("Using slower file reader ...");
      formatStr <- substr(fileName, nchar(fileName)-2, nchar(fileName))
      if(formatStr == "txt"){
        dat <- try(read.table(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
      }else{ # note, read.csv is more than read.table with sep=","
        dat <- try(read.csv(fileName, header=TRUE, comment.char = "", check.names=F, as.is=T));
      }
  }
  return(dat);
}

#' Adds an error message
#'@description The error message will be printed in all cases.
#'Used in higher functions.
#'@param msg Error message to print
#'@export
AddErrMsg <- function(msg){
  current.msg <<- msg;
  print(msg);
}

# this is to be validated (data.table), not done yet, very promising need to validate
# https://github.com/joey711/phyloseq/issues/517
fast_tax_glom <- function(physeq, taxrank, NArm=TRUE, bad_empty=c(NA, "", " ", "\t")){

	# Make a vector from the taxonomic data.
	CN  <- which(rank_names(physeq) %in% taxrank);
  tax.lvls <- rank_names(physeq)[1:CN];

	tax <- as(access(physeq, "tax_table"), "matrix")[, CN];
	if(NArm){ # if NArm is TRUE, remove the empty, white-space, NA values from
    keep_species <- names(tax)[ !(tax %in% bad_empty) ]
    physeq <- prune_taxa(keep_species, physeq)
	}

  taxS4_dt <- data.table(physeq@tax_table,keep.rownames = TRUE,key = "rn");
  otuS4_dt <- as.data.table(physeq@otu_table,keep.rownames = TRUE,key = "rn");
  tax_dt  <- as.data.table(lapply(taxS4_dt, unclass));
  otu_dt  <- as.data.table(lapply(otuS4_dt, unclass));
  otu_nms <- setdiff(names(otu_dt), 'rn');
  otab_dt <- merge(otu_dt, tax_dt, by = "rn");

  setkeyv(otab_dt, tax.lvls);
  o_wm <- copy(otab_dt);
  o_wm[, taxa_sums := Reduce(`+`, .SD), .SDc = otu_nms];
  o_wm[, taxa_max   := max(taxa_sums),   by = key(o_wm)];
  o_wm[, c(otu_nms) := lapply(.SD, sum), by = key(o_wm), .SDc = otu_nms];
  o_wm2 <- o_wm[taxa_max == taxa_sums, .SD, .SDc = c('rn')];

        # TODO here to recreate phyloseq

        # now recreate the new phyloseq object
        #physeqNew <- phyloseq(otu_table(countsNew, taxa_are_rows = TRUE),
        #                tax_table(TaxTableNew),
        #                sample_data(physeq));

	# "Empty" the values to the right of the rank, using NA_character_.
	if( CN < length(rank_names(physeqNew)) ){
		badcolumns <- (CN+1):length(rank_names(physeqNew))
		tax_table(physeqNew)[, badcolumns] <- NA_character_
	}
	return(physeqNew)
}

# same as fast_tax_glom_first, but cached, no repeat computing
fast_tax_glom_mem <- function(physeq, taxrank){
   # set up the cache for reuse (i.e updating colors)
   if(!exists("FastTaxGlom_mem")){
        require("memoise");
        FastTaxGlom_mem<<- memoise(fast_tax_glom_first); # the cache will be empty after 5 min
    }
    res <- FastTaxGlom_mem(physeq, taxrank);
    return(res);
}

fast_tax_glom_first <- function(physeq, taxrank){

  # setup data. we are going to glob the OTU table based on the Class Taxonomy
  CN  <- which(rank_names(physeq) %in% taxrank);
  if(length(CN) == 0){
      print("name not in the rank!");
      return(NULL);
  }
  tax <- as(access(physeq, "tax_table"), "matrix")[, 1:CN, drop=FALSE];
  tax <- apply(tax, 1, function(i){paste(i, sep=";_;", collapse=";_;")});
  # using Map-Reduce/vectorized
  otab2 <- data.frame(otu_table(physeq))
  taxdf <- data.frame(tax);
  otab2 <- merge(otab2, taxdf, by = "row.names");
  row.names(otab2) <- otab2$Row.names
  otab2 <- otab2[ , 2:ncol(otab2)]
  otab2 <- condenseOTUs(otab2,"tax");
  otab2<-otu_table(otab2,taxa_are_rows = T);
  colnames(otab2)<-sample_names(physeq);

  # check if this is meaningful 
  if(sum(is.na(otab2))/length(otab2) > 0.7){
    AddErrMsg("More than 70% values are missing at this level!");
    return(NULL);
  }

  if(length(phy_tree(physeq,errorIfNULL = FALSE))==0){
    phy_data<-merge_phyloseq(otab2,tax_table(physeq),sample_data(physeq));
  }else{
    phy_data<-merge_phyloseq(otab2,tax_table(physeq),sample_data(physeq),phy_tree(physeq));
  }
  return(phy_data)
}

#'  Helper Function that will take the OTU Table (as a data.frame())
#'  to which column representing taxa has been added.
#'  There is some munging to handle rownames and
#'  the tax column but the basic idea is there.
# convert tax col to numeric
# needed to allow colSums to work
condenseOTUs <- function(otutable, splitcol) {

  otutable[[splitcol]] <- as.numeric(factor(otutable[[splitcol]], labels=c(1:length(unique(otutable[[splitcol]])))))

  # split apply, combine
  splits <- split(otutable, otutable[[splitcol]])
  summed <- Map(colSums, splits)
  summeddf <- Reduce(rbind, summed, init = NULL)

  #get rownames.
  # this needs to be changed to get most abundant OTU if thats what is used elsewhere
  # rightnow it returns the first OTU in the group
  newrownames <- Map(function(x){rownames(x)[[1]]}, splits)

  #add back rowname and remove tax column
  rownames(summeddf) <- newrownames
  summeddf[, !colnames(summeddf) %in% c(splitcol)]
}

# need to return consistent color assignments for the same taxa
#'@export
GetSeriesColors <- function(taxa=NULL){
  if(!exists("pie.cols")){
    InitPieColors();
  }
    
  if(is.null(taxa)){
    taxa <- as.character(piedata$variable);
  }

  col.len <- length(pie.cols);
  taxa.len <- length(taxa);
  
  while(taxa.len > col.len){
    pie.cols <- c(pie.cols, pie.cols);
    col.len <- length(pie.cols);
  }

  # need to see if already assigned in the previous analysis
  hit.inx <- taxa %in% names(pie.cols);
    
  if(sum(!hit.inx) > 0){
    taxa1 <- taxa[!hit.inx];
    # get unused colors index in the original pie.cols
    hit.inx2 <- which(!(names(pie.cols) %in% taxa));
    # assign them to next in line
    free.inx <- hit.inx2[1:length(taxa1)];
    names(pie.cols)[free.inx] <- taxa1;
  }
    
  pie.cols <<- pie.cols;
  return(paste(pie.cols[taxa], collapse=","));
}

InitPieColors <- function(){
  pie.cols <- "ff6347,ffff00,ee82ee,f4a460,2e8b57,4682b4,9acd32,ffe4c4,8a2be2,a52a2a,deb887,5f9ea0,7fff00,d2691e,ff7f50,7fffd4,6495ed,fff8dc,dc143c,00ffff,00ff7f,00008b,b8860b,a9a9a9,006400,bdb76b,556b2f,9932cc,e9967a,8fbc8f,483d8b,2f4f4f,00ced1,9400d3,ff1493,00bfff,696969,1e90ff,b22222,fffaf0,228b22,ff00ff,dcdcdc,f8f8ff,ffd700,daa520,808080,008000,adff2f,f0fff0,ff69b4,cd5c5c,4b0082,fffff0,f0e68c,e6e6fa,fff0f5,7cfc00,fffacd,add8e6,f08080,e0ffff,fafad2,d3d3d3,90ee90,ffb6c1,ffa07a,20b2aa,87cefa,778899,b0c4de,ffffe0,00ff00,32cd32,faf0e6,800000,66cdaa,0000cd,ba55d3,9370d8,3cb371,7b68ee,00fa9a,48d1cc,c71585,191970,f5fffa,ffe4e1,ffe4b5,ffdead,fdf5e6,6b8e23,ffa500,ff4500,da70d6,eee8aa,98fb98,afeeee,d87093,ffefd5,ffdab9,cd853f,ffc0cb,dda0dd,b0e0e6,ff0000,bc8f8f,4169e1,8b4513,fa8072,fff5ee,a0522d,c0c0c0,87ceeb,6a5acd,708090,fffafa,d2b48c,008080,d8bfd8,40e0d0,f5deb3,ffffff,f5f5f5";
  pie.cols <- strsplit(pie.cols,",")[[1]];
  names(pie.cols) <- paste("nm", 1:length(pie.cols), sep="");
  pie.cols <<- pie.cols;
}

# utils to remove from
# within, leading and trailing spaces
ClearStrings<-function(query){
  # kill multiple white space
  query <- gsub(" +"," ",query);
  # remove leading and trailing space
  query <- sub("^[[:space:]]*(.*?)[[:space:]]*$", "\\1", query, perl=TRUE);
  return(query);
}

cleanMem <- function(n=10) { for (i in 1:n) gc() }

# based on phyloseq post: https://github.com/joey711/shiny-phyloseq/blob/master/panels/paneldoc/Transform.md
clr_transform <- function(x, base=2){
  x <- log((x / gm_mean(x)), base)
  x[!is.finite(x) | is.na(x)] <- 0.0
  return(x)
}

gm_mean <- function(x, na.rm=TRUE){
  # The geometric mean, with some error-protection bits.
  exp(sum(log(x[x > 0 & !is.na(x)]), na.rm=na.rm) / length(x))
}

# generate Latex table
#'@import xtable
GetSigTable<-function(mat, method){
  
  load_xtable();

  if(!isEmptyMatrix(mat)){ # test if empty
    cap<-"Important features identified by";
        
    if(nrow(mat)>50){
      smat<-as.matrix(mat[1:50,]); # only print top 50 if too many
      colnames(smat)<-colnames(mat); # make sure column names are also copied
      mat<-smat;
      cap<-"Top 50 features identified by";
    }
        
    # change the rowname to first column
    col1<-rownames(mat);
    cname<-colnames(mat);
    cname<-c("Features", cname);
    mat<-cbind(col1, mat);
    rownames(mat)<-NULL;
    colnames(mat)<-cname;
    print(xtable(mat, caption=paste(cap, method)), ,caption.placement="top", size="\\scriptsize");
  }else{
    print(paste("No significant features were found using the given threshold for", method));
  }
}

# test if a sig table matrix is empty
isEmptyMatrix<-function(mat){
  if(is.null(mat) | length(mat)==0){
    return(TRUE);
  }
  if(nrow(mat)==0 | ncol(mat)==0){
    return(TRUE);
  }
  if(is.na(mat[1,1])){
    return(TRUE);
  }
  return(FALSE);
}

#######################################
###########Color-palettes##############
#######################################

custom_col42<-  c("#781156","#A51876","#D21E96","#E43FAD","#EA6CC0","#F098D3","#114578","#185EA5","#1E78D2","#3F91E4",
                "#6CABEA","#98C4F0","#117878","#18A5A5","#3FE4E4","#6CEAEA","#98F0F0", "#117845","#18A55E","#1ED278",
                "#3FE491","#6CEAAB","#98F0C4","#787811","#A5A518","#D2D21E","#E4E43F","#EAEA6C","#F0F098","#F7F7C5",
                "#784511","#A55E18","#D2781E","#E4913F","#EAAB6C","#F0C498","#781122","#A5182F","#D21E2C","#E43F5B",
                "#EA6C81","#F098A7");

custom_final28<-c("#771155", "#AA4488", "#EA6CC0", "#CC99BB", "#114477", "#4477AA","#1E78D2", "#77AADD", "#117777",
                "#44AAAA", "#3FE4E4", "#77CCCC", "#117744","#44AA77", "#1ED278", "#88CCAA", "#771122", "#AA4455",
                "#D21E2C","#DD7788","#777711", "#AAAA44", "#D2D21E", "#DDDD77","#774411", "#AA7744", "#D2781E",
                "#DDAA77");

custom_col21<-  c("#771155","#AA4488", "#CC99BB", "#114477", "#4477AA", "#77AADD", "#117777", "#44AAAA", "#77CCCC",
                "#117744", "#44AA77", "#88CCAA", "#777711", "#AAAA44", "#DDDD77", "#774411", "#AA7744", "#DDAA77",
                "#771122", "#AA4455", "#DD7788");
col_vector<-c("#7FC97F","#BEAED4","#FDC086","#FFFF99","#386CB0","#F0027F","#BF5B17","#666666","#A6CEE3","#6A3D9A",
              "#33A02C","#FB9A99","#E31A1C","#FDBF6F","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB",
              "#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7",
              "#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#FFFF99","#FBB4AE","#B3CDE3",
              "#984EA3","#FF7F00","#FFFF33","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2",
              "#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#666666","#1B9E77","#D95F02",
              "#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#CCCCCC","#E41A1C","#377EB8","#4DAF4A","#D9D9D9",
              "#BC80BD","#CCEBC5","#FFED6F","#B2DF8A","#FF7F00","#B15928","#CAB2D6","#1F78B4","#A65628");

ComputeColorGradientCorr <- function(nd.vec, centered=TRUE){
  
  load_rcolorbrewer();
  if(sum(nd.vec<0, na.rm=TRUE) > 0){
    centered <- T;
  }else{
    centered <- F;
  }
    
  color <- c(grDevices::colorRampPalette(c("#003366", "#add8e6"))(50), grDevices::colorRampPalette(c("#FF9DA6", "#FF7783", "#E32636", "#BD0313"))(50));
  breaks <- generate_breaks(nd.vec, length(color), center = centered);
  return(scale_vec_colours(nd.vec, col = color, breaks = breaks));
}


# utility method to get anova results
.do.anova<- function(cls, data, nonpar){
  if(nonpar){
    anova.res <- apply(as.matrix(data), 2, function(x){kruskal.test(x ~ cls)});
    res <- unlist(lapply(anova.res, function(x) {c(x$statistic, x$p.value)}));
    return(matrix(res, nrow=length(anova.res), byrow=T));
  }else{
    aov.res <- apply(as.matrix(data), 2, function(x){aov(x ~ cls)});
    anova.res <- lapply(aov.res, anova);
    res <- unlist(lapply(anova.res, function(x) { c(x["F value"][1,], x["Pr(>F)"][1,])}));
    return(matrix(res, nrow=length(aov.res), byrow=T));
  }
}

# utility method to get p values
#'@import genefilter
.do.ttests<- function(cls, data, nonpar=F){

  data <- na.omit(data); 
  inx1 <- which(cls==levels(cls)[1]);
  inx2 <- which(cls==levels(cls)[2]); 
  
  if(nonpar){
    my.fun <- function(x){
      tmp <- try(wilcox.test(x[inx1], x[inx2]));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    };
  }else{
    my.fun <- function(x) {
      tmp <- try(t.test(x[inx1], x[inx2]));
      if(class(tmp) == "try-error") {
        return(c(NA, NA));
      }else{
        return(c(tmp$statistic, tmp$p.value));
      }
    };
  }
  res <- apply(data, 2, my.fun);
  return(t(res));
}

PerformUnivTests <- function(cls, data, nonpar){
    cls <- as.factor(cls);
    if(length(levels(cls)) > 2){
        return(.do.anova(cls, data, nonpar));
    }else{
        return(.do.ttests(cls, data, nonpar));
    }
}

## fast T-tests/F-tests using C++
PerformFastUnivTests <- function(data, cls, var.equal=TRUE){
    print("Performing fast univariate tests ....");

    # note, feature in rows for gene expression
    data <- t(as.matrix(data));
    cls <- as.factor(cls);
    if(length(levels(cls)) > 2){
        res <- try(rowcolFt(data, cls, var.equal = var.equal));
    }else{
        res <- try(rowcoltt(data, cls, FALSE, 1L, FALSE));
    }  

    if(class(res) == "try-error") {
        res <- cbind(NA, NA);
    }else{
        res <- cbind(res$statistic, res$p.value);
    }

    return(res);
}

###

# This function allows you to vectorise multiple `if` and `else if`
# statements. It is an R equivalent of the SQL `CASE WHEN` statement
#' https://github.com/s-fleck/lest/tree/master/R
case_when <- function(...) {
  formulas <- list(...)
  n <- length(formulas)
  
  if (n == 0) {
    stop("No cases provided")
  }
  
  query <- vector("list", n)
  value <- vector("list", n)
  
  for (i in seq_len(n)) {
    f <- formulas[[i]]
    if (!inherits(f, "formula") || length(f) != 3) {
      stop(sprintf(
        "Case %s (`%s`) must be a two-sided formula, not a %s",
        i,
        deparse_trunc(substitute(list(...))[[i + 1]]),
        typeof(f)
      ))
    }
    
    env <- environment(f)
    query[[i]] <- eval(f[[2]], env)
    
    if (!is.logical(query[[i]])) {
      stop(sprintf(
        "LHS of case %s (%s) must be a logical, not %s",
        i,
        backticks(deparse_trunc(f_lhs(f))),
        typeof(query[[i]])
      ))
    }
    
    value[[i]] <- eval(f[[3]], env)
  }
  
  lhs_lengths <- vapply(query, length, integer(1))
  rhs_lengths <- vapply(value, length, integer(1))
  all_lengths <- unique(c(lhs_lengths, rhs_lengths))
  
  if (length(all_lengths) <= 1) {
    m <- all_lengths[[1]]
  } else {
    non_atomic_lengths <- all_lengths[all_lengths != 1]
    m <- non_atomic_lengths[[1]]
    if (length(non_atomic_lengths) > 1) {
      inconsistent_lengths <- non_atomic_lengths[-1]
      lhs_problems <- lhs_lengths %in% inconsistent_lengths
      rhs_problems <- rhs_lengths %in% inconsistent_lengths
      
      bad_calls(
        formulas[lhs_problems | rhs_problems],
        inconsistent_lengths_message(inconsistent_lengths, m)
      )
    }
  }
  
  out <- value[[1]][rep(NA_integer_, m)]
  replaced <- rep(FALSE, m)
  
  for (i in seq_len(n)) {
    out <- replace_with(out, query[[i]] & !replaced, value[[i]], NULL)
    replaced <- replaced | (query[[i]] & !is.na(query[[i]]))
  }
  
  out
}

replace_with <- function (
  x,
  i,
  val,
  name,
  reason = NULL
){
  if (is.null(val)) {
    return(x)
  }
  
  assert_length_1_or_n(val, length(x), name, reason)
  assert_equal_type(val, x, name)
  assert_equal_class(val, x, name)
  
  i[is.na(i)] <- FALSE
  if (length(val) == 1L) {
    x[i] <- val
  } else {
    x[i] <- val[i]
  }
  x
}

assert_length_1_or_n <- function(
  x,
  n,
  header,
  reason
){
  chk <- check_length_1_or_n(x, n, header, reason)
  
  if (is.null(chk)){
    TRUE
  } else {
    stop(chk)
  }
}

check_length_1_or_n <- function(
  x,
  n,
  header,
  reason
){
  if (length(x) %in% c(1, n)){
    return(NULL)
  }
  
  
  if (is.null(reason))
    reason <- ""
  else
    reason <- paste0(" (", reason, ")")
  
  if (is.null(header))
    header <- ""
  else
    header <- paste0(header, " ")
  
  
  inconsistent_lengths_message(length(x), n, header = header, reason = reason)
}

inconsistent_lengths_message <- function(
  length_is,
  length_should,
  header = "",
  reason = ""
){
  if (length_should == 1) {
    sprintf("%smust be length 1%s, not %s", header, reason, paste(length_is, collapse = ", "))
  } else {
    sprintf("%smust be length %s%s or one, not %s", header, length_should, reason, paste(length_is, collapse = ", "))
  }
}

assert_equal_class <- function(
  x,
  template,
  header
){
  if (!is.object(x)) {
    return(TRUE)
    
  } else if (identical(class(x), class(template))) {
    return(TRUE)
    
  } else {
    
    if (is.null(header))
      header <- ""
    else
      header <- paste0(header, " ")
    
    
    stop(
      sprintf(
        "%smust be %s, not %s",
        header,
        paste(class(template), collapse = "/"),
        paste(class(x), collapse = "/")
      )
    )
  }
}

assert_equal_type <- function(
  x,
  template,
  header
){
  if (identical(typeof(x), typeof(template)))
    return(TRUE)
  
  if (is.null(header))
    header <- ""
  else
    header <- paste0(header, " ")
  
  stop(sprintf("%smust be type %s, not %s", header, typeof(template), typeof(x)))
}

# in public web, this is done by microservice
.perform.computing <- function(){
    dat.in <- qs::qread("dat.in.qs"); 
    dat.in$my.res <- dat.in$my.fun();
    qs::qsave(dat.in, file="dat.in.qs");    
}


fast.write <- function(dat, file, row.names=TRUE, quote="auto"){

    tryCatch(
        {
            if(is.data.frame(dat)){
                data.table::fwrite(dat, file, row.names=row.names);
            }else{
                write.csv(dat, file, row.names=row.names); 
            }
        }, error=function(e){
            print(e);
            write.csv(dat, file, row.names=row.names);   
        }, warning=function(w){
            print(w);
            write.csv(dat, file, row.names=row.names); 
        });
            
}

rowcolFt =  function(x, fac, var.equal, which = 1L) {
  
  if(!(which %in% c(1L, 2L)))
    stop(sQuote("which"), " must be 1L or 2L.")
  
  if(which==2L)
    x = t(x)

  if (typeof(x) == "integer")
      x[] <- as.numeric(x)

  sqr = function(x) x*x
  
  stopifnot(length(fac)==ncol(x), is.factor(fac), is.matrix(x))
  x   <- x[,!is.na(fac), drop=FALSE]
  fac <- fac[!is.na(fac)]

  ## Number of levels (groups)
  k <- nlevels(fac)

  ## xm: a nrow(x) x nlevels(fac) matrix with the means of each factor
  ## level
  xm <- matrix(
     sapply(levels(fac), function(fl) rowMeans(x[,which(fac==fl), drop=FALSE])),
     nrow = nrow(x),
     ncol = nlevels(fac))

  ## x1: a matrix of group means, with as many rows as x, columns correspond to groups 
  x1 <- xm[,fac, drop=FALSE]

  ## degree of freedom 1
  dff    <- k - 1

  if(var.equal){
    ## x0: a matrix of same size as x with overall means
    x0 <- matrix(rowMeans(x), ncol=ncol(x), nrow=nrow(x))
  
    ## degree of freedom 2
    dfr    <- ncol(x) - dff - 1

    ## mean sum of squares
    mssf   <- rowSums(sqr(x1 - x0)) / dff
    mssr   <- rowSums(sqr( x - x1)) / dfr

    ## F statistic
    fstat  <- mssf/mssr

  } else{

    ## a nrow(x) x nlevels(fac) matrix with the group size  of each factor
    ## level
    ni <- t(matrix(tapply(fac,fac,length),ncol=nrow(x),nrow=k))

    ## wi: a nrow(x) x nlevels(fac) matrix with the variance * group size of each factor
    ## level
    sss <- sqr(x-x1)
    x5 <- matrix(
       sapply(levels(fac), function(fl) rowSums(sss[,which(fac==fl), drop=FALSE])),
       nrow = nrow(sss),
       ncol = nlevels(fac))          
    wi <- ni*(ni-1) /x5

    ## u : Sum of wi
    u  <- rowSums(wi)

    ## F statistic
    MR <- rowSums(sqr((1 - wi/u)) * 1/(ni-1))*1/(sqr(k)-1)
    fsno <- 1/dff * rowSums(sqr(xm - rowSums(wi*xm)/u) * wi)
    fsdeno <- 1+ 2* (k-2)*MR
    fstat <- fsno/fsdeno

    ## degree of freedom 2: Vector with length nrow(x)
    dfr <- 1/(3 * MR)
  
  }
  
  res = data.frame(statistic = fstat,
                   p.value   = pf(fstat, dff, dfr, lower.tail=FALSE),
                   row.names = rownames(x))

  attr(res, "df") = c(dff=dff, dfr=dfr)
  return(res)
}

rowcoltt =  function(x, fac, tstatOnly, which, na.rm) {
    
  if(.on.public.web){
    dyn.load(.getDynLoadPath());
  }

  if (!missing(tstatOnly) && (!is.logical(tstatOnly) || is.na(tstatOnly)))
      stop(sQuote("tstatOnly"), " must be TRUE or FALSE.")
  
  f = checkfac(fac)
  if ((f$nrgrp > 2) || (f$nrgrp <= 0))
    stop("Number of groups is ", f$nrgrp, ", but must be >0 and <=2 for 'rowttests'.")

  if (typeof(x) == "integer")
      x[] <- as.numeric(x)

  cc = .Call("rowcolttests", x, f$fac, f$nrgrp, which-1L, na.rm)
    
  res = data.frame(statistic = cc$statistic,
                   dm        = cc$dm,
                   row.names = dimnames(x)[[which]])

  if (!tstatOnly)
    res = cbind(res, p.value = 2*pt(abs(res$statistic), cc$df, lower.tail=FALSE))

  attr(res, "df") = cc$df    
  return(res)
}

checkfac = function(fac) {

  if(is.numeric(fac)) {
    nrgrp = as.integer(max(fac, na.rm=TRUE)+1)
    fac   = as.integer(fac)
  }
  ## this must precede the factor test
  if(is.character(fac))
    fac = factor(fac)

  if (is.factor(fac)) {
    nrgrp = nlevels(fac)
    fac   = as.integer(as.integer(fac)-1)
  } 
  if(!is.integer(fac))
    stop("'fac' must be factor, character, numeric, or integer.")
  
  if(any(fac<0, na.rm=TRUE))
    stop("'fac' must not be negative.")
    
  return(list(fac=fac, nrgrp=nrgrp))
}

.getDynLoadPath <- function() {
    path = "../../rscripts/networkanalystr/src/MicrobiomeAnalyst.so";
    return(path)
}

# obtain a numeric matrix, exclude comments if any
.to.numeric.mat <- function(dat1){
  # now remove all comments in dat1
  # assign rownames after covert to matrix as data.frame does not allow duplicate names
  comments.inx <- grep("^#", dat1[,1]);
  if(sum(comments.inx) > 0){
    row.nms <- dat1[-comments.inx,1];
    dat1 <- dat1[-comments.inx,-1];
  }else{
    row.nms <- dat1[,1];
    dat1 <- dat1[,-1];
  }
  dimensions <- dim(dat1)
  col.nms <- colnames(dat1)
  dat1 <- sapply(dat1, as.numeric);
  dat1 <- matrix(data=dat1, ncol=dimensions[2], nrow=dimensions[1])
  rownames(dat1) <- row.nms;
  colnames(dat1) <- col.nms;
  return(dat1);
}