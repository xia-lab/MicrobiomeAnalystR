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
  
  if(.on.public.web){
    load_pryr();
  }
   sink(); # make sure print to screen
  print(mem_used());
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
  
  if(.on.public.web){
    load_rcolorbrewer();
  }

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
  msg.vec <<- c(msg.vec, msg);
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

fast_tax_glom_first <- function(physeq, taxrank){

  # setup data. we are going to glob the OTU table based on the Class Taxonomy
  CN  <- which( rank_names(physeq) %in% taxrank);
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
    
  if(length(phy_tree(physeq,errorIfNULL = FALSE))==0){
    phy_data<-merge_phyloslim(otab2,tax_table(physeq),sample_data(physeq));
  }else{
    phy_data<-merge_phyloslim(otab2,tax_table(physeq),sample_data(physeq),phy_tree(physeq));
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

# need to return consisten color assignments for the same taxa
GetSerieColors <- function(taxa=NULL){
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
  
  if(.on.public.web){
    load_xtable();
  }

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
