#'Function to merge user data with public data.
#'@description This function is used to merge microbiome data.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PrepareMergedData <- function(mbSetObj, metadata, keepfeat){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  set.seed(1315);
  metadata <<- metadata;
  data<- qs::qread("merged.data.qs");
  if(mbSetObj$module.type == "meta"){
      data <- subsetPhyloseqByDataset(mbSetObj, data);
  }

  if(keepfeat!="All features"){ # trim data
    feat_no<-round(ntaxa(data)*0.20);
    keepfeat<-" top 20 percent most abundant features";
    myTaxa <- names(sort(taxa_sums(data), decreasing = TRUE)[1:feat_no]);
    data<-prune_taxa(myTaxa,data);
  }
    
  otu.tab<-otu_table(data,taxa_are_rows =TRUE);
  taxa_names(otu.tab)->tax_nm;
  userrefdata <<-merge_phyloseq(data);
  sampleusrref<<-sample_data(data);
  feat.perc<<-keepfeat;
  
  if(.on.public.web){
    .set.mbSetObj(mbSetObj)
    return(rank_names(userrefdata));
  }else{
    return(.set.mbSetObj(mbSetObj))
  }
}

#'Function to perform reference data mapping.
#'@description This function performs reference data mapping.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import phyloseq
PerformRefDataMapping <- function(mbSetObj, refdataNm, taxo_type, sample_var, biome){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  #reading the reference data: (OTU abundance and tax info) and associated sample data file seperately.
  refdatafile_otu <- paste("/",refdataNm, "/", refdataNm, "_otu.rds", sep="");
  refdatafile_tax <- paste("/",refdataNm, "/", refdataNm,"_tax.rds", sep="");
  
  current.refset.otu <- .read.microbiomeanalyst.lib.rds(refdatafile_otu, "ppd", refdataNm);
  current.refset.tax <- .read.microbiomeanalyst.lib.rds(refdatafile_tax, "ppd", refdataNm);
  
  # create phyloseq or phyloslim object
  
  if(.on.public.web){
    OTU <- otu_table(current.refset.otu, taxa_are_rows = TRUE)
    TAX <- tax_table(current.refset.tax)
    current.refset <- phyloseq(OTU, TAX)
  }else{ # phyloseq should be loaded locally
    OTU <- otu_table(current.refset.otu, taxa_are_rows = TRUE)
    TAX <- tax_table(current.refset.tax)
    current.refset <- phyloseq(OTU, TAX)
  }
  
  #sample data
  refsmpldataNm <- paste("/", refdataNm, "/", refdataNm ,"_sampledata.csv", sep="");
  
  if(.on.public.web){
    refsmpldataloc <- paste(rpath, "libs/ppd",refsmpldataNm,sep="");
  }else{
    refsmpldataloc <- paste("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/libs/ppd", refsmpldataNm, sep="");
  }
  
  current.sample <- read.csv(refsmpldataloc,sep = "\t",header = T,row.names = 1);
  
  #reading user data
  data <- mbSetObj$dataSet$proc.phyobj;
  otu_no <- ntaxa(data);
  
  if(sample_var == "null"){
    sample_var <- colnames(sample_data(data))[1];
  }
  
  #since we use modified name for each OTU(to make it unique and convenient);in order to do mapping we need complete and original mapping label.
  #since our reference data has Greengenes OTU Ids as taxonomy identifier, so if data is already in this format, we will directly merge both of them(reference and users data;kepping all OTU's unique to user data)
  #if taxo_type is SILVA; first we have to map it with Greengenes OTU Ids from mapping file;#then we have to compare Greengenes OTU Ids OTUs between user and reference data
  if(taxo_type=="SILVA"){
    taxa_names(data) <- mbSetObj$dataSet$comp_taxnm;
    a <- taxa_names(data);
    #taxonomy mapping file
    
    if(.on.public.web){
      lib.path <- get.fun.lib.path("picrust12");
      otu.dic <<- qs::qread(paste0(picrust_path, "/greengenes_taxmap.qs"));
    }else{
      otu.dic <<- readRDS("https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/libs/picrust/greengenes_taxmap.rds");
    }
    
    #returns a vector of the positions of (first) matches of its first argument in user data(second argument).
    match_ind <- match(unique(otu.dic[ ,taxo_type]),a);
    #getting all the indices of user data that match(removing NA)
    match_ind1 <- match_ind[!is.na(match_ind)];
    
    if(length(match_ind1)==0){
      AddErrMsg("No SILVA labels match between user and taxonomy mapping file. Please check taxonomy labels ");
      return(0);
    }
    
    #getting only matched taxa in user data and then have to prune other taxa;
    match_taxa <- a[match_ind1];
    #now replacing the taxonomy names with Greengenes OTU Idss
    match_taxagg <- otu.dic[match(match_taxa,otu.dic[ ,taxo_type]),1];
    match_taxagg <- as.character(match_taxagg);
    
    #getting only matched taxa in user data and then have to prune other taxa;followed by replacing it with Greengenes OTU Ids;
    data = prune_taxa(match_taxa,data);
    taxa_names(data) <- match_taxagg;
  }
  
  # have to check whether 20% of OTUs common in user and reference data or not;
  taxa_ind <- match(taxa_names(data),taxa_names(current.refset));
  
  if(length(taxa_names(data)) > length(taxa_names(current.refset))){
    otu_nums <- length(taxa_names(current.refset));
  }else{
    otu_nums <- length(taxa_names(data));
  }
  
  if(length(which(taxa_ind != "NA")) < 0.20*otu_nums){
    err.vec <<- "";
    AddErrMsg(paste0(c("Less than 20 percent OTU match between user and reference data. ", length(which(taxa_ind != "NA")), "% match!"), collapse=" "));
    AddErrMsg("Please make sure that your uploaded features can be converted to Taxonomy ID (i.e. 182129, 4479944)");
    return(0);
  } else {
    msg <- paste("Dataset from",biome,"has been selected for comparison with user's data.")
  }
  
  #pruning reference dataset
  current.refset = prune_taxa(taxa_names(data), current.refset);
  
  if(taxo_type=="SILVA"){
    #for SILVA id we can compare only between common id that matches between user and reference data;pruning others OTUs in reference data.
    data <- prune_taxa(taxa_names(current.refset),data);
  }
  
  #merging both the user and reference phyloslim objects(sample data will be merged seperately)
  #storing taxonomic rank for users data
  #all reference data have kingdom column which can be removed;
  colnames(tax_table(current.refset))[c(1:7)] <- c("Kingdom","Phylum","Class","Order", "Family", "Genus","Species");
  tax_table(current.refset) <- tax_table(current.refset)[,-1];
  userdatarank <<- rank_names(current.refset);
  current.ref_userdata <- merge_phyloseq(otu_table(data),otu_table(current.refset),tax_table(data),tax_table(current.refset));
  #dummy variable for showing different shape for user and reference data
  sample_data(data) <-  sample_data(data)[,c(sample_var,"sample_id")];
  if(!"dataset" %in% colnames(sample_data(data)) ){
    sample_data(data) <- cbind(sample_data(data), dataset=rep(mbSetObj$dataSet$name,nrow(sample_data(data))));
  }
  current.sample$sample_id <- rownames(current.sample);
  current.sample$dataset <- rep("reference",nrow(current.sample));
  sam_data <- as.data.frame(sample_data(data),check.names=FALSE);
  #selecting primary variable data
  #user_sam<-sam_data[sample_var];
  #making the name of same variable same for reference sample data and then merging them
  colnames(current.sample) <- colnames(sam_data);
  
  current.ref_usersamdata <- rbind(current.sample,sam_data);
  current.ref_usersamdata$dataset <- as.factor(current.ref_usersamdata$dataset);    
  current.ref_usersamdata <- sample_data(current.ref_usersamdata);
  #storing the taxonomy rank from users data
  #final phyloslim object;(taxonomy table after merging will get distorted if taxa_ranks are not same in both user and reference data;but it's of no use)
  merged.data <- merge_phyloseq(current.ref_userdata,current.ref_usersamdata);
  sample_data(merged.data)$sample_id <- rownames(sample_data(merged.data));
  if(mbSetObj$module.type == "meta"){
    qs::qsave(merged.data, "merged.data.raw.qs");
  }
  
  #data filteration and transformation
  merged.data <- transform_sample_counts(merged.data, function(x) x / sum(x) );
  
  qs::qsave(merged.data, "merged.data.qs");
  
  current.msg <<- paste(msg, collapse=".");
  mbSetObj$dataSet$lib.msg <- current.msg;
  
  return(.set.mbSetObj(mbSetObj));
}

#'Function to create PCoA
#'@description This function creates data for plotting PCoA.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
#'@import vegan
PCoA3DAnal.16SRef <- function(mbSetObj, barplotNm, ordMeth, distName, taxrank, metadata, format="png", dpi=72){
  
  mbSetObj <- .get.mbSetObj(mbSetObj);

  suppressMessages(library(vegan));

  data <- userrefdata;

  if(taxrank!="OTU"){
    data <- fast_tax_glom_mem(data, taxrank);
    if(is.null(data)){
        AddErrMsg("Errors in projecting to the selected taxanomy level!");
        return(0);
    }
    nm <- as.character(tax_table(data)[,taxrank]);
    if(sum(is.na(nm))/length(nm) > 0.7){
      AddErrMsg("More than 70% values are missing at this level!");
      return(0);
    }
    #converting NA values to unassigned
    nm[is.na(nm)] <- "Not_Assigned";
    data1 <- as.matrix(otu_table(data));
    rownames(data1) <- nm;
    #all NA club together
    data1 <- rowsum(as.matrix(data1), rownames(data1));
  }
  
  GP.ord <- ordinate(data, ordMeth, distName);
  #creating 2D image for Report Generation
  barplotNm = paste(barplotNm, ".", format, sep="");
  mbSetObj$imgSet$ppd.2d<-barplotNm;
  
  Cairo::Cairo(file=barplotNm, width=720, height=500, type=format, bg="white",dpi=dpi);
  box = plot_ordination(data,GP.ord,color=metadata,shape="data");
  box$layers <- box$layers[-1];
  box=box+geom_point(size =4,alpha=0.8)+theme_bw();
  # used for area color for ellipse
  sam_data<-sample_data(data);
  clsLbl <- sam_data[[metadata]];
  box=box+ stat_ellipse(type="norm", linetype=2, geom = "polygon",alpha = 0.2, aes_string(fill = quo(clsLbl)), show.legend=FALSE);
  print(box);
  dev.off();
  
  # obtain variance explained
  sum.pca<-GP.ord;
  imp.pca<-sum.pca$values;
  imp.pca<-data.frame(imp.pca,check.names=FALSE);
  std.pca<-imp.pca[1,]; # eigen values
  var.pca<-imp.pca[,2]; # variance explained by each PC
  cum.pca<-imp.pca[5,]; # cummulated variance explained
  mbSetObj$analSet$topo.msg<-paste("Beta-diversity is performed using",feat.perc,"and",distName,"distance measure");
  mbSetObj$analSet$sum.pca<-append(sum.pca, list(std=std.pca, variance=var.pca, cum.var=cum.pca));
  mbSetObj$analSet$vis$feat.perc <- feat.perc
  mbSetObj$analSet$vis$distName <- distName
  return(.set.mbSetObj(mbSetObj))
}

#'Function to plot 3D score plot.
#'@description This function plots a 3D PCoA score plot.
#'@param mbSetObj Input the name of the mbSetObj.
#'@author Jeff Xia \email{jeff.xia@mcgill.ca}
#'McGill University, Canada
#'License: GNU GPL (>= 2)
#'@export
PlotUsrRefPCoA3DScore <- function(mbSetObj, imgName, format="json", inx1, inx2, inx3, variable){

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  pca <- mbSetObj$analSet$sum.pca;
  pca3d <- list();
  pca3d$score$axis <- paste("PC", 1:3, " (", 100*round(pca$variance[1:3], 3), "%)", sep="");
  coords <- data.frame(t(signif(pca$vectors[,c(inx1, inx2, inx3)], 5)),check.names=FALSE);
  colnames(coords) <- NULL;
  pca3d$score$xyz <- coords;
  pca3d$score$name <- sample_names(userrefdata);
  sam_data<-sample_data(userrefdata);
  cls<-as.character(sam_data[[variable]]);
  clsb<-as.character(sam_data[["dataset"]]);
  pca3d$score$facA <- cls;
  pca3d$score$facB <- clsb;
  variable <<- variable;
    
  # now set color for each group
  grp.num <- length(unique(cls)) + 1;
  cols <- 1:grp.num + 1;
  rgbcols <- col2rgb(cols);
  cols <- apply(rgbcols, 2, function(x){paste("rgb(", paste(x, collapse=","), ")", sep="")});
  pca3d$score$colors <- cols;

  mypos <- t(coords);
  colnames(mypos) <- paste("Dim", 1:3, sep="");

  pos.xyz <- mypos;
  rownames(pos.xyz) = pca3d$score$name;
  qs::qsave(pos.xyz, "score_pos_xyz.qs");

  json.obj <- rjson::toJSON(pca3d);
  #sink(imgName);
  #cat(json.obj);
  #sink();

  mbSetObj$analSet$pca <- pca3d;

  if(!exists("my.json.scatter")){
    .load.scripts.on.demand("utils_scatter3d.Rc");    
  }

  .set.mbSetObj(mbSetObj)
  my.json.scatter(mbSetObj, imgName, F);

  return(.set.mbSetObj(mbSetObj))
}
