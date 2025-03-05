##################################################
## R script for MicrobiomeAnalyst
## Functions for SparCC
##################################################

my.fast.spar <- function(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, output="network", opt="corr"){
  
  # only works for unix
  my.os <- tolower(GetOS());

  if(my.os == "mac" | my.os == "windows"){
    AddErrMsg("The built-in fastspar only works for unix system. Please visit https://github.com/scwatts/fastspar for more details.");
    return(0);
  }
  
  if(.on.public.web){
    spar_home <- paste0(rpath, "lib/fastspar/");
  }else if (file.exists("/home/jasmine/Downloads/fastspar/fastspar")){ #jas local
    spar_home <- "/home/jasmine/Downloads/fastspar/"
  }else{
    spar_home <- "https://www.microbiomeanalyst.ca/MicrobiomeAnalyst/resources/lib/fast_spar/"
  }
  
  path_fastspar <- paste(spar_home, "fastspar", sep="");
  path_fastspar_bs <- paste(spar_home, "fastspar_bootstrap", sep="");
  path_fastspar_pvals <- paste(spar_home, "fastspar_pvalues", sep="");
  
  # make sure they are executable
  my.cmd <- paste("chmod a+x", path_fastspar, path_fastspar_bs, path_fastspar_pvals);
system(my.cmd);
  
  mbSetObj <- .get.mbSetObj(mbSetObj);
  if(.on.public.web){
    load_igraph();
  }
  
  if(opt == "pattern"){
    data1 <- qs::qread("pattern_data.qs")
  }else{
    # first get OTU table to BIOM TSV format
    if(mbSetObj$module.type=="sdp"){
      data1 <- as.matrix(otu_table(mbSetObj$dataSet$norm.phyobj));
      taxrank <- "OTU";
    }else if(mbSetObj$module.type=="mdp"){
      if(!exists("phyloseq_objs")){
        phyloseq_objs <- qs::qread("phyloseq_objs.qs")
      }
      
      if(taxrank=="OTU"){
        data1 <- phyloseq_objs$count_tables$OTU
      }else{
        taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
        data1 <- phyloseq_objs$count_tables[[taxrank.inx]]
      }
    }
  }
  
  if(opt=="corr"){
    qs::qsave(data1, "sparcc_full_data.qs")  
  }
  
  if(!is.null(mbSetObj$dataSet$selected.grps) & opt == "corr"){
    data1 <- data1[,which(colnames(data1) %in% mbSetObj$dataSet$selected.grps)]
  }
  
  # filter data if necessary
  abund_filtered <- mbSetObj$dataSet$ab.filtered
  var_filtered <- mbSetObj$dataSet$var.filtered
  
  feat_num <- nrow(data1)
  
  # if number of feats is greater than 500, do filtering of poorly represented OTUs
  # first, do low count filtering 
  if(feat_num > 500 & !abund_filtered){
    minLen <- 0.2*ncol(data1); # filter out feats that do not have a minimum count of 4 in at least 20 percent of samples
    kept.inx <- apply(data1, MARGIN = 1,function(x) {sum(x >= 4) >= minLen});
    data1 <- data1[kept.inx, ];
  }
  
  # second, do low variance filtering
  if(feat_num > 500 & !var_filtered){
    filter.val <- apply(data1, 1, IQR, na.rm=T);
    rk <- rank(-filter.val, ties.method='random');
    remain <- rk < feat_num*(1-0.1); # filter out 10% of low variance feats based on IQR
    data1 <- data1[remain,];
  }
  
  # third, if still over 500 feats, rank and keep top 500
  if(feat_num > 500){
    filter.val <- apply(data1, 1, IQR, na.rm=T);
    rk <- rank(-filter.val, ties.method='random');
    remain <- rk < 500;
    data1 <- data1[remain,];
    current.msg <<- "Only the top 500 features are kept, ranked by their variance!"
  }
  
  # replace zeros with random numbers based on lowest non-zero count
  
  set.seed(12345)
  
  if(opt %in% c("corr", "feat")){
    lowest.count <- min(apply(data1, 1, FUN = function(x) {min(x[x > 0])}))
  }else{ #only pattern searching
    lowest.count <- min(apply(data1[-nrow(data1),], 1, FUN = function(x) {min(x[x > 0])}))
  }
  
  # replacements
  replacements <- apply(data1, 1, function(x) { (sample(1:lowest.count, size=sum(x==0), replace = F))/lowest.count} )
  
  # this works for rows
  replaced <- vector(length = nrow(data1), "list")
  
  for(i in 1:nrow(data1)){
    replaced[[i]] <- replace(as.vector(data1[i,]), as.vector(data1[i,]==0), replacements[[i]])
  }
  
  zero.output <- matrix(unlist(replaced), ncol = ncol(data1), byrow = TRUE)
  colnames(zero.output) <- colnames(data1)
  rownames(zero.output) <- rownames(data1)
  
  if(opt=="corr"){
    qs::qsave(zero.output, "sparcc_data.qs")
  }
  
  # add header
  new_col <- c("#OTU ID", colnames(zero.output))
  data <- cbind(rownames(zero.output), zero.output)
  colnames(data) <- new_col
  
  write.table(data, file="otu_table_corr.tsv", quote=FALSE, row.names = FALSE, sep="\t")
  
  otu_table <- "otu_table_corr.tsv"
  corr_output <- "sparcc_median_correlation.tsv"
  cov_output <- "sparcc_median_covariance.tsv"
  bootstrap_counts <- "bootstrap_counts"
  counts_prefix <- "bootstrap_counts/boot_data"
  bootstrap_correlation <- "bootstrap_correlation"
  corr_bs_output <- "${jnew}_median_correlation.tsv"
  cov_bs_output <- "${k}_median_covariance.tsv"
  corr_prefix <- "cor_boot_data_"
  test <- "text.txt"
  pval_output <- "sparcc_pvals.tsv"
  
  invisible(system(paste(path_fastspar, "--otu_table", otu_table, "--correlation", corr_output, "--covariance", cov_output, ";",
                         "mkdir", bootstrap_counts, ";",
                         path_fastspar_bs,  "--otu_table", otu_table, "--number", permNum, "--prefix", counts_prefix, ";",
                         "find", bootstrap_counts, ">>", test, ";",
                         "cat", test, "| while read line; do echo $line; j=$(basename ${line}); jnew=$(echo cor_${j}); k=$(echo cov_${j}); jnew=$(echo ${jnew} | sed 's/\\.tsv//'); k=$(echo ${k} | sed 's/\\.tsv//');", path_fastspar, "--otu_table ${line} --correlation",
                         corr_bs_output, "--covariance", cov_bs_output, "-i 5 ; done;",
                         path_fastspar_pvals, "--otu_table", otu_table, "--correlation", corr_output, "--prefix", corr_prefix, "--permutations", permNum, "--outfile", pval_output, ";",
                         "rm cor_boot_data*.tsv; rm cov_boot_data*.tsv;", "rm -rf", bootstrap_counts, "rm -rf text.txt"), intern=TRUE, ignore.stdout=TRUE))
  
  sparcc_results <- read.table(file = "sparcc_median_correlation.tsv", sep="\t", stringsAsFactors = FALSE)
  names <- sparcc_results[,1]
  sparcc_results_new <- as.matrix(sparcc_results[,-1])
  colnames(sparcc_results_new) <- rownames(sparcc_results_new) <- names
  
  sparcc_pvals <- read.table(file = "sparcc_pvals.tsv", sep="\t", stringsAsFactors = FALSE)
  sparcc_pvals_new <- as.matrix(sparcc_pvals[,-1])
  colnames(sparcc_pvals_new) <- rownames(sparcc_pvals_new) <- names
  
  if(output == "network"){
    sparcc_results_new[sparcc_results_new==0] <- .00001
    sparcc_corr <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(sparcc_results_new, weighted = TRUE))
    sparcc_pvals_new[sparcc_pvals_new==0] <- .00001
    sparcc_pvals <- igraph::as_data_frame(igraph::graph_from_adjacency_matrix(sparcc_pvals_new, weighted = TRUE))
    
    sparcc_combo <- cbind(sparcc_corr, sparcc_pvals[,3])
    colnames(sparcc_combo) <- c("Taxon1", "Taxon2", "Correlation", "P.Value")
    qs::qsave(sparcc_combo, "network_correlation.qs")
    
    sparcc_combo <- sparcc_combo[(abs(sparcc_combo[,3]) > corrCutoff & sparcc_combo[,4] < pvalCutoff),]
    fast.write(sparcc_combo, "correlation_table.csv", row.names = FALSE)
    .set.mbSetObj(mbSetObj)
    return(sparcc_combo)
  }else if(output=="heatmap"){
    .set.mbSetObj(mbSetObj)
    return(sparcc_results_new)
  }
}
