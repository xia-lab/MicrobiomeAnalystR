##################################################
## R script for MicrobiomeAnalyst
## Functions to compute correlation network
## Authors:  Jeff Xia, jeff.xia@mcgill.ca; Yao Lu, yao.lu5@mail.mcgill.ca
##################################################

my.corr.net <- function(mbSetObj, taxrank, cor.method="pearson", colorOpt="expr", 
                                      permNum=100, pvalCutoff=0.05, corrCutoff=0.3, abundOpt="mean", 
                                      corr.net.name, plotNet = FALSE, netType="static", netLayout="kk",
                                      netTextSize = 2.5){

  #save.image("mycorr.RData");

  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataSet$cor.method <- cor.method;
  mbSetObj$analSet$abund.opt <- abundOpt;
  current.msg <<- ""; 

  suppressMessages(library(ppcor));
  suppressMessages(library(igraph));
  
  if(cor.method == "sparcc"){
    if(.on.public.web){
      sparcc_results <- RunFastSpar_mem(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network")
    }else{
      sparcc_results <- RunFastSpar(mbSetObj, taxrank, permNum, pvalCutoff, corrCutoff, "network")
    }
    
    if(nrow(sparcc_results)==0){
      AddErrMsg("No correlations meet the p-value and correlation thresholds!")
      return(0)
    }else{
      mbSetObj$analSet$network_cor <- sparcc_results
    }
  }else if(cor.method == "secom_p1"){
  
  if(!exists("my.secom.anal")){ # public web on same user dir
    .load.scripts.on.demand("utils_secom.Rc");    
  }
 
   secomp1_results <- my.secom.anal(mbSetObj,taxrank,permNum,corrCutoff, pvalCutoff,"secomp1")
    # Check if analysis failed (returns 0) or has no results
    if(!is.data.frame(secomp1_results) || nrow(secomp1_results)==0){
      if(is.data.frame(secomp1_results)){
        AddErrMsg("No correlations meet the p-value and correlation thresholds!")
      }
      return(0)
    }else{
      mbSetObj$analSet$network_cor <- secomp1_results
    }

}else if(cor.method == "secom_p2"){
  
  if(!exists("my.secom.anal")){ # public web on same user dir
    .load.scripts.on.demand("utils_secom.Rc");    
  }

   secomp2_results <- my.secom.anal(mbSetObj,taxrank,permNum,corrCutoff, pvalCutoff,"secomp2")
   # Check if analysis failed (returns 0) or has no results
   if(!is.data.frame(secomp2_results) || nrow(secomp2_results)==0){
      if(is.data.frame(secomp2_results)){
        AddErrMsg("No correlations meet the p-value and correlation thresholds!")
      }
      return(0)
    }else{
      mbSetObj$analSet$network_cor <- secomp2_results
    }

}else if(cor.method == "secom_dist"){
  
  if(!exists("my.secom.anal")){ # public web on same user dir
    .load.scripts.on.demand("utils_secom.Rc");    
  }

   secomdis_results <- my.secom.anal(mbSetObj,taxrank,permNum,corrCutoff, pvalCutoff,"secomdis")
   # Check if analysis failed (returns 0) or has no results
   if(!is.data.frame(secomdis_results) || nrow(secomdis_results)==0){
      if(is.data.frame(secomdis_results)){
        AddErrMsg("No correlations meet the p-value and correlation thresholds!")
      }
      return(0)
    }else{
      mbSetObj$analSet$network_cor <- secomdis_results
    }

}else{
    if(!exists("phyloseq_objs")){
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
    }
 
    if(taxrank=="OTU"){
      data1 <- phyloseq_objs$count_tables$OTU
    }else{
      taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
      data1 <- phyloseq_objs$count_tables[[taxrank.inx]]
    }
  
    data <- t(data1);
    shadow_save(data, "network_cor_data.qs")
 
    data <- data[which(rownames(data) %in% mbSetObj$dataSet$selected.grps),]
   # data[data==0|is.na(data)] <- .00001

    if(ncol(data) > 1000){
      filter.val <- apply(data.matrix(data), 2, IQR, na.rm=T);
      rk <- rank(-filter.val, ties.method='random');
      data <- as.data.frame(data[,rk <=1000],check.names=FALSE);
    }
 
    mbSetObj$analSet$netcorr_data <- data;
    vars = data.frame(t(combn(colnames(data), 2)), stringsAsFactors = FALSE,check.names=FALSE)
    
    if(cor.method == "spearman"){
      cor.results <- do.call(rbind, mapply(SpearmanCorrFunc, vars[,1], vars[,2], MoreArgs = list(data=data), SIMPLIFY = FALSE))
    }else if(cor.method == "kendall"){
      cor.results <- do.call(rbind, mapply(KendallCorrFunc, vars[,1], vars[,2], MoreArgs = list(data=data), SIMPLIFY = FALSE))
    }else if(cor.method == "pearson"){
      cor.results <- do.call(rbind, mapply(PearsonCorrFunc, vars[,1], vars[,2], MoreArgs = list(data=data), SIMPLIFY = FALSE))
    }else{
      AddErrMsg("Invalid correlation method!")
      return(0)
    } 
    colnames(cor.results) <- c("Taxon1", "Taxon2", "Correlation", "P.value", "Statistic", "Method")
    shadow_save(cor.results, "network_correlation.qs")
    
    cor.results.filt <- cor.results[(abs(cor.results[,3]) > corrCutoff & cor.results[,4] < pvalCutoff),]
    cor.results.filt[,3] <- round(cor.results.filt[,3], digits=4)
    cor.results.filt[,4] <- round(cor.results.filt[,4], digits=4)
    mbSetObj$analSet$network_cor <- cor.results.filt;
  }
   fast.write(mbSetObj$analSet$network_cor, "correlation_table.csv", row.names = FALSE)
  mbSetObj$dataSet$corr.pval.cutoff <- pvalCutoff
  mbSetObj$analSet$corr_color_opt <- colorOpt;
 
  #network building only needed for web

  if(.on.public.web || plotNet == TRUE){
    all.taxons = unique(c(mbSetObj$analSet$network_cor[,1] , mbSetObj$analSet$network_cor[,2]))
    taxColNms = vector();
    
    if(taxrank == "OTU"){
      tax.tbl = as.data.frame(matrix(mbSetObj$dataSet$taxa_table[,1:7],ncol=7),check.names=FALSE)
      colnames(tax.tbl) = colnames(mbSetObj$dataSet$taxa_table[,1:7])
      colnames(tax.tbl)[7] = "OTU"
      tax.tbl[,7]=rownames(mbSetObj$dataSet$taxa_table)
      taxColNms = colnames(tax.tbl)
      tax.tbl = as.data.frame(tax.tbl[which(tax.tbl[,simpleCap(taxrank)] %in% all.taxons),],check.names=FALSE)
      colnames(tax.tbl) = taxColNms
    }else{
      tax.tbl = data.frame(mbSetObj$dataSet$taxa_table[,1:which( colnames(mbSetObj$dataSet$taxa_table)== simpleCap(taxrank))],check.names=FALSE)
      taxColNms = colnames(mbSetObj$dataSet$taxa_table[,1:which( colnames(mbSetObj$dataSet$taxa_table)== simpleCap(taxrank))])
        colnames(tax.tbl) = taxColNms
      tax.tbl = as.data.frame(tax.tbl[which(tax.tbl[,simpleCap(taxrank)] %in% all.taxons),],check.names=FALSE)
      colnames(tax.tbl) = taxColNms
 
    }
    
    inx = !duplicated(tax.tbl[,simpleCap(taxrank)])
    filt.taxa.table = data.frame(tax.tbl[inx,],check.names=FALSE)
    colnames(filt.taxa.table) = taxColNms 
    mbSetObj$analSet$filt.taxa.table = filt.taxa.table 
    .set.mbSetObj(mbSetObj);
    res <- SparccToNet(mbSetObj, corr.net.name, netType, netLayout, netTextSize);
    if(.on.public.web){
      if(res == 0){
        AddErrMsg(paste0("Errors during creating correlation network:", current.msg))
        return(0);
      }
    }else{
      if(length(res) == 1){
        AddErrMsg(paste0("Errors during creating correlation network:", current.msg))
        return(0);
      }
    }
  }
  
  if(cor.method=="sparcc"){
    method = "SparCC"
  }else if(cor.method=="secom_p1"){
    method = "SECOM (Pearson1)"; 
  }else if(cor.method=="secom_p2"){
    method = "SECOM (Pearson2)"; 
  }else if(cor.method=="secom_dist"){
    method = "SECOM (Distance)";  
  }else{
    method = paste(cor.method, "correlation"); 
  }
  
  # calculate microbial dysbiosis index
  feats_used <- unique(mbSetObj$analSet$network_cor[,1], mbSetObj$analSet$network_cor[,2]) 
  fc <- mbSetObj$analSet$diff_table
  clean_feats_used <- sub("^[^_]*__", "", unique(feats_used))
  match.inx <- which(fc$tax_name %in% clean_feats_used)
 
  if(cor.method=="sparcc"){
    abund_data <- qs::qread("sparcc_data.qs")
  }else if(cor.method %in% c("secom_p1","secom_p2","secom_dist")){
   abund_data <- qs::qread("secom_data.qs")
   }else{
    abund_data <- t(as.matrix(mbSetObj$analSet$netcorr_data))
  }
 
  abund.inx <- which(rownames(abund_data) %in% feats_used)
  abund_data.filt <- abund_data[abund.inx,]
  smpl <- data.frame(sample_data(mbSetObj$dataSet$proc.phyobj),check.names=FALSE);
 
  if(length(match.inx)!=0){
    
    fc.filt <- fc[match.inx,]
    
    if(length(mbSetObj$dataSet$comparison) > 2){ # only calculate mdi for 2 groups
      
      mbSetObj$analSet$mdi <- "MD-index: NA (only defined for two-group comparison)";
      
    }else{ # for only 2 groups
      if(length(unique(fc.filt$tax_name)) < length(fc.filt$tax_name)){
        # need to delete doubles in fc.filt
        keep <- c(TRUE, head(fc.filt$tax_name, -1) != tail(fc.filt$tax_name, -1))
        fc.filt <- fc.filt[keep,]
      }
      
      smpls.used <- rownames(smpl[smpl[mbSetObj$dataSet$meta]==fc.filt$treatment_1,,drop=FALSE])
      abund_data.filt <- abund_data.filt[,match(smpls.used, colnames(abund_data.filt))]
      
      log2fc <- fc.filt$log2_median_ratio
      names(log2fc) <- fc.filt$tax_name
      
      inc <- which(log2fc > 0) # total abundance increased in X
      dec <- which(log2fc < 0) # total abundance decreased in X
      
      if(length(inc)==0){
        mbSetObj$analSet$mdi <- "MD-index could not be calculated - no enriched taxa were increased.";
      }else if(length(dec)==0){
        mbSetObj$analSet$mdi <- "MD-index could not be calculated - no depleted taxa were detected."
      }else{
        abund <- c(sum(abund_data.filt[inc,]), sum(abund_data.filt[dec,]))
        abund[abund==0] <- 0.1
        
        MDI1 <- round(log(abund[1]/abund[2]), digits = 4) 
        first <- unique(fc.filt$treatment_1)
        second <- unique(fc.filt$treatment_2)
        groups <- paste(first, "/", second, sep="")
        
        mbSetObj$analSet$mdi <- paste(groups, "MD-index:", MDI1)
      }
    }
  }

  mbSetObj$analSet$cornet.taxa = taxrank;

  current.msg <<- paste(current.msg, method, "network analysis performed successfully!")
 
  return(.set.mbSetObj(mbSetObj));
}

