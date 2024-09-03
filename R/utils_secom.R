##################################################
## R script for MicrobiomeAnalyst
## Description:  
## Authors: Jeff Xia, jeff.xia@mcgill.ca; Yao Lu, yao.lu5@mail.mcgill.ca
###################################################

my.secom.anal<-function(mbSetObj,taxrank,R,corr_cut, max_p,mode,method='pearson',n_cl=2) {
  #print(corr_cut)
  #print(max_p)
  #print(mode)
  mbSetObj <- .get.mbSetObj(mbSetObj);
   if(!exists("phyloseq_prenorm_objs")){
     phyloseq_prenorm_objs <- qs::qread("phyloseq_prenorm_objs.qs")
   }
  
  feature_table <- phyloseq_prenorm_objs$count_tables[[taxrank]]
  meta_data <- mbSetObj$dataSet$sample_data
  if(!is.null(mbSetObj$dataSet$selected.grps)){
  feature_table <- feature_table[,which(colnames(feature_table) %in% mbSetObj$dataSet$selected.grps)]
  meta_data <- meta_data[rownames(meta_data) %in% mbSetObj$dataSet$selected.grps,]
}

  abn_list = abn_est(feature_table,meta_data, taxrank, pseudo=0,tax_keep=NULL)

  s_diff_hat = abn_list$s_diff_hat
  y_hat = abn_list$y_hat
  
  cl = parallel::makeCluster(n_cl)
  doParallel::registerDoParallel(cl)
 require("doRNG")

  if(mode=='secomdis'){
    res_corr = sparse_dist(mat = t(y_hat), wins_quant=c(0.05, 0.95), R=R, thresh_hard=0, max_p=max_p)
    parallel::stopCluster(cl)
    
    # To prevent FP from taxa with extremely small variances
      corr_s = cor(cbind(s_diff_hat, t(y_hat)),
                   use = "pairwise.complete.obs")[1, -1]
      fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
      fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
      fp_ind = (fp_ind1 * fp_ind2 == 1)
      diag(fp_ind) = FALSE
      res_corr$dcorr[fp_ind] = 0
      res_corr$dcorr_fl[fp_ind] = 0
      res_corr$dcorr_p[fp_ind] = 1
  
      res = reshape2::melt(res_corr$dcorr_fl)
      res$pval = reshape2::melt(res_corr$dcorr_p)$value
      res$method = "secom_dist"
  
  }else{
       
    if (method %in% c("pearson", "kendall", "spearman")) {
      res_corr = sparse_linear(mat = t(y_hat), wins_quant=c(0.05, 0.95), method, soft= FALSE,
                               thresh_len = 100, n_cv = 10, thresh_hard=0, max_p=max_p)
    } else {
      stop_txt = paste0("The correlation coefficien should be one of ",
                        "'pearson', 'kendall', 'spearman'")
      stop(stop_txt, call. = FALSE)
    }
    
    parallel::stopCluster(cl)
    
    # To prevent FP from taxa with extremely small variances
    corr_s = cor(cbind(s_diff_hat, t(y_hat)),
                 use = "pairwise.complete.obs")[1, -1]
    fp_ind1 = replicate(nrow(y_hat), corr_s > corr_cut)
    fp_ind2 = t(replicate(nrow(y_hat), corr_s > corr_cut))
    fp_ind = (fp_ind1 * fp_ind2 == 1)
    diag(fp_ind) = FALSE
    res_corr$corr[fp_ind] = 0
    res_corr$corr_th[fp_ind] = 0
    res_corr$corr_fl[fp_ind] = 0
    res_corr$corr_p[fp_ind] = 1
   
    
    if(mode=="secomp1"){
      res = reshape2::melt(res_corr$corr_th)
      res$pval = reshape2::melt(res_corr$corr_p)$value
      res$method = "secom_pearson1"
    }else if(mode=="secomp2"){
      res = reshape2::melt(res_corr$corr_fl)
      res$pval = reshape2::melt(res_corr$corr_p)$value
      res$method = "secom_pearson2"
    }
    
  }
   names(res) <- c("Taxon1", "Taxon2", "Correlation", "P.value", "Method")
   res<- res[which(res$Correlation !=0 & res$Taxon1 !=res$Taxon2),]
   res <- res[(abs(res[,3]) > corr_cut & res[,4] < max_p),]
 
   res[,3] <- round(res[,3], digits=4)
   res[,4] <- round(res[,4], digits=4)
   secom_data <- exp(1)^y_hat
   secom_data[is.na(secom_data)] <- 0
   qs::qsave(secom_data,"secom_data.qs")
   return(res);
  
}


###############function adpated from ANCOMBC package "https://github.com/FrederickHuangLin/ANCOMBC/"
# Bias-corrected abundance estimation
abn_est = function(feature_table,meta_data, tax_level, pseudo, prv_cut = 0.5, lib_cut = 1000,
                   tax_keep = NULL, samp_keep = NULL) {
  # Sampling fraction difference estimation

  O1 = feature_table
  if(nrow(feature_table)<50){
    curret.msg<<- paste0("The number of taxa used for estimating sample-specific biases is  ",nrow(feature_table),"  A large number of taxa (> 50) is required for the statistical consistency!")
  }
  s_diff_hat = s_diff_est(O1)
  samp_keep = names(s_diff_hat)
  # Discard taxa with prevalences < prv_cut
  
  prevalence = apply(feature_table, 1, function(x)
    sum(x != 0, na.rm = TRUE)/length(x[!is.na(x)]))
  tax_keep = which(prevalence >= prv_cut)
  

  if (length(tax_keep) > 0) {
    feature_table = feature_table[tax_keep, , drop = FALSE]
  } else {
    stop("No taxa remain under the current cutoff", call. = FALSE)
  }
  # Discard samples with library sizes < lib_cut
  if (is.null(samp_keep)) {
    lib_size = colSums(feature_table, na.rm = TRUE)
    samp_keep = which(lib_size >= lib_cut)
  }
  if (length(samp_keep) > 0){
    feature_table = feature_table[, samp_keep, drop = FALSE]
    meta_data = meta_data[samp_keep, , drop = FALSE]
  } else {
    stop("No samples remain under the current cutoff", call. = FALSE)
  }
  
 
  O2 = feature_table
  O2 = O2 + pseudo
  o = log(O2)
  o[is.infinite(o)] = NA
  n = ncol(o)
  d = nrow(o)
  
  # Bias-corrected abundance estimation
  y_hat = o - rowMeans(o, na.rm = TRUE) - t(replicate(d, s_diff_hat))
  
  abs_list = list(s_diff_hat = s_diff_hat, y_hat = y_hat)
  return(abs_list)
}


# Sampling fraction difference estimation
s_diff_est = function(feature_table) {
  if (nrow(feature_table) < 50) {
    warn_txt = sprintf(paste("The number of taxa used for estimating sample-specific biases is: ",
                             nrow(feature_table),
                             "A large number of taxa (> 50) is required for the statistical consistency",
                             sep = "\n"))
    warning(warn_txt, call. = FALSE)
  }
  o = log(feature_table)
  
  o[is.infinite(o)] = NA
  o_center = o - rowMeans(o, na.rm = TRUE)
  
  # Estimate weights
  wt = apply(o_center, 1, function(x) 1/var(x, na.rm = TRUE))
  o_center = o_center[is.finite(wt), ]
  wt = wt[is.finite(wt)]
  
  # Estimate sampling fraction difference
  s_diff_hat = apply(o_center, 2, function(x) {
    weighted.mean(x, wt, na.rm = TRUE)}
  )
  
  return(s_diff_hat)
}



sparse_linear = function(mat, wins_quant = c(0.05, 0.95), method, soft=FALSE, thresh_len=100,
                         n_cv=10, thresh_hard=0, max_p) {
  # Thresholding
  mat_thresh = function(mat, th, soft){
    mat_sign = sign(mat)
    mat_th = mat
    mat_th[abs(mat) <= th] = 0
    if (soft) {
      mat_th[abs(mat) > th] = abs(mat_th[abs(mat) > th]) - th
      mat_th = mat_th * mat_sign
    }
    return(mat_th)
  }
  
  # Threshold loss function
  thresh_loss = function(mat1, mat2, method, th, soft) {
    corr1 = cor(mat1, method = method, use = "pairwise.complete.obs")
    corr2 = cor(mat2, method = method, use = "pairwise.complete.obs")
    corr_diff = mat_thresh(corr1, th, soft) - corr2
    corr_diff[is.na(corr_diff)] = 0
    loss = norm(corr_diff, type = "F")
    return(loss)
  }
  
  # Filtering based on p-values
  p_filter = function(mat, mat_p, max_p){
    ind_p = mat_p
    ind_p[mat_p > max_p] = 0
    ind_p[mat_p <= max_p] = 1
    
    mat_filter = mat * ind_p
    return(mat_filter)
  }
  
  # Sort taxa
  sort_taxa = sort(colnames(mat))
  mat = mat[, sort_taxa]
  
  # Winsorization
  mat = apply(mat, 2, function(x)
    DescTools::Winsorize(x, val = quantile(x, probs = wins_quant, na.rm = TRUE)))
  
  # Co-occurrence matrix
  mat_occur = mat
  mat_occur[mat_occur != 0] = 1
  mat_occur[mat_occur == 0] = 0
  mat_occur[is.na(mat_occur)] = 0
  
  df_occur = as.data.frame(mat_occur) %>%
    tibble::rownames_to_column("sample_id") %>%
    tidyr::pivot_longer(cols = -.data$sample_id, names_to = "taxon",
                        values_to = "occur") %>%
    dplyr::filter(.data$occur == 1)
  
  mat_cooccur = crossprod(table(df_occur[, seq_len(2)]))
  diag(mat_cooccur) = colSums(mat_occur)
  
  if (any(mat_cooccur < 10)) {
    warn_txt = sprintf(paste("There are some pairs of taxa that have insufficient (< 10) overlapping samples",
                             "Proceed with caution since the point estimates for these pairs are unstable",
                             "For pairs of taxa with no overlapping samples, the point estimates will be replaced with 0s,",
                             "and the corresponding p-values will be replaced with 1s",
                             "Please check `mat_cooccur` for details about the co-occurrence pattern",
                             sep = "\n"))
    warning(warn_txt)
  }
  
  # Sample size for training and test sets
  n = dim(mat)[1]
  n1 = n - floor(n/log(n))
  n2 = n - n1
  d = dim(mat)[2]
  
  # Correlation matrix
  corr_list = suppressWarnings(Hmisc::rcorr(x = mat, type = method))
  corr = corr_list$r
  corr[mat_cooccur < 2] = 0
  corr[is.infinite(corr)] = 0
  
  # Cross-Validation
  max_thresh = max(abs(corr[corr != 1]), na.rm = TRUE)
  thresh_grid = seq(from = 0, to = max_thresh, length.out = thresh_len)

  loss_mat = foreach(i = seq_len(n_cv), .combine = rbind) %dorng% {
    index = sample(seq_len(n), size = n1, replace = FALSE)
    mat1 = mat[index,]
    mat2 = mat[-index,]
    loss = vapply(thresh_grid, FUN = thresh_loss,
                  mat1 = mat1, mat2 = mat2,
                  method = method, soft = soft,
                  FUN.VALUE = double(1))
  }
  
  # Correlation matrix after thresholding
  loss_vec = colMeans(loss_mat)
  thresh_opt = thresh_grid[which.min(loss_vec)]
  corr_th = mat_thresh(mat = corr, th = thresh_opt, soft = soft)
  corr_th = mat_thresh(mat = corr_th, th = thresh_hard, soft = FALSE)
  
  # Correlation matrix after filtering
  corr_p = corr_list$P
  diag(corr_p) = 0
  corr_p[mat_cooccur < 2] = 1
  corr_p[is.na(corr_p)] = 1
  corr_p[is.infinite(corr_p)] = 1
  corr_fl = p_filter(mat = corr, mat_p = corr_p, max_p = max_p)
  corr_fl = mat_thresh(mat = corr_fl, th = thresh_hard, soft = FALSE)
  
  # Output
  result = list(cv_error = loss_vec,
                thresh_grid = thresh_grid,
                thresh_opt = thresh_opt,
                mat_cooccur = mat_cooccur,
                corr = corr,
                corr_p = corr_p,
                corr_th = corr_th,
                corr_fl = corr_fl)
  return(result)
}



sparse_dist = function(mat, wins_quant= c(0.05, 0.95), R, thresh_hard = 0, max_p) {
    # Thresholding
    mat_thresh = function(mat, th, soft){
        mat_sign = sign(mat)
        mat_th = mat
        mat_th[abs(mat) <= th] = 0
        if (soft) {
            mat_th[abs(mat) > th] = abs(mat_th[abs(mat) > th]) - th
            mat_th = mat_th * mat_sign
        }
        return(mat_th)
    }

    # Filtering based on p-values
    p_filter = function(mat, mat_p, max_p){
        ind_p = mat_p
        ind_p[mat_p > max_p] = 0
        ind_p[mat_p <= max_p] = 1

        mat_filter = mat * ind_p
        return(mat_filter)
    }

    # Sort taxa
    sort_taxa = sort(colnames(mat))
    mat = mat[, sort_taxa]

    # Winsorization
    mat = apply(mat, 2, function(x)
        DescTools::Winsorize(x, val = quantile(x, probs = wins_quant, na.rm = TRUE)))

    # Co-occurrence matrix
    mat_occur = mat
    mat_occur[mat_occur != 0] = 1
    mat_occur[mat_occur == 0] = 0
    mat_occur[is.na(mat_occur)] = 0

    df_occur = as.data.frame(mat_occur) %>%
        tibble::rownames_to_column("sample_id") %>%
        tidyr::pivot_longer(cols = -.data$sample_id, names_to = "taxon",
                            values_to = "occur") %>%
        dplyr::filter(.data$occur == 1)

    mat_cooccur = crossprod(table(df_occur[, seq_len(2)]))
    diag(mat_cooccur) = colSums(mat_occur)

    if (any(mat_cooccur < 10)) {
        warn_txt = sprintf(paste("There are some pairs of taxa that have insufficient (< 10) overlapping samples",
                                 "Proceed with caution since the point estimates for these pairs are unstable",
                                 "For pairs of taxa with no overlapping samples, the point estimates will be replaced with 0s,",
                                 "and the corresponding p-values will be replaced with 1s",
                                 "Please check `mat_cooccur` for details about the co-occurrence pattern",
                                 sep = "\n"))
        warning(warn_txt)
    }

    # Calculation
    d = dim(mat)[2]
    taxanames = colnames(mat)
    comb = function(...) {
        mapply('rbind', ..., SIMPLIFY = FALSE)
    }

    idx = NULL
    dcorr_list = foreach(idx = seq_len(d - 1), .combine = 'comb', .multicombine = TRUE, .packages = "energy") %dorng% {
        dcorr_idx = rep(NA, d)
        p_val_idx = rep(NA, d)

        mat_x = mat[!is.na(mat[, idx]), ]
        x = mat_x[, idx]

        # Distance correlation
        dcorr_idx[(idx + 1):d] = apply(mat_x[, (idx + 1):d, drop = FALSE], 2,
                                       function(y) {
                                           z = x[!is.na(y)]
                                           y = y[!is.na(y)]
                                           dcor(z, y, index = 1.0)
                                           })

        # P-values
        p_val_idx[(idx + 1):d] = apply(mat_x[, (idx + 1):d, drop = FALSE], 2,
                                       function(y) {
                                           z = x[!is.na(y)]
                                           y = y[!is.na(y)]
                                           dcor.test(z, y, index = 1.0, R = R)$p.value
                                           })

        list(dcorr_idx, p_val_idx)
        }

    dcorr = rbind(dcorr_list[[1]], rep(NA, d))
    dcorr_p = rbind(dcorr_list[[2]], rep(NA, d))
    # Symmetrize the matrix
    dcorr[lower.tri(dcorr)] = t(dcorr)[lower.tri(dcorr)]
    diag(dcorr) = 1
    dcorr[mat_cooccur < 2] = 0
    dcorr[is.infinite(dcorr)] = 0
    dcorr_p[lower.tri(dcorr_p)] = t(dcorr_p)[lower.tri(dcorr_p)]
    diag(dcorr_p) = 0
    dcorr_p[mat_cooccur < 2] = 1
    dcorr_p[is.na(dcorr_p)] = 1
    dcorr_p[is.infinite(dcorr_p)] = 1
    dimnames(dcorr) = list(taxanames, taxanames)
    dimnames(dcorr_p) = list(taxanames, taxanames)

    dcorr_fl = p_filter(mat = dcorr, mat_p = dcorr_p, max_p = max_p)
    dcorr_fl = mat_thresh(mat = dcorr_fl, th = thresh_hard, soft = FALSE)

    # Output
    result = list(mat_cooccur = mat_cooccur,
                  dcorr = dcorr,
                  dcorr_p = dcorr_p,
                  dcorr_fl = dcorr_fl)
    return(result)
}