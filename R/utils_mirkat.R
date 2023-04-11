####Functions adapted from the MiRKAT package https://doi.org/10.1093/bioinformatics/btaa951


MiRKAT <- function(y, X = NULL, Ks, out_type = "C", 
                  method = "davies", omnibus = "permutation", 
                  nperm = 999, returnKRV = FALSE, returnR2 = FALSE){
  
  method <- match.arg(method, choices = c("davies", "permutation", "moment"))
  omnibus <- match.arg(tolower(omnibus), choices = c("permutation", "cauchy"))
  
  n = length(y)
  
  if (any(is.na(y))){
    ids = which(is.na(y))
    stop(paste("missing response for subject(s)", ids, ", please remove before proceeding \n")) 
  }
  
  if(is.null(X)==FALSE){
    if(NROW(X)!= length(y)) stop("Dimensions of X and y don't match.")
  }
  
  if (is.matrix(Ks)) {
    Ks = list(Ks)
  }
  
  if (is.list(Ks)) {  
    if((any(lapply(Ks, "nrow")!= n))|(any(lapply(Ks,  "ncol")!= n))){
      stop("Distance matrix need to be n x n, where n is the sample size \n ")    
    } 
    if (!is.list(Ks)) {
      stop("Distance needs to be a list of n x n matrices or a single n x n matrix \n")  
    }
  }
  
  
  if (!is.null(X)){
    if (any(is.na(X))){
      stop("NAs in  covariates X, please impute or remove subjects with missing covariate values") 
    }  
  }
  
  if (method == "moment" & n < 100 & out_type == "C"){
    
    warning("Continuous outcome: sample size < 100, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
  }
  if (method == "moment" & n < 200 & out_type == "D"){
    
    warning("Continuous outcome: sample size < 200, p-value using moment matching can be inaccurate at tails, davies or permutation is recommended")
  }
  
  if (!(out_type %in%  c("C", "D"))){
    stop("Only continuous and binary outcomes are supported by this function. Please choose out_type = \"C\" or \"D\", or use an alternative function for other outcome types.")
  }

  library("CompQuadForm");
  if(out_type  == "C"){
    re = MiRKAT_continuous(y, X = X, Ks = Ks, method = method, omnibus = omnibus, nperm = nperm, returnKRV = returnKRV, returnR2 = returnR2)  
  }
  
  if(out_type  == "D"){
    re = MiRKAT_binary(y, X = X, Ks = Ks, method = method, omnibus = omnibus, nperm = nperm, returnKRV = returnKRV, returnR2 = returnR2)  
  }
  
  return(re)
}

MiRKAT_continuous <- function(y, X = NULL, Ks, method, omnibus, 
                             nperm = 999, returnKRV = FALSE, returnR2 = FALSE){
  
  om <- substring(tolower(omnibus), 1, 1)
  n = length(y) 
  
  if (is.null(X)) {
    X1 <-  matrix(rep(1, length(y)), ncol=1)
  } else {
    X1 <- model.matrix(~. , as.data.frame(X))
  }
  
  KRVs = R2 = NULL 
  if (returnKRV | returnR2) {
    reskrv = scale(resid(lm(y ~ 1)))
    L = reskrv %*% t(reskrv)
    if (returnKRV) { KRVs <- unlist(lapply(Ks, FUN = function(k) calcKRVstat(k, L))) }
    if (returnR2) { R2 <- unlist(lapply(Ks, FUN = function(k) calcRsquared(k, L))) }
  }
  
  ## Prep RHS; take care of aliased variables and pivoting
  qX1 <- qr(X1)
  X1 <- X1[, qX1$pivot, drop=FALSE]
  X1 <- X1[, 1:qX1$rank, drop=FALSE]
  
  mod <- lm(y ~ X1-1)
  s2  = summary(mod)$s**2
  D0  = diag(n) 
  res = resid(mod)
  P0  = D0 - X1%*%solve(t(X1)%*%X1)%*%t(X1)
  px  = ncol(X1)
  
  ## Individual test statistics 
  Qs <- c() 
  for (i in 1:length(Ks)) {
    Qs <- c(Qs, getQ(K = Ks[[i]], res, s2))
  }
  
  ## Permuted test stats 
  if (method == "permutation" | (length(Ks) > 1 & om == "p")) {
    q_sim =  sapply(1:nperm, function(i) {
      ind <- sample(n)
      sapply(1:length(Ks),function(j){
        Q1 = as.numeric(res %*% Ks[[j]][ind, ind] %*% res/s2) 
      })
    })  
    q_sim = t(q_sim)
  }
  
  ## Individual p-values 
  if (method == "davies") {    
    lambda0 = lapply(Ks, getLambda_davies, P0)
    ps = rep(NA, length(Ks))
    for (i in 1:length(Ks)){
      ps[i] = getindivP_davies(Qs[[i]], lambda0[[i]], n, px)    
    }
    out_pvs = ps 
  } else if (method == "moment"){ 
    parm = sapply(Ks, getParamSat, P0, px)  
    p_sat = rep(NA, length(Ks))
    for (i in 1:length(Ks)){
      p_sat[i] = getSat(Qs[[i]], keppa_tlt = parm[,i]$keppa_tlt, niu_tlt = parm[,i]$niu_tlt)    
    }
    out_pvs = p_sat 
  } else if (method == substr("permutation", 1, nchar(method))) {
    # P = (# of statistics permutation >= observed test statistics + 1) / (#permutation + 1), 
    # thus the smallest p value will be  1/ (# permutation + 1).
    if (length(Ks) > 1) {
      Q_all = rbind(unlist(Qs), q_sim)
      p_all = 1 - (apply(Q_all, 2, rank) - 1)/(nperm + 1) 
      p_perm = p_all[1,]
      out_pvs = p_perm 
    } else {
      Q_all <- c(unlist(Qs), q_sim)
      p_all <- 1 - (rank(Q_all) - 1)/(nperm + 1)
      p_perm = p_all[1] 
      out_pvs <- p_perm
    }
    
  }
  
  ## name indiv p-values 
  names(out_pvs) = names(Ks) 
  
  if (length(out_pvs) == 1) {
    if (!returnKRV & !returnR2) {
      return(list(p_values = out_pvs))
    } else if (!returnKRV & returnR2) {
      return(list(p_values = out_pvs, R2 = R2))
    } else if (returnKRV & !returnR2) {
      return(list(p_values = out_pvs, KRV = KRVs))
    } else {
      return(list(p_values = out_pvs, KRV = KRVs, R2 = R2))    
    }
  }
  
  ## calculate omnibus p-value 
  if (om == "p") {
    Q_all = rbind(unlist(Qs), q_sim)
    p_all = 1 - (apply(Q_all, 2, rank)-1)/(nperm + 1) 
    p_perm = p_all[1,]
    minP_all= apply(p_all,1, min)
    p_final = rank(minP_all)[1]/(nperm  + 1)  
  } else if (om == "c") {
    cauchy.t <- sum(tan((0.5 - out_pvs)*pi))/length(out_pvs)
    p_final <- 1 - pcauchy(cauchy.t)
  } else {
    stop("I don't know that omnibus option. Please choose 'permutation' or 'Cauchy'.")
  }
  
  
  ## return all 
  if (is.null(KRVs) & is.null(R2)) {
    return(list(p_values = out_pvs , omnibus_p = p_final))
  } else if (is.null(KRVs) & !is.null(R2)) {
    return(list(p_values = out_pvs , omnibus_p = p_final, R2 = R2))
  } else if (!is.null(KRVs) & is.null(R2)) {
    return(list(p_values = out_pvs , omnibus_p = p_final, KRV = KRVs))
  } else {
    return(list(p_values = out_pvs , omnibus_p = p_final, KRV = KRVs, R2 = R2))    
  }
}

# KRV test statistic 
calcKRVstat <- function(K, L) {
  n = nrow(K) 
  I.n=diag(1,n)
  I.1=rep(1,n)
  H=I.n-I.1%*%t(I.1)/n
  K=H%*%K%*%H
  L=H%*%L%*%H
  A=K/tr(K%*%K)  ## standard-version of K
  W=L/tr(L%*%L)
  Fstar=tr(A%*%W)
  return(Fstar)
}

# trace of a matrix 
tr<- function(x){return(sum(diag(x))) }
calcRsquared <- function(K, L) {
  r1 <- cor(as.numeric(K), as.numeric(L)) 
  return(r1^2)
}

getQ <-  function(K, res, s2){    
  Q = c(1 / s2 * res %*% K %*% res)
}

getLambda_davies <-  function(K, P0){
  PKP = P0 %*% K %*% P0
  ee = eigen(PKP, symmetric = T)         
  lambda0 = ee$values[ee$values > 1e-10]
  return(lambda0)    
}

# Continuous outcome individual Davies p-values 
getindivP_davies <-  function(Q, lambda0, n, px){
  if (length(lambda0) >= n-px){ 
    # In fact, it can not go bigger than n-p because P0 has rank n-p
    lambda = c(lambda0 - c(Q)/(n-px))
    k = length(lambda)
    h = rep(1,k)
  }else{
    lambda = c(lambda0 - c(Q)/(n-px), c(-Q)/(n-px))
    k = length(lambda0)
    h = c(rep(1, k), n-px - k)
  }
  p_davies = davies(q = 0, lambda = lambda, h = h, acc = 1e-6)$Qq  
  p_davies = ifelse(p_davies < 0, 0, p_davies) 
  
  return(p_davies)  
}


D2K <- function(D){
    n <- nrow(D)
    centerM <- diag(n) - 1/n
    K <- -0.5*centerM %*% (D*D) %*% centerM
    eK <- eigen(K, symmetric=TRUE)
    K <- eK$vector %*% diag(pmax(0,eK$values)) %*% t(eK$vector)
    return(K)
  }