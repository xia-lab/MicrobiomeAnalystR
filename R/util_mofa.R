create_mofa_from_matrix <- function(data, groups = NULL) {
  
  # Make a dgCMatrix out of dgTMatrix
  if (all(sapply(data, function(x) is(x, "dgTMatrix")))) {
    data <- lapply(data, function(m) as(m, "dgCMatrix"))
  }
  
  # Set groups names
  if (is.null(groups)) {
    # message("No groups provided as argument... we assume that all samples are coming from the same group.\n")
    groups <- rep("group1", ncol(data[[1]]))
  }
  
  # Set views names
  if (is.null(names(data))) {
    default_views <- paste0("view_", seq_len(length(data)))
    message(paste0("View names are not specified in the data, using default: ", paste(default_views, collapse=", "), "\n"))
    names(data) <- default_views
  }
  views_names <- as.character(names(data))
  
  # Initialise MOFA object
  object <- new("MOFA")
  object@status <- "untrained"
  object@data <- .split_data_into_groups(data, groups)
  
  # groups_names <- as.character(unique(groups))
  groups_names <- names(object@data[[1]])
  
  # Set dimensionalities
  object@dimensions[["M"]] <- length(data)
  object@dimensions[["G"]] <- length(groups_names)
  object@dimensions[["D"]] <- sapply(data, nrow)
  object@dimensions[["N"]] <- sapply(groups_names, function(x) sum(groups == x))
  object@dimensions[["K"]] <- 0
  
  # Set features names
  for (m in seq_len(length(data))) {
    if (is.null(rownames(data[[m]]))) {
      warning(sprintf("Feature names are not specified for view %d, using default: feature1_v%d, feature2_v%d...", m, m, m))
      for (g in seq_len(length(object@data[[m]]))) {
        rownames(object@data[[m]][[g]]) <- paste0("feature_", seq_len(nrow(object@data[[m]][[g]])), "_v", m)
      }
    }
  }
  
  # Set samples names
  for (g in seq_len(object@dimensions[["G"]])) {
    if (is.null(colnames(object@data[[1]][[g]]))) {
      warning(sprintf("Sample names for group %d are not specified, using default: sample1_g%d, sample2_g%d,...", g, g, g))
      for (m in seq_len(object@dimensions[["M"]])) {
        colnames(object@data[[m]][[g]]) <- paste0("sample_", seq_len(ncol(object@data[[m]][[g]])), "_g", g)
      }
    }
  }
  
  # Set view names
  views_names(object) <- views_names
  
  # Set samples group names
  groups_names(object) <- groups_names

  # Create sample metadata
  object <- .create_samples_metadata(object)

  # Create features metadata
  object <- .create_features_metadata(object)

  # Rename duplicated features
  object <- .rename_duplicated_features(object)

  # Do quality control
  object <- .quality_control(object)

  return(object)
}


# (Hidden) function to split a list of matrices into a nested list of matrices
.split_data_into_groups <- function(data, groups) {
  group_indices <- split(seq_along(groups), factor(groups, exclude = character(0))) # factor call avoids dropping NA
  lapply(data, function(x) {
    lapply(group_indices, function(idx) {
      x[, idx, drop = FALSE]
    })
  })
}


# (Hidden) function to fill NAs for missing samples
.subset_augment<-function(mat, samp) {
  samp <- unique(samp)
  mat <- t(mat)
  aug_mat<-matrix(NA, ncol=ncol(mat), nrow=length(samp))
  aug_mat<-mat[match(samp,rownames(mat)),,drop=FALSE]
  rownames(aug_mat)<-samp
  colnames(aug_mat)<-colnames(mat)
  return(t(aug_mat))
}

.df_to_matrix <- function(x) {
  m <- as.matrix(x[,-1])
  rownames(m) <- x[[1]]
  if (ncol(m) == 1)
    colnames(m) <- colnames(x)[2:ncol(x)]
  m
}

.create_samples_metadata <- function(object) {
  # TO-DO: CHECK SAMPLE AND GROUP COLUMN IN PROVIDED METADATA
  foo <- lapply(object@data[[1]], colnames)
  tmp <- data.frame(
    sample = unname(unlist(foo)),
    group = unlist(lapply(names(foo), function(x) rep(x, length(foo[[x]])) )),
    stringsAsFactors = FALSE
  )
  if (.hasSlot(object, "samples_metadata") && (length(object@samples_metadata) > 0)) {
    object@samples_metadata <- cbind(tmp, object@samples_metadata[match(tmp$sample, rownames(object@samples_metadata)),, drop = FALSE])
  } else {
    object@samples_metadata <- tmp
  }
  return(object)
}

.create_features_metadata <- function(object) {
  tmp <- data.frame(
    feature = unname(unlist(lapply(object@data, function(x) rownames(x[[1]])))),
    view = unlist(lapply(seq_len(object@dimensions$M), function(x) rep(views_names(object)[[x]], object@dimensions$D[[x]]) )),
    stringsAsFactors = FALSE
  )
  if (.hasSlot(object, "features_metadata") && (length(object@features_metadata) > 0)) {
    object@features_metadata <- cbind(tmp, object@features_metadata[match(tmp$feature, rownames(object@features_metadata)),])
  } else {
    object@features_metadata <- tmp
  }
  return(object)
}

.rename_duplicated_features <- function(object) {
  feature_names <- unname(unlist(lapply(object@data, function(x) rownames(x[[1]]))))
  duplicated_names <- unique(feature_names[duplicated(feature_names)])
  if (length(duplicated_names)>0) 
    warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
  # Rename data matrices
  for (m in names(object@data)) {
    for (g in names(object@data[[m]])) {
      tmp <- which(rownames(object@data[[m]][[g]]) %in% duplicated_names)
      if (length(tmp)>0) {
        rownames(object@data[[m]][[g]])[tmp] <- paste(rownames(object@data[[m]][[g]])[tmp], m, sep="_")
      }
    }
  }
  
  # Rename features metadata
  tmp <- object@features_metadata[["feature"]] %in% duplicated_names
  object@features_metadata[tmp,"feature"] <- paste(object@features_metadata[tmp,"feature"], object@features_metadata[tmp,"view"], sep="_")
  return(object)
}

#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, data_options = data_opts)
get_default_data_options <- function(object) {
  
  # Define default data options
  data_options <- list(
    scale_views = FALSE,     # (logical) Scale views to unit variance?
    scale_groups = FALSE,    # (logical) Scale groups to unit variance?
    center_groups = TRUE,   # (logical) Center groups?
    use_float32 = TRUE       # (logical) Use float32 instead of float64 arrays to increase speed and memory usage
  )
  
  # Activate float32 arrays for large sample sizes  
  if (sum(object@dimensions$N)>1e5) {
    message("A lot of samples detected, using float32 arrays instead of float64 arrays to increase speed and memory usage. 
    You can modify this using the `data_options` argument of the `prepare_mofa` function.")
    data_options$use_float32 <- TRUE
  }
  
  # if data_options already exists, replace the default values but keep the additional ones
  if (length(object@data_options)>0)
    data_options <- modifyList(data_options, object@data_options)
  
  return(data_options)
}

.get_nodes_types <- function() {
  nodes_types <- list(
    multiview_nodes  = c("W", "AlphaW", "ThetaW"),
    multigroup_nodes = c("Z", "AlphaZ", "ThetaZ"),
    twodim_nodes     = c("Y", "Tau"),
    multivariate_singleview_node = "Sigma"
  )
}

get_default_model_options <- function(object) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!.hasSlot(object,"dimensions") | length(object@dimensions) == 0) 
    stop("dimensions of object need to be defined before getting the model options")
  if (.hasSlot(object,"data")) {
    if (length(object@data)==0) stop("data slot is empty")
  } else {
    stop("data slot not found")
  }
  
  # Guess likelihoods from the data
  # likelihoods <- .infer_likelihoods(object)
  likelihoods <- rep(x="gaussian", times=object@dimensions$M)
  names(likelihoods) <- views_names(object)
  
  # Define default model options
  model_options <- list(
    likelihoods = likelihoods,   # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    num_factors = 10,            # (numeric) initial number of latent factors
    spikeslab_factors = FALSE,   # (logical) Spike and Slab sparsity on the factors
    spikeslab_weights = FALSE,    # (logical) Spike and Slab sparsity on the weights
    ard_factors = FALSE,         # (logical) Group-wise ARD sparsity on the factors
    ard_weights = TRUE           # (logical) View-wise ARD sparsity on the weights
  )
  
  # (Heuristic) set the number of factors depending on the sample size
  N <- sum(object@dimensions$N)
  if (N<=25) {
    model_options$num_factors <- 5
  } else if (N>25 & N<=1e3) {
    model_options$num_factors <- 15
  } else if (N>1e3 & N<=1e4) {
    model_options$num_factors <- 20
  } else if (N>1e4) {
    model_options$num_factors <- 25
  }
  
  # Group-wise ARD sparsity on the factors only if there are multiple groups
  if (object@dimensions$G>1)
    model_options$ard_factors <- TRUE
  
  # if model_options already exist, replace the default values but keep the additional ones
  if (length(object@model_options)>0)
    model_options <- modifyList(model_options, object@model_options)
  
  return(model_options)
}

get_default_training_options <- function(object) {
  
  # Get default train options
  training_options <- list(
    maxiter = 1000,                # (numeric) Maximum number of iterations
    convergence_mode = 'fast',     # (string) Convergence mode based on change in the ELBO ("slow","medium","fast")
    drop_factor_threshold = -1,    # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,               # (logical) Verbosity
    startELBO = 1,                 # (numeric) First iteration to compute the ELBO
    freqELBO = 5,                  # (numeric) Frequency of ELBO calculation
    stochastic = FALSE,            # (logical) Do stochastic variational inference?
    gpu_mode = FALSE,              # (logical) Use GPU?
    gpu_device = NULL,             # (integer) Which GPU to use?
    seed = 42,                     # (numeric) random seed
    outfile = NULL,                # (string)  Output file name
    weight_views = FALSE,          # (logical) Weight the ELBO based on the number of features per view?
    save_interrupted = FALSE       # (logical) Save partially trained model when training is interrupted?
  )
  
  # if training_options already exist, replace the default values but keep the additional ones
  if (length(object@training_options)>0)
    training_options <- modifyList(training_options, object@training_options)
  
  return(training_options)
}

#' # Run the MOFA model
#' \dontrun{ MOFAmodel <- run_mofa(MOFAmodel, use_basilisk = TRUE) }
run_mofa <- function(object, outfile = NULL, save_data = TRUE) {
  object <<- object;
  outfile <<- outfile;
  save_data <<- save_data;
  #save.image("mofa.RData");
  if (file.exists(outfile)){
    message(paste0("Warning: Output file ", outfile, " already exists, it will be replaced"))
  }

  message("Connecting to the mofapy2 package using basilisk. Set 'use_basilisk' to FALSE if you prefer to manually set the python binary using 'reticulate'.")
  library('basilisk');
  mofa_env <- BasiliskEnvironment("mofa_env", pkgname="MOFA2", packages=.mofapy2_dependencies, pip = paste0("mofapy2==",.mofapy2_version));
  proc <- basiliskStart(mofa_env)
  on.exit(basiliskStop(proc))
  tmp <- basiliskRun(proc, function(object, outfile, save_data) {
      .run_mofa_reticulate(object, outfile, save_data)
  }, object=object, outfile=outfile, save_data=save_data);
  
  # Load the trained model
  object <- load_model(outfile)
  
  return(object)
}


.run_mofa_reticulate <- function(object, outfile, save_data) {
  
  # sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!requireNamespace("reticulate", quietly = TRUE)) {
    stop("Package \"reticulate\" is required but is not installed.", call. = FALSE)
  }
  library(reticulate);
  
  # Initiate reticulate
  mofa <- reticulate::import("mofapy2")
  
  # Call entry point
  mofa_entrypoint <- mofa$run.entry_point$entry_point()
  
  # Set data options
  mofa_entrypoint$set_data_options(
    scale_views = object@data_options$scale_views,
    scale_groups = object@data_options$scale_groups,
    center_groups = object@data_options$center_groups,
    use_float32 = object@data_options$use_float32
  )

  # Set samples metadata
  if (.hasSlot(object, "samples_metadata")) {
    mofa_entrypoint$data_opts$samples_metadata <- reticulate::r_to_py(lapply(object@data_options$groups,
                                                                 function(g) object@samples_metadata[object@samples_metadata$group == g,]))
  }

  # Set features metadata
  if (.hasSlot(object, "features_metadata")) {
    mofa_entrypoint$data_opts$features_metadata <- reticulate::r_to_py(unname(lapply(object@data_options$views,
                                                                         function(m) object@features_metadata[object@features_metadata$view == m,])))
  }

  # reticulate::r_to_py will convert a list with a single name to a string,
  # hence those are to be wrapped in `list()`
  maybe_list <- function(xs) {
	  if (length(xs) > 1) {
		  xs
	  } else {
		  list(xs)
	  }
  }
  
  # Set the data
  mofa_entrypoint$set_data_matrix(
    data = reticulate::r_to_py( unname(lapply(object@data, function(x) unname( lapply(x, function(y) reticulate::r_to_py(t(y)) ))) ) ),
    likelihoods = unname(object@model_options$likelihoods),
    views_names = reticulate::r_to_py(as.list(object@data_options$views)),
    groups_names = reticulate::r_to_py(as.list(object@data_options$groups)),
    samples_names = reticulate::r_to_py(lapply(unname(lapply(object@data[[1]], colnames)), maybe_list)),
    features_names = reticulate::r_to_py(lapply(unname(lapply(object@data, function(x) rownames(x[[1]]))), maybe_list))
  )

  # Set covariates
  if (.hasSlot(object, "covariates") && !is.null(object@covariates)) {
    sample_cov_to_py <- reticulate::r_to_py(unname(lapply(object@covariates, function(x) unname(reticulate::r_to_py(t(x))))))
    cov_names_2_py <- reticulate::r_to_py(covariates_names(object))
    mofa_entrypoint$set_covariates(sample_cov_to_py, cov_names_2_py)
  }
  
  # Set model options 
  mofa_entrypoint$set_model_options(
    factors     = object@model_options$num_factors,
    spikeslab_factors = object@model_options$spikeslab_factors, 
    spikeslab_weights = object@model_options$spikeslab_weights, 
    ard_factors       = object@model_options$ard_factors,
    ard_weights       = object@model_options$ard_weights 
  )
  
  # Set training options  
  mofa_entrypoint$set_train_options(
    iter             = object@training_options$maxiter,
    convergence_mode = object@training_options$convergence_mode,
    dropR2           = object@training_options$drop_factor_threshold,
    startELBO        = object@training_options$startELBO,
    freqELBO         = object@training_options$freqELBO,
    seed             = object@training_options$seed, 
    gpu_mode         = object@training_options$gpu_mode,
    gpu_device       = object@training_options$gpu_device,
    verbose          = object@training_options$verbose,
    outfile          = object@training_options$outfile,
    weight_views     = object@training_options$weight_views,
    save_interrupted = object@training_options$save_interrupted
  )
  
  
  # Set stochastic options
  if (object@training_options$stochastic) {
    mofa_entrypoint$set_stochastic_options(
      learning_rate    = object@stochastic_options$learning_rate,
      forgetting_rate  = object@stochastic_options$forgetting_rate,
      batch_size       = object@stochastic_options$batch_size,
      start_stochastic = object@stochastic_options$start_stochastic
    )
  }
  
  # Set mefisto options  
  if (.hasSlot(object, "covariates") && !is.null(object@covariates) & length(object@mefisto_options)>1) {
    warping_ref <- which(groups_names(object) == object@mefisto_options$warping_ref)
    mofa_entrypoint$set_smooth_options(
      scale_cov           = object@mefisto_options$scale_cov,
      start_opt           = as.integer(object@mefisto_options$start_opt),
      n_grid              = as.integer(object@mefisto_options$n_grid),
      opt_freq            = as.integer(object@mefisto_options$opt_freq),
      model_groups        = object@mefisto_options$model_groups,
      sparseGP            = object@mefisto_options$sparseGP,
      frac_inducing       = object@mefisto_options$frac_inducing,
      warping             = object@mefisto_options$warping,
      warping_freq        = as.integer(object@mefisto_options$warping_freq),
      warping_ref         = warping_ref-1, # 0-based python indexing
      warping_open_begin  = object@mefisto_options$warping_open_begin,
      warping_open_end    = object@mefisto_options$warping_open_end,
      warping_groups      = reticulate::r_to_py(object@mefisto_options$warping_groups)
    )
  }
  
  # Build the model
  mofa_entrypoint$build();
  
  # Run the model
  mofa_entrypoint$run();

  # Interpolate
  if (.hasSlot(object, "covariates") && !is.null(object@covariates) & length(object@mefisto_options)>1) {
    if(!is.null(object@mefisto_options$new_values)) {
      new_values <- object@mefisto_options$new_values
      if(is.null(dim(new_values))){
        new_values <- matrix(new_values, nrow = 1)
      }
      mofa_entrypoint$predict_factor(new_covariates = reticulate::r_to_py(t(new_values)))
    }
  }
  
  # Save the model output as an hdf5 file
  mofa_entrypoint$save(outfile, save_data = save_data);

}

get_factors <- function(object, groups = "all", factors = "all", scale = FALSE, as.data.frame = FALSE) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get groups
  groups <- .check_and_get_groups(object, groups)
  # Get factors
  factors <- .check_and_get_factors(object, factors)
  
  # Collect factors
  Z <- get_expectations(object, "Z", as.data.frame)
  if (as.data.frame) {
    Z <- Z[Z$factor%in%factors & Z$group%in%groups,]
    if (scale) Z$value <- Z$value/max(abs(Z$value),na.rm=TRUE)
  } else {
    Z <- lapply(Z[groups], function(z) z[,factors, drop=FALSE])
    if (scale) Z <- lapply(Z, function(x) x/max(abs(x)) )
    names(Z) <- groups
  }

  return(Z)
}

get_expectations <- function(object, variable, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(variable %in% names(object@expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  exp <- object@expectations[[variable]]

  # unlist single view nodes - single Sigma node across all groups using time warping
  if(variable == "Sigma")
    exp <- exp[[1]]
  
  # For memory and space efficiency, Y expectations are not saved to the model file when using only gaussian likelihoods.
  if (variable == "Y") {
    if ((length(object@expectations$Y) == 0) && all(object@model_options$likelihood == "gaussian")) {
      # message("Using training data slot as Y expectations since all the likelihoods are gaussian.")
      exp <- object@data
    }
  }
  
  # Convert to long data frame
  if (as.data.frame) {
    
    # Z node
    if (variable=="Z") {
      tmp <- reshape2::melt(exp, na.rm=TRUE)
      colnames(tmp) <- c("sample", "factor", "value", "group")
      tmp$sample <- as.character(tmp$sample)
      factor.cols <- c("sample", "factor", "group")
      factor.cols[factor.cols] <- lapply(factor.cols[factor.cols], factor)
    }
    
    # W node
    else if (variable=="W") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]], na.rm=TRUE)
        colnames(tmp) <- c("feature","factor","value")
        tmp$view <- m
        factor.cols <- c("view", "feature", "factor")
        tmp[factor.cols] <- lapply(tmp[factor.cols], factor)
        return(tmp)
      })
      tmp <- do.call(rbind.data.frame,tmp)
    }
    
    # Y node
    else if (variable=="Y") {
      tmp <- lapply(names(exp), function(m) {
        tmp <- lapply(names(exp[[m]]), function(g) {
          tmp <- reshape2::melt(exp[[m]][[g]], na.rm=TRUE)
          colnames(tmp) <- c("sample", "feature", "value")
          tmp$view <- m
          tmp$group <- g
          factor.cols <- c("view", "group", "feature", "factor")
          tmp[factor.cols] <- lapply(tmp[factor.cols], factor)
          return(tmp)
        })
      })
      tmp <- do.call(rbind, tmp)
    }
    
    exp <- tmp
  }
  return(exp)
}

.check_and_get_factors <- function(object, factors) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(factors)))
  if (is.numeric(factors)) {
    stopifnot(all(factors <= object@dimensions$K))
    factors_names(object)[factors] 
  } else {
    if (paste0(factors, collapse = "") == "all") { 
      factors_names(object)
    } else {
      stopifnot(all(factors %in% factors_names(object)))
      factors
    }
  }
  
}

.check_and_get_groups <- function(object, groups) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(groups)))
  if (is.numeric(groups)) {
    stopifnot(all(groups <= object@dimensions$G))
    groups_names(object)[groups] 
  } else {
    if (paste0(groups, collapse = "") == "all") { 
      groups_names(object)
    } else {
      stopifnot(all(groups %in% groups_names(object)))
      groups
    }
  }
}

get_weights <- function(object, views = "all", factors = "all", abs = FALSE, scale = FALSE, as.data.frame = FALSE) {

  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Get views
  views <- .check_and_get_views(object, views)
  factors <- .check_and_get_factors(object, factors)
  
  # Fetch weights
  weights <- get_expectations(object, "W", as.data.frame)
  
  if (as.data.frame) {
    weights <- weights[weights$view %in% views & weights$factor %in% factors, ]
    if (abs) weights$value <- abs(weights$value)
    if (scale) weights$value <- weights$value/max(abs(weights$value))
  } else {
    weights <- lapply(weights[views], function(x) x[,factors,drop=FALSE])
    if (abs) weights <- lapply(weights, abs)
    if (scale) weights <- lapply(weights, function(x) x/max(abs(x)) )
    names(weights) <- views
  }
  
  return(weights)
}

.check_and_get_views <- function(object, views, non_gaussian=TRUE) {
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(!any(duplicated(views)))
  if (is.numeric(views)) {
    stopifnot(all(views <= object@dimensions$M))
    views <- views_names(object)[views]
  } else {
    if (paste0(views, sep = "", collapse = "") == "all") { 
      views <- views_names(object)
    } else {
      stopifnot(all(views %in% views_names(object)))
    }
  }
  
  # Ignore non-gaussian views  
  if (isFALSE(non_gaussian)) {
    non_gaussian_views <- names(which(object@model_options$likelihoods!="gaussian"))
    views <- views[!views%in%non_gaussian_views]
  }
  
  return(views)
}


#' @importFrom stringi stri_enc_mark
.quality_control <- function(object, verbose = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Check views names
  if (verbose == TRUE) message("Checking views names...")
  stopifnot(!is.null(views_names(object)))
  stopifnot(!duplicated(views_names(object)))
  if (any(grepl("/", views_names(object)))) {
    stop("Some of the views names contain `/` symbol, which is not supported.
  This can be fixed e.g. with:
    views_names(object) <- gsub(\"/\", \"-\", views_names(object))")
  }
  
  # Check groups names
  if (verbose == TRUE) message("Checking groups names...")
  if (any(grepl("/", groups_names(object)))) {
    stop("Some of the groups names contain `/` symbol, which is not supported.
    This can be fixed e.g. with:
    groups_names(object) <- gsub(\"/\", \"-\", groups_names(object))")
  }
  stopifnot(!is.null(groups_names(object)))
  stopifnot(!duplicated(groups_names(object)))
  
  # Check samples names
  if (verbose == TRUE) message("Checking samples names...")
  stopifnot(!is.null(samples_names(object)))
  stopifnot(!duplicated(unlist(samples_names(object))))
  enc <- stringi::stri_enc_mark(unlist(samples_names(object)))
  if (any(enc!="ASCII")) {
    tmp <- unname(unlist(samples_names(object))[enc!="ASCII"])
    stop(sprintf("non-ascii characters detected in the following samples names, please rename them and run again create_mofa():\n- %s ", paste(tmp, collapse="\n- ")))
    print()
  }
  
  # Check features names
  if (verbose == TRUE) message("Checking features names...")
  stopifnot(!is.null(features_names(object)))
  stopifnot(!duplicated(unlist(features_names(object))))
  enc <- stringi::stri_enc_mark(unlist(features_names(object)))
  if (any(enc!="ASCII")) {
    tmp <- unname(unlist(features_names(object))[enc!="ASCII"])
    stop(sprintf("non-ascii characters detected in the following features names, please rename them and run again create_mofa():\n- %s ", paste(tmp, collapse="\n- ")))
    print()
  }
  
  # Check dimensionalities in the input data
  if (verbose == TRUE) message("Checking dimensions...")
  N <- object@dimensions$N
  D <- object@dimensions$D
  for (i in views_names(object)) {
    for (j in groups_names(object)) {
      stopifnot(ncol(object@data[[i]][[j]]) == N[[j]])
      stopifnot(nrow(object@data[[i]][[j]]) == D[[i]])
      stopifnot(length(colnames(object@data[[i]][[j]])) == N[[j]])
      stopifnot(length(rownames(object@data[[i]][[j]])) == D[[i]])
    }
  }
  
  # Check that there are no features with complete missing values (across all groups)
  if (object@status == "untrained" || object@data_options[["loaded"]]) {
      if (verbose == TRUE) message("Checking there are no features with complete missing values...")
      for (i in views_names(object)) {
        if (!(is(object@data[[i]][[1]], "dgCMatrix") || is(object@data[[i]][[1]], "dgTMatrix"))) {
          tmp <- as.data.frame(sapply(object@data[[i]], function(x) rowMeans(is.na(x)), simplify = TRUE))
          if (any(unlist(apply(tmp, 1, function(x) mean(x==1)))==1))
            warning("You have features which do not contain a single observation in any group, consider removing them...")
        }
      }
    }
    
  # check dimensionalities of sample_covariates 
  if (verbose == TRUE) message("Checking sample covariates...")
  if(.hasSlot(object, "covariates") && !is.null(object@covariates)){
    stopifnot(ncol(object@covariates) == sum(object@dimensions$N))
    stopifnot(nrow(object@covariates) == object@dimensions$C)
    stopifnot(all(unlist(samples_names(object)) == colnames(object@covariates)))
  }
  
  # Sanity checks that are exclusive for an untrained model  
  if (object@status == "untrained") {
    
    # Check features names
    if (verbose == TRUE) message("Checking features names...")
    tmp <- lapply(object@data, function(x) unique(lapply(x,rownames)))
    for (x in tmp) stopifnot(length(x)==1)
    for (x in tmp) if (any(duplicated(x[[1]]))) stop("There are duplicated features names within the same view. Please rename")
    all_names <- unname(unlist(tmp))
    duplicated_names <- unique(all_names[duplicated(all_names)])
    if (length(duplicated_names)>0) 
      warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
    for (i in names(object@data)) {
      for (j in names(object@data[[i]])) {
        tmp <- which(rownames(object@data[[i]][[j]]) %in% duplicated_names)
        if (length(tmp)>0) {
          rownames(object@data[[i]][[j]])[tmp] <- paste(rownames(object@data[[i]][[j]])[tmp], i, sep="_")
        }
      }
    }
    
  # Sanity checks that are exclusive for a trained model  
  } else if (object@status == "trained") {
    # Check expectations
    if (verbose == TRUE) message("Checking expectations...")
    stopifnot(all(c("W", "Z") %in% names(object@expectations)))
    # if(.hasSlot(object, "covariates") && !is.null(object@covariates)) stopifnot("Sigma" %in% names(object@expectations))
    stopifnot(all(sapply(object@expectations$W, is.matrix)))
    stopifnot(all(sapply(object@expectations$Z, is.matrix)))
    
    # Check for intercept factors
    if (object@data_options[["loaded"]]) { 
      if (verbose == TRUE) message("Checking for intercept factors...")
      if (!is.null(object@data)) {
        factors <- do.call("rbind",get_factors(object))
        r <- suppressWarnings( t(do.call('rbind', lapply(object@data, function(x) 
          abs(cor(colMeans(do.call("cbind",x),na.rm=TRUE),factors, use="pairwise.complete.obs"))
        ))) )
        intercept_factors <- which(rowSums(r>0.75)>0)
        if (length(intercept_factors)) {
            warning(sprintf("Factor(s) %s are strongly correlated with the total number of expressed features for at least one of your omics. Such factors appear when there are differences in the total 'levels' between your samples, *sometimes* because of poor normalisation in the preprocessing steps.\n",paste(intercept_factors,collapse=", ")))
        }
      }
    }
  
    # Check for correlated factors
    if (verbose == TRUE) message("Checking for highly correlated factors...")
    Z <- do.call("rbind",get_factors(object))
    op <- options(warn=-1) # suppress warnings
    
    noise <- matrix(rnorm(n=length(Z), mean=0, sd=1e-10), nrow(Z), ncol(Z))
    tmp <- cor(Z+noise); diag(tmp) <- NA
    options(op) # activate warnings again
    if (max(tmp,na.rm=TRUE)>0.5) {
      warning("The model contains highly correlated factors (see `plot_factor_cor(MOFAobject)`). \nWe recommend that you train the model with less factors and that you let it train for a longer time.\n")
    }
  
  }
  
  return(object)  
}


prepare_mofa <- function(object, data_options = NULL, model_options = NULL, 
                         training_options = NULL, stochastic_options = NULL, mefisto_options = NULL) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (any(object@dimensions$N<10) & !(.hasSlot(object, "covariates") && length(object@covariates)>=1)) warning("Some group(s) have less than 10 samples, MOFA will have little power to learn meaningful factors for these group(s)...")
  if (any(object@dimensions$D<15)) warning("Some view(s) have less than 15 features, MOFA will have little power to to learn meaningful factors for these view(s)....")
  if (any(object@dimensions$D>1e4)) warning("Some view(s) have a lot of features, it is recommended to perform a more stringent feature selection before creating the MOFA object....")
  if (length(object@samples_metadata)>0) { 
    stopifnot(c("sample","group") %in% colnames(object@samples_metadata))
  } else {
    stop("object@samples_metadata not found") 
  }
  if (length(object@features_metadata)>0) { 
    stopifnot(c("feature","view") %in% colnames(object@features_metadata))
  } else {
    stop("object@features_metadata not found") 
  }
  if (object@dimensions$G>1) {
    message("\n# Multi-group mode requested.")
    message("\nThis is an advanced option, if this is the first time that you are running MOFA, we suggest that you try do some exploration first without specifying groups. Two important remarks:")
    message("\n - The aim of the multi-group framework is to identify the sources of variability *within* the groups. If your aim is to find a factor that 'separates' the groups, you DO NOT want to use the multi-group framework. Please see the FAQ on the MOFA2 webpage.") 
    message("\n - It is important to account for the group effect before selecting highly variable features (HVFs). We suggest that either you calculate HVFs per group and then take the union, or regress out the group effect before HVF selection")
  }

  # Get data options
  message("Checking data options...")
  if (is.null(data_options)) {
    message("No data options specified, using default...")
    object@data_options <- get_default_data_options(object)
  } else {
    if (!is(data_options,"list") || !setequal(names(data_options), names(get_default_data_options(object)) ))
      stop("data_options are incorrectly specified, please read the documentation in get_default_data_options")
    object@data_options <- data_options
  }
  # if (any(nchar(unlist(samples_names(object)))>50))
  #   warning("Due to string size limitations in the HDF5 format, sample names will be trimmed to less than 50 characters")
  
  # Get training options
  if (is.null(training_options)) {
    message("No training options specified, using default...")
    object@training_options <- get_default_training_options(object)
  } else {
    message("Checking training options...")
    if (!is(training_options,"list") || !setequal(names(training_options), names(get_default_training_options(object)) ))
      stop("training_options are incorrectly specified, please read the documentation in get_default_training_options")
    object@training_options <- training_options
    
    if (object@training_options$maxiter<=100)
      warning("Maximum number of iterations is very small\n")
    if (object@training_options$startELBO<1) object@training_options$startELBO <- 1
    if (object@training_options$freqELBO<1) object@training_options$freqELBO <- 1
    if (!object@training_options$convergence_mode %in% c("fast","medium","slow")) 
      stop("Convergence mode has to be either 'fast', 'medium', or 'slow'")
  }
  
  # Get stochastic options
  if (is.null(stochastic_options)) {
    object@stochastic_options <- list()
  } else {
    if (isFALSE(object@training_options[["stochastic"]]))
      stop("stochastic_options have been provided but training_opts$stochastic is FALSE. If you want to use stochastic inference you have to set training_opts$stochastic = TRUE")
    # object@training_options$stochastic <- TRUE
  }
  
  if (object@training_options$stochastic) {
    message("Stochastic inference activated. Note that this is only recommended if you have a very large sample size (>1e4) and access to a GPU")
    
    if (is.null(stochastic_options)) {
      message("No stochastic options specified, using default...")
      object@stochastic_options <- get_default_stochastic_options(object)
    } else {
      object@training_options$stochastic <- TRUE
      message("Checking stochastic inference options...")
      if (!is(stochastic_options,"list") || !setequal(names(stochastic_options), names(get_default_stochastic_options(object)) ))
        stop("stochastic_options are incorrectly specified, please read the documentation in get_default_stochastic_options")
      
      if (!stochastic_options$batch_size %in% c(0.05,0.10,0.15,0.20,0.25,0.50))
        stop("Batch size has to be one of the following numeric values: 0.05, 0.10, 0.15, 0.20, 0.25, 0.50")
      if (stochastic_options$batch_size==1)
        warning("A batch size equal to 1 is equivalent to non-stochastic inference.")
      if (stochastic_options$learning_rate<=0 || stochastic_options$learning_rate>1)
        stop("The learning rate has to be a value between 0 and 1")
      if (stochastic_options$forgetting_rate<=0 || stochastic_options$forgetting_rate>1)
        stop("The forgetting rate has to be a value between 0 and 1")
      
      if (sum(object@dimensions$N)<1e4) warning("Stochastic inference is only recommended when you have a lot of samples (at least N>10,000))\n")
      
      object@stochastic_options <- stochastic_options
    }
  }
  
  # Get model options
  if (is.null(model_options)) {
    message("No model options specified, using default...")
    object@model_options <- get_default_model_options(object)
  } else {
    message("Checking model options...")
    if (!is(model_options,"list") || !setequal(names(model_options), names(get_default_model_options(object)) ))
      stop("model_options are incorrectly specified, please read the documentation in get_default_model_options")
    object@model_options <- model_options
  }
  if (object@model_options$num_factors > 50) warning("The number of factors is very large, training will be slow...")
  # if (!object@model_options$ard_weights) warning("model_options$ard_weights should always be set to TRUE")
  if (sum(object@dimensions$N) < 4 * object@model_options$num_factors) {
    warning(sprintf("The total number of samples is very small for learning %s factors.  
    Try to reduce the number of factors to obtain meaningful results. It should not exceed ~%s.",
                    object@model_options$num_factors, floor(min(object@dimensions$N/4))))
  }
  
  # Get mefisto covariates options
  if (.hasSlot(object, "covariates") && length(object@covariates)>=1) {
    if (is.null(mefisto_options)) {
        message("Covariates provided but no mefisto options specified, using default...")
        object@mefisto_options <- get_default_mefisto_options(object)
    } else {
      message("Checking inference options for mefisto covariates...")
      # message("mefisto covariates have been provided as prior information.")
      if (!is(mefisto_options,"list") || !setequal(names(mefisto_options), names(get_default_mefisto_options(object)) ))
        stop("mefisto_options are incorrectly specified, please read the documentation in get_default_mefisto_options")
      
      if (isTRUE(mefisto_options$sparseGP)) {
        if (object@dimensions[["N"]] < 1000) warning("Warning: sparseGPs should only be used when having a large sample size (>1e3)")
        if (isTRUE(mefisto_options$warping)) stop("Warping is not implemented in conjunction with sparseGPs")
      }
      
      # Check warping options
      if (isTRUE(mefisto_options$warping)) {
        stopifnot(object@dimensions[['G']] > 1) # check that multi-group is TRUE
        
        if (!is.null(mefisto_options$warping_ref)) {
          stopifnot(length(mefisto_options$warping_ref)==1)
          stopifnot(is.character(mefisto_options$warping_ref))
          stopifnot(mefisto_options$warping_ref %in% groups_names(object))
        }
        
        if (!is.null(mefisto_options$warping_groups)) {
          # check that warping groups are a partition of groups
          groups_ok <- sapply(unique(object@samples_metadata$group), function(g) {
            length(unique(mefisto_options$warping_groups[object@samples_metadata$group == g])) == 1
          })
          if (!all(groups_ok)) stop("Warping group assignment needs to be unique within each indiviudal group.")
        }
      }
      
      # Disable spike-slab on the factors
      if(isTRUE(model_options$spikeslab_factors)) {
        print("Spike-and-Slab sparsity prior on the factors is not available when using MEFISTO, setting to False")
        model_options$spikeslab_factors <- FALSE
      }
      
      # Disable stochastic inference
      if (isTRUE(model_options$stochastic)) {
        print("Stochastic inference is not available when using MEFISTO, setting to False")
        model_options$stochastic <- FALSE
        object@stochastic_options <- list()
      }
      
      # TO-DO: CHECKS ON MODEL_GROUPS
      
      object@mefisto_options <- mefisto_options
    }
    
  } else {
    object@mefisto_options <- list()
  }
  
  # Center the data
  # message("Centering the features (per group, this is a mandatory requirement)...")
  # for (m in views_names(object)) {
  #   if (model_options$likelihoods[[m]] == "gaussian") {
  #     for (g in groups_names(object)) {
  #       object@data[[m]][[g]] <- scale(object@data[[m]][[g]], center=T, scale=F)
  #     }
  #   }
  # }
  
  # Transform sparse matrices into dense ones
  # See https://github.com/rstudio/reticulate/issues/72
  for (m in views_names(object)) {
    for (g in groups_names(object)) {
      if (is(object@data[[m]][[g]], "dgCMatrix") || is(object@data[[m]][[g]], "dgTMatrix"))
        object@data[[m]][[g]] <- as(object@data[[m]][[g]], "matrix")
    }
  }
  
  return(object)
}



#' @title Get default training options
#' @name get_default_training_options
#' @description Function to obtain the default training options.
#' @param object an untrained \code{\link{MOFA}}
#' @details This function provides a default set of training options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step (see example), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' The training options are the following: \cr
#' \itemize{
#'  \item{\strong{maxiter}:}{ numeric value indicating the maximum number of iterations. 
#'  Default is 1000. Convergence is assessed using the ELBO statistic.}
#'  \item{\strong{drop_factor_threshold}:}{ numeric indicating the threshold on fraction of variance explained to consider a factor inactive and drop it from the model.
#'  For example, a value of 0.01 implies that factors explaining less than 1\% of variance (in each view) will be dropped. Default is -1 (no dropping of factors)}
#'  \item{\strong{convergence_mode}:}{ character indicating the convergence criteria, either "fast", "medium" or "slow", corresponding to 0.0005\%, 0.00005\% or 0.000005\% deltaELBO change. }
#'  \item{\strong{verbose}:}{ logical indicating whether to generate a verbose output.}
#'  \item{\strong{startELBO}:}{ integer indicating the first iteration to compute the ELBO (default is 1). }
#'  \item{\strong{freqELBO}:}{ integer indicating the first iteration to compute the ELBO (default is 1). }
#'  \item{\strong{stochastic}:}{ logical indicating whether to use stochastic variational inference (only required for very big data sets, default is \code{FALSE}).}
#'  \item{\strong{gpu_mode}:}{ logical indicating whether to use GPUs (see details).}
#'  \item{\strong{gpu_device}:}{ integer indicating which GPU to use.}
#'  \item{\strong{seed}:}{ numeric indicating the seed for reproducibility (default is 42).}
#' }
#' @return Returns a list with default training options
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Load default training options
#' train_opts <- get_default_training_options(MOFAmodel)
#' 
#' # Edit some of the training options
#' train_opts$convergence_mode <- "medium"
#' train_opts$startELBO <- 100
#' train_opts$seed <- 42
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, training_options = train_opts)
get_default_training_options <- function(object) {
  
  # Get default train options
  training_options <- list(
    maxiter = 1000,                # (numeric) Maximum number of iterations
    convergence_mode = 'fast',     # (string) Convergence mode based on change in the ELBO ("slow","medium","fast")
    drop_factor_threshold = -1,    # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,               # (logical) Verbosity
    startELBO = 1,                 # (numeric) First iteration to compute the ELBO
    freqELBO = 5,                  # (numeric) Frequency of ELBO calculation
    stochastic = FALSE,            # (logical) Do stochastic variational inference?
    gpu_mode = FALSE,              # (logical) Use GPU?
    gpu_device = NULL,             # (integer) Which GPU to use?
    seed = 42,                     # (numeric) random seed
    outfile = NULL,                # (string)  Output file name
    weight_views = FALSE,          # (logical) Weight the ELBO based on the number of features per view?
    save_interrupted = FALSE       # (logical) Save partially trained model when training is interrupted?
  )
  
  # if training_options already exist, replace the default values but keep the additional ones
  if (length(object@training_options)>0)
    training_options <- modifyList(training_options, object@training_options)
  
  return(training_options)
}


#' @title Get default data options
#' @name get_default_data_options
#' @description Function to obtain the default data options.
#' @param object an untrained \code{\link{MOFA}} object
#' @details This function provides a default set of data options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step (see example), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' The data options are the following: \cr
#' \itemize{
#'  \item{\strong{scale_views}:}{ logical indicating whether to scale views to have the same unit variance. 
#'  As long as the scale differences between the views is not too high, this is not required. Default is FALSE.}
#'  \item{\strong{scale_groups}:}{ logical indicating whether to scale groups to have the same unit variance. 
#'  As long as the scale differences between the groups is not too high, this is not required. Default is FALSE.}
#'  \item{\strong{use_float32}:}{ logical indicating whether use float32 instead of float64 arrays to increase speed and memory usage. Default is FALSE.}
#'  }
#' @return Returns a list with the default data options.
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Load default data options
#' data_opts <- get_default_data_options(MOFAmodel)
#' 
#' # Edit some of the data options
#' data_opts$scale_views <- TRUE
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, data_options = data_opts)
get_default_data_options <- function(object) {
  
  # Define default data options
  data_options <- list(
    scale_views = FALSE,     # (logical) Scale views to unit variance?
    scale_groups = FALSE,    # (logical) Scale groups to unit variance?
    center_groups = TRUE,   # (logical) Center groups?
    use_float32 = TRUE       # (logical) Use float32 instead of float64 arrays to increase speed and memory usage
  )
  
  # Activate float32 arrays for large sample sizes  
  if (sum(object@dimensions$N)>1e5) {
    message("A lot of samples detected, using float32 arrays instead of float64 arrays to increase speed and memory usage. 
    You can modify this using the `data_options` argument of the `prepare_mofa` function.")
    data_options$use_float32 <- TRUE
  }
  
  # if data_options already exists, replace the default values but keep the additional ones
  if (length(object@data_options)>0)
    data_options <- modifyList(data_options, object@data_options)
  
  return(data_options)
}

#' @title Get default model options
#' @name get_default_model_options
#' @description Function to obtain the default model options.
#' @param object an untrained \code{\link{MOFA}} object
#' @details This function provides a default set of model options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step (see example), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' The model options are the following: \cr
#' \itemize{
#'  \item{\strong{likelihoods}:}{ character vector with data likelihoods per view: 
#'  'gaussian' for continuous data (Default for all views), 'bernoulli' for binary data and 'poisson' for count data.}
#'  \item{\strong{num_factors}:}{ numeric value indicating the (initial) number of factors. Default is 15.}
#'  \item{\strong{spikeslab_factors}:}{ logical indicating whether to use spike and slab sparsity on the factors (Default is FALSE)}
#'  \item{\strong{spikeslab_weights}:}{ logical indicating whether to use spike and slab sparsity on the weights (Default is TRUE)}
#'  \item{\strong{ard_factors}:}{ logical indicating whether to use ARD sparsity on the factors (Default is TRUE only if using multiple groups)}
#'  \item{\strong{ard_weights}:}{ logical indicating whether to use ARD sparsity on the weights (Default is TRUE)}
#'  }
#' @return Returns a list with the default model options.
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # Load default model options
#' model_opts <- get_default_model_options(MOFAmodel)
#' 
#' # Edit some of the model options
#' model_opts$num_factors <- 10
#' model_opts$spikeslab_weights <- FALSE
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, model_options = model_opts)
get_default_model_options <- function(object) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  if (!.hasSlot(object,"dimensions") | length(object@dimensions) == 0) 
    stop("dimensions of object need to be defined before getting the model options")
  if (.hasSlot(object,"data")) {
    if (length(object@data)==0) stop("data slot is empty")
  } else {
    stop("data slot not found")
  }
  
  # Guess likelihoods from the data
  # likelihoods <- .infer_likelihoods(object)
  likelihoods <- rep(x="gaussian", times=object@dimensions$M)
  names(likelihoods) <- views_names(object)
  
  # Define default model options
  model_options <- list(
    likelihoods = likelihoods,   # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    num_factors = 10,            # (numeric) initial number of latent factors
    spikeslab_factors = FALSE,   # (logical) Spike and Slab sparsity on the factors
    spikeslab_weights = FALSE,    # (logical) Spike and Slab sparsity on the weights
    ard_factors = FALSE,         # (logical) Group-wise ARD sparsity on the factors
    ard_weights = TRUE           # (logical) View-wise ARD sparsity on the weights
  )
  
  # (Heuristic) set the number of factors depending on the sample size
  N <- sum(object@dimensions$N)
  if (N<=25) {
    model_options$num_factors <- 5
  } else if (N>25 & N<=1e3) {
    model_options$num_factors <- 15
  } else if (N>1e3 & N<=1e4) {
    model_options$num_factors <- 20
  } else if (N>1e4) {
    model_options$num_factors <- 25
  }
  
  # Group-wise ARD sparsity on the factors only if there are multiple groups
  if (object@dimensions$G>1)
    model_options$ard_factors <- TRUE
  
  # if model_options already exist, replace the default values but keep the additional ones
  if (length(object@model_options)>0)
    model_options <- modifyList(model_options, object@model_options)
  
  return(model_options)
}


#' @title Get default stochastic options
#' @name get_default_stochastic_options
#' @description Function to obtain the default options for stochastic variational inference.
#' @param object an untrained \code{\link{MOFA}}
#' @details This function provides a default set of stochastic inference options that can be modified and passed to the \code{\link{MOFA}} object
#' in the \code{\link{prepare_mofa}} step), i.e. after creating a \code{\link{MOFA}} object
#'  (using \code{\link{create_mofa}}) and before starting the training (using \code{\link{run_mofa}})
#' These options are only relevant when activating stochastic inference in training_options (see example).
#' The stochastic inference options are the following: \cr
#' \itemize{
#'  \item{\strong{batch_size}:}{ numeric value indicating the batch size (as a fraction)}. 
#'  Default is 0.5 (half of the data set).
#'  \item{\strong{learning_rate}:}{ numeric value indicating the learning rate. }
#'  Default is 1.0
#'  \item{\strong{forgetting_rate}:}{ numeric indicating the forgetting rate.}
#'  Default is 0.5
#'  \item{\strong{start_stochastic}:}{ integer indicating the first iteration to start stochastic inference}
#'  Default is 1
#'  }
#' @return Returns a list with default options
#' @importFrom utils modifyList
#' @export
#' @examples
#' # Using an existing simulated data with two groups and two views
#' file <- system.file("extdata", "test_data.RData", package = "MOFA2")
#' 
#' # Load data dt (in data.frame format)
#' load(file) 
#' 
#' # Create the MOFA object
#' MOFAmodel <- create_mofa(dt)
#' 
#' # activate stochastic inference in training options
#' train_opts <- get_default_training_options(MOFAmodel)
#' train_opts$stochastic <- TRUE
#' 
#' # Load default stochastic options
#' stochastic_opts <- get_default_stochastic_options(MOFAmodel)
#' 
#' # Edit some of the stochastic options
#' stochastic_opts$learning_rate <- 0.75
#' stochastic_opts$batch_size <- 0.25
#' 
#' # Prepare the MOFA object
#' MOFAmodel <- prepare_mofa(MOFAmodel, 
#'   training_options = train_opts,
#'   stochastic_options = stochastic_opts
#' )
#' 
get_default_stochastic_options <- function(object) {
  
  # Get default stochastic options
  stochastic_options <- list(
    batch_size = 0.5,        # Batch size (as a fraction)
    learning_rate = 1.0,     # Starting learning rate
    forgetting_rate = 0.5,   # Forgetting rate
    start_stochastic = 1     # First iteration to start stochastic inference
  )
  
  # if stochastic_options already exist, replace the default values but keep the additional ones
  if (length(object@stochastic_options)>0)
    stochastic_options <- modifyList(stochastic_options, object@stochastic_options)
  
  return(stochastic_options)
}

############################################
## Functions to load a trained MOFA model ##
############################################

#' @title Load a trained MOFA
#' @name load_model
#' @description Method to load a trained MOFA \cr
#' The training of mofa is done using a Python framework, and the model output is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the mofa Python framework
#' @param sort_factors logical indicating whether factors should be sorted by variance explained (default is TRUE)
#' @param on_disk logical indicating whether to work from memory (FALSE) or disk (TRUE). \cr
#' This should be set to TRUE when the training data is so big that cannot fit into memory. \cr
#' On-disk operations are performed using the \code{\link{HDF5Array}} and \code{\link{DelayedArray}} framework.
#' @param load_data logical indicating whether to load the training data (default is TRUE, it can be memory expensive)
#' @param remove_outliers logical indicating whether to mask outlier values.
#' @param remove_inactive_factors logical indicating whether to remove inactive factors from the model.
# #' @param remove_intercept_factors logical indicating whether to remove intercept factors for non-Gaussian views.
#' @param verbose logical indicating whether to print verbose output (default is FALSE)
#' @param load_interpol_Z (MEFISTO) logical indicating whether to load predictions for factor values based on latent processed (only
#'  relevant for models trained with covariates and Gaussian processes, where prediction was enabled)
#' @return a \code{\link{MOFA}} model
#' @importFrom rhdf5 h5read h5ls
#' @importFrom HDF5Array HDF5ArraySeed
#' @importFrom DelayedArray DelayedArray
#' @importFrom dplyr bind_rows
#' @export
#' @examples
#' #' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)

load_model <- function(file, sort_factors = TRUE, on_disk = FALSE, load_data = TRUE,
                       remove_outliers = FALSE, remove_inactive_factors = TRUE, verbose = FALSE,
                       load_interpol_Z = FALSE) {

  # Create new MOFAodel object
  object <- new("MOFA")
  object@status <- "trained"
  
  # Set on_disk option
  if (on_disk) { 
    object@on_disk <- TRUE 
  } else { 
      object@on_disk <- FALSE 
  }
  library('rhdf5');
  # Get groups and data set names from the hdf5 file object
  h5ls.out <- h5ls(file, datasetinfo = FALSE)
  
  ########################
  ## Load training data ##
  ########################

  # Load names
  if ("views" %in% h5ls.out$name) {
    view_names <- as.character( h5read(file, "views")[[1]] )
    group_names <- as.character( h5read(file, "groups")[[1]] )
    feature_names <- h5read(file, "features")[view_names]
    sample_names  <- h5read(file, "samples")[group_names] 
  } else {  # for old models
    feature_names <- h5read(file, "features")
    sample_names  <- h5read(file, "samples")
    view_names <- names(feature_names)
    group_names <- names(sample_names)
    h5ls.out <- h5ls.out[grep("variance_explained", h5ls.out$name, invert = TRUE),]
  }
  if("covariates" %in%  h5ls.out$name){
    covariate_names <- as.character( h5read(file, "covariates")[[1]])
  } else {
    covariate_names <- NULL
  }

  # Load training data (as nested list of matrices)
  data <- list(); intercepts <- list()
  if (load_data && "data"%in%h5ls.out$name) {
    
    object@data_options[["loaded"]] <- TRUE
    if (verbose) message("Loading data...")
    
    for (m in view_names) {
      data[[m]] <- list()
      intercepts[[m]] <- list()
      for (g in group_names) {
        if (on_disk) {
          # as DelayedArrays
          data[[m]][[g]] <- DelayedArray::DelayedArray( HDF5ArraySeed(file, name = sprintf("data/%s/%s", m, g) ) )
        } else {
          # as matrices
          data[[m]][[g]] <- h5read(file, sprintf("data/%s/%s", m, g) )
          tryCatch(intercepts[[m]][[g]] <- as.numeric( h5read(file, sprintf("intercepts/%s/%s", m, g) ) ), error = function(e) { NULL })
        }
        # Replace NaN by NA
        data[[m]][[g]][is.nan(data[[m]][[g]])] <- NA # this realised into memory, TO FIX
      }
    }
    
  # Create empty training data (as nested list of empty matrices, with the correct dimensions)
  } else {
    
    object@data_options[["loaded"]] <- FALSE
    
    for (m in view_names) {
      data[[m]] <- list()
      for (g in group_names) {
        data[[m]][[g]] <- .create_matrix_placeholder(rownames = feature_names[[m]], colnames = sample_names[[g]])
      }
    }
  }

  object@data <- data
  object@intercepts <- intercepts


  # Load metadata if any
  if ("samples_metadata" %in% h5ls.out$name) {
    object@samples_metadata <- bind_rows(lapply(group_names, function(g) as.data.frame(h5read(file, sprintf("samples_metadata/%s", g)))))
  }
  if ("features_metadata" %in% h5ls.out$name) {
    object@features_metadata <- bind_rows(lapply(view_names, function(m) as.data.frame(h5read(file, sprintf("features_metadata/%s", m)))))
  }
  
  ############################
  ## Load sample covariates ##
  ############################
  
  if (any(grepl("cov_samples", h5ls.out$group))){
    covariates <- list()
    for (g in group_names) {
      if (on_disk) {
        # as DelayedArrays
        covariates[[g]] <- DelayedArray::DelayedArray( HDF5ArraySeed(file, name = sprintf("cov_samples/%s", g) ) )
      } else {
        # as matrices
        covariates[[g]] <- h5read(file, sprintf("cov_samples/%s", g) )
      }    
    }
  } else covariates <- NULL
  object@covariates <- covariates

  if (any(grepl("cov_samples_transformed", h5ls.out$group))){
    covariates_warped <- list()
    for (g in group_names) {
      if (on_disk) {
        # as DelayedArrays
        covariates_warped[[g]] <- DelayedArray::DelayedArray( HDF5ArraySeed(file, name = sprintf("cov_samples_transformed/%s", g) ) )
      } else {
        # as matrices
        covariates_warped[[g]] <- h5read(file, sprintf("cov_samples_transformed/%s", g) )
      }    
    }
  } else covariates_warped <- NULL
  object@covariates_warped <- covariates_warped
  
  #######################
  ## Load interpolated factor values ##
  #######################
  
  interpolated_Z <- list()
  if (isTRUE(load_interpol_Z)) {
    
    if (isTRUE(verbose)) message("Loading interpolated factor values...")
    
    for (g in group_names) {
      interpolated_Z[[g]] <- list()
      if (on_disk) {
        # as DelayedArrays
        # interpolated_Z[[g]] <- DelayedArray::DelayedArray( HDF5ArraySeed(file, name = sprintf("Z_predictions/%s", g) ) )
      } else {
        # as matrices
        tryCatch( {
          interpolated_Z[[g]][["mean"]] <- h5read(file, sprintf("Z_predictions/%s/mean", g) )
        }, error = function(x) { print("Predicitions of Z not found, not loading it...") })
        tryCatch( {
          interpolated_Z[[g]][["variance"]] <- h5read(file, sprintf("Z_predictions/%s/variance", g) )
        }, error = function(x) { print("Variance of predictions of Z not found, not loading it...") })
        tryCatch( {
          interpolated_Z[[g]][["new_values"]] <- h5read(file, "Z_predictions/new_values")
        }, error = function(x) { print("New values of Z not found, not loading it...") })
      }
    }
  }
  object@interpolated_Z <- interpolated_Z
  
  #######################
  ## Load expectations ##
  #######################

  expectations <- list()
  node_names <- h5ls.out[h5ls.out$group=="/expectations","name"]

  if (verbose) message(paste0("Loading expectations for ", length(node_names), " nodes..."))

  if ("AlphaW" %in% node_names)
    expectations[["AlphaW"]] <- h5read(file, "expectations/AlphaW")[view_names]
  if ("AlphaZ" %in% node_names)
    expectations[["AlphaZ"]] <- h5read(file, "expectations/AlphaZ")[group_names]
  if ("Sigma" %in% node_names)
    expectations[["Sigma"]] <- h5read(file, "expectations/Sigma")
  if ("Z" %in% node_names)
    expectations[["Z"]] <- h5read(file, "expectations/Z")[group_names]
  if ("W" %in% node_names)
    expectations[["W"]] <- h5read(file, "expectations/W")[view_names]
  if ("ThetaW" %in% node_names)
    expectations[["ThetaW"]] <- h5read(file, "expectations/ThetaW")[view_names]
  if ("ThetaZ" %in% node_names)
    expectations[["ThetaZ"]] <- h5read(file, "expectations/ThetaZ")[group_names]
  # if ("Tau" %in% node_names)
  #   expectations[["Tau"]] <- h5read(file, "expectations/Tau")
  
  object@expectations <- expectations

  
  ########################
  ## Load model options ##
  ########################

  if (verbose) message("Loading model options...")

  tryCatch( {
    object@model_options <- as.list(h5read(file, 'model_options', read.attributes = TRUE))
  }, error = function(x) { print("Model options not found, not loading it...") })

  # Convert True/False strings to logical values
  for (i in names(object@model_options)) {
    if (object@model_options[i] == "False" || object@model_options[i] == "True") {
      object@model_options[i] <- as.logical(object@model_options[i])
    } else {
      object@model_options[i] <- object@model_options[i]
    }
  }

  ##########################################
  ## Load training options and statistics ##
  ##########################################

  if (verbose) message("Loading training options and statistics...")

  # Load training options
  if (length(object@training_options) == 0) {
    tryCatch( {
      object@training_options <- as.list(h5read(file, 'training_opts', read.attributes = TRUE))
    }, error = function(x) { print("Training opts not found, not loading it...") })
  }

  # Load training statistics
  tryCatch( {
    object@training_stats <- h5read(file, 'training_stats', read.attributes = TRUE)
    object@training_stats <- h5read(file, 'training_stats', read.attributes = TRUE)
  }, error = function(x) { print("Training stats not found, not loading it...") })

  #############################
  ## Load covariates options ##
  #############################
  
  if (any(grepl("cov_samples", h5ls.out$group))) { 
    if (isTRUE(verbose)) message("Loading covariates options...")
    tryCatch( {
      object@mefisto_options <- as.list(h5read(file, 'smooth_opts', read.attributes = TRUE))
    }, error = function(x) { print("Covariates options not found, not loading it...") })
    
    # Convert True/False strings to logical values
    for (i in names(object@mefisto_options)) {
      if (object@mefisto_options[i] == "False" | object@mefisto_options[i] == "True") {
        object@mefisto_options[i] <- as.logical(object@mefisto_options[i])
      } else {
        object@mefisto_options[i] <- object@mefisto_options[i]
      }
    }
    
  }
  
  
    
  #######################################
  ## Load variance explained estimates ##
  #######################################
  
  if ("variance_explained" %in% h5ls.out$name) {
    r2_list <- list(
      r2_total = h5read(file, "variance_explained/r2_total")[group_names],
      r2_per_factor = h5read(file, "variance_explained/r2_per_factor")[group_names]
    )
    object@cache[["variance_explained"]] <- r2_list
  }
  
  # Hack to fix the problems where variance explained values range from 0 to 1 (%)
  if (max(sapply(object@cache$variance_explained$r2_total,max,na.rm=TRUE),na.rm=TRUE)<1) {
    for (m in 1:length(view_names)) {
      for (g in 1:length(group_names)) {
        object@cache$variance_explained$r2_total[[g]][[m]] <- 100 * object@cache$variance_explained$r2_total[[g]][[m]]
        object@cache$variance_explained$r2_per_factor[[g]][,m] <- 100 * object@cache$variance_explained$r2_per_factor[[g]][,m]
      }
    }
  }
  
  ##############################
  ## Specify dimensionalities ##
  ##############################
  
  # Specify dimensionality of the data
  object@dimensions[["M"]] <- length(data)                            # number of views
  object@dimensions[["G"]] <- length(data[[1]])                       # number of groups
  object@dimensions[["N"]] <- sapply(data[[1]], ncol)                 # number of samples (per group)
  object@dimensions[["D"]] <- sapply(data, function(e) nrow(e[[1]]))  # number of features (per view)
  object@dimensions[["C"]] <- nrow(covariates[[1]])                        # number of covariates
  object@dimensions[["K"]] <- ncol(object@expectations$Z[[1]])        # number of factors
  
  # Assign sample and feature names (slow for large matrices)
  if (verbose) message("Assigning names to the different dimensions...")

  # Create default features names if they are null
  if (is.null(feature_names)) {
    print("Features names not found, generating default: feature1_view1, ..., featureD_viewM")
    feature_names <- lapply(seq_len(object@dimensions[["M"]]),
                            function(m) sprintf("feature%d_view_&d", as.character(seq_len(object@dimensions[["D"]][m])), m))
  } else {
    # Check duplicated features names
    all_names <- unname(unlist(feature_names))
    duplicated_names <- unique(all_names[duplicated(all_names)])
    if (length(duplicated_names)>0) 
      warning("There are duplicated features names across different views. We will add the suffix *_view* only for those features 
            Example: if you have both TP53 in mRNA and mutation data it will be renamed to TP53_mRNA, TP53_mutation")
    for (m in names(feature_names)) {
      tmp <- which(feature_names[[m]] %in% duplicated_names)
      if (length(tmp)>0) feature_names[[m]][tmp] <- paste(feature_names[[m]][tmp], m, sep="_")
    }
  }
  features_names(object) <- feature_names
  
  # Create default samples names if they are null
  if (is.null(sample_names)) {
    print("Samples names not found, generating default: sample1, ..., sampleN")
    sample_names <- lapply(object@dimensions[["N"]], function(n) paste0("sample", as.character(seq_len(n))))
  }
  samples_names(object) <- sample_names

  # Add covariates names
  if(!is.null(object@covariates)){
    # Create default covariates names if they are null
    if (is.null(covariate_names)) {
      print("Covariate names not found, generating default: covariate1, ..., covariateC")
      covariate_names <- paste0("sample", as.character(seq_len(object@dimensions[["C"]])))
    }
    covariates_names(object) <- covariate_names
  }
  
  # Set views names
  if (is.null(names(object@data))) {
    print("Views names not found, generating default: view1, ..., viewM")
    view_names <- paste0("view", as.character(seq_len(object@dimensions[["M"]])))
  }
  views_names(object) <- view_names
  
  # Set groups names
  if (is.null(names(object@data[[1]]))) {
    print("Groups names not found, generating default: group1, ..., groupG")
    group_names <- paste0("group", as.character(seq_len(object@dimensions[["G"]])))
  }
  groups_names(object) <- group_names
  
  # Set factors names
  factors_names(object)  <- paste0("Factor", as.character(seq_len(object@dimensions[["K"]])))
  
  ###################
  ## Parse factors ##
  ###################
  
  # Calculate variance explained estimates per factor
  if (is.null(object@cache[["variance_explained"]])) {
    object@cache[["variance_explained"]] <- calculate_variance_explained(object)
  } 
  
  # Remove inactive factors
  if (remove_inactive_factors) {
    r2 <- rowSums(do.call('cbind', lapply(object@cache[["variance_explained"]]$r2_per_factor, rowSums, na.rm=TRUE)))
    var.threshold <- 0.0001
    if (all(r2 < var.threshold)) {
      warning(sprintf("All %s factors were found to explain little or no variance so remove_inactive_factors option has been disabled.", length(r2)))
    } else if (any(r2 < var.threshold)) {
      object <- subset_factors(object, which(r2>=var.threshold), recalculate_variance_explained=FALSE)
      message(sprintf("%s factors were found to explain no variance and they were removed for downstream analysis. You can disable this option by setting load_model(..., remove_inactive_factors = FALSE)", sum(r2 < var.threshold)))
    }
  }
  
  # [Done in mofapy2] Sort factors by total variance explained
  if (sort_factors && object@dimensions$K>1) {

    # Sanity checks
    if (verbose) message("Re-ordering factors by their variance explained...")

    # Calculate variance explained per factor across all views
    r2 <- rowSums(sapply(object@cache[["variance_explained"]]$r2_per_factor, function(e) rowSums(e, na.rm = TRUE)))
    order_factors <- c(names(r2)[order(r2, decreasing = TRUE)])

    # re-order factors
    object <- subset_factors(object, order_factors)
  }

  # Mask outliers
  if (remove_outliers) {
    if (verbose) message("Removing outliers...")
    object <- .detect_outliers(object)
  }
  
  # Mask intercepts for non-Gaussian data
  if (any(object@model_options$likelihoods!="gaussian")) {
    for (m in names(which(object@model_options$likelihoods!="gaussian"))) {
      for (g in names(object@intercepts[[m]])) {
        object@intercepts[[m]][[g]] <- NA
      }
    }
  }

  ######################
  ## Quality controls ##
  ######################

  if (verbose) message("Doing quality control...")
  object <- .quality_control(object, verbose = verbose)
  
  return(object)
}


################################
## Functions to do subsetting ##
################################

#' @title Subset groups
#' @name subset_groups
#' @description Method to subset (or sort) groups
#' @param object a \code{\link{MOFA}} object.
#' @param groups character vector with the groups names, numeric vector with the groups indices
#' or logical vector with the groups to be kept as TRUE.
#' @return A \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Subset the first group
#' model <- subset_groups(model, groups = 1)
subset_groups <- function(object, groups) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(groups) <= object@dimensions[["G"]])
  
  # Define groups    
  groups <- .check_and_get_groups(object, groups)
  
  # Subset expectations
  if (length(object@expectations)>0) {
    if ("Z" %in% names(object@expectations) & length(object@expectations$Z)>0)
      object@expectations$Z <- object@expectations$Z[groups]
    if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0)
      object@expectations$Y <- sapply(object@expectations$Y, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }
  
  # Subset data
  if (length(object@data)>0) {
    object@data <- sapply(object@data, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }

  # Subset imputed data
  if (length(object@imputed_data)>0) {
    object@imputed_data <- sapply(object@imputed_data, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }
  
  # Subset intercepts
  if (length(object@intercepts[[1]])>0) {
    object@intercepts <- sapply(object@intercepts, function(x) x[groups], simplify = FALSE, USE.NAMES = TRUE) 
  }
  
  # Update dimensionality
  object@dimensions[["G"]] <- length(groups)
  object@dimensions[["N"]] <- object@dimensions[["N"]][groups]
  
  # Subset sample metadata
  stopifnot(groups%in%unique(object@samples_metadata$group))
  object@samples_metadata <- object@samples_metadata[object@samples_metadata$group %in% groups,]
  object@samples_metadata$group <- factor(object@samples_metadata$group, levels=groups)
  
  # Re-order samples
  samples <- unname(unlist(lapply(object@data[[1]],colnames)))
  object@samples_metadata <- object@samples_metadata[match(samples, object@samples_metadata$sample),]
  
  # Sanity checks
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@data[[1]],colnames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Z,rownames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Y,colnames)))
  
  # Update groups names
  # groups_names(object) <- groups # don't need to run this
  object@data_options$groups <- groups
  
  # Subset variance explained
  object@cache[["variance_explained"]]$r2_per_factor <- object@cache[["variance_explained"]]$r2_per_factor[groups]
  object@cache[["variance_explained"]]$r2_total <- object@cache[["variance_explained"]]$r2_total[groups]
  
  return(object)
}


#' @title Subset views
#' @name subset_views
#' @description Method to subset (or sort) views
#' @param object a \code{\link{MOFA}} object.
#' @param views character vector with the views names, numeric vector with the views indices,
#' or logical vector with the views to be kept as TRUE.
#' @return A \code{\link{MOFA}} object
#' @export
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Subset the first view
#' model <- subset_views(model, views = 1)
subset_views <- function(object, views) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(views) <= object@dimensions[["M"]])
  # warning("Removing views a posteriori is fine for an exploratory analysis, but you should removing them before training!")
  
  # Define views
  views <- .check_and_get_views(object, views)
  
  # Subset relevant slots
  if (length(object@expectations)>0) {
    object@expectations$W <- object@expectations$W[views]
    object@expectations$Y <- object@expectations$Y[views]
  }

  # Subset data
  if (length(object@data)>0) {
    object@data <- object@data[views]
  }
  
  # Subset imputed data
  if (length(object@imputed_data)>0) {
    object@imputed_data <- object@imputed_data[views]
  }
  
  # Subset intercepts
  if (length(object@intercepts[[1]])>0) {
    object@intercepts <- object@intercepts[views]
  }
  
  # Subset feature metadata
  if (length(object@features_metadata)>0) {
    object@features_metadata <- object@features_metadata[object@features_metadata$view %in% views,]
  }
  
  # Subset likelihoods
  object@model_options$likelihoods <- object@model_options$likelihoods[views]
  # Update dimensionality
  object@dimensions[["M"]] <- length(views)
  object@dimensions[["D"]] <- object@dimensions[["D"]][views]
  
  # Update view names
  object@data_options$views <- views
  
  # Subset variance explained
  if ((methods::.hasSlot(object, "cache")) && ("variance_explained" %in% names(object@cache))) {
    object@cache[["variance_explained"]]$r2_per_factor <- lapply(object@cache[["variance_explained"]]$r2_per_factor, function(x) x[,views,drop=FALSE])
    object@cache[["variance_explained"]]$r2_total <- lapply(object@cache[["variance_explained"]]$r2_total, function(x) x[views])
  }
  
  return(object)
}


#' @title Subset factors
#' @name subset_factors
#' @description Method to subset (or sort) factors
#' @param object a \code{\link{MOFA}} object.
#' @param factors character vector with the factor names, or numeric vector with the index of the factors.
#' @param recalculate_variance_explained logical indicating whether to recalculate variance explained values. Default is \code{TRUE}.
#' @export
#' @return A \code{\link{MOFA}} object
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # Subset factors 1 to 3
#' model <- subset_factors(model, factors = 1:3)
subset_factors <- function(object, factors, recalculate_variance_explained = TRUE) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(factors) <= object@dimensions[["K"]])
  
  # Define factors
  factors <- .check_and_get_factors(object, factors)
  
  # Subset expectations
  nodes_with_factors <- list(nodes = c("Z", "W", "Sigma", "AlphaZ", "AlphaW", "ThetaZ", "ThetaW"), axes = c(2, 2, 1, 0, 0, 0, 0))
  stopifnot(all(nodes_with_factors$axes %in% c(0, 1, 2)))

  if (length(object@expectations)>0) {
    for (i in seq_len(length(nodes_with_factors$nodes))) {
      node <- nodes_with_factors$nodes[i]
      axis <- nodes_with_factors$axes[i]
      if (node %in% names(object@expectations)) {
        if(node != "Sigma") {
        if (axis == 1) {
          object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[factors,,drop=FALSE], simplify = FALSE, USE.NAMES = TRUE)
        } else if (axis == 2) {
          object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[,factors,drop=FALSE], simplify = FALSE, USE.NAMES = TRUE)
        } else {
          object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[factors], simplify = FALSE, USE.NAMES = TRUE)
        }
        } else {
          if (axis == 1) {
            object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[factors,,,drop=FALSE], simplify = FALSE, USE.NAMES = TRUE)
          } else if (axis == 2) {
            object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[,factors,,drop=FALSE], simplify = FALSE, USE.NAMES = TRUE)
          } else if (axis == 3) {
            object@expectations[[node]] <- sapply(object@expectations[[node]], function(x) x[,,factors,drop=FALSE], simplify = FALSE, USE.NAMES = TRUE)
          } 
      }
    }
    }
  }
  
  # Subset interpolations
  if(length(object@interpolated_Z) > 0) {
    object@interpolated_Z <- lapply(object@interpolated_Z, function(g) { 
      if(!is.null(g$mean)) {
        m <- g$mean[factors, , drop = FALSE]
      }
      if(!is.null(g$variance)) {
        v <- g$variance[factors, , drop = FALSE]
      }
      list(mean = m, variance = v, new_values = g$new_values)
    })
  }
  
  # Remove total variance explained estimates  
  if (length(factors)<object@dimensions[["K"]]) {
    object@cache[["variance_explained"]]$r2_total <- lapply(object@cache[["variance_explained"]]$r2_per_factor, colSums)
    if (recalculate_variance_explained) {
      message("Recalculating total variance explained values (r2_total)...")
      object@cache[["variance_explained"]]$r2_total <- calculate_variance_explained(object)[["r2_total"]]
    }
  }
  
  # # Relalculate total variance explained estimates (not valid for non-orthogonal factors)
  # if (length(factors) < object@dimensions[["K"]]) {
  #   object@cache[["variance_explained"]]$r2_total <- lapply(object@cache[["variance_explained"]]$r2_per_factor, colSums)
  
  # Subset per-factor variance explained estimates
  if ((methods::.hasSlot(object, "cache")) && ("variance_explained" %in% names(object@cache))) {
    object@cache[["variance_explained"]]$r2_per_factor <- lapply(object@cache[["variance_explained"]]$r2_per_factor, function(x) x[factors,,drop=FALSE])
  }
  
  # Subset lengthscales per factor
  if (!is.null(object@training_stats$structural_sig)) {
    object@training_stats$structural_sig <- object@training_stats$structural_sig[factors,drop=FALSE]
  }
  if (!is.null(object@training_stats$length_scales)) {
    object@training_stats$length_scales <- object@training_stats$length_scales[factors,drop=FALSE]
  }
  if (!is.null(object@training_stats$scales)) {
    object@training_stats$scales <- object@training_stats$scales[factors,drop=FALSE]
  }
  
  # Update dimensionality
  object@dimensions[["K"]] <- length(factors)
  
  # Update factor names
  factors_names(object) <- paste0("Factor", as.character(seq_len(object@dimensions[["K"]])))
  
  
  return(object)
}



#' @title Subset samples
#' @name subset_samples
#' @description Method to subset (or sort) samples
#' @param object a \code{\link{MOFA}} object.
#' @param samples character vector with the sample names or numeric vector with the sample indices.
#' @export
#' @return A \code{\link{MOFA}} object
#' @examples
#' # Using an existing trained model on simulated data
#' file <- system.file("extdata", "model.hdf5", package = "MOFA2")
#' model <- load_model(file)
#' 
#' # (TO-DO) Remove a specific sample from the model (an outlier)
subset_samples <- function(object, samples) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  
  # Define samples
  samples <- .check_and_get_samples(object, samples)
  
  # Check if an entire group needs to be removed
  groups <- as.character(unique(object@samples_metadata[match(samples, object@samples_metadata$sample),]$group))
  if (length(groups)<length(groups_names(object))) object <- subset_groups(object, groups)
  
  # Subset data and expectations
  groups <- groups_names(object)
  tmp <- lapply(groups, function(g) samples_names(object)[[g]][samples_names(object)[[g]] %in% samples])
  names(tmp) <- groups
  
  for (g in groups) {
    samples_g <- tmp[[g]]
    
    # Subset expectations
    if (length(object@expectations)>0) {
      if ("Z" %in% names(object@expectations) & length(object@expectations$Z)>0) {
        object@expectations$Z[[g]] <- object@expectations$Z[[g]][samples_g,, drop=FALSE]
      }
      
      if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0) {
        for (m in views_names(object)) {
          object@expectations$Y[[m]][[g]] <- object@expectations$Y[[m]][[g]][,samples_g,drop=FALSE]
        }  
      }
      if ("Tau" %in% names(object@expectations) & length(object@expectations$Tau)>0) {
        for (m in views_names(object)) {
          object@expectations$Tau[[m]][[g]] <- object@expectations$Tau[[m]][[g]][samples_g, , drop=FALSE]
        }  
      }
    if(g == groups[1]) {# only one Sigma node
      if ("Sigma" %in% names(object@expectations) & length(object@expectations$Sigma)>0) {
        samples <- unique(unlist(tmp)) # TODO - make Sigma live on covariate level or expand to group-level
        object@expectations$Sigma[[1]] <- object@expectations$Sigma[[1]][,samples,samples,drop=FALSE]
      }
    }
    }

    # Subset data
    if (length(object@data)>0) { 
      for (m in views_names(object)) {
        object@data[[m]][[g]] <- object@data[[m]][[g]][,samples_g,drop=FALSE]
      }
    }
    
    # Subset imputed data
    if (length(object@imputed_data)>0) { 
      for (m in views_names(object)) {
        object@imputed_data[[m]][[g]] <- object@imputed_data[[m]][[g]][,samples_g,drop=FALSE]
      }
    }
    
  }

  # Subset sample metadata
  object@samples_metadata <- object@samples_metadata[object@samples_metadata$sample %in% samples,]
  
  # Update dimensionality
  object@dimensions[["N"]] <- sapply(tmp, length)
  
  # Update sample names
  samples_names(object) <- tmp
  
  # Sanity checks
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@data[[1]],colnames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Z,rownames)))
  stopifnot(object@samples_metadata$sample == unlist(lapply(object@expectations$Y,colnames)))
  
  # Remove variance explained estimates  
  warning("After subsetting the samples the variance explained estimates are not valid anymore, removing them...")
  object@cache[["variance_explained"]] <- NULL
  
  return(object)
}


#' @title Subset features
#' @name subset_features
#' @description Method to subset (or sort) features
#' @param object a \code{\link{MOFA}} object.
#' @param view character vector with the view name or integer with the view index
#' @param features character vector with the sample names, numeric vector with the feature indices 
#' or logical vector with the samples to be kept as TRUE.
#' @return A \code{\link{MOFA}} object
#' @export

subset_features <- function(object, view, features) {
  
  # Sanity checks
  if (!is(object, "MOFA")) stop("'object' has to be an instance of MOFA")
  stopifnot(length(features) <= sapply(object@dimensions[["D"]], sum))
  warning("Removing features a posteriori is fine for an exploratory analysis, but we recommend removing them before training!")

  if (is.numeric(view)) view <- views_names(object)[view]
  stopifnot(all(view %in% views_names(object)))  

  # Define features
  if (is.character(features)) {
    stopifnot(all(features %in% features_names(object)[[view]]))
  } else {
    features <- features_names(object)[[view]][features]
  }
  
  # Subset relevant slots
  if (length(object@expectations)>0) {
    if ("W" %in% names(object@expectations) & length(object@expectations$W)>0)
      object@expectations$W <- lapply(object@expectations$W, function(x) x[features,, drop=FALSE])
    if ("Y" %in% names(object@expectations) & length(object@expectations$Y)>0)
      object@expectations$Y[[view]] <- lapply(object@expectations$Y[[view]], function(x) x[features,])

  if (length(object@data)>0)
    object@data <- lapply(object@data, function(x) sapply(x, function(y) y[features,], simplify = FALSE, USE.NAMES = TRUE))

  if (length(object@expectations)>0)
  object@intercepts <- lapply(object@intercepts, function(x) sapply(x, function(y) y[features], simplify = FALSE, USE.NAMES = TRUE))

  if (length(object@imputed_data) != 0) {
    stop()
    # object@imputed_data <- lapply(object@imputed_data, function(x) sapply(x, function(y) y[,samples], simplify = FALSE, USE.NAMES = TRUE)) 
    }
  }
  # Update dimensionality
  object@dimensions[["D"]][[view]] <- length(features)
  
  # Update features names
  features_names(object)[[view]] <- features
  
  # Remove variance explained estimates  
  warning("After subsetting the features the variance explained estimates are not valid anymore, removing them...")
  object@cache[["variance_explained"]] <- NULL
  
  return(object)
}
