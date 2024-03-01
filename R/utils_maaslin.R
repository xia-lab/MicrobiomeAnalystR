
.perform.my.maaslin <- function(
    mbSetObj,
    analysis.var,
    is.norm = "false",
    comp = NULL,
    ref = NULL,
    block = "NA",
    taxrank = "NA",
    model = "LM",
    imgNm = "NA",
    thresh = 0.05){ 
  require(dplyr);
  require(R.utils); 

  mbSetObj <- .get.mbSetObj(mbSetObj);
  
  if (!exists('adj.vec')) {
    adj.bool <- F;
  } else {
    if (length(adj.vec) > 0) {
      adj.bool <- T;
    } else {
      adj.bool <- F;
    }
  }

  if(length(adj.vec) == 1){
    if(adj.vec == "") {adj.bool <- F}
  }

  thresh <- as.numeric(thresh);
  adj.vars <- adj.vec;
  
  if(!adj.bool){
    adj.vec <- "NA"; #for report purposes
  }

  if (mbSetObj$module.type=="sdp"){
    taxrank<-"OTU";
    if (is.norm == "false"){
      phyloseq_objs <- qs::qread("phyloseq_prenorm_objs.qs")
      input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
    } else {
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
      input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
    }
  } else {
    if (is.norm == "false"){
      phyloseq_objs <- qs::qread("phyloseq_prenorm_objs.qs")
      
      if (taxrank=="OTU"){
        input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
      } else {
        taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
        input.data <- phyloseq_objs$count_tables[[taxrank.inx]] %>% as.data.frame()
      } 
    } else {
      phyloseq_objs <- qs::qread("phyloseq_objs.qs")
      if (taxrank=="OTU") {
        input.data <- phyloseq_objs$count_tables$OTU %>% as.data.frame()
      } else {
        taxrank.inx <- which(names(phyloseq_objs$count_tables) %in% taxrank)
        input.data <- phyloseq_objs$count_tables[[taxrank.inx]] %>% as.data.frame()
      } 
    }
  }
  
  meta.nms <- colnames(mbSetObj$dataSet$sample_data)
  input.meta <- mbSetObj$dataSet$sample_data@.Data %>% as.data.frame()
  colnames(input.meta) <- meta.nms
  rownames(input.meta) <- input.meta$sample_id
  if(adj.bool){
    fixed.effects <- c(analysis.var, adj.vars)
    fixed.types <- mbSetObj$dataSet$meta.types[names(mbSetObj$dataSet$meta.types) %in% fixed.effects]
    fixed.types <- fixed.types[match(fixed.effects, names(fixed.types))]
  } else { # to do still
    fixed.effects <- analysis.var
    fixed.types <- mbSetObj$dataSet$meta.types[names(mbSetObj$dataSet$meta.types) == analysis.var]
  }
  
  analysis.type <- fixed.types[fixed.effects == analysis.var]
  disc.effects <- fixed.effects[fixed.types == "disc"]
  mbSetObj$analSet$adj.bool = adj.bool
  mbSetObj$analSet$block = block
  mbSetObj$analSet$analysis.type = analysis.type
  mbSetObj$analSet$disc.effects = disc.effects
  # build refs vector (may need to add for blocking too)
  if(length(disc.effects) > 0){
    if(analysis.type == "disc"){
      refs <- paste0(analysis.var, ",", ref)
      if(length(disc.effects) > 1){
        for(i in c(2:length(disc.effects))){
          ref.temp <- paste0(disc.effects[i], ",", levels(unlist(c(input.meta[,disc.effects[i]])))[1])
          refs <- c(refs, ref.temp)
        }
      }
    } else {
      refs <- c()
      if(length(disc.effects) > 1){
        for(i in c(1:length(disc.effects))){
          ref.temp <- paste0(disc.effects[i], ",", levels(unlist(c(input.meta[,disc.effects[i]])))[1])
          refs <- c(refs, ref.temp)
        }
      }
    }
  }

  # MaAslin does not require samples or orders to exactly match - it takes care of this
  # set normalized/transformation parameters
  if(is.norm == "false"){
    norm.method = "TSS"
    trans.method = "LOG"
  } else {
    norm.method = "NONE"
    trans.method = "NONE"
  }
if(model =="CPLM"){
  trans.method = "NONE"
}else  if(model=="NEGBIN"|model=="ZINB"){
    norm.method = "NONE"
    trans.method = "NONE"
}

  #record params for report
  mbSetObj$paramSet$cov <- list(
    primaryMeta = analysis.var,
    covariates = adj.vec,
    block = block,
    taxrank = taxrank,
    model = model,
    comparison = comp,
    reference = ref,
    p.lvl = thresh
  )
 
###set adjust paprameter
   if((!adj.bool) & (block == "NA")){
    maaslin.para.adj<<-0
  } else {
    refs2 <- refs[grep(paste0(analysis.var, ","), refs)];
    
   maaslin.para.noadj<<- list(
       input_data = input.data, 
       input_metadata = input.meta, 
      fixed_effects = c(analysis.var),
      reference = c(refs2),
       max_significance = 0.05,
      min_abundance = 0.0,
     min_prevalence = 0.0,
      min_variance = 0.0,
      normalization = norm.method,
     transform = trans.method);
  }

.set.mbSetObj(mbSetObj)
  if(block == "NA"){
    if(length(disc.effects) > 0){ # case: discrete variables, no blocking factor
     maaslin.para<<- list(input_data = input.data, 
        input_metadata = input.meta, 
        fixed_effects = c(fixed.effects),
        reference = c(refs),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method)
      return(1)
    } else { # case: no discrete variables, no blocking factor
    maaslin.para<<- list(input_data = input.data, 
        fixed_effects = c(fixed.effects),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method)
      return(1)
    }
  } else { # case: discrete variables, blocking factor (blocking factor must be discrete)
  maaslin.para <<-list(check= list(input_data = input.data[1,], 
      input_metadata = input.meta, 
      fixed_effects = c(fixed.effects),
      random_effects = c(block),
      reference = c(refs),
      max_significance = 0.05,
      min_abundance = 0.0,
      min_prevalence = 0.0,
      min_variance = 0.0,
      normalization = norm.method,
      transform = trans.method),
     test=list(input_data = input.data, 
        input_metadata = input.meta, 
        fixed_effects = c(fixed.effects),
        random_effects = c(block),
        reference = c(refs),
        max_significance = 0.05,
        min_abundance = 0.0,
        min_prevalence = 0.0,
        min_variance = 0.0,
        normalization = norm.method,
        transform = trans.method)
     )
    return(2)  
  }

}


############ use RserveMicro to perform MaasLin2
.prepare.maaslin2<-function(case,input_data,
                            input_metadata,
                            min_abundance = 0.0,
                            min_prevalence = 0.1,
                            min_variance = 0.0,
                            normalization = "TSS",
                            transform = "LOG",
                            analysis_method = model,
                            max_significance = 0.25,
                            random_effects = NULL,
                            fixed_effects = NULL,
                            correction = "BH",
                            standardize = TRUE,
                            cores = 1,
                            reference = NULL){
  require('data.table')
  require('dplyr')
if(case==1){
  input_data = maaslin.para$input_data
  if(exists("input_metadata",where = maaslin.para)){
    input_metadata = maaslin.para$input_metadata
    fixed_effects = maaslin.para$fixed_effects
  }
  reference = maaslin.para$reference
  max_significance = maaslin.para$max_significance
  min_abundance = maaslin.para$min_abundance
  min_prevalence = maaslin.para$min_prevalence
  min_variance = maaslin.para$min_variance
  normalization = maaslin.para$normalization
  transform = maaslin.para$transform
  
}else if(case==2){
  
  input_data = input_data
  input_metadata = maaslin.para$check$input_metadata
  fixed_effects = maaslin.para$check$fixed_effects
  random_effects = maaslin.para$check$random_effects
  reference = maaslin.para$check$reference
  max_significance = maaslin.para$check$max_significance
  min_abundance = maaslin.para$check$min_abundance
  min_prevalence = maaslin.para$check$min_prevalence
  min_variance = maaslin.para$check$min_variance
  normalization = maaslin.para$check$normalization
  transform = maaslin.para$check$transform
  
}else if(case==3){
  
  input_data = input_data
  input_metadata = maaslin.para$test$input_metadata
  fixed_effects = maaslin.para$test$fixed_effects
  random_effects = maaslin.para$test$random_effects
  reference = maaslin.para$test$reference
  max_significance = maaslin.para$test$max_significance
  min_abundance = maaslin.para$test$min_abundance
  min_prevalence = maaslin.para$test$min_prevalence
  min_variance = maaslin.para$test$min_variance
  normalization = maaslin.para$test$normalization
  transform = maaslin.para$test$transform
  
}else if(case==4){
  
  input_data = maaslin.para.noadj$input_data 
  input_metadata = maaslin.para.noadj$input_metadata
  fixed_effects = maaslin.para.noadj$fixed_effects
  reference = maaslin.para.noadj$reference
  max_significance = maaslin.para.noadj$max_significance
  min_abundance = maaslin.para.noadj$min_abundance
  min_prevalence = maaslin.para.noadj$min_prevalence
  min_variance = maaslin.para.noadj$min_variance
  normalization = maaslin.para.noadj$normalization
  transform = maaslin.para.noadj$transform
}

  # Allow for lower case variables
  normalization <- toupper(normalization);
  transform <- toupper(transform);
  analysis_method <- toupper(analysis_method);
  correction <- toupper(correction);
  
  #################################################################
  # Read in the data and metadata, create output folder, init log #
  #################################################################
  # if a character string then this is a file name, else it 
  # is a data frame
  if (is.character(input_data)) {
    data <- data.frame(data.table::fread(input_data, header = TRUE, sep = "\t"), row.names = 1);
    if (nrow(data) == 1) {
      # read again to get row name
      data <- read.table(input_data, header = TRUE, row.names = 1);
    }
  } else {
    data <- input_data;
  }
  if (is.character(input_metadata)) {
    metadata <- data.frame(data.table::fread(input_metadata, header = TRUE, sep = "\t"), row.names = 1);
    if (nrow(metadata) == 1) {
      metadata <- read.table(input_metadata, header = TRUE, row.names = 1);
    }
  } else {
    metadata <- input_metadata;
  }
  
  ###############################################################
  # Determine orientation of data in input and reorder to match #
  ###############################################################
  
  samples_row_row <- intersect(rownames(data), rownames(metadata))
  if (length(samples_row_row) > 0) {
    # this is the expected formatting so do not modify data frames
  } else {
    samples_column_row <- intersect(colnames(data), rownames(metadata));
    
    if (length(samples_column_row) == 0) {
      # modify possibly included special chars in sample names in metadata
      rownames(metadata) <- make.names(rownames(metadata));
      samples_column_row <- intersect(colnames(data), rownames(metadata));
    }
    
    if (length(samples_column_row) > 0) {
      # transpose data frame so samples are rows
      data <- as.data.frame(t(data));
    } else {
      samples_column_column <- intersect(colnames(data), colnames(metadata));
      if (length(samples_column_column) > 0) {
        data <- as.data.frame(t(data));
        metadata <- type.convert(as.data.frame(t(metadata)));
      } else {
        samples_row_column <- intersect(rownames(data), colnames(metadata));
        
        if (length(samples_row_column) == 0) {
          # modify possibly included special chars in sample names in data
          rownames(data) <- make.names(rownames(data));
          samples_row_column <- intersect(rownames(data), colnames(metadata));
        }
        
        if (length(samples_row_column) > 0) {
          metadata <- type.convert(as.data.frame(t(metadata)));
        } else {
          stop()
        }
      }
    }
  }
  
  # replace unexpected characters in feature names
  colnames(data) <- make.names(colnames(data))
  
  # check for samples without metadata
  extra_feature_samples <- setdiff(rownames(data), rownames(metadata))
  # check for metadata samples without features
  extra_metadata_samples <- setdiff(rownames(metadata), rownames(data))
  
  # get a set of the samples with both metadata and features
  intersect_samples <- intersect(rownames(data), rownames(metadata))
  
  # now order both data and metadata with the same sample ordering
  data <- data[intersect_samples, , drop = FALSE]
  metadata <- metadata[intersect_samples, , drop = FALSE]
  
  ###########################################
  # Compute the formula based on user input #
  ###########################################
  
  random_effects_formula <- NULL
  # use all metadata if no fixed effects are provided
  if (is.null(fixed_effects)) {
    fixed_effects <- colnames(metadata)
  } else {
    fixed_effects <- unlist(strsplit(fixed_effects, ",", fixed = TRUE))
    # remove any fixed effects not found in metadata names
    to_remove <- setdiff(fixed_effects, colnames(metadata))
    if (length(to_remove) > 0)
      fixed_effects <- setdiff(fixed_effects, to_remove)
    if (length(fixed_effects) == 0) {
      stop()
    }
  }
  
  if (!is.null(random_effects)) {
    random_effects <- unlist(strsplit(random_effects, ",", fixed = TRUE))
    
    # subtract random effects from fixed effects
    fixed_effects <- setdiff(fixed_effects, random_effects)
    
    # remove any random effects not found in metadata
    to_remove <- setdiff(random_effects, colnames(metadata))
    
    if (length(to_remove) > 0)
      random_effects <- setdiff(random_effects, to_remove)
    
    # create formula
    if (length(random_effects) > 0) {
      random_effects_formula_text <- paste("expr ~ (1 | ",
                                           paste(random_effects, ")", sep = '', collapse = " + (1 | "), sep = '');
      
      random_effects_formula <- tryCatch(as.formula(random_effects_formula_text),
                                         error = function(e)
                                           stop(paste("Invalid formula for random effects: ",
                                                      random_effects_formula_text)));
    }
  }
  
  # reduce metadata to only include fixed/random effects in formula
  effects_names <- union(fixed_effects, random_effects)
  metadata <- metadata[, effects_names, drop = FALSE]
  
  # create the fixed effects formula text
  formula_text <- paste("expr ~ ", paste(fixed_effects, collapse = " + "));
  formula <- tryCatch(as.formula(formula_text),
                      error = function(e)
                        stop(paste("Invalid formula.",
                                   "Please provide a different formula: ",
                                   formula_text)));
  #########################################################
  # Filter data based on min abundance and min prevalence #
  #########################################################
  
  # use ordered factor for variables with more than two levels
  # find variables with more than two levels
  if (is.null(reference)) {reference <- ","}
  
  for ( i in colnames(metadata) ) {
    mlevels <- unique(na.omit(metadata[,i]));
    numeric_levels <- grep('^-?[0-9.]+[eE+-]?', mlevels, value = T);
    if (all(c( ( length(mlevels[! (mlevels %in% c("UNK"))]) > 1 ),  # modification to allow setting reference when only two classes in metadata
         ( i %in% fixed_effects ),
         ( length(numeric_levels) == 0)))) {
      split_reference <- unlist(strsplit(reference, "[,;]"));
      if (! i %in% split_reference ) {
        stop(paste("Please provide the reference for the variable '",
                   i, "' which includes more than 2 levels: ",
                   paste(as.character(mlevels), collapse=", "), ".", sep=""));
      } else {
        ref <- split_reference[match(i,split_reference)+1];
        other_levels <- as.character(mlevels)[! as.character(mlevels) == ref];
        metadata[,i] <- factor(metadata[,i], levels=c(ref, other_levels));
      }
    }
  }       
  
  unfiltered_data <- data;
  unfiltered_metadata <- metadata;
  
  # require at least total samples * min prevalence values 
  # for each feature to be greater than min abundance
  total_samples <- nrow(unfiltered_data);
  min_samples <- total_samples * min_prevalence;
  
  # Filter by abundance using zero as value for NAs
  data_zeros <- unfiltered_data;
  data_zeros[is.na(data_zeros)] <- 0;
  filtered_data <- unfiltered_data[, colSums(data_zeros > min_abundance) > min_samples, drop = FALSE];
  total_filtered_features <- ncol(unfiltered_data) - ncol(filtered_data);
  filtered_feature_names <- setdiff(names(unfiltered_data), names(filtered_data));
  
  #################################
  # Filter data based on variance #
  #################################
  
  sds <- apply(filtered_data, 2, sd)
  variance_filtered_data <- filtered_data[, which(sds > min_variance), drop = FALSE]
  variance_filtered_features <- ncol(filtered_data) - ncol(variance_filtered_data)
  variance_filtered_feature_names <- setdiff(names(filtered_data), names(variance_filtered_data))
  filtered_data <- variance_filtered_data
  
  ######################
  # Normalize features #
  ######################
  
  filtered_data_norm <- normalizeFeatures(filtered_data, normalization = normalization)

  ################################
  # Standardize metadata, if set #
  ################################
  
  if (standardize) {
    metadata <- metadata %>% dplyr::mutate_if(is.numeric, scale)
  }
  
  ############################
  # Transform and run method #
  ############################
  
  # transform features
  filtered_data_norm_transformed <- transformFeatures(filtered_data_norm, transformation = transform)
  
  dat.in <- list(features=filtered_data_norm_transformed,
                 metadata=metadata, 
                 model=analysis_method, 
                 formula=formula, 
                 random_effects_formula=random_effects_formula,
                 correction=correction,
                 cores=cores);

  my.fun <- function(){
    features=dat.in$features
    metadata=dat.in$metadata
    model=dat.in$model
    formula=dat.in$formula
    random_effects_formula=dat.in$random_effects_formula
    correction=dat.in$correction
    cores=dat.in$cores
    
    if (is.null(formula))
      formula <-
        as.formula(paste(
          "expr ~ ", 
          paste(colnames(metadata), 
                collapse = "+")))
    
    if (!(is.null(random_effects_formula))) {
      formula <-
        paste(
          '. ~', 
          paste(all.vars(formula)[-1], collapse = ' + '), 
          '.', 
          sep = ' + ')
      formula <- update(random_effects_formula, formula)
    }
    
    #############################################################
    # Determine the function and summary for the model selected #
    #############################################################
    
    ################
    # Linear Model #
    ################
    
    if (model == "LM") {
      if (is.null(random_effects_formula)) {
        model_function <-
          function(formula, data, na.action) {
            return(glm(
              formula,
              data = data,
              family = 'gaussian',
              na.action = na.action
            ))
          }
        summary_function <- function(fit) {
          lm_summary <- summary(fit)$coefficients
          para <- as.data.frame(lm_summary)[-1, -3]
          para$name <- rownames(lm_summary)[-1]
          return(para)
        }
      } else {
        ranef_function <- lme4::ranef
        model_function <-
          function(formula, data, na.action) {
            return(lmerTest::lmer(
              formula, 
              data = data, 
              na.action = na.action))
          }
        summary_function <- function(fit) {
          lm_summary <- coef(summary(fit))
          para <- as.data.frame(lm_summary)[-1, -c(3:4)]
          para$name <- rownames(lm_summary)[-1]
          return(para)
        }
      }
    }
      ####################
        # Compound Poisson #
        ####################
      
        if (model == "CPLM") {
            if (is.null(random_effects_formula)) {
              model_function <- cplm::cpglm
              summary_function <- function(fit) {
                cplm_out <-
                  capture.output(
                    cplm_summary <- cplm::summary(fit)$coefficients)
                    para <- as.data.frame(cplm_summary)[-1, -3]
                    para$name <- rownames(cplm_summary)[-1]
                    logging::logdebug(
                        "Summary output\n%s", 
                        paste(cplm_out, collapse = "\n"))
                    return(para)
                    }
              } else {
                ranef_function <- glmmTMB::ranef
                model_function <-
                    function(formula, data, na.action) {
                        return(glmmTMB::glmmTMB(
                            formula,
                            data = data,
                            family=glmmTMB::tweedie(link = "log"),
                            ziformula = ~0,
                            na.action = na.action
                        ))
                    }
                summary_function <- function(fit) {
                  glmmTMB_summary <- coef(summary(fit))
                  para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
                  para$name <- rownames(glmmTMB_summary$cond)[-1]
                  return(para)
                }
              }
          }

#####################
        # Negative Binomial #
        #####################
        
        if (model == "NEGBIN") {
            if (is.null(random_effects_formula)) {
                model_function <- MASS::glm.nb
                summary_function <- function(fit) {
                  glm_summary <- summary(fit)$coefficients
                  para <- as.data.frame(glm_summary)[-1, -3]
                  para$name <- rownames(glm_summary)[-1]
                  return(para)
                  }
            } else {
              ranef_function <- glmmTMB::ranef
              model_function <-
                    function(formula, data, na.action) {
                        return(glmmTMB::glmmTMB(
                            formula,
                            data = data,
                            family=glmmTMB::nbinom2(link = "log"),
                            ziformula = ~0,
                            na.action = na.action
                        ))
                    }
                summary_function <- function(fit) {
                  glmmTMB_summary <- coef(summary(fit))
                  para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
                  para$name <- rownames(glmmTMB_summary$cond)[-1]
                  return(para)
                }
            }
          }
        
        ###################################
        # Zero-inflated Negative Binomial #
        ###################################
        
        if (model == "ZINB") {
          if (is.null(random_effects_formula)) {
            model_function <-
              function(formula, data, na.action) {
                return(pscl::zeroinfl(
                  formula,
                  data = data,
                  dist = "negbin",
                  na.action = na.action))
                }
            summary_function <- function(fit) {
              pscl_summary <- summary(fit)$coefficients$count
              rmidx <- grepl("Intercept",rownames(pscl_summary))|grepl("Log",rownames(pscl_summary))
              para <- as.data.frame(pscl_summary)[!rmidx, -3]
              para$name <- rownames(pscl_summary)[!rmidx] 
	            return(para)
	            }
            } else {
              ranef_function <- glmmTMB::ranef
              model_function <-
                function(formula, data, na.action) {
                  return(glmmTMB::glmmTMB(
                    formula,
                    data = data,
                    family=glmmTMB::nbinom2(link = "log"),
                    ziformula = ~1,
                    na.action = na.action))
                  }
              summary_function <- function(fit) {
                glmmTMB_summary <- coef(summary(fit))
                para <- as.data.frame(glmmTMB_summary$cond)[-1, -3]
                para$name <- rownames(glmmTMB_summary$cond)[-1]
                return(para)
              }
            }
    }
    #######################################
    # Init cluster for parallel computing #
    #######################################
    
    cluster <- NULL
    if (cores > 1)
    {
      cluster <- parallel::makeCluster(cores)
    }
    
    ##############################
    # Apply per-feature modeling #
    ##############################

    outputs <-lapply(seq_len(ncol(features)), function(x) {
        # Extract Features One by One
        featuresVector <- features[, x]
        
        dat_sub <- data.frame(expr = as.numeric(featuresVector), metadata) 
         fit <- tryCatch({
          fit1 <-
            model_function(
              formula, 
              data = dat_sub, 
              na.action = na.exclude)
        }, error = function(err) {
          fit1 <-
            try({
              model_function(
                formula, 
                data = dat_sub, 
                na.action = na.exclude)
            })
          return(fit1)
        })

        # Gather Output
        output <- list()
        if (all(!inherits(fit, "try-error"))) {
          output$para <- summary_function(fit)
          if (!(is.null(random_effects_formula))) {
            l <- ranef_function(fit)
            d<-as.vector(unlist(l))
            names(d)<-unlist(lapply(l, row.names))
            output$ranef<-d
          }
        }  else {
          output$para <-
            as.data.frame(matrix(NA, 
                                 nrow = ncol(metadata), ncol = 3))
          output$para$name <- colnames(metadata)
          if (!(is.null(random_effects_formula))) output$ranef <- NA
        }
 
        colnames(output$para) <- c('coef', 'stderr' , 'pval', 'name')
        output$para$feature <- colnames(features)[x]
        return(output)
      })
 
nafeat =which(!unlist(lapply(outputs,function(x) any(is.na(x[["para"]]$pval)|any(is.nan(x[["para"]]$pval))))))
outputs <- outputs[nafeat] 
    # stop the cluster
    if (!is.null(cluster))
      parallel::stopCluster(cluster)
    
    # bind the results for each feature
    paras <-
      do.call(rbind, lapply(outputs, function(x) {
        return(x$para)
      }))
    
    if (!(is.null(random_effects_formula))) {
      ranef <-
        do.call(rbind, lapply(outputs, function(x) {
          return(x$ranef)
        }))
      row.names(ranef) <- colnames(features) 
    }

nafeat = paras$feature[is.na(paras$pval)|is.nan(paras$pval)]
paras = paras[!paras$feature %in% nafeat,]

 
    ################################
    # Apply correction to p-values #
    ################################
    
    paras$qval <- as.numeric(p.adjust(paras$pval, method = correction))
    
    #####################################################
    # Determine the metadata names from the model names #
    #####################################################
    
    metadata_names <- colnames(metadata)
    # order the metadata names by decreasing length
    metadata_names_ordered <-
      metadata_names[order(
        nchar(metadata_names), decreasing = TRUE)]
    # find the metadata name based on the match 
    # to the beginning of the string
    extract_metadata_name <- function(name) {
      return(metadata_names_ordered[mapply(
        startsWith, 
        name, 
        metadata_names_ordered)][1])
    }
    paras$metadata <- unlist(lapply(paras$name, extract_metadata_name))

    # compute the value as the model contrast minus metadata
    paras$value <-
      mapply(function(x, y) {
        if (x == y)
          x
        else
          gsub(x, "", y)
      }, paras$metadata, paras$name)
 
    ##############################
    # Sort by decreasing q-value #
    ##############################
    
    paras <- paras[order(paras$qval, decreasing = FALSE), ]
    paras <-
      dplyr::select(
        paras,
        c('feature', 'metadata', 'value'),
        dplyr::everything())
    
    rownames(paras)<-NULL
    
    if (!(is.null(random_effects_formula))) {
      return(list("results" = paras, "ranef" = ranef))
    } else {
      return(list("results" = paras))
    }
  }
  dat.in <- list(features=filtered_data_norm_transformed,
                 metadata=metadata, 
                 model=analysis_method, 
                 formula=formula, 
                 random_effects_formula=random_effects_formula,
                 correction=correction,
                 cores=cores,
                 my.fun=my.fun,
                filtered_data_norm=filtered_data_norm);
  qs::qsave(dat.in, file="dat.in.qs");
  return(1);
}

.save.maaslin.res <- function(mbSetObj,case){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  dat.in <- qs::qread("dat.in.qs"); 
  filtered_data_norm <- dat.in$filtered_data_norm

  fit_data <- dat.in$my.res 
  fit_data$results$N <- apply(fit_data$results, 1, FUN = function(x)
    length(filtered_data_norm[, x[1]]));
  
  fit_data$results$N.not.zero <- apply(fit_data$results, 1, FUN = function(x)
    length(which(filtered_data_norm[, x[1]] > 0)));
  
  fit_data$input$data <- dat.in$filtered_data_norm_transformed
  fit_data$input$metadata <- dat.in$metadata
  fit_data$input$analysis_method <- dat.in$model
  fit_data$input$formula <- dat.in$formula
  fit_data$input$random_effects_formula <- dat.in$random_effects_formula
  fit_data$input$correction <- dat.in$correction
  if(case==4){
    mbSetObj$analSet$maaslin.noadj<-fit_data
  }else if(case !=2){
    mbSetObj$analSet$maaslin<-fit_data
  }
  .set.mbSetObj(mbSetObj)
  if(case==1|case==4){
      return(1)
  }else if(case==2){
    check.rank <- capture.output(fit_data, type=c("message"));
    
    if((length(grep("rank deficient", check.rank)) + length(grep("singular", check.rank))) > 0){
      # often random effects model matrix are rank deficient - check this way and return 
      # feedback that the experimental design does not support using a blocking factor.
      return(-2)
    }else{
      return(3)
    }
  }
}

PostProcessMaaslin <- function(mbSetObj,analysis.var,comp=NULL, thresh = 0.05,taxrank,is.norm,imgNm){

    mbSetObj <- .get.mbSetObj(mbSetObj);
    input.data<-maaslin.para$input_data
    res <- mbSetObj$analSet$maaslin$results
    inds <- !(res$feature %in% rownames(input.data)); # rownames that are all integers have "X" appended to front
    res$feature[inds] <- substring(res$feature[inds], 2);

    # get unadjusted results
    if((!mbSetObj$analSet$adj.bool) & (mbSetObj$analSet$block == "NA")){
        res.noadj <- res;
    } else {
        res.noadj <- mbSetObj$analSet$maaslin.noadj$results;
    }

    inds <- !(res.noadj$feature %in% rownames(input.data)) # rownames that are all integers have "X" appended to front
    res.noadj$feature[inds] <- substring(res.noadj$feature[inds], 2)

    # filter results to get only ones related to analysis var
    res <- res[res$metadata == analysis.var, ];

    # make res pretty
    res$coef <- signif(res$coef, digits = 3);
    res$stderr <- signif(res$stderr, digits = 3);
    res$pval <- signif(res$pval, digits = 3);
    res$qval <- signif(res$qval, digits = 3);

    if(  mbSetObj$analSet$analysis.type == "disc"){
        res <- res[res$value == comp, ];
        res.noadj <- res.noadj[res.noadj$value == comp, ];
        rownames(res) <- res$feature;
        rownames(res.noadj) <- res.noadj$feature;
        res <- res[ ,c("coef", "stderr", "pval", "qval")];
        res.noadj <- res.noadj[ ,c("coef", "stderr", "pval", "qval")];
        colnames(res) <- c("Log2FC", "St.Error", "P-value", "FDR");
        colnames(res.noadj) <- c("Log2FC", "St. Error", "P-value", "FDR");
    } else {
        rownames(res) <- res$feature;
        rownames(res.noadj) <- res.noadj$feature;
        res <- res[ ,c("coef", "stderr", "pval", "qval")];
        res.noadj <- res.noadj[ ,c("coef", "stderr", "pval", "qval")];
        colnames(res) <- c("Coefficient", "St.Error", "P-value", "FDR");
        colnames(res.noadj) <- c("Coefficient", "St.Error", "P-value", "FDR");
    }

    # write out/save results
    fileName <- "multifac_output.csv";
    fast.write(res, file = fileName);
    thresh <- as.numeric(thresh);
    # put results in mbSetObj, learn pattern of analysis set
    sigfeat <- rownames(res)[res$FDR < thresh];
    sig.count <- length(sigfeat);
    if(sig.count == 0){
        current.msg <<- "No significant features were identified using the given p value cutoff.";
    }else{
        current.msg <<- paste("A total of", sig.count, "significant features were identified!");
    }

    # process data for individual feature boxplot
    taxrank_boxplot <- taxrank;
    claslbl_boxplot <- as.factor(sample_data(mbSetObj$dataSet$norm.phyobj)[[analysis.var]]);
    nm_boxplot <- rownames(input.data);
    dat3t_boxplot <- as.data.frame(t(input.data),check.names=FALSE);
    colnames(dat3t_boxplot) <- nm_boxplot; 
    box_data <- dat3t_boxplot;
    box_data$class <- claslbl_boxplot;
    box_data$norm <- is.norm;


    # for graphical summary
    adj.mat <- res[, c("P-value", "FDR")]
    noadj.mat <- res.noadj[, c("P-value", "FDR")]

    colnames(adj.mat) <- c("pval.adj", "fdr.adj")
    colnames(noadj.mat) <- c("pval.no", "fdr.no")

    both.mat <- merge(adj.mat, noadj.mat, by = "row.names")
    both.mat$pval.adj <- -log10(both.mat$pval.adj)
    both.mat$fdr.adj <- -log10(both.mat$fdr.adj)
    both.mat$pval.no <- -log10(both.mat$pval.no)
    both.mat$fdr.no <- -log10(both.mat$fdr.no)

    rownames(both.mat) = both.mat[,1]

    # for plotting adjp vs p
    jsonNm <- gsub(".png", ".json", imgNm);
    jsonObj <- RJSONIO::toJSON(both.mat);
    sink(jsonNm);
    cat(jsonObj);
    sink();

    mbSetObj$analSet$cov.mat <- both.mat; 
    mbSetObj$analSet$multiboxdata <- box_data;
    mbSetObj$analSet$sig.count <- sig.count;
    mbSetObj$analSet$cov <- list();
    mbSetObj$analSet$cov$resTable <- mbSetObj$analSet$resTable <- res;
    mbSetObj$analSet$maas.resnoadj <- res.noadj;

    return(.set.mbSetObj(mbSetObj))
}


PlotCovariateMapMas <- function(mbSetObj, theme="default", imgName="NA", format="png", dpi=72, interactive=F){
  save.image("cov.RData");
  mbSetObj <- .get.mbSetObj(mbSetObj);
  both.mat <- mbSetObj$analSet$cov.mat
  both.mat <- both.mat[order(-both.mat[,"pval.adj"]),]
  logp_val <- -log10(mbSetObj$paramSet$cov$p.lvl)
  topFeature <- 5;
  if(nrow(both.mat) < topFeature){
    topFeature <- nrow(both.mat);
  }
  
  fileName <- paste0(imgName, ".", format);
  mbSetObj$imgSet$covAdj <- fileName;
  
  width <- 8;
  height <- 8.18;
  
  library(plotly)
  threshold <- logp_val               
  
  both.mat$category <- with(both.mat, case_when(
    pval.no > threshold & pval.adj > threshold ~ "Significant in both",
    pval.no > threshold & pval.adj <= threshold ~ "Significant in pval.no only",
    pval.adj > threshold & pval.no <= threshold ~ "Significant in pval.adj only",
    TRUE ~ "Non-significant"
  ))
  
  # Define a list or data frame mapping categories to properties
  category_properties <- data.frame(
    category = c("Significant in both", "Significant in pval.no only", 
                 "Significant in pval.adj only", "Non-significant"),
    color = c('#6699CC', '#94C973', '#E2808A', 'grey'),
    name = c("Significant", "Non-Sig. after adjustment", 
             "Sig. after adjustment", "Non-Significant")
  )
  
  p <- ggplot(both.mat, aes(x = pval.no, y = pval.adj, color = category, text = paste("Feature:", Row.names, 
                                                                                               "<br>Adjusted Pval:", signif(10^(-pval.adj), 4), 
                                                                                               "<br>Non-adjusted Pval:", signif(10^(-pval.no), 4)))) +
    geom_point(alpha = 0.5) +
    scale_color_manual(values = setNames(category_properties$color, category_properties$category), name="") +
    labs(x = "-log10(P-value): no covariate adjustment", y = "-log10(P-value): adjusted") +
    theme_minimal() +
    theme(legend.title = element_blank())
  

    Cairo::Cairo(file = fileName, unit="in", dpi=dpi, width=width, height=height, type=format);    
    print(p)
    dev.off()
  

  if(interactive){
    library(plotly);
    ggp_build <- layout(ggplotly(p,width = 800, height = 600, tooltip = c("text")), autosize = FALSE, margin = mbSetObj$imgSet$margin.config)
    .set.mbSetObj(mbSetObj)
    return(ggp_build);
  }else{
    return(.set.mbSetObj(mbSetObj));
  }
}
