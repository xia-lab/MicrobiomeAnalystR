################################################
###raw 16s data processing
### Author: Zhiqiang Pang, Yao Lu
########


MessageOutput <- function(msg, eol = "\n"){
  msg <- gsub("\"|\"|\\\\|\\[[0-9]\\] ", "", msg)
  write.table(msg, file = "seq_process_details.txt", quote = F, row.names = F, col.names = F, append = T, eol = eol)
}

PerformSeqCheck <- function(home_dir = ""){
  require(dada2)
  
  path <- home_dir
  list.files(path)
  fnFs <- sort(list.files(paste0(path, "/upload"), pattern="_R1.fastq", full.names = TRUE))
  fnRs <- sort(list.files(paste0(path, "/upload"), pattern="_R2.fastq", full.names = TRUE))
  
  if(length(fnRs)>=2){
    p1 <- plotQualityProfile(c(fnFs[1:2], fnRs[1:2]))
  } else if(length(fnRs)==1){
    p1 <- plotQualityProfile(c(fnFs[1], fnRs[1]))
  } else {
    p1 <- plotQualityProfile(fnFs[1:2])
  }
  qs::qsave(p1, "diagnotics_plot_src.qs");
  Cairo::Cairo(1600, 1250,file = paste0("diagnotics.png"),dpi = 180,bg = "white")
  print(p1)
  dev.off()
  return(1)
}

PerformSeqImport <- function(home_dir = ""){
  write.table(1.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  MessageOutput(paste0("<b>Loading R pacakge dada2, version: ", packageVersion("dada2"), " ...</b>"), " ")
  require(dada2)
  require(ggplot2)
  MessageOutput("OK, loaded!")
  if(!exists("dataObj")){
    dataObj <<- list()
  }
  ####preprocess using paired0-sequence as example
  
  path <- home_dir
  #list.files(path)
  fnFs <- sort(list.files(paste0(path, "/upload"), pattern="_R1.fastq", full.names = TRUE))
  fnRs <- sort(list.files(paste0(path, "/upload"), pattern="_R2.fastq", full.names = TRUE))

  # Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
  #sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
  sample.names <- unique(gsub("_R1.fastq|_R2.fastq|_R1.fq|_R2.fq", "", basename(fnFs)))
  MessageOutput("Step 1: Perform Data Import and Quality Check...")
  write.table(3.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  # Place filtered files in filtered/ subdirectory
  filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
  filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
  names(filtFs) <- sample.names
  names(filtRs) <- sample.names
  if(length(fnRs) == 0){
    dataObj$files <- data.frame(fnFs = fnFs,
                                filtFs = filtFs)
  } else {
    dataObj$files <- data.frame(fnFs = fnFs,
                                fnRs = fnRs,
                                filtFs = filtFs, 
                                filtRs = filtRs)
  }

  dataObj$sample.names <- sample.names;
  dataObj <<- dataObj;
  MessageOutput("OK, done!")
  write.table(4.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  return(1)
}

PerformSeqProcessing <- function(){
  MessageOutput("Step 2: Perform Sequencing data processing... ")
  write.table(5.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  
  params <- dataObj$params
  df_files <- dataObj$files;
  fnFs <- df_files$fnFs
  filtFs <- df_files$filtFs
  fnRs <- df_files$fnRs
  filtRs <- df_files$filtRs
  
  if(file.exists("/Users/lzy/Documents/examples_data_microbiomeanalyst/")){
    load.lib <- "/Users/lzy/Documents/examples_data_microbiomeanalyst/demo_16s/"
  }else{
    load.lib <-  "/home/glassfish/projects/microbiome_tax/"
  } 

  if(!exists("dataObj")){
    stop("SERVER ERROR <--- CHECK HERE")
  }
  write.table(6.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
 
  if(params$is_paired){
    truncLen<-c(params$forward_len, params$reverse_len);
    maxEE <- c(params$forward_maxEE, params$reverse_maxEE);
  } else {
    truncLen<-c(params$forward_len);
    maxEE <- c(params$forward_maxEE);
  }
  maxN <- params$maxN;
  truncQ <- params$trunQ;
  rm.phix <- params$rm_phix;
  trimLeft <- params$trimLeft;
  trimRight <- params$trimRight;
  write.table(8.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  dataObj$sample.names -> sample.names;
   MessageOutput("Start Sequencing data filtering and triming2 ... ")

  if(params$is_paired){
    out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=truncLen, 
                 trimLeft = trimLeft, trimRight = trimRight,
                 maxN=maxN, maxEE=maxEE, truncQ=truncQ, rm.phix=rm.phix,
                 compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
  } else {
    out <- filterAndTrim(fnFs, filtFs, truncLen=truncLen, 
                         trimLeft = trimLeft, trimRight = trimRight,
                         maxN=maxN, maxEE=maxEE, truncQ=truncQ, rm.phix=rm.phix,
                         compress=TRUE, multithread=TRUE) # On Windows set multithread=FALSE
  }
   MessageOutput("Start Sequencing data filtering and triming2 ... ")
  if(all(out[,"reads.out"]==0)){
    MessageOutput("ERROR! No reads passed the filter. Please revisit your filtering parameters.")
    stop("ERROR")
  }
  write.table(15.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  MessageOutput("OK, done!")
  ####Dereplicate
  if(params$is_paired){
    msg_d1 <- capture.output(derepF1 <- derepFastq(filtFs, verbose=TRUE), type = "message")
    sapply(msg_d1, MessageOutput)
    msg_d2 <- capture.output(derepR1 <- derepFastq(filtRs, verbose=TRUE), type = "message")
    sapply(msg_d2, MessageOutput)
  } else {
    msg_d1 <- capture.output(derepF1 <- derepFastq(filtFs, verbose=TRUE), type = "message")
    sapply(msg_d1, MessageOutput)
  }
  
  MessageOutput("OK, done!")
  write.table(20.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  ##Learn the error rates
  MessageOutput("Step 4: Perform Sequencing data dereplicating ... ")
  if(params$is_paired){
    msg_d3 <- capture.output(errF <- learnErrors(derepF1, multithread=TRUE)) 
    msg_d4 <- capture.output(errR <- learnErrors(derepR1, multithread=TRUE))
    sapply(msg_d3, MessageOutput)
    sapply(msg_d4, MessageOutput)
  } else {
    msg_d3 <- capture.output(errF <- learnErrors(derepF1, multithread=TRUE))
    sapply(msg_d3, MessageOutput)
  }

  write.table(30.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  MessageOutput("Plotting error frequnecy bar plots ... ")
  if(params$is_paired){
    p1 <- plotErrors(errF, nominalQ=TRUE)
    Cairo::Cairo(file = "error_images_f.png",  unit="in", dpi=72, width=12, height=10, type="png", bg="white");
    print(p1)
    dev.off();
    p2 <- plotErrors(errR, nominalQ=TRUE)
    Cairo::Cairo(file = "error_images_r.png",  unit="in", dpi=72, width=12, height=10, type="png", bg="white");
    print(p2)
    dev.off();
  } else {
    p0 <- plotErrors(errF, nominalQ=TRUE)
    Cairo::Cairo(file = "error_images_f.png",  unit="in", dpi=72, width=12, height=10, type="png", bg="white");
    print(p0)
    dev.off();
  }

  MessageOutput("OK, done!")
  
  ###Denoise
  MessageOutput("Step 5: Perform Sequencing data denoising ... ")
  if(params$is_paired){
    msg_d5 <- capture.output(dadaFs <- dada(derepF1, err=errF, multithread=TRUE))
    msg_d6 <- capture.output(dadaRs <- dada(derepR1, err=errR, multithread=TRUE))
    write.table(50.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
    sapply(msg_d5, MessageOutput)
    sapply(msg_d6, MessageOutput)
  } else {
    msg_d5 <- capture.output(dadaFs <- dada(derepF1, err=errF, multithread=TRUE))
    write.table(50.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
    sapply(msg_d5, MessageOutput)
  }

  MessageOutput("OK, done!")
  
  ####merge
  if(params$is_paired){
    MessageOutput("Step 6: Perform Sequencing data merging ... ")
    msg_d7 <- capture.output(mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE), type = "message")
    seqtab <- makeSequenceTable(mergers)
    sapply(msg_d7, MessageOutput)
    MessageOutput("OK, done!")
    write.table(60.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)

  } else {
    # code to processs single-ended data
    MessageOutput("Step 6: Perform Sequencing data generating ... ")
    seqtab <- makeSequenceTable(dadaFs);
    MessageOutput("OK, done!")
    write.table(60.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  }
  
  ####
  #####Remove chimeras
  MessageOutput("Step 7: Perform Sequencing chimeras removal ... ")
  msg_d8 <- capture.output(seqtab.nochim <- removeBimeraDenovo(seqtab,
                                                               method="consensus", 
                                                               multithread=TRUE, 
                                                               verbose=TRUE), type = "message")
  sapply(msg_d8, MessageOutput)
  MessageOutput("OK, done!")
  write.table(70.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  ##### assign taxonomy
  taxadb = params$taxadb
  MessageOutput("Step 8: Perform Sequencing taxonomy assignment ... ")
  if(taxadb=="gg"){
   taxa <- assignTaxonomy(seqtab.nochim, paste0(load.lib,"gg_13_8_train_set_97.fa.gz"),
                         multithread=TRUE)
  write.table(80.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)

  }else if(taxadb=="silva_fg"){
taxa <- assignTaxonomy(seqtab.nochim, paste0(load.lib,"silva_132.18s.99_rep_set.dada2.fa.gz"),
                         multithread=TRUE)
  write.table(80.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
   }else if(taxadb=="unite"){
  taxa <- assignTaxonomy(seqtab.nochim, paste0(load.lib,"sh_general_release_dynamic_10.10.2017.fasta"),
                         multithread=TRUE)
  write.table(80.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  }else{
  if(taxadb=="silva_py"){
     dbpath1 = "silva_nr99_v138.1_train_set.fa.gz"  
     dbpath2 = "silva_species_assignment_v138.1.fa.gz"  
   }else if(taxadb=="rdp_py"){
     dbpath1 = "rdp_train_set_18.fa.gz"  
     dbpath2 = "rdp_species_assignment_18.fa.gz"  
  }

  taxa <- assignTaxonomy(seqtab.nochim, paste0(load.lib,dbpath1),
                         multithread=TRUE)
  write.table(80.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  taxa <- addSpecies(taxa,  paste0(load.lib,dbpath2))

}
 MessageOutput("OK, done!")
  getN <- function(x) sum(getUniques(x))
  if(params$is_paired){
    track <- cbind(out, sapply(dadaFs, getN), sapply(mergers, getN), rowSums(seqtab), rowSums(seqtab.nochim));
    colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim");
  } else {
    track <- cbind(out, sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim));
    colnames(track) <- c("input", "filtered", "denoised", "tabled", "nonchim");
  }
  rownames(track) <- sample.names
  
  ##### Finalize results
  write.table(90.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  dataObj$res <- list();
  dataObj$res$seqtab <- t(seqtab);
  dataObj$res$seqtab.nochim <- t(seqtab.nochim);
  dataObj$res$track <- track;
  dataObj$res$taxa <- taxa;
  dataObj <<- dataObj;
  save(dataObj, file = "dataObj_completed.rda")
  write.table(100.0, file = "log_progress.txt", quote = F, row.names = F, col.names = F, append = F)
  MessageOutput("Everything has been finished successfully!")
  return(1)
}

### Functional utilities

checkParing <- function(){
  allUPFiles <- list.files("upload/", pattern = "fastq|fq|fastq.gz");
  if(length(allUPFiles) == 0){
    return(0)
  }
  #print(allUPFiles)
  allNames <- vapply(allUPFiles, function(x){
    gsub("_R1.fastq|_R2.fastq|_R1.fq|_R2.fq|_R1.fastq.gz|_R2.fastq.gz", "", x)
  }, FUN.VALUE = character(1L))
  allNames <- unname(allNames)
  if(length(unique(allNames)) == length(allNames)/2){
    return(2)
  }
  # return 0 means empty; return 1 means single; return 2 means paired.
  return(1)
}

GetSanityVec <- function(paired){
  # this function is checking the sainty of all fastq files
  allUPFiles <- list.files("upload/", pattern = "fastq|fq", full.names = TRUE);
  Sanity_vec <- vector();
  # 1. is real fastq file or not
  valid_vec <- vapply(allUPFiles, validate_fastq, FUN.VALUE = logical(1L));
  # 2. check file size (MB)
  size_vec <- vapply(allUPFiles, file.size, FUN.VALUE = double(1L))/(1024^2);
  size_vec <- signif(size_vec, 2)
  # 3. calculate how many reads in each file
  read_vec <- vapply(allUPFiles, calculate_Reads, FUN.VALUE = integer(1L));
  # 4. Pairing the name if they are paired-ended (paired == 2)
  if(file.exists("meta_info.qs")){
    meta_df <- qs::qread("meta_info.qs")
  } else if(file.exists("meta_its.qs")){
    meta_df <-qs::qread("meta_its.qs")
  }else{
    meta_df <- data.frame()
  }
  if(paired == 2){
    allNames <- vapply(allUPFiles, function(x){
      gsub("_R1.fastq|_R2.fastq|_R1.fq|_R2.fq", "", x)
    }, FUN.VALUE = character(1L));
    allNames <- unname(allNames);
    allNames <- unique(basename(allNames));
    
    Sanity_vec <- sapply(allNames, function(x){
      nmfp <- paste0(x, "_R1");
      nmrp <- paste0(x, "_R2");
      idxf <- grepl(nmfp, allUPFiles);
      idxr <- grepl(nmrp, allUPFiles);
      
      nmf <- basename(allUPFiles[idxf]);
      nmr <- basename(allUPFiles[idxr]);
      sizef <- unname(size_vec[idxf]);
      sizer <- unname(size_vec[idxr]);
      readf <- unname(read_vec[idxf]);
      readr <- unname(read_vec[idxr]);
      valdf <- unname(valid_vec[idxf]);
      valdr <- unname(valid_vec[idxr]);
      if(nrow(meta_df)!=0){
        group <- meta_df[grep(nmfp, meta_df[,1]),2];
      } else {
        group <- "NA";
      }
      paste(nmf, nmr, sizef, sizer, readf, readr, valdf, valdr, group, sep = " || ")
    })
  } else {
    allNames <- unique(basename(allUPFiles));
    Sanity_vec <- sapply(allNames, function(x){
      idxf <- grepl(x, allUPFiles);
      
      nmf <- basename(allUPFiles[idxf]);
      nmf <- gsub(".fastq", "", nmf)
      sizef <- unname(size_vec[idxf]);
      readf <- unname(read_vec[idxf]);
      valdf <- unname(valid_vec[idxf]);
      if(nrow(meta_df)!=0){
        group <- meta_df[grep(nmf, meta_df[,1]),2];
      } else {
        group <- "NA";
      }
      paste(nmf, "x" , sizef, "x" , readf, "x" , valdf, "x" , group, sep = " || ")
    })
  }
  return(Sanity_vec)
}

validate_fastq <- function(file){
  bool1 <- bool2 <- bool3 <- bool4 <- FALSE;
  con <- file(file,"r")
  lines <- readLines(con,n=4)
  bool1 <- grepl("^@[A-Za-z]+", lines[1]);
  bool2 <- grepl("^[A T C G N]*$", lines[2]);
  bool3 <- grepl("\\+", lines[3])
  bool4 <- grepl("[A-Za-z 0-9]+", lines[4])
  close(con)
  return(bool1 & bool2 & bool3 & bool4)
}

calculate_Reads <- function(file){
  return(as.integer(ceiling(R.utils::countLines.default(file)/4)))
}

read_seq_len <- function(){
  
  allUPFiles <- list.files("upload/", pattern = "fastq|fq", full.names = TRUE);
  # select first and last file to check the read len
  ffst <- allUPFiles[1];
  flst <- tail(allUPFiles,1)[1]
  # read top 12 lines (3 reads to as the sample)
  con <- file(ffst,"r")
  lines <- readLines(con,n=12)[c(2,6,10)]
  num1 <- unname(vapply(lines, nchar, FUN.VALUE = integer(1L)))
  close(con)
  
  con <- file(flst,"r")
  lines <- readLines(con,n=12)[c(2,6,10)]
  num2 <- unname(vapply(lines, nchar, FUN.VALUE = integer(1L)))
  close(con)
  return(ceiling(mean(c(num1, num2))))
}

composeParamFun <- function(trimLeft, trimRight, 
                            forward_len, reverse_len, 
                            maxN, 
                            forward_maxEE, reverse_maxEE,
                            trunQ,
                            minQ,
                            rm_phix,
                            is_paired,taxabd){
  if(!exists("dataObj")){
    dataObj <<- list()
  }
  funs <- params <- list()
  
  # Params section
  params$trimLeft <- trimLeft;
  params$trimRight <- trimRight;
  params$forward_len <- forward_len;
  params$reverse_len <- reverse_len;
  params$maxN <- maxN;
  params$forward_maxEE <- forward_maxEE;
  params$reverse_maxEE <- reverse_maxEE;
  params$trunQ <- trunQ;
  params$minQ <- minQ;
  params$is_paired <- is_paired;
  if(rm_phix == 0){
    params$rm_phix <- FALSE;
  } else {
    params$rm_phix <- TRUE;
  }
  params$taxadb <- taxabd;
  dataObj$params <- params;
  
  # meta section
  if(file.exists("meta_info.qs")){
    meta_df <- qs::qread("meta_info.qs")
  } else {
    meta_df <- NA;
  }
  dataObj$meta_info <- meta_df;
  
  # function section
  funs$MessageOutput <- MessageOutput;
  funs$PerformSeqCheck <- PerformSeqCheck;
  funs$PerformSeqImport <- PerformSeqImport;
  funs$PerformSeqProcessing <- PerformSeqProcessing;
  
  dataObj$funs <- funs;
  dataObj <<- dataObj;
  save(dataObj, file = "dataObj_param.rda")
  return(1);
}

sweaveRScript4exec <- function(users.path){
  str <- paste0('library(dada2)');
  
  # Set working dir & funcs to be used
  str <- paste0(str, ";\n", "setwd(\'",users.path,"\')");
  str <- paste0(str, ";\n", "load('dataObj_param.rda')");
  str <- paste0(str, ";\n", "dataObj <<- dataObj");
  str <- paste0(str, ";\n", "MessageOutput <- dataObj[['funs']][['MessageOutput']]");
  str <- paste0(str, ";\n", "PerformSeqCheck <- dataObj[['funs']][['PerformSeqCheck']]");
  str <- paste0(str, ";\n", "PerformSeqImport <- dataObj[['funs']][['PerformSeqImport']]");
  str <- paste0(str, ";\n", "PerformSeqProcessing <- dataObj[['funs']][['PerformSeqProcessing']]");
  
  ## Construct the exec pipeline
  str <- paste0(str, ';\n',  "PerformSeqImport('.')")
  str <- paste0(str, ';\n',  "PerformSeqProcessing()")
  
  # sink command for running
  sink("ExecuteRaw16Seq.R");
  cat(str);
  sink();
  return(1);
}

sweaveBash4exec <- function(users.path){

  ## Prepare Configuration script for slurm running
  conf_inf <- "#!/bin/bash\n#\n#SBATCH --job-name=16S_Processing\n#\n#SBATCH --ntasks=1\n#SBATCH --time=600:00\n#SBATCH --mem-per-cpu=5G\n#SBATCH --cpus-per-task=2\n"
  
  ## Prepare R script for running
  # need to require("dada2")
  str <- paste0('library(dada2)');
  
  # Set working dir & funcs to be used
  str <- paste0(str, ";\n", "setwd(\'",users.path,"\')");
  str <- paste0(str, ";\n", "load('dataObj_param.rda')");
  str <- paste0(str, ";\n", "dataObj <<- dataObj");
  str <- paste0(str, ";\n", "MessageOutput <- dataObj[['funs']][['MessageOutput']]");
  str <- paste0(str, ";\n", "PerformSeqCheck <- dataObj[['funs']][['PerformSeqCheck']]");
  str <- paste0(str, ";\n", "PerformSeqImport <- dataObj[['funs']][['PerformSeqImport']]");
  str <- paste0(str, ";\n", "PerformSeqProcessing <- dataObj[['funs']][['PerformSeqProcessing']]");
  
  ## Construct the exec pipeline
  str <- paste0(str, ';\n',  "PerformSeqImport('.')")
  str <- paste0(str, ';\n',  "PerformSeqProcessing()")
  
  # sink command for running
  sink("ExecuteRaw16Seq.sh");
  
  cat(conf_inf);
  cat(paste0("\nsrun R -e \"\n", str, "\n\""));
  
  sink();
  
  return(1)
}

sweaveBash4execPro <- function(users.path, isfromGoogle = FALSE, source_path){

  ## Prepare Configuration script for slurm running
  conf_inf <- paste0("#!/bin/bash\n#\n#SBATCH --job-name=16S_Processing\n#\n#SBATCH --ntasks=2\n#SBATCH --time=600:00\n#SBATCH --mem-per-cpu=5G\n#SBATCH --cpus-per-task=2\n#SBATCH --output=", users.path, "/seq_process_details.txt\n")
  
  ## to check if need to download files from google drive at the begining
  if(isfromGoogle){
    download_script <- paste0("Rscript --vanilla ", source_path, "/XiaLabPro/R/download_googledrive.R ", users.path)
  } else {
    download_script <- "";
  }
  
  
  ## Prepare R script for running
  # need to require("dada2")
  str <- paste0('library(dada2)');
  
  # Set working dir & funcs to be used
  str <- paste0(str, ";\n", "setwd(\'",users.path,"\')");
  str <- paste0(str, ";\n", "load('dataObj_param.rda')");
  str <- paste0(str, ";\n", "dataObj <<- dataObj");
  str <- paste0(str, ";\n", "MessageOutput <- dataObj[['funs']][['MessageOutput']]");
  str <- paste0(str, ";\n", "PerformSeqCheck <- dataObj[['funs']][['PerformSeqCheck']]");
  str <- paste0(str, ";\n", "PerformSeqImport <- dataObj[['funs']][['PerformSeqImport']]");
  str <- paste0(str, ";\n", "PerformSeqProcessing <- dataObj[['funs']][['PerformSeqProcessing']]");
  
  ## Construct the exec pipeline
  str <- paste0(str, ';\n',  "PerformSeqImport('.')")
  str <- paste0(str, ';\n',  "PerformSeqProcessing()")
  
  # sink command for running
  sink("ExecuteRaw16Seq.sh");
  
  cat(conf_inf);
  if(download_script!=""){
    cat("\n\n", download_script, "\n\n")
  }
  cat(paste0("\nR -e \"\n", str, "\n\""));
  
  sink();
  
  return(1)
}


ReadRawMeta<-function(fileName){
  if(grepl(".txt", fileName, fixed=T)){
    tbl=read.table(fileName,header=TRUE, stringsAsFactors = F, sep = "\t");
  }else if(grepl(".csv", fileName, fixed=T)){
    tbl = read.csv(fileName,header=TRUE, stringsAsFactors = F);
  }else{
    print("something very wrong : ===> wrongfiletype")
  }
  
  rawFileNms <-as.vector(tbl[,1])
  rawClassNms <-as.vector(tbl[,2])
  #rawFileNms <- sapply(strsplit(rawFileNms, "\\."), function(x) paste0(head(x,-1), collapse="."));
  rawFileNms <- sapply(rawFileNms, function(x){sub(".zip$|.fastq.gz$|.fq.gz$", "", x)})
  clsTable = table(rawClassNms)
  #check replicate number
  clsTypes = names(table(rawClassNms))
  
  if(length(rawFileNms)==0 | length(rawClassNms)==0){
    return(0);
  }
  
  meta_df <- data.frame(fileNms = rawFileNms,
                        class = rawClassNms);
  qs::qsave(meta_df, file = "meta_info.qs");

  rawFileNms<<-rawFileNms;
  rawClassNms<<-rawClassNms;
  return(1);
}

#' @noRd
GetRawFileNms <- function(){
  return(rawFileNms)
}

#' @noRd
GetRawClassNms <- function(){
  return(rawClassNms)
}

summaryRes <- function(mbSetObj = NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  load("dataObj_completed.rda");
  mbSetObj$dataObj <- dataObj;
  sampleMsg <- OTU_total_Msg <- Uni_taxa_Msg <- NonChim_Msg <- DeNoise_Msg <- "";
  if(length(dataObj[["res"]])>0){
    sampleMsg <- paste0("This job contains ", length(dataObj[["sample.names"]]), " samples.");
    OTU_total_Msg <- paste0("Total of ", nrow(dataObj[["res"]][["seqtab"]]), " OTUs and ", nrow(dataObj[["res"]][["seqtab.nochim"]]), " non-chimeric OTUs found.");
    Uni_taxa_Msg <- paste0(length(unique(dataObj[["res"]][["taxa"]][,2])), " phyla, ",
                           length(unique(dataObj[["res"]][["taxa"]][,3])), " classes, ",
                           length(unique(dataObj[["res"]][["taxa"]][,4])), " orders, ",
                           length(unique(dataObj[["res"]][["taxa"]][,5])), " families, ",
                           length(unique(dataObj[["res"]][["taxa"]][,6])), " genera and ",
                           length(unique(dataObj[["res"]][["taxa"]][,7])), " species have been found.");

    NonChim_Msg <- paste0(sum(dataObj[["res"]][["track"]][,ncol(dataObj[["res"]][["track"]])]), 
                          " (", signif(sum(dataObj[["res"]][["track"]][,ncol(dataObj[["res"]][["track"]])])/sum(dataObj[["res"]][["track"]][,1])*100,4), "%)",
                          " non-chimeirc OTUs found from all files.");
    DeNoise_Msg <- paste0(sum(dataObj[["res"]][["track"]][,3]), 
                          " (", signif(sum(dataObj[["res"]][["track"]][,3])/sum(dataObj[["res"]][["track"]][,1])*100,4), "%)",
                          " OTUs found from all files after de-noising.");
    
  }
  mbSetObj[["analSet"]]$summsg <- c(
    sampleMsg,
    OTU_total_Msg,
    NonChim_Msg,
    DeNoise_Msg,
    Uni_taxa_Msg)

  if (.on.public.web) {
    .set.mbSetObj(mbSetObj)
    return(c(
      sampleMsg,
      OTU_total_Msg,
      NonChim_Msg,
      DeNoise_Msg,
      Uni_taxa_Msg
    ))
  } else {
    print("Here is the summary of this processing: ")
    print(c(
      sampleMsg,
      OTU_total_Msg,
      NonChim_Msg,
      DeNoise_Msg,
      Uni_taxa_Msg
    ))
    return(.set.mbSetObj(mbSetObj))
  }
}

exportSampleNMs <- function(mbSetObj = NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataObj -> dataObj;
  if(is.null(dataObj)){
    load("dataObj_completed.rda");
    mbSetObj$dataObj <- dataObj;
  } else {
    dataObj <- mbSetObj$dataObj;
  }
  
  if (.on.public.web) {
    .set.mbSetObj(mbSetObj)
    return(dataObj[["sample.names"]]);
  } else {
    print(dataObj[["sample.names"]]);
    return(.set.mbSetObj(mbSetObj))
  }
}

exportSampleTrackTable <- function(mbSetObj = NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataObj -> dataObj;
  if(is.null(dataObj)){
    load("dataObj_completed.rda");
    mbSetObj$dataObj <- dataObj;
  } else {
    dataObj <- mbSetObj$dataObj;
  }
  if (.on.public.web) {
    .set.mbSetObj(mbSetObj)
    return(dataObj[["res"]][["track"]]);
  } else {
    print(dataObj[["res"]][["track"]]);
    return(.set.mbSetObj(mbSetObj))
  }
  
}

exportOTUtaxaInfo <- function(mbSetObj = NA, taxaLevel = 1){
  # taxaLevel: 0, OTU; 1, kingdom; 2, phylum; 3. class; 4, order; 5, family; 6, genera; 7, species
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataObj -> dataObj;
  if(is.null(dataObj)){
    load("dataObj_completed.rda");
    mbSetObj$dataObj <- dataObj;
  } else {
    dataObj <- mbSetObj$dataObj;
  }
  resv <- vector();
  if(taxaLevel == 0){
    resv <- rownames(dataObj[["res"]][["taxa"]])
  } else if(taxaLevel == 1) {
    resv <- dataObj[["res"]][["taxa"]][,1]
  } else if(taxaLevel == 2) {
    resv <- dataObj[["res"]][["taxa"]][,2]
  } else if(taxaLevel == 3) {
    resv <- dataObj[["res"]][["taxa"]][,3]
  } else if(taxaLevel == 4) {
    resv <- dataObj[["res"]][["taxa"]][,4]
  } else if(taxaLevel == 5) {
    resv <- dataObj[["res"]][["taxa"]][,5]
  } else if(taxaLevel == 6) {
    resv <- dataObj[["res"]][["taxa"]][,6]
  } else if(taxaLevel == 7) {
    resv <- dataObj[["res"]][["taxa"]][,7]
  } 
  resv <- unname(resv);
  resv[is.na(resv)] <- "NA"
  if (.on.public.web) {
    .set.mbSetObj(mbSetObj)
    return(resv);
  } else {
    print(resv);
    return(.set.mbSetObj(mbSetObj))
  }
}

generateResFigures <- function(mbSetObj = NA){
  mbSetObj <- .get.mbSetObj(mbSetObj);
  mbSetObj$dataObj -> dataObj;
  if(is.null(dataObj)){
    load("dataObj_completed.rda");
    mbSetObj$dataObj <- dataObj;
  } else {
    dataObj <- mbSetObj$dataObj;
  }
  require(phyloseq);
  require(ggplot2);
  seqtab.nochim <- cbind(rownames(dataObj[["res"]][["seqtab.nochim"]]), dataObj[["res"]][["seqtab.nochim"]]);
  colnames(seqtab.nochim)[1] <- "#NAME"
  colnames(seqtab.nochim)<-gsub("_F_filt.fastq.gz","", colnames(seqtab.nochim))
  write.table(seqtab.nochim, file = "microbiomeAnalyst_asv_seq.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  samples.out <- colnames(dataObj[["res"]][["seqtab.nochim"]])
  
  Samples <- sapply(strsplit(samples.out, "_F_filt.fastq.gz"), `[`, 1)
  colnames(dataObj[["res"]][["seqtab.nochim"]]) <- Samples
  samdf <- data.frame(Samples=Samples)
  meta_dt <- dataObj[["meta_info"]]
  grps <- sapply(Samples, function(x){
    unique(meta_dt[grep(x, meta_dt[,1]),2])[1]
  });
  samdf$Groups <- grps;
  rownames(samdf) <- Samples;
 
  taxa <- cbind(rownames(dataObj[["res"]][["taxa"]]),dataObj[["res"]][["taxa"]]) 
  colnames(taxa)[1] <- "#TAXONOMY"
  write.table(taxa, file = "microbiomeAnalyst_taxonomy_annotation.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  ## merge seqtab.nochim and taxa into OTU abundance table
  taxa_nms <- apply(dataObj[["res"]][["taxa"]], 1, FUN = function(x){
    if(is.na(x[6])){x[6] <- "uncultured"}
    if(is.na(x[7])){x[7] <- "uncultured bacterium"}
    paste(x, collapse = "; ")
  })
 
  otu_table <- cbind(NAME = taxa_nms, dataObj[["res"]][["seqtab.nochim"]])
  colnames(otu_table)[1] <- "#NAME"
  dataObj[["res"]][["otu_table"]] <- otu_table
  write.table(otu_table, file = "microbiomeAnalyst_16s_abund.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  ## save well-formatted meta data file
  meta_idx <- sapply(Samples, function(x){
    grep(x, meta_dt[,1])[1]
  });
  meta_dtx <- cbind(NAME = Samples, meta_dt[meta_idx, -1])
  colnames(meta_dtx) <- colnames(meta_dt)
  colnames(meta_dtx)[1] <- "#NAME"
  write.table(meta_dtx, file = "microbiomeAnalyst_16s_meta.txt", sep = "\t", row.names = FALSE, quote = FALSE)
  
  ps <- phyloseq(otu_table(dataObj[["res"]][["seqtab.nochim"]], taxa_are_rows=TRUE), 
                 sample_data(samdf), 
                 tax_table(dataObj[["res"]][["taxa"]]))
  mbSetObj$dataObj <- dataObj;
  mbSetObj$dataObj$ps_obj <- ps;
  # plot quick library size view
  libSizeQuickView(dataObj);
  
  # plot_richness
  p1 <- plot_richness(ps, measures=c("Shannon", "Simpson"), color="Groups")
  
  Cairo::Cairo(840, 500,file = paste0("richness.png"),dpi = 90,bg = "white")
  print(p1)
  dev.off()
  
  Cairo::CairoPDF(width = 8, height = 5,file = paste0("richness.pdf"))
  print(p1)
  dev.off()
  
  # plot_ordination
  ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
  ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
  
  p2 <- plot_ordination(ps.prop, ord.nmds.bray, color="Groups", title="Bray NMDS")
  Cairo::Cairo(840, 500,file = paste0("nmds.png"),dpi = 90,bg = "white")
  print(p2)
  dev.off()
  
  Cairo::CairoPDF(width = 8, height = 5,file = paste0("nmds.pdf"))
  print(p2)
  dev.off()
  
  # plot_bar
  top20 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:20]
  ps.top20 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
  ps.top20 <- prune_taxa(top20, ps.top20)
  p3 <- plot_bar(ps.top20, fill="Family") + facet_wrap(~Groups, scales="free_x")
  
  Cairo::Cairo(840, 500,file = paste0("bars.png"),dpi = 90,bg = "white")
  print(p3)
  dev.off()
  
  Cairo::CairoPDF(width = 8, height = 5,file = paste0("bars.pdf"))
  print(p3)
  dev.off()

  if (.on.public.web) {
    .set.mbSetObj(mbSetObj)
    return(1);
  } else {
    return(.set.mbSetObj(mbSetObj))
  }
}

libSizeQuickView <- function(dataObj = NULL, format = "png", dpi = 72){

  if(is.null(dataObj)){
    load("dataObj_completed.rda");
  } else {
    save(dataObj, file = "dataObj_completed.rda")
  }

  {
    otu_table <- dataObj[["res"]][["otu_table"]]
    rownames(otu_table) <- NULL
    rownms <- otu_table[,1]
    otu_table <- otu_table[,-1]
    otu_table <- apply(otu_table, 2, as.numeric)
    rownames(otu_table) <- rownms
    #data.proc <- qs::qread("data.proc");
    data_bef <- data.matrix(otu_table);
    smpl.sums <- colSums(data_bef);
    names(smpl.sums) <- colnames(data_bef);
    smpl.sums <- sort(smpl.sums);
  }
  
  smpl.sums <- rev(smpl.sums);
  vip.nms <- names(smpl.sums);
  names(smpl.sums) <- NULL;
  vip.nms <- substr(vip.nms, 1, 16);
  
  myH <- ncol(data_bef)*25 + 50;
  if(format == "png"){
      myW = 840
  } else {
      myW = 840
      myH = myH
  }
  
  # plot libsize_quickview
  Cairo::Cairo(file=paste0("libsize_quickview.", format), width=myW, height=myH, type=format, bg="white",dpi=dpi);
  xlim.ext <- GetExtendRange(smpl.sums, 10);
  par(mar=c(4,10,2,2));
  dotchart(smpl.sums, col="forestgreen", xlim=xlim.ext, pch=19, xlab="Read Counts", main="Library Raw Size Overview");
  mtext(side=2, at=1:length(vip.nms), vip.nms, las=2, line=1);
  text(x=smpl.sums,y=1:length(smpl.sums),labels= round(smpl.sums), col="blue", pos=4, xpd=T);
  dev.off();
  
}

plotSeqDiagnotics <- function(imageName = "diagnotics", format = "png", dpi = 120){
    if(file.exists("diagnotics_plot_src.qs")){
        p1 <- qs::qread("diagnotics_plot_src.qs")
        Cairo::Cairo(1600, 1250,file = paste0(imageName, ".",format), dpi = dpi, bg = "white", type=format)
        print(p1)
        dev.off()            
    }
    return(1)
}

### Testing code
#setwd("/home/qiang/Downloads/debug_dir/16s/MiSeq_SOP")
#PerformSeqCheck(home_dir = "/home/qiang/Downloads/debug_dir/16s/MiSeq_SOP")
#PerformSeqProcessing (dataObj, 
#                      truncLen=c(240,180), 
#                      maxN=0, 
#                      maxEE=c(2,2), 
#                      truncQ=2, 
#                      rm.phix=TRUE,
#                      compress=TRUE)
