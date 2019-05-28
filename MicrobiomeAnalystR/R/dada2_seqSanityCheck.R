##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#check number of paired files and reads quality
seqSanityCheck <- function(setParametersRes = setParametersRes, # results from setParameters
                           seq_folder = seq_folder, #dir of raw reads
                           file_pattern_F = file_pattern_F, # without compressed file extension, 1.fq, R1.fastq, L001_R1.fastq, must end with fq or fastq, do not use gz or bz2
                           file_pattern_R = file_pattern_R, # without compressed file extension, 2.fq, R2.fastq, L001_R2.fastq, must end with fq or fastq, do not use gz or bz2
                           plot_format = "pdf", #pdf, tiff, png, ...
                           ...){
  #tic("seqSanityCheck");
  #get R packages;
  require(dada2); require(R.utils);
  
  #get data from setParameters;
  file_compressed = setParametersRes$file_compressed;
  
  num_fF <- length(list.files(file.path(seq_folder),
                              pattern = file_pattern_F))
  print(paste0("detect ", num_fF," forward files"));
  
  num_fR <- length(list.files(file.path(seq_folder),
                              pattern = file_pattern_R));
  print(paste0("detect ", num_fR, " reverse files"));
  
  if(num_fF != num_fR){
    stop("Error: number of forward files are not equal to number of reverse files,
         please make sure all files are paired");
  };
  
  #decompress files
  if (file_compressed){
    fnFRCs <- sort(list.files(seq_folder, 
                              pattern = paste0(file_pattern_F, ".+", 
                                               "|",
                                               file_pattern_R, ".+"),
                              full.names = TRUE));
    print("decompress files");
    if (all(grepl("\\.gz$", fnFRCs))){
      lapply(fnFRCs,
             function(x){
               R.utils::gunzip(x, remove = FALSE)});
    } else if (all(grepl("\\.bz2$", fnFRCs))){
      lapply(fnFRCs,
             function(x){
               R.utils::bunzip2(x, remove = FALSE)})
    } else {
      stop("Please make sure the file format is either gz or bz2");
    };
  }
  
  #get files sorted
  print("sorting files");
  fnFs <- sort(list.files(seq_folder, 
                          pattern = paste0(file_pattern_F, "$"), 
                          full.names = TRUE));
  fnRs <- sort(list.files(seq_folder, 
                          pattern = paste0(file_pattern_R, "$"), 
                          full.names = TRUE));
  
  snF <- sapply(strsplit(basename(fnFs), file_pattern_F),
                `[`, 1);
  snR <- sapply(strsplit(basename(fnRs), file_pattern_R), 
                `[`, 1);
  ifelse(identical(snF, snR),
         sn <- snF,
         stop("Error: files are not paired"));
  print("all files are paired");
  
  print("plot reads quality graph");
  .plotQualityProfileLoop(fnFs, fnRs, sn = sn, plot_format, fd = seq_folder);
  
  #global variable;
  #seq_folder <<- seq_folder;
  #return objects
  seqSanityCheckResults <- list(fnFs = fnFs,
                                fnRs = fnRs,
                                sn = sn);
  #print(tac());
  return(seqSanityCheckResults);
}

