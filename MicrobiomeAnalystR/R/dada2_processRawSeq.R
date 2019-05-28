##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#conduct reads cleaning, error correction, reads merging
processRawSeq <- function(setParametersRes = setParametersRes, # results from setParameters
                          seqSanityCheckRes = seqSanityCheckRes, # results from seqSanityCheck
                          reads_trim_length_F_R = c(240, 160), #retained reads length for forward and reverse, 
                          #depending on the quailty of reads, the quailty graph can be obtained by seqSanityCheck results;
                          plot_format = "pdf", #pdf, tiff, png, ...
                          ... ){
  #tic("processRawSeq");
  #get data from setParameters and seqSanityCheck;
  file_compressed = setParametersRes$file_compressed;
  OS_is_windows = setParametersRes$OS_is_windows;
  sn = seqSanityCheckRes$sn;
  fnFs = seqSanityCheckRes$fnFs;
  fnRs = seqSanityCheckRes$fnRs;
  
  #filter and trim;
  print("creat filtered dir");
  ifelse(!dir.exists("filtered"),
         dir.create("filtered"),
         FALSE);#create dir;
  print("assign filtered files");
  filtFs <- file.path("filtered", paste0(sn, "_F_filt.fastq.gz"));
  filtRs <- file.path("filtered", paste0(sn, "_R_filt.fastq.gz"));
  
  print("name the filtered files");
  names(filtFs) <- sn;
  names(filtRs) <- sn;
  
  print("filter files");
  out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, 
                       truncLen = reads_trim_length_F_R,
                       maxN = 0, 
                       maxEE = c(2,2), 
                       truncQ = 2, 
                       rm.phix = TRUE,
                       compress = TRUE, 
                       multithread = !OS_is_windows); #on windows set multithread = FALSE;
  print("write filtered files matrix");
  write.table(out, 
              file = file.path("filtered", "filtered_data.txt"),
              sep = "\t",
              row.names = FALSE);
  print("plot quality score of clean reads");
  .plotQualityProfileLoop(filtFs, filtRs, sn, plot_format, fd = "filtered");
  
  #remove uncompressed files
  ##########################
  if(file_compressed){
    print("remove uncompressed files");
    file.remove(fnFs);
    file.remove(fnRs);
  }
  
  #learn the error rates
  #######################
  print("learn errors");
  errF <- learnErrors(filtFs, multithread = !OS_is_windows);
  errR <- learnErrors(filtRs, multithread = !OS_is_windows);
  
  print("plot error graph");
  plotErrors(errF, nominalQ = TRUE);
  ggsave(file.path("filtered",
                   paste0("error_plot_forward_reads.", plot_format)));
  plotErrors(errR, nominalQ = TRUE);
  ggsave(file.path("filtered",
                   paste0("error_plot_reverse_reads.", plot_format)));
  
  #sample inference#########
  ##########################
  print("error correction");
  dadaFs <- dada(filtFs, err = errF, multithread = !OS_is_windows);
  dadaRs <- dada(filtRs, err = errF, multithread = !OS_is_windows);
  ####consider the pooling samples will increase sensitivity to sequence variants 
  #with low requences in multiple samples;
  #dada(..., pool=TRUE);
  
  #merge paired reads##########
  ##############################
  print("merge reads")
  mergers <- mergePairs(dadaFs, 
                        filtFs, 
                        dadaRs, 
                        filtRs,
                        verbose = TRUE); #unoverlaped seq will be discarded
  
  #write mergered sequences
  print("inspect and create mergered_sesquence file")
  ifelse(!dir.exists("mergered_sequence"),
         dir.create("mergered_sequence"),
         FALSE);
  print("write mergered sequences");
  for(i in names(mergers)){
    write.table(mergers[[i]],
                file = file.path("mergered_sequence", paste0(i, "_mergered_sequences.txt")),
                row.names = FALSE,
                sep = "\t",
                quote = FALSE);
  }
  
  
  
  processRawSeqRes <- list(mergers = mergers,
                           out = out,
                           dadaFs = dadaFs,
                           dadaFs = dadaFs,
                           sn = sn);
  #print(tac());
  return(processRawSeqRes);
}

