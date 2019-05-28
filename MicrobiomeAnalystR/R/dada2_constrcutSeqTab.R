##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#generation of sequence table and chimera removal
constructSeqTab <- function(setParametersRes = setParametersRes, # results from setParameters
                            processRawSeqRes = processRawSeqRes,#results from processRawSeqRes
                            ...){
  #tic("constructSeqTab");
  #get data from processRawSeq
  OS_is_windows = setParametersRes$OS_is_windows;
  mergers = processRawSeqRes$mergers;
  out = processRawSeqRes$out;
  dadaFs = processRawSeqRes$dadaFs;
  dadaRs = processRawSeqRes$dadaFs;
  sn = processRawSeqRes$sn;
  #construct sequence table###########
  ###################################
  print("construct sequence table");
  seqtab <- makeSequenceTable(mergers);
  print("inspect and create sequence table file");
  ifelse(!dir.exists("seqence_table"),
         dir.create("sequence_table"),
         FALSE);
  print("write sequence abundance table");
  write.table(cbind.data.frame("Sample" = row.names(seqtab),
                               seqtab),
              file = file.path("sequence_table", "sequence_abundance_table.txt"),
              sep = "\t",
              quote = FALSE);
  print("write sequence length frequence table");
  seqFreq <- as.data.frame(table(nchar(getSequences(seqtab))));
  names(seqFreq) <- c("Length", "Frequence");
  write.table(seqFreq,
              file = file.path("sequence_table", "sequence_length_freqence.txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE);
  
  ##remove chimeras
  #####################
  print("remove chimeras");
  seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", 
                                      multithread = !OS_is_windows, 
                                      verbose = TRUE);
  print("write sequence table without chimera");
  write.table(cbind.data.frame("Sample" = row.names(seqtab.nochim),
                               seqtab.nochim),
              file = file.path("sequence_table", "sequence_abundance_table_without_chimera.txt"),
              sep = "\t",
              quote = FALSE);
  
  print("write chimera frequence");
  chimaera_table <- as.data.frame((1 - rowSums(seqtab.nochim)/rowSums(seqtab))*100);
  names(chimaera_table) <- "Frequence (%)";
  chimaera_table$Sample <- row.names(chimaera_table);
  write.table(chimaera_table,
              file = file.path("sequence_table", "chimera_frequence.txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE);
  
  #track reads through the pipeline;
  ##################################
  print("track reads");
  getN <- function(x) sum(getUniques(x))
  track <- cbind(out, 
                 sapply(dadaFs, getN), 
                 sapply(dadaRs, getN), 
                 sapply(mergers, getN), 
                 rowSums(seqtab.nochim));
  
  # If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
  colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
  rownames(track) <- sn;
  print("saving track reads");
  write.table(cbind.data.frame("Sample" = row.names(track),
                               track),
              file = file.path("sequence_table", "track_reads.txt"),
              sep = "\t",
              quote = FALSE);
  constrcutSeqTabRes <- list(seqtab.nochim = seqtab.nochim);
  #print(tac());
  return(constrcutSeqTabRes);
}
