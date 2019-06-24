##################################################
## R script for MicrobiomeAnalyst
## Description: Data/resource management functions
## Author: Jeff Xia, jeff.xia@mcgill.ca
###################################################

#set parameters for dada2

setParameters <- function(file_compressed = TRUE, # logic TURE or FALSE, compressed format gz, bz2,
                          OS_is_windows = FALSE, # operation system is windows, linux or mac;
                          ... ){
  parameters <- list(file_compressed = file_compressed,
                     OS_is_windows = OS_is_windows);
  return(parameters);
}

#a loop function for plotting reads quanlity graph.
.plotQualityProfileLoop <- function(f_F = f_F, #forward reads file
                                    f_R = f_R,#reverse reads file
                                    sn = sn, # sample name
                                    plot_format = plot_format, # pdf, tiff, ....
                                    fd = fd ){ # destion file
  fnFRs <- vector();
  for (i in seq(length(sn))){
    print(i);
    fnFRs[2*i -1] <- f_F[i];
    fnFRs[2*i] <- f_R[i];
  };
  n_plot_floor <- floor(length(sn) / 12);
  n_plot_remainder <- length(sn) %% 12;
  if (length(sn) <= 12){
    plotQualityProfile(fnFRs) -> p1; print(p1);
    ggsave(file.path(fd,
                     paste0("sequence_quality.", plot_format)),
           width = 20, height = 10);
  } else {
    for (i in seq(n_plot_floor)){
      plotQualityProfile(fnFRs[(i*12 - 11):(i*12)]) -> p2; print(p2);
      ggsave(file.path(fd,
                       paste0("sequence_quality_", i, ".", plot_format)),
             width = 20, height = 10)
    };
    if (n_plot_remainder > 0){
      plotQualityProfile(fnFRs[(n_plot_floor*12 + 1) : length(sn)]) -> p3; print(p3);
      ggsave(file.path(fd,
                       paste0("sequence_quality_", n_plot_floor+1, ".", plot_format)),
             width = 20, height = 10);
    };
  };
}

#check number of paired files and reads quality
seqSanityCheck <- function(setParametersRes = setParametersRes, # results from setParameters
                           seq_folder = seq_folder, #dir of raw reads
                           file_pattern_F = file_pattern_F, # without compressed file extension, .1.fq, R1.fastq, _L001_R1.fastq, must end with fq or fastq, do not use gz or bz2
                           file_pattern_R = file_pattern_R, # without compressed file extension, .2.fq, R2.fastq, _L001_R2.fastq, must end with fq or fastq, do not use gz or bz2
                           plot_format = "pdf", #pdf, tiff, png, ...
                           ...){
  #tic("seqSanityCheck");
  #get R packages;
  require(dada2); require(R.utils); #for everything

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
  write.table(cbind.data.frame("#SAMPLE" = row.names(seqtab),
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

  print("write chimera frequence");
  chimaera_table <- as.data.frame((1 - rowSums(seqtab.nochim)/rowSums(seqtab))*100);
  names(chimaera_table) <- "Frequence (%)";
  chimaera_table$Sample <- row.names(chimaera_table);
  write.table(chimaera_table,
              file = file.path("sequence_table", "chimera_frequence.txt"),
              row.names = FALSE,
              sep = "\t",
              quote = FALSE);

  print("write sequence table without chimera");
  write.table(cbind.data.frame("#NAME" = row.names(seqtab.nochim),
                               seqtab.nochim),
              file = file.path("sequence_table", "sequence_abundance_table_without_chimera.txt"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE);


  #transpose otu table for MicrobiomAnalyst submission
  seqtab.nochim_t <- as.data.frame(t(seqtab.nochim));
  write.table(cbind.data.frame("#NAME" = row.names(seqtab.nochim_t),
                               seqtab.nochim_t),
              file = file.path("sequence_table", "sequence_abundance_table_without_chimera_for_submission_to_MicrobiomeAnalyst.txt"),
              sep = "\t",
              row.names = FALSE,
              quote = FALSE);


  #extact sequences and convert to fasta for submission to Picrust2
  seq_table <- data.frame("ID" = paste0(">", colnames(seqtab.nochim)),
                          "seq" = colnames(seqtab.nochim));
  seq_table[] <- lapply(seq_table, as.character)

  seq_list <- list();
  for(i in 1:nrow(seq_table)){
    j <- 2*i;
    seq_list[[(j - 1)]] <- seq_table[i, 1];
    seq_list[[j]] <- seq_table[i, 2];
  }

  seq_list <- do.call(rbind.data.frame, seq_list)
  names(seq_list) <- NULL;
  write.table(seq_list,
              file = file.path("sequence_table",
                               "sequence_without_chimera_for_submission_to_Picrust2.fasta"),
              quote = FALSE,
              row.names = FALSE);

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

#download tax assignment database and assign tax infor
assignTax <- function(constructSeqTabRes = constructSeqTabRes, #results from constructSeqTab
                      ref_db = "silva",
                      ...){
  #tic("assignTax_time");
  #get data from constructSeqTabRes
  OS_is_windows = setParametersRes$OS_is_windows;
  seqtab.nochim = constructSeqTabRes$seqtab.nochim;

  #assign taxonomy########
  #########################
  print("check and create tax file");
  ifelse(!file.exists("tax"),
         dir.create("tax"),
         FALSE);

  print("downlod database");
  if (ref_db == "rdp"){
    # ref_db <- "rdp_train_set_16.fa.gz";
    ifelse(!file.exists("tax/rdp_train_set_16.fa.gz"),
           download.file(url = "https://zenodo.org/record/801828/files/rdp_train_set_16.fa.gz?download=1",
                         destfile = file.path("tax", "rdp_train_set_16.fa.gz"),
                         method = "auto"),
           FALSE);
    ifelse(!file.exists("tax/rdp_species_assignment_16.fa.gz"),
           download.file(url = "https://zenodo.org/record/801828/files/rdp_species_assignment_16.fa.gz?download=1",
                         destfile = file.path("tax", "rdp_species_assignment_16.fa.gz"),
                         method = "auto"),
           FALSE);
    print("assign taxa infor against RDP");
    taxa <- assignTaxonomy(seqtab.nochim,
                           file.path("tax", "rdp_train_set_16.fa.gz"),
                           multithread = !OS_is_windows);
    taxa <- addSpecies(taxa,
                       file.path("tax", "rdp_species_assignment_16.fa.gz"));

  } else if (ref_db == "silva"){
    # ref_db_tax <- "silva_nr_v132_train_set.fa.gz";
    # ref_db_spe <- "silva_species_assignment_v132.fa.gz";
    ifelse(!file.exists("tax/silva_nr_v132_train_set.fa.gz"),
           download.file(url = "https://zenodo.org/record/1172783/files/silva_nr_v132_train_set.fa.gz?download=1",
                         destfile = file.path("tax", "silva_nr_v132_train_set.fa.gz"),
                         method = "auto"),
           FALSE);
    ifelse(!file.exists("tax/silva_species_assignment_v132.fa.gz"),
           download.file(url = 'https://zenodo.org/record/1172783/files/silva_species_assignment_v132.fa.gz?download=1',
                         destfile = file.path("tax", "silva_species_assignment_v132.fa.gz"),
                         method = 'auto'),
           FALSE);
    print("assign taxa infor against silva");
    taxa <- assignTaxonomy(seqtab.nochim,
                           file.path("tax", "silva_nr_v132_train_set.fa.gz"),
                           multithread = !OS_is_windows);

    taxa <- addSpecies(taxa,
                       file.path("tax", "silva_species_assignment_v132.fa.gz"));
  } else if (ref_db == "greengenes") {
    #ref_db <- "gg_13_8_train_set_97.fa.gz";
    ifelse(!dir.exists("tax/gg_13_8_train_set_97.fa.gz"),
           download.file(url = "https://zenodo.org/record/158955/files/gg_13_8_train_set_97.fa.gz?download=1",
                         destfile = file.path("tax", "gg_13_8_train_set_97.fa.gz"),
                         method = "auto"),
           FALSE);
    print("assign taxa infor against greengenes");
    taxa <- assignTaxonomy(seqtab.nochim,
                           file.path("tax", "gg_13_8_train_set_97.fa.gz"),
                           multithread = !OS_is_windows);
  } else {
    stop("please specify a database!")
  };
  # taxa2 <- as.data.frame(taxa);
  # taxa2 <- within(taxa2,
  #                 "NAME" <-  paste(Kingdom,
  #                                  Phylum,
  #                                  Class,
  #                                  Order,
  #                                  Family,
  #                                  Genus,
  #                                  Species,
  #                                  sep = "; "));
  print("write taxa table");
  # write.table(cbind.data.frame("#NAME" = taxa2$NAME,
  #                              t(seqtab.nochim)),
  #             file = file.path("tax",
  #                              paste0("taxa_table_against_", ref_db,
  #                                     "_submit_to_MicrobiomeAnalyst.txt")),
  #             row.names = FALSE,
  #             quote = FALSE,
  #             sep = "\t");
  write.table(cbind.data.frame("#TAXONOMY" = row.names(taxa),
                         taxa),
              file = file.path("tax",
                               paste0("taxa_table_against_", ref_db, ".txt")),
              row.names = FALSE,
              quote = FALSE,
              sep = "\t");
  write.table(taxa,
              file = file.path("tax",
                               paste0("taxa_table_against_", ref_db, "_without_id.txt")),
              row.names = FALSE,
              quote = FALSE,
              sep = "\t");
  #inspect the taxonomic assignments;
  # taxa.print <- taxa # Removing sequence rownames for display only
  # rownames(taxa.print) <- NUSLL;
  #assignTaxRes <- list(taxa); return(assignTaxRes);
  #print(toc())
}

#contruction of plylogenetic tree, could be very slow
constructPhyloTree <- function(constrcutSeqTabRes = constrcutSeqTabRes, #results from constrcutSeqTab
                               ...){
  #tic("constructPhyloTree")
  require(DECIPHER);
  require(phangorn);

  #get data;
  seqtab.nochim = constructSeqTabRes$seqtab.nochim;

  print("get sequences");
  seqs <- getSequences(seqtab.nochim);
  names(seqs) <- seqs;
  print("do alignment");
  alignment <- AlignSeqs(DNAStringSet(seqs),
                         anchor = NA,
                         verbose = TRUE);
  print("construct tree");
  phangAlign <- phyDat(as(alignment, "matrix"), type = "DNA")
  dm <- dist.ml(phangAlign)
  treeNJ <- NJ(dm) # Note, tip order != sequence order
  fit = pml(treeNJ, data = phangAlign)
  fitGTR <- update(fit, k = 4, inv = 0.2)
  fitGTR <- optim.pml(fitGTR, model = "GTR",
                      optInv = TRUE,
                      optGamma = TRUE,
                      rearrangement = "stochastic",
                      control = pml.control(trace = 0));
  ifelse(!dir.exists("phylogenetic_tree"),
         dir.create("phylogenetic_tree"),
         FALSE);
  write.tree(fitGTR$tree,
             file = file.path("phylogenetic_tree", "phylotree.tre"));
  #detach("package:phangorn", unload=TRUE);
  #print(toc());
}
