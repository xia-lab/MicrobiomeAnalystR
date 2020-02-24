fastpR <- function(fastp_seq_clean,
                   raw_seq_dir,
                   n_threads = 8,
                   ...){
    pacman::p_load(crayon, parallel)
    if (n_threads > detectCores()){
        n_threads <- detectCores();
    };
    cat(green("starting reads cleaning"), "\n");
    cat(blue(format(Sys.time(), usetz = TRUE)), "\n");
    #create dir for trimmomatic
    ifelse(!dir.exists(fastp_seq_clean),
           dir.create(fastp_seq_clean),
           FALSE);
    smpname1 <- gsub("_R1.fastq.gz", "", 
                     list.files(raw_seq_dir, pattern = "_R1\\.fastq\\.gz"));
    cat(blue("detected", length(smpname1), "forward reads files"), "\n");
    smpname2 <- gsub("_R2.fastq.gz", "", 
                     list.files(raw_seq_dir, pattern = "_R2\\.fastq\\.gz"));
    cat(blue("detected", length(smpname2), "reverse reads files"), "\n");
    
    if(!identical(smpname1, smpname2)){
        stop(red("forward and reverse reads files are not matched,
                 please make sure they are same"))
    };
    smpname <- smpname1;
    cat(blue("detected", length(smpname), "forward and reverse reads files"), "\n");
    rf1 <- file.path(raw_seq_dir, paste0(smpname, "_R1.fastq.gz"));
    rf2 <- file.path(raw_seq_dir, paste0(smpname, "_R2.fastq.gz"));
    rfp1 <- file.path(fastp_seq_clean, paste0(smpname, "_R1.trim.paired.fastq.gz"));
    rfp2 <- file.path(fastp_seq_clean, paste0(smpname, "_R2.trim.paired.fastq.gz"));
    #rfup1 <- file.path(fastp_seq_clean, paste0(smpname, "_R1.trim.unpaired.fastq.gz"));
    #rfup2 <- file.path(fastp_seq_clean, paste0(smpname, "_R2.trim.unpaired.fastq.gz"));
    
    #run fastp
    lapply(1:length(smpname),
           FUN = function(x){
               cat(yellow("processing sample", smpname[x], 
                          "(", x, "/", length(smpname), "|",
                          round(x * 100 / length(smpname), 4),  "% )",
                          format(Sys.time(), usetz = TRUE)), "\n");
               system(paste("fastp",
                            "-w", n_threads,
                            "-c",
                            "-i", rf1[x], 
                            "-I", rf2[x],
                            "-o", rfp1[x], 
                            "-O", rfp2[x],
                            sep = " "));
           });
    cat(green("reads cleaning ended"), "\n");
    cat(blue(format(Sys.time(), usetz = TRUE)), "\n");
    }

rmHostCon <- function(host_free,
                      fastp_seq_clean,
                      host_genome.fa,
                      host_genome,
                      n_threads = 8,
                      ...){
    pacman::p_load(crayon, parallel)
    if (n_threads > detectCores()){
        n_threads <- detectCores();
    };
    cat(green("starting reads cleaning"), "\n");
    cat(blue(format(Sys.time(), usetz = TRUE)), "\n");
    #create dir for host free
    ifelse(!dir.exists(host_free),
           dir.create(host_free),
           FALSE);
    
    smpname1 <- gsub("_R1.trim.paired.fastq.gz", "",
                     list.files(fastp_seq_clean, 
                                pattern = "_R1\\.trim\\.paired\\.fastq\\.gz"));
    cat(blue("detected", length(smpname1), "forward reads files"), "\n");
    smpname2 <- gsub("_R2.trim.paired.fastq.gz", "",
                     list.files(fastp_seq_clean, 
                                pattern = "_R2\\.trim\\.paired\\.fastq\\.gz"));
    cat(blue("detected", length(smpname2), "reverse reads files"), "\n");
    
    if(!identical(smpname1, smpname2)){
        stop(red("forward and reverse reads files are not matched,
                     please make sure they are same"))
    };
    smpname <- smpname1;
    rfpi1 <- file.path(fastp_seq_clean, paste0(smpname, "_R1.trim.paired.fastq.gz"));
    rfpi2 <- file.path(fastp_seq_clean, paste0(smpname, "_R2.trim.paired.fastq.gz"));
    rfpo1 <- file.path(host_free, paste0(smpname, "_R1.trim.paired.clean.fastq.gz"));
    rfpo2 <- file.path(host_free, paste0(smpname, "_R2.trim.paired.clean.fastq.gz"));
    
    cat(blue("starting removing host contaminant using bbmap", format(Sys.time(), usetz = TRUE)), "\n");

    lapply(1:length(smpname),
           FUN = function(x){
               cat(yellow("starting reads mapping for sample", smpname[x],
                          "(", x, "/", length(smpname), "|",
                          round(x * 100 / length(smpname), 4), "% )",
                          format(Sys.time(), usetz = TRUE)), "\n");
               system(paste0("bbmap.sh quickmatch fast",
                             " ref=", host_genome.fa,
                             " path=", file.path(host_free, host_genome),
                             " in1=", rfpi1[x],
                             " in2=", rfpi2[x],
                             " outu1=", rfpo1[x],
                             " outu2=", rfpo2[x],
                             " threads=", n_threads));
           });
}

kaijuMap <- function(kaiju_map,
                     host_free, #output from bbmap
                     nodes.dmp,
                     kaiju_db.fmi,
                     run_mode = "greedy", # either mem or greedy;
                     n_threads = 8){
    pacman::p_load(crayon, parallel);
    #create dir for kaiju_map
    ifelse(!dir.exists(kaiju_map),
           dir.create(kaiju_map),
           FALSE);
    #check proteins_db exists or not
    if(dir.exists(host_free)){
            smpname1 <- gsub("_R1.trim.paired.clean.fastq.gz", "",
                             list.files(host_free, 
                                        pattern = "_R1\\.trim\\.paired\\.clean\\.fastq\\.gz"));
            cat(blue("detected", length(smpname1), "forward reads files"), "\n");
            smpname2 <- gsub("_R2.trim.paired.clean.fastq.gz", "",
                             list.files(host_free, 
                                   pattern = "_R2\\.trim\\.paired\\.clean\\.fastq\\.gz"));
            cat(blue("detected", length(smpname2), "reverse reads files"), "\n");
            
            if(!identical(smpname1, smpname2)){
                stop(red("forward and reverse reads files are not matched,
                         please make sure they are same"))
            };
            smpname <- smpname1;
            rfp1 <- file.path(host_free, paste0(smpname, "_R1.trim.paired.clean.fastq.gz"));
            rfp2 <- file.path(host_free, paste0(smpname, "_R2.trim.paired.clean.fastq.gz"));
            } else {
                stop(red("please check host_free dir!!"))
            }
  rfo <- file.path(kaiju_map, paste0(smpname, ".txt"));
    
    if (!file.exists(nodes.dmp) | !file.exists(kaiju_db.fmi)) {
        stop(red("Please provide nodes.dmp or kaiju_db.fmi files"));
    };
    
    cat(green("detected", length(smpname), "samples for Kaiju mapping"), "\n");
    cat(blue("starting reads Kaiju classification", 
             format(Sys.time(), usetz = TRUE)), "\n");
    lapply(1:length(smpname),
           FUN = function(x){
               cat(yellow("starting reads mapping for sample", smpname[x],
                          "(", x, "/", length(smpname), "|",
                          round(x * 100 / length(smpname), 4), "% )",
                          format(Sys.time(), usetz = TRUE)), "\n");
               system(paste("kaiju",
                            "-a", run_mode,
                            "-z", n_threads,
                            "-t", nodes.dmp,
                            "-f", kaiju_db.fmi,
                            "-i", rfp1[x],
                            "-j", rfp2[x],
                            "-o", rfo[x],
                            sep = " "));
             })
    cat(blue("reads mapping ended", 
             format(Sys.time(), usetz = TRUE)), "\n");
}

kaijuTaxa <- function(kaiju_taxa,
                      kaiju_map,
                      nodes.dmp,
                      names.dmp,
                      ...){
  pacman::p_load(crayon);
  #create dir for kaiju_map
  ifelse(!dir.exists(kaiju_taxa),
         dir.create(kaiju_taxa),
         FALSE);
  #check proteins_db exists or not
  if(dir.exists(kaiju_map)){
    smpname <- gsub("\\.txt", "", list.files(kaiju_map, pattern = "\\.txt"));
    cat(blue("detected", length(smpname), "files"), "\n");
  } else {
    stop(red("please check kaiju_map dir!!"))
  }
  
  rfi <- file.path(kaiju_map, paste0(smpname, ".txt"));
  rfo <- file.path(kaiju_taxa, paste0(smpname, "_summary.tsv"));
  
  if (!file.exists(nodes.dmp)) {
    stop(red("Please provide nodes.dmp file!"));
  };
  
  cat(green("detected", length(smpname), "samples for Kaiju taxa"), "\n");
  cat(blue("starting extracting Kaiju taxa", 
           format(Sys.time(), usetz = TRUE)), "\n");
  lapply(1:length(smpname),
         FUN = function(x){
           cat(yellow("starting sample", smpname[x],
                      "(", x, "/", length(smpname), "|",
                      round(x * 100 / length(smpname), 4), "% )",
                      format(Sys.time(), usetz = TRUE)), "\n");
           system(paste("kaiju2table",
                        "-t", nodes.dmp,
                        "-n", names.dmp,
                        "-r", "species",
                        "-l", "superkingdom,phylum,class,order,family,genus,species",
                        "-o", rfo[x],
                        rfi[x],
                        sep = " "));
         })
  cat(blue("reads mapping ended", 
           format(Sys.time(), usetz = TRUE)), "\n");
}

getOTUTable <- function(otu_dir,
                        kaiju_taxa,
                        n_threads = 8,
                        ...){
  pacman::p_load(crayon, data.table, tidyverse, dplyr, magrittr, plyr);
  #create dir for kaiju_map
  ifelse(!dir.exists(otu_dir),
         dir.create(otu_dir),
         FALSE);
  #check kaiju_taxa exists or not
  if(dir.exists(kaiju_taxa)){
    smpname <- gsub("_summary\\.tsv", "", 
                    list.files(kaiju_taxa, pattern = "_summary\\.tsv"));
    cat(blue("detected", length(smpname), "kaiju_taxa files"), "\n");
  } else {
    stop(red("please check kaiju_taxa dir!!"));
  }
  rfi <- file.path(kaiju_taxa, paste0(smpname, "_summary.tsv"));
  
  otu_table <- lapply(1:length(smpname),
                      FUN = function(x){
                        fread(file = rfi[x],
                              select = c(3, 4, 5),
                              header = TRUE,
                              quote = "",
                              nThread = n_threads) %>% 
                          mutate(taxon_id = as.character(taxon_id)) %>% 
                          select(taxon_id, taxon_name, everything()) %>% 
                          dplyr::rename(!!smpname[x] := "reads")
                      })
  
  
  join_all(otu_table,
           type = "full") -> otu_table;
  otu_table %<>% 
    replace(is.na(.), 0);
  write.table(otu_table,
              file.path(otu_dir, "kaiju_otu_table.txt"),
              row.names = FALSE, quote = F, sep = "\t")
}