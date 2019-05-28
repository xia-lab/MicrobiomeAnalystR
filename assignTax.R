assignTax <- function(constructSeqTabRes = constructSeqTabRes,
                       ref_db = "silva",
                       ...){
  tic("assignTax_time");
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
  
  print("write taxa table");
  write.table(data.frame("ID" = row.names(taxa),
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
  # rownames(taxa.print) <- NULL;
  #assignTaxRes <- list(taxa); return(assignTaxRes);
  print(toc())
}
