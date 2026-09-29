################################################
### Raw shotgun metagenomics data processing
### Species-resolved taxonomic and functional profiles with Meteor2
### (Ghozlane et al., Microbiome 2025; https://github.com/metagenopolis/meteor)
### Runs on the user's computer; the web server takes the resulting tables.
###
### Plan
### 1. CheckMeteor / ListMeteorCatalogues / DownloadMeteorCatalogue: tools and
###    reference catalogues (done).
### 2. RunMeteorProfiling: fastq -> mapping -> profile -> merge, restartable (done);
###    strain = TRUE keeps the filtered alignments that step 4 needs.
### 3. BuildCarrierTables: which species carry each resistance gene and CAZyme family
###    in each sample, for per-species host attribution (skeleton).
### 4. RunMeteorStrain: strain calls and pairwise strain distances (skeleton).
### 5. ReadMeteorResults / PackMeteorResults: read the merged tables, zip them for
###    upload (done). The zip picks up the carrier tables of step 3 automatically.
### Tested with Meteor2 2.0.22 and bowtie2 2.5.5 (dog gut catalogue, 3 samples).
################################################

#'Locate Meteor2 and bowtie2
#'@description Finds the meteor executable (argument, option
#'"MicrobiomeAnalystR.meteor", or PATH) and bowtie2, and reports the Meteor2 version.
#'Install with conda (conda create -n meteor -c conda-forge -c bioconda meteor)
#'or pip (pip install meteor) plus bowtie2.
#'@param meteor.bin Path to the meteor executable; NULL to search.
#'@return A list with the meteor path, the bowtie2 path and the Meteor2 version.
#'@export
CheckMeteor <- function(meteor.bin = NULL){
  if(is.null(meteor.bin)){
    meteor.bin <- getOption("MicrobiomeAnalystR.meteor", Sys.which("meteor"));
  }
  if(!nzchar(meteor.bin) || !file.exists(meteor.bin)){
    stop("meteor was not found. Install Meteor2 (conda: 'conda create -n meteor -c conda-forge -c bioconda meteor'; ",
         "pip: 'pip install meteor') and give its path as meteor.bin or options(MicrobiomeAnalystR.meteor = ...).");
  }
  # bowtie2 is looked up next to meteor first (conda or venv bin folder), then on PATH
  bt2 <- file.path(dirname(meteor.bin), "bowtie2");
  if(!file.exists(bt2)){
    bt2 <- Sys.which("bowtie2");
  }
  if(!nzchar(bt2)){
    stop("bowtie2 was not found on the PATH; Meteor2 needs it for read mapping.");
  }
  ver <- suppressWarnings(system2(meteor.bin, "--version", stdout = TRUE, stderr = TRUE));
  ver <- sub(".*version\\s*", "", ver[length(ver)]);
  return(list(meteor = unname(meteor.bin), bowtie2 = unname(bt2), version = ver));
}

#'List the Meteor2 reference catalogues
#'@description Returns the catalogue names the installed Meteor2 can download
#'(human gut, oral and skin; mouse, rat, pig, chicken, dog, cat and rabbit gut).
#'@param meteor.bin Path to the meteor executable; NULL to search.
#'@return A character vector of catalogue names.
#'@export
ListMeteorCatalogues <- function(meteor.bin = NULL){
  tools <- CheckMeteor(meteor.bin);
  hlp <- suppressWarnings(system2(tools$meteor, c("download", "-h"), stdout = TRUE, stderr = TRUE));
  hlp <- paste(hlp, collapse = " ");
  nms <- regmatches(hlp, regexpr("\\{[^}]+\\}", hlp));
  if(length(nms) == 0){
    stop("Could not read the catalogue list from 'meteor download -h'.");
  }
  return(strsplit(gsub("[{} ]", "", nms), ",")[[1]]);
}

#'Download a Meteor2 reference catalogue
#'@description Downloads and unpacks a gene catalogue (md5 checked). The full
#'catalogue is needed for functional profiles (KEGG modules, CAZymes, antibiotic
#'resistance); the light catalogue (fast = TRUE) gives species profiles only.
#'@param name Catalogue name, see ListMeteorCatalogues().
#'@param ref.dir Folder that receives the catalogue.
#'@param fast TRUE for the light (species only) catalogue.
#'@param meteor.bin Path to the meteor executable; NULL to search.
#'@return The path of the catalogue folder.
#'@export
DownloadMeteorCatalogue <- function(name, ref.dir, fast = FALSE, meteor.bin = NULL){
  tools <- CheckMeteor(meteor.bin);
  dir.create(ref.dir, showWarnings = FALSE, recursive = TRUE);
  args <- c("download", "-i", name, "-c", "-o", shQuote(ref.dir));
  if(fast){
    args <- c(args, "--fast");
  }
  .run.meteor(tools, args, file.path(ref.dir, paste0("download_", name, ".log")));
  cat.dir <- file.path(ref.dir, if(fast) paste0(name, "_taxo") else name);
  if(!dir.exists(cat.dir)){
    stop("The catalogue folder ", cat.dir, " was not created.");
  }
  return(normalizePath(cat.dir));
}

#'Profile raw shotgun reads with Meteor2
#'@description Imports the FASTQ files, maps each sample to the catalogue, computes
#'species (MSP) and functional profiles and merges all samples into one set of tables.
#'Reads should be quality-filtered and host reads removed beforehand. Samples already
#'processed in out.dir are skipped, so an interrupted run can be restarted.
#'@param fastq.dir Folder with the FASTQ files (.fastq/.fq, optionally .gz/.bz2/.xz).
#'Paired files end with _R1/_R2, .R1/.R2, _1/_2 or .1/.2 before the extension.
#'@param catalogue.dir Catalogue folder returned by DownloadMeteorCatalogue().
#'@param out.dir Output folder; the merged tables are written to out.dir/merged.
#'@param paired TRUE for paired-end files.
#'@param threads Number of threads for bowtie2.
#'@param normalization "coverage" (default), "fpkm" or "raw".
#'@param completeness Fraction of a module's KOs a species must carry for the module
#'to be counted as present in that species (Meteor2 default 0.9).
#'@param sample.mask Optional regular expression that extracts the sample name from a
#'file name, to group several runs of one library.
#'@param min.occurrence Keep species detected in at least this many samples.
#'@param strain TRUE to keep the filtered alignments needed by RunMeteorStrain().
#'@param meteor.bin Path to the meteor executable; NULL to search.
#'@return The path of the folder with the merged tables (invisibly).
#'@export
RunMeteorProfiling <- function(fastq.dir, catalogue.dir, out.dir, paired = TRUE, threads = 4,
                               normalization = c("coverage", "fpkm", "raw"), completeness = 0.9,
                               sample.mask = NULL, min.occurrence = 1, strain = FALSE, meteor.bin = NULL){
  normalization <- match.arg(normalization);
  tools <- CheckMeteor(meteor.bin);
  if(length(list.files(catalogue.dir, pattern = "_reference\\.json$")) == 0){
    stop("No *_reference.json in ", catalogue.dir, "; give the catalogue folder returned by DownloadMeteorCatalogue().");
  }
  dir.create(out.dir, showWarnings = FALSE, recursive = TRUE);
  out.dir <- normalizePath(out.dir);
  log.file <- file.path(out.dir, "meteor.log");
  fq.dir <- file.path(out.dir, "fastq");
  map.dir <- file.path(out.dir, "mapping");
  prof.dir <- file.path(out.dir, "profiles");
  merge.dir <- file.path(out.dir, "merged");

  # meteor fastq does not import into an existing folder, so a restarted run reuses it
  if(dir.exists(fq.dir)){
    message("Meteor2 ", tools$version, ": FASTQ files already imported (delete ", fq.dir, " to import again)");
  }else{
    message("Meteor2 ", tools$version, ": importing FASTQ files");
    args <- c("fastq", "-i", shQuote(normalizePath(fastq.dir)), "-o", shQuote(fq.dir));
    if(paired){
      args <- c(args, "-p");
    }
    if(!is.null(sample.mask)){
      args <- c(args, "-m", shQuote(sample.mask));
    }
    .run.meteor(tools, args, log.file);
  }
  samples <- sort(basename(list.dirs(fq.dir, recursive = FALSE)));
  if(length(samples) == 0){
    stop("No FASTQ files were recognised in ", fastq.dir, ".");
  }

  for(i in seq_along(samples)){
    s <- samples[i];
    if(length(list.files(file.path(prof.dir, s), pattern = "_census_stage_2\\.json$")) > 0){
      message(sprintf("[%d/%d] %s: already profiled", i, length(samples), s));
      next;
    }
    if(length(list.files(file.path(map.dir, s), pattern = "_census_stage_1\\.json$")) == 0){
      message(sprintf("[%d/%d] %s: mapping", i, length(samples), s));
      .run.meteor(tools, c("mapping", "-i", shQuote(file.path(fq.dir, s)), "-r", shQuote(catalogue.dir),
                           "-o", shQuote(map.dir), "-t", threads, if(strain) "--kf"), log.file);
    }
    message(sprintf("[%d/%d] %s: profiling", i, length(samples), s));
    .run.meteor(tools, c("profile", "-i", shQuote(file.path(map.dir, s)), "-r", shQuote(catalogue.dir),
                         "-o", shQuote(prof.dir), "-n", normalization, "--completeness", completeness), log.file);
  }

  message("Merging ", length(samples), " samples");
  dir.create(merge.dir, showWarnings = FALSE);
  .run.meteor(tools, c("merge", "-i", shQuote(prof.dir), "-r", shQuote(catalogue.dir), "-o", shQuote(merge.dir),
                       "-n", min.occurrence), log.file);
  message("Done: merged tables in ", merge.dir);
  return(invisible(merge.dir));
}

#'Build per-species carrier tables (not implemented yet)
#'@description Planned: for each sample, which species (MSPs) carry each resistance
#'gene and CAZyme family, used by MicrobiomeAnalyst for per-species host attribution.
#'Meteor2's merged tables do not give this: "_as_msp_sum" sums the abundance of the
#'species that carry an annotation without naming them.
#'@param out.dir Output folder of RunMeteorProfiling().
#'@param catalogue.dir Catalogue folder returned by DownloadMeteorCatalogue().
#'@param dbs Annotation databases to attribute.
#'@param prefix File prefix of the merged tables (default "output").
#'@return The paths of the carrier tables (invisibly).
#'@export
BuildCarrierTables <- function(out.dir, catalogue.dir, dbs = c("mustard", "resfinder", "resfinderfg", "dbcan"),
                               prefix = "output"){
  # Plan (inputs checked on the Meteor2 2.0.22 dog gut catalogue):
  # 1. Read <catalogue>_reference.json for the database folder and the file names of
  #    "msp" (msp_name, gene_id, gene_category) and of each db (gene_id, annotation).
  # 2. For each sample, read profiles/<sample>/<sample>_genes.tsv.xz (gene_id, gene_length,
  #    value) and keep genes with value > 0.
  # 3. Join detected genes with the MSP membership and the db annotation; an MSP carries an
  #    annotation in a sample when at least one of its annotated genes is detected there.
  # 4. Write merged/<prefix>_<db>_carriers.tsv: msp_name, annotation, one 0/1 column per
  #    sample (rows with no carrier in any sample dropped).
  # 5. Write merged/<prefix>_<db>_msp_share.tsv: per sample, the share of the db's gene
  #    signal (sum of value) in genes that belong to an MSP; the rest (often mobile
  #    elements) cannot be assigned to a host.
  # 6. Call it by default at the end of RunMeteorProfiling() for a full catalogue.
  stop("BuildCarrierTables() is not implemented yet.");
}

#'Call strains and strain distances with Meteor2 (not implemented yet)
#'@description Planned: runs 'meteor strain' for each sample and 'meteor tree' on all
#'of them, for strain sharing and persistence in MicrobiomeAnalyst. Needs freebayes and
#'a profiling run with strain = TRUE (Meteor2 keeps the filtered alignments only then).
#'@param out.dir Output folder of RunMeteorProfiling(strain = TRUE).
#'@param catalogue.dir Catalogue folder returned by DownloadMeteorCatalogue().
#'@param threads Number of threads.
#'@param meteor.bin Path to the meteor executable; NULL to search.
#'@return The path of the folder with the strain tables (invisibly).
#'@export
RunMeteorStrain <- function(out.dir, catalogue.dir, threads = 4, meteor.bin = NULL){
  # Plan:
  # 1. CheckMeteor(), plus freebayes next to meteor or on the PATH.
  # 2. Check that mapping/<sample>/<sample>.cram exists for each sample; if not, say that
  #    RunMeteorProfiling(strain = TRUE) must map the samples again.
  # 3. For each sample: meteor strain -i mapping/<sample> -r <catalogue> -o strain,
  #    skipping samples already done (restartable, as in RunMeteorProfiling()).
  # 4. meteor tree -i strain -o tree -t <threads>: per-species pairwise comparison tables
  #    (sample1, sample2, overlap statistics, distance, distance_category), distance
  #    matrices and trees; confirm the file names on the first run.
  # 5. Zip the comparison tables with a manifest for upload (strain tables are separate
  #    from the PackMeteorResults() zip).
  # Pairs are classified from 'distance' on the server; 'distance_category' differs
  # between Meteor2 versions and is not used.
  stop("RunMeteorStrain() is not implemented yet.");
}

#'Read merged Meteor2 tables
#'@description Reads the tables written by RunMeteorProfiling() (or by 'meteor merge').
#'@param merged.dir Folder with the merged tables.
#'@param prefix File prefix used by 'meteor merge' (default "output").
#'@return A list: species (MSP x sample matrix), taxonomy, modules (module x sample
#'matrix), module.info, completeness (long table: MSP, module, sample, completeness),
#'functions (named list of feature x sample matrices, e.g. kegg_as_msp_sum),
#'function.info (KO definitions, antimicrobial class of each resistance gene) and report.
#'Functional parts are NULL for a light catalogue.
#'@export
ReadMeteorResults <- function(merged.dir, prefix = "output"){
  f <- function(x) file.path(merged.dir, paste0(prefix, "_", x, ".tsv"));
  rd <- function(x) if(file.exists(f(x))) utils::read.delim(f(x), check.names = FALSE, quote = "", comment.char = "") else NULL;
  as.mat <- function(df){
    if(is.null(df)) return(NULL);
    m <- as.matrix(df[, -1, drop = FALSE]);
    rownames(m) <- df[[1]];
    return(m);
  }
  sp <- rd("msp");
  if(is.null(sp)){
    stop("No ", basename(f("msp")), " in ", merged.dir, ".");
  }
  cp <- rd("modules_completeness");
  if(!is.null(cp)){
    smp <- setdiff(colnames(cp), c("msp_name", "mod_id"));
    cp <- data.frame(msp_name = rep(cp$msp_name, length(smp)), mod_id = rep(cp$mod_id, length(smp)),
                     sample = rep(smp, each = nrow(cp)), completeness = unlist(cp[smp], use.names = FALSE));
    cp <- cp[!is.na(cp$completeness), ];
    rownames(cp) <- NULL;
  }
  fun.nms <- sub(paste0("^", prefix, "_(.*)\\.tsv$"), "\\1",
                 list.files(merged.dir, pattern = paste0("^", prefix, "_.*_as_(msp|genes)_sum\\.tsv$")));
  fun <- lapply(fun.nms, function(x) as.mat(rd(x)));
  names(fun) <- fun.nms;
  info <- list(kegg = rd("kegg_as_genes_sum_description"), mustard = rd("mustard_as_genes_sum_antimicrobial"));
  info <- info[!sapply(info, is.null)];
  return(list(species = as.mat(sp),
              taxonomy = rd("msp_taxonomy"),
              modules = as.mat(rd("modules")),
              module.info = rd("modules_definition"),
              completeness = cp,
              functions = if(length(fun) > 0) fun else NULL,
              function.info = if(length(info) > 0) info else NULL,
              report = rd("report")));
}

#'Pack merged Meteor2 tables for upload
#'@description Writes a zip file with the merged species, taxonomy, module, module
#'completeness and function tables (gene tables are left out) for upload to
#'MicrobiomeAnalyst.
#'@param merged.dir Folder with the merged tables.
#'@param zip.file Path of the zip file to write.
#'@param prefix File prefix used by 'meteor merge' (default "output").
#'@return The path of the zip file (invisibly).
#'@export
PackMeteorResults <- function(merged.dir, zip.file, prefix = "output"){
  fls <- list.files(merged.dir, pattern = paste0("^", prefix, "_.*\\.tsv$"));
  fls <- fls[!grepl(paste0("^", prefix, "_(genes|raw)\\.tsv$"), fls)];
  if(!any(fls == paste0(prefix, "_msp.tsv"))){
    stop("No ", prefix, "_msp.tsv in ", merged.dir, ".");
  }
  zip.file <- file.path(normalizePath(dirname(zip.file)), basename(zip.file));
  old.wd <- setwd(merged.dir);
  on.exit(setwd(old.wd));
  if(file.exists(zip.file)){
    file.remove(zip.file);
  }
  utils::zip(zip.file, fls, flags = "-q -j");
  message("Wrote ", zip.file, " (", length(fls), " tables)");
  return(invisible(zip.file));
}

# Run one meteor command; the output is appended to log.file and the tail is shown on failure
.run.meteor <- function(tools, args, log.file){
  # meteor calls bowtie2 by name, so its folder goes first on the PATH
  old.path <- Sys.getenv("PATH");
  on.exit(Sys.setenv(PATH = old.path));
  Sys.setenv(PATH = paste(c(unique(dirname(c(tools$bowtie2, tools$meteor))), old.path), collapse = .Platform$path.sep));
  cat("\n$ meteor", args, "\n", file = log.file, append = TRUE);
  tmp <- tempfile(fileext = ".log");
  status <- system2(tools$meteor, args, stdout = tmp, stderr = tmp);
  out <- readLines(tmp, warn = FALSE);
  cat(out, sep = "\n", file = log.file, append = TRUE);
  unlink(tmp);
  if(!identical(as.integer(status), 0L)){
    stop("meteor ", args[1], " failed (exit ", status, "):\n", paste(utils::tail(out, 15), collapse = "\n"),
         "\nFull log: ", log.file);
  }
  return(invisible(out));
}
