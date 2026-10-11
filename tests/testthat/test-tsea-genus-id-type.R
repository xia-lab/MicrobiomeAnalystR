context("Taxon set enrichment: a genus list is mapped by CrossReferencing / SpeciesMappingExact as a first-class level")

# The levels of a taxon list are genus, species, strain and mixed (the pages and the library variants run in that order). The name-mapping
# table holds the names of every rank in one column (taxa), so a genus name is mapped the way a species or a strain name is. The function is
# read from the SOURCE file, with the name table (the package reads microbe_db_new.rds) and the synonyms stubbed.

env <- new.env()
suppressMessages(sys.source(testthat::test_path("../../R/tsea_utils.R"), envir = env))

db <- data.frame(
  taxa      = c("Lactobacillus", "Lactobacillus acidophilus", "Bacteroides", "Escherichia coli K-12"),
  canonical = c("Lactobacillus", "Lactobacillus acidophilus", "Bacteroides", "Escherichia coli K-12"),
  species   = c("", "Lactobacillus acidophilus", "", "Escherichia coli"),
  genus     = c("Lactobacillus", "Lactobacillus", "Bacteroides", "Escherichia"),
  NCBITAX   = c("1578", "1579", "816", "83333"),
  GOLDMAPID = c("", "", "", ""),
  GOLD_ID   = c("", "", "", ""),
  stringsAsFactors = FALSE)
assign(".read.microbiomeanalyst.lib.rds", function(...) db, envir = env)
assign(".canonical.taxon.name", function(x) x, envir = env)

genus.names <- c("Lactobacillus", "bacteroides", "Nonexistus")
assign("name.map", list(), envir = env)     # CrossReferencing starts it; the mapping fills it

test_that("a genus list is mapped: its names are found in the table, an unknown name is not", {
  res <- env$SpeciesMappingExact(genus.names, "genus")
  expect_equal(nrow(res), 3L)
  expect_equal(as.character(res$Match), c("Lactobacillus", "Bacteroides", ""))
  expect_equal(as.character(res$NCBI_Taxonomy_ID[1]), "<a href=https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id=1578 target='_blank'>1578</a>")
  expect_equal(unname(env$name.map$match.state), c(1, 1, 0))
})

test_that("a genus list is mapped exactly as the same names are as a mixed list", {
  as.genus <- env$SpeciesMappingExact(genus.names, "genus")
  as.mixed <- env$SpeciesMappingExact(genus.names, "mixed")
  expect_identical(as.genus, as.mixed)
})

test_that("CrossReferencing keeps the genus ID type it was given", {
  mbSet <- list(dataSet = list(species = genus.names), analSet = list())
  assign(".get.mbSetObj", function(x) mbSet, envir = env)
  assign(".set.mbSetObj", function(x) { mbSet <<- x; invisible(x) }, envir = env)
  out <- env$CrossReferencing(NA, "genus")
  expect_equal(out$dataSet$q.type, "genus")
  expect_equal(as.character(out$analSet$resTable$Match), c("Lactobacillus", "Bacteroides", ""))
})
