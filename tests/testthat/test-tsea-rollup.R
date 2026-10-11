context("Taxon set enrichment: a list with species and strain names is rolled up to genus for a genus-level library")

# A genus-level library holds genus names only, so a list that names species and strains matched almost nothing of it (the Complete ranked
# list example: 36 names, 5 of them genera). For such a library the list is rolled up to its genera: one row per genus, a genus keeping the
# score of its strongest member; the background (the taxa measured) is rolled up the same way; the result says it. The functions are read
# from the SOURCE file with the name table stubbed.
skip_if_not_installed("fgsea")

env <- new.env()
suppressMessages(sys.source(testthat::test_path("../../R/tsea_utils.R"), envir = env))
msgs <- character(0)
assign("AddErrMsg", function(msg) msgs <<- c(msgs, msg), envir = env)
assign("fast.write", function(...) invisible(), envir = env)
assign("ExportResultMatArrow", function(...) invisible(), envir = env)
assign(".canonical.taxon.name", function(x) x, envir = env)
assign("run_func_via_microservice", function(func, args = list(), timeout_sec = 60, allow_fallback = TRUE) do.call(func, args), envir = env)

lib <- list(
  S1 = c("Bacteroides", "Prevotella", "Roseburia", "Faecalibacterium"),
  S2 = c("Lactobacillus", "Bifidobacterium", "Akkermansia", "Escherichia"),
  S3 = c("Bacteroides", "Lactobacillus", "Fusobacterium", "Clostridium"),
  S4 = c("Streptococcus", "Veillonella", "Roseburia", "Dialister"),
  S5 = c("Escherichia", "Klebsiella", "Enterococcus", "Fusobacterium"))
assign("current.mset", lib, envir = globalenv())
assign("current.msetlib", data.frame(name = names(lib), reference = "r", stringsAsFactors = FALSE), envir = globalenv())

mixed <- c("Bacteroides fragilis", "Bacteroides ovatus", "Prevotella copri", "Roseburia", "Faecalibacterium prausnitzii",
           "Lactobacillus reuteri", "Lactobacillus gasseri", "Bifidobacterium longum NT5", "Escherichia coli K-12",
           "Fusobacterium nucleatum", "Clostridium sp. CAG:127", "Streptococcus mutans", "Lachnospiraceae")
scores <- c(3, -4, 2, 1.5, 1, 0.5, 0.4, 0.3, -0.2, -3.5, -0.1, 0.05, 0.01)

mbset <- function(tset.type, species = mixed, sc = scores, tset.used = tset.type) {
  list(dataSet = list(species = species, scores = sc, q.type = "taxa", tset.type = tset.type, tset.used = tset.used,
                      tsea.input = "complete", tsea.cutoff.by = "top", tsea.cutoff = 0.3),
       analSet = list())
}
state <- new.env()
use <- function(mb) {
  state$mb <- mb
  assign(".get.mbSetObj", function(x) state$mb, envir = env)
  assign(".set.mbSetObj", function(x) { state$mb <- x; invisible(x) }, envir = env)
  invisible(NULL)
}
name.map <- function(mb) data.frame(query = mb$dataSet$species, Strain = mb$dataSet$species, stringsAsFactors = FALSE)

test_that("a species or a strain name is rolled up to its genus, any other name is left as it is", {
  g <- env$.tsea.to.level(c("Bacteroides fragilis", "Escherichia coli K-12", "[Eubacterium] siraeum", "Candidatus Accumulibacter phosphatis",
                            "Clostridium sp. CAG:127", "Roseburia", "Prevotella 9", "Lachnospiraceae", NA), "genus")
  expect_equal(g, c("Bacteroides", "Escherichia", "Eubacterium", "Accumulibacter", "Clostridium", "Roseburia", "Prevotella 9", "Lachnospiraceae", NA))
})

test_that("a strain name is rolled up to its species for a species-level library; a species, a genus and a family stay", {
  sp <- env$.tsea.to.level(c("Escherichia coli K-12", "Streptococcus phocae subsp. salmonis", "[Eubacterium] siraeum ATCC 29", "Bacteroides fragilis",
                             "Clostridium sp. CAG:127", "Roseburia", "Lachnospiraceae", NA), "species")
  expect_equal(sp, c("Escherichia coli", "Streptococcus phocae", "Eubacterium siraeum", "Bacteroides fragilis", "Clostridium sp. CAG:127",
                     "Roseburia", "Lachnospiraceae", NA))
  expect_equal(env$.tsea.to.level("Escherichia coli K-12", "strain"), "Escherichia coli K-12", info = "a strain-level library takes the names as they are")
})

test_that("a species-level library rolls the strains of a list up to their species, not to their genus", {
  mb <- mbset("host_int_species", species = c("Escherichia coli K-12", "Escherichia coli O157", "Bacteroides fragilis", "Roseburia"),
              sc = c(1, -2, 0.5, 0.1)); use(mb)
  r <- env$.tsea.roll.up(mb, name.map(mb))
  expect_equal(r$level, "species")
  expect_equal(r$nm.map$Strain, c("Escherichia coli", "Escherichia coli", "Bacteroides fragilis", "Roseburia"))
  expect_equal(r$set$dataSet$scores[1], -2)
  expect_equal(r$from, 4)
  expect_equal(r$to, 3)
})

test_that("a genus-level library rolls the list up: one row per genus, the strongest score, the other names left alone", {
  mb <- mbset("host_int_genus"); use(mb)
  r <- env$.tsea.roll.up(mb, name.map(mb))
  expect_false(is.null(r))
  expect_equal(r$level, "genus")
  expect_equal(r$from, 13)
  expect_equal(r$to, 11)
  first <- !duplicated(r$nm.map$Strain)
  expect_equal(sum(first), 11)
  expect_true("Lachnospiraceae" %in% r$nm.map$Strain, info = "a family has no genus to roll up to and stays")
  sc <- setNames(r$set$dataSet$scores[first], r$nm.map$Strain[first])
  expect_equal(sc[["Bacteroides"]], -4, info = "the member with the largest size of score keeps its score, sign included")
  expect_equal(sc[["Lactobacillus"]], 0.5)
  expect_equal(mb$dataSet$scores, scores, info = "the list as typed is not changed: the same run can be repeated against another library")
})

test_that("a list of p values keeps the smallest p of a genus, and a list's p and FDR columns too", {
  mb <- mbset("host_int_genus", species = c("A x", "A y", "B z"), sc = c(0.2, 0.01, 0.5)); use(mb)
  mb$dataSet$p.values <- c(0.2, 0.01, 0.5); mb$dataSet$fdr.values <- c(0.4, 0.05, 0.6)
  r <- env$.tsea.roll.up(mb, name.map(mb))
  first <- !duplicated(r$nm.map$Strain)
  expect_equal(r$set$dataSet$scores[first], c(0.01, 0.5))
  expect_equal(r$set$dataSet$p.values[first], c(0.01, 0.5))
  expect_equal(r$set$dataSet$fdr.values[first], c(0.05, 0.6))
})

test_that("nothing is rolled up for a library of the mixed level, a species library, or a list that has nothing to roll up", {
  mb <- mbset("host_int"); use(mb)
  expect_null(env$.tsea.roll.up(mb, name.map(mb)))
  mb <- mbset("host_int_strain"); use(mb)
  expect_null(env$.tsea.roll.up(mb, name.map(mb)))
  mb <- mbset("host_int_species", species = c("Bacteroides fragilis", "Prevotella copri", "Roseburia"), sc = c(1, 2, 3)); use(mb)
  expect_null(env$.tsea.roll.up(mb, name.map(mb)))
  mb <- mbset("host_int_genus", species = c("Bacteroides", "Prevotella", "Roseburia"), sc = c(1, 2, 3)); use(mb)
  expect_null(env$.tsea.roll.up(mb, name.map(mb)))
})

test_that("a list tested against the mixed library is read as typed", {
  assign("current.mset", list(M1 = c("Bacteroides fragilis", "Prevotella copri", "Roseburia", "Fusobacterium nucleatum"),
                              M2 = c("Lactobacillus reuteri", "Bifidobacterium longum NT5", "Escherichia coli K-12", "Lachnospiraceae"),
                              M3 = c("Bacteroides ovatus", "Prevotella copri", "Lactobacillus gasseri", "Streptococcus mutans")), envir = globalenv())
  on.exit(assign("current.mset", lib, envir = globalenv()))
  mb <- mbset("host_int"); use(mb)
  msgs <<- character(0)
  r <- env$CalculateHyperScore(NA)
  expect_true(is.list(r), info = paste(msgs, collapse = " | "))
  expect_null(state$mb$analSet$ora.background$rollup.from)
  expect_equal(state$mb$analSet$ora.background$list.size, 13)
})

test_that("the over-representation of a mixed list against a genus library tests its genera, and the record says it was rolled up", {
  mb <- mbset("host_int_genus"); use(mb)
  msgs <<- character(0)
  r <- env$CalculateHyperScore(NA)
  expect_true(is.list(r), info = paste(msgs, collapse = " | "))
  bg <- state$mb$analSet$ora.background
  expect_equal(bg$rollup.from, 13)
  expect_equal(bg$rollup.to, 11)
  expect_equal(bg$rollup.level, "genus")
  expect_true(bg$rollup.scored)
  expect_equal(bg$list.size, 11, info = "the list is now its 11 rows: ten genera and a family")
  expect_equal(bg$background, 10, info = "the taxa of the list that the library holds: the ten genera, not the family")
  expect_equal(bg$query.input, 3, info = "the top 30 percent of the 11, by the size of the score: Bacteroides, Fusobacterium, Prevotella")
  expect_true(all(unlist(state$mb$analSet$ora.hits) %in% c("Bacteroides", "Fusobacterium", "Prevotella")))
  expect_true("Bacteroides" %in% state$mb$analSet$ora.hits[["S3"]])
})

test_that("the taxa measured are rolled up with the list: a background of species counts as its genera", {
  mb <- mbset("host_int_genus", species = c("Bacteroides fragilis", "Prevotella copri", "Fusobacterium nucleatum"), sc = NULL)
  mb$dataSet$tsea.input <- NULL; use(mb)
  measured <- c("Bacteroides ovatus", "Prevotella copri", "Roseburia intestinalis", "Fusobacterium nucleatum", "Lactobacillus reuteri",
                "Faecalibacterium prausnitzii", "Escherichia coli", "Clostridium difficile")
  r <- env$CalculateHyperScore(NA, background = measured)
  expect_true(is.list(r), info = paste(msgs, collapse = " | "))
  bg <- state$mb$analSet$ora.background
  expect_equal(bg$background, 8, info = "the eight measured names are eight genera in the library; a species name matched no genus before")
  expect_equal(bg$rollup.from, 3)
  expect_equal(bg$rollup.to, 3)
})

test_that("the rank-based test of a mixed list against a genus library ranks its genera", {
  mb <- mbset("host_int_genus"); use(mb)
  msgs <<- character(0)
  r <- env$CalculateGseaScore(NA)
  expect_true(is.list(r), info = paste(msgs, collapse = " | "))
  expect_equal(state$mb$analSet$ora.background$rollup.to, 11)
  members <- unique(unlist(state$mb$analSet$ora.tested))
  expect_true(all(members %in% unique(unlist(lib))), info = "every member tested is a genus")
})

test_that("the level a list is read at: its finest level, the share of each, and whether its names are of one level", {
  lv <- function(nms) { mb <- list(dataSet = list(species = nms), analSet = list()); use(mb); env$GetTseaListLevel(NA) }
  g <- lv(c("Bacteroides", "Prevotella", "Roseburia"))
  expect_equal(g[1:2], c("genus", "genus"))
  expect_match(g[3], "3 genus", fixed = TRUE)
  sp <- lv(c("Bacteroides fragilis", "Prevotella copri", "Roseburia intestinalis"))
  expect_equal(sp[1:2], c("species", "species"))
  st <- lv(c("Escherichia coli K-12", "Escherichia coli O157", "Bacteroides fragilis NT5"))
  expect_equal(st[1:2], c("strain", "strain"))
  mixed.gs <- lv(c("Bacteroides", "Prevotella", "Roseburia", "Bacteroides fragilis", "Prevotella copri", "Roseburia intestinalis"))
  expect_equal(mixed.gs[1:2], c("species", "mixed"), info = "genus and species names: the finest level is species")
  mixed.all <- lv(c(mixed))
  expect_equal(mixed.all[1:2], c("strain", "mixed"))
  expect_match(mixed.all[3], "genus", fixed = TRUE)
  expect_match(mixed.all[3], "species", fixed = TRUE)
  expect_match(mixed.all[3], "strain", fixed = TRUE)
  fam <- lv(c("Lachnospiraceae", "Ruminococcaceae", "Bacteroidaceae"))
  expect_equal(fam[1:2], c("higher", "higher"))
  expect_equal(lv(character(0))[1:2], c("mixed", "mixed"), info = "nothing to read: the list is offered every level")
})
