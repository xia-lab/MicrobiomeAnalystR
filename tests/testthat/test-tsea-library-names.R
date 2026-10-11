context("Taxon set enrichment: a set whose name is the text NA keeps its name when a library is loaded")

# The libraries are read as tables, and the reader makes the text "NA" a missing value: a set named NA lost its name, answered to no name in
# the result (a row called NA, a set view that said the library held 0 members of it). The host-genetics library had one, a gene symbol
# that was missing. The set keeps the text.
env <- new.env()
suppressMessages(sys.source(testthat::test_path("../../R/tsea_utils.R"), envir = env))
dir <- tempfile("tsea_lib"); dir.create(file.path(dir, "libs", "tsea"), recursive = TRUE)
writeLines(c("name,member,reference",
             "AARS,Bacteroides; Prevotella,r1",
             "NA,Faecalibacterium; Lachnospira; Eubacterium,r2",
             "KCNIP4,Capnocytophaga; Lactococcus; Blautia,r3"), file.path(dir, "libs", "tsea", "tsea_host_int.csv"))
assign(".on.public.web", TRUE, envir = env)
assign("rpath", paste0(dir, "/"), envir = env)
assign(".readDataTable", function(p) utils::read.csv(p, stringsAsFactors = FALSE, check.names = FALSE), envir = env)
assign(".canonical.taxon.name", function(x) x, envir = env)
assign("AddErrMsg", function(msg) invisible(), envir = env)

test_that("a set named NA is loaded with that name, with its members", {
  expect_true(env$.tsea.load.lib("host_int"))
  sets <- get("current.mset", envir = globalenv())
  expect_false(anyNA(names(sets)), info = "no set is left without a name")
  expect_equal(names(sets), c("AARS", "NA", "KCNIP4"))
  expect_equal(sets[["NA"]], c("Faecalibacterium", "Lachnospira", "Eubacterium"))
  lib <- get("current.msetlib", envir = globalenv())
  expect_equal(lib[[1]], c("AARS", "NA", "KCNIP4"))
})
