context("Taxon set enrichment: the page's own network sizes its nodes as the enrichment network viewer does")

# The network on the results page and the one in the enrichment network viewer are the same sets. The viewer sizes a set by its number of hits
# (log10(hits + 1)) spread over 8 to 20, so the sets with the most and the fewest hits differ in size whatever the counts are. The page's network
# used log10(hits + 2) * 10, which gives 4.8 for one hit and 6.0 for three: nearly the same size for every set of a result with few hits.
skip_if_not_installed("igraph")
skip_if_not_installed("reshape2")
suppressMessages(library(igraph))
env <- new.env()
suppressMessages(sys.source(testthat::test_path("../../R/tsea_utils.R"), envir = env))

lib <- list(A = c("t1", "t2", "t3"), B = c("t1", "t2", "t4"), C = c("t1", "t2", "t5"), D = c("t1", "t3", "t6"))
assign("current.mset", lib, envir = globalenv())

sizes.of <- function(hits) {
  d <- tempfile("tsea_ov"); dir.create(d); old <- setwd(d); on.exit(setwd(old))
  env$PlotEnrichNet.Overview(hits, setNames(seq(0.01, by = 0.01, length.out = length(hits)), names(hits)))
  nodes <- jsonlite::fromJSON("tsea_network.json", simplifyVector = FALSE)$nodes
  setNames(vapply(nodes, function(n) n$size, 1), vapply(nodes, function(n) n$id, ""))
}

test_that("the set with the most hits is the largest and the one with the fewest the smallest, over the viewer's range", {
  s <- sizes.of(c(A = 1, B = 2, C = 6, D = 1))
  expect_equal(s[["A"]], 8, info = "the fewest hits: the smallest size of the viewer")
  expect_equal(s[["C"]], 20, info = "the most hits: the largest size of the viewer")
  expect_equal(s[["D"]], s[["A"]], info = "the same number of hits, the same size")
  expect_true(s[["B"]] > 8 && s[["B"]] < 20)
  expect_equal(s[["B"]], 8 + (log10(3) - log10(2)) / (log10(7) - log10(2)) * 12, tolerance = 1e-4, info = "log10(hits + 1), as the viewer")
})

test_that("a result whose sets all have the same number of hits draws them all at one size", {
  s <- sizes.of(c(A = 2, B = 2, C = 2, D = 2))
  expect_equal(unname(s), rep(12, 4))
})
