context("Taxon set enrichment: the enrichment network viewer's data (PrepareEnrichNetViewer, UpdateNetworkLayout)")

# The enrichment network viewer (the one MetaboAnalyst opens for its enrichment results) draws a JSON with two networks: the set
# network (nodes, edges) and the bipartite network of the sets and the taxa that hit them (bnodes, bedges), switched by the viewer's own
# View option. PrepareEnrichNetViewer writes that JSON, in the shape MetaboAnalyst's builder writes it, from a taxon set result
# (analSet$ora.mat and ora.hits: the over-representation and the rank-based result alike). The functions are read from the SOURCE file.
skip_if_not_installed("igraph")
skip_if_not_installed("jsonlite")
suppressMessages(library(igraph))

env <- new.env()
suppressMessages(sys.source(testthat::test_path("../../R/tsea_utils.R"), envir = env))
msgs <- character(0)
assign("AddErrMsg", function(msg) msgs <<- c(msgs, msg), envir = env)
assign("fast.write", function(dat, file, row.names = TRUE, quote = "auto") utils::write.csv(dat, file, row.names = row.names), envir = env)

# the taxon set library (what current.mset holds): the first two sets share most of their taxa, the third and fourth share their
# taxa, the fifth stands alone. Names carry the characters that break hand-made JSON.
lib <- list(
  "Set A" = c("Bacteroides fragilis", "Bacteroides ovatus", "Prevotella copri", "Escherichia coli", "Akkermansia muciniphila"),
  "Set B" = c("Bacteroides fragilis", "Bacteroides ovatus", "Prevotella copri", "Klebsiella pneumoniae"),
  "Set \"C\" & 5'-D" = c("Faecalibacterium prausnitzii", "Roseburia intestinalis", "Blautia obeum", "Eubacterium rectale"),
  "Set D" = c("Faecalibacterium prausnitzii", "Roseburia intestinalis", "Blautia obeum", "Clostridium sp. CAG:127"),
  "Set E" = c("Lactobacillus reuteri", "Lactobacillus gasseri", "Bifidobacterium longum")
)
hits <- list(
  "Set A" = c("Bacteroides fragilis", "Prevotella copri", "Escherichia coli"),
  "Set B" = c("Bacteroides fragilis", "Prevotella copri"),
  "Set \"C\" & 5'-D" = c("Faecalibacterium prausnitzii", "Roseburia intestinalis"),
  "Set D" = c("Faecalibacterium prausnitzii", "Roseburia intestinalis", "Blautia obeum"),
  "Set E" = c("Lactobacillus reuteri", "Bifidobacterium longum")
)
sets <- names(lib)
cols <- c("total", "expected", "hits", "Raw p", "Holm p", "FDR")
mat <- function(sets, hits, p, expected = rep(1.1, length(sets))) {
  m <- cbind(total = lengths(lib[sets]), expected = expected, hits = lengths(hits[sets]), "Raw p" = p, "Holm p" = pmin(1, p * 2), FDR = pmin(1, p * 1.5))
  rownames(m) <- sets
  m[order(m[, "Raw p"]), , drop = FALSE]
}

build <- function(m, set.hits = hits, set.lib = lib, then = NULL, tested = NULL) {
  d <- tempfile("tsea_enr"); dir.create(d)
  old <- setwd(d); on.exit(setwd(old))
  assign("current.mset", set.lib, envir = globalenv())
  assign(".get.mbSetObj", function(mbSetObj = NA) list(analSet = list(ora.mat = m, ora.hits = set.hits[rownames(m)], ora.tested = tested)), envir = env)
  msgs <<- character(0)
  res <- env$PrepareEnrichNetViewer(NA)
  f <- "enrichNet_tsea.json"
  out <- list(res = res, msg = tail(msgs, 1),
              json = if (file.exists(f)) jsonlite::fromJSON(f, simplifyVector = FALSE),
              csv = if (file.exists("enrichment_tsea.csv")) utils::read.csv("enrichment_tsea.csv", check.names = FALSE, row.names = 1))
  if (!is.null(then)) out$then <- then()
  out
}
ids <- function(nodes) vapply(nodes, function(n) n$id, "")
need.fn <- function(nm) {
  f <- get0(nm, envir = env, inherits = FALSE)
  if (!is.function(f)) stop(nm, " is not defined in R/tsea_utils.R")
  f
}

test_that("the viewer gets the two networks, in the shape MetaboAnalyst's builder writes them", {
  need.fn("PrepareEnrichNetViewer")
  r <- build(mat(sets, hits, c(0.001, 0.02, 0.04, 0.3, 0.01)))
  expect_equal(r$res, 1)
  d <- r$json
  expect_true(all(c("nodes", "edges", "bnodes", "bedges", "enr", "id", "sigCmpds", "hits", "analType", "pwType", "org",
                    "backgroundColor", "naviString") %in% names(d)))
  expect_equal(d$analType, "tsea")
  expect_equal(d$pwType, "taxon")
  expect_equal(d$naviString, "EnrichNetwork")

  # the set network: one node per set with what the viewer reads, an edge between the sets that share many taxa
  expect_equal(sort(ids(d$nodes)), sort(sets))
  for (n in d$nodes) {
    expect_true(all(c("id", "pwId", "label", "size", "true_size", "molType", "colorb", "colorw", "pval", "hits", "posx", "posy") %in% names(n)))
    expect_equal(n$molType, "pathway")
    expect_true(grepl("^#[0-9A-Fa-f]{6}", n$colorb) && grepl("^#[0-9A-Fa-f]{6}", n$colorw))
    expect_true(is.finite(n$posx) && is.finite(n$posy))
  }
  pairs <- vapply(d$edges, function(e) paste(sort(c(e$source, e$target)), collapse = " | "), "")
  expect_true(paste(sort(c("Set A", "Set B")), collapse = " | ") %in% pairs)
  expect_true(paste(sort(c("Set \"C\" & 5'-D", "Set D")), collapse = " | ") %in% pairs)
  expect_false(any(grepl("Set E", pairs)), info = "a set that shares nothing has no edge")
  for (e in d$edges) expect_true(!is.null(e$id) && !is.null(e$source) && !is.null(e$target) && e$source %in% sets && e$target %in% sets)

  # the bipartite network: the sets plus ALL the taxa in them (each set's members as the library holds them, the same members the set network
  # joins its sets by), a line from a set to each member; the taxa that hit are marked
  expect_equal(sort(ids(Filter(function(n) n$molType == "pathway", d$bnodes))), sort(sets))
  taxa <- ids(Filter(function(n) n$molType == "compound", d$bnodes))
  expect_equal(sort(taxa), sort(unique(unlist(lib))), info = "every member of the sets is drawn, not only the taxa that hit")
  expect_equal(length(d$bedges), sum(lengths(lib)))
  for (e in d$bedges) {
    expect_true(!is.null(e$id) && !is.null(e$source) && !is.null(e$target))
    expect_true(e$source %in% lib[[e$target]], info = "an edge joins a taxon to a set that holds it")
  }
  for (n in d$bnodes) expect_true(all(c("id", "label", "size", "true_size", "molType", "colorb", "colorw", "posx", "posy") %in% names(n)))
  for (n in Filter(function(n) n$molType == "compound", d$bnodes)) {
    expect_equal(n$isHit == 1, n$id %in% unique(unlist(hits)), info = paste("hit mark of", n$id))
    expect_true(n$tested == 1, info = "with no record of what was tested, no member is called untested")
  }

  # the viewer's table and the sets' hits
  expect_equal(unlist(d$id), rownames(mat(sets, hits, c(0.001, 0.02, 0.04, 0.3, 0.01))))
  expect_true(all(c("Hits", "Raw_p", "FDR") %in% names(d$enr[[1]])))
  expect_equal(sort(names(d$hits)), sort(sets))
  # the result table the viewer offers to download
  expect_equal(sort(rownames(r$csv)), sort(sets))
})

test_that("names with quotes, an ampersand and a prime survive", {
  need.fn("PrepareEnrichNetViewer")
  d <- build(mat(sets, hits, c(0.001, 0.02, 0.04, 0.3, 0.01)))$json
  expect_true("Set \"C\" & 5'-D" %in% ids(d$nodes))
  expect_true("Set \"C\" & 5'-D" %in% ids(d$bnodes))
})

test_that("a rank-based result (the leading-edge taxa, NES in the expected column) gives both networks too", {
  need.fn("PrepareEnrichNetViewer")
  r <- build(mat(sets, hits, c(0.002, 0.01, 0.03, 0.2, 0.05), expected = c(2.1, -1.8, 1.6, -1.2, 1.9)))
  expect_equal(r$res, 1)
  expect_equal(length(r$json$nodes), 5)
  expect_equal(length(r$json$bedges), sum(lengths(lib)))
})

test_that("sets that share no taxa give an empty edge list, and the bipartite network is still there", {
  need.fn("PrepareEnrichNetViewer")
  lib2 <- list(A = c("t1", "t2", "t3"), B = c("t4", "t5", "t6"), C = c("t7", "t8", "t9"))
  hits2 <- list(A = c("t1", "t2"), B = c("t4"), C = c("t7", "t8", "t9"))
  m <- cbind(total = rep(3, 3), expected = rep(1, 3), hits = lengths(hits2), "Raw p" = c(0.01, 0.02, 0.03), "Holm p" = rep(1, 3), FDR = rep(0.1, 3))
  rownames(m) <- names(lib2)
  r <- build(m, hits2, lib2)
  expect_equal(r$res, 1)
  expect_equal(length(r$json$nodes), 3)
  expect_equal(length(r$json$edges), 0)
  expect_equal(length(r$json$bedges), 9)
})

test_that("sets with one p value and the same number of hits still get a colour", {
  need.fn("PrepareEnrichNetViewer")
  lib2 <- list(A = c("t1", "t2", "t3"), B = c("t4", "t5", "t6"), C = c("t1", "t4", "t7"))
  hits2 <- list(A = c("t1", "t2"), B = c("t4", "t5"), C = c("t1", "t4"))
  m <- cbind(total = rep(3, 3), expected = rep(1, 3), hits = rep(2, 3), "Raw p" = rep(0.01, 3), "Holm p" = rep(1, 3), FDR = rep(0.1, 3))
  rownames(m) <- names(lib2)
  r <- build(m, hits2, lib2)
  expect_equal(r$res, 1)
  for (n in c(r$json$nodes, r$json$bnodes)) expect_true(grepl("^#[0-9A-Fa-f]{6}", n$colorb))
})

test_that("like the page's own network, only the top 50 sets are drawn", {
  need.fn("PrepareEnrichNetViewer")
  n <- 60
  nm <- sprintf("Set %02d", 1:n)
  lib2 <- setNames(lapply(1:n, function(k) sprintf("t%03d", (k * 3 + 0:5) %% 90)), nm)
  hits2 <- lapply(lib2, head, 3)
  m <- cbind(total = rep(6, n), expected = rep(1, n), hits = rep(3, n), "Raw p" = seq(0.001, 0.5, length.out = n), "Holm p" = rep(1, n), FDR = rep(0.1, n))
  rownames(m) <- nm
  r <- build(m, hits2, lib2)
  expect_equal(r$res, 1)
  expect_equal(sort(ids(r$json$nodes)), nm[1:50])
})

test_that("two sets that share taxa but no hit are joined through those taxa in the view with the taxa too", {
  need.fn("PrepareEnrichNetViewer")
  lib2 <- list(A = c("t1", "t2", "t3"), B = c("t1", "t4", "t5"))
  hits2 <- list(A = c("t2"), B = c("t4"))
  m <- cbind(total = c(3, 3), expected = c(1, 1), hits = c(1, 1), "Raw p" = c(0.01, 0.02), "Holm p" = c(1, 1), FDR = c(0.1, 0.1))
  rownames(m) <- names(lib2)
  d <- build(m, hits2, lib2)$json
  to.t1 <- vapply(Filter(function(e) e$source == "t1", d$bedges), function(e) e$target, "")
  expect_equal(sort(to.t1), c("A", "B"), info = "the taxon the two sets share is one node, and both sets reach it")
  expect_equal(length(d$edges), 1, info = "the set network joins the same two sets")
})

test_that("a member the study did not test is marked as not tested, a hit as a hit", {
  need.fn("PrepareEnrichNetViewer")
  lib2 <- list(A = c("t1", "t2", "t3", "t4"), B = c("t1", "t5", "t6", "t7"))
  hits2 <- list(A = c("t2"), B = c("t1"))
  tested2 <- list(A = c("t1", "t2"), B = c("t1", "t5"))
  m <- cbind(total = c(2, 2), expected = c(0.5, 0.5), hits = c(1, 1), "Raw p" = c(0.01, 0.02), "Holm p" = c(1, 1), FDR = c(0.1, 0.1))
  rownames(m) <- names(lib2)
  d <- build(m, hits2, lib2, tested = tested2)$json
  flag <- function(field) setNames(vapply(Filter(function(n) n$molType == "compound", d$bnodes), function(n) as.integer(n[[field]]), 1L),
                                    ids(Filter(function(n) n$molType == "compound", d$bnodes)))
  expect_equal(unname(flag("tested")[c("t1", "t2", "t3", "t4", "t5", "t6", "t7")]), c(1, 1, 0, 0, 1, 0, 0))
  expect_equal(unname(flag("isHit")[c("t1", "t2", "t3", "t5")]), c(1, 1, 0, 0))
})

test_that("a set's view says the library holds every member, and how many of them were tested", {
  need.fn("GetHTMLMetSet")
  assign("current.mset", list(S = c("a", "b", "c", "d", "e", "f")), envir = globalenv())
  assign("current.msetlib", data.frame(name = "S", reference = "ref", stringsAsFactors = FALSE), envir = globalenv())
  assign(".get.mbSetObj", function(mbSetObj = NA) list(analSet = list(ora.hits = list(S = "a"), ora.tested = list(S = c("a", "b", "c")))), envir = env)
  out <- env$GetHTMLMetSet(NA, "S")
  expect_equal(out[3], "ref")
  expect_match(out[2], "library holds 6 members of this set", fixed = TRUE)
  expect_match(out[2], "3 of them are in the list tested here", fixed = TRUE)
  expect_match(out[2], "not tested", fixed = TRUE)
  expect_match(out[2], "<font color=\"red\"><b>a</b></font>; b; c; <font color=\"grey\"", fixed = TRUE)
  # a result kept from before the tested members were recorded: the count of the library is still said, nothing is called untested
  assign(".get.mbSetObj", function(mbSetObj = NA) list(analSet = list(ora.hits = list(S = "a"))), envir = env)
  old <- env$GetHTMLMetSet(NA, "S")
  expect_match(old[2], "library holds 6 members of this set.", fixed = TRUE)
  expect_false(grepl("grey", old[2], fixed = TRUE))
})

test_that("a single enriched set cannot make a network: nothing is half-written and the message says taxon set", {
  need.fn("PrepareEnrichNetViewer")
  r <- build(mat(sets[1], hits, 0.01))
  expect_equal(r$res, 0)
  expect_match(r$msg, "Fewer than two")
  expect_match(r$msg, "taxon sets")
  expect_null(r$json)
})

test_that("a result with no table gives 0 and says so", {
  need.fn("PrepareEnrichNetViewer")
  assign(".get.mbSetObj", function(mbSetObj = NA) list(analSet = list()), envir = env)
  msgs <<- character(0)
  d <- tempfile("tsea_enr"); dir.create(d); old <- setwd(d); on.exit(setwd(old))
  expect_equal(env$PrepareEnrichNetViewer(NA), 0)
  expect_match(tail(msgs, 1), "No taxon set enrichment result")
})

test_that("the viewer's Layout option: positions for the network on screen, for either view", {
  need.fn("PrepareEnrichNetViewer"); need.fn("UpdateNetworkLayout")
  nodes.of <- function(view, algo) {
    build(mat(sets, hits, c(0.001, 0.02, 0.04, 0.3, 0.01)), then = function() {
      f <- env$UpdateNetworkLayout(algo, "layout1.json", view)
      expect_equal(f, "layout1.json")
      jsonlite::fromJSON("layout1.json", simplifyVector = FALSE)$nodes
    })$then
  }
  for (algo in c("Default", "FrR", "circle", "random", "lgl", "gopt")) {
    s <- nodes.of("enrichNet", algo)
    expect_equal(sort(ids(s)), sort(sets), info = algo)
    for (n in s) expect_true(is.finite(n$x) && is.finite(n$y), info = algo)
  }
  b <- nodes.of("enrichNet_bipartite", "Default")
  expect_equal(length(b), length(sets) + length(unique(unlist(lib))))
  expect_true(all(sets %in% ids(b)))
})
