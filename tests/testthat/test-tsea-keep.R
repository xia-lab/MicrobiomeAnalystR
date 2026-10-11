context("Taxon set enrichment: a complete list's two results are kept apart, with the members each was tested with")

# A complete list is tested two ways, and the page shows one at a time: the result in state is kept under a name (KeepTseaResult) and put back
# when its tab is chosen (ShowTseaResult). The members a set was tested with (ora.tested) belong to the result: the set's view marks the others
# as not tested, and the enrichment network viewer does too. They were not kept, so the first result was shown with the second one's members
# (a set that only the first reported had none, and its view said nothing of what was tested).
env <- new.env()
suppressMessages(sys.source(testthat::test_path("../../R/tsea_utils.R"), envir = env))
state <- new.env()
assign(".get.mbSetObj", function(x) state$mb, envir = env)
assign(".set.mbSetObj", function(x) { state$mb <- x; invisible(x) }, envir = env)

result <- function(sets, tested.of, tag) {
  m <- cbind(total = 2, expected = 0.1, hits = 1, "Raw p" = 0.1, "Holm p" = 1, FDR = 1)[rep(1, length(sets)), , drop = FALSE]
  rownames(m) <- sets
  list(ora.mat = m, ora.hits = setNames(lapply(sets, function(s) "h"), sets), ora.tested = setNames(lapply(sets, tested.of), sets),
       ora.background = list(mode = tag))
}

test_that("each result comes back with its own tested members, and a set only the first one reported keeps them", {
  top <- result(c("KCNIP4", "Other"), function(s) c("t1", "t2"), "top")
  rank <- result("Other", function(s) c("r1", "r2", "r3"), "rank")
  state$mb <- list(analSet = top)
  env$KeepTseaResult(NA, "top")
  state$mb$analSet[names(rank)] <- rank          # the second test overwrites the result in state
  env$KeepTseaResult(NA, "rank")

  env$ShowTseaResult(NA, "top")
  a <- state$mb$analSet
  expect_equal(rownames(a$ora.mat), c("KCNIP4", "Other"))
  expect_equal(a$ora.tested[["KCNIP4"]], c("t1", "t2"), info = "the set the rank-based result did not report is read with the members it was tested with")
  expect_equal(a$ora.tested[["Other"]], c("t1", "t2"), info = "and a set both reported is read with its own result's members, not the other's")
  expect_equal(a$ora.background$mode, "top")

  env$ShowTseaResult(NA, "rank")
  b <- state$mb$analSet
  expect_equal(rownames(b$ora.mat), "Other")
  expect_equal(b$ora.tested[["Other"]], c("r1", "r2", "r3"))
  expect_equal(b$ora.background$mode, "rank")
})

test_that("a kept result without tested members does not borrow another result's", {
  old <- result("S", function(s) "x", "old")
  old$ora.tested <- NULL
  state$mb <- list(analSet = c(old, list(tsea.keep = list(old = old[c("ora.mat", "ora.hits", "ora.background")]))))
  state$mb$analSet$ora.tested <- list(S = c("another", "result's", "members"))
  env$ShowTseaResult(NA, "old")
  expect_null(state$mb$analSet$ora.tested, info = "no record of what was tested is better than another result's")
})

test_that("a result that was not kept says so", {
  state$mb <- list(analSet = list())
  expect_error(env$ShowTseaResult(NA, "top"), "no such result")
})
