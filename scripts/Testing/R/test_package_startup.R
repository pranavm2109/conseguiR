#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
})

source("R/zzz.R")

package_startup_uses_lymphoma_graph_env <- function() {
  old_opts <- options(
    conseguiR.python = NULL,
    conseguiR.conda_env = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  .onLoad(NULL, NULL)
  python_path <- getOption("conseguiR.python")
  expect_null(python_path)
}

main <- function() {
  test_that("package startup avoids configuring Python implicitly", {
    package_startup_uses_lymphoma_graph_env()
  })
}

if (sys.nframe() == 0) {
  main()
}
