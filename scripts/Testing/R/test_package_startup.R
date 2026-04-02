#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
})

source("R/zzz.R")

package_startup_uses_lymphoma_graph_env <- function() {
  .onLoad(NULL, NULL)
  python_path <- getOption("conseguiR.python")
  if (is.null(python_path) || !nzchar(python_path)) {
    skip("Package startup did not configure a Python interpreter path.")
  }

  expect_true(file.exists(python_path))
  expect_true(grepl("lymphoma_graph_env", python_path, fixed = TRUE))
}

main <- function() {
  test_that("package startup configures the lymphoma_graph_env Python interpreter", {
    package_startup_uses_lymphoma_graph_env()
  })
}

if (sys.nframe() == 0) {
  main()
}
