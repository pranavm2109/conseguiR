#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
})

source("R/zzz.R")

package_startup_initializes_only_pkg_root <- function() {
  old_opts <- options(
    conseguiR.conda_env = NULL
  )
  on.exit(options(old_opts), add = TRUE)

  .conseguiR_state$pkg_root <- NULL
  .onLoad(NULL, NULL)
  expect_null(.conseguiR_state$pkg_root)
}

main <- function() {
  test_that("package startup avoids runtime side effects", {
    package_startup_initializes_only_pkg_root()
  })
}

if (sys.nframe() == 0) {
  main()
}
