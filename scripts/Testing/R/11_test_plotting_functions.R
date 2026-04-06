#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
})

source("R/user_api.R")

plot_test_dir <- "data/processed/test_outputs/plotting"
dir.create(plot_test_dir, recursive = TRUE, showWarnings = FALSE)

test_score_plotting_works <- function() {
  germline_path <- "data/processed/germline_gene_scores.tsv"
  if (!file.exists(germline_path)) {
    skip("Real germline score table is not available.")
  }

  output_path <- file.path(plot_test_dir, "germline_gene_scores_plot.pdf")
  bundle <- plot_scores(
    table = data.table::fread(germline_path),
    plot_file_path = output_path,
    plot_type = "top_bar",
    top_n = 20L,
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
}

test_diffusion_plotting_works <- function() {
  diffusion_path <- "data/processed/gene_reg_graph_diffusion_all_genes.tsv"
  if (!file.exists(diffusion_path)) {
    skip("Real diffusion table is not available.")
  }

  output_path <- file.path(plot_test_dir, "diffusion_plot.pdf")
  bundle <- plot_diffusion(
    table = data.table::fread(diffusion_path),
    plot_file_path = output_path,
    plot_type = "ranked_points",
    top_n = 50L,
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
}

test_score_plotting_requires_output_path_when_saving <- function() {
  expect_error(
    plot_scores(
      table = data.frame(feature_id = c("A", "B"), score = c(1, 2)),
      save_plot = TRUE
    ),
    "`plot_file_path` must be provided",
    fixed = TRUE
  )
}

test_diffusion_plotting_requires_output_path_when_saving <- function() {
  expect_error(
    plot_diffusion(
      table = data.frame(gene_name = c("A", "B"), prize = c(1, 2)),
      save_plot = TRUE
    ),
    "`plot_file_path` must be provided",
    fixed = TRUE
  )
}

main <- function() {
  test_that("plot_scores works on a real score table", {
    test_score_plotting_works()
  })

  test_that("plot_diffusion works on a real diffusion table", {
    test_diffusion_plotting_works()
  })

  test_that("plot_scores validates save arguments", {
    test_score_plotting_requires_output_path_when_saving()
  })

  test_that("plot_diffusion validates save arguments", {
    test_diffusion_plotting_requires_output_path_when_saving()
  })
}

if (sys.nframe() == 0) {
  main()
}
