#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
  library(data.table)
})

source("R/user_api.R")

plot_test_dir <- "data/processed/test_outputs/plotting"
dir.create(plot_test_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(plot_test_dir, c(
  "diffusion_plot.pdf",
  "germline_gene_scores_plot.pdf"
)))

test_rank_plotting_works <- function() {
  germline_path <- "data/processed/germline_gene_scores.tsv"
  if (!file.exists(germline_path)) {
    skip("Real germline score table is not available.")
  }

  output_path <- file.path(plot_test_dir, "germline_gene_scores_rank_plot.pdf")
  bundle <- plot_scores(
    table = fread(germline_path),
    plot_file_path = output_path,
    test_tail = "one_tailed",
    label_features = c("MYC", "BCL2"),
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
}

test_reg_rank_plotting_works <- function() {
  reg_path <- "data/processed/germline_reg_scores.tsv"
  if (!file.exists(reg_path)) {
    skip("Real germline regulatory score table is not available.")
  }

  output_path <- file.path(plot_test_dir, "germline_reg_scores_rank_plot.pdf")
  bundle <- plot_scores(
    table = fread(reg_path),
    which = "reg_scores",
    plot_file_path = output_path,
    test_tail = "one_tailed",
    label_features = c("MYC", "BCL2"),
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
  expect_true(any(bundle$objects$plot_data$should_label))
  expect_true(all(bundle$objects$plot_data$matched_label[bundle$objects$plot_data$should_label] %in% c("MYC", "BCL2")))
}

test_volcano_plotting_works <- function() {
  tbl <- data.table(
    gene_id = c("TP53", "MYC", "BCL2", "PAX5"),
    p_value = c(1e-8, 2e-4, 0.02, 0.3),
    zstat = c(5.8, 3.6, -2.4, 1.0)
  )

  output_path <- file.path(plot_test_dir, "somatic_gene_scores_volcano_plot.pdf")
  bundle <- plot_scores(
    table = tbl,
    plot_file_path = output_path,
    test_tail = "two_tailed",
    label_features = c("TP53", "MYC"),
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
}

test_plotting_requires_output_path_when_saving <- function() {
  expect_error(
    plot_scores(
      table = data.frame(feature_id = c("A", "B"), zstat = c(1, 2)),
      save_plot = TRUE
    ),
    "`plot_file_path` must be provided",
    fixed = TRUE
  )
}

main <- function() {
  test_that("plot_scores makes a rank plot on a real one-tailed table", {
    test_rank_plotting_works()
  })

  test_that("plot_scores makes a germline regulatory rank plot with backend gene labels", {
    test_reg_rank_plotting_works()
  })

  test_that("plot_scores makes a volcano plot on a two-tailed table", {
    test_volcano_plotting_works()
  })

  test_that("plot_scores validates save arguments", {
    test_plotting_requires_output_path_when_saving()
  })
}

if (sys.nframe() == 0) {
  main()
}
