#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(testthat)
  library(data.table)
})

source("R/zzz.R")
source("R/user_api.R")
.onLoad(libname = ".", pkgname = "conseguiR")

plot_test_dir <- file.path(tempdir(), "conseguiR_test_outputs", "plotting")
dir.create(plot_test_dir, recursive = TRUE, showWarnings = FALSE)
unlink(file.path(plot_test_dir, c(
  "diffusion_plot.pdf",
  "germline_gene_scores_plot.pdf"
)))

test_rank_plotting_works <- function() {
  tbl <- data.table(
    gene_id = c("MYC", "BCL2", "PAX5", "IRF4"),
    zstat = c(5.1, 3.4, 1.2, -0.7)
  )

  output_path <- file.path(plot_test_dir, "germline_gene_scores_rank_plot.pdf")
  bundle <- plot_scores(
    table = tbl,
    plot_file_path = output_path,
    test_tail = "one_tailed",
    label_features = c("MYC", "BCL2"),
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
}

test_reg_rank_plotting_works <- function() {
  tbl <- data.table(
    reg_elem_id = c("EH38E0080197", "EH38E2084302", "EH38E3951312", "EH38E2776544"),
    zstat = c(4.3, 2.8, 1.1, -0.4)
  )

  output_path <- file.path(plot_test_dir, "germline_reg_scores_rank_plot.pdf")
  bundle <- plot_scores(
    table = tbl,
    which = "reg_scores",
    plot_file_path = output_path,
    test_tail = "one_tailed",
    label_features = c("EH38E0080197", "EH38E2084302"),
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
  expect_true(any(bundle$objects$plot_data$should_label))
  expect_true(all(bundle$objects$plot_data$feature_id_plot[bundle$objects$plot_data$should_label] %in% c("EH38E0080197", "EH38E2084302")))
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
    test_tail = "one_tailed",
    plot_mode = "volcano",
    label_features = c("TP53", "MYC"),
    save_plot = TRUE
  )

  expect_s3_class(bundle, "conseguiR_bundle")
  expect_true(file.exists(output_path))
}

test_plot_scores_infers_semantics_from_bundle_metadata <- function() {
  score_bundle <- structure(
    list(
      bundle_type = "somatic_gene_scores",
      objects = list(
        gene_scores = data.table(
          gene_id = c("TP53", "MYC", "BCL2", "PAX5"),
          p_value = c(1e-8, 2e-4, 0.02, 0.3),
          zstat = c(5.8, 3.6, -2.4, 1.0)
        )
      ),
      output_paths = list(),
      config = list(
        test_tail = "one_tailed",
        default_plot_mode = "volcano"
      )
    ),
    class = c("conseguiR_somatic_gene_scores_bundle", "conseguiR_bundle", "list")
  )

  plot_bundle <- plot_scores(
    scores = score_bundle,
    which = "gene_scores",
    save_plot = FALSE
  )

  expect_identical(plot_bundle$config$test_tail, "one_tailed")
  expect_identical(plot_bundle$config$plot_mode, "volcano")
  expect_identical(plot_bundle$config$requested_test_tail, "auto")
  expect_identical(plot_bundle$config$requested_plot_mode, "auto")
}

test_volcano_plotting_can_drop_tukey_outliers <- function() {
  tbl <- data.table(
    gene_id = c("KMT2D", "TP53", "MYC", "BCL2", "PAX5", "IRF4", "EZH2"),
    p_value = c(1e-80, 1e-8, 2e-4, 0.02, 0.03, 0.04, 0.06),
    zstat = c(18, 5.8, 3.6, 2.2, 2.0, 1.9, 1.5)
  )

  kept <- plot_scores(
    table = tbl,
    test_tail = "one_tailed",
    plot_mode = "volcano",
    drop_tukey_outliers = FALSE,
    save_plot = FALSE
  )

  trimmed <- plot_scores(
    table = tbl,
    test_tail = "one_tailed",
    plot_mode = "volcano",
    drop_tukey_outliers = TRUE,
    save_plot = FALSE
  )

  expect_equal(nrow(kept$objects$plot_data), nrow(tbl))
  expect_lt(nrow(trimmed$objects$plot_data), nrow(kept$objects$plot_data))
  expect_false("KMT2D" %in% trimmed$objects$plot_data$feature_id_plot)
}

test_volcano_plotting_is_unclipped_by_default <- function() {
  tbl <- data.table(
    gene_id = c("KMT2D", "TNFRSF14", "B2M"),
    p_value = c(1e-80, 1e-10, 1e-6),
    zstat = c(37.065788, 6.673909, 7.778698)
  )

  bundle <- plot_scores(
    table = tbl,
    test_tail = "one_tailed",
    plot_mode = "volcano",
    save_plot = FALSE
  )

  plot_dt <- bundle$objects$plot_data
  expect_equal(plot_dt$z_display, plot_dt$z_plot)
  expect_equal(plot_dt$neglog10_p_display, plot_dt$neglog10_p)
}

test_volcano_plotting_can_clip_display_when_requested <- function() {
  tbl <- data.table(
    gene_id = c("KMT2D", "TNFRSF14", "B2M"),
    p_value = c(1e-80, 1e-10, 1e-6),
    zstat = c(37.065788, 6.673909, 7.778698)
  )

  bundle <- plot_scores(
    table = tbl,
    test_tail = "one_tailed",
    plot_mode = "volcano",
    clip_extreme_display = TRUE,
    save_plot = FALSE
  )

  plot_dt <- bundle$objects$plot_data
  expect_lt(plot_dt[gene_id == "KMT2D", z_display], plot_dt[gene_id == "KMT2D", z_plot])
  expect_lte(max(plot_dt$neglog10_p_display), max(plot_dt$neglog10_p))
}

test_plotting_requires_output_path_when_saving <- function() {
  expect_error(
    plot_scores(
      table = data.frame(feature_id = c("A", "B"), zstat = c(1, 2)),
      save_plot = TRUE
    ),
    regexp = "`plot_file_path` must be provided",
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

  test_that("plot_scores can make a volcano plot independently of one-tailed score semantics", {
    test_volcano_plotting_works()
  })

  test_that("plot_scores can infer score semantics from bundle metadata", {
    test_plot_scores_infers_semantics_from_bundle_metadata()
  })

  test_that("plot_scores can remove extreme volcano outliers with Tukey filtering", {
    test_volcano_plotting_can_drop_tukey_outliers()
  })

  test_that("plot_scores keeps volcano coordinates truthful by default", {
    test_volcano_plotting_is_unclipped_by_default()
  })

  test_that("plot_scores can still clip volcano display when requested", {
    test_volcano_plotting_can_clip_display_when_requested()
  })

  test_that("plot_scores validates save arguments", {
    test_plotting_requires_output_path_when_saving()
  })
}

if (sys.nframe() == 0) {
  main()
}
