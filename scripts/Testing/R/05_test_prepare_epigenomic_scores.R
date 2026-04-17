#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
})

source("scripts/Internals/R/05_prepare_epigenomic_scores.R")

default_epigenomic_track_dir <- "data/raw/Testing"
default_reg_ref_path <- "data/processed/GRCh38-cCREs.loc"
default_broken_bigwig_path <- "data/raw/Testing/broken_signal_track.bw"
default_epigenomic_test_output_dir <- "data/processed/test_outputs/epigenomic"

select_epigenomic_test_bigwigs <- function(track_dir) {
  bw_files <- list_epigenomic_track_files(track_dir = track_dir)
  keep <- !grepl("(_BL_|_FL_|broken_signal_track)", basename(bw_files))
  bw_files[keep]
}

make_epigenomic_test_path <- function(stem, ext = "") {
  dir.create(default_epigenomic_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_epigenomic_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

make_epigenomic_live_fixture <- function(n_tracks = 3L, n_reg_elements = 500L) {
  bw_files <- list_epigenomic_track_files(
    track_dir = default_epigenomic_track_dir,
    min_tracks = 3L
  )
  bw_files <- bw_files[!grepl("(_BL_|_FL_|broken_signal_track)", basename(bw_files))]
  bw_files <- bw_files[seq_len(min(n_tracks, length(bw_files)))]

  reg_dt <- fread(default_reg_ref_path, header = FALSE)
  reg_subset <- reg_dt[seq_len(min(n_reg_elements, nrow(reg_dt)))]
  reg_subset_path <- make_epigenomic_test_path("epigenomic_reg_subset", ".loc")
  fwrite(reg_subset, reg_subset_path, sep = "\t", col.names = FALSE)

  list(
    bw_files = bw_files,
    reg_ref_path = reg_subset_path
  )
}

test_run_epigenomic_reg_scoring_negative_too_few_tracks <- function() {
  fixture <- make_epigenomic_live_fixture(n_tracks = 2L, n_reg_elements = 100L)

  expect_error(
    run_epigenomic_reg_scoring(
      reg_ref_path = fixture$reg_ref_path,
      bw_files = fixture$bw_files,
      min_tracks = 3L,
      drop_mhc = TRUE,
      transform = "log1p",
      return_diagnostics = FALSE
    ),
    expected = "At least 3 epigenomic bigWig tracks are required",
    fixed = TRUE
  )
}

test_list_epigenomic_track_files_live <- function() {
  bw_files <- select_epigenomic_test_bigwigs(default_epigenomic_track_dir)

  expect_true(length(bw_files) > 0L)
  expect_true(all(file.exists(bw_files)))
  expect_false(any(grepl("broken_signal_track", basename(bw_files), fixed = TRUE)))
  invisible(bw_files)
}

test_validate_epigenomic_tracks_negative_broken_bigwig <- function() {
  expect_error(
    validate_epigenomic_tracks(default_broken_bigwig_path),
    expected = "Failed to validate bigWig file",
    fixed = TRUE
  )
}

test_load_regulatory_elements_for_epigenomic_scores_live <- function() {
  reg_gr <- load_regulatory_elements_for_epigenomic_scores(
    reg_ref_path = default_reg_ref_path,
    drop_mhc = TRUE
  )

  expect_s4_class(reg_gr, "GRanges")
  expect_gt(length(reg_gr), 0L)
  expect_true("reg_elem_id" %in% names(S4Vectors::mcols(reg_gr)))
  invisible(reg_gr)
}

test_run_epigenomic_reg_scoring_live <- function(print_scores = TRUE) {
  fixture <- make_epigenomic_live_fixture()

  result <- run_epigenomic_reg_scoring(
    reg_ref_path = fixture$reg_ref_path,
    bw_files = fixture$bw_files,
    drop_mhc = TRUE,
    transform = "log1p",
    return_diagnostics = FALSE
  )

  expect_true(is.data.frame(result))
  expect_identical(names(result), c("reg_elem_id", "zstat"))
  expect_gt(nrow(result), 0L)

  if (isTRUE(print_scores)) {
    message("Live epigenomic regulatory-element scores:")
    print(result[1:10])
  }

  invisible(result)
}

main <- function() {
  test_that("epigenomic track listing works on the testing directory", {
    test_list_epigenomic_track_files_live()
  })

  test_that("broken bigWig files fail validation clearly", {
    test_validate_epigenomic_tracks_negative_broken_bigwig()
  })

  test_that("regulatory elements load for epigenomic scoring", {
    test_load_regulatory_elements_for_epigenomic_scores_live()
  })

  test_that("epigenomic scoring fails clearly when fewer than three bigWigs are supplied", {
    test_run_epigenomic_reg_scoring_negative_too_few_tracks()
  })

  test_that("live epigenomic regulatory scoring runs on the testing tracks", {
    test_run_epigenomic_reg_scoring_live(print_scores = TRUE)
  })
}

if (sys.nframe() == 0) {
  main()
}
