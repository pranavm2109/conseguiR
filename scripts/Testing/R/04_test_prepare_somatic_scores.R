#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
  library(GenomicRanges)
})

source("scripts/Internals/R/04_prepare_somatic_scores.R")

default_somatic_path <- "data/raw/Testing/2026-01-09_no_CLL_lymph_only_pcawg_maf_tcga_order_hg38.maf"
default_reg_ref_path <- "data/raw/Testing/2026-01-26_UCSC_all_unfiltered_reg_elements.loc"
default_fishhook_covariate_path <- "data/raw/Testing/2026-01-26_all_reg_elems_sample_level_mut_frac_comparison_bet_only_memory_b_normal_and_non_cll_malig_b_cells.rds"
default_dndscv_refdb <- "data/raw/Testing/RefCDS_human_GRCh38.p12.rda"

make_dummy_maf <- function() {
  data.table(
    Tumor_Sample_Barcode = c("S1", "S1", "S2"),
    Chromosome = c("1", "1", "2"),
    Start_Position = c(1001L, 1050L, 2001L),
    End_Position = c(1001L, 1050L, 2001L),
    Reference_Allele = c("A", "G", "C"),
    Tumor_Seq_Allele2 = c("T", "A", "T"),
    Hugo_Symbol = c("TP53", "KRAS", "BRAF")
  )
}

find_dndscv_refdb <- function() {
  candidates <- c(
    default_dndscv_refdb,
    "data/raw/Testing/RefCDS_human_GRCh38.p12.rda",
    "data/raw/dndscv/hg38_refcds.rda",
    "data/raw/dndscv/RefCDS_hg38.rda",
    "data/raw/Reference/dndscv/hg38_refcds.rda",
    "data/raw/Reference/dndscv/RefCDS_hg38.rda"
  )

  candidates <- unique(candidates[!is.na(candidates)])
  hits <- candidates[file.exists(candidates)]

  if (length(hits) == 0L) {
    return(NULL)
  }

  hits[[1]]
}

get_negative_somatic_test_cases <- function() {
  list(
    list(
      path = "data/raw/Testing/somatic_maf_missing_sample.tsv",
      label = "missing sample identifier",
      expected_text = "sample identifier"
    ),
    list(
      path = "data/raw/Testing/somatic_maf_missing_chromosome.tsv",
      label = "missing chromosome",
      expected_text = "chromosome"
    ),
    list(
      path = "data/raw/Testing/somatic_maf_missing_start.tsv",
      label = "missing start position",
      expected_text = "start position"
    ),
    list(
      path = "data/raw/Testing/somatic_maf_missing_end.tsv",
      label = "missing end position",
      expected_text = "end position"
    ),
    list(
      path = "data/raw/Testing/somatic_maf_missing_ref.tsv",
      label = "missing reference allele",
      expected_text = "reference allele"
    ),
    list(
      path = "data/raw/Testing/somatic_maf_missing_alt.tsv",
      label = "missing alternate allele",
      expected_text = "alternate allele"
    )
  )
}

test_run_dndscv_gene_scoring_negative_missing_refdb <- function() {
  skip_if_not_installed("dndscv")

  maf <- make_dummy_maf()

  expect_error(
    run_dndscv_gene_scoring(maf = maf, refdb = NULL),
    expected = "`refdb` is required for dndscv scoring",
    fixed = TRUE
  )
}

test_run_dndscv_gene_scoring_negative_bad_refdb <- function() {
  skip_if_not_installed("dndscv")

  maf <- make_dummy_maf()

  expect_error(
    run_dndscv_gene_scoring(maf = maf, refdb = "data/raw/Testing/not_a_real_refdb.rda"),
    expected = "dndscv `refdb` file does not exist",
    fixed = TRUE
  )
}

test_run_fishhook_reg_scoring_negative_missing_covariate_column <- function() {
  skip_if_not_installed("fishHook")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

  maf <- make_dummy_maf()
  bad_covariates <- data.table(
    reg_elem_id = c("REG1", "REG2"),
    some_other_column = c(0.1, 0.2)
  )

  expect_error(
    run_fishhook_reg_scoring(
      maf = maf,
      reg_ref_path = "data/raw/Testing/reg_elements_valid.loc",
      fishhook_covariate_data = bad_covariates
    ),
    expected = "fishHook covariate data is missing required columns",
    fixed = TRUE
  )
}

test_run_fishhook_reg_scoring_negative_bad_reg_reference <- function() {
  skip_if_not_installed("fishHook")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

  maf <- make_dummy_maf()

  expect_error(
    run_fishhook_reg_scoring(
      maf = maf,
      reg_ref_path = "data/raw/Testing/does_not_exist.loc"
    ),
    expected = "Regulatory element reference file does not exist",
    fixed = TRUE
  )
}

test_prepare_somatic_inputs_from_dummy_maf <- function() {
  maf <- make_dummy_maf()

  maf_validated <- validate_somatic_maf(maf)
  dndscv_ready <- prepare_dndscv_input(maf)
  fishhook_ready <- prepare_fishhook_input(maf)
  fishhook_events <- make_fishhook_event_granges(maf)
  fishhook_hypotheses <- make_fishhook_hypothesis_granges("data/raw/Testing/reg_elements_valid.loc")

  expect_identical(names(maf_validated), c("sample_id", "chromosome", "start_position", "end_position", "ref", "alt"))
  expect_identical(names(dndscv_ready), c("sampleID", "chr", "pos", "ref", "mut"))
  expect_identical(
    names(fishhook_ready),
    c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Tumor_Seq_Allele2")
  )
  expect_s4_class(fishhook_events, "GRanges")
  expect_s4_class(fishhook_hypotheses, "GRanges")
  expect_equal(length(fishhook_events), nrow(maf))
  expect_gt(length(fishhook_hypotheses), 0L)
}

test_validate_somatic_maf_negative_case <- function(path, label, expected_text) {
  maf_raw <- fread(path)

  expect_error(
    validate_somatic_maf(maf_raw),
    expected = expected_text,
    fixed = TRUE,
    info = label
  )
}

test_extract_dndscv_gene_scores <- function() {
  dndscv_mock <- data.table(
    gene_name = c("TP53", "KRAS"),
    pallsubs_cv = c(1e-8, 0.02),
    wall_cv = c(3.2, 1.4)
  )

  scores <- extract_dndscv_gene_scores(dndscv_mock)

  expect_identical(names(scores), c("gene_id", "p_value", "zstat"))
  expect_identical(as.character(scores$gene_id[[1]]), "TP53")
  expect_true(all(scores$p_value > 0))
  expect_true(all(scores$zstat > 0))
  invisible(scores)
}

test_extract_fishhook_reg_scores <- function() {
  fishhook_mock <- data.table(
    reg_elem_id = c("GH01J000013", "GH01J000021"),
    p = c(1e-4, 0.03),
    effectsize = c(2.1, 1.2)
  )

  scores <- extract_fishhook_reg_scores(fishhook_mock)

  expect_identical(names(scores), c("reg_elem_id", "p_value", "zstat"))
  expect_identical(as.character(scores$reg_elem_id[[1]]), "GH01J000013")
  expect_true(all(scores$p_value > 0))
  expect_true(all(scores$zstat > 0))
  invisible(scores)
}

test_somatic_extreme_scores_are_capped <- function() {
  capped_from_p <- compute_signed_z_from_p(c(0, 1e-400, 1e-20))
  expect_true(all(is.finite(capped_from_p)))

  fishhook_mock <- data.table(
    reg_elem_id = c("GH01J000013", "GH01J000021"),
    p = c(1e-4, 0.03),
    zscore = c(Inf, -Inf)
  )

  scores <- extract_fishhook_reg_scores(fishhook_mock)
  expect_true(all(is.finite(scores$zstat)))
}

test_run_somatic_scoring_pipeline_with_mock_outputs <- function(print_scores = TRUE) {
  maf <- fread(default_somatic_path, nrows = 1000L)

  dndscv_mock <- data.table(
    gene_name = c("TP53", "KRAS"),
    pallsubs_cv = c(1e-8, 0.02),
    wall_cv = c(3.2, 1.4)
  )

  fishhook_mock <- data.table(
    reg_elem_id = c("GH01J000013", "GH01J000021"),
    p = c(1e-4, 0.03),
    effectsize = c(2.1, 1.2)
  )

  result <- run_somatic_scoring_pipeline(
    maf = maf,
    dndscv_result = dndscv_mock,
    fishhook_result = fishhook_mock
  )

  expect_true(is.data.frame(result$gene_scores))
  expect_true(is.data.frame(result$reg_scores))
  expect_identical(names(result$gene_scores), c("gene_id", "p_value", "zstat"))
  expect_identical(names(result$reg_scores), c("reg_elem_id", "p_value", "zstat"))
  expect_gt(nrow(result$gene_scores), 0L)
  expect_gt(nrow(result$reg_scores), 0L)

  if (isTRUE(print_scores)) {
    message("Final somatic gene scores from mocked dndscv statistics:")
    print(result$gene_scores)
    message("Final somatic regulatory-element scores from mocked fishHook statistics:")
    print(result$reg_scores)
  }

  invisible(result)
}

test_run_dndscv_gene_scoring_live <- function(refdb = find_dndscv_refdb()) {
  skip_if_not_installed("dndscv")
  skip_if(is.null(refdb), "No hg38 dndscv RefCDS file found for the real hg38 MAF.")

  maf <- fread(default_somatic_path, nrows = 1000L)
  result <- run_dndscv_gene_scoring(maf, refdb = refdb)

  expect_true(is.data.frame(result))
  expect_identical(names(result), c("gene_id", "p_value", "zstat"))
  expect_gt(nrow(result), 0L)
  message("Live dndscv gene scores:")
  print(result)
  invisible(result)
}

test_run_fishhook_reg_scoring_live <- function() {
  skip_if_not_installed("fishHook")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

  maf <- fread(default_somatic_path, nrows = 1000L)
  fishhook_covariate_data <- readRDS(default_fishhook_covariate_path)
  result <- run_fishhook_reg_scoring(
    maf = maf,
    reg_ref_path = default_reg_ref_path,
    fishhook_covariate_data = fishhook_covariate_data
  )

  expect_true(is.data.frame(result))
  expect_identical(names(result), c("reg_elem_id", "p_value", "zstat"))
  expect_gt(nrow(result), 0L)
  message("Live fishHook regulatory-element scores:")
  print(result)
  invisible(result)
}

test_run_somatic_scoring_pipeline_live <- function(
  print_scores = TRUE,
  refdb = find_dndscv_refdb()
) {
  skip_if_not_installed("dndscv")
  skip_if_not_installed("fishHook")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")
  skip_if(is.null(refdb), "No hg38 dndscv RefCDS file found for the real hg38 MAF.")

  maf <- fread(default_somatic_path, nrows = 1000L)
  fishhook_covariate_data <- readRDS(default_fishhook_covariate_path)

  result <- run_somatic_scoring_pipeline(
    maf = maf,
    reg_ref_path = default_reg_ref_path,
    fishhook_covariate_data = fishhook_covariate_data,
    refdb = refdb
  )

  expect_true(is.data.frame(result$gene_scores))
  expect_true(is.data.frame(result$reg_scores))
  expect_identical(names(result$gene_scores), c("gene_id", "p_value", "zstat"))
  expect_identical(names(result$reg_scores), c("reg_elem_id", "p_value", "zstat"))
  expect_gt(nrow(result$gene_scores), 0L)
  expect_gt(nrow(result$reg_scores), 0L)

  if (isTRUE(print_scores)) {
    message("Final somatic gene scores from live dndscv:")
    print(result$gene_scores)
    message("Final somatic regulatory-element scores from live fishHook:")
    print(result$reg_scores)
  }

  invisible(result)
}

main <- function() {
  test_that("dummy MAF can be converted into somatic scoring inputs", {
    test_prepare_somatic_inputs_from_dummy_maf()
  })

  test_that("somatic validators report clear errors for malformed MAF files", {
    for (case in get_negative_somatic_test_cases()) {
      test_validate_somatic_maf_negative_case(
        path = case$path,
        label = case$label,
        expected_text = case$expected_text
      )
    }
  })

  test_that("dndscv runner reports clear errors for missing or bad refdb inputs", {
    test_run_dndscv_gene_scoring_negative_missing_refdb()
    test_run_dndscv_gene_scoring_negative_bad_refdb()
  })

  test_that("fishHook runner reports clear errors for bad regulatory or covariate inputs", {
    test_run_fishhook_reg_scoring_negative_missing_covariate_column()
    test_run_fishhook_reg_scoring_negative_bad_reg_reference()
  })

  test_that("somatic score extraction works for dndscv genes", {
    test_extract_dndscv_gene_scores()
  })

  test_that("somatic score extraction works for fishHook regulatory elements", {
    test_extract_fishhook_reg_scores()
  })

  test_that("somatic extreme z-scores are capped to finite values", {
    test_somatic_extreme_scores_are_capped()
  })

  test_that("somatic scoring pipeline combines mocked gene and regulatory scores", {
    test_run_somatic_scoring_pipeline_with_mock_outputs(print_scores = TRUE)
  })

  test_that("live dndscv gene scoring runs when dndscv is installed", {
    test_run_dndscv_gene_scoring_live()
  })

  test_that("live fishHook regulatory scoring runs when fishHook is installed", {
    test_run_fishhook_reg_scoring_live()
  })

  test_that("live somatic scoring pipeline runs on the real MAF when dndscv and fishHook are installed", {
    test_run_somatic_scoring_pipeline_live(print_scores = TRUE)
  })
}

if (sys.nframe() == 0) {
  main()
}
