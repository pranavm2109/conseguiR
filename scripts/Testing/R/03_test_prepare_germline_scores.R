#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
})

source("scripts/Internals/R/03_prepare_germline_scores.R")

default_gwas_path <- "data/raw/Testing/34737426-GCST90043906-EFO_0000403.h.tsv"
default_gene_loc_path <- "data/raw/NCBI38/NCBI38.gene.loc"
default_reg_loc_path <- "data/processed/GRCh38-cCREs.loc"
default_reference_bfile <- "data/raw/g1000_eur/g1000_eur"
default_sample_size <- 456348L
default_germline_test_output_dir <- "data/processed/test_outputs/germline"
default_run_full_magma_tests <- identical(Sys.getenv("CONSEGUIR_RUN_FULL_MAGMA_TESTS", unset = "0"), "1")

make_germline_test_path <- function(stem, ext = "") {
  dir.create(default_germline_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_germline_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

make_tiny_gwas_fixture <- function() {
  data.table(
    hm_rsid = c("rs1", "rs2", "rs3"),
    hm_variant_id = c("1_100_A_G", "1_200_C_T", "2_300_G_A"),
    hm_chrom = c("1", "1", "2"),
    hm_pos = c(100L, 200L, 300L),
    p_value = c(0.01, 0.2, 0.05)
  )
}

find_reference_bfile <- function() {
  candidates <- c(
    default_reference_bfile,
    "data/raw/g1000_eur/g1000_eur",
    "data/raw/Reference/1000G_EUR_Phase3_plink/1000G.EUR.QC",
    "data/raw/Reference/g1000_eur/g1000_eur",
    "data/raw/LDREF/g1000_eur/g1000_eur",
    "tools/reference/g1000_eur/g1000_eur"
  )

  candidates <- unique(candidates[!is.na(candidates)])

  for (prefix in candidates) {
    if (is.null(prefix)) {
      next
    }

    required_files <- paste0(prefix, c(".bed", ".bim", ".fam"))
    if (all(file.exists(required_files))) {
      return(prefix)
    }
  }

  NULL
}

test_prepare_magma_gwas_cache_reuses_cached_outputs <- function(
  gwas_path = make_tiny_gwas_fixture()
) {
  output_prefix <- make_germline_test_path("magma_cached_inputs")

  first <- prepare_magma_gwas_cache(
    gwas_sumstats = gwas_path,
    cache_prefix = output_prefix,
    reuse_existing = FALSE
  )

  expect_true(file.exists(first$snp_loc_path))
  expect_true(file.exists(first$pval_path))
  expect_false(first$reused_existing)

  first_snp_mtime <- file.info(first$snp_loc_path)$mtime
  first_pval_mtime <- file.info(first$pval_path)$mtime

  Sys.sleep(1)

  second <- prepare_magma_gwas_cache(
    gwas_sumstats = gwas_path,
    cache_prefix = output_prefix,
    reuse_existing = TRUE
  )

  expect_true(second$reused_existing)
  expect_equal(file.info(second$snp_loc_path)$mtime, first_snp_mtime)
  expect_equal(file.info(second$pval_path)$mtime, first_pval_mtime)

  invisible(second)
}

test_run_magma_step1_annotation_reuses_existing_annotation <- function() {
  output_prefix <- make_germline_test_path("magma_step1_reuse")
  snp_loc_path <- paste0(output_prefix, ".snp_loc.tsv")
  pval_path <- paste0(output_prefix, ".pval.tsv")
  annot_path <- paste0(output_prefix, ".genes.annot")

  tiny_gwas <- data.table(
    hm_rsid = c("rs1", "rs2"),
    hm_variant_id = c("1_100_A_G", "1_200_C_T"),
    hm_chrom = c("1", "1"),
    hm_pos = c(100L, 200L),
    p_value = c(0.01, 0.2)
  )

  fwrite(data.table(V1 = c("rs1", "rs2"), V2 = c("1", "1"), V3 = c(100L, 200L)), snp_loc_path, sep = "\t", col.names = FALSE)
  fwrite(data.table(SNP = c("rs1", "rs2"), P = c(0.01, 0.2)), pval_path, sep = "\t")
  writeLines(c("GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs", "TEST\t1\t1\t1000\t2\trs1,rs2"), annot_path)

  result <- run_magma_step1_annotation(
    gwas_sumstats = tiny_gwas,
    gene_loc_path = default_gene_loc_path,
    output_prefix = output_prefix,
    reuse_prepared_inputs = TRUE,
    reuse_existing_annotation = TRUE
  )

  expect_true(isTRUE(result$reused_existing_annotation))
  expect_true(file.exists(result$annot_path))
  expect_null(result$command)
  invisible(result)
}

test_run_magma_step2_gene_analysis_reuses_existing_output <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  output_prefix <- make_germline_test_path("magma_step2_reuse")
  gene_annot_path <- make_germline_test_path("magma_gene_annot", ".genes.annot")
  pval_path <- make_germline_test_path("magma_pval", ".tsv")
  genes_out_path <- paste0(output_prefix, ".genes.out")

  writeLines(c("GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs", "TEST\t1\t1\t1000\t2\trs1,rs2"), gene_annot_path)
  fwrite(data.table(SNP = c("rs1", "rs2"), P = c(0.01, 0.2)), pval_path, sep = "\t")
  fwrite(data.table(GENE = "TEST", GENE_NAME = "TEST", ZSTAT = 2.1, P = 0.03), genes_out_path, sep = "\t")

  result <- run_magma_step2_gene_analysis(
    gene_annot_path = gene_annot_path,
    pval_path = pval_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = default_sample_size,
    reuse_existing_analysis = TRUE
  )

  expect_true(isTRUE(result$reused_existing_analysis))
  expect_true(file.exists(result$genes_out_path))
  expect_null(result$command)
  invisible(result)
}

test_run_magma_step1_annotation <- function(
  gwas_path = default_gwas_path,
  gene_loc_path = default_gene_loc_path,
  nrows = 50000L
) {
  output_prefix <- make_germline_test_path("magma_step1_test")

  result <- run_magma_step1_annotation(
    gwas_sumstats = gwas_path,
    gene_loc_path = gene_loc_path,
    output_prefix = output_prefix
  )

  message("MAGMA command used: ", result$command)

  expected_files <- c(
    result$snp_loc_path,
    result$pval_path,
    result$annot_path
  )

  if (!all(file.exists(expected_files))) {
    stop(
      "MAGMA step 1 test failed: not all expected output files were written. Missing: ",
      paste(expected_files[!file.exists(expected_files)], collapse = ", ")
    )
  }

  if (file.info(result$snp_loc_path)$size <= 0) {
    stop("MAGMA step 1 test failed: snp-loc file is empty.")
  }
  if (file.info(result$pval_path)$size <= 0) {
    stop("MAGMA step 1 test failed: pval file is empty.")
  }
  if (file.info(result$annot_path)$size <= 0) {
    stop("MAGMA step 1 test failed: .genes.annot file is empty.")
  }

  snp_loc_dt <- fread(result$snp_loc_path, header = FALSE)
  pval_dt <- fread(result$pval_path, header = TRUE)
  annot_preview <- readLines(result$annot_path, n = 5L)

  if (ncol(snp_loc_dt) != 3L) {
    stop("MAGMA step 1 test failed: snp-loc file does not have exactly 3 columns.")
  }
  if (!identical(names(pval_dt), c("SNP", "P"))) {
    stop("MAGMA step 1 test failed: pval file headers are not exactly SNP and P.")
  }
  if (any(is.na(snp_loc_dt$V1) | snp_loc_dt$V1 == "")) {
    stop("MAGMA step 1 test failed: written snp-loc file has missing SNP identifiers.")
  }
  if (any(is.na(snp_loc_dt$V2) | snp_loc_dt$V2 == "")) {
    stop("MAGMA step 1 test failed: written snp-loc file has missing chromosomes.")
  }
  if (any(is.na(snp_loc_dt$V3))) {
    stop("MAGMA step 1 test failed: written snp-loc file has missing positions.")
  }
  if (any(is.na(pval_dt$SNP) | pval_dt$SNP == "")) {
    stop("MAGMA step 1 test failed: written pval file has missing SNP identifiers.")
  }
  if (any(is.na(pval_dt$P))) {
    stop("MAGMA step 1 test failed: written pval file has missing p-values.")
  }
  if (length(annot_preview) == 0L) {
    stop("MAGMA step 1 test failed: .genes.annot preview is empty.")
  }

  message("MAGMA step 1 annotation test passed.")
  message("Output prefix: ", output_prefix)
  invisible(result)
}

test_extract_magma_zstat <- function() {
  genes_out_path <- make_germline_test_path("magma_genes_out", ".genes.out")
  zstat_out_path <- make_germline_test_path("magma_zstat", ".tsv")

  genes_out_dt <- data.table(
    GENE = c("1", "2"),
    GENE_NAME = c("GENEA", "GENEB"),
    NSNPS = c(12L, 8L),
    NPARAM = c(10L, 7L),
    ZSTAT = c(2.5, -1.75),
    P = c(0.0124, 0.0801)
  )

  fwrite(genes_out_dt, genes_out_path, sep = "\t")

  zstat_dt <- extract_magma_zstat(
    genes_out_path = genes_out_path,
    output_path = zstat_out_path
  )

  if (!file.exists(zstat_out_path) || file.info(zstat_out_path)$size <= 0) {
    stop("ZSTAT extraction test failed: output file was not written.")
  }
  if (!identical(names(zstat_dt), c("gene_id", "zstat"))) {
    stop("ZSTAT extraction test failed: unexpected output columns.")
  }
  if (!identical(as.character(zstat_dt$gene_id[[1]]), "GENEA")) {
    stop("ZSTAT extraction test failed: gene_id was not extracted correctly.")
  }
  if (!isTRUE(all.equal(as.numeric(zstat_dt$zstat[[1]]), 2.5))) {
    stop("ZSTAT extraction test failed: ZSTAT value was not extracted correctly.")
  }

  message("MAGMA ZSTAT extraction test passed.")
  invisible(zstat_dt)
}

test_extract_magma_feature_zstat_regulatory <- function() {
  genes_out_path <- make_germline_test_path("magma_reg_out", ".genes.out")

  genes_out_dt <- data.table(
    GENE = c("EH38E0080197", "EH38E2084302"),
    GENE_NAME = c("EH38E0080197", "EH38E2084302"),
    ZSTAT = c(1.25, -0.88),
    P = c(0.211, 0.379)
  )

  fwrite(genes_out_dt, genes_out_path, sep = "\t")

  zstat_dt <- extract_magma_feature_zstat(
    genes_out_path = genes_out_path,
    feature_type = "regulatory_element"
  )

  if (!identical(names(zstat_dt), c("feature_id", "zstat"))) {
    stop("Regulatory ZSTAT extraction test failed: unexpected output columns.")
  }
  if (!identical(as.character(zstat_dt$feature_id[[1]]), "EH38E0080197")) {
    stop("Regulatory ZSTAT extraction test failed: feature_id was not extracted correctly.")
  }

  message("MAGMA regulatory-element ZSTAT extraction test passed.")
  invisible(zstat_dt)
}

test_run_magma_step2_gene_analysis <- function(
  reference_bfile = find_reference_bfile(),
  gwas_path = default_gwas_path,
  feature_loc_path = default_gene_loc_path,
  feature_type = "gene",
  nrows = 50000L
) {
  if (is.null(reference_bfile)) {
    message("Skipping MAGMA step 2 test: no PLINK reference bfile found in the repository.")
    return(invisible(NULL))
  }

  step1_prefix <- make_germline_test_path("magma_step2_prereq")
  step1_result <- run_magma_step1_annotation(
    gwas_sumstats = gwas_path,
    gene_loc_path = feature_loc_path,
    output_prefix = step1_prefix
  )

  output_prefix <- make_germline_test_path("magma_step2_test")

  result <- run_magma_step2_gene_analysis(
    gene_annot_path = step1_result$annot_path,
    pval_path = step1_result$pval_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = default_sample_size
  )

  if (!file.exists(result$genes_out_path) || file.info(result$genes_out_path)$size <= 0) {
    stop("MAGMA step 2 test failed: .genes.out was not created.")
  }

  genes_out_dt <- fread(result$genes_out_path)
  if (!all(c("GENE", "ZSTAT") %in% names(genes_out_dt))) {
    stop("MAGMA step 2 test failed: .genes.out is missing GENE or ZSTAT.")
  }

  message("MAGMA step 2 ", feature_type, " analysis test passed.")
  invisible(result)
}

test_run_magma_feature_scoring_pipeline <- function(
  reference_bfile = find_reference_bfile(),
  feature_loc_path = default_gene_loc_path,
  feature_type = "gene",
  print_scores = TRUE
) {
  if (is.null(reference_bfile)) {
    message("Skipping full MAGMA pipeline test: no PLINK reference bfile found in the repository.")
    return(invisible(NULL))
  }

  output_prefix <- make_germline_test_path("magma_pipeline_test")

  result <- run_magma_feature_scoring_pipeline(
    gwas_sumstats = default_gwas_path,
    feature_loc_path = feature_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    feature_type = feature_type,
    sample_size = default_sample_size
  )

  if (!file.exists(result$zstat_output_path) || file.info(result$zstat_output_path)$size <= 0) {
    stop("Full MAGMA pipeline test failed: zstat output file was not written.")
  }
  if (!is.data.frame(result$zstat) || nrow(result$zstat) == 0L) {
    stop("Full MAGMA pipeline test failed: extracted zstat table is empty.")
  }
  if (!identical(names(result$zstat), c("feature_id", "zstat"))) {
    stop("Full MAGMA pipeline test failed: extracted zstat table is missing expected columns.")
  }

  if (isTRUE(print_scores)) {
    message("Final ", feature_type, " germline scores (MAGMA ZSTAT):")
    print(result$zstat)
  }

  message("Full MAGMA ", feature_type, " scoring pipeline test passed.")
  invisible(result)
}

run_all_germline_tests <- function(print_scores = TRUE) {
  test_prepare_magma_gwas_cache_reuses_cached_outputs()
  test_run_magma_step1_annotation_reuses_existing_annotation()
  test_run_magma_step2_gene_analysis_reuses_existing_output()
  test_extract_magma_zstat()
  test_extract_magma_feature_zstat_regulatory()
}

run_full_germline_integration_tests <- function(print_scores = TRUE) {
  test_run_magma_step1_annotation()
  test_run_magma_step2_gene_analysis(
    feature_loc_path = default_gene_loc_path,
    feature_type = "gene"
  )
  test_run_magma_step2_gene_analysis(
    feature_loc_path = default_reg_loc_path,
    feature_type = "regulatory_element"
  )

  gene_result <- test_run_magma_feature_scoring_pipeline(
    feature_loc_path = default_gene_loc_path,
    feature_type = "gene",
    print_scores = print_scores
  )
  reg_result <- test_run_magma_feature_scoring_pipeline(
    feature_loc_path = default_reg_loc_path,
    feature_type = "regulatory_element",
    print_scores = print_scores
  )

  invisible(list(
    gene = gene_result,
    regulatory_element = reg_result
  ))
}

main <- function() {
  test_that("MAGMA germline smoke tests pass", {
    expect_no_error(run_all_germline_tests(print_scores = FALSE))
  })

  if (isTRUE(default_run_full_magma_tests)) {
    test_that("MAGMA germline scoring works for genes and regulatory elements", {
      results <- run_full_germline_integration_tests(print_scores = TRUE)

      expect_true(is.data.frame(results$gene$zstat))
      expect_true(is.data.frame(results$regulatory_element$zstat))
      expect_gt(nrow(results$gene$zstat), 0L)
      expect_gt(nrow(results$regulatory_element$zstat), 0L)
      expect_identical(names(results$gene$zstat), c("feature_id", "zstat"))
      expect_identical(names(results$regulatory_element$zstat), c("feature_id", "zstat"))
    })
  } else {
    message(
      "Skipping full live MAGMA integration tests. Set CONSEGUIR_RUN_FULL_MAGMA_TESTS=1 to enable them."
    )
  }
}

if (sys.nframe() == 0) {
  main()
}
