#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
})

source("scripts/Internals/R/03_prepare_germline_scores.R")

default_gwas_path <- "data/raw/Testing/34737426-GCST90043906-EFO_0000403.h.tsv"
default_gene_loc_path <- "data/raw/NCBI38/NCBI38.gene.loc"
default_reg_loc_path <- "data/raw/GeneHancer/2026-01-26_UCSC_all_unfiltered_reg_elements.loc"
default_reference_bfile <- "data/raw/g1000_eur/g1000_eur"
default_sample_size <- 456348L

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

test_run_magma_step1_annotation <- function(
  gwas_path = default_gwas_path,
  gene_loc_path = default_gene_loc_path,
  nrows = 50000L
) {
  message("Reading GWAS file for MAGMA step 1 test: ", gwas_path)
  gwas_raw <- fread(gwas_path, nrows = nrows)

  output_prefix <- tempfile(pattern = "magma_step1_test_")

  result <- run_magma_step1_annotation(
    gwas_sumstats = gwas_raw,
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
  genes_out_path <- tempfile(pattern = "magma_genes_out_", fileext = ".genes.out")
  zstat_out_path <- tempfile(pattern = "magma_zstat_", fileext = ".tsv")

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
  genes_out_path <- tempfile(pattern = "magma_reg_out_", fileext = ".genes.out")

  genes_out_dt <- data.table(
    GENE = c("GH01J000013", "GH01J000021"),
    GENE_NAME = c("GH01J000013", "GH01J000021"),
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
  if (!identical(as.character(zstat_dt$feature_id[[1]]), "GH01J000013")) {
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

  message("Reading GWAS file for MAGMA step 2 test: ", gwas_path)
  gwas_raw <- fread(gwas_path, nrows = nrows)

  step1_prefix <- tempfile(pattern = "magma_step2_prereq_")
  step1_result <- run_magma_step1_annotation(
    gwas_sumstats = gwas_raw,
    gene_loc_path = feature_loc_path,
    output_prefix = step1_prefix
  )

  output_prefix <- tempfile(pattern = "magma_step2_test_")

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

  gwas_raw <- fread(default_gwas_path, nrows = 50000L)
  output_prefix <- tempfile(pattern = "magma_pipeline_test_")

  result <- run_magma_feature_scoring_pipeline(
    gwas_sumstats = gwas_raw,
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
  test_run_magma_step1_annotation()
  test_extract_magma_zstat()
  test_extract_magma_feature_zstat_regulatory()
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
  test_that("MAGMA germline scoring works for genes and regulatory elements", {
    results <- run_all_germline_tests(print_scores = TRUE)

    expect_true(is.data.frame(results$gene$zstat))
    expect_true(is.data.frame(results$regulatory_element$zstat))
    expect_gt(nrow(results$gene$zstat), 0L)
    expect_gt(nrow(results$regulatory_element$zstat), 0L)
    expect_identical(names(results$gene$zstat), c("feature_id", "zstat"))
    expect_identical(names(results$regulatory_element$zstat), c("feature_id", "zstat"))
  })
}

if (sys.nframe() == 0) {
  main()
}
