#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

source("scripts/Internals/R/00_harmonise_and_validate_inputs.R")

default_gwas_path <- "data/raw/Testing/34737426-GCST90043906-EFO_0000403.h.tsv"
negative_test_cases <- list(
  list(
    path = "data/raw/Testing/gwas_missing_variant_identifier.tsv",
    label = "missing variant identifier",
    expected_text = "GWAS SNP identifiers"
  ),
  list(
    path = "data/raw/Testing/gwas_missing_chromosome.tsv",
    label = "missing chromosome",
    expected_text = "chromosome"
  ),
  list(
    path = "data/raw/Testing/gwas_missing_base_pair_location.tsv",
    label = "missing base-pair position",
    expected_text = "base-pair position"
  ),
  list(
    path = "data/raw/Testing/gwas_missing_p_value.tsv",
    label = "missing p-value",
    expected_text = "p-value"
  )
)

default_somatic_path <- "data/raw/Testing/2026-01-09_no_CLL_lymph_only_pcawg_maf_tcga_order_hg38.maf"
negative_somatic_test_cases <- list(
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

default_reg_ref_path <- "data/raw/Testing/reg_elements_valid.loc"
default_bw_files <- c(
  "data/raw/Testing/SRR1020514_DLBCL_P265_H3K27ac_ChIPseq.bw",
  "data/raw/Testing/SRR1020516_DLBCL_P286_H3K27ac_ChIPseq.bw",
  "data/raw/Testing/SRR1020518_DLBCL_P397_H3K27ac_ChIPseq.bw"
)
default_test_output_dir <- file.path(tempdir(), "conseguiR_test_outputs", "inputs")

negative_epigenomic_test_cases <- list(
  list(
    mode = "reg_ref",
    path = "data/raw/Testing/reg_elements_missing_end.loc",
    label = "missing regulatory end column",
    expected_text = "at least 4 columns"
  ),
  list(
    mode = "reg_ref",
    path = "data/raw/Testing/reg_elements_bad_coordinates.loc",
    label = "bad regulatory coordinates",
    expected_text = "no usable rows"
  ),
  list(
    mode = "bw_files",
    files = c(
      "data/raw/Testing/SRR1020514_DLBCL_P265_H3K27ac_ChIPseq.bw",
      "data/raw/Testing/SRR1020516_DLBCL_P286_H3K27ac_ChIPseq.bw"
    ),
    label = "fewer than three bigwigs",
    expected_text = "at least three bigWig files"
  ),
  list(
    mode = "bw_files",
    files = c(
      "data/raw/Testing/SRR1020514_DLBCL_P265_H3K27ac_ChIPseq.bw",
      "data/raw/Testing/SRR1020516_DLBCL_P286_H3K27ac_ChIPseq.bw",
      "data/raw/Testing/broken_signal_track.bw"
    ),
    label = "broken bigwig file",
    expected_text = "Could not read bigWig header"
  )
)

test_validate_gwas_sumstats <- function(path = default_gwas_path) {
  message("Reading GWAS file: ", path)
  gwas_raw <- fread(path)

  message("Running validate_gwas_sumstats()")
  gwas_validated <- validate_gwas_sumstats(gwas_raw)
  magma_ready <- prepare_magma_input(gwas_raw)
  dir.create(default_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  snp_loc_outfile <- file.path(default_test_output_dir, "validation_test_snp_loc.tsv")
  pval_outfile <- file.path(default_test_output_dir, "validation_test_pval.tsv")
  write_magma_input_files(gwas_raw, snp_loc_path = snp_loc_outfile, pval_path = pval_outfile)
  snp_loc_written <- fread(snp_loc_outfile, header = FALSE)
  pval_written <- fread(pval_outfile, header = TRUE)

  required_cols <- c("variant_id", "chromosome", "base_pair_location", "p_value")
  missing_cols <- setdiff(required_cols, names(gwas_validated))
  snp_loc_required <- c("SNP", "CHR", "POS")
  pval_required <- c("SNP", "P")

  if (length(missing_cols) > 0) {
    stop("Test failed. Missing required output columns: ", paste(missing_cols, collapse = ", "))
  }
  if (length(setdiff(snp_loc_required, names(magma_ready$snp_loc))) > 0) {
    stop("Test failed. MAGMA snp_loc output is missing required columns.")
  }
  if (length(setdiff(pval_required, names(magma_ready$pval))) > 0) {
    stop("Test failed. MAGMA pval output is missing required columns.")
  }
  if (!identical(names(magma_ready$snp_loc), snp_loc_required)) {
    stop("Test failed. MAGMA snp_loc columns are not in the expected order: SNP, CHR, POS.")
  }
  if (!identical(names(magma_ready$pval), pval_required)) {
    stop("Test failed. MAGMA pval columns are not in the expected order: SNP, P.")
  }
  if (!file.exists(snp_loc_outfile) || !file.exists(pval_outfile)) {
    stop("Test failed. MAGMA output files were not written.")
  }
  if (ncol(snp_loc_written) != 3) {
    stop("Test failed. Written MAGMA snp_loc file does not have exactly 3 columns.")
  }
  if (!identical(names(pval_written), pval_required)) {
    stop("Test failed. Written MAGMA pval file headers are not SNP and P.")
  }
  if (!identical(as.character(magma_ready$snp_loc$SNP[[1]]), as.character(gwas_validated$variant_id[[1]]))) {
    stop("Test failed. MAGMA snp_loc SNP column does not match validated variant_id values.")
  }
  if (!identical(as.character(magma_ready$pval$SNP[[1]]), as.character(gwas_validated$variant_id[[1]]))) {
    stop("Test failed. MAGMA pval SNP column does not match validated variant_id values.")
  }
  if (!isTRUE(all.equal(as.numeric(magma_ready$pval$P[[1]]), as.numeric(gwas_validated$p_value[[1]])))) {
    stop("Test failed. MAGMA pval P column does not match validated p_value values.")
  }

  message("Test passed.")
  message("Rows: ", nrow(gwas_validated))
  message("Columns: ", paste(names(gwas_validated), collapse = ", "))
  print(head(gwas_validated))

  invisible(list(gwas = gwas_validated, magma = magma_ready))
}

test_validate_gwas_sumstats_negative_case <- function(path, label, expected_text) {
  message("Reading negative test file: ", path)
  gwas_raw <- fread(path)

  result <- tryCatch(
    {
      validate_gwas_sumstats(gwas_raw)
      NULL
    },
    error = function(e) e
  )

  if (is.null(result)) {
    stop("Negative test failed for ", label, ": validation unexpectedly succeeded.")
  }

  if (!grepl(expected_text, conditionMessage(result), fixed = TRUE)) {
    stop(
      "Negative test failed for ", label,
      ": expected error text containing '", expected_text,
      "', got '", conditionMessage(result), "'."
    )
  }

  message("Negative test passed for: ", label)
  invisible(result)
}

test_validate_somatic_maf <- function(path = default_somatic_path) {
  message("Reading somatic MAF file: ", path)
  maf_raw <- fread(path)

  message("Running validate_somatic_maf()")
  maf_validated <- validate_somatic_maf(maf_raw)
  dndscv_ready <- prepare_dndscv_input(maf_raw)
  fishhook_ready <- prepare_fishhook_input(maf_raw)

  somatic_required <- c("sample_id", "chromosome", "start_position", "end_position", "ref", "alt")
  dndscv_required <- c("sampleID", "chr", "pos", "ref", "mut")
  fishhook_required <- c("Tumor_Sample_Barcode", "Chromosome", "Start_Position", "End_Position")

  if (length(setdiff(somatic_required, names(maf_validated))) > 0) {
    stop("Somatic validation test failed: validated output is missing canonical columns.")
  }
  if (length(setdiff(dndscv_required, names(dndscv_ready))) > 0) {
    stop("Somatic validation test failed: dndscv output is missing required columns.")
  }
  if (length(setdiff(fishhook_required, names(fishhook_ready))) > 0) {
    stop("Somatic validation test failed: fishhook output is missing required columns.")
  }
  if (!identical(names(dndscv_ready), dndscv_required)) {
    stop("Somatic validation test failed: dndscv output columns are not in the expected order.")
  }
  if (!identical(names(fishhook_ready)[1:4], fishhook_required)) {
    stop("Somatic validation test failed: fishhook output leading columns are not in the expected order.")
  }
  if (!identical(as.character(dndscv_ready$sampleID[[1]]), as.character(maf_validated$sample_id[[1]]))) {
    stop("Somatic validation test failed: dndscv sampleID does not match validated sample_id.")
  }
  if (!identical(as.character(dndscv_ready$chr[[1]]), as.character(maf_validated$chromosome[[1]]))) {
    stop("Somatic validation test failed: dndscv chr does not match validated chromosome.")
  }
  if (!identical(as.integer(dndscv_ready$pos[[1]]), as.integer(maf_validated$start_position[[1]]))) {
    stop("Somatic validation test failed: dndscv pos does not match validated start_position.")
  }
  if (!identical(as.character(dndscv_ready$ref[[1]]), as.character(maf_validated$ref[[1]]))) {
    stop("Somatic validation test failed: dndscv ref does not match validated ref.")
  }
  if (!identical(as.character(dndscv_ready$mut[[1]]), as.character(maf_validated$alt[[1]]))) {
    stop("Somatic validation test failed: dndscv mut does not match validated alt.")
  }
  if (!identical(as.character(fishhook_ready$Tumor_Sample_Barcode[[1]]), as.character(maf_validated$sample_id[[1]]))) {
    stop("Somatic validation test failed: fishhook Tumor_Sample_Barcode does not match validated sample_id.")
  }
  if (!identical(as.character(fishhook_ready$Chromosome[[1]]), as.character(maf_validated$chromosome[[1]]))) {
    stop("Somatic validation test failed: fishhook Chromosome does not match validated chromosome.")
  }
  if (!identical(as.integer(fishhook_ready$Start_Position[[1]]), as.integer(maf_validated$start_position[[1]]))) {
    stop("Somatic validation test failed: fishhook Start_Position does not match validated start_position.")
  }
  if (!identical(as.integer(fishhook_ready$End_Position[[1]]), as.integer(maf_validated$end_position[[1]]))) {
    stop("Somatic validation test failed: fishhook End_Position does not match validated end_position.")
  }

  message("Somatic positive test passed.")
  message("Canonical columns: ", paste(names(maf_validated), collapse = ", "))
  invisible(list(maf = maf_validated, dndscv = dndscv_ready, fishhook = fishhook_ready))
}

test_validate_somatic_maf_negative_case <- function(path, label, expected_text) {
  message("Reading negative somatic test file: ", path)
  maf_raw <- fread(path)

  result <- tryCatch(
    {
      validate_somatic_maf(maf_raw)
      NULL
    },
    error = function(e) e
  )

  if (is.null(result)) {
    stop("Negative somatic test failed for ", label, ": validation unexpectedly succeeded.")
  }

  if (!grepl(expected_text, conditionMessage(result), fixed = TRUE)) {
    stop(
      "Negative somatic test failed for ", label,
      ": expected error text containing '", expected_text,
      "', got '", conditionMessage(result), "'."
    )
  }

  message("Negative somatic test passed for: ", label)
  invisible(result)
}

test_validate_epigenomic_inputs <- function(reg_ref_path = default_reg_ref_path, bw_files = default_bw_files) {
  message("Reading regulatory element reference: ", reg_ref_path)
  reg_gr <- validate_regulatory_element_reference(reg_ref_path)

  message("Running validate_epigenomic_bigwigs()")
  bw_summary <- validate_epigenomic_bigwigs(bw_files = bw_files, reg_gr = reg_gr)
  input_bundle <- validate_epigenomic_inputs(bw_files = bw_files, reg_ref_path = reg_ref_path)

  if (length(reg_gr) == 0) {
    stop("Epigenomic validation test failed: regulatory reference produced zero intervals.")
  }
  if (nrow(bw_summary) < 3) {
    stop("Epigenomic validation test failed: expected validation summary for at least 3 bigWig files.")
  }
  if (!all(c("file", "n_common_seqlevels", "n_test_intervals") %in% names(bw_summary))) {
    stop("Epigenomic validation test failed: bigWig summary is missing expected columns.")
  }
  if (!inherits(input_bundle$reg_gr, "GRanges")) {
    stop("Epigenomic validation test failed: bundled regulatory reference is not a GRanges.")
  }
  if (!is.data.frame(input_bundle$bigwig_summary)) {
    stop("Epigenomic validation test failed: bundled bigWig summary is not tabular.")
  }

  message("Epigenomic positive test passed.")
  message("Validated bigWig files: ", nrow(bw_summary))
  invisible(list(reg_gr = reg_gr, bigwig_summary = bw_summary, bundle = input_bundle))
}

test_validate_epigenomic_negative_case <- function(case, reg_ref_path = default_reg_ref_path, bw_files = default_bw_files) {
  result <- tryCatch(
    {
      if (identical(case$mode, "reg_ref")) {
        validate_regulatory_element_reference(case$path)
      } else if (identical(case$mode, "bw_files")) {
        reg_gr <- validate_regulatory_element_reference(reg_ref_path)
        validate_epigenomic_bigwigs(bw_files = case$files, reg_gr = reg_gr)
      } else {
        stop("Unknown epigenomic negative test mode: ", case$mode)
      }
      NULL
    },
    error = function(e) e
  )

  if (is.null(result)) {
    stop("Negative epigenomic test failed for ", case$label, ": validation unexpectedly succeeded.")
  }

  if (!grepl(case$expected_text, conditionMessage(result), fixed = TRUE)) {
    stop(
      "Negative epigenomic test failed for ", case$label,
      ": expected error text containing '", case$expected_text,
      "', got '", conditionMessage(result), "'."
    )
  }

  message("Negative epigenomic test passed for: ", case$label)
  invisible(result)
}

main <- function() {
  test_validate_gwas_sumstats()

  for (case in negative_test_cases) {
    test_validate_gwas_sumstats_negative_case(
      path = case$path,
      label = case$label,
      expected_text = case$expected_text
    )
  }

  test_validate_somatic_maf()

  for (case in negative_somatic_test_cases) {
    test_validate_somatic_maf_negative_case(
      path = case$path,
      label = case$label,
      expected_text = case$expected_text
    )
  }

  test_validate_epigenomic_inputs()

  for (case in negative_epigenomic_test_cases) {
    test_validate_epigenomic_negative_case(case = case)
  }

  message("All harmonisation and validation tests passed.")
}

if (sys.nframe() == 0) {
  main()
}
