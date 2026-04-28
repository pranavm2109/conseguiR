#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
})

source("R/zzz.R")
source("R/backend_resources.R")
source("R/user_api.R")
source("scripts/Internals/R/03_prepare_germline_scores.R")

.onLoad(libname = ".", pkgname = "conseguiR")

default_gwas_path <- "data/raw/Testing/34737426-GCST90043906-EFO_0000403.h.tsv"
default_gene_loc_path <- "data/raw/NCBI38/NCBI38.gene.loc"
default_reg_loc_path <- if (exists(".conseguiR_default_reg_loc_path", inherits = TRUE)) {
  .conseguiR_default_reg_loc_path()
} else {
  "data/processed/GRCh38-cCREs.loc"
}
default_reference_bfile <- "data/raw/g1000_eur/g1000_eur"
default_sample_size <- 456348L
default_germline_test_output_dir <- "data/processed/test_outputs/germline"
default_run_full_magma_tests <- identical(Sys.getenv("CONSEGUIR_RUN_FULL_MAGMA_TESTS", unset = "0"), "1")

expected_germline_wrapper_formals <- c(
  "gwas_sumstats", "reference_bfile", "sample_size", "sample_size_col",
  "magma_path", "magma_gwas_cache_prefix", "reuse_existing_gwas_cache",
  "reuse_existing_annotation", "reuse_existing_analysis", "keep_intermediates",
  "annotation_window", "filter_path", "ignore_strand", "nonhuman",
  "annotate_modifiers", "step1_general_args",
  "step1_extra_args", "gene_model", "gene_model_modifiers", "genes_only", "pval_use",
  "pval_snp_id", "pval_pval",
  "pval_duplicate", "bfile_synonyms", "bfile_synonym_dup",
  "gene_settings", "batch", "seed", "big_data",
  "step2_general_args",
  "step2_extra_args", "step1_args", "step2_args", "verbose"
)

expected_prepare_germline_formals <- c(
  "gwas_sumstats", "reference_bfile", "gene_output_prefix", "reg_output_prefix",
  "magma_gwas_cache_prefix", "gene_sample_size", "gene_sample_size_col",
  "reg_sample_size", "reg_sample_size_col", "gene_step1_args",
  "gene_step2_args", "reg_step1_args", "reg_step2_args", "shared_args",
  "verbose"
)

make_germline_test_path <- function(stem, ext = "") {
  dir.create(default_germline_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_germline_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

print_magma_coverage_summary <- function() {
  message("MAGMA coverage matrix in this script:")
  message("  Step 1 annotation parameters asserted:")
  message("    annotation_window, filter_path, ignore_strand, nonhuman, annotate_modifiers, step1_general_args, step1_extra_args")
  message("  Step 2 gene-analysis parameters asserted:")
  message("    gene_model, gene_model_modifiers, genes_only, sample_size, sample_size_col, pval_use, pval_snp_id, pval_pval, pval_duplicate, bfile_synonyms, bfile_synonym_dup, gene_settings, batch, seed, big_data, step2_general_args, step2_extra_args")
  message("  Wrapper paths asserted:")
  message("    direct gene wrapper, direct regulatory wrapper, prepare_germline_scores gene_step1_args/gene_step2_args/reg_step1_args/reg_step2_args/shared_args")
}

make_tiny_gwas_fixture <- function() {
  data.table(
    hm_rsid = c("rs1", "rs2", "rs3"),
    hm_variant_id = c("1_100_A_G", "1_200_C_T", "2_300_G_A"),
    hm_chrom = c("1", "1", "2"),
    hm_pos = c(100L, 200L, 300L),
    p_value = c(0.01, 0.2, 0.05),
    N = c(100L, 200L, 300L)
  )
}

make_fake_magma_executable <- function() {
  script_path <- make_germline_test_path("fake_magma", ".sh")
  args_log_path <- make_germline_test_path("fake_magma_args", ".txt")

  writeLines(c(
    "#!/usr/bin/env bash",
    "set -euo pipefail",
    paste0("printf '%s\\n' \"$@\" > ", shQuote(args_log_path)),
    "out=''",
    "annot=0",
    "for ((i=1; i<=$#; i++)); do",
    "  arg=\"${!i}\"",
    "  if [[ \"$arg\" == \"--annotate\" ]]; then annot=1; fi",
    "  if [[ \"$arg\" == \"--out\" ]]; then",
    "    j=$((i+1))",
    "    out=\"${!j}\"",
    "  fi",
    "done",
    "if [[ -z \"$out\" ]]; then",
    "  echo 'Missing --out argument' >&2",
    "  exit 1",
    "fi",
    "if [[ \"$annot\" == \"1\" ]]; then",
    "  printf 'GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs\nTEST\t1\t1\t1000\t2\trs1,rs2\n' > \"${out}.genes.annot\"",
    "else",
    "  printf 'GENE\tGENE_NAME\tZSTAT\tP\nTEST\tTEST\t2.1\t0.03\n' > \"${out}.genes.out\"",
    "fi"
  ), script_path)

  Sys.chmod(script_path, mode = "0755")

  list(
    path = script_path,
    args_log_path = args_log_path
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

test_magma_step1_argument_forwarding <- function() {
  fake <- make_fake_magma_executable()
  filter_path <- make_germline_test_path("magma_filter", ".txt")
  writeLines("TEST", filter_path)
  output_prefix <- make_germline_test_path("magma_step1_args")

  result <- run_magma_step1_annotation(
    gwas_sumstats = make_tiny_gwas_fixture(),
    gene_loc_path = default_gene_loc_path,
    output_prefix = output_prefix,
    magma_path = fake$path,
    annotation_window = c(35, 10),
    filter_path = filter_path,
    ignore_strand = TRUE,
    nonhuman = TRUE,
    extra_args = c("--step1-extra", "foo")
  )

  args <- readLines(fake$args_log_path)
  expect_true(file.exists(result$annot_path))
  expect_true("--annotate" %in% args)
  expect_true("window=35,10" %in% args)
  expect_true(paste0("filter=", filter_path) %in% args)
  expect_true("ignore-strand" %in% args)
  expect_true("nonhuman" %in% args)
  expect_true("--step1-extra" %in% args)
  expect_true("foo" %in% args)
}

test_magma_step2_argument_forwarding <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_step2_args")
  gene_annot_path <- make_germline_test_path("magma_step2_gene_annot", ".genes.annot")
  pval_path <- make_germline_test_path("magma_step2_pval", ".tsv")
  synonyms_path <- paste0(reference_bfile, ".synonyms")

  writeLines(c("GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs", "TEST\t1\t1\t1000\t2\trs1,rs2"), gene_annot_path)
  fwrite(data.table(SNP = c("rs1", "rs2"), P = c(0.01, 0.2)), pval_path, sep = "\t")

  result <- run_magma_step2_gene_analysis(
    gene_annot_path = gene_annot_path,
    pval_path = pval_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = default_sample_size,
    magma_path = fake$path,
    gene_model = "snp-wise=mean",
    genes_only = TRUE,
    pval_use = c("SNP", "P"),
    pval_duplicate = "drop",
    bfile_synonyms = if (file.exists(synonyms_path)) synonyms_path else NULL,
    bfile_synonym_dup = if (file.exists(synonyms_path)) "drop" else NULL,
    extra_args = c("--step2-extra", "bar")
  )

  args <- readLines(fake$args_log_path)
  expect_true(file.exists(result$genes_out_path))
  expect_true("--bfile" %in% args)
  expect_true(reference_bfile %in% args)
  expect_true("--gene-annot" %in% args)
  expect_true(gene_annot_path %in% args)
  expect_true("--pval" %in% args)
  expect_true(pval_path %in% args)
  expect_true("use=SNP,P" %in% args)
  expect_true(paste0("N=", default_sample_size) %in% args)
  expect_true("duplicate=drop" %in% args)
  if (file.exists(synonyms_path)) {
    expect_true(paste0("synonyms=", synonyms_path) %in% args)
    expect_true("synonym-dup=drop" %in% args)
  }
  expect_true("--gene-model" %in% args)
  expect_true("snp-wise=mean" %in% args)
  expect_true("--genes-only" %in% args)
  expect_true("--step2-extra" %in% args)
  expect_true("bar" %in% args)
}

test_magma_step2_sample_size_col_forwarding <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_step2_ncol_args")
  gene_annot_path <- make_germline_test_path("magma_step2_ncol_gene_annot", ".genes.annot")
  pval_path <- make_germline_test_path("magma_step2_ncol_pval", ".tsv")

  writeLines(c("GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs", "TEST\t1\t1\t1000\t2\trs1,rs2"), gene_annot_path)
  fwrite(data.table(SNP = c("rs1", "rs2"), P = c(0.01, 0.2), N = c(100L, 200L)), pval_path, sep = "\t")

  result <- run_magma_step2_gene_analysis(
    gene_annot_path = gene_annot_path,
    pval_path = pval_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size_col = "N",
    magma_path = fake$path,
    genes_only = FALSE
  )

  args <- readLines(fake$args_log_path)
  expect_true(file.exists(result$genes_out_path))
  expect_true("ncol=N" %in% args)
  expect_false("--genes-only" %in% args)
}

test_magma_step1_general_arg_forwarding <- function() {
  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_step1_general_args")

  result <- run_magma_step1_annotation(
    gwas_sumstats = make_tiny_gwas_fixture(),
    gene_loc_path = default_gene_loc_path,
    output_prefix = output_prefix,
    magma_path = fake$path,
    step1_args = list(
      annotation_window = c(50, 20),
      annotate_modifiers = c("custom-mod"),
      general_args = list(
        test_flag = TRUE,
        test_value = "alpha"
      ),
      extra_args = c("--tail-flag", "omega")
    )
  )

  args <- readLines(fake$args_log_path)
  expect_true(file.exists(result$annot_path))
  expect_true("window=50,20" %in% args)
  expect_true("custom-mod" %in% args)
  expect_true("--test-flag" %in% args)
  expect_true("--test-value" %in% args)
  expect_true("alpha" %in% args)
  expect_true("--tail-flag" %in% args)
  expect_true("omega" %in% args)
}

test_magma_step2_general_arg_forwarding <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_step2_general_args")
  gene_annot_path <- make_germline_test_path("magma_step2_general_gene_annot", ".genes.annot")
  pval_path <- make_germline_test_path("magma_step2_general_pval", ".tsv")

  writeLines(c("GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs", "TEST\t1\t1\t1000\t2\trs1,rs2"), gene_annot_path)
  fwrite(data.table(SNP = c("rs1", "rs2"), P = c(0.01, 0.2), N = c(100L, 200L)), pval_path, sep = "\t")

  result <- run_magma_step2_gene_analysis(
    gene_annot_path = gene_annot_path,
    pval_path = pval_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    magma_path = fake$path,
    step2_args = list(
      pval = list(
        use = c("SNP", "P"),
        ncol = "N",
        duplicate = "drop"
      ),
      bfile = list(
        synonyms = "synonyms.tsv",
        synonym_dup = "drop"
      ),
      gene_model = "snp-wise=mean",
      genes_only = FALSE,
      general_args = list(
        test_flag = TRUE,
        test_value = "beta"
      ),
      extra_args = c("--step2-tail", "gamma")
    )
  )

  args <- readLines(fake$args_log_path)
  expect_true(file.exists(result$genes_out_path))
  expect_true("use=SNP,P" %in% args)
  expect_true("ncol=N" %in% args)
  expect_true("duplicate=drop" %in% args)
  expect_true("synonyms=synonyms.tsv" %in% args)
  expect_true("synonym-dup=drop" %in% args)
  expect_true("--gene-model" %in% args)
  expect_true("snp-wise=mean" %in% args)
  expect_false("--genes-only" %in% args)
  expect_true("--test-flag" %in% args)
  expect_true("--test-value" %in% args)
  expect_true("beta" %in% args)
  expect_true("--step2-tail" %in% args)
  expect_true("gamma" %in% args)
}

test_magma_step2_explicit_pval_modifier_forwarding <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_step2_explicit_pval_args")
  gene_annot_path <- make_germline_test_path("magma_step2_explicit_pval_gene_annot", ".genes.annot")
  pval_path <- make_germline_test_path("magma_step2_explicit_pval", ".tsv")

  writeLines(c("GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs", "TEST\t1\t1\t1000\t2\trs1,rs2"), gene_annot_path)
  fwrite(data.table(rsid = c("rs1", "rs2"), pv = c(0.01, 0.2), NX = c(100L, 200L)), pval_path, sep = "\t")

  result <- run_magma_step2_gene_analysis(
    gene_annot_path = gene_annot_path,
    pval_path = pval_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    magma_path = fake$path,
    pval_snp_id = "rsid",
    pval_pval = "pv",
    sample_size = c(1000L, 1100L, 1200L),
    genes_only = FALSE
  )

  args <- readLines(fake$args_log_path)
  expect_true(file.exists(result$genes_out_path))
  expect_true("snp-id=rsid" %in% args)
  expect_true("pval=pv" %in% args)
  expect_true("N=1000,1100,1200" %in% args)
  expect_false("--genes-only" %in% args)
}

test_magma_step2_gene_settings_and_seed_forwarding <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_step2_gene_settings_args")
  gene_annot_path <- make_germline_test_path("magma_step2_gene_settings_gene_annot", ".genes.annot")
  pval_path <- make_germline_test_path("magma_step2_gene_settings_pval", ".tsv")
  snp_include_path <- make_germline_test_path("magma_gene_settings_snp_include", ".txt")
  snp_exclude_path <- make_germline_test_path("magma_gene_settings_snp_exclude", ".txt")
  indiv_include_path <- make_germline_test_path("magma_gene_settings_indiv_include", ".txt")
  indiv_exclude_path <- make_germline_test_path("magma_gene_settings_indiv_exclude", ".txt")

  writeLines(c("GENE\tCHR\tSTART\tSTOP\tNSNPS\tSNPs", "TEST\t1\t1\t1000\t2\trs1,rs2"), gene_annot_path)
  fwrite(data.table(SNP = c("rs1", "rs2"), P = c(0.01, 0.2), N = c(100L, 200L)), pval_path, sep = "\t")
  writeLines(c("rs1", "rs2"), snp_include_path)
  writeLines("rs3", snp_exclude_path)
  writeLines("F1 I1", indiv_include_path)
  writeLines("F2 I2", indiv_exclude_path)

  result <- run_magma_step2_gene_analysis(
    gene_annot_path = gene_annot_path,
    pval_path = pval_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    magma_path = fake$path,
    sample_size_col = "N",
    gene_model = "multi",
    gene_model_modifiers = "multi-show-all",
    gene_settings = list(
      snp_min_maf = 0.01,
      snp_min_mac = 2,
      snp_max_maf = 0.49,
      snp_max_mac = 100,
      snp_max_miss = 0.1,
      snp_diff = 1e-06,
      snp_include = snp_include_path,
      snp_exclude = snp_exclude_path,
      indiv_include = indiv_include_path,
      indiv_exclude = indiv_exclude_path,
      prune = 0.95,
      prune_prop = 0.5,
      prune_count = 20,
      fixed_permp = 5000,
      adap_permp = c(1e6, 25),
      min_perm = 1000
    ),
    batch = c(7, 20),
    seed = 12345,
    big_data = FALSE,
    genes_only = FALSE
  )

  args <- readLines(fake$args_log_path)
  expect_true(file.exists(result$genes_out_path))
  expect_true("--gene-model" %in% args)
  expect_true("multi" %in% args)
  expect_true("multi-show-all" %in% args)
  expect_true("--gene-settings" %in% args)
  expect_true("snp-min-maf=0.01" %in% args)
  expect_true("snp-min-mac=2" %in% args)
  expect_true("snp-max-maf=0.49" %in% args)
  expect_true("snp-max-mac=100" %in% args)
  expect_true("snp-max-miss=0.1" %in% args)
  expect_true("snp-diff=1e-06" %in% args || "snp-diff=0.000001" %in% args)
  expect_true(paste0("snp-include=", snp_include_path) %in% args)
  expect_true(paste0("snp-exclude=", snp_exclude_path) %in% args)
  expect_true(paste0("indiv-include=", indiv_include_path) %in% args)
  expect_true(paste0("indiv-exclude=", indiv_exclude_path) %in% args)
  expect_true("prune=0.95" %in% args)
  expect_true("prune-prop=0.5" %in% args)
  expect_true("prune-count=20" %in% args)
  expect_true("fixed-permp=5000" %in% args)
  expect_true("adap-permp=1e+06,25" %in% args || "adap-permp=1000000,25" %in% args)
  expect_true("min-perm=1000" %in% args)
  expect_true("--batch" %in% args)
  expect_true("7" %in% args)
  expect_true("20" %in% args)
  expect_true("--seed" %in% args)
  expect_true("12345" %in% args)
  expect_true("--big-data" %in% args)
  expect_true("off" %in% args)
}

test_prepare_germline_scores_forwards_all_stage_lists <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  gwas_dt <- make_tiny_gwas_fixture()
  output_prefix <- make_germline_test_path("magma_prepare_germline")

  result <- prepare_germline_scores(
    gwas_sumstats = gwas_dt,
    reference_bfile = reference_bfile,
    gene_sample_size = default_sample_size,
    reg_sample_size_col = "N",
    gene_step1_args = list(
      annotation_window = c(35, 10),
      ignore_strand = TRUE,
      nonhuman = TRUE,
      annotate_modifiers = "gene-mod",
      general_args = list(
        gene_step1_flag = TRUE
      ),
      extra_args = c("--gene-step1-extra", "foo")
    ),
    gene_step2_args = list(
      gene_model = "snp-wise=mean",
      genes_only = TRUE,
      pval = list(
        use = c("SNP", "P"),
        duplicate = "drop"
      ),
      general_args = list(
        gene_step2_flag = TRUE
      ),
      extra_args = c("--gene-step2-extra", "bar")
    ),
    reg_step1_args = list(
      annotation_window = c(0, 0),
      ignore_strand = FALSE,
      nonhuman = FALSE,
      general_args = list(
        reg_step1_flag = TRUE
      ),
      extra_args = c("--reg-step1-extra", "alpha")
    ),
    reg_step2_args = list(
      genes_only = FALSE,
      pval = list(
        snp_id = "SNP",
        pval = "P",
        ncol = "N"
      ),
      general_args = list(
        reg_step2_flag = TRUE
      ),
      extra_args = c("--reg-step2-extra", "beta")
    ),
    shared_args = list(
      magma_path = fake$path
    ),
    verbose = FALSE
  )

  expect_s3_class(result, "conseguiR_bundle")
  expect_true(is.data.frame(result$gene_scores))
  expect_true(is.data.frame(result$reg_scores))
  expect_s3_class(result$gene_result, "conseguiR_bundle")
  expect_s3_class(result$reg_result, "conseguiR_bundle")

  gene_step1 <- result$gene_result$pipeline$step1$args
  gene_step2 <- result$gene_result$pipeline$step2$args
  reg_step1 <- result$reg_result$pipeline$step1$args
  reg_step2 <- result$reg_result$pipeline$step2$args

  expect_true("window=35,10" %in% gene_step1)
  expect_true("ignore-strand" %in% gene_step1)
  expect_true("nonhuman" %in% gene_step1)
  expect_true("gene-mod" %in% gene_step1)
  expect_true("--gene-step1-flag" %in% gene_step1)
  expect_true("--gene-step1-extra" %in% gene_step1)
  expect_true("foo" %in% gene_step1)

  expect_true("--gene-model" %in% gene_step2)
  expect_true("snp-wise=mean" %in% gene_step2)
  expect_true("use=SNP,P" %in% gene_step2)
  expect_true(paste0("N=", default_sample_size) %in% gene_step2)
  expect_true("duplicate=drop" %in% gene_step2)
  expect_true("--genes-only" %in% gene_step2)
  expect_true("--gene-step2-flag" %in% gene_step2)
  expect_true("--gene-step2-extra" %in% gene_step2)
  expect_true("bar" %in% gene_step2)

  expect_true("window=0,0" %in% reg_step1)
  expect_true("--reg-step1-flag" %in% reg_step1)
  expect_true("--reg-step1-extra" %in% reg_step1)
  expect_true("alpha" %in% reg_step1)

  expect_true("ncol=N" %in% reg_step2)
  expect_true("snp-id=SNP" %in% reg_step2)
  expect_true("pval=P" %in% reg_step2)
  expect_false("--genes-only" %in% reg_step2)
  expect_true("--reg-step2-flag" %in% reg_step2)
  expect_true("--reg-step2-extra" %in% reg_step2)
  expect_true("beta" %in% reg_step2)
}

test_germline_wrapper_surface_matches_supported_api <- function() {
  gene_formals <- names(formals(run_germline_gene_scoring))
  reg_formals <- names(formals(run_germline_regulatory_scoring))
  prepare_formals <- names(formals(prepare_germline_scores))

  expect_true(all(expected_germline_wrapper_formals %in% gene_formals))
  expect_true(all(expected_germline_wrapper_formals %in% reg_formals))
  expect_true(all(expected_prepare_germline_formals %in% prepare_formals))
}

test_germline_wrapper_respects_output_and_cache_controls <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_wrapper_output_cache")
  cache_prefix <- make_germline_test_path("magma_wrapper_cache")

  first <- run_germline_gene_scoring(
    gwas_sumstats = make_tiny_gwas_fixture(),
    gene_loc_path = default_gene_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = default_sample_size,
    magma_path = fake$path,
    magma_gwas_cache_prefix = cache_prefix,
    keep_intermediates = TRUE,
    reuse_existing_gwas_cache = FALSE,
    reuse_existing_annotation = FALSE,
    reuse_existing_analysis = FALSE
  )

  expect_true(file.exists(paste0(output_prefix, ".zstat.tsv")))
  expect_true(file.exists(paste0(cache_prefix, ".snp_loc.tsv")))
  expect_true(file.exists(paste0(cache_prefix, ".pval.tsv")))

  second <- run_germline_gene_scoring(
    gwas_sumstats = make_tiny_gwas_fixture(),
    gene_loc_path = default_gene_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = default_sample_size,
    magma_path = fake$path,
    magma_gwas_cache_prefix = cache_prefix,
    keep_intermediates = TRUE,
    reuse_existing_gwas_cache = TRUE,
    reuse_existing_annotation = TRUE,
    reuse_existing_analysis = TRUE
  )

  expect_true(isTRUE(second$pipeline$step1$reused_existing_annotation))
  expect_true(isTRUE(second$pipeline$step2$reused_existing_analysis))
  expect_true(file.exists(first$output_paths$gene_scores_path))
}

test_germline_wrapper_forwards_explicit_step_args <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_wrapper_gene")
  filter_path <- make_germline_test_path("magma_wrapper_gene_filter", ".txt")
  writeLines(c("rs1", "rs2"), filter_path)
  result <- run_germline_gene_scoring(
    gwas_sumstats = make_tiny_gwas_fixture(),
    gene_loc_path = default_gene_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = default_sample_size,
    magma_path = fake$path,
    annotation_window = c(35, 10),
    filter_path = filter_path,
    ignore_strand = TRUE,
    nonhuman = TRUE,
    annotate_modifiers = "wrapper-gene-mod",
    step1_general_args = list(wrapper_step1_flag = TRUE),
    step1_extra_args = c("--wrapper-step1-extra", "foo"),
    gene_model = "snp-wise=mean",
    gene_model_modifiers = "multi-show-all",
    genes_only = TRUE,
    pval_use = c("SNP", "P"),
    pval_snp_id = "SNP",
    pval_pval = "P",
    pval_duplicate = "drop",
    gene_settings = list(snp_min_maf = 0.01, snp_max_miss = 0.1),
    batch = c(2, "chr"),
    seed = 123,
    big_data = TRUE,
    step2_general_args = list(wrapper_step2_flag = TRUE),
    step2_extra_args = c("--wrapper-step2-extra", "bar")
  )

  expect_true("window=35,10" %in% result$pipeline$step1$args)
  expect_true(paste0("filter=", filter_path) %in% result$pipeline$step1$args)
  expect_true("ignore-strand" %in% result$pipeline$step1$args)
  expect_true("nonhuman" %in% result$pipeline$step1$args)
  expect_true("wrapper-gene-mod" %in% result$pipeline$step1$args)
  expect_true("--wrapper-step1-flag" %in% result$pipeline$step1$args)
  expect_true("--wrapper-step1-extra" %in% result$pipeline$step1$args)
  expect_true("foo" %in% result$pipeline$step1$args)
  expect_true("--gene-model" %in% result$pipeline$step2$args)
  expect_true("snp-wise=mean" %in% result$pipeline$step2$args)
  expect_true("use=SNP,P" %in% result$pipeline$step2$args)
  expect_true("snp-id=SNP" %in% result$pipeline$step2$args)
  expect_true("pval=P" %in% result$pipeline$step2$args)
  expect_true("duplicate=drop" %in% result$pipeline$step2$args)
  expect_true("multi-show-all" %in% result$pipeline$step2$args)
  expect_true("--gene-settings" %in% result$pipeline$step2$args)
  expect_true("snp-min-maf=0.01" %in% result$pipeline$step2$args)
  expect_true("snp-max-miss=0.1" %in% result$pipeline$step2$args)
  expect_true("--batch" %in% result$pipeline$step2$args)
  expect_true("2" %in% result$pipeline$step2$args)
  expect_true("chr" %in% result$pipeline$step2$args)
  expect_true("--seed" %in% result$pipeline$step2$args)
  expect_true("123" %in% result$pipeline$step2$args)
  expect_true("--big-data" %in% result$pipeline$step2$args)
  expect_true("--genes-only" %in% result$pipeline$step2$args)
  expect_true("--wrapper-step2-flag" %in% result$pipeline$step2$args)
  expect_true("--wrapper-step2-extra" %in% result$pipeline$step2$args)
  expect_true("bar" %in% result$pipeline$step2$args)
}

test_regulatory_wrapper_forwards_explicit_step_args <- function(
  reference_bfile = find_reference_bfile()
) {
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }

  fake <- make_fake_magma_executable()
  output_prefix <- make_germline_test_path("magma_wrapper_reg")
  filter_path <- make_germline_test_path("magma_wrapper_reg_filter", ".txt")
  writeLines(c("rs1", "rs2"), filter_path)
  result <- run_germline_regulatory_scoring(
    gwas_sumstats = make_tiny_gwas_fixture(),
    reg_loc_path = default_reg_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size_col = "N",
    magma_path = fake$path,
    annotation_window = c(0, 0),
    filter_path = filter_path,
    ignore_strand = TRUE,
    nonhuman = TRUE,
    annotate_modifiers = "wrapper-reg-mod",
    step1_general_args = list(wrapper_reg_step1_flag = TRUE),
    step1_extra_args = c("--wrapper-reg-step1-extra", "alpha"),
    gene_model = "multi",
    gene_model_modifiers = "multi-show-all",
    genes_only = FALSE,
    pval_use = c("SNP", "P"),
    pval_snp_id = "SNP",
    pval_pval = "P",
    pval_duplicate = "drop",
    bfile_synonyms = "reg.synonyms",
    bfile_synonym_dup = "drop",
    gene_settings = list(
      snp_min_maf = 0.01,
      snp_max_miss = 0.1
    ),
    batch = c("X", "chr"),
    seed = 321,
    big_data = FALSE,
    step2_general_args = list(wrapper_reg_step2_flag = TRUE),
    step2_extra_args = c("--wrapper-reg-step2-extra", "beta")
  )

  expect_true("window=0,0" %in% result$pipeline$step1$args)
  expect_true(paste0("filter=", filter_path) %in% result$pipeline$step1$args)
  expect_true("ignore-strand" %in% result$pipeline$step1$args)
  expect_true("nonhuman" %in% result$pipeline$step1$args)
  expect_true("wrapper-reg-mod" %in% result$pipeline$step1$args)
  expect_true("--wrapper-reg-step1-flag" %in% result$pipeline$step1$args)
  expect_true("--wrapper-reg-step1-extra" %in% result$pipeline$step1$args)
  expect_true("alpha" %in% result$pipeline$step1$args)
  expect_true("--gene-model" %in% result$pipeline$step2$args)
  expect_true("multi" %in% result$pipeline$step2$args)
  expect_true("multi-show-all" %in% result$pipeline$step2$args)
  expect_true("use=SNP,P" %in% result$pipeline$step2$args)
  expect_true("ncol=N" %in% result$pipeline$step2$args)
  expect_true("snp-id=SNP" %in% result$pipeline$step2$args)
  expect_true("pval=P" %in% result$pipeline$step2$args)
  expect_true("duplicate=drop" %in% result$pipeline$step2$args)
  expect_true("synonyms=reg.synonyms" %in% result$pipeline$step2$args)
  expect_true("synonym-dup=drop" %in% result$pipeline$step2$args)
  expect_true("--gene-settings" %in% result$pipeline$step2$args)
  expect_true("snp-min-maf=0.01" %in% result$pipeline$step2$args)
  expect_true("snp-max-miss=0.1" %in% result$pipeline$step2$args)
  expect_true("--batch" %in% result$pipeline$step2$args)
  expect_true("X" %in% result$pipeline$step2$args)
  expect_true("chr" %in% result$pipeline$step2$args)
  expect_true("--seed" %in% result$pipeline$step2$args)
  expect_true("321" %in% result$pipeline$step2$args)
  expect_true("--big-data" %in% result$pipeline$step2$args)
  expect_true("off" %in% result$pipeline$step2$args)
  expect_false("--genes-only" %in% result$pipeline$step2$args)
  expect_true("--wrapper-reg-step2-flag" %in% result$pipeline$step2$args)
  expect_true("--wrapper-reg-step2-extra" %in% result$pipeline$step2$args)
  expect_true("beta" %in% result$pipeline$step2$args)
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
  print_magma_coverage_summary()

  test_that("MAGMA germline smoke tests pass", {
    expect_no_error(run_all_germline_tests(print_scores = FALSE))
  })

  test_that("MAGMA step 1 forwards supported wrapper arguments", {
    test_magma_step1_argument_forwarding()
  })

  test_that("MAGMA step 2 forwards supported wrapper arguments", {
    test_magma_step2_argument_forwarding()
  })

  test_that("MAGMA step 2 forwards sample_size_col and genes_only controls", {
    test_magma_step2_sample_size_col_forwarding()
  })

  test_that("MAGMA step 2 forwards explicit snp-id, pval, and multi-value N modifiers", {
    test_magma_step2_explicit_pval_modifier_forwarding()
  })

  test_that("MAGMA step 2 forwards gene-settings, seed, big-data, and gene-model modifiers", {
    test_magma_step2_gene_settings_and_seed_forwarding()
  })

  test_that("MAGMA step 1 forwards generalized step-list arguments", {
    test_magma_step1_general_arg_forwarding()
  })

  test_that("MAGMA step 2 forwards generalized step-list arguments", {
    test_magma_step2_general_arg_forwarding()
  })

  test_that("germline gene wrapper forwards explicit step arguments", {
    test_germline_wrapper_forwards_explicit_step_args()
  })

  test_that("germline regulatory wrapper forwards explicit step arguments", {
    test_regulatory_wrapper_forwards_explicit_step_args()
  })

  test_that("prepare_germline_scores forwards all four MAGMA stage bundles", {
    test_prepare_germline_scores_forwards_all_stage_lists()
  })

  test_that("germline wrapper surface exposes the supported MAGMA API", {
    test_germline_wrapper_surface_matches_supported_api()
  })

  test_that("germline wrapper respects output and cache controls", {
    test_germline_wrapper_respects_output_and_cache_controls()
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
