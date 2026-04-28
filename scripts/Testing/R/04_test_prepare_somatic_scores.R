#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
  library(GenomicRanges)
})

source("R/zzz.R")
source("R/backend_resources.R")
source("R/user_api.R")
source("scripts/Internals/R/04_prepare_somatic_scores.R")

.onLoad(libname = ".", pkgname = "conseguiR")

expected_dndscv_formals <- c(
  "maf", "refdb", "output_path", "cv", "max_muts_per_gene_per_sample",
  "max_coding_muts_per_sample", "gene_list", "sm", "kc", "use_indel_sites",
  "min_indels", "maxcovs", "constrain_wnon_wspl", "outp", "numcode",
  "outmats", "mingenecovs", "onesided", "dc", "dndscv_args", "verbose"
)

expected_fishhook_constructor_formals <- c(
  "maf", "reg_ref_path", "output_path", "eligible_gr", "fishhook_covariates",
  "fishhook_covariate_data", "idcol", "constructor_out_path",
  "constructor_use_local_mut_density", "constructor_local_mut_density_bin",
  "constructor_mc_cores", "constructor_na_rm", "constructor_pad",
  "constructor_max_slice", "constructor_ff_chunk", "constructor_max_chunk",
  "constructor_idcap", "constructor_weight_events", "constructor_nb"
)

expected_fishhook_score_formals <- c(
  "score_sets", "score_model", "score_return_model", "score_nb", "score_iter",
  "score_subsample", "score_seed", "score_verbose", "score_mc_cores",
  "score_p_randomized", "score_class_return", "fishhook_args", "verbose"
)

expected_prepare_somatic_formals <- c(
  "maf", "refdb", "reg_ref_path", "gene_output_path", "reg_output_path",
  "gene_cv", "gene_max_muts_per_gene_per_sample",
  "gene_max_coding_muts_per_sample", "gene_list", "sm", "kc",
  "use_indel_sites", "min_indels", "maxcovs", "constrain_wnon_wspl",
  "outp", "numcode", "outmats", "mingenecovs", "onesided", "dc",
  "dndscv_args", "eligible_gr", "fishhook_covariates",
  "fishhook_covariate_data", "fishhook_idcol", "constructor_out_path",
  "constructor_use_local_mut_density", "constructor_local_mut_density_bin",
  "constructor_mc_cores", "constructor_na_rm", "constructor_pad",
  "constructor_max_slice", "constructor_ff_chunk", "constructor_max_chunk",
  "constructor_idcap", "constructor_weight_events", "constructor_nb",
  "score_sets", "score_model", "score_return_model", "score_nb",
  "score_iter", "score_subsample", "score_seed", "score_verbose",
  "score_mc_cores", "score_p_randomized", "score_class_return",
  "fishhook_args", "verbose"
)

default_somatic_path <- "data/raw/Testing/2026-01-09_no_CLL_lymph_only_pcawg_maf_tcga_order_hg38.maf"
default_reg_ref_path <- .conseguiR_default_reg_loc_path()
default_fishhook_covariate_path <- "data/raw/Testing/2026-01-26_all_reg_elems_sample_level_mut_frac_comparison_bet_only_memory_b_normal_and_non_cll_malig_b_cells.rds"
default_dndscv_refdb <- "data/raw/Testing/RefCDS_human_GRCh38.p12.rda"
default_somatic_test_output_dir <- "data/processed/test_outputs/somatic"

make_somatic_test_path <- function(stem, ext = "") {
  dir.create(default_somatic_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_somatic_test_output_dir,
    paste0(
      stem, "_",
      format(Sys.time(), "%Y%m%d%H%M%S"), "_",
      sprintf("%06d", sample.int(999999L, 1L)),
      ext
    )
  )
}

print_somatic_coverage_summary <- function() {
  message("Somatic coverage matrix in this script:")
  message("  dndscv parameters asserted:")
  message("    refdb, gene_list, sm, kc, cv, max_muts_per_gene_per_sample, max_coding_muts_per_sample, use_indel_sites, min_indels, maxcovs, constrain_wnon_wspl, outp, numcode, outmats, mingenecovs, onesided, dc, dndscv_args, verbose")
  message("  fishHook constructor parameters asserted:")
  message("    eligible_gr, fishhook_covariates, fishhook_covariate_data, idcol, constructor_out_path, constructor_use_local_mut_density, constructor_local_mut_density_bin, constructor_mc_cores, constructor_na_rm, constructor_pad, constructor_max_slice, constructor_ff_chunk, constructor_max_chunk, constructor_idcap, constructor_weight_events, constructor_nb, verbose")
  message("  fishHook score parameters asserted:")
  message("    score_sets, score_model, score_return_model, score_nb, score_iter, score_subsample, score_seed, score_verbose, score_mc_cores, score_p_randomized, score_class_return, fishhook_args, verbose")
  message("  Wrapper paths asserted:")
  message("    run_somatic_gene_scoring, run_somatic_regulatory_scoring, prepare_somatic_scores")
}

fishhook_covariate_overlap_count <- function(reg_ref_path = default_reg_ref_path,
                                             covariate_path = default_fishhook_covariate_path) {
  if (!file.exists(reg_ref_path) || !file.exists(covariate_path)) {
    return(0L)
  }

  reg_gr <- validate_regulatory_element_reference(reg_ref_path)
  covariate_data <- readRDS(covariate_path)
  cov_dt <- data.table(covariate_data)

  if (!"reg_elem_id" %in% names(S4Vectors::mcols(reg_gr)) || !"reg_elem_id" %in% names(cov_dt)) {
    return(0L)
  }

  reg_ids <- as.character(S4Vectors::mcols(reg_gr)$reg_elem_id)
  cov_ids <- as.character(cov_dt$reg_elem_id)
  sum(reg_ids %in% cov_ids)
}

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

with_traced_namespace_function <- function(fun_name, pkg, tracer_expr, code) {
  trace(
    fun_name,
    where = asNamespace(pkg),
    tracer = tracer_expr,
    print = FALSE
  )
  on.exit(untrace(fun_name, where = asNamespace(pkg)), add = TRUE)
  force(code)
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
    regexp = "`refdb` is required for dndscv scoring",
    fixed = TRUE
  )
}

test_run_dndscv_gene_scoring_negative_bad_refdb <- function() {
  skip_if_not_installed("dndscv")

  maf <- make_dummy_maf()

  expect_error(
    run_dndscv_gene_scoring(maf = maf, refdb = "data/raw/Testing/not_a_real_refdb.rda"),
    regexp = "dndscv `refdb` file does not exist",
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
    regexp = "fishHook covariate data is missing required columns",
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
    regexp = "Regulatory element reference file does not exist",
    fixed = TRUE
  )
}

test_run_dndscv_gene_scoring_forwards_supported_args <- function() {
  skip_if_not_installed("dndscv")

  refdb <- find_dndscv_refdb()
  if (is.null(refdb)) {
    skip("No dndscv refdb found in the repository.")
  }

  captured <- new.env(parent = emptyenv())
  cv_mat <- matrix(
    c(1, 2, 3, 4),
    nrow = 2L,
    dimnames = list(c("TP53", "KRAS"), c("cov1", "cov2"))
  )

  expect_error(
    with_traced_namespace_function(
      "dndscv",
      "dndscv",
      bquote({
        assign("args", as.list(environment()), envir = .(captured))
        stop("captured dndscv call")
      }),
      run_somatic_gene_scoring(
        maf = make_dummy_maf(),
        refdb = refdb,
        max_muts_per_gene_per_sample = 9L,
        max_coding_muts_per_sample = 1234L,
        gene_list = c("TP53", "KRAS"),
        sm = "192r_3w",
        kc = "cgc81",
        cv = cv_mat,
        use_indel_sites = FALSE,
        min_indels = 1L,
        maxcovs = 7L,
        constrain_wnon_wspl = FALSE,
        outp = 1L,
        numcode = 1L,
        outmats = TRUE,
        mingenecovs = 100L,
        onesided = TRUE,
        dc = c(TP53 = 1, KRAS = 2),
        dndscv_args = list()
      )
    ),
    regexp = "captured dndscv call",
    fixed = TRUE
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_true(is.data.frame(args$mutations))
  expect_identical(names(args$mutations), c("sampleID", "chr", "pos", "ref", "mut"))
  expect_identical(args$refdb, refdb)
  expect_true(is.matrix(args$cv))
  expect_identical(dim(args$cv), c(2L, 2L))
  expect_identical(colnames(args$cv), c("cov1", "cov2"))
  expect_equal(args$max_muts_per_gene_per_sample, 9L)
  expect_equal(args$max_coding_muts_per_sample, 1234L)
  expect_identical(args$sm, "192r_3w")
  expect_identical(args$kc, "cgc81")
  expect_identical(args$use_indel_sites, FALSE)
  expect_equal(args$min_indels, 1L)
  expect_equal(args$maxcovs, 7L)
  expect_identical(args$constrain_wnon_wspl, FALSE)
  expect_equal(args$outp, 1L)
  expect_equal(args$numcode, 1L)
  expect_identical(args$outmats, TRUE)
  expect_equal(args$mingenecovs, 100L)
  expect_identical(args$onesided, TRUE)
  expect_equal(as.character(args$gene_list), c("TP53", "KRAS"))
  expect_true(is.numeric(args$dc))
  expect_equal(unname(args$dc), c(1, 2))
}

test_run_fishhook_reg_scoring_forwards_constructor_args <- function() {
  skip_if_not_installed("fishHook")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

  captured <- new.env(parent = emptyenv())
  eligible_gr <- make_fishhook_eligible_hg38()
  cov_gr_full <- make_fishhook_hypothesis_granges("data/raw/Testing/reg_elements_valid.loc")
  cov_gr <- cov_gr_full[seq_len(min(10L, length(cov_gr_full)))]
  cov_obj <- fishHook::Cov(data = cov_gr, name = "dummy_covariate")

  expect_error(
    with_traced_namespace_function(
      "Fish",
      "fishHook",
      bquote({
        assign("args", as.list(environment()), envir = .(captured))
        stop("captured Fish constructor")
      }),
      run_somatic_regulatory_scoring(
        maf = make_dummy_maf(),
        reg_ref_path = "data/raw/Testing/reg_elements_valid.loc",
        eligible_gr = eligible_gr,
        fishhook_covariates = list(cov_obj),
        idcol = "Tumor_Sample_Barcode",
        constructor_out_path = tempdir(),
        constructor_use_local_mut_density = TRUE,
        constructor_local_mut_density_bin = 500000,
        constructor_mc_cores = 3L,
        constructor_na_rm = FALSE,
        constructor_pad = 25,
        constructor_max_slice = 5000L,
        constructor_ff_chunk = 25000L,
        constructor_max_chunk = 500000L,
        constructor_idcap = 2,
        constructor_weight_events = TRUE,
        constructor_nb = FALSE,
        verbose = TRUE
      )
    ),
    regexp = "captured Fish constructor",
    fixed = TRUE
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_s4_class(args$hypotheses, "GRanges")
  expect_s4_class(args$events, "GRanges")
  expect_s4_class(args$eligible, "GRanges")
  expect_identical(args$idcol, "Tumor_Sample_Barcode")
  expect_true(is.list(args$covariates))
  expect_length(args$covariates, 1L)
  expect_identical(args$out.path, tempdir())
  expect_identical(args$use_local_mut_density, TRUE)
  expect_equal(args$local_mut_density_bin, 500000)
  expect_equal(args$mc.cores, 3L)
  expect_identical(args$na.rm, FALSE)
  expect_equal(args$pad, 25)
  expect_identical(args$verbose, TRUE)
  expect_equal(args$max.slice, 5000L)
  expect_equal(args$ff.chunk, 25000L)
  expect_equal(args$max.chunk, 500000L)
  expect_equal(args$idcap, 2)
  expect_identical(args$weightEvents, TRUE)
  expect_identical(args$nb, FALSE)
}

test_run_fishhook_reg_scoring_forwards_score_args <- function() {
  skip_if_not_installed("fishHook")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

  captured <- new.env(parent = emptyenv())
  score_sets <- list(dummy_set = c("EH38E0080197", "EH38E2084302"))

  expect_error(
    with_traced_namespace_function(
      "score.hypotheses",
      "fishHook",
      bquote({
        assign("args", as.list(environment()), envir = .(captured))
        stop("captured fishHook score")
      }),
      run_somatic_regulatory_scoring(
        maf = make_dummy_maf(),
        reg_ref_path = "data/raw/Testing/reg_elements_valid.loc",
        score_model = "negbin",
        score_return_model = FALSE,
        score_nb = FALSE,
        score_iter = 5L,
        score_subsample = 250L,
        score_sets = score_sets,
        score_seed = 99L,
        score_verbose = FALSE,
        score_mc_cores = 2L,
        score_p_randomized = FALSE,
        score_class_return = FALSE,
        fishhook_args = list()
      )
    ),
    regexp = "captured fishHook score",
    fixed = TRUE
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_identical(args$model, "negbin")
  expect_identical(args$return.model, FALSE)
  expect_identical(args$nb, FALSE)
  expect_equal(args$iter, 5L)
  expect_equal(args$subsample, 250L)
  expect_identical(args$sets, score_sets)
  expect_equal(args$seed, 99L)
  expect_identical(args$verbose, FALSE)
  expect_equal(args$mc.cores, 2L)
  expect_identical(args$p.randomized, FALSE)
  expect_identical(args$classReturn, FALSE)
}

test_somatic_wrapper_surfaces_match_supported_api <- function() {
  gene_formals <- names(formals(run_somatic_gene_scoring))
  reg_formals <- names(formals(run_somatic_regulatory_scoring))
  prepare_formals <- names(formals(prepare_somatic_scores))

  expect_true(all(expected_dndscv_formals %in% gene_formals))
  expect_true(all(expected_fishhook_constructor_formals %in% reg_formals))
  expect_true(all(expected_fishhook_score_formals %in% reg_formals))
  expect_true(all(expected_prepare_somatic_formals %in% prepare_formals))
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
    regexp = expected_text,
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

test_extract_dndscv_gene_scores_onesided <- function() {
  dndscv_mock <- data.table(
    gene_name = c("TP53", "FBXW7"),
    ppos_cv = c(1e-6, 0.8),
    pneg_cv = c(0.9, 1e-4),
    wall_cv = c(2.4, 0.6)
  )

  scores <- extract_dndscv_gene_scores(dndscv_mock)

  expect_true(scores[gene_id == "TP53", zstat] > 0)
  expect_true(scores[gene_id == "FBXW7", zstat] < 0)
  expect_equal(scores[gene_id == "TP53", p_value], 1e-6)
  expect_equal(scores[gene_id == "FBXW7", p_value], 1e-4)
}

test_extract_fishhook_reg_scores <- function() {
  fishhook_mock <- data.table(
    reg_elem_id = c("EH38E0080197", "EH38E2084302"),
    p = c(1e-4, 0.03),
    effectsize = c(2.1, 1.2)
  )

  scores <- extract_fishhook_reg_scores(fishhook_mock)

  expect_identical(names(scores), c("reg_elem_id", "p_value", "zstat"))
  expect_identical(as.character(scores$reg_elem_id[[1]]), "EH38E0080197")
  expect_true(all(scores$p_value > 0))
  expect_true(all(scores$zstat > 0))
  invisible(scores)
}

test_extract_fishhook_reg_scores_directional <- function() {
  fishhook_mock <- data.table(
    reg_elem_id = c("EH38E0080197", "EH38E2084302"),
    p = c(1e-4, 0.9),
    p.neg = c(0.9, 1e-5),
    effectsize = c(2.1, 0.6)
  )

  scores <- extract_fishhook_reg_scores(fishhook_mock)

  expect_true(scores[reg_elem_id == "EH38E0080197", zstat] > 0)
  expect_true(scores[reg_elem_id == "EH38E2084302", zstat] < 0)
  expect_equal(scores[reg_elem_id == "EH38E0080197", p_value], 1e-4)
  expect_equal(scores[reg_elem_id == "EH38E2084302", p_value], 0.9)
}

test_somatic_extreme_scores_are_capped <- function() {
  capped_from_p <- compute_signed_z_from_p(c(0, 1e-400, 1e-20))
  expect_true(all(is.finite(capped_from_p)))

  fishhook_mock <- data.table(
    reg_elem_id = c("EH38E0080197", "EH38E2084302"),
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
    reg_elem_id = c("EH38E0080197", "EH38E2084302"),
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
  output_path <- make_somatic_test_path("dndscv_live_scores", ".tsv")
  result <- run_somatic_gene_scoring(maf, refdb = refdb, output_path = output_path)

  expect_s3_class(result, "conseguiR_bundle")
  expect_true(is.data.frame(result$gene_scores))
  expect_identical(names(result$gene_scores), c("gene_id", "p_value", "zstat"))
  expect_gt(nrow(result$gene_scores), 0L)
  expect_true(file.exists(output_path))
  message("Live dndscv gene scores:")
  print(result$gene_scores)
  invisible(result)
}

test_run_fishhook_reg_scoring_live <- function() {
  skip_if_not_installed("fishHook")
  skip_if_not_installed("BSgenome.Hsapiens.UCSC.hg38")

  overlap_n <- fishhook_covariate_overlap_count()
  skip_if(
    overlap_n == 0L,
    "fishHook live covariate file does not overlap the current regulatory universe."
  )

  maf <- fread(default_somatic_path, nrows = 1000L)
  fishhook_covariate_data <- readRDS(default_fishhook_covariate_path)
  output_path <- make_somatic_test_path("fishhook_live_scores", ".tsv")
  result <- run_somatic_regulatory_scoring(
    maf = maf,
    reg_ref_path = default_reg_ref_path,
    fishhook_covariate_data = fishhook_covariate_data,
    output_path = output_path
  )

  expect_s3_class(result, "conseguiR_bundle")
  expect_true(is.data.frame(result$reg_scores))
  expect_identical(names(result$reg_scores), c("reg_elem_id", "p_value", "zstat"))
  expect_gt(nrow(result$reg_scores), 0L)
  expect_true(file.exists(output_path))
  message("Live fishHook regulatory-element scores:")
  print(result$reg_scores)
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
  overlap_n <- fishhook_covariate_overlap_count()
  skip_if(
    overlap_n == 0L,
    "fishHook live covariate file does not overlap the current regulatory universe."
  )

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
  print_somatic_coverage_summary()

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

  test_that("dndscv runner forwards supported passthrough arguments", {
    test_run_dndscv_gene_scoring_forwards_supported_args()
  })

  test_that("fishHook runner reports clear errors for bad regulatory or covariate inputs", {
    test_run_fishhook_reg_scoring_negative_missing_covariate_column()
    test_run_fishhook_reg_scoring_negative_bad_reg_reference()
  })

  test_that("fishHook runner forwards supported constructor arguments", {
    test_run_fishhook_reg_scoring_forwards_constructor_args()
  })

  test_that("fishHook runner forwards supported score arguments", {
    test_run_fishhook_reg_scoring_forwards_score_args()
  })

  test_that("somatic wrapper surfaces expose the supported dndscv and fishHook APIs", {
    test_somatic_wrapper_surfaces_match_supported_api()
  })

  test_that("somatic score extraction works for dndscv genes", {
    test_extract_dndscv_gene_scores()
    test_extract_dndscv_gene_scores_onesided()
  })

  test_that("somatic score extraction works for fishHook regulatory elements", {
    test_extract_fishhook_reg_scores()
    test_extract_fishhook_reg_scores_directional()
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
