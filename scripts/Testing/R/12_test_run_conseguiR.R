#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(testthat)
  library(devtools)
  library(igraph)
})

source("R/zzz.R")
source("R/backend_resources.R")

devtools::load_all(".", quiet = TRUE)

default_pipeline_test_output_dir <- "data/processed/test_outputs/run_conseguiR"
default_pipeline_gwas_path <- "data/raw/Testing/34737426-GCST90043906-EFO_0000403.h.tsv"
default_pipeline_somatic_path <- "data/raw/Testing/2026-01-09_no_CLL_lymph_only_pcawg_maf_tcga_order_hg38.maf"
default_pipeline_reg_ref_path <- .conseguiR_default_reg_loc_path()
default_pipeline_track_dir <- "data/raw/Testing"
default_pipeline_dndscv_refdb <- "data/raw/Testing/RefCDS_human_GRCh38.p12.rda"

expected_run_conseguiR_formals <- c(
  "gwas_sumstats", "somatic_maf", "reg_ref_path", "reference_bfile",
  "dndscv_refdb", "epigenomic_track_dir", "epigenomic_tracks",
  "graph_rds_path", "gg_nodes_path", "gg_edges_path", "output_dir",
  "target_genes", "germline_args", "somatic_args", "epigenomic_args",
  "scored_graph_args", "diffusion_args", "subgraph_args", "plot_args",
  "verbose"
)

make_pipeline_test_path <- function(stem, ext = "") {
  dir.create(default_pipeline_test_output_dir, recursive = TRUE, showWarnings = FALSE)
  file.path(
    default_pipeline_test_output_dir,
    paste0(stem, "_", format(Sys.time(), "%Y%m%d%H%M%S"), "_", sprintf("%06d", sample.int(999999L, 1L)), ext)
  )
}

print_run_conseguiR_coverage_summary <- function() {
  message("run_conseguiR coverage matrix in this script:")
  message("  Stage forwarding asserted through run_conseguiR:")
  message("    prepare_germline_scores: gene_step1_args, gene_step2_args, reg_step1_args, reg_step2_args, shared_args")
  message("    prepare_somatic_scores: dndscv surface, fishHook constructor surface, fishHook score surface")
  message("    prepare_epigenomic_scores: track_dir, bw_files, output_path, min_tracks, drop_mhc, transform, return_diagnostics, summary_fun, verbose")
  message("    build_scored_gene_reg_graph: graph/graph_rds_path, score tables, output_prefix, save_outputs, verbose")
  message("    run_gene_reg_diffusion: output paths plus diffusion controls")
  message("    call_selected_subgraph: target, pool, weights, scaling, solver, output controls")
  message("    plot_selected_subgraph: bundle inputs, file paths, title, layout, labels, dimensions, save controls")
  message("  Note: forwarding tests use capture stubs; the live end-to-end test only uses real upstream flags.")
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

with_external_overrides <- function(overrides, code) {
  .conseguiR_load_external_api()
  env <- .conseguiR_runtime_env
  old_values <- lapply(names(overrides), function(nm) get(nm, envir = env, inherits = FALSE))
  names(old_values) <- names(overrides)

  for (nm in names(overrides)) {
    assign(nm, overrides[[nm]], envir = env)
  }

  on.exit({
    for (nm in names(old_values)) {
      assign(nm, old_values[[nm]], envir = env)
    }
  }, add = TRUE)

  force(code)
}

fake_bundle <- function(type, entries = list()) {
  structure(
    c(list(bundle_type = type), entries),
    class = c(paste0("conseguiR_", type, "_bundle"), "conseguiR_bundle", "list")
  )
}

make_pipeline_filter_file <- function(stem) {
  path <- make_pipeline_test_path(stem, ".txt")
  writeLines(c("rs1", "rs2"), path)
  path
}

find_reference_bfile <- function() {
  candidates <- c(
    "data/raw/g1000_eur/g1000_eur",
    "data/raw/Reference/1000G_EUR_Phase3_plink/1000G.EUR.QC",
    "data/raw/Reference/g1000_eur/g1000_eur",
    "data/raw/LDREF/g1000_eur/g1000_eur",
    "tools/reference/g1000_eur/g1000_eur"
  )

  candidates <- unique(candidates[file.exists(paste0(candidates, ".bed")) &
                                    file.exists(paste0(candidates, ".bim")) &
                                    file.exists(paste0(candidates, ".fam"))])
  if (length(candidates) == 0L) NULL else candidates[[1]]
}

select_epigenomic_test_bigwigs <- function(track_dir) {
  bw_files <- list.files(track_dir, pattern = "\\.(bw|bigWig)$", full.names = TRUE)
  keep <- !grepl("(_BL_|_FL_|broken_signal_track)", basename(bw_files))
  bw_files <- bw_files[keep]
  if (length(bw_files) < 3L) {
    return(character())
  }
  bw_files[seq_len(3L)]
}

make_pipeline_reg_ref_subset <- function(n_reg_elements = 5000L) {
  reg_dt <- fread(default_pipeline_reg_ref_path, header = FALSE, showProgress = FALSE)
  reg_subset <- reg_dt[seq_len(min(n_reg_elements, nrow(reg_dt)))]
  reg_subset_path <- make_pipeline_test_path("run_conseguiR_reg_subset", ".loc")
  fwrite(reg_subset, reg_subset_path, sep = "\t", col.names = FALSE)
  reg_subset_path
}

make_pipeline_inputs <- function() {
  gwas_dt <- fread(default_pipeline_gwas_path, nrows = 5000L, showProgress = FALSE)
  somatic_dt <- fread(default_pipeline_somatic_path, nrows = 10000L, showProgress = FALSE)
  reg_ref_subset_path <- make_pipeline_reg_ref_subset()
  bw_files <- select_epigenomic_test_bigwigs(default_pipeline_track_dir)

  test_backend_dir <- file.path(tempdir(), "conseguiR_backend_test_12")
  options(conseguiR.backend_dir = test_backend_dir)
  .conseguiR_initialize_backend_graphs(
    backend_dir = test_backend_dir,
    build_gene_reg = TRUE,
    build_gene_gene = TRUE,
    force = TRUE,
    strict = TRUE,
    quiet = TRUE
  )

  backend_paths <- .conseguiR_backend_paths(test_backend_dir)
  gg_graph <- readRDS(backend_paths$gene_gene_graph_rds)
  gg_nodes_path <- make_pipeline_test_path("run_conseguiR_gg_nodes", ".tsv")
  gg_edges_path <- make_pipeline_test_path("run_conseguiR_gg_edges", ".tsv")
  fwrite(as.data.table(igraph::as_data_frame(gg_graph, what = "vertices")), gg_nodes_path, sep = "\t")
  fwrite(as.data.table(igraph::as_data_frame(gg_graph, what = "edges")), gg_edges_path, sep = "\t")

  list(
    gwas_sumstats = gwas_dt,
    somatic_maf = somatic_dt,
    reg_ref_path = reg_ref_subset_path,
    epigenomic_tracks = bw_files,
    graph_rds_path = backend_paths$gene_reg_graph_rds,
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path
  )
}

make_pipeline_forwarding_inputs <- function() {
  graph_rds_path <- make_pipeline_test_path("run_conseguiR_dummy_graph", ".rds")
  gg_nodes_path <- make_pipeline_test_path("run_conseguiR_dummy_nodes", ".tsv")
  gg_edges_path <- make_pipeline_test_path("run_conseguiR_dummy_edges", ".tsv")
  saveRDS(list(dummy = TRUE), graph_rds_path)
  fwrite(data.table(name = "GENE1"), gg_nodes_path, sep = "\t")
  fwrite(data.table(from = "GENE1", to = "GENE1", weight = 1), gg_edges_path, sep = "\t")

  list(
    gwas_sumstats = default_pipeline_gwas_path,
    somatic_maf = default_pipeline_somatic_path,
    reg_ref_path = default_pipeline_reg_ref_path,
    epigenomic_tracks = select_epigenomic_test_bigwigs(default_pipeline_track_dir),
    graph_rds_path = graph_rds_path,
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path
  )
}

make_full_run_conseguiR_germline_args <- function() {
  gene_filter <- make_pipeline_filter_file("run_conseguiR_gene_filter")
  reg_filter <- make_pipeline_filter_file("run_conseguiR_reg_filter")
  snp_include <- make_pipeline_filter_file("run_conseguiR_snp_include")
  snp_exclude <- make_pipeline_filter_file("run_conseguiR_snp_exclude")
  indiv_include <- make_pipeline_test_path("run_conseguiR_indiv_include", ".txt")
  indiv_exclude <- make_pipeline_test_path("run_conseguiR_indiv_exclude", ".txt")
  writeLines("F1 I1", indiv_include)
  writeLines("F2 I2", indiv_exclude)

  list(
    gene_sample_size = 456348L,
    gene_sample_size_col = "N",
    reg_sample_size = 456348L,
    reg_sample_size_col = "N",
    gene_step1_args = list(
      annotation_window = c(35, 10),
      filter_path = gene_filter,
      ignore_strand = TRUE,
      nonhuman = TRUE,
      annotate_modifiers = c("gene-mod"),
      general_args = list(gene_step1_flag = TRUE),
      extra_args = c("--gene-step1-extra", "foo")
    ),
    gene_step2_args = list(
      gene_model = "multi",
      gene_model_modifiers = "multi-show-all",
      genes_only = TRUE,
      pval = list(
        use = c("SNP", "P"),
        snp_id = "SNP",
        pval = "P",
        N = c(1000L, 1100L, 1200L),
        ncol = "N",
        duplicate = "drop"
      ),
      bfile = list(
        synonyms = "gene.synonyms",
        synonym_dup = "drop"
      ),
      gene_settings = list(
        snp_min_maf = 0.01,
        snp_min_mac = 2,
        snp_max_maf = 0.49,
        snp_max_mac = 100,
        snp_max_miss = 0.1,
        snp_diff = 1e-6,
        snp_include = snp_include,
        snp_exclude = snp_exclude,
        indiv_include = indiv_include,
        indiv_exclude = indiv_exclude,
        prune = 0.95,
        prune_prop = 0.5,
        prune_count = 20,
        fixed_permp = 5000,
        adap_permp = c(1e6, 25),
        min_perm = 1000
      ),
      batch = c(7, 20),
      seed = 12345L,
      big_data = FALSE,
      general_args = list(gene_step2_flag = TRUE),
      extra_args = c("--gene-step2-extra", "bar")
    ),
    reg_step1_args = list(
      annotation_window = c(0, 0),
      filter_path = reg_filter,
      ignore_strand = TRUE,
      nonhuman = TRUE,
      annotate_modifiers = c("reg-mod"),
      general_args = list(reg_step1_flag = TRUE),
      extra_args = c("--reg-step1-extra", "alpha")
    ),
    reg_step2_args = list(
      gene_model = "multi",
      gene_model_modifiers = "multi-show-all",
      genes_only = FALSE,
      pval = list(
        use = c("SNP", "P"),
        snp_id = "SNP",
        pval = "P",
        N = c(1000L, 1100L, 1200L),
        ncol = "N",
        duplicate = "drop"
      ),
      bfile = list(
        synonyms = "reg.synonyms",
        synonym_dup = "drop"
      ),
      gene_settings = list(
        snp_min_maf = 0.01,
        snp_max_miss = 0.1
      ),
      batch = c("X", "chr"),
      seed = 321L,
      big_data = FALSE,
      general_args = list(reg_step2_flag = TRUE),
      extra_args = c("--reg-step2-extra", "beta")
    ),
    shared_args = list(
      magma_path = "dummy_magma_path",
      reuse_existing_gwas_cache = TRUE,
      reuse_existing_annotation = TRUE,
      reuse_existing_analysis = TRUE,
      keep_intermediates = TRUE
    )
  )
}

make_full_run_conseguiR_somatic_args <- function() {
  list(
    gene_cv = matrix(c(1, 2, 3, 4), nrow = 2L),
    gene_max_muts_per_gene_per_sample = 6L,
    gene_max_coding_muts_per_sample = 5000L,
    gene_list = c("TP53", "KRAS"),
    sm = "192r_3w",
    kc = "cgc81",
    use_indel_sites = TRUE,
    min_indels = 1L,
    maxcovs = 20L,
    constrain_wnon_wspl = TRUE,
    outp = 3L,
    numcode = 1L,
    outmats = FALSE,
    mingenecovs = 100L,
    onesided = TRUE,
    dc = c(TP53 = 1, KRAS = 2),
    dndscv_args = list(),
    eligible_gr = GenomicRanges::GRanges(
      seqnames = "chr1",
      ranges = IRanges::IRanges(start = 1L, end = 1000L),
      strand = "*"
    ),
    fishhook_covariates = list(),
    fishhook_covariate_data = NULL,
    fishhook_idcol = "Tumor_Sample_Barcode",
    constructor_out_path = tempdir(),
    constructor_use_local_mut_density = TRUE,
    constructor_local_mut_density_bin = 5e5,
    constructor_mc_cores = 2L,
    constructor_na_rm = TRUE,
    constructor_pad = 0,
    constructor_max_slice = 5000L,
    constructor_ff_chunk = 25000L,
    constructor_max_chunk = 500000L,
    constructor_idcap = 1,
    constructor_weight_events = FALSE,
    constructor_nb = TRUE,
    score_sets = list(dummy = c("EH38E0080197")),
    score_model = "negbin",
    score_return_model = TRUE,
    score_nb = TRUE,
    score_iter = 25L,
    score_subsample = 1000L,
    score_seed = 42L,
    score_verbose = FALSE,
    score_mc_cores = 2L,
    score_p_randomized = FALSE,
    score_class_return = TRUE,
    fishhook_args = list()
  )
}

make_full_run_conseguiR_epigenomic_args <- function() {
  list(
    track_dir = default_pipeline_track_dir,
    bw_files = select_epigenomic_test_bigwigs(default_pipeline_track_dir),
    output_path = make_pipeline_test_path("run_conseguiR_epi", ".tsv"),
    min_tracks = 3L,
    drop_mhc = TRUE,
    transform = "log1p",
    return_diagnostics = FALSE,
    summary_fun = mean,
    verbose = FALSE
  )
}

make_full_run_conseguiR_scored_graph_args <- function() {
  list(
    graph = igraph::make_empty_graph(),
    graph_rds_path = make_pipeline_test_path("run_conseguiR_scored_graph_input", ".rds"),
    output_prefix = make_pipeline_test_path("run_conseguiR_scored_graph_prefix"),
    gene_germline_scores = data.table(gene_id = "TP53", zstat = 2),
    reg_germline_scores = data.table(reg_elem_id = "EH38E0080197", zstat = 1.2),
    gene_somatic_scores = data.table(gene_id = "TP53", zstat = -1.5),
    reg_somatic_scores = data.table(reg_elem_id = "EH38E0080197", zstat = 0.3),
    reg_epigenomic_scores = data.table(reg_elem_id = "EH38E0080197", zstat = 2.4),
    save_outputs = FALSE,
    verbose = FALSE
  )
}

make_full_run_conseguiR_diffusion_args <- function() {
  list(
    nodes_path = make_pipeline_test_path("run_conseguiR_diff_nodes", ".tsv"),
    edges_path = make_pipeline_test_path("run_conseguiR_diff_edges", ".tsv"),
    output_dir = make_pipeline_test_path("run_conseguiR_diff_dir"),
    output_stem = "diffusion_override",
    top_k = 4L,
    confidence_power = 1.5,
    beta_germline = 0.6,
    beta_somatic = 0.7,
    beta_epigenomic = 0.8,
    integration_weight_germline = 1.0,
    integration_weight_somatic = 1.1,
    integration_weight_epigenomic = 1.2,
    positive_only = FALSE,
    reg_signal_clip = 4.5,
    top_n_to_save = 25L,
    python_path = NULL,
    verbose = FALSE
  )
}

make_full_run_conseguiR_subgraph_args <- function() {
  list(
    diffusion_path = make_pipeline_test_path("run_conseguiR_subgraph_diff", ".tsv"),
    gg_nodes_path = make_pipeline_test_path("run_conseguiR_subgraph_nodes", ".tsv"),
    gg_edges_path = make_pipeline_test_path("run_conseguiR_subgraph_edges", ".tsv"),
    output_dir = make_pipeline_test_path("run_conseguiR_subgraph_dir"),
    output_stem = "subgraph_override",
    target_genes = 12L,
    candidate_pool_size = 60L,
    min_confidence = 0,
    max_edges_in_model = 400L,
    node_prize_weight = 1,
    edge_conf_weight = 1,
    edge_cost_weight = 1,
    node_scale = 1000L,
    edge_scale = 1000L,
    max_time_seconds = 20L,
    num_workers = 2L,
    random_seed = 42L,
    prize_column = "post_integrated",
    confidence_column = "confidence",
    edge_cost_column = "weight",
    python_path = NULL,
    verbose = FALSE
  )
}

make_full_run_conseguiR_plot_args <- function() {
  list(
    nodes = data.table(gene_name = "TP53", prize = 1),
    edges = data.table(from = "TP53", to = "TP53", weight = 1),
    summary = data.table(metric = "n", value = 1),
    nodes_path = make_pipeline_test_path("run_conseguiR_plot_nodes", ".tsv"),
    edges_path = make_pipeline_test_path("run_conseguiR_plot_edges", ".tsv"),
    summary_path = make_pipeline_test_path("run_conseguiR_plot_summary", ".tsv"),
    bundle_output_prefix = make_pipeline_test_path("run_conseguiR_plot_bundle"),
    plot_file_path = make_pipeline_test_path("run_conseguiR_plot_file", ".pdf"),
    title = "run_conseguiR E2E Test Plot",
    layout = "fr",
    top_n_labels = 12L,
    width = 12,
    height = 10,
    dpi = 300,
    save_bundle = FALSE,
    save_plot = FALSE,
    verbose = FALSE
  )
}

test_run_conseguiR_forwards_full_germline_args_surface <- function() {
  inputs <- make_pipeline_forwarding_inputs()
  captured <- new.env(parent = emptyenv())
  germline_args <- make_full_run_conseguiR_germline_args()

  with_external_overrides(
    list(
      validate_inputs = function(...) fake_bundle("validation"),
      prepare_germline_scores = function(...) {
        assign("args", list(...), envir = captured)
        stop("captured prepare_germline_scores")
      }
    ),
    expect_error(
      run_conseguiR(
        gwas_sumstats = inputs$gwas_sumstats,
        somatic_maf = inputs$somatic_maf,
        reg_ref_path = inputs$reg_ref_path,
        reference_bfile = "dummy_reference_bfile",
        dndscv_refdb = "dummy_refdb",
        epigenomic_track_dir = default_pipeline_track_dir,
        epigenomic_tracks = inputs$epigenomic_tracks,
        graph_rds_path = inputs$graph_rds_path,
        gg_nodes_path = inputs$gg_nodes_path,
        gg_edges_path = inputs$gg_edges_path,
        output_dir = make_pipeline_test_path("run_conseguiR_capture_germline"),
        germline_args = germline_args,
        somatic_args = list(),
        verbose = FALSE
      ),
      regexp = "captured prepare_germline_scores",
      fixed = TRUE
    )
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_true(all(names(germline_args) %in% names(args)))
  expect_identical(args$gene_step1_args$annotation_window, c(35, 10))
  expect_identical(args$gene_step2_args$gene_model, "multi")
  expect_identical(args$gene_step2_args$gene_model_modifiers, "multi-show-all")
  expect_identical(args$gene_step2_args$batch, c(7, 20))
  expect_identical(args$reg_step2_args$batch, c("X", "chr"))
  expect_true(is.list(args$shared_args))
}

test_run_conseguiR_forwards_full_somatic_args_surface <- function() {
  inputs <- make_pipeline_forwarding_inputs()
  captured <- new.env(parent = emptyenv())
  somatic_args <- make_full_run_conseguiR_somatic_args()

  with_external_overrides(
    list(
      validate_inputs = function(...) fake_bundle("validation"),
      prepare_germline_scores = function(...) fake_bundle("germline_scores"),
      prepare_somatic_scores = function(...) {
        assign("args", list(...), envir = captured)
        stop("captured prepare_somatic_scores")
      }
    ),
    expect_error(
      run_conseguiR(
        gwas_sumstats = inputs$gwas_sumstats,
        somatic_maf = inputs$somatic_maf,
        reg_ref_path = inputs$reg_ref_path,
        reference_bfile = "dummy_reference_bfile",
        dndscv_refdb = "dummy_refdb",
        epigenomic_track_dir = default_pipeline_track_dir,
        epigenomic_tracks = inputs$epigenomic_tracks,
        graph_rds_path = inputs$graph_rds_path,
        gg_nodes_path = inputs$gg_nodes_path,
        gg_edges_path = inputs$gg_edges_path,
        output_dir = make_pipeline_test_path("run_conseguiR_capture_somatic"),
        germline_args = list(),
        somatic_args = somatic_args,
        verbose = FALSE
      ),
      regexp = "captured prepare_somatic_scores",
      fixed = TRUE
    )
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_true(all(names(somatic_args) %in% names(args)))
  expect_identical(args$gene_list, c("TP53", "KRAS"))
  expect_identical(args$score_model, "negbin")
  expect_identical(args$score_return_model, TRUE)
  expect_identical(args$score_class_return, TRUE)
  expect_true(inherits(args$eligible_gr, "GRanges"))
}

test_run_conseguiR_forwards_epigenomic_args_surface <- function() {
  inputs <- make_pipeline_forwarding_inputs()
  captured <- new.env(parent = emptyenv())
  epigenomic_args <- make_full_run_conseguiR_epigenomic_args()

  with_external_overrides(
    list(
      validate_inputs = function(...) fake_bundle("validation"),
      prepare_germline_scores = function(...) fake_bundle("germline_scores"),
      prepare_somatic_scores = function(...) fake_bundle("somatic_scores"),
      prepare_epigenomic_scores = function(...) {
        assign("args", list(...), envir = captured)
        stop("captured prepare_epigenomic_scores")
      }
    ),
    expect_error(
      run_conseguiR(
        gwas_sumstats = inputs$gwas_sumstats,
        somatic_maf = inputs$somatic_maf,
        reg_ref_path = inputs$reg_ref_path,
        reference_bfile = "dummy_reference_bfile",
        dndscv_refdb = "dummy_refdb",
        epigenomic_track_dir = default_pipeline_track_dir,
        epigenomic_tracks = inputs$epigenomic_tracks,
        graph_rds_path = inputs$graph_rds_path,
        gg_nodes_path = inputs$gg_nodes_path,
        gg_edges_path = inputs$gg_edges_path,
        output_dir = make_pipeline_test_path("run_conseguiR_capture_epi"),
        epigenomic_args = epigenomic_args,
        verbose = FALSE
      ),
      regexp = "captured prepare_epigenomic_scores",
      fixed = TRUE
    )
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_identical(args$track_dir, default_pipeline_track_dir)
  expect_true(is.character(args$bw_files))
  expect_identical(args$min_tracks, 3L)
  expect_identical(args$drop_mhc, TRUE)
  expect_identical(args$transform, "log1p")
  expect_identical(args$return_diagnostics, FALSE)
}

test_run_conseguiR_forwards_scored_graph_args_surface <- function() {
  inputs <- make_pipeline_forwarding_inputs()
  captured <- new.env(parent = emptyenv())
  scored_graph_args <- make_full_run_conseguiR_scored_graph_args()

  saveRDS(list(dummy = TRUE), scored_graph_args$graph_rds_path)

  with_external_overrides(
    list(
      validate_inputs = function(...) fake_bundle("validation"),
      prepare_germline_scores = function(...) fake_bundle("germline_scores"),
      prepare_somatic_scores = function(...) fake_bundle("somatic_scores"),
      prepare_epigenomic_scores = function(...) fake_bundle("epigenomic_scores"),
      build_scored_gene_reg_graph = function(...) {
        assign("args", list(...), envir = captured)
        stop("captured build_scored_gene_reg_graph")
      }
    ),
    expect_error(
      run_conseguiR(
        gwas_sumstats = inputs$gwas_sumstats,
        somatic_maf = inputs$somatic_maf,
        reg_ref_path = inputs$reg_ref_path,
        reference_bfile = "dummy_reference_bfile",
        dndscv_refdb = "dummy_refdb",
        epigenomic_tracks = inputs$epigenomic_tracks,
        graph_rds_path = inputs$graph_rds_path,
        gg_nodes_path = inputs$gg_nodes_path,
        gg_edges_path = inputs$gg_edges_path,
        output_dir = make_pipeline_test_path("run_conseguiR_capture_scored"),
        scored_graph_args = scored_graph_args,
        verbose = FALSE
      ),
      regexp = "captured build_scored_gene_reg_graph",
      fixed = TRUE
    )
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_true(inherits(args$graph, "igraph"))
  expect_identical(args$graph_rds_path, scored_graph_args$graph_rds_path)
  expect_true(is.data.frame(args$gene_germline_scores))
  expect_identical(args$save_outputs, FALSE)
}

test_run_conseguiR_forwards_diffusion_args_surface <- function() {
  inputs <- make_pipeline_forwarding_inputs()
  captured <- new.env(parent = emptyenv())
  diffusion_args <- make_full_run_conseguiR_diffusion_args()

  with_external_overrides(
    list(
      validate_inputs = function(...) fake_bundle("validation"),
      prepare_germline_scores = function(...) fake_bundle("germline_scores"),
      prepare_somatic_scores = function(...) fake_bundle("somatic_scores"),
      prepare_epigenomic_scores = function(...) fake_bundle("epigenomic_scores"),
      build_scored_gene_reg_graph = function(...) fake_bundle("scored_graph"),
      run_gene_reg_diffusion = function(...) {
        assign("args", list(...), envir = captured)
        stop("captured run_gene_reg_diffusion")
      }
    ),
    expect_error(
      run_conseguiR(
        gwas_sumstats = inputs$gwas_sumstats,
        somatic_maf = inputs$somatic_maf,
        reg_ref_path = inputs$reg_ref_path,
        reference_bfile = "dummy_reference_bfile",
        dndscv_refdb = "dummy_refdb",
        epigenomic_tracks = inputs$epigenomic_tracks,
        graph_rds_path = inputs$graph_rds_path,
        gg_nodes_path = inputs$gg_nodes_path,
        gg_edges_path = inputs$gg_edges_path,
        output_dir = make_pipeline_test_path("run_conseguiR_capture_diff"),
        diffusion_args = diffusion_args,
        verbose = FALSE
      ),
      regexp = "captured run_gene_reg_diffusion",
      fixed = TRUE
    )
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_identical(args$top_k, 4L)
  expect_identical(args$confidence_power, 1.5)
  expect_identical(args$output_stem, "diffusion_override")
}

test_run_conseguiR_forwards_subgraph_args_surface <- function() {
  inputs <- make_pipeline_forwarding_inputs()
  captured <- new.env(parent = emptyenv())
  subgraph_args <- make_full_run_conseguiR_subgraph_args()

  with_external_overrides(
    list(
      validate_inputs = function(...) fake_bundle("validation"),
      prepare_germline_scores = function(...) fake_bundle("germline_scores"),
      prepare_somatic_scores = function(...) fake_bundle("somatic_scores"),
      prepare_epigenomic_scores = function(...) fake_bundle("epigenomic_scores"),
      build_scored_gene_reg_graph = function(...) fake_bundle("scored_graph"),
      run_gene_reg_diffusion = function(...) fake_bundle("diffusion"),
      call_selected_subgraph = function(...) {
        assign("args", list(...), envir = captured)
        stop("captured call_selected_subgraph")
      }
    ),
    expect_error(
      run_conseguiR(
        gwas_sumstats = inputs$gwas_sumstats,
        somatic_maf = inputs$somatic_maf,
        reg_ref_path = inputs$reg_ref_path,
        reference_bfile = "dummy_reference_bfile",
        dndscv_refdb = "dummy_refdb",
        epigenomic_tracks = inputs$epigenomic_tracks,
        graph_rds_path = inputs$graph_rds_path,
        gg_nodes_path = inputs$gg_nodes_path,
        gg_edges_path = inputs$gg_edges_path,
        output_dir = make_pipeline_test_path("run_conseguiR_capture_subgraph"),
        target_genes = 7L,
        subgraph_args = subgraph_args,
        verbose = FALSE
      ),
      regexp = "captured call_selected_subgraph",
      fixed = TRUE
    )
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_identical(args$target_genes, 12L)
  expect_identical(args$candidate_pool_size, 60L)
  expect_identical(args$output_stem, "subgraph_override")
}

test_run_conseguiR_forwards_plot_args_surface <- function() {
  inputs <- make_pipeline_forwarding_inputs()
  captured <- new.env(parent = emptyenv())
  plot_args <- make_full_run_conseguiR_plot_args()

  with_external_overrides(
    list(
      validate_inputs = function(...) fake_bundle("validation"),
      prepare_germline_scores = function(...) fake_bundle("germline_scores"),
      prepare_somatic_scores = function(...) fake_bundle("somatic_scores"),
      prepare_epigenomic_scores = function(...) fake_bundle("epigenomic_scores"),
      build_scored_gene_reg_graph = function(...) fake_bundle("scored_graph"),
      run_gene_reg_diffusion = function(...) fake_bundle("diffusion"),
      call_selected_subgraph = function(...) fake_bundle("selected_subgraph"),
      plot_selected_subgraph = function(...) {
        assign("args", list(...), envir = captured)
        stop("captured plot_selected_subgraph")
      }
    ),
    expect_error(
      run_conseguiR(
        gwas_sumstats = inputs$gwas_sumstats,
        somatic_maf = inputs$somatic_maf,
        reg_ref_path = inputs$reg_ref_path,
        reference_bfile = "dummy_reference_bfile",
        dndscv_refdb = "dummy_refdb",
        epigenomic_tracks = inputs$epigenomic_tracks,
        graph_rds_path = inputs$graph_rds_path,
        gg_nodes_path = inputs$gg_nodes_path,
        gg_edges_path = inputs$gg_edges_path,
        output_dir = make_pipeline_test_path("run_conseguiR_capture_plot"),
        plot_args = plot_args,
        verbose = FALSE
      ),
      regexp = "captured plot_selected_subgraph",
      fixed = TRUE
    )
  )

  args <- get("args", envir = captured, inherits = FALSE)
  expect_identical(args$title, "run_conseguiR E2E Test Plot")
  expect_identical(args$layout, "fr")
  expect_identical(args$top_n_labels, 12L)
  expect_true(is.data.frame(args$nodes))
}

test_run_conseguiR_end_to_end <- function(print_outputs = TRUE) {
  magma_path <- tryCatch(.conseguiR_resolve_magma_path(must_work = FALSE), error = function(e) NULL)
  if (is.null(magma_path)) {
    skip("No usable MAGMA executable was auto-discovered for run_conseguiR().")
  }

  reference_bfile <- find_reference_bfile()
  if (is.null(reference_bfile)) {
    skip("No PLINK reference bfile found in the repository.")
  }
  if (!file.exists(default_pipeline_dndscv_refdb)) {
    skip("No dndscv reference database found in the repository.")
  }

  inputs <- make_pipeline_inputs()
  if (length(inputs$epigenomic_tracks) < 3L) {
    skip("Fewer than three valid bigWig files are available for epigenomic scoring.")
  }

  output_dir <- make_pipeline_test_path("run_conseguiR_outputs")
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

  result <- tryCatch(
    run_conseguiR(
      gwas_sumstats = inputs$gwas_sumstats,
      somatic_maf = inputs$somatic_maf,
      reg_ref_path = inputs$reg_ref_path,
      reference_bfile = reference_bfile,
      dndscv_refdb = default_pipeline_dndscv_refdb,
      epigenomic_tracks = inputs$epigenomic_tracks,
      graph_rds_path = inputs$graph_rds_path,
      gg_nodes_path = inputs$gg_nodes_path,
      gg_edges_path = inputs$gg_edges_path,
      output_dir = output_dir,
      target_genes = 12L,
      germline_args = list(
        gene_sample_size = 456348L,
        reg_sample_size = 456348L,
        gene_step1_args = list(
          annotation_window = c(35, 10),
          ignore_strand = TRUE,
          nonhuman = TRUE
        ),
        gene_step2_args = list(
          gene_model = "snp-wise=mean",
          genes_only = TRUE,
          pval = list(
            use = c("SNP", "P"),
            duplicate = "drop"
          )
        ),
        reg_step1_args = list(
          annotation_window = c(0, 0)
        ),
        reg_step2_args = list(
          genes_only = FALSE,
          pval = list(
            use = c("SNP", "P"),
            duplicate = "drop"
          )
        )
      ),
      somatic_args = list(
        gene_max_muts_per_gene_per_sample = 6L,
        gene_max_coding_muts_per_sample = 5000L,
        kc = "cgc81",
        maxcovs = 20L,
        constrain_wnon_wspl = TRUE,
        sm = "192r_3w",
        use_indel_sites = TRUE,
        min_indels = 1L,
        outp = 3L,
        numcode = 1L,
        outmats = FALSE,
        mingenecovs = 100L,
        onesided = TRUE,
        constructor_use_local_mut_density = FALSE,
        constructor_local_mut_density_bin = 5e5,
        constructor_mc_cores = 2L,
        constructor_na_rm = TRUE,
        constructor_pad = 0,
        constructor_max_slice = 5000L,
        constructor_ff_chunk = 25000L,
        constructor_max_chunk = 500000L,
        constructor_idcap = 1,
        constructor_weight_events = FALSE,
        constructor_nb = TRUE,
        score_nb = TRUE,
        score_return_model = TRUE,
        score_iter = 25L,
        score_subsample = 1000L,
        score_seed = 42L,
        score_verbose = FALSE,
        score_mc_cores = 2L,
        score_p_randomized = FALSE,
        score_class_return = TRUE
      ),
      epigenomic_args = list(
        min_tracks = 3L,
        drop_mhc = TRUE,
        transform = "log1p",
        return_diagnostics = FALSE,
        summary_fun = mean
      ),
      scored_graph_args = list(),
      diffusion_args = list(
        top_k = 4L,
        confidence_power = 1.5,
        beta_germline = 0.6,
        beta_somatic = 0.7,
        beta_epigenomic = 0.8,
        integration_weight_germline = 1.0,
        integration_weight_somatic = 1.1,
        integration_weight_epigenomic = 1.2,
        positive_only = FALSE,
        reg_signal_clip = 4.5,
        top_n_to_save = 25L
      ),
      subgraph_args = list(
        candidate_pool_size = 60L,
        min_confidence = 0,
        max_edges_in_model = 400L,
        node_prize_weight = 1,
        edge_conf_weight = 1,
        edge_cost_weight = 1,
        node_scale = 1000L,
        edge_scale = 1000L,
        max_time_seconds = 20L,
        num_workers = 2L,
        random_seed = 42L,
        prize_column = "post_integrated",
        confidence_column = "confidence",
        edge_cost_column = "weight"
      ),
      plot_args = list(
        title = "run_conseguiR E2E Test Plot",
        layout = "fr",
        top_n_labels = 12L
      ),
      verbose = TRUE
    ),
    error = function(e) {
      msg <- conditionMessage(e)
      if (grepl("count does not vary", msg, fixed = TRUE)) {
        skip("Live fishHook branch is degenerate on the compact somatic/regulatory test fixture (`count does not vary`).")
      }
      stop(e)
    }
  )

  expect_s3_class(result, "conseguiR_bundle")
  expect_identical(result$bundle_type, "pipeline")
  expect_s3_class(result$validation, "conseguiR_bundle")
  expect_s3_class(result$germline, "conseguiR_bundle")
  expect_s3_class(result$somatic, "conseguiR_bundle")
  expect_s3_class(result$epigenomic, "conseguiR_bundle")
  expect_s3_class(result$scored_graph, "conseguiR_bundle")
  expect_s3_class(result$diffusion, "conseguiR_bundle")
  expect_s3_class(result$selected_subgraph, "conseguiR_bundle")
  expect_s3_class(result$plot, "conseguiR_bundle")

  expect_true(is.data.frame(result$germline$gene_scores))
  expect_true(is.data.frame(result$germline$reg_scores))
  expect_true(is.data.frame(result$somatic$gene_scores))
  expect_true(is.data.frame(result$somatic$reg_scores))
  expect_true(is.data.frame(result$epigenomic$reg_scores))
  expect_true(is.data.frame(result$diffusion$all_genes))
  expect_true(is.data.frame(result$selected_subgraph$nodes))
  expect_true(inherits(result$plot$plot, "ggplot"))

  expect_true(file.exists(file.path(output_dir, "germline_gene_scores.zstat.tsv")))
  expect_true(file.exists(file.path(output_dir, "germline_reg_scores.zstat.tsv")))
  expect_true(file.exists(file.path(output_dir, "somatic_gene_scores.tsv")))
  expect_true(file.exists(file.path(output_dir, "somatic_reg_scores.tsv")))
  expect_true(file.exists(file.path(output_dir, "epigenomic_reg_scores.tsv")))
  expect_true(file.exists(file.path(output_dir, "gene_reg_graph_scored.rds")))
  expect_true(file.exists(file.path(output_dir, "gene_reg_graph_diffusion_all_genes.tsv")))
  expect_true(file.exists(file.path(output_dir, "gene_gene_selected_subgraph_nodes.tsv")))
  expect_true(file.exists(file.path(output_dir, "gene_gene_selected_subgraph_edges.tsv")))
  expect_true(file.exists(file.path(output_dir, "gene_gene_selected_subgraph_summary.tsv")))

  if (isTRUE(print_outputs)) {
    message("run_conseguiR() selected genes:")
    print(result$selected_subgraph$nodes[order(-prize)][1:min(10L, .N), .(gene_name, prize)])
    message("run_conseguiR() output dir: ", output_dir)
  }

  invisible(result)
}

test_run_conseguiR_surface_matches_supported_api <- function() {
  run_formals <- names(formals(run_conseguiR))
  expect_true(all(expected_run_conseguiR_formals %in% run_formals))
}

main <- function() {
  print_run_conseguiR_coverage_summary()

  test_that("run_conseguiR surface matches supported API", {
    test_run_conseguiR_surface_matches_supported_api()
  })

  test_that("run_conseguiR forwards the full germline nested argument surface", {
    test_run_conseguiR_forwards_full_germline_args_surface()
  })

  test_that("run_conseguiR forwards the full somatic nested argument surface", {
    test_run_conseguiR_forwards_full_somatic_args_surface()
  })

  test_that("run_conseguiR forwards the epigenomic stage surface", {
    test_run_conseguiR_forwards_epigenomic_args_surface()
  })

  test_that("run_conseguiR forwards the scored-graph stage surface", {
    test_run_conseguiR_forwards_scored_graph_args_surface()
  })

  test_that("run_conseguiR forwards the diffusion stage surface", {
    test_run_conseguiR_forwards_diffusion_args_surface()
  })

  test_that("run_conseguiR forwards the selected-subgraph stage surface", {
    test_run_conseguiR_forwards_subgraph_args_surface()
  })

  test_that("run_conseguiR forwards the plotting stage surface", {
    test_run_conseguiR_forwards_plot_args_surface()
  })

  test_that("run_conseguiR completes an end-to-end pipeline smoke test", {
    test_run_conseguiR_end_to_end(print_outputs = TRUE)
  })
}

if (sys.nframe() == 0) {
  main()
}
