#' @keywords internal
.conseguiR_runtime_env <- new.env(parent = baseenv())

#' @keywords internal
.conseguiR_runtime_file <- function(relpath) {
  candidates <- c(
    if (!is.null(.conseguiR_pkg_root)) file.path(.conseguiR_pkg_root, relpath),
    file.path(getwd(), relpath)
  )

  candidates <- unique(candidates[!is.na(candidates)])
  existing <- candidates[file.exists(candidates)]

  if (length(existing) == 0L) {
    stop("Could not locate required runtime file: ", relpath)
  }

  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_load_external_api <- function() {
  if (isTRUE(get0(".loaded", envir = .conseguiR_runtime_env, inherits = FALSE))) {
    return(invisible())
  }

  api_path <- .conseguiR_runtime_file("scripts/Externals/R/00_user_callable_functions.R")
  sys.source(api_path, envir = .conseguiR_runtime_env)
  assign(".loaded", TRUE, envir = .conseguiR_runtime_env)
  invisible()
}

#' @keywords internal
.conseguiR_external_fun <- function(name) {
  .conseguiR_load_external_api()

  if (!exists(name, envir = .conseguiR_runtime_env, inherits = FALSE)) {
    stop("External API function is not available: ", name)
  }

  get(name, envir = .conseguiR_runtime_env, inherits = FALSE)
}

#' Validate raw conseguiR inputs
#'
#' Validates raw GWAS, somatic, regulatory-reference, and epigenomic inputs
#' using the package's internal validation layer.
#'
#' @inheritParams .conseguiR_external_fun
#' @param gwas_sumstats GWAS summary statistics path or table.
#' @param somatic_maf Somatic MAF path or table.
#' @param reg_ref_path Regulatory-element reference path.
#' @param epigenomic_tracks Optional vector of bigWig paths.
#' @param epigenomic_track_dir Optional directory containing bigWig tracks.
#' @param epigenomic_exclude_patterns Patterns used to exclude bigWigs when
#'   discovering tracks from `epigenomic_track_dir`.
#'
#' @return A validation bundle containing validated objects and config.
#' @export
validate_inputs <- function(
  gwas_sumstats = NULL,
  somatic_maf = NULL,
  reg_ref_path = NULL,
  epigenomic_tracks = NULL,
  epigenomic_track_dir = NULL,
  epigenomic_exclude_patterns = c("_BL_", "_FL_")
) {
  .conseguiR_external_fun("validate_inputs")(
    gwas_sumstats = gwas_sumstats,
    somatic_maf = somatic_maf,
    reg_ref_path = reg_ref_path,
    epigenomic_tracks = epigenomic_tracks,
    epigenomic_track_dir = epigenomic_track_dir,
    epigenomic_exclude_patterns = epigenomic_exclude_patterns
  )
}

#' Initialize backend graph resources
#'
#' Builds the package's unscored backend graph resources when they are missing
#' and the required raw graph resources are available.
#'
#' @inheritParams initialize_backend_graphs
#'
#' @return A backend-initialization result describing the backend directory,
#'   output paths, and graph build status.
#' @export
initialize_backend_graphs <- function(
  backend_dir = NULL,
  build_gene_reg = TRUE,
  build_gene_gene = TRUE,
  force = FALSE,
  strict = TRUE,
  quiet = FALSE
) {
  .conseguiR_initialize_backend_graphs(
    backend_dir = backend_dir,
    build_gene_reg = build_gene_reg,
    build_gene_gene = build_gene_gene,
    force = force,
    strict = strict,
    quiet = quiet
  )
}

#' Run MAGMA germline gene scoring
#'
#' This wrapper exposes both MAGMA stages for the gene-level run:
#'
#' - step 1: annotation
#' - step 2: gene analysis
#'
#' Pass stage-specific MAGMA arguments through `step1_args` and `step2_args`.
#'
#' @inheritParams validate_inputs
#' @param gene_loc_path Gene location file for MAGMA step 1.
#' @param reference_bfile PLINK reference prefix for MAGMA step 2.
#' @param output_prefix Output prefix for saved artifacts.
#' @param sample_size Fixed sample size for MAGMA.
#' @param sample_size_col Optional sample size column name.
#' @param magma_path Path to the MAGMA executable.
#' @param magma_gwas_cache_prefix Optional shared MAGMA GWAS cache prefix.
#' @param reuse_existing_gwas_cache Whether to reuse the shared MAGMA cache.
#' @param reuse_existing_annotation Whether to reuse an existing MAGMA
#'   annotation output.
#' @param reuse_existing_analysis Whether to reuse an existing MAGMA gene
#'   analysis output.
#' @param keep_intermediates Whether to keep intermediate MAGMA files.
#' @param step1_args Named list of MAGMA step 1 arguments. Supported entries
#'   include `annotation_window`, `filter_path`, `ignore_strand`, `nonhuman`,
#'   and `extra_args`.
#' @param step2_args Named list of MAGMA step 2 arguments. Supported entries
#'   include `gene_model`, `genes_only`, `pval_use`, `pval_duplicate`,
#'   `bfile_synonyms`, `bfile_synonym_dup`, and `extra_args`.
#'
#' @return A germline gene score bundle.
#' @export
run_germline_gene_scoring <- function(
  gwas_sumstats,
  gene_loc_path,
  reference_bfile,
  output_prefix = "data/processed/germline_gene_scores",
  sample_size = NULL,
  sample_size_col = NULL,
  magma_path = "tools/magma_v1/magma",
  magma_gwas_cache_prefix = NULL,
  reuse_existing_gwas_cache = TRUE,
  reuse_existing_annotation = FALSE,
  reuse_existing_analysis = FALSE,
  keep_intermediates = FALSE,
  step1_args = list(),
  step2_args = list()
) {
  .conseguiR_external_fun("run_germline_gene_scoring")(
    gwas_sumstats = gwas_sumstats,
    gene_loc_path = gene_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = sample_size,
    sample_size_col = sample_size_col,
    magma_path = magma_path,
    magma_gwas_cache_prefix = magma_gwas_cache_prefix,
    reuse_existing_gwas_cache = reuse_existing_gwas_cache,
    reuse_existing_annotation = reuse_existing_annotation,
    reuse_existing_analysis = reuse_existing_analysis,
    keep_intermediates = keep_intermediates,
    step1_args = step1_args,
    step2_args = step2_args
  )
}

#' Run MAGMA germline regulatory scoring
#'
#' This wrapper exposes both MAGMA stages for the regulatory-element run:
#'
#' - step 1: annotation
#' - step 2: gene analysis
#'
#' Pass stage-specific MAGMA arguments through `step1_args` and `step2_args`.
#'
#' @inheritParams run_germline_gene_scoring
#' @param reg_loc_path Regulatory-element location file for MAGMA step 1.
#'
#' @return A germline regulatory score bundle.
#' @export
run_germline_regulatory_scoring <- function(
  gwas_sumstats,
  reg_loc_path,
  reference_bfile,
  output_prefix = "data/processed/germline_reg_scores",
  sample_size = NULL,
  sample_size_col = NULL,
  magma_path = "tools/magma_v1/magma",
  magma_gwas_cache_prefix = NULL,
  reuse_existing_gwas_cache = TRUE,
  reuse_existing_annotation = FALSE,
  reuse_existing_analysis = FALSE,
  keep_intermediates = FALSE,
  step1_args = list(),
  step2_args = list()
) {
  .conseguiR_external_fun("run_germline_regulatory_scoring")(
    gwas_sumstats = gwas_sumstats,
    reg_loc_path = reg_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    sample_size = sample_size,
    sample_size_col = sample_size_col,
    magma_path = magma_path,
    magma_gwas_cache_prefix = magma_gwas_cache_prefix,
    reuse_existing_gwas_cache = reuse_existing_gwas_cache,
    reuse_existing_annotation = reuse_existing_annotation,
    reuse_existing_analysis = reuse_existing_analysis,
    keep_intermediates = keep_intermediates,
    step1_args = step1_args,
    step2_args = step2_args
  )
}

#' Prepare germline scores for genes and regulatory elements
#'
#' This wrapper orchestrates two separate MAGMA runs:
#'
#' - one gene-level run with its own MAGMA step 1 and step 2 settings
#' - one regulatory-level run with its own MAGMA step 1 and step 2 settings
#'
#' In other words, MAGMA stage customization is exposed twice:
#'
#' - `gene_step1_args` and `gene_step2_args`
#' - `reg_step1_args` and `reg_step2_args`
#'
#' @inheritParams run_germline_gene_scoring
#' @param reg_loc_path Regulatory-element location file.
#' @param gene_output_prefix Output prefix for gene-level germline scores.
#' @param reg_output_prefix Output prefix for regulatory germline scores.
#' @param gene_sample_size Fixed sample size for the gene run.
#' @param gene_sample_size_col Optional sample size column for the gene run.
#' @param reg_sample_size Fixed sample size for the regulatory run.
#' @param reg_sample_size_col Optional sample size column for the regulatory run.
#' @param gene_step1_args Named list of gene MAGMA step 1 arguments.
#' @param gene_step2_args Named list of gene MAGMA step 2 arguments.
#' @param reg_step1_args Named list of regulatory MAGMA step 1 arguments.
#' @param reg_step2_args Named list of regulatory MAGMA step 2 arguments.
#' @param shared_args Named list of arguments passed to both runs.
#'
#' @return A germline score bundle with gene and regulatory score tables.
#' @export
prepare_germline_scores <- function(
  gwas_sumstats,
  reference_bfile,
  gene_loc_path = "data/raw/NCBI38/NCBI38.gene.loc",
  reg_loc_path = "data/raw/GeneHancer/2026-01-26_UCSC_all_unfiltered_reg_elements.loc",
  gene_output_prefix = "data/processed/germline_gene_scores",
  reg_output_prefix = "data/processed/germline_reg_scores",
  magma_gwas_cache_prefix = "data/processed/magma_shared_gwas_cache",
  gene_sample_size = NULL,
  gene_sample_size_col = NULL,
  reg_sample_size = NULL,
  reg_sample_size_col = NULL,
  gene_step1_args = list(),
  gene_step2_args = list(),
  reg_step1_args = list(),
  reg_step2_args = list(),
  shared_args = list()
) {
  .conseguiR_external_fun("prepare_germline_scores")(
    gwas_sumstats = gwas_sumstats,
    reference_bfile = reference_bfile,
    gene_loc_path = gene_loc_path,
    reg_loc_path = reg_loc_path,
    gene_output_prefix = gene_output_prefix,
    reg_output_prefix = reg_output_prefix,
    magma_gwas_cache_prefix = magma_gwas_cache_prefix,
    gene_sample_size = gene_sample_size,
    gene_sample_size_col = gene_sample_size_col,
    reg_sample_size = reg_sample_size,
    reg_sample_size_col = reg_sample_size_col,
    gene_step1_args = gene_step1_args,
    gene_step2_args = gene_step2_args,
    reg_step1_args = reg_step1_args,
    reg_step2_args = reg_step2_args,
    shared_args = shared_args
  )
}

#' Run dndscv somatic gene scoring
#'
#' This wrapper exposes the core dndscv settings used by the current pipeline
#' and leaves less common controls available through `dndscv_args`.
#'
#' @param maf Somatic MAF path or table.
#' @param refdb dndscv reference database path.
#' @param output_path Optional output path for saved scores.
#' @param cv Optional dndscv covariates.
#' @param max_muts_per_gene_per_sample dndscv gene-level mutation cap.
#' @param max_coding_muts_per_sample dndscv sample-level coding mutation cap.
#' @param dndscv_args Named list of additional dndscv arguments.
#'
#' @return A somatic gene score bundle.
#' @export
run_somatic_gene_scoring <- function(
  maf,
  refdb,
  output_path = "data/processed/somatic_gene_scores.tsv",
  cv = NULL,
  max_muts_per_gene_per_sample = 6L,
  max_coding_muts_per_sample = 5000L,
  dndscv_args = list()
) {
  .conseguiR_external_fun("run_somatic_gene_scoring")(
    maf = maf,
    refdb = refdb,
    output_path = output_path,
    cv = cv,
    max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,
    max_coding_muts_per_sample = max_coding_muts_per_sample,
    dndscv_args = dndscv_args
  )
}

#' Run fishHook somatic regulatory scoring
#'
#' This wrapper exposes the main fishHook configuration surface used by the
#' package workflow and leaves less common model tuning available through
#' `fishhook_args`.
#'
#' @inheritParams run_somatic_gene_scoring
#' @param reg_ref_path Regulatory-element reference path.
#' @param eligible_gr Optional fishHook eligible territory.
#' @param fishhook_covariates Optional fishHook covariate objects/specifications.
#' @param fishhook_covariate_data Optional tabular covariate data.
#' @param idcol Sample identifier column for fishHook.
#' @param fishhook_args Named list of additional fishHook scoring arguments.
#'
#' @return A somatic regulatory score bundle.
#' @export
run_somatic_regulatory_scoring <- function(
  maf,
  reg_ref_path,
  output_path = "data/processed/somatic_reg_scores.tsv",
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  idcol = "Tumor_Sample_Barcode",
  fishhook_args = list()
) {
  .conseguiR_external_fun("run_somatic_regulatory_scoring")(
    maf = maf,
    reg_ref_path = reg_ref_path,
    output_path = output_path,
    eligible_gr = eligible_gr,
    fishhook_covariates = fishhook_covariates,
    fishhook_covariate_data = fishhook_covariate_data,
    idcol = idcol,
    fishhook_args = fishhook_args
  )
}

#' Prepare somatic scores for genes and regulatory elements
#'
#' This wrapper orchestrates two distinct somatic scoring branches:
#'
#' - a gene-level `dndscv` run
#' - a regulatory `fishHook` run
#'
#' The dndscv and fishHook settings are kept separate so users can tune each
#' modeling framework independently.
#'
#' @inheritParams run_somatic_gene_scoring
#' @param reg_ref_path Regulatory-element reference path.
#' @param gene_output_path Output path for somatic gene scores.
#' @param reg_output_path Output path for somatic regulatory scores.
#' @param gene_cv Optional dndscv covariates.
#' @param gene_max_muts_per_gene_per_sample dndscv gene-level mutation cap.
#' @param gene_max_coding_muts_per_sample dndscv sample-level coding mutation cap.
#' @param eligible_gr Optional fishHook eligible territory.
#' @param fishhook_covariates Optional fishHook covariate objects/specifications.
#' @param fishhook_covariate_data Optional tabular covariate data.
#' @param fishhook_idcol Sample identifier column for fishHook.
#' @param fishhook_args Named list of additional fishHook arguments.
#'
#' @return A somatic score bundle with gene and regulatory score tables.
#' @export
prepare_somatic_scores <- function(
  maf,
  refdb,
  reg_ref_path,
  gene_output_path = "data/processed/somatic_gene_scores.tsv",
  reg_output_path = "data/processed/somatic_reg_scores.tsv",
  gene_cv = NULL,
  gene_max_muts_per_gene_per_sample = 6L,
  gene_max_coding_muts_per_sample = 5000L,
  dndscv_args = list(),
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  fishhook_idcol = "Tumor_Sample_Barcode",
  fishhook_args = list()
) {
  .conseguiR_external_fun("prepare_somatic_scores")(
    maf = maf,
    refdb = refdb,
    reg_ref_path = reg_ref_path,
    gene_output_path = gene_output_path,
    reg_output_path = reg_output_path,
    gene_cv = gene_cv,
    gene_max_muts_per_gene_per_sample = gene_max_muts_per_gene_per_sample,
    gene_max_coding_muts_per_sample = gene_max_coding_muts_per_sample,
    dndscv_args = dndscv_args,
    eligible_gr = eligible_gr,
    fishhook_covariates = fishhook_covariates,
    fishhook_covariate_data = fishhook_covariate_data,
    fishhook_idcol = fishhook_idcol,
    fishhook_args = fishhook_args
  )
}

#' Prepare regulatory epigenomic scores
#'
#' @param reg_ref_path Regulatory-element reference path.
#' @param track_dir Optional directory containing bigWig tracks.
#' @param bw_files Optional explicit vector of bigWig paths.
#' @param output_path Optional output path for saved scores.
#' @param exclude_patterns Patterns used to exclude tracks discovered from a
#'   directory.
#' @param min_tracks Minimum number of tracks required.
#' @param drop_mhc Whether to exclude extended MHC elements.
#' @param transform Signal transformation method.
#' @param return_diagnostics Whether to return the diagnostic table and raw
#'   matrix alongside z scores.
#' @param summary_fun Summary function applied to each bigWig over each
#'   regulatory element.
#'
#' @return An epigenomic score bundle.
#' @export
prepare_epigenomic_scores <- function(
  reg_ref_path,
  track_dir = NULL,
  bw_files = NULL,
  output_path = "data/processed/epigenomic_reg_scores.tsv",
  exclude_patterns = c("_BL_", "_FL_"),
  min_tracks = 3L,
  drop_mhc = TRUE,
  transform = "log1p",
  return_diagnostics = TRUE,
  summary_fun = mean
) {
  .conseguiR_external_fun("prepare_epigenomic_scores")(
    reg_ref_path = reg_ref_path,
    track_dir = track_dir,
    bw_files = bw_files,
    output_path = output_path,
    exclude_patterns = exclude_patterns,
    min_tracks = min_tracks,
    drop_mhc = drop_mhc,
    transform = transform,
    return_diagnostics = return_diagnostics,
    summary_fun = summary_fun
  )
}

#' Build a scored gene-regulatory graph
#'
#' @param graph Optional in-memory no-score `igraph`.
#' @param graph_rds_path Path to the backend no-score gene-reg graph.
#' @param output_prefix Output prefix for the saved scored graph artifacts.
#' @param germline_scores Optional germline score bundle.
#' @param somatic_scores Optional somatic score bundle.
#' @param epigenomic_scores Optional epigenomic score bundle.
#' @param gene_germline_scores Optional explicit gene germline score table.
#' @param reg_germline_scores Optional explicit regulatory germline score table.
#' @param gene_somatic_scores Optional explicit gene somatic score table.
#' @param reg_somatic_scores Optional explicit regulatory somatic score table.
#' @param reg_epigenomic_scores Optional explicit regulatory epigenomic score
#'   table.
#' @param save_outputs Whether to save the scored graph to disk.
#'
#' @return A scored gene-reg graph bundle containing the graph, nodes, and edges.
#' @export
build_scored_gene_reg_graph <- function(
  graph = NULL,
  graph_rds_path = NULL,
  output_prefix = "data/processed/gene_reg_graph_scored",
  germline_scores = NULL,
  somatic_scores = NULL,
  epigenomic_scores = NULL,
  gene_germline_scores = NULL,
  reg_germline_scores = NULL,
  gene_somatic_scores = NULL,
  reg_somatic_scores = NULL,
  reg_epigenomic_scores = NULL,
  save_outputs = TRUE
) {
  if (is.null(graph_rds_path)) {
    initialize_backend_graphs(strict = FALSE, quiet = TRUE)
    graph_rds_path <- .conseguiR_backend_paths()$gene_reg_graph_rds
  }

  .conseguiR_external_fun("build_scored_gene_reg_graph")(
    graph = graph,
    graph_rds_path = graph_rds_path,
    output_prefix = output_prefix,
    germline_scores = germline_scores,
    somatic_scores = somatic_scores,
    epigenomic_scores = epigenomic_scores,
    gene_germline_scores = gene_germline_scores,
    reg_germline_scores = reg_germline_scores,
    gene_somatic_scores = gene_somatic_scores,
    reg_somatic_scores = reg_somatic_scores,
    reg_epigenomic_scores = reg_epigenomic_scores,
    save_outputs = save_outputs
  )
}

#' Run diffusion on a scored gene-regulatory graph
#'
#' This wrapper exposes the main hyperparameters used by the current
#' regulatory-to-gene diffusion implementation.
#'
#' @param scored_graph Optional scored graph bundle.
#' @param nodes_path Optional explicit scored node table path.
#' @param edges_path Optional explicit scored edge table path.
#' @param output_dir Output directory for diffusion artifacts.
#' @param output_stem Output stem for diffusion artifacts.
#' @param top_k Diffusion top-k parameter.
#' @param confidence_power Edge-confidence exponent.
#' @param beta_germline Germline contribution weight.
#' @param beta_somatic Somatic contribution weight.
#' @param beta_epigenomic Epigenomic contribution weight.
#' @param positive_only Whether to restrict to positive regulatory signal.
#' @param reg_signal_clip Regulatory signal clip value.
#' @param top_n_to_save Number of top genes to save separately.
#' @param python_path Optional explicit Python interpreter path.
#'
#' @return A diffusion bundle containing full and top-gene diffusion tables.
#' @export
run_gene_reg_diffusion <- function(
  scored_graph = NULL,
  nodes_path = NULL,
  edges_path = NULL,
  output_dir = "data/processed",
  output_stem = "gene_reg_graph_diffusion",
  top_k = 3L,
  confidence_power = 2.0,
  beta_germline = 0.5,
  beta_somatic = 0.5,
  beta_epigenomic = 0.7,
  positive_only = FALSE,
  reg_signal_clip = 5.0,
  top_n_to_save = 50L,
  python_path = NULL
) {
  .conseguiR_external_fun("run_gene_reg_diffusion")(
    scored_graph = scored_graph,
    nodes_path = nodes_path,
    edges_path = edges_path,
    output_dir = output_dir,
    output_stem = output_stem,
    top_k = top_k,
    confidence_power = confidence_power,
    beta_germline = beta_germline,
    beta_somatic = beta_somatic,
    beta_epigenomic = beta_epigenomic,
    positive_only = positive_only,
    reg_signal_clip = reg_signal_clip,
    top_n_to_save = top_n_to_save,
    python_path = python_path
  )
}

#' Call the selected subgraph from diffusion results
#'
#' This wrapper exposes the main candidate-selection, objective-weighting, and
#' solver controls used by the current subgraph-calling implementation.
#'
#' @param diffusion Optional diffusion bundle.
#' @param diffusion_path Optional explicit diffusion table path.
#' @param gg_nodes_path Gene-gene node table path.
#' @param gg_edges_path Gene-gene edge table path.
#' @param output_dir Output directory for selected-subgraph artifacts.
#' @param output_stem Output stem for selected-subgraph artifacts.
#' @param target_genes Requested number of genes in the selected subgraph.
#' @param candidate_pool_size Candidate pool size for the solver.
#' @param min_confidence Minimum edge confidence allowed into the solver.
#' @param max_edges_in_model Maximum number of edges in the optimization model.
#' @param node_prize_weight Node prize weight.
#' @param edge_conf_weight Edge confidence weight.
#' @param edge_cost_weight Edge cost weight.
#' @param node_scale Integer scaling factor for node prizes.
#' @param edge_scale Integer scaling factor for edge terms.
#' @param max_time_seconds Solver time limit in seconds.
#' @param num_workers Solver worker count.
#' @param random_seed Solver random seed.
#' @param prize_column Diffusion prize column to optimize.
#' @param confidence_column Gene-gene confidence column.
#' @param edge_cost_column Gene-gene edge-cost column.
#' @param python_path Optional explicit Python interpreter path.
#'
#' @return A selected subgraph bundle.
#' @export
call_selected_subgraph <- function(
  diffusion = NULL,
  diffusion_path = NULL,
  gg_nodes_path = NULL,
  gg_edges_path = NULL,
  output_dir = "data/processed",
  output_stem = "gene_gene_selected_subgraph",
  target_genes = 50L,
  candidate_pool_size = 400L,
  min_confidence = 0,
  max_edges_in_model = 12000L,
  node_prize_weight = 1,
  edge_conf_weight = 1,
  edge_cost_weight = 1,
  node_scale = 1000L,
  edge_scale = 1000L,
  max_time_seconds = 600L,
  num_workers = 8L,
  random_seed = 42L,
  prize_column = "post_norm",
  confidence_column = "confidence",
  edge_cost_column = "weight",
  python_path = NULL
) {
  if (is.null(gg_nodes_path) || is.null(gg_edges_path)) {
    initialize_backend_graphs(strict = FALSE, quiet = TRUE)
    backend_paths <- .conseguiR_backend_paths()
    gg_nodes_path <- gg_nodes_path %||% backend_paths$gene_gene_graph_nodes
    gg_edges_path <- gg_edges_path %||% backend_paths$gene_gene_graph_edges
  }

  .conseguiR_external_fun("call_selected_subgraph")(
    diffusion = diffusion,
    diffusion_path = diffusion_path,
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path,
    output_dir = output_dir,
    output_stem = output_stem,
    target_genes = target_genes,
    candidate_pool_size = candidate_pool_size,
    min_confidence = min_confidence,
    max_edges_in_model = max_edges_in_model,
    node_prize_weight = node_prize_weight,
    edge_conf_weight = edge_conf_weight,
    edge_cost_weight = edge_cost_weight,
    node_scale = node_scale,
    edge_scale = edge_scale,
    max_time_seconds = max_time_seconds,
    num_workers = num_workers,
    random_seed = random_seed,
    prize_column = prize_column,
    confidence_column = confidence_column,
    edge_cost_column = edge_cost_column,
    python_path = python_path
  )
}

#' Plot a selected subgraph and build a visualization bundle
#'
#' This wrapper exposes the main input-resolution, saving, layout, labelling,
#' and rendering controls used by the current plotting implementation.
#'
#' @param selected_subgraph Optional selected subgraph bundle.
#' @param nodes Optional selected-subgraph node table.
#' @param edges Optional selected-subgraph edge table.
#' @param summary Optional selected-subgraph summary table.
#' @param nodes_path Optional node table path.
#' @param edges_path Optional edge table path.
#' @param summary_path Optional summary table path.
#' @param bundle_output_prefix Output prefix for the saved visualization bundle.
#' @param plot_file_path Optional output path for the saved figure.
#' @param title Plot title.
#' @param layout Graph layout name.
#' @param top_n_labels Number of node labels to show.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_bundle Whether to save the visualization bundle.
#' @param save_plot Whether to save the figure.
#'
#' @return A plot bundle containing the ggplot object and visualization bundle.
#' @export
plot_selected_subgraph <- function(
  selected_subgraph = NULL,
  nodes = NULL,
  edges = NULL,
  summary = NULL,
  nodes_path = NULL,
  edges_path = NULL,
  summary_path = NULL,
  bundle_output_prefix = "data/processed/gene_gene_selected_subgraph_plot_bundle",
  plot_file_path = NULL,
  title = "conseguiR Selected Gene Subgraph",
  layout = "fr",
  top_n_labels = Inf,
  width = 12,
  height = 10,
  dpi = 300,
  save_bundle = TRUE,
  save_plot = !is.null(plot_file_path)
) {
  .conseguiR_external_fun("plot_selected_subgraph")(
    selected_subgraph = selected_subgraph,
    nodes = nodes,
    edges = edges,
    summary = summary,
    nodes_path = nodes_path,
    edges_path = edges_path,
    summary_path = summary_path,
    bundle_output_prefix = bundle_output_prefix,
    plot_file_path = plot_file_path,
    title = title,
    layout = layout,
    top_n_labels = top_n_labels,
    width = width,
    height = height,
    dpi = dpi,
    save_bundle = save_bundle,
    save_plot = save_plot
  )
}

#' Run the full conseguiR pipeline end to end
#'
#' @param gwas_sumstats GWAS summary statistics path or table.
#' @param somatic_maf Somatic MAF path or table.
#' @param reg_ref_path Regulatory-element reference path.
#' @param reference_bfile PLINK reference prefix for MAGMA.
#' @param dndscv_refdb dndscv reference database path.
#' @param epigenomic_track_dir Optional directory containing epigenomic tracks.
#' @param epigenomic_tracks Optional explicit vector of bigWig paths.
#' @param gene_loc_path Gene location file path.
#' @param reg_loc_path Regulatory-element location file path.
#' @param graph_rds_path Backend no-score gene-reg graph path.
#' @param gg_nodes_path Gene-gene node table path.
#' @param gg_edges_path Gene-gene edge table path.
#' @param output_dir Output directory for pipeline artifacts.
#' @param target_genes Requested selected-subgraph size.
#' @param germline_args Named list of overrides passed to
#'   `prepare_germline_scores()`.
#' @param somatic_args Named list of overrides passed to
#'   `prepare_somatic_scores()`.
#' @param epigenomic_args Named list of overrides passed to
#'   `prepare_epigenomic_scores()`.
#' @param scored_graph_args Named list of overrides passed to
#'   `build_scored_gene_reg_graph()`.
#' @param diffusion_args Named list of overrides passed to
#'   `run_gene_reg_diffusion()`.
#' @param subgraph_args Named list of overrides passed to
#'   `call_selected_subgraph()`.
#' @param plot_args Named list of overrides passed to
#'   `plot_selected_subgraph()`.
#'
#' @return A pipeline bundle containing all stage bundles.
#' @export
run_conseguiR <- function(
  gwas_sumstats,
  somatic_maf,
  reg_ref_path,
  reference_bfile,
  dndscv_refdb,
  epigenomic_track_dir = NULL,
  epigenomic_tracks = NULL,
  gene_loc_path = "data/raw/NCBI38/NCBI38.gene.loc",
  reg_loc_path = "data/raw/GeneHancer/2026-01-26_UCSC_all_unfiltered_reg_elements.loc",
  graph_rds_path = NULL,
  gg_nodes_path = NULL,
  gg_edges_path = NULL,
  output_dir = "data/processed",
  target_genes = 50L,
  germline_args = list(),
  somatic_args = list(),
  epigenomic_args = list(),
  scored_graph_args = list(),
  diffusion_args = list(),
  subgraph_args = list(),
  plot_args = list()
) {
  initialize_backend_graphs(strict = FALSE, quiet = TRUE)
  backend_paths <- .conseguiR_backend_paths()
  graph_rds_path <- graph_rds_path %||% backend_paths$gene_reg_graph_rds
  gg_nodes_path <- gg_nodes_path %||% backend_paths$gene_gene_graph_nodes
  gg_edges_path <- gg_edges_path %||% backend_paths$gene_gene_graph_edges

  .conseguiR_external_fun("run_conseguiR")(
    gwas_sumstats = gwas_sumstats,
    somatic_maf = somatic_maf,
    reg_ref_path = reg_ref_path,
    reference_bfile = reference_bfile,
    dndscv_refdb = dndscv_refdb,
    epigenomic_track_dir = epigenomic_track_dir,
    epigenomic_tracks = epigenomic_tracks,
    gene_loc_path = gene_loc_path,
    reg_loc_path = reg_loc_path,
    graph_rds_path = graph_rds_path,
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path,
    output_dir = output_dir,
    target_genes = target_genes,
    germline_args = germline_args,
    somatic_args = somatic_args,
    epigenomic_args = epigenomic_args,
    scored_graph_args = scored_graph_args,
    diffusion_args = diffusion_args,
    subgraph_args = subgraph_args,
    plot_args = plot_args
  )
}
