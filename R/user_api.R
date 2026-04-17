#' @keywords internal
.conseguiR_runtime_env <- new.env(parent = globalenv())

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
.conseguiR_runtime_file <- function(relpath) {
  pkg_root <- .conseguiR_state$pkg_root
  candidates <- c(
    if (!is.null(pkg_root)) file.path(pkg_root, relpath),
    file.path(getwd(), relpath),
    {
      pkg_path <- system.file(relpath, package = "conseguiR")
      if (nzchar(pkg_path)) pkg_path else NA_character_
    }
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

#' @keywords internal
.conseguiR_call_external <- function(name, args) {
  do.call(.conseguiR_external_fun(name), args)
}

#' @keywords internal
.conseguiR_default_gene_gene_paths <- function(
  gg_nodes_path = NULL,
  gg_edges_path = NULL
) {
  if (is.null(gg_nodes_path) || is.null(gg_edges_path)) {
    initialize_backend_graphs(strict = FALSE, quiet = TRUE)
    backend_paths <- .conseguiR_backend_paths()
    gg_nodes_path <- gg_nodes_path %||% backend_paths$gene_gene_graph_nodes
    gg_edges_path <- gg_edges_path %||% backend_paths$gene_gene_graph_edges
  }

  list(
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path
  )
}

#' @keywords internal
.conseguiR_pipeline_paths <- function(
  graph_rds_path = NULL,
  gg_nodes_path = NULL,
  gg_edges_path = NULL
) {
  initialize_backend_graphs(strict = FALSE, quiet = TRUE)
  backend_paths <- .conseguiR_backend_paths()
  gene_gene_paths <- .conseguiR_default_gene_gene_paths(
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path
  )

  list(
    graph_rds_path = graph_rds_path %||% backend_paths$gene_reg_graph_rds,
    gg_nodes_path = gene_gene_paths$gg_nodes_path,
    gg_edges_path = gene_gene_paths$gg_edges_path
  )
}

#' @keywords internal
.conseguiR_pipeline_args <- function(
  gwas_sumstats,
  somatic_maf,
  reg_ref_path,
  reference_bfile,
  dndscv_refdb,
  epigenomic_track_dir,
  epigenomic_tracks,
  paths,
  output_dir,
  target_genes,
  germline_args,
  somatic_args,
  epigenomic_args,
  scored_graph_args,
  diffusion_args,
  subgraph_args,
  plot_args,
  verbose
) {
  list(
    gwas_sumstats = gwas_sumstats,
    somatic_maf = somatic_maf,
    reg_ref_path = reg_ref_path,
    reference_bfile = reference_bfile,
    dndscv_refdb = dndscv_refdb,
    epigenomic_track_dir = epigenomic_track_dir,
    epigenomic_tracks = epigenomic_tracks,
    graph_rds_path = paths$graph_rds_path,
    gg_nodes_path = paths$gg_nodes_path,
    gg_edges_path = paths$gg_edges_path,
    output_dir = output_dir,
    target_genes = target_genes,
    germline_args = germline_args,
    somatic_args = somatic_args,
    epigenomic_args = epigenomic_args,
    scored_graph_args = scored_graph_args,
    diffusion_args = diffusion_args,
    subgraph_args = subgraph_args,
    plot_args = plot_args,
    verbose = verbose
  )
}

#' Validate raw conseguiR inputs
#'
#' Validates raw GWAS, somatic, regulatory-reference, and epigenomic inputs
#' using the package's internal validation layer.
#'
#' @param gwas_sumstats GWAS summary statistics path or table.
#' @param somatic_maf Somatic MAF path or table.
#' @param reg_ref_path Regulatory-element reference path.
#' @param epigenomic_tracks Optional vector of bigWig paths.
#' @param epigenomic_track_dir Optional directory containing bigWig tracks.
#' @param verbose Logical scalar. If `TRUE`, show progress bars and stage
#'   messages when available.
#'
#' @details
#' `validate_inputs()` is intentionally a lightweight sanity check, not a full
#' dry-run of the downstream pipeline. The goal is to answer: "are these inputs
#' structured well enough for conseguiR to start?" rather than "can every
#' downstream tool finish successfully on this machine right now?".
#'
#' Input formatting rules:
#'
#' - `gwas_sumstats` can be either a file path or a data frame/data.table. When
#'   a table is supplied, it should contain columns that can be mapped to GWAS
#'   identifiers, chromosome, position, and p-value, for example
#'   `hm_variant_id`, `hm_chrom`, `hm_pos`, `p_value`, or the canonical columns
#'   `variant_id`, `chromosome`, `base_pair_location`, `p_value`.
#' - A practical GWAS minimum is therefore: one SNP identifier column, one
#'   chromosome column, one base-pair position column, and one p-value column.
#'   If harmonized columns such as `hm_rsid` and `hm_pos` are available, they
#'   are preferred.
#' - `somatic_maf` can be either a file path or a data frame/data.table. The
#'   table should contain sample identifier, chromosome, start, end, reference,
#'   and alternate allele columns in standard MAF-style names or names that the
#'   internal validator can harmonize.
#' - A practical somatic minimum is a MAF-like table with columns equivalent to
#'   `Tumor_Sample_Barcode`, `Chromosome`, `Start_Position`,
#'   `End_Position`, `Reference_Allele`, and `Tumor_Seq_Allele2`.
#' - `reg_ref_path` must be a path to a tab-delimited regulatory-element file.
#'   The current pipeline expects at least four columns corresponding to
#'   regulatory-element identifier, chromosome, start, and end. Extra columns
#'   are allowed and ignored by validation.
#' - `epigenomic_tracks` must be a character vector of full file paths to
#'   `.bw`/`.bigWig` files, for example
#'   `c(\"sample1.bw\", \"sample2.bw\", \"sample3.bw\")`.
#' - `epigenomic_track_dir` must be a single directory path containing bigWig
#'   files. Use either `epigenomic_tracks` or `epigenomic_track_dir`.
#'
#' @examples
#' validate_inputs()
#'
#' @return A validation bundle containing validated objects and config.
#' @export
validate_inputs <- function(
  gwas_sumstats = NULL,
  somatic_maf = NULL,
  reg_ref_path = NULL,
  epigenomic_tracks = NULL,
  epigenomic_track_dir = NULL,
  verbose = FALSE
) {
  .conseguiR_external_fun("validate_inputs")(
    gwas_sumstats = gwas_sumstats,
    somatic_maf = somatic_maf,
    reg_ref_path = reg_ref_path,
    epigenomic_tracks = epigenomic_tracks,
    epigenomic_track_dir = epigenomic_track_dir,
    verbose = verbose
  )
}

#' Initialize backend graph resources
#'
#' Builds the package's unscored backend graph resources when they are missing
#' and the required raw graph resources are available.
#'
#' @param backend_dir Optional backend cache directory. When `NULL`,
#'   `conseguiR` uses its default backend cache location.
#' @param build_gene_reg Logical scalar. If `TRUE`, ensure the gene-regulatory
#'   backend graph resources are available.
#' @param build_gene_gene Logical scalar. If `TRUE`, ensure the gene-gene
#'   backend graph resources are available.
#' @param force Logical scalar. If `TRUE`, rebuild or reseed backend resources
#'   even when cached outputs already exist.
#' @param strict Logical scalar. If `TRUE`, error when required backend inputs
#'   are unavailable. If `FALSE`, return a best-effort status object instead.
#' @param quiet Logical scalar. If `TRUE`, suppress non-essential backend
#'   initialization messages.
#' @param verbose Logical scalar. If `TRUE`, print backend initialization
#'   status messages.
#'
#' @examples
#' initialize_backend_graphs()
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
  quiet = FALSE,
  verbose = FALSE
) {
  .conseguiR_initialize_backend_graphs(
    backend_dir = backend_dir,
    build_gene_reg = build_gene_reg,
    build_gene_gene = build_gene_gene,
    force = force,
    strict = strict,
    quiet = if (isTRUE(verbose)) FALSE else quiet
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
#' If `gene_loc_path` is `NULL`, the package will try to use the backend-owned
#' gene location resource shipped with `conseguiR`.
#'
#' @inheritParams validate_inputs
#' @param gene_loc_path Optional gene location file for MAGMA step 1. When
#'   `NULL`, `conseguiR` uses its backend gene location resource.
#' @param reference_bfile PLINK reference prefix for MAGMA step 2.
#' @param output_prefix Output prefix for saved artifacts.
#' @param sample_size Fixed sample size for MAGMA.
#' @param sample_size_col Optional sample size column name.
#' @param magma_path Optional path to the MAGMA executable. When `NULL`,
#'   `conseguiR` searches in this order: `options(conseguiR.magma_path = ...)`,
#'   `Sys.getenv("CONSEGUIR_MAGMA_PATH")`, then `magma` on `PATH`.
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
#' @param verbose Logical scalar. If `TRUE`, show a progress bar and MAGMA step
#'   output.
#'
#' @details
#' Exact formatting:
#'
#' - `gwas_sumstats`: either a file path or a data frame/data.table with GWAS
#'   columns.
#' - `gene_loc_path`: a single path to a MAGMA-compatible gene location file.
#' - `reference_bfile`: the shared PLINK reference prefix, typed without file
#'   suffixes, for example `\"/path/to/g1000_eur/g1000_eur\"`.
#' - `output_prefix`: a single path prefix, not a directory. For example
#'   `\"results/germline_gene\"` will yield MAGMA files such as
#'   `results/germline_gene.genes.out`.
#' - `sample_size`: use this when the full GWAS should be analyzed with one
#'   fixed sample size.
#' - `sample_size_col`: use this when the GWAS table already contains a sample
#'   size column that MAGMA should read per row instead.
#' - `sample_size` and `sample_size_col` are alternatives. In a typical run you
#'   use one or the other, not both.
#' - `step1_args` and `step2_args`: named lists. For example:
#'
#' `step1_args = list(
#'   annotation_window = c(35, 10),
#'   filter_path = NULL,
#'   ignore_strand = FALSE,
#'   nonhuman = FALSE,
#'   extra_args = character()
#' )`
#'
#' `step2_args = list(
#'   gene_model = \"snp-wise=mean\",
#'   genes_only = TRUE,
#'   pval_use = c(\"SNP\", \"P\"),
#'   pval_duplicate = \"drop\",
#'   bfile_synonyms = NULL,
#'   bfile_synonym_dup = NULL,
#'   extra_args = character()
#' )`
#'
#' The most common wrapper-level MAGMA inputs usually look like this:
#'
#' - `reference_bfile = "/path/to/g1000_eur/g1000_eur"`
#' - `sample_size = 456348`
#' - `sample_size_col = "n_total"` when the GWAS already contains a usable
#'   sample-size column
#'
#' In plain language:
#'
#' - MAGMA step 1 (`step1_args`) controls how SNPs are assigned to features.
#'   The most important knob is usually `annotation_window`, which expands the
#'   feature boundaries upstream and downstream before MAGMA annotates SNPs.
#' - MAGMA step 2 (`step2_args`) controls how annotated SNP-level signal is
#'   collapsed into one feature-level statistic. The most important knobs are
#'   usually `gene_model`, `pval_use`, and `pval_duplicate`.
#' - `pval_use = c("SNP", "P")` tells MAGMA which columns in the p-value file
#'   correspond to SNP IDs and p-values.
#' - `pval_duplicate` follows MAGMA's own duplicate-SNP handling, for example
#'   `"drop"`, `"first"`, `"last"`, or `"error"`.
#' - `genes_only = TRUE` is a practical default when you want the standard
#'   feature-level MAGMA result rather than broader auxiliary outputs.
#' - `gene_model = "snp-wise=mean"` is a common starting point because it asks
#'   MAGMA to summarize SNP-level signal across the assigned feature rather than
#'   relying on a single top SNP.
#'
#' MAGMA manual:
#' \url{https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf}
#'
#' @examples
#' names(formals(run_germline_gene_scoring))
#'
#' @return A germline gene score bundle.
#' @export
run_germline_gene_scoring <- function(
  gwas_sumstats,
  gene_loc_path = NULL,
  reference_bfile,
  output_prefix = NULL,
  sample_size = NULL,
  sample_size_col = NULL,
  magma_path = NULL,
  magma_gwas_cache_prefix = NULL,
  reuse_existing_gwas_cache = TRUE,
  reuse_existing_annotation = FALSE,
  reuse_existing_analysis = FALSE,
  keep_intermediates = FALSE,
  step1_args = list(),
  step2_args = list(),
  verbose = FALSE
) {
  gene_loc_path <- gene_loc_path %||% .conseguiR_default_gene_loc_path()
  if (is.null(gene_loc_path)) {
    stop("No gene location resource was provided and no backend gene location resource could be found.")
  }

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
    step2_args = step2_args,
    verbose = verbose
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
#' If `reg_loc_path` is `NULL`, the package will try to use the backend-owned
#' regulatory location resource shipped with `conseguiR`.
#'
#' @inheritParams run_germline_gene_scoring
#' @param reg_loc_path Optional regulatory-element location file for MAGMA step
#'   1. When `NULL`, `conseguiR` uses its backend regulatory location resource.
#' @param verbose Logical scalar. If `TRUE`, show a progress bar and MAGMA step
#'   output.
#'
#' @details
#' `reg_loc_path` must point to a MAGMA-compatible regulatory-element location
#' file. The stage argument lists use the same format as
#' `run_germline_gene_scoring()`: `step1_args = list(...)` for annotation-stage
#' settings and `step2_args = list(...)` for gene-analysis-stage settings.
#'
#' In practice this function behaves like the gene-level MAGMA wrapper, except
#' that the "features" being annotated and scored are regulatory elements
#' rather than genes. A common pattern is to use a narrower
#' `annotation_window` here than for genes, because regulatory intervals are
#' already localized genomic objects.
#'
#' MAGMA manual:
#' \url{https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf}
#'
#' @examples
#' names(formals(run_germline_regulatory_scoring))
#'
#' @return A germline regulatory score bundle.
#' @export
run_germline_regulatory_scoring <- function(
  gwas_sumstats,
  reg_loc_path = NULL,
  reference_bfile,
  output_prefix = NULL,
  sample_size = NULL,
  sample_size_col = NULL,
  magma_path = NULL,
  magma_gwas_cache_prefix = NULL,
  reuse_existing_gwas_cache = TRUE,
  reuse_existing_annotation = FALSE,
  reuse_existing_analysis = FALSE,
  keep_intermediates = FALSE,
  step1_args = list(),
  step2_args = list(),
  verbose = FALSE
) {
  reg_loc_path <- reg_loc_path %||% .conseguiR_default_reg_loc_path()
  if (is.null(reg_loc_path)) {
    stop("No regulatory location resource was provided and no backend regulatory location resource could be found.")
  }

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
    step2_args = step2_args,
    verbose = verbose
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
#' @param verbose Logical scalar. If `TRUE`, show progress bars and MAGMA step
#'   output for both the gene and regulatory branches.
#'
#' @details
#' `prepare_germline_scores()` orchestrates two MAGMA runs, one for genes and
#' one for regulatory elements. The main tuning lists are:
#'
#' - `gene_step1_args`: MAGMA annotation-stage settings for the gene run.
#' - `gene_step2_args`: MAGMA gene-analysis-stage settings for the gene run.
#' - `reg_step1_args`: MAGMA annotation-stage settings for the regulatory run.
#' - `reg_step2_args`: MAGMA gene-analysis-stage settings for the regulatory
#'   run.
#' - `shared_args`: additional wrapper-level arguments passed to both runs.
#'
#' Typical entries for the step-specific lists include:
#'
#' - step 1 lists: `annotation_window`, `filter_path`, `ignore_strand`,
#'   `nonhuman`, `extra_args`
#' - step 2 lists: `gene_model`, `genes_only`, `pval_use`, `pval_duplicate`,
#'   `bfile_synonyms`, `bfile_synonym_dup`, `extra_args`
#'
#' Example:
#'
#' `prepare_germline_scores(
#'   gwas_sumstats = gwas_path,
#'   reference_bfile = \"/path/to/g1000_eur/g1000_eur\",
#'   gene_step1_args = list(annotation_window = c(35, 10)),
#'   gene_step2_args = list(gene_model = \"snp-wise=mean\", pval_use = c(\"SNP\", \"P\")),
#'   reg_step1_args = list(annotation_window = c(0, 0)),
#'   reg_step2_args = list(gene_model = \"snp-wise=mean\", pval_use = c(\"SNP\", \"P\"))
#' )`
#'
#' Practical decision rules:
#'
#' - use the same `reference_bfile` for both the gene and regulatory runs
#' - use `gene_sample_size` / `gene_sample_size_col` for the gene run and
#'   `reg_sample_size` / `reg_sample_size_col` for the regulatory run when the
#'   two branches need different MAGMA p-value inputs
#' - keep `gene_step1_args` and `reg_step1_args` separate when you want a wider
#'   annotation window for genes but a tighter one for already-localized
#'   regulatory elements
#' - keep `shared_args` for wrapper-level settings you truly want to send to
#'   both branches rather than duplicating them by hand
#'
#' MAGMA manual:
#' \url{https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf}
#'
#' @examples
#' names(formals(prepare_germline_scores))
#'
#' @return A germline score bundle with gene and regulatory score tables.
#' @export
prepare_germline_scores <- function(
  gwas_sumstats,
  reference_bfile,
  gene_output_prefix = NULL,
  reg_output_prefix = NULL,
  magma_gwas_cache_prefix = NULL,
  gene_sample_size = NULL,
  gene_sample_size_col = NULL,
  reg_sample_size = NULL,
  reg_sample_size_col = NULL,
  gene_step1_args = list(),
  gene_step2_args = list(),
  reg_step1_args = list(),
  reg_step2_args = list(),
  shared_args = list(),
  verbose = FALSE
) {
  .conseguiR_external_fun("prepare_germline_scores")(
    gwas_sumstats = gwas_sumstats,
    reference_bfile = reference_bfile,
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
    shared_args = shared_args,
    verbose = verbose
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
#' @param verbose Logical scalar. If `TRUE`, show stage messages and dndscv
#'   output.
#'
#' @details
#' Exact formatting:
#'
#' - `maf`: either a file path or a data frame/data.table containing somatic
#'   mutation records.
#' - `refdb`: a single path to a dndscv reference database `.rda` file.
#' - `cv`: optional covariate structure passed through to `dndscv()`. In most
#'   user workflows this is `NULL`, but advanced users can pass a covariate
#'   matrix or object that matches dndscv's expectations.
#' - `dndscv_args`: a named list of additional dndscv arguments, for example
#'   `list(sm = \"192r_3w\", kc = \"cgc81\")`.
#' - `maf` should already be MAF-like before it reaches this function. In
#'   practice that means sample, chromosome, start, end, reference allele, and
#'   alternate allele information should already be available under standard
#'   MAF-style column names or names the package validator can harmonize.
#'
#' In plain language:
#'
#' - `max_muts_per_gene_per_sample` caps how many mutations from one sample can
#'   contribute to one gene before dndscv filters them down.
#' - `max_coding_muts_per_sample` caps extremely hypermutated samples in the
#'   coding analysis.
#' - the most common advanced dndscv passthrough is `sm`, which controls the
#'   substitution model, for example `\"192r_3w\"`.
#' - if you supply `cv`, it should already be formatted the way `dndscv()`
#'   expects it. `conseguiR` does not reshape arbitrary covariate tables into a
#'   dndscv-ready object for you.
#' - `refdb` is not just any annotation file. It must be a dndscv-compatible
#'   reference `.rda` resource built for the same genome build as the MAF you
#'   are scoring.
#'
#' Minimal examples:
#'
#' `dndscv_args = list(
#'   sm = \"192r_3w\",
#'   kc = \"cgc81\"
#' )`
#'
#' `cv = NULL`
#'
#' If you do use `cv`, think of it as an advanced dndscv-native object rather
#' than a generic spreadsheet of covariates. The package passes it through; it
#' does not convert a plain data frame into the structure expected by dndscv.
#'
#' dndscv documentation:
#' \url{https://rdrr.io/github/im3sanger/dndscv/man/dndscv.html}
#'
#' @examples
#' names(formals(run_somatic_gene_scoring))
#'
#' @return A somatic gene score bundle.
#' @export
run_somatic_gene_scoring <- function(
  maf,
  refdb,
  output_path = NULL,
  cv = NULL,
  max_muts_per_gene_per_sample = 6L,
  max_coding_muts_per_sample = 5000L,
  dndscv_args = list(),
  verbose = FALSE
) {
  .conseguiR_external_fun("run_somatic_gene_scoring")(
    maf = maf,
    refdb = refdb,
    output_path = output_path,
    cv = cv,
    max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,
    max_coding_muts_per_sample = max_coding_muts_per_sample,
    dndscv_args = dndscv_args,
    verbose = verbose
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
#' @param verbose Logical scalar. If `TRUE`, show stage messages and fishHook
#'   output.
#'
#' @details
#' Exact formatting:
#'
#' - `reg_ref_path`: a path to the regulatory-element reference file.
#' - `eligible_gr`: a `GRanges` object or `NULL`.
#' - `fishhook_covariates`: the covariate specification object expected by the
#'   internal fishHook runner, or `NULL`.
#' - `fishhook_covariate_data`: a data frame/data.table containing one row per
#'   regulatory element and the covariate columns needed by fishHook. In
#'   practice this should include a regulatory-element identifier column plus
#'   the numeric or categorical covariates you want fishHook to model.
#' - `idcol`: the sample identifier column used by the regulatory somatic input
#'   handed to fishHook. It should match the sample identifier field in the MAF
#'   after harmonization.
#' - `fishhook_args`: a named list of additional fishHook arguments.
#' - `reg_ref_path` should resolve to the same regulatory universe you want to
#'   score. In practice this means the regulatory-element IDs in
#'   `fishhook_covariate_data` should correspond to the IDs in that reference.
#'
#' In plain language:
#'
#' - `eligible_gr` defines the territory in which fishHook is allowed to place
#'   mutations when estimating the background model.
#' - `fishhook_covariate_data` is where users typically supply replication,
#'   accessibility, mappability, or GC-like covariates if they have them.
#' - if you do not have a custom territory or covariates yet, `eligible_gr =
#'   NULL` and `fishhook_covariate_data = NULL` are reasonable starting points.
#' - `idcol` should match the sample identifier column name used in the somatic
#'   mutation table after harmonization. In most MAF-based workflows the
#'   default `Tumor_Sample_Barcode` is the right choice.
#'
#' Minimal covariate-data example:
#'
#' `covariate_dt = data.frame(
#'   reg_elem_id = c(\"GH01J000001\", \"GH01J000002\"),
#'   accessibility = c(1.2, 0.4),
#'   replication_timing = c(0.7, -0.1),
#'   gc_content = c(0.44, 0.51)
#' )`
#'
#' `fishhook_covariates = NULL`
#'
#' `fishhook_args = list()`
#'
#' A good first-pass mental model is:
#'
#' - `eligible_gr` says where mutations are allowed to occur in the background
#'   model
#' - `fishhook_covariate_data` says which feature-level covariates explain that
#'   background
#' - `fishhook_covariates` is only needed when you already have a fishHook-side
#'   covariate specification you want to pass through directly
#'
#' fishHook documentation:
#' \url{https://rdrr.io/github/mskilab/fish.hook/man/FishHook.html}
#'
#' fishHook tutorial:
#' \url{https://mskilab.com/fishHook/tutorial.html}
#'
#' @examples
#' names(formals(run_somatic_regulatory_scoring))
#'
#' @return A somatic regulatory score bundle.
#' @export
run_somatic_regulatory_scoring <- function(
  maf,
  reg_ref_path,
  output_path = NULL,
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  idcol = "Tumor_Sample_Barcode",
  fishhook_args = list(),
  verbose = FALSE
) {
  .conseguiR_external_fun("run_somatic_regulatory_scoring")(
    maf = maf,
    reg_ref_path = reg_ref_path,
    output_path = output_path,
    eligible_gr = eligible_gr,
    fishhook_covariates = fishhook_covariates,
    fishhook_covariate_data = fishhook_covariate_data,
    idcol = idcol,
    fishhook_args = fishhook_args,
    verbose = verbose
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
#' @param verbose Logical scalar. If `TRUE`, show progress bars and tool output
#'   for both the dndscv and fishHook branches.
#'
#' @details
#' `prepare_somatic_scores()` combines two independent somatic branches:
#'
#' - gene-level dndscv controls: `gene_cv`,
#'   `gene_max_muts_per_gene_per_sample`,
#'   `gene_max_coding_muts_per_sample`, `dndscv_args`
#' - regulatory-level fishHook controls: `eligible_gr`,
#'   `fishhook_covariates`, `fishhook_covariate_data`, `fishhook_idcol`,
#'   `fishhook_args`
#'
#' A good mental model is:
#'
#' - `prepare_somatic_scores()` does not fit one joint model
#' - it runs a gene-oriented dndscv branch and a regulatory-element-oriented
#'   fishHook branch separately
#' - then it returns both score tables together as one somatic bundle
#'
#' Formatting notes for covariate-bearing inputs:
#'
#' - `gene_cv` should already be formatted for direct use by `dndscv()`
#' - `fishhook_covariate_data` should contain one row per regulatory element
#'   plus the covariate columns you want fishHook to use
#' - `fishhook_covariates` should already be a fishHook-ready specification if
#'   you provide it
#' - `fishhook_idcol` should match the sample identifier field used by the MAF
#'   after harmonization, typically `Tumor_Sample_Barcode`
#'
#' In practice, a simple starting configuration is often:
#'
#' - `gene_cv = NULL`
#' - `dndscv_args = list(sm = "192r_3w")`
#' - `eligible_gr = NULL`
#' - `fishhook_covariate_data = NULL` or a one-row-per-reg-element table
#' - `fishhook_args = list()`
#'
#' Example:
#'
#' `prepare_somatic_scores(
#'   maf = maf_path,
#'   refdb = dndscv_refdb,
#'   reg_ref_path = reg_ref_path,
#'   dndscv_args = list(sm = \"192r_3w\"),
#'   fishhook_args = list()
#' )`
#'
#' dndscv documentation:
#' \url{https://rdrr.io/github/im3sanger/dndscv/man/dndscv.html}
#'
#' fishHook documentation:
#' \url{https://rdrr.io/github/mskilab/fish.hook/man/FishHook.html}
#'
#' fishHook tutorial:
#' \url{https://mskilab.com/fishHook/tutorial.html}
#'
#' @examples
#' names(formals(prepare_somatic_scores))
#'
#' @return A somatic score bundle with gene and regulatory score tables.
#' @export
prepare_somatic_scores <- function(
  maf,
  refdb,
  reg_ref_path,
  gene_output_path = NULL,
  reg_output_path = NULL,
  gene_cv = NULL,
  gene_max_muts_per_gene_per_sample = 6L,
  gene_max_coding_muts_per_sample = 5000L,
  dndscv_args = list(),
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  fishhook_idcol = "Tumor_Sample_Barcode",
  fishhook_args = list(),
  verbose = FALSE
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
    fishhook_args = fishhook_args,
    verbose = verbose
  )
}

#' Prepare regulatory epigenomic scores
#'
#' @param reg_ref_path Regulatory-element reference path.
#' @param track_dir Optional directory containing bigWig tracks.
#' @param bw_files Optional explicit vector of bigWig paths.
#' @param output_path Optional output path for saved scores.
#' @param min_tracks Minimum number of tracks required.
#' @param drop_mhc Whether to exclude extended MHC elements.
#' @param transform Signal transformation method.
#' @param return_diagnostics Whether to return the diagnostic table and raw
#'   matrix alongside z scores.
#' @param summary_fun Summary function applied to each bigWig over each
#'   regulatory element.
#' @param verbose Logical scalar. If `TRUE`, show a progress bar and stage
#'   messages during epigenomic scoring.
#'
#' @details
#' Exact formatting:
#'
#' - `track_dir`: a single directory path containing bigWig files.
#' - `bw_files`: a character vector of bigWig paths. Supply at least three
#'   tracks, for example
#'   `c(\"sample1.bw\", \"sample2.bw\", \"sample3.bw\")`.
#' - `summary_fun`: a function object such as `mean` or `max`.
#'
#' The current epigenomic score is not a generic "activity" score. Instead,
#' conseguiR quantifies each regulatory element across the supplied tracks,
#' summarizes each track over each element with `summary_fun`, and then scores
#' elements by how variable they are across tracks. The intuition is that
#' highly variable elements are more context-specific, whereas uniformly active
#' elements are more housekeeping-like.
#'
#' In practice:
#'
#' - `transform = "log1p"` is a reasonable default for signal stabilization
#' - `summary_fun = mean` is a reasonable default for average signal over each
#'   regulatory element
#' - at least three tracks are required so that cross-track variability is
#'   meaningful
#'
#' @examples
#' names(formals(prepare_epigenomic_scores))
#'
#' @return An epigenomic score bundle.
#' @export
prepare_epigenomic_scores <- function(
  reg_ref_path,
  track_dir = NULL,
  bw_files = NULL,
  output_path = NULL,
  min_tracks = 3L,
  drop_mhc = TRUE,
  transform = "log1p",
  return_diagnostics = TRUE,
  summary_fun = mean,
  verbose = FALSE
) {
  .conseguiR_external_fun("prepare_epigenomic_scores")(
    reg_ref_path = reg_ref_path,
    track_dir = track_dir,
    bw_files = bw_files,
    output_path = output_path,
    min_tracks = min_tracks,
    drop_mhc = drop_mhc,
    transform = transform,
    return_diagnostics = return_diagnostics,
    summary_fun = summary_fun,
    verbose = verbose
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
#' @param verbose Logical scalar. If `TRUE`, show stage messages while building
#'   the scored graph.
#'
#' @details
#' Score-table formatting:
#'
#' - gene-level score tables must contain a gene identifier column and a `zstat`
#'   column
#' - regulatory-level score tables must contain a regulatory-element identifier
#'   column and a `zstat` column
#' - the safest explicit formats are:
#'   - gene tables with columns like `gene_id`, `zstat`
#'   - regulatory tables with columns like `reg_elem_id`, `zstat`
#'
#' You can either pass full score bundles from earlier stages or pass explicit
#' score tables directly through `gene_germline_scores`,
#' `reg_germline_scores`, `gene_somatic_scores`, `reg_somatic_scores`, and
#' `reg_epigenomic_scores`.
#'
#' @examples
#' names(formals(build_scored_gene_reg_graph))
#'
#' toy_graph <- igraph::graph_from_data_frame(
#'   d = data.frame(
#'     from = "GH01J000001",
#'     to = "TP53",
#'     confidence = 0.9
#'   ),
#'   vertices = data.frame(
#'     name = c("TP53", "GH01J000001"),
#'     node_id = c("TP53", "GH01J000001"),
#'     node_type = c("gene", "reg")
#'   ),
#'   directed = TRUE
#' )
#'
#' scored_graph <- build_scored_gene_reg_graph(
#'   graph = toy_graph,
#'   graph_rds_path = tempfile(fileext = ".rds"),
#'   gene_germline_scores = data.frame(gene_id = "TP53", zstat = 2),
#'   reg_germline_scores = data.frame(reg_elem_id = "GH01J000001", zstat = 1.2),
#'   gene_somatic_scores = data.frame(gene_id = "TP53", zstat = -1.5),
#'   reg_somatic_scores = data.frame(reg_elem_id = "GH01J000001", zstat = 0.3),
#'   reg_epigenomic_scores = data.frame(reg_elem_id = "GH01J000001", zstat = 2.4),
#'   save_outputs = FALSE
#' )
#'
#' names(scored_graph$objects)
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
  save_outputs = FALSE,
  verbose = FALSE
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
    save_outputs = save_outputs,
    verbose = verbose
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
#' @param python_path Deprecated and ignored. Python is managed internally via
#'   `basilisk`.
#' @param verbose Logical scalar. If `TRUE`, show stage messages while running
#'   diffusion.
#'
#' @details
#' Exact formatting:
#'
#' - `scored_graph`: the bundle returned by `build_scored_gene_reg_graph()`
#' - `nodes_path` and `edges_path`: paths to the scored node and edge tables if
#'   you are not passing `scored_graph`
#' - `top_k`: a positive integer
#' - `confidence_power`, `beta_germline`, `beta_somatic`,
#'   `beta_epigenomic`, `reg_signal_clip`: numeric scalars
#' - `positive_only`: a single logical value
#'
#' In plain language:
#'
#' - `top_k` controls how many regulatory neighbors contribute most strongly to
#'   each gene's update step
#' - `beta_germline`, `beta_somatic`, and `beta_epigenomic` control how much
#'   each modality contributes to the regulatory signal that is propagated
#' - `confidence_power` upweights or downweights high-confidence regulatory
#'   links relative to weaker ones
#' - `positive_only = TRUE` suppresses negative regulatory signal before
#'   diffusion
#' - `reg_signal_clip` caps extreme regulatory signal before propagation so one
#'   extreme feature does not dominate the update
#'
#' @examples
#' names(formals(run_gene_reg_diffusion))
#'
#' @return A diffusion bundle containing full and top-gene diffusion tables.
#' @export
run_gene_reg_diffusion <- function(
  scored_graph = NULL,
  nodes_path = NULL,
  edges_path = NULL,
  output_dir = NULL,
  output_stem = NULL,
  top_k = 3L,
  confidence_power = 2.0,
  beta_germline = 0.5,
  beta_somatic = 0.5,
  beta_epigenomic = 0.7,
  positive_only = FALSE,
  reg_signal_clip = 5.0,
  top_n_to_save = 50L,
  python_path = NULL,
  verbose = FALSE
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
    python_path = python_path,
    verbose = verbose
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
#' @param python_path Deprecated and ignored. Python is managed internally via
#'   `basilisk`.
#' @param verbose Logical scalar. If `TRUE`, show stage messages while calling
#'   the selected subgraph.
#'
#' @details
#' Exact formatting:
#'
#' - `diffusion`: the bundle returned by `run_gene_reg_diffusion()`
#' - `diffusion_path`: path to the full diffusion results table if `diffusion`
#'   is not supplied
#' - `target_genes`, `candidate_pool_size`, `max_edges_in_model`,
#'   `max_time_seconds`, `num_workers`, `random_seed`: integer scalars
#' - `node_prize_weight`, `edge_conf_weight`, `edge_cost_weight`: numeric
#'   scalars
#' - `prize_column`, `confidence_column`, `edge_cost_column`: single column
#'   names present in the input tables
#'
#' Practical column examples:
#'
#' - `prize_column = "post_norm"` means optimize using the post-diffusion gene
#'   score column
#' - `confidence_column = "confidence"` means use the confidence column from
#'   the backend gene-gene edge table
#' - `edge_cost_column = "weight"` means use the edge-cost/penalty column from
#'   the backend gene-gene edge table
#'
#' In plain language:
#'
#' - `target_genes` is the size of the final selected subgraph you want back
#' - `candidate_pool_size` is the number of top candidate genes that are handed
#'   to the optimization stage before the final smaller subgraph is chosen
#' - `node_prize_weight` rewards high-scoring genes from diffusion
#' - `edge_conf_weight` rewards keeping confident gene-gene edges
#' - `edge_cost_weight` penalizes expensive edges
#' - `max_edges_in_model` and `max_time_seconds` are practical solver-budget
#'   controls for difficult graphs
#'
#' @examples
#' names(formals(call_selected_subgraph))
#'
#' @return A selected subgraph bundle.
#' @export
call_selected_subgraph <- function(
  diffusion = NULL,
  diffusion_path = NULL,
  gg_nodes_path = NULL,
  gg_edges_path = NULL,
  output_dir = NULL,
  output_stem = NULL,
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
  python_path = NULL,
  verbose = FALSE
) {
  backend_paths <- .conseguiR_default_gene_gene_paths(
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path
  )

  .conseguiR_call_external("call_selected_subgraph", list(
    diffusion = diffusion,
    diffusion_path = diffusion_path,
    gg_nodes_path = backend_paths$gg_nodes_path,
    gg_edges_path = backend_paths$gg_edges_path,
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
    python_path = python_path,
    verbose = verbose
  ))
}

#' Plot score tables from scoring or diffusion outputs
#'
#' Creates a rank plot for one-tailed outputs and a volcano plot for two-tailed
#' outputs.
#'
#' @param scores Optional bundle returned by a scoring wrapper or
#'   `run_gene_reg_diffusion()`.
#' @param table Optional explicit table.
#' @param which Optional bundle component name such as `gene_scores`,
#'   `reg_scores`, `all_genes`, or `top_genes`.
#' @param plot_file_path Optional output path for the saved figure.
#' @param test_tail One of `auto`, `one_tailed`, or `two_tailed`.
#' @param feature_column Optional explicit feature-label column.
#' @param z_column Z-score column name.
#' @param p_value_column Optional explicit p-value column for volcano plots.
#' @param label_features Optional character vector of features to label.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages while building
#'   the plot.
#'
#' @details
#' Exact formatting:
#' - `scores`: a bundle returned by `run_germline_gene_scoring()`,
#'   `run_germline_regulatory_scoring()`, `prepare_germline_scores()`,
#'   `run_somatic_gene_scoring()`, `run_somatic_regulatory_scoring()`,
#'   `prepare_somatic_scores()`, `prepare_epigenomic_scores()`, or
#'   `run_gene_reg_diffusion()`
#' - `table`: a data frame/data.table containing at least a feature column and
#'   a `zstat`-like numeric column
#' - `which`: when `scores` contains multiple tables, use values like
#'   `"gene_scores"`, `"reg_scores"`, `"all_genes"`, or `"top_genes"`
#' - `test_tail = "two_tailed"` expects a usable p-value column in the table
#' - `plot_file_path`: a single file path such as `"scores_plot.pdf"`
#'
#' @examples
#' plot_scores(
#'   table = data.frame(feature_id = c("A", "B"), zstat = c(1, -1)),
#'   feature_column = "feature_id",
#'   z_column = "zstat",
#'   save_plot = FALSE
#' )
#'
#' @return A plot bundle containing the ggplot object and plotting data.
#' @export
plot_scores <- function(
  scores = NULL,
  table = NULL,
  which = NULL,
  plot_file_path = NULL,
  test_tail = "auto",
  feature_column = NULL,
  z_column = "zstat",
  p_value_column = NULL,
  label_features = NULL,
  title = "conseguiR Scores",
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_scores")(
    scores = scores,
    table = table,
    which = which,
    plot_file_path = plot_file_path,
    test_tail = test_tail,
    feature_column = feature_column,
    z_column = z_column,
    p_value_column = p_value_column,
    label_features = label_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot germline gene signal before or after diffusion
#'
#' @param germline_scores Optional germline score bundle.
#' @param scored_graph Optional scored gene-reg graph bundle.
#' @param diffusion Optional diffusion bundle.
#' @param nodes_path Optional explicit scored node-table path.
#' @param diffusion_path Optional explicit diffusion table path.
#' @param stage One of `"pre"` or `"post"`.
#' @param plot_file_path Optional output path for the saved figure.
#' @param label_features Optional gene symbols to highlight and label.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @examples
#' bundle <- list(gene_scores = data.frame(gene_id = c("A", "B"), zstat = c(1, -1)))
#' plot_germline_gene_scores(germline_scores = bundle, stage = "pre", save_plot = FALSE)
#'
#' @return A plot bundle containing the ggplot object and plotting data.
#' @export
plot_germline_gene_scores <- function(
  germline_scores = NULL,
  scored_graph = NULL,
  diffusion = NULL,
  nodes_path = NULL,
  diffusion_path = NULL,
  stage = c("pre", "post"),
  plot_file_path = NULL,
  label_features = NULL,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_germline_gene_scores")(
    germline_scores = germline_scores,
    scored_graph = scored_graph,
    diffusion = diffusion,
    nodes_path = nodes_path,
    diffusion_path = diffusion_path,
    stage = stage,
    plot_file_path = plot_file_path,
    label_features = label_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot germline regulatory signal
#'
#' @param germline_scores Optional germline score bundle.
#' @param scored_graph Optional scored gene-reg graph bundle.
#' @param nodes_path Optional explicit scored node-table path.
#' @param plot_file_path Optional output path for the saved figure.
#' @param label_features Optional gene symbols to highlight and label.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @examples
#' bundle <- list(reg_scores = data.frame(reg_elem_id = c("r1", "r2"), zstat = c(1, -1)))
#' plot_germline_reg_scores(germline_scores = bundle, save_plot = FALSE)
#'
#' @return A plot bundle containing the ggplot object and plotting data.
#' @export
plot_germline_reg_scores <- function(
  germline_scores = NULL,
  scored_graph = NULL,
  nodes_path = NULL,
  plot_file_path = NULL,
  label_features = NULL,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_germline_reg_scores")(
    germline_scores = germline_scores,
    scored_graph = scored_graph,
    nodes_path = nodes_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot somatic gene signal before or after diffusion
#'
#' @param somatic_scores Optional somatic score bundle.
#' @param scored_graph Optional scored gene-reg graph bundle.
#' @param diffusion Optional diffusion bundle.
#' @param nodes_path Optional explicit scored node-table path.
#' @param diffusion_path Optional explicit diffusion table path.
#' @param stage One of `"pre"` or `"post"`.
#' @param plot_file_path Optional output path for the saved figure.
#' @param label_features Optional gene symbols to highlight and label.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @examples
#' bundle <- list(
#'   gene_scores = data.frame(gene_id = c("A", "B"), zstat = c(2, -2), p_value = c(0.01, 0.02))
#' )
#' plot_somatic_gene_scores(somatic_scores = bundle, stage = "pre", save_plot = FALSE)
#'
#' @return A plot bundle containing the ggplot object and plotting data.
#' @export
plot_somatic_gene_scores <- function(
  somatic_scores = NULL,
  scored_graph = NULL,
  diffusion = NULL,
  nodes_path = NULL,
  diffusion_path = NULL,
  stage = c("pre", "post"),
  plot_file_path = NULL,
  label_features = NULL,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_somatic_gene_scores")(
    somatic_scores = somatic_scores,
    scored_graph = scored_graph,
    diffusion = diffusion,
    nodes_path = nodes_path,
    diffusion_path = diffusion_path,
    stage = stage,
    plot_file_path = plot_file_path,
    label_features = label_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot somatic regulatory signal
#'
#' @param somatic_scores Optional somatic score bundle.
#' @param scored_graph Optional scored gene-reg graph bundle.
#' @param nodes_path Optional explicit scored node-table path.
#' @param plot_file_path Optional output path for the saved figure.
#' @param label_features Optional gene symbols to highlight and label.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @examples
#' bundle <- list(
#'   reg_scores = data.frame(reg_elem_id = c("r1", "r2"), zstat = c(2, -2), p_value = c(0.01, 0.02))
#' )
#' plot_somatic_reg_scores(somatic_scores = bundle, save_plot = FALSE)
#'
#' @return A plot bundle containing the ggplot object and plotting data.
#' @export
plot_somatic_reg_scores <- function(
  somatic_scores = NULL,
  scored_graph = NULL,
  nodes_path = NULL,
  plot_file_path = NULL,
  label_features = NULL,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_somatic_reg_scores")(
    somatic_scores = somatic_scores,
    scored_graph = scored_graph,
    nodes_path = nodes_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot post-diffusion epigenomic gene signal
#'
#' @param diffusion Optional diffusion bundle.
#' @param diffusion_path Optional explicit diffusion table path.
#' @param plot_file_path Optional output path for the saved figure.
#' @param label_features Optional gene symbols to highlight and label.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @examples
#' diffusion_bundle <- list(
#'   all_genes = data.frame(
#'     gene_name = c("A", "B"),
#'     post_epigenomic = c(1.5, 0.5),
#'     post_norm = c(2, 1)
#'   )
#' )
#' plot_epigenomic_gene_scores(diffusion = diffusion_bundle, save_plot = FALSE)
#'
#' @return A plot bundle containing the ggplot object and plotting data.
#' @export
plot_epigenomic_gene_scores <- function(
  diffusion = NULL,
  diffusion_path = NULL,
  plot_file_path = NULL,
  label_features = NULL,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_epigenomic_gene_scores")(
    diffusion = diffusion,
    diffusion_path = diffusion_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot epigenomic regulatory signal
#'
#' @param epigenomic_scores Optional epigenomic score bundle.
#' @param scored_graph Optional scored gene-reg graph bundle.
#' @param nodes_path Optional explicit scored node-table path.
#' @param plot_file_path Optional output path for the saved figure.
#' @param label_features Optional features to highlight and label.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @examples
#' bundle <- list(reg_scores = data.frame(reg_elem_id = c("r1", "r2"), zstat = c(1, -1)))
#' plot_epigenomic_reg_scores(epigenomic_scores = bundle, save_plot = FALSE)
#'
#' @return A plot bundle containing the ggplot object and plotting data.
#' @export
plot_epigenomic_reg_scores <- function(
  epigenomic_scores = NULL,
  scored_graph = NULL,
  nodes_path = NULL,
  plot_file_path = NULL,
  label_features = NULL,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_epigenomic_reg_scores")(
    epigenomic_scores = epigenomic_scores,
    scored_graph = scored_graph,
    nodes_path = nodes_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot a locus-centered multimodal context panel
#'
#' Creates a genome-track-style panel over a user-specified locus. Genes are
#' drawn as rectangles, regulatory elements as circles, and regulatory links as
#' arcs. The panel combines post-diffusion gene scores with locus-matched
#' regulatory scores from the scored graph.
#'
#' @param chromosome Locus chromosome, for example `"8"` or `"chr8"`.
#' @param start Locus start coordinate.
#' @param end Locus end coordinate.
#' @param postdiff_gene_reg_graph Optional merged post-diffusion gene-reg graph
#'   object or `.rds` path containing node and edge tables.
#' @param scored_graph Optional scored gene-reg graph bundle.
#' @param diffusion Optional diffusion bundle.
#' @param selected_subgraph Optional selected-subgraph bundle used to highlight
#'   genes and filter displayed regulatory arcs.
#' @param nodes_path Optional explicit scored node-table path.
#' @param edges_path Optional explicit scored edge-table path.
#' @param diffusion_path Optional explicit diffusion table path.
#' @param selected_nodes_path Optional explicit selected-subgraph node-table
#'   path.
#' @param label_features Optional gene symbols to label. For regulatory-element
#'   labels, the plot uses the top-scoring regulatory element for each requested
#'   gene label within the locus.
#' @param gwas_sumstats Optional GWAS summary statistics path or table used for
#'   locus SNP labeling.
#' @param label_top_gwas_snp Logical scalar. If `TRUE`, label the top GWAS SNP
#'   inside the top germline regulatory element in the locus.
#' @param rsid_pmid Optional cached rsID-to-PMID evidence table with at least
#'   `rsid` and `pmid` columns.
#' @param label_top_lit_snps Integer count of literature-backed SNPs to label.
#'   SNPs are restricted to regulatory elements in the locus. If no
#'   literature-backed SNPs survive the lookup/filtering steps, `conseguiR`
#'   falls back to top GWAS SNPs from the top germline regulatory elements.
#' @param pmid_query Deprecated. This argument is currently ignored; locus SNP
#'   labeling now uses dbSNP citation support directly.
#' @param pmid_page_size Maximum number of PMIDs per rsID retained from the
#'   dbSNP citation lookup.
#' @param plot_file_path Optional output path for the saved figure.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @details
#' Exact formatting:
#' - `chromosome`: a single chromosome string such as `"8"` or `"chr8"`
#' - `start`, `end`: scalar genomic coordinates defining the locus
#' - `postdiff_gene_reg_graph`: an optional merged post-diffusion graph object
#'   containing node and edge tables
#' - `scored_graph`: the bundle returned by `build_scored_gene_reg_graph()`
#' - `diffusion`: the bundle returned by `run_gene_reg_diffusion()`
#' - `selected_subgraph`: the bundle returned by `call_selected_subgraph()`
#' - `label_features`: a character vector of gene symbols such as
#'   `c("MYC", "BCL2", "BCL6")`
#' - `rsid_pmid`: an optional cached literature-support table with at least
#'   `rsid` and `pmid` columns
#' - `label_top_lit_snps`: number of literature-backed SNPs to label before
#'   falling back to top GWAS SNPs in the top germline regulatory elements
#'
#' Track semantics:
#' - the top three rows show regulatory-element somatic, epigenomic, and
#'   germline z-scores
#' - the `Reg elements` row shows regulatory elements colored by their combined
#'   pre-diffusion norm
#' - the bottom gene row shows post-diffusion `conseguiR` scores for genes in
#'   the locus
#' - thin black diagonal curves connect regulatory elements to genes
#' - locus SNP labels prefer literature-backed SNPs when available and
#'   otherwise fall back to top GWAS SNPs in the top germline regulatory
#'   elements
#'
#' @examples
#' names(formals(plot_locus_context))
#'
#' @return A plot bundle containing the ggplot object and locus plotting data.
#' @export
plot_locus_context <- function(
  chromosome,
  start,
  end,
  postdiff_gene_reg_graph = NULL,
  scored_graph = NULL,
  diffusion = NULL,
  selected_subgraph = NULL,
  nodes_path = NULL,
  edges_path = NULL,
  diffusion_path = NULL,
  selected_nodes_path = NULL,
  label_features = NULL,
  gwas_sumstats = NULL,
  label_top_gwas_snp = FALSE,
  rsid_pmid = NULL,
  label_top_lit_snps = 0L,
  pmid_query = NULL,
  pmid_page_size = 1000L,
  plot_file_path = NULL,
  title = NULL,
  width = 14,
  height = 9,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_call_external("plot_locus_context", list(
    chromosome = chromosome,
    start = start,
    end = end,
    postdiff_gene_reg_graph = postdiff_gene_reg_graph,
    scored_graph = scored_graph,
    diffusion = diffusion,
    selected_subgraph = selected_subgraph,
    nodes_path = nodes_path,
    edges_path = edges_path,
    diffusion_path = diffusion_path,
    selected_nodes_path = selected_nodes_path,
    label_features = label_features,
    gwas_sumstats = gwas_sumstats,
    label_top_gwas_snp = label_top_gwas_snp,
    rsid_pmid = rsid_pmid,
    label_top_lit_snps = label_top_lit_snps,
    pmid_query = pmid_query,
    pmid_page_size = pmid_page_size,
    plot_file_path = plot_file_path,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  ))
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
#' @param verbose Logical scalar. If `TRUE`, show stage messages while building
#'   the plot bundle.
#'
#' @details
#' Exact formatting:
#'
#' - `selected_subgraph`: the bundle returned by `call_selected_subgraph()`
#' - `nodes`, `edges`, `summary`: data frames/data.tables containing the
#'   selected subgraph outputs
#' - `nodes_path`, `edges_path`, `summary_path`: file paths to those tables if
#'   the in-memory objects are not supplied
#' - `plot_file_path`: a single file path ending in `.pdf`, `.png`, or another
#'   graphics device extension supported by the plotting helper
#' - `layout`: a single layout keyword such as `\"fr\"`
#'
#' @examples
#' names(formals(plot_selected_subgraph))
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
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
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
    save_plot = save_plot,
    verbose = verbose
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
#' @param graph_rds_path Backend no-score gene-reg graph path.
#' @param gg_nodes_path Gene-gene node table path.
#' @param gg_edges_path Gene-gene edge table path.
#' @param output_dir Optional output directory for saved pipeline artifacts.
#'   When `NULL`, `run_conseguiR()` runs in object-first mode and returns the
#'   stage bundles without treating disk output as the default interface.
#' @param target_genes Requested selected-subgraph size.
#' @param germline_args Named list of overrides passed to
#'   `prepare_germline_scores()`. This list may contain both gene- and
#'   regulatory-run settings, for example `gene_sample_size`,
#'   `reg_sample_size`, `gene_step1_args`, `gene_step2_args`,
#'   `reg_step1_args`, `reg_step2_args`, and `shared_args`.
#' @param somatic_args Named list of overrides passed to
#'   `prepare_somatic_scores()`. This list may contain both dndscv and fishHook
#'   settings, for example `gene_cv`,
#'   `gene_max_muts_per_gene_per_sample`,
#'   `gene_max_coding_muts_per_sample`, `dndscv_args`, `eligible_gr`,
#'   `fishhook_covariates`, `fishhook_covariate_data`, `fishhook_idcol`, and
#'   `fishhook_args`.
#' @param epigenomic_args Named list of overrides passed to
#'   `prepare_epigenomic_scores()`. Typical entries include `track_dir`,
#'   `bw_files`, `min_tracks`, `drop_mhc`, `transform`,
#'   `return_diagnostics`, and `summary_fun`.
#' @param scored_graph_args Named list of overrides passed to
#'   `build_scored_gene_reg_graph()`.
#' @param diffusion_args Named list of overrides passed to
#'   `run_gene_reg_diffusion()`.
#' @param subgraph_args Named list of overrides passed to
#'   `call_selected_subgraph()`.
#' @param plot_args Named list of overrides passed to
#'   `plot_selected_subgraph()`.
#' @param verbose Logical scalar. If `TRUE`, show progress bars and stage
#'   output across the pipeline.
#'
#' @details
#' `run_conseguiR()` is a thin orchestration wrapper. Most tuning is passed
#' through stage-specific named lists:
#'
#' - `germline_args` controls the two MAGMA branches inside
#'   `prepare_germline_scores()`
#' - `somatic_args` controls the dndscv and fishHook branches inside
#'   `prepare_somatic_scores()`
#' - `epigenomic_args` controls `prepare_epigenomic_scores()`
#' - `scored_graph_args`, `diffusion_args`, `subgraph_args`, and `plot_args`
#'   are forwarded to their corresponding downstream stages
#'
#' Formatting notes:
#'
#' - each `*_args` input should be a named list
#' - nested MAGMA settings belong inside `gene_step1_args`, `gene_step2_args`,
#'   `reg_step1_args`, and `reg_step2_args`
#' - fishHook covariate-bearing inputs still follow the same rules as
#'   `prepare_somatic_scores()`: one row per regulatory element in
#'   `fishhook_covariate_data`, and a preformatted fishHook-ready object in
#'   `fishhook_covariates` if supplied
#'
#' The gene and regulatory location resources are backend-managed by the
#' package and are not user-facing arguments in this high-level wrapper.
#'
#' `conseguiR` now uses an object-first design. In plain terms, that means the
#' compute stages return R objects/bundles by default, and file writing is
#' optional rather than being the main API. If you supply `output_dir`,
#' `run_conseguiR()` will save stage artifacts there; if you leave
#' `output_dir = NULL`, it will still run the full workflow and return the
#' resulting objects in memory.
#'
#' Exact list formatting:
#'
#' `germline_args = list(
#'   gene_sample_size = 456348,
#'   reg_sample_size = 456348,
#'   gene_step1_args = list(annotation_window = c(35, 10)),
#'   gene_step2_args = list(gene_model = \"snp-wise=mean\", pval_use = c(\"SNP\", \"P\")),
#'   reg_step1_args = list(annotation_window = c(0, 0)),
#'   reg_step2_args = list(gene_model = \"snp-wise=mean\", pval_use = c(\"SNP\", \"P\"))
#' )`
#'
#' `somatic_args = list(
#'   gene_cv = NULL,
#'   dndscv_args = list(sm = \"192r_3w\"),
#'   eligible_gr = NULL,
#'   fishhook_covariate_data = covariate_dt,
#'   fishhook_args = list()
#' )`
#'
#' `epigenomic_args = list(
#'   bw_files = c(\"track1.bw\", \"track2.bw\", \"track3.bw\"),
#'   min_tracks = 3L,
#'   transform = \"log1p\"
#' )`
#'
#' @examples
#' names(formals(run_conseguiR))
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
  graph_rds_path = NULL,
  gg_nodes_path = NULL,
  gg_edges_path = NULL,
  output_dir = NULL,
  target_genes = 50L,
  germline_args = list(),
  somatic_args = list(),
  epigenomic_args = list(),
  scored_graph_args = list(),
  diffusion_args = list(),
  subgraph_args = list(),
  plot_args = list(),
  verbose = FALSE
) {
  paths <- .conseguiR_pipeline_paths(
    graph_rds_path = graph_rds_path,
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path
  )

  .conseguiR_call_external("run_conseguiR", .conseguiR_pipeline_args(
    gwas_sumstats = gwas_sumstats,
    somatic_maf = somatic_maf,
    reg_ref_path = reg_ref_path,
    reference_bfile = reference_bfile,
    dndscv_refdb = dndscv_refdb,
    epigenomic_track_dir = epigenomic_track_dir,
    epigenomic_tracks = epigenomic_tracks,
    paths = paths,
    output_dir = output_dir,
    target_genes = target_genes,
    germline_args = germline_args,
    somatic_args = somatic_args,
    epigenomic_args = epigenomic_args,
    scored_graph_args = scored_graph_args,
    diffusion_args = diffusion_args,
    subgraph_args = subgraph_args,
    plot_args = plot_args,
    verbose = verbose
  ))
}
