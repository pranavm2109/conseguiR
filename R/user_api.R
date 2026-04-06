#' @keywords internal
.conseguiR_runtime_env <- new.env(parent = baseenv())

#' @keywords internal
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

#' @keywords internal
.conseguiR_runtime_file <- function(relpath) {
  pkg_root <- if (exists(".conseguiR_pkg_root", inherits = FALSE)) .conseguiR_pkg_root else NULL
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
#' @param verbose Logical scalar. If `TRUE`, show progress bars and stage
#'   messages when available.
#'
#' @details
#' Input formatting rules:
#'
#' - `gwas_sumstats` can be either a file path or a data frame/data.table. When
#'   a table is supplied, it should contain columns that can be mapped to GWAS
#'   identifiers, chromosome, position, and p-value, for example
#'   `hm_variant_id`, `hm_chrom`, `hm_pos`, `p_value`, or the canonical columns
#'   `variant_id`, `chromosome`, `base_pair_location`, `p_value`.
#' - `somatic_maf` can be either a file path or a data frame/data.table. The
#'   table should contain sample identifier, chromosome, start, end, reference,
#'   and alternate allele columns in standard MAF-style names or names that the
#'   internal validator can harmonize.
#' - `reg_ref_path` must be a path to a tab-delimited regulatory-element file.
#'   The current pipeline expects at least four columns corresponding to
#'   chromosome, start, end, and regulatory-element identifier.
#' - `epigenomic_tracks` must be a character vector of full file paths to
#'   `.bw`/`.bigWig` files, for example
#'   `c(\"sample1.bw\", \"sample2.bw\", \"sample3.bw\")`.
#' - `epigenomic_track_dir` must be a single directory path containing bigWig
#'   files. Use either `epigenomic_tracks` or `epigenomic_track_dir`.
#'
#' @examples
#' validate_inputs()
#'
#' \dontrun{
#' validate_inputs(
#'   gwas_sumstats = "study_gwas.tsv",
#'   somatic_maf = "study_somatic.maf",
#'   reg_ref_path = "regulatory_elements.loc",
#'   epigenomic_tracks = c("track1.bw", "track2.bw", "track3.bw")
#' )
#' }
#'
#' @return A validation bundle containing validated objects and config.
#' @export
validate_inputs <- function(
  gwas_sumstats = NULL,
  somatic_maf = NULL,
  reg_ref_path = NULL,
  epigenomic_tracks = NULL,
  epigenomic_track_dir = NULL,
  epigenomic_exclude_patterns = c("_BL_", "_FL_"),
  verbose = FALSE
) {
  .conseguiR_external_fun("validate_inputs")(
    gwas_sumstats = gwas_sumstats,
    somatic_maf = somatic_maf,
    reg_ref_path = reg_ref_path,
    epigenomic_tracks = epigenomic_tracks,
    epigenomic_track_dir = epigenomic_track_dir,
    epigenomic_exclude_patterns = epigenomic_exclude_patterns,
    verbose = verbose
  )
}

#' Initialize backend graph resources
#'
#' Builds the package's unscored backend graph resources when they are missing
#' and the required raw graph resources are available.
#'
#' @inheritParams initialize_backend_graphs
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
#' MAGMA manual:
#' \url{https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf}
#'
#' @examples
#' \dontrun{
#' run_germline_gene_scoring(
#'   gwas_sumstats = "study_gwas.tsv",
#'   reference_bfile = "/path/to/g1000_eur/g1000_eur",
#'   sample_size = 456348,
#'   step1_args = list(annotation_window = c(35, 10)),
#'   step2_args = list(
#'     gene_model = "snp-wise=mean",
#'     pval_use = c("SNP", "P")
#'   )
#' )
#' }
#'
#' @return A germline gene score bundle.
#' @export
run_germline_gene_scoring <- function(
  gwas_sumstats,
  gene_loc_path = NULL,
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
#' MAGMA manual:
#' \url{https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf}
#'
#' @examples
#' \dontrun{
#' run_germline_regulatory_scoring(
#'   gwas_sumstats = "study_gwas.tsv",
#'   reference_bfile = "/path/to/g1000_eur/g1000_eur",
#'   sample_size = 456348,
#'   step1_args = list(annotation_window = c(0, 0)),
#'   step2_args = list(
#'     gene_model = "snp-wise=mean",
#'     pval_use = c("SNP", "P")
#'   )
#' )
#' }
#'
#' @return A germline regulatory score bundle.
#' @export
run_germline_regulatory_scoring <- function(
  gwas_sumstats,
  reg_loc_path = NULL,
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
#' MAGMA manual:
#' \url{https://ibg.colorado.edu/cdrom2021/Day10-posthuma/magma_session/manual_v1.09a.pdf}
#'
#' @examples
#' \dontrun{
#' prepare_germline_scores(
#'   gwas_sumstats = "study_gwas.tsv",
#'   reference_bfile = "/path/to/g1000_eur/g1000_eur",
#'   gene_step1_args = list(annotation_window = c(35, 10)),
#'   gene_step2_args = list(gene_model = "snp-wise=mean", pval_use = c("SNP", "P")),
#'   reg_step1_args = list(annotation_window = c(0, 0)),
#'   reg_step2_args = list(gene_model = "snp-wise=mean", pval_use = c("SNP", "P"))
#' )
#' }
#'
#' @return A germline score bundle with gene and regulatory score tables.
#' @export
prepare_germline_scores <- function(
  gwas_sumstats,
  reference_bfile,
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
#' - `dndscv_args`: a named list of additional dndscv arguments, for example
#'   `list(sm = \"192r_3w\", kc = \"cgc81\")`.
#'
#' dndscv documentation:
#' \url{https://rdrr.io/github/im3sanger/dndscv/man/dndscv.html}
#'
#' @examples
#' \dontrun{
#' run_somatic_gene_scoring(
#'   maf = "study_somatic.maf",
#'   refdb = "RefCDS_human_GRCh38.rda",
#'   dndscv_args = list(sm = "192r_3w")
#' )
#' }
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
#'   regulatory element and the covariate columns needed by fishHook.
#' - `fishhook_args`: a named list of additional fishHook arguments.
#'
#' fishHook documentation:
#' \url{https://rdrr.io/github/mskilab/fish.hook/man/FishHook.html}
#'
#' fishHook tutorial:
#' \url{https://mskilab.com/fishHook/tutorial.html}
#'
#' @examples
#' \dontrun{
#' run_somatic_regulatory_scoring(
#'   maf = "study_somatic.maf",
#'   reg_ref_path = "regulatory_elements.loc",
#'   fishhook_covariate_data = covariate_dt,
#'   fishhook_args = list()
#' )
#' }
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
#' \dontrun{
#' prepare_somatic_scores(
#'   maf = "study_somatic.maf",
#'   refdb = "RefCDS_human_GRCh38.rda",
#'   reg_ref_path = "regulatory_elements.loc",
#'   dndscv_args = list(sm = "192r_3w"),
#'   fishhook_covariate_data = covariate_dt
#' )
#' }
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
#' @param exclude_patterns Patterns used to exclude tracks discovered from a
#'   directory.
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
#' - `exclude_patterns`: a character vector of regex fragments used to exclude
#'   files discovered from `track_dir`.
#' - `summary_fun`: a function object such as `mean` or `max`.
#'
#' @examples
#' \dontrun{
#' prepare_epigenomic_scores(
#'   reg_ref_path = "regulatory_elements.loc",
#'   bw_files = c("track1.bw", "track2.bw", "track3.bw"),
#'   min_tracks = 3L,
#'   transform = "log1p"
#' )
#' }
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
  summary_fun = mean,
  verbose = FALSE
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
#'
#' You can either pass full score bundles from earlier stages or pass explicit
#' score tables directly through `gene_germline_scores`,
#' `reg_germline_scores`, `gene_somatic_scores`, `reg_somatic_scores`, and
#' `reg_epigenomic_scores`.
#'
#' @examples
#' \dontrun{
#' build_scored_gene_reg_graph(
#'   gene_germline_scores = "germline_gene_scores.tsv",
#'   reg_germline_scores = "germline_reg_scores.tsv",
#'   gene_somatic_scores = "somatic_gene_scores.tsv",
#'   reg_somatic_scores = "somatic_reg_scores.tsv",
#'   reg_epigenomic_scores = "epigenomic_reg_scores.tsv"
#' )
#' }
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
  save_outputs = TRUE,
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
#' @param python_path Optional explicit Python interpreter path.
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
#'
#' @examples
#' \dontrun{
#' scored_graph <- build_scored_gene_reg_graph(
#'   gene_germline_scores = "germline_gene_scores.tsv",
#'   reg_germline_scores = "germline_reg_scores.tsv",
#'   gene_somatic_scores = "somatic_gene_scores.tsv",
#'   reg_somatic_scores = "somatic_reg_scores.tsv",
#'   reg_epigenomic_scores = "epigenomic_reg_scores.tsv"
#' )
#'
#' run_gene_reg_diffusion(
#'   scored_graph = scored_graph,
#'   top_k = 3L,
#'   beta_epigenomic = 0.7
#' )
#' }
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
#' @param python_path Optional explicit Python interpreter path.
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
#' @examples
#' \dontrun{
#' diffusion <- run_gene_reg_diffusion(
#'   nodes_path = "gene_reg_graph_scored_nodes.tsv.gz",
#'   edges_path = "gene_reg_graph_scored_edges.tsv.gz"
#' )
#'
#' call_selected_subgraph(
#'   diffusion = diffusion,
#'   target_genes = 50L
#' )
#' }
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
  python_path = NULL,
  verbose = FALSE
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
    python_path = python_path,
    verbose = verbose
  )
}

#' Plot score tables from germline, somatic, or epigenomic bundles
#'
#' Creates a default score visualization from a score bundle or an explicit
#' score table.
#'
#' @param scores Optional score bundle returned by a scoring wrapper.
#' @param table Optional explicit score table.
#' @param which Optional bundle component name such as `gene_scores` or
#'   `reg_scores`.
#' @param plot_file_path Optional output path for the saved figure.
#' @param plot_type Plot type. Supported values are `ranked_points`,
#'   `top_bar`, and `histogram`.
#' @param top_n Number of top-ranked features used by ranked/bar plots.
#' @param feature_column Optional explicit feature-label column.
#' @param value_column Optional explicit numeric score column.
#' @param highlight_features Optional character vector of features to
#'   highlight.
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
#'   `prepare_somatic_scores()`, or `prepare_epigenomic_scores()`
#' - `table`: a data frame/data.table containing one feature label column and
#'   one numeric score column
#' - `which`: when `scores` contains multiple tables, use values like
#'   `"gene_scores"` or `"reg_scores"`
#' - `plot_file_path`: a single file path such as `"scores_plot.pdf"`
#'
#' @examples
#' \dontrun{
#' germline <- prepare_germline_scores(
#'   gwas_sumstats = "study_gwas.tsv",
#'   reference_bfile = "/path/to/g1000_eur/g1000_eur",
#'   gene_sample_size = 456348
#' )
#'
#' plot_scores(
#'   scores = germline,
#'   which = "gene_scores",
#'   plot_file_path = "germline_gene_scores.pdf",
#'   plot_type = "top_bar",
#'   top_n = 25L
#' )
#' }
#'
#' @return A plot bundle containing the ggplot object.
#' @export
plot_scores <- function(
  scores = NULL,
  table = NULL,
  which = NULL,
  plot_file_path = NULL,
  plot_type = "ranked_points",
  top_n = 25L,
  feature_column = NULL,
  value_column = NULL,
  highlight_features = NULL,
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
    plot_type = plot_type,
    top_n = top_n,
    feature_column = feature_column,
    value_column = value_column,
    highlight_features = highlight_features,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
  )
}

#' Plot diffusion results
#'
#' Creates a default diffusion visualization from a diffusion bundle or an
#' explicit diffusion table.
#'
#' @param diffusion Optional diffusion bundle returned by
#'   `run_gene_reg_diffusion()`.
#' @param table Optional explicit diffusion table.
#' @param which Which diffusion table to use: `all_genes` or `top_genes`.
#' @param plot_file_path Optional output path for the saved figure.
#' @param plot_type Plot type. Supported values are `ranked_points`,
#'   `top_bar`, and `histogram`.
#' @param top_n Number of top-ranked genes used by ranked/bar plots.
#' @param gene_column Optional explicit gene-label column.
#' @param score_column Optional explicit numeric diffusion-score column.
#' @param highlight_genes Optional character vector of genes to highlight.
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
#' - `diffusion`: the bundle returned by `run_gene_reg_diffusion()`
#' - `table`: a data frame/data.table with a gene-label column and a numeric
#'   diffusion score column
#' - `which`: either `"all_genes"` or `"top_genes"`
#' - `plot_file_path`: a single file path such as `"diffusion_plot.pdf"`
#'
#' @examples
#' \dontrun{
#' diffusion <- run_gene_reg_diffusion(
#'   nodes_path = "gene_reg_graph_scored_nodes.tsv.gz",
#'   edges_path = "gene_reg_graph_scored_edges.tsv.gz"
#' )
#'
#' plot_diffusion(
#'   diffusion = diffusion,
#'   which = "top_genes",
#'   plot_file_path = "diffusion_top_genes.pdf",
#'   highlight_genes = c("MYC", "PAX5", "BCL2")
#' )
#' }
#'
#' @return A plot bundle containing the ggplot object.
#' @export
plot_diffusion <- function(
  diffusion = NULL,
  table = NULL,
  which = "all_genes",
  plot_file_path = NULL,
  plot_type = "ranked_points",
  top_n = 50L,
  gene_column = NULL,
  score_column = NULL,
  highlight_genes = NULL,
  title = "conseguiR Diffusion Scores",
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = FALSE
) {
  .conseguiR_external_fun("plot_diffusion")(
    diffusion = diffusion,
    table = table,
    which = which,
    plot_file_path = plot_file_path,
    plot_type = plot_type,
    top_n = top_n,
    gene_column = gene_column,
    score_column = score_column,
    highlight_genes = highlight_genes,
    title = title,
    width = width,
    height = height,
    dpi = dpi,
    save_plot = save_plot,
    verbose = verbose
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
#' \dontrun{
#' selected_subgraph <- call_selected_subgraph(
#'   diffusion_path = "gene_reg_graph_diffusion_all_genes.tsv",
#'   target_genes = 50L
#' )
#'
#' plot_selected_subgraph(
#'   selected_subgraph = selected_subgraph,
#'   plot_file_path = "selected_subgraph.pdf",
#'   title = "Selected Gene Subgraph"
#' )
#' }
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
#' @param output_dir Output directory for pipeline artifacts.
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
#'   `bw_files`, `exclude_patterns`, `min_tracks`, `drop_mhc`, `transform`,
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
#' The gene and regulatory location resources are backend-managed by the
#' package and are not user-facing arguments in this high-level wrapper.
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
#'   exclude_patterns = c(\"_BL_\", \"_FL_\"),
#'   min_tracks = 3L,
#'   transform = \"log1p\"
#' )`
#'
#' @examples
#' \dontrun{
#' run_conseguiR(
#'   gwas_sumstats = "study_gwas.tsv",
#'   somatic_maf = "study_somatic.maf",
#'   reg_ref_path = "regulatory_elements.loc",
#'   reference_bfile = "/path/to/g1000_eur/g1000_eur",
#'   dndscv_refdb = "RefCDS_human_GRCh38.rda",
#'   epigenomic_tracks = c("track1.bw", "track2.bw", "track3.bw"),
#'   target_genes = 50L,
#'   germline_args = list(
#'     gene_sample_size = 456348,
#'     reg_sample_size = 456348,
#'     gene_step1_args = list(annotation_window = c(35, 10)),
#'     gene_step2_args = list(gene_model = "snp-wise=mean", pval_use = c("SNP", "P")),
#'     reg_step1_args = list(annotation_window = c(0, 0)),
#'     reg_step2_args = list(gene_model = "snp-wise=mean", pval_use = c("SNP", "P"))
#'   ),
#'   somatic_args = list(
#'     dndscv_args = list(sm = "192r_3w"),
#'     fishhook_covariate_data = covariate_dt
#'   ),
#'   epigenomic_args = list(
#'     bw_files = c("track1.bw", "track2.bw", "track3.bw"),
#'     min_tracks = 3L
#'   )
#' )
#' }
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
  output_dir = "data/processed",
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
    plot_args = plot_args,
    verbose = verbose
  )
}
