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
  old_wd <- getwd()
  on.exit(setwd(old_wd), add = TRUE)

  pkg_root <- .conseguiR_state$pkg_root
  if (is.character(pkg_root) && length(pkg_root) == 1L && nzchar(pkg_root) && dir.exists(pkg_root)) {
    setwd(pkg_root)
  }

  sys.source(api_path, envir = .conseguiR_runtime_env)
  assign(".loaded", TRUE, envir = .conseguiR_runtime_env)
  invisible()
}

#' @keywords internal
.conseguiR_reset_runtime_env <- function() {
  loaded_names <- ls(envir = .conseguiR_runtime_env, all.names = TRUE)
  if (length(loaded_names) > 0L) {
    rm(list = loaded_names, envir = .conseguiR_runtime_env)
  }

  .conseguiR_state$basilisk_env <- NULL
  .conseguiR_state$basilisk_env_pkgname <- NULL
  .conseguiR_state$basilisk_status <- NULL

  invisible()
}

#' Reload the sourced conseguiR runtime helpers
#'
#' Refreshes the mutable runtime environment used by the sourced external API
#' path. This is mainly useful during development, vignette rendering from a
#' live checkout, or after syncing updated source files into an existing R
#' session. Installed-package users typically do not need to call this.
#'
#' @param pkg_root Optional repository root to use when resolving sourced
#'   runtime files. When omitted, `conseguiR` reuses the currently recorded
#'   package root if available.
#' @param rebind Logical scalar. If `TRUE`, refresh convenience objects such as
#'   `copy`, `fread`, and `fwrite` inside the runtime environment after
#'   reloading.
#'
#' @return Invisibly returns the runtime environment.
#' @export
reload_conseguiR_runtime <- function(pkg_root = NULL, rebind = TRUE) {
  if (!is.null(pkg_root)) {
    pkg_root <- normalizePath(pkg_root, winslash = "/", mustWork = TRUE)
    .conseguiR_state$pkg_root <- pkg_root
  }

  .conseguiR_reset_runtime_env()
  .conseguiR_load_external_api()

  if (isTRUE(rebind)) {
    if (requireNamespace("data.table", quietly = TRUE)) {
      assign("copy", data.table::copy, envir = .conseguiR_runtime_env)
      assign("as.data.table", data.table::as.data.table, envir = .conseguiR_runtime_env)
      assign("fread", data.table::fread, envir = .conseguiR_runtime_env)
      assign("fwrite", data.table::fwrite, envir = .conseguiR_runtime_env)
    }

    assign(
      "deep_copy_object",
      function(x) unserialize(serialize(x, NULL)),
      envir = .conseguiR_runtime_env
    )
  }

  invisible(.conseguiR_runtime_env)
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
  candidate_pool_size,
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
    candidate_pool_size = candidate_pool_size,
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
#' @param reg_ref_path Optional regulatory-element reference path. When `NULL`,
#'   `conseguiR` uses its backend-owned ENCODE-derived regulatory universe.
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
  verbose = TRUE
) {
  if (is.null(reg_ref_path) && (!is.null(epigenomic_tracks) || !is.null(epigenomic_track_dir))) {
    reg_ref_path <- .conseguiR_default_reg_loc_path()
  }

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
#' Materializes packaged backend graph seeds into a writable backend directory
#' for downstream stages.
#'
#' @param backend_dir Optional backend directory. When `NULL`, `conseguiR`
#'   uses its default backend resource location.
#' @param build_gene_reg Logical scalar. If `TRUE`, ensure the gene-regulatory
#'   backend graph resources are available.
#' @param build_gene_gene Logical scalar. If `TRUE`, ensure the gene-gene
#'   backend graph resources are available.
#' @param force Logical scalar. If `TRUE`, re-materialize backend resources even
#'   when previous outputs already exist.
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
#'   output paths, and graph materialization status.
#' @export
initialize_backend_graphs <- function(
  backend_dir = NULL,
  build_gene_reg = TRUE,
  build_gene_gene = TRUE,
  force = FALSE,
  strict = TRUE,
  quiet = FALSE,
  verbose = TRUE
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

#' Resolve the backend gene-location resource path
#'
#' Returns the MAGMA-compatible gene location file that `conseguiR` resolves
#' from its packaged/backend resources.
#'
#' In normal installed-package use this resource should always be available.
#' `strict = FALSE` is mainly useful for advanced fallback workflows,
#' development checkouts, or partial environments where you want to handle a
#' missing resource yourself.
#'
#' @param strict Logical scalar. If `TRUE` (recommended default), error when
#'   the gene location resource cannot be resolved. If `FALSE`, return `NULL`
#'   instead.
#'
#' @examples
#' get_backend_gene_loc_path()
#'
#' @return A normalized absolute path to the backend gene location file, or
#'   `NULL` when `strict = FALSE` and the resource cannot be resolved.
#' @export
get_backend_gene_loc_path <- function(strict = TRUE) {
  path <- .conseguiR_default_gene_loc_path()
  if (is.null(path) && isTRUE(strict)) {
    stop(
      "Could not resolve the backend gene location resource. ",
      "Try checking that the packaged backend resources are installed."
    )
  }

  path
}

#' Resolve the backend regulatory-location resource path
#'
#' Returns the MAGMA-compatible regulatory-element location file that
#' `conseguiR` resolves from its packaged/backend resources.
#'
#' In normal installed-package use this resource should always be available.
#' `strict = FALSE` is mainly useful for advanced fallback workflows,
#' development checkouts, or partial environments where you want to handle a
#' missing resource yourself.
#'
#' @param strict Logical scalar. If `TRUE` (recommended default), error when
#'   the regulatory-element location resource cannot be resolved. If `FALSE`,
#'   return `NULL` instead.
#'
#' @examples
#' get_backend_reg_loc_path()
#'
#' @return A normalized absolute path to the backend regulatory-element
#'   location file, or `NULL` when `strict = FALSE` and the resource cannot be
#'   resolved.
#' @export
get_backend_reg_loc_path <- function(strict = TRUE) {
  path <- .conseguiR_default_reg_loc_path()
  if (is.null(path) && isTRUE(strict)) {
    stop(
      "Could not resolve the backend regulatory-element location resource. ",
      "Try checking that the packaged backend resources are installed."
    )
  }

  path
}

#' Manually set the MAGMA executable path
#'
#' Stores a validated absolute path to the MAGMA executable for the current R
#' session. This is a fallback for cases where `conseguiR` cannot auto-discover
#' MAGMA on package load.
#'
#' @param path Absolute path to the MAGMA executable.
#'
#' @examples
#' \dontrun{
#' setMAGMAPath("/absolute/path/to/magma")
#' }
#'
#' @return The normalized MAGMA path, invisibly.
#' @export
setMAGMAPath <- function(path) {
  resolved <- .conseguiR_resolve_magma_path(magma_path = path, must_work = TRUE)
  options(conseguiR.magma_path = resolved)
  .conseguiR_state$magma_path <- resolved
  invisible(resolved)
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
#' @param keep_intermediates Whether to keep intermediate MAGMA files.
#' @param annotation_window Optional MAGMA step 1 window passed as
#'   `window=before,after`.
#' @param filter_path Optional MAGMA step 1 filter file.
#' @param ignore_strand Whether to add `ignore-strand` in MAGMA step 1.
#' @param nonhuman Whether to add `nonhuman` in MAGMA step 1.
#' @param annotate_modifiers Optional additional MAGMA step 1 annotation
#'   modifiers appended after `--annotate`.
#' @param step1_general_args Named list converted into additional MAGMA step 1
#'   flags.
#' @param step1_extra_args Character vector appended verbatim to the MAGMA step
#'   1 command.
#' @param gene_model Optional MAGMA step 2 gene model.
#' @param gene_model_modifiers Optional additional modifiers appended after the
#'   MAGMA step 2 gene model value.
#' @param genes_only Whether to add `--genes-only` in MAGMA step 2.
#' @param pval_use Optional two-element character vector naming the SNP-ID and
#'   p-value columns supplied to MAGMA via `use=...`.
#' @param pval_snp_id Optional explicit SNP-ID column name supplied to MAGMA
#'   via `snp-id=...`.
#' @param pval_pval Optional explicit p-value column name supplied to MAGMA via
#'   `pval=...`.
#' @param pval_duplicate Optional MAGMA duplicate-SNP handling mode.
#' @param bfile_synonyms Optional PLINK synonym file supplied to MAGMA via
#'   `synonyms=...`.
#' @param bfile_synonym_dup Optional MAGMA synonym duplicate handling mode.
#' @param gene_settings Named list of MAGMA `--gene-settings` modifiers, for
#'   example SNP filters, pruning controls, and permutation controls.
#' @param batch Optional MAGMA `--batch` value.
#' @param seed Optional MAGMA `--seed` value.
#' @param big_data Optional MAGMA `--big-data` switch. Use `TRUE` for `on` and
#'   `FALSE` for `off`.
#' @param step2_general_args Named list converted into additional MAGMA step 2
#'   flags.
#' @param step2_extra_args Character vector appended verbatim to the MAGMA step
#'   2 command.
#' @param step1_args Named list of MAGMA step 1 arguments. Supported entries
#'   include `annotation_window`, `filter_path`, `ignore_strand`, `nonhuman`,
#'   `annotate_modifiers`, `general_args`, and `extra_args`.
#' @param step2_args Named list of MAGMA step 2 arguments. Supported entries
#'   include `gene_model`, `gene_model_modifiers`, `genes_only`, nested
#'   `pval = list(...)`, nested `bfile = list(...)`, `gene_settings`, `batch`,
#'   `seed`, `big_data`, `general_args`, and `extra_args`.
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
#'   annotation_window = c(15, 15),
#'   filter_path = NULL,
#'   ignore_strand = FALSE,
#'   nonhuman = FALSE,
#'   annotate_modifiers = NULL,
#'   general_args = list(),
#'   extra_args = character()
#' )`
#'
#' `step2_args = list(
#'   gene_model = \"snp-wise=mean\",
#'   gene_model_modifiers = NULL,
#'   genes_only = TRUE,
#'   pval = list(use = c(\"SNP\", \"P\"), duplicate = \"drop\"),
#'   bfile = list(synonyms = NULL, synonym_dup = NULL),
#'   gene_settings = list(),
#'   batch = NULL,
#'   seed = NULL,
#'   big_data = NULL,
#'   general_args = list(),
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
#'   The usual controls are `annotation_window`, `filter_path`,
#'   `ignore_strand`, and `nonhuman`.
#' - MAGMA step 2 (`step2_args`) controls how SNP-level p-values are read and
#'   how feature-level statistics are fit. The usual controls are
#'   `gene_model`, `genes_only`, the `pval` list, the `bfile` list,
#'   `gene_settings`, `batch`, `seed`, and `big_data`.
#' - `pval_use = c("SNP", "P")` is the compact MAGMA syntax for the SNP-ID and
#'   p-value columns; `pval_snp_id` and `pval_pval` expose the explicit
#'   `snp-id=...` and `pval=...` forms.
#' - `sample_size` and `sample_size_col` correspond to MAGMA's fixed `N=...`
#'   and column-driven `ncol=...` p-value modes.
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
  keep_intermediates = FALSE,
  annotation_window = c(15, 15),
  filter_path = NULL,
  ignore_strand = FALSE,
  nonhuman = FALSE,
  annotate_modifiers = NULL,
  step1_general_args = list(),
  step1_extra_args = character(),
  gene_model = NULL,
  gene_model_modifiers = NULL,
  genes_only = TRUE,
  pval_use = NULL,
  pval_snp_id = NULL,
  pval_pval = NULL,
  pval_duplicate = NULL,
  bfile_synonyms = NULL,
  bfile_synonym_dup = NULL,
  gene_settings = list(),
  batch = NULL,
  seed = NULL,
  big_data = NULL,
  step2_general_args = list(),
  step2_extra_args = character(),
  step1_args = list(),
  step2_args = list(),
  verbose = TRUE
) {
  gene_loc_path <- gene_loc_path %||% .conseguiR_default_gene_loc_path()
  if (is.null(gene_loc_path)) {
    stop("No gene location resource was provided and no backend gene location resource could be found.")
  }

  args <- list(
    gwas_sumstats = gwas_sumstats,
    gene_loc_path = gene_loc_path,
    reference_bfile = reference_bfile,
    verbose = verbose
  )

  if (!missing(output_prefix)) args$output_prefix <- output_prefix
  if (!missing(sample_size)) args$sample_size <- sample_size
  if (!missing(sample_size_col)) args$sample_size_col <- sample_size_col
  if (!missing(magma_path)) args$magma_path <- magma_path
  if (!missing(keep_intermediates)) args$keep_intermediates <- keep_intermediates
  if (!missing(annotation_window)) args$annotation_window <- annotation_window
  if (!missing(filter_path)) args$filter_path <- filter_path
  if (!missing(ignore_strand)) args$ignore_strand <- ignore_strand
  if (!missing(nonhuman)) args$nonhuman <- nonhuman
  if (!missing(annotate_modifiers)) args$annotate_modifiers <- annotate_modifiers
  if (!missing(step1_general_args)) args$step1_general_args <- step1_general_args
  if (!missing(step1_extra_args)) args$step1_extra_args <- step1_extra_args
  if (!missing(gene_model)) args$gene_model <- gene_model
  if (!missing(gene_model_modifiers)) args$gene_model_modifiers <- gene_model_modifiers
  if (!missing(genes_only)) args$genes_only <- genes_only
  if (!missing(pval_use)) args$pval_use <- pval_use
  if (!missing(pval_snp_id)) args$pval_snp_id <- pval_snp_id
  if (!missing(pval_pval)) args$pval_pval <- pval_pval
  if (!missing(pval_duplicate)) args$pval_duplicate <- pval_duplicate
  if (!missing(bfile_synonyms)) args$bfile_synonyms <- bfile_synonyms
  if (!missing(bfile_synonym_dup)) args$bfile_synonym_dup <- bfile_synonym_dup
  if (!missing(gene_settings)) args$gene_settings <- gene_settings
  if (!missing(batch)) args$batch <- batch
  if (!missing(seed)) args$seed <- seed
  if (!missing(big_data)) args$big_data <- big_data
  if (!missing(step2_general_args)) args$step2_general_args <- step2_general_args
  if (!missing(step2_extra_args)) args$step2_extra_args <- step2_extra_args
  if (!missing(step1_args)) args$step1_args <- step1_args
  if (!missing(step2_args)) args$step2_args <- step2_args

  do.call(.conseguiR_external_fun("run_germline_gene_scoring"), args)
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
#' rather than genes.
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
  keep_intermediates = FALSE,
  annotation_window = c(15, 15),
  filter_path = NULL,
  ignore_strand = FALSE,
  nonhuman = FALSE,
  annotate_modifiers = NULL,
  step1_general_args = list(),
  step1_extra_args = character(),
  gene_model = NULL,
  gene_model_modifiers = NULL,
  genes_only = TRUE,
  pval_use = NULL,
  pval_snp_id = NULL,
  pval_pval = NULL,
  pval_duplicate = NULL,
  bfile_synonyms = NULL,
  bfile_synonym_dup = NULL,
  gene_settings = list(),
  batch = NULL,
  seed = NULL,
  big_data = NULL,
  step2_general_args = list(),
  step2_extra_args = character(),
  step1_args = list(),
  step2_args = list(),
  verbose = TRUE
) {
  reg_loc_path <- reg_loc_path %||% .conseguiR_default_reg_loc_path()
  if (is.null(reg_loc_path)) {
    stop("No regulatory location resource was provided and no backend regulatory location resource could be found.")
  }

  args <- list(
    gwas_sumstats = gwas_sumstats,
    reg_loc_path = reg_loc_path,
    reference_bfile = reference_bfile,
    verbose = verbose
  )

  if (!missing(output_prefix)) args$output_prefix <- output_prefix
  if (!missing(sample_size)) args$sample_size <- sample_size
  if (!missing(sample_size_col)) args$sample_size_col <- sample_size_col
  if (!missing(magma_path)) args$magma_path <- magma_path
  if (!missing(keep_intermediates)) args$keep_intermediates <- keep_intermediates
  if (!missing(annotation_window)) args$annotation_window <- annotation_window
  if (!missing(filter_path)) args$filter_path <- filter_path
  if (!missing(ignore_strand)) args$ignore_strand <- ignore_strand
  if (!missing(nonhuman)) args$nonhuman <- nonhuman
  if (!missing(annotate_modifiers)) args$annotate_modifiers <- annotate_modifiers
  if (!missing(step1_general_args)) args$step1_general_args <- step1_general_args
  if (!missing(step1_extra_args)) args$step1_extra_args <- step1_extra_args
  if (!missing(gene_model)) args$gene_model <- gene_model
  if (!missing(gene_model_modifiers)) args$gene_model_modifiers <- gene_model_modifiers
  if (!missing(genes_only)) args$genes_only <- genes_only
  if (!missing(pval_use)) args$pval_use <- pval_use
  if (!missing(pval_snp_id)) args$pval_snp_id <- pval_snp_id
  if (!missing(pval_pval)) args$pval_pval <- pval_pval
  if (!missing(pval_duplicate)) args$pval_duplicate <- pval_duplicate
  if (!missing(bfile_synonyms)) args$bfile_synonyms <- bfile_synonyms
  if (!missing(bfile_synonym_dup)) args$bfile_synonym_dup <- bfile_synonym_dup
  if (!missing(gene_settings)) args$gene_settings <- gene_settings
  if (!missing(batch)) args$batch <- batch
  if (!missing(seed)) args$seed <- seed
  if (!missing(big_data)) args$big_data <- big_data
  if (!missing(step2_general_args)) args$step2_general_args <- step2_general_args
  if (!missing(step2_extra_args)) args$step2_extra_args <- step2_extra_args
  if (!missing(step1_args)) args$step1_args <- step1_args
  if (!missing(step2_args)) args$step2_args <- step2_args

  do.call(.conseguiR_external_fun("run_germline_regulatory_scoring"), args)
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
#'   `nonhuman`, `annotate_modifiers`, `general_args`, `extra_args`
#' - step 2 lists: `gene_model`, `gene_model_modifiers`, `genes_only`,
#'   `pval = list(...)`, `bfile = list(...)`, `gene_settings`, `batch`,
#'   `seed`, `big_data`, `general_args`, `extra_args`
#'
#' Example:
#'
#' `prepare_germline_scores(
#'   gwas_sumstats = gwas_path,
#'   reference_bfile = \"/path/to/g1000_eur/g1000_eur\",
#'   gene_step1_args = list(annotation_window = c(15, 15)),
#'   gene_step2_args = list(gene_model = \"snp-wise=mean\", pval_use = c(\"SNP\", \"P\")),
#'   reg_step1_args = list(annotation_window = c(15, 15)),
#'   reg_step2_args = list(gene_model = \"snp-wise=mean\", pval_use = c(\"SNP\", \"P\"))
#' )`
#'
#' Practical decision rules:
#'
#' - use the same `reference_bfile` for both the gene and regulatory runs
#' - use `gene_sample_size` / `gene_sample_size_col` for the gene run and
#'   `reg_sample_size` / `reg_sample_size_col` for the regulatory run when the
#'   two branches need different MAGMA p-value inputs
#' - keep `gene_step1_args` and `reg_step1_args` separate when you want
#'   different annotation behavior between genes and regulatory elements
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
  gene_sample_size = NULL,
  gene_sample_size_col = NULL,
  reg_sample_size = NULL,
  reg_sample_size_col = NULL,
  gene_step1_args = list(annotation_window = c(15, 15)),
  gene_step2_args = list(),
  reg_step1_args = list(annotation_window = c(15, 15)),
  reg_step2_args = list(),
  shared_args = list(),
  verbose = TRUE
) {
  .conseguiR_external_fun("prepare_germline_scores")(
    gwas_sumstats = gwas_sumstats,
    reference_bfile = reference_bfile,
    gene_output_prefix = gene_output_prefix,
    reg_output_prefix = reg_output_prefix,
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
#' This wrapper exposes the relevant `dndscv()` parameter surface used by the
#' package, while still allowing additional passthrough through `dndscv_args`.
#'
#' @param maf Somatic MAF path or table.
#' @param refdb dndscv reference database path.
#' @param output_path Optional output path for saved scores.
#' @param cv Optional dndscv covariates.
#' @param max_muts_per_gene_per_sample dndscv gene-level mutation cap.
#' @param max_coding_muts_per_sample dndscv sample-level coding mutation cap.
#' @param gene_list Optional vector of genes passed to `dndscv()`.
#' @param sm dndscv substitution model.
#' @param kc dndscv known-cancer-gene setting.
#' @param use_indel_sites Whether dndscv should use indel sites.
#' @param min_indels Minimum number of indels for dndscv indel modeling.
#' @param maxcovs Maximum number of covariates passed to dndscv.
#' @param constrain_wnon_wspl Whether to constrain `wnon == wspl` in dndscv.
#' @param outp dndscv output level.
#' @param numcode Genetic code used by dndscv.
#' @param outmats Whether dndscv should return count/exposure matrices.
#' @param mingenecovs Minimum number of genes for dndscv covariate modeling.
#' @param onesided dndscv one-sided testing switch. Defaults to `TRUE` so the
#'   returned somatic gene scores follow the directional enrichment-oriented
#'   interpretation used throughout the package.
#' @param dc Optional dndscv duplex-coverage vector.
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
#' - by default, `conseguiR` asks `dndscv()` for one-sided positive and
#'   negative selection tests when that API surface is available. In practice
#'   this means the package prefers directional columns such as `ppos_cv` and
#'   `pneg_cv` for score extraction.
#' - if those one-sided columns are not present, `conseguiR` falls back to the
#'   older signed two-sided extraction path so older `dndscv` outputs still
#'   work.
#' - if you supply `cv`, it should already be formatted the way `dndscv()`
#'   expects it. `conseguiR` does not reshape arbitrary covariate tables into a
#'   dndscv-ready object for you.
#' - `refdb` is not just any annotation file. It must be a dndscv-compatible
#'   reference `.rda` resource built for the same genome build as the MAF you
#'   are scoring.
#'
#' Minimal examples:
#'
#' `dndscv_args = list()`
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
  gene_list = NULL,
  sm = "192r_3w",
  kc = "cgc81",
  use_indel_sites = TRUE,
  min_indels = 5L,
  maxcovs = 20L,
  constrain_wnon_wspl = TRUE,
  outp = 3L,
  numcode = 1L,
  outmats = FALSE,
  mingenecovs = 500L,
  onesided = TRUE,
  dc = NULL,
  dndscv_args = list(),
  verbose = TRUE
) {
  args <- list(
    maf = maf,
    refdb = refdb,
    verbose = verbose
  )

  if (!missing(output_path)) args$output_path <- output_path
  if (!missing(cv)) args$cv <- cv
  if (!missing(max_muts_per_gene_per_sample)) args$max_muts_per_gene_per_sample <- max_muts_per_gene_per_sample
  if (!missing(max_coding_muts_per_sample)) args$max_coding_muts_per_sample <- max_coding_muts_per_sample
  if (!missing(gene_list)) args$gene_list <- gene_list
  if (!missing(sm)) args$sm <- sm
  if (!missing(kc)) args$kc <- kc
  if (!missing(use_indel_sites)) args$use_indel_sites <- use_indel_sites
  if (!missing(min_indels)) args$min_indels <- min_indels
  if (!missing(maxcovs)) args$maxcovs <- maxcovs
  if (!missing(constrain_wnon_wspl)) args$constrain_wnon_wspl <- constrain_wnon_wspl
  if (!missing(outp)) args$outp <- outp
  if (!missing(numcode)) args$numcode <- numcode
  if (!missing(outmats)) args$outmats <- outmats
  if (!missing(mingenecovs)) args$mingenecovs <- mingenecovs
  if (!missing(onesided)) args$onesided <- onesided
  if (!missing(dc)) args$dc <- dc
  if (!missing(dndscv_args)) args$dndscv_args <- dndscv_args

  do.call(.conseguiR_external_fun("run_somatic_gene_scoring"), args)
}

#' Run fishHook somatic regulatory scoring
#'
#' This wrapper exposes the relevant fishHook constructor and
#' `score.hypotheses()` parameter surfaces used by the package, while still
#' allowing additional passthrough through `fishhook_args`.
#'
#' @inheritParams run_somatic_gene_scoring
#' @param reg_ref_path Regulatory-element reference path.
#' @param eligible_gr Optional fishHook eligible territory.
#' @param fishhook_covariates Optional fishHook covariate objects/specifications.
#' @param fishhook_covariate_data Optional tabular covariate data.
#' @param idcol Sample identifier column for fishHook.
#' @param constructor_out_path Optional fishHook constructor output path.
#' @param constructor_use_local_mut_density Whether fishHook should add the
#'   local-mutation-density constructor branch.
#' @param constructor_local_mut_density_bin Bin size for fishHook local mutation
#'   density.
#' @param constructor_mc_cores fishHook constructor parallel worker count.
#' @param constructor_na_rm Whether fishHook constructor covariates drop `NA`.
#' @param constructor_pad fishHook constructor padding.
#' @param constructor_max_slice fishHook constructor `max.slice`.
#' @param constructor_ff_chunk fishHook constructor `ff.chunk`.
#' @param constructor_max_chunk fishHook constructor `max.chunk`.
#' @param constructor_idcap fishHook constructor `idcap`.
#' @param constructor_weight_events Whether fishHook constructor should weight
#'   events.
#' @param constructor_nb fishHook constructor negative-binomial switch.
#' @param score_sets Optional sets passed to `fishHook::score.hypotheses()`.
#' @param score_model Optional score model passed to fishHook.
#' @param score_return_model Whether to request the fishHook fitted model.
#' @param score_nb Optional score-side negative-binomial switch.
#' @param score_iter Optional score-side maximum iteration count.
#' @param score_subsample Optional fishHook score subsample size.
#' @param score_seed Optional fishHook score seed.
#' @param score_verbose Optional fishHook score verbosity override.
#' @param score_mc_cores Optional fishHook score parallel worker count.
#' @param score_p_randomized Optional fishHook randomized-p setting.
#' @param score_class_return Whether fishHook should return the classed result
#'   object.
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
#' - when fishHook returns directional significance columns, `conseguiR`
#'   extracts regulatory scores from enrichment-side `p` and depletion-side
#'   `p.neg` directly rather than reconstructing direction later from a generic
#'   p-value.
#' - if you do not have a custom territory or covariates yet, `eligible_gr =
#'   NULL` and `fishhook_covariate_data = NULL` are reasonable starting points.
#' - `idcol` should match the sample identifier column name used in the somatic
#'   mutation table after harmonization. In most MAF-based workflows the
#'   default `Tumor_Sample_Barcode` is the right choice.
#'
#' Minimal covariate-data example:
#'
#' `covariate_dt = data.frame(
#'   reg_elem_id = c(\"EH38E0080197\", \"EH38E2084302\"),
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
  reg_ref_path = NULL,
  output_path = NULL,
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  idcol = "Tumor_Sample_Barcode",
  constructor_out_path = NULL,
  constructor_use_local_mut_density = FALSE,
  constructor_local_mut_density_bin = 1e6,
  constructor_mc_cores = 1L,
  constructor_na_rm = TRUE,
  constructor_pad = 0,
  constructor_max_slice = 1e5,
  constructor_ff_chunk = 1e6,
  constructor_max_chunk = 1e12,
  constructor_idcap = 1,
  constructor_weight_events = FALSE,
  constructor_nb = TRUE,
  score_sets = NULL,
  score_model = NULL,
  score_return_model = TRUE,
  score_nb = NULL,
  score_iter = NULL,
  score_subsample = NULL,
  score_seed = NULL,
  score_verbose = NULL,
  score_mc_cores = NULL,
  score_p_randomized = NULL,
  score_class_return = TRUE,
  fishhook_args = list(),
  verbose = TRUE
) {
  reg_ref_path <- reg_ref_path %||% .conseguiR_default_reg_loc_path()
  if (is.null(reg_ref_path)) {
    stop("No regulatory reference was provided and no backend ENCODE-derived regulatory resource could be found.")
  }

  args <- list(
    maf = maf,
    reg_ref_path = reg_ref_path,
    verbose = verbose
  )

  if (!missing(output_path)) args$output_path <- output_path
  if (!missing(eligible_gr)) args$eligible_gr <- eligible_gr
  if (!missing(fishhook_covariates)) args$fishhook_covariates <- fishhook_covariates
  if (!missing(fishhook_covariate_data)) args$fishhook_covariate_data <- fishhook_covariate_data
  if (!missing(idcol)) args$idcol <- idcol
  if (!missing(constructor_out_path)) args$constructor_out_path <- constructor_out_path
  if (!missing(constructor_use_local_mut_density)) args$constructor_use_local_mut_density <- constructor_use_local_mut_density
  if (!missing(constructor_local_mut_density_bin)) args$constructor_local_mut_density_bin <- constructor_local_mut_density_bin
  if (!missing(constructor_mc_cores)) args$constructor_mc_cores <- constructor_mc_cores
  if (!missing(constructor_na_rm)) args$constructor_na_rm <- constructor_na_rm
  if (!missing(constructor_pad)) args$constructor_pad <- constructor_pad
  if (!missing(constructor_max_slice)) args$constructor_max_slice <- constructor_max_slice
  if (!missing(constructor_ff_chunk)) args$constructor_ff_chunk <- constructor_ff_chunk
  if (!missing(constructor_max_chunk)) args$constructor_max_chunk <- constructor_max_chunk
  if (!missing(constructor_idcap)) args$constructor_idcap <- constructor_idcap
  if (!missing(constructor_weight_events)) args$constructor_weight_events <- constructor_weight_events
  if (!missing(constructor_nb)) args$constructor_nb <- constructor_nb
  if (!missing(score_sets)) args$score_sets <- score_sets
  if (!missing(score_model)) args$score_model <- score_model
  if (!missing(score_return_model)) args$score_return_model <- score_return_model
  if (!missing(score_nb)) args$score_nb <- score_nb
  if (!missing(score_iter)) args$score_iter <- score_iter
  if (!missing(score_subsample)) args$score_subsample <- score_subsample
  if (!missing(score_seed)) args$score_seed <- score_seed
  if (!missing(score_verbose)) args$score_verbose <- score_verbose
  if (!missing(score_mc_cores)) args$score_mc_cores <- score_mc_cores
  if (!missing(score_p_randomized)) args$score_p_randomized <- score_p_randomized
  if (!missing(score_class_return)) args$score_class_return <- score_class_return
  if (!missing(fishhook_args)) args$fishhook_args <- fishhook_args

  do.call(.conseguiR_external_fun("run_somatic_regulatory_scoring"), args)
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
#' @param reg_ref_path Optional regulatory-element reference path. When `NULL`,
#'   `conseguiR` uses its backend-owned ENCODE-derived regulatory universe.
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
#'   `gene_max_coding_muts_per_sample`, `gene_list`, `sm`, `kc`,
#'   `use_indel_sites`, `min_indels`, `maxcovs`, `constrain_wnon_wspl`,
#'   `outp`, `numcode`, `outmats`, `mingenecovs`, `onesided`, `dc`,
#'   `dndscv_args`
#' - regulatory-level fishHook controls: `eligible_gr`,
#'   `fishhook_covariates`, `fishhook_covariate_data`, `fishhook_idcol`,
#'   constructor-side controls, score-side controls, and `fishhook_args`
#'
#' A good mental model is:
#'
#' - `prepare_somatic_scores()` does not fit one joint model
#' - it runs a gene-oriented dndscv branch and a regulatory-element-oriented
#'   fishHook branch separately
#' - then it returns both score tables together as one somatic bundle
#' - in the current package defaults, dndscv score extraction is one-sided when
#'   supported by the installed dndscv version, while fishHook extraction uses
#'   directional `p`/`p.neg` output when available
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
  reg_ref_path = NULL,
  gene_output_path = NULL,
  reg_output_path = NULL,
  gene_cv = NULL,
  gene_max_muts_per_gene_per_sample = 6L,
  gene_max_coding_muts_per_sample = 5000L,
  gene_list = NULL,
  sm = "192r_3w",
  kc = "cgc81",
  use_indel_sites = TRUE,
  min_indels = 5L,
  maxcovs = 20L,
  constrain_wnon_wspl = TRUE,
  outp = 3L,
  numcode = 1L,
  outmats = FALSE,
  mingenecovs = 500L,
  onesided = TRUE,
  dc = NULL,
  dndscv_args = list(),
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  fishhook_idcol = "Tumor_Sample_Barcode",
  constructor_out_path = NULL,
  constructor_use_local_mut_density = FALSE,
  constructor_local_mut_density_bin = 1e6,
  constructor_mc_cores = 1L,
  constructor_na_rm = TRUE,
  constructor_pad = 0,
  constructor_max_slice = 1e5,
  constructor_ff_chunk = 1e6,
  constructor_max_chunk = 1e12,
  constructor_idcap = 1,
  constructor_weight_events = FALSE,
  constructor_nb = TRUE,
  score_sets = NULL,
  score_model = NULL,
  score_return_model = TRUE,
  score_nb = NULL,
  score_iter = NULL,
  score_subsample = NULL,
  score_seed = NULL,
  score_verbose = NULL,
  score_mc_cores = NULL,
  score_p_randomized = NULL,
  score_class_return = TRUE,
  fishhook_args = list(),
  verbose = TRUE
) {
  reg_ref_path <- reg_ref_path %||% .conseguiR_default_reg_loc_path()
  if (is.null(reg_ref_path)) {
    stop("No regulatory reference was provided and no backend ENCODE-derived regulatory resource could be found.")
  }

  args <- list(
    maf = maf,
    refdb = refdb,
    reg_ref_path = reg_ref_path,
    verbose = verbose
  )

  if (!missing(gene_output_path)) args$gene_output_path <- gene_output_path
  if (!missing(reg_output_path)) args$reg_output_path <- reg_output_path
  if (!missing(gene_cv)) args$gene_cv <- gene_cv
  if (!missing(gene_max_muts_per_gene_per_sample)) args$gene_max_muts_per_gene_per_sample <- gene_max_muts_per_gene_per_sample
  if (!missing(gene_max_coding_muts_per_sample)) args$gene_max_coding_muts_per_sample <- gene_max_coding_muts_per_sample
  if (!missing(gene_list)) args$gene_list <- gene_list
  if (!missing(sm)) args$sm <- sm
  if (!missing(kc)) args$kc <- kc
  if (!missing(use_indel_sites)) args$use_indel_sites <- use_indel_sites
  if (!missing(min_indels)) args$min_indels <- min_indels
  if (!missing(maxcovs)) args$maxcovs <- maxcovs
  if (!missing(constrain_wnon_wspl)) args$constrain_wnon_wspl <- constrain_wnon_wspl
  if (!missing(outp)) args$outp <- outp
  if (!missing(numcode)) args$numcode <- numcode
  if (!missing(outmats)) args$outmats <- outmats
  if (!missing(mingenecovs)) args$mingenecovs <- mingenecovs
  if (!missing(onesided)) args$onesided <- onesided
  if (!missing(dc)) args$dc <- dc
  if (!missing(dndscv_args)) args$dndscv_args <- dndscv_args
  if (!missing(eligible_gr)) args$eligible_gr <- eligible_gr
  if (!missing(fishhook_covariates)) args$fishhook_covariates <- fishhook_covariates
  if (!missing(fishhook_covariate_data)) args$fishhook_covariate_data <- fishhook_covariate_data
  if (!missing(fishhook_idcol)) args$fishhook_idcol <- fishhook_idcol
  if (!missing(constructor_out_path)) args$constructor_out_path <- constructor_out_path
  if (!missing(constructor_use_local_mut_density)) args$constructor_use_local_mut_density <- constructor_use_local_mut_density
  if (!missing(constructor_local_mut_density_bin)) args$constructor_local_mut_density_bin <- constructor_local_mut_density_bin
  if (!missing(constructor_mc_cores)) args$constructor_mc_cores <- constructor_mc_cores
  if (!missing(constructor_na_rm)) args$constructor_na_rm <- constructor_na_rm
  if (!missing(constructor_pad)) args$constructor_pad <- constructor_pad
  if (!missing(constructor_max_slice)) args$constructor_max_slice <- constructor_max_slice
  if (!missing(constructor_ff_chunk)) args$constructor_ff_chunk <- constructor_ff_chunk
  if (!missing(constructor_max_chunk)) args$constructor_max_chunk <- constructor_max_chunk
  if (!missing(constructor_idcap)) args$constructor_idcap <- constructor_idcap
  if (!missing(constructor_weight_events)) args$constructor_weight_events <- constructor_weight_events
  if (!missing(constructor_nb)) args$constructor_nb <- constructor_nb
  if (!missing(score_sets)) args$score_sets <- score_sets
  if (!missing(score_model)) args$score_model <- score_model
  if (!missing(score_return_model)) args$score_return_model <- score_return_model
  if (!missing(score_nb)) args$score_nb <- score_nb
  if (!missing(score_iter)) args$score_iter <- score_iter
  if (!missing(score_subsample)) args$score_subsample <- score_subsample
  if (!missing(score_seed)) args$score_seed <- score_seed
  if (!missing(score_verbose)) args$score_verbose <- score_verbose
  if (!missing(score_mc_cores)) args$score_mc_cores <- score_mc_cores
  if (!missing(score_p_randomized)) args$score_p_randomized <- score_p_randomized
  if (!missing(score_class_return)) args$score_class_return <- score_class_return
  if (!missing(fishhook_args)) args$fishhook_args <- fishhook_args

  do.call(.conseguiR_external_fun("prepare_somatic_scores"), args)
}

#' Prepare regulatory epigenomic scores
#'
#' @param reg_ref_path Optional regulatory-element reference path. When `NULL`,
#'   `conseguiR` uses its backend-owned ENCODE-derived regulatory universe.
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
  reg_ref_path = NULL,
  track_dir = NULL,
  bw_files = NULL,
  output_path = NULL,
  min_tracks = 3L,
  drop_mhc = TRUE,
  transform = "log1p",
  return_diagnostics = TRUE,
  summary_fun = mean,
  verbose = TRUE
) {
  reg_ref_path <- reg_ref_path %||% .conseguiR_default_reg_loc_path()
  if (is.null(reg_ref_path)) {
    stop("No regulatory reference was provided and no backend ENCODE-derived regulatory resource could be found.")
  }

  args <- list(
    reg_ref_path = reg_ref_path,
    verbose = verbose
  )

  if (!missing(track_dir)) args$track_dir <- track_dir
  if (!missing(bw_files)) args$bw_files <- bw_files
  if (!missing(output_path)) args$output_path <- output_path
  if (!missing(min_tracks)) args$min_tracks <- min_tracks
  if (!missing(drop_mhc)) args$drop_mhc <- drop_mhc
  if (!missing(transform)) args$transform <- transform
  if (!missing(return_diagnostics)) args$return_diagnostics <- return_diagnostics
  if (!missing(summary_fun)) args$summary_fun <- summary_fun

  do.call(.conseguiR_external_fun("prepare_epigenomic_scores"), args)
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
#'     from = "EH38E0080197",
#'     to = "TP53",
#'     confidence = 0.9
#'   ),
#'   vertices = data.frame(
#'     name = c("TP53", "EH38E0080197"),
#'     node_id = c("TP53", "EH38E0080197"),
#'     node_type = c("gene", "reg")
#'   ),
#'   directed = TRUE
#' )
#'
#' scored_graph <- build_scored_gene_reg_graph(
#'   graph = toy_graph,
#'   graph_rds_path = tempfile(fileext = ".rds"),
#'   gene_germline_scores = data.frame(gene_id = "TP53", zstat = 2),
#'   reg_germline_scores = data.frame(reg_elem_id = "EH38E0080197", zstat = 1.2),
#'   gene_somatic_scores = data.frame(gene_id = "TP53", zstat = -1.5),
#'   reg_somatic_scores = data.frame(reg_elem_id = "EH38E0080197", zstat = 0.3),
#'   reg_epigenomic_scores = data.frame(reg_elem_id = "EH38E0080197", zstat = 2.4),
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
  verbose = TRUE
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
#' @param integration_weight_germline Germline weight used in the signed
#'   cross-modality integration score.
#' @param integration_weight_somatic Somatic weight used in the signed
#'   cross-modality integration score.
#' @param integration_weight_epigenomic Epigenomic weight used in the signed
#'   cross-modality integration score.
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
#'   `beta_epigenomic`, `integration_weight_germline`,
#'   `integration_weight_somatic`, `integration_weight_epigenomic`,
#'   `reg_signal_clip`: numeric scalars
#' - `positive_only`: a single logical value
#'
#' In plain language:
#'
#' - `top_k` controls how many regulatory neighbors contribute most strongly to
#'   each gene's update step
#' - `beta_germline`, `beta_somatic`, and `beta_epigenomic` control how much
#'   each modality contributes to the regulatory signal that is propagated
#' - `integration_weight_germline`, `integration_weight_somatic`, and
#'   `integration_weight_epigenomic` control the weighted signed Stouffer-style
#'   combination used to integrate post-diffusion modality scores into the main
#'   ranking columns
#' - `confidence_power` upweights or downweights high-confidence regulatory
#'   links relative to weaker ones
#' - `positive_only = TRUE` suppresses negative regulatory signal before
#'   diffusion
#' - `reg_signal_clip` caps extreme regulatory signal before propagation so one
#'   extreme feature does not dominate the update
#' - the output diffusion table keeps the legacy Euclidean norm
#'   (`prediff_norm`, `post_norm`) and a nonnegative magnitude-style
#'   vulnerability summary (`prediff_vulnerability`, `post_vulnerability`) for
#'   auditing, but the main ranking columns (`prediff_rank`, `post_rank`) are
#'   based on signed integrated scores (`prediff_integrated`,
#'   `post_integrated`) computed with a weighted signed Stouffer-style
#'   combination so negative modality contributions are penalized rather than
#'   clipped away
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
  integration_weight_germline = 1.0,
  integration_weight_somatic = 1.0,
  integration_weight_epigenomic = 1.0,
  positive_only = FALSE,
  reg_signal_clip = 5.0,
  top_n_to_save = 50L,
  python_path = NULL,
  verbose = TRUE
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
    integration_weight_germline = integration_weight_germline,
    integration_weight_somatic = integration_weight_somatic,
    integration_weight_epigenomic = integration_weight_epigenomic,
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
#' @param candidate_pool_size Candidate pool size for the solver. Must be at
#'   least as large as `target_genes`, and must not exceed the number of
#'   available diffusion-ranked genes. Larger values expand the solver search
#'   space and can increase runtime.
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
#' - `prize_column = "post_integrated"` means optimize using the signed
#'   post-diffusion integrated score, which rewards positive modality support
#'   and penalizes negative modality support
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
#' - `candidate_pool_size` must be at least as large as `target_genes`
#' - `candidate_pool_size` must not exceed the number of available
#'   diffusion-ranked genes
#' - larger `candidate_pool_size` values increase the solver search space and
#'   can increase runtime
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
  prize_column = "post_integrated",
  confidence_column = "confidence",
  edge_cost_column = "weight",
  python_path = NULL,
  verbose = TRUE
) {
  backend_paths <- .conseguiR_default_gene_gene_paths(
    gg_nodes_path = gg_nodes_path,
    gg_edges_path = gg_edges_path
  )
  args <- as.list(environment())
  args$gg_nodes_path <- backend_paths$gg_nodes_path
  args$gg_edges_path <- backend_paths$gg_edges_path
  args$backend_paths <- NULL

  .conseguiR_call_external("call_selected_subgraph", args)
}

#' Plot score tables from scoring or diffusion outputs
#'
#' Creates either a rank plot or a volcano plot for score tables.
#'
#' @param scores Optional bundle returned by a scoring wrapper or
#'   `run_gene_reg_diffusion()`.
#' @param table Optional explicit table.
#' @param which Optional bundle component name such as `gene_scores`,
#'   `reg_scores`, `all_genes`, or `top_genes`.
#' @param plot_file_path Optional output path for the saved figure.
#' @param test_tail One of `auto`, `one_tailed`, or `two_tailed`. This
#'   describes the score interpretation, not the plot geometry. When
#'   `scores` is a conseguiR bundle and `test_tail = "auto"`, the plotting
#'   helper uses the score semantics stored in that bundle.
#' @param plot_mode One of `auto`, `rank`, or `volcano`. Use this to choose the
#'   figure geometry explicitly when needed. In practice, use `plot_mode = "rank"`
#'   for MAGMA-style germline rankings and `plot_mode = "volcano"` for somatic
#'   plots where you want z-scores against `-log10(p)`.
#' @param feature_column Optional explicit feature-label column.
#' @param z_column Z-score column name.
#' @param p_value_column Optional explicit p-value column for volcano plots.
#' @param drop_tukey_outliers Logical scalar. If `TRUE`, remove extreme
#'   volcano-plot outliers using Tukey's upper-fence rule on `-log10(p)`.
#' @param clip_extreme_display Logical scalar. If `TRUE`, cap the displayed
#'   volcano axes for readability. If `FALSE`, plot the actual z-scores and
#'   `-log10(p)` values without display clipping.
#' @param label_features Optional character vector of features to label.
#' @param label_max_per_feature Positive integer. For any requested label, keep
#'   at most this many plotted labels, prioritizing the largest z-scores.
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
#' - when `scores` is a conseguiR bundle, `test_tail = "auto"` uses the bundle
#'   metadata instead of making you restate the score semantics manually
#' - when `table` is supplied directly, `test_tail` remains an interpretation
#'   hint because the raw table does not carry bundle metadata
#' - `plot_mode`, not `test_tail`, is what decides whether the figure is drawn
#'   as a rank plot or a volcano plot
#' - `drop_tukey_outliers = TRUE` trims extreme volcano-plot y-axis outliers
#'   before plotting, which can help when one feature is far more significant
#'   than the rest
#' - `plot_file_path`: a single file path such as `"scores_plot.pdf"`
#'
#' @examples
#' plot_scores(
#'   table = data.frame(feature_id = c("A", "B"), zstat = c(1, -1)),
#'   feature_column = "feature_id",
#'   z_column = "zstat",
#'   plot_mode = "rank",
#'   save_plot = FALSE
#' )
#'
#' plot_scores(
#'   table = data.frame(
#'     feature_id = c("A", "B"),
#'     zstat = c(3.2, 1.4),
#'     p_value = c(1e-5, 0.03)
#'   ),
#'   feature_column = "feature_id",
#'   z_column = "zstat",
#'   p_value_column = "p_value",
#'   plot_mode = "volcano",
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
  plot_mode = "auto",
  feature_column = NULL,
  z_column = "zstat",
  p_value_column = NULL,
  drop_tukey_outliers = FALSE,
  clip_extreme_display = FALSE,
  label_features = NULL,
  label_max_per_feature = 1L,
  title = "conseguiR Scores",
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
) {
  .conseguiR_external_fun("plot_scores")(
    scores = scores,
    table = table,
    which = which,
    plot_file_path = plot_file_path,
    test_tail = test_tail,
    plot_mode = plot_mode,
    feature_column = feature_column,
    z_column = z_column,
    p_value_column = p_value_column,
    drop_tukey_outliers = drop_tukey_outliers,
    clip_extreme_display = clip_extreme_display,
    label_features = label_features,
    label_max_per_feature = label_max_per_feature,
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
#' @param label_max_per_feature Positive integer limiting the number of labels
#'   retained per requested feature label.
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
  label_max_per_feature = 1L,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
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
    label_max_per_feature = label_max_per_feature,
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
#' @param label_max_per_feature Positive integer limiting the number of labels
#'   retained per requested feature label.
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
  label_max_per_feature = 1L,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
) {
  .conseguiR_external_fun("plot_germline_reg_scores")(
    germline_scores = germline_scores,
    scored_graph = scored_graph,
    nodes_path = nodes_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    label_max_per_feature = label_max_per_feature,
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
#' @param label_max_per_feature Positive integer limiting the number of labels
#'   retained per requested feature label.
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
  label_max_per_feature = 1L,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
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
    label_max_per_feature = label_max_per_feature,
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
#' @param label_max_per_feature Positive integer limiting the number of labels
#'   retained per requested feature label.
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
  label_max_per_feature = 1L,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
) {
  .conseguiR_external_fun("plot_somatic_reg_scores")(
    somatic_scores = somatic_scores,
    scored_graph = scored_graph,
    nodes_path = nodes_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    label_max_per_feature = label_max_per_feature,
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
#' @param label_max_per_feature Positive integer limiting the number of labels
#'   retained per requested feature label.
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
  label_max_per_feature = 1L,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
) {
  .conseguiR_external_fun("plot_epigenomic_gene_scores")(
    diffusion = diffusion,
    diffusion_path = diffusion_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    label_max_per_feature = label_max_per_feature,
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
#' @param label_max_per_feature Positive integer limiting the number of labels
#'   retained per requested feature label.
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
  label_max_per_feature = 1L,
  title = NULL,
  width = 10,
  height = 7,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
) {
  .conseguiR_external_fun("plot_epigenomic_reg_scores")(
    epigenomic_scores = epigenomic_scores,
    scored_graph = scored_graph,
    nodes_path = nodes_path,
    plot_file_path = plot_file_path,
    label_features = label_features,
    label_max_per_feature = label_max_per_feature,
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
#' @param label_top_gwas_snp Logical scalar or non-negative integer. Use
#'   `FALSE` or `0` to disable GWAS SNP labels, `TRUE` for one top GWAS SNP,
#'   or a positive integer to label that many top GWAS SNPs in the plotted
#'   locus.
#' @param rsid_pmid Optional rsID-to-PMID evidence table with at least
#'   `rsid` and `pmid` columns.
#' @param label_top_lit_snps Integer count of literature-backed SNPs to label.
#'   SNPs are restricted to regulatory elements in the locus. If no
#'   literature-backed SNPs survive the lookup/filtering steps, `conseguiR`
#'   falls back to top GWAS SNPs from the top germline regulatory elements.
#' @param pmid_query Deprecated. This argument is currently ignored; locus SNP
#'   labeling now uses dbSNP citation support directly.
#' @param pmid_page_size Maximum number of PMIDs retained per queried entity
#'   during built-in literature lookups. For locus SNP labels this applies to
#'   the dbSNP rsID citation lookup; for validated locus plots it also bounds
#'   the regulatory-element literature screening queries.
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
#' - `rsid_pmid`: an optional literature-support table with at least
#'   `rsid` and `pmid` columns
#' - `label_top_lit_snps`: number of literature-backed SNPs to label before
#'   falling back to top GWAS SNPs in the top germline regulatory elements
#'
#' Track semantics:
#' - the top three rows show regulatory-element somatic, epigenomic, and
#'   germline input scores from the scored graph
#' - the `Reg elements` row shows regulatory elements colored by their combined
#'   pre-diffusion norm
#' - the bottom gene row shows post-diffusion `conseguiR` scores for genes in
#'   the locus
#' - thin black diagonal curves connect regulatory elements to genes
#' - locus SNP labels prefer literature-backed SNPs when available and
#'   otherwise fall back to the top GWAS SNP in the plotted locus
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
  verbose = TRUE
) {
  .conseguiR_call_external("plot_locus_context", as.list(environment()))
}

#' Plot a validated locus-centered multimodal context panel
#'
#' Creates a paper-oriented locus panel that keeps only regulatory elements in
#' the requested window that have literature support for the regulatory
#' elements themselves. This is a citation-filtered sibling of
#' [plot_locus_context()] intended to reduce clutter and focus attention on
#' biologically validated regulatory elements and their linked genes.
#'
#' @inheritParams plot_locus_context
#' @param strict_gene_filter Logical scalar. In ordinary use, keep this as
#'   `TRUE` so the gene row is restricted to genes linked to
#'   literature-supported regulatory elements. Set to `FALSE` to keep all
#'   in-window genes while still filtering the regulatory-element rows to
#'   literature-supported elements.
#'
#' @details
#' Track semantics:
#' - the top three rows show somatic, epigenomic, and germline regulatory input
#'   scores, but only for literature-supported regulatory elements in the locus
#' - the `Reg elements` row shows the same literature-supported regulatory
#'   elements colored by their combined score
#' - the bottom gene row shows linked post-diffusion gene scores
#' - locus SNP labels can still use top GWAS SNPs or dbSNP-backed
#'   literature-supported SNPs
#'
#' Regulatory-element support is assessed by querying NCBI literature search for
#' the regulatory elements themselves using their interval coordinates,
#' identifiers, and linked-gene labels when available. This function errors
#' when no literature-supported regulatory elements are found in the requested
#' interval.
#'
#' @examples
#' names(formals(plot_validated_locus_context))
#'
#' @return A plot bundle containing the ggplot object and validated locus
#'   plotting data.
#' @export
plot_validated_locus_context <- function(
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
  strict_gene_filter = TRUE,
  plot_file_path = NULL,
  title = NULL,
  width = 14,
  height = 9,
  dpi = 300,
  save_plot = !is.null(plot_file_path),
  verbose = TRUE
) {
  .conseguiR_call_external("plot_validated_locus_context", as.list(environment()))
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
  verbose = TRUE
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
#'   When `NULL`, `run_conseguiR()` returns the stage bundles directly in
#'   memory without treating disk output as the default interface.
#' @param target_genes Requested selected-subgraph size.
#' @param candidate_pool_size Number of top diffusion-ranked candidate genes
#'   handed to the selected-subgraph solver. Must be at least as large as
#'   `target_genes`, and must not exceed the number of candidate genes
#'   available after diffusion. Larger values expand the solver search space
#'   and can increase runtime.
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
#' - `candidate_pool_size` controls how many top diffusion-ranked genes are
#'   considered by the selected-subgraph optimization stage before the final
#'   `target_genes` solution is chosen
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
#' The regulatory-element universe used across germline, somatic, and
#' epigenomic scoring is backend-managed by the package by default. In the
#' common workflow, users do not need to provide `reg_ref_path` manually.
#'
#' `conseguiR` returns R objects and bundles by default, and file writing is
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
#'   gene_step1_args = list(annotation_window = c(15, 15), ignore_strand = TRUE),
#'   gene_step2_args = list(
#'     gene_model = \"snp-wise=mean\",
#'     pval = list(use = c(\"SNP\", \"P\"), duplicate = \"drop\")
#'   ),
#'   reg_step1_args = list(annotation_window = c(15, 15)),
#'   reg_step2_args = list(
#'     pval = list(use = c(\"SNP\", \"P\"), duplicate = \"drop\")
#'   ),
#'   shared_args = list()
#' )`
#'
#' `somatic_args = list(
#'   gene_cv = NULL,
#'   sm = \"192r_3w\",
#'   kc = \"cgc81\",
#'   eligible_gr = NULL,
#'   fishhook_covariate_data = covariate_dt,
#'   constructor_use_local_mut_density = FALSE,
#'   score_iter = 25L,
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
  reg_ref_path = NULL,
  reference_bfile,
  dndscv_refdb,
  epigenomic_track_dir = NULL,
  epigenomic_tracks = NULL,
  graph_rds_path = NULL,
  gg_nodes_path = NULL,
  gg_edges_path = NULL,
  output_dir = NULL,
  target_genes = 50L,
  candidate_pool_size = 400L,
  germline_args = list(),
  somatic_args = list(),
  epigenomic_args = list(),
  scored_graph_args = list(),
  diffusion_args = list(),
  subgraph_args = list(),
  plot_args = list(),
  verbose = TRUE
) {
  reg_ref_path <- reg_ref_path %||% .conseguiR_default_reg_loc_path()
  if (is.null(reg_ref_path)) {
    stop("No regulatory reference was provided and no backend ENCODE-derived regulatory resource could be found.")
  }

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
    candidate_pool_size = candidate_pool_size,
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
