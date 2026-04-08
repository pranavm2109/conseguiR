#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

conseguiR_runtime_file <- function(relpath) {
  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(candidate)
  }

  pkg_path <- system.file(relpath, package = "conseguiR")
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }

  stop("Could not locate required runtime file: ", relpath)
}

sys.source(conseguiR_runtime_file("scripts/Internals/R/00_harmonise_and_validate_inputs.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/03_prepare_germline_scores.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/04_prepare_somatic_scores.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/05_prepare_epigenomic_scores.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/06_impose_scores_on_gene_reg_graph_nodes.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/07_run_diffusion_on_gene_reg_graph.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/08_call_subgraph.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/09_create_selected_subgraph_visualisation_bundle.R"), envir = environment())
sys.source(conseguiR_runtime_file("scripts/Internals/R/10_plot_stage_outputs.R"), envir = environment())

internal_run_gene_reg_diffusion <- run_gene_reg_diffusion

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

as_list_or_empty <- function(x) {
  if (is.null(x)) {
    return(list())
  }

  if (!is.list(x)) {
    stop("Expected a list of arguments, but received: ", class(x)[[1]])
  }

  x
}

read_table_if_path <- function(x, ...) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.character(x) && length(x) == 1L) {
    return(as.data.table(fread(x, ...)))
  }

  as.data.table(x)
}

ensure_output_parent <- function(path) {
  dir.create(dirname(path), recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

write_bundle_table <- function(table, path = NULL) {
  if (is.null(path)) {
    return(NULL)
  }

  ensure_output_parent(path)
  fwrite(as.data.table(table), path, sep = "\t")
  path
}

make_temp_output_prefix <- function(stem) {
  file.path(
    tempdir(),
    paste0(
      "conseguiR_",
      stem,
      "_",
      format(Sys.time(), "%Y%m%d%H%M%S"),
      "_",
      sprintf("%06d", sample.int(999999L, 1L))
    )
  )
}

make_temp_output_dir <- function(stem) {
  path <- file.path(
    tempdir(),
    paste0(
      "conseguiR_",
      stem,
      "_",
      format(Sys.time(), "%Y%m%d%H%M%S"),
      "_",
      sprintf("%06d", sample.int(999999L, 1L))
    )
  )
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  path
}

materialize_table_path <- function(table, path = NULL, stem = "table") {
  if (is.null(table)) {
    return(path)
  }

  path <- path %||% paste0(make_temp_output_prefix(stem), ".tsv")
  ensure_output_parent(path)
  fwrite(as.data.table(table), path, sep = "\t")
  path
}

new_bundle <- function(type, objects = list(), output_paths = list(), config = list()) {
  structure(
    c(
      list(
        bundle_type = type,
        objects = objects,
        output_paths = output_paths,
        config = config
      ),
      objects
    ),
    class = c(paste0("conseguiR_", type, "_bundle"), "conseguiR_bundle", "list")
  )
}

resolve_bundle_component <- function(x, preferred_name) {
  if (is.null(x)) {
    return(NULL)
  }

  if (is.data.frame(x) || data.table::is.data.table(x)) {
    return(as.data.table(x))
  }

  if (is.list(x) && preferred_name %in% names(x)) {
    return(as.data.table(x[[preferred_name]]))
  }

  if (is.list(x) && "objects" %in% names(x) && preferred_name %in% names(x$objects)) {
    return(as.data.table(x$objects[[preferred_name]]))
  }

  NULL
}

resolve_output_path <- function(bundle, preferred_name) {
  if (is.null(bundle) || !is.list(bundle)) {
    return(NULL)
  }

  if ("output_paths" %in% names(bundle) && preferred_name %in% names(bundle$output_paths)) {
    return(bundle$output_paths[[preferred_name]])
  }

  if (preferred_name %in% names(bundle)) {
    return(bundle[[preferred_name]])
  }

  NULL
}

#' Run a conseguiR external function with merged argument lists
#'
#' @param fun Function to call.
#' @param base_args Named list of base arguments.
#' @param extra_args Named list of additional arguments.
#'
#' @return The result of `fun`.
run_with_args <- function(fun, base_args = list(), extra_args = list()) {
  do.call(fun, c(base_args, extra_args))
}

verbose_message <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(...)
  }
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
#' @return A validation bundle containing validated objects and config.
validate_inputs <- function(
  gwas_sumstats = NULL,
  somatic_maf = NULL,
  reg_ref_path = NULL,
  epigenomic_tracks = NULL,
  epigenomic_track_dir = NULL,
  verbose = FALSE
) {
  verbose_message(verbose, "Validating inputs...")
  objects <- list()

  if (!is.null(gwas_sumstats)) {
    verbose_message(verbose, "Validating GWAS summary statistics...")
    objects$gwas <- if (is.character(gwas_sumstats) && length(gwas_sumstats) == 1L) {
      validate_gwas_sumstats_path(
        sumstats_path = gwas_sumstats,
        show_progress = isTRUE(verbose)
      )
    } else {
      validate_gwas_sumstats(gwas_sumstats)
    }
  }

  if (!is.null(somatic_maf)) {
    verbose_message(verbose, "Validating somatic MAF...")
    objects$somatic_maf <- if (is.character(somatic_maf) && length(somatic_maf) == 1L) {
      validate_somatic_maf_path(
        maf_path = somatic_maf,
        show_progress = isTRUE(verbose)
      )
    } else {
      validate_somatic_maf(somatic_maf)
    }
  }

  if (!is.null(reg_ref_path)) {
    verbose_message(verbose, "Validating regulatory-element reference...")
    objects$regulatory_elements <- validate_regulatory_element_reference(
      reg_ref_path,
      nrows = 1000L
    )
  }

  if (!is.null(epigenomic_tracks) || !is.null(epigenomic_track_dir)) {
    verbose_message(verbose, "Validating epigenomic tracks...")
    bw_files <- epigenomic_tracks
    if (is.null(bw_files)) {
      bw_files <- list_epigenomic_track_files(track_dir = epigenomic_track_dir)
    }

    if (is.null(reg_ref_path)) {
      stop("`reg_ref_path` is required when validating epigenomic tracks.")
    }

    objects$epigenomic <- validate_epigenomic_inputs(
      bw_files = bw_files,
      reg_ref_path = reg_ref_path,
      verbose = verbose
    )
  }

  new_bundle(
    type = "validation",
    objects = objects,
    config = list(
      has_gwas = !is.null(gwas_sumstats),
      has_somatic_maf = !is.null(somatic_maf),
      has_reg_ref = !is.null(reg_ref_path),
      has_epigenomic = !is.null(epigenomic_tracks) || !is.null(epigenomic_track_dir)
    )
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
#' User-facing wrapper around the internal MAGMA gene-scoring pipeline.
#'
#' @param gwas_sumstats GWAS summary statistics path or table.
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
run_germline_gene_scoring <- function(
  gwas_sumstats,
  gene_loc_path,
  reference_bfile,
  output_prefix = NULL,
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
  step1_args <- as_list_or_empty(step1_args)
  step2_args <- as_list_or_empty(step2_args)
  requested_output_prefix <- output_prefix
  output_prefix <- output_prefix %||% make_temp_output_prefix("germline_gene_scores")
  zstat_output_path <- if (is.null(requested_output_prefix)) NULL else paste0(output_prefix, ".zstat.tsv")

  result <- run_with_args(
    run_magma_feature_scoring_pipeline,
    base_args = list(
      gwas_sumstats = gwas_sumstats,
      feature_loc_path = gene_loc_path,
      reference_bfile = reference_bfile,
      output_prefix = output_prefix,
      feature_type = "gene",
      sample_size = sample_size,
      sample_size_col = sample_size_col,
      magma_path = magma_path,
      magma_gwas_cache_prefix = magma_gwas_cache_prefix,
      reuse_existing_gwas_cache = reuse_existing_gwas_cache,
      reuse_existing_annotation = reuse_existing_annotation,
      reuse_existing_analysis = reuse_existing_analysis,
      keep_intermediates = keep_intermediates,
      verbose = verbose,
      annotation_window = step1_args$annotation_window %||% NULL,
      filter_path = step1_args$filter_path %||% NULL,
      ignore_strand = step1_args$ignore_strand %||% FALSE,
      nonhuman = step1_args$nonhuman %||% FALSE,
      step1_extra_args = step1_args$extra_args %||% character(),
      gene_model = step2_args$gene_model %||% NULL,
      genes_only = step2_args$genes_only %||% TRUE,
      pval_use = step2_args$pval_use %||% NULL,
      pval_duplicate = step2_args$pval_duplicate %||% NULL,
      bfile_synonyms = step2_args$bfile_synonyms %||% NULL,
      bfile_synonym_dup = step2_args$bfile_synonym_dup %||% NULL,
      step2_extra_args = step2_args$extra_args %||% character(),
      zstat_output_path = zstat_output_path
    )
  )

  new_bundle(
    type = "germline_gene_scores",
    objects = list(
      gene_scores = result$zstat,
      pipeline = result
    ),
    output_paths = list(
      gene_scores_path = zstat_output_path
    ),
    config = list(
      output_prefix = requested_output_prefix,
      sample_size = sample_size,
      sample_size_col = sample_size_col,
      step1_args = step1_args,
      step2_args = step2_args
    )
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
#' User-facing wrapper around the internal MAGMA regulatory scoring pipeline.
#'
#' @param gwas_sumstats GWAS summary statistics path or table.
#' @param reg_loc_path Regulatory-element location file for MAGMA step 1.
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
#' @param step1_args Named list of MAGMA step 1 arguments.
#' @param step2_args Named list of MAGMA step 2 arguments.
#'
#' @return A germline regulatory score bundle.
run_germline_regulatory_scoring <- function(
  gwas_sumstats,
  reg_loc_path,
  reference_bfile,
  output_prefix = NULL,
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
  step1_args <- as_list_or_empty(step1_args)
  step2_args <- as_list_or_empty(step2_args)
  requested_output_prefix <- output_prefix
  output_prefix <- output_prefix %||% make_temp_output_prefix("germline_reg_scores")
  zstat_output_path <- if (is.null(requested_output_prefix)) NULL else paste0(output_prefix, ".zstat.tsv")

  result <- run_with_args(
    run_magma_feature_scoring_pipeline,
    base_args = list(
      gwas_sumstats = gwas_sumstats,
      feature_loc_path = reg_loc_path,
      reference_bfile = reference_bfile,
      output_prefix = output_prefix,
      feature_type = "regulatory_element",
      sample_size = sample_size,
      sample_size_col = sample_size_col,
      magma_path = magma_path,
      magma_gwas_cache_prefix = magma_gwas_cache_prefix,
      reuse_existing_gwas_cache = reuse_existing_gwas_cache,
      reuse_existing_annotation = reuse_existing_annotation,
      reuse_existing_analysis = reuse_existing_analysis,
      keep_intermediates = keep_intermediates,
      verbose = verbose,
      annotation_window = step1_args$annotation_window %||% NULL,
      filter_path = step1_args$filter_path %||% NULL,
      ignore_strand = step1_args$ignore_strand %||% FALSE,
      nonhuman = step1_args$nonhuman %||% FALSE,
      step1_extra_args = step1_args$extra_args %||% character(),
      gene_model = step2_args$gene_model %||% NULL,
      genes_only = step2_args$genes_only %||% TRUE,
      pval_use = step2_args$pval_use %||% NULL,
      pval_duplicate = step2_args$pval_duplicate %||% NULL,
      bfile_synonyms = step2_args$bfile_synonyms %||% NULL,
      bfile_synonym_dup = step2_args$bfile_synonym_dup %||% NULL,
      step2_extra_args = step2_args$extra_args %||% character(),
      zstat_output_path = zstat_output_path
    )
  )

  new_bundle(
    type = "germline_regulatory_scores",
    objects = list(
      reg_scores = result$zstat,
      pipeline = result
    ),
    output_paths = list(
      reg_scores_path = zstat_output_path
    ),
    config = list(
      output_prefix = requested_output_prefix,
      sample_size = sample_size,
      sample_size_col = sample_size_col,
      step1_args = step1_args,
      step2_args = step2_args
    )
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
#' Runs both MAGMA germline scoring paths with separately customizable step 1
#' and step 2 argument lists for genes and regulatory elements.
#'
#' @param gwas_sumstats GWAS summary statistics path or table.
#' @param reference_bfile PLINK reference prefix for MAGMA step 2.
#' @param gene_output_prefix Output prefix for gene-level germline scores.
#' @param reg_output_prefix Output prefix for regulatory germline scores.
#' @param magma_gwas_cache_prefix Shared MAGMA GWAS cache prefix.
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
  shared_args <- as_list_or_empty(shared_args)
  verbose_message(verbose, "Preparing germline scores...")
  gene_loc_path <- if (exists(".conseguiR_default_gene_loc_path", inherits = TRUE)) .conseguiR_default_gene_loc_path() else "data/raw/NCBI38/NCBI38.gene.loc"
  reg_loc_path <- if (exists(".conseguiR_default_reg_loc_path", inherits = TRUE)) .conseguiR_default_reg_loc_path() else "data/raw/GeneHancer/2026-01-26_UCSC_all_unfiltered_reg_elements.loc"

  gene_bundle <- run_with_args(
    run_germline_gene_scoring,
    base_args = c(
      list(
        gwas_sumstats = gwas_sumstats,
        gene_loc_path = gene_loc_path,
        reference_bfile = reference_bfile,
        output_prefix = gene_output_prefix,
        sample_size = gene_sample_size,
        sample_size_col = gene_sample_size_col,
        magma_gwas_cache_prefix = magma_gwas_cache_prefix,
        verbose = verbose,
        step1_args = gene_step1_args,
        step2_args = gene_step2_args
      ),
      shared_args
    )
  )

  reg_bundle <- run_with_args(
    run_germline_regulatory_scoring,
    base_args = c(
      list(
        gwas_sumstats = gwas_sumstats,
        reg_loc_path = reg_loc_path,
        reference_bfile = reference_bfile,
        output_prefix = reg_output_prefix,
        sample_size = reg_sample_size %||% gene_sample_size,
        sample_size_col = reg_sample_size_col %||% gene_sample_size_col,
        magma_gwas_cache_prefix = magma_gwas_cache_prefix,
        verbose = verbose,
        step1_args = reg_step1_args,
        step2_args = reg_step2_args
      ),
      shared_args
    )
  )

  new_bundle(
    type = "germline_scores",
    objects = list(
      gene_scores = gene_bundle$gene_scores,
      reg_scores = reg_bundle$reg_scores,
      gene_result = gene_bundle,
      reg_result = reg_bundle
    ),
    output_paths = list(
      gene_scores_path = gene_bundle$output_paths$gene_scores_path,
      reg_scores_path = reg_bundle$output_paths$reg_scores_path,
      magma_gwas_cache_prefix = magma_gwas_cache_prefix
    ),
    config = list(
      gene_step1_args = gene_step1_args,
      gene_step2_args = gene_step2_args,
      reg_step1_args = reg_step1_args,
      reg_step2_args = reg_step2_args
    )
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
  dndscv_args <- as_list_or_empty(dndscv_args)
  verbose_message(verbose, "Preparing somatic gene scores...")
  maf <- read_table_if_path(maf, showProgress = FALSE)

  scores <- run_with_args(
    run_dndscv_gene_scoring,
    base_args = c(
      list(
        maf = maf,
        refdb = refdb,
        cv = cv,
        max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,
        max_coding_muts_per_sample = max_coding_muts_per_sample,
        verbose = verbose
      ),
      dndscv_args
    )
  )

  output_path <- write_bundle_table(scores, output_path)

  new_bundle(
    type = "somatic_gene_scores",
    objects = list(gene_scores = scores),
    output_paths = list(gene_scores_path = output_path),
    config = list(
      refdb = refdb,
      cv = cv,
      max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,
      max_coding_muts_per_sample = max_coding_muts_per_sample,
      dndscv_args = dndscv_args
    )
  )
}

#' Run fishHook somatic regulatory scoring
#'
#' This wrapper exposes the main fishHook configuration surface used by the
#' package workflow and leaves less common model tuning available through
#' `fishhook_args`.
#'
#' @param maf Somatic MAF path or table.
#' @param reg_ref_path Regulatory-element reference path.
#' @param output_path Optional output path for saved scores.
#' @param eligible_gr Optional fishHook eligible territory.
#' @param fishhook_covariates Optional fishHook covariate objects/specifications.
#' @param fishhook_covariate_data Optional tabular covariate data.
#' @param idcol Sample identifier column for fishHook.
#' @param fishhook_args Named list of additional fishHook scoring arguments.
#'
#' @return A somatic regulatory score bundle.
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
  fishhook_args <- as_list_or_empty(fishhook_args)
  verbose_message(verbose, "Preparing somatic regulatory scores...")
  maf <- read_table_if_path(maf, showProgress = FALSE)

  scores <- run_with_args(
    run_fishhook_reg_scoring,
    base_args = c(
      list(
        maf = maf,
        reg_ref_path = reg_ref_path,
        eligible_gr = eligible_gr,
        fishhook_covariates = fishhook_covariates,
        fishhook_covariate_data = fishhook_covariate_data,
        idcol = idcol,
        verbose = verbose
      ),
      fishhook_args
    )
  )

  output_path <- write_bundle_table(scores, output_path)

  new_bundle(
    type = "somatic_regulatory_scores",
    objects = list(reg_scores = scores),
    output_paths = list(reg_scores_path = output_path),
    config = list(
      reg_ref_path = reg_ref_path,
      idcol = idcol,
      fishhook_args = fishhook_args
    )
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
#' @param maf Somatic MAF path or table.
#' @param refdb dndscv reference database path.
#' @param reg_ref_path Regulatory-element reference path.
#' @param gene_output_path Output path for somatic gene scores.
#' @param reg_output_path Output path for somatic regulatory scores.
#' @param gene_cv Optional dndscv covariates.
#' @param gene_max_muts_per_gene_per_sample dndscv gene-level mutation cap.
#' @param gene_max_coding_muts_per_sample dndscv sample-level coding mutation cap.
#' @param dndscv_args Named list of additional dndscv arguments.
#' @param eligible_gr Optional fishHook eligible territory.
#' @param fishhook_covariates Optional fishHook covariate objects/specifications.
#' @param fishhook_covariate_data Optional tabular covariate data.
#' @param fishhook_idcol Sample identifier column for fishHook.
#' @param fishhook_args Named list of additional fishHook arguments.
#'
#' @return A somatic score bundle with gene and regulatory score tables.
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
  verbose_message(verbose, "Preparing somatic scores...")
  gene_bundle <- run_somatic_gene_scoring(
    maf = maf,
    refdb = refdb,
    output_path = gene_output_path,
    cv = gene_cv,
    max_muts_per_gene_per_sample = gene_max_muts_per_gene_per_sample,
    max_coding_muts_per_sample = gene_max_coding_muts_per_sample,
    dndscv_args = dndscv_args,
    verbose = verbose
  )

  reg_bundle <- run_somatic_regulatory_scoring(
    maf = maf,
    reg_ref_path = reg_ref_path,
    output_path = reg_output_path,
    eligible_gr = eligible_gr,
    fishhook_covariates = fishhook_covariates,
    fishhook_covariate_data = fishhook_covariate_data,
    idcol = fishhook_idcol,
    fishhook_args = fishhook_args,
    verbose = verbose
  )

  new_bundle(
    type = "somatic_scores",
    objects = list(
      gene_scores = gene_bundle$gene_scores,
      reg_scores = reg_bundle$reg_scores,
      gene_result = gene_bundle,
      reg_result = reg_bundle
    ),
    output_paths = list(
      gene_scores_path = gene_bundle$output_paths$gene_scores_path,
      reg_scores_path = reg_bundle$output_paths$reg_scores_path
    ),
    config = list(
      dndscv_args = dndscv_args,
      fishhook_args = fishhook_args
    )
  )
}

#' Prepare regulatory epigenomic scores
#'
#' Quantifies bigWig signal across regulatory elements and converts cross-track
#' variation into regulatory-element epigenomic z scores.
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
#'
#' @return An epigenomic score bundle.
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
  verbose_message(verbose, "Preparing epigenomic scores...")
  result <- run_epigenomic_reg_scoring(
    track_dir = track_dir,
    reg_ref_path = reg_ref_path,
    bw_files = bw_files,
    min_tracks = min_tracks,
    drop_mhc = drop_mhc,
    transform = transform,
    return_diagnostics = return_diagnostics,
    summary_fun = summary_fun,
    verbose = verbose
  )

  if (is.list(result) && "zscores" %in% names(result)) {
    reg_scores <- as.data.table(result$zscores)
    diagnostics <- result$diagnostics %||% NULL
  } else {
    reg_scores <- as.data.table(result)
    diagnostics <- NULL
  }

  output_path <- write_bundle_table(reg_scores, output_path)

  new_bundle(
    type = "epigenomic_scores",
    objects = list(
      reg_scores = reg_scores,
      diagnostics = diagnostics,
      raw_result = result
    ),
    output_paths = list(
      reg_scores_path = output_path
    ),
    config = list(
      reg_ref_path = reg_ref_path,
      track_dir = track_dir,
      bw_files = bw_files,
      min_tracks = min_tracks,
      transform = transform
    )
  )
}

#' Build a scored gene-regulatory graph
#'
#' Merges germline, somatic, and epigenomic score tables into the backend
#' no-score gene-regulatory graph.
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
#' @return A scored gene-reg graph bundle containing the graph, nodes, and
#'   edges.
build_scored_gene_reg_graph <- function(
  graph = NULL,
  graph_rds_path = "data/processed/gene_reg_graph_no_scores.rds",
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
  verbose_message(verbose, "Building scored gene-regulatory graph...")
  gene_germline_scores <- gene_germline_scores %||% resolve_bundle_component(germline_scores, "gene_scores")
  reg_germline_scores <- reg_germline_scores %||% resolve_bundle_component(germline_scores, "reg_scores")
  gene_somatic_scores <- gene_somatic_scores %||% resolve_bundle_component(somatic_scores, "gene_scores")
  reg_somatic_scores <- reg_somatic_scores %||% resolve_bundle_component(somatic_scores, "reg_scores")
  reg_epigenomic_scores <- reg_epigenomic_scores %||% resolve_bundle_component(epigenomic_scores, "reg_scores")

  result <- prepare_scored_gene_reg_graph(
    graph = graph,
    graph_rds_path = graph_rds_path,
    output_prefix = output_prefix,
    gene_somatic_scores = gene_somatic_scores,
    gene_germline_scores = gene_germline_scores,
    reg_somatic_scores = reg_somatic_scores,
    reg_germline_scores = reg_germline_scores,
    reg_epigenomic_scores = reg_epigenomic_scores,
    save_outputs = save_outputs
  )

  new_bundle(
    type = "scored_gene_reg_graph",
    objects = list(
      graph = result$graph,
      nodes = result$nodes,
      edges = result$edges
    ),
    output_paths = list(
      graph_rds_path = if (isTRUE(save_outputs)) paste0(output_prefix, ".rds") else NULL,
      nodes_path = if (isTRUE(save_outputs)) paste0(output_prefix, "_nodes.tsv.gz") else NULL,
      edges_path = if (isTRUE(save_outputs)) paste0(output_prefix, "_edges.tsv.gz") else NULL
    ),
    config = list(
      graph_rds_path = graph_rds_path,
      output_prefix = output_prefix
    )
  )
}

#' Run diffusion on a scored gene-regulatory graph
#'
#' User-facing wrapper for the Python-backed diffusion stage.
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
  verbose_message(verbose, "Running gene-regulatory diffusion...")
  requested_output_dir <- output_dir
  requested_output_stem <- output_stem
  nodes <- resolve_bundle_component(scored_graph, "nodes")
  edges <- resolve_bundle_component(scored_graph, "edges")
  nodes_path <- materialize_table_path(
    table = nodes,
    path = nodes_path %||% resolve_output_path(scored_graph, "nodes_path"),
    stem = "gene_reg_graph_scored_nodes"
  )
  edges_path <- materialize_table_path(
    table = edges,
    path = edges_path %||% resolve_output_path(scored_graph, "edges_path"),
    stem = "gene_reg_graph_scored_edges"
  )
  output_dir <- output_dir %||% make_temp_output_dir("gene_reg_diffusion")
  output_stem <- output_stem %||% "gene_reg_graph_diffusion"

  validate_scored_gene_reg_nodes(read_scored_gene_reg_nodes(nodes_path %||% default_diffusion_config$nodes_path))
  validate_scored_gene_reg_edges(read_scored_gene_reg_edges(edges_path %||% default_diffusion_config$edges_path))

  result <- internal_run_gene_reg_diffusion(
    nodes_path = nodes_path %||% default_diffusion_config$nodes_path,
    edges_path = edges_path %||% default_diffusion_config$edges_path,
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

  new_bundle(
    type = "diffusion",
    objects = list(
      all_genes = read_diffusion_results(result$output_paths$all_genes_path),
      top_genes = read_diffusion_results(result$output_paths$top_genes_path)
    ),
    output_paths = if (!is.null(requested_output_dir) || !is.null(requested_output_stem)) as.list(result$output_paths) else list(
      all_genes_path = NULL,
      top_genes_path = NULL
    ),
    config = c(
      result$config,
      list(
        output_dir = requested_output_dir,
        output_stem = requested_output_stem
      )
    )
  )
}

#' Call the selected subgraph from diffusion results
#'
#' User-facing wrapper for the Python-backed cardinality-constrained subgraph
#' selection stage.
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
  verbose_message(verbose, "Calling selected subgraph...")
  requested_output_dir <- output_dir
  requested_output_stem <- output_stem
  diffusion_table <- resolve_bundle_component(diffusion, "all_genes")
  diffusion_path <- materialize_table_path(
    table = diffusion_table,
    path = diffusion_path %||% resolve_output_path(diffusion, "all_genes_path"),
    stem = "gene_reg_graph_diffusion_all_genes"
  )
  if (is.null(gg_nodes_path) || is.null(gg_edges_path)) {
    initialize_backend_graphs(strict = FALSE, quiet = TRUE)
    backend_paths <- .conseguiR_backend_paths()
    gg_nodes_path <- gg_nodes_path %||% backend_paths$gene_gene_graph_nodes
    gg_edges_path <- gg_edges_path %||% backend_paths$gene_gene_graph_edges
  }
  output_dir <- output_dir %||% make_temp_output_dir("selected_subgraph")
  output_stem <- output_stem %||% "gene_gene_selected_subgraph"

  validate_diffusion_results(read_diffusion_results(diffusion_path %||% default_subgraph_config$diffusion_path), prize_column = prize_column)
  validate_gene_gene_nodes(read_gene_gene_nodes(gg_nodes_path))
  validate_gene_gene_edges(
    read_gene_gene_edges(gg_edges_path),
    confidence_column = confidence_column,
    edge_cost_column = edge_cost_column
  )

  result <- run_cardinality_subgraph_calling(
    diffusion_path = diffusion_path %||% default_subgraph_config$diffusion_path,
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

  new_bundle(
    type = "selected_subgraph",
    objects = list(
      nodes = read_selected_subgraph_nodes(result$output_paths$nodes_path),
      edges = read_selected_subgraph_edges(result$output_paths$edges_path),
      summary = read_selected_subgraph_summary(result$output_paths$summary_path)
    ),
    output_paths = if (!is.null(requested_output_dir) || !is.null(requested_output_stem)) as.list(result$output_paths) else list(
      nodes_path = NULL,
      edges_path = NULL,
      summary_path = NULL,
      graphml_path = NULL
    ),
    config = c(
      result$config,
      list(
        output_dir = requested_output_dir,
        output_stem = requested_output_stem
      )
    )
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
  verbose_message(verbose, "Preparing selected-subgraph plot...")
  nodes <- nodes %||% resolve_bundle_component(selected_subgraph, "nodes")
  edges <- edges %||% resolve_bundle_component(selected_subgraph, "edges")
  summary <- summary %||% resolve_bundle_component(selected_subgraph, "summary")
  nodes_path <- nodes_path %||% resolve_output_path(selected_subgraph, "nodes_path")
  edges_path <- edges_path %||% resolve_output_path(selected_subgraph, "edges_path")
  summary_path <- summary_path %||% resolve_output_path(selected_subgraph, "summary_path")

  bundle <- prepare_selected_subgraph_visualisation_bundle(
    nodes = nodes,
    edges = edges,
    summary = summary,
    nodes_path = nodes_path %||% default_selected_subgraph_plot_config$nodes_path,
    edges_path = edges_path %||% default_selected_subgraph_plot_config$edges_path,
    summary_path = summary_path %||% default_selected_subgraph_plot_config$summary_path
  )

  plot_obj <- create_selected_subgraph_plot(
    bundle = bundle,
    layout = layout,
    top_n_labels = top_n_labels,
    title = title,
    subtitle = NULL
  )

  if (isTRUE(save_bundle)) {
    save_selected_subgraph_visualisation_bundle(
      bundle = bundle,
      output_prefix = bundle_output_prefix
    )
  }

  if (isTRUE(save_plot)) {
    if (is.null(plot_file_path) || !nzchar(plot_file_path)) {
      stop("`plot_file_path` must be provided when `save_plot = TRUE`.")
    }

    save_selected_subgraph_plot(
      bundle = bundle,
      file_path = plot_file_path,
      width = width,
      height = height,
      dpi = dpi,
      layout = layout,
      top_n_labels = top_n_labels,
      title = title,
      subtitle = NULL
    )
  }

  new_bundle(
    type = "plot",
    objects = list(
      plot = plot_obj,
      bundle = bundle
    ),
    output_paths = list(
      bundle_rds_path = if (isTRUE(save_bundle)) paste0(bundle_output_prefix, ".rds") else NULL,
      bundle_nodes_path = if (isTRUE(save_bundle)) paste0(bundle_output_prefix, "_nodes.tsv.gz") else NULL,
      bundle_edges_path = if (isTRUE(save_bundle)) paste0(bundle_output_prefix, "_edges.tsv.gz") else NULL,
      bundle_summary_path = if (isTRUE(save_bundle)) paste0(bundle_output_prefix, "_summary.tsv") else NULL,
      plot_file_path = if (isTRUE(save_plot)) plot_file_path else NULL
    ),
    config = list(
      title = title,
      layout = layout,
      top_n_labels = top_n_labels
    )
  )
}

#' Plot score tables from germline, somatic, epigenomic, or diffusion bundles
#'
#' Creates a rank plot for one-tailed outputs and a volcano plot for two-tailed
#' outputs.
#'
#' @param scores Optional bundle containing a plottable score table.
#' @param table Optional explicit score table.
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
#'
#' @return A plot bundle containing the ggplot object and plotting data.
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
  verbose_message(verbose, "Preparing score plot...")

  plot_obj <- create_score_plot(
    bundle = scores,
    table = table,
    which = which,
    test_tail = test_tail,
    feature_column = feature_column,
    z_column = z_column,
    p_value_column = p_value_column,
    label_features = label_features,
    title = title
  )

  export_width <- width
  export_height <- height
  if (identical(plot_obj$plot_mode, "rank") && nrow(plot_obj$plot_data) > 100000L) {
    export_width <- max(export_width, 12)
    export_height <- max(export_height, 8)
  }

  if (isTRUE(save_plot)) {
    if (is.null(plot_file_path) || !nzchar(plot_file_path)) {
      stop("`plot_file_path` must be provided when `save_plot = TRUE`.")
    }

    save_score_plot(
      bundle = scores,
      table = table,
      file_path = plot_file_path,
      which = which,
      test_tail = test_tail,
      feature_column = feature_column,
      z_column = z_column,
      p_value_column = p_value_column,
      label_features = label_features,
      title = title,
      width = export_width,
      height = export_height,
      dpi = dpi
    )
  }

  invisible(new_bundle(
    type = "score_plot",
    objects = list(
      plot = plot_obj$plot,
      plot_data = plot_obj$plot_data
    ),
    output_paths = list(
      plot_file_path = if (isTRUE(save_plot)) plot_file_path else NULL
    ),
    config = list(
      which = which,
      test_tail = test_tail,
      feature_column = feature_column,
      z_column = z_column,
      p_value_column = p_value_column,
      plot_mode = plot_obj$plot_mode,
      export_width = export_width,
      export_height = export_height
    )
  ))
}

resolve_plot_nodes <- function(scored_graph = NULL, nodes_path = NULL) {
  nodes <- resolve_bundle_component(scored_graph, "nodes")
  if (!is.null(nodes)) {
    return(as.data.table(nodes))
  }

  nodes_path <- nodes_path %||% resolve_output_path(scored_graph, "nodes_path")
  nodes_path <- nodes_path %||% default_diffusion_config$nodes_path
  validate_scored_gene_reg_nodes(read_scored_gene_reg_nodes(nodes_path))
}

resolve_plot_diffusion <- function(diffusion = NULL, diffusion_path = NULL, which = "all_genes") {
  dt <- resolve_bundle_component(diffusion, which)
  if (!is.null(dt)) {
    return(data.table::as.data.table(dt))
  }

  diffusion_path <- diffusion_path %||% resolve_output_path(diffusion, paste0(which, "_path"))
  diffusion_path <- diffusion_path %||% resolve_output_path(diffusion, "all_genes_path")
  if (is.null(diffusion_path) || !nzchar(diffusion_path)) {
    stop("A diffusion bundle or diffusion_path is required for post-diffusion plotting.")
  }

  read_diffusion_results(diffusion_path)
}

resolve_plot_edges <- function(scored_graph = NULL, edges_path = NULL) {
  edges <- resolve_bundle_component(scored_graph, "edges")
  if (!is.null(edges)) {
    return(as.data.table(edges))
  }

  edges_path <- edges_path %||% resolve_output_path(scored_graph, "edges_path")
  edges_path <- edges_path %||% default_diffusion_config$edges_path
  validate_scored_gene_reg_edges(read_scored_gene_reg_edges(edges_path))
}

merge_postdiff_scores_into_nodes <- function(nodes_dt, diffusion_dt = NULL) {
  nodes_dt <- as.data.table(copy(nodes_dt))

  for (col in c("post_germline", "post_somatic", "post_epigenomic", "post_norm")) {
    if (!col %in% names(nodes_dt)) {
      nodes_dt[, (col) := NA_real_]
    }
  }

  if (is.null(diffusion_dt) || nrow(diffusion_dt) == 0L) {
    return(nodes_dt)
  }

  diffusion_dt <- as.data.table(copy(diffusion_dt))
  diffusion_dt[, gene_name := as.character(gene_name)]
  keep_cols <- intersect(c("gene_name", "post_germline", "post_somatic", "post_epigenomic", "post_norm"), names(diffusion_dt))
  if (!"gene_name" %in% keep_cols) {
    return(nodes_dt)
  }

  gene_map <- nodes_dt[node_type == "gene", .(node_id = as.character(node_id), gene_name = as.character(name))]
  gene_map <- merge(gene_map, diffusion_dt[, ..keep_cols], by = "gene_name", all.x = TRUE)

  for (col in setdiff(keep_cols, "gene_name")) {
    nodes_dt[match(gene_map$node_id, as.character(node_id)), (col) := safe_numeric(gene_map[[col]])]
  }

  nodes_dt
}

resolve_postdiff_gene_reg_graph <- function(
  postdiff_gene_reg_graph = NULL,
  scored_graph = NULL,
  diffusion = NULL,
  nodes_path = NULL,
  edges_path = NULL,
  diffusion_path = NULL
) {
  if (is.character(postdiff_gene_reg_graph) && length(postdiff_gene_reg_graph) == 1L && file.exists(postdiff_gene_reg_graph)) {
    postdiff_gene_reg_graph <- readRDS(postdiff_gene_reg_graph)
  }

  nodes_dt <- resolve_bundle_component(postdiff_gene_reg_graph, "nodes")
  edges_dt <- resolve_bundle_component(postdiff_gene_reg_graph, "edges")

  if (!is.null(nodes_dt) && !is.null(edges_dt)) {
    nodes_dt <- as.data.table(nodes_dt)
    edges_dt <- as.data.table(edges_dt)
    if (!"post_norm" %in% names(nodes_dt)) {
      diffusion_dt <- resolve_bundle_component(postdiff_gene_reg_graph, "all_genes")
      nodes_dt <- merge_postdiff_scores_into_nodes(nodes_dt, diffusion_dt)
    }
    return(list(nodes = nodes_dt, edges = edges_dt))
  }

  nodes_dt <- resolve_plot_nodes(scored_graph = scored_graph, nodes_path = nodes_path)
  edges_dt <- resolve_plot_edges(scored_graph = scored_graph, edges_path = edges_path)
  diffusion_dt <- resolve_plot_diffusion(diffusion = diffusion, diffusion_path = diffusion_path, which = "all_genes")

  list(
    nodes = merge_postdiff_scores_into_nodes(nodes_dt, diffusion_dt),
    edges = edges_dt
  )
}

default_score_table_path <- function(modality = c("germline", "somatic", "epigenomic"), target = c("gene", "reg")) {
  modality <- match.arg(modality)
  target <- match.arg(target)

  switch(
    paste(modality, target, sep = "_"),
    germline_gene = "data/processed/germline_gene_scores.tsv",
    germline_reg = "data/processed/germline_reg_scores.tsv",
    somatic_gene = "data/processed/somatic_gene_scores.tsv",
    somatic_reg = "data/processed/somatic_reg_scores.tsv",
    epigenomic_reg = "data/processed/epigenomic_reg_scores.tsv",
    NULL
  )
}

build_prediff_gene_table <- function(score_bundle = NULL, scored_graph = NULL, nodes_path = NULL, modality = c("germline", "somatic", "epigenomic")) {
  modality <- match.arg(modality)

  preferred_component <- if (identical(modality, "epigenomic")) NULL else "gene_scores"
  if (!is.null(preferred_component)) {
    dt <- resolve_bundle_component(score_bundle, preferred_component)
    if (!is.null(dt)) {
      return(list(
        table = data.table::as.data.table(dt),
        feature_column = intersect(c("gene_name", "gene_id", "feature_id"), names(dt))[[1]],
        z_column = "zstat",
        p_value_column = if ("p_value" %in% names(dt)) "p_value" else NULL
      ))
    }
  }

  default_path <- default_score_table_path(modality = modality, target = "gene")
  if (!is.null(default_path) && file.exists(default_path)) {
    dt <- data.table::as.data.table(data.table::fread(default_path, showProgress = FALSE))
    return(list(
      table = dt,
      feature_column = intersect(c("gene_name", "gene_id", "feature_id"), names(dt))[[1]],
      z_column = "zstat",
      p_value_column = if ("p_value" %in% names(dt)) "p_value" else NULL
    ))
  }

  nodes <- resolve_plot_nodes(scored_graph = scored_graph, nodes_path = nodes_path)
  score_col <- paste0(modality, "_score")
  if (!score_col %in% names(nodes)) {
    stop("Scored node table is missing expected column: ", score_col)
  }

  gene_nodes <- nodes[node_type == "gene"]
  gene_label_col <- if ("name" %in% names(gene_nodes)) "name" else "node_id"
  list(
    table = gene_nodes[, .(
      gene_name = as.character(get(gene_label_col)),
      zstat = as.numeric(get(score_col))
    )],
    feature_column = "gene_name",
    z_column = "zstat",
    p_value_column = NULL
  )
}

build_prediff_reg_table <- function(score_bundle = NULL, scored_graph = NULL, nodes_path = NULL, modality = c("germline", "somatic", "epigenomic")) {
  modality <- match.arg(modality)

  dt <- resolve_bundle_component(score_bundle, "reg_scores")
  if (!is.null(dt)) {
    feature_column <- intersect(c("reg_elem_id", "feature_id"), names(dt))[[1]]
    return(list(
      table = data.table::as.data.table(dt),
      feature_column = feature_column,
      z_column = "zstat",
      p_value_column = if ("p_value" %in% names(dt)) "p_value" else NULL
    ))
  }

  default_path <- default_score_table_path(modality = modality, target = "reg")
  if (!is.null(default_path) && file.exists(default_path)) {
    dt <- data.table::as.data.table(data.table::fread(default_path, showProgress = FALSE))
    feature_column <- intersect(c("reg_elem_id", "feature_id"), names(dt))[[1]]
    return(list(
      table = dt,
      feature_column = feature_column,
      z_column = "zstat",
      p_value_column = if ("p_value" %in% names(dt)) "p_value" else NULL
    ))
  }

  nodes <- resolve_plot_nodes(scored_graph = scored_graph, nodes_path = nodes_path)
  score_col <- paste0(modality, "_score")
  if (!score_col %in% names(nodes)) {
    stop("Scored node table is missing expected column: ", score_col)
  }

  reg_nodes <- nodes[node_type == "reg"]
  list(
    table = reg_nodes[, .(
      reg_elem_id = as.character(node_id),
      zstat = as.numeric(get(score_col))
    )],
    feature_column = "reg_elem_id",
    z_column = "zstat",
    p_value_column = NULL
  )
}

build_postdiff_gene_table <- function(diffusion = NULL, diffusion_path = NULL, modality = c("germline", "somatic", "epigenomic")) {
  modality <- match.arg(modality)
  dt <- resolve_plot_diffusion(diffusion = diffusion, diffusion_path = diffusion_path, which = "all_genes")
  z_column <- paste0("post_", modality)
  if (!z_column %in% names(dt)) {
    stop("Diffusion table is missing expected post-diffusion column: ", z_column)
  }

  feature_column <- if ("gene_name" %in% names(dt)) "gene_name" else "node_id"
  list(
    table = data.table::as.data.table(dt),
    feature_column = feature_column,
    z_column = z_column,
    p_value_column = NULL
  )
}

run_score_plot_spec <- function(spec, plot_file_path = NULL, test_tail = "one_tailed", label_features = NULL, title = "conseguiR Scores", width = 10, height = 7, dpi = 300, save_plot = !is.null(plot_file_path), verbose = FALSE) {
  plot_scores(
    table = spec$table,
    feature_column = spec$feature_column,
    z_column = spec$z_column,
    p_value_column = spec$p_value_column,
    plot_file_path = plot_file_path,
    test_tail = test_tail,
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
  stage <- match.arg(stage)
  spec <- if (identical(stage, "pre")) {
    build_prediff_gene_table(germline_scores, scored_graph, nodes_path, modality = "germline")
  } else {
    build_postdiff_gene_table(diffusion, diffusion_path, modality = "germline")
  }

  title <- title %||% if (identical(stage, "pre")) "Pre-diffusion Germline Gene Signal" else "Post-diffusion Germline Gene Signal"
  run_score_plot_spec(spec, plot_file_path, test_tail = "one_tailed", label_features = label_features, title = title, width = width, height = height, dpi = dpi, save_plot = save_plot, verbose = verbose)
}

#' Plot germline regulatory signal
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
  spec <- build_prediff_reg_table(germline_scores, scored_graph, nodes_path, modality = "germline")
  title <- title %||% "Germline Regulatory Signal"
  run_score_plot_spec(spec, plot_file_path, test_tail = "one_tailed", label_features = label_features, title = title, width = width, height = height, dpi = dpi, save_plot = save_plot, verbose = verbose)
}

#' Plot somatic gene signal before or after diffusion
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
  stage <- match.arg(stage)
  spec <- if (identical(stage, "pre")) {
    build_prediff_gene_table(somatic_scores, scored_graph, nodes_path, modality = "somatic")
  } else {
    build_postdiff_gene_table(diffusion, diffusion_path, modality = "somatic")
  }

  title <- title %||% if (identical(stage, "pre")) "Pre-diffusion Somatic Gene Signal" else "Post-diffusion Somatic Gene Signal"
  run_score_plot_spec(spec, plot_file_path, test_tail = "two_tailed", label_features = label_features, title = title, width = width, height = height, dpi = dpi, save_plot = save_plot, verbose = verbose)
}

#' Plot somatic regulatory signal
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
  spec <- build_prediff_reg_table(somatic_scores, scored_graph, nodes_path, modality = "somatic")
  title <- title %||% "Somatic Regulatory Signal"
  run_score_plot_spec(spec, plot_file_path, test_tail = "two_tailed", label_features = label_features, title = title, width = width, height = height, dpi = dpi, save_plot = save_plot, verbose = verbose)
}

#' Plot post-diffusion epigenomic gene signal
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
  spec <- build_postdiff_gene_table(diffusion, diffusion_path, modality = "epigenomic")
  title <- title %||% "Post-diffusion Epigenomic Gene Signal"
  run_score_plot_spec(spec, plot_file_path, test_tail = "one_tailed", label_features = label_features, title = title, width = width, height = height, dpi = dpi, save_plot = save_plot, verbose = verbose)
}

#' Plot epigenomic regulatory signal
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
  spec <- build_prediff_reg_table(epigenomic_scores, scored_graph, nodes_path, modality = "epigenomic")
  title <- title %||% "Epigenomic Regulatory Signal"
  run_score_plot_spec(spec, plot_file_path, test_tail = "one_tailed", label_features = label_features, title = title, width = width, height = height, dpi = dpi, save_plot = save_plot, verbose = verbose)
}

#' Run the full conseguiR pipeline end to end
#'
#' Orchestrates validation, score preparation, scored-graph construction,
#' diffusion, subgraph calling, and plotting.
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
  verbose_message(verbose, "Running end-to-end conseguiR pipeline...")
  requested_output_dir <- output_dir
  if (is.null(graph_rds_path) || is.null(gg_nodes_path) || is.null(gg_edges_path)) {
    initialize_backend_graphs(strict = FALSE, quiet = TRUE)
    backend_paths <- .conseguiR_backend_paths()
    graph_rds_path <- graph_rds_path %||% backend_paths$gene_reg_graph_rds
    gg_nodes_path <- gg_nodes_path %||% backend_paths$gene_gene_graph_nodes
    gg_edges_path <- gg_edges_path %||% backend_paths$gene_gene_graph_edges
  }
  validation <- validate_inputs(
    gwas_sumstats = gwas_sumstats,
    somatic_maf = somatic_maf,
    reg_ref_path = reg_ref_path,
    epigenomic_tracks = epigenomic_tracks,
    epigenomic_track_dir = epigenomic_track_dir,
    verbose = verbose
  )

  germline <- run_with_args(
    prepare_germline_scores,
    base_args = c(
      list(
        gwas_sumstats = gwas_sumstats,
        reference_bfile = reference_bfile,
        gene_output_prefix = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "germline_gene_scores") else NULL,
        reg_output_prefix = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "germline_reg_scores") else NULL,
        magma_gwas_cache_prefix = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "magma_shared_gwas_cache") else NULL,
        verbose = verbose
      ),
      as_list_or_empty(germline_args)
    )
  )

  somatic <- run_with_args(
    prepare_somatic_scores,
    base_args = c(
      list(
        maf = somatic_maf,
        refdb = dndscv_refdb,
        reg_ref_path = reg_ref_path,
        gene_output_path = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "somatic_gene_scores.tsv") else NULL,
        reg_output_path = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "somatic_reg_scores.tsv") else NULL,
        verbose = verbose
      ),
      as_list_or_empty(somatic_args)
    )
  )

  epigenomic <- run_with_args(
    prepare_epigenomic_scores,
    base_args = c(
      list(
        reg_ref_path = reg_ref_path,
        track_dir = epigenomic_track_dir,
        bw_files = epigenomic_tracks,
        output_path = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "epigenomic_reg_scores.tsv") else NULL,
        verbose = verbose
      ),
      as_list_or_empty(epigenomic_args)
    )
  )

  scored_graph <- run_with_args(
    build_scored_gene_reg_graph,
    base_args = c(
      list(
        graph_rds_path = graph_rds_path,
        output_prefix = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "gene_reg_graph_scored") else "data/processed/gene_reg_graph_scored",
        germline_scores = germline,
        somatic_scores = somatic,
        epigenomic_scores = epigenomic,
        save_outputs = !is.null(requested_output_dir),
        verbose = verbose
      ),
      as_list_or_empty(scored_graph_args)
    )
  )

  diffusion <- run_with_args(
    run_gene_reg_diffusion,
    base_args = c(
      list(
        scored_graph = scored_graph,
        output_dir = requested_output_dir,
        output_stem = if (!is.null(requested_output_dir)) "gene_reg_graph_diffusion" else NULL,
        verbose = verbose
      ),
      as_list_or_empty(diffusion_args)
    )
  )

  selected_subgraph <- run_with_args(
    call_selected_subgraph,
    base_args = c(
      list(
        diffusion = diffusion,
        gg_nodes_path = gg_nodes_path,
        gg_edges_path = gg_edges_path,
        output_dir = requested_output_dir,
        output_stem = if (!is.null(requested_output_dir)) "gene_gene_selected_subgraph" else NULL,
        target_genes = target_genes,
        verbose = verbose
      ),
      as_list_or_empty(subgraph_args)
    )
  )

  plot_bundle <- run_with_args(
    plot_selected_subgraph,
    base_args = c(
      list(
        selected_subgraph = selected_subgraph,
        bundle_output_prefix = if (!is.null(requested_output_dir)) file.path(requested_output_dir, "gene_gene_selected_subgraph_plot_bundle") else "data/processed/gene_gene_selected_subgraph_plot_bundle",
        plot_file_path = NULL,
        save_bundle = FALSE,
        save_plot = FALSE,
        verbose = verbose
      ),
      as_list_or_empty(plot_args)
    )
  )

  new_bundle(
    type = "pipeline",
    objects = list(
      validation = validation,
      germline = germline,
      somatic = somatic,
      epigenomic = epigenomic,
      scored_graph = scored_graph,
      diffusion = diffusion,
      selected_subgraph = selected_subgraph,
      plot = plot_bundle
    ),
    output_paths = list(
      output_dir = requested_output_dir
    ),
    config = list(
      target_genes = target_genes
    )
  )
}

#' Plot a locus-centered multimodal context panel
#'
#' @param chromosome Locus chromosome, for example `"8"` or `"chr8"`.
#' @param start Locus start coordinate.
#' @param end Locus end coordinate.
#' @param scored_graph Optional scored gene-reg graph bundle.
#' @param diffusion Optional diffusion bundle.
#' @param selected_subgraph Optional selected-subgraph bundle used to highlight
#'   genes and filter the arc layer.
#' @param nodes_path Optional explicit scored node-table path.
#' @param edges_path Optional explicit scored edge-table path.
#' @param diffusion_path Optional explicit diffusion table path.
#' @param selected_nodes_path Optional explicit selected-subgraph node-table
#'   path.
#' @param label_features Optional gene symbols to label.
#' @param gwas_sumstats Optional GWAS summary statistics path or table used to
#'   label the top SNP in the locus.
#' @param label_top_gwas_snp Logical scalar. If `TRUE`, label the top GWAS SNP
#'   in the window by rsID/variant ID.
#' @param rsid_pmid Optional rsID-to-PMID evidence table or path. Must contain
#'   at least `rsid` and `pmid` columns. When omitted, `conseguiR` queries
#'   LitVar in memory for SNP literature support.
#' @param label_top_lit_snps Integer count of literature-backed SNPs to label.
#'   SNPs are restricted to regulatory elements in the locus and ranked by GWAS
#'   significance, with PMID support used as an additional prioritization field.
#' @param pmid_query Optional disease query term such as `"DLBCL"` or
#'   `"lymphoma"`. When supplied, `conseguiR` first queries LitVar for rsID to
#'   PMID links and then filters those PMIDs using Europe PMC query results. If
#'   no literature-backed SNPs are found, the plot falls back to top GWAS SNPs
#'   inside the top regulatory elements by germline z-score.
#' @param pmid_page_size Maximum number of PMIDs per rsID retained from LitVar,
#'   and the Europe PMC page size used when filtering by `pmid_query`.
#' @param plot_file_path Optional output path for the saved figure.
#' @param title Plot title.
#' @param width Plot width in inches.
#' @param height Plot height in inches.
#' @param dpi Plot DPI.
#' @param save_plot Whether to save the figure.
#' @param verbose Logical scalar. If `TRUE`, show stage messages.
#'
#' @return A plot bundle containing the ggplot object and locus plotting data.
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
  verbose_message(verbose, "Preparing locus context plot...")

  graph_dt <- resolve_postdiff_gene_reg_graph(
    postdiff_gene_reg_graph = postdiff_gene_reg_graph,
    scored_graph = scored_graph,
    diffusion = diffusion,
    nodes_path = nodes_path,
    edges_path = edges_path,
    diffusion_path = diffusion_path
  )
  nodes_dt <- graph_dt$nodes
  edges_dt <- graph_dt$edges
  selected_genes <- resolve_locus_selected_genes(
    selected_subgraph = selected_subgraph,
    nodes_path = selected_nodes_path
  )

  plot_obj <- create_locus_context_plot(
    nodes = nodes_dt,
    edges = edges_dt,
    selected_genes = selected_genes,
    chromosome = chromosome,
    start = start,
    end = end,
    label_features = label_features,
    title = title,
    gwas_sumstats = gwas_sumstats,
    label_top_gwas_snp = label_top_gwas_snp,
    rsid_pmid = rsid_pmid,
    label_top_lit_snps = label_top_lit_snps,
    pmid_query = pmid_query,
    pmid_page_size = pmid_page_size,
    verbose = verbose
  )

  if (isTRUE(verbose) && identical(plot_obj$plot_data$snp_label_mode, "fallback")) {
    message("No dbSNP-backed SNPs remained after optional disease filtering; falling back to top GWAS SNPs inside the top regulatory elements by germline score.")
  }

  if (isTRUE(save_plot)) {
    if (is.null(plot_file_path) || !nzchar(plot_file_path)) {
      stop("`plot_file_path` must be provided when `save_plot = TRUE`.")
    }

    save_locus_context_plot(
      nodes = nodes_dt,
      edges = edges_dt,
      selected_genes = selected_genes,
      chromosome = chromosome,
      start = start,
      end = end,
      file_path = plot_file_path,
      label_features = label_features,
      title = title,
      gwas_sumstats = gwas_sumstats,
      label_top_gwas_snp = label_top_gwas_snp,
      rsid_pmid = rsid_pmid,
      label_top_lit_snps = label_top_lit_snps,
      pmid_query = pmid_query,
      pmid_page_size = pmid_page_size,
      width = width,
      height = height,
      dpi = dpi
    )
  }

  invisible(new_bundle(
    type = "locus_plot",
    objects = list(
      plot = plot_obj$plot,
      plot_data = plot_obj$plot_data
    ),
    output_paths = list(
      plot_file_path = if (isTRUE(save_plot)) plot_file_path else NULL
    ),
    config = list(
      chromosome = chromosome,
      start = start,
      end = end,
      title = title,
      label_top_lit_snps = label_top_lit_snps
    )
  ))
}
