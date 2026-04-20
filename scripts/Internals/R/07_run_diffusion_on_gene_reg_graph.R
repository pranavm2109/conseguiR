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

default_diffusion_config <- list(
  nodes_path = "data/processed/gene_reg_graph_scored_nodes.tsv.gz",
  edges_path = "data/processed/gene_reg_graph_scored_edges.tsv.gz",
  output_dir = "data/processed",
  output_stem = "gene_reg_graph_diffusion",
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
  top_n_to_save = 50L
)

## Python execution is now delegated to the shared basilisk helper sourced by
## the user-facing external API runtime.

read_scored_gene_reg_nodes <- function(path = default_diffusion_config$nodes_path) {
  if (!file.exists(path)) {
    stop("Scored gene-reg node file does not exist: ", path)
  }

  as.data.table(fread(path, sep = "\t", showProgress = FALSE))
}

read_scored_gene_reg_edges <- function(path = default_diffusion_config$edges_path) {
  if (!file.exists(path)) {
    stop("Scored gene-reg edge file does not exist: ", path)
  }

  as.data.table(fread(path, sep = "\t", showProgress = FALSE))
}

validate_scored_gene_reg_nodes <- function(nodes) {
  dt <- as.data.table(nodes)
  required_cols <- c("node_id", "node_type", "somatic_score", "germline_score", "epigenomic_score")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Scored gene-reg node table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  bad_types <- setdiff(unique(dt$node_type), c("gene", "reg"))
  if (length(bad_types) > 0L) {
    stop("Unsupported node_type values found in scored gene-reg nodes: ", paste(bad_types, collapse = ", "))
  }

  if (any(is.na(dt$node_id)) || any(duplicated(dt$node_id))) {
    stop("Scored gene-reg node table must contain unique, non-missing node_id values.")
  }

  for (col in c("somatic_score", "germline_score", "epigenomic_score")) {
    if (any(is.na(dt[[col]]))) {
      stop("Score column ", col, " contains missing values.")
    }
  }

  dt
}

validate_scored_gene_reg_edges <- function(edges) {
  dt <- as.data.table(edges)
  required_cols <- c("from", "to", "confidence")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Scored gene-reg edge table is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (any(is.na(dt$from)) || any(is.na(dt$to))) {
    stop("Scored gene-reg edge table contains missing endpoint identifiers.")
  }

  if (any(is.na(dt$confidence))) {
    stop("Scored gene-reg edge table contains missing confidence values.")
  }

  dt
}

run_gene_reg_diffusion <- function(
  nodes_path = default_diffusion_config$nodes_path,
  edges_path = default_diffusion_config$edges_path,
  output_dir = default_diffusion_config$output_dir,
  output_stem = default_diffusion_config$output_stem,
  top_k = default_diffusion_config$top_k,
  confidence_power = default_diffusion_config$confidence_power,
  beta_germline = default_diffusion_config$beta_germline,
  beta_somatic = default_diffusion_config$beta_somatic,
  beta_epigenomic = default_diffusion_config$beta_epigenomic,
  integration_weight_germline = default_diffusion_config$integration_weight_germline,
  integration_weight_somatic = default_diffusion_config$integration_weight_somatic,
  integration_weight_epigenomic = default_diffusion_config$integration_weight_epigenomic,
  positive_only = default_diffusion_config$positive_only,
  reg_signal_clip = default_diffusion_config$reg_signal_clip,
  top_n_to_save = default_diffusion_config$top_n_to_save,
  python_path = NULL
) {
  if (!is.null(python_path) && nzchar(python_path)) {
    warning(
      "`python_path` is deprecated and ignored. ",
      "conseguiR now runs diffusion inside a managed basilisk environment.",
      call. = FALSE
    )
  }

  config <- list(
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
    top_n_to_save = top_n_to_save
  )

  output_paths <- .conseguiR_basilisk_run_python_module(
    script_path = conseguiR_runtime_file("scripts/Internals/Python/07_run_diffusion_on_gene_reg_graph.py"),
    module_name = "conseguiR_step7_r",
    function_name = "run_gene_reg_diffusion",
    config = config,
    config_class_name = "DiffusionConfig"
  )

  list(
    config = config,
    output_paths = output_paths
  )
}
