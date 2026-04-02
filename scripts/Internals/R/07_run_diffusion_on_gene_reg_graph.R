#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

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
  positive_only = FALSE,
  reg_signal_clip = 5.0,
  top_n_to_save = 50L
)

find_python <- function() {
  configured <- getOption("conseguiR.python", NULL)
  if (!is.null(configured) && nzchar(configured) && file.exists(configured)) {
    return(configured)
  }

  conda <- Sys.which("conda")
  if (nzchar(conda)) {
    out <- tryCatch(
      system2(
        conda,
        args = c("run", "-n", "lymphoma_graph_env", "python", "-c", shQuote("import sys; print(sys.executable)")),
        stdout = TRUE,
        stderr = FALSE
      ),
      error = function(e) character()
    )
    out <- trimws(out)
    if (length(out) == 1L && nzchar(out) && file.exists(out)) {
      return(out)
    }
  }

  if (identical(Sys.getenv("CONDA_DEFAULT_ENV"), "lymphoma_graph_env")) {
    prefix <- Sys.getenv("CONDA_PREFIX")
    if (nzchar(prefix)) {
      python_path <- file.path(prefix, "bin", "python")
      if (file.exists(python_path)) {
        return(python_path)
      }
    }
  }

  python <- Sys.which("python")
  if (nzchar(python)) {
    return(python)
  }

  python3 <- Sys.which("python3")
  if (nzchar(python3)) {
    return(python3)
  }

  stop("Could not find a Python interpreter on PATH. Install Python 3 or make it available as 'python3' or 'python'.")
}

python_module_execute <- function(script_path, module_name, function_name, config, python_path = NULL) {
  if (!file.exists(script_path)) {
    stop("Python script not found: ", script_path)
  }

  if (is.null(python_path) || python_path == "") {
    python_path <- find_python()
  }

  if (!file.exists(python_path)) {
    stop("Python interpreter not found at: ", python_path)
  }

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("The 'jsonlite' package is required for Python-backed wrappers. Install it with install.packages('jsonlite').")
  }

  payload <- list(
    script_path = normalizePath(script_path, mustWork = TRUE),
    module_name = module_name,
    function_name = function_name,
    config = config
  )

  json_text <- jsonlite::toJSON(payload, auto_unbox = TRUE)
  python_code <- c(
    "import importlib.util, json, sys",
    "from pathlib import Path",
    "payload = json.loads(sys.stdin.read())",
    "script_path = Path(payload['script_path'])",
    "spec = importlib.util.spec_from_file_location(payload['module_name'], str(script_path))",
    "module = importlib.util.module_from_spec(spec)",
    "assert spec.loader is not None",
    "sys.modules[spec.name] = module",
    "spec.loader.exec_module(module)",
    "config = getattr(module, 'DiffusionConfig')(**payload['config'])",
    "result = getattr(module, payload['function_name'])(config=config)",
    "sys.stdout.write(json.dumps(result['output_paths']))"
  )

  temp_script <- tempfile(fileext = ".py")
  stdout_path <- tempfile(fileext = ".json")
  stderr_path <- tempfile(fileext = ".log")
  writeLines(python_code, temp_script)
  on.exit(unlink(c(temp_script, stdout_path, stderr_path), force = TRUE), add = TRUE)

  status <- system2(
    command = python_path,
    args = c(temp_script),
    input = json_text,
    stdout = stdout_path,
    stderr = stderr_path
  )

  stderr_lines <- if (file.exists(stderr_path)) readLines(stderr_path, warn = FALSE) else character()
  stdout_lines <- if (file.exists(stdout_path)) readLines(stdout_path, warn = FALSE) else character()
  if (!is.null(status) && status != 0L) {
    stop(
      "Python wrapper failed with exit status ", status, ".\n",
      paste(c(stderr_lines, stdout_lines), collapse = "\n")
    )
  }

  output_text <- paste(stdout_lines, collapse = "\n")
  if (!nzchar(trimws(output_text))) {
    stop(
      "Python wrapper produced no JSON output.",
      if (length(stderr_lines) > 0L) paste0("\n", paste(stderr_lines, collapse = "\n")) else ""
    )
  }
  jsonlite::fromJSON(output_text)
}

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
  positive_only = default_diffusion_config$positive_only,
  reg_signal_clip = default_diffusion_config$reg_signal_clip,
  top_n_to_save = default_diffusion_config$top_n_to_save,
  python_path = NULL
) {
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
    positive_only = positive_only,
    reg_signal_clip = reg_signal_clip,
    top_n_to_save = top_n_to_save
  )

  output_paths <- python_module_execute(
    script_path = "scripts/Internals/Python/07_run_diffusion_on_gene_reg_graph.py",
    module_name = "conseguiR_step7_r",
    function_name = "run_gene_reg_diffusion",
    config = config,
    python_path = python_path
  )

  list(
    config = config,
    output_paths = output_paths
  )
}
