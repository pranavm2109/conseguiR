#' @keywords internal
.conseguiR_backend_dir <- function(create = FALSE) {
  opt_dir <- getOption("conseguiR.backend_dir")
  if (is.character(opt_dir) && length(opt_dir) == 1L && nzchar(opt_dir)) {
    if (isTRUE(create)) {
      dir.create(opt_dir, recursive = TRUE, showWarnings = FALSE)
    }
    return(normalizePath(opt_dir, winslash = "/", mustWork = FALSE))
  }

  candidates <- c(
    if (!is.null(.conseguiR_pkg_root)) file.path(.conseguiR_pkg_root, "data", "processed"),
    file.path(getwd(), "data", "processed"),
    file.path(tools::R_user_dir("conseguiR", which = "cache"), "backend")
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])

  existing <- candidates[file.exists(candidates)]
  chosen <- if (length(existing) > 0L) existing[[1]] else candidates[[length(candidates)]]

  if (isTRUE(create)) {
    dir.create(chosen, recursive = TRUE, showWarnings = FALSE)
  }

  normalizePath(chosen, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.conseguiR_backend_paths <- function(backend_dir = .conseguiR_backend_dir(create = TRUE)) {
  list(
    backend_dir = backend_dir,
    gene_reg_graph_rds = file.path(backend_dir, "gene_reg_graph_no_scores.rds"),
    gene_reg_graph_nodes = file.path(backend_dir, "gene_reg_graph_no_scores_nodes.tsv.gz"),
    gene_reg_graph_edges = file.path(backend_dir, "gene_reg_graph_no_scores_edges.tsv.gz"),
    gene_gene_graph_rds = file.path(backend_dir, "gene_gene_graph.rds"),
    gene_gene_graph_nodes = file.path(backend_dir, "gene_gene_graph_nodes.tsv.gz"),
    gene_gene_graph_edges = file.path(backend_dir, "gene_gene_graph_edges.tsv.gz")
  )
}

#' @keywords internal
.conseguiR_data_candidates <- function(relpath) {
  unique(c(
    if (!is.null(.conseguiR_pkg_root)) file.path(.conseguiR_pkg_root, relpath),
    file.path(getwd(), relpath)
  ))
}

#' @keywords internal
.conseguiR_find_data_path <- function(relpath) {
  candidates <- .conseguiR_data_candidates(relpath)
  existing <- candidates[file.exists(candidates)]
  if (length(existing) == 0L) {
    return(NULL)
  }
  normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_internal_script_env <- function(relpath) {
  env <- new.env(parent = baseenv())
  sys.source(.conseguiR_runtime_file(relpath), envir = env)
  env
}

#' @keywords internal
.conseguiR_graph_files_exist <- function(paths) {
  all(file.exists(unlist(paths, use.names = FALSE)))
}

.conseguiR_initialize_backend_graphs <- function(
  backend_dir = NULL,
  build_gene_reg = TRUE,
  build_gene_gene = TRUE,
  force = FALSE,
  strict = TRUE,
  quiet = FALSE
) {
  if (!is.null(backend_dir)) {
    options(conseguiR.backend_dir = backend_dir)
  }

  backend_dir <- .conseguiR_backend_dir(create = TRUE)
  options(conseguiR.backend_dir = backend_dir)
  paths <- .conseguiR_backend_paths(backend_dir)

  status <- list()

  if (isTRUE(build_gene_reg)) {
    target_paths <- list(
      paths$gene_reg_graph_rds,
      paths$gene_reg_graph_nodes,
      paths$gene_reg_graph_edges
    )

    if (!isTRUE(force) && .conseguiR_graph_files_exist(target_paths)) {
      status$gene_reg <- "reused"
    } else {
      interactions_path <- .conseguiR_find_data_path("data/raw/GeneHancer/gh_interactions_hg38_primary_assembly")
      reg_elements_path <- .conseguiR_find_data_path("data/raw/GeneHancer/gh_reg_elements_hg38_primary_assembly")

      if (is.null(interactions_path) || is.null(reg_elements_path)) {
        msg <- paste(
          "Unable to initialize the backend gene-regulatory graph because the",
          "required GeneHancer resources were not found."
        )
        if (isTRUE(strict)) stop(msg)
        status$gene_reg <- "skipped"
      } else {
        env <- .conseguiR_internal_script_env("scripts/Internals/R/01_prepare_gene_reg_graph.R")
        config <- env$default_config
        config$interactions_path <- interactions_path
        config$reg_elements_path <- reg_elements_path
        config$output_prefix <- file.path(backend_dir, "gene_reg_graph_no_scores")
        env$prepare_gene_reg_graph(config = config)
        status$gene_reg <- "built"
      }
    }
  }

  if (isTRUE(build_gene_gene)) {
    target_paths <- list(
      paths$gene_gene_graph_rds,
      paths$gene_gene_graph_nodes,
      paths$gene_gene_graph_edges
    )

    if (!isTRUE(force) && .conseguiR_graph_files_exist(target_paths)) {
      status$gene_gene <- "reused"
    } else {
      links_path <- .conseguiR_find_data_path("data/raw/STRING/9606.protein.links.v12.0.txt")
      info_path <- .conseguiR_find_data_path("data/raw/STRING/9606.protein.info.v12.0.txt")

      if (is.null(links_path) || is.null(info_path)) {
        msg <- paste(
          "Unable to initialize the backend gene-gene graph because the",
          "required STRING resources were not found."
        )
        if (isTRUE(strict)) stop(msg)
        status$gene_gene <- "skipped"
      } else {
        env <- .conseguiR_internal_script_env("scripts/Internals/R/02_prepare_gene_gene_graph.R")
        config <- env$default_config
        config$protein_links_path <- links_path
        config$protein_info_path <- info_path
        config$output_prefix <- file.path(backend_dir, "gene_gene_graph")
        env$prepare_gene_gene_graph(config = config)
        status$gene_gene <- "built"
      }
    }
  }

  if (!isTRUE(quiet)) {
    message("Backend graph initialization complete.")
  }

  structure(
    list(
      backend_dir = backend_dir,
      output_paths = paths,
      status = status
    ),
    class = c("conseguiR_backend_init", "list")
  )
}
