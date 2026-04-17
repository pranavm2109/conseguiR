.conseguiR_backend_cache <- new.env(parent = emptyenv())

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
    if (!is.null(.conseguiR_state$pkg_root)) file.path(.conseguiR_state$pkg_root, "data", "processed"),
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
  cache_key <- paste0("backend_paths::", normalizePath(backend_dir, winslash = "/", mustWork = FALSE))
  if (exists(cache_key, envir = .conseguiR_backend_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .conseguiR_backend_cache, inherits = FALSE))
  }

  paths <- list(
    backend_dir = backend_dir,
    gene_reg_graph_rds = file.path(backend_dir, "gene_reg_graph_no_scores.rds"),
    gene_reg_graph_nodes = file.path(backend_dir, "gene_reg_graph_no_scores_nodes.tsv.gz"),
    gene_reg_graph_edges = file.path(backend_dir, "gene_reg_graph_no_scores_edges.tsv.gz"),
    gene_gene_graph_rds = file.path(backend_dir, "gene_gene_graph.rds"),
    gene_gene_graph_nodes = file.path(backend_dir, "gene_gene_graph_nodes.tsv.gz"),
    gene_gene_graph_edges = file.path(backend_dir, "gene_gene_graph_edges.tsv.gz")
  )
  assign(cache_key, paths, envir = .conseguiR_backend_cache)
  paths
}

#' @keywords internal
.conseguiR_existing_graph_paths <- function(paths) {
  flat <- unlist(paths, use.names = FALSE)
  flat[file.exists(flat)]
}

#' @keywords internal
.conseguiR_backend_resource_cache_dir <- function(create = FALSE) {
  path <- file.path(
    tools::R_user_dir("conseguiR", which = "cache"),
    "backend_resources"
  )
  if (isTRUE(create)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
  }
  normalizePath(path, winslash = "/", mustWork = FALSE)
}

#' @keywords internal
.conseguiR_backend_resource_dirs <- function(include_cache = TRUE) {
  option_dir <- getOption("conseguiR.backend_resource_dir")
  if (!is.null(option_dir)) {
    option_dir <- trimws(as.character(option_dir))
    option_dir <- option_dir[nzchar(option_dir)]
  }

  installed <- tryCatch(
    system.file("extdata", "backend", package = "conseguiR"),
    error = function(e) ""
  )
  candidates <- c(
    if (isTRUE(include_cache)) .conseguiR_backend_resource_cache_dir(),
    option_dir,
    installed,
    if (!is.null(.conseguiR_state$pkg_root)) {
      file.path(.conseguiR_state$pkg_root, "extdata", "backend")
    },
    if (!is.null(.conseguiR_state$pkg_root)) {
      file.path(.conseguiR_state$pkg_root, "inst", "extdata", "backend")
    },
    file.path(getwd(), "extdata", "backend"),
    file.path(getwd(), "inst", "extdata", "backend")
  )
  candidates <- unique(candidates[!is.na(candidates) & nzchar(candidates)])
  candidates[dir.exists(candidates)]
}

#' @keywords internal
.conseguiR_backend_seed_dir <- function() {
  cache_key <- "backend_seed_dir"
  if (exists(cache_key, envir = .conseguiR_backend_cache, inherits = FALSE)) {
    return(get(cache_key, envir = .conseguiR_backend_cache, inherits = FALSE))
  }

  existing <- .conseguiR_backend_resource_dirs(include_cache = TRUE)
  seed_dir <- if (length(existing) > 0L) {
    normalizePath(existing[[1]], winslash = "/", mustWork = TRUE)
  } else {
    NULL
  }
  assign(cache_key, seed_dir, envir = .conseguiR_backend_cache)
  seed_dir
}

#' @keywords internal
.conseguiR_backend_resource_candidates <- function(filename) {
  dirs <- .conseguiR_backend_resource_dirs(include_cache = TRUE)
  if (length(dirs) == 0L) {
    return(character())
  }

  candidates <- file.path(dirs, filename)
  unique(candidates[file.exists(candidates)])
}

#' @keywords internal
.conseguiR_backend_resource_path <- function(filename) {
  candidates <- .conseguiR_backend_resource_candidates(filename)
  if (length(candidates) == 0L) {
    return(NULL)
  }

  normalizePath(candidates[[1]], winslash = "/", mustWork = TRUE)
}

#' @keywords internal
.conseguiR_prefer_repo_relative_path <- function(relpath) {
  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(relpath)
  }

  NULL
}

#' @keywords internal
.conseguiR_default_gene_loc_path <- function() {
  path <- .conseguiR_prefer_repo_relative_path("data/raw/NCBI38/NCBI38.gene.loc")
  if (!is.null(path)) {
    return(path)
  }
  path <- .conseguiR_find_data_path("data/raw/NCBI38/NCBI38.gene.loc")
  if (!is.null(path)) {
    return(path)
  }

  .conseguiR_backend_resource_path("NCBI38.gene.loc")
}

#' @keywords internal
.conseguiR_default_reg_loc_path <- function() {
  path <- .conseguiR_prefer_repo_relative_path("data/raw/GeneHancer/2026-01-26_UCSC_all_unfiltered_reg_elements.loc")
  if (!is.null(path)) {
    return(path)
  }
  path <- .conseguiR_find_data_path("data/raw/GeneHancer/2026-01-26_UCSC_all_unfiltered_reg_elements.loc")
  if (!is.null(path)) {
    return(path)
  }

  .conseguiR_backend_resource_path("2026-01-26_UCSC_all_unfiltered_reg_elements.loc")
}

#' @keywords internal
.conseguiR_seed_backend_graph <- function(kind = c("gene_reg", "gene_gene"), backend_dir) {
  kind <- match.arg(kind)
  seed_dir <- .conseguiR_backend_seed_dir()
  if (is.null(seed_dir)) {
    return(FALSE)
  }

  files <- switch(kind,
    gene_reg = c(
      "gene_reg_graph_no_scores_nodes.tsv.gz",
      "gene_reg_graph_no_scores_edges.tsv.gz",
      "gene_reg_graph_no_scores.rds"
    ),
    gene_gene = c(
      "gene_gene_graph_nodes.tsv.gz",
      "gene_gene_graph_edges.tsv.gz",
      "gene_gene_graph.rds"
    )
  )

  src <- file.path(seed_dir, files)
  existing_src <- src[file.exists(src)]
  if (length(existing_src) < 2L) {
    return(FALSE)
  }

  dest <- file.path(backend_dir, basename(existing_src))
  dir.create(backend_dir, recursive = TRUE, showWarnings = FALSE)
  ok <- file.copy(existing_src, dest, overwrite = TRUE)
  all(ok)
}

#' @keywords internal
.conseguiR_data_candidates <- function(relpath) {
  unique(c(
    if (!is.null(.conseguiR_state$pkg_root)) file.path(.conseguiR_state$pkg_root, relpath),
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

#' @keywords internal
.conseguiR_backend_init_cache_key <- function(
  backend_dir,
  build_gene_reg,
  build_gene_gene,
  force,
  strict
) {
  paste(
    "backend_init",
    normalizePath(backend_dir, winslash = "/", mustWork = FALSE),
    build_gene_reg,
    build_gene_gene,
    force,
    strict,
    sep = "::"
  )
}

#' @keywords internal
.conseguiR_cached_backend_init <- function(
  init_cache_key,
  paths,
  build_gene_reg,
  build_gene_gene,
  quiet
) {
  if (!exists(init_cache_key, envir = .conseguiR_backend_cache, inherits = FALSE)) {
    return(NULL)
  }

  cached <- get(init_cache_key, envir = .conseguiR_backend_cache, inherits = FALSE)
  required_paths <- c(
    if (isTRUE(build_gene_reg)) {
      .conseguiR_existing_graph_paths(
        paths[c("gene_reg_graph_nodes", "gene_reg_graph_edges")]
      )
    },
    if (isTRUE(build_gene_gene)) {
      .conseguiR_existing_graph_paths(
        paths[c("gene_gene_graph_nodes", "gene_gene_graph_edges")]
      )
    }
  )

  if (
    is.null(cached) ||
      length(required_paths) == 0L ||
      !.conseguiR_graph_files_exist(required_paths)
  ) {
    return(NULL)
  }

  if (!isTRUE(quiet)) {
    message("Backend graph initialization complete.")
  }

  cached
}

#' @keywords internal
.conseguiR_initialize_gene_reg_backend <- function(
  paths,
  backend_dir,
  force,
  strict
) {
  target_paths <- list(
    paths$gene_reg_graph_nodes,
    paths$gene_reg_graph_edges
  )

  if (!isTRUE(force) && .conseguiR_graph_files_exist(target_paths)) {
    return("reused")
  }

  if (.conseguiR_seed_backend_graph("gene_reg", backend_dir = backend_dir)) {
    return("seeded")
  }

  interactions_path <- .conseguiR_find_data_path(
    "data/raw/GeneHancer/gh_interactions_hg38_primary_assembly"
  )
  reg_elements_path <- .conseguiR_find_data_path(
    "data/raw/GeneHancer/gh_reg_elements_hg38_primary_assembly"
  )

  if (is.null(interactions_path) || is.null(reg_elements_path)) {
    msg <- paste(
      "Unable to initialize the backend gene-regulatory graph because the",
      "required GeneHancer resources were not found."
    )
    if (isTRUE(strict)) {
      stop(msg)
    }
    return("skipped")
  }

  env <- .conseguiR_internal_script_env(
    "scripts/Internals/R/01_prepare_gene_reg_graph.R"
  )
  config <- env$default_config
  config$interactions_path <- interactions_path
  config$reg_elements_path <- reg_elements_path
  config$output_prefix <- file.path(backend_dir, "gene_reg_graph_no_scores")
  env$prepare_gene_reg_graph(config = config)
  "built"
}

#' @keywords internal
.conseguiR_initialize_gene_gene_backend <- function(
  paths,
  backend_dir,
  force,
  strict
) {
  target_paths <- list(
    paths$gene_gene_graph_nodes,
    paths$gene_gene_graph_edges
  )

  if (!isTRUE(force) && .conseguiR_graph_files_exist(target_paths)) {
    return("reused")
  }

  if (.conseguiR_seed_backend_graph("gene_gene", backend_dir = backend_dir)) {
    return("seeded")
  }

  links_path <- .conseguiR_find_data_path(
    "data/raw/STRING/9606.protein.links.v12.0.txt"
  )
  info_path <- .conseguiR_find_data_path(
    "data/raw/STRING/9606.protein.info.v12.0.txt"
  )

  if (is.null(links_path) || is.null(info_path)) {
    msg <- paste(
      "Unable to initialize the backend gene-gene graph because the",
      "required STRING resources were not found."
    )
    if (isTRUE(strict)) {
      stop(msg)
    }
    return("skipped")
  }

  env <- .conseguiR_internal_script_env(
    "scripts/Internals/R/02_prepare_gene_gene_graph.R"
  )
  config <- env$default_config
  config$protein_links_path <- links_path
  config$protein_info_path <- info_path
  config$output_prefix <- file.path(backend_dir, "gene_gene_graph")
  env$prepare_gene_gene_graph(config = config)
  "built"
}

#' @keywords internal
.conseguiR_backend_init_result <- function(
  backend_dir,
  paths,
  status,
  init_cache_key,
  force
) {
  result <- structure(
    list(
      backend_dir = backend_dir,
      output_paths = paths,
      status = status
    ),
    class = c("conseguiR_backend_init", "list")
  )

  if (!isTRUE(force)) {
    assign(init_cache_key, result, envir = .conseguiR_backend_cache)
  }

  result
}

#' @keywords internal
.conseguiR_backend_init_context <- function(
  backend_dir,
  build_gene_reg,
  build_gene_gene,
  force,
  strict
) {
  if (!is.null(backend_dir)) {
    options(conseguiR.backend_dir = backend_dir)
  }

  resolved_backend_dir <- .conseguiR_backend_dir(create = TRUE)
  options(conseguiR.backend_dir = resolved_backend_dir)
  paths <- .conseguiR_backend_paths(resolved_backend_dir)
  init_cache_key <- .conseguiR_backend_init_cache_key(
    backend_dir = resolved_backend_dir,
    build_gene_reg = build_gene_reg,
    build_gene_gene = build_gene_gene,
    force = force,
    strict = strict
  )

  list(
    backend_dir = resolved_backend_dir,
    paths = paths,
    init_cache_key = init_cache_key
  )
}

#' @keywords internal
.conseguiR_backend_init_status <- function(
  paths,
  backend_dir,
  build_gene_reg,
  build_gene_gene,
  force,
  strict
) {
  status <- list()

  if (isTRUE(build_gene_reg)) {
    status$gene_reg <- .conseguiR_initialize_gene_reg_backend(
      paths = paths,
      backend_dir = backend_dir,
      force = force,
      strict = strict
    )
  }

  if (isTRUE(build_gene_gene)) {
    status$gene_gene <- .conseguiR_initialize_gene_gene_backend(
      paths = paths,
      backend_dir = backend_dir,
      force = force,
      strict = strict
    )
  }

  status
}

.conseguiR_initialize_backend_graphs <- function(
  backend_dir = NULL,
  build_gene_reg = TRUE,
  build_gene_gene = TRUE,
  force = FALSE,
  strict = TRUE,
  quiet = FALSE
) {
  context <- .conseguiR_backend_init_context(
    backend_dir = backend_dir,
    build_gene_reg = build_gene_reg,
    build_gene_gene = build_gene_gene,
    force = force,
    strict = strict
  )
  paths <- context$paths
  backend_dir <- context$backend_dir
  init_cache_key <- context$init_cache_key

  if (!isTRUE(force)) {
    cached <- .conseguiR_cached_backend_init(
      init_cache_key = init_cache_key,
      paths = paths,
      build_gene_reg = build_gene_reg,
      build_gene_gene = build_gene_gene,
      quiet = quiet
    )
    if (!is.null(cached)) {
      return(cached)
    }
  }

  status <- .conseguiR_backend_init_status(
    paths = paths,
    backend_dir = backend_dir,
    build_gene_reg = build_gene_reg,
    build_gene_gene = build_gene_gene,
    force = force,
    strict = strict
  )

  if (!isTRUE(quiet)) {
    message("Backend graph initialization complete.")
  }

  .conseguiR_backend_init_result(
    backend_dir = backend_dir,
    paths = paths,
    status = status,
    init_cache_key = init_cache_key,
    force = force
  )
}
