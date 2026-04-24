#' @keywords internal
.conseguiR_runtime_requirements <- function() {
  list(
    core_r_packages = c(
      "BSgenome.Hsapiens.UCSC.hg38",
      "basilisk",
      "data.table",
      "digest",
      "dndscv",
      "fishHook",
      "gUtils",
      "ggrepel",
      "ggnewscale",
      "ggplot2",
      "igraph",
      "matrixStats",
      "reticulate",
      "rtracklayer",
      "scales"
    ),
    optional_r_packages = c("RCy3", "tidygraph", "org.Hs.eg.db"),
    python_modules = c("ortools", "pandas", "numpy")
  )
}

#' @keywords internal
.conseguiR_magma_resolution_note <- function() {
  paste(
    "Install MAGMA separately and provide it via `magma_path`,",
    "`options(conseguiR.magma_path = \"/path/to/magma\")`,",
    "the `CONSEGUIR_MAGMA_PATH` environment variable, or by making `magma`",
    "available on your system PATH."
  )
}

#' @keywords internal
.conseguiR_normalize_magma_candidate <- function(path, source_label) {
  if (
    is.null(path) ||
      !is.character(path) ||
      length(path) != 1L ||
      !nzchar(path)
  ) {
    stop(
      source_label,
      " must contain a single non-empty character string."
    )
  }
  path
}

#' @keywords internal
.conseguiR_magma_candidates <- function(magma_path = NULL) {
  explicit_path <- if (!is.null(magma_path)) {
    .conseguiR_normalize_magma_candidate(
      magma_path,
      "`magma_path`"
    )
  } else {
    NULL
  }
  option_path <- getOption("conseguiR.magma_path", NULL)
  if (!is.null(option_path)) {
    option_path <- .conseguiR_normalize_magma_candidate(
      option_path,
      "`options(conseguiR.magma_path = ...)`"
    )
  }
  env_path <- Sys.getenv("CONSEGUIR_MAGMA_PATH", unset = "")
  if (!nzchar(env_path)) {
    env_path <- NULL
  }
  path_path <- Sys.which("magma")
  if (!nzchar(path_path)) {
    path_path <- NULL
  }
  autodiscovered_paths <- tryCatch(
    .conseguiR_magma_autodiscovery_candidates(),
    error = function(e) character()
  )

  unique(Filter(
    Negate(is.null),
    c(explicit_path, option_path, env_path, path_path, autodiscovered_paths)
  ))
}

#' @keywords internal
.conseguiR_find_usable_magma_path <- function(candidates) {
  for (candidate in candidates) {
    if (!file.exists(candidate)) {
      next
    }
    if (isTRUE(file.info(candidate)$isdir)) {
      next
    }
    if (file.access(candidate, mode = 1) == 0) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
  }

  NULL
}

#' @keywords internal
.conseguiR_resolve_magma_path <- function(magma_path = NULL, must_work = TRUE) {
  candidates <- .conseguiR_magma_candidates(magma_path)
  if (length(candidates) == 0L) {
    if (isTRUE(must_work)) {
      stop(
        "Could not locate a MAGMA executable. ",
        .conseguiR_magma_resolution_note(),
        call. = FALSE
      )
    }
    return(NULL)
  }

  resolved <- .conseguiR_find_usable_magma_path(candidates)
  if (!is.null(resolved)) {
    return(resolved)
  }

  if (isTRUE(must_work)) {
    stop(
      "Found MAGMA candidate path(s), but none were usable executables: ",
      paste(shQuote(candidates), collapse = ", "),
      ". ",
      .conseguiR_magma_resolution_note(),
      call. = FALSE
    )
  }

  NULL
}

#' @keywords internal
.conseguiR_named_namespace_check <- function(packages) {
  stats::setNames(
    vapply(
      packages,
      requireNamespace,
      quietly = TRUE,
      FUN.VALUE = logical(1)
    ),
    packages
  )
}

#' @keywords internal
.conseguiR_runtime_status_message <- function(status) {
  paste0(
    "conseguiR runtime check: ",
    if (isTRUE(status$core_ok)) "core OK" else "core incomplete",
    if (!isTRUE(status$python_stage_ok)) {
      " (Python-backed stages unavailable)"
    } else {
      " (full pipeline OK)"
    }
  )
}

#' Check the conseguiR runtime environment
#'
#' Checks whether the expected core R packages and optional external runtimes
#' are available for the current `conseguiR` session.
#'
#' The package baseline is considered usable when the core R dependencies and
#' MAGMA are available. Python is reported separately because it is only needed
#' for the diffusion and selected-subgraph stages.
#'
#' @param quiet Whether to suppress the summary message.
#'
#' @examples
#' status <- check_conseguiR_runtime(quiet = TRUE)
#' is.list(status)
#'
#' @return A named list describing the current runtime status.
#' @export
check_conseguiR_runtime <- function(quiet = FALSE) {
  req <- .conseguiR_runtime_requirements()

  core_r_packages <- .conseguiR_named_namespace_check(req$core_r_packages)
  optional_r_packages <- .conseguiR_named_namespace_check(
    req$optional_r_packages
  )

  py_status <- .conseguiR_basilisk_python_status(req$python_modules)
  python_path <- py_status$python_path
  python_ok <- isTRUE(py_status$python_ok)
  python_modules <- py_status$python_modules

  magma_path <- tryCatch(
    .conseguiR_resolve_magma_path(must_work = FALSE),
    error = function(...) NULL
  )
  magma_ok <- !is.null(magma_path) &&
    file.exists(magma_path) &&
    file.access(magma_path, mode = 1) == 0
  python_stage_ok <- isTRUE(python_ok) && all(python_modules)
  core_ok <- isTRUE(magma_ok) && all(core_r_packages)

  status <- list(
    conda_env = NA_character_,
    python_path = python_path,
    python_ok = python_ok,
    core_r_packages = core_r_packages,
    optional_r_packages = optional_r_packages,
    python_modules = python_modules,
    python_error = py_status$error,
    magma_path = magma_path,
    magma_ok = magma_ok,
    core_ok = core_ok,
    python_stage_ok = python_stage_ok,
    full_pipeline_ok = core_ok && python_stage_ok,
    ok = core_ok
  )

  options(conseguiR.runtime_status = status)

  if (!isTRUE(quiet)) {
    message(.conseguiR_runtime_status_message(status))
  }

  status
}
