#' @keywords internal
.conseguiR_runtime_requirements <- function() {
  list(
    core_r_packages = c(
      "dndscv",
      "fishHook",
      "rtracklayer",
      "BSgenome.Hsapiens.UCSC.hg38",
      "data.table",
      "igraph",
      "ggplot2",
      "ggrepel",
      "basilisk",
      "reticulate"
    ),
    optional_r_packages = character(),
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
.conseguiR_resolve_magma_path <- function(magma_path = NULL, must_work = TRUE) {
  explicit_path <- NULL
  if (!is.null(magma_path)) {
    if (!is.character(magma_path) || length(magma_path) != 1L || !nzchar(magma_path)) {
      stop("`magma_path` must be NULL or a single non-empty character string.")
    }
    explicit_path <- magma_path
  }

  option_path <- getOption("conseguiR.magma_path", NULL)
  if (!is.null(option_path) && (!is.character(option_path) || length(option_path) != 1L || !nzchar(option_path))) {
    stop("`options(conseguiR.magma_path = ...)` must contain a single non-empty character string.")
  }

  env_path <- Sys.getenv("CONSEGUIR_MAGMA_PATH", unset = "")
  if (!nzchar(env_path)) {
    env_path <- NULL
  }

  path_path <- Sys.which("magma")
  if (!nzchar(path_path)) {
    path_path <- NULL
  }

  dev_relpath <- "tools/magma_v1/magma"
  dev_path <- NULL
  pkg_root <- if (exists(".conseguiR_pkg_root", inherits = FALSE)) .conseguiR_pkg_root else NULL
  dev_candidates <- c(
    if (!is.null(pkg_root)) file.path(pkg_root, dev_relpath),
    file.path(getwd(), dev_relpath)
  )
  dev_candidates <- unique(dev_candidates[nzchar(dev_candidates)])
  existing_dev <- dev_candidates[file.exists(dev_candidates)]
  if (length(existing_dev) > 0L) {
    dev_path <- existing_dev[[1L]]
  }

  candidates <- unique(Filter(Negate(is.null), c(explicit_path, option_path, env_path, path_path, dev_path)))

  if (length(candidates) == 0L) {
    if (isTRUE(must_work)) {
      stop("Could not locate a MAGMA executable. ", .conseguiR_magma_resolution_note(), call. = FALSE)
    }
    return(NULL)
  }

  for (candidate in candidates) {
    candidate_exists <- file.exists(candidate)
    if (!candidate_exists) {
      next
    }
    if (isTRUE(file.info(candidate)$isdir)) {
      next
    }
    if (file.access(candidate, mode = 1) == 0) {
      return(normalizePath(candidate, winslash = "/", mustWork = TRUE))
    }
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
#' @return A named list describing the current runtime status.
#' @export
check_conseguiR_runtime <- function(quiet = FALSE) {
  req <- .conseguiR_runtime_requirements()

  core_r_packages <- setNames(
    vapply(req$core_r_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
    req$core_r_packages
  )

  optional_r_packages <- setNames(
    vapply(req$optional_r_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
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
  magma_ok <- !is.null(magma_path) && file.exists(magma_path) && file.access(magma_path, mode = 1) == 0
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
    message(
      "conseguiR runtime check: ",
      if (isTRUE(status$core_ok)) "core OK" else "core incomplete",
      if (!isTRUE(status$python_stage_ok)) " (Python-backed stages unavailable)" else " (full pipeline OK)"
    )
  }

  status
}
