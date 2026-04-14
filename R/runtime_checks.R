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
    python_modules = c("ortools", "pandas", "numpy"),
    magma_relpath = "tools/magma_v1/magma"
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

  magma_path <- .conseguiR_find_data_path(req$magma_relpath)
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
