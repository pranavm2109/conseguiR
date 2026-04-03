#' @keywords internal
.conseguiR_runtime_requirements <- function() {
  list(
    r_packages = c(
      "dndscv",
      "fishHook",
      "rtracklayer",
      "BSgenome.Hsapiens.UCSC.hg38",
      "data.table",
      "igraph",
      "ggplot2",
      "ggrepel",
      "reticulate"
    ),
    python_modules = c("ortools", "pandas", "numpy"),
    magma_relpath = "tools/magma_v1/magma"
  )
}

#' Check the conseguiR runtime environment
#'
#' Checks whether the expected R packages, Python interpreter, Python modules,
#' and MAGMA executable are available for the current `conseguiR` runtime.
#'
#' @param quiet Whether to suppress the summary message.
#'
#' @return A named list describing the current runtime status.
#' @export
check_conseguiR_runtime <- function(quiet = FALSE) {
  req <- .conseguiR_runtime_requirements()

  r_packages <- setNames(
    vapply(req$r_packages, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1)),
    req$r_packages
  )

  python_path <- getOption("conseguiR.python")
  if (is.character(python_path)) {
    python_path <- trimws(python_path)
    python_path <- python_path[nzchar(python_path)]
    python_path <- python_path[file.exists(python_path) & !dir.exists(python_path)]
    if (length(python_path) > 0L) {
      python_path <- python_path[[1]]
    } else {
      python_path <- NULL
    }
  } else {
    python_path <- NULL
  }

  python_ok <- is.character(python_path) && length(python_path) == 1L &&
    nzchar(python_path) && file.exists(python_path)

  python_modules <- setNames(rep(FALSE, length(req$python_modules)), req$python_modules)

  if (isTRUE(python_ok)) {
    py_mods <- paste(sprintf("'%s'", req$python_modules), collapse = ", ")
    module_expr <- sprintf(
      "import importlib.util as u, json; mods=[%s]; print(json.dumps({m:(u.find_spec(m) is not None) for m in mods}))",
      py_mods
    )

    out <- tryCatch(
      system2(
        command = python_path,
        args = c("-c", shQuote(module_expr)),
        stdout = TRUE,
        stderr = TRUE
      ),
      error = function(e) NULL
    )

    if (!is.null(out) && length(out) > 0L) {
      line <- tail(out, 1)
      parsed <- tryCatch(jsonlite::fromJSON(line), error = function(e) NULL)
      if (!is.null(parsed)) {
        python_modules[names(parsed)] <- as.logical(parsed)
      }
    }
  }

  magma_path <- .conseguiR_find_data_path(req$magma_relpath)
  magma_ok <- !is.null(magma_path) && file.exists(magma_path) && file.access(magma_path, mode = 1) == 0

  status <- list(
    conda_env = getOption("conseguiR.conda_env"),
    python_path = python_path,
    python_ok = python_ok,
    r_packages = r_packages,
    python_modules = python_modules,
    magma_path = magma_path,
    magma_ok = magma_ok,
    ok = isTRUE(python_ok) && isTRUE(magma_ok) && all(r_packages) && all(python_modules)
  )

  options(conseguiR.runtime_status = status)

  if (!isTRUE(quiet)) {
    message(
      "conseguiR runtime check: ",
      if (isTRUE(status$ok)) "OK" else "incomplete"
    )
  }

  status
}
