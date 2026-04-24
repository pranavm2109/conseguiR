# Package startup support for optional Python-backed stages.
# This file is intentionally not exported and should avoid expensive side effects
# at load time so package startup remains lightweight and predictable.

.conseguiR_state <- new.env(parent = emptyenv())
.conseguiR_state$pkg_root <- NULL
.conseguiR_state$basilisk_status <- NULL
.conseguiR_state$magma_path <- NULL
.conseguiR_default_conda_env <- function() {
  env <- getOption("conseguiR.conda_env", Sys.getenv("CONSEGUIR_CONDA_ENV", unset = ""))
  env <- trimws(as.character(env %||% ""))
  if (!nzchar(env)) {
    return(NULL)
  }

  env
}

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

consequIR_valid_python_bin <- function(path) {
  is.character(path) && length(path) == 1L && nzchar(path) && file.exists(path) && !dir.exists(path)
}

consequIR_python_from_conda_run <- function(env = .conseguiR_default_conda_env()) {
  if (is.null(env) || !nzchar(env)) {
    return(NULL)
  }

  if (!nzchar(Sys.which("conda"))) {
    return(NULL)
  }

  out <- tryCatch(
    system2(
      "conda",
      args = c("run", "-n", env, "python", "-c", shQuote("import sys; print(sys.executable)")),
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) NULL
  )

  if (is.null(out) || length(out) == 0L) {
    return(NULL)
  }

  out <- trimws(out)
  out <- out[nzchar(out)]
  out <- out[file.exists(out) & !dir.exists(out)]

  if (length(out) == 0L) {
    return(NULL)
  }

  python_path <- out[[1]]
  if (consequIR_valid_python_bin(python_path)) {
    python_path
  } else {
    NULL
  }
}

consequIR_python_from_reticulate <- function(env = .conseguiR_default_conda_env()) {
  if (is.null(env) || !nzchar(env)) {
    return(NULL)
  }

  if (!requireNamespace("reticulate", quietly = TRUE)) {
    return(NULL)
  }

  envs <- tryCatch(reticulate::conda_list(), error = function(e) NULL)
  if (is.null(envs) || !("name" %in% names(envs))) {
    return(NULL)
  }

  row <- envs[envs$name == env, , drop = FALSE]
  if (nrow(row) != 1L || is.na(row$python)) {
    return(NULL)
  }

  python_path <- as.character(row$python)
  if (consequIR_valid_python_bin(python_path)) {
    python_path
  } else {
    NULL
  }
}

consequIR_python_from_conda_prefix <- function() {
  prefix <- Sys.getenv("CONDA_PREFIX")
  if (!nzchar(prefix) || !dir.exists(prefix)) {
    return(NULL)
  }

  python_path <- file.path(prefix, "bin", "python")
  if (consequIR_valid_python_bin(python_path)) {
    python_path
  } else {
    NULL
  }
}

consequIR_find_python <- function() {
  configured_env <- .conseguiR_default_conda_env()
  current_env <- Sys.getenv("CONDA_DEFAULT_ENV")
  if ((!is.null(configured_env) && identical(current_env, configured_env)) || nzchar(Sys.getenv("CONDA_PREFIX"))) {
    python_path <- consequIR_python_from_conda_prefix()
    if (!is.null(python_path)) {
      return(python_path)
    }

    python_path <- Sys.which("python3")
    if (!nzchar(python_path)) {
      python_path <- Sys.which("python")
    }
    if (!is.null(python_path) && consequIR_valid_python_bin(python_path)) {
      return(python_path)
    }
  }

  python_path <- consequIR_python_from_conda_run(configured_env)
  if (!is.null(python_path)) {
    return(python_path)
  }

  python_path <- consequIR_python_from_reticulate(configured_env)
  if (!is.null(python_path)) {
    return(python_path)
  }

  NULL
}

.conseguiR_magma_autodiscovery_candidates <- function() {
  roots <- unique(Filter(
    function(x) is.character(x) && length(x) == 1L && nzchar(x) && dir.exists(x),
    c(
      .conseguiR_state$pkg_root,
      getwd(),
      Sys.getenv("HOME", unset = "")
    )
  ))

  patterns <- c(
    "magma",
    "bin/magma",
    "tools/magma",
    "tools/magma/bin/magma",
    "tools/magma_v1/magma",
    "tools/magma_v1.1/magma",
    "tools/magma_v1.10/magma",
    "tools/magma_linux/magma",
    "tools/magma_v1_linux/magma",
    "tools/magma_v1.1_linux/magma",
    "tools/magma_v1.10_linux/magma"
  )

  candidates <- unlist(
    lapply(
      roots,
      function(root) {
        unique(c(
          file.path(root, patterns),
          Sys.glob(file.path(root, "tools", "magma*", "magma"))
        ))
      }
    ),
    use.names = FALSE
  )

  candidates <- unique(candidates[file.exists(candidates) & !dir.exists(candidates)])
  as.character(candidates)
}

.conseguiR_autoconfigure_magma_path <- function() {
  existing <- getOption("conseguiR.magma_path", NULL)
  if (is.character(existing) && length(existing) == 1L && nzchar(existing)) {
    .conseguiR_state$magma_path <- existing
    return(invisible(existing))
  }

  env_path <- Sys.getenv("CONSEGUIR_MAGMA_PATH", unset = "")
  if (nzchar(env_path) && file.exists(env_path) && !dir.exists(env_path) && file.access(env_path, mode = 1) == 0) {
    options(conseguiR.magma_path = normalizePath(env_path, winslash = "/", mustWork = TRUE))
    .conseguiR_state$magma_path <- getOption("conseguiR.magma_path")
    return(invisible(.conseguiR_state$magma_path))
  }

  path_path <- Sys.which("magma")
  if (nzchar(path_path) && file.exists(path_path) && !dir.exists(path_path) && file.access(path_path, mode = 1) == 0) {
    options(conseguiR.magma_path = normalizePath(path_path, winslash = "/", mustWork = TRUE))
    .conseguiR_state$magma_path <- getOption("conseguiR.magma_path")
    return(invisible(.conseguiR_state$magma_path))
  }

  candidates <- .conseguiR_magma_autodiscovery_candidates()
  usable <- candidates[file.access(candidates, mode = 1) == 0]
  if (length(usable) > 0L) {
    options(conseguiR.magma_path = normalizePath(usable[[1]], winslash = "/", mustWork = TRUE))
    .conseguiR_state$magma_path <- getOption("conseguiR.magma_path")
    return(invisible(.conseguiR_state$magma_path))
  }

  invisible(NULL)
}

.onLoad <- function(libname, pkgname) {
  if (!is.null(libname) && !is.null(pkgname)) {
    .conseguiR_state$pkg_root <- file.path(libname, pkgname)
  }

  .conseguiR_autoconfigure_magma_path()

  invisible()
}

.onAttach <- function(libname, pkgname) {
  warm_python <- getOption(
    "conseguiR.warm_python_on_attach",
    Sys.getenv("CONSEGUIR_WARM_PYTHON_ON_ATTACH", unset = "true")
  )
  warm_python <- tolower(trimws(as.character(warm_python %||% "true")))
  if (!(warm_python %in% c("true", "1", "yes", "y"))) {
    return(invisible())
  }

  warmup <- tryCatch(
    .conseguiR_basilisk_warmup(),
    error = function(e) e
  )

  if (inherits(warmup, "error")) {
    .conseguiR_state$basilisk_status <- list(
      python_path = NULL,
      python_ok = FALSE,
      python_modules = stats::setNames(logical(), character()),
      error = conditionMessage(warmup)
    )
    return(invisible())
  }

  invisible()
}
