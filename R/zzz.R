# Package startup support for optional Python-backed stages.
# This file is intentionally not exported and should avoid expensive side effects
# at load time so package startup remains lightweight and predictable.

.conseguiR_state <- new.env(parent = emptyenv())
.conseguiR_state$pkg_root <- NULL
.conseguiR_state$basilisk_status <- NULL
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
  if (!is.null(configured_env) && identical(current_env, configured_env)) {
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

.onLoad <- function(libname, pkgname) {
  if (!is.null(libname) && !is.null(pkgname)) {
    .conseguiR_state$pkg_root <- file.path(libname, pkgname)
  }

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
