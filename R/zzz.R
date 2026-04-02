# Package startup support for the lymphoma_graph_env Python environment.
# This file is intentionally not exported; it only sets package-level options on load.

lymphoma_graph_env_name <- "lymphoma_graph_env"

consequIR_valid_python_bin <- function(path) {
  is.character(path) && length(path) == 1L && nzchar(path) && file.exists(path) && !dir.exists(path)
}

consequIR_python_from_conda_run <- function(env = lymphoma_graph_env_name) {
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

  if (is.null(out) || length(out) != 1L) {
    return(NULL)
  }

  python_path <- trimws(out)
  if (consequIR_valid_python_bin(python_path)) {
    python_path
  } else {
    NULL
  }
}

consequIR_python_from_reticulate <- function(env = lymphoma_graph_env_name) {
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

consequIR_prepend_path <- function(prefix) {
  if (!nzchar(prefix) || !dir.exists(prefix)) {
    return(invisible(NULL))
  }

  env_bin <- file.path(prefix, "bin")
  if (!nzchar(env_bin) || !dir.exists(env_bin)) {
    return(invisible(NULL))
  }

  current_path <- Sys.getenv("PATH")
  path_sep <- .Platform$path.sep
  path_entries <- unlist(strsplit(current_path, path_sep, fixed = TRUE), use.names = FALSE)
  if (env_bin %in% path_entries) {
    return(invisible(NULL))
  }

  Sys.setenv(PATH = paste(c(env_bin, current_path), collapse = path_sep))
  invisible(NULL)
}

consequIR_find_python <- function() {
  current_env <- Sys.getenv("CONDA_DEFAULT_ENV")
  if (identical(current_env, lymphoma_graph_env_name)) {
    python_path <- consequIR_python_from_conda_prefix()
    if (!is.null(python_path)) {
      return(python_path)
    }

    python_path <- Sys.which("python3")
    if (!nzchar(python_path)) {
      python_path <- Sys.which("python")
    }
    if (!is.null(python_path) && consequiR_valid_python_bin(python_path)) {
      return(python_path)
    }
  }

  python_path <- consequIR_python_from_conda_run(lymphoma_graph_env_name)
  if (!is.null(python_path)) {
    return(python_path)
  }

  python_path <- consequIR_python_from_reticulate(lymphoma_graph_env_name)
  if (!is.null(python_path)) {
    return(python_path)
  }

  NULL
}

.onLoad <- function(libname, pkgname) {
  if (!is.null(getOption("conseguiR.python"))) {
    return(invisible())
  }

  python_path <- consequIR_find_python()
  if (!is.null(python_path)) {
    options(
      conseguiR.python = python_path,
      conseguiR.conda_env = lymphoma_graph_env_name
    )
  }

  conda_prefix <- Sys.getenv("CONDA_PREFIX")
  if (nzchar(conda_prefix) && identical(Sys.getenv("CONDA_DEFAULT_ENV"), lymphoma_graph_env_name)) {
    consequIR_prepend_path(conda_prefix)
  }

  invisible()
}
