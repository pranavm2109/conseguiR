## Shared basilisk-backed Python helpers for the installed/runtime wrapper path.

.conseguiR_basilisk_pkgname <- function() {
  if (requireNamespace("conseguiR", quietly = TRUE)) {
    return("conseguiR")
  }
  "base"
}

.conseguiR_basilisk_env <- local({
  env <- NULL
  env_pkgname <- NULL

  function() {
    current_pkgname <- .conseguiR_basilisk_pkgname()
    if (is.null(env) || !identical(env_pkgname, current_pkgname)) {
      env_args <- list(
        envname = "conseguiR_py",
        pkgname = current_pkgname,
        packages = c(
          "python=3.11",
          "numpy=1.26.4",
          "pandas=2.2.3"
        ),
        pip = c(
          "ortools==9.10.4067"
        )
      )
      env <<- do.call(basilisk::BasiliskEnvironment, env_args)
      env_pkgname <<- current_pkgname
    }

    env
  }
})

.conseguiR_direct_python_path <- function() {
  candidates <- character()

  configured_python <- Sys.getenv("CONSEGUIR_PYTHON", unset = "")
  if (nzchar(configured_python)) {
    candidates <- c(candidates, configured_python)
  }

  reticulate_python <- Sys.getenv("RETICULATE_PYTHON", unset = "")
  if (nzchar(reticulate_python)) {
    candidates <- c(candidates, reticulate_python)
  }

  conda_prefix <- Sys.getenv("CONDA_PREFIX", unset = "")
  if (nzchar(conda_prefix)) {
    candidates <- c(candidates, file.path(conda_prefix, "bin", "python"))
  }

  if (exists("consequIR_find_python", mode = "function")) {
    discovered <- tryCatch(consequIR_find_python(), error = function(e) NULL)
    if (is.character(discovered) && length(discovered) == 1L && nzchar(discovered)) {
      candidates <- c(candidates, discovered)
    }
  }

  python3_path <- Sys.which("python3")
  if (nzchar(python3_path)) {
    candidates <- c(candidates, python3_path)
  }

  python_path <- Sys.which("python")
  if (nzchar(python_path)) {
    candidates <- c(candidates, python_path)
  }

  candidates <- unique(candidates[nzchar(candidates)])
  valid <- Filter(
    function(path) is.character(path) && length(path) == 1L && file.exists(path) && !dir.exists(path),
    candidates
  )

  if (length(valid) < 1L) {
    return(NULL)
  }

  valid[[1]]
}

.conseguiR_run_python_module_direct <- function(
  script_path,
  module_name,
  function_name,
  config,
  config_class_name,
  result_field = "output_paths",
  python_path = .conseguiR_direct_python_path()
) {
  if (is.null(python_path) || !nzchar(python_path)) {
    stop("No configured Python interpreter is available for direct conseguiR execution.")
  }

  reticulate::use_python(python_path, required = TRUE)

  importlib <- reticulate::import("importlib.util", convert = FALSE)
  sys <- reticulate::import("sys", convert = FALSE)
  spec <- importlib$spec_from_file_location(module_name, script_path)
  module <- importlib$module_from_spec(spec)
  sys$modules[[module_name]] <- module
  spec$loader$exec_module(module)

  config_class <- module[[config_class_name]]
  run_fun <- module[[function_name]]
  config_obj <- do.call(config_class, config)
  result <- run_fun(config = config_obj)

  reticulate::py_to_r(result[[result_field]])
}

.conseguiR_basilisk_run_python_module <- function(
  script_path,
  module_name,
  function_name,
  config,
  config_class_name,
  result_field = "output_paths"
) {
  if (!file.exists(script_path)) {
    stop("Python script not found: ", script_path)
  }

  direct_python <- .conseguiR_direct_python_path()
  if (!is.null(direct_python)) {
    direct_result <- tryCatch(
      .conseguiR_run_python_module_direct(
        script_path = script_path,
        module_name = module_name,
        function_name = function_name,
        config = config,
        config_class_name = config_class_name,
        result_field = result_field,
        python_path = direct_python
      ),
      error = function(e) e
    )

    if (!inherits(direct_result, "error")) {
      return(direct_result)
    }
  }

  tryCatch(
    basilisk::basiliskRun(
      env = .conseguiR_basilisk_env(),
      fun = function(script_path, module_name, function_name, config, config_class_name, result_field) {
        importlib <- reticulate::import("importlib.util", convert = FALSE)
        sys <- reticulate::import("sys", convert = FALSE)
        spec <- importlib$spec_from_file_location(module_name, script_path)
        module <- importlib$module_from_spec(spec)
        sys$modules[[module_name]] <- module
        spec$loader$exec_module(module)

        config_class <- module[[config_class_name]]
        run_fun <- module[[function_name]]
        config_obj <- do.call(config_class, config)
        result <- run_fun(config = config_obj)

        reticulate::py_to_r(result[[result_field]])
      },
      script_path = normalizePath(script_path, winslash = "/", mustWork = TRUE),
      module_name = module_name,
      function_name = function_name,
      config = config,
      config_class_name = config_class_name,
      result_field = result_field
    ),
    error = function(e) {
      stop(
        "Python-backed conseguiR stage failed inside the managed basilisk environment. ",
        "This usually means the environment could not be provisioned ",
        "or one of the required Python modules is unavailable.\n",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
}
