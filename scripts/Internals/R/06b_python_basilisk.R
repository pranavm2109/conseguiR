## Shared basilisk-backed Python helpers for the installed/runtime wrapper path.

.conseguiR_basilisk_env <- local({
  env <- NULL

  function() {
    if (is.null(env)) {
      env <<- basilisk::BasiliskEnvironment(
        envname = "conseguiR_py",
        pkgname = "conseguiR",
        packages = c(
          "python=3.11",
          "numpy=1.26.4",
          "pandas=2.2.3"
        ),
        pip = c(
          "ortools==9.10.4067"
        )
      )
    }

    env
  }
})

.conseguiR_direct_python_path <- function() {
  if (!exists("consequIR_find_python", mode = "function")) {
    return(NULL)
  }

  path <- tryCatch(consequIR_find_python(), error = function(e) NULL)
  if (is.null(path) || !is.character(path) || length(path) != 1L || !nzchar(path)) {
    return(NULL)
  }

  path
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
