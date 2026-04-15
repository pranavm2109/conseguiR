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
