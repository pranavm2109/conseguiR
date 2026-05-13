#' @keywords internal
.conseguiR_basilisk_pkgname <- function() {
    if (requireNamespace("conseguiR", quietly = TRUE)) {
        return("conseguiR")
    }
    "base"
}

#' @keywords internal
.conseguiR_basilisk_env <- local({
    function() {
        env <- .conseguiR_state$basilisk_env
        current_pkgname <- .conseguiR_basilisk_pkgname()
        cached_pkgname <- .conseguiR_state$basilisk_env_pkgname %||% NULL
        if (is.null(env) || !identical(cached_pkgname, current_pkgname)) {
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
            env <- do.call(basilisk::BasiliskEnvironment, env_args)
            .conseguiR_state$basilisk_env <- env
            .conseguiR_state$basilisk_env_pkgname <- current_pkgname
        }

        env
    }
})

#' @keywords internal
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

#' @keywords internal
.conseguiR_user_requested_direct_python <- function() {
    any(nzchar(c(
        Sys.getenv("CONSEGUIR_PYTHON", unset = ""),
        Sys.getenv("RETICULATE_PYTHON", unset = "")
    )))
}

#' @keywords internal
.conseguiR_write_python_config_tsv <- function(config, path) {
    stopifnot(is.list(config), is.character(path), length(path) == 1L)

    keys <- names(config)
    if (is.null(keys) || any(!nzchar(keys))) {
        stop("Python config must be a named list.")
    }

    infer_type <- function(x) {
        if (is.null(x)) {
            return("none")
        }
        if (is.logical(x)) {
            return("bool")
        }
        if (is.integer(x)) {
            return("int")
        }
        if (is.numeric(x)) {
            if (length(x) == 1L && is.finite(x) && identical(as.numeric(as.integer(x)), as.numeric(x))) {
                return("int")
            }
            return("float")
        }
        "str"
    }

    encode_value <- function(x, type_name) {
        if (is.null(x)) {
            return("")
        }
        if (type_name == "bool") {
            return(ifelse(isTRUE(x), "true", "false"))
        }
        as.character(x)
    }

    rows <- lapply(keys, function(key) {
        value <- config[[key]]
        type_name <- infer_type(value)
        data.frame(
            key = key,
            value = encode_value(value, type_name),
            type = type_name,
            stringsAsFactors = FALSE
        )
    })

    utils::write.table(
        do.call(rbind, rows),
        file = path,
        sep = "\t",
        quote = TRUE,
        row.names = FALSE,
        col.names = TRUE
    )

    path
}

#' @keywords internal
.conseguiR_read_python_output_tsv <- function(path) {
    dt <- utils::read.delim(
        path,
        sep = "\t",
        header = TRUE,
        stringsAsFactors = FALSE,
        check.names = FALSE
    )

    if (!all(c("key", "value") %in% names(dt))) {
        stop("Python output mapping file is malformed: ", path)
    }

    stats::setNames(as.list(dt$value), dt$key)
}

#' @keywords internal
.conseguiR_run_python_module_subprocess <- function(
  script_path,
  config,
  python_path
) {
    if (is.null(python_path) || !nzchar(python_path)) {
        stop("No configured Python interpreter is available for direct conseguiR execution.")
    }

    config_path <- tempfile("conseguiR_py_config_", fileext = ".tsv")
    output_path <- tempfile("conseguiR_py_output_", fileext = ".tsv")
    on.exit(unlink(c(config_path, output_path), force = TRUE), add = TRUE)

    .conseguiR_write_python_config_tsv(config, config_path)

    out <- tryCatch(
        system2(
            python_path,
            args = c(
                normalizePath(script_path, winslash = "/", mustWork = TRUE),
                "--config-tsv", config_path,
                "--output-paths-tsv", output_path
            ),
            stdout = TRUE,
            stderr = TRUE
        ),
        error = function(e) e
    )

    if (inherits(out, "error")) {
        stop(conditionMessage(out), call. = FALSE)
    }

    status <- attr(out, "status")
    if (!is.null(status) && status != 0L) {
        stop(
            "Direct Python subprocess execution failed.\n",
            paste(out, collapse = "\n"),
            call. = FALSE
        )
    }

    if (!file.exists(output_path)) {
        stop(
            "Direct Python subprocess did not produce the expected output mapping file.\n",
            paste(out, collapse = "\n"),
            call. = FALSE
        )
    }

    .conseguiR_read_python_output_tsv(output_path)
}

#' @keywords internal
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

#' @keywords internal
.conseguiR_basilisk_python_status <- function(modules = c("numpy", "pandas", "ortools")) {
    status <- list(
        python_path = NULL,
        python_ok = FALSE,
        python_modules = stats::setNames(rep(FALSE, length(modules)), modules),
        error = NULL,
        source = "basilisk"
    )

    direct_python <- .conseguiR_direct_python_path()
    if (!is.null(direct_python)) {
        probe <- tryCatch(
            {
                reticulate::use_python(direct_python, required = TRUE)
                importlib <- reticulate::import("importlib.util", convert = FALSE)
                found <- vapply(
                    modules,
                    function(mod) !is.null(importlib$find_spec(mod)),
                    logical(1)
                )
                list(
                    python_path = direct_python,
                    python_modules = found
                )
            },
            error = function(e) e
        )

        if (!inherits(probe, "error")) {
            status$python_path <- probe$python_path
            status$python_modules[names(probe$python_modules)] <- as.logical(probe$python_modules)
            status$python_ok <- is.character(status$python_path) && length(status$python_path) == 1L && nzchar(status$python_path)
            status$source <- "direct"
            return(status)
        }
    }

    probe <- tryCatch(
        basilisk::basiliskRun(
            env = .conseguiR_basilisk_env(),
            fun = function(modules) {
                importlib <- reticulate::import("importlib.util", convert = FALSE)
                sys <- reticulate::import("sys", convert = FALSE)
                found <- vapply(
                    modules,
                    function(mod) !is.null(importlib$find_spec(mod)),
                    logical(1)
                )

                list(
                    python_path = reticulate::py_to_r(sys$executable),
                    python_modules = found
                )
            },
            modules = modules
        ),
        error = function(e) e
    )

    if (inherits(probe, "error")) {
        status$error <- conditionMessage(probe)
        return(status)
    }

    status$python_path <- probe$python_path
    status$python_modules[names(probe$python_modules)] <- as.logical(probe$python_modules)
    status$python_ok <- is.character(status$python_path) && length(status$python_path) == 1L && nzchar(status$python_path)
    status
}

#' @keywords internal
.conseguiR_basilisk_warmup <- function(modules = c("numpy", "pandas", "ortools")) {
    status <- .conseguiR_basilisk_python_status(modules = modules)
    .conseguiR_state$basilisk_status <- status
    status
}

#' @keywords internal
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
            .conseguiR_run_python_module_subprocess(
                script_path = script_path,
                config = config,
                python_path = direct_python
            ),
            error = function(e) e
        )

        if (!inherits(direct_result, "error")) {
            return(direct_result)
        }

        if (isTRUE(.conseguiR_user_requested_direct_python())) {
            stop(
                "Python-backed conseguiR stage failed while using the explicitly configured direct Python interpreter.\n",
                conditionMessage(direct_result),
                call. = FALSE
            )
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
