#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

conseguiR_runtime_file <- function(relpath) {
  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(candidate)
  }

  pkg_path <- system.file(relpath, package = "conseguiR")
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }

  stop("Could not locate required runtime file: ", relpath)
}

sys.source(conseguiR_runtime_file("scripts/Internals/R/00_harmonise_and_validate_inputs.R"), envir = environment())

conseguiR_verbose_message <- function(verbose, ...) {
  if (isTRUE(verbose)) {
    message(...)
  }
}

conseguiR_verbose_cat <- function(lines, verbose) {
  if (isTRUE(verbose) && length(lines) > 0L) {
    cat(paste(lines, collapse = "\n"), "\n")
  }
}

conseguiR_magma_resolution_note <- function() {
  paste(
    "Install MAGMA separately and provide it via `magma_path`,",
    "`options(conseguiR.magma_path = \"/path/to/magma\")`,",
    "the `CONSEGUIR_MAGMA_PATH` environment variable, or by making `magma`",
    "available on your system PATH."
  )
}

autodiscovered_magma_candidates <- function() {
  roots <- unique(Filter(
    function(x) is.character(x) && length(x) == 1L && nzchar(x) && dir.exists(x),
    c(getwd(), Sys.getenv("HOME", unset = ""))
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

  unique(candidates[file.exists(candidates) & !dir.exists(candidates)])
}

resolve_magma_path <- function(magma_path = NULL, must_work = TRUE) {
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

  candidates <- unique(Filter(Negate(is.null), c(
    explicit_path,
    option_path,
    env_path,
    path_path,
    autodiscovered_magma_candidates()
  )))

  if (length(candidates) == 0L) {
    if (isTRUE(must_work)) {
      stop("Could not locate a MAGMA executable. ", conseguiR_magma_resolution_note(), call. = FALSE)
    }
    return(NULL)
  }

  for (candidate in candidates) {
    if (!file.exists(candidate)) {
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
      conseguiR_magma_resolution_note(),
      call. = FALSE
    )
  }

  NULL
}

default_magma_gwas_select_columns <- unique(c(
  "hm_rsid",
  "hm_variant_id",
  "hm_chrom",
  "hm_pos",
  "p_value",
  "rs_id",
  "rsid",
  "variant_id",
  "chromosome",
  "base_pair_location",
  "p",
  "P",
  "pval"
))

read_gwas_sumstats_for_magma <- function(
  gwas_sumstats,
  select_columns = default_magma_gwas_select_columns
) {
  if (is.character(gwas_sumstats) && length(gwas_sumstats) == 1L) {
    if (!file.exists(gwas_sumstats)) {
      stop("GWAS summary statistics file does not exist: ", gwas_sumstats)
    }

    header_dt <- fread(gwas_sumstats, nrows = 0L, showProgress = FALSE)
    available_columns <- intersect(select_columns, names(header_dt))

    if (length(available_columns) == 0L) {
      stop(
        "GWAS summary statistics file does not contain any recognized MAGMA columns. ",
        "Available columns were: ", paste(names(header_dt), collapse = ", ")
      )
    }

    return(fread(
      gwas_sumstats,
      select = available_columns,
      showProgress = FALSE
    ))
  }

  as.data.table(gwas_sumstats)
}

prepare_magma_input_files <- function(
  gwas_sumstats,
  output_prefix,
  snp_loc_path = paste0(output_prefix, ".snp_loc.tsv"),
  pval_path = paste0(output_prefix, ".pval.tsv"),
  reuse_existing = TRUE,
  verbose = FALSE
) {
  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(snp_loc_path), recursive = TRUE, showWarnings = FALSE)
  dir.create(dirname(pval_path), recursive = TRUE, showWarnings = FALSE)

  existing_ready <- all(file.exists(c(snp_loc_path, pval_path))) &&
    file.info(snp_loc_path)$size > 0 &&
    file.info(pval_path)$size > 0

  if (isTRUE(reuse_existing) && existing_ready) {
    conseguiR_verbose_message(verbose, "Reusing existing MAGMA GWAS cache.")
    return(list(
      snp_loc_path = snp_loc_path,
      pval_path = pval_path,
      reused_existing = TRUE
    ))
  }

  conseguiR_verbose_message(verbose, "Reading GWAS summary statistics for MAGMA preprocessing...")
  gwas_loaded <- read_gwas_sumstats_for_magma(gwas_sumstats)
  conseguiR_verbose_message(verbose, "Validating GWAS columns required for MAGMA...")
  gwas_validated <- validate_gwas_sumstats(gwas_loaded)
  conseguiR_verbose_message(verbose, "Constructing MAGMA SNP-location and p-value tables...")
  magma_inputs <- prepare_magma_input(gwas_validated)

  if (nrow(magma_inputs$snp_loc) == 0L) {
    stop("Prepared MAGMA snp-loc table is empty.")
  }
  if (nrow(magma_inputs$pval) == 0L) {
    stop("Prepared MAGMA pval table is empty.")
  }
  if (any(is.na(magma_inputs$snp_loc$SNP) | magma_inputs$snp_loc$SNP == "")) {
    stop("Prepared MAGMA snp-loc table contains missing SNP identifiers.")
  }
  if (any(is.na(magma_inputs$snp_loc$CHR) | magma_inputs$snp_loc$CHR == "")) {
    stop("Prepared MAGMA snp-loc table contains missing chromosomes.")
  }
  if (any(is.na(magma_inputs$snp_loc$POS))) {
    stop("Prepared MAGMA snp-loc table contains missing positions.")
  }
  if (any(is.na(magma_inputs$pval$SNP) | magma_inputs$pval$SNP == "")) {
    stop("Prepared MAGMA pval table contains missing SNP identifiers.")
  }
  if (any(is.na(magma_inputs$pval$P))) {
    stop("Prepared MAGMA pval table contains missing p-values.")
  }

  conseguiR_verbose_message(verbose, "Writing MAGMA SNP-location cache...")
  fwrite(magma_inputs$snp_loc, snp_loc_path, sep = "\t", col.names = FALSE)
  conseguiR_verbose_message(verbose, "Writing MAGMA p-value cache...")
  fwrite(magma_inputs$pval, pval_path, sep = "\t", col.names = TRUE)

  if (!file.exists(snp_loc_path) || file.info(snp_loc_path)$size <= 0) {
    stop("Failed to write a usable MAGMA snp-loc file: ", snp_loc_path)
  }
  if (!file.exists(pval_path) || file.info(pval_path)$size <= 0) {
    stop("Failed to write a usable MAGMA pval file: ", pval_path)
  }

  list(
    snp_loc_path = snp_loc_path,
    pval_path = pval_path,
    reused_existing = FALSE
  )
}

prepare_magma_gwas_cache <- function(
  gwas_sumstats,
  cache_prefix,
  snp_loc_path = paste0(cache_prefix, ".snp_loc.tsv"),
  pval_path = paste0(cache_prefix, ".pval.tsv"),
  reuse_existing = TRUE,
  verbose = FALSE
) {
  prepare_magma_input_files(
    gwas_sumstats = gwas_sumstats,
    output_prefix = cache_prefix,
    snp_loc_path = snp_loc_path,
    pval_path = pval_path,
    reuse_existing = reuse_existing,
    verbose = verbose
  )
}

run_magma_step1_annotation <- function(
  gwas_sumstats,
  gene_loc_path,
  output_prefix,
  magma_path = NULL,
  snp_loc_path = paste0(output_prefix, ".snp_loc.tsv"),
  pval_path = paste0(output_prefix, ".pval.tsv"),
  annotation_window = NULL,
  filter_path = NULL,
  ignore_strand = FALSE,
  nonhuman = FALSE,
  reuse_prepared_inputs = TRUE,
  reuse_existing_annotation = FALSE,
  extra_args = character(),
  verbose = FALSE
) {
  magma_path <- resolve_magma_path(magma_path)

  if (!file.exists(gene_loc_path)) {
    stop("Gene/regulatory location file does not exist: ", gene_loc_path)
  }

  prepared_inputs <- prepare_magma_input_files(
    gwas_sumstats = gwas_sumstats,
    output_prefix = output_prefix,
    snp_loc_path = snp_loc_path,
    pval_path = pval_path,
    reuse_existing = reuse_prepared_inputs
  )

  snp_loc_path <- prepared_inputs$snp_loc_path
  pval_path <- prepared_inputs$pval_path
  annot_path <- paste0(output_prefix, ".genes.annot")

  existing_annotation_ready <- file.exists(annot_path) &&
    file.info(annot_path)$size > 0

  if (isTRUE(reuse_existing_annotation) && existing_annotation_ready) {
    return(list(
      status = 0L,
      stdout_stderr = "Reused existing MAGMA annotation output.",
      command = NULL,
      output_prefix = output_prefix,
      snp_loc_path = snp_loc_path,
      pval_path = pval_path,
      gene_loc_path = gene_loc_path,
      annot_path = annot_path,
      magma_path = magma_path,
      args = NULL,
      reused_existing_inputs = prepared_inputs$reused_existing,
      reused_existing_annotation = TRUE
    ))
  }

  annotate_flag <- "--annotate"
  annotate_modifiers <- character()

  if (!is.null(annotation_window)) {
    if (length(annotation_window) == 1L) {
      annotate_modifiers <- c(annotate_modifiers, paste0("window=", annotation_window))
    } else if (length(annotation_window) == 2L) {
      annotate_modifiers <- c(annotate_modifiers, paste0("window=", annotation_window[[1]], ",", annotation_window[[2]]))
    } else {
      stop("annotation_window must be NULL, length 1, or length 2.")
    }
  }

  if (!is.null(filter_path)) {
    if (!file.exists(filter_path)) {
      stop("filter_path does not exist: ", filter_path)
    }
    annotate_modifiers <- c(annotate_modifiers, paste0("filter=", filter_path))
  }

  if (isTRUE(ignore_strand)) {
    annotate_modifiers <- c(annotate_modifiers, "ignore-strand")
  }

  if (isTRUE(nonhuman)) {
    annotate_modifiers <- c(annotate_modifiers, "nonhuman")
  }

  if (length(annotate_modifiers) > 0) {
    annotate_flag <- paste0("--annotate ", paste(annotate_modifiers, collapse = " "))
  }

  args <- c(
    strsplit(annotate_flag, " ", fixed = TRUE)[[1]],
    "--snp-loc", snp_loc_path,
    "--gene-loc", gene_loc_path,
    "--out", output_prefix,
    extra_args
  )

  command_string <- paste(c(shQuote(magma_path), shQuote(args)), collapse = " ")
  conseguiR_verbose_message(verbose, "Running MAGMA step 1 annotation command:")
  conseguiR_verbose_message(verbose, command_string)

  status <- system2(
    command = magma_path,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )

  exit_status <- attr(status, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }

  if (!identical(exit_status, 0L)) {
    stop(
      "MAGMA step 1 annotation failed with status ", exit_status, ".\n",
      paste(status, collapse = "\n")
    )
  }

  conseguiR_verbose_cat(status, verbose)

  list(
    status = exit_status,
    stdout_stderr = status,
    command = command_string,
    output_prefix = output_prefix,
    snp_loc_path = snp_loc_path,
    pval_path = pval_path,
    gene_loc_path = gene_loc_path,
    annot_path = annot_path,
    magma_path = magma_path,
    args = args,
    reused_existing_inputs = prepared_inputs$reused_existing,
    reused_existing_annotation = FALSE
  )
}

run_magma_step2_gene_analysis <- function(
  gene_annot_path,
  pval_path,
  reference_bfile,
  output_prefix,
  sample_size = NULL,
  sample_size_col = NULL,
  magma_path = NULL,
  gene_model = NULL,
  genes_only = TRUE,
  pval_use = c("SNP", "P"),
  pval_duplicate = NULL,
  bfile_synonyms = NULL,
  bfile_synonym_dup = NULL,
  reuse_existing_analysis = FALSE,
  extra_args = character(),
  verbose = FALSE
) {
  magma_path <- resolve_magma_path(magma_path)
  if (!file.exists(gene_annot_path)) {
    stop("MAGMA gene annotation file does not exist: ", gene_annot_path)
  }
  if (!file.exists(pval_path)) {
    stop("MAGMA p-value file does not exist: ", pval_path)
  }
  if (!file.exists(paste0(reference_bfile, ".bed"))) {
    stop("Reference bfile is missing .bed: ", paste0(reference_bfile, ".bed"))
  }
  if (!file.exists(paste0(reference_bfile, ".bim"))) {
    stop("Reference bfile is missing .bim: ", paste0(reference_bfile, ".bim"))
  }
  if (!file.exists(paste0(reference_bfile, ".fam"))) {
    stop("Reference bfile is missing .fam: ", paste0(reference_bfile, ".fam"))
  }

  if (is.null(sample_size) && is.null(sample_size_col)) {
    stop("Provide either `sample_size` or `sample_size_col` for MAGMA step 2.")
  }
  if (!is.null(sample_size) && !is.null(sample_size_col)) {
    stop("Provide only one of `sample_size` or `sample_size_col`, not both.")
  }

  dir.create(dirname(output_prefix), recursive = TRUE, showWarnings = FALSE)
  genes_out_path <- paste0(output_prefix, ".genes.out")
  genes_raw_path <- paste0(output_prefix, ".genes.raw")

  existing_analysis_ready <- file.exists(genes_out_path) &&
    file.info(genes_out_path)$size > 0

  if (isTRUE(reuse_existing_analysis) && existing_analysis_ready) {
    return(list(
      status = 0L,
      stdout_stderr = "Reused existing MAGMA gene analysis output.",
      command = NULL,
      output_prefix = output_prefix,
      gene_annot_path = gene_annot_path,
      pval_path = pval_path,
      reference_bfile = reference_bfile,
      genes_out_path = genes_out_path,
      genes_raw_path = genes_raw_path,
      magma_path = magma_path,
      args = NULL,
      reused_existing_analysis = TRUE
    ))
  }

  pval_modifier <- character()
  if (!is.null(pval_use)) {
    if (length(pval_use) != 2L) {
      stop("`pval_use` must have length 2: SNP column and p-value column.")
    }
    pval_modifier <- c(pval_modifier, paste0("use=", paste(pval_use, collapse = ",")))
  }
  if (!is.null(sample_size)) {
    pval_modifier <- c(pval_modifier, paste0("N=", sample_size))
  }
  if (!is.null(sample_size_col)) {
    pval_modifier <- c(pval_modifier, paste0("ncol=", sample_size_col))
  }
  if (!is.null(pval_duplicate)) {
    pval_modifier <- c(pval_modifier, paste0("duplicate=", pval_duplicate))
  }

  bfile_modifier <- character()
  if (!is.null(bfile_synonyms)) {
    bfile_modifier <- c(bfile_modifier, paste0("synonyms=", bfile_synonyms))
  }
  if (!is.null(bfile_synonym_dup)) {
    bfile_modifier <- c(bfile_modifier, paste0("synonym-dup=", bfile_synonym_dup))
  }

  args <- c(
    "--bfile", reference_bfile, bfile_modifier,
    "--gene-annot", gene_annot_path,
    "--pval", pval_path,
    pval_modifier,
    "--out", output_prefix
  )

  if (!is.null(gene_model)) {
    args <- c(args, "--gene-model", gene_model)
  }

  if (isTRUE(genes_only)) {
    args <- c(args, "--genes-only")
  }

  args <- c(args, extra_args)

  command_string <- paste(c(shQuote(magma_path), shQuote(args)), collapse = " ")
  conseguiR_verbose_message(verbose, "Running MAGMA step 2 gene analysis command:")
  conseguiR_verbose_message(verbose, command_string)

  status <- system2(
    command = magma_path,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )

  exit_status <- attr(status, "status")
  if (is.null(exit_status)) {
    exit_status <- 0L
  }

  if (!identical(exit_status, 0L)) {
    stop(
      "MAGMA step 2 gene analysis failed with status ", exit_status, ".\n",
      paste(status, collapse = "\n")
    )
  }

  conseguiR_verbose_cat(status, verbose)

  if (!file.exists(genes_out_path) || file.info(genes_out_path)$size <= 0) {
    stop("MAGMA step 2 did not produce a usable .genes.out file: ", genes_out_path)
  }

  list(
    status = exit_status,
    stdout_stderr = status,
    command = command_string,
    output_prefix = output_prefix,
    gene_annot_path = gene_annot_path,
    pval_path = pval_path,
    reference_bfile = reference_bfile,
    genes_out_path = genes_out_path,
    genes_raw_path = genes_raw_path,
    magma_path = magma_path,
    args = args,
    reused_existing_analysis = FALSE
  )
}

extract_magma_feature_zstat <- function(
  genes_out_path,
  feature_type = c("gene", "regulatory_element"),
  output_path = NULL,
  cleanup_paths = NULL
) {
  feature_type <- match.arg(feature_type)

  if (!file.exists(genes_out_path)) {
    stop("MAGMA gene output file does not exist: ", genes_out_path)
  }

  genes_dt <- fread(genes_out_path)

  required_cols <- c("GENE", "ZSTAT")
  missing_cols <- setdiff(required_cols, names(genes_dt))
  if (length(missing_cols) > 0) {
    stop(
      "MAGMA gene output is missing required columns: ",
      paste(missing_cols, collapse = ", ")
    )
  }

  id_col <- if ("GENE_NAME" %in% names(genes_dt)) "GENE_NAME" else "GENE"
  out <- genes_dt[, .(
    feature_id = as.character(get(id_col)),
    zstat = as.numeric(ZSTAT)
  )]

  if (!is.null(output_path)) {
    dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)
    fwrite(out, output_path, sep = "\t")
  }

  if (!is.null(cleanup_paths) && length(cleanup_paths) > 0L) {
    cleanup_paths <- unique(cleanup_paths[file.exists(cleanup_paths)])
    if (length(cleanup_paths) > 0L) {
      unlink(cleanup_paths, force = TRUE)
    }
  }

  out
}

extract_magma_zstat <- function(
  genes_out_path,
  output_path = NULL,
  cleanup_paths = NULL
) {
  out <- extract_magma_feature_zstat(
    genes_out_path = genes_out_path,
    feature_type = "gene",
    output_path = output_path,
    cleanup_paths = cleanup_paths
  )

  setnames(out, "feature_id", "gene_id")

  if (!is.null(output_path)) {
    fwrite(out, output_path, sep = "\t")
  }

  out
}

run_magma_feature_scoring_pipeline <- function(
  gwas_sumstats,
  feature_loc_path,
  reference_bfile,
  output_prefix,
  feature_type = c("gene", "regulatory_element"),
  sample_size = NULL,
  sample_size_col = NULL,
  magma_path = NULL,
  annotation_window = NULL,
  filter_path = NULL,
  ignore_strand = FALSE,
  nonhuman = FALSE,
  gene_model = NULL,
  genes_only = TRUE,
  pval_use = NULL,
  pval_duplicate = NULL,
  bfile_synonyms = NULL,
  bfile_synonym_dup = NULL,
  magma_gwas_cache_prefix = NULL,
  reuse_existing_gwas_cache = TRUE,
  reuse_existing_annotation = FALSE,
  reuse_existing_analysis = FALSE,
  step1_extra_args = character(),
  step2_extra_args = character(),
  keep_intermediates = FALSE,
  zstat_output_path = paste0(output_prefix, ".zstat.tsv"),
  verbose = FALSE
) {
  feature_type <- match.arg(feature_type)
  pb <- if (isTRUE(verbose)) utils::txtProgressBar(min = 0, max = 4, style = 3) else NULL
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)
  step1_prefix <- paste0(output_prefix, ".step1")
  step2_prefix <- paste0(output_prefix, ".step2")
  shared_snp_loc_path <- NULL
  shared_pval_path <- NULL

  if (!is.null(magma_gwas_cache_prefix)) {
    if (!is.null(pb)) utils::setTxtProgressBar(pb, 1)
    conseguiR_verbose_message(verbose, "Preparing shared MAGMA GWAS cache...")
    gwas_cache <- prepare_magma_gwas_cache(
      gwas_sumstats = gwas_sumstats,
      cache_prefix = magma_gwas_cache_prefix,
      reuse_existing = reuse_existing_gwas_cache,
      verbose = verbose
    )
    shared_snp_loc_path <- gwas_cache$snp_loc_path
    shared_pval_path <- gwas_cache$pval_path
  }

  if (is.null(shared_snp_loc_path)) {
    shared_snp_loc_path <- paste0(step1_prefix, ".snp_loc.tsv")
  }
  if (is.null(shared_pval_path)) {
    shared_pval_path <- paste0(step1_prefix, ".pval.tsv")
  }

  if (is.null(magma_gwas_cache_prefix)) {
    if (!is.null(pb)) utils::setTxtProgressBar(pb, 1)
    conseguiR_verbose_message(verbose, "Preparing MAGMA inputs for step 1...")
  }

  step1 <- run_magma_step1_annotation(
    gwas_sumstats = gwas_sumstats,
    gene_loc_path = feature_loc_path,
    output_prefix = step1_prefix,
    magma_path = magma_path,
    snp_loc_path = shared_snp_loc_path,
    pval_path = shared_pval_path,
    annotation_window = annotation_window,
    filter_path = filter_path,
    ignore_strand = ignore_strand,
    nonhuman = nonhuman,
    reuse_prepared_inputs = TRUE,
    reuse_existing_annotation = reuse_existing_annotation,
    extra_args = step1_extra_args,
    verbose = verbose
  )
  if (!is.null(pb)) utils::setTxtProgressBar(pb, 2)
  conseguiR_verbose_message(verbose, "MAGMA step 1 complete.")

  step2 <- run_magma_step2_gene_analysis(
    gene_annot_path = step1$annot_path,
    pval_path = step1$pval_path,
    reference_bfile = reference_bfile,
    output_prefix = step2_prefix,
    sample_size = sample_size,
    sample_size_col = sample_size_col,
    magma_path = magma_path,
    gene_model = gene_model,
    genes_only = genes_only,
    pval_use = pval_use,
    pval_duplicate = pval_duplicate,
    bfile_synonyms = bfile_synonyms,
    bfile_synonym_dup = bfile_synonym_dup,
    reuse_existing_analysis = reuse_existing_analysis,
    extra_args = step2_extra_args,
    verbose = verbose
  )
  if (!is.null(pb)) utils::setTxtProgressBar(pb, 3)
  conseguiR_verbose_message(verbose, "MAGMA step 2 complete.")

  cleanup_paths <- NULL
  if (!isTRUE(keep_intermediates)) {
    cleanup_paths <- c(
      if (is.null(magma_gwas_cache_prefix)) step1$snp_loc_path,
      if (is.null(magma_gwas_cache_prefix)) step1$pval_path,
      step1$annot_path,
      paste0(step1$output_prefix, ".log"),
      paste0(step1$output_prefix, ".log.tmp"),
      step2$genes_out_path,
      step2$genes_raw_path,
      paste0(step2$output_prefix, ".log"),
      paste0(step2$output_prefix, ".log.tmp")
    )
  }

  zstat <- extract_magma_feature_zstat(
    genes_out_path = step2$genes_out_path,
    feature_type = feature_type,
    output_path = zstat_output_path,
    cleanup_paths = cleanup_paths
  )

  if (!is.null(pb)) utils::setTxtProgressBar(pb, 4)
  conseguiR_verbose_message(verbose, "MAGMA z-score extraction complete.")

  list(
    step1 = step1,
    step2 = step2,
    zstat = zstat,
    zstat_output_path = zstat_output_path,
    feature_type = feature_type,
    magma_gwas_cache_prefix = magma_gwas_cache_prefix
  )
}

run_magma_gene_scoring_pipeline <- function(
  gwas_sumstats,
  gene_loc_path,
  reference_bfile,
  output_prefix,
  sample_size = NULL,
  sample_size_col = NULL,
  magma_path = NULL,
  annotation_window = NULL,
  filter_path = NULL,
  ignore_strand = FALSE,
  nonhuman = FALSE,
  gene_model = NULL,
  genes_only = TRUE,
  pval_use = NULL,
  pval_duplicate = NULL,
  bfile_synonyms = NULL,
  bfile_synonym_dup = NULL,
  magma_gwas_cache_prefix = NULL,
  reuse_existing_gwas_cache = TRUE,
  reuse_existing_annotation = FALSE,
  reuse_existing_analysis = FALSE,
  step1_extra_args = character(),
  step2_extra_args = character(),
  keep_intermediates = FALSE,
  zstat_output_path = paste0(output_prefix, ".zstat.tsv")
) {
  run_magma_feature_scoring_pipeline(
    gwas_sumstats = gwas_sumstats,
    feature_loc_path = gene_loc_path,
    reference_bfile = reference_bfile,
    output_prefix = output_prefix,
    feature_type = "gene",
    sample_size = sample_size,
    sample_size_col = sample_size_col,
    magma_path = magma_path,
    annotation_window = annotation_window,
    filter_path = filter_path,
    ignore_strand = ignore_strand,
    nonhuman = nonhuman,
    gene_model = gene_model,
    genes_only = genes_only,
    pval_use = pval_use,
    pval_duplicate = pval_duplicate,
    bfile_synonyms = bfile_synonyms,
    bfile_synonym_dup = bfile_synonym_dup,
    magma_gwas_cache_prefix = magma_gwas_cache_prefix,
    reuse_existing_gwas_cache = reuse_existing_gwas_cache,
    reuse_existing_annotation = reuse_existing_annotation,
    reuse_existing_analysis = reuse_existing_analysis,
    step1_extra_args = step1_extra_args,
    step2_extra_args = step2_extra_args,
    keep_intermediates = keep_intermediates,
    zstat_output_path = zstat_output_path
  )
}
