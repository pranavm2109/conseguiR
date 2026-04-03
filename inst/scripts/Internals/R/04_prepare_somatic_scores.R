#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
})

conseguiR_runtime_file <- function(relpath) {
  pkg_path <- system.file(relpath, package = "conseguiR")
  if (nzchar(pkg_path) && file.exists(pkg_path)) {
    return(pkg_path)
  }

  candidate <- file.path(getwd(), relpath)
  if (file.exists(candidate)) {
    return(candidate)
  }

  stop("Could not locate required runtime file: ", relpath)
}

sys.source(conseguiR_runtime_file("scripts/Internals/R/00_harmonise_and_validate_inputs.R"), envir = environment())

# Somatic scoring design:
# - gene-level scores come from dndscv
# - regulatory-element-level scores come from fishHook
# - final output should be a simple score table with identifier and zstat

compute_signed_z_from_p <- function(p_value, effect_direction = NULL, min_p = 1e-300) {
  p_value <- as.numeric(p_value)
  p_value[p_value <= 0] <- min_p
  p_value[p_value > 1] <- NA_real_

  zstat <- stats::qnorm(1 - (p_value / 2))

  if (!is.null(effect_direction)) {
    zstat <- zstat * sign(as.numeric(effect_direction))
  }

  zstat
}

detect_dndscv_refdb_chr_style <- function(refdb) {
  if (is.null(refdb) || !is.character(refdb) || length(refdb) != 1L || !file.exists(refdb)) {
    return(NULL)
  }

  env <- new.env(parent = emptyenv())
  load(refdb, envir = env)

  if (!exists("RefCDS", envir = env, inherits = FALSE)) {
    return(NULL)
  }

  refcds <- get("RefCDS", envir = env, inherits = FALSE)
  if (!is.list(refcds) || length(refcds) == 0L) {
    return(NULL)
  }

  chr_values <- vapply(
    refcds,
    FUN.VALUE = character(1),
    FUN = function(x) {
      chr <- x$chr %||% NA_character_
      chr[[1]]
    }
  )
  chr_values <- chr_values[!is.na(chr_values) & nzchar(chr_values)]

  if (length(chr_values) == 0L) {
    return(NULL)
  }

  if (all(grepl("^chr", chr_values))) {
    return("UCSC")
  }

  "NCBI"
}

harmonize_dndscv_chr_style <- function(mutation_dt, refdb) {
  style <- detect_dndscv_refdb_chr_style(refdb)

  if (is.null(style)) {
    return(mutation_dt)
  }

  out <- copy(mutation_dt)

  if (identical(style, "UCSC")) {
    out[, chr := ifelse(grepl("^chr", chr), chr, paste0("chr", chr))]
  } else {
    out[, chr := sub("^chr", "", chr)]
  }

  out
}

extract_dndscv_gene_scores <- function(
  dndscv_result,
  gene_col = c("gene_name", "gene", "symbol"),
  p_col = c("pallsubs_cv", "pglobal_cv", "pmis_cv"),
  effect_col = c("wall_cv", "wmis_cv")
) {
  dt <- as.data.table(dndscv_result)

  gene_col <- pick_first_existing_column(dt, gene_col, "dndscv gene identifier")
  p_col <- pick_first_existing_column(dt, p_col, "dndscv p-value")
  effect_col <- intersect(effect_col, names(dt))
  effect_col <- if (length(effect_col) > 0) effect_col[[1]] else NULL

  out <- dt[, .(
    gene_id = as.character(get(gene_col)),
    zstat = compute_signed_z_from_p(
      p_value = get(p_col),
      effect_direction = if (!is.null(effect_col)) get(effect_col) - 1 else NULL
    )
  )]

  out <- unique(out[!is.na(gene_id) & gene_id != "" & !is.na(zstat)])
  out[, gene_id := toupper(gene_id)]
  out
}

extract_fishhook_reg_scores <- function(
  fishhook_result,
  reg_col = c("id", "reg_elem_id", "element_id", "name"),
  p_col = c("p", "pvalue", "p_value", "p.val"),
  effect_col = c("effectsize", "zscore", "z", "coef", "beta", "observed_over_expected")
) {
  dt <- as.data.table(fishhook_result)

  reg_col <- pick_first_existing_column(dt, reg_col, "fishHook regulatory element identifier")
  p_col <- pick_first_existing_column(dt, p_col, "fishHook p-value")
  effect_col <- intersect(effect_col, names(dt))
  effect_col <- if (length(effect_col) > 0) effect_col[[1]] else NULL

  out <- dt[, .(
    reg_elem_id = as.character(get(reg_col)),
    zstat = if (!is.null(effect_col) && effect_col %in% c("zscore", "z")) {
      as.numeric(get(effect_col))
    } else {
      compute_signed_z_from_p(
        p_value = get(p_col),
        effect_direction = if (!is.null(effect_col)) get(effect_col) else NULL
      )
    }
  )]

  out <- unique(out[!is.na(reg_elem_id) & reg_elem_id != "" & !is.na(zstat)])
  out
}

make_fishhook_event_granges <- function(maf) {
  dt <- prepare_fishhook_input(maf)

  gr <- GRanges(
    seqnames = dt$Chromosome,
    ranges = IRanges(start = dt$Start_Position, end = dt$End_Position),
    strand = "*"
  )

  mcols(gr)$Tumor_Sample_Barcode <- dt$Tumor_Sample_Barcode
  mcols(gr)$Reference_Allele <- dt$Reference_Allele
  mcols(gr)$Tumor_Seq_Allele2 <- dt$Tumor_Seq_Allele2
  seqlevelsStyle(gr) <- "UCSC"
  gr
}

make_fishhook_hypothesis_granges <- function(reg_ref_path) {
  reg_gr <- validate_regulatory_element_reference(reg_ref_path)

  if (!"reg_elem_id" %in% names(mcols(reg_gr))) {
    stop("Regulatory element reference must contain `reg_elem_id` metadata.")
  }

  seqlevelsStyle(reg_gr) <- "UCSC"
  reg_gr[, "reg_elem_id"]
}

make_fishhook_eligible_hg38 <- function() {
  if (!requireNamespace("BSgenome.Hsapiens.UCSC.hg38", quietly = TRUE)) {
    stop("Package `BSgenome.Hsapiens.UCSC.hg38` is required to construct the default hg38 eligible territory.")
  }

  si <- seqinfo(BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38)
  si <- keepStandardChromosomes(si, pruning.mode = "coarse")

  eligible_gr <- GRanges(
    seqnames = seqnames(si),
    ranges = IRanges(start = 1L, end = seqlengths(si)),
    strand = "*"
  )

  seqinfo(eligible_gr) <- si
  genome(eligible_gr) <- "hg38"
  eligible_gr
}

make_fishhook_covariates_from_data <- function(covariate_data, reg_gr) {
  cov_dt <- as.data.table(covariate_data)

  required_cols <- c("reg_elem_id", "n_samples_mut_norm")
  missing_cols <- setdiff(required_cols, names(cov_dt))
  if (length(missing_cols) > 0) {
    stop("fishHook covariate data is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  reg_ids <- as.character(mcols(reg_gr)$reg_elem_id)
  cov_dt[, reg_elem_id := as.character(reg_elem_id)]
  cov_dt <- unique(cov_dt, by = "reg_elem_id")

  keep_idx <- reg_ids %in% cov_dt$reg_elem_id

  if (!any(keep_idx)) {
    stop("No regulatory elements in the hypotheses could be matched to covariate data.")
  }

  reg_gr <- reg_gr[keep_idx]
  reg_ids <- as.character(mcols(reg_gr)$reg_elem_id)
  match_idx <- match(reg_ids, cov_dt$reg_elem_id)

  reg_lengths <- width(reg_gr)
  aid_per_base <- cov_dt$n_samples_mut_norm[match_idx] / reg_lengths

  mcols(reg_gr)$reg_elem_length <- as.numeric(reg_lengths)
  mcols(reg_gr)$aid_per_base <- as.numeric(aid_per_base)

  cov_len_obj <- fishHook::Cov(
    data = reg_gr,
    field = "reg_elem_length",
    name = "RE_length",
    na.rm = TRUE
  )
  cov_aid_obj <- fishHook::Cov(
    data = reg_gr,
    field = "aid_per_base",
    name = "AID_per_base",
    na.rm = TRUE
  )

  list(
    covariates = c(cov_len_obj, cov_aid_obj),
    hypotheses = reg_gr
  )
}

instantiate_fishhook_covariates <- function(covariates) {
  if (is.null(covariates) || length(covariates) == 0L) {
    return(NULL)
  }

  if (!is.list(covariates)) {
    stop("`covariates` must be a list of fishHook Cov objects or covariate specifications.")
  }

  out <- vector("list", length(covariates))

  for (i in seq_along(covariates)) {
    cov_i <- covariates[[i]]

    if (inherits(cov_i, "Covariate")) {
      out[[i]] <- cov_i
      next
    }

    if (!is.list(cov_i) || is.null(cov_i$data)) {
      stop("Each covariate specification must be a list with at least a `data` element.")
    }

    cov_name <- cov_i$name %||% paste0("covariate_", i)
    cov_field <- cov_i$field %||% NULL

    if (!inherits(cov_i$data, "GRanges")) {
      stop("Each covariate `data` element must be a GRanges object.")
    }

    if (is.null(cov_field)) {
      out[[i]] <- fishHook::Cov(data = cov_i$data, name = cov_name)
    } else {
      out[[i]] <- fishHook::Cov(data = cov_i$data, field = cov_field, name = cov_name)
    }
  }

  out
}

run_dndscv_gene_scoring <- function(
  maf,
  refdb,
  cv = NULL,
  max_muts_per_gene_per_sample = 6L,
  max_coding_muts_per_sample = 5000L,
  ...
) {
  maf_ready <- prepare_dndscv_input(maf)

  if (!requireNamespace("dndscv", quietly = TRUE)) {
    stop("Package `dndscv` is not installed. Install it or pass a precomputed result.")
  }

  if (missing(refdb) || is.null(refdb) || !nzchar(refdb)) {
    stop("`refdb` is required for dndscv scoring. Provide the hg38 RefCDS path used for this MAF.")
  }

  if (!file.exists(refdb)) {
    stop("dndscv `refdb` file does not exist: ", refdb)
  }

  maf_ready <- maf_ready[!chr %in% c("MT", "M", "chrM", "chrMT")]
  maf_ready <- harmonize_dndscv_chr_style(maf_ready, refdb = refdb)

  dndscv_args <- list(
    mutations = as.data.frame(maf_ready),
    cv = cv,
    refdb = refdb,
    max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,
    max_coding_muts_per_sample = max_coding_muts_per_sample,
    ...
  )

  dndscv_fit <- do.call(dndscv::dndscv, dndscv_args)

  if (!"sel_cv" %in% names(dndscv_fit)) {
    stop("dndscv result does not contain `sel_cv`.")
  }

  extract_dndscv_gene_scores(dndscv_fit$sel_cv)
}

run_fishhook_reg_scoring <- function(
  maf,
  reg_ref_path,
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  idcol = "Tumor_Sample_Barcode",
  ...
) {
  event_gr <- make_fishhook_event_granges(maf)
  hypothesis_gr <- make_fishhook_hypothesis_granges(reg_ref_path)

  if (!requireNamespace("fishHook", quietly = TRUE)) {
    stop("Package `fishHook` is not installed. Install it or pass a precomputed result.")
  }

  if (!requireNamespace("gUtils", quietly = TRUE)) {
    stop("Package `gUtils` is required for fishHook scoring in this environment.")
  }

  if (!"package:gUtils" %in% search()) {
    suppressPackageStartupMessages(
      library("gUtils", character.only = TRUE)
    )
  }

  if (is.null(eligible_gr)) {
    eligible_gr <- make_fishhook_eligible_hg38()
  }

  if (!inherits(eligible_gr, "GRanges")) {
    stop("`eligible_gr` must be a GRanges object.")
  }

  seqlevelsStyle(event_gr) <- "UCSC"
  seqlevelsStyle(hypothesis_gr) <- "UCSC"
  seqlevelsStyle(eligible_gr) <- "UCSC"

  common_hyp_seqlevels <- intersect(seqlevels(eligible_gr), seqlevels(hypothesis_gr))
  common_evt_seqlevels <- intersect(seqlevels(eligible_gr), seqlevels(event_gr))

  seqlevels(hypothesis_gr) <- common_hyp_seqlevels
  seqinfo(hypothesis_gr) <- seqinfo(eligible_gr)[common_hyp_seqlevels]
  seqlevels(event_gr) <- common_evt_seqlevels
  seqinfo(event_gr) <- seqinfo(eligible_gr)[common_evt_seqlevels]

  genome(hypothesis_gr) <- "hg38"
  genome(event_gr) <- "hg38"
  genome(eligible_gr) <- "hg38"

  if (is.null(fishhook_covariates) && !is.null(fishhook_covariate_data)) {
    cov_bundle <- make_fishhook_covariates_from_data(
      covariate_data = fishhook_covariate_data,
      reg_gr = hypothesis_gr
    )
    hypothesis_gr <- cov_bundle$hypotheses
    covs <- cov_bundle$covariates
  } else {
    covs <- instantiate_fishhook_covariates(fishhook_covariates)
  }

  fish <- fishHook::Fish(
    hypotheses = hypothesis_gr,
    events = event_gr,
    eligible = eligible_gr,
    covariates = covs,
    idcol = idcol
  )

  fish$score(...)

  if (is.null(fish$res)) {
    stop("fishHook did not produce a `$res` result table.")
  }

  res_dt <- as.data.table(fish$res)
  if (!"reg_elem_id" %in% names(res_dt) && "reg_elem_id" %in% names(mcols(hypothesis_gr))) {
    res_dt[, reg_elem_id := as.character(mcols(hypothesis_gr)$reg_elem_id)]
  }

  extract_fishhook_reg_scores(res_dt)
}

run_somatic_scoring_pipeline <- function(
  maf,
  reg_ref_path = NULL,
  dndscv_result = NULL,
  fishhook_result = NULL,
  refdb = NULL,
  cv = NULL,
  max_muts_per_gene_per_sample = 6L,
  max_coding_muts_per_sample = 5000L,
  eligible_gr = NULL,
  fishhook_covariates = NULL,
  fishhook_covariate_data = NULL,
  idcol = "Tumor_Sample_Barcode",
  ...
) {
  gene_scores <- if (!is.null(dndscv_result)) {
    extract_dndscv_gene_scores(dndscv_result)
  } else {
    run_dndscv_gene_scoring(
      maf = maf,
      refdb = refdb,
      cv = cv,
      max_muts_per_gene_per_sample = max_muts_per_gene_per_sample,
      max_coding_muts_per_sample = max_coding_muts_per_sample,
      ...
    )
  }

  reg_scores <- NULL
  if (!is.null(fishhook_result)) {
    reg_scores <- extract_fishhook_reg_scores(fishhook_result)
  } else if (!is.null(reg_ref_path)) {
    reg_scores <- run_fishhook_reg_scoring(
      maf = maf,
      reg_ref_path = reg_ref_path,
      eligible_gr = eligible_gr,
      fishhook_covariates = fishhook_covariates,
      fishhook_covariate_data = fishhook_covariate_data,
      idcol = idcol,
      ...
    )
  }

  list(
    gene_scores = gene_scores,
    reg_scores = reg_scores
  )
}
