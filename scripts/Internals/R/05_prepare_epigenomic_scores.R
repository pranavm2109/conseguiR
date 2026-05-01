#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
  library(matrixStats)
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

# Epigenomic scoring design:
# - regulatory-element scores come from cross-track variation in bigWig signal
# - signal is quantified per regulatory element in each track
# - final output should be a simple score table with reg_elem_id and zstat

make_extended_mhc_gr <- function(
  seqname = "chr6",
  start = 25726063L,
  end = 33400644L
) {
  GRanges(
    seqnames = seqname,
    ranges = IRanges(start = start, end = end),
    strand = "*"
  )
}

drop_extended_mhc_overlaps <- function(gr, mhc_gr = make_extended_mhc_gr()) {
  if (length(gr) == 0L) {
    return(gr)
  }

  mhc_gr <- keepSeqlevels(mhc_gr, intersect(seqlevels(mhc_gr), seqlevels(gr)), pruning.mode = "coarse")
  if (length(mhc_gr) == 0L) {
    return(gr)
  }

  seqlevelsStyle(mhc_gr) <- seqlevelsStyle(gr)[[1]]
  gr[!overlapsAny(gr, mhc_gr, ignore.strand = TRUE)]
}

zscore_vec <- function(x) {
  x <- as.numeric(x)
  mu <- mean(x, na.rm = TRUE)
  sig <- stats::sd(x, na.rm = TRUE)

  if (!is.finite(sig) || sig == 0) {
    return(rep(0, length(x)))
  }

  (x - mu) / sig
}

deep_copy_object <- function(x) {
  unserialize(serialize(x, NULL))
}

list_epigenomic_track_files <- function(
  track_dir,
  pattern = "\\.bw$",
  exclude_patterns = NULL,
  min_tracks = 3L
) {
  if (!dir.exists(track_dir)) {
    stop("Epigenomic track directory does not exist: ", track_dir)
  }

  bw_files <- list.files(track_dir, pattern = pattern, full.names = TRUE)

  if (length(exclude_patterns) > 0L) {
    exclude_regex <- paste(exclude_patterns, collapse = "|")
    bw_files <- bw_files[!grepl(exclude_regex, basename(bw_files))]
  }

  if (length(bw_files) == 0L) {
    stop("No epigenomic bigWig files were found after filtering in: ", track_dir)
  }

  bw_files <- validate_epigenomic_tracks(bw_files)

  if (length(bw_files) < min_tracks) {
    stop("At least ", min_tracks, " epigenomic bigWig tracks are required; found ", length(bw_files), ".")
  }

  bw_files
}

load_regulatory_elements_for_epigenomic_scores <- function(
  reg_ref_path,
  drop_mhc = TRUE
) {
  reg_gr <- validate_regulatory_element_reference(reg_ref_path)

  if (!"reg_elem_id" %in% names(mcols(reg_gr))) {
    stop("Regulatory element reference must contain `reg_elem_id` metadata.")
  }

  seqlevelsStyle(reg_gr) <- "UCSC"

  if (isTRUE(drop_mhc)) {
    reg_gr <- drop_extended_mhc_overlaps(reg_gr)
  }

  reg_gr
}

harmonize_regulatory_seqlevels_to_track <- function(reg_gr, bw_file) {
  bw_seqinfo <- seqinfo(BigWigFile(bw_file))

  seqlevelsStyle(reg_gr) <- "UCSC"
  common_seqlevels <- intersect(seqlevels(reg_gr), seqlevels(bw_seqinfo))

  if (length(common_seqlevels) == 0L) {
    stop("No shared seqlevels between regulatory elements and bigWig file: ", basename(bw_file))
  }

  reg_gr <- keepSeqlevels(reg_gr, common_seqlevels, pruning.mode = "coarse")
  seqinfo(reg_gr) <- bw_seqinfo[common_seqlevels]
  reg_gr
}

extract_reg_signal_from_bigwig <- function(bw_file, reg_gr, summary_fun = mean) {
  signal_list <- tryCatch(
    {
      import(bw_file, which = reg_gr, as = "NumericList")
    },
    error = function(e) {
      stop("Failed to import bigWig file `", bw_file, "`: ", conditionMessage(e))
    }
  )

  signal_vec <- vapply(
    signal_list,
    FUN.VALUE = numeric(1),
    FUN = function(x) {
      if (length(x) == 0L) {
        return(0)
      }

      x <- as.numeric(x)
      x <- x[is.finite(x)]

      if (length(x) == 0L) {
        return(0)
      }

      summary_fun(x)
    }
  )

  if (length(signal_vec) != length(reg_gr)) {
    stop("Signal vector length does not match regulatory-element count for: ", basename(bw_file))
  }

  signal_vec
}

quantify_epigenomic_signal_matrix <- function(
  bw_files,
  reg_gr,
  summary_fun = mean
) {
  signal_mat <- vapply(
    bw_files,
    FUN.VALUE = numeric(length(reg_gr)),
    FUN = function(bw_file) {
      extract_reg_signal_from_bigwig(
        bw_file = bw_file,
        reg_gr = reg_gr,
        summary_fun = summary_fun
      )
    }
  )

  signal_mat <- as.matrix(signal_mat)
  colnames(signal_mat) <- basename(bw_files)
  signal_mat
}

compute_epigenomic_reg_scores <- function(
  signal_mat,
  reg_gr,
  transform = c("log1p", "none")
) {
  transform <- match.arg(transform)

  if (!is.matrix(signal_mat)) {
    signal_mat <- as.matrix(signal_mat)
  }

  score_mat <- switch(
    transform,
    log1p = log1p(signal_mat),
    none = signal_mat
  )

  sd_signal <- rowSds(score_mat, na.rm = TRUE)
  zstat <- zscore_vec(sd_signal)

  data.table(
    reg_elem_id = as.character(mcols(reg_gr)$reg_elem_id),
    zstat = zstat
  )[!is.na(reg_elem_id) & reg_elem_id != ""]
}

compute_epigenomic_reg_diagnostics <- function(
  signal_mat,
  reg_gr,
  transform = c("log1p", "none")
) {
  transform <- match.arg(transform)

  if (!is.matrix(signal_mat)) {
    signal_mat <- as.matrix(signal_mat)
  }

  score_mat <- switch(
    transform,
    log1p = log1p(signal_mat),
    none = signal_mat
  )

  mean_signal <- rowMeans(score_mat, na.rm = TRUE)
  sd_signal <- rowSds(score_mat, na.rm = TRUE)
  var_signal <- rowVars(score_mat, na.rm = TRUE)

  out <- data.table(
    reg_elem_id = as.character(mcols(reg_gr)$reg_elem_id),
    mean_signal_across_tracks = mean_signal,
    sd_signal_across_tracks = sd_signal,
    var_signal_across_tracks = var_signal,
    zstat = zscore_vec(sd_signal)
  )

  if ("reg_elem_name" %in% names(mcols(reg_gr))) {
    out[, reg_elem_name := as.character(mcols(reg_gr)$reg_elem_name)]
    setcolorder(
      out,
      c(
        "reg_elem_id",
        "reg_elem_name",
        "mean_signal_across_tracks",
        "sd_signal_across_tracks",
        "var_signal_across_tracks",
        "zstat"
      )
    )
  }

  setorder(out, -zstat)
  out
}

run_epigenomic_reg_scoring <- function(
  track_dir = NULL,
  reg_ref_path,
  bw_files = NULL,
  exclude_patterns = NULL,
  min_tracks = 3L,
  drop_mhc = TRUE,
  transform = c("log1p", "none"),
  return_diagnostics = FALSE,
  summary_fun = mean,
  verbose = FALSE
) {
  transform <- match.arg(transform)
  pb <- if (isTRUE(verbose)) utils::txtProgressBar(min = 0, max = 5, style = 3) else NULL
  on.exit(if (!is.null(pb)) close(pb), add = TRUE)

  if (is.null(bw_files)) {
    conseguiR_verbose_message(verbose, "Discovering epigenomic tracks...")
    bw_files <- list_epigenomic_track_files(
      track_dir = track_dir,
      exclude_patterns = exclude_patterns,
      min_tracks = min_tracks
    )
  } else {
    conseguiR_verbose_message(verbose, "Validating supplied epigenomic tracks...")
    bw_files <- validate_epigenomic_tracks(bw_files)
    if (length(bw_files) < min_tracks) {
      stop("At least ", min_tracks, " epigenomic bigWig tracks are required; found ", length(bw_files), ".")
    }
  }
  if (!is.null(pb)) utils::setTxtProgressBar(pb, 1)

  conseguiR_verbose_message(verbose, "Loading regulatory elements for epigenomic scoring...")
  reg_gr <- load_regulatory_elements_for_epigenomic_scores(
    reg_ref_path = reg_ref_path,
    drop_mhc = drop_mhc
  )
  if (!is.null(pb)) utils::setTxtProgressBar(pb, 2)

  conseguiR_verbose_message(verbose, "Harmonizing regulatory elements to track seqlevels...")
  reg_gr <- harmonize_regulatory_seqlevels_to_track(
    reg_gr = reg_gr,
    bw_file = bw_files[[1]]
  )
  if (!is.null(pb)) utils::setTxtProgressBar(pb, 3)

  conseguiR_verbose_message(verbose, "Quantifying bigWig signal over regulatory elements...")
  signal_mat <- quantify_epigenomic_signal_matrix(
    bw_files = bw_files,
    reg_gr = reg_gr,
    summary_fun = summary_fun
  )
  if (!is.null(pb)) utils::setTxtProgressBar(pb, 4)

  conseguiR_verbose_message(verbose, "Computing epigenomic z-scores...")
  zscores <- compute_epigenomic_reg_scores(
    signal_mat = signal_mat,
    reg_gr = reg_gr,
    transform = transform
  )
  if (!is.null(pb)) utils::setTxtProgressBar(pb, 5)

  if (!isTRUE(return_diagnostics)) {
    return(zscores)
  }

  diagnostics <- compute_epigenomic_reg_diagnostics(
    signal_mat = signal_mat,
    reg_gr = reg_gr,
    transform = transform
  )

  list(
    zscores = zscores,
    diagnostics = diagnostics,
    signal_matrix = deep_copy_object(signal_mat),
    reg_gr = deep_copy_object(reg_gr),
    bw_files = bw_files
  )
}
