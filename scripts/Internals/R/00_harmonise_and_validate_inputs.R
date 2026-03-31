#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(GenomicRanges)
  library(IRanges)
  library(GenomeInfoDb)
  library(rtracklayer)
})

# Minimal harmonization helpers for early conseguiR scripts.
#
# GWAS output columns:
# - variant_id: taken from rs_id, rsid, hm_rsid, or variant_id
# - chromosome: taken from chromosome, hm_chrom, chr, or CHR
# - base_pair_location: taken from base_pair_location, hm_pos, bp, pos, or BP
# - p_value: taken from p_value, p, P, or pval
#
# Somatic MAF output columns:
# - sample_id
# - chromosome
# - start_position
# - end_position
# - ref
# - alt

ensure_data_table <- function(x) {
  if (!is.data.frame(x)) {
    stop("Expected a data.frame or data.table.")
  }

  as.data.table(x)
}

normalize_colnames <- function(x) {
  dt <- copy(ensure_data_table(x))
  clean_names <- names(dt)
  clean_names <- gsub("\\.", "_", clean_names)
  clean_names <- gsub("\\s+", "_", clean_names)
  clean_names <- gsub("-", "_", clean_names, fixed = TRUE)
  setnames(dt, clean_names)
  dt
}

pick_first_existing_column <- function(dt, candidates, field_name) {
  hit <- intersect(candidates, names(dt))

  if (length(hit) == 0) {
    stop(
      "Could not find a column for ",
      field_name,
      ". Tried: ",
      paste(candidates, collapse = ", ")
    )
  }

  hit[[1]]
}

validate_gwas_sumstats <- function(sumstats) {
  dt <- normalize_colnames(sumstats)

  id_col <- pick_first_existing_column(
    dt,
    c("rs_id", "rsid", "hm_rsid", "variant_id"),
    "variant identifier"
  )
  chr_col <- pick_first_existing_column(
    dt,
    c("chromosome", "hm_chrom", "chr", "CHR"),
    "chromosome"
  )
  pos_col <- pick_first_existing_column(
    dt,
    c("base_pair_location", "hm_pos", "bp", "pos", "BP"),
    "base-pair position"
  )
  p_col <- pick_first_existing_column(
    dt,
    c("p_value", "p", "P", "pval"),
    "p-value"
  )

  out <- unique(dt[, .(
    variant_id = as.character(get(id_col)),
    chromosome = as.character(get(chr_col)),
    base_pair_location = as.integer(get(pos_col)),
    p_value = as.numeric(get(p_col))
  )])

  out[, variant_id := trimws(variant_id)]
  out[, chromosome := trimws(chromosome)]
  out <- out[
    !is.na(variant_id) & variant_id != "" &
    !is.na(chromosome) & chromosome != "" &
    !is.na(base_pair_location) &
    !is.na(p_value)
  ]

  if (nrow(out) == 0L) {
    stop("GWAS summary statistics have no usable rows after validation.")
  }

  attr(out, "conseguiR_input_type") <- "gwas_sumstats_magma_minimal"
  out
}

prepare_magma_input <- function(sumstats) {
  dt <- validate_gwas_sumstats(sumstats)

  snp_loc <- unique(dt[, .(
    SNP = variant_id,
    CHR = chromosome,
    POS = base_pair_location
  )])

  pval <- unique(dt[, .(
    SNP = variant_id,
    P = p_value
  )])

  attr(snp_loc, "conseguiR_input_type") <- "magma_snp_loc"
  attr(pval, "conseguiR_input_type") <- "magma_pval"

  list(
    snp_loc = snp_loc,
    pval = pval
  )
}

write_magma_input_files <- function(sumstats, snp_loc_path, pval_path) {
  magma_input <- prepare_magma_input(sumstats)

  fwrite(magma_input$snp_loc, snp_loc_path, sep = "\t", col.names = FALSE)
  fwrite(magma_input$pval, pval_path, sep = "\t", col.names = TRUE)

  invisible(list(
    snp_loc_path = snp_loc_path,
    pval_path = pval_path
  ))
}

validate_somatic_maf <- function(maf) {
  dt <- normalize_colnames(maf)

  sample_col <- pick_first_existing_column(
    dt,
    c("Tumor_Sample_Barcode", "tumor_sample_barcode", "sampleID", "sample_id", "idcol"),
    "sample identifier"
  )
  chr_col <- pick_first_existing_column(
    dt,
    c("Chromosome", "chromosome", "chr", "CHR"),
    "chromosome"
  )
  start_col <- pick_first_existing_column(
    dt,
    c("Start_Position", "Start_position", "start_position", "Start", "start"),
    "start position"
  )
  end_col <- pick_first_existing_column(
    dt,
    c("End_Position", "End_position", "end_position", "End", "end"),
    "end position"
  )
  ref_col <- pick_first_existing_column(
    dt,
    c("Reference_Allele", "reference_allele", "ref"),
    "reference allele"
  )
  alt_col <- pick_first_existing_column(
    dt,
    c("Tumor_Seq_Allele2", "tumor_seq_allele2", "mut", "alt"),
    "alternate allele"
  )

  out <- unique(dt[, .(
    sample_id = as.character(get(sample_col)),
    chromosome = as.character(get(chr_col)),
    start_position = as.integer(get(start_col)),
    end_position = as.integer(get(end_col)),
    ref = as.character(get(ref_col)),
    alt = as.character(get(alt_col))
  )])

  attr(out, "conseguiR_input_type") <- "somatic_maf_minimal"
  out
}

prepare_dndscv_input <- function(maf) {
  dt <- validate_somatic_maf(maf)

  out <- dt[, .(
    sampleID = sample_id,
    chr = chromosome,
    pos = start_position,
    ref = ref,
    mut = alt
  )]

  attr(out, "conseguiR_input_type") <- "dndscv_input"
  out
}

prepare_fishhook_input <- function(maf) {
  dt <- validate_somatic_maf(maf)

  out <- dt[, .(
    Tumor_Sample_Barcode = sample_id,
    Chromosome = chromosome,
    Start_Position = start_position,
    End_Position = end_position,
    Reference_Allele = ref,
    Tumor_Seq_Allele2 = alt
  )]

  attr(out, "conseguiR_input_type") <- "fishhook_input"
  out
}

validate_regulatory_element_reference <- function(reg_ref_path) {
  if (!file.exists(reg_ref_path)) {
    stop("Regulatory element reference file does not exist: ", reg_ref_path)
  }

  reg_dt <- fread(reg_ref_path, header = FALSE)

  if (ncol(reg_dt) < 4) {
    stop("Regulatory element reference must have at least 4 columns: reg_elem_id, chrom, start, end.")
  }

  base_names <- c("reg_elem_id", "chrom", "start_raw", "end_raw", "strand", "reg_elem_name")
  setnames(reg_dt, names(reg_dt)[seq_len(min(length(base_names), ncol(reg_dt)))], base_names[seq_len(min(length(base_names), ncol(reg_dt)))])

  reg_dt[, seqname := ifelse(grepl("^chr", chrom), chrom, paste0("chr", chrom))]
  reg_dt[, start := as.integer(start_raw)]
  reg_dt[, end := as.integer(end_raw)]

  has_zero_start <- any(reg_dt$start == 0L, na.rm = TRUE)
  if (has_zero_start) {
    reg_dt[, start := start + 1L]
  }

  reg_dt <- reg_dt[
    !is.na(reg_elem_id) & reg_elem_id != "" &
    !is.na(seqname) & seqname != "" &
    !is.na(start) & !is.na(end) &
    start > 0L & end >= start
  ]

  if (nrow(reg_dt) == 0) {
    stop("Regulatory element reference has no usable rows after cleaning.")
  }

  reg_gr <- GRanges(
    seqnames = reg_dt$seqname,
    ranges = IRanges(start = reg_dt$start, end = reg_dt$end),
    strand = "*"
  )

  mcols(reg_gr)$reg_elem_id <- as.character(reg_dt$reg_elem_id)
  if ("reg_elem_name" %in% names(reg_dt)) {
    mcols(reg_gr)$reg_elem_name <- as.character(reg_dt$reg_elem_name)
  }

  attr(reg_gr, "conseguiR_input_type") <- "regulatory_element_reference"
  reg_gr
}

validate_epigenomic_bigwigs <- function(bw_files, reg_gr, exclude_patterns = c("_BL_", "_FL_")) {
  if (length(bw_files) == 1L && dir.exists(bw_files)) {
    bw_files <- list.files(bw_files, pattern = "\\.(bw|bigWig|bigwig)$", full.names = TRUE)
  }

  if (length(exclude_patterns) > 0) {
    exclude_regex <- paste(exclude_patterns, collapse = "|")
    bw_files <- bw_files[!grepl(exclude_regex, basename(bw_files))]
  }

  if (length(bw_files) < 3) {
    stop("Need at least three bigWig files after filtering.")
  }

  missing_files <- bw_files[!file.exists(bw_files)]
  if (length(missing_files) > 0) {
    stop("These bigWig files do not exist: ", paste(missing_files, collapse = ", "))
  }

  reg_gr_test <- reg_gr
  file_summaries <- vector("list", length(bw_files))

  for (i in seq_along(bw_files)) {
    bw_path <- bw_files[[i]]

    bw_info <- tryCatch(
      seqinfo(BigWigFile(bw_path)),
      error = function(e) {
        stop("Could not read bigWig header for ", bw_path, ": ", conditionMessage(e))
      }
    )

    seqlevelsStyle(reg_gr_test) <- "UCSC"
    common_seqlevels <- intersect(seqlevels(reg_gr_test), seqlevels(bw_info))

    if (length(common_seqlevels) == 0) {
      stop("No overlapping seqlevels between regulatory elements and bigWig file: ", bw_path)
    }

    reg_subset <- keepSeqlevels(reg_gr_test, common_seqlevels, pruning.mode = "coarse")
    if (length(reg_subset) == 0) {
      stop("No regulatory elements remain after seqlevel harmonization for bigWig file: ", bw_path)
    }

    test_gr <- reg_subset[seq_len(min(10L, length(reg_subset)))]

    imported <- tryCatch(
      import(bw_path, which = test_gr, as = "NumericList"),
      error = function(e) {
        stop("Failed to import signal from bigWig file ", bw_path, ": ", conditionMessage(e))
      }
    )

    if (length(imported) != length(test_gr)) {
      stop("Imported signal length does not match queried regulatory elements for bigWig file: ", bw_path)
    }

    file_summaries[[i]] <- data.table(
      file = bw_path,
      n_common_seqlevels = length(common_seqlevels),
      n_test_intervals = length(test_gr)
    )
  }

  out <- rbindlist(file_summaries)
  attr(out, "conseguiR_input_type") <- "epigenomic_bigwig_validation"
  out
}

validate_epigenomic_inputs <- function(bw_files, reg_ref_path, exclude_patterns = c("_BL_", "_FL_")) {
  reg_gr <- validate_regulatory_element_reference(reg_ref_path)
  bw_summary <- validate_epigenomic_bigwigs(
    bw_files = bw_files,
    reg_gr = reg_gr,
    exclude_patterns = exclude_patterns
  )

  list(
    reg_gr = reg_gr,
    bigwig_summary = bw_summary
  )
}
