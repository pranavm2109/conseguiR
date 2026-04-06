#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(devtools)
  library(data.table)
})

load_all(".")

# Manual end-to-end verbose rerun script.
# Run from the repo root in an interactive R session:
# source("scripts/Testing/R/99_manual_verbose_full_rerun.R")

repo_root <- getwd()

relpath <- function(...) {
  file.path(...)
}

cfg <- list(
  clean_outputs = TRUE,
  target_genes = 50L,
  lymphoma_driver_labels = c(
    "MYC", "BCL2", "BCL6", "TP53", "EZH2",
    "KMT2D", "CREBBP", "EP300", "CARD11", "CD79B",
    "PIM1", "BTG2", "SOCS1", "XPO1", "IRF4"
  ),
  inputs = list(
    gwas_sumstats = relpath("data/raw/Testing/34737426-GCST90043906-EFO_0000403.h.tsv"),
    somatic_maf = relpath("data/raw/Testing/2026-01-09_no_CLL_lymph_only_pcawg_maf_tcga_order_hg38.maf"),
    reg_ref_path = relpath("data/raw/Testing/2026-01-26_UCSC_all_unfiltered_reg_elements.loc"),
    gene_loc_path = relpath("data/raw/NCBI38/NCBI38.gene.loc"),
    reg_loc_path = relpath("data/raw/GeneHancer/2026-01-26_UCSC_all_unfiltered_reg_elements.loc"),
    reference_bfile = relpath("data/raw/g1000_eur/g1000_eur"),
    dndscv_refdb = relpath("data/raw/Testing/RefCDS_human_GRCh38.p12.rda"),
    epigenomic_track_dir = relpath("data/raw/Testing"),
    fishhook_covariate_rds = relpath(
      "data/raw/Testing/2026-01-26_all_reg_elems_sample_level_mut_frac_comparison_bet_only_memory_b_normal_and_non_cll_malig_b_cells.rds"
    )
  ),
  output = list(
    score_dir = relpath("data/processed/scores"),
    graph_dir = relpath("data/processed/graphs"),
    diffusion_dir = relpath("data/processed/diffusion"),
    subgraph_dir = relpath("data/processed/subgraphs"),
    figure_dir = relpath("data/processed/figures")
  )
)

ensure_dir <- function(path) {
  dir.create(path, recursive = TRUE, showWarnings = FALSE)
  invisible(path)
}

list_generated_files <- function(path) {
  files <- list.files(path, all.files = FALSE, full.names = TRUE, recursive = TRUE)
  files <- files[!grepl("/\\._", files, fixed = FALSE)]
  sort(files[file.exists(files)])
}

clean_dir <- function(path) {
  if (!dir.exists(path)) {
    dir.create(path, recursive = TRUE, showWarnings = FALSE)
    return(invisible(path))
  }

  old <- list.files(path, all.files = TRUE, full.names = TRUE, no.. = TRUE)
  if (length(old) > 0L) {
    unlink(old, recursive = TRUE, force = TRUE)
  }

  invisible(path)
}

message_block <- function(...) {
  cat("\n", paste(rep("=", 80), collapse = ""), "\n", sep = "")
  cat(..., "\n", sep = "")
  cat(paste(rep("=", 80), collapse = ""), "\n\n", sep = "")
}

prepare_manual_magma_gwas_file <- function(raw_gwas_path, output_path) {
  message("Reading GWAS summary statistics for MAGMA preprocessing...")
  dt <- data.table::fread(
    raw_gwas_path,
    select = c("hm_rsid", "hm_variant_id", "hm_chrom", "hm_pos", "p_value"),
    showProgress = TRUE
  )

  before_n <- nrow(dt)
  dt <- dt[
    !is.na(hm_chrom) & trimws(as.character(hm_chrom)) != "" &
      !is.na(hm_pos) &
      !is.na(p_value) &
      (
        (!is.na(hm_rsid) & trimws(hm_rsid) != "") |
          (!is.na(hm_variant_id) & trimws(hm_variant_id) != "")
      )
  ]
  dropped_n <- before_n - nrow(dt)
  if (dropped_n > 0L) {
    message("Dropped ", dropped_n, " malformed GWAS row(s) before MAGMA preprocessing.")
  }

  # Prefer harmonized rsIDs when available; otherwise fall back to harmonized variant IDs.
  dt[, magma_snp_id := data.table::fifelse(
    !is.na(hm_rsid) & trimws(hm_rsid) != "",
    hm_rsid,
    hm_variant_id
  )]
  dt <- dt[!is.na(magma_snp_id) & trimws(magma_snp_id) != ""]

  message("Writing minimized MAGMA GWAS file...")
  data.table::fwrite(dt, output_path, sep = "\t")

  if (!file.exists(output_path) || file.info(output_path)$size <= 0) {
    stop("Failed to write minimized MAGMA GWAS file: ", output_path)
  }

  output_path
}

prepare_manual_magma_cache <- function(minimized_gwas_path, cache_prefix) {
  message("Reading minimized GWAS file for MAGMA cache generation...")
  dt <- data.table::fread(minimized_gwas_path, showProgress = TRUE)

  required_cols <- c("hm_variant_id", "hm_chrom", "hm_pos", "p_value")
  missing_cols <- setdiff(required_cols, names(dt))
  if (length(missing_cols) > 0L) {
    stop("Minimized GWAS file is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  before_n <- nrow(dt)
  dt <- dt[
    !is.na(hm_chrom) & trimws(as.character(hm_chrom)) != "" &
      !is.na(hm_pos) &
      !is.na(p_value) &
      (
        (!is.na(hm_rsid) & trimws(hm_rsid) != "") |
          (!is.na(hm_variant_id) & trimws(hm_variant_id) != "")
      )
  ]
  dropped_n <- before_n - nrow(dt)
  if (dropped_n > 0L) {
    message("Dropped ", dropped_n, " malformed minimized GWAS row(s) before MAGMA cache generation.")
  }

  dt[, magma_snp_id := data.table::fifelse(
    !is.na(hm_rsid) & trimws(hm_rsid) != "",
    hm_rsid,
    hm_variant_id
  )]
  dt <- dt[!is.na(magma_snp_id) & trimws(magma_snp_id) != ""]

  snp_loc_path <- paste0(cache_prefix, ".snp_loc.tsv")
  pval_path <- paste0(cache_prefix, ".pval.tsv")

  message("Building MAGMA SNP-location table...")
  snp_loc <- unique(dt[, .(
    SNP = magma_snp_id,
    CHR = hm_chrom,
    POS = hm_pos
  )])

  message("Building MAGMA p-value table...")
  pval <- unique(dt[, .(
    SNP = magma_snp_id,
    P = p_value
  )])

  message("Writing MAGMA cache files...")
  data.table::fwrite(snp_loc, snp_loc_path, sep = "\t", col.names = FALSE)
  data.table::fwrite(pval, pval_path, sep = "\t", col.names = TRUE)

  if (!file.exists(snp_loc_path) || file.info(snp_loc_path)$size <= 0) {
    stop("Failed to write MAGMA snp-loc cache file: ", snp_loc_path)
  }
  if (!file.exists(pval_path) || file.info(pval_path)$size <= 0) {
    stop("Failed to write MAGMA p-value cache file: ", pval_path)
  }

  invisible(list(
    snp_loc_path = snp_loc_path,
    pval_path = pval_path
  ))
}

stop_if_missing <- function(paths) {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0L) {
    stop("Missing required input(s):\n", paste(missing, collapse = "\n"))
  }
}

input_files_to_check <- c(
  cfg$inputs$gwas_sumstats,
  cfg$inputs$somatic_maf,
  cfg$inputs$reg_ref_path,
  cfg$inputs$gene_loc_path,
  cfg$inputs$reg_loc_path,
  cfg$inputs$dndscv_refdb,
  cfg$inputs$fishhook_covariate_rds
)

stop_if_missing(input_files_to_check)

reference_bfile_required <- paste0(cfg$inputs$reference_bfile, c(".bed", ".bim", ".fam"))
stop_if_missing(reference_bfile_required)

if (!dir.exists(cfg$inputs$epigenomic_track_dir)) {
  stop("Epigenomic track directory does not exist:\n", cfg$inputs$epigenomic_track_dir)
}

if (isTRUE(cfg$clean_outputs)) {
  message_block("Cleaning score/graph/diffusion/subgraph/figure outputs")
  clean_dir(cfg$output$score_dir)
  clean_dir(cfg$output$graph_dir)
  clean_dir(cfg$output$diffusion_dir)
  clean_dir(cfg$output$subgraph_dir)
  clean_dir(cfg$output$figure_dir)
} else {
  lapply(cfg$output, ensure_dir)
}

message_block("Initializing backend graphs")
initialize_backend_graphs(verbose = TRUE)

fishhook_covariate_data <- readRDS(cfg$inputs$fishhook_covariate_rds)

message_block("Preparing minimized GWAS file for MAGMA")
gwas_magma_input <- prepare_manual_magma_gwas_file(
  raw_gwas_path = cfg$inputs$gwas_sumstats,
  output_path = file.path(cfg$output$score_dir, "gwas_magma_minimal.tsv")
)

message_block("Preparing shared MAGMA GWAS cache")
prepare_manual_magma_cache(
  minimized_gwas_path = gwas_magma_input,
  cache_prefix = file.path(cfg$output$score_dir, "magma_shared_gwas_cache")
)

message_block("Running germline gene scoring")
germline_gene <- run_germline_gene_scoring(
  gwas_sumstats = gwas_magma_input,
  gene_loc_path = cfg$inputs$gene_loc_path,
  reference_bfile = cfg$inputs$reference_bfile,
  output_prefix = file.path(cfg$output$score_dir, "germline_gene_scores"),
  sample_size = 456348L,
  magma_gwas_cache_prefix = file.path(cfg$output$score_dir, "magma_shared_gwas_cache"),
  reuse_existing_gwas_cache = TRUE,
  verbose = TRUE
)

message_block("Running germline regulatory scoring")
germline_reg <- run_germline_regulatory_scoring(
  gwas_sumstats = gwas_magma_input,
  reg_loc_path = cfg$inputs$reg_loc_path,
  reference_bfile = cfg$inputs$reference_bfile,
  output_prefix = file.path(cfg$output$score_dir, "germline_reg_scores"),
  sample_size = 456348L,
  magma_gwas_cache_prefix = file.path(cfg$output$score_dir, "magma_shared_gwas_cache"),
  reuse_existing_gwas_cache = TRUE,
  verbose = TRUE
)

germline_scores <- list(
  gene_scores = germline_gene$gene_scores,
  reg_scores = germline_reg$reg_scores
)

message_block("Running somatic scoring")
somatic_scores <- prepare_somatic_scores(
  maf = cfg$inputs$somatic_maf,
  refdb = cfg$inputs$dndscv_refdb,
  reg_ref_path = cfg$inputs$reg_ref_path,
  gene_output_path = file.path(cfg$output$score_dir, "somatic_gene_scores.tsv"),
  reg_output_path = file.path(cfg$output$score_dir, "somatic_reg_scores.tsv"),
  fishhook_covariate_data = fishhook_covariate_data,
  verbose = TRUE
)

message_block("Running epigenomic scoring")
epigenomic_scores <- prepare_epigenomic_scores(
  reg_ref_path = cfg$inputs$reg_ref_path,
  track_dir = cfg$inputs$epigenomic_track_dir,
  output_path = file.path(cfg$output$score_dir, "epigenomic_reg_scores.tsv"),
  exclude_patterns = c("_BL_", "_FL_", "broken_signal_track"),
  verbose = TRUE
)

message_block("Building scored gene-reg graph")
scored_graph <- build_scored_gene_reg_graph(
  graph_rds_path = relpath("data/processed/gene_reg_graph_no_scores.rds"),
  germline_scores = germline_scores,
  somatic_scores = somatic_scores,
  epigenomic_scores = epigenomic_scores,
  output_prefix = file.path(cfg$output$graph_dir, "gene_reg_graph_scored"),
  verbose = TRUE
)

message("Scored graph node counts:")
print(scored_graph$nodes[, .N, by = node_type][order(node_type)])

message_block("Running gene-reg diffusion")
diffusion <- run_gene_reg_diffusion(
  scored_graph = scored_graph,
  output_dir = cfg$output$diffusion_dir,
  output_stem = "gene_reg_graph_diffusion",
  verbose = TRUE
)

message_block("Calling selected subgraph")
selected_subgraph <- call_selected_subgraph(
  diffusion = diffusion,
  gg_nodes_path = relpath("data/processed/gene_gene_graph_nodes.tsv.gz"),
  gg_edges_path = relpath("data/processed/gene_gene_graph_edges.tsv.gz"),
  output_dir = cfg$output$subgraph_dir,
  output_stem = "gene_gene_selected_subgraph",
  target_genes = cfg$target_genes,
  verbose = TRUE
)

message_block("Plotting pre- and post-diffusion figures")
plot_germline_gene_scores(
  germline_scores = germline_scores,
  stage = "pre",
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "germline_genes_pre_rank.pdf"),
  verbose = TRUE
)

plot_germline_gene_scores(
  germline_scores = germline_scores,
  diffusion = diffusion,
  stage = "post",
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "germline_genes_post_rank.pdf"),
  verbose = TRUE
)

plot_germline_reg_scores(
  germline_scores = germline_scores,
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "germline_regs_rank.pdf"),
  verbose = TRUE
)

plot_somatic_gene_scores(
  somatic_scores = somatic_scores,
  stage = "pre",
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "somatic_genes_pre_volcano.pdf"),
  verbose = TRUE
)

plot_somatic_gene_scores(
  somatic_scores = somatic_scores,
  diffusion = diffusion,
  stage = "post",
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "somatic_genes_post_volcano.pdf"),
  verbose = TRUE
)

plot_somatic_reg_scores(
  somatic_scores = somatic_scores,
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "somatic_regs_volcano.pdf"),
  verbose = TRUE
)

plot_epigenomic_reg_scores(
  epigenomic_scores = epigenomic_scores,
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "epigenomic_regs.pdf"),
  verbose = TRUE
)

plot_epigenomic_gene_scores(
  diffusion = diffusion,
  label_features = cfg$lymphoma_driver_labels,
  plot_file_path = file.path(cfg$output$figure_dir, "epigenomic_genes_post.pdf"),
  verbose = TRUE
)

plot_selected_subgraph(
  selected_subgraph = selected_subgraph,
  bundle_output_prefix = file.path(cfg$output$subgraph_dir, "gene_gene_selected_subgraph_plot_bundle"),
  plot_file_path = file.path(cfg$output$figure_dir, "gene_gene_selected_subgraph_plot.pdf"),
  verbose = TRUE
)

message_block("Run complete")
message("Generated files:")
print(list_generated_files(cfg$output$score_dir))
print(list_generated_files(cfg$output$graph_dir))
print(list_generated_files(cfg$output$diffusion_dir))
print(list_generated_files(cfg$output$subgraph_dir))
print(list_generated_files(cfg$output$figure_dir))
