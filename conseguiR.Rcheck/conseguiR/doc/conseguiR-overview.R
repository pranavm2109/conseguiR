## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  message = FALSE,
  warning = FALSE
)

## ----libraries----------------------------------------------------------------
library(conseguiR)
library(data.table)
library(igraph)

## ----validation-inputs--------------------------------------------------------
toy_gwas <- data.table(
  variant_id = c("rs1", "rs2", "rs3"),
  chromosome = c("1", "1", "1"),
  base_pair_location = c(101L, 220L, 340L),
  p_value = c(0.05, 1e-4, 0.2)
)

toy_somatic <- data.table(
  sample_id = c("S1", "S2", "S2"),
  chromosome = c("1", "1", "1"),
  start_position = c(105L, 230L, 360L),
  end_position = c(105L, 230L, 360L),
  ref = c("A", "G", "C"),
  alt = c("T", "A", "T")
)

validation <- validate_inputs(
  gwas_sumstats = toy_gwas,
  somatic_maf = toy_somatic,
  verbose = FALSE
)

validation$bundle_type

## ----toy-graph----------------------------------------------------------------
toy_vertices <- data.frame(
  name = c("MYC", "BCL2", "GH01J000001", "GH01J000002"),
  node_id = c("MYC", "BCL2", "GH01J000001", "GH01J000002"),
  node_type = c("gene", "gene", "reg", "reg"),
  stringsAsFactors = FALSE
)

toy_edges <- data.frame(
  from = c("GH01J000001", "GH01J000001", "GH01J000002"),
  to = c("MYC", "BCL2", "MYC"),
  confidence = c(0.9, 0.6, 0.8),
  stringsAsFactors = FALSE
)

toy_graph <- graph_from_data_frame(
  d = toy_edges,
  vertices = toy_vertices,
  directed = TRUE
)

## ----toy-scores---------------------------------------------------------------
toy_gene_germline <- data.table(
  gene_id = c("MYC", "BCL2"),
  zstat = c(3.2, 1.4)
)

toy_reg_germline <- data.table(
  reg_elem_id = c("GH01J000001", "GH01J000002"),
  zstat = c(2.6, 1.1)
)

toy_gene_somatic <- data.table(
  gene_id = c("MYC", "BCL2"),
  zstat = c(1.8, -0.9),
  p_value = c(0.02, 0.4)
)

toy_reg_somatic <- data.table(
  reg_elem_id = c("GH01J000001", "GH01J000002"),
  zstat = c(0.7, -0.4),
  p_value = c(0.3, 0.7)
)

toy_reg_epigenomic <- data.table(
  reg_elem_id = c("GH01J000001", "GH01J000002"),
  zstat = c(2.1, 0.3)
)

## ----scored-graph-------------------------------------------------------------
scored_graph <- build_scored_gene_reg_graph(
  graph = toy_graph,
  gene_germline_scores = toy_gene_germline,
  reg_germline_scores = toy_reg_germline,
  gene_somatic_scores = toy_gene_somatic,
  reg_somatic_scores = toy_reg_somatic,
  reg_epigenomic_scores = toy_reg_epigenomic,
  save_outputs = FALSE,
  verbose = FALSE
)

scored_graph$bundle_type

## ----scored-nodes-------------------------------------------------------------
scored_graph$objects$nodes[, .(
  node_id,
  node_type,
  germline_score,
  somatic_score,
  epigenomic_score
)]

## ----germline-plot------------------------------------------------------------
germline_gene_plot <- plot_germline_gene_scores(
  scored_graph = scored_graph,
  stage = "pre",
  label_features = c("MYC", "BCL2"),
  save_plot = FALSE
)

germline_gene_plot$plot

## ----epigenomic-plot----------------------------------------------------------
epigenomic_reg_plot <- plot_epigenomic_reg_scores(
  scored_graph = scored_graph,
  label_features = c("GH01J000001"),
  save_plot = FALSE
)

epigenomic_reg_plot$plot

## ----runtime-check------------------------------------------------------------
runtime_status <- check_conseguiR_runtime(quiet = TRUE)

list(
  core_ok = runtime_status$core_ok,
  python_stage_ok = runtime_status$python_stage_ok,
  magma_ok = runtime_status$magma_ok
)

## ----full-pipeline-shape, eval=FALSE------------------------------------------
# germline_scores <- prepare_germline_scores(...)
# somatic_scores <- prepare_somatic_scores(...)
# epigenomic_scores <- prepare_epigenomic_scores(...)
# 
# scored_graph <- build_scored_gene_reg_graph(
#   germline_scores = germline_scores,
#   somatic_scores = somatic_scores,
#   epigenomic_scores = epigenomic_scores
# )
# 
# diffusion <- run_gene_reg_diffusion(scored_graph = scored_graph)
# selected_subgraph <- call_selected_subgraph(diffusion = diffusion)

## ----session-info-------------------------------------------------------------
sessionInfo()

