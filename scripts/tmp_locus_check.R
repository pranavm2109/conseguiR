#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
})

repo_root <- getwd()
runtime_env <- globalenv()

sys.source(file.path(repo_root, "R", "zzz.R"), envir = runtime_env)
sys.source(file.path(repo_root, "R", "backend_resources.R"), envir = runtime_env)
sys.source(file.path(repo_root, "R", "user_api.R"), envir = runtime_env)

.conseguiR_state$pkg_root <- repo_root
setwd(repo_root)
initialize_backend_graphs(strict = FALSE, quiet = TRUE)
.conseguiR_load_external_api()

assign("copy", data.table::copy, envir = .conseguiR_runtime_env)
assign("as.data.table", data.table::as.data.table, envir = .conseguiR_runtime_env)
assign("fread", data.table::fread, envir = .conseguiR_runtime_env)
assign("fwrite", data.table::fwrite, envir = .conseguiR_runtime_env)

prepare_epigenomic_scores <- get("prepare_epigenomic_scores", envir = .conseguiR_runtime_env)
build_scored_gene_reg_graph <- get("build_scored_gene_reg_graph", envir = .conseguiR_runtime_env)
plot_locus_context <- get("plot_locus_context", envir = .conseguiR_runtime_env)

gene_scores <- fread(
  file.path(repo_root, "data/processed/test_outputs/germline/magma_zstat_20260427125250_983097.tsv"),
  showProgress = FALSE
)

ccre <- fread(.conseguiR_default_reg_loc_path(), header = FALSE, showProgress = FALSE)
ccre_sub <- ccre[V2 == 7 & V3 <= 106500000 & V4 >= 106100000]

outdir <- file.path(repo_root, "data/processed/test_outputs/plotting")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

ccre_path <- file.path(outdir, "ccre_nampt_subset.loc")
fwrite(ccre_sub, ccre_path, sep = "\t", col.names = FALSE)

epi <- prepare_epigenomic_scores(
  reg_ref_path = ccre_path,
  bw_files = c(
    file.path(repo_root, "data/raw/Testing/SRR1020514_DLBCL_P265_H3K27ac_ChIPseq.bw"),
    file.path(repo_root, "data/raw/Testing/SRR1020516_DLBCL_P286_H3K27ac_ChIPseq.bw"),
    file.path(repo_root, "data/raw/Testing/SRR1020518_DLBCL_P397_H3K27ac_ChIPseq.bw")
  ),
  verbose = FALSE
)

sg <- build_scored_gene_reg_graph(
  graph_rds_path = .conseguiR_backend_paths()$gene_reg_graph_rds,
  gene_germline_scores = gene_scores,
  reg_epigenomic_scores = epi$reg_scores,
  save_outputs = FALSE,
  verbose = FALSE
)

post <- list(objects = list(nodes = sg$objects$nodes, edges = sg$objects$edges))

plt <- plot_locus_context(
  chromosome = "7",
  start = 106100000,
  end = 106500000,
  postdiff_gene_reg_graph = post,
  gwas_sumstats = NULL,
  label_top_gwas_snp = FALSE,
  title = "NAMPT locus context",
  save_plot = FALSE,
  verbose = FALSE
)

cat("plot class", paste(class(plt$objects$plot), collapse = ", "), "\n")
cat(
  "tracks", nrow(plt$objects$plot_data$tracks),
  "features", nrow(plt$objects$plot_data$features),
  "links", nrow(plt$objects$plot_data$links), "\n"
)
print(plt$objects$plot_data$locus)
