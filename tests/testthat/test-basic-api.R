test_that("validate_inputs returns a validation bundle for in-memory tables", {
  gwas <- data.table::data.table(
    variant_id = c("rs1", "rs2"),
    chromosome = c("1", "1"),
    base_pair_location = c(101L, 202L),
    p_value = c(0.05, 1e-3)
  )

  somatic <- data.table::data.table(
    sample_id = c("S1", "S2"),
    chromosome = c("1", "1"),
    start_position = c(101L, 202L),
    end_position = c(101L, 202L),
    ref = c("A", "G"),
    alt = c("T", "A")
  )

  res <- conseguiR::validate_inputs(
    gwas_sumstats = gwas,
    somatic_maf = somatic
  )

  expect_s3_class(res, "conseguiR_bundle")
  expect_identical(res$bundle_type, "validation")
  expect_true(is.data.frame(res$objects$gwas))
  expect_true(is.data.frame(res$objects$somatic_maf))
})

test_that("runtime checks return a structured status object", {
  status <- conseguiR::check_conseguiR_runtime(quiet = TRUE)

  expect_true(is.list(status))
  expect_true("ok" %in% names(status))
  expect_true("core_r_packages" %in% names(status))
  expect_true("python_stage_ok" %in% names(status))
  expect_true("python_modules" %in% names(status))
})

test_that("run_conseguiR exposes candidate_pool_size on the public API", {
  run_formals <- formals(conseguiR::run_conseguiR)

  expect_true("candidate_pool_size" %in% names(run_formals))
  expect_equal(as.integer(run_formals$candidate_pool_size), 400L)
})

test_that("validated locus plotting is exposed on the public API", {
  plot_formals <- formals(conseguiR::plot_validated_locus_context)

  expect_true(is.function(conseguiR::plot_validated_locus_context))
  expect_true("strict_gene_filter" %in% names(plot_formals))
  expect_identical(plot_formals$strict_gene_filter, TRUE)
  expect_true("label_top_gwas_snp" %in% names(plot_formals))
  expect_true("label_top_lit_snps" %in% names(plot_formals))
})

test_that("internal pipeline arg assembly forwards candidate_pool_size", {
  args <- conseguiR:::.conseguiR_pipeline_args(
    gwas_sumstats = "gwas.tsv",
    somatic_maf = "somatic.maf",
    reg_ref_path = "regs.loc",
    reference_bfile = "ref_prefix",
    dndscv_refdb = "refdb.rda",
    epigenomic_track_dir = "tracks",
    epigenomic_tracks = c("a.bw", "b.bw"),
    paths = list(
      graph_rds_path = "graph.rds",
      gg_nodes_path = "gg_nodes.tsv",
      gg_edges_path = "gg_edges.tsv"
    ),
    output_dir = "out",
    target_genes = 25L,
    candidate_pool_size = 125L,
    germline_args = list(alpha = 1),
    somatic_args = list(beta = 2),
    epigenomic_args = list(gamma = 3),
    scored_graph_args = list(delta = 4),
    diffusion_args = list(epsilon = 5),
    subgraph_args = list(zeta = 6),
    plot_args = list(eta = 7),
    verbose = FALSE
  )

  expect_identical(args$candidate_pool_size, 125L)
  expect_identical(args$target_genes, 25L)
})
