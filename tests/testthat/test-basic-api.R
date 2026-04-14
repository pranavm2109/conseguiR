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

  res <- validate_inputs(
    gwas_sumstats = gwas,
    somatic_maf = somatic
  )

  expect_s3_class(res, "conseguiR_bundle")
  expect_identical(res$bundle_type, "validation")
  expect_true(is.data.frame(res$objects$gwas))
  expect_true(is.data.frame(res$objects$somatic_maf))
})

test_that("runtime checks return a structured status object", {
  status <- check_conseguiR_runtime(quiet = TRUE)

  expect_true(is.list(status))
  expect_true("ok" %in% names(status))
  expect_true("core_r_packages" %in% names(status))
  expect_true("python_stage_ok" %in% names(status))
  expect_true("python_modules" %in% names(status))
})
