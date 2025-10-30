test_that("endo_load_demo returns expected structure", {
  demo <- endo_load_demo()
  expect_type(demo, "list")
  expect_true(all(c("counts", "pheno", "annot") %in% names(demo)))
})

test_that("endo_load_annotation(minimal=TRUE) returns a data.frame when available", {
  ann <- endo_load_annotation(minimal = TRUE)
  expect_true(is.null(ann) || is.data.frame(ann))
})

test_that("endo_load_gse201926(sample_only=TRUE) returns demo-like structure", {
  demo <- endo_load_gse201926(sample_only = TRUE)
  expect_type(demo, "list")
  expect_true(all(c("counts", "pheno", "annot") %in% names(demo)))
})


