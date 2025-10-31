test_that("plotEndometrialPCA returns ggplot object with expected layers", {
  data(gse201926_sample)
  
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  # Test with pheno
  p <- plotEndometrialPCA(mat_t, pheno = gse201926_sample$pheno)
  
  expect_s3_class(p, "ggplot")
  expect_true("GeomPoint" %in% class(p$layers[[1]]$geom))
  expect_no_error(print(p))
})

test_that("plotEndometrialPCA works without pheno", {
  data(gse201926_sample)
  
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  p <- plotEndometrialPCA(mat_t)
  
  expect_s3_class(p, "ggplot")
  # expect_no_error(print(p))
})

test_that("plotEndometrialPCA is deterministic with fixed seed", {
  data(gse201926_sample)
  
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  set.seed(123)
  p1 <- plotEndometrialPCA(mat_t, pheno = gse201926_sample$pheno)
  
  set.seed(123)
  p2 <- plotEndometrialPCA(mat_t, pheno = gse201926_sample$pheno)
  
  # PCA data should be identical
  expect_equal(p1$data, p2$data)
})

test_that("plotEndometrialLibsize returns ggplot object", {
  data(gse201926_sample)
  
  p <- plotEndometrialLibsize(gse201926_sample$counts, pheno = gse201926_sample$pheno)
  
  expect_s3_class(p, "ggplot")
  expect_no_error(print(p))
})

test_that("plotEndometrialLibsize works without pheno", {
  data(gse201926_sample)
  
  p <- plotEndometrialLibsize(gse201926_sample$counts)
  
  expect_s3_class(p, "ggplot")
  expect_no_error(print(p))
})

test_that("plotEndometrialZeros returns ggplot object", {
  data(gse201926_sample)
  
  p_gene <- plotEndometrialZeros(gse201926_sample$counts, by = "gene")
  p_sample <- plotEndometrialZeros(gse201926_sample$counts, by = "sample")
  
  expect_s3_class(p_gene, "ggplot")
  expect_s3_class(p_sample, "ggplot")
  expect_no_error(print(p_gene))
  expect_no_error(print(p_sample))
})

test_that("plotEndometrialZeros handles both by options", {
  data(gse201926_sample)
  
  expect_no_error(plotEndometrialZeros(gse201926_sample$counts, by = "gene"))
  expect_no_error(plotEndometrialZeros(gse201926_sample$counts, by = "sample"))
})

test_that("all plots can be generated from gse201926_sample without errors", {
  data(gse201926_sample)
  
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  expect_no_error({
    p1 <- plotEndometrialPCA(mat_t, pheno = gse201926_sample$pheno)
    p2 <- plotEndometrialLibsize(gse201926_sample$counts, pheno = gse201926_sample$pheno)
    p3 <- plotEndometrialZeros(gse201926_sample$counts, by = "gene")
    p4 <- plotEndometrialZeros(gse201926_sample$counts, by = "sample")
  })
  
  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")
  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

