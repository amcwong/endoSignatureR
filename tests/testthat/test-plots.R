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

test_that("plotEndometrialMA returns ggplot object with expected layers", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  p <- plotEndometrialMA(de_table)

  expect_s3_class(p, "ggplot")
  expect_true("GeomPoint" %in% class(p$layers[[1]]$geom))
  expect_no_error(print(p))
})

test_that("plotEndometrialMA works with highlight_genes", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Highlight top 5 genes by absolute log2FC
  top_genes <- head(de_table$gene_id[order(abs(de_table$log2FC), decreasing = TRUE)], 5)

  p <- plotEndometrialMA(de_table, highlight_genes = top_genes)

  expect_s3_class(p, "ggplot")
  expect_no_error(print(p))
})

test_that("plotEndometrialVolcano returns ggplot object with expected layers", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  p <- plotEndometrialVolcano(de_table)

  expect_s3_class(p, "ggplot")
  expect_true("GeomPoint" %in% class(p$layers[[1]]$geom))
  expect_no_error(print(p))
})

test_that("plotEndometrialVolcano works with highlight_genes", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Highlight top 5 genes by FDR
  top_genes <- head(de_table$gene_id[order(de_table$FDR)], 5)

  p <- plotEndometrialVolcano(de_table, highlight_genes = top_genes)

  expect_s3_class(p, "ggplot")
  expect_no_error(print(p))
})

test_that("plotEndometrialMA and plotEndometrialVolcano handle edge cases", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Test with different thresholds
  expect_no_error({
    p1 <- plotEndometrialMA(de_table, fdr_threshold = 0.01, log2fc_threshold = 0.5)
    p2 <- plotEndometrialVolcano(de_table, fdr_threshold = 0.01, log2fc_threshold = 0.5)
  })

  expect_s3_class(p1, "ggplot")
  expect_s3_class(p2, "ggplot")

  # Test with empty highlight_genes (should work)
  expect_no_error({
    p3 <- plotEndometrialMA(de_table, highlight_genes = character(0))
    p4 <- plotEndometrialVolcano(de_table, highlight_genes = character(0))
  })

  expect_s3_class(p3, "ggplot")
  expect_s3_class(p4, "ggplot")
})

test_that("MA and Volcano plots can be generated from gse201926_sample DE results", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  expect_no_error({
    p_ma <- plotEndometrialMA(de_table)
    p_volcano <- plotEndometrialVolcano(de_table)
  })

  expect_s3_class(p_ma, "ggplot")
  expect_s3_class(p_volcano, "ggplot")
})
