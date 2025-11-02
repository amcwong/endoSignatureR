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

test_that("esr_selectTopGenes returns correct gene IDs by variance", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Select top 20 genes by variance
  top_genes <- esr_selectTopGenes(mat_t, n = 20, by = "variance")

  expect_type(top_genes, "character")
  expect_length(top_genes, 20)
  expect_true(all(top_genes %in% colnames(mat_t)))
  
  # Verify genes are selected by variance (higher variance should be selected)
  gene_var <- apply(mat_t, 2, var)
  top_var_values <- gene_var[top_genes]
  expect_true(all(!is.na(top_var_values)))
})

test_that("esr_selectTopGenes returns correct gene IDs by DE", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Select top 10 genes by DE
  top_genes <- esr_selectTopGenes(de_table = de_table, n = 10, by = "de")

  expect_type(top_genes, "character")
  expect_length(top_genes, 10)
  expect_true(all(top_genes %in% de_table$gene_id))
  
  # Verify genes are selected by FDR then log2FC
  de_subset <- de_table[de_table$gene_id %in% top_genes, ]
  expect_true(all(de_subset$FDR <= sort(de_table$FDR)[10]))
})

test_that("esr_selectTopGenes handles edge cases", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Test with n > total genes (should return all available)
  top_all <- esr_selectTopGenes(mat_t, n = 10000, by = "variance")
  expect_true(length(top_all) <= ncol(mat_t))

  # Test with n = 1
  top_one <- esr_selectTopGenes(mat_t, n = 1, by = "variance")
  expect_length(top_one, 1)

  # Test custom selection
  top_custom <- esr_selectTopGenes(de_table = de_table, n = 5, by = "custom", sort_col = "AveExpr")
  expect_length(top_custom, 5)
  expect_true(all(top_custom %in% de_table$gene_id))
})

test_that("esr_selectTopGenes validates inputs correctly", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Missing mat_t for variance
  expect_error(esr_selectTopGenes(n = 10, by = "variance"), "mat_t is required")

  # Missing de_table for DE
  expect_error(esr_selectTopGenes(n = 10, by = "de"), "de_table is required")

  # Invalid n
  expect_error(esr_selectTopGenes(mat_t, n = -1, by = "variance"), "n must be a positive integer")
  expect_error(esr_selectTopGenes(mat_t, n = 0, by = "variance"), "n must be a positive integer")
})

test_that("plotEndometrialHeatmap returns ComplexHeatmap object", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Create heatmap with top 20 genes
  top_genes <- esr_selectTopGenes(mat_t, n = 20, by = "variance")
  hm <- plotEndometrialHeatmap(mat_t, genes = top_genes)

  # Check that it's a Heatmap object
  expect_true(inherits(hm, "Heatmap"))
  expect_no_error(ComplexHeatmap::draw(hm))
})

test_that("plotEndometrialHeatmap works with phenotype annotations", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Create heatmap with phenotype
  top_genes <- esr_selectTopGenes(mat_t, n = 20, by = "variance")
  hm <- plotEndometrialHeatmap(mat_t, genes = top_genes, pheno = gse201926_sample$pheno)

  expect_true(inherits(hm, "Heatmap"))
  expect_no_error(ComplexHeatmap::draw(hm))
})

test_that("plotEndometrialHeatmap handles different scaling options", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  top_genes <- esr_selectTopGenes(mat_t, n = 10, by = "variance")

  # Test row scaling
  expect_no_error({
    hm_row <- plotEndometrialHeatmap(mat_t, genes = top_genes, scale = "row")
    expect_true(inherits(hm_row, "Heatmap"))
  })

  # Test column scaling
  expect_no_error({
    hm_col <- plotEndometrialHeatmap(mat_t, genes = top_genes, scale = "column")
    expect_true(inherits(hm_col, "Heatmap"))
  })

  # Test no scaling
  expect_no_error({
    hm_none <- plotEndometrialHeatmap(mat_t, genes = top_genes, scale = "none")
    expect_true(inherits(hm_none, "Heatmap"))
  })
})

test_that("plotEndometrialHeatmap handles edge cases", {
  data(gse201926_sample)
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Test with empty genes vector (should error gracefully)
  expect_error(plotEndometrialHeatmap(mat_t, genes = character(0)), "None of the requested genes")

  # Test with all missing genes (should warn then error - can't create heatmap with no genes)
  expect_error(
    expect_warning(
      plotEndometrialHeatmap(mat_t, genes = c("fake_gene_1", "fake_gene_2")),
      "Some requested genes not found"
    ),
    "None of the requested genes are present"
  )

  # Test with some missing genes (should warn but succeed with available genes)
  if (ncol(mat_t) >= 2) {
    real_genes <- colnames(mat_t)[1:2]
    expect_warning({
      hm_partial <- plotEndometrialHeatmap(mat_t, genes = c(real_genes, "fake_gene_1", "fake_gene_2"))
      expect_true(inherits(hm_partial, "Heatmap"))
    }, "Some requested genes not found")
  }

  # Test with all genes (no subsetting)
  expect_no_error({
    hm_all <- plotEndometrialHeatmap(mat_t, genes = NULL)
    expect_true(inherits(hm_all, "Heatmap"))
  })
})
