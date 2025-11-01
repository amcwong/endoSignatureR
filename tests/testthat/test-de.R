test_that("esr_analyzeDifferentialExpression accepts valid inputs and returns expected structure", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Check structure
  expect_true(is.data.frame(de_table))
  expect_true(nrow(de_table) > 0)

  # Check required columns
  required_cols <- c("gene_id", "log2FC", "pvalue", "FDR", "AveExpr", "t")
  expect_true(all(required_cols %in% names(de_table)))

  # Check column types
  expect_true(is.character(de_table$gene_id))
  expect_true(is.numeric(de_table$log2FC))
  expect_true(is.numeric(de_table$pvalue))
  expect_true(is.numeric(de_table$FDR))
  expect_true(is.numeric(de_table$AveExpr))
  expect_true(is.numeric(de_table$t))
})

test_that("esr_analyzeDifferentialExpression validates input dimensions and ID matching", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Test with valid inputs (should work)
  expect_no_error({
    de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
  })

  # Test with mismatched sample IDs (should handle gracefully by matching common samples)
  pheno_mismatch <- gse201926_sample$pheno
  pheno_mismatch$sample_id[1] <- "nonexistent_sample"

  expect_no_error({
    de_table <- esr_analyzeDifferentialExpression(mat_t, pheno_mismatch)
  })

  # Should still have results for matching samples
  expect_true(nrow(de_table) > 0)
})

test_that("esr_analyzeDifferentialExpression produces deterministic output with fixed seed", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Set seed and run
  set.seed(123)
  de_table1 <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno, seed = 123)

  # Set seed again and run
  set.seed(123)
  de_table2 <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno, seed = 123)

  # Results should be identical
  expect_equal(de_table1, de_table2)
})

test_that("esr_analyzeDifferentialExpression handles edge cases", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Test with missing required columns (should error)
  pheno_bad <- gse201926_sample$pheno[, c("sample_id"), drop = FALSE]

  expect_error({
    esr_analyzeDifferentialExpression(mat_t, pheno_bad)
  })

  # Test with single group (should error)
  pheno_single <- gse201926_sample$pheno
  pheno_single$group <- "PS"

  expect_error({
    esr_analyzeDifferentialExpression(mat_t, pheno_single)
  })
})

test_that("esr_analyzeDifferentialExpression verifies sorting by FDR then absolute log2FC", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)

  # Check sorting: FDR ascending, then absolute log2FC descending
  # Use vectorized check instead of loop to avoid thousands of expectations
  if (nrow(de_table) > 1) {
    # Compare each row with the next
    fdr_next <- de_table$FDR[-1]
    fdr_curr <- de_table$FDR[-nrow(de_table)]
    log2fc_abs_next <- abs(de_table$log2FC[-1])
    log2fc_abs_curr <- abs(de_table$log2FC[-nrow(de_table)])

    # Either FDR increases, or if FDR is equal, absolute log2FC decreases or stays same
    sorting_ok <- fdr_curr < fdr_next |
      (fdr_curr == fdr_next & log2fc_abs_curr >= log2fc_abs_next)

    # Single vectorized expectation instead of thousands
    expect_true(all(sorting_ok),
      info = paste(
        "Table is not properly sorted. First violation at row",
        which.min(sorting_ok) + 1
      )
    )
  }
})

test_that("esr_analyzeDifferentialExpression works with gse201926_sample", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Should complete without errors
  expect_no_error({
    de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
  })

  # Should have reasonable number of genes (at least some)
  expect_true(nrow(de_table) > 0)
  expect_true(nrow(de_table) <= ncol(mat_t)) # Can't have more genes than input

  # Check that FDR values are in valid range [0, 1]
  expect_true(all(de_table$FDR >= 0 & de_table$FDR <= 1))

  # Check that p-values are in valid range [0, 1]
  expect_true(all(de_table$pvalue >= 0 & de_table$pvalue <= 1))
})
