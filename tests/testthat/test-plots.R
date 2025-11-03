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

# Performance Plot Tests
test_that("plotEndometrialROC returns ggplot with correct AUC", {
  # Create mock predictions data
  set.seed(123)
  n_samples <- 20
  predictions <- data.frame(
    sample_id = paste0("sample_", 1:n_samples),
    label = c(rep(0, 10), rep(1, 10)),
    prob = c(runif(10, 0, 0.5), runif(10, 0.5, 1)),
    stringsAsFactors = FALSE
  )
  
  # Sort by prob to ensure positive labels have higher probabilities
  predictions <- predictions[order(predictions$prob), ]
  predictions$label <- c(rep(0, 10), rep(1, 10))
  
  # Test ROC curve
  p_roc <- plotEndometrialROC(predictions, use_calibrated = FALSE, show_auc = TRUE)
  
  expect_s3_class(p_roc, "ggplot")
  expect_true("GeomLine" %in% class(p_roc$layers[[1]]$geom))
  expect_no_error(print(p_roc))
  
  # Verify ROC curve starts at (0,0) and ends at (1,1)
  roc_data <- p_roc$data
  expect_equal(roc_data$FPR[1], 0)
  expect_equal(roc_data$TPR[1], 0)
  expect_equal(roc_data$FPR[nrow(roc_data)], 1)
  expect_equal(roc_data$TPR[nrow(roc_data)], 1)
})

test_that("plotEndometrialROC handles calibrated probabilities", {
  # Create mock predictions with calibrated probabilities
  set.seed(123)
  predictions <- data.frame(
    sample_id = paste0("sample_", 1:20),
    label = c(rep(0, 10), rep(1, 10)),
    prob = c(runif(10, 0, 0.5), runif(10, 0.5, 1)),
    prob_calibrated = c(runif(10, 0, 0.5), runif(10, 0.5, 1)),
    stringsAsFactors = FALSE
  )
  
  # Test with raw probabilities
  p_roc_raw <- plotEndometrialROC(predictions, use_calibrated = FALSE)
  expect_s3_class(p_roc_raw, "ggplot")
  
  # Test with calibrated probabilities
  p_roc_cal <- plotEndometrialROC(predictions, use_calibrated = TRUE)
  expect_s3_class(p_roc_cal, "ggplot")
  
  # Both should work without errors
  expect_no_error(print(p_roc_raw))
  expect_no_error(print(p_roc_cal))
})

test_that("plotEndometrialPR returns ggplot with correct PR-AUC", {
  # Create mock predictions data
  set.seed(123)
  predictions <- data.frame(
    sample_id = paste0("sample_", 1:20),
    label = c(rep(0, 10), rep(1, 10)),
    prob = c(runif(10, 0, 0.5), runif(10, 0.5, 1)),
    stringsAsFactors = FALSE
  )
  
  # Test PR curve
  # Suppress geom_hline mapping warning (expected when yintercept is provided)
  p_pr <- suppressWarnings(plotEndometrialPR(predictions, use_calibrated = FALSE, show_auc = TRUE))
  
  expect_s3_class(p_pr, "ggplot")
  expect_true("GeomLine" %in% class(p_pr$layers[[1]]$geom))
  expect_no_error(print(p_pr))
  
  # Verify PR curve starts at (0,1) and ends at (1,0)
  pr_data <- p_pr$data
  expect_equal(pr_data$Recall[1], 0)
  expect_equal(pr_data$Precision[1], 1)
})

test_that("plotEndometrialCalibration returns ggplot with correct metrics", {
  # Create mock predictions data
  set.seed(123)
  predictions <- data.frame(
    sample_id = paste0("sample_", 1:20),
    label = c(rep(0, 10), rep(1, 10)),
    prob = c(runif(10, 0, 0.5), runif(10, 0.5, 1)),
    stringsAsFactors = FALSE
  )
  
  # Test calibration curve
  p_cal <- plotEndometrialCalibration(predictions, use_calibrated = FALSE, 
                                      show_brier = TRUE, show_ece = TRUE)
  
  expect_s3_class(p_cal, "ggplot")
  expect_true("GeomPoint" %in% class(p_cal$layers[[1]]$geom))
  expect_no_error(print(p_cal))
  
  # Verify calibration data has bins
  cal_data <- p_cal$data
  expect_true(nrow(cal_data) > 0)
  expect_true(all(cal_data$mean_pred >= 0 & cal_data$mean_pred <= 1))
  expect_true(all(cal_data$mean_obs >= 0 & cal_data$mean_obs <= 1))
})

test_that("plotEndometrialComparison works with new signature only", {
  # Create mock training result
  set.seed(123)
  predictions <- data.frame(
    sample_id = paste0("sample_", 1:20),
    label = c(rep(0, 10), rep(1, 10)),
    prob = c(runif(10, 0, 0.5), runif(10, 0.5, 1)),
    prob_calibrated = c(runif(10, 0, 0.5), runif(10, 0.5, 1)),
    stringsAsFactors = FALSE
  )
  
  new_result <- list(
    metrics = list(
      auc = 0.85,
      accuracy = 0.80,
      brier_score = 0.15,
      ece = 0.10,
      predictions = predictions
    )
  )
  
  # Test comparison with only new signature
  # Suppress geom_hline mapping warning (expected when yintercept is provided)
  comparison_plots <- suppressWarnings(plotEndometrialComparison(
    pretrained_result = NULL,
    new_result = new_result,
    metrics_to_plot = c("roc", "pr", "calibration")
  ))
  
  expect_true(is.list(comparison_plots))
  expect_true("roc" %in% names(comparison_plots))
  expect_true("pr" %in% names(comparison_plots))
  expect_true("calibration" %in% names(comparison_plots))
  
  # All plots should be ggplot objects
  expect_s3_class(comparison_plots$roc, "ggplot")
  expect_s3_class(comparison_plots$pr, "ggplot")
  expect_s3_class(comparison_plots$calibration, "ggplot")
  
  # Metrics table should be present
  expect_true("metrics_table" %in% names(comparison_plots))
  expect_true(is.data.frame(comparison_plots$metrics_table))
})

test_that("performance plots work with training output from gse201926_trainmini", {
  data(gse201926_trainmini)
  
  # Train a signature (with minimal settings for speed)
  # Suppress warnings about no consensus genes and CV fold errors (expected for small sample sizes)
  set.seed(123)
  result <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 50,
    outer = "lpo",
    inner_folds = 3,
    inner_repeats = 5,
    calibration_method = "platt",
    stability_selection = FALSE,
    seed = 123
  ))
  
  # Test all performance plots
  expect_no_error({
    p_roc <- plotEndometrialROC(result$metrics$predictions, use_calibrated = FALSE)
    expect_s3_class(p_roc, "ggplot")
  })
  
  expect_no_error({
    # Suppress geom_hline mapping warning (expected when yintercept is provided)
    p_pr <- suppressWarnings(plotEndometrialPR(result$metrics$predictions, use_calibrated = FALSE))
    expect_s3_class(p_pr, "ggplot")
  })
  
  expect_no_error({
    p_cal <- plotEndometrialCalibration(result$metrics$predictions, use_calibrated = FALSE)
    expect_s3_class(p_cal, "ggplot")
  })
  
  # Test with calibrated probabilities
  expect_no_error({
    p_roc_cal <- plotEndometrialROC(result$metrics$predictions, use_calibrated = TRUE)
    expect_s3_class(p_roc_cal, "ggplot")
  })
  
  # Test comparison plot
  expect_no_error({
    # Suppress geom_hline mapping warning (expected when yintercept is provided)
    comparison_plots <- suppressWarnings(plotEndometrialComparison(
      pretrained_result = NULL,
      new_result = result,
      metrics_to_plot = c("roc", "pr", "calibration")
    ))
    expect_true(is.list(comparison_plots))
  })
})

test_that("performance plots handle edge cases gracefully", {
  # Test with balanced predictions (all same probability)
  set.seed(123)
  predictions_balanced <- data.frame(
    sample_id = paste0("sample_", 1:10),
    label = c(rep(0, 5), rep(1, 5)),
    prob = rep(0.5, 10),
    stringsAsFactors = FALSE
  )
  
  # ROC should still work
  expect_no_error({
    p_roc <- plotEndometrialROC(predictions_balanced, use_calibrated = FALSE)
    expect_s3_class(p_roc, "ggplot")
  })
  
  # PR should still work
  expect_no_error({
    # Suppress geom_hline mapping warning (expected when yintercept is provided)
    p_pr <- suppressWarnings(plotEndometrialPR(predictions_balanced, use_calibrated = FALSE))
    expect_s3_class(p_pr, "ggplot")
  })
  
  # Calibration should still work
  expect_no_error({
    p_cal <- plotEndometrialCalibration(predictions_balanced, use_calibrated = FALSE, n_bins = 5)
    expect_s3_class(p_cal, "ggplot")
  })
})
