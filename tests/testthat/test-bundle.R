test_that("esr_computeQCMetrics returns expected structure", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  qc_metrics <- esr_computeQCMetrics(
    counts = gse201926_sample$counts,
    mat_t = mat_t,
    pheno = gse201926_sample$pheno
  )

  # Check structure
  expect_true(is.data.frame(qc_metrics))
  expect_true(nrow(qc_metrics) > 0)

  # Check required columns
  required_cols <- c("sample_id", "library_size", "pct_zeros", "n_genes")
  expect_true(all(required_cols %in% names(qc_metrics)))

  # Check column types
  expect_true(is.character(qc_metrics$sample_id))
  expect_true(is.numeric(qc_metrics$library_size))
  expect_true(is.numeric(qc_metrics$pct_zeros))
  expect_true(is.numeric(qc_metrics$n_genes))

  # Check additional columns when mat_t provided
  if ("genes_after_filter" %in% names(qc_metrics)) {
    expect_true(is.numeric(qc_metrics$genes_after_filter))
  }
})

test_that("esr_computeQCMetrics computes correct values", {
  data(gse201926_sample)

  # Calculate library sizes manually
  lib_sizes_manual <- colSums(gse201926_sample$counts, na.rm = TRUE)

  qc_metrics <- esr_computeQCMetrics(counts = gse201926_sample$counts)

  # Library sizes should match
  expect_equal(qc_metrics$library_size, lib_sizes_manual, ignore_attr = TRUE)

  # Percentage zeros should be reasonable (between 0 and 100)
  expect_true(all(qc_metrics$pct_zeros >= 0 & qc_metrics$pct_zeros <= 100))

  # Number of genes should be reasonable
  expect_true(all(qc_metrics$n_genes >= 0))
  expect_true(all(qc_metrics$n_genes <= nrow(gse201926_sample$counts)))
})

test_that("esr_computeQCMetrics works with pheno", {
  data(gse201926_sample)

  qc_metrics <- esr_computeQCMetrics(
    counts = gse201926_sample$counts,
    pheno = gse201926_sample$pheno
  )

  # Should have group column if pheno provided
  if ("group" %in% names(qc_metrics)) {
    expect_true(all(qc_metrics$group %in% c("PS", "PIS")))
  }
})

test_that("esr_computeQCMetrics handles missing mat_t", {
  data(gse201926_sample)

  # Should work without mat_t
  expect_no_error({
    qc_metrics <- esr_computeQCMetrics(counts = gse201926_sample$counts)
  })

  # Should have required columns
  expect_true("sample_id" %in% names(qc_metrics))
  expect_true("library_size" %in% names(qc_metrics))
})

test_that("esr_createAnalysisBundle returns expected structure", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
  top_genes <- esr_selectTopGenes(de_table = de_table, n = 10, by = "de")
  qc_metrics <- esr_computeQCMetrics(
    counts = gse201926_sample$counts,
    mat_t = mat_t,
    pheno = gse201926_sample$pheno
  )

  bundle <- esr_createAnalysisBundle(
    counts_t = mat_t,
    de_table = de_table,
    selected_genes = top_genes,
    qc_metrics = qc_metrics,
    pheno = gse201926_sample$pheno,
    annot = gse201926_sample$annot
  )

  # Check structure
  expect_true(is.list(bundle))
  expect_true("esr_analysis_bundle" %in% class(bundle))

  # Check required component
  expect_true("counts_t" %in% names(bundle))
  expect_true(is.matrix(bundle$counts_t))

  # Check optional components
  expect_true("de_table" %in% names(bundle))
  expect_true("selected_genes" %in% names(bundle))
  expect_true("qc_metrics" %in% names(bundle))
  expect_true("pheno" %in% names(bundle))
  expect_true("annot" %in% names(bundle))
})

test_that("esr_createAnalysisBundle validates component alignment", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
  qc_metrics <- esr_computeQCMetrics(
    counts = gse201926_sample$counts,
    mat_t = mat_t,
    pheno = gse201926_sample$pheno
  )

  # Should work with aligned components
  expect_no_error({
    bundle <- esr_createAnalysisBundle(
      counts_t = mat_t,
      de_table = de_table,
      qc_metrics = qc_metrics,
      pheno = gse201926_sample$pheno
    )
  })

  # Sample IDs should align
  sample_ids_counts <- rownames(bundle$counts_t)
  sample_ids_qc <- bundle$qc_metrics$sample_id
  sample_ids_pheno <- bundle$pheno$sample_id

  expect_true(length(intersect(sample_ids_counts, sample_ids_qc)) > 0)
  expect_true(length(intersect(sample_ids_counts, sample_ids_pheno)) > 0)
})

test_that("esr_createAnalysisBundle handles missing optional components", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Should work with only counts_t
  expect_no_error({
    bundle <- esr_createAnalysisBundle(counts_t = mat_t)
  })

  # Optional components should be NULL
  expect_null(bundle$de_table)
  expect_null(bundle$selected_genes)
  expect_null(bundle$qc_metrics)
})

test_that("esr_createAnalysisBundle warns on misaligned components", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Create misaligned qc_metrics
  qc_metrics_bad <- data.frame(
    sample_id = c("bad1", "bad2", "bad3"),
    library_size = c(1000, 2000, 3000),
    pct_zeros = c(10, 20, 30),
    n_genes = c(100, 200, 300),
    stringsAsFactors = FALSE
  )

  # Should warn but not error
  expect_warning(
    {
      bundle <- esr_createAnalysisBundle(
        counts_t = mat_t,
        qc_metrics = qc_metrics_bad
      )
    },
    "No matching sample IDs"
  )
})

test_that("complete Mode 3 workflow produces valid bundle", {
  data(gse201926_sample)

  # Set seed for reproducibility
  set.seed(123)

  # Run complete Mode 3 workflow
  # 1. Validation
  validation_result <- esr_validateEndometrial(
    gse201926_sample$counts,
    gse201926_sample$pheno,
    annot = gse201926_sample$annot
  )

  # 2. Transformation
  mat_t <- esr_transform_log1p_cpm(validation_result$X)

  # 3. QC metrics
  qc_metrics <- esr_computeQCMetrics(
    counts = validation_result$X,
    mat_t = mat_t,
    pheno = validation_result$pheno
  )

  # 4. DE analysis
  de_table <- esr_analyzeDifferentialExpression(
    mat_t,
    validation_result$pheno,
    seed = 123
  )

  # 5. Gene selection
  selected_genes <- esr_selectTopGenes(
    de_table = de_table,
    n = 10,
    by = "de"
  )

  # 6. Create bundle
  bundle <- esr_createAnalysisBundle(
    counts_t = mat_t,
    de_table = de_table,
    selected_genes = selected_genes,
    qc_metrics = qc_metrics,
    pheno = validation_result$pheno,
    annot = validation_result$annot,
    validation_issues = validation_result$issues
  )

  # Verify bundle structure
  expect_true(is.list(bundle))
  expect_true("counts_t" %in% names(bundle))
  expect_true("de_table" %in% names(bundle))
  expect_true("selected_genes" %in% names(bundle))
  expect_true("qc_metrics" %in% names(bundle))

  # Verify component types
  expect_true(is.matrix(bundle$counts_t))
  expect_true(is.data.frame(bundle$de_table))
  expect_true(is.character(bundle$selected_genes))
  expect_true(is.data.frame(bundle$qc_metrics))
})

test_that("Mode 3 full pipeline with visualization: validate → transform → QC → DE → visualize → bundle", {
  data(gse201926_sample)

  # Set seed for reproducibility
  set.seed(123)

  # 1. Validation
  validation_result <- esr_validateEndometrial(
    gse201926_sample$counts,
    gse201926_sample$pheno,
    annot = gse201926_sample$annot
  )

  # 2. Transformation
  mat_t <- esr_transform_log1p_cpm(validation_result$X)

  # 3. QC metrics
  qc_metrics <- esr_computeQCMetrics(
    counts = validation_result$X,
    mat_t = mat_t,
    pheno = validation_result$pheno
  )

  # 4. DE analysis
  de_table <- esr_analyzeDifferentialExpression(
    mat_t,
    validation_result$pheno,
    seed = 123
  )

  # 5. Gene selection
  selected_genes <- esr_selectTopGenes(
    de_table = de_table,
    n = 10,
    by = "de"
  )

  # 6. Visualize (all Mode 3 plots)
  p1 <- plotEndometrialPCA(mat_t, pheno = validation_result$pheno)
  p2 <- plotEndometrialMA(de_table)
  p3 <- plotEndometrialVolcano(de_table)
  p4 <- plotEndometrialHeatmap(
    mat_t,
    genes = selected_genes,
    pheno = validation_result$pheno
  )

  # 7. Create bundle
  bundle <- esr_createAnalysisBundle(
    counts_t = mat_t,
    de_table = de_table,
    selected_genes = selected_genes,
    qc_metrics = qc_metrics,
    pheno = validation_result$pheno,
    annot = validation_result$annot
  )

  # Verify all steps completed successfully
  expect_true(is.matrix(mat_t))
  expect_true(is.data.frame(de_table))
  expect_true(is.character(selected_genes))
  expect_true(is.data.frame(qc_metrics))
  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))
  expect_true(inherits(p3, "ggplot"))
  expect_true(inherits(p4, "Heatmap"))
  expect_true(is.list(bundle))
  expect_true("esr_analysis_bundle" %in% class(bundle))
})

test_that("bundle components are reproducible with fixed seed", {
  data(gse201926_sample)

  # Run workflow twice with same seed
  set.seed(123)
  mat_t1 <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table1 <- esr_analyzeDifferentialExpression(mat_t1, gse201926_sample$pheno, seed = 123)
  selected_genes1 <- esr_selectTopGenes(de_table = de_table1, n = 10, by = "de")
  qc_metrics1 <- esr_computeQCMetrics(counts = gse201926_sample$counts, mat_t = mat_t1)
  bundle1 <- esr_createAnalysisBundle(
    counts_t = mat_t1,
    de_table = de_table1,
    selected_genes = selected_genes1,
    qc_metrics = qc_metrics1
  )

  set.seed(123)
  mat_t2 <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table2 <- esr_analyzeDifferentialExpression(mat_t2, gse201926_sample$pheno, seed = 123)
  selected_genes2 <- esr_selectTopGenes(de_table = de_table2, n = 10, by = "de")
  qc_metrics2 <- esr_computeQCMetrics(counts = gse201926_sample$counts, mat_t = mat_t2)
  bundle2 <- esr_createAnalysisBundle(
    counts_t = mat_t2,
    de_table = de_table2,
    selected_genes = selected_genes2,
    qc_metrics = qc_metrics2
  )

  # Bundle components should be identical
  expect_equal(bundle1$counts_t, bundle2$counts_t)
  expect_equal(bundle1$de_table, bundle2$de_table)
  expect_equal(bundle1$selected_genes, bundle2$selected_genes)
  expect_equal(bundle1$qc_metrics, bundle2$qc_metrics)
})

test_that("bundle works with all plot functions", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
  selected_genes <- esr_selectTopGenes(de_table = de_table, n = 10, by = "de")
  qc_metrics <- esr_computeQCMetrics(
    counts = gse201926_sample$counts,
    mat_t = mat_t,
    pheno = gse201926_sample$pheno
  )

  bundle <- esr_createAnalysisBundle(
    counts_t = mat_t,
    de_table = de_table,
    selected_genes = selected_genes,
    qc_metrics = qc_metrics,
    pheno = gse201926_sample$pheno
  )

  # Test all plot functions work with bundle components
  expect_no_error({
    p1 <- plotEndometrialPCA(bundle$counts_t, pheno = bundle$pheno)
  })

  expect_no_error({
    p2 <- plotEndometrialMA(bundle$de_table)
  })

  expect_no_error({
    p3 <- plotEndometrialVolcano(bundle$de_table)
  })

  expect_no_error({
    p4 <- plotEndometrialHeatmap(
      bundle$counts_t,
      genes = bundle$selected_genes,
      pheno = bundle$pheno
    )
  })

  # Verify plots are created
  expect_true(inherits(p1, "ggplot"))
  expect_true(inherits(p2, "ggplot"))
  expect_true(inherits(p3, "ggplot"))
  expect_true(inherits(p4, "Heatmap"))
})

test_that("esr_createAnalysisBundle validates required inputs", {
  data(gse201926_sample)

  # Missing counts_t should error
  expect_error(
    {
      esr_createAnalysisBundle()
    },
    "counts_t is required"
  )
})

test_that("esr_createAnalysisBundle handles edge cases", {
  data(gse201926_sample)

  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)

  # Empty de_table should be handled gracefully (warnings expected for misalignment)
  de_table_empty <- data.frame(
    gene_id = character(),
    log2FC = numeric(),
    FDR = numeric(),
    stringsAsFactors = FALSE
  )

  expect_no_error({
    suppressWarnings({
      bundle <- esr_createAnalysisBundle(
        counts_t = mat_t,
        de_table = de_table_empty
      )
    })
  })

  # Empty selected_genes should be handled gracefully (warnings expected for misalignment)
  expect_no_error({
    suppressWarnings({
      bundle <- esr_createAnalysisBundle(
        counts_t = mat_t,
        selected_genes = character(0)
      )
    })
  })
})

# [END]
