test_that("esr_transform_log1p_cpm produces deterministic output", {
  data(gse201926_sample)
  
  # Set seed for reproducibility
  set.seed(123)
  result1 <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  set.seed(123)
  result2 <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  expect_equal(result1, result2)
})

test_that("esr_transform_log1p_cpm returns correct dimensions (samples x genes)", {
  data(gse201926_sample)
  
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  # Original: genes x samples (200 x 12)
  # Transformed: samples x genes (12 x genes_after_filtering)
  expect_equal(nrow(mat_t), ncol(gse201926_sample$counts))
  expect_lte(ncol(mat_t), nrow(gse201926_sample$counts))
})

test_that("esr_transform_log1p_cpm filters genes based on CPM threshold", {
  data(gse201926_sample)
  
  # Use strict filtering
  mat_t_strict <- esr_transform_log1p_cpm(
    gse201926_sample$counts,
    cpm_min = 10,
    cpm_min_samples = 6
  )
  
  # Use lenient filtering
  mat_t_lenient <- esr_transform_log1p_cpm(
    gse201926_sample$counts,
    cpm_min = 0.1,
    cpm_min_samples = 1
  )
  
  # Strict filtering should have fewer genes
  expect_lte(ncol(mat_t_strict), ncol(mat_t_lenient))
})

test_that("esr_transform_log1p_cpm handles edge cases", {
  data(gse201926_sample)
  
  # Very low counts
  X_low <- matrix(0, nrow = 10, ncol = 5)
  colnames(X_low) <- paste0("S", 1:5)
  rownames(X_low) <- paste0("G", 1:10)
  
  mat_t <- esr_transform_log1p_cpm(X_low, cpm_min = 1, cpm_min_samples = 1)
  
  expect_true(is.matrix(mat_t))
  expect_equal(nrow(mat_t), ncol(X_low))
})

test_that("esr_transform_log1p_cpm produces non-negative values", {
  data(gse201926_sample)
  
  mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
  
  expect_true(all(mat_t >= 0, na.rm = TRUE))
})

test_that("esr_transform_log1p_cpm handles data.frame input", {
  data(gse201926_sample)
  
  # Convert to data.frame
  X_df <- as.data.frame(gse201926_sample$counts)
  
  mat_t <- esr_transform_log1p_cpm(X_df)
  
  expect_true(is.matrix(mat_t))
  expect_equal(nrow(mat_t), ncol(X_df))
})

