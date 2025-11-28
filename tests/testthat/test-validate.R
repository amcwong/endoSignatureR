test_that("esr_validateEndometrial accepts valid inputs and returns expected structure", {
  data(gse201926_sample)

  result <- esr_validateEndometrial(
    gse201926_sample$counts,
    gse201926_sample$pheno,
    annot = gse201926_sample$annot
  )

  expect_type(result, "list")
  expect_true(all(c("X", "pheno", "annot", "issues") %in% names(result)))
  expect_true(is.matrix(result$X))
  expect_true(is.data.frame(result$pheno))
  expect_true(is.data.frame(result$issues))
  expect_true(all(c("type", "message") %in% names(result$issues)))
})

test_that("esr_validateEndometrial detects mismatched sample IDs", {
  data(gse201926_sample)

  # Create mismatched pheno
  pheno_mismatched <- gse201926_sample$pheno
  pheno_mismatched$sample_id[1] <- "MISMATCHED_SAMPLE"

  result <- esr_validateEndometrial(
    gse201926_sample$counts,
    pheno_mismatched
  )

  # Should have warnings about mismatched IDs
  mismatch_issues <- result$issues[grepl("not found", result$issues$message), ]
  expect_gt(nrow(mismatch_issues), 0)
})

test_that("esr_validateEndometrial detects class imbalance", {
  data(gse201926_sample)

  # Create imbalanced pheno
  pheno_imbalanced <- gse201926_sample$pheno
  pheno_imbalanced$group <- c(rep("PIS", 2), rep("PS", 10))

  result <- esr_validateEndometrial(
    gse201926_sample$counts,
    pheno_imbalanced
  )

  # Should have warning about class imbalance
  imbalance_issues <- result$issues[grepl("imbalance", result$issues$message), ]
  expect_gt(nrow(imbalance_issues), 0)
})

test_that("esr_validateEndometrial handles missing annot gracefully", {
  data(gse201926_sample)

  result <- esr_validateEndometrial(
    gse201926_sample$counts,
    gse201926_sample$pheno,
    annot = NULL
  )

  expect_true(is.null(result$annot))
  expect_true(is.data.frame(result$issues))
})

test_that("esr_validateEndometrial validates required columns", {
  data(gse201926_sample)

  # Missing sample_id
  pheno_no_id <- gse201926_sample$pheno[, "group", drop = FALSE]

  result <- esr_validateEndometrial(
    gse201926_sample$counts,
    pheno_no_id
  )

  error_issues <- result$issues[result$issues$type == "error", ]
  expect_gt(nrow(error_issues), 0)
})

test_that("esr_validateEndometrial handles non-numeric X with error", {
  # Create non-numeric matrix
  X_char <- matrix(LETTERS[1:12], nrow = 2, ncol = 6)
  colnames(X_char) <- paste0("S", 1:6)

  pheno_test <- data.frame(
    sample_id = colnames(X_char),
    group = rep(c("PS", "PIS"), each = 3),
    stringsAsFactors = FALSE
  )

  result <- esr_validateEndometrial(X_char, pheno_test)

  error_issues <- result$issues[result$issues$type == "error", ]
  expect_gt(nrow(error_issues), 0)
})

# [END]
