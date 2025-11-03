# Test Export Functionality

test_that("esr_exportSignature exports CSV with correct structure", {
  # Load test data
  data(gse201926_trainmini)
  data(folds_demo)

  # Train signature
  # Suppress warnings about no consensus genes (expected for small sample sizes)
  result <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 100,
    outer_folds = NULL, # Use folds_demo
    seed = 123
  ))

  # Export to temp directory
  export_dir <- tempfile()
  dir.create(export_dir)
  on.exit(unlink(export_dir, recursive = TRUE), add = TRUE)

  # Export signature
  paths <- esr_exportSignature(
    signature = result$signature,
    result = result,
    dir = export_dir,
    formats = "csv"
  )

  # Verify CSV file exists
  csv_path <- file.path(export_dir, "endometrial_signature.csv")
  expect_true(file.exists(csv_path))

  # Read CSV back
  df <- readr::read_csv(csv_path, show_col_types = FALSE)

  # Verify expected columns
  expect_true("gene_id" %in% names(df))
  expect_true("coefficient" %in% names(df))
  expect_true("selection_frequency" %in% names(df))

  # Verify gene IDs match signature panel
  expect_equal(sort(df$gene_id[df$gene_id != "(Intercept)"]), sort(result$signature$panel))

  # Verify coefficients match (for non-intercept rows)
  # Remove names from data.frame column values for comparison
  sig_coefs <- result$signature$coefficients[match(df$gene_id[df$gene_id != "(Intercept)"], names(result$signature$coefficients))]
  expect_equal(as.numeric(df$coefficient[df$gene_id != "(Intercept)"]), as.numeric(sig_coefs), tolerance = 1e-6)

  # Verify intercept row exists if requested
  expect_true("(Intercept)" %in% df$gene_id)
  # Extract intercept value and compare (remove names)
  intercept_from_csv <- as.numeric(df$coefficient[df$gene_id == "(Intercept)"])
  intercept_from_sig <- as.numeric(result$signature$intercept)
  expect_equal(intercept_from_csv, intercept_from_sig, tolerance = 1e-6)

  # Verify selection frequencies match
  # Remove names from data.frame column values for comparison
  sig_freqs <- result$signature$selection_frequency[match(df$gene_id[df$gene_id != "(Intercept)"], names(result$signature$selection_frequency))]
  expect_equal(as.numeric(df$selection_frequency[df$gene_id != "(Intercept)"]), as.numeric(sig_freqs))
})

test_that("esr_exportSignature exports JSON with correct schema", {
  # Load test data
  data(gse201926_trainmini)
  data(folds_demo)

  # Train signature
  # Suppress warnings about no consensus genes (expected for small sample sizes)
  result <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 100,
    outer_folds = NULL, # Use folds_demo
    seed = 123
  ))

  # Export to temp directory
  export_dir <- tempfile()
  dir.create(export_dir)
  on.exit(unlink(export_dir, recursive = TRUE), add = TRUE)

  # Export recipe
  paths <- esr_exportSignature(
    signature = result$signature,
    result = result,
    dir = export_dir,
    formats = "json"
  )

  # Verify JSON file exists
  json_path <- file.path(export_dir, "endometrial_recipe.json")
  expect_true(file.exists(json_path))

  # Read JSON back
  json_data <- jsonlite::fromJSON(json_path)

  # Verify expected structure
  expect_true("preprocessing" %in% names(json_data))
  expect_true("training" %in% names(json_data))
  expect_true("signature" %in% names(json_data))
  expect_true("reproducibility" %in% names(json_data))

  # Verify preprocessing parameters
  expect_true("transform" %in% names(json_data$preprocessing))
  expect_true("cpm_min" %in% names(json_data$preprocessing))
  expect_true("top_k" %in% names(json_data$preprocessing))

  # Verify training parameters
  expect_true("lambda_rule" %in% names(json_data$training))
  expect_true("aggregation_method" %in% names(json_data$training))

  # Verify signature metadata
  expect_true("n_genes" %in% names(json_data$signature))
  expect_equal(json_data$signature$n_genes, length(result$signature$panel))
  # Compare intercept values (JSON stores as numeric, signature might have names)
  # Use larger tolerance for JSON comparisons due to JSON serialization precision loss
  # JSON format limits precision - values are rounded to ~5-6 significant digits
  # For intercept values around 0.05-0.1, absolute tolerance of 1e-3 (0.001) is appropriate
  intercept_json <- as.numeric(json_data$signature$intercept)
  intercept_sig <- as.numeric(result$signature$intercept)
  expect_equal(intercept_json, intercept_sig, tolerance = 1e-3)

  # Verify reproducibility info
  expect_true("package_version" %in% names(json_data$reproducibility))
  expect_true("r_version" %in% names(json_data$reproducibility))
})

test_that("esr_exportSignature exports model card with required sections", {
  # Load test data
  data(gse201926_trainmini)
  data(folds_demo)

  # Train signature
  # Suppress warnings about no consensus genes (expected for small sample sizes)
  result <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 100,
    outer_folds = NULL, # Use folds_demo
    seed = 123
  ))

  # Export to temp directory
  export_dir <- tempfile()
  dir.create(export_dir)
  on.exit(unlink(export_dir, recursive = TRUE), add = TRUE)

  # Export model card
  paths <- esr_exportSignature(
    signature = result$signature,
    result = result,
    dir = export_dir,
    formats = "md"
  )

  # Verify Markdown file exists
  md_path <- file.path(export_dir, "endometrial_model_card.md")
  expect_true(file.exists(md_path))

  # Read Markdown back
  md_lines <- readLines(md_path)

  # Verify required sections are present
  expect_true(any(grepl("^# Endometrial Signature Model Card$", md_lines)))
  expect_true(any(grepl("^## Model Details$", md_lines)))
  expect_true(any(grepl("^## Model Purpose$", md_lines)))
  expect_true(any(grepl("^## Training Data$", md_lines)))
  expect_true(any(grepl("^## Model Architecture$", md_lines)))
  expect_true(any(grepl("^## Performance Metrics$", md_lines)))
  expect_true(any(grepl("^## Limitations$", md_lines)))
  expect_true(any(grepl("^## Intended Use$", md_lines)))

  # Verify performance metrics match signature metrics
  if (!is.null(result$metrics)) {
    auc_line <- md_lines[grepl("AUC-ROC", md_lines)]
    expect_true(length(auc_line) > 0)
    # Extract AUC value from line (format: "- **AUC-ROC**: 0.xxx")
    if (length(auc_line) > 0 && !grepl("N/A", auc_line)) {
      auc_value <- as.numeric(gsub(".*: (.*)$", "\\1", auc_line))
      expect_equal(auc_value, result$metrics$auc, tolerance = 1e-2)
    }
  }
})

test_that("esr_exportSignature handles standalone signature", {
  # Load test data
  data(gse201926_trainmini)
  data(folds_demo)

  # Train signature
  # Suppress warnings about no consensus genes (expected for small sample sizes)
  result <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 100,
    outer_folds = NULL, # Use folds_demo
    seed = 123
  ))

  # Extract standalone signature (no result)
  signature <- result$signature

  # Export to temp directory
  export_dir <- tempfile()
  dir.create(export_dir)
  on.exit(unlink(export_dir, recursive = TRUE), add = TRUE)

  # Export standalone signature (CSV and JSON should work; model card will be minimal)
  paths <- esr_exportSignature(
    signature = signature,
    result = NULL, # No result
    dir = export_dir,
    formats = c("csv", "json", "md")
  )

  # Verify CSV exists
  csv_path <- file.path(export_dir, "endometrial_signature.csv")
  expect_true(file.exists(csv_path))

  # Verify JSON exists
  json_path <- file.path(export_dir, "endometrial_recipe.json")
  expect_true(file.exists(json_path))

  # Verify model card exists (will be minimal without result)
  md_path <- file.path(export_dir, "endometrial_model_card.md")
  expect_true(file.exists(md_path))

  # Verify model card indicates metrics not available
  md_lines <- readLines(md_path)
  expect_true(any(grepl("signature-only export", md_lines, ignore.case = TRUE)))
})

test_that("esr_exportSignature validates inputs and handles errors", {
  # Test invalid signature structure (missing required fields)
  invalid_sig <- list(panel = c("GENE1", "GENE2"))
  expect_error(
    esr_exportSignature(signature = invalid_sig, dir = tempdir()),
    "signature must include 'coefficients'"
  )

  # Test invalid signature structure (not a list)
  expect_error(
    esr_exportSignature(signature = "not a list", dir = tempdir()),
    "signature must be a list"
  )

  # Test invalid formats
  valid_sig <- list(
    panel = c("GENE1", "GENE2"),
    coefficients = c(GENE1 = 0.5, GENE2 = -0.3),
    intercept = 0.1,
    selection_frequency = c(GENE1 = 3L, GENE2 = 2L)
  )
  expect_error(
    esr_exportSignature(signature = valid_sig, dir = tempdir(), formats = "invalid"),
    "Invalid formats"
  )

  # Test NULL signature
  expect_error(
    esr_exportSignature(signature = NULL, dir = tempdir()),
    "signature must be provided"
  )
})

test_that("export is reproducible with fixed seed", {
  # Load test data
  data(gse201926_trainmini)
  data(folds_demo)

  # Train signature with fixed seed
  # Suppress warnings about no consensus genes (expected for small sample sizes)
  result1 <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 100,
    outer_folds = NULL,
    seed = 123
  ))

  result2 <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 100,
    outer_folds = NULL,
    seed = 123
  ))

  # Export both signatures
  export_dir1 <- tempfile()
  export_dir2 <- tempfile()
  dir.create(export_dir1)
  dir.create(export_dir2)
  on.exit(unlink(c(export_dir1, export_dir2), recursive = TRUE), add = TRUE)

  paths1 <- esr_exportSignature(
    signature = result1$signature,
    result = result1,
    dir = export_dir1
  )

  paths2 <- esr_exportSignature(
    signature = result2$signature,
    result = result2,
    dir = export_dir2
  )

  # Verify exported files have same content (same signature = same export)
  csv1 <- readr::read_csv(file.path(export_dir1, "endometrial_signature.csv"), show_col_types = FALSE)
  csv2 <- readr::read_csv(file.path(export_dir2, "endometrial_signature.csv"), show_col_types = FALSE)
  # Compare data frames (ignoring row names)
  expect_equal(csv1$gene_id, csv2$gene_id)
  expect_equal(as.numeric(csv1$coefficient), as.numeric(csv2$coefficient), tolerance = 1e-6)
  expect_equal(as.numeric(csv1$selection_frequency), as.numeric(csv2$selection_frequency))

  json1 <- jsonlite::fromJSON(file.path(export_dir1, "endometrial_recipe.json"))
  json2 <- jsonlite::fromJSON(file.path(export_dir2, "endometrial_recipe.json"))
  expect_equal(json1, json2)
})

test_that("esr_exportSignature includes bootstrap frequencies when available", {
  # Load test data
  data(gse201926_trainmini)
  data(folds_demo)

  # Train signature with stability selection
  # Suppress warnings about no consensus genes (expected for small sample sizes)
  result <- suppressWarnings(esr_trainEndometrialSignature(
    X = gse201926_trainmini$counts,
    pheno = gse201926_trainmini$pheno,
    top_k = 100,
    outer_folds = NULL,
    stability_selection = TRUE,
    stability_resamples = 50,
    seed = 123
  ))

  # Export to temp directory
  export_dir <- tempfile()
  dir.create(export_dir)
  on.exit(unlink(export_dir, recursive = TRUE), add = TRUE)

  # Export signature with stability
  paths <- esr_exportSignature(
    signature = result$signature,
    result = result,
    dir = export_dir,
    formats = "csv",
    include_stability = TRUE
  )

  # Read CSV back
  csv_path <- file.path(export_dir, "endometrial_signature.csv")
  df <- readr::read_csv(csv_path, show_col_types = FALSE)

  # Verify bootstrap_frequency column exists if stability was computed
  if (!is.null(result$stability$bootstrap_frequency) && length(result$stability$bootstrap_frequency) > 0) {
    expect_true("bootstrap_frequency" %in% names(df))
  # Verify bootstrap frequencies match (for genes in panel)
  # Remove names from data.frame column values for comparison
  sig_bootstrap <- result$stability$bootstrap_frequency[match(df$gene_id[df$gene_id != "(Intercept)"], names(result$stability$bootstrap_frequency))]
  expect_equal(as.numeric(df$bootstrap_frequency[df$gene_id != "(Intercept)"]), as.numeric(sig_bootstrap), tolerance = 1e-6)
  }
})

