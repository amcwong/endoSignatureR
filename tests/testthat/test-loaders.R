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

# Tests for user data loading functions
test_that("esr_loadCountsFromFile loads counts correctly", {
  counts_file <- system.file("extdata", "gse201926_raw_counts.tsv.gz", package = "endoSignatureR")
  skip_if_not(file.exists(counts_file), "Test data file not found")

  counts <- esr_loadCountsFromFile(counts_file, gene_id_col = "GeneID")
  expect_true(is.matrix(counts))
  expect_true(nrow(counts) > 0)
  expect_true(ncol(counts) > 0)
  expect_true(!is.null(rownames(counts)))
  expect_true(!is.null(colnames(counts)))
  expect_true(all(is.numeric(counts)))
})

test_that("esr_loadCountsFromFile handles missing file", {
  expect_error(
    esr_loadCountsFromFile("nonexistent_file.tsv"),
    "Counts file not found"
  )
})

test_that("esr_loadCountsFromFile handles missing gene_id_col", {
  counts_file <- system.file("extdata", "gse201926_raw_counts.tsv.gz", package = "endoSignatureR")
  skip_if_not(file.exists(counts_file), "Test data file not found")

  expect_error(
    esr_loadCountsFromFile(counts_file, gene_id_col = "NonExistentCol"),
    "Gene ID column"
  )
})

test_that("esr_loadPhenoFromFile loads phenotype from series matrix", {
  pheno_file <- system.file("extdata", "gse201926_series_matrix.txt.gz", package = "endoSignatureR")
  skip_if_not(file.exists(pheno_file), "Test data file not found")

  pheno <- esr_loadPhenoFromFile(pheno_file, format = "series_matrix")
  expect_true(is.data.frame(pheno))
  expect_true("sample_id" %in% names(pheno))
  expect_true("group" %in% names(pheno))
  expect_true(nrow(pheno) > 0)
})

test_that("esr_loadPhenoFromFile handles missing file", {
  expect_error(
    esr_loadPhenoFromFile("nonexistent_file.tsv"),
    "Phenotype file not found"
  )
})

test_that("esr_loadAnnotFromFile loads annotation correctly", {
  annot_file <- system.file("extdata", "gse201926_annotation.tsv.gz", package = "endoSignatureR")
  skip_if_not(file.exists(annot_file), "Test data file not found")

  annot <- esr_loadAnnotFromFile(annot_file, gene_id_col = "GeneID")
  expect_true(is.data.frame(annot))
  expect_true("GeneID" %in% names(annot))
  expect_true(nrow(annot) > 0)
})

test_that("esr_loadAnnotFromFile handles missing gene_id_col", {
  annot_file <- system.file("extdata", "gse201926_annotation.tsv.gz", package = "endoSignatureR")
  skip_if_not(file.exists(annot_file), "Test data file not found")

  expect_error(
    esr_loadAnnotFromFile(annot_file, gene_id_col = "NonExistentCol"),
    "Gene ID column"
  )
})

test_that("esr_loadFromFiles loads all files correctly", {
  counts_file <- system.file("extdata", "gse201926_raw_counts.tsv.gz", package = "endoSignatureR")
  pheno_file <- system.file("extdata", "gse201926_series_matrix.txt.gz", package = "endoSignatureR")
  annot_file <- system.file("extdata", "gse201926_annotation.tsv.gz", package = "endoSignatureR")
  skip_if_not(all(file.exists(c(counts_file, pheno_file, annot_file))), "Test data files not found")

  user_data <- esr_loadFromFiles(
    counts_file = counts_file,
    pheno_file = pheno_file,
    annot_file = annot_file,
    pheno_format = "series_matrix",
    validate = FALSE
  )

  expect_type(user_data, "list")
  expect_true("counts" %in% names(user_data))
  expect_true("pheno" %in% names(user_data))
  expect_true("annot" %in% names(user_data))
  expect_true(is.matrix(user_data$counts))
  expect_true(is.data.frame(user_data$pheno))
  expect_true(is.data.frame(user_data$annot))
})

test_that("esr_loadFromFiles handles unlabeled data (Mode 1)", {
  counts_file <- system.file("extdata", "gse201926_raw_counts.tsv.gz", package = "endoSignatureR")
  skip_if_not(file.exists(counts_file), "Test data file not found")

  user_data <- esr_loadFromFiles(
    counts_file = counts_file,
    pheno_file = NULL,
    annot_file = NULL,
    validate = FALSE
  )

  expect_type(user_data, "list")
  expect_true("counts" %in% names(user_data))
  expect_true(is.null(user_data$pheno))
  expect_true(is.null(user_data$annot))
})

test_that("esr_loadFromFiles aligns IDs correctly", {
  counts_file <- system.file("extdata", "gse201926_raw_counts.tsv.gz", package = "endoSignatureR")
  pheno_file <- system.file("extdata", "gse201926_series_matrix.txt.gz", package = "endoSignatureR")
  annot_file <- system.file("extdata", "gse201926_annotation.tsv.gz", package = "endoSignatureR")
  skip_if_not(all(file.exists(c(counts_file, pheno_file, annot_file))), "Test data files not found")

  user_data <- esr_loadFromFiles(
    counts_file = counts_file,
    pheno_file = pheno_file,
    annot_file = annot_file,
    pheno_format = "series_matrix",
    align_ids = TRUE,
    validate = FALSE
  )

  # Check that sample IDs align between counts and pheno
  if (!is.null(user_data$pheno)) {
    expect_equal(colnames(user_data$counts), user_data$pheno$sample_id)
  }

  # Check that gene IDs align between counts and annot
  if (!is.null(user_data$annot)) {
    expect_equal(rownames(user_data$counts), user_data$annot$GeneID)
  }
})

test_that("esr_loadFromFiles validates data when validate=TRUE", {
  counts_file <- system.file("extdata", "gse201926_raw_counts.tsv.gz", package = "endoSignatureR")
  pheno_file <- system.file("extdata", "gse201926_series_matrix.txt.gz", package = "endoSignatureR")
  annot_file <- system.file("extdata", "gse201926_annotation.tsv.gz", package = "endoSignatureR")
  skip_if_not(all(file.exists(c(counts_file, pheno_file, annot_file))), "Test data files not found")

  user_data <- esr_loadFromFiles(
    counts_file = counts_file,
    pheno_file = pheno_file,
    annot_file = annot_file,
    pheno_format = "series_matrix",
    validate = TRUE
  )

  expect_type(user_data, "list")
  expect_true("issues" %in% names(user_data))
  expect_true(is.data.frame(user_data$issues))
})
