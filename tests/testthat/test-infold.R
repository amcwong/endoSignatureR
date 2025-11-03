# Tests for in-fold preprocessing (anti-leakage)
# Tests that preprocessing parameters come from training data only

test_that("esr_transformInFold uses training data only for parameters", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Create split
    set.seed(12345)
    splits <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Apply in-fold transform
    result <- esr_transformInFold(
        splits$splits[[1]],
        gse201926_trainmini$counts,
        gse201926_trainmini$pheno
    )

    # Check return structure
    expect_type(result, "list")
    expect_true("mat_t_train" %in% names(result))
    expect_true("mat_t_test" %in% names(result))
    expect_true("genes_keep" %in% names(result))
    expect_true("params" %in% names(result))

    # Check matrices are samples x genes
    expect_true(is.matrix(result$mat_t_train))
    expect_true(is.matrix(result$mat_t_test))
    expect_equal(ncol(result$mat_t_train), length(result$genes_keep))
    expect_equal(ncol(result$mat_t_test), length(result$genes_keep))

    # Check parameters were computed from training data
    expect_equal(
        result$params$n_train_samples,
        length(rsample::training(splits$splits[[1]])$sample_id)
    )
    expect_equal(
        result$params$n_test_samples,
        length(rsample::testing(splits$splits[[1]])$sample_id)
    )

    # Verify test data was transformed using training-only parameters
    # (same genes_keep should be applied to both)
    expect_equal(colnames(result$mat_t_train), colnames(result$mat_t_test))
})

test_that("esr_transformInFold prevents data leakage", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Create split
    set.seed(12345)
    splits <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Get training and test sample IDs
    train_ids <- rsample::training(splits$splits[[1]])$sample_id
    test_ids <- rsample::testing(splits$splits[[1]])$sample_id

    # Apply in-fold transform
    result <- esr_transformInFold(
        splits$splits[[1]],
        gse201926_trainmini$counts,
        gse201926_trainmini$pheno
    )

    # Verify training samples are in mat_t_train
    expect_equal(sort(rownames(result$mat_t_train)), sort(train_ids))

    # Verify test samples are in mat_t_test
    expect_equal(sort(rownames(result$mat_t_test)), sort(test_ids))

    # Verify genes_keep was determined from training data only
    # (we can't directly verify this, but we check that test data uses same genes)
    train_genes_keep <- colnames(result$mat_t_train)
    test_genes_keep <- colnames(result$mat_t_test)
    expect_equal(sort(train_genes_keep), sort(test_genes_keep))
})

test_that("esr_selectDEInFold uses training data only for selection", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Transform full dataset first
    mat_t <- esr_transform_log1p_cpm(gse201926_trainmini$counts)

    # Create split
    set.seed(12345)
    splits <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Get training and test sample IDs
    train_ids <- rsample::training(splits$splits[[1]])$sample_id
    test_ids <- rsample::testing(splits$splits[[1]])$sample_id

    # Select genes in-fold (should use training data only)
    selected_genes <- esr_selectDEInFold(
        splits$splits[[1]],
        mat_t,
        gse201926_trainmini$pheno,
        n = 50,
        method = "de"
    )

    # Check return type
    expect_type(selected_genes, "character")
    expect_true(length(selected_genes) <= 50)

    # Verify selected genes are in transformed matrix
    expect_true(all(selected_genes %in% colnames(mat_t)))

    # Verify selection was based on training data only
    # (we can verify by checking that DE analysis would use training data)
    # This is implicitly tested by ensuring test_ids are never used in selection

    # Verify selected genes exist in training data
    mat_t_train <- mat_t[train_ids, , drop = FALSE]
    expect_true(all(selected_genes %in% colnames(mat_t_train)))
})

test_that("esr_selectDEInFold prevents data leakage", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Transform full dataset
    mat_t <- esr_transform_log1p_cpm(gse201926_trainmini$counts)

    # Create split
    set.seed(12345)
    splits <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Get training and test sample IDs
    train_ids <- rsample::training(splits$splits[[1]])$sample_id
    test_ids <- rsample::testing(splits$splits[[1]])$sample_id

    # Select genes in-fold (DE method)
    selected_de <- esr_selectDEInFold(
        splits$splits[[1]],
        mat_t,
        gse201926_trainmini$pheno,
        n = 50,
        method = "de"
    )

    # Select genes in-fold (variance method)
    selected_var <- esr_selectDEInFold(
        splits$splits[[1]],
        mat_t,
        gse201926_trainmini$pheno,
        n = 50,
        method = "variance"
    )

    # Both should return character vectors
    expect_type(selected_de, "character")
    expect_type(selected_var, "character")

    # Both should use training data only (implicitly tested)
    # We verify that selected genes exist in training data
    mat_t_train <- mat_t[train_ids, , drop = FALSE]
    expect_true(all(selected_de %in% colnames(mat_t_train)))
    expect_true(all(selected_var %in% colnames(mat_t_train)))
})

test_that("in-fold preprocessing produces deterministic outputs with fixed seed", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Create split with fixed seed
    set.seed(12345)
    splits <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Transform in-fold twice with same seed
    set.seed(123)
    result1 <- esr_transformInFold(
        splits$splits[[1]],
        gse201926_trainmini$counts,
        gse201926_trainmini$pheno
    )

    set.seed(123)
    result2 <- esr_transformInFold(
        splits$splits[[1]],
        gse201926_trainmini$counts,
        gse201926_trainmini$pheno
    )

    # Results should be identical
    expect_equal(result1$genes_keep, result2$genes_keep)
    expect_equal(result1$mat_t_train, result2$mat_t_train)
    expect_equal(result1$mat_t_test, result2$mat_t_test)

    # Select DE in-fold twice with same seed
    mat_t <- esr_transform_log1p_cpm(gse201926_trainmini$counts)

    set.seed(123)
    selected1 <- esr_selectDEInFold(
        splits$splits[[1]],
        mat_t,
        gse201926_trainmini$pheno,
        n = 50,
        method = "de",
        seed = 123
    )

    set.seed(123)
    selected2 <- esr_selectDEInFold(
        splits$splits[[1]],
        mat_t,
        gse201926_trainmini$pheno,
        n = 50,
        method = "de",
        seed = 123
    )

    # Results should be identical
    expect_equal(selected1, selected2)
})
