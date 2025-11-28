# Tests for CV split functionality
# Tests split determinism, structure, balance, and contracts

test_that("folds_demo loads correctly and has expected structure", {
    # This test will fail until folds_demo is created
    # Skip for now if dataset doesn't exist
    skip_if_not(exists("folds_demo", where = asNamespace("endoSignatureR"), mode = "any"))

    data(folds_demo)

    expect_type(folds_demo, "list")
    expect_true("outer_splits" %in% names(folds_demo))
    expect_true("inner_splits" %in% names(folds_demo))
    expect_true("outer_seed" %in% names(folds_demo))
    expect_true("inner_seed" %in% names(folds_demo))
    expect_true("n_outer_folds" %in% names(folds_demo))
    expect_type(folds_demo$outer_seed, "double") # R stores seeds as numeric (double)
    expect_type(folds_demo$inner_seed, "double") # R stores seeds as numeric (double)
})

test_that("splits are deterministic with fixed seed", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Create splits with same seed twice
    set.seed(12345)
    splits1 <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    set.seed(12345)
    splits2 <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Splits should be identical
    expect_equal(nrow(splits1), nrow(splits2))

    # Compare training/test indices for first fold
    train1 <- rsample::training(splits1$splits[[1]])$sample_id
    train2 <- rsample::training(splits2$splits[[1]])$sample_id
    expect_equal(sort(train1), sort(train2))

    test1 <- rsample::testing(splits1$splits[[1]])$sample_id
    test2 <- rsample::testing(splits2$splits[[1]])$sample_id
    expect_equal(sort(test1), sort(test2))
})

test_that("split artifacts maintain class balance", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Create stratified splits
    set.seed(12345)
    splits <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Check class balance in each fold
    for (i in seq_len(nrow(splits))) {
        train_data <- rsample::training(splits$splits[[i]])
        test_data <- rsample::testing(splits$splits[[i]])

        # Both training and test should have both classes
        expect_true(length(unique(train_data$group)) >= 2)
        expect_true(length(unique(test_data$group)) >= 2)

        # Training set should not have all samples from one class
        train_counts <- table(train_data$group)
        expect_false(any(train_counts == 0))

        # Test set should not have all samples from one class
        test_counts <- table(test_data$group)
        expect_false(any(test_counts == 0))
    }
})

test_that("split contracts are valid (can extract training/test indices)", {
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))

    library(rsample)
    data(gse201926_trainmini)

    # Create splits
    set.seed(12345)
    splits <- rsample::vfold_cv(
        data = gse201926_trainmini$pheno,
        v = 3,
        strata = "group",
        breaks = 2
    )

    # Test that we can extract training/test indices
    for (i in seq_len(nrow(splits))) {
        train_data <- rsample::training(splits$splits[[i]])
        test_data <- rsample::testing(splits$splits[[i]])

        # Should return data frames
        expect_true(is.data.frame(train_data))
        expect_true(is.data.frame(test_data))

        # Should have sample_id column
        expect_true("sample_id" %in% names(train_data))
        expect_true("sample_id" %in% names(test_data))

        # Training and test should be disjoint
        expect_equal(length(intersect(train_data$sample_id, test_data$sample_id)), 0)

        # Training and test should cover all samples
        all_samples <- union(train_data$sample_id, test_data$sample_id)
        expect_equal(sort(all_samples), sort(gse201926_trainmini$pheno$sample_id))
    }
})

test_that("folds_demo outer_splits structure is valid", {
    skip_if_not(exists("folds_demo", where = asNamespace("endoSignatureR"), mode = "any"))
    skip_if_not_installed("rsample")

    library(rsample)
    data(folds_demo)

    expect_true(inherits(folds_demo$outer_splits, "vfold_cv"))
    expect_true(nrow(folds_demo$outer_splits) > 0)

    # Check first outer split
    train_data <- rsample::training(folds_demo$outer_splits$splits[[1]])
    test_data <- rsample::testing(folds_demo$outer_splits$splits[[1]])

    expect_true(is.data.frame(train_data))
    expect_true(is.data.frame(test_data))
    expect_equal(length(intersect(train_data$sample_id, test_data$sample_id)), 0)
})

test_that("folds_demo inner_splits structure is valid", {
    skip_if_not(exists("folds_demo", where = asNamespace("endoSignatureR"), mode = "any"))
    skip_if_not_installed("rsample")

    library(rsample)
    data(folds_demo)

    # Check inner splits exist
    expect_true(is.list(folds_demo$inner_splits))
    expect_equal(length(folds_demo$inner_splits), nrow(folds_demo$outer_splits))

    # Check first inner split if available
    if (!is.null(folds_demo$inner_splits[[1]])) {
        expect_true(inherits(folds_demo$inner_splits[[1]], "vfold_cv"))
        expect_true(nrow(folds_demo$inner_splits[[1]]) > 0)

        # Inner split should use training data from outer split (anti-leakage)
        outer_train <- rsample::training(folds_demo$outer_splits$splits[[1]])$sample_id
        inner_train <- rsample::training(folds_demo$inner_splits[[1]]$splits[[1]])$sample_id

        # All inner training samples should be in outer training set
        expect_true(all(inner_train %in% outer_train))
    }
})

# [END]
