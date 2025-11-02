## code to prepare `folds_demo` dataset
## Demo CV splits with fixed seeds for reproducibility
## Used in vignettes and tests to avoid recomputing splits

suppressPackageStartupMessages({
    library(rsample)
    library(usethis)
    library(dplyr)
})

# Note: This script should be run AFTER gse201926_trainmini is created
# We'll load the training mini dataset to create splits

# Check if gse201926_trainmini exists
if (!exists("gse201926_trainmini")) {
    # Try to load it
    if (file.exists("data/gse201926_trainmini.rda")) {
        load("data/gse201926_trainmini.rda")
    } else {
        # Create it first by sourcing the build script
        message("gse201926_trainmini not found. Creating it first...")
        source("data-raw/build_gse201926_trainmini.R")
    }
}

# Extract phenotype for creating splits
pheno <- gse201926_trainmini$pheno

# Verify we have balanced classes
group_counts <- table(pheno$group)
message("Class distribution for splits: ", paste(names(group_counts), group_counts, sep = "=", collapse = ", "))

# Set fixed seeds for reproducibility
outer_seed <- 12345
inner_seed <- 67890

message("Creating outer CV splits with seed ", outer_seed, "...")

# For small dataset (12 samples), use Leave-Pair-Out (LPO) or stratified K-fold
# Since we have 6 PS + 6 PIS, we can use stratified K-fold with k=3 or k=4
# This gives us outer splits for cross-validation

# Create outer splits (K-fold with stratification)
# Use k=3 for 12 samples (3 folds × 4 samples each)
set.seed(outer_seed)
outer_splits <- rsample::vfold_cv(
    data = pheno,
    v = 3,
    strata = "group",
    breaks = 2  # Minimal strata breaks for small dataset
)

message("Created ", nrow(outer_splits), " outer CV splits")

# Create inner splits for each outer fold
# Inner splits will be used for hyperparameter tuning within each outer fold
message("Creating inner CV splits for each outer fold with seed ", inner_seed, "...")

inner_splits_list <- lapply(1:nrow(outer_splits), function(i) {
    # Get training data for this outer fold
    outer_train <- rsample::training(outer_splits$splits[[i]])
    outer_train_pheno <- outer_train
    
    # Create inner splits on training data only (anti-leakage)
    set.seed(inner_seed + i)  # Different seed for each outer fold
    
    # Use v=2 or v=3 for inner folds (smaller than outer since training set is smaller)
    # For 8 samples (after leaving 4 out), use v=2 (each fold has 4 samples)
    if (nrow(outer_train) >= 4) {
        inner_splits <- rsample::vfold_cv(
            data = outer_train_pheno,
            v = 2,
            strata = "group",
            breaks = 2
        )
    } else {
        # If too small, use LOOCV or skip inner splits
        inner_splits <- NULL
    }
    
    return(inner_splits)
})

# Store split artifacts with metadata
folds_demo <- list(
    outer_splits = outer_splits,
    inner_splits = inner_splits_list,
    outer_seed = as.integer(outer_seed),  # Ensure integer type
    inner_seed = as.integer(inner_seed),  # Ensure integer type
    n_outer_folds = nrow(outer_splits),
    n_inner_folds = if (!is.null(inner_splits_list[[1]])) nrow(inner_splits_list[[1]]) else NULL,
    split_type = "vfold_cv",
    strata = "group",
    created_date = Sys.Date()
)

message("Created folds_demo with:")
message("  - ", folds_demo$n_outer_folds, " outer CV splits")
message("  - ", if (!is.null(folds_demo$n_inner_folds)) folds_demo$n_inner_folds else "N/A", " inner CV splits per outer fold")
message("  - Outer seed: ", folds_demo$outer_seed)
message("  - Inner seed: ", folds_demo$inner_seed)

# Save as package data
usethis::use_data(folds_demo, overwrite = TRUE, compress = "xz")
message("Saved folds_demo to data/")

