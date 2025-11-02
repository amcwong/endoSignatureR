#' Train Endometrial Signature (Placeholder)
#'
#' Trains a PS vs PIS signature using a defined pipeline.
#'
#' @param X Matrix/data.frame of counts (genes x samples).
#' @param pheno Data.frame with labels.
#' @param ... Additional arguments (ignored in placeholder).
#'
#' @return A list containing a placeholder signature and metrics.
#' @export
esr_trainEndometrialSignature <- function(X, pheno, ...) {
  signature <- list(panel = character(), coefficients = numeric(), recipe = list())
  metrics <- list(auc = NA_real_)
  return(list(signature = signature, metrics = metrics))
}

#' Transform Counts In-Fold (Anti-leakage)
#'
#' Applies transform/filter within CV fold, computing parameters from training data only
#' and applying them to test data to prevent data leakage.
#'
#' @param split An `rsample` split object (e.g., from `rsample::vfold_cv()`).
#' @param counts Matrix/data.frame of raw counts (genes x samples) for full dataset.
#' @param pheno Data.frame with sample metadata; must include `sample_id`.
#' @param transform Character scalar; transformation method. Defaults to "log1p-cpm".
#' @param cpm_min Minimum CPM threshold for gene filtering. Defaults to 1.
#' @param cpm_min_samples Minimum number of samples that must meet CPM threshold. Defaults to 4.
#'
#' @return A list with elements:
#' \describe{
#'   \item{mat_t_train}{Transformed matrix (samples x genes) for training data}
#'   \item{mat_t_test}{Transformed matrix (samples x genes) for test data}
#'   \item{genes_keep}{Character vector of gene IDs kept after filtering}
#'   \item{params}{List with preprocessing parameters used (for reproducibility)}
#' }
#'
#' @details
#' This function enforces anti-leakage by:
#' - Computing CPM filtering parameters from training data only
#' - Applying training-based parameters to test data
#' - Ensuring test data never influences preprocessing decisions
#'
#' @examples
#' \dontrun{
#' library(rsample)
#' data(gse201926_trainmini)
#' set.seed(123)
#' splits <- vfold_cv(gse201926_trainmini$pheno, v = 3)
#' result <- esr_transformInFold(splits$splits[[1]], 
#'                               gse201926_trainmini$counts,
#'                               gse201926_trainmini$pheno)
#' dim(result$mat_t_train)
#' dim(result$mat_t_test)
#' }
#' @export
esr_transformInFold <- function(split, counts, pheno, 
                                transform = "log1p-cpm",
                                cpm_min = 1, 
                                cpm_min_samples = 4) {
    if (!requireNamespace("rsample", quietly = TRUE)) {
        stop("rsample package is required for in-fold preprocessing")
    }
    
    # Extract training and test sample IDs from split
    train_data <- rsample::training(split)
    test_data <- rsample::testing(split)
    train_sample_ids <- train_data$sample_id
    test_sample_ids <- test_data$sample_id
    
    # Match sample IDs to counts matrix column names
    counts_colnames <- colnames(counts)
    if (is.null(counts_colnames)) {
        stop("counts matrix must have column names (sample IDs)")
    }
    
    # Diagnostic information for debugging
    if (any(is.na(match(train_sample_ids, counts_colnames)))) {
        missing_train <- train_sample_ids[is.na(match(train_sample_ids, counts_colnames))]
        stop(
            "Some training sample IDs not found in counts column names.\n",
            "Training sample IDs from split: ", paste(head(train_sample_ids, 5), collapse = ", "), "...\n",
            "Counts column names: ", paste(head(counts_colnames, 5), collapse = ", "), "...\n",
            "Missing training IDs: ", paste(missing_train, collapse = ", "), "\n",
            "All training IDs: ", paste(train_sample_ids, collapse = ", "), "\n",
            "All counts columns: ", paste(counts_colnames, collapse = ", ")
        )
    }
    
    train_indices <- match(train_sample_ids, counts_colnames)
    test_indices <- match(test_sample_ids, counts_colnames)
    
    if (any(is.na(train_indices))) {
        missing_train <- train_sample_ids[is.na(train_indices)]
        stop(
            "Some training sample IDs not found in counts column names.\n",
            "Missing: ", paste(missing_train, collapse = ", "), "\n",
            "Available: ", paste(head(counts_colnames, 10), collapse = ", ")
        )
    }
    if (any(is.na(test_indices))) {
        missing_test <- test_sample_ids[is.na(test_indices)]
        stop(
            "Some test sample IDs not found in counts column names.\n",
            "Missing: ", paste(missing_test, collapse = ", "), "\n",
            "Available: ", paste(head(counts_colnames, 10), collapse = ", ")
        )
    }
    
    # Get training and test counts (genes x samples)
    counts_train <- counts[, train_indices, drop = FALSE]
    counts_test <- counts[, test_indices, drop = FALSE]
    
    # Transform training data (computes parameters)
    if (transform == "log1p-cpm") {
        # Apply transform on training data to compute parameters
        mat_t_train <- esr_transform_log1p_cpm(counts_train, 
                                                cpm_min = cpm_min, 
                                                cpm_min_samples = cpm_min_samples)
        
        # Get genes kept from training transform
        # We need to determine which genes were kept by comparing dimensions
        # The transform function filters internally, so we need to check
        # We'll compute CPM on training data to determine kept genes
        lib_sizes_train <- colSums(counts_train, na.rm = TRUE)
        cpm_train <- t(t(counts_train) / lib_sizes_train * 1e6)
        genes_keep <- rownames(counts_train)[rowSums(cpm_train >= cpm_min, na.rm = TRUE) >= cpm_min_samples]
        
        # Transform test data using training-only parameters
        # Apply same CPM filter and transformation
        counts_test_filtered <- counts_test[genes_keep, , drop = FALSE]
        lib_sizes_test <- colSums(counts_test_filtered, na.rm = TRUE)
        cpm_test <- t(t(counts_test_filtered) / lib_sizes_test * 1e6)
        mat_t_test <- log1p(cpm_test)
        mat_t_test <- t(mat_t_test)  # Transpose to samples x genes
        
        # Ensure test matrix has same genes as training (in same order)
        # mat_t_train columns are genes_keep (after filtering)
        # mat_t_test columns should match
        test_colnames <- colnames(mat_t_test)
        train_colnames <- colnames(mat_t_train)
        
        if (!setequal(test_colnames, train_colnames)) {
            # Align test to training gene order
            common_genes <- intersect(train_colnames, test_colnames)
            mat_t_test <- mat_t_test[, common_genes, drop = FALSE]
            mat_t_train <- mat_t_train[, common_genes, drop = FALSE]
        }
        
        # Store parameters for reproducibility
        params <- list(
            transform = transform,
            cpm_min = cpm_min,
            cpm_min_samples = cpm_min_samples,
            n_genes_kept = length(genes_keep),
            n_train_samples = ncol(counts_train),
            n_test_samples = ncol(counts_test)
        )
        
    } else {
        stop("Only 'log1p-cpm' transform is currently supported")
    }
    
    return(list(
        mat_t_train = mat_t_train,
        mat_t_test = mat_t_test,
        genes_keep = genes_keep,
        params = params
    ))
}

#' Select Top Genes by DE In-Fold (Anti-leakage)
#'
#' Selects top-K genes via differential expression analysis within CV fold,
#' using training data only to prevent data leakage.
#'
#' @param split An `rsample` split object (e.g., from `rsample::vfold_cv()`).
#' @param mat_t A samples x genes transformed matrix for full dataset.
#' @param pheno Data.frame with sample metadata; must include `sample_id` and group column.
#' @param group_col Character scalar; name of phenotype column with PS/PIS labels. Defaults to "group".
#' @param n Integer; number of genes to select. Defaults to 50.
#' @param method Character scalar; selection method: "de" (differential expression) or "variance". Defaults to "de".
#' @param seed Integer; random seed for reproducibility. Defaults to 123.
#'
#' @return A character vector of selected gene IDs (length `n` or fewer).
#'
#' @details
#' This function enforces anti-leakage by:
#' - Performing DE analysis on training data only
#' - Selecting genes based on training-only statistics
#' - Never using test/validation data for gene selection
#'
#' @examples
#' \dontrun{
#' library(rsample)
#' data(gse201926_trainmini)
#' mat_t <- esr_transform_log1p_cpm(gse201926_trainmini$counts)
#' set.seed(123)
#' splits <- vfold_cv(gse201926_trainmini$pheno, v = 3)
#' selected <- esr_selectDEInFold(splits$splits[[1]], 
#'                                mat_t,
#'                                gse201926_trainmini$pheno,
#'                                n = 50)
#' length(selected)
#' }
#' @export
esr_selectDEInFold <- function(split, mat_t, pheno, 
                               group_col = "group",
                               n = 50,
                               method = c("de", "variance"),
                               seed = 123) {
    if (!requireNamespace("rsample", quietly = TRUE)) {
        stop("rsample package is required for in-fold preprocessing")
    }
    
    method <- match.arg(method)
    
    # Extract training sample IDs from split
    train_data <- rsample::training(split)
    train_sample_ids <- train_data$sample_id
    
    # Match sample IDs to mat_t row names
    mat_t_rownames <- rownames(mat_t)
    if (is.null(mat_t_rownames)) {
        stop("mat_t matrix must have row names (sample IDs)")
    }
    
    # Diagnostic information for debugging
    train_indices <- match(train_sample_ids, mat_t_rownames)
    if (any(is.na(train_indices))) {
        missing_train <- train_sample_ids[is.na(train_indices)]
        stop(
            "Some training sample IDs not found in mat_t row names.\n",
            "Training sample IDs from split: ", paste(head(train_sample_ids, 5), collapse = ", "), "...\n",
            "mat_t row names: ", paste(head(mat_t_rownames, 5), collapse = ", "), "...\n",
            "Missing training IDs: ", paste(missing_train, collapse = ", "), "\n",
            "All training IDs: ", paste(train_sample_ids, collapse = ", "), "\n",
            "All mat_t rows: ", paste(mat_t_rownames, collapse = ", ")
        )
    }
    
    # Get training data only
    mat_t_train <- mat_t[train_indices, , drop = FALSE]
    pheno_train <- pheno[pheno$sample_id %in% train_sample_ids, , drop = FALSE]
    
    # Ensure ordering matches
    pheno_train <- pheno_train[match(rownames(mat_t_train), pheno_train$sample_id), , drop = FALSE]
    
    if (method == "de") {
        # Perform DE analysis on training data only
        de_table <- esr_analyzeDifferentialExpression(
            mat_t_train,
            pheno_train,
            group_col = group_col,
            seed = seed
        )
        
        # Select top-K genes from DE table
        selected_genes <- esr_selectTopGenes(
            de_table = de_table,
            n = n,
            by = "de"
        )
    } else if (method == "variance") {
        # Select top-K genes by variance on training data only
        selected_genes <- esr_selectTopGenes(
            mat_t = mat_t_train,
            n = n,
            by = "variance"
        )
    }
    
    return(selected_genes)
}

# [END]


