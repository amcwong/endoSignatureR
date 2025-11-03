#' Train Endometrial Signature with Nested Cross-Validation
#'
#' Trains a PS vs PIS signature using nested cross-validation with glmnet (LASSO/Elastic Net).
#' Uses in-fold preprocessing and coefficient aggregation (Option 2) for stability with small n.
#'
#' @param X Matrix/data.frame of raw counts (genes x samples). Must have column names (sample IDs).
#' @param pheno Data.frame with sample metadata; must include `sample_id` and `group` columns.
#' @param transform Character scalar; transformation method. Defaults to "log1p-cpm".
#' @param cpm_min Numeric; minimum CPM threshold for gene filtering. Defaults to 1.
#' @param cpm_min_samples Integer; minimum number of samples that must meet CPM threshold. Defaults to 4.
#' @param top_k Integer; number of top genes to select via DE analysis. Defaults to 300.
#' @param outer Character scalar; outer CV type. "kfold" for K-fold CV or "lpo" for Leave-Pair-Out. Defaults to "kfold".
#' @param outer_folds Integer; number of outer CV folds (if outer = "kfold"). If NULL, uses folds_demo structure. Defaults to NULL.
#' @param inner_folds Integer; number of inner CV folds for λ tuning. Defaults to 5.
#' @param inner_repeats Integer; number of inner CV repeats. Defaults to 10.
#' @param lambda_rule Character scalar; λ selection rule. "1se" (1 standard error rule, sparser) or "min" (minimum CV error). Defaults to "1se".
#' @param min_folds Integer; minimum number of outer folds that must select a gene for consensus. Defaults to 2.
#' @param aggregation_method Character scalar; method for aggregating coefficients. "mean" or "median". Defaults to "mean".
#' @param seed Integer; random seed for reproducibility. Defaults to 123.
#' @param ... Additional arguments (currently unused; future: alpha for elastic net).
#'
#' @return A list with elements:
#' \describe{
#'   \item{signature}{List containing:
#'     \describe{
#'       \item{panel}{Character vector of consensus gene IDs with non-zero aggregated coefficients}
#'       \item{coefficients}{Named numeric vector of aggregated coefficients for panel genes}
#'       \item{intercept}{Numeric scalar; aggregated logistic regression intercept}
#'       \item{selection_frequency}{Named integer vector; number of folds that selected each gene}
#'       \item{fold_coefficients}{List of per-fold coefficients (one element per outer fold)}
#'       \item{recipe}{List of preprocessing parameters (aggregated across folds)}
#'     }
#'   }
#'   \item{metrics}{List containing:
#'     \describe{
#'       \item{auc}{Numeric scalar; AUC from outer CV predictions (aggregated)}
#'       \item{accuracy}{Numeric scalar; classification accuracy from outer CV (threshold 0.5)}
#'       \item{predictions}{Data.frame with columns: sample_id, fold, prob, label, pred}
#'     }
#'   }
#'   \item{splits}{List containing outer/inner splits used (or reference to folds_demo)}
#'   \item{seeds}{List of seeds used for reproducibility}
#'   \item{cv_results}{List of inner CV results (λ tuning per outer fold)}
#'   \item{aggregation}{List containing aggregation metadata (method, min_folds, consensus_genes)}
#' }
#'
#' @details
#' This function implements nested cross-validation with:
#' - **Outer CV**: Model evaluation on held-out test sets
#' - **Inner CV**: Hyperparameter tuning (λ selection) within each outer fold
#' - **In-fold preprocessing**: Transforms/filters and gene selection within folds only (anti-leakage)
#' - **Coefficient aggregation**: Option 2 - aggregates coefficients across outer folds for stability
#'
#' For small n (n=12) and weak signals, coefficient aggregation (Option 2) provides:
#' - Lower variance by averaging out fold-specific noise
#' - Higher stability by requiring consensus across folds
#' - Better generalization compared to single-fold signatures
#'
#' @examples
#' \dontrun{
#' library(rsample)
#' data(gse201926_trainmini)
#' data(folds_demo)
#'
#' # Train signature with default parameters
#' result <- esr_trainEndometrialSignature(
#'   X = gse201926_trainmini$counts,
#'   pheno = gse201926_trainmini$pheno,
#'   top_k = 100,
#'   outer_folds = NULL, # Use folds_demo
#'   seed = 123
#' )
#'
#' # Extract signature
#' signature <- result$signature
#' length(signature$panel)
#' head(signature$coefficients)
#'
#' # Check performance
#' result$metrics$auc
#' result$metrics$accuracy
#' }
#' @export
esr_trainEndometrialSignature <- function(X, pheno,
                                          transform = "log1p-cpm",
                                          cpm_min = 1,
                                          cpm_min_samples = 4,
                                          top_k = 300,
                                          outer = c("kfold", "lpo"),
                                          outer_folds = NULL,
                                          inner_folds = 5,
                                          inner_repeats = 10,
                                          lambda_rule = c("1se", "min"),
                                          min_folds = 2,
                                          aggregation_method = c("mean", "median"),
                                          seed = 123,
                                          ...) {
  # Check required packages
  if (!requireNamespace("glmnet", quietly = TRUE)) {
    stop("glmnet package is required for training")
  }
  if (!requireNamespace("rsample", quietly = TRUE)) {
    stop("rsample package is required for cross-validation")
  }

  # Match arguments
  transform <- match.arg(transform)
  outer <- match.arg(outer)
  lambda_rule <- match.arg(lambda_rule)
  aggregation_method <- match.arg(aggregation_method)

  # Set seed for reproducibility
  set.seed(seed)

  # Validate inputs
  if (!is.matrix(X) && !is.data.frame(X)) {
    stop("X must be a matrix or data.frame")
  }
  X <- as.matrix(X)

  if (is.null(colnames(X))) {
    stop("X must have column names (sample IDs)")
  }

  if (!is.data.frame(pheno)) {
    stop("pheno must be a data.frame")
  }

  if (!"sample_id" %in% names(pheno)) {
    stop("pheno must contain a 'sample_id' column")
  }

  if (!"group" %in% names(pheno)) {
    stop("pheno must contain a 'group' column")
  }

  # Validate group levels
  groups <- unique(pheno$group)
  if (length(groups) != 2) {
    stop("pheno$group must contain exactly 2 unique values (PS and PIS)")
  }

  # Match sample IDs
  sample_ids <- colnames(X)
  pheno_samples <- pheno$sample_id

  common_samples <- intersect(sample_ids, pheno_samples)
  if (length(common_samples) == 0) {
    stop("No matching sample IDs found between X column names and pheno$sample_id")
  }

  if (length(common_samples) < length(sample_ids)) {
    warning("Some samples in X are missing from pheno. Using only common samples.")
    X <- X[, common_samples, drop = FALSE]
  }

  # Subset and reorder pheno
  pheno_matched <- pheno[pheno$sample_id %in% common_samples, , drop = FALSE]
  pheno_matched <- pheno_matched[match(common_samples, pheno_matched$sample_id), , drop = FALSE]

  # Verify class balance
  group_counts <- table(pheno_matched$group)
  if (any(group_counts < 2)) {
    stop("Each group must have at least 2 samples for cross-validation")
  }

  # Create or use outer CV splits
  folds_data <- NULL
  if (is.null(outer_folds)) {
    # Try to use folds_demo if available
    tryCatch(
      {
        # Try to load folds_demo using data()
        data(folds_demo, package = "endoSignatureR", envir = environment())
        if (exists("folds_demo", envir = environment())) {
          folds_data <- get("folds_demo", envir = environment())
        }
      },
      error = function(e) {
        # If loading fails, folds_data remains NULL
      }
    )

    if (!is.null(folds_data) && "outer_splits" %in% names(folds_data)) {
      outer_splits <- folds_data$outer_splits
      n_outer_folds <- nrow(outer_splits)
      message("Using folds_demo for outer CV splits (", n_outer_folds, " folds)")
    } else {
      # Create outer splits
      set.seed(seed)
      outer_splits <- rsample::vfold_cv(
        data = pheno_matched,
        v = 3,
        strata = "group",
        breaks = 2
      )
      n_outer_folds <- nrow(outer_splits)
    }
  } else {
    # Use specified number of folds
    set.seed(seed)
    outer_splits <- rsample::vfold_cv(
      data = pheno_matched,
      v = outer_folds,
      strata = "group",
      breaks = 2
    )
    n_outer_folds <- outer_folds
  }

  # Store seeds
  seeds <- list(
    main = seed,
    outer = if (!is.null(folds_data) && "outer_seed" %in% names(folds_data)) folds_data$outer_seed else seed,
    inner = if (!is.null(folds_data) && "inner_seed" %in% names(folds_data)) folds_data$inner_seed else seed + 1
  )

  # Initialize storage for outer loop
  outer_predictions <- list()
  outer_signatures <- list()
  inner_cv_results <- list()
  fold_recipes <- list()

  # Outer CV loop
  for (fold_idx in seq_len(n_outer_folds)) {
    tryCatch(
      {
        outer_split <- outer_splits$splits[[fold_idx]]
        train_pheno <- rsample::training(outer_split)
        test_pheno <- rsample::testing(outer_split)

        test_sample_ids <- test_pheno$sample_id

        # Apply in-fold preprocessing (transform + filter)
        transform_result <- esr_transformInFold(
          split = outer_split,
          counts = X,
          pheno = pheno_matched,
          transform = transform,
          cpm_min = cpm_min,
          cpm_min_samples = cpm_min_samples
        )

        mat_t_train <- transform_result$mat_t_train
        mat_t_test <- transform_result$mat_t_test
        genes_keep <- transform_result$genes_keep

        # Select top-K genes via DE (in-fold, anti-leakage)
        selected_genes <- esr_selectDEInFold(
          split = outer_split,
          mat_t = rbind(mat_t_train, mat_t_test), # Full transformed matrix
          pheno = pheno_matched,
          group_col = "group",
          n = top_k,
          method = "de",
          seed = seed + fold_idx
        )

        # Subset to selected genes
        mat_t_train_selected <- mat_t_train[, selected_genes, drop = FALSE]
        mat_t_test_selected <- mat_t_test[, selected_genes, drop = FALSE]

        # Store recipe for this fold
        fold_recipes[[fold_idx]] <- list(
          transform = transform,
          cpm_min = cpm_min,
          cpm_min_samples = cpm_min_samples,
          top_k = top_k,
          n_genes_kept = length(genes_keep),
          n_genes_selected = length(selected_genes),
          genes_keep = genes_keep,
          selected_genes = selected_genes
        )

        # Prepare labels for glmnet (convert PS/PIS to 0/1)
        train_labels <- as.integer(train_pheno$group == groups[2]) # PIS = 1, PS = 0
        test_labels <- as.integer(test_pheno$group == groups[2])

        # Inner CV for λ tuning
        set.seed(seeds$inner + fold_idx)

        # Create inner splits on outer training data
        inner_train_pheno <- train_pheno
        if (nrow(inner_train_pheno) >= 4) {
          inner_splits <- rsample::vfold_cv(
            data = inner_train_pheno,
            v = inner_folds,
            strata = "group",
            breaks = 2
          )
        } else {
          # Too small for inner CV, skip inner loop and use cv.glmnet directly
          inner_splits <- NULL
        }

        # Use cv.glmnet for λ tuning (handles inner CV internally)
        if (!is.null(inner_splits)) {
          # Manual inner CV if needed, but cv.glmnet is simpler
          cv_fit <- glmnet::cv.glmnet(
            x = mat_t_train_selected,
            y = train_labels,
            family = "binomial",
            nfolds = inner_folds,
            alpha = 1.0, # LASSO
            standardize = TRUE,
            type.measure = "auc"
          )
        } else {
          # Use cv.glmnet directly (no manual inner CV)
          cv_fit <- glmnet::cv.glmnet(
            x = mat_t_train_selected,
            y = train_labels,
            family = "binomial",
            nfolds = min(inner_folds, nrow(mat_t_train_selected)),
            alpha = 1.0,
            standardize = TRUE,
            type.measure = "auc"
          )
        }

        # Select λ based on rule
        if (lambda_rule == "1se") {
          lambda_selected <- cv_fit$lambda.1se
        } else {
          lambda_selected <- cv_fit$lambda.min
        }

        # Store inner CV results
        inner_cv_results[[fold_idx]] <- list(
          lambda_min = cv_fit$lambda.min,
          lambda_1se = cv_fit$lambda.1se,
          lambda_selected = lambda_selected,
          lambda_rule = lambda_rule,
          cvm = cv_fit$cvm,
          cvsd = cv_fit$cvsd,
          nzero = cv_fit$nzero
        )

        # Train final model on full outer training data with selected λ
        final_fit <- glmnet::glmnet(
          x = mat_t_train_selected,
          y = train_labels,
          family = "binomial",
          alpha = 1.0,
          standardize = TRUE,
          lambda = lambda_selected
        )

        # Extract signature from this fold
        # Use stats::coef() generic (glmnet provides a method)
        coef_fit <- coef(final_fit, s = lambda_selected)
        coef_vec <- as.numeric(coef_fit)
        names(coef_vec) <- rownames(coef_fit)

        # Remove intercept (stored separately)
        intercept_fold <- coef_vec[1]
        coef_fold <- coef_vec[-1]

        # Get non-zero coefficients (signature genes)
        non_zero_genes <- names(coef_fold)[coef_fold != 0]
        coef_fold <- coef_fold[non_zero_genes]

        # Store signature for this fold
        outer_signatures[[fold_idx]] <- list(
          panel = non_zero_genes,
          coefficients = coef_fold,
          intercept = intercept_fold,
          n_genes = length(non_zero_genes)
        )

        # Predict on outer test data
        # Use stats::predict() generic (glmnet provides a method)
        predictions_prob <- predict(final_fit,
          newx = mat_t_test_selected,
          s = lambda_selected,
          type = "response"
        )
        predictions_prob <- as.numeric(predictions_prob)
        predictions_binary <- as.integer(predictions_prob >= 0.5)

        # Store predictions
        outer_predictions[[fold_idx]] <- data.frame(
          sample_id = test_sample_ids,
          fold = fold_idx,
          prob = predictions_prob,
          label = test_labels,
          pred = predictions_binary,
          stringsAsFactors = FALSE
        )
      },
      error = function(e) {
        warning("Error in outer CV fold ", fold_idx, ": ", conditionMessage(e))
        # Store empty predictions for this fold to maintain structure
        outer_predictions[[fold_idx]] <- data.frame(
          sample_id = character(0),
          fold = integer(0),
          prob = numeric(0),
          label = integer(0),
          pred = integer(0),
          stringsAsFactors = FALSE
        )
        # Store empty signature for this fold to maintain structure
        outer_signatures[[fold_idx]] <- list(
          panel = character(0),
          coefficients = numeric(0),
          intercept = 0.0,
          n_genes = 0L
        )
        # Store empty inner CV results
        inner_cv_results[[fold_idx]] <- list(
          lambda_min = NA_real_,
          lambda_1se = NA_real_,
          lambda_selected = NA_real_,
          lambda_rule = lambda_rule,
          error = conditionMessage(e)
        )
      }
    )
  }

  # Aggregate predictions across all outer folds
  # Handle empty list case to avoid NULL return
  if (length(outer_predictions) == 0) {
    all_predictions <- data.frame(
      sample_id = character(0),
      fold = integer(0),
      prob = numeric(0),
      label = integer(0),
      pred = integer(0),
      stringsAsFactors = FALSE
    )
  } else {
    all_predictions <- do.call(rbind, outer_predictions)
  }

  # Compute performance metrics
  # Simple AUC computation (trapezoidal rule)
  compute_auc <- function(labels, probs) {
    if (length(labels) == 0 || length(probs) == 0) {
      return(NA_real_)
    }
    if (length(unique(labels)) != 2) {
      return(NA_real_)
    }

    # Sort by probabilities (descending)
    ord <- order(probs, decreasing = TRUE)
    labels_sorted <- labels[ord]

    # Count true positives and false positives at each threshold
    n_pos <- sum(labels_sorted == 1)
    n_neg <- sum(labels_sorted == 0)

    if (n_pos == 0 || n_neg == 0) {
      return(NA_real_)
    }

    # Compute TPR and FPR
    cum_tp <- cumsum(labels_sorted == 1)
    cum_fp <- cumsum(labels_sorted == 0)

    tpr <- cum_tp / n_pos
    fpr <- cum_fp / n_neg

    # AUC using trapezoidal rule
    auc <- sum(diff(c(0, fpr, 1)) * (c(0, tpr) + c(tpr, 1)) / 2)

    return(auc)
  }

  # Compute metrics safely, handling empty predictions
  if (nrow(all_predictions) == 0) {
    auc_value <- NA_real_
    accuracy_value <- NA_real_
  } else {
    auc_value <- compute_auc(all_predictions$label, all_predictions$prob)
    accuracy_value <- mean(all_predictions$pred == all_predictions$label)
  }

  # Aggregate signatures across outer folds (Option 2)
  # Handle case where all folds failed (outer_signatures is empty)
  if (length(outer_signatures) == 0) {
    # All folds failed - create empty signature structure
    all_genes <- character(0)
    selection_freq <- integer(0)
    consensus_genes <- character(0)
    aggregated_coef <- numeric(0)
    aggregated_intercept <- 0.0
    selection_freq_subset <- integer(0)
    warning("All outer CV folds failed - returning empty signature. Check for errors in preprocessing or model training.")
  } else {
    all_genes <- unique(unlist(lapply(outer_signatures, function(sig) sig$panel)))

    # Count selection frequency (as integer)
    if (length(all_genes) > 0) {
      selection_freq <- rep(0L, length(all_genes))
      names(selection_freq) <- all_genes
      storage.mode(selection_freq) <- "integer"

      for (sig in outer_signatures) {
        selected_in_fold <- sig$panel
        if (length(selected_in_fold) > 0) {
          selection_freq[selected_in_fold] <- selection_freq[selected_in_fold] + 1
        }
      }

      # Filter to consensus genes (selected in ≥min_folds folds)
      consensus_genes <- names(selection_freq)[selection_freq >= min_folds]
    } else {
      selection_freq <- integer(0)
      consensus_genes <- character(0)
    }

    # Aggregate coefficients for consensus genes
    if (length(consensus_genes) > 0) {
      aggregated_coef <- numeric(length(consensus_genes))
      names(aggregated_coef) <- consensus_genes

      for (gene in consensus_genes) {
        coef_values <- numeric()
        for (sig in outer_signatures) {
          if (gene %in% sig$panel && length(sig$coefficients) > 0) {
            coef_values <- c(coef_values, sig$coefficients[gene])
          }
        }

        if (aggregation_method == "mean") {
          aggregated_coef[gene] <- mean(coef_values, na.rm = TRUE)
        } else {
          aggregated_coef[gene] <- median(coef_values, na.rm = TRUE)
        }
      }

      # Aggregate intercept
      intercepts <- sapply(outer_signatures, function(sig) if (length(sig$intercept) > 0) sig$intercept else 0.0)
      if (aggregation_method == "mean") {
        aggregated_intercept <- mean(intercepts, na.rm = TRUE)
      } else {
        aggregated_intercept <- median(intercepts, na.rm = TRUE)
      }

      # Ensure selection_frequency is integer type
      selection_freq_subset <- as.integer(selection_freq[consensus_genes])
      names(selection_freq_subset) <- consensus_genes
    } else {
      # No consensus genes - warn and use genes from fold with most genes (fallback)
      fold_sizes <- sapply(outer_signatures, function(sig) if (length(sig$panel) > 0) length(sig$panel) else 0)
      if (length(fold_sizes) > 0 && max(fold_sizes) > 0) {
        best_fold <- which.max(fold_sizes)
        if (best_fold > 0 && best_fold <= length(outer_signatures)) {
          consensus_genes <- outer_signatures[[best_fold]]$panel
          aggregated_coef <- outer_signatures[[best_fold]]$coefficients
          aggregated_intercept <- if (length(outer_signatures[[best_fold]]$intercept) > 0) outer_signatures[[best_fold]]$intercept else 0.0
          selection_freq_subset <- rep(1L, length(consensus_genes))
          names(selection_freq_subset) <- consensus_genes
          if (length(consensus_genes) == 0) {
            selection_freq_subset <- integer(0)
          }
          warning("No consensus genes found (selected in ≥", min_folds, " folds). Using genes from fold ", best_fold, " (", length(consensus_genes), " genes).")
        } else {
          consensus_genes <- character(0)
          aggregated_coef <- numeric(0)
          aggregated_intercept <- 0.0
          selection_freq_subset <- integer(0)
        }
      } else {
        # All folds have empty signatures
        consensus_genes <- character(0)
        aggregated_coef <- numeric(0)
        aggregated_intercept <- 0.0
        selection_freq_subset <- integer(0)
        warning("No consensus genes found and all folds have empty signatures. Returning empty signature.")
      }
    }
  }

  # Store per-fold coefficients for debugging
  fold_coefficients <- lapply(outer_signatures, function(sig) sig$coefficients)
  names(fold_coefficients) <- paste0("fold_", seq_along(fold_coefficients))

  # Aggregate recipes (use most common settings)
  aggregated_recipe <- list(
    transform = transform,
    cpm_min = cpm_min,
    cpm_min_samples = cpm_min_samples,
    top_k = top_k,
    n_folds = n_outer_folds,
    recipes = fold_recipes
  )

  # Build return structure
  result <- list(
    signature = list(
      panel = consensus_genes,
      coefficients = aggregated_coef,
      intercept = aggregated_intercept,
      selection_frequency = selection_freq_subset,
      fold_coefficients = fold_coefficients,
      recipe = aggregated_recipe
    ),
    metrics = list(
      auc = auc_value,
      accuracy = accuracy_value,
      predictions = all_predictions
    ),
    splits = list(
      outer_splits = outer_splits,
      n_outer_folds = n_outer_folds
    ),
    seeds = seeds,
    cv_results = inner_cv_results,
    aggregation = list(
      method = aggregation_method,
      min_folds = min_folds,
      consensus_genes = consensus_genes,
      n_consensus_genes = length(consensus_genes),
      selection_frequencies = selection_freq
    )
  )

  return(result)
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
#' result <- esr_transformInFold(
#'   splits$splits[[1]],
#'   gse201926_trainmini$counts,
#'   gse201926_trainmini$pheno
#' )
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

  # Clean sample IDs (remove quotes and whitespace)
  train_sample_ids <- trimws(gsub('^"|"$', "", as.character(train_sample_ids)))
  test_sample_ids <- trimws(gsub('^"|"$', "", as.character(test_sample_ids)))

  # Match sample IDs to counts matrix column names
  counts_colnames <- colnames(counts)
  if (is.null(counts_colnames)) {
    stop("counts matrix must have column names (sample IDs)")
  }

  # Clean counts column names (remove quotes and whitespace)
  counts_colnames_clean <- trimws(gsub('^"|"$', "", as.character(counts_colnames)))

  # Diagnostic information for debugging
  if (any(is.na(match(train_sample_ids, counts_colnames_clean)))) {
    missing_train <- train_sample_ids[is.na(match(train_sample_ids, counts_colnames_clean))]
    stop(
      "Some training sample IDs not found in counts column names.\n",
      "Training sample IDs from split: ", paste(head(train_sample_ids, 5), collapse = ", "), "...\n",
      "Counts column names: ", paste(head(counts_colnames_clean, 5), collapse = ", "), "...\n",
      "Missing training IDs: ", paste(missing_train, collapse = ", "), "\n",
      "All training IDs: ", paste(train_sample_ids, collapse = ", "), "\n",
      "All counts columns: ", paste(counts_colnames_clean, collapse = ", ")
    )
  }

  train_indices <- match(train_sample_ids, counts_colnames_clean)
  test_indices <- match(test_sample_ids, counts_colnames_clean)

  if (any(is.na(train_indices))) {
    missing_train <- train_sample_ids[is.na(train_indices)]
    stop(
      "Some training sample IDs not found in counts column names.\n",
      "Missing: ", paste(missing_train, collapse = ", "), "\n",
      "Available: ", paste(head(counts_colnames_clean, 10), collapse = ", ")
    )
  }
  if (any(is.na(test_indices))) {
    missing_test <- test_sample_ids[is.na(test_indices)]
    stop(
      "Some test sample IDs not found in counts column names.\n",
      "Missing: ", paste(missing_test, collapse = ", "), "\n",
      "Available: ", paste(head(counts_colnames_clean, 10), collapse = ", ")
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
      cpm_min_samples = cpm_min_samples
    )

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
    mat_t_test <- t(mat_t_test) # Transpose to samples x genes

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
#'   mat_t,
#'   gse201926_trainmini$pheno,
#'   n = 50
#' )
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

  # Clean sample IDs (remove quotes and whitespace)
  train_sample_ids <- trimws(gsub('^"|"$', "", as.character(train_sample_ids)))

  # Match sample IDs to mat_t row names
  mat_t_rownames <- rownames(mat_t)
  if (is.null(mat_t_rownames)) {
    stop("mat_t matrix must have row names (sample IDs)")
  }

  # Clean mat_t rownames (remove quotes and whitespace)
  mat_t_rownames_clean <- trimws(gsub('^"|"$', "", as.character(mat_t_rownames)))

  # Diagnostic information for debugging
  train_indices <- match(train_sample_ids, mat_t_rownames_clean)
  if (any(is.na(train_indices))) {
    missing_train <- train_sample_ids[is.na(train_indices)]
    stop(
      "Some training sample IDs not found in mat_t row names.\n",
      "Training sample IDs from split: ", paste(head(train_sample_ids, 5), collapse = ", "), "...\n",
      "mat_t row names: ", paste(head(mat_t_rownames_clean, 5), collapse = ", "), "...\n",
      "Missing training IDs: ", paste(missing_train, collapse = ", "), "\n",
      "All training IDs: ", paste(train_sample_ids, collapse = ", "), "\n",
      "All mat_t rows: ", paste(mat_t_rownames_clean, collapse = ", ")
    )
  }

  # Get training data only (use original rownames for subsetting)
  # Match cleaned IDs back to original rownames
  train_indices_original <- which(mat_t_rownames_clean %in% train_sample_ids)
  mat_t_train <- mat_t[train_indices_original, , drop = FALSE]

  # Clean pheno sample_id for matching
  pheno_sample_id_clean <- trimws(gsub('^"|"$', "", as.character(pheno$sample_id)))
  pheno_train <- pheno[pheno_sample_id_clean %in% train_sample_ids, , drop = FALSE]

  # Ensure ordering matches
  mat_t_train_rownames_clean <- trimws(gsub('^"|"$', "", as.character(rownames(mat_t_train))))
  pheno_train_sample_id_clean <- trimws(gsub('^"|"$', "", as.character(pheno_train$sample_id)))
  pheno_train <- pheno_train[match(mat_t_train_rownames_clean, pheno_train_sample_id_clean), , drop = FALSE]

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
