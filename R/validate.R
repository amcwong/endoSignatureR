#' Validate Endometrial Data Structures
#'
#' Performs basic schema checks for endometrial bulk RNA-seq inputs including
#' ID matching, class imbalance detection, and annotation validation.
#'
#' @param X A matrix/data.frame of gene expression (genes x samples).
#' @param pheno A data.frame with sample metadata; must include `sample_id` and `group` columns by default.
#' @param annot Optional data.frame with gene annotations.
#' @param label_col Character scalar; name of the phenotype column containing PS/PIS labels. Defaults to "group".
#'
#' @return A list with elements:
#' \describe{
#'   \item{X}{The input expression matrix (cleaned/converted to matrix)}
#'   \item{pheno}{The input phenotype data.frame}
#'   \item{annot}{The input annotation data.frame (or NULL)}
#'   \item{issues}{A data.frame with columns `type` and `message` containing detected validation issues}
#' }
#'
#' @import limma
#' @examples
#' data(gse201926_sample)
#' result <- esr_validateEndometrial(
#'     gse201926_sample$counts,
#'     gse201926_sample$pheno,
#'     annot = gse201926_sample$annot
#' )
#' result$issues
#' @export
esr_validateEndometrial <- function(X, pheno, annot = NULL, label_col = "group") {
    issues <- data.frame(type = character(), message = character(), stringsAsFactors = FALSE)

    # Convert X to matrix if needed
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }

    # Validate X is numeric
    if (!is.numeric(X)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "Expression matrix X must contain numeric values",
            stringsAsFactors = FALSE
        ))
        return(list(X = X, pheno = pheno, annot = annot, issues = issues))
    }

    # Validate X dimensions
    if (nrow(X) == 0 || ncol(X) == 0) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "Expression matrix X must have at least one gene and one sample",
            stringsAsFactors = FALSE
        ))
        return(list(X = X, pheno = pheno, annot = annot, issues = issues))
    }

    # Validate pheno structure
    if (!is.data.frame(pheno)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "pheno must be a data.frame",
            stringsAsFactors = FALSE
        ))
        return(list(X = X, pheno = pheno, annot = annot, issues = issues))
    }

    # Check required columns in pheno
    if (!"sample_id" %in% names(pheno)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "pheno must contain a 'sample_id' column",
            stringsAsFactors = FALSE
        ))
    }

    if (!label_col %in% names(pheno)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = paste0("pheno must contain a '", label_col, "' column"),
            stringsAsFactors = FALSE
        ))
    }

    # Check ID matching between X and pheno
    if ("sample_id" %in% names(pheno) && !is.null(colnames(X))) {
        X_samples <- colnames(X)
        pheno_samples <- pheno$sample_id
        missing_in_pheno <- setdiff(X_samples, pheno_samples)
        missing_in_X <- setdiff(pheno_samples, X_samples)

        if (length(missing_in_pheno) > 0) {
            issues <- rbind(issues, data.frame(
                type = "warning",
                message = paste0("Samples in X not found in pheno: ", paste(missing_in_pheno, collapse = ", ")),
                stringsAsFactors = FALSE
            ))
        }

        if (length(missing_in_X) > 0) {
            issues <- rbind(issues, data.frame(
                type = "warning",
                message = paste0("Samples in pheno not found in X: ", paste(missing_in_X, collapse = ", ")),
                stringsAsFactors = FALSE
            ))
        }

        # Calculate mapping rate
        if (length(X_samples) > 0 && length(pheno_samples) > 0) {
            mapped <- length(intersect(X_samples, pheno_samples))
            mapping_rate <- mapped / max(length(X_samples), length(pheno_samples))
            if (mapping_rate < 0.8) {
                issues <- rbind(issues, data.frame(
                    type = "warning",
                    message = paste0("Low ID mapping rate: ", round(mapping_rate * 100, 1), "% of samples matched"),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    # Check class balance
    if (label_col %in% names(pheno)) {
        labels <- pheno[[label_col]]
        label_table <- table(labels)

        if (length(label_table) >= 2) {
            proportions <- prop.table(label_table)
            min_prop <- min(proportions)
            max_prop <- max(proportions)

            if (min_prop < 0.2) {
                issues <- rbind(issues, data.frame(
                    type = "warning",
                    message = paste0(
                        "Class imbalance detected: smallest class is ",
                        round(min_prop * 100, 1), "% of samples"
                    ),
                    stringsAsFactors = FALSE
                ))
            }

            if (max_prop > 0.8) {
                issues <- rbind(issues, data.frame(
                    type = "warning",
                    message = paste0(
                        "Class imbalance detected: largest class is ",
                        round(max_prop * 100, 1), "% of samples"
                    ),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }

    # Validate annot structure if provided
    if (!is.null(annot)) {
        if (!is.data.frame(annot)) {
            issues <- rbind(issues, data.frame(
                type = "warning",
                message = "annot should be a data.frame; ignoring",
                stringsAsFactors = FALSE
            ))
        } else if (!is.null(rownames(X))) {
            # Check if gene IDs match
            annot_gene_ids <- annot$GeneID
            if (!is.null(annot_gene_ids)) {
                X_genes <- rownames(X)
                matched_genes <- length(intersect(X_genes, annot_gene_ids))
                if (matched_genes == 0) {
                    issues <- rbind(issues, data.frame(
                        type = "warning",
                        message = "No matching gene IDs found between X rownames and annot$GeneID",
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }

    result <- list(X = X, pheno = pheno, annot = annot, issues = issues)
    return(result)
}

#' Transform Counts to log1p-CPM
#'
#' Converts raw counts to counts per million (CPM) with log1p transformation,
#' optionally filtering low-expressed genes.
#'
#' @param X A matrix/data.frame of raw counts (genes x samples).
#' @param cpm_min Minimum CPM threshold for gene filtering. Defaults to 1.
#' @param cpm_min_samples Minimum number of samples that must meet CPM threshold. Defaults to 4.
#'
#' @return A transformed matrix (samples x genes) ready for PCA/visualization.
#'   Rows are samples, columns are genes. Genes failing CPM filter are removed.
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' dim(mat_t) # samples x genes
#' @export
esr_transform_log1p_cpm <- function(X, cpm_min = 1, cpm_min_samples = 4) {
    # Convert to matrix if needed
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }

    # Calculate library sizes (column sums)
    lib_sizes <- colSums(X, na.rm = TRUE)

    # Calculate CPM (counts per million)
    cpm <- t(t(X) / lib_sizes * 1e6)

    # Filter genes with CPM < cpm_min in fewer than cpm_min_samples
    genes_keep <- rowSums(cpm >= cpm_min, na.rm = TRUE) >= cpm_min_samples

    # Subset and apply log1p transformation
    cpm_filtered <- cpm[genes_keep, , drop = FALSE]
    mat_t <- log1p(cpm_filtered)

    # Transpose to samples x genes format
    mat_t <- t(mat_t)

    return(mat_t)
}

#' Analyze Differential Expression
#'
#' Performs differential expression analysis using limma on transformed data
#' with a design matrix for PS vs PIS comparison.
#'
#' @param mat_t A samples x genes transformed matrix (e.g., from `esr_transform_log1p_cpm()`).
#' @param pheno A data.frame with sample metadata; must include `sample_id` and group column.
#' @param group_col Character scalar; name of the phenotype column containing PS/PIS labels. Defaults to "group".
#' @param transform Character scalar; type of transform used (informational only, not used in analysis). Defaults to "log1p-cpm".
#' @param contrast Character vector or NULL; contrasts to extract. If NULL, defaults to comparing last level vs first level of group. Defaults to NULL.
#' @param fdr_method Character scalar; method for FDR adjustment. Defaults to "BH" (Benjamini-Hochberg).
#' @param seed Integer; random seed for reproducibility. Defaults to 123.
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{gene_id}{Character; gene identifiers (rownames from mat_t)}
#'   \item{log2FC}{Numeric; log2 fold change}
#'   \item{pvalue}{Numeric; unadjusted p-value}
#'   \item{FDR}{Numeric; FDR-adjusted p-value}
#'   \item{AveExpr}{Numeric; mean expression across samples}
#'   \item{t}{Numeric; moderated t-statistic}
#' }
#' Results are sorted by FDR (ascending), then by absolute log2FC (descending).
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
#' head(de_table)
#' @importFrom utils head
#' @export
esr_analyzeDifferentialExpression <- function(mat_t, pheno, group_col = "group",
                                              transform = "log1p-cpm",
                                              contrast = NULL,
                                              fdr_method = "BH",
                                              seed = 123) {
    if (!requireNamespace("limma", quietly = TRUE)) {
        stop("limma package is required for differential expression analysis")
    }

    # Set seed for reproducibility
    set.seed(seed)

    # Validate mat_t is a matrix
    if (!is.matrix(mat_t)) {
        mat_t <- as.matrix(mat_t)
    }

    # Validate pheno structure
    if (!is.data.frame(pheno)) {
        stop("pheno must be a data.frame")
    }

    # Check required columns
    if (!"sample_id" %in% names(pheno)) {
        stop("pheno must contain a 'sample_id' column")
    }

    if (!group_col %in% names(pheno)) {
        stop(paste0("pheno must contain a '", group_col, "' column"))
    }

    # Check ID matching
    mat_samples <- rownames(mat_t)
    pheno_samples <- pheno$sample_id

    if (is.null(mat_samples)) {
        stop("mat_t must have rownames (sample IDs)")
    }

    # Match samples
    common_samples <- intersect(mat_samples, pheno_samples)

    if (length(common_samples) == 0) {
        stop("No matching sample IDs found between mat_t rownames and pheno$sample_id")
    }

    # Subset to common samples
    mat_t_matched <- mat_t[common_samples, , drop = FALSE]
    pheno_matched <- pheno[pheno$sample_id %in% common_samples, , drop = FALSE]

    # Ensure pheno is ordered to match mat_t
    pheno_matched <- pheno_matched[match(common_samples, pheno_matched$sample_id), , drop = FALSE]

    # Check sample size per group
    group_factor <- factor(pheno_matched[[group_col]])
    group_counts <- table(group_factor)

    if (length(group_counts) < 2) {
        stop("At least 2 groups are required for differential expression analysis")
    }

    if (any(group_counts < 2)) {
        warning("Some groups have fewer than 2 samples; DE analysis may be unreliable")
    }

    if (any(group_counts < 4)) {
        warning("Some groups have fewer than 4 samples; DE analysis has limited power")
    }

    # Create design matrix
    design <- stats::model.matrix(~group_factor)
    colnames(design)[1] <- "(Intercept)"

    # Fit linear model
    fit <- limma::lmFit(t(mat_t_matched), design) # limma expects genes x samples
    fit <- limma::eBayes(fit)

    # Extract contrasts if specified, otherwise use default (last vs first level)
    if (!is.null(contrast)) {
        # User-specified contrasts
        cont_matrix <- limma::makeContrasts(contrasts = contrast, levels = design)
        fit2 <- limma::contrasts.fit(fit, cont_matrix)
        fit2 <- limma::eBayes(fit2)

        # Extract top table
        top_table <- limma::topTable(fit2, number = Inf, adjust.method = fdr_method, sort.by = "none")
    } else {
        # Default: compare last level vs first level
        # Extract top table for second coefficient (group comparison)
        top_table <- limma::topTable(fit, coef = 2, number = Inf, adjust.method = fdr_method, sort.by = "none")
    }

    # Format output
    de_table <- data.frame(
        gene_id = rownames(top_table),
        log2FC = top_table$logFC,
        pvalue = top_table$P.Value,
        FDR = top_table$adj.P.Val,
        AveExpr = top_table$AveExpr,
        t = top_table$t,
        stringsAsFactors = FALSE
    )

    # Sort by FDR (ascending), then by absolute log2FC (descending)
    de_table <- de_table[order(de_table$FDR, -abs(de_table$log2FC)), ]
    rownames(de_table) <- NULL

    return(de_table)
}

#' Compute Confusion Matrix and Performance Metrics
#'
#' Computes confusion matrix and performance metrics for binary classification predictions.
#'
#' @param predictions A data.frame from `esr_classifyEndometrial()` with `probability` column,
#'   or a list with `predictions` element.
#' @param labels A factor or character vector of true labels (PS/PIS or 0/1).
#' @param threshold Decision threshold for positive class. Defaults to 0.5.
#' @return A list containing:
#' \describe{
#'   \item{confusion_matrix}{2x2 matrix with predicted (rows) and actual (columns) labels}
#'   \item{metrics}{Data.frame with performance metrics (accuracy, sensitivity, specificity, PPV, NPV, F1)}
#'   \item{threshold}{Numeric; threshold used}
#' }
#'
#' @examples
#' \dontrun{
#' data(gse201926_sample)
#' signature <- esr_loadPretrainedSignature()
#' predictions <- esr_classifyEndometrial(
#'     X_new = gse201926_sample$counts,
#'     signature = signature,
#'     threshold = 0.5
#' )
#' labels <- gse201926_sample$pheno$group
#' cm <- esr_computeConfusionMatrix(predictions, labels, threshold = 0.5)
#' cm$confusion_matrix
#' cm$metrics
#' }
#'
#' @export
esr_computeConfusionMatrix <- function(predictions, labels, threshold = 0.5) {
    # Validate inputs
    if (is.null(predictions)) {
        stop("predictions must be provided")
    }

    if (is.null(labels)) {
        stop("labels must be provided")
    }

    # Handle list result (if alerts included)
    if (is.list(predictions) && "predictions" %in% names(predictions)) {
        predictions_df <- predictions$predictions
    } else if (is.data.frame(predictions)) {
        predictions_df <- predictions
    } else {
        stop("predictions must be a data.frame or list with predictions element")
    }

    # Validate required columns
    if (!"probability" %in% names(predictions_df)) {
        stop("predictions must contain 'probability' column")
    }

    # Check lengths match
    if (nrow(predictions_df) != length(labels)) {
        stop("predictions and labels must have same length")
    }

    # Validate threshold
    if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
        stop("threshold must be numeric between 0 and 1")
    }

    # Extract probabilities
    probs <- predictions_df$probability

    # Convert labels to factor (PS/PIS)
    if (is.factor(labels)) {
        labels_char <- as.character(labels)
    } else {
        labels_char <- as.character(labels)
    }

    # Handle binary labels (0/1)
    labels_binary <- ifelse(labels_char == "PS" | labels_char == "0" | labels_char == 0, 0, 1)
    labels_factor <- factor(ifelse(labels_binary == 0, "PS", "PIS"), levels = c("PS", "PIS"))

    # Apply threshold to get binary predictions
    preds_binary <- ifelse(probs >= threshold, 1, 0)
    preds_factor <- factor(ifelse(preds_binary == 0, "PS", "PIS"), levels = c("PS", "PIS"))

    # Compute confusion matrix
    cm <- table(Predicted = preds_factor, Actual = labels_factor)

    # Extract values
    tn <- cm["PS", "PS"]
    fp <- cm["PIS", "PS"]
    fn <- cm["PS", "PIS"]
    tp <- cm["PIS", "PIS"]

    # Handle NA values (if no samples in a category)
    tn <- if (is.na(tn)) 0 else tn
    fp <- if (is.na(fp)) 0 else fp
    fn <- if (is.na(fn)) 0 else fn
    tp <- if (is.na(tp)) 0 else tp

    # Compute performance metrics
    total <- tp + tn + fp + fn

    accuracy <- if (total > 0) (tp + tn) / total else NA_real_
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    f1 <- if ((ppv + sensitivity) > 0) 2 * (ppv * sensitivity) / (ppv + sensitivity) else NA_real_

    # Create metrics data.frame
    metrics_df <- data.frame(
        metric = c(
            "Accuracy", "Sensitivity (Recall)", "Specificity",
            "PPV (Precision)", "NPV", "F1 Score"
        ),
        value = c(accuracy, sensitivity, specificity, ppv, npv, f1),
        stringsAsFactors = FALSE
    )

    # Return results
    return(list(
        confusion_matrix = cm,
        metrics = metrics_df,
        threshold = threshold
    ))
}

#' Compare Performance Across Multiple Thresholds
#'
#' Compares performance metrics across multiple thresholds for threshold selection.
#'
#' @param predictions A data.frame from `esr_classifyEndometrial()` with `probability` column,
#'   or a list with `predictions` element.
#' @param labels A factor or character vector of true labels (PS/PIS or 0/1).
#' @param thresholds Numeric vector of thresholds to compare. Defaults to seq(0.1, 0.9, 0.1).
#' @return A data.frame with columns: threshold, accuracy, sensitivity, specificity, ppv, npv, f1
#'
#' @examples
#' \dontrun{
#' data(gse201926_sample)
#' signature <- esr_loadPretrainedSignature()
#' predictions <- esr_classifyEndometrial(
#'     X_new = gse201926_sample$counts,
#'     signature = signature,
#'     threshold = 0.5
#' )
#' labels <- gse201926_sample$pheno$group
#' comparison <- esr_compareThresholds(predictions, labels, thresholds = seq(0.3, 0.7, 0.1))
#' print(comparison)
#' }
#'
#' @export
esr_compareThresholds <- function(predictions, labels, thresholds = seq(0.1, 0.9, 0.1)) {
    # Validate inputs
    if (is.null(predictions)) {
        stop("predictions must be provided")
    }

    if (is.null(labels)) {
        stop("labels must be provided")
    }

    if (!is.numeric(thresholds) || any(thresholds < 0 | thresholds > 1)) {
        stop("thresholds must be numeric values between 0 and 1")
    }

    # Compute metrics for each threshold
    results_list <- lapply(thresholds, function(thresh) {
        cm_result <- esr_computeConfusionMatrix(predictions, labels, threshold = thresh)
        metrics <- cm_result$metrics

        # Extract values
        accuracy <- metrics$value[metrics$metric == "Accuracy"]
        sensitivity <- metrics$value[metrics$metric == "Sensitivity (Recall)"]
        specificity <- metrics$value[metrics$metric == "Specificity"]
        ppv <- metrics$value[metrics$metric == "PPV (Precision)"]
        npv <- metrics$value[metrics$metric == "NPV"]
        f1 <- metrics$value[metrics$metric == "F1 Score"]

        return(data.frame(
            threshold = thresh,
            accuracy = accuracy,
            sensitivity = sensitivity,
            specificity = specificity,
            ppv = ppv,
            npv = npv,
            f1 = f1,
            stringsAsFactors = FALSE
        ))
    })

    # Combine results
    comparison_df <- do.call(rbind, results_list)

    return(comparison_df)
}



# Compute Confusion Matrix and Performance Metrics
#
# Computes confusion matrix and performance metrics for binary classification predictions.
#
# @param predictions A data.frame from `esr_classifyEndometrial()` with `probability` column,
#   or a list with `predictions` element.
# @param labels A factor or character vector of true labels (PS/PIS or 0/1).
# @param threshold Decision threshold for positive class. Defaults to 0.5.
# @return A list containing:
#   \describe{
#     \item{confusion_matrix}{2x2 matrix with predicted (rows) and actual (columns) labels}
#     \item{metrics}{Data.frame with performance metrics (accuracy, sensitivity, specificity, PPV, NPV, F1)}
#     \item{threshold}{Numeric; threshold used}
#   }
#
# @examples
# \dontrun{
# data(gse201926_sample)
# signature <- esr_loadPretrainedSignature()
# predictions <- esr_classifyEndometrial(
#   X_new = gse201926_sample$counts,
#   signature = signature,
#   threshold = 0.5
# )
# labels <- gse201926_sample$pheno$group
# cm <- esr_computeConfusionMatrix(predictions, labels, threshold = 0.5)
# cm$confusion_matrix
# cm$metrics
# }
#
# @export
esr_computeConfusionMatrix <- function(predictions, labels, threshold = 0.5) {
    # Validate inputs
    if (is.null(predictions)) {
        stop("predictions must be provided")
    }

    if (is.null(labels)) {
        stop("labels must be provided")
    }

    # Handle list result (if alerts included)
    if (is.list(predictions) && "predictions" %in% names(predictions)) {
        predictions_df <- predictions$predictions
    } else if (is.data.frame(predictions)) {
        predictions_df <- predictions
    } else {
        stop("predictions must be a data.frame or list with predictions element")
    }

    # Validate required columns
    if (!"probability" %in% names(predictions_df)) {
        stop("predictions must contain 'probability' column")
    }

    # Check lengths match
    if (nrow(predictions_df) != length(labels)) {
        stop("predictions and labels must have same length")
    }

    # Validate threshold
    if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
        stop("threshold must be numeric between 0 and 1")
    }

    # Extract probabilities
    probs <- predictions_df$probability

    # Convert labels to factor (PS/PIS)
    if (is.factor(labels)) {
        labels_char <- as.character(labels)
    } else {
        labels_char <- as.character(labels)
    }

    # Handle binary labels (0/1)
    labels_binary <- ifelse(labels_char == "PS" | labels_char == "0" | labels_char == 0, 0, 1)
    labels_factor <- factor(ifelse(labels_binary == 0, "PS", "PIS"), levels = c("PS", "PIS"))

    # Apply threshold to get binary predictions
    preds_binary <- ifelse(probs >= threshold, 1, 0)
    preds_factor <- factor(ifelse(preds_binary == 0, "PS", "PIS"), levels = c("PS", "PIS"))

    # Compute confusion matrix
    cm <- table(Predicted = preds_factor, Actual = labels_factor)

    # Extract values
    tn <- cm["PS", "PS"]
    fp <- cm["PIS", "PS"]
    fn <- cm["PS", "PIS"]
    tp <- cm["PIS", "PIS"]

    # Handle NA values (if no samples in a category)
    tn <- if (is.na(tn)) 0 else tn
    fp <- if (is.na(fp)) 0 else fp
    fn <- if (is.na(fn)) 0 else fn
    tp <- if (is.na(tp)) 0 else tp

    # Compute performance metrics
    total <- tp + tn + fp + fn

    accuracy <- if (total > 0) (tp + tn) / total else NA_real_
    sensitivity <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
    specificity <- if ((tn + fp) > 0) tn / (tn + fp) else NA_real_
    ppv <- if ((tp + fp) > 0) tp / (tp + fp) else NA_real_
    npv <- if ((tn + fn) > 0) tn / (tn + fn) else NA_real_
    f1 <- if ((ppv + sensitivity) > 0) 2 * (ppv * sensitivity) / (ppv + sensitivity) else NA_real_

    # Create metrics data.frame
    metrics_df <- data.frame(
        metric = c(
            "Accuracy", "Sensitivity (Recall)", "Specificity",
            "PPV (Precision)", "NPV", "F1 Score"
        ),
        value = c(accuracy, sensitivity, specificity, ppv, npv, f1),
        stringsAsFactors = FALSE
    )

    # Return results
    return(list(
        confusion_matrix = cm,
        metrics = metrics_df,
        threshold = threshold
    ))
}

# Compare Performance Across Multiple Thresholds
#
# Compares performance metrics across multiple thresholds for threshold selection.
#
# @param predictions A data.frame from `esr_classifyEndometrial()` with `probability` column,
#   or a list with `predictions` element.
# @param labels A factor or character vector of true labels (PS/PIS or 0/1).
# @param thresholds Numeric vector of thresholds to compare. Defaults to seq(0.1, 0.9, 0.1).
# @return A data.frame with columns: threshold, accuracy, sensitivity, specificity, ppv, npv, f1
#
# @examples
# \dontrun{
# data(gse201926_sample)
# signature <- esr_loadPretrainedSignature()
# predictions <- esr_classifyEndometrial(
#   X_new = gse201926_sample$counts,
#   signature = signature,
#   threshold = 0.5
# )
# labels <- gse201926_sample$pheno$group
# comparison <- esr_compareThresholds(predictions, labels, thresholds = seq(0.3, 0.7, 0.1))
# print(comparison)
# }
#
# @export
esr_compareThresholds <- function(predictions, labels, thresholds = seq(0.1, 0.9, 0.1)) {
    # Validate inputs
    if (is.null(predictions)) {
        stop("predictions must be provided")
    }

    if (is.null(labels)) {
        stop("labels must be provided")
    }

    if (!is.numeric(thresholds) || any(thresholds < 0 | thresholds > 1)) {
        stop("thresholds must be numeric values between 0 and 1")
    }

    # Compute metrics for each threshold
    results_list <- lapply(thresholds, function(thresh) {
        cm_result <- esr_computeConfusionMatrix(predictions, labels, threshold = thresh)
        metrics <- cm_result$metrics

        # Extract values
        accuracy <- metrics$value[metrics$metric == "Accuracy"]
        sensitivity <- metrics$value[metrics$metric == "Sensitivity (Recall)"]
        specificity <- metrics$value[metrics$metric == "Specificity"]
        ppv <- metrics$value[metrics$metric == "PPV (Precision)"]
        npv <- metrics$value[metrics$metric == "NPV"]
        f1 <- metrics$value[metrics$metric == "F1 Score"]

        return(data.frame(
            threshold = thresh,
            accuracy = accuracy,
            sensitivity = sensitivity,
            specificity = specificity,
            ppv = ppv,
            npv = npv,
            f1 = f1,
            stringsAsFactors = FALSE
        ))
    })

    # Combine results
    comparison_df <- do.call(rbind, results_list)

    return(comparison_df)
}

# [END]
