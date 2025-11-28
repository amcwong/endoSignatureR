#' Export Endometrial Signature
#'
#' Exports a trained signature to reproducible artifact files: CSV (signature panel and coefficients),
#' JSON (preprocessing recipe and training parameters), and Markdown (model card with provenance and performance).
#'
#' @param signature List containing signature object (from `result$signature` or standalone).
#'   Must include `panel`, `coefficients`, and `intercept` at minimum.
#' @param result Optional list from `esr_trainEndometrialSignature()`. If provided, extracts
#'   signature, metrics, calibration, stability, and other metadata for complete export.
#' @param dir Character; output directory. Defaults to "export".
#' @param formats Character vector; which formats to export: "csv", "json", "md", or "all".
#'   Defaults to c("csv", "json", "md").
#' @param include_intercept Logical; include intercept in CSV export (as separate row).
#'   Defaults to TRUE.
#' @param include_stability Logical; include bootstrap frequencies in CSV if available.
#'   Defaults to TRUE.
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly, a character vector of paths to exported files.
#'
#' @details
#' This function exports three artifact files:
#' - **endometrial_signature.csv**: Gene panel with coefficients, selection frequencies, and optional bootstrap frequencies
#' - **endometrial_recipe.json**: Preprocessing recipe, training parameters, signature metadata, and reproducibility information
#' - **endometrial_model_card.md**: Model documentation with provenance, performance metrics, limitations, and intended use
#'
#' If `result` is provided, the export includes full metadata (performance metrics, calibration info, stability info).
#' If only `signature` is provided, the export includes signature information only (model card will be minimal).
#'
#' @import readr
#' @import jsonlite
#' @examples
#' \dontrun{
#' library(rsample)
#' data(gse201926_trainmini)
#'
#' # Train signature
#' result <- esr_trainEndometrialSignature(
#'   X = gse201926_trainmini$counts,
#'   pheno = gse201926_trainmini$pheno,
#'   top_k = 100,
#'   seed = 123
#' )
#'
#' # Export all formats
#' paths <- esr_exportSignature(
#'   signature = result$signature,
#'   result = result,
#'   dir = tempdir()
#' )
#'
#' # Export only CSV
#' esr_exportSignature(
#'   signature = result$signature,
#'   result = result,
#'   formats = "csv",
#'   dir = tempdir()
#' )
#'
#' # Export standalone signature (no result)
#' esr_exportSignature(
#'   signature = result$signature,
#'   dir = tempdir()
#' )
#' }
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
#'
#' Ooms, J. (2014). The jsonlite package: A practical and consistent
#' mapping between JSON data and R objects. arXiv:1403.2805.
#' <https://arxiv.org/abs/1403.2805>.
#' @export
esr_exportSignature <- function(signature, result = NULL, dir = "export",
                                formats = c("csv", "json", "md"),
                                include_intercept = TRUE,
                                include_stability = TRUE,
                                ...) {
  # Validate inputs
  if (is.null(signature)) {
    stop("signature must be provided")
  }
  if (!is.list(signature)) {
    stop("signature must be a list")
  }
  if (!("panel" %in% names(signature))) {
    stop("signature must include 'panel' element")
  }
  if (!("coefficients" %in% names(signature))) {
    stop("signature must include 'coefficients' element")
  }
  if (!("intercept" %in% names(signature))) {
    stop("signature must include 'intercept' element")
  }

  # Normalize formats
  if ("all" %in% formats) {
    formats <- c("csv", "json", "md")
  }
  valid_formats <- c("csv", "json", "md")
  invalid_formats <- setdiff(formats, valid_formats)
  if (length(invalid_formats) > 0) {
    stop(
      "Invalid formats: ", paste(invalid_formats, collapse = ", "),
      ". Valid formats: ", paste(valid_formats, collapse = ", ")
    )
  }

  # Create directory if doesn't exist
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  # Extract data from result if provided
  if (!is.null(result)) {
    metrics <- result$metrics
    calibration <- result$calibration
    stability <- result$stability
    splits <- result$splits
    seeds <- result$seeds
    aggregation <- result$aggregation
  } else {
    metrics <- NULL
    calibration <- NULL
    stability <- NULL
    splits <- NULL
    seeds <- NULL
    aggregation <- NULL
  }

  # Export paths
  paths <- character(0)

  # Export CSV
  if ("csv" %in% formats) {
    csv_path <- export_signature_csv(
      signature = signature,
      stability = stability,
      dir = dir,
      include_intercept = include_intercept,
      include_stability = include_stability
    )
    paths <- c(paths, csv_path)
  }

  # Export JSON
  if ("json" %in% formats) {
    json_path <- export_signature_json(
      signature = signature,
      metrics = metrics,
      calibration = calibration,
      stability = stability,
      splits = splits,
      seeds = seeds,
      aggregation = aggregation,
      dir = dir
    )
    paths <- c(paths, json_path)
  }

  # Export Model Card
  if ("md" %in% formats) {
    md_path <- export_model_card(
      signature = signature,
      metrics = metrics,
      calibration = calibration,
      stability = stability,
      splits = splits,
      seeds = seeds,
      aggregation = aggregation,
      dir = dir
    )
    paths <- c(paths, md_path)
  }

  return(invisible(paths))
}

#' Export Signature to CSV
#'
#' Exports signature panel and coefficients to CSV format.
#'
#' @param signature Signature object.
#' @param stability Stability object (optional).
#' @param dir Output directory.
#' @param include_intercept Logical; include intercept row.
#' @param include_stability Logical; include bootstrap frequencies.
#' @return Path to exported CSV file.
#' @keywords internal
export_signature_csv <- function(signature, stability = NULL, dir = "export",
                                 include_intercept = TRUE,
                                 include_stability = TRUE) {
  # Extract gene data
  panel <- signature$panel
  coefficients <- signature$coefficients
  selection_frequency <- signature$selection_frequency

  # Match coefficients and frequencies to panel
  # Convert to numeric to remove names
  coef_vec <- as.numeric(coefficients[match(panel, names(coefficients))])
  # Convert selection_frequency to integer (it represents counts)
  freq_vec <- as.integer(selection_frequency[match(panel, names(selection_frequency))])

  # Extract bootstrap frequencies if available
  bootstrap_freq_vec <- NULL
  if (!is.null(stability) && include_stability && "bootstrap_frequency" %in% names(stability)) {
    bootstrap_freq <- stability$bootstrap_frequency
    bootstrap_freq_vec <- bootstrap_freq[match(panel, names(bootstrap_freq))]
    # Convert to numeric, handle missing values
    bootstrap_freq_vec <- as.numeric(bootstrap_freq_vec)
  }

  # Build data frame
  if (!is.null(bootstrap_freq_vec)) {
    signature_tbl <- data.frame(
      gene_id = panel,
      coefficient = coef_vec,
      selection_frequency = freq_vec,
      bootstrap_frequency = bootstrap_freq_vec,
      stringsAsFactors = FALSE
    )
  } else {
    signature_tbl <- data.frame(
      gene_id = panel,
      coefficient = coef_vec,
      selection_frequency = freq_vec,
      stringsAsFactors = FALSE
    )
  }

  # Sort by absolute coefficient magnitude (descending)
  signature_tbl <- signature_tbl[order(abs(signature_tbl$coefficient), decreasing = TRUE), , drop = FALSE]

  # Add intercept row if requested
  if (include_intercept && !is.null(signature$intercept)) {
    # Unwrap intercept if it's a named scalar
    intercept_value <- as.numeric(signature$intercept)
    intercept_row <- data.frame(
      gene_id = "(Intercept)",
      coefficient = intercept_value,
      selection_frequency = NA_integer_,
      stringsAsFactors = FALSE
    )
    if (!is.null(bootstrap_freq_vec)) {
      intercept_row$bootstrap_frequency <- NA_real_
    }
    # Intercept at the end
    signature_tbl <- rbind(signature_tbl, intercept_row)
  }

  # Write CSV
  csv_path <- file.path(dir, "endometrial_signature.csv")
  readr::write_csv(signature_tbl, csv_path)

  return(csv_path)
}

#' Export Recipe to JSON
#'
#' Exports preprocessing recipe and training parameters to JSON format.
#'
#' @param signature Signature object.
#' @param metrics Metrics object (optional).
#' @param calibration Calibration object (optional).
#' @param stability Stability object (optional).
#' @param splits Splits object (optional).
#' @param seeds Seeds object (optional).
#' @param aggregation Aggregation object (optional).
#' @param dir Output directory.
#' @return Path to exported JSON file.
#' @keywords internal
export_signature_json <- function(signature, metrics = NULL, calibration = NULL,
                                  stability = NULL, splits = NULL, seeds = NULL,
                                  aggregation = NULL, dir = "export") {
  # Extract preprocessing parameters
  recipe <- signature$recipe
  preprocessing <- list()
  if (!is.null(recipe)) {
    preprocessing$transform <- recipe$transform %||% "log1p-cpm"
    preprocessing$cpm_min <- recipe$cpm_min %||% 1.0
    preprocessing$cpm_min_samples <- recipe$cpm_min_samples %||% 4
    preprocessing$top_k <- recipe$top_k %||% 300
  } else {
    preprocessing$transform <- "log1p-cpm"
    preprocessing$cpm_min <- 1.0
    preprocessing$cpm_min_samples <- 4
    preprocessing$top_k <- 300
  }

  # Extract training parameters (from recipe or defaults)
  training <- list()
  training$outer <- if (!is.null(splits) && "outer_splits" %in% names(splits)) {
    if (is.list(splits$outer_splits) && length(splits$outer_splits) > 0) {
      # Try to infer from split structure
      if (length(splits$outer_splits) <= 3) "lpo" else "kfold"
    } else {
      "kfold"
    }
  } else {
    "kfold"
  }
  training$outer_folds <- if (!is.null(splits) && "n_outer_folds" %in% names(splits)) {
    splits$n_outer_folds
  } else {
    NULL
  }
  training$inner_folds <- recipe$inner_folds %||% 5
  training$inner_repeats <- recipe$inner_repeats %||% 10
  training$lambda_rule <- recipe$lambda_rule %||% "1se"
  training$min_folds <- if (!is.null(aggregation) && "min_folds" %in% names(aggregation)) {
    aggregation$min_folds
  } else {
    2
  }
  training$aggregation_method <- if (!is.null(aggregation) && "method" %in% names(aggregation)) {
    aggregation$method
  } else {
    "mean"
  }
  training$calibration_method <- if (!is.null(calibration) && "method" %in% names(calibration)) {
    calibration$method
  } else {
    "platt"
  }
  training$stability_selection <- if (!is.null(stability) && !is.null(stability$bootstrap_frequency)) {
    length(stability$bootstrap_frequency) > 0
  } else {
    FALSE
  }
  training$stability_resamples <- if (training$stability_selection && !is.null(recipe$stability_resamples)) {
    recipe$stability_resamples
  } else {
    100
  }

  # Extract signature metadata
  signature_meta <- list()
  signature_meta$n_genes <- length(signature$panel)
  # Unwrap intercept if it's a named scalar
  intercept_value <- signature$intercept %||% 0.0
  signature_meta$intercept <- as.numeric(intercept_value)
  signature_meta$gene_namespace <- recipe$gene_namespace %||% "ensembl_gene_id"

  # Extract reproducibility info
  reproducibility <- list()
  reproducibility$seed <- if (!is.null(seeds) && "main" %in% names(seeds)) {
    seeds$main
  } else if (!is.null(seeds) && length(seeds) > 0) {
    seeds[[1]]
  } else {
    NULL
  }
  reproducibility$package_version <- paste0("endoSignatureR_", utils::packageVersion("endoSignatureR"))
  reproducibility$r_version <- R.version.string

  # Build JSON structure
  json_data <- list(
    preprocessing = preprocessing,
    training = training,
    signature = signature_meta,
    reproducibility = reproducibility
  )

  # Write JSON
  json_path <- file.path(dir, "endometrial_recipe.json")
  jsonlite::write_json(
    json_data,
    json_path,
    pretty = TRUE,
    auto_unbox = TRUE
  )

  return(json_path)
}

#' Export Model Card to Markdown
#'
#' Exports model card in Markdown format documenting signature provenance, performance, limitations, and intended use.
#'
#' @param signature Signature object.
#' @param metrics Metrics object (optional).
#' @param calibration Calibration object (optional).
#' @param stability Stability object (optional).
#' @param splits Splits object (optional).
#' @param seeds Seeds object (optional).
#' @param aggregation Aggregation object (optional).
#' @param dir Output directory.
#' @return Path to exported Markdown file.
#' @keywords internal
export_model_card <- function(signature, metrics = NULL, calibration = NULL,
                              stability = NULL, splits = NULL, seeds = NULL,
                              aggregation = NULL, dir = "export") {
  # Get current date
  date_str <- format(Sys.Date(), "%Y-%m-%d")

  # Extract package version
  pkg_version <- paste0("endoSignatureR v", utils::packageVersion("endoSignatureR"))
  r_version <- R.version.string

  # Build model card
  lines <- c(
    "# Endometrial Signature Model Card",
    "",
    "## Model Details",
    paste0("- **Model Name**: Endometrial PS vs PIS Signature"),
    paste0("- **Version**: 1.0.0"),
    paste0("- **Date**: ", date_str),
    paste0("- **Package**: ", pkg_version),
    paste0("- **R Version**: ", r_version),
    ""
  )

  # Model Purpose
  lines <- c(
    lines,
    "## Model Purpose",
    "This signature classifies endometrial bulk RNA-seq samples as Proliferative Secretory (PS) or Proliferative Inflammatory Secretory (PIS) using LASSO logistic regression with nested cross-validation.",
    ""
  )

  # Training Data (minimal if no metrics)
  if (!is.null(metrics) && "predictions" %in% names(metrics)) {
    preds <- metrics$predictions
    n_samples <- if (is.data.frame(preds)) nrow(preds) else 0
    # Try to infer class balance from predictions
    if (is.data.frame(preds) && "label" %in% names(preds)) {
      n_ps <- sum(preds$label == 0, na.rm = TRUE)
      n_pis <- sum(preds$label == 1, na.rm = TRUE)
      class_balance <- paste0("Balanced (", n_ps, " PS, ", n_pis, " PIS)")
    } else {
      class_balance <- "Unknown"
    }
  } else {
    n_samples <- "Unknown"
    class_balance <- "Unknown"
  }

  lines <- c(
    lines,
    "## Training Data",
    paste0("- **Samples**: ", n_samples),
    paste0("- **Genes**: ", length(signature$panel), " (signature panel)"),
    paste0("- **Class Balance**: ", class_balance),
    ""
  )

  # Model Architecture
  recipe <- signature$recipe
  transform_method <- recipe$transform %||% "log1p-cpm"
  top_k <- recipe$top_k %||% 300
  cpm_min <- recipe$cpm_min %||% 1
  cpm_min_samples <- recipe$cpm_min_samples %||% 4
  lambda_rule <- recipe$lambda_rule %||% "1se"
  inner_folds <- recipe$inner_folds %||% 5
  inner_repeats <- recipe$inner_repeats %||% 10
  n_outer_folds <- if (!is.null(splits) && "n_outer_folds" %in% names(splits)) splits$n_outer_folds else "Unknown"
  aggregation_method <- if (!is.null(aggregation) && "method" %in% names(aggregation)) aggregation$method else "mean"
  min_folds <- if (!is.null(aggregation) && "min_folds" %in% names(aggregation)) aggregation$min_folds else 2

  lines <- c(
    lines,
    "## Model Architecture",
    paste0("- **Algorithm**: LASSO logistic regression (glmnet)"),
    paste0("- **Preprocessing**: ", transform_method, " transformation, CPM filtering (min=", cpm_min, ", min_samples=", cpm_min_samples, "), DE-based top-K selection (K=", top_k, ")"),
    paste0("- **CV Structure**: Nested CV (outer: ", n_outer_folds, "-fold, inner: ", inner_folds, "-fold with ", inner_repeats, " repeats)"),
    paste0("- **Lambda Selection**: ", lambda_rule, " rule (", if (lambda_rule == "1se") "1 standard error rule for sparsity" else "minimum CV error", ")"),
    paste0("- **Coefficient Aggregation**: ", aggregation_method, " aggregation across outer folds (Option 2)"),
    paste0("- **Consensus Genes**: Genes selected in >=", min_folds, " outer folds"),
    ""
  )

  # Performance Metrics
  if (!is.null(metrics)) {
    auc <- metrics$auc %||% NA_real_
    accuracy <- metrics$accuracy %||% NA_real_
    brier_score <- metrics$brier_score %||% NA_real_
    ece <- metrics$ece %||% NA_real_

    lines <- c(
      lines,
      "## Performance Metrics",
      paste0("- **AUC-ROC**: ", if (is.na(auc)) "N/A" else round(auc, 3)),
      paste0("- **Accuracy**: ", if (is.na(accuracy)) "N/A" else round(accuracy, 3)),
      paste0("- **Brier Score**: ", if (is.na(brier_score)) "N/A" else round(brier_score, 3)),
      paste0("- **Expected Calibration Error (ECE)**: ", if (is.na(ece)) "N/A" else round(ece, 3)),
      paste0("- **Signature Size**: ", length(signature$panel), " genes"),
      ""
    )
  } else {
    lines <- c(
      lines,
      "## Performance Metrics",
      paste0("- **Signature Size**: ", length(signature$panel), " genes"),
      "- Performance metrics not available (signature-only export)",
      ""
    )
  }

  # Calibration
  if (!is.null(calibration) && "method" %in% names(calibration)) {
    cal_method <- calibration$method
    lines <- c(
      lines,
      "## Calibration",
      paste0("- **Method**: ", cal_method, if (cal_method == "platt") " scaling" else if (cal_method == "isotonic") " regression" else ""),
      "- **Calibrated Probabilities**: Included in predictions",
      ""
    )
  } else {
    lines <- c(
      lines,
      "## Calibration",
      "- **Method**: Not available",
      ""
    )
  }

  # Stability Selection
  if (!is.null(stability) && !is.null(stability$bootstrap_frequency) && length(stability$bootstrap_frequency) > 0) {
    n_resamples <- recipe$stability_resamples %||% 100
    threshold <- recipe$stability_threshold %||% 0.7
    n_stable <- if (!is.null(stability$stable_genes)) length(stability$stable_genes) else 0

    lines <- c(
      lines,
      "## Stability Selection",
      paste0("- **Bootstrap Resamples**: ", n_resamples),
      paste0("- **Stability Threshold**: ", threshold),
      paste0("- **Stable Genes**: ", n_stable),
      ""
    )
  } else {
    lines <- c(
      lines,
      "## Stability Selection",
      "- **Status**: Not performed",
      ""
    )
  }

  # Limitations
  lines <- c(
    lines,
    "## Limitations",
    "- **Tissue Specificity**: Trained on endometrial tissue only; not portable across tissues",
    "- **Small n Uncertainty**: Limited sample size introduces uncertainty; avoid clinical claims",
    "- **Binary Only**: Supports only PS vs PIS classification; multiclass/continuous outcomes out of scope",
    "- **Batch Effects**: Users should provide batch information and use `~ batch + group` designs",
    ""
  )

  # Intended Use
  lines <- c(
    lines,
    "## Intended Use",
    "- **Use Case**: Classification of endometrial bulk RNA-seq samples into PS vs PIS states",
    "- **Target Audience**: Endometrial researchers, pathologists, clinicians",
    "- **Not Intended For**: Non-endometrial tissues, clinical decision-making without validation",
    ""
  )

  # Training Parameters
  main_seed <- if (!is.null(seeds) && "main" %in% names(seeds)) {
    seeds$main
  } else if (!is.null(seeds) && length(seeds) > 0) {
    seeds[[1]]
  } else {
    "Not recorded"
  }

  # Build training parameters lines
  training_param_lines <- c(
    "## Training Parameters",
    paste0("- **Outer CV**: ", n_outer_folds, "-fold (stratified)"),
    paste0("- **Inner CV**: ", inner_folds, "-fold with ", inner_repeats, " repeats"),
    paste0("- **Lambda Rule**: ", lambda_rule),
    paste0("- **Top-K Selection**: ", top_k, " genes"),
    paste0("- **Aggregation Method**: ", aggregation_method)
  )
  # Add calibration method if available
  if (!is.null(calibration) && "method" %in% names(calibration)) {
    training_param_lines <- c(training_param_lines, paste0("- **Calibration Method**: ", calibration$method))
  }
  training_param_lines <- c(training_param_lines, paste0("- **Seed**: ", main_seed), "")
  lines <- c(lines, training_param_lines)

  # Reproducibility
  outer_seed <- if (!is.null(seeds) && "outer" %in% names(seeds)) seeds$outer else NULL
  inner_seed <- if (!is.null(seeds) && "inner" %in% names(seeds)) seeds$inner else NULL

  lines <- c(
    lines,
    "## Reproducibility",
    paste0("- **Seeds**: ", if (!is.null(outer_seed) && !is.null(inner_seed)) {
      paste0("Outer CV seed=", outer_seed, ", Inner CV seed=", inner_seed)
    } else {
      paste0("Main seed=", main_seed)
    }),
    "- **CV Folds**: Pre-computed or generated deterministically",
    paste0("- **Package Versions**: ", pkg_version, ", glmnet, limma, ..."),
    ""
  )

  # References
  lines <- c(
    lines,
    "## References",
    "- LASSO: Tibshirani (1996), \"Regression shrinkage and selection via the lasso\"",
    "- Nested CV: Varma & Simon (2006), \"Bias in error estimation when using cross-validation for model selection\"",
    ""
  )

  # Write Markdown
  md_path <- file.path(dir, "endometrial_model_card.md")
  writeLines(lines, md_path)

  return(md_path)
}

#' Export Stability Frequencies to CSV
#'
#' Exports bootstrap stability frequencies to CSV format.
#'
#' @param stability Stability object containing bootstrap_frequency.
#' @param dir Output directory.
#' @return Path to exported CSV file.
#' @keywords internal
export_stability_csv <- function(stability, dir = "export") {
  # Validate inputs
  if (is.null(stability)) {
    stop("stability must be provided")
  }
  if (!is.list(stability)) {
    stop("stability must be a list")
  }
  if (!("bootstrap_frequency" %in% names(stability))) {
    stop("stability must include 'bootstrap_frequency' element")
  }

  # Extract bootstrap frequencies
  bootstrap_freq <- stability$bootstrap_frequency

  # Create data frame
  stability_df <- data.frame(
    gene_id = names(bootstrap_freq),
    bootstrap_frequency = as.numeric(bootstrap_freq),
    stringsAsFactors = FALSE
  )

  # Sort by frequency (descending)
  stability_df <- stability_df[order(stability_df$bootstrap_frequency, decreasing = TRUE), ]

  # Write CSV
  csv_path <- file.path(dir, "endometrial_stability.csv")
  readr::write_csv(stability_df, csv_path)

  return(csv_path)
}

# Internal helper for exporting predictions to CSV
export_predictions_csv <- function(predictions, dir = "export",
                                   filename = "endometrial_predictions.csv") {
  # Validate inputs
  if (is.null(predictions)) {
    stop("predictions must be provided")
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
  required_cols <- c("sample", "score", "probability", "prediction")
  missing_cols <- setdiff(required_cols, names(predictions_df))
  if (length(missing_cols) > 0) {
    stop(paste0("predictions must contain columns: ", paste(missing_cols, collapse = ", ")))
  }

  # Create directory if doesn't exist
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  # Prepare data for export
  # Select columns in order: sample, score, probability, prediction, confidence_lower, confidence_upper
  export_cols <- c("sample", "score", "probability", "prediction")
  if ("confidence_lower" %in% names(predictions_df)) {
    export_cols <- c(export_cols, "confidence_lower")
  }
  if ("confidence_upper" %in% names(predictions_df)) {
    export_cols <- c(export_cols, "confidence_upper")
  }

  # Extract columns
  export_df <- predictions_df[, export_cols, drop = FALSE]

  # Format predictions as character (PS/PIS)
  export_df$prediction <- as.character(export_df$prediction)

  # Round numeric columns to appropriate precision
  export_df$score <- round(export_df$score, 6)
  export_df$probability <- round(export_df$probability, 6)
  if ("confidence_lower" %in% names(export_df)) {
    export_df$confidence_lower <- round(export_df$confidence_lower, 6)
  }
  if ("confidence_upper" %in% names(export_df)) {
    export_df$confidence_upper <- round(export_df$confidence_upper, 6)
  }

  # Write CSV
  csv_path <- file.path(dir, filename)
  readr::write_csv(export_df, csv_path)

  return(csv_path)
}

#' Export Predictions to CSV
#'
#' Exports predictions from `esr_classifyEndometrial()` to CSV format.
#'
#' @param predictions A data.frame from `esr_classifyEndometrial()` or a list with `predictions` element.
#' @param dir Character; output directory. Defaults to "export".
#' @param formats Character vector; which formats to export (currently only "csv").
#'   Defaults to "csv".
#' @param filename Character; output filename. Defaults to "endometrial_predictions.csv".
#' @param ... Additional arguments (currently unused).
#'
#' @return Invisibly, a character vector of paths to exported files.
#'
#' @details
#' This function exports predictions to CSV format with columns:
#' - sample: Sample IDs
#' - score: Raw signature scores
#' - probability: Calibrated or sigmoid-transformed probabilities
#' - prediction: Binary predictions (PS/PIS)
#' - confidence_lower: Lower confidence bound (if available)
#' - confidence_upper: Upper confidence bound (if available)
#'
#' @examples
#' \dontrun{
#' # Load signature and classify samples
#' signature <- esr_loadPretrainedSignature()
#' data(gse201926_sample)
#' predictions <- esr_classifyEndometrial(
#'   X_new = gse201926_sample$counts,
#'   signature = signature
#' )
#'
#' # Export predictions
#' export_path <- esr_exportPredictions(
#'   predictions = predictions,
#'   dir = tempdir()
#' )
#'
#' # Read exported file
#' read.csv(export_path)
#' }
#'
#' @references
#' R Core Team (2025). R: A language and environment for statistical
#' computing. R Foundation for Statistical Computing, Vienna, Austria.
#' <https://www.R-project.org/>.
#' @export
esr_exportPredictions <- function(predictions, dir = "export",
                                  formats = "csv", filename = NULL, ...) {
  # Validate inputs
  if (is.null(predictions)) {
    stop("predictions must be provided")
  }

  # Normalize formats
  if ("all" %in% formats) {
    formats <- "csv"
  }
  valid_formats <- c("csv")
  invalid_formats <- setdiff(formats, valid_formats)
  if (length(invalid_formats) > 0) {
    stop(
      "Invalid formats: ", paste(invalid_formats, collapse = ", "),
      ". Valid formats: ", paste(valid_formats, collapse = ", ")
    )
  }

  # Create directory if doesn't exist
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  # Default filename
  if (is.null(filename)) {
    filename <- "endometrial_predictions.csv"
  }

  # Export paths
  paths <- character(0)

  # Export CSV
  if ("csv" %in% formats) {
    csv_path <- export_predictions_csv(
      predictions = predictions,
      dir = dir,
      filename = filename
    )
    paths <- c(paths, csv_path)
  }

  return(invisible(paths))
}


# Internal helper for exporting predictions to CSV
export_predictions_csv <- function(predictions, dir = "export",
                                   filename = "endometrial_predictions.csv") {
  # Validate inputs
  if (is.null(predictions)) {
    stop("predictions must be provided")
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
  required_cols <- c("sample", "score", "probability", "prediction")
  missing_cols <- setdiff(required_cols, names(predictions_df))
  if (length(missing_cols) > 0) {
    stop(paste0("predictions must contain columns: ", paste(missing_cols, collapse = ", ")))
  }

  # Create directory if doesn't exist
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  # Prepare data for export
  # Select columns in order: sample, score, probability, prediction, confidence_lower, confidence_upper
  export_cols <- c("sample", "score", "probability", "prediction")
  if ("confidence_lower" %in% names(predictions_df)) {
    export_cols <- c(export_cols, "confidence_lower")
  }
  if ("confidence_upper" %in% names(predictions_df)) {
    export_cols <- c(export_cols, "confidence_upper")
  }

  # Extract columns
  export_df <- predictions_df[, export_cols, drop = FALSE]

  # Format predictions as character (PS/PIS)
  export_df$prediction <- as.character(export_df$prediction)

  # Round numeric columns to appropriate precision
  export_df$score <- round(export_df$score, 6)
  export_df$probability <- round(export_df$probability, 6)
  if ("confidence_lower" %in% names(export_df)) {
    export_df$confidence_lower <- round(export_df$confidence_lower, 6)
  }
  if ("confidence_upper" %in% names(export_df)) {
    export_df$confidence_upper <- round(export_df$confidence_upper, 6)
  }

  # Write CSV
  csv_path <- file.path(dir, filename)
  readr::write_csv(export_df, csv_path)

  return(csv_path)
}

# Public export function for predictions
esr_exportPredictions <- function(predictions, dir = "export",
                                  formats = "csv", filename = NULL, ...) {
  # Validate inputs
  if (is.null(predictions)) {
    stop("predictions must be provided")
  }

  # Normalize formats
  if ("all" %in% formats) {
    formats <- "csv"
  }
  valid_formats <- c("csv")
  invalid_formats <- setdiff(formats, valid_formats)
  if (length(invalid_formats) > 0) {
    stop(
      "Invalid formats: ", paste(invalid_formats, collapse = ", "),
      ". Valid formats: ", paste(valid_formats, collapse = ", ")
    )
  }

  # Create directory if doesn't exist
  if (!dir.exists(dir)) {
    dir.create(dir, recursive = TRUE)
  }

  # Default filename
  if (is.null(filename)) {
    filename <- "endometrial_predictions.csv"
  }

  # Export paths
  paths <- character(0)

  # Export CSV
  if ("csv" %in% formats) {
    csv_path <- export_predictions_csv(
      predictions = predictions,
      dir = dir,
      filename = filename
    )
    paths <- c(paths, csv_path)
  }

  return(invisible(paths))
}


# Helper function for NULL coalescing
`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

# [END]
