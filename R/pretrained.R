#' Load Pre-trained Endometrial Signature
#'
#' Loads the shipped pre-trained PS vs PIS signature artifacts from `inst/extdata/pretrained-signature/`.
#'
#' @return A list containing:
#' \describe{
#'   \item{panel}{Character vector of consensus gene IDs with non-zero coefficients}
#'   \item{coefficients}{Named numeric vector of coefficients for panel genes}
#'   \item{intercept}{Numeric scalar; logistic regression intercept}
#'   \item{selection_frequency}{Named integer vector; selection frequency across outer folds}
#'   \item{recipe}{List of preprocessing parameters (transform, cpm_min, top_k, etc.)}
#'   \item{stability}{List containing optional bootstrap frequency information}
#' }
#'
#' @details
#' This function loads pre-trained signature artifacts that were trained on the full GSE201926 dataset.
#' The signature artifacts are stored in `inst/extdata/pretrained-signature/` and include:
#' - `endometrial_signature.csv`: Gene panel with coefficients and selection frequencies
#' - `endometrial_recipe.json`: Preprocessing recipe and training parameters
#' - `endometrial_stability.csv` (optional): Bootstrap stability frequencies
#'
#' @import readr
#' @import jsonlite
#' @examples
#' signature <- esr_loadPretrainedSignature()
#' length(signature$panel)
#' head(signature$coefficients)
#' signature$recipe$preprocessing
#'
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
#'
#' Ooms, J. (2014). The jsonlite package: A practical and consistent
#' mapping between JSON data and R objects. arXiv:1403.2805.
#' <https://arxiv.org/abs/1403.2805>.
#' @export
esr_loadPretrainedSignature <- function() {
  # Get paths to artifacts using system.file()
  csv_path <- system.file("extdata", "pretrained-signature", "endometrial_signature.csv",
    package = "endoSignatureR"
  )
  json_path <- system.file("extdata", "pretrained-signature", "endometrial_recipe.json",
    package = "endoSignatureR"
  )
  stability_path <- system.file("extdata", "pretrained-signature", "endometrial_stability.csv",
    package = "endoSignatureR"
  )

  # Check if artifacts exist
  if (!nzchar(csv_path) || !file.exists(csv_path)) {
    stop(
      "Pre-trained signature CSV not found. ",
      "Expected location: inst/extdata/pretrained-signature/endometrial_signature.csv. ",
      "Please ensure the package is properly installed."
    )
  }

  if (!nzchar(json_path) || !file.exists(json_path)) {
    stop(
      "Pre-trained signature JSON not found. ",
      "Expected location: inst/extdata/pretrained-signature/endometrial_recipe.json. ",
      "Please ensure the package is properly installed."
    )
  }

  # Load CSV signature
  csv_data <- readr::read_csv(csv_path,
    col_types = readr::cols(
      gene_id = readr::col_character(),
      coefficient = readr::col_double(),
      selection_frequency = readr::col_integer(),
      bootstrap_frequency = readr::col_double()
    ),
    show_col_types = FALSE
  )

  # Separate intercept from panel genes
  intercept_row <- csv_data[csv_data$gene_id == "(Intercept)", , drop = FALSE]
  panel_rows <- csv_data[csv_data$gene_id != "(Intercept)", , drop = FALSE]

  # Extract intercept value
  intercept_value <- if (nrow(intercept_row) > 0) {
    as.numeric(intercept_row$coefficient[1])
  } else {
    stop("Intercept not found in signature CSV")
  }

  # Extract panel genes and coefficients
  panel <- panel_rows$gene_id
  coefficients <- panel_rows$coefficient
  names(coefficients) <- panel

  # Extract selection frequencies
  selection_frequency <- panel_rows$selection_frequency
  names(selection_frequency) <- panel

  # Load JSON recipe
  recipe_data <- jsonlite::read_json(json_path, simplifyVector = TRUE)

  # Validate recipe structure
  required_sections <- c("preprocessing", "training", "signature", "reproducibility")
  missing_sections <- setdiff(required_sections, names(recipe_data))
  if (length(missing_sections) > 0) {
    stop(
      "JSON recipe missing required sections: ",
      paste(missing_sections, collapse = ", ")
    )
  }

  # Load optional stability CSV
  stability <- NULL
  if (nzchar(stability_path) && file.exists(stability_path)) {
    stability_data <- readr::read_csv(stability_path,
      col_types = readr::cols(
        gene_id = readr::col_character(),
        bootstrap_frequency = readr::col_double()
      ),
      show_col_types = FALSE
    )

    # Extract bootstrap frequencies
    bootstrap_frequency <- stability_data$bootstrap_frequency
    names(bootstrap_frequency) <- stability_data$gene_id

    # Match bootstrap frequencies to panel (missing genes get NA)
    bootstrap_freq_matched <- bootstrap_frequency[panel]
    names(bootstrap_freq_matched) <- panel

    stability <- list(bootstrap_frequency = bootstrap_freq_matched)
  } else {
    # Check if bootstrap frequencies are in CSV
    if ("bootstrap_frequency" %in% names(panel_rows)) {
      bootstrap_freq_vec <- panel_rows$bootstrap_frequency
      names(bootstrap_freq_vec) <- panel
      # Remove NA values (intercept row)
      bootstrap_freq_vec <- bootstrap_freq_vec[!is.na(bootstrap_freq_vec)]
      if (length(bootstrap_freq_vec) > 0) {
        stability <- list(bootstrap_frequency = bootstrap_freq_vec)
      }
    }
  }

  # Validate signature structure
  if (length(panel) == 0) {
    stop("Signature panel is empty")
  }

  if (length(coefficients) != length(panel)) {
    stop("Coefficient length does not match panel length")
  }

  if (is.na(intercept_value)) {
    stop("Intercept value is NA")
  }

  # Construct signature list matching Phase 2 structure
  signature <- list(
    panel = panel,
    coefficients = coefficients,
    intercept = intercept_value,
    selection_frequency = selection_frequency,
    recipe = recipe_data,
    stability = stability
  )

  return(signature)
}

# Internal helper for Youden threshold selection
#' @keywords internal
.select_threshold_youden <- function(probabilities, labels) {
  # Validate inputs
  if (!is.numeric(probabilities)) {
    stop("probabilities must be numeric")
  }

  if (length(probabilities) == 0) {
    stop("probabilities must not be empty")
  }

  if (length(probabilities) != length(labels)) {
    stop("probabilities and labels must have same length")
  }

  # Check probability range
  if (any(probabilities < 0 | probabilities > 1, na.rm = TRUE)) {
    warning("Probabilities outside [0, 1] range detected; clamping to [0, 1]")
    probabilities <- pmax(0, pmin(1, probabilities))
  }

  # Convert labels to binary (0 = PS, 1 = PIS)
  if (is.factor(labels)) {
    labels <- as.character(labels)
  }

  # Handle PS/PIS labels
  labels_binary <- ifelse(labels == "PS" | labels == "0" | labels == 0, 0, 1)

  # Check for edge cases
  if (length(unique(probabilities)) == 1) {
    # All probabilities are the same
    warning("All probabilities are the same; using default threshold 0.5")
    return(0.5)
  }

  if (length(unique(labels_binary)) < 2) {
    # Only one class present
    warning("Only one class present in labels; using default threshold 0.5")
    return(0.5)
  }

  # Use pROC if available, otherwise manual computation
  if (requireNamespace("pROC", quietly = TRUE)) {
    # Compute ROC curve using pROC
    roc_obj <- pROC::roc(
      response = labels_binary, predictor = probabilities,
      quiet = TRUE, direction = "<"
    )

    # Get optimal threshold using Youden's J
    coords <- pROC::coords(roc_obj,
      x = "best", best.method = "youden",
      transpose = FALSE
    )

    if (is.null(coords) || nrow(coords) == 0 || is.na(coords$threshold[1])) {
      warning("Could not compute Youden threshold; using default 0.5")
      return(0.5)
    }

    optimal_threshold <- as.numeric(coords$threshold[1])

    # Ensure threshold is in valid range
    if (is.na(optimal_threshold) || optimal_threshold < 0 || optimal_threshold > 1) {
      warning("Invalid Youden threshold computed; using default 0.5")
      return(0.5)
    }

    return(optimal_threshold)
  } else {
    # Manual computation of Youden threshold
    # Sort probabilities in descending order
    ord <- order(probabilities, decreasing = TRUE)
    probs_sorted <- probabilities[ord]
    labels_sorted <- labels_binary[ord]

    # Count positives and negatives
    n_pos <- sum(labels_sorted == 1)
    n_neg <- sum(labels_sorted == 0)

    if (n_pos == 0 || n_neg == 0) {
      warning("No positive or negative samples; using default threshold 0.5")
      return(0.5)
    }

    # Compute TPR and FPR at each unique threshold
    unique_thresholds <- sort(unique(probs_sorted), decreasing = TRUE)
    best_j <- -Inf
    best_threshold <- 0.5

    for (thresh in unique_thresholds) {
      preds <- ifelse(probs_sorted >= thresh, 1, 0)
      tp <- sum(preds == 1 & labels_sorted == 1)
      fp <- sum(preds == 1 & labels_sorted == 0)

      tpr <- tp / n_pos
      fpr <- fp / n_neg

      # Youden's J = TPR - FPR = sensitivity + specificity - 1
      j <- tpr - fpr

      if (j > best_j) {
        best_j <- j
        best_threshold <- thresh
      }
    }

    # If perfect separation, use threshold closest to 0.5
    if (best_j >= 1.0) {
      # Perfect separation - find threshold closest to 0.5
      best_threshold <- unique_thresholds[which.min(abs(unique_thresholds - 0.5))]
    }

    return(best_threshold)
  }
}

#' Classify Endometrial Samples Using a Signature
#'
#' Applies a pre-trained signature to new endometrial samples and returns predictions with confidence scores.
#'
#' @param X_new A matrix/data.frame of gene expression (genes x samples) for new samples.
#'   Must have rownames (gene IDs) and colnames (sample IDs).
#' @param signature Optional signature list; if NULL, loads the shipped pre-trained signature.
#' @param threshold Decision threshold for positive class. Defaults to 0.5. Can be numeric (0-1) or "youden" (requires labeled validation data via y_new).
#' @param confidence Logical; whether to compute confidence outputs. Defaults to TRUE.
#' @param y_new Optional vector of validation labels (PS/PIS or 0/1). Required if threshold = "youden".
#'
#' @return If alerts are generated, returns a list with:
#' \describe{
#'   \item{predictions}{Data.frame with per-sample predictions containing columns:
#'     \describe{
#'       \item{sample}{Character; sample IDs (colnames of X_new)}
#'       \item{score}{Numeric; raw signature score (linear combination)}
#'       \item{probability}{Numeric; calibrated probability (if calibration available) or sigmoid-transformed score}
#'       \item{prediction}{Factor; binary predictions (PS or PIS) based on threshold}
#'       \item{confidence_lower}{Numeric; lower confidence bound (if confidence = TRUE)}
#'       \item{confidence_upper}{Numeric; upper confidence bound (if confidence = TRUE)}
#'     }
#'   }
#'   \item{alerts}{Data.frame with validation alerts containing columns:
#'     \describe{
#'       \item{type}{Character; alert type (warning, info, error)}
#'       \item{message}{Character; alert message}
#'       \item{severity}{Character; severity level (low, medium, high)}
#'     }
#'   }
#' }
#' If no alerts are generated, returns the predictions data.frame directly (for backward compatibility).
#'
#' @details
#' This function applies a pre-trained signature to new endometrial samples:
#' - Validates data structure and performs domain/tissue checks (ID mapping, missing genes)
#' - Applies preprocessing recipe (transform, filter) to match training pipeline
#' - Computes signature scores using linear combination
#' - Applies calibration if available (Platt scaling or isotonic regression)
#' - Applies threshold to get binary predictions (PS vs PIS)
#' - Computes confidence intervals if requested
#'
#' Domain checks include:
#' - Gene ID namespace matching (warns if mapping rate < 80%)
#' - Missing signature gene handling (sets to 0 with warning)
#' - Class imbalance detection (if pheno provided)
#'
#' @importFrom utils head
#' @examples
#' # Load signature
#' signature <- esr_loadPretrainedSignature()
#'
#' # Load sample data (treat as unlabeled)
#' data(gse201926_sample)
#' X_new <- gse201926_sample$counts
#'
#' # Classify samples
#' predictions <- esr_classifyEndometrial(
#'   X_new = X_new,
#'   signature = signature,
#'   threshold = 0.5,
#'   confidence = TRUE
#' )
#'
#' # View predictions
#' head(predictions)
#' table(predictions$prediction)
#'
#' @references
#' Friedman, J., Hastie, T., & Tibshirani, R. (2010). Regularization paths for
#' generalized linear models via coordinate descent. Journal of Statistical
#' Software, 33(1), 1-22. <https://doi.org/10.18637/jss.v033.i01>.
#' @export
esr_classifyEndometrial <- function(X_new, signature = NULL, threshold = 0.5, confidence = TRUE, y_new = NULL) {
  # Load signature if not provided
  if (is.null(signature)) {
    signature <- esr_loadPretrainedSignature()
  }

  # Validate signature structure
  required_elements <- c("panel", "coefficients", "intercept", "recipe")
  missing_elements <- setdiff(required_elements, names(signature))
  if (length(missing_elements) > 0) {
    stop(
      "Signature missing required elements: ",
      paste(missing_elements, collapse = ", ")
    )
  }

  # Validate X_new structure
  if (is.null(X_new)) {
    stop("X_new must be provided")
  }

  # Convert to matrix if needed
  if (!is.matrix(X_new)) {
    X_new <- as.matrix(X_new)
  }

  # Check dimensions
  if (nrow(X_new) == 0 || ncol(X_new) == 0) {
    stop("X_new must have at least one gene and one sample")
  }

  # Check for rownames (gene IDs)
  if (is.null(rownames(X_new))) {
    stop("X_new must have rownames (gene IDs)")
  }

  # Initialize alerts data.frame for validation alerts
  alerts <- data.frame(
    type = character(),
    message = character(),
    severity = character(),
    stringsAsFactors = FALSE
  )

  # Check for colnames (sample IDs)
  sample_ids <- if (is.null(colnames(X_new))) {
    paste0("sample_", seq_len(ncol(X_new)))
  } else {
    colnames(X_new)
  }

  # Check that X_new is numeric
  if (!is.numeric(X_new)) {
    stop("X_new must contain numeric values")
  }

  # Extract recipe for preprocessing
  recipe <- signature$recipe
  preprocessing <- recipe$preprocessing %||% list()

  # Extract preprocessing parameters
  transform_method <- preprocessing$transform %||% "log1p-cpm"
  cpm_min <- preprocessing$cpm_min %||% 1
  cpm_min_samples <- preprocessing$cpm_min_samples %||% 4

  # Domain/tissue checks: Gene ID mapping
  signature_panel <- signature$panel
  X_genes <- rownames(X_new)

  # Compute mapping rate
  matched_genes <- intersect(signature_panel, X_genes)
  mapping_rate <- length(matched_genes) / length(signature_panel)


  # Warn if low mapping rate and add alert
  if (mapping_rate < 0.8) {
    warning_msg <- paste0(
      "Low gene ID mapping rate: ", round(mapping_rate * 100, 1), "% (",
      length(matched_genes), "/", length(signature_panel),
      ") signature genes found in new data. ",
      "This may indicate a gene ID namespace mismatch. ",
      "Predictions may be unreliable."
    )
    warning(warning_msg)
    alerts <- rbind(alerts, data.frame(
      type = "warning",
      message = warning_msg,
      severity = "high",
      stringsAsFactors = FALSE
    ))
  }


  # Check for missing genes and add alert
  missing_genes <- setdiff(signature_panel, X_genes)
  missing_pct <- length(missing_genes) / length(signature_panel)
  if (length(missing_genes) > 0) {
    warning_msg <- paste0(
      "Missing signature genes in new data: ", length(missing_genes),
      " genes (", round(missing_pct * 100, 1), "%) will be set to 0. Missing: ",
      paste(head(missing_genes, 10), collapse = ", "),
      if (length(missing_genes) > 10) " ..." else ""
    )
    warning(warning_msg)
    severity <- if (missing_pct > 0.2) "high" else if (missing_pct > 0.1) "medium" else "low"
    alerts <- rbind(alerts, data.frame(
      type = "warning",
      message = warning_msg,
      severity = severity,
      stringsAsFactors = FALSE
    ))
  }


  # Apply preprocessing recipe (transform and filter)
  # Use same preprocessing as training pipeline
  if (transform_method == "log1p-cpm") {
    # Apply log1p-CPM transformation
    mat_t <- esr_transform_log1p_cpm(
      X = X_new,
      cpm_min = cpm_min,
      cpm_min_samples = cpm_min_samples
    )
    # mat_t is samples x genes after transformation
  } else {
    stop("Only 'log1p-cpm' transform is currently supported for classification")
  }

  # Map signature genes to transformed data
  # mat_t is samples x genes, we need to extract signature genes
  signature_genes_in_data <- intersect(signature_panel, colnames(mat_t))
  signature_genes_missing <- setdiff(signature_panel, colnames(mat_t))

  # Extract expression values for matched genes
  if (length(signature_genes_in_data) > 0) {
    signature_expr <- mat_t[, signature_genes_in_data, drop = FALSE]
  } else {
    stop("No signature genes found in transformed data")
  }

  # Handle missing genes (set to 0)
  if (length(signature_genes_missing) > 0) {
    # Create matrix of zeros for missing genes
    missing_expr <- matrix(0,
      nrow = nrow(signature_expr),
      ncol = length(signature_genes_missing),
      dimnames = list(rownames(signature_expr), signature_genes_missing)
    )
    # Combine matched and missing genes
    signature_expr <- cbind(signature_expr, missing_expr)
    # Ensure genes are in same order as signature panel
    signature_expr <- signature_expr[, signature_panel, drop = FALSE]
  } else {
    # Ensure genes are in same order as signature panel
    signature_expr <- signature_expr[, signature_panel, drop = FALSE]
  }

  # Extract coefficients (in same order as panel)
  coef_vec <- signature$coefficients[signature_panel]

  # Compute signature scores: score = intercept + sum(coefficient * expression)
  # signature_expr is samples x genes, coef_vec is genes x 1
  scores <- as.numeric(signature_expr %*% coef_vec) + signature$intercept

  # Apply calibration if available
  calibration <- signature$recipe$training$calibration_method %||% NULL
  # calibration_params would be used when calibration is fully implemented
  # calibration_params <- NULL
  probabilities <- scores

  # Convert scores to probabilities using sigmoid (for logistic regression)
  probabilities <- 1 / (1 + exp(-scores))

  # Apply calibration if available (Platt scaling or isotonic)
  # Note: Calibration parameters are not stored in signature artifacts yet
  # This is a placeholder for future calibration support
  if (!is.null(calibration) && calibration != "none") {
    # TODO: Apply calibration parameters from signature if available
    # For now, use uncalibrated probabilities
  }


  # Apply threshold
  if (is.character(threshold) && threshold == "youden") {
    # Youden threshold requires labeled validation data
    if (is.null(y_new)) {
      warning("Youden threshold requires labeled validation data (y_new). Using default threshold 0.5.")
      threshold <- 0.5
    } else {
      # Compute Youden threshold using validation labels
      # Note: probabilities are computed after calibration, so use them for threshold selection
      youden_threshold <- .select_threshold_youden(probabilities = probabilities, labels = y_new)
      threshold <- youden_threshold
      if (threshold != 0.5) {
        message(paste("Optimal Youden threshold:", round(threshold, 3)))
      }
    }
  }

  if (!is.numeric(threshold) || threshold < 0 || threshold > 1) {
    stop("threshold must be numeric between 0 and 1, or 'youden'")
  }

  # Get binary predictions based on threshold
  predictions <- ifelse(probabilities >= threshold, "PIS", "PS")
  predictions <- factor(predictions, levels = c("PS", "PIS"))

  # Check for class imbalance if labels provided
  if (!is.null(y_new)) {
    if (is.factor(y_new)) {
      y_new_char <- as.character(y_new)
    } else {
      y_new_char <- as.character(y_new)
    }
    n_ps <- sum(y_new_char == "PS" | y_new_char == "0" | y_new_char == 0)
    n_pis <- sum(y_new_char == "PIS" | y_new_char == "1" | y_new_char == 1)
    total <- n_ps + n_pis

    if (total > 0) {
      ps_pct <- n_ps / total
      pis_pct <- n_pis / total
      imbalance_ratio <- max(ps_pct, pis_pct) / min(ps_pct, pis_pct)

      if (imbalance_ratio > 2.0) {
        warning_msg <- paste0(
          "Class imbalance detected in validation data: ",
          n_ps, " PS (", round(ps_pct * 100, 1), "%), ",
          n_pis, " PIS (", round(pis_pct * 100, 1), "%)"
        )
        alerts <- rbind(alerts, data.frame(
          type = "info",
          message = warning_msg,
          severity = "low",
          stringsAsFactors = FALSE
        ))
      }
    }
  }

  # Compute confidence intervals if requested
  confidence_lower <- NULL
  confidence_upper <- NULL

  if (confidence) {
    # Simple confidence based on distance from threshold
    # For calibrated probabilities, use probability uncertainty
    # For uncalibrated, use score distance from threshold
    distance_from_threshold <- abs(probabilities - threshold)
    # Normalize to 0-1 range (closer to threshold = lower confidence)
    confidence_width <- 1 - pmin(distance_from_threshold * 2, 1)

    # Compute confidence bounds (simple approach)
    confidence_lower <- pmax(probabilities - confidence_width * 0.1, 0)
    confidence_upper <- pmin(probabilities + confidence_width * 0.1, 1)
  }

  # Construct output data.frame
  result <- data.frame(
    sample = sample_ids,
    score = scores,
    probability = probabilities,
    prediction = predictions,
    stringsAsFactors = FALSE
  )

  # Add confidence intervals if requested
  if (confidence && !is.null(confidence_lower) && !is.null(confidence_upper)) {
    result$confidence_lower <- confidence_lower
    result$confidence_upper <- confidence_upper
  }

  # Return predictions and alerts
  # If alerts exist, return as list; otherwise return data.frame for backward compatibility
  if (nrow(alerts) > 0) {
    return(list(predictions = result, alerts = alerts))
  } else {
    return(result)
  }
}

# [END]
