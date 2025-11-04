# Endometrial Signature Model Card

## Model Details
- **Model Name**: Endometrial PS vs PIS Signature
- **Version**: 1.0.0
- **Date**: 2025-11-03
- **Package**: endoSignatureR v0.0.0.9000
- **R Version**: R version 4.5.1 (2025-06-13)

## Model Purpose
This signature classifies endometrial bulk RNA-seq samples as Proliferative Secretory (PS) or Proliferative Inflammatory Secretory (PIS) using LASSO logistic regression with nested cross-validation.

## Training Data
- **Samples**: 12
- **Genes**: 3 (signature panel)
- **Class Balance**: Balanced (6 PS, 6 PIS)

## Model Architecture
- **Algorithm**: LASSO logistic regression (glmnet)
- **Preprocessing**: log1p-cpm transformation, CPM filtering (min=1, min_samples=4), DE-based top-K selection (K=300)
- **CV Structure**: Nested CV (outer: 3-fold, inner: 5-fold with 10 repeats)
- **Lambda Selection**: 1se rule (1 standard error rule for sparsity)
- **Coefficient Aggregation**: mean aggregation across outer folds (Option 2)
- **Consensus Genes**: Genes selected in ≥2 outer folds

## Performance Metrics
- **AUC-ROC**: 0.639
- **Accuracy**: 0.75
- **Brier Score**: 0.303
- **Expected Calibration Error (ECE)**: 0.316
- **Signature Size**: 3 genes

## Calibration
- **Method**: platt scaling
- **Calibrated Probabilities**: Included in predictions

## Stability Selection
- **Bootstrap Resamples**: 100
- **Stability Threshold**: 0.7
- **Stable Genes**: 0

## Limitations
- **Tissue Specificity**: Trained on endometrial tissue only; not portable across tissues
- **Small n Uncertainty**: Limited sample size introduces uncertainty; avoid clinical claims
- **Binary Only**: Supports only PS vs PIS classification; multiclass/continuous outcomes out of scope
- **Batch Effects**: Users should provide batch information and use `~ batch + group` designs

## Intended Use
- **Use Case**: Classification of endometrial bulk RNA-seq samples into PS vs PIS states
- **Target Audience**: Endometrial researchers, pathologists, clinicians
- **Not Intended For**: Non-endometrial tissues, clinical decision-making without validation

## Training Parameters
- **Outer CV**: 3-fold (stratified)
- **Inner CV**: 5-fold with 10 repeats
- **Lambda Rule**: 1se
- **Top-K Selection**: 300 genes
- **Aggregation Method**: mean
- **Calibration Method**: platt
- **Seed**: 12345

## Reproducibility
- **Seeds**: Outer CV seed=12345, Inner CV seed=67890
- **CV Folds**: Pre-computed or generated deterministically
- **Package Versions**: endoSignatureR v0.0.0.9000, glmnet, limma, ...

## References
- LASSO: Tibshirani (1996), "Regression shrinkage and selection via the lasso"
- Nested CV: Varma & Simon (2006), "Bias in error estimation when using cross-validation for model selection"

