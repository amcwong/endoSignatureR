# Pre-trained Signature Artifacts

This directory contains the pre-trained PS vs PIS signature artifacts trained on the full GSE201926 endometrial lesion RNA-seq dataset.

## Artifacts

- **`endometrial_signature.csv`**: Gene panel with coefficients, selection frequencies, and optional bootstrap frequencies
- **`endometrial_recipe.json`**: Preprocessing recipe, training parameters, signature metadata, and reproducibility information
- **`endometrial_model_card.md`**: Model documentation with provenance, performance metrics, limitations, and intended use
- **`endometrial_stability.csv`** (optional): Bootstrap stability frequencies for signature genes

## Training Parameters

The pre-trained signature was trained using the Phase 2 pipeline (`esr_trainEndometrialSignature()`) with the following parameters:

- **Dataset**: Full GSE201926 (all 12 samples, all genes after filtering)
- **Transform**: log1p-CPM transformation
- **CV Structure**: Leave-Pair-Out (LPO) for small n (n=12) with balanced classes
- **Lambda Selection**: 1se rule (1 standard error rule for sparsity)
- **Aggregation**: Mean aggregation across outer folds
- **Calibration**: Platt scaling for probability calibration
- **Stability Selection**: Bootstrap resampling (500 resamples)

## Access

These artifacts can be accessed via `system.file()`:

```r
# CSV signature
csv_path <- system.file("extdata", "pretrained-signature", "endometrial_signature.csv", package = "endoSignatureR")

# JSON recipe
json_path <- system.file("extdata", "pretrained-signature", "endometrial_recipe.json", package = "endoSignatureR")

# Model card
md_path <- system.file("extdata", "pretrained-signature", "endometrial_model_card.md", package = "endoSignatureR")

# Stability CSV (if available)
stability_path <- system.file("extdata", "pretrained-signature", "endometrial_stability.csv", package = "endoSignatureR")
```

## Usage

The pre-trained signature can be loaded using `esr_loadPretrainedSignature()` (Phase 3.2) and applied to new endometrial samples using `esr_classifyEndometrial()`.

For more information, see the package vignette: `vignette("endoSignatureR-intro", package = "endoSignatureR")`


