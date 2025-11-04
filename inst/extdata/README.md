# EndoSignatureR extdata — GSE201926 Renamed Sources

This directory contains gzipped source files from GEO/NCBI used by the package. Files were lightly renamed for consistency inside the package.

## Filename mappings (original → packaged)

- GSE201926_raw_counts_GRCh38.p13_NCBI.tsv → gse201926_raw_counts.tsv.gz
- GSE201926_norm_counts_TPM_GRCh38.p13_NCBI.tsv → gse201926_tpm_counts.tsv.gz
- Human.GRCh38.p13.annot.tsv → gse201926_annotation.tsv.gz
- GSE201926_series_matrix.txt → gse201926_series_matrix.txt.gz

Notes

- Files may exist as `.gz`; code should accept either `.gz` or plain.
- Access paths in code via:
  ```r
  system.file("extdata", "gse201926_raw_counts.tsv", package = "endoSignatureR")
  ```
  Prefer `.gz` if present.

## How these are used

- Full workflows: loaders read these sources directly (see `R/data-loading.R`).
- Vignettes/examples: use compact `.rda` datasets in `data/` for speed and offline knitting.

## Provenance and reproducibility

- Builder scripts in `data-raw/` construct minimal `.rda` datasets from these sources.
- See `validation/YYYYMMDD-xxxx-validation.Rmd` for steps to gzip, rebuild, and verify.
