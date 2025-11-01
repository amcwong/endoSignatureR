## code to prepare `gse201926_sample` dataset

library(readr)
library(dplyr)

# Use endo_parse_metadata from package if available, otherwise source the R file
# This avoids conflicts when package is loaded via devtools::load_all()
if (exists("endo_parse_metadata", envir = asNamespace("endoSignatureR"), mode = "function")) {
    # Package is loaded, use namespace function
    endo_parse_metadata <- get("endo_parse_metadata", envir = asNamespace("endoSignatureR"))
} else if (!exists("endo_parse_metadata", envir = .GlobalEnv)) {
    # Package not loaded and function not in global env, source the R file
    source("R/data-loading.R")
}

# Load full raw counts matrix from packaged extdata (supports .gz)
pick_ext <- function(name) {
    gz <- file.path("inst", "extdata", paste0(name, ".gz"))
    plain <- file.path("inst", "extdata", name)
    if (file.exists(gz)) {
        return(gz)
    }
    if (file.exists(plain)) {
        return(plain)
    }
    stop("Missing extdata file: ", name)
}

# Load full raw counts matrix
raw_counts_path <- pick_ext("gse201926_raw_counts.tsv")
raw_counts <- readr::read_tsv(raw_counts_path, show_col_types = FALSE)

# Load real phenotype data from series matrix using package parser
series_matrix_path <- pick_ext("gse201926_series_matrix.txt")
pheno_data <- endo_parse_metadata(series_matrix_path)

# Create full counts matrix (all genes, no filtering)
counts_matrix <- as.matrix(raw_counts[, -1])
rownames(counts_matrix) <- raw_counts$GeneID

# Set column names to match phenotype sample_ids exactly
# The raw counts file column names should correspond to phenotype sample_ids
# We need to align them so they match for downstream analysis
raw_colnames <- colnames(raw_counts)[-1]

# Clean sample IDs from both sources for matching
raw_colnames_clean <- trimws(gsub('^"|"$', "", raw_colnames))
pheno_sample_ids_clean <- trimws(gsub('^"|"$', "", pheno_data$sample_id))

# Match counts matrix columns to phenotype sample_ids
# First, check if they contain the same IDs (order may differ)
if (setequal(raw_colnames_clean, pheno_sample_ids_clean)) {
    # Reorder counts matrix columns to match phenotype sample_id order
    # Create mapping from cleaned raw colnames to phenotype IDs
    matched_indices <- match(pheno_sample_ids_clean, raw_colnames_clean)
    if (any(is.na(matched_indices))) {
        stop("Could not match all phenotype sample_ids to counts matrix columns")
    }

    # Reorder and rename columns
    counts_matrix <- counts_matrix[, matched_indices, drop = FALSE]
    colnames(counts_matrix) <- pheno_data$sample_id # Use phenotype IDs exactly
} else {
    # IDs don't match - try to match what we can
    common_ids <- intersect(raw_colnames_clean, pheno_sample_ids_clean)
    if (length(common_ids) == 0) {
        stop("No matching sample IDs between counts matrix and phenotype data")
    }

    # Match what we can
    matched_indices <- match(common_ids, raw_colnames_clean)
    counts_matrix <- counts_matrix[, matched_indices, drop = FALSE]
    colnames(counts_matrix) <- pheno_data$sample_id[pheno_data$sample_id %in% common_ids]

    warning(
        "Only ", length(common_ids), " out of ", length(pheno_sample_ids_clean),
        " sample IDs matched between counts and phenotype"
    )
}

# Load full gene annotation (no filtering)
annot_path <- pick_ext("gse201926_annotation.tsv")
annot <- readr::read_tsv(annot_path, show_col_types = FALSE)

# Create final data object with full dataset
gse201926_sample <- list(
    counts = counts_matrix,
    pheno = pheno_data,
    annot = annot
)

# Save as package data
usethis::use_data(gse201926_sample, overwrite = TRUE, compress = "xz")
