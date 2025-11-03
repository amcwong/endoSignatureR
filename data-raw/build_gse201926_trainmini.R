## code to prepare `gse201926_trainmini` dataset
## Training subset: 800-1500 genes × 12 samples (all samples, balanced labels)
## Derived from full GSE201926 data in inst/extdata

suppressPackageStartupMessages({
    library(readr)
    library(dplyr)
    library(usethis)
})

# Use endo_load_gse201926 from package namespace (avoid creating global copies)
# Always prefer package namespace to avoid conflicts with devtools::load_all()
if (exists("endo_load_gse201926", envir = asNamespace("endoSignatureR"), mode = "function")) {
    # Use function from package namespace directly (no global assignment)
    load_fn <- get("endo_load_gse201926", envir = asNamespace("endoSignatureR"))
} else {
    # If package not loaded, source into local environment
    source_env <- new.env()
    source("R/data-loading.R", local = source_env)
    # Use function from local environment (no global assignment)
    load_fn <- get("endo_load_gse201926", envir = source_env)
}

# Load full dataset from inst/extdata
message("Loading full GSE201926 dataset from inst/extdata...")
full_data <- load_fn(sample_only = FALSE)

# Clean up: remove any global copies that might mask package namespace
# This prevents conflicts when running devtools::document() or devtools::load_all()
if (exists("endo_load_gse201926", envir = .GlobalEnv, inherits = FALSE)) {
    rm(list = "endo_load_gse201926", envir = .GlobalEnv)
}

# Extract components
counts_full <- full_data$counts
pheno_full <- full_data$pheno
annot_full <- full_data$annot

# Convert counts to matrix (genes x samples)
# Handle both data.frame and matrix input
if (is.data.frame(counts_full)) {
    # If data.frame, GeneID column becomes rownames, other columns are samples
    if ("GeneID" %in% names(counts_full)) {
        gene_ids <- counts_full$GeneID
        # Remove GeneID column - only keep sample columns
        sample_cols <- names(counts_full)[!names(counts_full) %in% "GeneID"]
        counts_matrix <- as.matrix(counts_full[, sample_cols, drop = FALSE])
        rownames(counts_matrix) <- gene_ids
        # Set column names to sample columns (clean them)
        colnames(counts_matrix) <- trimws(gsub('^"|"$', "", sample_cols))
    } else {
        # No GeneID column, check if first column should be rownames
        if (is.null(rownames(counts_full)) && ncol(counts_full) > 0) {
            # Try using first column as rownames if it looks like gene IDs
            first_col <- counts_full[, 1, drop = TRUE]
            if (is.character(first_col) || is.numeric(first_col)) {
                gene_ids <- first_col
                counts_matrix <- as.matrix(counts_full[, -1, drop = FALSE])
                rownames(counts_matrix) <- gene_ids
                colnames(counts_matrix) <- trimws(gsub('^"|"$', "", names(counts_full)[-1]))
            } else {
                stop("Cannot determine gene IDs: counts data.frame has no GeneID column and first column is not suitable for rownames")
            }
        } else if (!is.null(rownames(counts_full))) {
            counts_matrix <- as.matrix(counts_full)
            colnames(counts_matrix) <- trimws(gsub('^"|"$', "", colnames(counts_matrix)))
        } else {
            stop("Cannot determine gene IDs: counts data.frame has no GeneID column and no rownames")
        }
    }
} else if (is.matrix(counts_full)) {
    counts_matrix <- counts_full
    # If matrix has GeneID in column names, remove it
    if (!is.null(colnames(counts_matrix)) && "GeneID" %in% colnames(counts_matrix)) {
        geneid_col_idx <- which(colnames(counts_matrix) == "GeneID")
        if (length(geneid_col_idx) == 1 && is.null(rownames(counts_matrix))) {
            # Use GeneID column as rownames, then remove it
            rownames(counts_matrix) <- counts_matrix[, geneid_col_idx]
            counts_matrix <- counts_matrix[, -geneid_col_idx, drop = FALSE]
        } else {
            # Just remove the GeneID column
            counts_matrix <- counts_matrix[, -geneid_col_idx, drop = FALSE]
        }
    }
    # Clean column names
    if (!is.null(colnames(counts_matrix))) {
        colnames(counts_matrix) <- trimws(gsub('^"|"$', "", colnames(counts_matrix)))
    }
} else {
    stop("counts must be a matrix or data.frame")
}

# Ensure we have rownames (GeneID) and no GeneID in column names
if (is.null(rownames(counts_matrix))) {
    stop("Counts matrix must have rownames (GeneID)")
}

# Final check: remove GeneID from column names if it somehow got there
if (!is.null(colnames(counts_matrix)) && "GeneID" %in% colnames(counts_matrix)) {
    geneid_idx <- which(colnames(counts_matrix) == "GeneID")
    counts_matrix <- counts_matrix[, -geneid_idx, drop = FALSE]
    warning("Removed GeneID from column names - GeneID should only be rownames")
}

# Clean phenotype sample_id values (remove quotes if present)
# The series matrix parser may include quotes as part of the string
pheno_sample_ids_clean <- trimws(gsub('^"|"$', "", pheno_full$sample_id))

# Ensure counts matrix column names match phenotype sample_id
# This is critical for rsample splits to work correctly
counts_colnames <- colnames(counts_matrix)

# Clean counts column names (remove quotes and trim whitespace)
counts_colnames_clean <- trimws(gsub('^"|"$', "", counts_colnames))

# Remove "GeneID" from column names if present (it should be rownames, not a column)
# This should have been handled above, but double-check
if (!is.null(colnames(counts_matrix)) && "GeneID" %in% colnames(counts_matrix)) {
    geneid_idx <- which(colnames(counts_matrix) == "GeneID")
    counts_matrix <- counts_matrix[, -geneid_idx, drop = FALSE]
    message("Removed GeneID from column names (GeneID should only be rownames)")
}

if (is.null(counts_colnames) || length(counts_colnames) == 0) {
    # No column names - set them from phenotype
    if (length(pheno_sample_ids_clean) == ncol(counts_matrix)) {
        colnames(counts_matrix) <- pheno_sample_ids_clean
        message("Set counts column names from phenotype sample_id (cleaned)")
    } else {
        stop(
            "Cannot set column names: phenotype has ", length(pheno_sample_ids_clean),
            " samples but counts has ", ncol(counts_matrix), " columns"
        )
    }
} else {
    # Column names exist - clean and align with phenotype
    # Match cleaned phenotype IDs to cleaned column names
    matched_indices <- match(pheno_sample_ids_clean, counts_colnames_clean)

    if (any(is.na(matched_indices))) {
        # Some IDs don't match - try without cleaning first
        matched_indices_raw <- match(pheno_full$sample_id, counts_colnames)
        if (!any(is.na(matched_indices_raw))) {
            # Raw match works - use raw IDs
            matched_indices <- matched_indices_raw
            pheno_sample_ids_clean <- pheno_full$sample_id
            counts_colnames_clean <- counts_colnames
        } else {
            # Neither works - detailed error
            missing_ids <- pheno_sample_ids_clean[is.na(matched_indices)]
            stop(
                "Cannot align counts column names with phenotype sample_id.\n",
                "Missing sample IDs: ", paste(head(missing_ids, 5), collapse = ", "), "\n",
                "Available column names: ", paste(head(counts_colnames_clean, 5), collapse = ", "), "\n",
                "Phenotype sample_id (raw): ", paste(head(pheno_full$sample_id, 5), collapse = ", "), "\n",
                "Counts columns (raw): ", paste(head(counts_colnames, 5), collapse = ", ")
            )
        }
    }

    # Reorder and rename to match phenotype exactly (use cleaned IDs)
    counts_matrix <- counts_matrix[, matched_indices, drop = FALSE]
    colnames(counts_matrix) <- pheno_sample_ids_clean
    message("Aligned counts column names with phenotype sample_id (cleaned quotes)")

    # Also update phenotype to use cleaned sample_id
    pheno_full$sample_id <- pheno_sample_ids_clean
}

counts_full <- counts_matrix

# Use all 12 samples (already balanced: 6 PS + 6 PIS)
# For training mini dataset, we subset genes to 800-1500 for faster nested CV
message("Selecting subset of genes for training mini dataset...")

# Strategy: Select most variable genes (or protein-coding if available)
# Set seed for reproducibility
set.seed(12345)

# Calculate gene variance across samples
gene_vars <- apply(counts_full, 1, var, na.rm = TRUE)
gene_vars <- sort(gene_vars, decreasing = TRUE)

# Select top variable genes (target: 800-1500 genes)
# Use 1200 as target (middle of range) for good representation
target_n_genes <- 1200
n_available <- length(gene_vars)
n_select <- min(target_n_genes, n_available)

selected_genes <- names(gene_vars)[1:n_select]

# Alternative: If annotation has GeneType, prefer protein-coding genes
if (!is.null(annot_full) && "GeneType" %in% names(annot_full)) {
    pc_genes <- annot_full$GeneID[annot_full$GeneType == "protein_coding"]
    pc_genes <- intersect(pc_genes, rownames(counts_full))

    if (length(pc_genes) >= 800) {
        # Use protein-coding genes, sorted by variance
        pc_vars <- gene_vars[names(gene_vars) %in% pc_genes]
        pc_vars <- sort(pc_vars, decreasing = TRUE)
        n_pc_select <- min(target_n_genes, length(pc_vars))
        selected_genes <- names(pc_vars)[1:n_pc_select]
        message("Selected ", length(selected_genes), " protein-coding genes (most variable)")
    } else {
        message("Selected ", length(selected_genes), " most variable genes (protein-coding filter insufficient)")
    }
} else {
    message("Selected ", length(selected_genes), " most variable genes")
}

# Subset counts to selected genes
counts_trainmini <- counts_full[selected_genes, , drop = FALSE]

# Subset annotation to selected genes
if (!is.null(annot_full) && "GeneID" %in% names(annot_full)) {
    annot_trainmini <- annot_full[annot_full$GeneID %in% selected_genes, , drop = FALSE]
} else {
    annot_trainmini <- annot_full
}

# Verify class balance in phenotype
table(pheno_full$group)
message("Class distribution: ", paste(table(pheno_full$group), collapse = " vs "))

# Create final training mini dataset
gse201926_trainmini <- list(
    counts = counts_trainmini,
    pheno = pheno_full,
    annot = annot_trainmini
)

message("Created gse201926_trainmini:")
message("  - Counts: ", nrow(counts_trainmini), " genes × ", ncol(counts_trainmini), " samples")
message("  - Pheno: ", nrow(pheno_full), " samples (", sum(pheno_full$group == "PS"), " PS, ", sum(pheno_full$group == "PIS"), " PIS)")
message("  - Annot: ", nrow(annot_trainmini), " genes")

# Save as package data
usethis::use_data(gse201926_trainmini, overwrite = TRUE, compress = "xz")
message("Saved gse201926_trainmini to data/")
