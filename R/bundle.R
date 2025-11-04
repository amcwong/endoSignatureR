#' Compute QC Metrics for Endometrial Data
#'
#' Computes quality control metrics including library sizes, percentage of zeros,
#' and filtering statistics for endometrial RNA-seq data.
#'
#' @param counts A matrix/data.frame of raw counts (genes × samples).
#' @param mat_t Optional transformed matrix (samples × genes). If provided, used
#'   to compute filtering statistics.
#' @param pheno Optional data.frame with sample metadata; used for per-group
#'   statistics if provided.
#' @param group_col Character scalar; name of the phenotype column for grouping.
#'   Defaults to "group".
#'
#' @return A data.frame with columns:
#' \describe{
#'   \item{sample_id}{Character; sample identifiers}
#'   \item{library_size}{Numeric; total read counts per sample}
#'   \item{pct_zeros}{Numeric; percentage of zero counts per sample}
#'   \item{n_genes}{Numeric; number of genes (non-zero) per sample}
#' }
#' If `mat_t` is provided, additional columns:
#' \describe{
#'   \item{genes_after_filter}{Numeric; number of genes after filtering}
#' }
#' If `pheno` is provided, the data.frame includes the group column.
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' qc_metrics <- esr_computeQCMetrics(
#'   counts = gse201926_sample$counts,
#'   mat_t = mat_t,
#'   pheno = gse201926_sample$pheno
#' )
#' head(qc_metrics)
#' @importFrom utils head
#' @export
esr_computeQCMetrics <- function(counts, mat_t = NULL, pheno = NULL, group_col = "group") {
    # Convert counts to matrix if needed
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }

    # Validate counts dimensions
    if (nrow(counts) == 0 || ncol(counts) == 0) {
        stop("counts must have at least one gene and one sample")
    }

    # Get sample IDs from column names or create them
    if (!is.null(colnames(counts))) {
        sample_ids <- colnames(counts)
    } else {
        sample_ids <- paste0("sample_", seq_len(ncol(counts)))
    }

    # Calculate library sizes (column sums)
    lib_sizes <- colSums(counts, na.rm = TRUE)

    # Calculate percentage of zeros per sample
    n_zeros <- colSums(counts == 0, na.rm = TRUE)
    n_total <- nrow(counts)
    pct_zeros <- (n_zeros / n_total) * 100

    # Calculate number of genes (non-zero) per sample
    n_genes <- n_total - n_zeros

    # Create base QC metrics data.frame
    qc_metrics <- data.frame(
        sample_id = sample_ids,
        library_size = lib_sizes,
        pct_zeros = pct_zeros,
        n_genes = n_genes,
        stringsAsFactors = FALSE
    )

    # Add filtering statistics if mat_t provided
    if (!is.null(mat_t)) {
        if (!is.matrix(mat_t) && !is.data.frame(mat_t)) {
            stop("mat_t must be a matrix or data.frame")
        }
        if (!is.matrix(mat_t)) {
            mat_t <- as.matrix(mat_t)
        }

        # Get genes after filtering (columns in mat_t)
        genes_after_filter <- ncol(mat_t)
        qc_metrics$genes_after_filter <- genes_after_filter

        # Also compute per-sample filtering stats if sample IDs match
        if (!is.null(rownames(mat_t))) {
            # Match samples between counts and mat_t
            common_samples <- intersect(sample_ids, rownames(mat_t))
            if (length(common_samples) > 0) {
                # For matched samples, genes_after_filter is the same
                # This is already added above as a constant
                # Could compute per-sample if needed, but typically filtering is gene-level
            }
        }
    }

    # Add group information if pheno provided
    if (!is.null(pheno)) {
        if (!is.data.frame(pheno)) {
            stop("pheno must be a data.frame")
        }

        # Check if sample_id column exists
        if ("sample_id" %in% names(pheno)) {
            # Merge with pheno to add group
            qc_metrics <- merge(qc_metrics, pheno[, c("sample_id", group_col)], 
                                by = "sample_id", all.x = TRUE)
            # Rename group column for consistency
            if (group_col != "group" && group_col %in% names(qc_metrics)) {
                qc_metrics$group <- qc_metrics[[group_col]]
            }
        } else if (group_col %in% names(pheno)) {
            # If no sample_id but group_col exists, try to match by position
            if (nrow(pheno) == nrow(qc_metrics)) {
                qc_metrics$group <- pheno[[group_col]]
            } else {
                warning("pheno row count doesn't match qc_metrics; skipping group merge")
            }
        }
    }

    # Sort by sample_id for consistency
    qc_metrics <- qc_metrics[order(qc_metrics$sample_id), ]
    rownames(qc_metrics) <- NULL

    return(qc_metrics)
}

#' Create Analysis Bundle for Mode 3 Workflow
#'
#' Constructs an analysis bundle that packages all Mode 3 (Standalone Visualization
#' & Analysis) outputs together. The bundle provides a unified structure for
#' passing data between functions, exporting results, and integrating with Shiny.
#'
#' @param counts_t A transformed matrix (samples × genes). Required.
#' @param de_table Optional data.frame from `esr_analyzeDifferentialExpression()`
#'   with columns: `gene_id`, `log2FC`, `FDR`, etc.
#' @param selected_genes Optional character vector of selected gene IDs.
#' @param qc_metrics Optional data.frame from `esr_computeQCMetrics()` with QC
#'   metrics per sample.
#' @param pheno Optional data.frame with sample metadata.
#' @param annot Optional data.frame with gene annotations.
#' @param transform_params Optional list of transformation parameters used.
#' @param de_params Optional list of DE analysis parameters used.
#' @param validation_issues Optional data.frame with validation issues (from
#'   `esr_validateEndometrial()`).
#' @param ... Additional optional fields to include in bundle.
#'
#' @return A list (analysis bundle) with components:
#' \describe{
#'   \item{counts_t}{Transformed matrix (samples × genes)}
#'   \item{de_table}{DE table (data.frame or NULL)}
#'   \item{selected_genes}{Selected gene IDs (character vector or NULL)}
#'   \item{qc_metrics}{QC metrics (data.frame or NULL)}
#'   \item{pheno}{Sample metadata (data.frame or NULL)}
#'   \item{annot}{Gene annotations (data.frame or NULL)}
#'   \item{transform_params}{Transformation parameters (list or NULL)}
#'   \item{de_params}{DE analysis parameters (list or NULL)}
#'   \item{validation_issues}{Validation issues (data.frame or NULL)}
#' }
#' Additional fields may be included via `...`.
#'
#' @details
#' The bundle validates component alignment:
#' - Sample IDs: `counts_t` rownames, `qc_metrics$sample_id`, and `pheno$sample_id`
#'   must align if all provided.
#' - Gene IDs: `counts_t` colnames, `de_table$gene_id`, and `selected_genes` must
#'   align (subject to filtering).
#'
#' Warnings are issued for mismatches, but bundle creation proceeds.
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
#' top_genes <- esr_selectTopGenes(de_table = de_table, n = 10, by = "de")
#' qc_metrics <- esr_computeQCMetrics(
#'   counts = gse201926_sample$counts,
#'   mat_t = mat_t,
#'   pheno = gse201926_sample$pheno
#' )
#'
#' bundle <- esr_createAnalysisBundle(
#'   counts_t = mat_t,
#'   de_table = de_table,
#'   selected_genes = top_genes,
#'   qc_metrics = qc_metrics,
#'   pheno = gse201926_sample$pheno,
#'   annot = gse201926_sample$annot
#' )
#'
#' # Access bundle components
#' bundle$counts_t
#' bundle$de_table
#' bundle$selected_genes
#' @importFrom utils head
#' @export
esr_createAnalysisBundle <- function(counts_t, de_table = NULL, selected_genes = NULL,
                                     qc_metrics = NULL, pheno = NULL, annot = NULL,
                                     transform_params = NULL, de_params = NULL,
                                     validation_issues = NULL, ...) {
    # Validate required input
    if (missing(counts_t) || is.null(counts_t)) {
        stop("counts_t is required")
    }

    # Convert counts_t to matrix if needed
    if (!is.matrix(counts_t)) {
        if (is.data.frame(counts_t)) {
            counts_t <- as.matrix(counts_t)
        } else {
            stop("counts_t must be a matrix or data.frame")
        }
    }

    # Validate counts_t dimensions
    if (nrow(counts_t) == 0 || ncol(counts_t) == 0) {
        stop("counts_t must have at least one sample and one gene")
    }

    # Get sample IDs from counts_t rownames
    sample_ids_counts <- rownames(counts_t)
    if (is.null(sample_ids_counts)) {
        sample_ids_counts <- paste0("sample_", seq_len(nrow(counts_t)))
        rownames(counts_t) <- sample_ids_counts
    }

    # Get gene IDs from counts_t colnames
    gene_ids_counts <- colnames(counts_t)
    if (is.null(gene_ids_counts)) {
        gene_ids_counts <- paste0("gene_", seq_len(ncol(counts_t)))
        colnames(counts_t) <- gene_ids_counts
    }

    # Validate and align sample IDs
    if (!is.null(qc_metrics)) {
        if (!is.data.frame(qc_metrics)) {
            stop("qc_metrics must be a data.frame")
        }
        if ("sample_id" %in% names(qc_metrics)) {
            qc_sample_ids <- qc_metrics$sample_id
            common_samples_qc <- intersect(sample_ids_counts, qc_sample_ids)
            if (length(common_samples_qc) == 0) {
                warning("No matching sample IDs found between counts_t rownames and qc_metrics$sample_id")
            } else if (length(common_samples_qc) < length(sample_ids_counts) || 
                       length(common_samples_qc) < length(qc_sample_ids)) {
                warning("Some sample IDs don't match between counts_t and qc_metrics")
            }
        }
    }

    if (!is.null(pheno)) {
        if (!is.data.frame(pheno)) {
            stop("pheno must be a data.frame")
        }
        if ("sample_id" %in% names(pheno)) {
            pheno_sample_ids <- pheno$sample_id
            common_samples_pheno <- intersect(sample_ids_counts, pheno_sample_ids)
            if (length(common_samples_pheno) == 0) {
                warning("No matching sample IDs found between counts_t rownames and pheno$sample_id")
            } else if (length(common_samples_pheno) < length(sample_ids_counts) || 
                       length(common_samples_pheno) < length(pheno_sample_ids)) {
                warning("Some sample IDs don't match between counts_t and pheno")
            }
        }
    }

    # Validate and align gene IDs
    if (!is.null(de_table)) {
        if (!is.data.frame(de_table)) {
            stop("de_table must be a data.frame")
        }
        if ("gene_id" %in% names(de_table)) {
            de_gene_ids <- de_table$gene_id
            common_genes_de <- intersect(gene_ids_counts, de_gene_ids)
            if (length(common_genes_de) == 0) {
                warning("No matching gene IDs found between counts_t colnames and de_table$gene_id")
            } else if (length(common_genes_de) < length(de_gene_ids)) {
                # It's OK if de_table has more genes (before filtering), but warn if none match
                if (length(common_genes_de) == 0) {
                    warning("No genes from de_table match counts_t colnames (may indicate filtering)")
                }
            }
        } else {
            warning("de_table provided but doesn't contain 'gene_id' column")
        }
    }

    if (!is.null(selected_genes)) {
        if (!is.character(selected_genes)) {
            stop("selected_genes must be a character vector")
        }
        common_genes_selected <- intersect(gene_ids_counts, selected_genes)
        if (length(common_genes_selected) == 0) {
            warning("No matching gene IDs found between counts_t colnames and selected_genes")
        } else if (length(common_genes_selected) < length(selected_genes)) {
            missing_genes <- setdiff(selected_genes, gene_ids_counts)
            warning(paste0("Some selected genes not found in counts_t: ", 
                          paste(head(missing_genes, 5), collapse = ", "),
                          if (length(missing_genes) > 5) " ..." else ""))
        }
    }

    # Create bundle list
    bundle <- list(
        counts_t = counts_t,
        de_table = de_table,
        selected_genes = selected_genes,
        qc_metrics = qc_metrics,
        pheno = pheno,
        annot = annot,
        transform_params = transform_params,
        de_params = de_params,
        validation_issues = validation_issues
    )

    # Add additional fields from ...
    additional_fields <- list(...)
    if (length(additional_fields) > 0) {
        bundle <- c(bundle, additional_fields)
    }

    # Add class for potential S3 methods in future
    class(bundle) <- c("esr_analysis_bundle", "list")

    return(bundle)
}

# [END]

