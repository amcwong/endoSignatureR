#' Plot Endometrial PCA
#'
#' Produces a PCA plot for exploratory analysis with group coloring.
#'
#' @param mat_t A samples x genes transformed matrix.
#' @param pheno Optional data.frame with sample metadata containing a `group` column for coloring.
#' @param group_col Character scalar; name of the phenotype column for group coloring. Defaults to "group".
#'
#' @return A ggplot object showing PC1 vs PC2 colored by group.
#'
#' @importFrom utils head
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' plotEndometrialPCA(mat_t, pheno = gse201926_sample$pheno)
#' @export
plotEndometrialPCA <- function(mat_t, pheno = NULL, group_col = "group") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Compute PCA
    set.seed(123) # For reproducibility
    pca_result <- stats::prcomp(mat_t, center = TRUE, scale. = TRUE)

    # Extract PC1 and PC2
    pca_df <- data.frame(
        PC1 = pca_result$x[, 1],
        PC2 = pca_result$x[, 2],
        sample_id = rownames(pca_result$x),
        stringsAsFactors = FALSE
    )

    # Add group information if provided
    if (!is.null(pheno) && group_col %in% names(pheno)) {
        if ("sample_id" %in% names(pheno)) {
            pca_df <- merge(pca_df, pheno[, c("sample_id", group_col)], by = "sample_id", all.x = TRUE)
            pca_df$group <- pca_df[[group_col]]
        } else {
            pca_df$group <- factor(1)
        }
    } else {
        pca_df$group <- factor(1)
    }

    # Calculate variance explained
    var_explained <- summary(pca_result)$importance[2, c(1, 2)]

    # Count samples
    n_samples <- nrow(pca_df)
    n_groups <- length(unique(pca_df$group))

    # Create plot with better point distinction
    p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = .data$PC1, y = .data$PC2, color = .data$group, shape = .data$group)) +
        ggplot2::geom_point(size = 4, alpha = 0.8, stroke = 1.2) +
        ggplot2::scale_shape_manual(values = c("PIS" = 16, "PS" = 17, "1" = 16)) + # Different shapes for groups
        ggplot2::labs(
            x = paste0("PC1 (", round(var_explained[1] * 100, 1), "%)"),
            y = paste0("PC2 (", round(var_explained[2] * 100, 1), "%)"),
            title = paste0("PCA Plot (n = ", n_samples, " samples)"),
            color = "Group",
            shape = "Group"
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        # Add sample count annotation
        ggplot2::annotate(
            "text",
            x = Inf, y = Inf,
            label = paste0("n = ", n_samples, " samples\n", n_groups, " groups"),
            hjust = 1.1, vjust = 1.5,
            size = 3.5,
            color = "gray40"
        )

    return(p)
}

#' Plot Endometrial Library Sizes
#'
#' Visualizes library size (total read counts) distribution across samples.
#'
#' @param counts A matrix/data.frame of raw counts (genes x samples).
#' @param pheno Optional data.frame with sample metadata containing a `group` column for grouping.
#' @param group_col Character scalar; name of the phenotype column for grouping. Defaults to "group".
#'
#' @return A ggplot object showing library size distribution (histogram or boxplot by group).
#'
#' @examples
#' data(gse201926_sample)
#' plotEndometrialLibsize(gse201926_sample$counts, pheno = gse201926_sample$pheno)
#' @export
plotEndometrialLibsize <- function(counts, pheno = NULL, group_col = "group") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Convert to matrix if needed
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }

    # Calculate library sizes (column sums)
    lib_sizes <- colSums(counts, na.rm = TRUE)

    # Create data frame
    lib_df <- data.frame(
        sample_id = names(lib_sizes),
        library_size = lib_sizes,
        stringsAsFactors = FALSE
    )

    # Add group information if provided
    if (!is.null(pheno) && group_col %in% names(pheno)) {
        if ("sample_id" %in% names(pheno)) {
            lib_df <- merge(lib_df, pheno[, c("sample_id", group_col)], by = "sample_id", all.x = TRUE)
            lib_df$group <- lib_df[[group_col]]
        } else {
            lib_df$group <- factor(1)
        }
    } else {
        lib_df$group <- factor(1)
    }

    # Count samples
    n_samples <- nrow(lib_df)
    n_groups <- if (length(unique(lib_df$group)) > 1) length(unique(lib_df$group)) else 1

    # Create plot
    if (length(unique(lib_df$group)) > 1) {
        # Boxplot if grouping available
        p <- ggplot2::ggplot(lib_df, ggplot2::aes(x = .data$group, y = .data$library_size, fill = .data$group)) +
            ggplot2::geom_boxplot(alpha = 0.7, outlier.shape = NA) +
            ggplot2::geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.6, color = "black") +
            ggplot2::labs(
                x = "Group",
                y = "Library Size (total counts)",
                title = paste0("Library Size Distribution by Group (n = ", n_samples, " samples)"),
                fill = "Group"
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(
                plot.title = ggplot2::element_text(hjust = 0.5),
                legend.position = "right"
            ) +
            # Add sample count annotation
            ggplot2::annotate(
                "text",
                x = Inf, y = Inf,
                label = paste0("n = ", n_samples, " samples"),
                hjust = 1.1, vjust = 1.5,
                size = 3.5,
                color = "gray40"
            )
    } else {
        # Histogram if no grouping
        p <- ggplot2::ggplot(lib_df, ggplot2::aes(x = .data$library_size)) +
            ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
            ggplot2::labs(
                x = "Library Size (total counts)",
                y = "Frequency",
                title = paste0("Library Size Distribution (n = ", n_samples, " samples)")
            ) +
            ggplot2::theme_bw() +
            ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
            # Add sample count annotation
            ggplot2::annotate(
                "text",
                x = Inf, y = Inf,
                label = paste0("n = ", n_samples, " samples"),
                hjust = 1.1, vjust = 1.5,
                size = 3.5,
                color = "gray40"
            )
    }

    return(p)
}

#' Plot Endometrial Zero Counts
#'
#' Visualizes the percentage of zero counts per gene or per sample.
#'
#' @param counts A matrix/data.frame of raw counts (genes x samples).
#' @param by Character scalar; "gene" (default) or "sample" to compute zeros per gene or per sample.
#'
#' @return A ggplot object showing zero count percentage distribution.
#'
#' @examples
#' data(gse201926_sample)
#' plotEndometrialZeros(gse201926_sample$counts, by = "gene")
#' plotEndometrialZeros(gse201926_sample$counts, by = "sample")
#' @export
plotEndometrialZeros <- function(counts, by = c("gene", "sample")) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    by <- match.arg(by)

    # Convert to matrix if needed
    if (!is.matrix(counts)) {
        counts <- as.matrix(counts)
    }

    # Calculate zero percentages
    if (by == "gene") {
        # Zeros per gene (across samples)
        n_samples <- ncol(counts)
        zeros_per_gene <- rowSums(counts == 0, na.rm = TRUE)
        pct_zeros <- (zeros_per_gene / n_samples) * 100

        zero_df <- data.frame(
            gene_id = names(pct_zeros),
            pct_zeros = pct_zeros,
            stringsAsFactors = FALSE
        )

        x_label <- "Percentage of Samples with Zero Counts"
        title_text <- "Zero Count Distribution per Gene"
    } else {
        # Zeros per sample (across genes)
        n_genes <- nrow(counts)
        zeros_per_sample <- colSums(counts == 0, na.rm = TRUE)
        pct_zeros <- (zeros_per_sample / n_genes) * 100

        zero_df <- data.frame(
            sample_id = names(pct_zeros),
            pct_zeros = pct_zeros,
            stringsAsFactors = FALSE
        )

        x_label <- "Percentage of Genes with Zero Counts"
        title_text <- "Zero Count Distribution per Sample"
    }

    # Count samples or genes
    n_items <- nrow(zero_df)
    item_type <- ifelse(by == "gene", "genes", "samples")

    # Calculate summary statistics for annotation
    mean_zeros <- round(mean(pct_zeros), 2)
    min_zeros <- round(min(pct_zeros), 2)
    max_zeros <- round(max(pct_zeros), 2)
    perfect_coverage <- sum(pct_zeros == 0)

    # Create plot with clearer percentage display
    p <- ggplot2::ggplot(zero_df, ggplot2::aes(x = .data$pct_zeros)) +
        ggplot2::geom_histogram(bins = 30, fill = "steelblue", alpha = 0.7, color = "black") +
        ggplot2::scale_x_continuous(
            labels = function(x) paste0(x, "%"),
            breaks = scales::pretty_breaks(n = 6)
        ) +
        ggplot2::labs(
            x = paste0(x_label, " (%)"),
            y = "Frequency",
            title = paste0(title_text, " (n = ", n_items, " ", item_type, ")")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5)) +
        # Add count annotation
        ggplot2::annotate(
            "text",
            x = Inf, y = Inf,
            label = paste0("n = ", n_items, " ", item_type),
            hjust = 1.1, vjust = 1.5,
            size = 3.5,
            color = "gray40"
        ) +
        # Add summary statistics annotation
        ggplot2::annotate(
            "text",
            x = Inf, y = -Inf,
            label = paste0(
                "Range: ", min_zeros, "% - ", max_zeros, "%\n",
                "Mean: ", mean_zeros, "%",
                if (by == "sample" && perfect_coverage > 0) {
                    paste0("\nPerfect coverage (0%): ", perfect_coverage, " samples")
                } else {
                    ""
                }
            ),
            hjust = 1.1, vjust = -0.5,
            size = 3,
            color = "gray30"
        )

    return(p)
}

#' Plot Endometrial MA Plot
#'
#' Visualizes log2 fold change vs mean expression (MA plot) for differential expression results.
#'
#' @param de_table A data.frame from `esr_analyzeDifferentialExpression()` with columns: `gene_id`, `log2FC`, `AveExpr`, `FDR`.
#' @param fdr_threshold Numeric; FDR threshold for significance. Defaults to 0.05.
#' @param log2fc_threshold Numeric; log2 fold change threshold for large effect. Defaults to 1.
#' @param highlight_genes Optional character vector of gene IDs to highlight (e.g., with labels).
#' @param annot_col Optional character scalar; column name in de_table for gene annotations (not yet used).
#'
#' @return A ggplot object showing MA plot with significance coloring.
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
#' plotEndometrialMA(de_table)
#' @export
plotEndometrialMA <- function(de_table, fdr_threshold = 0.05,
                              log2fc_threshold = 1,
                              highlight_genes = NULL,
                              annot_col = NULL) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Validate required columns
    required_cols <- c("gene_id", "log2FC", "AveExpr", "FDR")
    missing_cols <- setdiff(required_cols, names(de_table))
    if (length(missing_cols) > 0) {
        stop(paste0("de_table must contain columns: ", paste(missing_cols, collapse = ", ")))
    }

    # Validate thresholds
    if (!is.numeric(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold > 1) {
        stop("fdr_threshold must be a numeric value between 0 and 1")
    }

    if (!is.numeric(log2fc_threshold) || log2fc_threshold < 0) {
        stop("log2fc_threshold must be a non-negative numeric value")
    }

    # Create significance categories
    de_table$significant <- with(de_table, {
        ifelse(FDR < fdr_threshold & abs(log2FC) >= log2fc_threshold, "Significant & Large FC",
            ifelse(FDR < fdr_threshold & abs(log2FC) < log2fc_threshold, "Significant Only",
                "Not Significant"
            )
        )
    })

    de_table$significant <- factor(de_table$significant,
        levels = c("Significant & Large FC", "Significant Only", "Not Significant")
    )

    # Count genes by category
    n_sig_large <- sum(de_table$significant == "Significant & Large FC")
    n_sig_only <- sum(de_table$significant == "Significant Only")
    n_total <- nrow(de_table)

    # Color palette (color-blind friendly)
    colors_map <- c(
        "Significant & Large FC" = "#d62728", # red
        "Significant Only" = "#ff7f0e", # orange
        "Not Significant" = "#7f7f7f" # gray
    )

    # Create plot
    p <- ggplot2::ggplot(de_table, ggplot2::aes(x = .data$AveExpr, y = .data$log2FC, color = .data$significant)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.5) +
        ggplot2::geom_hline(yintercept = log2fc_threshold, linetype = "dashed", color = "gray40", linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = -log2fc_threshold, linetype = "dashed", color = "gray40", linewidth = 0.5) +
        ggplot2::geom_hline(yintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
        ggplot2::scale_color_manual(
            name = "Significance",
            values = colors_map,
            labels = c(
                paste0("Significant & Large FC (n=", n_sig_large, ")"),
                paste0("Significant Only (n=", n_sig_only, ")"),
                paste0("Not Significant (n=", n_total - n_sig_large - n_sig_only, ")")
            )
        ) +
        ggplot2::labs(
            x = "Mean Expression (log2)",
            y = "Log2 Fold Change",
            title = paste0("MA Plot (n = ", n_total, " genes, FDR < ", fdr_threshold, ", |log2FC| >= ", log2fc_threshold, ")")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        # Add annotation
        ggplot2::annotate(
            "text",
            x = Inf, y = Inf,
            label = paste0("n = ", n_total, " genes"),
            hjust = 1.1, vjust = 1.5,
            size = 3.5,
            color = "gray40"
        )

    # Optionally highlight specific genes
    if (!is.null(highlight_genes)) {
        de_table_highlight <- de_table[de_table$gene_id %in% highlight_genes, ]
        if (nrow(de_table_highlight) > 0) {
            p <- p + ggplot2::geom_point(
                data = de_table_highlight,
                ggplot2::aes(x = .data$AveExpr, y = .data$log2FC),
                color = "black",
                shape = 21,
                fill = "yellow",
                size = 3,
                stroke = 1.5
            )
        }
    }

    return(p)
}

#' Plot Endometrial Volcano Plot
#'
#' Visualizes -log10(p-value) vs log2 fold change (Volcano plot) for differential expression results.
#'
#' @param de_table A data.frame from `esr_analyzeDifferentialExpression()` with columns: `gene_id`, `log2FC`, `FDR`.
#' @param fdr_threshold Numeric; FDR threshold for significance. Defaults to 0.05.
#' @param log2fc_threshold Numeric; log2 fold change threshold for large effect. Defaults to 1.
#' @param highlight_genes Optional character vector of gene IDs to highlight (e.g., with labels).
#' @param pvalue_col Character scalar; column name for p-values to plot. Defaults to "FDR".
#' @param log2fc_col Character scalar; column name for log2 fold change. Defaults to "log2FC".
#'
#' @return A ggplot object showing Volcano plot with significance coloring.
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
#' plotEndometrialVolcano(de_table)
#' @export
plotEndometrialVolcano <- function(de_table, fdr_threshold = 0.05,
                                   log2fc_threshold = 1,
                                   highlight_genes = NULL,
                                   pvalue_col = "FDR",
                                   log2fc_col = "log2FC") {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Validate required columns
    required_cols <- c("gene_id", log2fc_col, pvalue_col)
    missing_cols <- setdiff(required_cols, names(de_table))
    if (length(missing_cols) > 0) {
        stop(paste0("de_table must contain columns: ", paste(missing_cols, collapse = ", ")))
    }

    # Validate thresholds
    if (!is.numeric(fdr_threshold) || fdr_threshold <= 0 || fdr_threshold > 1) {
        stop("fdr_threshold must be a numeric value between 0 and 1")
    }

    if (!is.numeric(log2fc_threshold) || log2fc_threshold < 0) {
        stop("log2fc_threshold must be a non-negative numeric value")
    }

    # Extract columns (handle potential column name issues)
    de_table$log2FC_plot <- de_table[[log2fc_col]]
    de_table$pvalue_plot <- de_table[[pvalue_col]]

    # Calculate -log10(p-value)
    de_table$neg_log10_pval <- -log10(de_table$pvalue_plot + 1e-300) # Add small epsilon to avoid log(0)

    # Create significance categories
    de_table$significant <- with(de_table, {
        ifelse(pvalue_plot < fdr_threshold & abs(log2FC_plot) >= log2fc_threshold, "Significant & Large FC",
            ifelse(pvalue_plot < fdr_threshold & abs(log2FC_plot) < log2fc_threshold, "Significant Only",
                ifelse(pvalue_plot >= fdr_threshold & abs(log2FC_plot) >= log2fc_threshold, "Large FC Only",
                    "Not Significant"
                )
            )
        )
    })

    de_table$significant <- factor(de_table$significant,
        levels = c("Significant & Large FC", "Significant Only", "Large FC Only", "Not Significant")
    )

    # Count genes by category
    n_sig_large <- sum(de_table$significant == "Significant & Large FC")
    n_sig_only <- sum(de_table$significant == "Significant Only")
    n_large_fc <- sum(de_table$significant == "Large FC Only")
    n_total <- nrow(de_table)

    # Color palette (color-blind friendly)
    colors_map <- c(
        "Significant & Large FC" = "#d62728", # red
        "Significant Only" = "#ff7f0e", # orange
        "Large FC Only" = "#2ca02c", # green
        "Not Significant" = "#7f7f7f" # gray
    )

    # Create plot
    p <- ggplot2::ggplot(de_table, ggplot2::aes(x = .data$log2FC_plot, y = .data$neg_log10_pval, color = .data$significant)) +
        ggplot2::geom_point(alpha = 0.6, size = 1.5) +
        ggplot2::geom_vline(xintercept = log2fc_threshold, linetype = "dashed", color = "gray40", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = -log2fc_threshold, linetype = "dashed", color = "gray40", linewidth = 0.5) +
        ggplot2::geom_vline(xintercept = 0, linetype = "solid", color = "black", linewidth = 0.3) +
        ggplot2::geom_hline(yintercept = -log10(fdr_threshold), linetype = "dashed", color = "gray40", linewidth = 0.5) +
        ggplot2::scale_color_manual(
            name = "Significance",
            values = colors_map,
            labels = c(
                paste0("Significant & Large FC (n=", n_sig_large, ")"),
                paste0("Significant Only (n=", n_sig_only, ")"),
                paste0("Large FC Only (n=", n_large_fc, ")"),
                paste0("Not Significant (n=", n_total - n_sig_large - n_sig_only - n_large_fc, ")")
            )
        ) +
        ggplot2::labs(
            x = "Log2 Fold Change",
            y = paste0("-Log10(", pvalue_col, ")"),
            title = paste0("Volcano Plot (n = ", n_total, " genes, ", pvalue_col, " < ", fdr_threshold, ", |log2FC| >= ", log2fc_threshold, ")")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        # Add quadrants annotation
        ggplot2::annotate(
            "text",
            x = -Inf, y = Inf,
            label = "Down-regulated\nSignificant",
            hjust = -0.1, vjust = 1.5,
            size = 3,
            color = "gray50"
        ) +
        ggplot2::annotate(
            "text",
            x = Inf, y = Inf,
            label = "Up-regulated\nSignificant",
            hjust = 1.1, vjust = 1.5,
            size = 3,
            color = "gray50"
        ) +
        # Add count annotation
        ggplot2::annotate(
            "text",
            x = Inf, y = -Inf,
            label = paste0("n = ", n_total, " genes"),
            hjust = 1.1, vjust = -0.5,
            size = 3.5,
            color = "gray40"
        )

    # Optionally highlight specific genes
    if (!is.null(highlight_genes)) {
        de_table_highlight <- de_table[de_table$gene_id %in% highlight_genes, ]
        if (nrow(de_table_highlight) > 0) {
            p <- p + ggplot2::geom_point(
                data = de_table_highlight,
                ggplot2::aes(x = .data$log2FC_plot, y = .data$neg_log10_pval),
                color = "black",
                shape = 21,
                fill = "yellow",
                size = 3,
                stroke = 1.5
            )
        }
    }

    return(p)
}

#' Select Top Genes by Variance, DE, or Custom Criteria
#'
#' Selects top-K genes based on specified criteria (variance, differential expression, or custom sorting).
#'
#' @param mat_t A samples × genes transformed matrix (e.g., from `esr_transform_log1p_cpm()`). Required for variance-based selection.
#' @param de_table Optional data.frame from `esr_analyzeDifferentialExpression()` with columns: `gene_id`, `log2FC`, `FDR`. Required for DE-based selection.
#' @param n Integer; number of genes to select. Defaults to 50.
#' @param by Character scalar; selection method: "variance" (default), "de", or "custom".
#' @param sort_col Character scalar; column name in `de_table` for custom sorting (only used when `by = "custom"`). Defaults to NULL.
#'
#' @return Character vector of gene IDs (length `n`, or fewer if fewer genes available).
#'
#' @details
#' Selection methods:
#' - `by = "variance"`: Selects genes with highest variance across samples (useful for exploratory heatmaps). Requires `mat_t`.
#' - `by = "de"`: Selects top-K genes from DE table (sorted by FDR ascending, then by absolute log2FC descending). Requires `de_table`.
#' - `by = "custom"`: Uses `sort_col` to sort genes (ascending order). Requires `de_table` and `sort_col`.
#'
#' For variance-based selection, ties in variance are broken deterministically by gene ID (alphabetical order).
#' For DE-based selection, genes are first sorted by FDR (ascending), then by absolute log2FC (descending).
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#'
#' # Select top 20 genes by variance
#' top_var <- esr_selectTopGenes(mat_t, n = 20, by = "variance")
#' head(top_var)
#'
#' # Select top 10 genes by DE
#' de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_sample$pheno)
#' top_de <- esr_selectTopGenes(de_table = de_table, n = 10, by = "de")
#' head(top_de)
#' @importFrom stats var
#' @export
esr_selectTopGenes <- function(mat_t = NULL, de_table = NULL, n = 50,
                               by = c("variance", "de", "custom"), sort_col = NULL) {
    by <- match.arg(by)

    # Validate inputs based on selection method
    if (by == "variance") {
        if (is.null(mat_t)) {
            stop("mat_t is required for variance-based selection")
        }
        if (!is.matrix(mat_t) && !is.data.frame(mat_t)) {
            stop("mat_t must be a matrix or data.frame")
        }
        # Convert to matrix if needed
        if (!is.matrix(mat_t)) {
            mat_t <- as.matrix(mat_t)
        }
    } else if (by == "de" || by == "custom") {
        if (is.null(de_table)) {
            stop("de_table is required for DE-based or custom selection")
        }
        if (!is.data.frame(de_table)) {
            stop("de_table must be a data.frame")
        }
        if (!"gene_id" %in% names(de_table)) {
            stop("de_table must contain a 'gene_id' column")
        }
    }

    # Validate n
    if (!is.numeric(n) || n <= 0 || length(n) != 1) {
        stop("n must be a positive integer")
    }
    n <- as.integer(n)

    # Perform selection based on method
    if (by == "variance") {
        # Calculate variance per gene (across samples)
        # mat_t is samples × genes, so we transpose to get genes × samples
        gene_var <- apply(mat_t, 2, function(x) var(x, na.rm = TRUE))

        # Get gene IDs (column names of mat_t)
        gene_ids <- colnames(mat_t)

        if (is.null(gene_ids)) {
            stop("mat_t must have column names (gene IDs)")
        }

        # Sort by variance (descending), then by gene ID (ascending) for deterministic ties
        var_df <- data.frame(
            gene_id = gene_ids,
            variance = gene_var,
            stringsAsFactors = FALSE
        )
        var_df <- var_df[order(-var_df$variance, var_df$gene_id), ]

        # Select top n
        n_available <- nrow(var_df)
        n_select <- min(n, n_available)
        selected_genes <- var_df$gene_id[1:n_select]
    } else if (by == "de") {
        # Validate required columns
        required_cols <- c("gene_id", "FDR", "log2FC")
        missing_cols <- setdiff(required_cols, names(de_table))
        if (length(missing_cols) > 0) {
            stop(paste0("de_table must contain columns: ", paste(missing_cols, collapse = ", ")))
        }

        # Sort by FDR (ascending), then by absolute log2FC (descending)
        de_df <- data.frame(
            gene_id = de_table$gene_id,
            FDR = de_table$FDR,
            abs_log2FC = abs(de_table$log2FC),
            stringsAsFactors = FALSE
        )
        de_df <- de_df[order(de_df$FDR, -de_df$abs_log2FC), ]

        # Select top n
        n_available <- nrow(de_df)
        n_select <- min(n, n_available)
        selected_genes <- de_df$gene_id[1:n_select]
    } else if (by == "custom") {
        if (is.null(sort_col)) {
            stop("sort_col must be specified for custom selection")
        }
        if (!sort_col %in% names(de_table)) {
            stop(paste0("sort_col '", sort_col, "' not found in de_table"))
        }

        # Sort by sort_col (ascending)
        de_df <- data.frame(
            gene_id = de_table$gene_id,
            sort_value = de_table[[sort_col]],
            stringsAsFactors = FALSE
        )
        de_df <- de_df[order(de_df$sort_value), ]

        # Select top n
        n_available <- nrow(de_df)
        n_select <- min(n, n_available)
        selected_genes <- de_df$gene_id[1:n_select]
    }

    return(selected_genes)
}

#' Plot Endometrial Heatmap
#'
#' Creates a heatmap visualization of gene expression patterns across samples using ComplexHeatmap.
#'
#' @param mat_t A samples × genes transformed matrix (e.g., from `esr_transform_log1p_cpm()`).
#' @param genes Optional character vector of gene IDs to plot. If NULL, plots all genes. For large matrices, recommend selecting subset (e.g., using `esr_selectTopGenes()`).
#' @param pheno Optional data.frame with sample metadata for column annotations.
#' @param group_col Character scalar; name of the phenotype column for group coloring. Defaults to "group".
#' @param annot_cols Optional character vector of additional annotation column names from pheno to display.
#' @param cluster_rows Logical; whether to cluster genes (rows). Defaults to TRUE.
#' @param cluster_cols Logical; whether to cluster samples (columns). Defaults to TRUE.
#' @param show_row_names Logical; whether to show gene names (rows). Defaults to FALSE (typically set to FALSE for large heatmaps).
#' @param show_col_names Logical; whether to show sample names (columns). Defaults to TRUE.
#' @param scale Character scalar; scaling method: "row" (z-score per gene), "column" (z-score per sample), or "none" (no scaling). Defaults to "row".
#' @param color_palette Optional color palette function (default: color-blind friendly).
#' @param ... Additional arguments passed to `ComplexHeatmap::Heatmap()`.
#'
#' @return A ComplexHeatmap Heatmap object (use `ComplexHeatmap::draw()` to display).
#'
#' @details
#' This function creates publication-quality heatmaps using ComplexHeatmap. For large datasets (e.g., full 39,368 genes),
#' it is recommended to select a subset of genes (e.g., top-50 by variance or top-K by DE) to ensure readable output.
#'
#' Column annotations (e.g., group labels) can be added via the `pheno` parameter.
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#'
#' # Create heatmap with top 20 genes by variance
#' top_genes <- esr_selectTopGenes(mat_t, n = 20, by = "variance")
#' p <- plotEndometrialHeatmap(mat_t, genes = top_genes, pheno = gse201926_sample$pheno)
#' ComplexHeatmap::draw(p)
#'
#' # Create heatmap with phenotype annotations
#' p <- plotEndometrialHeatmap(mat_t,
#'     genes = top_genes, pheno = gse201926_sample$pheno,
#'     show_row_names = TRUE, scale = "row"
#' )
#' ComplexHeatmap::draw(p)
#' @importFrom utils head
#' @export
plotEndometrialHeatmap <- function(mat_t, genes = NULL, pheno = NULL,
                                   group_col = "group", annot_cols = NULL,
                                   cluster_rows = TRUE, cluster_cols = TRUE,
                                   show_row_names = FALSE, show_col_names = TRUE,
                                   scale = c("row", "column", "none"),
                                   color_palette = NULL, ...) {
    if (!requireNamespace("ComplexHeatmap", quietly = TRUE)) {
        stop("ComplexHeatmap package is required for heatmap plotting")
    }

    scale <- match.arg(scale)

    # Validate mat_t
    if (!is.matrix(mat_t) && !is.data.frame(mat_t)) {
        stop("mat_t must be a matrix or data.frame")
    }
    if (!is.matrix(mat_t)) {
        mat_t <- as.matrix(mat_t)
    }

    # Check that mat_t is samples × genes (rows = samples, cols = genes)
    if (is.null(rownames(mat_t))) {
        stop("mat_t must have rownames (sample IDs)")
    }
    if (is.null(colnames(mat_t))) {
        stop("mat_t must have colnames (gene IDs)")
    }

    # Subset to selected genes if specified
    if (!is.null(genes)) {
        if (!is.character(genes)) {
            stop("genes must be a character vector of gene IDs")
        }

        # Check which genes are available
        available_genes <- intersect(genes, colnames(mat_t))
        missing_genes <- setdiff(genes, colnames(mat_t))

        if (length(missing_genes) > 0) {
            warning(paste0(
                "Some requested genes not found in mat_t: ",
                paste(head(missing_genes, 5), collapse = ", "),
                if (length(missing_genes) > 5) " ..." else ""
            ))
        }

        if (length(available_genes) == 0) {
            stop("None of the requested genes are present in mat_t")
        }

        # Subset to available genes
        mat_t <- mat_t[, available_genes, drop = FALSE]
        genes <- available_genes # Update to actual available genes
    }

    # Prepare data for heatmap (genes × samples for ComplexHeatmap)
    # mat_t is samples × genes, need to transpose to genes × samples
    mat_heatmap <- t(mat_t)

    # Prepare column annotations if pheno provided
    col_annot <- NULL
    if (!is.null(pheno)) {
        if (!is.data.frame(pheno)) {
            stop("pheno must be a data.frame")
        }

        # Check sample ID matching
        sample_ids <- rownames(mat_t)
        if (!"sample_id" %in% names(pheno)) {
            # Try to match by rownames if sample_id not present
            if (nrow(pheno) == length(sample_ids) &&
                all(rownames(pheno) %in% sample_ids)) {
                pheno$sample_id <- rownames(pheno)
            } else {
                stop("pheno must contain a 'sample_id' column or rownames must match mat_t rownames")
            }
        }

        # Match samples
        pheno_matched <- pheno[pheno$sample_id %in% sample_ids, , drop = FALSE]
        if (nrow(pheno_matched) == 0) {
            warning("No matching samples found between mat_t and pheno; skipping annotations")
        } else {
            # Order pheno to match mat_t column order (samples)
            pheno_matched <- pheno_matched[match(sample_ids, pheno_matched$sample_id), , drop = FALSE]

            # Prepare annotation columns
            annot_list <- list()

            # Add group annotation if available
            if (group_col %in% names(pheno_matched)) {
                annot_list[[group_col]] <- pheno_matched[[group_col]]
            }

            # Add additional annotations if specified
            if (!is.null(annot_cols)) {
                for (col in annot_cols) {
                    if (col %in% names(pheno_matched)) {
                        annot_list[[col]] <- pheno_matched[[col]]
                    } else {
                        warning(paste0("Annotation column '", col, "' not found in pheno; skipping"))
                    }
                }
            }

            # Create column annotation if we have any annotations
            if (length(annot_list) > 0) {
                col_annot <- ComplexHeatmap::HeatmapAnnotation(
                    df = data.frame(annot_list, stringsAsFactors = FALSE),
                    col = list(group = c("PS" = "#2ca02c", "PIS" = "#d62728")) # Color-blind friendly colors
                )
            }
        }
    }

    # Set up color palette if not provided
    if (is.null(color_palette)) {
        # Default color-blind friendly palette (blue-white-red)
        color_palette <- circlize::colorRamp2(
            breaks = c(-2, 0, 2),
            colors = c("#0571b0", "#ffffff", "#ca0020")
        )
    }

    # Set up scaling
    if (scale == "row") {
        # Z-score per gene (row-wise)
        mat_heatmap_scaled <- t(scale(t(mat_heatmap)))
    } else if (scale == "column") {
        # Z-score per sample (column-wise)
        mat_heatmap_scaled <- scale(mat_heatmap)
    } else {
        # No scaling
        mat_heatmap_scaled <- mat_heatmap
    }

    # Create heatmap
    hm <- ComplexHeatmap::Heatmap(
        matrix = mat_heatmap_scaled,
        name = if (scale == "none") "Expression" else paste0("Z-score (", scale, ")"),
        cluster_rows = cluster_rows,
        cluster_columns = cluster_cols,
        show_row_names = show_row_names,
        show_column_names = show_col_names,
        col = color_palette,
        top_annotation = col_annot,
        column_title = paste0(
            "Heatmap (", ncol(mat_heatmap_scaled), " samples, ",
            nrow(mat_heatmap_scaled), " genes)"
        ),
        ...
    )

    return(hm)
}

#' Plot Endometrial ROC Curve
#'
#' Creates a Receiver Operating Characteristic (ROC) curve to visualize signature performance.
#'
#' @param predictions A data.frame with columns: `label` (0/1 binary labels), `prob` (raw probabilities), optionally `prob_calibrated` (calibrated probabilities).
#' @param use_calibrated Logical; use calibrated probabilities if available. Defaults to FALSE.
#' @param show_auc Logical; display AUC in plot annotation. Defaults to TRUE.
#' @param show_ci Logical; display confidence intervals (requires pROC package). Defaults to FALSE.
#' @param color_palette Optional color palette (default: color-blind friendly).
#' @param ... Additional arguments passed to `ggplot2::geom_line()`.
#'
#' @return A ggplot object showing ROC curve (FPR on x-axis, TPR on y-axis).
#'
#' @details
#' This function creates ROC curves for signature performance visualization. The ROC curve shows
#' the trade-off between true positive rate (TPR) and false positive rate (FPR) across different
#' probability thresholds.
#'
#' AUC (Area Under Curve) represents the probability that the classifier will rank a randomly
#' chosen positive instance higher than a randomly chosen negative instance. AUC = 1.0 indicates
#' perfect classification, AUC = 0.5 indicates random classification.
#'
#' @examples
#' \dontrun{
#' # After training a signature
#' result <- esr_trainEndometrialSignature(X, pheno)
#'
#' # Plot ROC curve with raw probabilities
#' p_roc_raw <- plotEndometrialROC(result$metrics$predictions, use_calibrated = FALSE)
#'
#' # Plot ROC curve with calibrated probabilities
#' p_roc_cal <- plotEndometrialROC(result$metrics$predictions, use_calibrated = TRUE)
#' }
#'
#' @export
plotEndometrialROC <- function(predictions, use_calibrated = FALSE, show_auc = TRUE,
                               show_ci = FALSE, color_palette = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Validate predictions
    if (!is.data.frame(predictions)) {
        stop("predictions must be a data.frame")
    }

    required_cols <- c("label", "prob")
    missing_cols <- setdiff(required_cols, names(predictions))
    if (length(missing_cols) > 0) {
        stop(paste0("predictions must contain columns: ", paste(missing_cols, collapse = ", ")))
    }

    # Select probability column
    if (use_calibrated && "prob_calibrated" %in% names(predictions)) {
        prob_col <- "prob_calibrated"
        prob_label <- "Calibrated Probability"
    } else {
        prob_col <- "prob"
        prob_label <- "Raw Probability"
        if (use_calibrated) {
            warning("Calibrated probabilities not available; using raw probabilities")
        }
    }

    probs <- predictions[[prob_col]]
    labels <- predictions$label

    # Validate data
    if (length(probs) == 0 || length(labels) == 0) {
        stop("predictions must contain at least one observation")
    }

    if (length(unique(labels)) != 2) {
        stop("predictions must contain exactly 2 unique label values (0 and 1)")
    }

    # Convert labels to numeric if needed
    if (!is.numeric(labels)) {
        labels <- as.numeric(as.factor(labels)) - 1
    }

    # Compute ROC curve
    # Sort by probabilities (descending)
    ord <- order(probs, decreasing = TRUE)
    probs_sorted <- probs[ord]
    labels_sorted <- labels[ord]

    # Count positives and negatives
    n_pos <- sum(labels_sorted == 1)
    n_neg <- sum(labels_sorted == 0)

    if (n_pos == 0 || n_neg == 0) {
        stop("predictions must contain at least one positive and one negative label")
    }

    # Compute TPR and FPR at each threshold
    cum_tp <- cumsum(labels_sorted == 1)
    cum_fp <- cumsum(labels_sorted == 0)

    tpr <- cum_tp / n_pos
    fpr <- cum_fp / n_neg

    # Add (0,0) and (1,1) for complete curve
    tpr <- c(0, tpr, 1)
    fpr <- c(0, fpr, 1)

    # Compute AUC using trapezoidal rule
    auc_value <- sum(diff(fpr) * (tpr[-1] + tpr[-length(tpr)]) / 2)

    # Create data frame for plotting
    roc_df <- data.frame(
        FPR = fpr,
        TPR = tpr,
        stringsAsFactors = FALSE
    )

    # Set up color palette
    if (is.null(color_palette)) {
        color_line <- "#2ca02c" # Color-blind friendly green
    } else {
        color_line <- color_palette
    }

    # Create plot
    p <- ggplot2::ggplot(roc_df, ggplot2::aes(x = .data$FPR, y = .data$TPR)) +
        ggplot2::geom_line(color = color_line, linewidth = 1.2, ...) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
        ggplot2::labs(
            x = "False Positive Rate",
            y = "True Positive Rate",
            title = paste0("ROC Curve (", prob_label, ")")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))

    # Add AUC annotation
    if (show_auc) {
        p <- p + ggplot2::annotate(
            "text",
            x = 0.6, y = 0.2,
            label = paste0("AUC = ", round(auc_value, 3)),
            size = 4.5,
            color = "black",
            fontface = "bold"
        )
    }

    # Add confidence intervals if requested (requires pROC)
    if (show_ci) {
        if (requireNamespace("pROC", quietly = TRUE)) {
            tryCatch(
                {
                    roc_obj <- pROC::roc(labels, probs, quiet = TRUE)
                    ci_obj <- pROC::ci.se(roc_obj, specificities = seq(0, 1, 0.01))
                    # Add CI to plot (simplified - could be more sophisticated)
                    # For now, just add note that CI is computed
                    auc_ci <- pROC::ci.auc(roc_obj, quiet = TRUE)
                    if (show_auc) {
                        p <- p + ggplot2::annotate(
                            "text",
                            x = 0.6, y = 0.15,
                            label = paste0("95% CI: [", round(auc_ci[1], 3), ", ", round(auc_ci[3], 3), "]"),
                            size = 3.5,
                            color = "gray40"
                        )
                    }
                },
                error = function(e) {
                    warning("Could not compute confidence intervals: ", conditionMessage(e))
                }
            )
        } else {
            warning("pROC package is required for confidence intervals. Install with: install.packages('pROC')")
        }
    }

    return(p)
}

#' Plot Endometrial Precision-Recall Curve
#'
#' Creates a Precision-Recall (PR) curve to visualize signature performance, especially useful for imbalanced datasets.
#'
#' @param predictions A data.frame with columns: `label` (0/1 binary labels), `prob` (raw probabilities), optionally `prob_calibrated` (calibrated probabilities).
#' @param use_calibrated Logical; use calibrated probabilities if available. Defaults to FALSE.
#' @param show_auc Logical; display PR-AUC in plot annotation. Defaults to TRUE.
#' @param color_palette Optional color palette (default: color-blind friendly).
#' @param ... Additional arguments passed to `ggplot2::geom_line()`.
#'
#' @return A ggplot object showing PR curve (Recall on x-axis, Precision on y-axis).
#'
#' @details
#' This function creates Precision-Recall curves for signature performance visualization. PR curves
#' are especially useful for imbalanced datasets where the positive class is rare.
#'
#' PR-AUC (Area Under Precision-Recall Curve) represents the average precision across all recall
#' values. Higher PR-AUC indicates better performance. The baseline (random classifier) is equal to
#' the prevalence of the positive class.
#'
#' @examples
#' \dontrun{
#' # After training a signature
#' result <- esr_trainEndometrialSignature(X, pheno)
#'
#' # Plot PR curve with raw probabilities
#' p_pr_raw <- plotEndometrialPR(result$metrics$predictions, use_calibrated = FALSE)
#'
#' # Plot PR curve with calibrated probabilities
#' p_pr_cal <- plotEndometrialPR(result$metrics$predictions, use_calibrated = TRUE)
#' }
#'
#' @export
plotEndometrialPR <- function(predictions, use_calibrated = FALSE, show_auc = TRUE,
                              color_palette = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Validate predictions
    if (!is.data.frame(predictions)) {
        stop("predictions must be a data.frame")
    }

    required_cols <- c("label", "prob")
    missing_cols <- setdiff(required_cols, names(predictions))
    if (length(missing_cols) > 0) {
        stop(paste0("predictions must contain columns: ", paste(missing_cols, collapse = ", ")))
    }

    # Select probability column
    if (use_calibrated && "prob_calibrated" %in% names(predictions)) {
        prob_col <- "prob_calibrated"
        prob_label <- "Calibrated Probability"
    } else {
        prob_col <- "prob"
        prob_label <- "Raw Probability"
        if (use_calibrated) {
            warning("Calibrated probabilities not available; using raw probabilities")
        }
    }

    probs <- predictions[[prob_col]]
    labels <- predictions$label

    # Validate data
    if (length(probs) == 0 || length(labels) == 0) {
        stop("predictions must contain at least one observation")
    }

    if (length(unique(labels)) != 2) {
        stop("predictions must contain exactly 2 unique label values (0 and 1)")
    }

    # Convert labels to numeric if needed
    if (!is.numeric(labels)) {
        labels <- as.numeric(as.factor(labels)) - 1
    }

    # Compute PR curve
    # Sort by probabilities (descending)
    ord <- order(probs, decreasing = TRUE)
    probs_sorted <- probs[ord]
    labels_sorted <- labels[ord]

    # Count positives and negatives
    n_pos <- sum(labels_sorted == 1)
    n_total <- length(labels_sorted)

    if (n_pos == 0) {
        stop("predictions must contain at least one positive label")
    }

    # Compute precision and recall at each threshold
    cum_tp <- cumsum(labels_sorted == 1)
    cum_all <- 1:n_total

    precision <- cum_tp / cum_all
    recall <- cum_tp / n_pos

    # Add (0,1) and (1,0) for complete curve (handle edge cases)
    precision <- c(1, precision)
    recall <- c(0, recall)

    # Handle NaN precision (when denominator is 0)
    precision[is.nan(precision)] <- 0

    # Compute PR-AUC using trapezoidal rule
    pr_auc_value <- sum(diff(recall) * (precision[-1] + precision[-length(precision)]) / 2)

    # Compute baseline (prevalence of positive class)
    prevalence <- n_pos / n_total

    # Create data frame for plotting
    pr_df <- data.frame(
        Recall = recall,
        Precision = precision,
        stringsAsFactors = FALSE
    )

    # Set up color palette
    if (is.null(color_palette)) {
        color_line <- "#d62728" # Color-blind friendly red
    } else {
        color_line <- color_palette
    }

    # Create plot
    p <- ggplot2::ggplot(pr_df, ggplot2::aes(x = .data$Recall, y = .data$Precision)) +
        ggplot2::geom_line(color = color_line, linewidth = 1.2, ...) +
        ggplot2::geom_hline(
            yintercept = prevalence, linetype = "dashed", color = "gray40", linewidth = 0.8
        ) +
        ggplot2::labs(
            x = "Recall",
            y = "Precision",
            title = paste0("Precision-Recall Curve (", prob_label, ")")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))

    # Add PR-AUC annotation
    if (show_auc) {
        p <- p + ggplot2::annotate(
            "text",
            x = 0.6, y = 0.2,
            label = paste0("PR-AUC = ", round(pr_auc_value, 3)),
            size = 4.5,
            color = "black",
            fontface = "bold"
        ) +
            ggplot2::annotate(
                "text",
                x = 0.6, y = 0.15,
                label = paste0("Baseline = ", round(prevalence, 3)),
                size = 3.5,
                color = "gray40"
            )
    }

    return(p)
}

#' Plot Endometrial Calibration Curve
#'
#' Creates a calibration curve to visualize how well predicted probabilities match observed frequencies.
#'
#' @param predictions A data.frame with columns: `label` (0/1 binary labels), `prob` (raw probabilities), optionally `prob_calibrated` (calibrated probabilities).
#' @param use_calibrated Logical; use calibrated probabilities if available. Defaults to FALSE.
#' @param n_bins Integer; number of bins for calibration plot. Defaults to 10.
#' @param show_brier Logical; display Brier Score in annotation. Defaults to TRUE.
#' @param show_ece Logical; display Expected Calibration Error (ECE) in annotation. Defaults to TRUE.
#' @param color_palette Optional color palette (default: color-blind friendly).
#' @param ... Additional arguments passed to `ggplot2::geom_point()` or `ggplot2::geom_bar()`.
#'
#' @return A ggplot object showing calibration curve (predicted probability on x-axis, observed frequency on y-axis).
#'
#' @details
#' This function creates calibration curves to assess probability calibration quality. A well-calibrated
#' model should have predicted probabilities that match observed frequencies (points close to the
#' diagonal y = x).
#'
#' Brier Score measures the mean squared error between predicted probabilities and binary outcomes.
#' Lower Brier Score indicates better calibration (perfect calibration has Brier Score = 0).
#'
#' ECE (Expected Calibration Error) measures the weighted average of absolute difference between
#' predicted probabilities and observed frequencies within bins. Lower ECE indicates better calibration.
#'
#' @examples
#' \dontrun{
#' # After training a signature
#' result <- esr_trainEndometrialSignature(X, pheno, calibration_method = "platt")
#'
#' # Plot calibration curve with raw probabilities
#' p_cal_raw <- plotEndometrialCalibration(result$metrics$predictions, use_calibrated = FALSE)
#'
#' # Plot calibration curve with calibrated probabilities
#' p_cal_cal <- plotEndometrialCalibration(result$metrics$predictions, use_calibrated = TRUE)
#' }
#'
#' @export
plotEndometrialCalibration <- function(predictions, use_calibrated = FALSE, n_bins = 10,
                                       show_brier = TRUE, show_ece = TRUE,
                                       color_palette = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Validate predictions
    if (!is.data.frame(predictions)) {
        stop("predictions must be a data.frame")
    }

    required_cols <- c("label", "prob")
    missing_cols <- setdiff(required_cols, names(predictions))
    if (length(missing_cols) > 0) {
        stop(paste0("predictions must contain columns: ", paste(missing_cols, collapse = ", ")))
    }

    # Validate n_bins
    if (!is.numeric(n_bins) || n_bins <= 0 || length(n_bins) != 1) {
        stop("n_bins must be a positive integer")
    }
    n_bins <- as.integer(n_bins)

    # Select probability column
    if (use_calibrated && "prob_calibrated" %in% names(predictions)) {
        prob_col <- "prob_calibrated"
        prob_label <- "Calibrated Probability"
    } else {
        prob_col <- "prob"
        prob_label <- "Raw Probability"
        if (use_calibrated) {
            warning("Calibrated probabilities not available; using raw probabilities")
        }
    }

    probs <- predictions[[prob_col]]
    labels <- predictions$label

    # Validate data
    if (length(probs) == 0 || length(labels) == 0) {
        stop("predictions must contain at least one observation")
    }

    if (length(unique(labels)) != 2) {
        stop("predictions must contain exactly 2 unique label values (0 and 1)")
    }

    # Convert labels to numeric if needed
    if (!is.numeric(labels)) {
        labels <- as.numeric(as.factor(labels)) - 1
    }

    # Bin probabilities
    bin_edges <- seq(0, 1, length.out = n_bins + 1)
    bin_edges[1] <- 0
    bin_edges[length(bin_edges)] <- 1.0001 # Ensure max prob is included
    bin_ids <- cut(probs, breaks = bin_edges, include.lowest = TRUE, right = FALSE)

    # Compute calibration statistics per bin
    bin_stats <- data.frame(
        bin_id = levels(bin_ids),
        bin_mid = (bin_edges[-length(bin_edges)] + bin_edges[-1]) / 2,
        n_obs = as.numeric(table(bin_ids)),
        mean_pred = tapply(probs, bin_ids, mean, na.rm = TRUE),
        mean_obs = tapply(labels, bin_ids, mean, na.rm = TRUE),
        stringsAsFactors = FALSE
    )

    # Remove bins with no observations
    bin_stats <- bin_stats[bin_stats$n_obs > 0, ]

    if (nrow(bin_stats) == 0) {
        stop("No valid bins found. Check that probabilities are in [0, 1] range.")
    }

    # Compute Brier Score
    brier_score <- mean((probs - labels)^2, na.rm = TRUE)

    # Compute ECE (Expected Calibration Error)
    ece_value <- sum(bin_stats$n_obs / length(probs) * abs(bin_stats$mean_pred - bin_stats$mean_obs), na.rm = TRUE)

    # Set up color palette
    if (is.null(color_palette)) {
        color_point <- "#1f77b4" # Color-blind friendly blue
        color_line <- "#ff7f0e" # Color-blind friendly orange
    } else {
        color_point <- color_palette
        color_line <- color_palette
    }

    # Create plot
    p <- ggplot2::ggplot(bin_stats, ggplot2::aes(x = .data$mean_pred, y = .data$mean_obs)) +
        ggplot2::geom_point(size = 3, color = color_point, ...) +
        ggplot2::geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
        ggplot2::labs(
            x = "Mean Predicted Probability",
            y = "Mean Observed Frequency",
            title = paste0("Calibration Curve (", prob_label, ")")
        ) +
        ggplot2::theme_bw() +
        ggplot2::theme(
            plot.title = ggplot2::element_text(hjust = 0.5),
            legend.position = "right"
        ) +
        ggplot2::coord_fixed(ratio = 1, xlim = c(0, 1), ylim = c(0, 1))

    # Add Brier Score and ECE annotations
    if (show_brier || show_ece) {
        annotation_text <- character(0)
        if (show_brier) {
            annotation_text <- c(annotation_text, paste0("Brier Score = ", round(brier_score, 3)))
        }
        if (show_ece) {
            annotation_text <- c(annotation_text, paste0("ECE = ", round(ece_value, 3)))
        }

        p <- p + ggplot2::annotate(
            "text",
            x = 0.05, y = 0.95,
            label = paste(annotation_text, collapse = "\n"),
            size = 4,
            color = "black",
            fontface = "bold",
            hjust = 0,
            vjust = 1
        )
    }

    return(p)
}

#' Plot Endometrial Signature Comparison
#'
#' Creates comparison plots to compare pre-trained vs new signature performance (ROC, PR, and Calibration curves).
#'
#' @param pretrained_result Optional list from pre-trained signature application (Phase 3; default: NULL). If NULL, plots only new signature.
#' @param new_result List from `esr_trainEndometrialSignature()` with `metrics` and `predictions` components.
#' @param metrics_to_plot Character vector; which plots to include: "roc", "pr", "calibration", or "all". Defaults to c("roc", "pr", "calibration").
#' @param show_metrics_table Logical; display metrics comparison table. Defaults to TRUE.
#' @param color_palette Optional color palette (default: color-blind friendly).
#' @param ... Additional arguments passed to plot functions.
#'
#' @return A list of ggplot objects (or a single plot if only one metric requested). If `show_metrics_table=TRUE`, includes a metrics table.
#'
#' @details
#' This function creates side-by-side or overlay comparison plots for signature performance.
#' It compares ROC curves, PR curves, and calibration curves between pre-trained and new signatures.
#'
#' If `pretrained_result` is NULL, only the new signature is plotted (single signature mode).
#' If `pretrained_result` is provided, both signatures are plotted for comparison.
#'
#' The metrics table (if `show_metrics_table=TRUE`) compares AUC, accuracy, Brier Score, and ECE
#' between signatures.
#'
#' @examples
#' \dontrun{
#' # After training a signature
#' new_result <- esr_trainEndometrialSignature(X, pheno)
#'
#' # Compare with pre-trained signature (if available)
#' pretrained_result <- esr_loadPretrainedSignature()
#'
#' # Generate comparison plots
#' comparison_plots <- plotEndometrialComparison(
#'     pretrained_result = pretrained_result,
#'     new_result = new_result,
#'     metrics_to_plot = "all"
#' )
#'
#' # Plot only new signature (no comparison)
#' new_plots <- plotEndometrialComparison(
#'     pretrained_result = NULL,
#'     new_result = new_result
#' )
#' }
#'
#' @export
plotEndometrialComparison <- function(pretrained_result = NULL, new_result,
                                      metrics_to_plot = c("roc", "pr", "calibration"),
                                      show_metrics_table = TRUE,
                                      color_palette = NULL, ...) {
    if (!requireNamespace("ggplot2", quietly = TRUE)) {
        stop("ggplot2 package is required for plotting")
    }

    # Validate new_result
    if (!is.list(new_result)) {
        stop("new_result must be a list from esr_trainEndometrialSignature()")
    }

    if (!"metrics" %in% names(new_result)) {
        stop("new_result must contain 'metrics' component")
    }

    if (!"predictions" %in% names(new_result$metrics)) {
        stop("new_result$metrics must contain 'predictions' data.frame")
    }

    new_predictions <- new_result$metrics$predictions

    # Validate pretrained_result if provided
    pretrained_predictions <- NULL
    if (!is.null(pretrained_result)) {
        if (!is.list(pretrained_result)) {
            stop("pretrained_result must be a list or NULL")
        }
        # Check if it has predictions (structure may vary by Phase 3 implementation)
        if ("predictions" %in% names(pretrained_result)) {
            pretrained_predictions <- pretrained_result$predictions
        } else if ("metrics" %in% names(pretrained_result) && "predictions" %in% names(pretrained_result$metrics)) {
            pretrained_predictions <- pretrained_result$metrics$predictions
        } else {
            warning("pretrained_result does not contain predictions; plotting only new signature")
            pretrained_result <- NULL
        }
    }

    # Match metrics_to_plot
    valid_metrics <- c("roc", "pr", "calibration", "all")
    if (any(!metrics_to_plot %in% valid_metrics)) {
        stop(paste0("metrics_to_plot must be one of: ", paste(valid_metrics, collapse = ", ")))
    }

    if ("all" %in% metrics_to_plot) {
        metrics_to_plot <- c("roc", "pr", "calibration")
    }

    # Set up color palette
    if (is.null(color_palette)) {
        color_new <- "#2ca02c" # Color-blind friendly green
        color_pretrained <- "#d62728" # Color-blind friendly red
    } else {
        color_new <- color_palette
        color_pretrained <- color_palette
    }

    # Create plots
    plots_list <- list()

    # ROC curve
    if ("roc" %in% metrics_to_plot) {
        p_roc_new <- plotEndometrialROC(new_predictions,
            use_calibrated = FALSE,
            show_auc = TRUE, color_palette = color_new, ...
        )
        if (!is.null(pretrained_predictions)) {
            p_roc_pretrained <- plotEndometrialROC(pretrained_predictions,
                use_calibrated = FALSE,
                show_auc = TRUE, color_palette = color_pretrained, ...
            )
            # Overlay both curves (simplified - could use patchwork for better layout)
            p_roc_new <- p_roc_new + ggplot2::geom_line(
                data = p_roc_pretrained$data,
                ggplot2::aes(x = .data$FPR, y = .data$TPR),
                color = color_pretrained, linewidth = 1.2, ...
            )
        }
        plots_list[["roc"]] <- p_roc_new
    }

    # PR curve
    if ("pr" %in% metrics_to_plot) {
        p_pr_new <- plotEndometrialPR(new_predictions,
            use_calibrated = FALSE,
            show_auc = TRUE, color_palette = color_new, ...
        )
        if (!is.null(pretrained_predictions)) {
            p_pr_pretrained <- plotEndometrialPR(pretrained_predictions,
                use_calibrated = FALSE,
                show_auc = TRUE, color_palette = color_pretrained, ...
            )
            # Overlay both curves
            p_pr_new <- p_pr_new + ggplot2::geom_line(
                data = p_pr_pretrained$data,
                ggplot2::aes(x = .data$Recall, y = .data$Precision),
                color = color_pretrained, linewidth = 1.2, ...
            )
        }
        plots_list[["pr"]] <- p_pr_new
    }

    # Calibration curve
    if ("calibration" %in% metrics_to_plot) {
        p_cal_new <- plotEndometrialCalibration(new_predictions,
            use_calibrated = FALSE,
            show_brier = TRUE, show_ece = TRUE,
            color_palette = color_new, ...
        )
        if (!is.null(pretrained_predictions)) {
            p_cal_pretrained <- plotEndometrialCalibration(pretrained_predictions,
                use_calibrated = FALSE,
                show_brier = TRUE, show_ece = TRUE,
                color_palette = color_pretrained, ...
            )
            # Overlay both curves
            p_cal_new <- p_cal_new + ggplot2::geom_point(
                data = p_cal_pretrained$data,
                ggplot2::aes(x = .data$mean_pred, y = .data$mean_obs),
                color = color_pretrained, size = 3, shape = 17, ...
            )
        }
        plots_list[["calibration"]] <- p_cal_new
    }

    # Create metrics table if requested
    metrics_table <- NULL
    if (show_metrics_table) {
        new_metrics <- new_result$metrics
        metrics_table <- data.frame(
            Metric = c("AUC", "Accuracy", "Brier Score", "ECE"),
            New_Signature = c(
                ifelse("auc" %in% names(new_metrics), round(new_metrics$auc, 3), NA),
                ifelse("accuracy" %in% names(new_metrics), round(new_metrics$accuracy, 3), NA),
                ifelse("brier_score" %in% names(new_metrics), round(new_metrics$brier_score, 3), NA),
                ifelse("ece" %in% names(new_metrics), round(new_metrics$ece, 3), NA)
            ),
            stringsAsFactors = FALSE
        )

        if (!is.null(pretrained_result)) {
            pretrained_metrics <- NULL
            if ("metrics" %in% names(pretrained_result)) {
                pretrained_metrics <- pretrained_result$metrics
            } else if ("auc" %in% names(pretrained_result)) {
                pretrained_metrics <- pretrained_result
            }

            if (!is.null(pretrained_metrics)) {
                metrics_table$Pretrained_Signature <- c(
                    ifelse("auc" %in% names(pretrained_metrics), round(pretrained_metrics$auc, 3), NA),
                    ifelse("accuracy" %in% names(pretrained_metrics), round(pretrained_metrics$accuracy, 3), NA),
                    ifelse("brier_score" %in% names(pretrained_metrics), round(pretrained_metrics$brier_score, 3), NA),
                    ifelse("ece" %in% names(pretrained_metrics), round(pretrained_metrics$ece, 3), NA)
                )
            }
        }

        plots_list[["metrics_table"]] <- metrics_table
    }

    # Return single plot if only one requested, otherwise list
    if (length(plots_list) == 1 && !show_metrics_table) {
        return(plots_list[[1]])
    } else {
        return(plots_list)
    }
}

# [END]
