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
    p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, color = group, shape = group)) +
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
        p <- ggplot2::ggplot(lib_df, ggplot2::aes(x = group, y = library_size, fill = group)) +
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
        p <- ggplot2::ggplot(lib_df, ggplot2::aes(x = library_size)) +
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
    p <- ggplot2::ggplot(zero_df, ggplot2::aes(x = pct_zeros)) +
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

# [END]
