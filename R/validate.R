#' Validate Endometrial Data Structures
#'
#' Performs basic schema checks for endometrial bulk RNA-seq inputs including
#' ID matching, class imbalance detection, and annotation validation.
#'
#' @param X A matrix/data.frame of gene expression (genes x samples).
#' @param pheno A data.frame with sample metadata; must include `sample_id` and `group` columns by default.
#' @param annot Optional data.frame with gene annotations.
#' @param label_col Character scalar; name of the phenotype column containing PS/PIS labels. Defaults to "group".
#'
#' @return A list with elements:
#' \describe{
#'   \item{X}{The input expression matrix (cleaned/converted to matrix)}
#'   \item{pheno}{The input phenotype data.frame}
#'   \item{annot}{The input annotation data.frame (or NULL)}
#'   \item{issues}{A data.frame with columns `type` and `message` containing detected validation issues}
#' }
#'
#' @examples
#' data(gse201926_sample)
#' result <- esr_validateEndometrial(
#'   gse201926_sample$counts,
#'   gse201926_sample$pheno,
#'   annot = gse201926_sample$annot
#' )
#' result$issues
#' @export
esr_validateEndometrial <- function(X, pheno, annot = NULL, label_col = "group") {
    issues <- data.frame(type = character(), message = character(), stringsAsFactors = FALSE)
    
    # Convert X to matrix if needed
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    
    # Validate X is numeric
    if (!is.numeric(X)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "Expression matrix X must contain numeric values",
            stringsAsFactors = FALSE
        ))
        return(list(X = X, pheno = pheno, annot = annot, issues = issues))
    }
    
    # Validate X dimensions
    if (nrow(X) == 0 || ncol(X) == 0) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "Expression matrix X must have at least one gene and one sample",
            stringsAsFactors = FALSE
        ))
        return(list(X = X, pheno = pheno, annot = annot, issues = issues))
    }
    
    # Validate pheno structure
    if (!is.data.frame(pheno)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "pheno must be a data.frame",
            stringsAsFactors = FALSE
        ))
        return(list(X = X, pheno = pheno, annot = annot, issues = issues))
    }
    
    # Check required columns in pheno
    if (!"sample_id" %in% names(pheno)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = "pheno must contain a 'sample_id' column",
            stringsAsFactors = FALSE
        ))
    }
    
    if (!label_col %in% names(pheno)) {
        issues <- rbind(issues, data.frame(
            type = "error",
            message = paste0("pheno must contain a '", label_col, "' column"),
            stringsAsFactors = FALSE
        ))
    }
    
    # Check ID matching between X and pheno
    if ("sample_id" %in% names(pheno) && !is.null(colnames(X))) {
        X_samples <- colnames(X)
        pheno_samples <- pheno$sample_id
        missing_in_pheno <- setdiff(X_samples, pheno_samples)
        missing_in_X <- setdiff(pheno_samples, X_samples)
        
        if (length(missing_in_pheno) > 0) {
            issues <- rbind(issues, data.frame(
                type = "warning",
                message = paste0("Samples in X not found in pheno: ", paste(missing_in_pheno, collapse = ", ")),
                stringsAsFactors = FALSE
            ))
        }
        
        if (length(missing_in_X) > 0) {
            issues <- rbind(issues, data.frame(
                type = "warning",
                message = paste0("Samples in pheno not found in X: ", paste(missing_in_X, collapse = ", ")),
                stringsAsFactors = FALSE
            ))
        }
        
        # Calculate mapping rate
        if (length(X_samples) > 0 && length(pheno_samples) > 0) {
            mapped <- length(intersect(X_samples, pheno_samples))
            mapping_rate <- mapped / max(length(X_samples), length(pheno_samples))
            if (mapping_rate < 0.8) {
                issues <- rbind(issues, data.frame(
                    type = "warning",
                    message = paste0("Low ID mapping rate: ", round(mapping_rate * 100, 1), "% of samples matched"),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
    
    # Check class balance
    if (label_col %in% names(pheno)) {
        labels <- pheno[[label_col]]
        label_table <- table(labels)
        
        if (length(label_table) >= 2) {
            proportions <- prop.table(label_table)
            min_prop <- min(proportions)
            max_prop <- max(proportions)
            
            if (min_prop < 0.2) {
                issues <- rbind(issues, data.frame(
                    type = "warning",
                    message = paste0("Class imbalance detected: smallest class is ", 
                                   round(min_prop * 100, 1), "% of samples"),
                    stringsAsFactors = FALSE
                ))
            }
            
            if (max_prop > 0.8) {
                issues <- rbind(issues, data.frame(
                    type = "warning",
                    message = paste0("Class imbalance detected: largest class is ", 
                                   round(max_prop * 100, 1), "% of samples"),
                    stringsAsFactors = FALSE
                ))
            }
        }
    }
    
    # Validate annot structure if provided
    if (!is.null(annot)) {
        if (!is.data.frame(annot)) {
            issues <- rbind(issues, data.frame(
                type = "warning",
                message = "annot should be a data.frame; ignoring",
                stringsAsFactors = FALSE
            ))
        } else if (!is.null(rownames(X))) {
            # Check if gene IDs match
            annot_gene_ids <- annot$GeneID
            if (!is.null(annot_gene_ids)) {
                X_genes <- rownames(X)
                matched_genes <- length(intersect(X_genes, annot_gene_ids))
                if (matched_genes == 0) {
                    issues <- rbind(issues, data.frame(
                        type = "warning",
                        message = "No matching gene IDs found between X rownames and annot$GeneID",
                        stringsAsFactors = FALSE
                    ))
                }
            }
        }
    }
    
    result <- list(X = X, pheno = pheno, annot = annot, issues = issues)
    return(result)
}

#' Transform Counts to log1p-CPM
#'
#' Converts raw counts to counts per million (CPM) with log1p transformation,
#' optionally filtering low-expressed genes.
#'
#' @param X A matrix/data.frame of raw counts (genes x samples).
#' @param cpm_min Minimum CPM threshold for gene filtering. Defaults to 1.
#' @param cpm_min_samples Minimum number of samples that must meet CPM threshold. Defaults to 4.
#'
#' @return A transformed matrix (samples x genes) ready for PCA/visualization.
#'   Rows are samples, columns are genes. Genes failing CPM filter are removed.
#'
#' @examples
#' data(gse201926_sample)
#' mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
#' dim(mat_t)  # samples x genes
#' @export
esr_transform_log1p_cpm <- function(X, cpm_min = 1, cpm_min_samples = 4) {
    # Convert to matrix if needed
    if (!is.matrix(X)) {
        X <- as.matrix(X)
    }
    
    # Calculate library sizes (column sums)
    lib_sizes <- colSums(X, na.rm = TRUE)
    
    # Calculate CPM (counts per million)
    cpm <- t(t(X) / lib_sizes * 1e6)
    
    # Filter genes with CPM < cpm_min in fewer than cpm_min_samples
    genes_keep <- rowSums(cpm >= cpm_min, na.rm = TRUE) >= cpm_min_samples
    
    # Subset and apply log1p transformation
    cpm_filtered <- cpm[genes_keep, , drop = FALSE]
    mat_t <- log1p(cpm_filtered)
    
    # Transpose to samples x genes format
    mat_t <- t(mat_t)
    
    return(mat_t)
}

# [END]
