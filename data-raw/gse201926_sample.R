## code to prepare `gse201926_sample` dataset

library(readr)
library(dplyr)

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

raw_counts_path <- pick_ext("gse201926_raw_counts.tsv")
raw_counts <- readr::read_tsv(raw_counts_path, show_col_types = FALSE)

# Select top 200 most variable genes
gene_vars <- apply(raw_counts[, -1], 1, var)
top_genes <- order(gene_vars, decreasing = TRUE)[1:200]

# Create subset counts matrix
counts_matrix <- as.matrix(raw_counts[top_genes, -1])
rownames(counts_matrix) <- raw_counts$GeneID[top_genes]
colnames(counts_matrix) <- colnames(raw_counts)[-1]

# Create phenotype data
pheno_data <- data.frame(
    sample_id = c(
        "GSM6081173", "GSM6081174", "GSM6081175", "GSM6081176", "GSM6081177", "GSM6081178",
        "GSM6081179", "GSM6081180", "GSM6081181", "GSM6081182", "GSM6081183", "GSM6081184"
    ),
    group = c(rep("PIS", 6), rep("PS", 6)),
    stringsAsFactors = FALSE
)

# Load gene annotation and subset from packaged extdata
annot_path <- pick_ext("gse201926_annotation.tsv")
annot <- readr::read_tsv(annot_path, show_col_types = FALSE)
annot_subset <- annot[annot$GeneID %in% raw_counts$GeneID[top_genes], ]

# Create final data object
gse201926_sample <- list(
    counts = counts_matrix,
    pheno = pheno_data,
    annot = annot_subset
)

# Save as package data
usethis::use_data(gse201926_sample, overwrite = TRUE)
