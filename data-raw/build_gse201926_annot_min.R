# Build minimal annotation dataset for vignettes/examples

suppressPackageStartupMessages({
    library(readr)
    library(usethis)
    library(dplyr)
})

pick_extdata <- function(name) {
    # Prefer local development path if running inside package root
    local_ext <- file.path("inst", "extdata")
    local_gz <- file.path(local_ext, paste0(name, ".gz"))
    local_plain <- file.path(local_ext, name)
    if (file.exists(local_gz)) {
        return(normalizePath(local_gz))
    }
    if (file.exists(local_plain)) {
        return(normalizePath(local_plain))
    }

    # Fallback to installed package path
    gz <- system.file("extdata", paste0(name, ".gz"), package = "endoSignatureR")
    plain <- system.file("extdata", name, package = "endoSignatureR")
    if (nzchar(gz) && file.exists(gz)) {
        return(gz)
    }
    if (nzchar(plain) && file.exists(plain)) {
        return(plain)
    }
    stop("File not found in inst/extdata: ", name)
}

annot_path <- pick_extdata("gse201926_annotation.tsv")

message("Reading annotation from: ", annot_path)
annot_df <- readr::read_tsv(annot_path, show_col_types = FALSE)

# Heuristic mapping to minimal columns; adapt if column names differ
candidate_cols <- c(
    gene_id = "GeneID|gene_id|ensembl_gene_id",
    symbol = "Symbol|symbol|hgnc_symbol",
    biotype = "GeneType|gene_biotype|biotype"
)

pick_first_match <- function(pattern, df) {
    cols <- grep(pattern, names(df), ignore.case = TRUE, value = TRUE)
    if (length(cols) == 0) {
        return(NA_character_)
    }
    cols[[1]]
}

gene_col <- pick_first_match(candidate_cols[["gene_id"]], annot_df)
sym_col <- pick_first_match(candidate_cols[["symbol"]], annot_df)
bio_col <- pick_first_match(candidate_cols[["biotype"]], annot_df)

stopifnot(!is.na(gene_col))

gse201926_annot_min <- annot_df %>%
    transmute(
        gene_id = .data[[gene_col]],
        symbol = if (!is.na(sym_col)) .data[[sym_col]] else NA_character_,
        biotype = if (!is.na(bio_col)) .data[[bio_col]] else NA_character_
    )

message("Saving data object: gse201926_annot_min (", nrow(gse201926_annot_min), " rows)")
usethis::use_data(gse201926_annot_min, overwrite = TRUE, compress = "xz")

# [END]
