# Build minimal phenotype dataset for vignettes/examples

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

series_path <- pick_extdata("gse201926_series_matrix.txt")

# Reuse the package's parser if available
if (exists("endo_parse_metadata")) {
    pheno <- endo_parse_metadata(series_path)
} else {
    lines <- readLines(series_path)
    data_start <- which(grepl("^!series_matrix_table_begin", lines))
    header_line <- lines[data_start + 1]
    sample_ids <- unlist(strsplit(header_line, "\t"))[-1]
    groups <- ifelse(grepl("GSM608117[3-8]", sample_ids), "PIS", "PS")
    title_line_idx <- which(grepl("^!Sample_title", lines))
    titles <- if (length(title_line_idx) > 0) unlist(strsplit(lines[title_line_idx[1]], "\t"))[-1] else rep(NA_character_, length(sample_ids))
    pheno <- data.frame(sample_id = sample_ids, group = groups, title = titles, stringsAsFactors = FALSE)
}

# Optional batch detection (placeholder: none provided in series matrix)
gse201926_pheno_min <- pheno %>%
    transmute(sample_id = .data$sample_id, group = .data$group, batch = NA_character_)

message("Saving data object: gse201926_pheno_min (", nrow(gse201926_pheno_min), " rows)")
usethis::use_data(gse201926_pheno_min, overwrite = TRUE, compress = "xz")
