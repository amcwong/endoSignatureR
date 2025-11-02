#' GSE201926 external data files (inst/extdata)
#'
#' The following files were copied into `inst/extdata/` and lightly renamed
#' for consistency. Original filenames (left) → package filenames (right):
#'
#' - `GSE201926_raw_counts_GRCh38.p13_NCBI.tsv` → `gse201926_raw_counts.tsv`
#' - `GSE201926_norm_counts_TPM_GRCh38.p13_NCBI.tsv` → `gse201926_tpm_counts.tsv`
#' - `Human.GRCh38.p13.annot.tsv` → `gse201926_annotation.tsv`
#' - `GSE201926_series_matrix.txt` → `gse201926_series_matrix.txt` (unchanged)
#'
#' These files are accessed via `system.file("extdata", <name>, package = "endoSignatureR")`.
#'
#' Load GSE201926 sample dataset
#'
#' Loads the bundled subset of the GSE201926 endometrial lesion RNA-seq dataset
#' containing 200 most variable genes across 12 samples (6 PS + 6 PIS).
#'
#' @return List containing:
#' \describe{
#'   \item{counts}{Matrix of raw read counts (200 genes × 12 samples)}
#'   \item{pheno}{Data frame with sample metadata (sample_id, group)}
#'   \item{annot}{Data frame with gene annotation}
#' }
#'
#' @examples
#' demo_data <- endo_load_demo()
#' head(demo_data$counts[, 1:3])
#' demo_data$pheno
#'
#' @export
endo_load_demo <- function() {
    utils::data("gse201926_sample", package = "endoSignatureR", envir = environment())
    return(get("gse201926_sample", envir = environment()))
}

#' Load full GSE201926 dataset
#'
#' @param sample_only Logical; if TRUE, return the bundled demo subset
#'   (`gse201926_sample`) instead of reading full matrices from `inst/extdata`.
#'
#' Loads the complete GSE201926 endometrial lesion RNA-seq dataset
#' containing all 39,368 genes across 12 samples (6 PS + 6 PIS).
#'
#' @return List containing:
#' \describe{
#'   \item{counts}{Data frame of raw read counts (39,368 genes × 12 samples)}
#'   \item{pheno}{Data frame with sample metadata (sample_id, group)}
#'   \item{annot}{Data frame with gene annotation}
#' }
#'
#' @examples
#' full_data <- endo_load_gse201926()
#' dim(full_data$counts)
#' head(full_data$pheno)
#'
#' @export
endo_load_gse201926 <- function(sample_only = FALSE) {
    if (isTRUE(sample_only)) {
        utils::data("gse201926_sample", package = "endoSignatureR", envir = environment())
        return(get("gse201926_sample", envir = environment()))
    }

    # Helper to pick .gz if present, otherwise plain
    pick_extdata <- function(name) {
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

    # Load raw counts
    counts_file <- pick_extdata("gse201926_raw_counts.tsv")
    raw_counts <- readr::read_tsv(counts_file, show_col_types = FALSE)

    # Load phenotype data
    pheno_file <- pick_extdata("gse201926_series_matrix.txt")
    pheno <- endo_parse_metadata(pheno_file)

    # Load annotation (full)
    annot_file <- pick_extdata("gse201926_annotation.tsv")
    annot <- readr::read_tsv(annot_file, show_col_types = FALSE)

    return(list(
        counts = raw_counts,
        pheno = pheno,
        annot = annot
    ))
}

#' Parse GSE201926 series matrix metadata
#'
#' Parses the GEO series matrix file to extract sample metadata including
#' sample IDs and group labels (PS vs PIS).
#'
#' @param file_path Path to series matrix file
#' @return Data frame with sample metadata containing:
#' \describe{
#'   \item{sample_id}{Character vector of sample IDs (GSM6081173-1184)}
#'   \item{group}{Character vector of group labels (PS/PIS)}
#'   \item{title}{Character vector of sample titles}
#' }
#'
#' @examples
#' # Handle gzipped or plain file
#' pheno_file_gz <- system.file("extdata", "gse201926_series_matrix.txt.gz", package = "endoSignatureR")
#' if (nzchar(pheno_file_gz) && file.exists(pheno_file_gz)) {
#'     pheno <- endo_parse_metadata(pheno_file_gz)
#'     head(pheno)
#' }
#'
#' @export
endo_parse_metadata <- function(file_path) {
    # Check if file exists (gzipped or plain)
    if (!file.exists(file_path)) {
        stop("Series matrix file not found: ", file_path)
    }

    # Read the series matrix file (readLines handles .gz automatically)
    lines <- readLines(file_path)

    # Find the data table section
    data_start <- which(grepl("^!series_matrix_table_begin", lines))
    data_end <- which(grepl("^!series_matrix_table_end", lines))

    if (length(data_start) == 0 || length(data_end) == 0) {
        stop("Could not find data table section in series matrix file")
    }

    # Extract sample IDs from header
    header_line <- lines[data_start + 1]
    sample_ids <- unlist(strsplit(header_line, "\t"))[-1] # Remove first empty element
    
    # Clean sample IDs (remove quotes and trim whitespace)
    # The series matrix file may include quotes as part of the string
    sample_ids <- trimws(gsub('^"|"$', "", sample_ids))

    # Extract sample titles
    title_line_idx <- which(grepl("^!Sample_title", lines))
    if (length(title_line_idx) > 0) {
        title_line <- lines[title_line_idx[1]]
        titles <- unlist(strsplit(title_line, "\t"))[-1]
        # Clean titles as well
        titles <- trimws(gsub('^"|"$', "", titles))
    } else {
        titles <- rep(NA_character_, length(sample_ids))
    }

    # Determine group labels based on sample IDs
    # PIS samples: GSM6081173-1178, PS samples: GSM6081179-1184
    groups <- ifelse(grepl("GSM608117[3-8]", sample_ids), "PIS", "PS")

    # Create phenotype data frame
    pheno_data <- data.frame(
        sample_id = sample_ids,
        group = groups,
        title = titles,
        stringsAsFactors = FALSE
    )

    return(pheno_data)
}

#' Load gene annotation data
#'
#' @param minimal Logical; if TRUE, return the bundled minimal annotation
#'   (`gse201926_annot_min`) when available; otherwise read the full
#'   annotation from `inst/extdata`.
#' Loads the gene annotation file for GSE201926 containing Ensembl gene IDs,
#' gene symbols, descriptions, and other metadata.
#'
#' @return Data frame with gene annotation containing:
#' \describe{
#'   \item{GeneID}{Ensembl gene ID}
#'   \item{Symbol}{Gene symbol}
#'   \item{Description}{Gene description}
#'   \item{GeneType}{Type of gene (protein-coding, pseudo, etc.)}
#'   \item{EnsemblGeneID}{Full Ensembl gene ID}
#'   \item{ChrAcc, ChrStart, ChrStop}{Chromosomal coordinates}
#'   \item{GOFunction, GOProcess, GOComponent}{Gene Ontology terms}
#' }
#'
#' @examples
#' annot <- endo_load_annotation()
#' # Safely preview a few columns if present
#' cols <- intersect(c("GeneID", "Symbol", "Description"), names(annot))
#' head(annot[, cols, drop = FALSE])
#'
#' @export
endo_load_annotation <- function(minimal = TRUE) {
    if (isTRUE(minimal)) {
        # Prefer packaged minimal annotation if available
        if (exists("gse201926_annot_min", inherits = FALSE)) {
            return(get("gse201926_annot_min"))
        }
        if (exists("gse201926_annot_min", where = .GlobalEnv)) {
            return(get("gse201926_annot_min", envir = .GlobalEnv))
        }
    }

    pick_extdata <- function(name) {
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

    annot_file <- pick_extdata("gse201926_annotation.tsv")
    return(readr::read_tsv(annot_file, show_col_types = FALSE))
}

#' Load TPM normalized expression data
#'
#' Loads the TPM (Transcripts Per Million) normalized expression matrix
#' for GSE201926, useful for exploratory analysis and visualization.
#'
#' @return Data frame with TPM normalized expression values
#' (39,368 genes × 12 samples)
#'
#' @examples
#' tpm_data <- endo_load_tpm()
#' dim(tpm_data)
#' head(tpm_data[, 1:3])
#'
#' @export
endo_load_tpm <- function() {
    pick_extdata <- function(name) {
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
    tpm_file <- pick_extdata("gse201926_tpm_counts.tsv")
    return(readr::read_tsv(tpm_file, show_col_types = FALSE))
}
