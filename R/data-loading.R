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
#' @references
#' R Core Team (2025). R: A language and environment for statistical
#' computing. R Foundation for Statistical Computing, Vienna, Austria.
#' <https://www.R-project.org/>.
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
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
#'
#' Edgar, R., Domrachev, M., & Lash, A. E. (2002). Gene Expression
#' Omnibus: NCBI gene expression and hybridization array data repository.
#' Nucleic Acids Research, 30(1), 207-210.
#' <https://doi.org/10.1093/nar/30.1.207>.
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
#' pheno_file_gz <- system.file("extdata",
#'     "gse201926_series_matrix.txt.gz",
#'     package = "endoSignatureR"
#' )
#' if (nzchar(pheno_file_gz) && file.exists(pheno_file_gz)) {
#'     pheno <- endo_parse_metadata(pheno_file_gz)
#'     head(pheno)
#' }
#'
#' @references
#' Edgar, R., Domrachev, M., & Lash, A. E. (2002). Gene Expression
#' Omnibus: NCBI gene expression and hybridization array data repository.
#' Nucleic Acids Research, 30(1), 207-210.
#' <https://doi.org/10.1093/nar/30.1.207>.
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
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
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
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
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

#' Load Counts Data from File
#'
#' Loads raw counts data from TSV or CSV file and converts it to the matrix format
#' expected by package functions (genes × samples).
#'
#' @param file_path Path to counts file (TSV or CSV). May be compressed (.gz).
#' @param gene_id_col Character scalar; name of the gene ID column. Defaults to "GeneID".
#' @param sample_id_col Optional character scalar; name of sample ID column. If NULL, uses column names after gene_id_col.
#' @param transpose Logical; if TRUE, input is samples × genes and will be transposed. Defaults to FALSE.
#' @param ... Additional arguments passed to `readr::read_tsv()` or `readr::read_csv()`.
#'
#' @return Matrix with genes as rows (rownames = GeneIDs) and samples as columns (colnames = sample IDs).
#'
#' @details
#' This function loads raw counts data from user files and converts it to the format expected
#' by package functions. The function:
#' - Auto-detects TSV vs CSV format based on file extension
#' - Handles compressed files (.gz) automatically
#' - Converts data to numeric matrix format
#' - Sets rownames to GeneIDs and colnames to sample IDs
#' - Cleans sample IDs (removes quotes, trims whitespace)
#'
#' @examples
#' \dontrun{
#' # Load counts from TSV file
#' counts <- esr_loadCountsFromFile("my_counts.tsv")
#' dim(counts)
#' head(rownames(counts)) # Gene IDs
#' head(colnames(counts)) # Sample IDs
#'
#' # Load from compressed CSV file
#' counts <- esr_loadCountsFromFile("my_counts.csv.gz", gene_id_col = "EnsemblID")
#' }
#'
#' @import readr
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
#' @export
esr_loadCountsFromFile <- function(file_path, gene_id_col = "GeneID",
                                   sample_id_col = NULL, transpose = FALSE, ...) {
    # Check if file exists
    if (!file.exists(file_path)) {
        stop("Counts file not found: ", file_path)
    }

    # Auto-detect file format (TSV or CSV)
    file_ext <- tools::file_ext(gsub("\\.gz$", "", file_path))
    is_tsv <- file_ext %in% c("tsv", "txt", "") || file_ext == ""

    # Read file
    if (is_tsv) {
        counts_df <- readr::read_tsv(file_path, show_col_types = FALSE, ...)
    } else {
        counts_df <- readr::read_csv(file_path, show_col_types = FALSE, ...)
    }

    # Check that gene_id_col exists
    if (!gene_id_col %in% names(counts_df)) {
        stop(paste0("Gene ID column '", gene_id_col, "' not found in counts file. Available columns: ", paste(names(counts_df), collapse = ", ")))
    }

    # Extract gene IDs and sample IDs
    gene_ids <- counts_df[[gene_id_col]]
    sample_cols <- setdiff(names(counts_df), gene_id_col)

    # Convert to matrix
    counts_matrix <- as.matrix(counts_df[, sample_cols, drop = FALSE])

    # Transpose if needed
    if (transpose) {
        counts_matrix <- t(counts_matrix)
    }

    # Set rownames to gene IDs
    rownames(counts_matrix) <- gene_ids

    # Set colnames to sample IDs (clean if needed)
    if (is.null(sample_id_col)) {
        sample_ids <- sample_cols
    } else {
        if (!sample_id_col %in% names(counts_df)) {
            stop(paste0("Sample ID column '", sample_id_col, "' not found in counts file"))
        }
        sample_ids <- counts_df[[sample_id_col]]
    }

    # Clean sample IDs (remove quotes, trim whitespace)
    sample_ids <- trimws(gsub('^"|"$', "", sample_ids))
    colnames(counts_matrix) <- sample_ids

    # Validate that all values are numeric
    if (!is.numeric(counts_matrix)) {
        warning("Counts matrix contains non-numeric values. Converting to numeric may cause issues.")
        counts_matrix <- apply(counts_matrix, 2, as.numeric)
    }

    # Check for missing values
    if (any(is.na(counts_matrix))) {
        warning("Counts matrix contains missing values (NA). This may cause issues in downstream analysis.")
    }

    return(counts_matrix)
}

#' Load Phenotype/Metadata from File
#'
#' Loads phenotype/metadata from series matrix file (GEO format) or simple TSV/CSV file
#' and converts it to the format expected by package functions.
#'
#' @param file_path Path to phenotype file (series matrix, TSV, or CSV). May be compressed (.gz).
#' @param sample_id_col Optional character scalar; name of sample ID column. Auto-detected if NULL.
#' @param group_col Optional character scalar; name of group column. Auto-detected if NULL.
#' @param format Character scalar; file format: "auto" (default, auto-detect), "series_matrix", "tsv", or "csv".
#' @param ... Additional arguments passed to `readr::read_tsv()` or `readr::read_csv()`.
#'
#' @return Data frame with required columns `sample_id` and `group`, plus any additional metadata columns.
#'
#' @details
#' This function loads phenotype/metadata from user files and converts it to the format expected
#' by package functions. The function:
#' - Auto-detects series matrix format vs simple TSV/CSV
#' - Handles compressed files (.gz) automatically
#' - Auto-detects column names if not specified
#' - Maps columns to expected format (`sample_id`, `group`)
#' - Validates required columns exist
#'
#' For series matrix format, uses `endo_parse_metadata()` internally.
#' For simple TSV/CSV, auto-detects common column name patterns:
#' - Sample ID: `sample_id`, `Sample_ID`, `sample`, `Sample`, `GSM_ID`, etc.
#' - Group: `group`, `Group`, `label`, `Label`, `class`, `Class`, etc.
#'
#' @examples
#' \dontrun{
#' # Load from series matrix file (GEO format)
#' pheno <- esr_loadPhenoFromFile("my_series_matrix.txt", format = "series_matrix")
#'
#' # Load from simple TSV file
#' pheno <- esr_loadPhenoFromFile("my_metadata.tsv")
#'
#' # Load with explicit column names
#' pheno <- esr_loadPhenoFromFile("my_metadata.csv",
#'     sample_id_col = "SampleID",
#'     group_col = "Class"
#' )
#' }
#'
#' @references
#' Edgar, R., Domrachev, M., & Lash, A. E. (2002). Gene Expression
#' Omnibus: NCBI gene expression and hybridization array data repository.
#' Nucleic Acids Research, 30(1), 207-210.
#' <https://doi.org/10.1093/nar/30.1.207>.
#'
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
#' @export
esr_loadPhenoFromFile <- function(file_path, sample_id_col = NULL, group_col = NULL,
                                  format = c("auto", "series_matrix", "tsv", "csv"), ...) {
    format <- match.arg(format)

    # Check if file exists
    if (!file.exists(file_path)) {
        stop("Phenotype file not found: ", file_path)
    }

    # Auto-detect format if needed
    if (format == "auto") {
        # Check if it's a series matrix file
        lines <- readLines(file_path, n = 10)
        if (any(grepl("!series_matrix_table_begin", lines))) {
            format <- "series_matrix"
        } else {
            # Auto-detect TSV vs CSV
            file_ext <- tools::file_ext(gsub("\\.gz$", "", file_path))
            if (file_ext %in% c("tsv", "txt", "")) {
                format <- "tsv"
            } else {
                format <- "csv"
            }
        }
    }

    # Handle series matrix format
    if (format == "series_matrix") {
        pheno <- endo_parse_metadata(file_path)
        # Ensure required columns exist
        if (!"sample_id" %in% names(pheno)) {
            stop("Series matrix parsing failed to produce 'sample_id' column")
        }
        if (!"group" %in% names(pheno)) {
            stop("Series matrix parsing failed to produce 'group' column")
        }
        return(pheno)
    }

    # Handle simple TSV/CSV format
    is_tsv <- (format == "tsv")
    if (is_tsv) {
        pheno_df <- readr::read_tsv(file_path, show_col_types = FALSE, ...)
    } else {
        pheno_df <- readr::read_csv(file_path, show_col_types = FALSE, ...)
    }

    # Auto-detect column names if not specified
    if (is.null(sample_id_col)) {
        # Try common sample ID column names
        sample_id_patterns <- c("sample_id", "Sample_ID", "sample", "Sample", "GSM_ID", "GSM", "sample_id", "sampleID")
        sample_id_col <- NULL
        for (pattern in sample_id_patterns) {
            if (pattern %in% names(pheno_df)) {
                sample_id_col <- pattern
                break
            }
        }
        if (is.null(sample_id_col)) {
            # Use first column as fallback
            sample_id_col <- names(pheno_df)[1]
            warning(paste0("Sample ID column not found. Using first column '", sample_id_col, "' as sample_id"))
        }
    }

    if (is.null(group_col)) {
        # Try common group column names
        group_patterns <- c("group", "Group", "label", "Label", "class", "Class", "phenotype", "Phenotype")
        group_col <- NULL
        for (pattern in group_patterns) {
            if (pattern %in% names(pheno_df)) {
                group_col <- pattern
                break
            }
        }
        if (is.null(group_col)) {
            stop("Group column not found. Please specify 'group_col' parameter or ensure file contains a column named 'group', 'Group', 'label', or 'Label'")
        }
    }

    # Validate columns exist
    if (!sample_id_col %in% names(pheno_df)) {
        stop(paste0("Sample ID column '", sample_id_col, "' not found in phenotype file. Available columns: ", paste(names(pheno_df), collapse = ", ")))
    }

    if (!group_col %in% names(pheno_df)) {
        stop(paste0("Group column '", group_col, "' not found in phenotype file. Available columns: ", paste(names(pheno_df), collapse = ", ")))
    }

    # Create output data frame with required columns
    pheno <- data.frame(
        sample_id = trimws(gsub('^"|"$', "", as.character(pheno_df[[sample_id_col]]))),
        group = trimws(gsub('^"|"$', "", as.character(pheno_df[[group_col]]))),
        stringsAsFactors = FALSE
    )

    # Add any additional columns
    other_cols <- setdiff(names(pheno_df), c(sample_id_col, group_col))
    if (length(other_cols) > 0) {
        for (col in other_cols) {
            pheno[[col]] <- pheno_df[[col]]
        }
    }

    # Validate group labels (warn if not PS/PIS)
    unique_groups <- unique(pheno$group)
    if (!all(unique_groups %in% c("PS", "PIS"))) {
        warning(paste0("Group labels found: ", paste(unique_groups, collapse = ", "), ". Expected 'PS' or 'PIS'. Please verify group labels are correct."))
    }

    return(pheno)
}

#' Load Annotation Data from File
#'
#' Loads gene annotation from TSV or CSV file and converts it to the format expected
#' by package functions.
#'
#' @param file_path Path to annotation file (TSV or CSV). May be compressed (.gz).
#' @param gene_id_col Character scalar; name of the gene ID column. Defaults to "GeneID".
#' @param ... Additional arguments passed to `readr::read_tsv()` or `readr::read_csv()`.
#'
#' @return Data frame with gene annotation. Must contain `GeneID` column matching counts rownames.
#'
#' @details
#' This function loads gene annotation from user files. The function:
#' - Auto-detects TSV vs CSV format based on file extension
#' - Handles compressed files (.gz) automatically
#' - Validates GeneID column exists
#' - Preserves all annotation columns
#'
#' @examples
#' \dontrun{
#' # Load annotation from TSV file
#' annot <- esr_loadAnnotFromFile("my_annotation.tsv")
#' head(annot[, c("GeneID", "Symbol", "Description")])
#'
#' # Load with different GeneID column name
#' annot <- esr_loadAnnotFromFile("my_annotation.csv", gene_id_col = "EnsemblID")
#' }
#'
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
#' @export
esr_loadAnnotFromFile <- function(file_path, gene_id_col = "GeneID", ...) {
    # Check if file exists
    if (!file.exists(file_path)) {
        stop("Annotation file not found: ", file_path)
    }

    # Auto-detect file format (TSV or CSV)
    file_ext <- tools::file_ext(gsub("\\.gz$", "", file_path))
    is_tsv <- file_ext %in% c("tsv", "txt", "") || file_ext == ""

    # Read file
    if (is_tsv) {
        annot_df <- readr::read_tsv(file_path, show_col_types = FALSE, ...)
    } else {
        annot_df <- readr::read_csv(file_path, show_col_types = FALSE, ...)
    }

    # Check that gene_id_col exists
    if (!gene_id_col %in% names(annot_df)) {
        stop(paste0("Gene ID column '", gene_id_col, "' not found in annotation file. Available columns: ", paste(names(annot_df), collapse = ", ")))
    }

    # Rename gene_id_col to GeneID for consistency
    if (gene_id_col != "GeneID") {
        names(annot_df)[names(annot_df) == gene_id_col] <- "GeneID"
    }

    return(annot_df)
}

#' Load Data from Files (High-Level Convenience Function)
#'
#' Loads counts, phenotype, and annotation files and returns them in the standardized format
#' expected by package functions. This function handles ID alignment between files automatically.
#'
#' @param counts_file Path to counts file (TSV or CSV). Required.
#' @param pheno_file Path to phenotype file (series matrix, TSV, or CSV). Optional (NULL for unlabeled data).
#' @param annot_file Path to annotation file (TSV or CSV). Optional.
#' @param counts_gene_id_col Character scalar; gene ID column name in counts file. Defaults to "GeneID".
#' @param pheno_sample_id_col Optional character scalar; sample ID column name in phenotype file. Auto-detected if NULL.
#' @param pheno_group_col Optional character scalar; group column name in phenotype file. Auto-detected if NULL.
#' @param annot_gene_id_col Character scalar; gene ID column name in annotation file. Defaults to "GeneID".
#' @param pheno_format Character scalar; phenotype file format: "auto" (default, auto-detect), "series_matrix", "tsv", or "csv".
#' @param align_ids Logical; if TRUE (default), align sample/gene IDs between files automatically. If FALSE, return files as-is.
#' @param validate Logical; if TRUE (default), run `esr_validateEndometrial()` on loaded data. If FALSE, skip validation.
#' @param ... Additional arguments passed to file reading functions.
#'
#' @return List containing:
#' \describe{
#'   \item{counts}{Matrix with genes as rows, samples as columns}
#'   \item{pheno}{Data frame with `sample_id` and `group` columns (or NULL if pheno_file not provided)}
#'   \item{annot}{Data frame with `GeneID` column (or NULL if annot_file not provided)}
#'   \item{issues}{Data frame with validation issues (if validate = TRUE)}
#' }
#'
#' @details
#' This is a high-level convenience function that loads all data files and aligns IDs automatically.
#' The function:
#' - Loads counts, phenotype, and annotation files using individual loading functions
#' - Aligns sample IDs between counts and pheno (reorders columns to match)
#' - Aligns gene IDs between counts and annot (reorders rows to match)
#' - Optionally validates data structure using `esr_validateEndometrial()`
#' - Returns data in format matching `gse201926_sample` structure
#'
#' @examples
#' \dontrun{
#' # Load unlabeled data (Mode 1)
#' user_data <- esr_loadFromFiles(
#'     counts_file = "my_counts.tsv",
#'     pheno_file = NULL,
#'     validate = FALSE
#' )
#'
#' # Load labeled data (Mode 2)
#' user_data <- esr_loadFromFiles(
#'     counts_file = "my_counts.tsv",
#'     pheno_file = "my_metadata.tsv",
#'     annot_file = "my_annotation.tsv",
#'     validate = TRUE
#' )
#'
#' # Check validation issues
#' if (nrow(user_data$issues) > 0) {
#'     print(user_data$issues)
#' }
#' }
#'
#' @references
#' Wickham, H., Hester, J., & Bryan, J. (2024). readr: Read rectangular
#' text data. R package. <https://CRAN.R-project.org/package=readr>.
#'
#' Edgar, R., Domrachev, M., & Lash, A. E. (2002). Gene Expression
#' Omnibus: NCBI gene expression and hybridization array data repository.
#' Nucleic Acids Research, 30(1), 207-210.
#' <https://doi.org/10.1093/nar/30.1.207>.
#' @export
esr_loadFromFiles <- function(counts_file, pheno_file = NULL, annot_file = NULL,
                              counts_gene_id_col = "GeneID",
                              pheno_sample_id_col = NULL, pheno_group_col = NULL,
                              annot_gene_id_col = "GeneID",
                              pheno_format = c("auto", "series_matrix", "tsv", "csv"),
                              align_ids = TRUE, validate = TRUE, ...) {
    pheno_format <- match.arg(pheno_format)

    # Load counts file
    counts <- esr_loadCountsFromFile(
        file_path = counts_file,
        gene_id_col = counts_gene_id_col,
        ...
    )

    # Load phenotype file (if provided)
    pheno <- NULL
    if (!is.null(pheno_file)) {
        pheno <- esr_loadPhenoFromFile(
            file_path = pheno_file,
            sample_id_col = pheno_sample_id_col,
            group_col = pheno_group_col,
            format = pheno_format,
            ...
        )
    }

    # Load annotation file (if provided)
    annot <- NULL
    if (!is.null(annot_file)) {
        annot <- esr_loadAnnotFromFile(
            file_path = annot_file,
            gene_id_col = annot_gene_id_col,
            ...
        )
    }

    # Align IDs if requested
    if (align_ids) {
        # Align sample IDs between counts and pheno
        if (!is.null(pheno)) {
            counts_sample_ids <- colnames(counts)
            pheno_sample_ids <- pheno$sample_id

            # Check for matching sample IDs
            common_samples <- intersect(counts_sample_ids, pheno_sample_ids)
            if (length(common_samples) == 0) {
                warning("No matching sample IDs found between counts and phenotype. Sample IDs may not align correctly.")
            } else {
                # Reorder counts columns to match pheno order
                if (length(common_samples) < length(counts_sample_ids)) {
                    warning(paste0("Only ", length(common_samples), " out of ", length(counts_sample_ids), " sample IDs match between counts and phenotype. Subsetting to common samples."))
                    counts <- counts[, common_samples, drop = FALSE]
                }
                # Reorder pheno to match counts order
                pheno <- pheno[match(common_samples, pheno$sample_id), , drop = FALSE]
                # Reorder counts to match pheno order
                counts <- counts[, pheno$sample_id, drop = FALSE]
            }
        }

        # Align gene IDs between counts and annot
        if (!is.null(annot)) {
            counts_gene_ids <- rownames(counts)
            annot_gene_ids <- as.character(annot$GeneID) # Convert to character for matching

            # Check for matching gene IDs
            common_genes <- intersect(counts_gene_ids, annot_gene_ids)
            if (length(common_genes) == 0) {
                warning("No matching gene IDs found between counts and annotation. Gene IDs may not align correctly.")
            } else {
                # Subset counts to common genes (preserve counts order)
                if (length(common_genes) < length(counts_gene_ids)) {
                    warning(paste0("Only ", length(common_genes), " out of ", length(counts_gene_ids), " gene IDs match between counts and annotation. Subsetting to common genes."))
                }
                counts <- counts[common_genes, , drop = FALSE]

                # Convert annot$GeneID to character for matching
                annot$GeneID <- as.character(annot$GeneID)

                # Reorder annot to match counts order (use rownames(counts) as the ordering)
                annot <- annot[match(rownames(counts), annot$GeneID), , drop = FALSE]
            }
        }
    }

    # Validate data structure if requested
    issues <- NULL
    if (validate) {
        validated <- esr_validateEndometrial(
            X = counts,
            pheno = pheno,
            annot = annot
        )
        issues <- validated$issues
        # Use validated data (may have been cleaned)
        counts <- validated$X
        pheno <- validated$pheno
        annot <- validated$annot
    }

    # Return data in standard format
    result <- list(
        counts = counts,
        pheno = pheno,
        annot = annot
    )

    # Add issues if validation was performed
    if (!is.null(issues)) {
        result$issues <- issues
    }

    return(result)
}
