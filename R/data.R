#' GSE201926 Sample Dataset
#'
#' The complete GSE201926 endometrial lesion RNA-seq dataset containing
#' all genes across 12 samples (6 progestin-sensitive (PS)
#' and 6 progestin-insensitive (PIS) patients).
#'
#' @format A list with three components:
#' \describe{
#'   \item{counts}{Matrix of raw read counts (39,368 genes × 12 samples).
#'     Rows are genes (numeric GeneIDs matching the raw counts file), columns are samples (GSM IDs).}
#'   \item{pheno}{Data frame with sample metadata containing:
#'     \describe{
#'       \item{sample_id}{Character vector of sample IDs (GSM6081173-1184)}
#'       \item{group}{Character vector of group labels (PS/PIS)}
#'       \item{title}{Character vector of sample titles from GEO}
#'     }
#'   }
#'   \item{annot}{Data frame with gene annotation containing:
#'     \describe{
#'       \item{GeneID}{Numeric gene ID (matches counts rownames)}
#'       \item{Symbol}{Gene symbol}
#'       \item{Description}{Gene description}
#'       \item{GeneType}{Type of gene (protein-coding, pseudo, etc.)}
#'       \item{EnsemblGeneID}{Full Ensembl gene ID (Ensembl format, e.g., ENSG00000000003)}
#'       \item{ChrAcc, ChrStart, ChrStop}{Chromosomal coordinates}
#'       \item{GOFunction, GOProcess, GOComponent}{Gene Ontology terms}
#'     }
#'   }
#' }
#'
#' @details
#' This dataset contains the full GSE201926 dataset loaded from package
#' extdata files, providing the complete gene expression matrix and annotation
#' for realistic analysis workflows in examples, tests, and vignettes.
#'
#' The original study investigated progestin sensitivity in endometrial lesions
#' using RNA-seq. Samples were collected from patients with endometrial
#' endometrioid carcinoma or endometrial atypical hyperplasia, and classified
#' as either progestin-sensitive (PS) or progestin-insensitive (PIS) based
#' on clinical response.
#'
#' Source files used to generate this dataset (original names → package names):
#' - `GSE201926_raw_counts_GRCh38.p13_NCBI.tsv` → `inst/extdata/gse201926_raw_counts.tsv`
#' - `Human.GRCh38.p13.annot.tsv` → `inst/extdata/gse201926_annotation.tsv`
#' - `GSE201926_series_matrix.txt` → `inst/extdata/gse201926_series_matrix.txt`
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201926}
#'
#' @references
#' The original GSE201926 dataset was published as part of research on
#' progestin sensitivity in endometrial lesions. For full details, see
#' the GEO database entry.
#'
#' @examples
#' # Load the sample dataset
#' data(gse201926_sample)
#'
#' # Examine the data structure
#' str(gse201926_sample)
#'
#' # Check dimensions
#' dim(gse201926_sample$counts)
#'
#' # View sample metadata
#' gse201926_sample$pheno
#'
#' # View gene annotation
#' head(gse201926_sample$annot[, c("GeneID", "Symbol", "Description")])
#'
#' # Check group distribution
#' table(gse201926_sample$pheno$group)
#'
"gse201926_sample"

#' GSE201926 Minimal Annotation
#'
#' A minimal gene annotation table for GSE201926 suitable for examples and
#' vignettes. Contains only essential fields for mapping and display.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{gene_id}{Ensembl gene ID}
#'   \item{symbol}{HGNC gene symbol}
#'   \item{biotype}{Gene biotype (e.g., protein_coding)}
#' }
#'
#' @details
#' Built from `inst/extdata/gse201926_annotation.tsv(.gz)` in `data-raw/build_gse201926_annot_min.R`.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201926}
#'
#' @examples
#' # data(gse201926_annot_min) # loaded on package attach if included
#' head(gse201926_annot_min)
"gse201926_annot_min"

#' GSE201926 Minimal Phenotype
#'
#' A minimal phenotype/metadata table for GSE201926 suitable for examples and
#' vignettes. Contains only sample identifiers, group labels, and optional batch.
#'
#' @format A data frame with columns:
#' \describe{
#'   \item{sample_id}{GSM sample identifier}
#'   \item{group}{PS or PIS}
#'   \item{batch}{Optional batch label if available}
#' }
#'
#' @details
#' Built from `inst/extdata/gse201926_series_matrix.txt(.gz)` in `data-raw/build_gse201926_pheno_min.R`.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201926}
#'
#' @examples
#' # data(gse201926_pheno_min)
#' head(gse201926_pheno_min)
"gse201926_pheno_min"
