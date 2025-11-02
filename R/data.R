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
#' **Data Structure:**
#' - `counts`: Matrix with numeric GeneIDs as rownames (not Ensembl format).
#'   Column names (sample IDs) are aligned with `pheno$sample_id` for seamless
#'   downstream analysis.
#' - `pheno`: Contains real sample metadata from the GEO series matrix,
#'   including `sample_id`, `group` (PS/PIS), and `title` columns.
#' - `annot`: Full annotation matching all genes in the counts matrix.
#'   Contains both numeric `GeneID` (matches counts rownames) and
#'   `EnsemblGeneID` (Ensembl format) columns.
#'
#' **Sample ID Alignment:**
#' The counts matrix column names are aligned with phenotype `sample_id` to
#' ensure compatibility with functions like `esr_analyzeDifferentialExpression()`
#' and other analysis workflows. Both use the same sample IDs from the series
#' matrix file.
#'
#' The original study investigated progestin sensitivity in endometrial lesions
#' using RNA-seq. Samples were collected from patients with endometrial
#' endometrioid carcinoma or endometrial atypical hyperplasia, and classified
#' as either progestin-sensitive (PS) or progestin-insensitive (PIS) based
#' on clinical response.
#'
#' **Source Files** (original names → package names):
#' - `GSE201926_raw_counts_GRCh38.p13_NCBI.tsv` → `inst/extdata/gse201926_raw_counts.tsv`
#' - `Human.GRCh38.p13.annot.tsv` → `inst/extdata/gse201926_annotation.tsv`
#' - `GSE201926_series_matrix.txt` → `inst/extdata/gse201926_series_matrix.txt`
#'
#' **Builder Script:** `data-raw/gse201926_sample.R` loads the full dataset from
#' these extdata files, aligns sample IDs between counts and phenotype, and
#' saves the complete dataset to `data/gse201926_sample.rda`.
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

#' GSE201926 Training Mini Dataset
#'
#' A medium-sized training subset of the GSE201926 endometrial lesion RNA-seq dataset
#' containing 800-1500 genes across 12 samples (6 PS + 6 PIS), designed for fast nested
#' cross-validation in vignettes and tests (nested CV <60s).
#'
#' @format A list with three components:
#' \describe{
#'   \item{counts}{Matrix of raw read counts (800-1500 genes × 12 samples).
#'     Rows are genes (numeric GeneIDs), columns are samples (GSM IDs).
#'     Genes selected based on variance (preferring protein-coding if available).}
#'   \item{pheno}{Data frame with sample metadata containing:
#'     \describe{
#'       \item{sample_id}{Character vector of sample IDs (GSM6081173-1184)}
#'       \item{group}{Character vector of group labels (PS/PIS)}
#'       \item{title}{Character vector of sample titles from GEO}
#'     }
#'   }
#'   \item{annot}{Data frame with gene annotation for selected genes containing:
#'     \describe{
#'       \item{GeneID}{Numeric gene ID (matches counts rownames)}
#'       \item{Symbol}{Gene symbol}
#'       \item{Description}{Gene description}
#'       \item{GeneType}{Type of gene (protein-coding, pseudo, etc.)}
#'       \item{EnsemblGeneID}{Full Ensembl gene ID (Ensembl format, e.g., ENSG00000000003)}
#'       \item{Other annotation columns}{ChrAcc, ChrStart, ChrStop, GO terms, etc.}
#'     }
#'   }
#' }
#'
#' @details
#' This dataset is a subset of the full GSE201926 dataset, designed specifically for
#' Mode 2 (Signature Validation) training workflows. The dataset balances:
#'
#' **Computational Efficiency**: 800-1500 genes × 12 samples enables nested CV in <60s
#' **Realistic Training Scenario**: Maintains class balance (6 PS + 6 PIS)
#' **Gene Selection**: Selected based on variance (most variable genes, preferring
#' protein-coding if available in annotation)
#'
#' **Data Structure:**
#' - `counts`: Matrix with numeric GeneIDs as rownames (not Ensembl format).
#'   Column names (sample IDs) are aligned with `pheno$sample_id` for seamless
#'   downstream analysis.
#' - `pheno`: Contains real sample metadata from the GEO series matrix,
#'   including `sample_id`, `group` (PS/PIS), and `title` columns.
#' - `annot`: Subset of full annotation matching selected genes only.
#'   Contains both numeric `GeneID` (matches counts rownames) and
#'   `EnsemblGeneID` (Ensembl format) columns.
#'
#' **Sample ID Alignment:**
#' The counts matrix column names are aligned with phenotype `sample_id` to
#' ensure compatibility with functions like `esr_trainEndometrialSignature()` and
#' other training workflows.
#'
#' **Builder Script:** `data-raw/build_gse201926_trainmini.R` loads the full dataset
#' from `inst/extdata` files, selects subset of genes (800-1500 based on variance),
#' and saves the training mini dataset to `data/gse201926_trainmini.rda`.
#'
#' @source \url{https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE201926}
#'
#' @references
#' The original GSE201926 dataset was published as part of research on
#' progestin sensitivity in endometrial lesions. For full details, see
#' the GEO database entry.
#'
#' @examples
#' # Load the training mini dataset
#' data(gse201926_trainmini)
#'
#' # Examine the data structure
#' str(gse201926_trainmini)
#'
#' # Check dimensions
#' dim(gse201926_trainmini$counts)
#'
#' # View sample metadata
#' gse201926_trainmini$pheno
#'
#' # Verify class balance
#' table(gse201926_trainmini$pheno$group)
#'
#' # View gene annotation
#' head(gse201926_trainmini$annot[, c("GeneID", "Symbol", "Description")])
"gse201926_trainmini"

#' Demo CV Splits Dataset
#'
#' Pre-computed cross-validation splits with fixed seeds for reproducible examples
#' in vignettes and tests. Used to avoid recomputing splits during demonstrations.
#'
#' @format A list with elements:
#' \describe{
#'   \item{outer_splits}{An `rsample` `vfold_cv` object containing outer CV splits.
#'     For `gse201926_trainmini`, this typically contains 3 outer folds with stratified
#'     sampling to maintain class balance.}
#'   \item{inner_splits}{A list of `rsample` `vfold_cv` objects, one for each outer fold.
#'     Each inner split contains CV folds for hyperparameter tuning within the outer
#'     training set (anti-leakage). Typically 2 inner folds per outer fold.}
#'   \item{outer_seed}{Integer; fixed random seed used for creating outer splits.
#'     Defaults to 12345.}
#'   \item{inner_seed}{Integer; fixed random seed used for creating inner splits.
#'     Defaults to 67890.}
#'   \item{n_outer_folds}{Integer; number of outer CV folds (typically 3).}
#'   \item{n_inner_folds}{Integer; number of inner CV folds per outer fold (typically 2).}
#'   \item{split_type}{Character; type of split used ("vfold_cv" for K-fold cross-validation).}
#'   \item{strata}{Character; column name used for stratified sampling ("group" for PS/PIS balance).}
#'   \item{created_date}{Date; date when folds were created.}
#' }
#'
#' @details
#' This dataset stores pre-computed CV splits for `gse201926_trainmini` to enable:
#'
#' **Reproducibility**: Fixed seeds ensure identical splits across runs
#' **Speed**: Pre-computed splits avoid recomputation in vignettes/tests
#' **Determinism**: Same seed always produces same splits
#' **Anti-leakage**: Inner splits are created only on training data (nested CV)
#'
#' **Usage:**
#' - Load folds: `data(folds_demo)`
#' - Access outer splits: `folds_demo$outer_splits`
#' - Access inner splits: `folds_demo$inner_splits[[1]]` (for first outer fold)
#' - Get training/test data: `rsample::training(folds_demo$outer_splits$splits[[1]])`
#'
#' **Split Structure:**
#' - **Outer splits**: Used for model selection/validation
#'   - 3-fold stratified CV for `gse201926_trainmini` (12 samples)
#'   - Each fold maintains class balance (PS/PIS)
#' - **Inner splits**: Used for hyperparameter tuning within each outer fold
#'   - 2-fold CV for each outer training set (8 samples after leaving 4 out)
#'   - Ensures test data never influences tuning decisions (anti-leakage)
#'
#' **Builder Script:** `data-raw/build_folds_demo.R` creates splits using `rsample::vfold_cv()`
#' with fixed seeds and saves to `data/folds_demo.rda`. Script must be run after
#' `gse201926_trainmini` is created.
#'
#' @source Created from `gse201926_trainmini` using `rsample::vfold_cv()` with fixed seeds.
#'
#' @references
#' See `vignette("mode2-signature-validation", package = "endoSignatureR")` for
#' examples of using `folds_demo` in training workflows.
#'
#' @examples
#' # Load demo folds
#' data(folds_demo)
#'
#' # Check structure
#' str(folds_demo)
#'
#' # Access outer splits
#' folds_demo$outer_splits
#'
#' # Get training/test data for first outer fold
#' library(rsample)
#' train_data <- training(folds_demo$outer_splits$splits[[1]])
#' test_data <- testing(folds_demo$outer_splits$splits[[1]])
#'
#' # Check class balance in first fold
#' table(train_data$group)
#' table(test_data$group)
#'
#' # Access inner splits for first outer fold
#' if (!is.null(folds_demo$inner_splits[[1]])) {
#'     inner_fold_1 <- folds_demo$inner_splits[[1]]
#'     head(inner_fold_1)
#' }
"folds_demo"
