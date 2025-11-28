#' Launch EndoSignatureR Shiny Application
#'
#' Launches the EndoSignatureR Shiny application, which provides a three-mode
#' workflow for analyzing endometrial bulk RNA-seq datasets:
#'
#' @details
#' \itemize{
#' \item \strong{Mode 1: Rapid Classification} - Apply pre-trained signature to
#' unlabeled samples for PS vs PIS predictions with confidence intervals. Results
#' include predictions table and signature visualization plots (coefficient lollipop,
#' stability bars).
#' \item \strong{Mode 2: Signature Validation} - Train and validate a new
#' signature on labeled cohorts using best practices (coming soon)
#' \item \strong{Mode 3: Visualization & Analysis} - Perform QC, exploratory
#' analysis, and differential expression visualization on endometrial data
#' }
#'
#' The Shiny application provides an interactive interface for:
#' \itemize{
#' \item Uploading endometrial RNA-seq data (counts, phenotype, annotation)
#' \item Loading demo data from the bundled \code{gse201926_sample} dataset
#' \item Performing quality control and exploratory data analysis (library size,
#' zeros, PCA)
#' \item Running differential expression analysis
#' \item Visualizing results (MA plots, volcano plots, heatmaps)
#' \item Exporting analysis results and plots
#' \item Applying pre-trained signature for classification (Mode 1)
#' \item Viewing signature visualization plots after classification (Mode 1)
#' }
#'
#' Mode 1 (Rapid Classification) and Mode 3 (Visualization & Analysis) are fully
#' implemented. Mode 2 will be available in a future release.
#'
#' @return A Shiny app object. The app will open in the default web browser.
#'
#' @examples
#' \dontrun{
#' # Launch the Shiny application
#' runEndoSignatureR()
#' }
#'
#' @references
#' Chang, W., Cheng, J., Allaire, J., Sievert, C., Schloerke, B., Xie, Y., ... &
#' Allen, J. (2024). shiny: Web Application Framework for R. R package version
#' 1.8.0. https://CRAN.R-project.org/package=shiny
#'
#' Wickham, H. (2016). ggplot2: Elegant Graphics for Data Analysis.
#' Springer-Verlag New York. ISBN 978-3-319-24277-4.
#'
#' Gu, Z., Eils, R., & Schlesner, M. (2016). Complex heatmaps reveal patterns
#' and correlations in multidimensional genomic data. Bioinformatics, 32(18),
#' 2847-2849.
#'
#' Ritchie, M. E., Phipson, B., Wu, D., Hu, Y., Law, C. W., Shi, W., & Smyth, G. K.
#' (2015). limma powers differential expression analyses for RNA-sequencing and
#' microarray studies. Nucleic Acids Research, 43(7), e47.
#'
#' @importFrom shiny shinyApp
#' @export
runEndoSignatureR <- function() {
    app_file <- system.file(
      "shiny-scripts", "app.r", package = "endoSignatureR"
    )

    if (!file.exists(app_file)) {
        stop("Shiny app file not found. Expected location: ", app_file)
    }

    # Source the app file in a new environment to get ui and server
    app_env <- new.env()
    source(app_file, local = app_env)

    # Get ui and server from the environment
    if (exists("ui", envir = app_env) && exists("server", envir = app_env)) {
        ui <- get("ui", envir = app_env)
        server <- get("server", envir = app_env)
        return(shinyApp(ui = ui, server = server))
    } else {
        stop("Could not find ui and server objects in app.r file")
    }
}

# [END]
