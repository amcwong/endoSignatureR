# EndoSignatureR Shiny Application
# Three-Mode Workflow: Mode 1 (Rapid Classification), Mode 2 (Signature Validation), Mode 3 (Visualization & Analysis)

library(shiny)

# Define UI
ui <- fluidPage(
    titlePanel("EndoSignatureR: Endometrial RNA-seq Analysis"),

    # Main tabset with three modes
    tabsetPanel(
        id = "modeTabs",

        # Landing Page (Welcome/Getting Started)
        tabPanel(
            "Welcome",
            value = "welcome",
            div(
                style = "max-width: 900px; margin: 0 auto; padding: 20px;",
                h2("Welcome to EndoSignatureR", style = "text-align: center; color: #2c3e50; margin-bottom: 30px;"),
                p("EndoSignatureR is an R package and Shiny application designed to make it easy for researchers to analyze and compare small endometrial bulk RNA-seq datasets using a standardized, reproducible workflow.", style = "font-size: 16px; line-height: 1.6; text-align: center; margin-bottom: 40px;"),
                h3("What the App Does", style = "color: #34495e; margin-top: 40px;"),
                p("EndoSignatureR provides three complementary workflows for endometrial RNA-seq analysis:"),
                tags$ul(
                    tags$li(strong("Rapid Classification (Mode 1):"), " Apply a pre-trained PS vs PIS signature to unlabeled endometrial samples to obtain predictions with confidence intervals."),
                    tags$li(strong("Signature Validation (Mode 2):"), " Train and validate a new signature on labeled cohorts using best practices (nested cross-validation, LASSO, calibration) and compare to the pre-trained signature."),
                    tags$li(strong("Visualization & Analysis (Mode 3):"), " Perform standalone quality control, exploratory analysis, and differential expression visualization on endometrial data without requiring training or classification.")
                ),
                br(),
                h3("What It Can Do", style = "color: #34495e; margin-top: 40px;"),
                p("The three-mode workflow allows you to:"),
                tags$ul(
                    tags$li("Classify unlabeled endometrial samples using a pre-trained signature"),
                    tags$li("Train and validate signatures on your own labeled data"),
                    tags$li("Explore and visualize endometrial RNA-seq data with comprehensive QC and DE analysis"),
                    tags$li("Export results in multiple formats for further analysis or reporting"),
                    tags$li("Compare signature performance and visualize results with publication-ready plots")
                ),
                br(),
                h3("How to Use It", style = "color: #34495e; margin-top: 40px;"),
                p("Getting started is straightforward:"),
                tags$ul(
                    tags$li(strong("Mode 1 and Mode 3"), " can be run with just your data - no prerequisites required. Simply upload your data and start analyzing."),
                    tags$li(strong("Some functionality in Mode 1"), " requires Mode 2 completion. For example, using a user-trained signature (instead of the pre-trained signature) requires first training a signature in Mode 2."),
                    tags$li(strong("Additional information"), " is available in dropdown menus next to export options, providing ideas about data exploration and detailed descriptions of what each export contains.")
                ),
                br(),
                p(em("Tip: Each export option has an information dropdown (ℹ️) that explains what's in the export, what you can do with it, and format specifications.")),
                br(),
                h3("Get Started", style = "color: #34495e; margin-top: 40px; text-align: center;"),
                div(
                    style = "text-align: center; margin: 30px 0;",
                    actionButton("navToMode1", "Go to Mode 1: Rapid Classification",
                        class = "btn-primary btn-lg",
                        style = "margin: 10px; padding: 15px 30px; font-size: 16px;"
                    ),
                    br(),
                    actionButton("navToMode2", "Go to Mode 2: Signature Validation",
                        class = "btn-primary btn-lg",
                        style = "margin: 10px; padding: 15px 30px; font-size: 16px;"
                    ),
                    br(),
                    actionButton("navToMode3", "Go to Mode 3: Visualization & Analysis",
                        class = "btn-primary btn-lg",
                        style = "margin: 10px; padding: 15px 30px; font-size: 16px;"
                    )
                ),
                br(),
                hr(),
                p("For detailed information about our methodology, visit the ", strong("Methodology"), " tab.", style = "text-align: center; color: #7f8c8d;")
            )
        ),

        # Mode 1: Rapid Classification
        tabPanel(
            "Mode 1: Rapid Classification",
            sidebarLayout(
                sidebarPanel(
                    h3("Mode 1: Rapid Classification"),
                    p("Apply pre-trained signature to unlabeled endometrial samples for PS vs PIS predictions."),
                    br(),

                    # Demo data section
                    h4("Demo Data"),
                    p("Load or download unlabeled demo dataset (gse201926_sample, counts only):"),
                    actionButton("loadDemoMode1", "Load Demo Data", class = "btn-primary"),
                    br(), br(),
                    downloadButton("downloadDemoMode1", "Download Demo Data", class = "btn-secondary"),
                    br(), br(),
                    hr(),

                    # Data upload section
                    h4("Upload Your Data"),
                    p("Upload unlabeled endometrial RNA-seq data (counts matrix only, no phenotype labels):"),
                    fileInput("countsFileMode1", "Counts Matrix (TSV/CSV)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "Genes × Samples"
                    ),
                    fileInput("annotFileMode1", "Annotation Data (TSV/CSV, optional)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "Gene annotations"
                    ),
                    br(),

                    # Data status
                    verbatimTextOutput("dataStatusMode1"),
                    br(),

                    # Signature selection
                    h4("Signature Selection"),
                    p("Select which signature to use for classification:"),
                    uiOutput("signatureSelectionUI"),
                    br(),

                    # Upload signature files (for user-provided signatures)
                    h4("Upload Signature Files (Optional)"),
                    p("If you downloaded signature files from Mode 2, you can upload them here:"),
                    fileInput("signatureCSVFile", "Signature CSV File",
                        accept = c(".csv", ".tsv"),
                        placeholder = "endometrial_signature.csv"
                    ),
                    fileInput("signatureJSONFile", "Signature JSON File",
                        accept = c(".json"),
                        placeholder = "endometrial_recipe.json"
                    ),
                    fileInput("signatureStabilityFile", "Stability CSV File (Optional)",
                        accept = c(".csv", ".tsv"),
                        placeholder = "endometrial_stability.csv"
                    ),
                    br(),

                    # Signature information
                    h4("Signature Information"),
                    verbatimTextOutput("signatureInfo"),
                    br(),

                    # Classification controls
                    h4("Classification"),
                    p("Apply selected signature to classify samples:"),
                    numericInput("thresholdMode1", "Classification Threshold",
                        value = 0.5, min = 0, max = 1, step = 0.05
                    ),
                    checkboxInput("confidenceMode1", "Include Confidence Intervals", value = TRUE),
                    br(),
                    actionButton("runClassification", "Classify Samples", class = "btn-success"),
                    br(), br(),

                    # Status messages
                    verbatimTextOutput("statusMessageMode1")
                ),
                mainPanel(
                    tabsetPanel(
                        id = "mode1Tabs",

                        # Results Tab
                        tabPanel(
                            "Results",
                            h3("Classification Results"),
                            verbatimTextOutput("signatureSourceInfo"),
                            br(),
                            conditionalPanel(
                                condition = "output.classificationRun == false",
                                p("Click 'Classify Samples' in the sidebar to run classification.")
                            ),
                            conditionalPanel(
                                condition = "output.classificationRun == true",
                                h4("Summary"),
                                verbatimTextOutput("classificationSummary"),
                                br(),
                                h4("Predictions Table"),
                                p("Download full results using the button below."),
                                div(
                                    style = "display: flex; align-items: center; gap: 10px; margin-bottom: 10px;",
                                    downloadButton("downloadPredictions", "Download Predictions (CSV)"),
                                    tags$details(
                                        tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                                        div(
                                            style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                            h5("What's in this export:"),
                                            p("The CSV file contains classification results with the following columns:"),
                                            tags$ul(
                                                tags$li(strong("sample_id:"), " Unique identifier for each sample"),
                                                tags$li(strong("predicted_class:"), " Predicted class (PS or PIS)"),
                                                tags$li(strong("probability:"), " Predicted probability for the assigned class"),
                                                tags$li(strong("confidence_lower:"), " Lower bound of confidence interval (if confidence intervals were requested)"),
                                                tags$li(strong("confidence_upper:"), " Upper bound of confidence interval (if confidence intervals were requested)")
                                            ),
                                            br(),
                                            h5("What you can do with this export:"),
                                            tags$ul(
                                                tags$li("Use predictions for downstream analysis in R or other tools"),
                                                tags$li("Create custom visualizations or reports"),
                                                tags$li("Share results with collaborators"),
                                                tags$li("Integrate predictions into clinical or research workflows"),
                                                tags$li("Compare predictions across different signatures or thresholds")
                                            ),
                                            br(),
                                            h5("Format specifications:"),
                                            p("CSV format (comma-separated values), UTF-8 encoding. Can be opened in Excel, R, Python, or any text editor.")
                                        )
                                    )
                                ),
                                br(), br(),
                                tableOutput("predictionsTable"),
                                br(), br(),
                                hr(),
                                h4("Signature Visualization"),
                                verbatimTextOutput("signatureVisualizationInfo"),
                                br(),
                                h5("Coefficient Lollipop Plot"),
                                p("Shows the coefficients (weights) for each gene in the signature panel."),
                                plotOutput("coefLollipopPlot", height = "600px"),
                                downloadButton("downloadCoefLollipop", "Download Plot"),
                                br(), br(),
                                conditionalPanel(
                                    condition = "output.stabilityAvailable == true",
                                    h5("Stability Bars Plot"),
                                    p("Shows the selection frequency of genes across cross-validation folds."),
                                    plotOutput("stabilityBarsPlot", height = "600px"),
                                    downloadButton("downloadStabilityBars", "Download Plot")
                                ),
                                conditionalPanel(
                                    condition = "output.stabilityAvailable == false",
                                    h5("Stability Bars Plot"),
                                    p("Stability information is not available for the pre-trained signature.")
                                )
                            )
                        )
                    )
                )
            )
        ),

        # Mode 2: Signature Validation
        tabPanel(
            "Mode 2: Signature Validation",
            sidebarLayout(
                sidebarPanel(
                    h3("Mode 2: Signature Validation"),
                    p("Train a new signature using nested cross-validation and compare to pre-trained signature."),
                    br(),

                    # Demo data section
                    h4("Demo Data"),
                    p("Load or download labeled demo dataset (gse201926_trainmini):"),
                    actionButton("loadDemoMode2", "Load Demo Data", class = "btn-primary"),
                    br(), br(),
                    downloadButton("downloadDemoMode2", "Download Demo Data", class = "btn-secondary"),
                    br(), br(),
                    hr(),

                    # Data upload section
                    h4("Upload Your Data"),
                    p("Upload labeled endometrial RNA-seq data (counts + phenotype with PS/PIS labels):"),
                    fileInput("countsFileMode2", "Counts Matrix (TSV/CSV)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "Genes × Samples"
                    ),
                    fileInput("phenoFileMode2", "Phenotype Metadata (TSV/CSV, required)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "sample_id, group columns"
                    ),
                    fileInput("annotFileMode2", "Annotation Data (TSV/CSV, optional)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "Gene annotations"
                    ),
                    br(),

                    # Data status
                    verbatimTextOutput("dataStatusMode2"),
                    br(),

                    # Training parameters section
                    h4("Training Parameters"),
                    p("Adjust training parameters (defaults are recommended):"),
                    br(),

                    # Transform method
                    radioButtons("transformMode2", "Transform Method",
                        choices = list(
                            "log1p-CPM" = "log1p-cpm",
                            "VST" = "vst"
                        ),
                        selected = "log1p-cpm"
                    ),

                    # CPM thresholds
                    numericInput("cpmMinMode2", "CPM Minimum Threshold",
                        value = 1, min = 0, max = 10, step = 0.1
                    ),
                    numericInput("cpmMinSamplesMode2", "CPM Min Samples",
                        value = 4, min = 1, max = 20, step = 1
                    ),

                    # Top-K genes
                    numericInput("topKMode2", "Top-K Genes",
                        value = 300, min = 50, max = 1000, step = 50
                    ),

                    # CV settings
                    radioButtons("outerMode2", "Outer CV Method",
                        choices = list(
                            "K-fold" = "kfold",
                            "Leave-Pair-Out" = "lpo"
                        ),
                        selected = "kfold"
                    ),
                    numericInput("outerFoldsMode2", "Outer Folds (NULL = auto)",
                        value = NULL, min = 2, max = 10, step = 1
                    ),
                    numericInput("innerFoldsMode2", "Inner Folds",
                        value = 5, min = 3, max = 10, step = 1
                    ),
                    numericInput("innerRepeatsMode2", "Inner Repeats",
                        value = 10, min = 1, max = 20, step = 1
                    ),

                    # Lambda rule
                    radioButtons("lambdaRuleMode2", "Lambda Selection Rule",
                        choices = list(
                            "1-SE (recommended)" = "1se",
                            "Minimum" = "min"
                        ),
                        selected = "1se"
                    ),

                    # Calibration method
                    radioButtons("calibrationMode2", "Calibration Method",
                        choices = list(
                            "Platt Scaling" = "platt",
                            "Isotonic Regression" = "isotonic",
                            "None" = "none"
                        ),
                        selected = "platt"
                    ),

                    # Stability selection
                    checkboxInput("stabilityMode2", "Stability Selection (auto-enabled for n >= 20)", value = FALSE),
                    numericInput("stabilityResamplesMode2", "Stability Resamples",
                        value = 100, min = 50, max = 500, step = 50
                    ),
                    p(
                        style = "font-size: 12px; color: #666; font-style: italic;",
                        "Note: Stability Selection is required to generate the stability bars plot."
                    ),

                    # Seed
                    numericInput("seedMode2", "Random Seed",
                        value = 123, min = 1, max = 9999, step = 1
                    ),
                    br(), hr(), br(),

                    # Training action button
                    h4("Training"),
                    p("Click to train signature using nested cross-validation:"),
                    actionButton("runTraining", "Train Signature", class = "btn-success"),
                    br(), br(),

                    # Status messages
                    verbatimTextOutput("statusMessageMode2")
                ),
                mainPanel(
                    tabsetPanel(
                        id = "mode2Tabs",

                        # Overview Tab
                        tabPanel(
                            "Overview",
                            h3("Training Summary"),
                            br(),
                            conditionalPanel(
                                condition = "output.trainingRun == false",
                                p("Click 'Train Signature' in the sidebar to run training.")
                            ),
                            conditionalPanel(
                                condition = "output.trainingRun == true",
                                h4("Training Summary"),
                                verbatimTextOutput("trainingSummary"),
                                br(),
                                h4("Performance Metrics"),
                                verbatimTextOutput("performanceMetrics"),
                                br(),
                                h4("Signature Panel Information"),
                                verbatimTextOutput("signaturePanelInfo")
                            )
                        ),

                        # Performance Tab
                        tabPanel(
                            "Performance",
                            h3("Performance Plots"),
                            br(),
                            conditionalPanel(
                                condition = "output.trainingRun == false",
                                p("Please run training first to see performance plots.")
                            ),
                            conditionalPanel(
                                condition = "output.trainingRun == true",
                                h4("ROC Curve"),
                                radioButtons("rocCalibratedMode2", "Use Calibrated Predictions",
                                    choices = list("Yes" = "TRUE", "No" = "FALSE"),
                                    selected = "TRUE",
                                    inline = TRUE
                                ),
                                plotOutput("rocPlotMode2", height = "500px"),
                                downloadButton("downloadROCMode2", "Download Plot"),
                                br(), br(), hr(), br(),
                                h4("Precision-Recall Curve"),
                                radioButtons("prCalibratedMode2", "Use Calibrated Predictions",
                                    choices = list("Yes" = "TRUE", "No" = "FALSE"),
                                    selected = "TRUE",
                                    inline = TRUE
                                ),
                                plotOutput("prPlotMode2", height = "500px"),
                                downloadButton("downloadPRMode2", "Download Plot"),
                                br(), br(), hr(), br(),
                                h4("Calibration Plot"),
                                radioButtons("calCalibratedMode2", "Use Calibrated Predictions",
                                    choices = list("Yes" = "TRUE", "No" = "FALSE"),
                                    selected = "TRUE",
                                    inline = TRUE
                                ),
                                plotOutput("calibrationPlotMode2", height = "500px"),
                                downloadButton("downloadCalibrationMode2", "Download Plot")
                            )
                        ),

                        # Comparison Tab
                        tabPanel(
                            "Comparison",
                            h3("Pre-trained vs New Signature Comparison"),
                            br(),
                            conditionalPanel(
                                condition = "output.trainingRun == false",
                                p("Please run training first to see comparison plots.")
                            ),
                            conditionalPanel(
                                condition = "output.trainingRun == true && output.pretrainedAvailable == false",
                                p("Pre-trained signature not available. Comparison plots require both pre-trained and new signatures.")
                            ),
                            conditionalPanel(
                                condition = "output.trainingRun == true && output.pretrainedAvailable == true",
                                h4("Comparison Plots"),
                                p("Side-by-side comparison of pre-trained and newly trained signatures:"),
                                plotOutput("comparisonPlotMode2", height = "800px"),
                                downloadButton("downloadComparisonMode2", "Download Plot"),
                                br(), br(),
                                h4("Comparison Metrics"),
                                verbatimTextOutput("comparisonMetrics")
                            )
                        ),

                        # Signature Panel Tab
                        tabPanel(
                            "Signature Panel",
                            h3("Signature Visualization"),
                            br(),
                            conditionalPanel(
                                condition = "output.trainingRun == false",
                                p("Please run training first to see signature plots.")
                            ),
                            conditionalPanel(
                                condition = "output.trainingRun == true",
                                h4("Coefficient Lollipop Plot"),
                                p("Shows the coefficients (weights) for each gene in the signature panel."),
                                plotOutput("coefLollipopPlotMode2", height = "600px"),
                                downloadButton("downloadCoefLollipopMode2", "Download Plot"),
                                br(), br(), hr(), br(),
                                conditionalPanel(
                                    condition = "output.stabilityAvailableMode2 == true",
                                    h4("Stability Bars Plot"),
                                    p("Shows the selection frequency of genes across cross-validation folds."),
                                    plotOutput("stabilityBarsPlotMode2", height = "600px"),
                                    downloadButton("downloadStabilityBarsMode2", "Download Plot")
                                ),
                                conditionalPanel(
                                    condition = "output.stabilityAvailableMode2 == false",
                                    h4("Stability Bars Plot"),
                                    p("Stability information is not available for this signature.")
                                )
                            )
                        ),

                        # Export Tab
                        tabPanel(
                            "Export",
                            h3("Export Results"),
                            br(),
                            conditionalPanel(
                                condition = "output.trainingRun == false",
                                p("Please run training first to export results.")
                            ),
                            conditionalPanel(
                                condition = "output.trainingRun == true",
                                h4("Signature Artifacts"),
                                p("Export signature in standard formats:"),
                                div(
                                    style = "margin-bottom: 15px;",
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                                        downloadButton("downloadSignatureCSVMode2", "Download Signature (CSV)"),
                                        tags$details(
                                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                                            div(
                                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                                h5("What's in this export:"),
                                                p("The CSV file contains the gene signature panel with:"),
                                                tags$ul(
                                                    tags$li(strong("gene_id:"), " Gene identifier"),
                                                    tags$li(strong("coefficient:"), " LASSO coefficient (weight) for each gene"),
                                                    tags$li(strong("intercept:"), " Model intercept (in separate row or metadata)"),
                                                    tags$li(strong("metadata:"), " Training parameters and preprocessing information")
                                                ),
                                                br(),
                                                h5("What you can do with this export:"),
                                                tags$ul(
                                                    tags$li("Use the signature in other tools or pipelines"),
                                                    tags$li("Share signature with collaborators"),
                                                    tags$li("Version control your trained signatures"),
                                                    tags$li("Compare signatures across different training runs"),
                                                    tags$li("Load signature in Mode 1 for classification")
                                                ),
                                                br(),
                                                h5("Format specifications:"),
                                                p("CSV format (comma-separated values), UTF-8 encoding. Contains gene panel with coefficients and metadata.")
                                            )
                                        )
                                    )
                                ),
                                br(),
                                div(
                                    style = "margin-bottom: 15px;",
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                                        downloadButton("downloadSignatureJSONMode2", "Download Signature (JSON)"),
                                        tags$details(
                                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                                            div(
                                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                                h5("What's in this export:"),
                                                p("The JSON file contains the complete signature recipe with:"),
                                                tags$ul(
                                                    tags$li(strong("signature:"), " Gene panel with coefficients"),
                                                    tags$li(strong("recipe:"), " Complete preprocessing parameters (transform method, CPM thresholds, Top-K selection)"),
                                                    tags$li(strong("metadata:"), " Training information (date, parameters, performance metrics)")
                                                ),
                                                br(),
                                                h5("What you can do with this export:"),
                                                tags$ul(
                                                    tags$li("Reproduce exact preprocessing steps in other pipelines"),
                                                    tags$li("Share complete signature recipe with collaborators"),
                                                    tags$li("Version control signature with full metadata"),
                                                    tags$li("Load signature in Mode 1 with complete preprocessing information")
                                                ),
                                                br(),
                                                h5("Format specifications:"),
                                                p("JSON format (JavaScript Object Notation), UTF-8 encoding. Contains structured data with signature and recipe information.")
                                            )
                                        )
                                    )
                                ),
                                br(),
                                div(
                                    style = "margin-bottom: 15px;",
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                                        downloadButton("downloadModelCardMode2", "Download Model Card (MD)"),
                                        tags$details(
                                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                                            div(
                                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                                h5("What's in this export:"),
                                                p("The Markdown file contains comprehensive model documentation with:"),
                                                tags$ul(
                                                    tags$li(strong("Model information:"), " Signature details, training date, data source"),
                                                    tags$li(strong("Performance metrics:"), " AUC, accuracy, calibration metrics"),
                                                    tags$li(strong("Training parameters:"), " All hyperparameters and settings used"),
                                                    tags$li(strong("Statistical assumptions:"), " Assumptions and limitations"),
                                                    tags$li(strong("Reproducibility information:"), " Seeds, fold specifications, session info")
                                                ),
                                                br(),
                                                h5("What you can do with this export:"),
                                                tags$ul(
                                                    tags$li("Share model information with collaborators or reviewers"),
                                                    tags$li("Document methodology for publications"),
                                                    tags$li("Track model versions and training configurations"),
                                                    tags$li("Ensure reproducibility and transparency")
                                                ),
                                                br(),
                                                h5("Format specifications:"),
                                                p("Markdown format (.md), UTF-8 encoding. Can be viewed in any text editor or rendered as HTML. Contains formatted documentation sections.")
                                            )
                                        )
                                    )
                                ),
                                br(), hr(), br(),
                                h4("Stability Data"),
                                p("Export bootstrap stability frequencies (if available):"),
                                div(
                                    style = "margin-bottom: 15px;",
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                                        downloadButton("downloadStabilityCSVMode2", "Download Stability CSV"),
                                        tags$details(
                                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                                            div(
                                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                                h5("What's in this export:"),
                                                p("The CSV file contains gene selection frequencies across resamples with:"),
                                                tags$ul(
                                                    tags$li(strong("gene_id:"), " Gene identifier"),
                                                    tags$li(strong("selection_frequency:"), " Frequency (0-1) of gene selection across bootstrap resamples"),
                                                    tags$li(strong("in_signature:"), " Whether gene is in final signature panel")
                                                ),
                                                br(),
                                                h5("What you can do with this export:"),
                                                tags$ul(
                                                    tags$li("Analyze signature stability and identify core genes"),
                                                    tags$li("Compare gene selection frequencies across different training runs"),
                                                    tags$li("Identify highly stable genes for further investigation"),
                                                    tags$li("Assess signature robustness and reproducibility")
                                                ),
                                                br(),
                                                h5("Format specifications:"),
                                                p("CSV format (comma-separated values), UTF-8 encoding. Contains gene stability metrics from bootstrap resampling.")
                                            )
                                        )
                                    )
                                ),
                                br(), hr(), br(),
                                h4("Predictions"),
                                p("Export predictions from cross-validation:"),
                                div(
                                    style = "margin-bottom: 15px;",
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                                        downloadButton("downloadPredictionsMode2", "Download Predictions (CSV)"),
                                        tags$details(
                                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                                            div(
                                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                                h5("What's in this export:"),
                                                p("The CSV file contains cross-validation predictions with:"),
                                                tags$ul(
                                                    tags$li(strong("sample_id:"), " Sample identifier"),
                                                    tags$li(strong("true_label:"), " True class label (PS or PIS)"),
                                                    tags$li(strong("predicted_class:"), " Predicted class from outer CV fold"),
                                                    tags$li(strong("probability:"), " Predicted probability"),
                                                    tags$li(strong("fold:"), " Outer CV fold identifier")
                                                ),
                                                br(),
                                                h5("What you can do with this export:"),
                                                tags$ul(
                                                    tags$li("Create custom performance visualizations"),
                                                    tags$li("Analyze per-sample prediction confidence"),
                                                    tags$li("Identify misclassified samples for further investigation"),
                                                    tags$li("Compare predictions across different models or parameters")
                                                ),
                                                br(),
                                                h5("Format specifications:"),
                                                p("CSV format (comma-separated values), UTF-8 encoding. Contains out-of-fold predictions from nested cross-validation.")
                                            )
                                        )
                                    )
                                ),
                                br(), hr(), br(),
                                h4("Full Results"),
                                p("Export complete training result (RDS format):"),
                                div(
                                    style = "margin-bottom: 15px;",
                                    div(
                                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                                        downloadButton("downloadFullResultsMode2", "Download Full Results (RDS)"),
                                        tags$details(
                                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                                            div(
                                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                                h5("What's in this export:"),
                                                p("The RDS file contains the complete R object with all training results:"),
                                                tags$ul(
                                                    tags$li(strong("signature:"), " Trained signature object"),
                                                    tags$li(strong("performance:"), " Performance metrics and plots"),
                                                    tags$li(strong("predictions:"), " Cross-validation predictions"),
                                                    tags$li(strong("stability:"), " Stability selection results"),
                                                    tags$li(strong("metadata:"), " Training parameters, seeds, fold specifications"),
                                                    tags$li(strong("session_info:"), " R session information for reproducibility")
                                                ),
                                                br(),
                                                h5("What you can do with this export:"),
                                                tags$ul(
                                                    tags$li("Load complete results in R for further analysis"),
                                                    tags$li("Reproduce all plots and metrics"),
                                                    tags$li("Extract specific components for custom analysis"),
                                                    tags$li("Share complete analysis state with collaborators"),
                                                    tags$li("Archive training results for reproducibility")
                                                ),
                                                br(),
                                                h5("Format specifications:"),
                                                p("RDS format (R Data Serialization), binary format. Can only be loaded in R using readRDS(). Contains complete R object with all training artifacts.")
                                            )
                                        )
                                    )
                                )
                            )
                        )
                    )
                )
            )
        ),

        # Mode 3: Visualization & Analysis
        tabPanel(
            "Mode 3: Visualization & Analysis",
            sidebarLayout(
                sidebarPanel(
                    h3("Mode 3: Visualization & Analysis"),
                    p("Upload endometrial RNA-seq data or use demo data to perform quality control,
        exploratory analysis, and differential expression visualization."),
                    br(),

                    # Demo data section
                    h4("Demo Data"),
                    p("Load or download the bundled demo dataset (gse201926_sample):"),
                    actionButton("loadDemoMode3", "Load Demo Data", class = "btn-primary"),
                    br(), br(),
                    downloadButton("downloadDemoMode3", "Download Demo Data", class = "btn-secondary"),
                    br(), br(),
                    hr(),

                    # Data upload section
                    h4("Upload Your Data"),
                    p("Upload your own endometrial RNA-seq data files:"),
                    fileInput("countsFileMode3", "Counts Matrix (TSV/CSV)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "Genes × Samples"
                    ),
                    fileInput("phenoFileMode3", "Phenotype Metadata (TSV/CSV)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "sample_id, group columns"
                    ),
                    fileInput("annotFileMode3", "Annotation Data (TSV/CSV, optional)",
                        accept = c(".tsv", ".csv", ".txt"),
                        placeholder = "Gene annotations"
                    ),
                    br(),

                    # Data status
                    verbatimTextOutput("dataStatusMode3"),
                    br(),

                    # Status messages
                    verbatimTextOutput("statusMessageMode3")
                ),
                mainPanel(
                    uiOutput("mode3TabsUI")
                )
            )
        ),

        # Methodology Information Page
        tabPanel(
            "Methodology",
            div(
                style = "max-width: 1000px; margin: 0 auto; padding: 20px;",
                h2("Methodology: Small Sample Signature Extraction", style = "text-align: center; color: #2c3e50; margin-bottom: 30px;"),
                h3("Introduction", style = "color: #34495e; margin-top: 40px;"),
                p("EndoSignatureR uses a carefully designed approach for extracting gene signatures from small endometrial RNA-seq datasets. This methodology addresses the unique challenges of working with small sample sizes (small n) and high-dimensional data (large p, where p >> n)."),
                br(),
                h3("How We Approached Small Sample Signature Extraction", style = "color: #34495e; margin-top: 40px;"),
                p("The challenge of small sample signature extraction involves:"),
                tags$ul(
                    tags$li("High-dimensional data: Thousands of genes (p) but few samples (n), where p >> n"),
                    tags$li("Risk of overfitting: Standard methods can overfit to small training sets"),
                    tags$li("Need for generalization: Signatures must work on new, unseen data"),
                    tags$li("Biological relevance: Selected genes should be biologically meaningful")
                ),
                p("Our approach combines several best practices to address these challenges:"),
                tags$ul(
                    tags$li("Nested cross-validation to estimate true generalization performance"),
                    tags$li("LASSO regularization to select sparse, interpretable gene sets"),
                    tags$li("In-fold preprocessing to prevent data leakage"),
                    tags$li("Differential expression screening to focus on biologically relevant genes"),
                    tags$li("Stability selection to identify robust gene signatures"),
                    tags$li("Probability calibration to improve prediction reliability")
                ),
                br(),
                h3("Key Methodological Choices", style = "color: #34495e; margin-top: 40px;"),
                h4("1. Nested Cross-Validation", style = "color: #7f8c8d; margin-top: 20px;"),
                p("Nested cross-validation (CV) uses two levels of data splitting:"),
                tags$ul(
                    tags$li(strong("Outer CV:"), " Evaluates generalization performance. Each fold serves as a test set, with the remaining data used for training."),
                    tags$li(strong("Inner CV:"), " Tunes hyperparameters (e.g., LASSO λ) using only training data from the outer fold."),
                    tags$li(strong("Why it matters:"), " Provides unbiased estimates of model performance on new data, crucial for small datasets where a single train/test split would be unreliable.")
                ),
                br(),
                h4("2. LASSO Regularization", style = "color: #7f8c8d; margin-top: 20px;"),
                p("LASSO (Least Absolute Shrinkage and Selection Operator) logistic regression:"),
                tags$ul(
                    tags$li(strong("Sparsity:"), " L1 penalty drives most coefficients to zero, automatically selecting a small subset of genes."),
                    tags$li(strong("Interpretability:"), " Sparse models are easier to interpret and validate biologically."),
                    tags$li(strong("Regularization:"), " Prevents overfitting by penalizing large coefficients."),
                    tags$li(strong("Why it matters:"), " In p >> n scenarios, LASSO helps identify the most important genes while avoiding overfitting.")
                ),
                br(),
                h4("3. In-Fold Preprocessing", style = "color: #7f8c8d; margin-top: 20px;"),
                p("All preprocessing operations occur inside cross-validation folds:"),
                tags$ul(
                    tags$li(strong("Transformation:"), " Normalization parameters (e.g., library size factors) computed from training data only."),
                    tags$li(strong("Filtering:"), " CPM thresholds and gene filtering based on training data only."),
                    tags$li(strong("Gene selection:"), " Top-K or DE-based selection performed separately in each fold."),
                    tags$li(strong("Why it matters:"), " Prevents data leakage - test data never influences preprocessing decisions, ensuring valid performance estimates.")
                ),
                br(),
                h4("4. Stability Selection", style = "color: #7f8c8d; margin-top: 20px;"),
                p("Stability selection assesses how consistently genes are selected across resamples:"),
                tags$ul(
                    tags$li(strong("Bootstrap resampling:"), " Multiple resamples of the training data."),
                    tags$li(strong("Selection frequency:"), " Frequency (0-1) that each gene is selected across resamples."),
                    tags$li(strong("Robust genes:"), " Genes selected frequently are more stable and likely to generalize."),
                    tags$li(strong("Why it matters:"), " Identifies core signature genes that are robust to data variation, important for small datasets where results can be sensitive to specific samples.")
                ),
                br(),
                h4("5. Differential Expression Screening", style = "color: #7f8c8d; margin-top: 20px;"),
                p("Top-K gene selection based on differential expression statistics:"),
                tags$ul(
                    tags$li(strong("DE statistics:"), " Genes ranked by fold change, p-values, or FDR."),
                    tags$li(strong("Top-K selection:"), " Select top K genes (e.g., 300) before LASSO training."),
                    tags$li(strong("Biological relevance:"), " Focuses on genes that show meaningful differences between classes."),
                    tags$li(strong("Why it matters:"), " In p >> n scenarios, DE screening improves stability by focusing on biologically relevant features and reducing noise.")
                ),
                br(),
                h3("Rationale for Small-n Considerations", style = "color: #34495e; margin-top: 40px;"),
                h4("Leave-Pair-Out Cross-Validation", style = "color: #7f8c8d; margin-top: 20px;"),
                p("For very small, balanced datasets, we use leave-pair-out (LPO) CV:"),
                tags$ul(
                    tags$li(strong("Method:"), " Each test set contains one sample from each class (1 PS + 1 PIS)."),
                    tags$li(strong("Maximizes training data:"), " Leaves maximum samples for training while maintaining valid test sets."),
                    tags$li(strong("Class balance:"), " Ensures test sets are balanced, important for small datasets."),
                    tags$li(strong("When used:"), " For tiny balanced cohorts where standard k-fold CV would leave too few training samples.")
                ),
                br(),
                h4("Anti-Leakage Guards", style = "color: #7f8c8d; margin-top: 20px;"),
                p("Strict separation between training and test data:"),
                tags$ul(
                    tags$li(strong("No test data in preprocessing:"), " All transformations, filters, and selections use only training data."),
                    tags$li(strong("Fold-specific preprocessing:"), " Each CV fold has its own preprocessing parameters."),
                    tags$li(strong("Validation:"), " Test data is only used for final evaluation, never for model development."),
                    tags$li(strong("Why it matters:"), " Prevents optimistic bias in performance estimates, critical for small datasets where even small leaks can dramatically inflate performance.")
                ),
                br(),
                h4("Probability Calibration", style = "color: #7f8c8d; margin-top: 20px;"),
                p("Calibration ensures predicted probabilities reflect true class frequencies:"),
                tags$ul(
                    tags$li(strong("Platt scaling or isotonic regression:"), " Adjusts raw model probabilities."),
                    tags$li(strong("Reliability:"), " Calibrated probabilities are more reliable for decision-making."),
                    tags$li(strong("Why it matters:"), " Small datasets can produce poorly calibrated probabilities; calibration improves reliability for clinical or research applications.")
                ),
                br(),
                h3("Statistical Assumptions and Limitations", style = "color: #34495e; margin-top: 40px;"),
                h4("Data Type Assumptions", style = "color: #7f8c8d; margin-top: 20px;"),
                tags$ul(
                    tags$li("Bulk RNA-seq gene counts from endometrial tissue (not single-cell, not microarray, not other tissues)"),
                    tags$li("The pre-trained signature is specific to endometrial tissue and is not portable across tissues or divergent processing pipelines")
                ),
                br(),
                h4("Design Assumptions", style = "color: #7f8c8d; margin-top: 20px;"),
                tags$ul(
                    tags$li("Binary labels (PS vs PIS) with independent biological replicates"),
                    tags$li("The package supports only binary classification; multiclass or continuous outcomes are out of scope for v1")
                ),
                br(),
                h4("Modeling Assumptions", style = "color: #7f8c8d; margin-top: 20px;"),
                tags$ul(
                    tags$li("Linear log-odds relationship (logistic regression) with L1 penalty (LASSO)"),
                    tags$li("The L1 penalty drives most coefficients to zero, yielding an interpretable gene signature")
                ),
                br(),
                h4("Resampling Assumptions", style = "color: #7f8c8d; margin-top: 20px;"),
                tags$ul(
                    tags$li("Nested cross-validation is used as a proxy for generalization in small-n datasets"),
                    tags$li("Leave-pair-out (LPO) cross-validation is used when n is tiny and class-balanced to maximize training data while maintaining valid test sets")
                ),
                br(),
                h4("Limitations", style = "color: #7f8c8d; margin-top: 20px;"),
                tags$ul(
                    tags$li(strong("Tissue specificity:"), " Trained signatures are not portable across tissues or divergent processing pipelines"),
                    tags$li(strong("Small-n uncertainty:"), " Results should be interpreted with caution; report confidence intervals, permutation p-values, and stability metrics"),
                    tags$li(strong("Binary only:"), " Multiclass or continuous outcomes are out of scope for v1"),
                    tags$li(strong("Batch effects:"), " Users should provide batch fields and use appropriate designs (~ batch + group) or pre-correct for batch effects")
                ),
                br(),
                h3("References and Further Reading", style = "color: #34495e; margin-top: 40px;"),
                p("For more detailed information, please refer to:"),
                tags$ul(
                    tags$li("Package vignettes: ", code("browseVignettes('endoSignatureR')"), " for workflow-specific tutorials"),
                    tags$li("README: ", code("?endoSignatureR"), " for package overview and installation"),
                    tags$li("Function documentation: Individual function help pages (e.g., ", code("?esr_trainEndometrialSignature"), ") for detailed parameter descriptions")
                ),
                p("Key methodological references:"),
                tags$ul(
                    tags$li("LASSO: Friedman, Hastie & Tibshirani (2010). Regularization paths for generalized linear models via coordinate descent."),
                    tags$li("Nested CV: Varma & Simon (2006). Bias in error estimation when using cross-validation for model selection."),
                    tags$li("Stability Selection: Meinshausen & Bühlmann (2010). Stability selection."),
                    tags$li("Differential Expression: Ritchie et al. (2015). limma powers differential expression analyses.")
                ),
                br(),
                hr(),
                p(em("This methodology is designed specifically for small endometrial RNA-seq datasets. Results should be interpreted in the context of the statistical assumptions and limitations outlined above."), style = "text-align: center; color: #7f8c8d; margin-top: 30px;")
            )
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    # Increase maximum upload size to 500MB (default is 5MB)
    options(shiny.maxRequestSize = 500 * 1024^2)

    # Navigation observers for landing page buttons
    observeEvent(input$navToMode1, {
        updateTabsetPanel(session, "modeTabs", selected = "Mode 1: Rapid Classification")
    })

    observeEvent(input$navToMode2, {
        updateTabsetPanel(session, "modeTabs", selected = "Mode 2: Signature Validation")
    })

    observeEvent(input$navToMode3, {
        updateTabsetPanel(session, "modeTabs", selected = "Mode 3: Visualization & Analysis")
    })

    # Reactive values to store data
    values <- reactiveValues(
        # Mode 1 data
        counts_mode1 = NULL,
        annot_mode1 = NULL,
        predictions = NULL,
        signature = NULL,
        selectedSignature = NULL, # Currently selected signature (pre-trained or user-trained)
        mode1_completed = FALSE,
        classificationRun = FALSE,

        # Mode 2 data
        counts_mode2 = NULL,
        pheno_mode2 = NULL,
        annot_mode2 = NULL,
        trainingResult = NULL,
        mode2_completed = FALSE,
        trainingRun = FALSE,

        # Mode 3 data
        counts = NULL,
        pheno = NULL,
        annot = NULL,
        counts_t = NULL,
        qc_metrics = NULL,
        de_table = NULL,
        selected_genes = NULL,
        bundle = NULL,
        dataLoaded = FALSE,
        deRun = FALSE,
        # Plot generation flags
        maPlotGenerated = FALSE,
        volcanoPlotGenerated = FALSE,
        heatmapPlotGenerated = FALSE,
        # Plot parameters
        maPlotParams = list(fdr_threshold = 0.05, log2fc_threshold = 1, point_size = 1.5, point_alpha = 0.6, legend_position = "right", theme = "bw"),
        volcanoPlotParams = list(fdr_threshold = 0.05, log2fc_threshold = 1, point_size = 1.5, point_alpha = 0.6, legend_position = "right", theme = "bw"),
        heatmapPlotParams = list(n_genes = 50, scale = "row", show_row_names = FALSE),
        # QC/EDA plot parameters
        libsizePlotParams = list(point_size = 2, point_alpha = 0.6, bins = 30, legend_position = "right", theme = "bw"),
        zerosPlotParams = list(bins = 30, theme = "bw"),
        pcaPlotParams = list(point_size = 4, point_alpha = 0.8, legend_position = "right", theme = "bw")
    )

    # Load pre-trained signature on app startup
    observe({
        tryCatch(
            {
                values$signature <- endoSignatureR::esr_loadPretrainedSignature()
                # Initialize selected signature to pre-trained (will be overridden if user-trained available)
                if (is.null(values$selectedSignature)) {
                    values$selectedSignature <- values$signature
                }
            },
            error = function(e) {
                # Signature loading will fail gracefully if files not found
                # User will see error message in signature info output
            }
        )
    })

    # Load signature from uploaded files
    observeEvent(
        {
            input$signatureCSVFile
            input$signatureJSONFile
            input$signatureStabilityFile
        },
        {
            # Require CSV and JSON files
            if (is.null(input$signatureCSVFile) || is.null(input$signatureJSONFile)) {
                return()
            }

            tryCatch(
                {
                    # Load CSV signature (similar to esr_loadPretrainedSignature logic)
                    csv_data <- readr::read_csv(
                        input$signatureCSVFile$datapath,
                        col_types = readr::cols(
                            gene_id = readr::col_character(),
                            coefficient = readr::col_double(),
                            selection_frequency = readr::col_integer(),
                            bootstrap_frequency = readr::col_double()
                        ),
                        show_col_types = FALSE
                    )

                    # Separate intercept from panel genes
                    intercept_row <- csv_data[csv_data$gene_id == "(Intercept)", , drop = FALSE]
                    panel_rows <- csv_data[csv_data$gene_id != "(Intercept)", , drop = FALSE]

                    # Extract intercept value
                    intercept_value <- if (nrow(intercept_row) > 0) {
                        as.numeric(intercept_row$coefficient[1])
                    } else {
                        stop("Intercept not found in signature CSV")
                    }

                    # Extract panel genes and coefficients
                    panel <- panel_rows$gene_id
                    coefficients <- panel_rows$coefficient
                    names(coefficients) <- panel

                    # Extract selection frequencies
                    selection_frequency <- panel_rows$selection_frequency
                    names(selection_frequency) <- panel

                    # Load JSON recipe
                    recipe_data <- jsonlite::read_json(input$signatureJSONFile$datapath, simplifyVector = TRUE)

                    # Validate recipe structure
                    required_sections <- c("preprocessing", "training", "signature", "reproducibility")
                    missing_sections <- setdiff(required_sections, names(recipe_data))
                    if (length(missing_sections) > 0) {
                        stop(
                            "JSON recipe missing required sections: ",
                            paste(missing_sections, collapse = ", ")
                        )
                    }

                    # Load optional stability CSV
                    stability <- NULL
                    if (!is.null(input$signatureStabilityFile)) {
                        stability_data <- readr::read_csv(
                            input$signatureStabilityFile$datapath,
                            col_types = readr::cols(
                                gene_id = readr::col_character(),
                                bootstrap_frequency = readr::col_double()
                            ),
                            show_col_types = FALSE
                        )

                        # Extract bootstrap frequencies
                        bootstrap_frequency <- stability_data$bootstrap_frequency
                        names(bootstrap_frequency) <- stability_data$gene_id

                        # Match bootstrap frequencies to panel
                        bootstrap_freq_matched <- bootstrap_frequency[panel]
                        names(bootstrap_freq_matched) <- panel
                        bootstrap_freq_matched <- bootstrap_freq_matched[!is.na(bootstrap_freq_matched)]

                        if (length(bootstrap_freq_matched) > 0) {
                            stability <- list(bootstrap_frequency = bootstrap_freq_matched)
                        }
                    } else {
                        # Check if bootstrap frequencies are in CSV
                        if ("bootstrap_frequency" %in% names(panel_rows)) {
                            bootstrap_freq_vec <- panel_rows$bootstrap_frequency
                            names(bootstrap_freq_vec) <- panel
                            bootstrap_freq_vec <- bootstrap_freq_vec[!is.na(bootstrap_freq_vec)]
                            if (length(bootstrap_freq_vec) > 0) {
                                stability <- list(bootstrap_frequency = bootstrap_freq_vec)
                            }
                        }
                    }

                    # Validate signature structure
                    if (length(panel) == 0) {
                        stop("Signature panel is empty")
                    }

                    if (length(coefficients) != length(panel)) {
                        stop("Coefficient length does not match panel length")
                    }

                    if (is.na(intercept_value)) {
                        stop("Intercept value is NA")
                    }

                    # Construct signature list
                    values$uploadedSignature <- list(
                        panel = panel,
                        coefficients = coefficients,
                        intercept = intercept_value,
                        selection_frequency = selection_frequency,
                        recipe = recipe_data,
                        stability = stability
                    )

                    # Update selected signature if "uploaded" is selected
                    if (!is.null(input$signatureSelection) && input$signatureSelection == "uploaded") {
                        values$selectedSignature <- values$uploadedSignature
                    } else if (is.null(values$selectedSignature)) {
                        # If no selection yet, default to uploaded if it's the only option
                        values$selectedSignature <- values$uploadedSignature
                    }

                    output$statusMessageMode1 <- renderText("Signature files loaded successfully!")
                },
                error = function(e) {
                    output$statusMessageMode1 <- renderText(paste("Error loading signature files:", e$message))
                }
            )
        }
    )

    # Update selected signature when user changes selection
    observeEvent(input$signatureSelection, {
        if (input$signatureSelection == "pretrained") {
            values$selectedSignature <- values$signature
            # Clear predictions if signature changes after classification
            if (values$classificationRun) {
                values$predictions <- NULL
                values$classificationRun <- FALSE
            }
        } else if (input$signatureSelection == "user-trained") {
            if (!is.null(values$trainingResult) && !is.null(values$trainingResult$signature)) {
                values$selectedSignature <- values$trainingResult$signature
                # Clear predictions if signature changes after classification
                if (values$classificationRun) {
                    values$predictions <- NULL
                    values$classificationRun <- FALSE
                }
            } else {
                # User-trained signature not available, revert to pre-trained
                updateRadioButtons(session, "signatureSelection", selected = "pretrained")
                values$selectedSignature <- values$signature
            }
        } else if (input$signatureSelection == "uploaded") {
            if (!is.null(values$uploadedSignature)) {
                values$selectedSignature <- values$uploadedSignature
                # Clear predictions if signature changes after classification
                if (values$classificationRun) {
                    values$predictions <- NULL
                    values$classificationRun <- FALSE
                }
            } else {
                # Uploaded signature not available, revert to pre-trained
                updateRadioButtons(session, "signatureSelection", selected = "pretrained")
                values$selectedSignature <- values$signature
            }
        }
    })

    # Also update selected signature when Mode 2 training completes
    observeEvent(values$mode2_completed, {
        if (values$mode2_completed && !is.null(values$trainingResult)) {
            # Default to user-trained signature if it just became available
            if (is.null(input$signatureSelection) || input$signatureSelection == "pretrained") {
                # Auto-select user-trained as default
                values$selectedSignature <- values$trainingResult$signature
            } else if (input$signatureSelection == "user-trained") {
                # Update if user has already selected user-trained
                values$selectedSignature <- values$trainingResult$signature
            }
        }
    })

    # Display signature information (for Mode 1) - shows selected signature
    output$signatureInfo <- renderText({
        if (is.null(values$selectedSignature)) {
            return("Signature not available. Please ensure package is properly installed, train a signature in Mode 2, or upload signature files.")
        }

        # Check if signatureSelection input exists and has a valid value
        sig_selection <- if (is.null(input$signatureSelection) || length(input$signatureSelection) == 0) {
            "unknown"
        } else {
            input$signatureSelection
        }

        sig <- values$selectedSignature
        sig_type <- if (sig_selection == "pretrained") {
            "Pre-trained"
        } else if (sig_selection == "user-trained") {
            "User-trained"
        } else if (sig_selection == "uploaded") {
            "Uploaded"
        } else {
            "Selected"
        }

        info_text <- paste(
            sig_type, "signature loaded successfully!\n",
            "Panel size:", length(sig$panel), "genes\n",
            "Preprocessing:", sig$recipe$preprocessing$transform %||% "log1p-cpm"
        )

        if (sig_selection == "pretrained") {
            info_text <- paste(info_text, "\nTrained on: GSE201926 dataset")
        } else if (sig_selection == "user-trained") {
            info_text <- paste(info_text, "\nTrained on: Your data (Mode 2)")
        } else if (sig_selection == "uploaded") {
            info_text <- paste(info_text, "\nSource: Uploaded signature files")
        }

        return(info_text)
    })

    # Display signature source information in Results tab
    output$signatureSourceInfo <- renderText({
        if (!values$classificationRun || is.null(values$selectedSignature)) {
            return("")
        }

        sig_type <- if (input$signatureSelection == "pretrained") {
            "Pre-trained"
        } else if (input$signatureSelection == "user-trained") {
            "User-trained"
        } else if (input$signatureSelection == "uploaded") {
            "Uploaded"
        } else {
            "Selected"
        }
        paste("Using", sig_type, "signature for classification")
    })

    # Display signature visualization info
    output$signatureVisualizationInfo <- renderText({
        if (!values$classificationRun || is.null(values$selectedSignature)) {
            return("")
        }

        sig_type <- if (input$signatureSelection == "pretrained") {
            "pre-trained"
        } else if (input$signatureSelection == "user-trained") {
            "user-trained"
        } else if (input$signatureSelection == "uploaded") {
            "uploaded"
        } else {
            "selected"
        }
        paste("Visualization of the", sig_type, "signature used for classification.")
    })

    # Display signature information (for Signature Plots tab)
    output$signatureInfoMode3 <- renderText({
        if (is.null(values$signature)) {
            return("Pre-trained signature not available.")
        }
        paste(
            "Pre-trained signature loaded successfully!\n",
            "Panel size:", length(values$signature$panel), "genes\n",
            "Preprocessing:", values$signature$recipe$preprocessing$transform %||% "log1p-cpm\n",
            "Trained on: GSE201926 dataset"
        )
    })

    # Output flags for conditional panels
    output$classificationRun <- reactive({
        values$classificationRun
    })
    outputOptions(output, "classificationRun", suspendWhenHidden = FALSE)

    output$mode1Completed <- reactive({
        values$mode1_completed
    })
    outputOptions(output, "mode1Completed", suspendWhenHidden = FALSE)

    output$stabilityAvailable <- reactive({
        if (is.null(values$signature)) {
            return(FALSE)
        }
        !is.null(values$signature$stability) || !is.null(values$signature$selection_frequency)
    })
    outputOptions(output, "stabilityAvailable", suspendWhenHidden = FALSE)

    output$deRun <- reactive({
        values$deRun
    })
    outputOptions(output, "deRun", suspendWhenHidden = FALSE)

    output$maPlotGenerated <- reactive({
        values$maPlotGenerated
    })
    outputOptions(output, "maPlotGenerated", suspendWhenHidden = FALSE)

    output$volcanoPlotGenerated <- reactive({
        values$volcanoPlotGenerated
    })
    outputOptions(output, "volcanoPlotGenerated", suspendWhenHidden = FALSE)

    output$heatmapPlotGenerated <- reactive({
        values$heatmapPlotGenerated
    })
    outputOptions(output, "heatmapPlotGenerated", suspendWhenHidden = FALSE)

    # Mode 2 output flags
    output$trainingRun <- reactive({
        values$trainingRun
    })
    outputOptions(output, "trainingRun", suspendWhenHidden = FALSE)

    output$pretrainedAvailable <- reactive({
        !is.null(values$signature)
    })
    outputOptions(output, "pretrainedAvailable", suspendWhenHidden = FALSE)

    output$stabilityAvailableMode2 <- reactive({
        if (is.null(values$trainingResult)) {
            return(FALSE)
        }
        !is.null(values$trainingResult$signature$selection_frequency) ||
            !is.null(values$trainingResult$stability)
    })
    outputOptions(output, "stabilityAvailableMode2", suspendWhenHidden = FALSE)

    # Mode 2 completion flag for conditional UI
    output$mode2Completed <- reactive({
        values$mode2_completed
    })
    outputOptions(output, "mode2Completed", suspendWhenHidden = FALSE)

    # Render signature selection UI (conditional on Mode 2 completion and uploaded signature)
    output$signatureSelectionUI <- renderUI({
        choices_list <- list("Pre-trained Signature" = "pretrained")

        # Add user-trained option if Mode 2 training is completed
        if (values$mode2_completed && !is.null(values$trainingResult)) {
            choices_list[["User-Trained Signature (from Mode 2)"]] <- "user-trained"
        }

        # Add uploaded signature option if signature files have been uploaded
        if (!is.null(values$uploadedSignature)) {
            choices_list[["Uploaded Signature"]] <- "uploaded"
        }

        # Determine default selection: prefer user-trained if available, otherwise pre-trained
        default_selection <- if (values$mode2_completed && !is.null(values$trainingResult)) {
            "user-trained" # Default to user-trained if available
        } else if (!is.null(values$uploadedSignature)) {
            "uploaded" # Default to uploaded if available
        } else {
            "pretrained" # Fall back to pre-trained
        }

        # Preserve current selection if valid, otherwise use default
        current_selection <- if (exists("input$signatureSelection") && !is.null(input$signatureSelection)) {
            if (input$signatureSelection %in% unlist(choices_list)) {
                input$signatureSelection
            } else {
                default_selection
            }
        } else {
            default_selection
        }

        tagList(
            radioButtons("signatureSelection", "Signature to Use",
                choices = choices_list,
                selected = current_selection
            ),
            if (!values$mode2_completed && is.null(values$trainingResult) && is.null(values$uploadedSignature)) {
                p(strong("Note:"), " User-trained signature option will appear after training a signature in Mode 2, or upload signature files above.",
                    style = "color: #666; font-size: 0.9em;"
                )
            }
        )
    })

    # Render Mode 3 tabs (QC/EDA, DE, Export)
    output$mode3TabsUI <- renderUI({
        tabs <- list(
            # QC/EDA Tab
            tabPanel(
                "QC/EDA",
                h3("Quality Control and Exploratory Data Analysis"),
                br(),
                h4("Library Size"),
                fluidRow(
                    column(
                        6,
                        plotOutput("libsizePlot", height = "400px"),
                        downloadButton("downloadLibsize", "Download Plot")
                    ),
                    column(
                        6,
                        h5("Styling Options"),
                        numericInput("libsizePointSize", "Point Size",
                            value = 2, min = 0.5, max = 10, step = 0.5
                        ),
                        sliderInput("libsizePointAlpha", "Point Transparency",
                            value = 0.6, min = 0, max = 1, step = 0.1
                        ),
                        numericInput("libsizeBins", "Histogram Bins",
                            value = 30, min = 10, max = 100, step = 5
                        ),
                        selectInput("libsizeTheme", "Theme",
                            choices = list(
                                "Black & White" = "bw",
                                "Classic" = "classic",
                                "Minimal" = "minimal",
                                "Light" = "light",
                                "Dark" = "dark"
                            ),
                            selected = "bw"
                        ),
                        selectInput("libsizeLegendPos", "Legend Position",
                            choices = list(
                                "Right" = "right",
                                "Left" = "left",
                                "Top" = "top",
                                "Bottom" = "bottom",
                                "None" = "none"
                            ),
                            selected = "right"
                        ),
                        actionButton("updateLibsize", "Update Plot", class = "btn-primary")
                    )
                ),
                br(), br(), hr(), br(),
                h4("Percentage of Zeros"),
                fluidRow(
                    column(
                        6,
                        plotOutput("zerosPlot", height = "400px"),
                        downloadButton("downloadZeros", "Download Plot")
                    ),
                    column(
                        6,
                        h5("Styling Options"),
                        numericInput("zerosBins", "Histogram Bins",
                            value = 30, min = 10, max = 100, step = 5
                        ),
                        selectInput("zerosTheme", "Theme",
                            choices = list(
                                "Black & White" = "bw",
                                "Classic" = "classic",
                                "Minimal" = "minimal",
                                "Light" = "light",
                                "Dark" = "dark"
                            ),
                            selected = "bw"
                        ),
                        actionButton("updateZeros", "Update Plot", class = "btn-primary")
                    )
                ),
                br(), br(), hr(), br(),
                h4("Principal Component Analysis (PCA)"),
                fluidRow(
                    column(
                        6,
                        plotOutput("pcaPlot", height = "500px"),
                        downloadButton("downloadPCA", "Download Plot")
                    ),
                    column(
                        6,
                        h5("Styling Options"),
                        numericInput("pcaPointSize", "Point Size",
                            value = 4, min = 1, max = 10, step = 0.5
                        ),
                        sliderInput("pcaPointAlpha", "Point Transparency",
                            value = 0.8, min = 0, max = 1, step = 0.1
                        ),
                        selectInput("pcaTheme", "Theme",
                            choices = list(
                                "Black & White" = "bw",
                                "Classic" = "classic",
                                "Minimal" = "minimal",
                                "Light" = "light",
                                "Dark" = "dark"
                            ),
                            selected = "bw"
                        ),
                        selectInput("pcaLegendPos", "Legend Position",
                            choices = list(
                                "Right" = "right",
                                "Left" = "left",
                                "Top" = "top",
                                "Bottom" = "bottom",
                                "None" = "none"
                            ),
                            selected = "right"
                        ),
                        actionButton("updatePCA", "Update Plot", class = "btn-primary")
                    )
                )
            ),

            # DE Tab
            tabPanel(
                "Differential Expression",
                h3("Differential Expression Analysis"),
                br(),

                # DE Analysis Container at Top
                wellPanel(
                    h4("Run Differential Expression Analysis"),
                    p("Genes are selected by FDR (ascending), then by absolute log2FC (descending)."),
                    conditionalPanel(
                        condition = "output.dataLoaded == false",
                        p(strong("Please load data first in the sidebar."), style = "color: red;")
                    ),
                    conditionalPanel(
                        condition = "output.dataLoaded == true && output.deRun == false",
                        actionButton("runDE", "Run DE Analysis", class = "btn-success btn-lg"),
                        br(), br(),
                        verbatimTextOutput("deStatusMessage")
                    ),
                    conditionalPanel(
                        condition = "output.deRun == true",
                        p(strong("DE Analysis Complete!"), style = "color: green;"),
                        p("All plots have been generated with default settings. Use the sections below to customize and regenerate individual plots."),
                        verbatimTextOutput("deStatusMessage")
                    )
                ),
                br(), br(),
                conditionalPanel(
                    condition = "output.deRun == true",
                    # MA Plot Section
                    h4("MA Plot"),
                    conditionalPanel(
                        condition = "output.maPlotGenerated == true",
                        plotOutput("maPlot", height = "500px"),
                        downloadButton("downloadMA", "Download Plot"),
                        br(), br()
                    ),
                    tags$details(
                        tags$summary(h5("Customize MA Plot Settings", style = "cursor: pointer; padding: 8px; border: 2px solid #337ab7; border-radius: 4px; background-color: #f5f5f5;")),
                        wellPanel(
                            style = "border: 2px solid #337ab7; margin-top: 10px;",
                            p("Adjust thresholds and styling options, then click 'Update Plot' to regenerate with new settings."),
                            fluidRow(
                                column(
                                    6,
                                    numericInput("maFDR", "FDR Threshold",
                                        value = 0.05, min = 0.001, max = 1, step = 0.01
                                    ),
                                    numericInput("maLog2FC", "log2FC Threshold",
                                        value = 1, min = 0, max = 10, step = 0.1
                                    )
                                ),
                                column(
                                    6,
                                    numericInput("maPointSize", "Point Size",
                                        value = 1.5, min = 0.5, max = 5, step = 0.5
                                    ),
                                    sliderInput("maPointAlpha", "Point Transparency",
                                        value = 0.6, min = 0, max = 1, step = 0.1
                                    ),
                                    selectInput("maTheme", "Theme",
                                        choices = list(
                                            "Black & White" = "bw",
                                            "Classic" = "classic",
                                            "Minimal" = "minimal",
                                            "Light" = "light",
                                            "Dark" = "dark"
                                        ),
                                        selected = "bw"
                                    ),
                                    selectInput("maLegendPos", "Legend Position",
                                        choices = list(
                                            "Right" = "right",
                                            "Left" = "left",
                                            "Top" = "top",
                                            "Bottom" = "bottom",
                                            "None" = "none"
                                        ),
                                        selected = "right"
                                    )
                                )
                            ),
                            actionButton("generateMA", "Update MA Plot", class = "btn-primary")
                        )
                    ),
                    br(), br(), hr(), br(),

                    # Volcano Plot Section
                    h4("Volcano Plot"),
                    conditionalPanel(
                        condition = "output.volcanoPlotGenerated == true",
                        plotOutput("volcanoPlot", height = "500px"),
                        downloadButton("downloadVolcano", "Download Plot"),
                        br(), br()
                    ),
                    tags$details(
                        tags$summary(h5("Customize Volcano Plot Settings", style = "cursor: pointer; padding: 8px; border: 2px solid #337ab7; border-radius: 4px; background-color: #f5f5f5;")),
                        wellPanel(
                            style = "border: 2px solid #337ab7; margin-top: 10px;",
                            p("Adjust thresholds and styling options, then click 'Update Plot' to regenerate with new settings."),
                            fluidRow(
                                column(
                                    6,
                                    numericInput("volcanoFDR", "FDR Threshold",
                                        value = 0.05, min = 0.001, max = 1, step = 0.01
                                    ),
                                    numericInput("volcanoLog2FC", "log2FC Threshold",
                                        value = 1, min = 0, max = 10, step = 0.1
                                    )
                                ),
                                column(
                                    6,
                                    numericInput("volcanoPointSize", "Point Size",
                                        value = 1.5, min = 0.5, max = 5, step = 0.5
                                    ),
                                    sliderInput("volcanoPointAlpha", "Point Transparency",
                                        value = 0.6, min = 0, max = 1, step = 0.1
                                    ),
                                    selectInput("volcanoTheme", "Theme",
                                        choices = list(
                                            "Black & White" = "bw",
                                            "Classic" = "classic",
                                            "Minimal" = "minimal",
                                            "Light" = "light",
                                            "Dark" = "dark"
                                        ),
                                        selected = "bw"
                                    ),
                                    selectInput("volcanoLegendPos", "Legend Position",
                                        choices = list(
                                            "Right" = "right",
                                            "Left" = "left",
                                            "Top" = "top",
                                            "Bottom" = "bottom",
                                            "None" = "none"
                                        ),
                                        selected = "right"
                                    )
                                )
                            ),
                            actionButton("generateVolcano", "Update Volcano Plot", class = "btn-primary")
                        )
                    ),
                    br(), br(), hr(), br(),

                    # Heatmap Section
                    h4("Heatmap"),
                    conditionalPanel(
                        condition = "output.heatmapPlotGenerated == true",
                        plotOutput("heatmapPlot", height = "600px"),
                        downloadButton("downloadHeatmap", "Download Plot"),
                        br(), br()
                    ),
                    tags$details(
                        tags$summary(h5("Customize Heatmap Settings", style = "cursor: pointer; padding: 8px; border: 2px solid #337ab7; border-radius: 4px; background-color: #f5f5f5;")),
                        wellPanel(
                            style = "border: 2px solid #337ab7; margin-top: 10px;",
                            p("Adjust number of genes, scaling, and display options. Click 'Update Plot' to regenerate with new settings."),
                            fluidRow(
                                column(
                                    4,
                                    numericInput("heatmapGenes", "Number of Genes",
                                        value = 50, min = 10, max = 500, step = 10
                                    )
                                ),
                                column(
                                    4,
                                    selectInput("heatmapScale", "Scaling Method",
                                        choices = list(
                                            "Row (z-score per gene)" = "row",
                                            "Column (z-score per sample)" = "column",
                                            "None" = "none"
                                        ),
                                        selected = "row"
                                    )
                                ),
                                column(
                                    4,
                                    checkboxInput("heatmapShowRowNames", "Show Gene Names", value = FALSE),
                                    br(),
                                    actionButton("generateHeatmap", "Update Heatmap", class = "btn-primary")
                                )
                            )
                        )
                    ),
                    br(), br(), hr(), br(),

                    # DE Results Table
                    h4("DE Results Table"),
                    p("Showing top 50 genes. Use export to download full table."),
                    tableOutput("deTable")
                )
            )
        )

        # Add Export Tab
        tabs <- append(tabs, list(
            tabPanel(
                "Export",
                h3("Export Results"),
                br(),
                p("Download analysis results and plots:"),
                br(),
                h4("Analysis Bundle"),
                p("Complete analysis bundle (RDS format):"),
                div(
                    style = "margin-bottom: 15px;",
                    div(
                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                        downloadButton("downloadBundle", "Download Analysis Bundle"),
                        tags$details(
                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                            div(
                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                h5("What's in this export:"),
                                p("The RDS file contains a complete analysis bundle with:"),
                                tags$ul(
                                    tags$li(strong("counts_t:"), " Transformed/normalized expression matrix"),
                                    tags$li(strong("de_table:"), " Differential expression statistics table"),
                                    tags$li(strong("selected_genes:"), " List of selected genes (top-K or DE-based)"),
                                    tags$li(strong("qc_metrics:"), " Quality control metrics (library sizes, zero percentages, etc.)")
                                ),
                                br(),
                                h5("What you can do with this export:"),
                                tags$ul(
                                    tags$li("Load complete analysis state in R for further analysis"),
                                    tags$li("Share complete analysis with collaborators"),
                                    tags$li("Reproduce all plots and visualizations"),
                                    tags$li("Extract individual components for custom analysis"),
                                    tags$li("Archive analysis results for reproducibility")
                                ),
                                br(),
                                h5("Format specifications:"),
                                p("RDS format (R Data Serialization), binary format. Can only be loaded in R using readRDS(). Contains complete analysis bundle object.")
                            )
                        )
                    )
                ),
                br(), br(),
                h4("Individual Components"),
                div(
                    style = "margin-bottom: 15px;",
                    div(
                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                        downloadButton("downloadCountsT", "Download Transformed Counts (TSV)"),
                        tags$details(
                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                            div(
                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                h5("What's in this export:"),
                                p("The TSV file contains the normalized/transformed expression matrix with:"),
                                tags$ul(
                                    tags$li(strong("Rows:"), " Genes (gene identifiers in first column)"),
                                    tags$li(strong("Columns:"), " Samples (sample identifiers as column headers)"),
                                    tags$li(strong("Values:"), " Transformed expression values (log1p-CPM or VST)"),
                                    tags$li(strong("Transformation method:"), " Specified in preprocessing (log1p-CPM or VST)")
                                ),
                                br(),
                                h5("What you can do with this export:"),
                                tags$ul(
                                    tags$li("Use transformed counts in other analysis tools (R, Python, etc.)"),
                                    tags$li("Create custom visualizations or heatmaps"),
                                    tags$li("Perform additional statistical analyses"),
                                    tags$li("Share normalized data with collaborators"),
                                    tags$li("Use as input for other machine learning pipelines")
                                ),
                                br(),
                                h5("Format specifications:"),
                                p("TSV format (tab-separated values), UTF-8 encoding. Genes as rows, samples as columns. Can be opened in Excel, R, Python, or any text editor.")
                            )
                        )
                    )
                ),
                br(),
                div(
                    style = "margin-bottom: 15px;",
                    div(
                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                        downloadButton("downloadDETable", "Download DE Table (TSV)"),
                        tags$details(
                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                            div(
                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                h5("What's in this export:"),
                                p("The TSV file contains differential expression statistics with:"),
                                tags$ul(
                                    tags$li(strong("gene_id:"), " Gene identifier"),
                                    tags$li(strong("log2FC:"), " Log2 fold change between groups"),
                                    tags$li(strong("pvalue:"), " P-value from differential expression test"),
                                    tags$li(strong("FDR:"), " False discovery rate (adjusted p-value)"),
                                    tags$li(strong("AveExpr:"), " Average expression across all samples")
                                ),
                                br(),
                                h5("What you can do with this export:"),
                                tags$ul(
                                    tags$li("Identify significantly differentially expressed genes"),
                                    tags$li("Create custom volcano plots or MA plots"),
                                    tags$li("Perform gene set enrichment analysis (GSEA)"),
                                    tags$li("Filter genes by fold change or significance thresholds"),
                                    tags$li("Compare DE results across different analyses")
                                ),
                                br(),
                                h5("Format specifications:"),
                                p("TSV format (tab-separated values), UTF-8 encoding. One row per gene with DE statistics. Can be opened in Excel, R, Python, or any text editor.")
                            )
                        )
                    )
                ),
                br(),
                div(
                    style = "margin-bottom: 15px;",
                    div(
                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                        downloadButton("downloadQCMetrics", "Download QC Metrics (TSV)"),
                        tags$details(
                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                            div(
                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                h5("What's in this export:"),
                                p("The TSV file contains quality control metrics for each sample with:"),
                                tags$ul(
                                    tags$li(strong("sample_id:"), " Sample identifier"),
                                    tags$li(strong("library_size:"), " Total read count (library size)"),
                                    tags$li(strong("zero_percentage:"), " Percentage of genes with zero counts"),
                                    tags$li(strong("median_count:"), " Median count across all genes"),
                                    tags$li(strong("other_metrics:"), " Additional QC metrics as computed")
                                ),
                                br(),
                                h5("What you can do with this export:"),
                                tags$ul(
                                    tags$li("Assess data quality and identify problematic samples"),
                                    tags$li("Compare QC metrics across samples or batches"),
                                    tags$li("Create custom QC visualizations"),
                                    tags$li("Filter samples based on QC thresholds"),
                                    tags$li("Document data quality for publications")
                                ),
                                br(),
                                h5("Format specifications:"),
                                p("TSV format (tab-separated values), UTF-8 encoding. One row per sample with QC metrics. Can be opened in Excel, R, Python, or any text editor.")
                            )
                        )
                    )
                ),
                br(),
                div(
                    style = "margin-bottom: 15px;",
                    div(
                        style = "display: flex; align-items: center; gap: 10px; margin-bottom: 5px;",
                        downloadButton("downloadSelectedGenes", "Download Selected Genes (TSV)"),
                        tags$details(
                            tags$summary(style = "cursor: pointer; color: #3498db; font-weight: bold;", "ℹ️ More Information"),
                            div(
                                style = "margin-top: 10px; padding: 15px; background-color: #f8f9fa; border-left: 4px solid #3498db;",
                                h5("What's in this export:"),
                                p("The TSV file contains the list of selected genes with:"),
                                tags$ul(
                                    tags$li(strong("gene_id:"), " Gene identifier"),
                                    tags$li(strong("selection_criteria:"), " How gene was selected (top-K, DE-based, variance-based, etc.)"),
                                    tags$li(strong("rank:"), " Rank or score used for selection (if applicable)")
                                ),
                                br(),
                                h5("What you can do with this export:"),
                                tags$ul(
                                    tags$li("Use selected genes for downstream analysis"),
                                    tags$li("Create gene lists for pathway analysis"),
                                    tags$li("Compare selected genes across different analyses"),
                                    tags$li("Use as input for signature training or other methods"),
                                    tags$li("Share gene lists with collaborators")
                                ),
                                br(),
                                h5("Format specifications:"),
                                p("TSV format (tab-separated values), UTF-8 encoding. One row per selected gene. Can be opened in Excel, R, Python, or any text editor.")
                            )
                        )
                    )
                )
            )
        ))

        # Return tabsetPanel with all tabs
        do.call(tabsetPanel, c(list(id = "mode3Tabs"), tabs))
    })

    # ===== MODE 1: RAPID CLASSIFICATION =====

    # Load Mode 1 demo data
    observeEvent(input$loadDemoMode1, {
        tryCatch(
            {
                demo_data <- endoSignatureR::endo_load_demo()
                # For Mode 1, use only counts (unlabeled)
                values$counts_mode1 <- demo_data$counts
                values$annot_mode1 <- demo_data$annot

                output$dataStatusMode1 <- renderText({
                    paste(
                        "Demo data loaded successfully (unlabeled)!\n",
                        "Samples:", ncol(values$counts_mode1), "\n",
                        "Genes:", nrow(values$counts_mode1)
                    )
                })

                output$statusMessageMode1 <- renderText("Demo data loaded successfully!")
            },
            error = function(e) {
                output$statusMessageMode1 <- renderText(paste("Error loading demo data:", e$message))
            }
        )
    })

    # Download Mode 1 demo data
    output$downloadDemoMode1 <- downloadHandler(
        filename = function() {
            "gse201926_sample_unlabeled.zip"
        },
        content = function(file) {
            temp_dir <- tempdir()
            demo_data <- endoSignatureR::endo_load_demo()

            # Write counts (unlabeled)
            # Convert counts matrix to data.frame with GeneID column
            # This matches the format expected by esr_loadCountsFromFile
            counts_df <- data.frame(
                GeneID = rownames(demo_data$counts),
                as.data.frame(demo_data$counts),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            readr::write_tsv(
                counts_df,
                file.path(temp_dir, "gse201926_sample_unlabeled_counts.tsv")
            )

            # Write annot (optional)
            readr::write_tsv(
                demo_data$annot,
                file.path(temp_dir, "gse201926_sample_annot.tsv")
            )

            # Create zip file
            utils::zip(file,
                files = c(
                    file.path(temp_dir, "gse201926_sample_unlabeled_counts.tsv"),
                    file.path(temp_dir, "gse201926_sample_annot.tsv")
                ),
                flags = "-j"
            )
        }
    )

    # Load Mode 1 uploaded data
    observeEvent(
        {
            input$countsFileMode1
            input$annotFileMode1
        },
        {
            if (!is.null(input$countsFileMode1)) {
                tryCatch(
                    {
                        # Load counts
                        values$counts_mode1 <- endoSignatureR::esr_loadCountsFromFile(input$countsFileMode1$datapath)

                        # Load annot if provided
                        if (!is.null(input$annotFileMode1)) {
                            values$annot_mode1 <- endoSignatureR::esr_loadAnnotFromFile(input$annotFileMode1$datapath)
                        } else {
                            values$annot_mode1 <- NULL
                        }

                        # Validate data (pheno is optional for unlabeled data)
                        validation <- endoSignatureR::esr_validateEndometrial(
                            values$counts_mode1,
                            pheno = NULL, values$annot_mode1
                        )

                        values$counts_mode1 <- validation$X
                        values$annot_mode1 <- validation$annot

                        output$dataStatusMode1 <- renderText({
                            paste(
                                "Data loaded successfully!\n",
                                "Samples:", ncol(values$counts_mode1), "\n",
                                "Genes:", nrow(values$counts_mode1), "\n",
                                if (nrow(validation$issues) > 0) {
                                    paste("\nWarnings:", nrow(validation$issues))
                                } else {
                                    "\nNo validation issues"
                                }
                            )
                        })

                        output$statusMessageMode1 <- renderText("Data loaded and validated successfully!")
                    },
                    error = function(e) {
                        output$statusMessageMode1 <- renderText(paste("Error loading data:", e$message))
                    }
                )
            }
        }
    )

    # Run classification
    observeEvent(input$runClassification, {
        if (is.null(values$counts_mode1)) {
            output$statusMessageMode1 <- renderText("Please load data first!")
            return()
        }

        if (is.null(values$selectedSignature)) {
            sig_type <- if (input$signatureSelection == "pretrained") {
                "Pre-trained"
            } else if (input$signatureSelection == "user-trained") {
                "User-trained"
            } else if (input$signatureSelection == "uploaded") {
                "Uploaded"
            } else {
                "Selected"
            }
            output$statusMessageMode1 <- renderText(paste(sig_type, "signature not available. Please ensure package is properly installed, train a signature in Mode 2, or upload signature files."))
            return()
        }

        tryCatch(
            {
                output$statusMessageMode1 <- renderText("Running classification...")

                # Run classification with selected signature
                values$predictions <- endoSignatureR::esr_classifyEndometrial(
                    X_new = values$counts_mode1,
                    signature = values$selectedSignature,
                    threshold = input$thresholdMode1,
                    confidence = input$confidenceMode1
                )

                # Set completion flags
                values$classificationRun <- TRUE
                values$mode1_completed <- TRUE

                sig_type <- if (input$signatureSelection == "pretrained") {
                    "Pre-trained"
                } else if (input$signatureSelection == "user-trained") {
                    "User-trained"
                } else if (input$signatureSelection == "uploaded") {
                    "Uploaded"
                } else {
                    "Selected"
                }
                output$statusMessageMode1 <- renderText(paste("Classification completed successfully using", sig_type, "signature!"))

                # Automatically navigate to Results tab
                updateTabsetPanel(session, "mode1Tabs", selected = "Results")
            },
            error = function(e) {
                output$statusMessageMode1 <- renderText(paste("Error in classification:", e$message))
            }
        )
    })

    # Classification summary
    output$classificationSummary <- renderText({
        if (is.null(values$predictions)) {
            return(NULL)
        }

        n_ps <- sum(values$predictions$prediction == "PS", na.rm = TRUE)
        n_pis <- sum(values$predictions$prediction == "PIS", na.rm = TRUE)
        n_total <- nrow(values$predictions)

        summary_text <- paste(
            "Total samples:", n_total, "\n",
            "PS predictions:", n_ps, "\n",
            "PIS predictions:", n_pis
        )

        if ("confidence_interval_lower" %in% names(values$predictions)) {
            mean_conf <- mean(
                values$predictions$confidence_interval_upper - values$predictions$confidence_interval_lower,
                na.rm = TRUE
            )
            summary_text <- paste(
                summary_text, "\n",
                "Mean confidence interval width:", round(mean_conf, 3)
            )
        }

        return(summary_text)
    })

    # Predictions table
    output$predictionsTable <- renderTable(
        {
            if (is.null(values$predictions)) {
                return(NULL)
            }
            values$predictions
        },
        digits = 4
    )

    # Download predictions
    output$downloadPredictions <- downloadHandler(
        filename = "predictions.csv",
        content = function(file) {
            if (!is.null(values$predictions)) {
                readr::write_csv(values$predictions, file)
            }
        }
    )

    # ===== MODE 2: SIGNATURE VALIDATION =====

    # Load Mode 2 demo data
    observeEvent(input$loadDemoMode2, {
        tryCatch(
            {
                data(gse201926_trainmini, package = "endoSignatureR", envir = environment())
                demo_data <- get("gse201926_trainmini", envir = environment())
                values$counts_mode2 <- demo_data$counts
                values$pheno_mode2 <- demo_data$pheno
                values$annot_mode2 <- demo_data$annot

                # Validate data
                validation <- endoSignatureR::esr_validateEndometrial(
                    values$counts_mode2,
                    values$pheno_mode2,
                    values$annot_mode2
                )

                values$counts_mode2 <- validation$X
                values$pheno_mode2 <- validation$pheno
                values$annot_mode2 <- validation$annot

                output$dataStatusMode2 <- renderText({
                    paste(
                        "Demo data loaded successfully (labeled training data)!\n",
                        "Samples:", ncol(values$counts_mode2), "\n",
                        "Genes:", nrow(values$counts_mode2), "\n",
                        "Groups:", paste(unique(values$pheno_mode2$group), collapse = ", "), "\n",
                        if (nrow(validation$issues) > 0) {
                            paste("\nWarnings:", nrow(validation$issues))
                        } else {
                            "\nNo validation issues"
                        }
                    )
                })

                output$statusMessageMode2 <- renderText("Demo data loaded successfully!")
            },
            error = function(e) {
                output$statusMessageMode2 <- renderText(paste("Error loading demo data:", e$message))
            }
        )
    })

    # Download Mode 2 demo data
    output$downloadDemoMode2 <- downloadHandler(
        filename = function() {
            "gse201926_trainmini.zip"
        },
        content = function(file) {
            temp_dir <- tempdir()
            data(gse201926_trainmini, package = "endoSignatureR", envir = environment())
            demo_data <- get("gse201926_trainmini", envir = environment())

            # Convert counts matrix to data.frame with GeneID column
            counts_df <- data.frame(
                GeneID = rownames(demo_data$counts),
                as.data.frame(demo_data$counts),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            readr::write_tsv(
                counts_df,
                file.path(temp_dir, "gse201926_trainmini_counts.tsv")
            )
            readr::write_tsv(
                demo_data$pheno,
                file.path(temp_dir, "gse201926_trainmini_pheno.tsv")
            )
            readr::write_tsv(
                demo_data$annot,
                file.path(temp_dir, "gse201926_trainmini_annot.tsv")
            )

            utils::zip(file,
                files = c(
                    file.path(temp_dir, "gse201926_trainmini_counts.tsv"),
                    file.path(temp_dir, "gse201926_trainmini_pheno.tsv"),
                    file.path(temp_dir, "gse201926_trainmini_annot.tsv")
                ),
                flags = "-j"
            )
        }
    )

    # Load Mode 2 uploaded data
    observeEvent(
        {
            input$countsFileMode2
            input$phenoFileMode2
            input$annotFileMode2
        },
        {
            if (!is.null(input$countsFileMode2) && !is.null(input$phenoFileMode2)) {
                tryCatch(
                    {
                        values$counts_mode2 <- endoSignatureR::esr_loadCountsFromFile(input$countsFileMode2$datapath)
                        values$pheno_mode2 <- endoSignatureR::esr_loadPhenoFromFile(input$phenoFileMode2$datapath)

                        if (!is.null(input$annotFileMode2)) {
                            values$annot_mode2 <- endoSignatureR::esr_loadAnnotFromFile(input$annotFileMode2$datapath)
                        } else {
                            values$annot_mode2 <- NULL
                        }

                        # Validate data
                        validation <- endoSignatureR::esr_validateEndometrial(
                            values$counts_mode2,
                            values$pheno_mode2,
                            values$annot_mode2
                        )

                        values$counts_mode2 <- validation$X
                        values$pheno_mode2 <- validation$pheno
                        values$annot_mode2 <- validation$annot

                        output$dataStatusMode2 <- renderText({
                            paste(
                                "Data loaded successfully!\n",
                                "Samples:", ncol(values$counts_mode2), "\n",
                                "Genes:", nrow(values$counts_mode2), "\n",
                                "Groups:", paste(unique(values$pheno_mode2$group), collapse = ", "), "\n",
                                if (nrow(validation$issues) > 0) {
                                    paste("\nWarnings:", nrow(validation$issues))
                                } else {
                                    "\nNo validation issues"
                                }
                            )
                        })

                        output$statusMessageMode2 <- renderText("Data loaded and validated successfully!")
                    },
                    error = function(e) {
                        output$statusMessageMode2 <- renderText(paste("Error loading data:", e$message))
                    }
                )
            }
        }
    )

    # Run training
    observeEvent(input$runTraining, {
        if (is.null(values$counts_mode2) || is.null(values$pheno_mode2)) {
            output$statusMessageMode2 <- renderText("Please load data first (counts and pheno required)!")
            return()
        }

        tryCatch(
            {
                output$statusMessageMode2 <- renderText("Training in progress... This may take 30-60 seconds.")

                # Prepare parameters
                outer_folds <- if (is.null(input$outerFoldsMode2) || is.na(input$outerFoldsMode2)) {
                    NULL
                } else {
                    input$outerFoldsMode2
                }

                stability_selection <- if (input$stabilityMode2) {
                    TRUE
                } else {
                    NULL # Auto-detect based on sample size
                }

                # Run training with progress indicator
                withProgress(
                    message = "Training signature",
                    detail = "This may take 30-60 seconds...",
                    value = 0,
                    {
                        setProgress(0.1, detail = "Starting nested cross-validation...")

                        # Suppress expected warnings about consensus genes and sample sizes
                        result <- suppressWarnings(
                            endoSignatureR::esr_trainEndometrialSignature(
                                X = values$counts_mode2,
                                pheno = values$pheno_mode2,
                                transform = input$transformMode2,
                                cpm_min = input$cpmMinMode2,
                                cpm_min_samples = input$cpmMinSamplesMode2,
                                top_k = input$topKMode2,
                                outer = input$outerMode2,
                                outer_folds = outer_folds,
                                inner_folds = input$innerFoldsMode2,
                                inner_repeats = input$innerRepeatsMode2,
                                lambda_rule = input$lambdaRuleMode2,
                                calibration_method = input$calibrationMode2,
                                stability_selection = stability_selection,
                                stability_resamples = input$stabilityResamplesMode2,
                                seed = input$seedMode2
                            )
                        )

                        setProgress(1.0, detail = "Training complete!")
                    }
                )

                values$trainingResult <- result
                values$trainingRun <- TRUE
                values$mode2_completed <- TRUE

                # Check for small dataset warnings and provide informative message
                n_genes <- length(result$signature$panel)
                n_samples <- ncol(values$counts_mode2)

                status_msg <- "Training completed successfully!"
                if (n_genes <= 1 || n_samples < 10) {
                    status_msg <- paste0(
                        status_msg,
                        "\n\nNote: With small sample sizes (< 10 samples) or small gene sets, ",
                        "you may see: (1) No consensus genes found across folds (genes from a single fold will be used), ",
                        "and (2) Bootstrap stability frequencies will not be computed. ",
                        "These are expected outcomes with limited data."
                    )
                }

                output$statusMessageMode2 <- renderText(status_msg)

                # Automatically navigate to Overview tab
                updateTabsetPanel(session, "mode2Tabs", selected = "Overview")
            },
            error = function(e) {
                output$statusMessageMode2 <- renderText(paste("Error in training:", e$message))
            }
        )
    })

    # Training summary
    output$trainingSummary <- renderText({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }

        result <- values$trainingResult
        paste(
            "Samples used:", ncol(values$counts_mode2), "\n",
            "Genes after filtering:", length(result$signature$panel), "\n",
            "Signature panel size:", length(result$signature$panel), "genes\n",
            "Outer CV folds:", result$splits$n_outer_folds, "\n",
            "Transform method:", result$signature$recipe$preprocessing$transform %||% "log1p-cpm", "\n",
            "Calibration method:", result$calibration$method, "\n",
            "Lambda rule:", if (!is.null(result$cv_results) && length(result$cv_results) > 0) {
                "1se/min (varies by fold)"
            } else {
                "N/A"
            }
        )
    })

    # Performance metrics
    output$performanceMetrics <- renderText({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }

        metrics <- values$trainingResult$metrics
        paste(
            "AUC:", round(metrics$auc, 4), "\n",
            "Accuracy:", round(metrics$accuracy, 4), "\n",
            "Brier Score:", round(metrics$brier_score, 4), "\n",
            "ECE (Expected Calibration Error):", round(metrics$ece, 4)
        )
    })

    # Signature panel information
    output$signaturePanelInfo <- renderText({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }

        sig <- values$trainingResult$signature
        paste(
            "Number of genes in panel:", length(sig$panel), "\n",
            "Preprocessing method:", sig$recipe$preprocessing$transform %||% "log1p-cpm", "\n",
            "Calibration method:", values$trainingResult$calibration$method
        )
    })

    # ROC plot
    output$rocPlotMode2 <- renderPlot({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }
        use_calibrated <- input$rocCalibratedMode2 == "TRUE"
        endoSignatureR::plotEndometrialROC(
            values$trainingResult$metrics$predictions,
            use_calibrated = use_calibrated,
            show_auc = TRUE
        )
    })

    # PR plot
    output$prPlotMode2 <- renderPlot({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }
        use_calibrated <- input$prCalibratedMode2 == "TRUE"
        endoSignatureR::plotEndometrialPR(
            values$trainingResult$metrics$predictions,
            use_calibrated = use_calibrated,
            show_auc = TRUE
        )
    })

    # Calibration plot
    output$calibrationPlotMode2 <- renderPlot({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }
        use_calibrated <- input$calCalibratedMode2 == "TRUE"
        endoSignatureR::plotEndometrialCalibration(
            values$trainingResult$metrics$predictions,
            use_calibrated = use_calibrated
        )
    })

    # Comparison plot
    output$comparisonPlotMode2 <- renderPlot({
        if (is.null(values$trainingResult) || is.null(values$signature)) {
            return(NULL)
        }

        # For comparison, we need predictions from pre-trained signature
        # Since we don't have pre-trained predictions readily available,
        # we'll show the new signature plots with a note
        # In a full implementation, we would apply pre-trained signature to the same data
        tryCatch(
            {
                # Try to create comparison plot
                # Note: plotEndometrialComparison expects pretrained_result and new_result
                # We may need to adapt this based on available data
                endoSignatureR::plotEndometrialComparison(
                    pretrained_result = NULL, # Would need pre-trained predictions
                    new_result = values$trainingResult
                )
            },
            error = function(e) {
                # If comparison fails, show new signature plots
                endoSignatureR::plotEndometrialROC(
                    values$trainingResult$metrics$predictions,
                    use_calibrated = TRUE,
                    show_auc = TRUE
                )
            }
        )
    })

    # Comparison metrics
    output$comparisonMetrics <- renderText({
        if (is.null(values$trainingResult) || is.null(values$signature)) {
            return("Pre-trained signature not available for comparison.")
        }

        # Try to compute comparison metrics
        tryCatch(
            {
                # Note: esr_compareSignatures may be a placeholder
                # This would need pre-trained predictions to work properly
                comparison <- endoSignatureR::esr_compareSignatures(
                    pretrained_result = NULL, # Would need pre-trained predictions
                    new_result = values$trainingResult
                )
                paste(
                    "Comparison metrics:\n",
                    "Delta AUC:", if (!is.na(comparison$delta_auc)) round(comparison$delta_auc, 4) else "N/A"
                )
            },
            error = function(e) {
                "Comparison metrics require pre-trained signature predictions on the same data."
            }
        )
    })

    # Coefficient lollipop plot
    output$coefLollipopPlotMode2 <- renderPlot({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialCoefLollipop(
            signature = values$trainingResult$signature,
            annot = values$annot_mode2
        )
    })

    # Stability bars plot
    output$stabilityBarsPlotMode2 <- renderPlot({
        if (is.null(values$trainingResult)) {
            return(NULL)
        }
        if (is.null(values$trainingResult$signature$selection_frequency) &&
            is.null(values$trainingResult$stability)) {
            return(NULL)
        }

        # Create a signature object that includes stability info if available
        # The function expects signature$stability$bootstrap_frequency, but in trainingResult
        # bootstrap_frequency is at trainingResult$stability$bootstrap_frequency
        sig_for_plot <- values$trainingResult$signature

        # Check if bootstrap_frequency exists in trainingResult$stability
        has_bootstrap <- !is.null(values$trainingResult$stability) &&
            !is.null(values$trainingResult$stability$bootstrap_frequency)

        # If bootstrap exists, add it to signature$stability for the function
        if (has_bootstrap) {
            sig_for_plot$stability <- list(
                bootstrap_frequency = values$trainingResult$stability$bootstrap_frequency
            )
            freq_type <- "bootstrap"
        } else {
            # Use selection frequency (which is already in signature$selection_frequency)
            freq_type <- "selection"
        }

        # Try to create the plot with error handling
        tryCatch(
            {
                # Suppress warnings about bootstrap_frequency not available
                suppressWarnings({
                    endoSignatureR::plotEndometrialStabilityBars(
                        signature = sig_for_plot,
                        annot = values$annot_mode2,
                        frequency_type = freq_type
                    )
                })
            },
            error = function(e) {
                # If no valid frequencies found, show a message plot instead
                if (grepl("No valid frequencies found", e$message, fixed = TRUE)) {
                    ggplot2::ggplot() +
                        ggplot2::annotate("text",
                            x = 0.5, y = 0.5,
                            label = "Stability frequency data not available.\nAll frequency values are missing or invalid.",
                            size = 5, hjust = 0.5, vjust = 0.5
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Stability Bars Plot")
                } else if (grepl("bootstrap_frequency not available", e$message, fixed = TRUE)) {
                    # Handle case where bootstrap frequency is requested but not available
                    ggplot2::ggplot() +
                        ggplot2::annotate("text",
                            x = 0.5, y = 0.5,
                            label = "Bootstrap frequency data not available.\nUsing selection frequency instead.\nNote: Small sample sizes (< 10 samples) prevent bootstrap resampling.",
                            size = 4, hjust = 0.5, vjust = 0.5
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Stability Bars Plot")
                } else {
                    # For other errors, show error message
                    ggplot2::ggplot() +
                        ggplot2::annotate("text",
                            x = 0.5, y = 0.5,
                            label = paste("Error creating stability plot:\n", e$message),
                            size = 4, hjust = 0.5, vjust = 0.5
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Stability Bars Plot - Error")
                }
            }
        )
    })

    # Download handlers for Mode 2 plots
    output$downloadROCMode2 <- downloadHandler(
        filename = "roc_plot.png",
        content = function(file) {
            use_calibrated <- input$rocCalibratedMode2 == "TRUE"
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialROC(
                    values$trainingResult$metrics$predictions,
                    use_calibrated = use_calibrated,
                    show_auc = TRUE
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadPRMode2 <- downloadHandler(
        filename = "pr_plot.png",
        content = function(file) {
            use_calibrated <- input$prCalibratedMode2 == "TRUE"
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialPR(
                    values$trainingResult$metrics$predictions,
                    use_calibrated = use_calibrated,
                    show_auc = TRUE
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadCalibrationMode2 <- downloadHandler(
        filename = "calibration_plot.png",
        content = function(file) {
            use_calibrated <- input$calCalibratedMode2 == "TRUE"
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialCalibration(
                    values$trainingResult$metrics$predictions,
                    use_calibrated = use_calibrated
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadComparisonMode2 <- downloadHandler(
        filename = "comparison_plot.png",
        content = function(file) {
            if (is.null(values$trainingResult) || is.null(values$signature)) {
                return()
            }
            # Recreate the comparison plot
            p <- tryCatch(
                {
                    endoSignatureR::plotEndometrialComparison(
                        pretrained_result = NULL,
                        new_result = values$trainingResult
                    )
                },
                error = function(e) {
                    endoSignatureR::plotEndometrialROC(
                        values$trainingResult$metrics$predictions,
                        use_calibrated = TRUE,
                        show_auc = TRUE
                    )
                }
            )
            ggplot2::ggsave(file, plot = p, width = 12, height = 10, dpi = 300)
        }
    )

    output$downloadCoefLollipopMode2 <- downloadHandler(
        filename = "coef_lollipop_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialCoefLollipop(
                    signature = values$trainingResult$signature,
                    annot = values$annot_mode2
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadStabilityBarsMode2 <- downloadHandler(
        filename = "stability_bars_plot.png",
        content = function(file) {
            # Create a signature object that includes stability info if available
            sig_for_plot <- values$trainingResult$signature

            # Check if bootstrap_frequency exists in trainingResult$stability
            has_bootstrap <- !is.null(values$trainingResult$stability) &&
                !is.null(values$trainingResult$stability$bootstrap_frequency)

            # If bootstrap exists, add it to signature$stability for the function
            if (has_bootstrap) {
                sig_for_plot$stability <- list(
                    bootstrap_frequency = values$trainingResult$stability$bootstrap_frequency
                )
                freq_type <- "bootstrap"
            } else {
                # Use selection frequency (which is already in signature$selection_frequency)
                freq_type <- "selection"
            }

            # Try to create the plot with error handling
            p <- tryCatch(
                {
                    # Suppress warnings about bootstrap_frequency not available
                    suppressWarnings({
                        endoSignatureR::plotEndometrialStabilityBars(
                            signature = sig_for_plot,
                            annot = values$annot_mode2,
                            frequency_type = freq_type
                        )
                    })
                },
                error = function(e) {
                    # If no valid frequencies found, show a message plot instead
                    if (grepl("No valid frequencies found", e$message, fixed = TRUE)) {
                        ggplot2::ggplot() +
                            ggplot2::annotate("text",
                                x = 0.5, y = 0.5,
                                label = "Stability frequency data not available.\nAll frequency values are missing or invalid. Note: Stability Selection is required to generate the stability bars plot.",
                                size = 5, hjust = 0.5, vjust = 0.5
                            ) +
                            ggplot2::theme_void() +
                            ggplot2::labs(title = "Stability Bars Plot")
                    } else if (grepl("bootstrap_frequency not available", e$message, fixed = TRUE)) {
                        # Handle case where bootstrap frequency is requested but not available
                        ggplot2::ggplot() +
                            ggplot2::annotate("text",
                                x = 0.5, y = 0.5,
                                label = "Bootstrap frequency data not available.\nUsing selection frequency instead.\nNote: Small sample sizes (< 10 samples) prevent bootstrap resampling.",
                                size = 4, hjust = 0.5, vjust = 0.5
                            ) +
                            ggplot2::theme_void() +
                            ggplot2::labs(title = "Stability Bars Plot")
                    } else {
                        # For other errors, show error message
                        ggplot2::ggplot() +
                            ggplot2::annotate("text",
                                x = 0.5, y = 0.5,
                                label = paste("Error creating stability plot:\n", e$message),
                                size = 4, hjust = 0.5, vjust = 0.5
                            ) +
                            ggplot2::theme_void() +
                            ggplot2::labs(title = "Stability Bars Plot - Error")
                    }
                }
            )

            ggplot2::ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
        }
    )

    # Export handlers for Mode 2
    output$downloadSignatureCSVMode2 <- downloadHandler(
        filename = "endometrial_signature.csv",
        content = function(file) {
            if (!is.null(values$trainingResult)) {
                # Export only CSV format
                endoSignatureR::esr_exportSignature(
                    signature = values$trainingResult$signature,
                    dir = tempdir(),
                    result = values$trainingResult,
                    formats = "csv"
                )
                # Copy the generated file
                csv_file <- file.path(tempdir(), "endometrial_signature.csv")
                if (file.exists(csv_file)) {
                    file.copy(csv_file, file, overwrite = TRUE)
                } else {
                    stop("CSV file was not created. Please try again.")
                }
            } else {
                stop("No training result available. Please train a signature first.")
            }
        }
    )

    output$downloadSignatureJSONMode2 <- downloadHandler(
        filename = "endometrial_recipe.json",
        content = function(file) {
            if (!is.null(values$trainingResult)) {
                # Export only JSON format
                endoSignatureR::esr_exportSignature(
                    signature = values$trainingResult$signature,
                    dir = tempdir(),
                    result = values$trainingResult,
                    formats = "json"
                )
                # Copy the generated file (export function creates endometrial_recipe.json)
                json_file <- file.path(tempdir(), "endometrial_recipe.json")
                if (file.exists(json_file)) {
                    file.copy(json_file, file, overwrite = TRUE)
                } else {
                    stop("JSON file was not created. Please try again.")
                }
            } else {
                stop("No training result available. Please train a signature first.")
            }
        }
    )

    output$downloadModelCardMode2 <- downloadHandler(
        filename = "endometrial_model_card.md",
        content = function(file) {
            if (!is.null(values$trainingResult)) {
                # Export only Markdown format
                endoSignatureR::esr_exportSignature(
                    signature = values$trainingResult$signature,
                    dir = tempdir(),
                    result = values$trainingResult,
                    formats = "md"
                )
                # Copy the generated file
                md_file <- file.path(tempdir(), "endometrial_model_card.md")
                if (file.exists(md_file)) {
                    file.copy(md_file, file, overwrite = TRUE)
                } else {
                    stop("Model card file was not created. Please try again.")
                }
            } else {
                stop("No training result available. Please train a signature first.")
            }
        }
    )

    output$downloadStabilityCSVMode2 <- downloadHandler(
        filename = "endometrial_stability.csv",
        content = function(file) {
            if (is.null(values$trainingResult)) {
                # Write a CSV with error message instead of stopping
                error_df <- data.frame(
                    note = "No training result available. Please train a signature first.",
                    stringsAsFactors = FALSE
                )
                readr::write_csv(error_df, file)
                return()
            }

            # Check if bootstrap frequency exists
            has_bootstrap <- !is.null(values$trainingResult$stability) &&
                !is.null(values$trainingResult$stability$bootstrap_frequency) &&
                length(values$trainingResult$stability$bootstrap_frequency) > 0

            if (!has_bootstrap) {
                # Write a CSV with message explaining why data is not available
                message_df <- data.frame(
                    note = "Stability frequency data is not available. Bootstrap resampling requires at least 10 samples. The stability CSV will only be generated when bootstrap frequencies are computed during training.",
                    stringsAsFactors = FALSE
                )
                readr::write_csv(message_df, file)
                return()
            }

            # Extract bootstrap frequencies
            bootstrap_freq <- values$trainingResult$stability$bootstrap_frequency

            # Create data frame
            stability_df <- data.frame(
                gene_id = names(bootstrap_freq),
                bootstrap_frequency = as.numeric(bootstrap_freq),
                stringsAsFactors = FALSE
            )

            # Sort by frequency (descending)
            stability_df <- stability_df[order(stability_df$bootstrap_frequency, decreasing = TRUE), ]

            # Write CSV
            readr::write_csv(stability_df, file)
        }
    )

    output$downloadPredictionsMode2 <- downloadHandler(
        filename = "predictions.csv",
        content = function(file) {
            if (!is.null(values$trainingResult) && !is.null(values$trainingResult$metrics$predictions)) {
                readr::write_csv(values$trainingResult$metrics$predictions, file)
            }
        }
    )

    output$downloadFullResultsMode2 <- downloadHandler(
        filename = "training_results.rds",
        content = function(file) {
            if (!is.null(values$trainingResult)) {
                saveRDS(values$trainingResult, file)
            }
        }
    )

    # ===== MODE 3: VISUALIZATION & ANALYSIS =====

    # Load Mode 3 demo data
    observeEvent(input$loadDemoMode3, {
        tryCatch(
            {
                demo_data <- endoSignatureR::endo_load_demo()
                values$counts <- demo_data$counts
                values$pheno <- demo_data$pheno
                values$annot <- demo_data$annot
                values$dataLoaded <- TRUE

                # Transform counts
                values$counts_t <- endoSignatureR::esr_transform_log1p_cpm(values$counts)

                # Compute QC metrics
                values$qc_metrics <- endoSignatureR::esr_computeQCMetrics(
                    counts = values$counts,
                    mat_t = values$counts_t,
                    pheno = values$pheno
                )

                output$dataStatusMode3 <- renderText({
                    paste(
                        "Demo data loaded successfully!\n",
                        "Samples:", ncol(values$counts), "\n",
                        "Genes:", nrow(values$counts)
                    )
                })

                output$statusMessageMode3 <- renderText("Demo data loaded successfully!")
            },
            error = function(e) {
                output$statusMessageMode3 <- renderText(paste("Error loading demo data:", e$message))
            }
        )
    })

    # Download Mode 3 demo data
    output$downloadDemoMode3 <- downloadHandler(
        filename = function() {
            "gse201926_sample_demo.zip"
        },
        content = function(file) {
            temp_dir <- tempdir()
            demo_data <- endoSignatureR::endo_load_demo()

            # Convert counts matrix to data.frame with GeneID column
            # This matches the format expected by esr_loadCountsFromFile
            counts_df <- data.frame(
                GeneID = rownames(demo_data$counts),
                as.data.frame(demo_data$counts),
                stringsAsFactors = FALSE,
                check.names = FALSE
            )
            readr::write_tsv(
                counts_df,
                file.path(temp_dir, "gse201926_sample_counts.tsv")
            )
            readr::write_tsv(
                demo_data$pheno,
                file.path(temp_dir, "gse201926_sample_pheno.tsv")
            )
            readr::write_tsv(
                demo_data$annot,
                file.path(temp_dir, "gse201926_sample_annot.tsv")
            )

            utils::zip(file,
                files = c(
                    file.path(temp_dir, "gse201926_sample_counts.tsv"),
                    file.path(temp_dir, "gse201926_sample_pheno.tsv"),
                    file.path(temp_dir, "gse201926_sample_annot.tsv")
                ),
                flags = "-j"
            )
        }
    )

    # Load Mode 3 uploaded data
    observeEvent(
        {
            input$countsFileMode3
            input$phenoFileMode3
            input$annotFileMode3
        },
        {
            if (!is.null(input$countsFileMode3) && !is.null(input$phenoFileMode3)) {
                tryCatch(
                    {
                        values$counts <- endoSignatureR::esr_loadCountsFromFile(input$countsFileMode3$datapath)
                        values$pheno <- endoSignatureR::esr_loadPhenoFromFile(input$phenoFileMode3$datapath)

                        if (!is.null(input$annotFileMode3)) {
                            values$annot <- endoSignatureR::esr_loadAnnotFromFile(input$annotFileMode3$datapath)
                        } else {
                            values$annot <- NULL
                        }

                        validation <- endoSignatureR::esr_validateEndometrial(
                            values$counts, values$pheno, values$annot
                        )

                        values$counts <- validation$X
                        values$pheno <- validation$pheno
                        values$annot <- validation$annot

                        values$counts_t <- endoSignatureR::esr_transform_log1p_cpm(values$counts)

                        values$qc_metrics <- endoSignatureR::esr_computeQCMetrics(
                            counts = values$counts,
                            mat_t = values$counts_t,
                            pheno = values$pheno
                        )

                        values$dataLoaded <- TRUE

                        output$dataStatusMode3 <- renderText({
                            paste(
                                "Data loaded successfully!\n",
                                "Samples:", ncol(values$counts), "\n",
                                "Genes:", nrow(values$counts), "\n",
                                if (nrow(validation$issues) > 0) {
                                    paste("\nWarnings:", nrow(validation$issues))
                                } else {
                                    "\nNo validation issues"
                                }
                            )
                        })

                        output$statusMessageMode3 <- renderText("Data loaded and validated successfully!")
                    },
                    error = function(e) {
                        output$statusMessageMode3 <- renderText(paste("Error loading data:", e$message))
                    }
                )
            }
        }
    )

    # Output flag for data loaded
    output$dataLoaded <- reactive({
        values$dataLoaded
    })
    outputOptions(output, "dataLoaded", suspendWhenHidden = FALSE)

    # Run DE analysis
    observeEvent(input$runDE, {
        if (!values$dataLoaded) {
            output$deStatusMessage <- renderText("Error: Please load data first!")
            return()
        }

        tryCatch(
            {
                output$deStatusMessage <- renderText("Running differential expression analysis...")

                # Suppress expected warnings about small sample sizes in groups
                values$de_table <- suppressWarnings(
                    endoSignatureR::esr_analyzeDifferentialExpression(
                        mat_t = values$counts_t,
                        pheno = values$pheno
                    )
                )

                # Select top genes (default 50 for initial selection)
                # Users can adjust number when generating heatmap
                values$selected_genes <- endoSignatureR::esr_selectTopGenes(
                    de_table = values$de_table,
                    n = 50,
                    by = "de"
                )

                values$bundle <- endoSignatureR::esr_createAnalysisBundle(
                    counts_t = values$counts_t,
                    de_table = values$de_table,
                    selected_genes = values$selected_genes,
                    qc_metrics = values$qc_metrics,
                    pheno = values$pheno,
                    annot = values$annot
                )

                values$deRun <- TRUE

                # Generate all plots automatically with default settings
                # Replace entire lists to trigger reactivity
                values$maPlotParams <- list(
                    fdr_threshold = 0.05,
                    log2fc_threshold = 1,
                    point_size = 1.5,
                    point_alpha = 0.6,
                    legend_position = "right",
                    theme = "bw"
                )
                values$maPlotGenerated <- TRUE

                values$volcanoPlotParams <- list(
                    fdr_threshold = 0.05,
                    log2fc_threshold = 1,
                    point_size = 1.5,
                    point_alpha = 0.6,
                    legend_position = "right",
                    theme = "bw"
                )
                values$volcanoPlotGenerated <- TRUE

                values$heatmapPlotParams <- list(
                    n_genes = 50,
                    scale = "row",
                    show_row_names = FALSE
                )
                values$heatmapPlotGenerated <- TRUE

                output$statusMessageMode3 <- renderText("Differential expression analysis completed!")
                output$deStatusMessage <- renderText({
                    paste(
                        "DE Analysis completed successfully!\n",
                        "Total genes analyzed:", nrow(values$de_table), "\n",
                        "Top", length(values$selected_genes), "genes selected for heatmap.\n",
                        "All plots generated with default settings."
                    )
                })
            },
            error = function(e) {
                output$deStatusMessage <- renderText(paste("Error in DE analysis:", e$message))
                output$statusMessageMode3 <- renderText(paste("Error in DE analysis:", e$message))
            }
        )
    })

    # QC/EDA Plots
    # Update Libsize Plot
    observeEvent(input$updateLibsize, {
        if (!values$dataLoaded) {
            return()
        }
        # Replace entire list to trigger reactivity
        values$libsizePlotParams <- list(
            point_size = input$libsizePointSize,
            point_alpha = input$libsizePointAlpha,
            bins = input$libsizeBins,
            theme = input$libsizeTheme,
            legend_position = input$libsizeLegendPos
        )
    })

    output$libsizePlot <- renderPlot({
        if (!values$dataLoaded) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialLibsize(
            values$counts, values$pheno,
            point_size = values$libsizePlotParams$point_size,
            point_alpha = values$libsizePlotParams$point_alpha,
            bins = values$libsizePlotParams$bins,
            legend_position = values$libsizePlotParams$legend_position,
            theme = values$libsizePlotParams$theme
        )
    })

    # Update Zeros Plot
    observeEvent(input$updateZeros, {
        if (!values$dataLoaded) {
            return()
        }
        # Replace entire list to trigger reactivity
        values$zerosPlotParams <- list(
            bins = input$zerosBins,
            theme = input$zerosTheme
        )
    })

    output$zerosPlot <- renderPlot({
        if (!values$dataLoaded) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialZeros(
            values$counts,
            by = "sample",
            bins = values$zerosPlotParams$bins,
            theme = values$zerosPlotParams$theme
        )
    })

    # Update PCA Plot
    observeEvent(input$updatePCA, {
        if (!values$dataLoaded) {
            return()
        }
        # Replace entire list to trigger reactivity
        values$pcaPlotParams <- list(
            point_size = input$pcaPointSize,
            point_alpha = input$pcaPointAlpha,
            theme = input$pcaTheme,
            legend_position = input$pcaLegendPos
        )
    })

    output$pcaPlot <- renderPlot({
        if (!values$dataLoaded) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialPCA(
            values$counts_t, values$pheno,
            point_size = values$pcaPlotParams$point_size,
            point_alpha = values$pcaPlotParams$point_alpha,
            legend_position = values$pcaPlotParams$legend_position,
            theme = values$pcaPlotParams$theme
        )
    })

    # Generate MA Plot (regenerate with custom settings)
    observeEvent(input$generateMA, {
        if (!values$deRun || is.null(values$de_table)) {
            return()
        }
        # Replace entire list to trigger reactivity
        values$maPlotParams <- list(
            fdr_threshold = input$maFDR,
            log2fc_threshold = input$maLog2FC,
            point_size = input$maPointSize,
            point_alpha = input$maPointAlpha,
            legend_position = input$maLegendPos,
            theme = input$maTheme
        )
        values$maPlotGenerated <- TRUE
    })

    output$maPlot <- renderPlot({
        if (!values$maPlotGenerated || is.null(values$de_table)) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialMA(
            values$de_table,
            fdr_threshold = values$maPlotParams$fdr_threshold,
            log2fc_threshold = values$maPlotParams$log2fc_threshold,
            point_size = values$maPlotParams$point_size,
            point_alpha = values$maPlotParams$point_alpha,
            legend_position = values$maPlotParams$legend_position,
            theme = values$maPlotParams$theme
        )
    })

    # Generate Volcano Plot (regenerate with custom settings)
    observeEvent(input$generateVolcano, {
        if (!values$deRun || is.null(values$de_table)) {
            return()
        }
        # Replace entire list to trigger reactivity
        values$volcanoPlotParams <- list(
            fdr_threshold = input$volcanoFDR,
            log2fc_threshold = input$volcanoLog2FC,
            point_size = input$volcanoPointSize,
            point_alpha = input$volcanoPointAlpha,
            legend_position = input$volcanoLegendPos,
            theme = input$volcanoTheme
        )
        values$volcanoPlotGenerated <- TRUE
    })

    output$volcanoPlot <- renderPlot({
        if (!values$volcanoPlotGenerated || is.null(values$de_table)) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialVolcano(
            values$de_table,
            fdr_threshold = values$volcanoPlotParams$fdr_threshold,
            log2fc_threshold = values$volcanoPlotParams$log2fc_threshold,
            point_size = values$volcanoPlotParams$point_size,
            point_alpha = values$volcanoPlotParams$point_alpha,
            legend_position = values$volcanoPlotParams$legend_position,
            theme = values$volcanoPlotParams$theme
        )
    })

    # Generate Heatmap (regenerate with custom settings)
    observeEvent(input$generateHeatmap, {
        if (!values$deRun || is.null(values$de_table)) {
            return()
        }
        # Update selected genes based on new number
        values$selected_genes <- endoSignatureR::esr_selectTopGenes(
            de_table = values$de_table,
            n = input$heatmapGenes,
            by = "de"
        )
        # Replace entire list to trigger reactivity
        values$heatmapPlotParams <- list(
            n_genes = input$heatmapGenes,
            scale = input$heatmapScale,
            show_row_names = input$heatmapShowRowNames
        )
        values$heatmapPlotGenerated <- TRUE
    })

    output$heatmapPlot <- renderPlot({
        if (!values$heatmapPlotGenerated || is.null(values$selected_genes)) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialHeatmap(
            mat_t = values$counts_t,
            genes = values$selected_genes,
            pheno = values$pheno,
            scale = values$heatmapPlotParams$scale,
            show_row_names = values$heatmapPlotParams$show_row_names
        )
    })

    # DE Table
    output$deTable <- renderTable(
        {
            if (!values$deRun || is.null(values$de_table)) {
                return(NULL)
            }
            head(values$de_table, 50)
        },
        digits = 4
    )

    # Signature Plots (shown in Mode 1 Results tab after classification)
    output$coefLollipopPlot <- renderPlot({
        if (is.null(values$selectedSignature) || !values$classificationRun) {
            return(NULL)
        }
        # Use appropriate annotation based on signature type
        annot_to_use <- if (input$signatureSelection == "pretrained") {
            values$annot_mode1
        } else if (input$signatureSelection == "user-trained") {
            values$annot_mode2
        } else {
            # For uploaded signature, try annot_mode1 first, then annot_mode2, then NULL
            if (!is.null(values$annot_mode1)) {
                values$annot_mode1
            } else if (!is.null(values$annot_mode2)) {
                values$annot_mode2
            } else {
                NULL
            }
        }
        endoSignatureR::plotEndometrialCoefLollipop(
            signature = values$selectedSignature,
            annot = annot_to_use
        )
    })

    output$stabilityBarsPlot <- renderPlot({
        if (is.null(values$selectedSignature) || !values$classificationRun) {
            return(NULL)
        }

        # For user-trained signature, check trainingResult$stability structure
        if (input$signatureSelection == "user-trained" && !is.null(values$trainingResult)) {
            # User-trained signature: check both signature$selection_frequency and trainingResult$stability
            has_stability <- !is.null(values$trainingResult$stability)
            has_selection_freq <- !is.null(values$selectedSignature$selection_frequency)
            if (!has_stability && !has_selection_freq) {
                return(NULL)
            }
        } else {
            # Pre-trained signature: check signature$stability and signature$selection_frequency
            if (is.null(values$selectedSignature$stability) &&
                is.null(values$selectedSignature$selection_frequency)) {
                return(NULL)
            }
        }

        # Use appropriate annotation based on signature type
        annot_to_use <- if (input$signatureSelection == "pretrained") {
            values$annot_mode1
        } else if (input$signatureSelection == "user-trained") {
            values$annot_mode2
        } else {
            # For uploaded signature, try annot_mode1 first, then annot_mode2, then NULL
            if (!is.null(values$annot_mode1)) {
                values$annot_mode1
            } else if (!is.null(values$annot_mode2)) {
                values$annot_mode2
            } else {
                NULL
            }
        }

        # Prepare signature object for plotting
        sig_for_plot <- values$selectedSignature

        # For user-trained signature, add stability info if available
        if (input$signatureSelection == "user-trained" && !is.null(values$trainingResult$stability)) {
            if (!is.null(values$trainingResult$stability$bootstrap_frequency)) {
                sig_for_plot$stability <- list(
                    bootstrap_frequency = values$trainingResult$stability$bootstrap_frequency
                )
            }
        }
        # For uploaded signature, stability info is already in the signature object

        # Check which frequency type is available
        has_bootstrap <- !is.null(sig_for_plot$stability) &&
            !is.null(sig_for_plot$stability$bootstrap_frequency)
        has_selection <- !is.null(sig_for_plot$selection_frequency)

        # Use selection frequency if bootstrap not available
        # Default to selection if bootstrap not available (has_selection should be TRUE if we got here)
        freq_type <- if (has_bootstrap) "bootstrap" else if (has_selection) "selection" else "selection"

        # Try to create the plot with error handling
        tryCatch(
            {
                # Suppress warnings about bootstrap_frequency not available
                suppressWarnings({
                    endoSignatureR::plotEndometrialStabilityBars(
                        signature = sig_for_plot,
                        annot = annot_to_use,
                        frequency_type = freq_type
                    )
                })
            },
            error = function(e) {
                # If no valid frequencies found, show a message plot instead
                if (grepl("No valid frequencies found", e$message, fixed = TRUE)) {
                    ggplot2::ggplot() +
                        ggplot2::annotate("text",
                            x = 0.5, y = 0.5,
                            label = "Stability frequency data not available.\nAll frequency values are missing or invalid.",
                            size = 5, hjust = 0.5, vjust = 0.5
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Stability Bars Plot")
                } else if (grepl("bootstrap_frequency not available", e$message, fixed = TRUE)) {
                    # Handle case where bootstrap frequency is requested but not available
                    ggplot2::ggplot() +
                        ggplot2::annotate("text",
                            x = 0.5, y = 0.5,
                            label = "Bootstrap frequency data not available.\nUsing selection frequency instead.\nNote: Small sample sizes (< 10 samples) prevent bootstrap resampling.",
                            size = 4, hjust = 0.5, vjust = 0.5
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Stability Bars Plot")
                } else {
                    # For other errors, show error message
                    ggplot2::ggplot() +
                        ggplot2::annotate("text",
                            x = 0.5, y = 0.5,
                            label = paste("Error creating stability plot:\n", e$message),
                            size = 4, hjust = 0.5, vjust = 0.5
                        ) +
                        ggplot2::theme_void() +
                        ggplot2::labs(title = "Stability Bars Plot - Error")
                }
            }
        )
    })

    # Download handlers for plots
    output$downloadLibsize <- downloadHandler(
        filename = "libsize_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialLibsize(
                    values$counts, values$pheno,
                    point_size = values$libsizePlotParams$point_size,
                    point_alpha = values$libsizePlotParams$point_alpha,
                    bins = values$libsizePlotParams$bins,
                    legend_position = values$libsizePlotParams$legend_position,
                    theme = values$libsizePlotParams$theme
                ),
                width = 10, height = 6, dpi = 300
            )
        }
    )

    output$downloadZeros <- downloadHandler(
        filename = "zeros_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialZeros(
                    values$counts,
                    by = "sample",
                    bins = values$zerosPlotParams$bins,
                    theme = values$zerosPlotParams$theme
                ),
                width = 10, height = 6, dpi = 300
            )
        }
    )

    output$downloadPCA <- downloadHandler(
        filename = "pca_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialPCA(
                    values$counts_t, values$pheno,
                    point_size = values$pcaPlotParams$point_size,
                    point_alpha = values$pcaPlotParams$point_alpha,
                    legend_position = values$pcaPlotParams$legend_position,
                    theme = values$pcaPlotParams$theme
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadMA <- downloadHandler(
        filename = "ma_plot.png",
        content = function(file) {
            if (!values$maPlotGenerated || is.null(values$de_table)) {
                return()
            }
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialMA(
                    values$de_table,
                    fdr_threshold = values$maPlotParams$fdr_threshold,
                    log2fc_threshold = values$maPlotParams$log2fc_threshold,
                    point_size = values$maPlotParams$point_size,
                    point_alpha = values$maPlotParams$point_alpha,
                    legend_position = values$maPlotParams$legend_position,
                    theme = values$maPlotParams$theme
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadVolcano <- downloadHandler(
        filename = "volcano_plot.png",
        content = function(file) {
            if (!values$volcanoPlotGenerated || is.null(values$de_table)) {
                return()
            }
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialVolcano(
                    values$de_table,
                    fdr_threshold = values$volcanoPlotParams$fdr_threshold,
                    log2fc_threshold = values$volcanoPlotParams$log2fc_threshold,
                    point_size = values$volcanoPlotParams$point_size,
                    point_alpha = values$volcanoPlotParams$point_alpha,
                    legend_position = values$volcanoPlotParams$legend_position,
                    theme = values$volcanoPlotParams$theme
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadHeatmap <- downloadHandler(
        filename = "heatmap_plot.png",
        content = function(file) {
            if (!values$heatmapPlotGenerated || is.null(values$selected_genes)) {
                return()
            }
            png(file, width = 12, height = 10, units = "in", res = 300)
            ComplexHeatmap::draw(endoSignatureR::plotEndometrialHeatmap(
                mat_t = values$counts_t,
                genes = values$selected_genes,
                pheno = values$pheno,
                scale = values$heatmapPlotParams$scale,
                show_row_names = values$heatmapPlotParams$show_row_names
            ))
            dev.off()
        }
    )

    output$downloadCoefLollipop <- downloadHandler(
        filename = "coef_lollipop_plot.png",
        content = function(file) {
            # Use appropriate annotation based on signature type
            annot_to_use <- if (input$signatureSelection == "pretrained") {
                values$annot_mode1
            } else if (input$signatureSelection == "user-trained") {
                values$annot_mode2
            } else {
                # For uploaded signature, try annot_mode1 first, then annot_mode2, then NULL
                if (!is.null(values$annot_mode1)) {
                    values$annot_mode1
                } else if (!is.null(values$annot_mode2)) {
                    values$annot_mode2
                } else {
                    NULL
                }
            }
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialCoefLollipop(
                    signature = values$selectedSignature,
                    annot = annot_to_use
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadStabilityBars <- downloadHandler(
        filename = "stability_bars_plot.png",
        content = function(file) {
            # Use appropriate annotation based on signature type
            annot_to_use <- if (input$signatureSelection == "pretrained") {
                values$annot_mode1
            } else if (input$signatureSelection == "user-trained") {
                values$annot_mode2
            } else {
                # For uploaded signature, try annot_mode1 first, then annot_mode2, then NULL
                if (!is.null(values$annot_mode1)) {
                    values$annot_mode1
                } else if (!is.null(values$annot_mode2)) {
                    values$annot_mode2
                } else {
                    NULL
                }
            }

            # Prepare signature object for plotting
            sig_for_plot <- values$selectedSignature

            # For user-trained signature, add stability info if available
            if (input$signatureSelection == "user-trained" && !is.null(values$trainingResult$stability)) {
                if (!is.null(values$trainingResult$stability$bootstrap_frequency)) {
                    sig_for_plot$stability <- list(
                        bootstrap_frequency = values$trainingResult$stability$bootstrap_frequency
                    )
                }
            }
            # For uploaded signature, stability info is already in the signature object

            # Check which frequency type is available
            has_bootstrap <- !is.null(sig_for_plot$stability) &&
                !is.null(sig_for_plot$stability$bootstrap_frequency)
            has_selection <- !is.null(sig_for_plot$selection_frequency)

            # Use selection frequency if bootstrap not available
            # Default to selection if bootstrap not available (has_selection should be TRUE if we got here)
            freq_type <- if (has_bootstrap) "bootstrap" else if (has_selection) "selection" else "selection"

            # Try to create the plot with error handling
            p <- tryCatch(
                {
                    # Suppress warnings about bootstrap_frequency not available
                    suppressWarnings({
                        endoSignatureR::plotEndometrialStabilityBars(
                            signature = sig_for_plot,
                            annot = annot_to_use,
                            frequency_type = freq_type
                        )
                    })
                },
                error = function(e) {
                    # If no valid frequencies found, show a message plot instead
                    if (grepl("No valid frequencies found", e$message, fixed = TRUE)) {
                        ggplot2::ggplot() +
                            ggplot2::annotate("text",
                                x = 0.5, y = 0.5,
                                label = "Stability frequency data not available.\nAll frequency values are missing or invalid.",
                                size = 5, hjust = 0.5, vjust = 0.5
                            ) +
                            ggplot2::theme_void() +
                            ggplot2::labs(title = "Stability Bars Plot")
                    } else if (grepl("bootstrap_frequency not available", e$message, fixed = TRUE)) {
                        # Handle case where bootstrap frequency is requested but not available
                        ggplot2::ggplot() +
                            ggplot2::annotate("text",
                                x = 0.5, y = 0.5,
                                label = "Bootstrap frequency data not available.\nUsing selection frequency instead.\nNote: Small sample sizes (< 10 samples) prevent bootstrap resampling.",
                                size = 4, hjust = 0.5, vjust = 0.5
                            ) +
                            ggplot2::theme_void() +
                            ggplot2::labs(title = "Stability Bars Plot")
                    } else {
                        # For other errors, show error message
                        ggplot2::ggplot() +
                            ggplot2::annotate("text",
                                x = 0.5, y = 0.5,
                                label = paste("Error creating stability plot:\n", e$message),
                                size = 4, hjust = 0.5, vjust = 0.5
                            ) +
                            ggplot2::theme_void() +
                            ggplot2::labs(title = "Stability Bars Plot - Error")
                    }
                }
            )

            ggplot2::ggsave(file, plot = p, width = 10, height = 8, dpi = 300)
        }
    )

    # Export handlers
    output$downloadBundle <- downloadHandler(
        filename = "analysis_bundle.rds",
        content = function(file) {
            if (is.null(values$bundle)) {
                values$bundle <- endoSignatureR::esr_createAnalysisBundle(
                    counts_t = values$counts_t,
                    de_table = values$de_table,
                    selected_genes = values$selected_genes,
                    qc_metrics = values$qc_metrics,
                    pheno = values$pheno,
                    annot = values$annot
                )
            }
            saveRDS(values$bundle, file)
        }
    )

    output$downloadCountsT <- downloadHandler(
        filename = "transformed_counts.tsv",
        content = function(file) {
            readr::write_tsv(as.data.frame(values$counts_t), file)
        }
    )

    output$downloadDETable <- downloadHandler(
        filename = "de_table.tsv",
        content = function(file) {
            if (!is.null(values$de_table)) {
                readr::write_tsv(values$de_table, file)
            }
        }
    )

    output$downloadQCMetrics <- downloadHandler(
        filename = "qc_metrics.tsv",
        content = function(file) {
            if (!is.null(values$qc_metrics)) {
                readr::write_tsv(values$qc_metrics, file)
            }
        }
    )

    output$downloadSelectedGenes <- downloadHandler(
        filename = "selected_genes.tsv",
        content = function(file) {
            if (!is.null(values$selected_genes)) {
                readr::write_tsv(data.frame(gene_id = values$selected_genes), file)
            }
        }
    )
}

# Note: ui and server are defined above and will be used by runEndoSignatureR()
# The shinyApp() call is made in runEndoSignatureR() function

# [END]
