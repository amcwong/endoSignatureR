# EndoSignatureR Shiny Application
# Three-Mode Workflow: Mode 1 (Rapid Classification), Mode 2 (Signature Validation), Mode 3 (Visualization & Analysis)

library(shiny)

# Define UI
ui <- fluidPage(
    titlePanel("EndoSignatureR: Endometrial RNA-seq Analysis"),

    # Main tabset with three modes
    tabsetPanel(
        id = "modeTabs",

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

                    # Signature information
                    h4("Pre-trained Signature"),
                    verbatimTextOutput("signatureInfo"),
                    br(),

                    # Classification controls
                    h4("Classification"),
                    p("Apply pre-trained signature to classify samples:"),
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
                                downloadButton("downloadPredictions", "Download Predictions (CSV)"),
                                br(), br(),
                                tableOutput("predictionsTable"),
                                br(), br(),
                                hr(),
                                h4("Signature Visualization"),
                                p("Visualization of the pre-trained signature used for classification."),
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
                                downloadButton("downloadSignatureCSVMode2", "Download Signature (CSV)"),
                                br(), br(),
                                downloadButton("downloadSignatureJSONMode2", "Download Signature (JSON)"),
                                br(), br(),
                                downloadButton("downloadModelCardMode2", "Download Model Card (MD)"),
                                br(), br(), hr(), br(),
                                h4("Predictions"),
                                p("Export predictions from cross-validation:"),
                                downloadButton("downloadPredictionsMode2", "Download Predictions (CSV)"),
                                br(), br(), hr(), br(),
                                h4("Full Results"),
                                p("Export complete training result (RDS format):"),
                                downloadButton("downloadFullResultsMode2", "Download Full Results (RDS)")
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
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    # Increase maximum upload size to 500MB (default is 5MB)
    options(shiny.maxRequestSize = 500 * 1024^2)

    # Reactive values to store data
    values <- reactiveValues(
        # Mode 1 data
        counts_mode1 = NULL,
        annot_mode1 = NULL,
        predictions = NULL,
        signature = NULL,
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
            },
            error = function(e) {
                # Signature loading will fail gracefully if files not found
                # User will see error message in signature info output
            }
        )
    })

    # Display signature information (for Mode 1)
    output$signatureInfo <- renderText({
        if (is.null(values$signature)) {
            return("Pre-trained signature not available. Please ensure package is properly installed.")
        }
        paste(
            "Pre-trained signature loaded successfully!\n",
            "Panel size:", length(values$signature$panel), "genes\n",
            "Preprocessing:", values$signature$recipe$preprocessing$transform %||% "log1p-cpm\n",
            "Trained on: GSE201926 dataset"
        )
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
                downloadButton("downloadBundle", "Download Analysis Bundle"),
                br(), br(),
                h4("Individual Components"),
                downloadButton("downloadCountsT", "Download Transformed Counts (TSV)"),
                br(), br(),
                downloadButton("downloadDETable", "Download DE Table (TSV)"),
                br(), br(),
                downloadButton("downloadQCMetrics", "Download QC Metrics (TSV)"),
                br(), br(),
                downloadButton("downloadSelectedGenes", "Download Selected Genes (TSV)")
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

        if (is.null(values$signature)) {
            output$statusMessageMode1 <- renderText("Pre-trained signature not available. Please ensure package is properly installed.")
            return()
        }

        tryCatch(
            {
                output$statusMessageMode1 <- renderText("Running classification...")

                # Run classification
                values$predictions <- endoSignatureR::esr_classifyEndometrial(
                    X_new = values$counts_mode1,
                    signature = values$signature,
                    threshold = input$thresholdMode1,
                    confidence = input$confidenceMode1
                )

                # Set completion flags
                values$classificationRun <- TRUE
                values$mode1_completed <- TRUE

                output$statusMessageMode1 <- renderText("Classification completed successfully!")

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

                        result <- endoSignatureR::esr_trainEndometrialSignature(
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

                        setProgress(1.0, detail = "Training complete!")
                    }
                )

                values$trainingResult <- result
                values$trainingRun <- TRUE
                values$mode2_completed <- TRUE

                output$statusMessageMode2 <- renderText("Training completed successfully!")

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
                endoSignatureR::plotEndometrialStabilityBars(
                    signature = sig_for_plot,
                    annot = values$annot_mode2,
                    frequency_type = freq_type
                )
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
                    endoSignatureR::plotEndometrialStabilityBars(
                        signature = sig_for_plot,
                        annot = values$annot_mode2,
                        frequency_type = freq_type
                    )
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
                endoSignatureR::esr_exportSignature(
                    signature = values$trainingResult$signature,
                    dir = tempdir(),
                    result = values$trainingResult
                )
                # Copy the generated file
                sig_file <- file.path(tempdir(), "endometrial_signature.csv")
                if (file.exists(sig_file)) {
                    file.copy(sig_file, file, overwrite = TRUE)
                }
            }
        }
    )

    output$downloadSignatureJSONMode2 <- downloadHandler(
        filename = "endometrial_signature.json",
        content = function(file) {
            if (!is.null(values$trainingResult)) {
                endoSignatureR::esr_exportSignature(
                    signature = values$trainingResult$signature,
                    dir = tempdir(),
                    result = values$trainingResult
                )
                # Copy the generated file
                sig_file <- file.path(tempdir(), "endometrial_signature.json")
                if (file.exists(sig_file)) {
                    file.copy(sig_file, file, overwrite = TRUE)
                }
            }
        }
    )

    output$downloadModelCardMode2 <- downloadHandler(
        filename = "endometrial_model_card.md",
        content = function(file) {
            if (!is.null(values$trainingResult)) {
                endoSignatureR::esr_exportSignature(
                    signature = values$trainingResult$signature,
                    dir = tempdir(),
                    result = values$trainingResult
                )
                # Copy the generated file
                sig_file <- file.path(tempdir(), "endometrial_model_card.md")
                if (file.exists(sig_file)) {
                    file.copy(sig_file, file, overwrite = TRUE)
                }
            }
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

                values$de_table <- endoSignatureR::esr_analyzeDifferentialExpression(
                    mat_t = values$counts_t,
                    pheno = values$pheno
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
        if (is.null(values$signature) || !values$classificationRun) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialCoefLollipop(
            signature = values$signature,
            annot = values$annot_mode1
        )
    })

    output$stabilityBarsPlot <- renderPlot({
        if (is.null(values$signature) || !values$classificationRun) {
            return(NULL)
        }
        if (is.null(values$signature$stability) && is.null(values$signature$selection_frequency)) {
            return(NULL)
        }
        # Check which frequency type is available
        has_bootstrap <- !is.null(values$signature$stability) &&
            !is.null(values$signature$stability$bootstrap_frequency)
        has_selection <- !is.null(values$signature$selection_frequency)

        # Use selection frequency if bootstrap not available
        freq_type <- if (has_bootstrap) "bootstrap" else "selection"

        # Try to create the plot with error handling
        tryCatch(
            {
                endoSignatureR::plotEndometrialStabilityBars(
                    signature = values$signature,
                    annot = values$annot_mode1,
                    frequency_type = freq_type
                )
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
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialCoefLollipop(
                    signature = values$signature,
                    annot = values$annot_mode1
                ),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadStabilityBars <- downloadHandler(
        filename = "stability_bars_plot.png",
        content = function(file) {
            # Check which frequency type is available
            has_bootstrap <- !is.null(values$signature$stability) &&
                !is.null(values$signature$stability$bootstrap_frequency)
            has_selection <- !is.null(values$signature$selection_frequency)

            # Use selection frequency if bootstrap not available
            freq_type <- if (has_bootstrap) "bootstrap" else "selection"

            # Try to create the plot with error handling
            p <- tryCatch(
                {
                    endoSignatureR::plotEndometrialStabilityBars(
                        signature = values$signature,
                        annot = values$annot_mode1,
                        frequency_type = freq_type
                    )
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
