# EndoSignatureR Shiny Application
# Mode 3: Standalone Visualization & Analysis

library(shiny)

# Define UI
ui <- fluidPage(
    titlePanel("EndoSignatureR: Endometrial RNA-seq Analysis"),

    # Sidebar layout
    sidebarLayout(
        sidebarPanel(
            h3("Mode 3: Visualization & Analysis"),
            p("Upload endometrial RNA-seq data or use demo data to perform quality control,
        exploratory analysis, and differential expression visualization."),
            br(),

            # Demo data section
            h4("Demo Data"),
            p("Load or download the bundled demo dataset (gse201926_sample):"),
            actionButton("loadDemo", "Load Demo Data", class = "btn-primary"),
            br(), br(),
            downloadButton("downloadDemo", "Download Demo Data", class = "btn-secondary"),
            br(), br(),
            hr(),

            # Data upload section
            h4("Upload Your Data"),
            p("Upload your own endometrial RNA-seq data files:"),
            fileInput("countsFile", "Counts Matrix (TSV/CSV)",
                accept = c(".tsv", ".csv", ".txt"),
                placeholder = "Genes × Samples"
            ),
            fileInput("phenoFile", "Phenotype Metadata (TSV/CSV)",
                accept = c(".tsv", ".csv", ".txt"),
                placeholder = "sample_id, group columns"
            ),
            fileInput("annotFile", "Annotation Data (TSV/CSV, optional)",
                accept = c(".tsv", ".csv", ".txt"),
                placeholder = "Gene annotations"
            ),
            br(),

            # Data status
            verbatimTextOutput("dataStatus"),
            br(),

            # DE analysis controls
            h4("Differential Expression"),
            p("Run differential expression analysis:"),
            numericInput("topGenes", "Number of top genes for heatmap",
                value = 50, min = 10, max = 500, step = 10
            ),
            p("Genes are selected by FDR (ascending), then by absolute log2FC (descending)."),
            actionButton("runDE", "Run DE Analysis", class = "btn-success"),
            br(), br(),

            # Status messages
            verbatimTextOutput("statusMessage")
        ),
        mainPanel(
            tabsetPanel(
                id = "mainTabs",

                # QC/EDA Tab
                tabPanel(
                    "QC/EDA",
                    h3("Quality Control and Exploratory Data Analysis"),
                    br(),
                    h4("Library Size"),
                    plotOutput("libsizePlot", height = "400px"),
                    downloadButton("downloadLibsize", "Download Plot"),
                    br(), br(),
                    h4("Percentage of Zeros"),
                    plotOutput("zerosPlot", height = "400px"),
                    downloadButton("downloadZeros", "Download Plot"),
                    br(), br(),
                    h4("Principal Component Analysis (PCA)"),
                    plotOutput("pcaPlot", height = "500px"),
                    downloadButton("downloadPCA", "Download Plot")
                ),

                # DE Tab
                tabPanel(
                    "Differential Expression",
                    h3("Differential Expression Analysis"),
                    br(),
                    conditionalPanel(
                        condition = "output.deRun == false",
                        p("Click 'Run DE Analysis' in the sidebar to perform differential expression analysis.")
                    ),
                    conditionalPanel(
                        condition = "output.deRun == true",
                        h4("MA Plot"),
                        plotOutput("maPlot", height = "500px"),
                        downloadButton("downloadMA", "Download Plot"),
                        br(), br(),
                        h4("Volcano Plot"),
                        plotOutput("volcanoPlot", height = "500px"),
                        downloadButton("downloadVolcano", "Download Plot"),
                        br(), br(),
                        h4("Heatmap"),
                        plotOutput("heatmapPlot", height = "600px"),
                        downloadButton("downloadHeatmap", "Download Plot"),
                        br(), br(),
                        h4("DE Results Table"),
                        p("Showing top 50 genes. Use export to download full table."),
                        tableOutput("deTable")
                    )
                ),

                # Export Tab
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
            )
        )
    )
)

# Define server logic
server <- function(input, output, session) {
    # Reactive values to store data
    values <- reactiveValues(
        counts = NULL,
        pheno = NULL,
        annot = NULL,
        counts_t = NULL,
        qc_metrics = NULL,
        de_table = NULL,
        selected_genes = NULL,
        bundle = NULL,
        dataLoaded = FALSE,
        deRun = FALSE
    )

    # Output flag for DE tab conditional panel
    output$deRun <- reactive({
        values$deRun
    })
    outputOptions(output, "deRun", suspendWhenHidden = FALSE)

    # Load demo data
    observeEvent(input$loadDemo, {
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

                output$dataStatus <- renderText({
                    paste(
                        "Demo data loaded successfully!\n",
                        "Samples:", ncol(values$counts), "\n",
                        "Genes:", nrow(values$counts)
                    )
                })

                output$statusMessage <- renderText("Demo data loaded successfully!")
            },
            error = function(e) {
                output$statusMessage <- renderText(paste("Error loading demo data:", e$message))
            }
        )
    })

    # Download demo data
    output$downloadDemo <- downloadHandler(
        filename = function() {
            "gse201926_sample_demo.zip"
        },
        content = function(file) {
            # Create temporary directory
            temp_dir <- tempdir()

            # Load demo data
            demo_data <- endoSignatureR::endo_load_demo()

            # Write files
            readr::write_tsv(
                as.data.frame(demo_data$counts),
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

            # Create zip file
            utils::zip(file,
                files = c(
                    file.path(temp_dir, "gse201926_sample_counts.tsv"),
                    file.path(temp_dir, "gse201926_sample_pheno.tsv"),
                    file.path(temp_dir, "gse201926_sample_annot.tsv")
                ),
                flags = "-j"
            ) # -j stores just the file names, not the paths
        }
    )

    # Load uploaded data
    observeEvent(
        {
            input$countsFile
            input$phenoFile
            input$annotFile
        },
        {
            if (!is.null(input$countsFile) && !is.null(input$phenoFile)) {
                tryCatch(
                    {
                        # Load counts
                        values$counts <- endoSignatureR::esr_loadCountsFromFile(input$countsFile$datapath)

                        # Load pheno
                        values$pheno <- endoSignatureR::esr_loadPhenoFromFile(input$phenoFile$datapath)

                        # Load annot if provided
                        if (!is.null(input$annotFile)) {
                            values$annot <- endoSignatureR::esr_loadAnnotFromFile(input$annotFile$datapath)
                        } else {
                            values$annot <- NULL
                        }

                        # Validate data
                        validation <- endoSignatureR::esr_validateEndometrial(
                            values$counts, values$pheno, values$annot
                        )

                        values$counts <- validation$X
                        values$pheno <- validation$pheno
                        values$annot <- validation$annot

                        # Transform counts
                        values$counts_t <- endoSignatureR::esr_transform_log1p_cpm(values$counts)

                        # Compute QC metrics
                        values$qc_metrics <- endoSignatureR::esr_computeQCMetrics(
                            counts = values$counts,
                            mat_t = values$counts_t,
                            pheno = values$pheno
                        )

                        values$dataLoaded <- TRUE

                        output$dataStatus <- renderText({
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

                        output$statusMessage <- renderText("Data loaded and validated successfully!")
                    },
                    error = function(e) {
                        output$statusMessage <- renderText(paste("Error loading data:", e$message))
                    }
                )
            }
        }
    )

    # Run DE analysis
    observeEvent(input$runDE, {
        if (!values$dataLoaded) {
            output$statusMessage <- renderText("Please load data first!")
            return()
        }

        tryCatch(
            {
                output$statusMessage <- renderText("Running differential expression analysis...")

                # Run DE analysis
                values$de_table <- endoSignatureR::esr_analyzeDifferentialExpression(
                    mat_t = values$counts_t,
                    pheno = values$pheno
                )

                # Select top genes (sorts by FDR then abs_log2FC)
                values$selected_genes <- endoSignatureR::esr_selectTopGenes(
                    de_table = values$de_table,
                    n = input$topGenes,
                    by = "de"
                )

                # Create analysis bundle
                values$bundle <- endoSignatureR::esr_createAnalysisBundle(
                    counts_t = values$counts_t,
                    de_table = values$de_table,
                    selected_genes = values$selected_genes,
                    qc_metrics = values$qc_metrics,
                    pheno = values$pheno,
                    annot = values$annot
                )

                values$deRun <- TRUE
                output$statusMessage <- renderText("Differential expression analysis completed!")
            },
            error = function(e) {
                output$statusMessage <- renderText(paste("Error in DE analysis:", e$message))
            }
        )
    })

    # QC/EDA Plots
    output$libsizePlot <- renderPlot({
        if (!values$dataLoaded) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialLibsize(values$counts, values$pheno)
    })

    output$zerosPlot <- renderPlot({
        if (!values$dataLoaded) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialZeros(values$counts, by = "sample")
    })

    output$pcaPlot <- renderPlot({
        if (!values$dataLoaded) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialPCA(values$counts_t, values$pheno)
    })

    # DE Plots
    output$maPlot <- renderPlot({
        if (!values$deRun || is.null(values$de_table)) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialMA(values$de_table)
    })

    output$volcanoPlot <- renderPlot({
        if (!values$deRun || is.null(values$de_table)) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialVolcano(values$de_table)
    })

    output$heatmapPlot <- renderPlot({
        if (!values$deRun || is.null(values$selected_genes)) {
            return(NULL)
        }
        endoSignatureR::plotEndometrialHeatmap(
            mat_t = values$counts_t,
            genes = values$selected_genes,
            pheno = values$pheno
        )
    })

    # DE Table
    output$deTable <- renderTable(
        {
            if (!values$deRun || is.null(values$de_table)) {
                return(NULL)
            }
            # Show top 50 rows
            head(values$de_table, 50)
        },
        digits = 4
    )

    # Download handlers for plots
    output$downloadLibsize <- downloadHandler(
        filename = "libsize_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialLibsize(values$counts, values$pheno),
                width = 10, height = 6, dpi = 300
            )
        }
    )

    output$downloadZeros <- downloadHandler(
        filename = "zeros_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialZeros(values$counts, by = "sample"),
                width = 10, height = 6, dpi = 300
            )
        }
    )

    output$downloadPCA <- downloadHandler(
        filename = "pca_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialPCA(values$counts_t, values$pheno),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadMA <- downloadHandler(
        filename = "ma_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialMA(values$de_table),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadVolcano <- downloadHandler(
        filename = "volcano_plot.png",
        content = function(file) {
            ggplot2::ggsave(file,
                plot = endoSignatureR::plotEndometrialVolcano(values$de_table),
                width = 10, height = 8, dpi = 300
            )
        }
    )

    output$downloadHeatmap <- downloadHandler(
        filename = "heatmap_plot.png",
        content = function(file) {
            png(file, width = 12, height = 10, units = "in", res = 300)
            ComplexHeatmap::draw(endoSignatureR::plotEndometrialHeatmap(
                mat_t = values$counts_t,
                genes = values$selected_genes,
                pheno = values$pheno
            ))
            dev.off()
        }
    )

    # Export handlers
    output$downloadBundle <- downloadHandler(
        filename = "analysis_bundle.rds",
        content = function(file) {
            if (is.null(values$bundle)) {
                # Create bundle if not exists
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
