# Integration Tests for Mode 2 (Signature Validation) Full Pipeline
# Tests complete workflow: data → validate → train → evaluate → visualize → export

test_that("Mode 2 full pipeline: data → validate → train → evaluate → visualize → export", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))
    
    library(glmnet)
    library(rsample)
    
    # 1. Load data
    data(gse201926_trainmini)
    
    # 2. Validate data
    validation <- esr_validateEndometrial(
        gse201926_trainmini$counts,
        gse201926_trainmini$pheno,
        annot = gse201926_trainmini$annot
    )
    
    # 3. Train signature
    set.seed(123)
    result <- suppressWarnings(esr_trainEndometrialSignature(
        X = validation$X,
        pheno = validation$pheno,
        annot = validation$annot,
        top_k = 50,
        outer_folds = NULL, # Use folds_demo
        seed = 123
    ))
    
    # 4. Visualize (performance plots)
    # Extract predictions (already has label and prob columns)
    pred_df <- result$metrics$predictions
    
    p1 <- plotEndometrialROC(pred_df, use_calibrated = FALSE)
    p2 <- plotEndometrialPR(pred_df, use_calibrated = FALSE)
    p3 <- plotEndometrialCalibration(pred_df, use_calibrated = FALSE)
    p4 <- plotEndometrialCoefLollipop(result$signature)
    
    # 5. Export signature
    temp_dir <- tempfile()
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))
    
    export_paths <- esr_exportSignature(
        signature = result$signature,
        result = result,
        dir = temp_dir,
        formats = c("csv", "json", "md")
    )
    
    # 6. Compare with pre-trained (if available)
    csv_path <- system.file("extdata", "pretrained-signature", 
                            "endometrial_signature.csv", package = "endoSignatureR")
    if (nchar(csv_path) > 0 && file.exists(csv_path)) {
        pretrained <- esr_loadPretrainedSignature()
        # Use plotEndometrialComparison for visualization comparison
        comparison_plot <- plotEndometrialComparison(
            pretrained_result = NULL, # Can't use pretrained_result directly, need metrics
            new_result = result
        )
        
        # Note: esr_compareSignatures is placeholder, but test that it runs
        # comparison <- esr_compareSignatures(pretrained, result)
    }
    
    # Verify all steps completed successfully
    expect_true(is.list(result))
    expect_true("signature" %in% names(result))
    expect_true("metrics" %in% names(result))
    expect_true(inherits(p1, "ggplot"))
    expect_true(inherits(p2, "ggplot"))
    expect_true(inherits(p3, "ggplot"))
    expect_true(inherits(p4, "ggplot"))
    # Verify all export files exist
    expect_true(is.character(export_paths))
    expect_true(all(file.exists(export_paths)))
})

# [END]
