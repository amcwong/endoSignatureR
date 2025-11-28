# Integration Tests for Mode 1 (Rapid Classification) Full Pipeline
# Tests complete workflow: data → signature → classify → visualize → export

test_that("Mode 1 full pipeline: data → signature → classify → visualize → export", {
    skip_if_not_installed("glmnet")
    
    # Skip if pretrained artifacts don't exist
    csv_path <- system.file("extdata", "pretrained-signature", 
                            "endometrial_signature.csv", package = "endoSignatureR")
    skip_if(
        nchar(csv_path) == 0 || !file.exists(csv_path),
        "Pretrained signature artifacts not found"
    )
    
    # 1. Load data (simulating data loading step)
    data(gse201926_sample)
    X_new <- gse201926_sample$counts
    
    # 2. Load signature
    signature <- esr_loadPretrainedSignature()
    
    # 3. Classify samples
    predictions <- esr_classifyEndometrial(
        X_new = X_new,
        signature = signature,
        threshold = 0.5,
        confidence = TRUE
    )
    
    # Extract predictions if wrapped in list
    if (is.list(predictions) && "predictions" %in% names(predictions)) {
        pred_df <- predictions$predictions
    } else {
        pred_df <- predictions
    }
    
    # 4. Visualize (signature plots)
    p1 <- plotEndometrialCoefLollipop(signature)
    p2 <- plotEndometrialStabilityBars(signature)
    
    # 5. Export
    temp_dir <- tempfile()
    dir.create(temp_dir)
    on.exit(unlink(temp_dir, recursive = TRUE))
    
    export_path <- esr_exportPredictions(
        predictions = pred_df,
        dir = temp_dir,
        formats = "csv"
    )
    
    # Verify all steps completed successfully
    expect_true(is.data.frame(pred_df))
    expect_true(nrow(pred_df) == ncol(X_new))
    expect_true("sample" %in% names(pred_df))
    expect_true("prediction" %in% names(pred_df))
    expect_true("probability" %in% names(pred_df))
    expect_true(inherits(p1, "ggplot"))
    expect_true(inherits(p2, "ggplot"))
    expect_true(file.exists(export_path))
})

