# Integration Tests for Cross-Mode Workflows
# Tests workflows that span multiple modes

test_that("Mode 2 trained signature works in Mode 1 classification", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))
    skip_if_not(exists("gse201926_sample", where = asNamespace("endoSignatureR"), mode = "any"))
    
    library(glmnet)
    library(rsample)
    
    # Train signature (Mode 2)
    data(gse201926_trainmini)
    set.seed(123)
    result <- suppressWarnings(esr_trainEndometrialSignature(
        X = gse201926_trainmini$counts,
        pheno = gse201926_trainmini$pheno,
        top_k = 50,
        outer_folds = NULL,
        seed = 123
    ))
    
    # Use for classification (Mode 1)
    data(gse201926_sample)
    predictions <- esr_classifyEndometrial(
        X_new = gse201926_sample$counts,
        signature = result$signature,
        threshold = 0.5,
        confidence = FALSE
    )
    
    # Extract predictions if wrapped
    if (is.list(predictions) && "predictions" %in% names(predictions)) {
        pred_df <- predictions$predictions
    } else {
        pred_df <- predictions
    }
    
    # Verify predictions are valid
    expect_true(is.data.frame(pred_df))
    expect_true(nrow(pred_df) == ncol(gse201926_sample$counts))
    expect_true("prediction" %in% names(pred_df))
    expect_true("probability" %in% names(pred_df))
})

test_that("Mode 1 predictions can be visualized with Mode 3 functions", {
    skip_if_not_installed("pROC")
    
    csv_path <- system.file("extdata", "pretrained-signature", 
                            "endometrial_signature.csv", package = "endoSignatureR")
    skip_if(
        nchar(csv_path) == 0 || !file.exists(csv_path),
        "Pretrained signature artifacts not found"
    )
    
    # Classify samples (Mode 1)
    data(gse201926_sample)
    signature <- esr_loadPretrainedSignature()
    predictions <- esr_classifyEndometrial(
        X_new = gse201926_sample$counts,
        signature = signature,
        threshold = 0.5
    )
    
    # Extract predictions
    if (is.list(predictions) && "predictions" %in% names(predictions)) {
        pred_df <- predictions$predictions
    } else {
        pred_df <- predictions
    }
    
    # Use Mode 3 visualization with predictions
    # Transform data for PCA
    mat_t <- esr_transform_log1p_cpm(gse201926_sample$counts)
    
    # Create pheno with predictions
    pheno_with_pred <- gse201926_sample$pheno
    pheno_with_pred$predicted_group <- pred_df$prediction
    
    # Visualize with Mode 3 functions
    p1 <- plotEndometrialPCA(mat_t, pheno = pheno_with_pred, group_col = "predicted_group")
    
    # Verify visualization works
    expect_true(inherits(p1, "ggplot"))
})

test_that("Mode 3 DE analysis can inform Mode 2 training", {
    skip_if_not_installed("glmnet")
    skip_if_not_installed("rsample")
    skip_if_not(exists("gse201926_trainmini", where = asNamespace("endoSignatureR"), mode = "any"))
    
    library(glmnet)
    library(rsample)
    
    # Perform DE analysis (Mode 3)
    data(gse201926_trainmini)
    mat_t <- esr_transform_log1p_cpm(gse201926_trainmini$counts)
    de_table <- esr_analyzeDifferentialExpression(mat_t, gse201926_trainmini$pheno)
    
    # Select top genes from DE (Mode 3)
    top_genes <- esr_selectTopGenes(de_table = de_table, n = 50, by = "de")
    
    # Use selected genes to inform training (Mode 2)
    # Note: esr_trainEndometrialSignature uses top_k parameter, but we can verify
    # that DE-selected genes are reasonable for training
    expect_true(length(top_genes) > 0)
    # Verify that selected genes exist in the DE table (gene_id column)
    expect_true(all(top_genes %in% de_table$gene_id))
    
    # Train signature (Mode 2) - top_k should select from similar gene set
    set.seed(123)
    result <- suppressWarnings(esr_trainEndometrialSignature(
        X = gse201926_trainmini$counts,
        pheno = gse201926_trainmini$pheno,
        top_k = 50, # Similar to DE selection
        outer_folds = NULL,
        seed = 123
    ))
    
    # Verify training completed
    expect_true(is.list(result))
    expect_true("signature" %in% names(result))
})

