# Test Pre-trained Artifact Validation

test_that("pretrained signature CSV has correct structure", {
    # Get pretrained artifact path (in pretrained-signature subfolder)
    csv_path <- system.file("extdata", "pretrained-signature", "endometrial_signature.csv", package = "endoSignatureR")

    # Skip test if artifact doesn't exist (e.g., during development)
    skip_if(
        nchar(csv_path) == 0 || !file.exists(csv_path),
        "Pretrained signature CSV not found in inst/extdata/pretrained-signature/"
    )

    # Read CSV with explicit column types to ensure proper parsing
    signature_csv <- readr::read_csv(
        csv_path,
        col_types = readr::cols(
            gene_id = readr::col_character(),
            coefficient = readr::col_double(),
            selection_frequency = readr::col_integer(),
            bootstrap_frequency = readr::col_double()
        )
    )

    # Verify expected columns
    expect_true("gene_id" %in% names(signature_csv))
    expect_true("coefficient" %in% names(signature_csv))
    expect_true("selection_frequency" %in% names(signature_csv))

    # Verify bootstrap_frequency is optional (may or may not be present)
    has_bootstrap <- "bootstrap_frequency" %in% names(signature_csv)

    # Verify data types (with explicit column types, these should be correct)
    expect_type(signature_csv$gene_id, "character")
    expect_type(signature_csv$coefficient, "double")
    expect_type(signature_csv$selection_frequency, "integer")

    if (has_bootstrap) {
        expect_type(signature_csv$bootstrap_frequency, "double")
    }

    # Verify non-empty signature
    non_intercept <- signature_csv$gene_id != "(Intercept)"
    expect_true(sum(non_intercept) > 0, "Signature should contain genes")

    # Verify coefficients are numeric (non-NA for genes)
    expect_true(all(!is.na(signature_csv$coefficient[non_intercept])))

    # Verify selection frequencies are non-negative integers
    expect_true(all(signature_csv$selection_frequency[non_intercept] >= 0))

    if (has_bootstrap) {
        # Verify bootstrap frequencies are between 0 and 1 (exclude NA values)
        bootstrap_non_intercept <- signature_csv$bootstrap_frequency[non_intercept]
        bootstrap_non_na <- bootstrap_non_intercept[!is.na(bootstrap_non_intercept)]
        if (length(bootstrap_non_na) > 0) {
            expect_true(all(bootstrap_non_na >= 0 & bootstrap_non_na <= 1),
                info = "Bootstrap frequencies should be between 0 and 1"
            )
        }
    }
})

test_that("pretrained signature JSON has correct schema", {
    # Get pretrained artifact path (in pretrained-signature subfolder)
    json_path <- system.file("extdata", "pretrained-signature", "endometrial_recipe.json", package = "endoSignatureR")

    # Skip test if artifact doesn't exist (e.g., during development)
    skip_if(
        nchar(json_path) == 0 || !file.exists(json_path),
        "Pretrained recipe JSON not found in inst/extdata/pretrained-signature/"
    )

    # Read JSON
    recipe_json <- jsonlite::fromJSON(json_path)

    # Verify JSON is a list
    expect_type(recipe_json, "list")

    # Verify required sections
    required_sections <- c("preprocessing", "training", "signature", "reproducibility")
    expect_true(all(required_sections %in% names(recipe_json)))

    # Verify preprocessing parameters
    expect_true("preprocessing" %in% names(recipe_json))
    expect_type(recipe_json$preprocessing, "list")
    expect_true("transform" %in% names(recipe_json$preprocessing))
    expect_true("cpm_min" %in% names(recipe_json$preprocessing))
    expect_true("top_k" %in% names(recipe_json$preprocessing))

    # Verify training parameters
    expect_true("training" %in% names(recipe_json))
    expect_type(recipe_json$training, "list")
    expect_true("lambda_rule" %in% names(recipe_json$training))
    expect_true("calibration_method" %in% names(recipe_json$training))

    # Verify signature metadata
    expect_true("signature" %in% names(recipe_json))
    expect_type(recipe_json$signature, "list")
    expect_true("n_genes" %in% names(recipe_json$signature))
    expect_true("intercept" %in% names(recipe_json$signature))
    expect_type(recipe_json$signature$n_genes, "integer")
    expect_type(recipe_json$signature$intercept, "double")

    # Verify reproducibility info
    expect_true("reproducibility" %in% names(recipe_json))
    expect_type(recipe_json$reproducibility, "list")
    expect_true("package_version" %in% names(recipe_json$reproducibility))
    expect_true("r_version" %in% names(recipe_json$reproducibility))
})

test_that("pretrained artifacts are accessible via system.file()", {
    # Test CSV access (in pretrained-signature subfolder)
    csv_path <- system.file("extdata", "pretrained-signature", "endometrial_signature.csv", package = "endoSignatureR")
    skip_if(
        nchar(csv_path) == 0 || !file.exists(csv_path),
        "Pretrained signature CSV not found"
    )
    expect_true(nchar(csv_path) > 0)
    expect_true(file.exists(csv_path))

    # Test JSON access (in pretrained-signature subfolder)
    json_path <- system.file("extdata", "pretrained-signature", "endometrial_recipe.json", package = "endoSignatureR")
    skip_if(
        nchar(json_path) == 0 || !file.exists(json_path),
        "Pretrained recipe JSON not found"
    )
    expect_true(nchar(json_path) > 0)
    expect_true(file.exists(json_path))

    # Test model card access (optional, in pretrained-signature subfolder)
    md_path <- system.file("extdata", "pretrained-signature", "endometrial_model_card.md", package = "endoSignatureR")
    if (nchar(md_path) > 0 && file.exists(md_path)) {
        expect_true(file.exists(md_path))
    }

    # Test stability CSV access (optional, in pretrained-signature subfolder)
    stability_path <- system.file("extdata", "pretrained-signature", "endometrial_stability.csv", package = "endoSignatureR")
    if (nchar(stability_path) > 0 && file.exists(stability_path)) {
        expect_true(file.exists(stability_path))
    }
})

test_that("pretrained artifacts match expected schema", {
    # Get artifact paths (in pretrained-signature subfolder)
    csv_path <- system.file("extdata", "pretrained-signature", "endometrial_signature.csv", package = "endoSignatureR")
    json_path <- system.file("extdata", "pretrained-signature", "endometrial_recipe.json", package = "endoSignatureR")

    # Skip test if artifacts don't exist
    skip_if(
        nchar(csv_path) == 0 || !file.exists(csv_path) ||
            nchar(json_path) == 0 || !file.exists(json_path),
        "Pretrained artifacts not found"
    )

    # Read CSV with explicit column types to ensure proper parsing
    signature_csv <- readr::read_csv(
        csv_path,
        col_types = readr::cols(
            gene_id = readr::col_character(),
            coefficient = readr::col_double(),
            selection_frequency = readr::col_integer(),
            bootstrap_frequency = readr::col_double()
        )
    )

    # Read JSON
    recipe_json <- jsonlite::fromJSON(json_path)

    # Verify n_genes in JSON matches number of genes in CSV (excluding intercept)
    n_genes_csv <- sum(signature_csv$gene_id != "(Intercept)")
    n_genes_json <- recipe_json$signature$n_genes

    expect_equal(n_genes_csv, n_genes_json,
        info = "Number of genes in CSV should match JSON"
    )

    # Verify intercept in JSON matches intercept in CSV (if present)
    intercept_row <- signature_csv[signature_csv$gene_id == "(Intercept)", ]
    if (nrow(intercept_row) > 0) {
        # Use more lenient tolerance for floating-point comparison (CSV vs JSON may have slight rounding differences)
        expect_equal(intercept_row$coefficient, recipe_json$signature$intercept,
            tolerance = 1e-4,
            info = "Intercept in CSV should match JSON (within tolerance)"
        )
    }

    # Verify transform method is valid
    valid_transforms <- c("log1p-cpm", "vst")
    expect_true(recipe_json$preprocessing$transform %in% valid_transforms,
        info = "Transform method should be valid"
    )

    # Verify lambda_rule is valid
    valid_lambda_rules <- c("1se", "min")
    expect_true(recipe_json$training$lambda_rule %in% valid_lambda_rules,
        info = "Lambda rule should be valid"
    )

    # Verify calibration method is valid
    valid_calibration <- c("platt", "isotonic", "none")
    expect_true(recipe_json$training$calibration_method %in% valid_calibration,
        info = "Calibration method should be valid"
    )
})

test_that("pretrained stability CSV has correct structure (if available)", {
    # Get stability artifact path (in pretrained-signature subfolder)
    stability_path <- system.file("extdata", "pretrained-signature", "endometrial_stability.csv", package = "endoSignatureR")

    # Skip test if stability CSV doesn't exist (optional)
    skip_if(
        nchar(stability_path) == 0 || !file.exists(stability_path),
        "Pretrained stability CSV not found (optional)"
    )

    # Read stability CSV with explicit column types to ensure gene_id is character
    stability_csv <- readr::read_csv(
        stability_path,
        col_types = readr::cols(
            gene_id = readr::col_character(),
            bootstrap_frequency = readr::col_double()
        )
    )

    # Verify expected columns
    expect_true("gene_id" %in% names(stability_csv))
    expect_true("bootstrap_frequency" %in% names(stability_csv))

    # Verify data types
    expect_type(stability_csv$gene_id, "character")
    expect_type(stability_csv$bootstrap_frequency, "double")

    # Verify bootstrap frequencies are between 0 and 1 (exclude NA values)
    bootstrap_non_na <- stability_csv$bootstrap_frequency[!is.na(stability_csv$bootstrap_frequency)]
    if (length(bootstrap_non_na) > 0) {
        expect_true(all(bootstrap_non_na >= 0 & bootstrap_non_na <= 1),
            info = "Bootstrap frequencies should be between 0 and 1"
        )
    }

    # Verify stability CSV is sorted by frequency (descending)
    if (nrow(stability_csv) > 1) {
        frequencies <- stability_csv$bootstrap_frequency
        expect_true(all(diff(frequencies) <= 0),
            info = "Stability CSV should be sorted by frequency (descending)"
        )
    }
})
