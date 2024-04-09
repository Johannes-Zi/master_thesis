# This script iterates over the elastic net models, that were trained during the
# regression step of STICHIT. It count the number of features (segments, that
# are putative regulatory elements) that were used for the model and count how
# many of those features have coefficients, that are 0 (and thus meaningless).

# Define Input parameters
input_file_path <- paste(
  "C:\\Users\\johan\\Desktop\\regession_output_example\\",
  "Elasticnet_Regression_Model_Segmentation_ENSG00000001167_10_Spearman.RData",
  sep = ""
)

# Echo path of imported file
cat(sprintf("File path for import:\n%s\n", input_file_path))

# Load the .RData file - accesible as 'elasticnet_model'
#load(file.choose())
load(input_file_path)

# Print imported df
print(View(elasticnet_model))
print(View(elasticnet_model$model$lambda.min))

print(elasticnet_model.model.lambda.min)

# Count the number of features used in the glmnet model
#num_features <- length(coef(your_model))

# Count the number of features with coefficients equal to 0
#num_zero_coeff_features <- sum(coef(your_model) == 0)

# Print the results
#cat("Number of features used in the model:", num_features, "\n")
#cat("Number of features with coefficients equal to 0:", num_zero_coeff_features, "\n")