# This script iterates over the elastic net models, that were trained during the
# regression step of STICHIT. It count the number of features (segments, that
# are putative regulatory elements) that were used for the model and count how
# many of those features have coefficients, that are 0 (and thus meaningless).

# Define Input parameters
input_file_path <- "C://User//johan//Desktop//regession_output_example//
Elasticnet_Regression_Model_Segmentation_ENSG00000001167_10_Pearson.RData"

print(sprintf("Imported file:", input_file_path))
#cat(input_file_path)

quit()

# Load the RData file
load(file_path)

# View the head of the dataframe
head(your_dataframe)
