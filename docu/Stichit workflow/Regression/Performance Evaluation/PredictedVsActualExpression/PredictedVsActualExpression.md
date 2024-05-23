To get an better overview what was happening in the quite often occurring cases, described and over-viewed in [[OnlyZeroCoeffOLS]], we looked into five specific gene models.

The combined dataset of the standard cross validation was then filtered for entries that have OLS with only zero coefficients and p_val and q_value respectively < 0,05. (LeaveOneOut could't be used, due to p and q val are not calculated in that type of cross validation).
The resulting entries were sorted descending based on their Pearson correlation and a subset of 5 [[OnlyZeroCoeffOLS_entries.csv_reduced| test samples]] was selected.

For those gene models we gathered next to the [[Performance_Overview.txt | Performance Evaluation]] the [[counts.matrix.norm.tsv | RNASeq expression data]], the gene specific models with their lambda.min feature coefficients, the gene specific feature segments that were used as input by the regression. (The collected files are listed below in the data sections together with example paths)

Additional to collecting model specific data for more insights, we performed compared the actual gene expression values of the initial dataset with the values that were predicted by the trained models based on the ATAC seq data. We did this to get and impression if the calculated correlations are driven by some outliers, or if the predictions are consistent within each other.

To do this, we had to load the gene specific ATAC Seq data over the samples together with the respective RNA Seq expression. Both is represented in the gene specific output of the segmentation.
The elastic nets were trained with normalized and scaled output of the segmentation. This means the ATAC Seq and the RNA Seq data has to be transformed into the same space as the data used for training.

The normalization is done by log transforming the data:
log(value +1); +1 to avoid log of zero

Finally, the data is being scaled. The values are standardized by subtracting the column mean and dividing by the column standard deviation. This results in data where each column has a mean of 0 and a standard deviation of 1.

# Observations
* the small RNA Seq expression variation in gene ENSG00000259056 and ENSG00000249092, leads in the dotplot to accommodations on specific x axis locations (original gene expression values)
* The segmentation of ENSG00000156232 seems to 
# Data

### input
```r
#' Input paths for the gene specific models

# Path to the segmentation output file
segmentation_path_1 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000125629_10/Segmentation_ENSG00000125629_10_Pearson.txt"
# Path to the elastic net model file
elnet_path_1 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000125629_10/Elasticnet_Regression_Model_Segmentation_ENSG00000125629_10_Pearson.RData"

segmentation_path_2 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000156232_10/Segmentation_ENSG00000156232_10_Pearson.txt"
elnet_path_2 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000156232_10/Elasticnet_Regression_Model_Segmentation_ENSG00000156232_10_Pearson.RData"

segmentation_path_3 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000249092_10/Segmentation_ENSG00000249092_10_Pearson.txt"
elnet_path_3 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000249092_10/Elasticnet_Regression_Model_Segmentation_ENSG00000249092_10_Pearson.RData"

segmentation_path_4 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000259056_10/Segmentation_ENSG00000259056_10_Pearson.txt"
elnet_path_4 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000259056_10/Elasticnet_Regression_Model_Segmentation_ENSG00000259056_10_Pearson.RData"

segmentation_path_5 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000211788_10/Segmentation_ENSG00000211788_10_Pearson.txt"
elnet_path_5 <- "C:/Users/johan/Desktop/local_master_thesis_data/PredictedVsActualExpression/ENSG00000211788_10/Elasticnet_Regression_Model_Segmentation_ENSG00000211788_10_Pearson.RData"
```

### output
```bash
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\PredictedVsActualExpression"
```

### additional
 * [[counts.matrix.norm.tsv | RNASeq expression data]]
```bash
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\raw_data\RNA\counts.matrix.norm.tsv"
```
* **gene specific models**
```bash
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\leaveOneOut_regression\output_data\regression_output.zip/*_Pearson.Rdata"

"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\standard_regression\output_data\regression_output.zip/*_Pearson.Rdata"
```
* **gene specific segmentation feature sets**
```
Output of the segmentation
eg. Segmentation_ENSG00000125629_10_Pearson.txt
```

# Code
https://github.com/Johannes-Zi/master_thesis/blob/main/regression/performance_evaluation/PredictedVsActualExpression/PredictedVsActualExpression.r

# Runs

| gene            | figure             | inital RNA Seq                                          | feature coefficients                            | normalized ATAC/ RNA Seq                        | normalized, scaled, centered ATAC/ RNA Seq           |
| --------------- | ------------------ | ------------------------------------------------------- | ----------------------------------------------- | ----------------------------------------------- | ---------------------------------------------------- |
| ENSG00000125629 | [[fig_2105242319]] | [[counts.matrix.norm.tsv_ENSG00000125629_10.txt]]       | [[feature_coefficients_ENSG00000125629_10.txt]] | [[segmentation_normalized_ENSG00000125629.tsv]] | [[segmentation_scaled_centered_ENSG00000125629.tsv]] |
| ENSG00000156232 | [[fig_2105242320]] | [[counts.matrix.norm.tsv_entry_ENSG00000156232_10.txt]] | [[feature_coefficients_ENSG00000156232_10.txt]] | [[segmentation_normalized_ENSG00000156232.tsv]] | [[segmentation_scaled_centered_ENSG00000156232.tsv]] |
| ENSG00000249092 | [[fig_2105242321]] | [[counts.matrix.norm.tsv_ENSG00000249092_10.txt]]       | [[feature_coefficients_ENSG00000249092_10.txt]] | [[segmentation_normalized_ENSG00000249092.tsv]] | [[segmentation_scaled_centered_ENSG00000249092.tsv]] |
| ENSG00000259056 | [[fig_2105242322]] | [[counts.matrix.norm.tsv_ENSG00000259056_10.txt]]       | [[feature_coefficients_ENSG00000259056_10.txt]] | [[segmentation_normalized_ENSG00000259056.tsv]] | [[segmentation_scaled_centered_ENSG00000259056.tsv]] |
| ENSG00000211788 | [[fig_2105242323]] | [[counts.matrix.norm.tsv_ENSG00000211788_10.txt]]       | [[feature_coefficients_ENSG00000211788_10.txt]] | [[segmentation_normalized_ENSG00000211788.tsv]] | [[segmentation_scaled_centered_ENSG00000211788.tsv]] |
