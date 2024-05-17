It occurred to us that there are OLS Models that have only zero coefficients, but their elastic net counterpart has normal coefficients. This is in contrast to the fact that both methods are based on  the same gene specific feature segment set, pre-selected at the en of the segmentation step.
## Description:
To check how often that occurs and what might be the reason for this, I gathered the data of the Pearson based segment feature preselection of the OLS and ElNet models and created a  [[OnlyZeroCoeffOLS_entries.csv | combined dataset]] with gene id specific entries.

The combined dataset of the standard cross validation was then filtered for entries that have OLS with only zero coefficients and p_val and q_value respectively < 0,05. (LeaveOneOut could't be used, due to p and q val are not calculated in that type of cross validation).
The resulting entries were sorted descending based on their Pearson correlation and a subset of 5 [[OnlyZeroCoeffOLS_entries.csv_reduced| test samples]] was selected.

Next, for the [[PredictedVsActualExpression]], the [[Performance_Overview.txt | Performance Evaluation]], for each text sample additional data was gathered in form of the [[counts.matrix.norm.tsv | RNASeq expression data]], the gene specific models with their lambda.min feature coefficients, the gene specific feature segments that were used as input by the regression.
## Observations/ Comments
 * If the RNA Seq counts are pretty small, it can occur that during the normalization, upstream of the ElNet training, the values can be transformed to negative representations within the "normalized values training space" and this can lead to negative coefficients.
 * 
## Data
### input
### Performance Evaluations
```bash
input_directory_path_LeaveOneOut <- 
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\leaveOneOut_regression/performance_evaluation/Performance_Overview.txt"

input_directory_path_standard <- 
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\standard_regression/performance_evaluation/Performance_Overview.txt"
```
## output
 * [[OnlyZeroCoeffOLS_entries.csv]]
```bash
# Path
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\OnlyZeroCoeffOLS\OnlyZeroCoeffOLS_entries.csv"
```
* **ElNet lambda.min feature coefficients**
```
# eg.
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\OnlyZeroCoeffOLS\ENSG00000125629_10\feature_coefficients_ENSG00000125629_10.txt"
```
## additional
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
## Code
https://github.com/Johannes-Zi/master_thesis/blob/main/regression/performance_evaluation/OnlyZeroCoeffOLS/OnlyZeroCoeffOLS.r

## Runs performed
### Combined dataset
[[OnlyZeroCoeffOLS_entries.csv]]