It occurred to us that there are OLS Models that have only zero coefficients, but their elastic net counterpart have normal coefficients. This is in contrast to the fact that both methods are based on  the same gene specific feature segment set, pre-selected at the en of the segmentation step.
## Description:
To check how often that occurs and what might be the reason for this, I gathered the data of the Pearson based segment feature preselection of the OLS and ElNet models and created a  [[OnlyZeroCoeffOLS_entries.csv | combined dataset]] with gene id specific entries.
## Observations/ Comments
 * If the RNA Seq counts are pretty small, it can occur that during the normalization, upstream of the ElNet training, the values can be transformed to negative representations within the "normalized values training space" and this can lead to negative coefficients.
## Data
### input
#### Performance Evaluations
```bash
input_directory_path_LeaveOneOut <- 
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\leaveOneOut_regression/performance_evaluation/Performance_Overview.txt"

input_directory_path_standard <- 
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\standard_regression/performance_evaluation/Performance_Overview.txt"
```
#### OLS Models
```bash
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\standard_regression\output_data\regression_output.zip/regression_output/*_Pearson.RData"
```
#### Elastic-net models
```bash
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\regression\leaveOneOut_regression\output_data\regression_output.zip/regression_output/*_Pearson.RData"
```
### output
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
## Code
https://github.com/Johannes-Zi/master_thesis/blob/main/regression/performance_evaluation/OnlyZeroCoeffOLS/OnlyZeroCoeffOLS.r

## Runs performed
only one run performed for the standard regression version - beacuse in LeaveOneOut no p and q values were calculated in the performance report and thus no filtering with that metrics could be done.
### Combined dataset
[[OnlyZeroCoeffOLS_entries.csv]]