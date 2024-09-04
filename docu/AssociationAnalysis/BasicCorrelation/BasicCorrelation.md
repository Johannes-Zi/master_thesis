(code: https://github.com/Johannes-Zi/master_thesis/blob/main/AssociationAnalysis/BasicCorrelation/BasicCorrelation.r )
(data: C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/BasicCorrelation/ )
# Idea
The idea is to determine if the by Stichit detected potentially regulatory active segments can be associated with specific clinical parameters. This is done by correlating the accessibility of each segment together with a measured clinical parameter (like blood pressure). This can be done by using each patient as data point and correlation the ATAC Seq values of the segment with the patient specifically measured clinical parameter.
The result of this analysis is to determine segments with a strong correlation to clinical parameter, which could mean that the accessibility of those segments could play a role in the mechanism behind pulmonary hypertension.

# Input
```
	# Import clinical metadata
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    
	# Import elnet segments atac data
	elnet_segments_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/CollapsedSegmentation/elnet_model_segments_df.tsv"

	# Import gene model specific outer cross validation performance estimations
	LeaveOneOutCV_performance_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.tsv"
```

# Workflow
## 1. Import data
- clinical metadata with patient specific parameters
- ATAC Seq values of by Stichit selected segments
## 2. Join datasets
* Extracts patient ids that are represented in the ATAC seq dataframe
* reduce the clinical metadata set to the ATAC patient ids
## 3. Cross validation based input filter
Filters out the gene specific models based on their performance in the LeaveOneOut outer cross validation performed for the gene specific elastic net models.

==Filtered out gene models with CV Pearson correlation < |0.4| ==

## 4. Spearman Correlations
Correlate the ATAC accessibility of a specific segment across all patients together with the corresponding patient specific values of a specific specific clinical parameter.
Those correlations were done for all segments and also each segment values were correlated against all clinical parameters.

The Spearman correlation was chosen because the data does not need to follow any distributions and thus does not need to be normalized. (laut Marcel)
## 5. Potential filter based on p and q values
This Filter was discussed, but is not applicable based on the calculated values (p close to evenly distributed and q values not smaller than 0.5)
# Output
The data is stored for each clinical parameter separately
```
# In local repo
C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/BasicCorrelation/combined_corr_cv_pear_04_thres/

# Second storage location in cloud
"C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data_correlations/combined_corr_cv_pear_04_thres/"
```
## Correlations
 * stored as  ```<param>_spearman_correlations.tsv``` in a directory with the respectively correlated clinical parameter
## Plots
* ``` p_value_hist.png  # histogram plot of the p_values of the correlations```
* ``` q_value_hist.png  # histogram plot of the p_values of the correlations```
* ``` spearman_corr_density.png  # density plot of the spearman_correlations ```
* ``` spearman_corr_density_smaller_than_005_Spearman.png  # density plot of the spearman_correlations over the p_value < 0.05 -  Filter the dataframe to include only rows where spearman_corr is below 0.05 ```
* ``` spearman_corr_density_p_value.png  # density plot of the spearman_correlations seperated by p_value < 0.05 ```
## Runlog

```
Reduce clinical metadata to the patient ids which are represented in the elnet segments atac data

Extracted represented patient ids in the elnet segments atac data:  26

Number of patients represented in the clinical metadata with RNA Seq ids:  52

Patients represented in the reduced clinical metadata:  26

  

Apply elnet outer cross validation based pearson correlation filtering

with threshold Pearson >= 0.4 or <= -0.4

Number of unique genes represented in LeaveOneOutCV_performance:  9170

Number of entries with negative Pearson correlation:  571

Number of genes after filtering:  5098

  

Filter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list.

Number of unique genes represented by elnet segments:  8987

Number of elnet segments before filtering:  33162

Number of elnet segments after filtering:  19692

  

Run the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over all clinical parameters.
```
