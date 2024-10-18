(script: https://github.com/Johannes-Zi/master_thesis/blob/main/AssociationAnalysis/CorrelationVisualization/visualize_top_correlations.r)

# Idea
Behind every geneid that was used as input for the DisGeNet or GO enrichment analysis there is a correlation between the ATAC accessibility of at least one segment and the values of a clinical parameter. This correlation is one of the top 250 correlations for the specific parameter. But this correlation can be driven by outliers, resulting in overestimation of the importance of specific genes in the gene set for the enrichment analysis.
To circumvent this problem all correlations of the top genes used for the analysis were visualized automatically to enable a manual curation of the most interesting results based on the correlation.

For this, the input gene sets of the enrichment analysis of the different clinical parameters were imported and combined to a single set of ==1720 genes== (to be able to check for each top correlation also the performance on other parameters).
Those gene ids were then used as filter to visualize only segments that are associated to those genes. This resulted in a set of ==6973 segments==. The correlations o all those segments across all clinical parameters were calculated and visualized.

The Code is an adapted version of the BasicCorrelation code that includes additionally the filter to reduce the input an thus the total correlations. And an additional Visualization for each calculated correlation.
It has to be noted that the visualized linear correlations were computed with the ggplot geom_smooth linear model and are not/ can be based on the Spearman correlation (because its rank based and also capturing nonlinear dependencies).

# Example plot
![[AssociationAnalysis/CorrelationVisualization/X0.711765807251866_ENSG00000182310_chr19.51690548.51690747.png]]

# Input Data
```
	# Import clinical metadata
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    
	# Import elnet segments atac data
	elnet_segments_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/CollapsedSegmentation/elnet_model_segments_df.tsv"

	# Import gene model specific outer cross validation performance estimations
	LeaveOneOutCV_performance_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.tsv"

	# Import filtered correlations gene ids (DisGeNet/ GO-enrichment input set)
	filtered_correlations_top_directory_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/go_enrichment/spear_thres_04_up_to_250/corr_based_gene_filtering/"
```
# Workflow
- Same workflow as [[BasicCorrelation]], but with an additional filter step (for the segments that belong to the top 1720 genes), placed after the step that filters out the gene models that performed bad in the outer cross validation of DisGeNet. (Additional filter to the Stichit outer cv filter and the filter of using only coding genes as input for the enrichment analysis)
- Additional Visualization of the clinical parameter + ATAC accessability with ggplot geom smooth linear model based linear correlations.

# Output
#### General output directories
```
# Cloud
C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\CorrelationVisualization\runs\v1\

# VSC
C:\Users\johan\VSCode_projects\bioinf_master\AssociationAnalysis\CorrelationVisualization\runs\v1
```
#### Parameter specific output
```
# Calculated Spearman correlations
<clinical_parameter>_spearman_correlations.tsv

# Visualized correlations
./plots/
X<correlation_value>_<gene_id>_<chromosomal_location>.png
# Example
X-0.1524323_ENSG00000091490_chr4.25861760.25861849.png
```

#### Logging
```
# In the respective run directory
run1.1.log
run1.2.log
```