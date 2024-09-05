(script: https://github.com/Johannes-Zi/master_thesis/blob/main/AssociationAnalysis/disgenet2r/run_disgenet2r.r)

# Idea
The idea is to use the gene-set of the top correlating (potentially regulatory) segments (ATAC vs clinical parameters), to determine if in there are subset of genes that are associated with specific disease and over represented in the geneset.
# Input
```
	# Import clinical metadata and extract column names
	clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"



```
# Workflow
## 1. Import Correlations
* Import clinical metadata to extract he column names
* Import clinical parameter specific segmental correlations
## 2. Filter Correlations / transform to genes
2.1 Filter out correlations below |corr| < 0.4
2.2 Extract list of genes that are represented by the segments
2.3 Extract genes that are represented by more than one segment
2.4 Use top XXX genes (based on segment correlations) for DisGeNet


Output
```
# Log file of the clinical parameter specific filtering
<clinical_param>_filtering.log

# Genes that are represented by more than one segment
<clinical_param>_multiple_represented_duplicated_genes.csv

# Segments corresponding to genes that are represented by morethan one segment
<clinical_param>_multiple_represented_genes_segments.csv

# Segments with the corresponding genes that passed all filters and can be used as starting point for the DisGeNet analysis
<clinical_param>_reduced_correlation_results.csv

# Additional insigths across alls clinical parameters can be found in the run directory
```

## 3. Run DisGeNet

2.3 Filter out non-coding genes
##### Input filters

# Output
```
	# Reduced correlation dataset
	output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/disgenet2r/runs/combined_corr_cv_pear_04_thres/spear_thres_04_up_to_200_coding_genes/corr_based_gene_filtering/"

	# 
```

# Runlog