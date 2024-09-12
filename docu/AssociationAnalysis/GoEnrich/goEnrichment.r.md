(script: https://github.com/Johannes-Zi/master_thesis/blob/main/AssociationAnalysis/GoEnrich/goEnrichment.r)

# Idea
Performed GO enrichment analysis with similar approach to the [[disgenet2r | DisGeNet analysis]]  by using the top correlations that resulted out of the correlation between clinical parameters vs segmental ATAC accessibility.

# Input Data
```
	# Import clinical metadata and extract column names
	clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
```
# Workflow
## 1. Import Correlations
Same as [[disgenet2r | DisGeNet analysis]] 

## 2. Filter Correlations / transform to genes
Same as [[disgenet2r | DisGeNet analysis]] - usage of the ==top 250 correlations per clinical parameter==

## 3. Run GO enrichment
1. The gene names are translated from ENSEBL to NCBI ENRZ IDS to fit as GoEnrich input
2. Additional information is collected for each gene (Alias, gene name, gene type, etc.)
3. The rate of failed mappings is calculated (now always zero because previous filter kick out exceptions)
4. Run GoEnrich - ==p and q value cutoff of 0.6 (because of much noise and small input gene set)==
5. Save whole DisGeNet results for each clinical parameter
# Output
#### General output directories
```
# Cloud
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\go_enrichment\spear_thres_04_up_to_250"

# VSC
"C:\Users\johan\VSCode_projects\bioinf_master\AssociationAnalysis\GoEnrich\runs\spear_thres_04_up_to_250"
```
#### Parameter specific output
```
# Input geneset with metadata used for GoEnrich run
<clinical_parameter>_GoEnrich_input_genes_metadata.csv

# All GoEnrich results for this clinical parameter
<clinical_parameter>_GoEnrich_results.csv
```
#### Combined GoEnrich results across all clinical parameters
```
# In the respective run directory
GoEnrichresults_0.4pqcutoff_uptotop250correlations_latest_version.csv
```

#### Logging
```
# In the respective run directory
run_091224.log
```