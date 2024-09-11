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
2.4 Filter out non-coding genes
2.5 Filter out genes that cant be translated from ENSEBL to ENTRZ ID
2.5 Use ==top 250 genes== (based on segment correlations) for DisGeNet


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

* The code iterates over the clinical parameters and uses for each clinical parameter the respective top 250 correlations (output of the previously describe filter) as input for the DisGeNet analysis.
*  ==The input gene list can contain multiple representations of the same gene, because there can be multiple ENTREZ ids for the same ENSEBL id - this is no problem because DisGeNet handles that. This is also the reason for less than 250 genes in the DisGeNet analysis==

#### Workflow
1. The gene names are translated from ENSEBL to NCBI ENRZ IDS to fit as DisGeNet input
2. Additional information is collected for Each gene (Alias, gene name, gene type, etc.)
3. The rate of failed mappings is calculated (now always zero because previous filter kick out exceptions)
4. Run DisGeNet - ==p and q value cutoff of 0.4 (because of much noise and small input gene set)==
5. Save whole DisGeNet results for clinical parameter
6. Extract result entries which contain the ==substrings "pulmo", "arter" or "hypertension"== and save them as separate output (potentially associated disease)
7. Count for each gene that occurs in the potentially associated disease  how often they are assigned to a disease and save those gene list as output 
# Output
#### General output directories
```
# General output directories
# Inside VSC

C:\Users\johan\VSCode_projects\bioinf_master\AssociationAnalysis\disgenet2r\runs\combined_corr_cv_pear_04_thres\spear_thres_04_up_to_250\DisGeNet_results/

# In Cloud
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\disgenet_data\spear_thres_04_up_to_250" 
```
#### Parameter specific output
```
# Input used for DisGeNet run
<clinical_parameter>_DisGeNet_input_genes_metadata.csv

# All DisGeNet results
<clinical_parameter>_DisGeNet_results.csv

# Filtered DisGeNet results that are disease that might be assotiated to pulmanory hypertension
<clinical_parameter>_filtered_DisGeNet_results.csv

# Gene occurance counts of the filtered DisGeNet output described above
<clinical_parameter>_filtered_DisGeNet_results_gene_abundancies.csv
```
#### Combined DisGeNet results (unfiltered) across all clinical parameters
```
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\disgenet_data\spear_thres_04_up_to_250\DisGeNETresults_0.4pqcutoff_uptotop250correlations_latest_version.csv"
```
# Runlog
```
# Logfile location
"C:\Users\johan\OneDrive\dateien_cloud\Master\Semester_4\Masterarbeit\data\pulmanory_hypertension\disgenet_data\spear_thres_04_up_to_250\run_disgenet2r.r_runlog_091124.log"

C:\Users\johan\VSCode_projects\bioinf_master\AssociationAnalysis\disgenet2r\runs\combined_corr_cv_pear_04_thres\spear_thres_04_up_to_250\run_disgenet2r.r_runlog_091124.log
```

# Visualization
