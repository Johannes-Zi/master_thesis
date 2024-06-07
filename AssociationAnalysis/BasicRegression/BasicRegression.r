
#' Import clinical metadata
#' 
import_clinical_meatadata <- function(file_path) {
  # Import clinical metadata
  clinical_metadata <- read.csv(file_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  return(clinical_metadata)
}

#'
#' Import elnet segments atac data
#' 
import_elnet_segments_atac_data <- function(file_path) {
  # Import segments atac data
  segments_atac_data <- read.csv(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(segments_atac_data)
}

#'
#' Extract represented patient ids in the elnet segments atac data
#' 
extract_represented_patient_ids <- function(elnet_segments_atac_data) {
  # Get List with all columnnames
  columnnames <- colnames(elnet_segments_atac_data)

  # Extract columnnames wihich end on .ATAC
  atac_columns <- grep(".ATAC", columnnames, value = TRUE)
  
  # Extract out of columnnames represented patient ids
  patient_ids <- gsub(".ATAC", "", atac_columns)

  return(patient_ids)
}

#'
#' Returns reduced clinical metadata df with only the patient ids which are represented in the elnet segments atac data
#' 
reduced_clinical_metadata <- function(clinical_metadata, patient_ids) {
  
  # Replace - seperator within the rna.final.ids column with dots
  clinical_metadata$rna.final.ids <- gsub("-", ".", clinical_metadata$rna.final.ids)

  # Reduce clinical metadata to only the patient ids which are represented in the elnet segments atac data
  reduced_clinical_metadata <- clinical_metadata[clinical_metadata$rna.final.ids %in% patient_ids, ]

  # Reset rownames
  rownames(reduced_clinical_metadata) <- 1:nrow(reduced_clinical_metadata)

  return(reduced_clinical_metadata)
}

#'
#' Function calculates Spearman correlation between a handed over column of the clinical metadata and the atac data
#' represented as rows in the elnet segments atac data
#'
calculate_spearman_correlation <- function(clinical_metadata_rna_ids, clinical_metadata_column, elnet_segments_atac_data) {
  # Create df with the handed over column of the clinical metadata ant the RNA ids
  clinical_metadata_column_df <- data.frame( clinical_metadata_rna_ids, clinical_metadata_column)
  
  # Sort rows by patient ids
  clinical_metadata_column_df <- clinical_metadata_column_df[order(clinical_metadata_column_df$clinical_metadata_rna_ids),]

  # Set rownames to the RNA ids
  rownames(clinical_metadata_column_df) <- clinical_metadata_column_df$clinical_metadata_rna_ids

  # Remove RNA ids from the df
  clinical_metadata_column_df <- clinical_metadata_column_df[-1]

  print(clinical_metadata_column_df)

  # Extract first row of the elnet segments atac data
  current_segment_atac <- elnet_segments_atac_data[1,]
  
  # Extract segment and gene_id
  current_segment <- current_segment_atac$segment
  current_gene_id <- current_segment_atac$gene_id

  # Remove segment and gene_id from the elnet segments atac data
  current_segment_atac <- current_segment_atac[-c(1,2)]

  # Transpose the elnet segments atac data
  current_segment_atac <- as.data.frame(t(current_segment_atac))

  # Remove trailing .ATAC from the rownames
  rownames(current_segment_atac) <- gsub(".ATAC", "", rownames(current_segment_atac))
  
  # Rename the columnname 1 to ATAC
  colnames(current_segment_atac) <- "ATAC"
  
  # Add column based on rownames
  current_segment_atac$ids <- rownames(current_segment_atac)

  # Sort dataframe by ids
  current_segment_atac <- current_segment_atac[order(current_segment_atac$ids),]
  
  # Remove ids column
  current_segment_atac <- current_segment_atac[-2]
  
  print(current_segment_atac)
  print(clinical_metadata_column_df)
  
  # Extract the numeric vector from clinical_metadata_column_df
  # Assuming it has only one column after removing the IDs
  clinical_vector <- as.numeric(clinical_metadata_column_df[, 1])

  # Extract the numeric vector from current_segment_atac
  atac_vector <- as.numeric(current_segment_atac$ATAC)

  # Calculate spearman correlation between the two dfs
  spearman_correlation <- cor.test(clinical_vector, atac_vector, method = "spearman")
  print(spearman_correlation)
}

# Import clinical metadata
clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)

# Import elnet segmeants atac data
elnet_segments_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/CollapsedSegmentation/elnet_model_segments_df.tsv"
elnet_segments_atac <- import_elnet_segments_atac_data(elnet_segments_path)

# Extract represented patient ids in the elnet segments atac data
patient_ids <- extract_represented_patient_ids(elnet_segments_atac)

# Reduce clinical metadata to only the patient ids which are represented in the elnet segments atac data
reduced_clinical_metadata <- reduced_clinical_metadata(clinical_metadata, patient_ids)
print(reduced_clinical_metadata)

# Print sorted and ; separated final rna ids
#cat(paste("Sorted and ; separated final rna ids:", paste(sort(reduced_clinical_metadata$`rna.final.ids`), collapse = ";")))
# Print sorted and ; separated final patient_ids
#cat(paste("Sorted and ; separated final patient_ids:", paste(sort(patient_ids), collapse = ";")))


# Example call for regresion
calculate_spearman_correlation(reduced_clinical_metadata$rna.final.ids,reduced_clinical_metadata$`mPAP`, elnet_segments_atac)
