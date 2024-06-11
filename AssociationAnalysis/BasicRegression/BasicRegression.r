
create_new_library <- function() {
new_library_path <- "C:/Users/johan/Desktop/local_master_thesis_data/temp_library/"
dir.create(new_library_path, recursive = TRUE)

# Add to PATH
.libPaths(c(new_library_path, .libPaths()))

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", lib = new_library_path)

BiocManager::install("qvalue", lib = new_library_path, force = TRUE)
}

#create_new_library()


# Create new library

#if (!require("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")

#BiocManager::install("qvalue")

library(qvalue)

#'
#' Import clinical metadata
import_clinical_meatadata <- function(file_path) {
  # Import clinical metadata
  clinical_metadata <- read.csv(file_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  return(clinical_metadata)
}

#'
#' Import elnet segments atac data
import_elnet_segments_atac_data <- function(file_path) {
  # Import segments atac data
  segments_atac_data <- read.csv(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  return(segments_atac_data)
}

#'
#' Extract represented patient ids in the elnet segments atac data
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
#' Create df with  spearman correlation across all elnet segments and the handed over clinical metadata column
create_spearman_correlations_df <- function(clinical_vector, elnet_segments_atac_data){

  # Initialize progress bar
  pb <- txtProgressBar(min = 0, max = nrow(elnet_segments_atac_data), style = 3)

  # Initialize empty df for the spearman correlations
  spearman_correlations_df <- data.frame(segment = character(), gene_id = character(), spearman_corr = numeric(),
                                         p_value = numeric(), stringsAsFactors = FALSE)

  # Iterate over the elnet segments atac data
  for (i in 1:nrow(elnet_segments_atac_data)){
    # Extract current row of the elnet segments atac data
    current_segment_atac <- elnet_segments_atac_data[i,]

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

    # Extract the numeric vector from current_segment_atac
    atac_vector <- as.numeric(current_segment_atac$ATAC)

    # Calculate spearman correlation between the two dfs
    spearman_correlation <- suppressWarnings(cor.test(clinical_vector, atac_vector, method = "spearman"))

    # Add new row to the spearman_correlations_df
    spearman_correlations_df <- rbind(spearman_correlations_df, data.frame(segment = current_segment, 
                                      gene_id = current_gene_id, spearman_corr = spearman_correlation$estimate,
                                      p_value = spearman_correlation$p.value, stringsAsFactors = FALSE))
    # Reset rownames
    rownames(spearman_correlations_df) <- 1:nrow(spearman_correlations_df)

    # Update progress bar
    setTxtProgressBar(pb, i)
  }

  # Close progress bar
  close(pb)

  # Add q_value column to the spearman_correlations_df
  spearman_correlations_df$q_value <- qvalue(spearman_correlations_df$p_value)$qvalues

  return(spearman_correlations_df)
}

#'
#' Function calculates Spearman correlation between a handed over column of the clinical metadata and the atac data
#' represented as rows in the elnet segments atac data
perform_spearman_correlation <- function(output_path, clinical_metadata_rna_ids, clinical_metadata_column, elnet_segments_atac_data) {
  
  # Extract the name of the clinical metadata column
  passed_paramter <- deparse(substitute(clinical_metadata_column))
  parts <- strsplit(passed_paramter, "\\$")[[1]]
  clinical_metadata_column_name <- if(length(parts) > 1) parts[2] else parts[1]

  # Create output direcotry based on the handed over output path and the clinical metadata column name
  output_dir <- paste(output_path, clinical_metadata_column_name, sep = "/")

  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  #
  # Extract clinical metadata column
  #

  # Create df with the handed over column of the clinical metadata ant the RNA ids
  clinical_metadata_column_df <- data.frame( clinical_metadata_rna_ids, clinical_metadata_column)
  
  # Sort rows by patient ids
  clinical_metadata_column_df <- clinical_metadata_column_df[order(clinical_metadata_column_df$clinical_metadata_rna_ids),]

  # Set rownames to the RNA ids
  rownames(clinical_metadata_column_df) <- clinical_metadata_column_df$clinical_metadata_rna_ids

  # Remove RNA ids from the df
  clinical_metadata_column_df <- clinical_metadata_column_df[-1]
  
  # Extract the numeric vector from clinical_metadata_column_df
  # Assuming it has only one column after removing the IDs
  clinical_vector <- as.numeric(clinical_metadata_column_df[, 1])


  #
  # Create Correlation df by looping over the elnet segments with their respective atac data 
  #
  spearman_correlations_df <- create_spearman_correlations_df(clinical_vector = clinical_vector, elnet_segments_atac_data = elnet_segments_atac_data)
  
  #
  # Create output files
  #
  
  # Save spearman_correlations_df as csv file with seperator ;
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "spearman_correlations.csv", sep = "_"), sep = "/")
  print(output_file_path)
  write.table(spearman_correlations_df, output_file_path, row.names = FALSE, sep = ";")
  
  return(spearman_correlations_df)
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

# Print sorted and ; separated final rna ids
#cat(paste("Sorted and ; separated final rna ids:", paste(sort(reduced_clinical_metadata$`rna.final.ids`), collapse = ";")))
# Print sorted and ; separated final patient_ids
#cat(paste("Sorted and ; separated final patient_ids:", paste(sort(patient_ids), collapse = ";")))


# create output_path based on cwd
output_path <- paste(getwd(), "AssociationAnalysis/BasicRegression/combined_correlations", sep = "/")

# Run Regression
spearman_correlations_df_final <- perform_spearman_correlation(output_path = output_path, reduced_clinical_metadata$rna.final.ids,reduced_clinical_metadata$`mPAP`, elnet_segments_atac)

nrow(spearman_correlations_df_final)

# Sort spearman_correlations_df_final by p_value
spearman_correlations_df_final <- spearman_correlations_df_final[order(spearman_correlations_df_final$p_value),]
head(spearman_correlations_df_final)
# Create filtered df, which only rows with p and q value below 0.05
filtered_spearman_correlations_df <- spearman_correlations_df_final[spearman_correlations_df_final$p_value < 0.05 & spearman_correlations_df_final$q_value < 0.05,]

head(filtered_spearman_correlations_df)
