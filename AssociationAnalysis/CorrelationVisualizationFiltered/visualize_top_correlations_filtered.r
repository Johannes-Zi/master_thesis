library(qvalue)
library(ggplot2)
library(crayon)
library(dplyr)


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
#' Load ENSEBL gene ids of top correlations between clinical parameteres and segmental ATAC valeus for all clinical 
#' parameters.
load_filtered_correlations_gene_ids <- function(filtered_correlations_top_directory_path) {
  # Extract directories in the filtered_correlations_top_directory_path
  directories <- list.dirs(filtered_correlations_top_directory_path, full.names = FALSE, recursive = FALSE)

  # Iterate over the directories and load the results.csv files
  filtered_correlations_gene_ids <- c()
  for (directory in directories) {
    # Load the results.csv file
    results_file_path <- paste(filtered_correlations_top_directory_path, paste(directory, paste(directory, "_reduced_correlation_results.csv", sep = ""), sep = "/"), sep = "")
    results_df <- read.csv(results_file_path, header = TRUE, sep = ",", stringsAsFactors = FALSE)

    # Extract the gene ids
    gene_ids <- results_df$gene_id

    # Add the gene ids to the filtered_correlations_gene_ids current vector
    filtered_correlations_gene_ids <- c(filtered_correlations_gene_ids, gene_ids)

    # Drop the duplicated gene ids
    filtered_correlations_gene_ids <- unique(filtered_correlations_gene_ids)
  }

  # Return the filtered_correlations_gene_ids vector
  return(filtered_correlations_gene_ids)
}


#' Function to read, sort, open, and close files in a directory because Windows sort is dumb and i need to use sort
#' by creation date
process_and_sort_files_in_directory <- function(directory_path) {
  # List all files in the directory
  filenames <- list.files(directory_path, full.names = TRUE)
  
  # Sort filenames by name
  sorted_filenames <- sort(filenames)
  
  # Loop through each file
  for (filename in sorted_filenames) {
    # Rename the file by adding a leading "X" to the filename
    new_filename <- file.path(dirname(filename), paste0("X", basename(filename)))

    # Copy the file to the new filename
    file.copy(filename, new_filename)
    
    # Delete the original file
    file.remove(filename)
  }
}

#'
#' Function to perform the permutation test by removing different random subsets
outlier_removal_permutation_test <- function(x, y, x_no_outliers, y_no_outliers, n_permutations = 1000) {
  # Observed difference between input dataset correlations
  cor_diff_observed <- abs(cor(x, y, method = "spearman") - cor(x_no_outliers, y_no_outliers, method = "spearman"))
  
  # Permutation test
  perm_diffs <- replicate(n_permutations, {
    # Identify original outliers
    original_outliers <- setdiff(x, x_no_outliers)

    # Permutation test uses the same number of outliers as the observed number of outliers
    # Otherwise the distribution would be driven by a different number of outliers,
    # which could lead to a over or underestimation of the p-value thus significance
    # of the outlier influence on the correlation
    number_of_outliers <- length(original_outliers)

    # Randomly select an outlier set of the selected size
    selected_outlier_indices <- sample(length(x), number_of_outliers)
    
    # Create a reduced dataset by removing the selected outliers
    reduced_x <- x[-selected_outlier_indices]
    reduced_y <- y[-selected_outlier_indices]
    
    # Recalculate correlations for the reduced dataset
    cor_diff <- cor(reduced_x, reduced_y, method = "spearman") - cor(x_no_outliers, y_no_outliers, method = "spearman")
    return(cor_diff)
  })
  
  # Calculate p-value based on permutation distribution
  permutated_differences <- perm_diffs

  # Calculate p-value
  p_value_perm <- mean(abs(permutated_differences) >= cor_diff_observed)

  if (!is.na(p_value_perm) && p_value_perm < 0.05) {
    print("##### significant influnce on correlation!!! #####")

    print("p_value_perm:")
    print(p_value_perm)

    print("cor_diff_observed:")
    print(cor_diff_observed)

    # Create a ggplot2 that shows the distribution of the differences
    ggplot(data = data.frame(differences = permutated_differences), aes(x = differences)) +
      geom_histogram(binwidth = 0.01, fill = "blue", color = "black") +
      geom_vline(xintercept = cor_diff_observed, color = "red") +
      labs(title = "Permutation test for outlier removal", x = "Difference in Spearman correlation", y = "Frequency") +
      theme_minimal()

    # Save the plot as a png file
    ggsave("permutation_test.png", width = 10, height = 10, dpi = 150, bg = "white")
  }

  return(p_value_perm)
}


# Iteriert über die segmente eine klinischen Parameters
#'
#' Create df with  spearman correlation between all elnet segments and the handed over clinical metadata column
create_spearman_correlations_df <- function(clinical_vector, elnet_segments_atac_data, plot_output_dir, clinical_parameter, patient_vector, condition_vector){

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

    # Calculate spearman correlation between the two vectors
    spearman_correlation <- suppressWarnings(cor.test(clinical_vector, atac_vector, method = "spearman"))

    # Update progress bar
    setTxtProgressBar(pb, i)

    # Create dataframe with the patient_id vector, the atac_vector and the clinical_vector
    data_df <- data.frame(patient_id = patient_vector, atac_vector = atac_vector, clinical_vector = clinical_vector, condition_vector = condition_vector)

    #
    # First filter - Skip segment if the spearman correlation is below abs(0.4)
    #

    # Bug blocker for the case that the correlation is NA
    if (is.na(spearman_correlation$estimate)) {
      next
    }

    # Skip to next segment if the amount of the spearman correlation is below 0.4
    if (abs(spearman_correlation$estimate) < 0.4) {
      next
    }

    #
    # Second filter - check if two third of the ATAC values are not zero
    #

    # Skip if more than one third of the ATAC values are zero
    if (sum(atac_vector == 0) > length(atac_vector) / 3) {
      next
    }

    #
    # Third filter - check if the difference between the healthy and disease group is significant
    #

    # Extract the atac_vector values wich have a cooresponding condition_vector value of cteph, pah, pah-lung
    pah_atac_vector <- data_df[data_df$condition_vector == "pah",]$atac_vector
    pah_lung_atac_vector <- data_df[data_df$condition_vector == "ph-lung",]$atac_vector
    cteph_atac_vector <- data_df[data_df$condition_vector == "cteph",]$atac_vector
    disease_atac_vector <- c(pah_atac_vector, pah_lung_atac_vector, cteph_atac_vector)

    # Extract healthy atac_vector
    healthy_atac_vector <- data_df[data_df$condition_vector == "healthy",]$atac_vector

    # Calculate paired t-test between the healthy and disease group
    wilcox_test <- wilcox.test(healthy_atac_vector, disease_atac_vector, exact = FALSE)

    # Skip to next segment if the p-value of the paired t-test is above 0.05
    if (wilcox_test$p.value > 0.05) {
      next
    }

    #
    # Fourth Filter - check if the correlation is driven by outliers
    #

    # Calculate IQR and outlier bounds
    Q1 <- quantile(atac_vector, 0.25)
    Q3 <- quantile(atac_vector, 0.75)
    IQR_value <- IQR(atac_vector)

    # Outlier thresholds
    lower_bound <- Q1 - 1.5 * IQR_value
    upper_bound <- Q3 + 1.5 * IQR_value

    # Identify outliers
    outliers <- atac_vector[atac_vector < lower_bound | atac_vector > upper_bound]

    # If there are outliers, perform permutation test to check if the influence of the outliers on the correlation 
    # is significant
    if (length(outliers) > 0) {
      # message("\nOutliers detected for geneid: ", current_gene_id, " Segment: ", current_segment)
      # message("Outliers: ")
      # print(outliers)
      # print("")
      # Determine the indices of the outliers
      outlier_indices <- which(atac_vector %in% outliers)
      
      # Remove outliers from atac_vector based on the indices
      atac_vector_no_outliers <- atac_vector[-outlier_indices]
      # Remove corresponding values from clinical_vector
      clinical_vector_no_outliers <- clinical_vector[-outlier_indices]

      # Calculate correlation of reduced datasets
      spearman_correlation_no_outliers <- suppressWarnings(cor.test(clinical_vector_no_outliers, atac_vector_no_outliers, method = "spearman"))

      p_value_perm <- outlier_removal_permutation_test(atac_vector, clinical_vector, atac_vector_no_outliers, clinical_vector_no_outliers)

      # message("Correlation without outliers: ", spearman_correlation_no_outliers$estimate)
      # message("Correlation with outliers: ", spearman_correlation$estimate)
      # message("P-value of permutation test: ", p_value_perm)

    } else {
      p_value_perm <- NA
    }

    # Skip to next segment if the p-value suggests that the corrleation is outlier driven and the correlation with the reduced dataset is below abs(0.4)
    if (!is.na(p_value_perm) && p_value_perm < 0.05) {
      # Skip only if the even though the p-value is significant the correlation is higher than |0.4|
      if (abs(spearman_correlation_no_outliers$estimate) < 0.4) {
        message("Geneid: ", current_gene_id, " Segment: ", current_segment, " is outlier driven and the outlier excluded correlation is below |0.4|")
        next
      }
    }

    # Add new row to the spearman_correlations_df if all filters passed
    spearman_correlations_df <- rbind(spearman_correlations_df, data.frame(segment = current_segment, 
                                      gene_id = current_gene_id, spearman_corr = spearman_correlation$estimate,
                                      p_value = spearman_correlation$p.value, stringsAsFactors = FALSE))
    # Reset rownames
    rownames(spearman_correlations_df) <- 1:nrow(spearman_correlations_df)

    # Create dotplot to visualize the correlation
    suppressWarnings(suppressMessages(({
      p <- ggplot(data = data_df, aes(x = atac_vector, y = clinical_vector, color = condition_vector)) +
        geom_point() +
        geom_smooth(aes(x = atac_vector, y = clinical_vector), method = "lm", se = TRUE, color = "black") +
        geom_smooth(method = "lm", se = FALSE) +
        labs(title = paste("Correlation between ATAC-seq signal and", clinical_parameter, "for", current_gene_id, "in",
                           current_segment), 
            x = "ATAC-seq signal",
            y = clinical_parameter) +
        theme_minimal()
    })))

    # Save plot as png with white background
    suppressWarnings(suppressMessages({
      plot_output_file_path_png <- paste(plot_output_dir, paste(paste(spearman_correlation$estimate, current_gene_id, current_segment, sep = "_"), ".png", sep = ""), sep = "/")
      ggsave(plot_output_file_path_png, plot = p, width = 10, height = 10, dpi = 150, bg = "white")
    }))
  }

  # Sort files for windows by opening them enabling sort by date option
  process_and_sort_files_in_directory(plot_output_dir)

  # Close progress bar
  close(pb)

  # Add q_value column to the spearman_correlations_df
  #spearman_correlations_df$q_value <- qvalue(p = spearman_correlations_df$p_value)$qvalues

  return(spearman_correlations_df)
}


# Erstellt directory und segment spezifische correlation visualisierungen
#'
#' Function calculates Spearman correlation between a handed over column of the clinical metadata and the atac data
#' represented as rows in the elnet segments atac data
perform_spearman_correlation <- function(output_path, clinical_metadata_rna_ids, clinical_metadata_column, elnet_segments_atac_data, clinical_metadata_column_name, condition_vector) {

  # Create output direcotry based on the handed over output path and the clinical metadata column name
  output_dir <- paste(output_path, clinical_metadata_column_name, sep = "/")

  # Create output directory if it does not exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Create plot output directory within the output directory if it does not exist
  plot_output_dir <- paste(output_dir, "plots/", sep = "/")
  if (!dir.exists(plot_output_dir)) {
    dir.create(plot_output_dir, recursive = TRUE)
  }

  #
  # Extract clinical metadata column
  #
  #

  # Combine data and create df with the handed over column of the clinical metadata and the RNA ids
  clinical_metadata_column_df <- data.frame(clinical_metadata_rna_ids, clinical_metadata_column)

  # Sort rows by patient ids
  clinical_metadata_column_df <- clinical_metadata_column_df[order(clinical_metadata_column_df$clinical_metadata_rna_ids),]

  # Set rownames to the RNA ids
  rownames(clinical_metadata_column_df) <- clinical_metadata_column_df$clinical_metadata_rna_ids

  # Create vector with the patient ids in frm of the rownames
  patient_vector <- rownames(clinical_metadata_column_df)

  # Remove RNA ids from the df
  clinical_metadata_column_df <- clinical_metadata_column_df[-1]

  # Extract the numeric vector from clinical_metadata_column_df
  # Assuming it has only one column after removing the IDs
  # Change comma to dot before making it numeric
  clinical_vector <- as.numeric(gsub(",", ".", clinical_metadata_column_df[, 1]))

  #
  # Create Correlation df by looping over the elnet segments with their respective atac data 
  #
  spearman_correlations_df <- create_spearman_correlations_df(clinical_vector = clinical_vector, 
                                elnet_segments_atac_data = elnet_segments_atac_data, plot_output_dir = plot_output_dir,
                                clinical_parameter = clinical_metadata_column_name, patient_vector = patient_vector, 
                                condition_vector = condition_vector)
  
  #
  # Create output files
  #
  
  # Save spearman_correlations_df as csv file with seperator ;
  output_file_path <- paste(output_dir, paste(clinical_metadata_column_name, "spearman_correlations.tsv", sep = "_"), sep = "/")
  write.table(spearman_correlations_df, output_file_path, row.names = FALSE, sep = ";")
  print(paste("->", output_file_path))

  return(spearman_correlations_df)
}


# Iteriert über die klinischen Parameter
#'
#' Iterate over all columns of the clinical metadata and perform the spearman correlations
iterate_over_clinical_metadata <- function(output_path, clinical_metadata, elnet_segments_atac_data) {

  # Initialize empty list to store the spearman_correlations_df_final of each clinical metadata column
  spearman_correlations_df_final_list <- list()

  # Iterate over columns of the clinical metadata and perform spearman correlation and save output files
  for (i in 4:ncol(clinical_metadata)){

    # Extract current column
    current_column <- clinical_metadata[,i]


    # Create temporary sorted df based on patient ids to hand over the condition vector in the rigth order
    temp_df <- data.frame(clinical_metadata$rna.final.ids, clinical_metadata$condition)
    temp_df <- temp_df[order(clinical_metadata$rna.final.ids),]
    condition_vector <- temp_df$clinical_metadata.condition

    # Perform spearman correlation
    message(paste(bold(cyan("\nPerform spearman correlations for:")), magenta(colnames(clinical_metadata)[i])))
    current_spearman_correlations_df_final <- perform_spearman_correlation(output_path = output_path, clinical_metadata$rna.final.ids, 
                                                                           current_column, elnet_segments_atac_data, clinical_metadata_column_name=colnames(clinical_metadata)[i], condition_vector)
    
    # Append current spearman_correlations_df_final to the list
    spearman_correlations_df_final_list[[i]] <- current_spearman_correlations_df_final
  }

  # Return list with all spearman_correlations_df_final
  return(spearman_correlations_df_final_list)
}


#'
#' Load the gene model specific outer cross validation performance estimations
load_performance_evaluation <- function(LeaveOneOutCV_performance_path) {
  # Load the data
  LeaveOneOutCV_performance <- read.table(LeaveOneOutCV_performance_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  
  # Add column with gene name based on the segment name
  LeaveOneOutCV_performance$gene_name <- sapply(strsplit(LeaveOneOutCV_performance$Sample_Name, "_"), function(x) x[2])

  # Drop rows of the column Sample_Name, if they do not contain the word Pearson
  LeaveOneOutCV_performance <- LeaveOneOutCV_performance[grepl("Pearson", LeaveOneOutCV_performance$Sample_Name),]

  return(LeaveOneOutCV_performance)
}




#' 
#' Main
#' 
if (TRUE) {
  #'
  #' Start logging
  #' 
  # Path to the log file
  #log_file_path <- paste(getwd(), "AssociationAnalysis/BasicCorrelation/combined_correlations", "log.txt", sep = "/")
  #sink(log_file_path, append = FALSE, split = TRUE)

  #'
  #' Import required data
  if (TRUE) {
    cat(paste(bold(cyan("\nImport required data\n"))))
    # Import clinical metadata
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)
    message(paste("Imported clinical metadata from: \n", clinical_metadata_path))

    # Import elnet segments atac data
    elnet_segments_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/CollapsedSegmentation/elnet_model_segments_df.tsv"
    elnet_segments_atac <- import_elnet_segments_atac_data(elnet_segments_path)
    message(paste("Imported elnet segments atac data from: \n", elnet_segments_path))

    # Import gene model specific outer cross validation performance estimations
    LeaveOneOutCV_performance_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/regression/leaveOneOut_regression/performance_evaluation/Performance_Overview.tsv"
    LeaveOneOutCV_performance <- load_performance_evaluation(LeaveOneOutCV_performance_path)  # Load the data
    message(paste("Imported gene model specific outer cross validation performance estimations from: \n", LeaveOneOutCV_performance_path))

    # Import filtered correlations gene ids (DisGeNet/ GO-enrichment input set)
    filtered_correlations_top_directory_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/go_enrichment/spear_thres_04_up_to_250/corr_based_gene_filtering/"
    filtered_correlations_gene_ids <- load_filtered_correlations_gene_ids(filtered_correlations_top_directory_path)
    message(paste("Imported filtered correlations gene ids from: \n", filtered_correlations_top_directory_path))
    message(paste("Total numer of unique gene ids extracted from all prefiltered clinical parameter gene sets: ", magenta(length(filtered_correlations_gene_ids))))
  }

  #'
  #' Reduce clinical metadata to the patient ids which are represented in the elnet segments atac data
  if (TRUE) {
    message(paste(bold(cyan("\nReduce clinical metadata to the patient ids which are represented in the elnet segments atac data"))))
    # Extract represented patient ids in the elnet segments atac data
    patient_ids <- extract_represented_patient_ids(elnet_segments_atac)
    message(paste("Extracted represented patient ids in the elnet segments atac data: ", magenta(length(patient_ids))))

    # Print how many unique patient ids are represented in the clinical metadata with RNA Seq ids
    message(paste("Number of patients represented in the clinical metadata with RNA Seq ids: ", magenta(length(unique(clinical_metadata$rna.final.ids)))))
    # Reduce clinical metadata to the patient ids which are represented in the elnet segments atac data
    reduced_clinical_metadata <- reduced_clinical_metadata(clinical_metadata, patient_ids)
    message(paste("Number of patients represented in the reduced clinical metadata: ", magenta(nrow(reduced_clinical_metadata))))
  }

  #'
  #' Create a list of genes which have a cross validation Pearson correlation above a defined threshold
  if (TRUE) {
    thresx <- 0.4 # Defines the threshold for the outer cross validation threshold to filter out the gene models below the threshold
    cat(paste(bold(cyan("\nApply elnet outer cross validation based pearson correlation filtering\n"))))
    cat(paste( "with threshold Pearson >=", red(thresx), "or <=", red(-thresx), "\n"))
    cat(paste("Number of unique genes represented in LeaveOneOutCV_performance: ", magenta(length(unique(LeaveOneOutCV_performance$gene_name)), "\n")))
    cat(paste("Number of entries with negative Pearson correlation: ", magenta(nrow(LeaveOneOutCV_performance[LeaveOneOutCV_performance$Pearson < 0,]), "\n")))
    # Filter out the genes with a amount of Pearson correlation below the threshold
    genes_with_cv_pearson_correlation_over_thresx <- LeaveOneOutCV_performance[LeaveOneOutCV_performance$Pearson >= thresx | LeaveOneOutCV_performance$Pearson <= -thresx, ]$gene_name
    cat(paste("Number of genes after filtering: ", magenta(length(genes_with_cv_pearson_correlation_over_thresx), "\n")))
  }

  #'
  #' Filter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list
  if (TRUE) {
    cat(paste(bold(cyan("\nFilter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list.\n"))))
    cat(paste("Number of unique genes represented by elnet segments: ", magenta(length(unique(elnet_segments_atac$gene_id)), "\n")))
    cat(paste("Number of elnet segments before filtering: ", magenta(nrow(elnet_segments_atac), "\n")))
    # Filter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list
    elnet_segments_atac_filtered <- elnet_segments_atac[elnet_segments_atac$gene_id %in% genes_with_cv_pearson_correlation_over_thresx,]
    cat(paste("Number of elnet segments after filtering: ", magenta(nrow(elnet_segments_atac_filtered), "\n")))
  }

  #'
  #' Filter the elnet segments atac data based on the filtered correlations gene ids (gene ids represented in the DisGeNet or GO-enrichment input set)
  #' for a reduced correlation analysis set that can be automatically visualized
  if (TRUE) {
    cat(paste(bold(cyan("\nFilter the elnet segments down to the genes represented in the GO enrichment and DisGeNet input fo a comprehensive visualization.\n"))))
    cat(paste("Number of unique genes represented by elnet segments: ", magenta(length(unique(elnet_segments_atac$gene_id)), "\n")))
    cat(paste("Numebr of genes represented in the merged GO enrichmen clinical parameter gene sets: ", magenta(length(filtered_correlations_gene_ids), "\n")))
    cat(paste("Number of elnet segments before filtering: ", magenta(nrow(elnet_segments_atac), "\n")))
    # Filter the elnet segments atac data based on the elnet outer cv Pearson filtered gene list
    elnet_segments_atac_filtered <- elnet_segments_atac[elnet_segments_atac$gene_id %in% filtered_correlations_gene_ids,]
    cat(paste("Number of elnet segments after filtering: ", magenta(nrow(elnet_segments_atac_filtered), "\n")))
  }

  #'
  #' Run the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over 
  #' all clinical parameters
  #' 
  #' Additional segment specific filter steps implemented within the loop
  if (TRUE) {
    cat(paste(bold(cyan("\nRun the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over all clinical parameters.\n"))))
    # Create output_path based on cwd
    output_path <- paste(getwd(), "AssociationAnalysis/CorrelationVisualizationFiltered/runs/v1", sep = "/")
    # Run Correlations with filtered elnet segments atac data
    list_with_all_correlation_dfs_filtered <- iterate_over_clinical_metadata(output_path = output_path, reduced_clinical_metadata, elnet_segments_atac_filtered)
  }

  # Stop logging
  #sink()
}
