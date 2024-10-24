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


# Iteriert über die segmente eine klinischen Parameters
#'
#' Create df with  spearman correlation between all elnet segments and the handed over clinical metadata column
create_spearman_correlations_df <- function(clinical_vector, elnet_segments_atac_data,
    plot_output_dir, clinical_parameter, patient_vector, condition_vector, gene_id_to_visualize, segment_to_visualize) {

  # Initialize progress bar
  #pb <- txtProgressBar(min = 0, max = nrow(elnet_segments_atac_data), style = 3)

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

    # Check if the current gene_id and the current segment are the ones which should be visualized
    if ((current_gene_id != gene_id_to_visualize) || (current_segment != segment_to_visualize)) {
      next
    }
    message(paste("Visualize spearman correlation of", current_gene_id, "at", current_segment))

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
    #setTxtProgressBar(pb, i)

    # Create dataframe with the patient_id vector, the atac_vector and the clinical_vector
    data_df <- data.frame(patient_id = patient_vector, atac_vector = atac_vector, clinical_vector = clinical_vector, condition_vector = condition_vector)

    # Create dotplot to visualize the correlation with larger font and bigger dots
    # Create dotplot to visualize the correlation with larger font and bigger dots
    suppressWarnings(suppressMessages(({
      p <- ggplot(data = data_df, aes(x = atac_vector, y = clinical_vector, color = condition_vector, shape = condition_vector)) +
        geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black", linetype = "dashed", size = 2.5, show.legend = FALSE) +  # Combined smoothing line with confidence interval
        geom_smooth(method = "lm", se = FALSE, size = 2.5, show.legend = FALSE) +  # Group-specific smoothing lines without confidence intervals
        geom_point(size = 6) +  # size of the dots
        labs(title = paste(""),
            x = "ATAC-Seq enhancer signal",
            y = "Systolic pulmonary artery pressure",
            color = "Condition",
            shape = "Condition") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
        theme(
          #plot.title = element_blank(),
          text = element_text(size = 20), 
          plot.title = element_text(size = 40, colour = "#161616", face = "bold"),  #  font size for the title
          axis.title = element_text(size = 40, face = "bold"),  #  font size for axis titles
          axis.text = element_text(size = 34),  #  font size for axis text
          axis.title.x = element_text(margin = margin(t = 15), colour = "#161616"),  #  space above the x-axis title
          axis.title.y = element_text(margin = margin(r = 15), colour = "#161616"),  # Ad space to the right of the y-axis title
          legend.position = c(0.8, 0.9),  # Position the legend inside the plot area
          legend.background = element_rect(fill = alpha('#ffffffed', 0.5)),  # Semi-transparent background for the legend
          legend.text = element_text(size = 38, colour = "#161616"),  #  the size of the legend text
          legend.title = element_text(size = 40, colour = "#161616", face = "bold"),  #  the size of the legend title
          legend.key.size = unit(1.5, "lines")  #  the size of the legend keys
        ) +
        scale_shape_manual(values = c(18, 15, 16, 17)) +  #  shapes for different groups
        scale_color_manual(values = c("#000000", "#009E73", "#0072B2", "#D55E00")) +  #  colors for different groups
        guides(
          color = guide_legend(override.aes = list(shape = c(18, 15, 16, 17))),  #  color and shape legends
          shape = guide_legend(override.aes = list(color = c("#000000", "#009E73", "#0072B2", "#D55E00")))  # Ensure shapes are shown in the legend
        )
    })))

    # Save plot as png with white background
    suppressWarnings(suppressMessages({
      plot_output_file_path_png <- paste(plot_output_dir, paste(paste(spearman_correlation$estimate, current_gene_id, current_segment, "dotplot", sep = "_"), ".svg", sep = ""), sep = "/")
      ggsave(plot_output_file_path_png, plot = p, width = 12, height = 10, dpi = 150, bg = "white")
    }))


    # Create boxplot to display the x-axis values of the different condition groups
    suppressWarnings(suppressMessages({
      boxplot <- ggplot(data = data_df, 
                        aes(x = factor(condition_vector, levels = c("healthy", "ph-lung", "pah", "cteph")),  # Manually set the order of the groups
                            y = atac_vector, 
                            color = condition_vector, 
                            fill = condition_vector)) +
        geom_boxplot(outlier.shape = NA, alpha = 0.5, width = 0.8, size=1.5, color = "#818181") +  # Add boxplot with semi-transparent fill and no outliers
        geom_jitter(shape = 16, position = position_jitter(height = 0, width = 0.15), size = 6) +  # Add jittered points for better visibility
        labs(title = paste(""),
            x = "Condition",
            y = "ATAC-Seq enhancer signal",
            color = "Condition",
            fill = "Condition") +
        theme_minimal() +
        theme(plot.title = element_text(hjust = 0.5)) +  # Center the title
        theme(
          text = element_text(size = 20),  # Increase font size for all text
          plot.title = element_text(size = 40, colour = "#161616", face = "bold"),  # Increase font size for the title
          axis.title = element_text(size = 40, face = "bold"),  # Increase font size for axis titles
          axis.text = element_text(size = 34),  # Increase font size for axis text
          axis.title.x = element_text(margin = margin(t = 15), colour = "#161616"),  # Add space above the x-axis title
          axis.title.y = element_text(margin = margin(r = 15), colour = "#161616"),  # Add space to the right of the y-axis title
          axis.text.x = element_text(angle = 35, vjust = 1, hjust = 1),
          legend.position = "none"  # Remove the legend
        ) +
        scale_color_manual(values = c("#000000", "#009E73", "#0072B2", "#D55E00")) +  # Specify colors for different groups
        scale_fill_manual(values = c("#000000", "#009E73", "#0072B2", "#D55E00"))  # Specify fill colors for different groups
    }))

    # Save boxplot as svg with white background
    suppressWarnings(suppressMessages({
      boxplot_output_file_path_svg <- paste(plot_output_dir, paste(paste(spearman_correlation$estimate, current_gene_id, current_segment, "boxplot", sep = "_"), ".svg", sep = ""), sep = "/")
      ggsave(boxplot_output_file_path_svg, plot = boxplot, width = 6, height = 12, dpi = 150, bg = "white")
    }))

    # Check if the current gene_id and the current segment are the ones which should be visualized and thus the loop can be stopped
    if (current_gene_id == gene_id_to_visualize && current_segment == segment_to_visualize) {
      break
    }
  }

  # Sort files for windows by opening them enabling sort by date option
  process_and_sort_files_in_directory(plot_output_dir)

  # Close progress bar
  #close(pb)

  # Add q_value column to the spearman_correlations_df
  #spearman_correlations_df$q_value <- qvalue(p = spearman_correlations_df$p_value)$qvalues

  return(spearman_correlations_df)
}


# Erstellt directory und segment spezifische correlation visualisierungen
#'
#' Function calculates Spearman correlation between a handed over column of the clinical metadata and the atac data
#' represented as rows in the elnet segments atac data
perform_spearman_correlation <- function(output_path, clinical_metadata_rna_ids, clinical_metadata_column, elnet_segments_atac_data, clinical_metadata_column_name, condition_vector, gene_id_to_visualize, segment_to_visualize) {

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

  # Create vector with the patient ids in form of the rownames
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
  spearman_correlations_df <- create_spearman_correlations_df(
    clinical_vector = clinical_vector,
    elnet_segments_atac_data = elnet_segments_atac_data,
    plot_output_dir = plot_output_dir,
    clinical_parameter = clinical_metadata_column_name,
    patient_vector = patient_vector, condition_vector = condition_vector,
    gene_id_to_visualize = gene_id_to_visualize,
    segment_to_visualize = segment_to_visualize)
  
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
iterate_over_clinical_metadata <- function(output_path, clinical_metadata, elnet_segments_atac_data, clinical_parameter_to_visualize, gene_id_to_visualize, segment_to_visualize) {
  
  # Initialize empty list to store the spearman_correlations_df_final of each clinical metadata column
  spearman_correlations_df_final_list <- list()

  # Iterate over columns of the clinical metadata and perform spearman correlation and save output files
  for (i in 4:ncol(clinical_metadata)){

    # Extract current column
    current_column <- clinical_metadata[,i]

    # Check if the current column is the column which should be visualized
    if (colnames(clinical_metadata)[i] != clinical_parameter_to_visualize) {
      next
    }
    print(paste("Perform spearman correlation for: ", colnames(clinical_metadata)[i]))

    # Create temporary sorted df based on patient ids to hand over the condition vector in the rigth order
    temp_df <- data.frame(clinical_metadata$rna.final.ids, clinical_metadata$condition)
    temp_df <- temp_df[order(clinical_metadata$rna.final.ids),]
    condition_vector <- temp_df$clinical_metadata.condition

    # Perform spearman correlation
    message(paste(bold(cyan("\nPerform spearman correlations for:")), magenta(colnames(clinical_metadata)[i])))
    current_spearman_correlations_df_final <- perform_spearman_correlation(output_path = output_path,
        clinical_metadata$rna.final.ids, current_column, elnet_segments_atac_data,
        clinical_metadata_column_name = colnames(clinical_metadata)[i], condition_vector, 
        gene_id_to_visualize = gene_id_to_visualize, segment_to_visualize = segment_to_visualize)
    
    # Append current spearman_correlations_df_final to the list
    spearman_correlations_df_final_list[[i]] <- current_spearman_correlations_df_final

    # Check if the current column is the column which should be visualized and thus the loop can be stopped
    if (colnames(clinical_metadata)[i] == clinical_parameter_to_visualize) {
      break
    }
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
  #' Create a list of genes which have a CV Pearson correlation above a defined threshold
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
  if (TRUE) {
    cat(paste(bold(cyan("\nRun the correlation analysis with the filtered elnet segments atac data together with the clinical metadata over all clinical parameters.\n"))))
    # Create output_path based on cwd
    output_path <- paste(getwd(), "AssociationAnalysis/CorrelationVisualization/runs/v2_single_plots", sep = "/")



    # Configuration for the plot that should be visualized
    clinical_parameter_to_visualize <- "sPAP.excl..ZVD"
    # gene_id_to_visualize <- "ENSG00000106991"
    # segment_to_visualize <- "chr9.127853598.127853857"
    gene_id_to_visualize <- "ENSG00000139567"
    segment_to_visualize <- "chr12.51907898.51908147"



    # Run Correlations with filtered elnet segments atac data
    list_with_all_correlation_dfs_filtered <- iterate_over_clinical_metadata(output_path = output_path, 
        reduced_clinical_metadata, elnet_segments_atac_filtered, 
        clinical_parameter_to_visualize = clinical_parameter_to_visualize, gene_id_to_visualize = gene_id_to_visualize,
        segment_to_visualize = segment_to_visualize)
  }

  # Stop logging
  #sink()
}
