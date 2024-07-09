#library(enrichplot)
#library(gprofiler2)
#library(org.Hs.eg.db)
#library(DOSE)
library(crayon)
library(ggplot2)
library(svglite)


#'
#' Script to run DisGeNET analysis for each clinical parameter
#' 
#' This script imports the clinical metadata and the correlation results of the correlation analysis.
#' The script then filters the correlation results to only include the top X correlations.
#' The script then translates the Ensembl gene ids to NCBI gene names and runs DisGeNET analysis for each clinical 
#' parameter.
#' The script then runs DisGeNET analysis for each clinical parameter and saves the results as a csv file.
#' 

# To be ignored
if (FALSE) {
  # if (!require("BiocManager", quietly = TRUE))
  #     install.packages("BiocManager")

  # BiocManager::install("clusterProfiler")
  #BiocManager::install("org.Hs.eg.db")
  #BiocManager::install("gprofiler2")
  #BiocManager::install("enrichplot")

  #install.packages("fastmap")
}

#'
#' Import clinical metadata
import_clinical_meatadata <- function(file_path) {
  # Import clinical metadata
  clinical_metadata <- read.csv(file_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  return(clinical_metadata)
}

#'
#' Import Regresion results in form of the gene specific correlations
import_correlation_results <- function(toplevel_path, clinical_parameter_names, corr_type){

  # Initialize list to store correlation results
  correlation_results_list <- list()

  message(paste(underline("Import correlation results and filter out rows with NA values")))
  message("Filtered out rows with NA values:")
  # Iterate over all clinical parameters
  for (clinical_parameter in clinical_parameter_names) {
    # Create input path
    correlation_results_path <- paste0(toplevel_path, clinical_parameter, "/", clinical_parameter, "_", corr_type, "_correlations.tsv")

    # Import correlation results
    correlation_results <- read.csv(correlation_results_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)


    # Filter out rows, which have a column entry with "NA"
    rows_before <- nrow(correlation_results)
    correlation_results <- correlation_results[!apply(correlation_results, 1, function(row) any(is.na(row))),]

    # Print how many rows were filtered out due to NA values
    rows_after <- nrow(correlation_results)
    rows_filtered <- rows_before - rows_after
    message(paste(rows_filtered, "Filtered out for", clinical_parameter))

    # Append correlation results to list
    correlation_results_list[[clinical_parameter]] <- correlation_results
  }
  message("\n")

  # Return list with all correlation results of each clinical parameter
  return(correlation_results_list)
}

#'
#' Creates List with reduced datasets with the top correlation of each clinical parameter
create_reduced_correlation_datasets <- function(correlation_results_list, correlation_threshold, n_top_correlations, 
                                                output_dir) {
  # Initialize list to store reduced correlation datasets
  reduced_correlation_dfs_list <- list()

  #' Overview dataframe with filter results of each clinical parameter
  #' Columns: clinical parameter, number of correlation results, number of
  #' correlation results above threshold, number of genes represented by more
  #' than one segment, number of genes in the final set
  filter_df <- data.frame(clinical_parameter = character(),
                          num_correlation_results = numeric(),
                          num_correlation_results_above_threshold = numeric(),
                          num_genes_more_than_one_segment = numeric(),
                          num_genes_final_set = numeric())

  # Dataframe to store all correlation results above the correlation threshold across all clinical parameters
  reduced_correlation_results_all_params <- data.frame()

  # Iterate over all clinical parameters and creae respective output files
  for (clinical_parameter in names(correlation_results_list)) {

    message(paste(underline(clinical_parameter)))
    # Extract current correlation results
    correlation_results <- correlation_results_list[[clinical_parameter]]
    message(paste("Number of correlation results: ", nrow(correlation_results)))

    # Transform the correlation results into a df
    correlation_results <- as.data.frame(correlation_results)

    # Add colum with the absolute value of the correlation
    correlation_results$absolute_spearman_corr <- abs(correlation_results$spearman_corr)

    # Extract rows with absolute correlation scores above the threshold 
    reduced_correlation_results <- correlation_results[correlation_results$absolute_spearman_corr >= 
                                                       correlation_threshold,]

    num_correlation_results_above_threshold <- nrow(reduced_correlation_results)
    message(paste("Number of correlation results above threshold: ", num_correlation_results_above_threshold))

    # Sort the reduced correlation results by the absolute correlation
    reduced_correlation_results <- reduced_correlation_results[order(reduced_correlation_results$absolute_spearman_corr,
                                                                     decreasing = TRUE),]


    # Append reduced correlation results to the overall reduced correlation results df with an additional column for the clinical parameter
    reduced_correlation_results_temp <- reduced_correlation_results
    reduced_correlation_results_temp$clinical_parameter <- clinical_parameter
    reduced_correlation_results_temp <- reduced_correlation_results_temp[,c(ncol(reduced_correlation_results_temp), 1:(ncol(reduced_correlation_results_temp)-1))]
    reduced_correlation_results_all_params <- rbind(reduced_correlation_results_all_params, reduced_correlation_results_temp)


    # Identify rows with duplicated genes (both first occurrence and subsequent duplicates)
    duplicated_genes <- duplicated(reduced_correlation_results$gene) | duplicated(reduced_correlation_results$gene, fromLast = TRUE)
    # Filter to only include duplicated genes
    duplicated_gene_rows <- reduced_correlation_results[duplicated_genes, ]
    # Sort duplicated gene rows by gene name abundace
    duplicated_gene_rows <- duplicated_gene_rows[order(duplicated_gene_rows$gene),]


    # Create parameter specific output directory using file.path for better path handling
    parameter_otput_dir <- file.path(output_dir, clinical_parameter)
    dir.create(parameter_otput_dir, recursive = TRUE, showWarnings = FALSE)

    # Save the duplicated gene rows as a csv file using file.path
    duplicated_genes_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter,
                                                                          "_multiple_represented_genes_segments.csv"))
    write.csv(duplicated_gene_rows, duplicated_genes_output_path, row.names = FALSE)

    # Save a list of the duplicated genes as csv file using file.path
    duplicated_genes_list_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter,
                                                                               "_multiple_represented_duplicated_genes.csv"))
    write.csv(data.frame(gene = unique(duplicated_gene_rows$gene)), duplicated_genes_list_output_path, row.names = FALSE)

    # Count the unique genes among the duplicates
    num_genes_more_than_one_segment <- length(unique(duplicated_gene_rows$gene))
    message(paste("Number of genes represented by more than one segment:", num_genes_more_than_one_segment))

    # Filter out duplication in the gene column - retrieve unique gene ids
    reduced_correlation_results <- reduced_correlation_results[!duplicated(reduced_correlation_results$gene),]

    # Create a reduced correlation df and use only up to the top n_top_correlations absolute correlations
    reduced_correlation_results <- head(reduced_correlation_results, n = min(nrow(reduced_correlation_results), 500))
    message(paste("Number of genes in the final set (unique + above threshold + max top n correlations): ",
                  nrow(reduced_correlation_results), "\n"))

    # Append reduced correlation results df to list
    reduced_correlation_dfs_list[[clinical_parameter]] <- reduced_correlation_results

    # Save the reduced correlation results as a csv file using file.path
    reduced_correlation_output_path <- file.path(parameter_otput_dir, paste0(clinical_parameter,
                                                                             "_reduced_correlation_results.csv"))
    write.csv(reduced_correlation_results, reduced_correlation_output_path, row.names = FALSE)

    # Append filter results to the overview df
    filter_df <- rbind(filter_df, data.frame(clinical_parameter = clinical_parameter,
                                             num_correlation_results = nrow(correlation_results),
                                             num_correlation_results_above_threshold = num_correlation_results_above_threshold,
                                             num_genes_more_than_one_segment = num_genes_more_than_one_segment,
                                             num_genes_final_set = nrow(reduced_correlation_results)))

    # Create log file path
    log_file_path <- file.path(parameter_otput_dir, paste0(clinical_parameter, "_filtering.log"))
    # Create log file with clinical paramter specific information, like the messages above
    sink(log_file_path)
    cat(paste(underline(clinical_parameter), "\n"))
    cat(paste("Number of correlation results: ", nrow(correlation_results), "\n"))
    cat(paste("Number of correlation results above threshold: ", num_correlation_results_above_threshold, "\n"))
    cat(paste("Number of genes represented by more than one segment:", num_genes_more_than_one_segment, "\n"))
    cat(paste("Number of genes in the final set (unique + above threshold + max top n correlations): ",
              nrow(reduced_correlation_results), "\n"))
    sink()
  }

  # Save the reduced correlation results of all clinical parameters as a csv file using file.path
  reduced_correlation_output_path <- file.path(output_dir, "corr_results_above_thres_all_params_combined.csv")
  write.csv(reduced_correlation_results_all_params, reduced_correlation_output_path, row.names = FALSE)

  # Count for each gene how many times it is represented in the reduced correlation results
  gene_counts <- table(reduced_correlation_results_all_params$gene)
  # Sort the gene counts in descending order
  gene_counts <- sort(gene_counts, decreasing = TRUE)
  # Save the gene counts as a csv file using file.path
  gene_counts_output_path <- file.path(output_dir, "number_of_segments_per_gene_across_all_params.csv")
  write.csv(data.frame(gene = names(gene_counts), count = as.numeric(gene_counts)), gene_counts_output_path, 
                       row.names = FALSE)
  # Create a ggplot histogram plot with the gene counts
  ggplot(data = data.frame(gene = names(gene_counts), count = as.numeric(gene_counts)), aes(x = count)) +
    geom_histogram(binwidth = 1) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Save the plot as a svg file using file.path
  ggsave(file.path(output_dir, "segments_per_gene_across_all_params.svg"))
  
  # Create a ggplot histogram plot with number of genes above threshold for each clinical parameter
  ggplot(filter_df, aes(x = clinical_parameter, y = num_correlation_results_above_threshold)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) 
  # Create the output path for the plot
  plot_output_path <- paste0(parameter_otput_dir, "num_genes_above_threshold_per_clinical_parameter.svg")
  # Save the plot as a svg file based on the parameter_otput_dir
  ggsave(plot_output_path)

  # Create a ggplot histogram plot with number of genes represented by more than one segment for each clinical parameter
  ggplot(filter_df, aes(x = clinical_parameter, y = num_genes_more_than_one_segment)) +
    geom_bar(stat = "identity") +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Create the output path for the plot
  plot_output_path <- paste0(parameter_otput_dir, "num_genes_more_than_one_segment_per_clinical_parameter.svg")
  # Save the plot as a svg file based on the parameter_otput_dir
  ggsave(plot_output_path)

  # Return list with reduced correlation datasets
  return(reduced_correlation_dfs_list)
}

#'
#'  Run DisGeNET analysis for each clinical parameter
run_digenet_analysis <- function(reduced_correlation_datasets_list){
  # Initialize list to store DisGeNET results
  disgenet_results_list <- list()

  # Initialize empty df to  store combined DisGeNET results
  combined_disgenet_results_df <- data.frame()

  # Extract clinical parameters
  clinical_parameters <- names(reduced_correlation_datasets_list)

  # Iterate over all clinical parameters
  for (i in 1:length(reduced_correlation_datasets_list)) {
    # Current clinical parameter name
    clinical_parameter <- clinical_parameters[i]

    message("\n##################")
    message(paste("current parameter:", cyan(clinical_parameter)))
    message("##################\n")

    # gene_names <- reduced_correlation_datasets_list[[clinical_parameter]]$gene
    # gene_names_str <- paste(gene_names, collapse = ", ")
    # message(cyan("\nGene list:"), gene_names_str)

    # Current reduced correlation results
    reduced_correlation_results <- reduced_correlation_datasets_list[[clinical_parameter]]

    # Extract gene names
    gene_names <- reduced_correlation_results$gene
    #print(gene_names)

    # Traslate the Ensemble ids to NCBI gene names
    translated_ids = bitr(gene_names, fromType="ENSEMBL", toType="ENTREZID", OrgDb="org.Hs.eg.db")

    # Extract the NCBI gene names
    ncbi_gene_names <- translated_ids$ENTREZID

    #print(length(ncbi_gene_names))

    #x <- gconvert(gene_names, organism = "hsapiens", target="ENTREZGENE_ACC")

    #ncbi_gene_names <- x$target
    #print(ncbi_gene_names)

    # Run DisGeNET analysis
    current_results <- enrichDGN(gene = ncbi_gene_names, qvalueCutoff  = 0.4, pvalueCutoff  = 0.4)
    #       pvalueCutoff  = 0.2,
    #       pAdjustMethod = "BH",
    #       minGSSize     = 5,
    #       maxGSSize     = 500,
    #       qvalueCutoff  = 0.2,
    #       readable      = FALSE)
    message(cyan("number of used gene ids used:"))
    message(paste(length(ncbi_gene_names), "\n"))

    print(current_results)

    # Append DisGeNET results to list
    disgenet_results_list[[clinical_parameter]] <- current_results

    # Extract result df
    current_result_df <- current_results@result

    # Remove rownames and use standard rownames
    rownames(current_result_df) <- NULL

    # Add cluster column and set to clinical parameter
    current_result_df$Cluster <- clinical_parameter

    # Move cluster column to first position
    current_result_df <- current_result_df[,c(ncol(current_result_df), 1:(ncol(current_result_df)-1))]

    # Append DisGeNET result df to combined df
    combined_disgenet_results_df <- rbind(combined_disgenet_results_df, current_result_df)
  }

  # Return list with DisGeNET results
  return(list(disgenet_results_list = disgenet_results_list, combined_disgenet_results_df = combined_disgenet_results_df))
}



#'
#' Main
if (TRUE) {

  #' Import clinical metadata and correlation results
  if (TRUE) {
    message(bold(cyan("Importing clinical metadata and correlation results")))
    # Import clinical metadata and extract column names
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    message(cyan("Import clinical metadata"))
    message(underline("From:\n"), clinical_metadata_path)
    clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)
    column_names <- colnames(clinical_metadata[4:ncol(clinical_metadata)])  # Extract column names
    message(underline("Number of imported clinical parameters:"), length(column_names))
    message(underline("Extracted clinical paramters:\n"), paste(column_names, collapse = ", "), "\n")

    # Import correlation results based on current cwd
    correlation_type <- "spearman"  # Defines the correlation type to be imported from the correlation results
    correlation_data_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/BasicCorrelation/combined_corr_cv_pear_04_thres/"
    correlation_results_list <- import_correlation_results(correlation_data_dir, column_names, corr_type = correlation_type)
    message(cyan("Import correlation results"))
    message(underline("From:\n"), correlation_data_dir)
    message(underline("Imported correlation type: \t"), correlation_type)
    # Concatenate and print the number of imported correlation results for each clinical parameter in a single line
    info_string <- paste(sapply(column_names, function(name) {
      paste(name, ":", nrow(correlation_results_list[[name]]))
    }), collapse = ", ")
    message(underline("Number of imported correlation results (genes) per clinical parameter: \n"), info_string)
    }

  # Create reduced correlation datasets and save results as output
  if (TRUE) {
    output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/disgenet2r/runs/combined_corr_cv_pear_04_thres/spear_thres_05_up_to_500/corr_based_gene_filtering/"
    correlation_threshold = 0.5  # Defines the correlation threshold to filter out low correlations
    n_top_correlations = 500  # Defines the number of top correlations to be included in the reduced correlation datasets
    message(bold(cyan("\nCreating reduced correlation datasets")))
    message(paste("Filtered out genes with correlation > ", magenta(correlation_threshold), "or <", magenta(-correlation_threshold)))
    message(paste("Included up to top", magenta(n_top_correlations), "correlations\n"))
    reduced_correlation_datasets_list <- create_reduced_correlation_datasets(correlation_results_list, 
    correlation_threshold = correlation_threshold , n_top_correlations = n_top_correlations, output_dir = output_dir)
  }

  # Run DisGeNET analysis
  if (FALSE) {
    DisGeNETresults <- run_digenet_analysis(reduced_correlation_datasets_list)
  }

  # Save the DisGeNET results df as csv
  if (FALSE) {
    write.csv(DisGeNETresults$combined_disgenet_results_df, "DisGeNETresults_0.4pqcutoff_uptotop225correlations.csv", row.names = FALSE)
  }
}