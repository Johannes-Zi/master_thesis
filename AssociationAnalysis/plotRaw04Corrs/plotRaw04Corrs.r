library(enrichplot)
library(clusterProfiler)
library(gprofiler2)
library(org.Hs.eg.db)
library(DOSE)
library(crayon)
library(ggplot2)
library(svglite)
library(dplyr)

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
create_reduced_correlation_datasets <- function(correlation_results_list, correlation_threshold, output_dir) {
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
    message(paste("Number of segmental correlations: ", nrow(correlation_results)))

    # Transform the correlation results into a df
    correlation_results <- as.data.frame(correlation_results)

    # Add colum with the absolute value of the correlation
    correlation_results$absolute_spearman_corr <- abs(correlation_results$spearman_corr)

    # Extract rows with absolute correlation scores above the threshold 
    reduced_correlation_results <- correlation_results[correlation_results$absolute_spearman_corr >= 
                                                       correlation_threshold,]

    num_correlation_results_above_threshold <- nrow(reduced_correlation_results)
    message(paste("Number of segmental correlations above threshold: ", num_correlation_results_above_threshold))

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

    # Filter out noncoding genes
    if (TRUE) {
      message(paste("Number of genes before filtering out noncoding genes: ", nrow(reduced_correlation_results)))
      # Extract gene names
      gene_names <- reduced_correlation_results$gene_id

      translated_GENETYPE_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENETYPE", OrgDb = "org.Hs.eg.db")))

      # Merge the gene types with the reduced correlation results
      reduced_correlation_results <- merge(reduced_correlation_results, translated_GENETYPE_ids, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)

      # Filter out noncoding genes
      reduced_correlation_results <- reduced_correlation_results[reduced_correlation_results$GENETYPE == "protein-coding",]

      # Remove the GENETYPE column from the dataframe - Replaces by NA
      reduced_correlation_results <- subset(reduced_correlation_results, select = -GENETYPE)  
      # Throw out entries with NA in gene type column
      reduced_correlation_results <- reduced_correlation_results[!is.na(reduced_correlation_results$gene_id),]
      # Reset rownames
      rownames(reduced_correlation_results) <- NULL

      message(paste("Number of genes after filtering out noncoding genes: ",nrow(reduced_correlation_results)))
    }

    # Filter out coding-genes that cant be translated to NCBI gene names
    if (TRUE) {
      message(paste("Number of coding genes before filtering out genes that cant be tranlated to ENTRZID: ", nrow(reduced_correlation_results)))
      # Extract gene names
      gene_names <- reduced_correlation_results$gene_id

      translated_ENTREZID_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")))

      # Merge the gene types with the reduced correlation results
      reduced_correlation_results <- merge(reduced_correlation_results, translated_ENTREZID_ids, by.x = "gene_id", by.y = "ENSEMBL", all.x = TRUE)

      # Filter out genes with NA in the ENTRZID column
      reduced_correlation_results <- reduced_correlation_results[!is.na(reduced_correlation_results$ENTREZID),]

      # Remove the GENETYPE column from the dataframe
      reduced_correlation_results <- subset(reduced_correlation_results, select = -ENTREZID)
      message(paste("Number of coding genes before filtering out genes that cant be tranlated to ENTRZID: ",nrow(reduced_correlation_results)))
    }
    message(paste("Number of genes in the final set (unique + above threshold + coding gene): ",
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
  
  print(head(filter_df))
  # Drop the row with clinical parameter "SMW"
  filter_df <- filter_df[filter_df$clinical_parameter != "SMW",]
  print("Paramter SMW dropped")

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
  # Creates a ggplot histogram plot with the gene counts
  ggplot(data = data.frame(gene = names(gene_counts), count = as.numeric(gene_counts)), aes(x = count)) +
    geom_histogram(binwidth = 1) +
    xlab("Number of gene specific segments across all clinical parameters") +
    ylab("Number of genes") +
    ggtitle("Gene count distribution of gene specific segment counts")
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
  # Save the plot as a svg file using file.path
  ggsave(file.path(output_dir, "segments_per_gene_across_all_params.svg"))

  # Creates a ggplot histogram plot with number of genes above threshold for each clinical parameter
  ggplot(filter_df, aes(x = clinical_parameter, y = num_correlation_results_above_threshold)) +
    geom_bar(stat = "identity") +
    xlab("Clinical parameter") +
    ylab("Number of genes above the correlation threshold") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1)) 
  # Creates the output path for the plot
  plot_output_path <- file.path(output_dir, "num_genes_above_threshold_per_clinical_parameter.svg")
  # Save the plot as a svg file based on the parameter_otput_dir
  ggsave(plot_output_path)

  # Creates a ggplot histogram plot with number of genes represented by more than one segment for each clinical parameter
  ggplot(filter_df, aes(x = clinical_parameter, y = num_genes_more_than_one_segment)) +
    geom_bar(stat = "identity") +
    xlab("Clinical parameter") +
    ylab("Number of genes represented by more than one segment") +
    theme(axis.text.x = element_text(angle = 50, hjust = 1))
  # Creates the output path for the plot
  plot_output_path <- file.path(output_dir, "num_genes_more_than_one_segment_per_clinical_parameter.svg")
  # Save the plot as a svg file based on the parameter_otput_dir
  ggsave(plot_output_path)

  # Returns list with reduced correlation datasets
  return(reduced_correlation_dfs_list)
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
    output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/plotRaw04Corrs/runs/11_01_25/"
    correlation_threshold = 0.4  # Defines the correlation threshold to filter out low correlations
    message(bold(cyan("\nCreating reduced correlation datasets")))
    message(paste("Filtered out segments with correlation > ", magenta(correlation_threshold), "or <", magenta(-correlation_threshold)))
    reduced_correlation_datasets_list <- create_reduced_correlation_datasets(correlation_results_list, 
    correlation_threshold = correlation_threshold , output_dir = output_dir)
  }
}