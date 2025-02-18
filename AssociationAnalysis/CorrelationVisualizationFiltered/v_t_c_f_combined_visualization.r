library(ggplot2)
library(dplyr)
library(tidyr)

# Function that loads the data from the input directories
load_data <- function(input_dir) {

  # Dataframe that hold all the input data
  combined_df <- data.frame()

  # Iterate over directories in input_dir
  for (dir in list.dirs(input_dir, full.names = FALSE, recursive = FALSE)) {
    # Set directory name as the name of the current clinical parameter
    clinical_parameter <- dir

    # Build the path to the input file
    input_file <- paste(input_dir, dir, "/", dir, "_spearman_correlations.tsv", sep = "")

    # Load the input tsv file to a data frame
    input_df <- read.table(input_file, header = TRUE, sep = ";")

    # Add the clinical parameter to the data frame
    input_df$clinical_parameter <- clinical_parameter

    # Add the data frame to the combined data frame
    combined_df <- rbind(combined_df, input_df)
  }
  # Rename parameters
  combined_df$clinical_parameter <- gsub("age", "age", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("\\bAT\\b", "pulmonary acceleration time", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("AT_ET", "AT / ET", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("AVDO2", "arteriovenous oxygen difference", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("BNP", "BNP - heart failure marker", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Cardiac.Index", "caridac performance index", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("CVP", "central venous pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("DLCO", "lung CO diffusion capacity", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("eGFR", "glomerular filtration rate", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("heart.rate", "heart rate", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("height", "height", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("HZV_Fick", "cardiac output (Fick principle)", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("HZV_Thermodil", "cardiac output (thermodilution)", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Kreatinin", "kreatinin", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("mPAP", "mean Pulmonary Artery Pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("NYHA", "classification of heart failure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PA_diastolic", "diastolic PAP", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PA_systolic", "systolic PAP", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Paradoxe_Septumbewegung", "ventricular septum movement", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PAWP", "pulmonary capillary occlusion pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Perikarderguss", "pericardium fluid accumulation", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("PVR", "pulmonary vascular resistance", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("RA_area", "right atrium size", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Rrsys", "systolic blood pressure at rest", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("RVEDD", "right ventricle size at diastole", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("\\bS\\b", "impaired ventricular contraction TAPSV", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("sPAP.excl..ZVD", "systolic PAP without CVP-est.", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("TAPSE", "impaired ventricular contraction TAPSE", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("Tiffeneau.Index.FEV1.VC", "airflow obstruction in lung", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("ZVD", "central venous pressure", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("VCI_diameter", "inferior vena cava diameter", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("ven_SO2", "venous oxygen saturation", combined_df$clinical_parameter)
  combined_df$clinical_parameter <- gsub("\\bFEV1\\b", "one second forced expiratory volume", combined_df$clinical_parameter)

  # Legende hinzufügen:
  # AT / ET: (pulmonary acceleration time) / (right ventricular ejection time)
  # measured by RHC link
  return(combined_df)
}


# Function that plots the number of df entries for each clinical parameter
plot_entries <- function(input_df, output_dir) {
  # Plot the number of entries for each clinical parameter
  ggplot(input_df, aes(x = clinical_parameter)) +
    geom_bar() +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 40, hjust = 1),  # Rotate x-axis labels by 45 degrees
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.background = element_rect(fill = "white", color = NA)  # Add white background to the plot
    ) +
    labs(title = "Number of top segments for each clinical parameter",
         x = "Clinical parameter",
         y = "Number of segments")

  # Save the plot to a file 
  ggsave(paste(output_dir, "number_of_entries_per_param.png", sep = ""))

  # Save the number of entries for each clinical parameter to a csv file
  write.table(as.data.frame(table(input_df$clinical_parameter)), paste(output_dir, "number_of_entries_per_param.csv", sep = ""), sep = ",", row.names = FALSE)
}

# Function that creates a barchart with the most abundand gene_id entries across all clinical parameters
plot_top_gene_ids <- function(input_df, output_dir) {
  # Create a data frame with the number of entries for each gene_id
  gene_id_counts <- as.data.frame(table(input_df$gene_id))

  # Sort the data frame by the number of entries
  gene_id_counts <- gene_id_counts[order(-gene_id_counts$Freq), ]

  # Reset the row names
  rownames(gene_id_counts) <- NULL
  
  # Ensure the Var1 column is a factor with levels ordered by Freq
  gene_id_counts$Var1 <- factor(gene_id_counts$Var1, levels = gene_id_counts$Var1[order(gene_id_counts$Freq, decreasing = TRUE)])

  # Plot the Freqency of each gene_id as a barchart
  ggplot(gene_id_counts[1:40, ], aes(y = Var1, x = Freq)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.background = element_rect(fill = "white", color = NA)  # Add white background to the plot
    ) +
    labs(
        title = "Top 10 most abundant gene_ids",
        y = "Gene ID",
        x = "Number of entries"
    )

  # Save the plot to a file
  ggsave(paste(output_dir, "top_gene_ids.png", sep = ""))


  # Calculate the distribution of frequencies
  frequency_distribution <- as.data.frame(table(gene_id_counts$Freq))
  colnames(frequency_distribution) <- c("Frequency", "Number_of_Genes")

  # Convert Frequency to numeric for proper plotting
  frequency_distribution$Frequency <- as.numeric(as.character(frequency_distribution$Frequency))

  # Plot the distribution of frequencies
  ggplot(frequency_distribution, aes(x = Frequency, y = Number_of_Genes)) +
    geom_bar(stat = "identity") +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.background = element_rect(fill = "white", color = NA),  # Add white background to the plot
    ) +
    labs(
      title = "Distribution of Gene occurrences",
      x = "number of clinical parameters a gene is represented as top correlation",
      y = "Number of Genes"
    )

  # Save the plot to a file
  ggsave(paste(output_dir, "gene_frequency_distribution.png", sep = ""))

  # Save the gene_id_counts to a file
  write.table(gene_id_counts, paste(output_dir, "gene_id_counts.tsv", sep = ""), sep = "\t", row.names = FALSE)
  
  # Print the gene_id_counts to the console comma separated as a single line
  print(paste(gene_id_counts$Var1, collapse = ","))
  }

# Function used below for the heatmap
reorder_percentage_matrix <- function(percentage_matrix) {
  # Use percentages between variables as distance
  dd <- as.dist((100 - percentage_matrix) / 100)
  hc <- hclust(dd)
  percentage_matrix <- percentage_matrix[hc$order, hc$order]
  return(percentage_matrix)
}

# Function that creates a heatmap of the gene_id that pairs of the clinical parameters share
clinical_parameter_gene_assotiations_heatmap <- function(input_df, output_dir) {

  # Create cross-tabulation of clinical parameters and gene IDs
  # Simply a list matrix, with one axis beein the clinical parameters and the other the gene IDs
  # If there is a representetion of a gene_id for a clinical parameter, the value is 1, otherwise 0
  cross_tab <- table(input_df$clinical_parameter, input_df$gene_id)

  # Calculate the co-occurrences of gene IDs for each pair of clinical parameters
  # by performing matrix multiplication between the cross-tabulation matrix and its transpose
  co_occurrences <- cross_tab %*% t(cross_tab)

  # Calculate percentages
  total_genes <- ncol(cross_tab)
  message("Total genes represented in heatmap dataframe: ", total_genes)
  percentage_matrix <- (co_occurrences / total_genes) * 100

  # Reorder the percentage matrix
  percentage_matrix <- reorder_percentage_matrix(percentage_matrix)

  # Convert matrix to data frame
  percentage_df <- as.data.frame(as.table(percentage_matrix))
  colnames(percentage_df) <- c("Clinical_Parameter_X", "Clinical_Parameter_Y", "Percentage")

  # Reverse the order of the y-axis by reordering the factor levels
  percentage_df$Clinical_Parameter_Y <- factor(percentage_df$Clinical_Parameter_Y, levels = rev(levels(percentage_df$Clinical_Parameter_Y)))

  # Set Percentage to NA for tiles above the diagonal
  percentage_df <- percentage_df %>%
    mutate(Percentage = ifelse(as.numeric(factor(Clinical_Parameter_X, levels = levels(Clinical_Parameter_Y))) < as.numeric(Clinical_Parameter_Y), NA, Percentage))

  # Create the heatmap
  ggplot(percentage_df, aes(x = Clinical_Parameter_X, y = Clinical_Parameter_Y, fill = Percentage)) +
    geom_tile() +
    scale_fill_gradient(low = "white", high = "#0007c4", na.value = "white") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, color = "black"),  # Rotate x-axis labels by 45 degrees
      axis.text.y = element_text(angle = 45, color = "black"),
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.background = element_rect(fill = "white", color = NA),
      legend.position = c(0.9, 0.8)
    ) +
    labs(title = "Co-occurrence of Gene IDs Across Clinical Parameters",
        x = NULL,
        y = NULL,
        fill = "Gene-set \noverlap\npercentage")

  # Save the plot to a file
  ggsave(paste(output_dir, "clinical_parameter_gene_associations_heatmap.png", sep = ""), dpi = 200, height = 8, width = 8) 
  
  }


# Function that exports thefiltered segments to a bed file
export_filtered_segments <- function(input_df, output_dir) {
  # Copy the input_df to a new data frame
  output_df <- input_df
  # Combine the columns gene_id, spearman_corr and p_value to a single column called name and seperate the values with a underscore
  output_df$geneSpearmanPval <- paste(output_df$gene_id, output_df$spearman_corr, output_df$p_value, sep = "_")

  # Remove the columns gene_id, spearman_corr and p_value
  output_df <- output_df[, !(names(output_df) %in% c("gene_id", "spearman_corr", "p_value", "clinical_parameter"))]
  
  # Separate the segment column into chrom, chromStart, and chromEnd
  output_df <- separate(output_df, segment, into = c("chr", "start", "end"), sep = "\\.")
  
  print(head(output_df))

  # Drop duplicates based on chrom, chromStart, and chromEnd
  output_df <- output_df %>% distinct(chr, start, end, .keep_all = TRUE)

  # Save the output_df to a bed file
  # Create the file path
  file_path <- paste(output_dir, "filtered_segments.bed", sep = "/")
  # Write the header with a #
  writeLines(paste("#", paste(colnames(output_df), collapse = "\t"), sep = ""), con = file_path)
  # Append the data to the file
  write.table(output_df, file = file_path, sep = "\t", row.names = FALSE, col.names = FALSE, quote = FALSE, append = TRUE)
}

if (TRUE) {
  # Path to directory with input files
  input_dir = "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/CorrelationVisualizationFiltered/runs/v1/"

  input_df <- load_data(input_dir)

  output_dir = "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/CorrelationVisualizationFiltered/visulizations/v1/"

  # Plot the number of segments for each clinical parameter
  plot_entries(input_df, output_dir)

  # Drop the clinical parameters that have no meaningful correlations
  clinical_params_to_drop <- c("age", "height", "weight", "NYHA", "ZVD", "heart.rate", "Paradoxe_Septumbewegung", "Perikarderguss", "SMW")

  # Drop the rows with the clinical parameters that have no meaningful correlations
  input_df <- input_df[!input_df$clinical_parameter %in% clinical_params_to_drop, ]
  # Reset the row names
  rownames(input_df) <- NULL

  # Plot the number of Gene IDs across all clinical parameters
  plot_top_gene_ids(input_df, output_dir)

  # Create reduced df version without the columns segement, spearman_corr and p_value
  reduced_df <- input_df[, c("clinical_parameter", "gene_id")]

  # Drop duplicated rows
  reduced_df <- unique(reduced_df)

  # Create a heatmap of the gene_id that pairs of the clinical parameters share
  clinical_parameter_gene_assotiations_heatmap(reduced_df, output_dir) 

  # Export the filtered segments to a bed file
  export_filtered_segments(input_df, output_dir)

}
