library(ggplot2)

# Function that loads the data from the input directories
load_data <- function(input_dir) {

  # Dataframe that hold all the input data
  combined_df <- data.frame()

  # Iterate over directories in input_dir
  for (dir in list.dirs(input_dir, full.names = FALSE, recursive = FALSE)) {
    # Set directory name as the name of the current clinical parameter
    clinical_parameter <- dir

    # Build the path to the input file
    input_file <- paste(input_dir, dir, '/', dir, '_spearman_correlations.tsv', sep = '')

    # Load the input tsv file to a data frame
    input_df <- read.table(input_file, header = TRUE, sep = ';')

    # Add the clinical parameter to the data frame
    input_df$clinical_parameter <- clinical_parameter

    # Add the data frame to the combined data frame
    combined_df <- rbind(combined_df, input_df)
  }
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
    labs(title = 'Number of top segments for each clinical parameter',
         x = 'Clinical parameter',
         y = 'Number of segments')

  # Save the plot to a file 
  ggsave(paste(output_dir, 'number_of_entries_per_param.png', sep = ''))
}

# Function that creates a barchart with the most abundand gene_id entries across all clinical parameters
plot_top_gene_ids <- function(input_df, output_dir) {
  # Create a data frame with the number of entries for each gene_id
  gene_id_counts <- as.data.frame(table(input_df$gene_id))

  # Sort the data frame by the number of entries
  gene_id_counts <- gene_id_counts[order(-gene_id_counts$Freq), ]

  # Reset the row names
  rownames(gene_id_counts) <- NULL

  print(gene_id_counts)
  
  # Ensure the Var1 column is a factor with levels ordered by Freq
  gene_id_counts$Var1 <- factor(gene_id_counts$Var1, levels = gene_id_counts$Var1[order(gene_id_counts$Freq, decreasing = TRUE)])

  # Plot the Freqency of each gene_id as a barchart
  ggplot(gene_id_counts[1:40, ], aes(y = Var1, x = Freq)) +
    geom_bar(stat = 'identity') +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.background = element_rect(fill = "white", color = NA)  # Add white background to the plot
    ) +
    labs(
        title = 'Top 10 most abundant gene_ids',
        y = 'Gene ID',
        x = 'Number of entries'
    )

  # Save the plot to a file
  ggsave(paste(output_dir, 'top_gene_ids.png', sep = ''))


  # Calculate the distribution of frequencies
  frequency_distribution <- as.data.frame(table(gene_id_counts$Freq))
  colnames(frequency_distribution) <- c("Frequency", "Number_of_Genes")

  # Convert Frequency to numeric for proper plotting
  frequency_distribution$Frequency <- as.numeric(as.character(frequency_distribution$Frequency))

  # Plot the distribution of frequencies
  ggplot(frequency_distribution, aes(x = Frequency, y = Number_of_Genes)) +
    geom_bar(stat = 'identity') +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5),  # Center the title
      plot.background = element_rect(fill = "white", color = NA),  # Add white background to the plot
    ) +
    labs(
      title = 'Distribution of Gene occurrences',
      x = 'number of clinical parameters a gene is represented as top correlation',
      y = 'Number of Genes'
    )

  # Save the plot to a file
  ggsave(paste(output_dir, 'gene_frequency_distribution.png', sep = ''))

  # Save the gene_id_counts to a file
  write.table(gene_id_counts, paste(output_dir, 'gene_id_counts.tsv', sep = ''), sep = '\t', row.names = FALSE)
  
  # Prit the gene_id_counts to the console comma separated as a single line
  print(paste(gene_id_counts$Var1, collapse = ','))
  }

if (TRUE) {
  # Path to directory with input files
  input_dir = 'C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/CorrelationVisualizationFiltered/runs/v1/'

  input_df <- load_data(input_dir)

  output_dir = 'C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/CorrelationVisualizationFiltered/visulizations/v1/'
  
  # Plot the number of segments for each clinical parameter
  plot_entries(input_df, output_dir)

  # Plot the number of Gene IDs across all clinical parameters
  plot_top_gene_ids(input_df, output_dir)

}
