library(ggplot2)
library(RColorBrewer)
library(rstudioapi)
library(dplyr)

# Import Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
head(df_standard)

preprocess_performance_evaluation_df <- function(input_df) {
  # This function preprocesses the performance evaluation dataframe by adapting the filename column and extracting gene names.
  #
  # Args:
  #   input_df: The input dataframe containing the performance evaluation data.
  #
  # Returns:
  #   The preprocessed dataframe with adapted filename column and extracted gene names.

  # Adapt the filename column
  filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")  # Split the Sample_Name column by "_"
  head(filename_parts)  # Print the first few filename parts

  # Extract the segmentation preselection correlation
  segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
  head(segmentation_preselection_corellation)  # Print the first few segmentation preselection correlations

  # Create a new column in the dataframe for the segmentation preselection correlation
  adapted_df <- input_df
  adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

  # Extract the gene names
  extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
  adapted_df$Sample_Name <- extracted_gene_names

  # Rename the Sample_Name column to gene_name
  adapted_df <- rename(adapted_df, gene_name = Sample_Name)

  return(adapted_df)
}


create_violin_plots <- function(df_LeaveOneOut, df_standard, target_column, output_directory) {
  # This function creates violin plots comparing the target column between two dataframes.
  #
  # Args:
  #   df_LeaveOneOut: The dataframe containing the LeaveOneOut performance evaluation data.
  #   df_standard: The dataframe containing the standard performance evaluation data.
  #   target_column: The column to compare between the dataframes.
  #   output_directory: The directory to save the generated plot.
  #
  # Returns:
  #   A list containing the generated violin plot and the combined dataframe.

  # Extract x axis names
  df_LeaveOneOut_name <- tail(strsplit(deparse(substitute(df_LeaveOneOut)), "_")[[1]], 1)
  df_standard_name <- tail(strsplit(deparse(substitute(df_standard)), "_")[[1]], 1)

  # Create new data frames that only contain the target column and the new CV_Type column
  df_LeaveOneOut_new <- data.frame(CV_Type = df_LeaveOneOut_name, target = df_LeaveOneOut[[target_column]])
  df_standard_new <- data.frame(CV_Type = df_standard_name, target = df_standard[[target_column]])
  print(head(df_LeaveOneOut_new))
  print(head(df_standard_new))

  # Combine the new data frames into one
  df_combined <- rbind(df_LeaveOneOut_new, df_standard_new)

  # Create violin plots
  violinplot <- ggplot(df_combined, aes(x = CV_Type, y = target, fill = CV_Type)) +
    geom_violin(trim = FALSE, width = 0.5) +
    labs(title = paste(target_column, "Comparison", sep = " "), y = target_column, x = "CV Type") +
    theme_gray() +
    theme(plot.title = element_text(hjust = 0.5, size = 16, face = "bold"),
          axis.title = element_text(size = 14, face = "bold"),
          axis.text = element_text(size = 12, face = "bold")) +
    scale_fill_brewer(palette = "Set2") + # Apply Brewer color palette
    theme(legend.position = "none") # Remove legend

  # Save the plot to the target directory
  ggsave(filename = paste(output_directory, "/", target_column, "_violinplot.png", sep = ""), plot = violinplot)
  return(violinplot)
}

output_path <- dirname(rstudioapi::getActiveDocumentContext()$path)
print(output_path)

# Preprocess the Performance_Evaluation dataframes
df_LeaveOneOut_preprocessed_1 <- process_regression_performance_evaluation_df(df_LeaveOneOut)
head(df_LeaveOneOut_preprocessed_1)
df_LeaveOneOut_preprocessed_2 <- df_LeaveOneOut_preprocessed_1[df_LeaveOneOut_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_LeaveOneOut_preprocessed_2)

df_standard_preprocessed_1 <- process_regression_performance_evaluation_df(df_standard)
head(df_standard_preprocessed_1)
df_standard_preprocessed_2 <- df_standard_preprocessed_1[df_standard_preprocessed_1$segmentation_preselection_corellation == "Pearson", ]
head(df_standard_preprocessed_2)

# Call the function with the dataframes as arguments
create_violin_plots(df_LeaveOneOut, df_standard, target_column = "Pearson", output_directory = output_path)
