library(dplyr)

# Import Performance_Overview files
input_directory_path_LeaveOneOut <- "C:/Users/johan/Desktop/LOneOCV_regression/performance_evaluation/Performance_Overview.txt"
#input_directory_path_standard <- "C:/Users/johan/Desktop/standard_regression/performance_evaluation/Performance_Overview.txt"

# Create dataframes
df_LeaveOneOut <- read.table(input_directory_path_LeaveOneOut, header = TRUE, sep = "\t")
#df_standard <- read.table(input_directory_path_standard, header = TRUE, sep = "\t")
head(df_LeaveOneOut)
#head(df_standard)

process_regression_performance_evaluation_df <- function(input_df) {
    # Adapt the filename column
    filename_parts <- strsplit(as.character(input_df$Sample_Name), "_")
    head(filename_parts)
    segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(substr(x[4], 1, nchar(x[4])-4), sep = "_"))
    head(segmentation_preselection_corellation)
    adapted_df <- input_df
    adapted_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

    extracted_gene_names <- sapply(filename_parts, function(x) paste(x[2], x[3], sep = "_"))
    adapted_df$Sample_Name <- extracted_gene_names

    adapted_df <- rename(adapted_df, gene_name = Sample_Name)

    return(adapted_df)
}

processed_performance_evaluation_df <- process_regression_performance_evaluation_df(df_LeaveOneOut)
head(processed_performance_evaluation_df)

# Import elastic net models
# Function to process each file
process_file <- function(file_path) {
  # Load the .RData file - accessible as 'elasticnet_model'
  load(file_path)

  # Get the filename without extension
  filename <- tools::file_path_sans_ext(basename(file_path))

  # Get ident of lambda min based on the list with all lambdas
  lambda_min_value <- elasticnet_model$model$lambda.min
  lambda_min_index <- which(elasticnet_model$model$lambda == lambda_min_value)

  # Get the corresponding nzero value
  nzero_value <- elasticnet_model$model$nzero[lambda_min_index]
  number_of_features <- elasticnet_model$model$glmnet.fit$dim[1]

  # Create a data frame with the results
  result <- data.frame(filename = filename,
                       number_of_features = number_of_features,
                       nzero_value = nzero_value)
  return(result)
}

process_directory <- function(directory_path, import_limit, file_pattern) {
  # Get the list of .RData files in the directory
  file_list <- list.files(directory_path, pattern = file_pattern,
                          full.names = TRUE)
  print(paste("Number of detected files:", length(file_list)))
  print(paste("Used filename pattern", file_pattern))

  # Limit the number of files to import
  if (length(file_list) > import_limit) {
    file_list <- file_list[1:import_limit]
    print(paste("Number of imported files limited to", import_limit))
  }

  # Initialize an empty data frame to store the results
  result_df <- data.frame(filename = character(),
                          number_of_features = numeric(),
                          nzero_value = numeric(),
                          stringsAsFactors = FALSE)

  # Process each file and append the results to the data frame
  for (i in seq_along(file_list)) {
    file_path <- file_list[i]
    result <- process_file(file_path)
    result_df <- rbind(result_df, result)
    
    # Print the current number of processed files every 50 files
    if (i %% 50 == 0) {
      print(paste("Processed", i, "files"))
    }
  }

  return(result_df)
}

# Define the input directory path
#input_directory_path <- "C:/Users/johan/Desktop/standard_regression/regression_output/regression_output/"
input_directory_path <- "C:/Users/johan/Desktop/LOneOCV_regression/regression_output/"
cat("imported directory path:\n", input_directory_path, "\n")

# Process the directory and get the results
output_df <- process_directory(input_directory_path, 20000, "\\Pearson.RData$")
head(output_df)
#output_df <- process_directory(input_directory_path, 20000, "\\Spearman.RData$")


process_regression_models_df <- function(input_df) {
    # Create a new dataframe using the filename as key - extracting only the filename
    new_df <- data.frame(filename = unique(input_df$filename))

    # Merge the new dataframe with the input_df dataframe
    merged_df <- merge(new_df, input_df, by = "filename", all.x = TRUE)

    # Adapt the filename column
    filename_parts <- strsplit(as.character(merged_df$filename), "_")
    segmentation_preselection_corellation <- sapply(filename_parts, function(x) paste(x[7], sep = "_"))
    merged_df$segmentation_preselection_corellation <- segmentation_preselection_corellation

    extracted_gene_names <- sapply(filename_parts, function(x) paste(x[5], x[6], sep = "_"))
    merged_df$filename <- extracted_gene_names

    merged_df <- rename(merged_df, gene_name = filename)

    return(merged_df)
}

# Call the function with the output_df dataframe
processed_regression_models_df <- process_regression_models_df(output_df)
head(processed_regression_models_df)




# Combine the two dataframes
combined_df <- merge(processed_performance_evaluation_df, processed_regression_models_df, by = c("gene_name", "segmentation_preselection_corellation"), all = TRUE)
head(combined_df)

combined_filtered_df <- combined_df[combined_df$segmentation_preselection_corellation == "Pearson", ]
head(combined_filtered_df)

# Filter the data frame to only show rows with missing values in the nzero_value or Pearson columns
missing_values_df <- combined_filtered_df[is.na(combined_filtered_df$nzero_value) | is.na(combined_filtered_df$Pearson), ]

# Remove rows with missing values in the nzero_value or Pearson columns
combined_filtered2_df <- combined_filtered_df[!is.na(combined_filtered_df$nzero_value) & !is.na(combined_filtered_df$Pearson), ]
head(combined_filtered2_df)

# Print the rows with missing values
print(missing_values_df)

library(ggplot2)

# Create the dotplot
dotplot <- ggplot(combined_filtered2_df, aes(x = nzero_value, y = Pearson)) +
    geom_point() +
    xlab("Number of Nonzero Values") +
    ylab("Pearson Value")

# Display the dotplot
dotplot
ggsave("dotplot.png", dotplot)

# Load the required packages
library(MASS)
library(ggplot2)
library(viridis)

# Define the function to get density
get_density <- function(x, y, ...) {
  dens <- MASS::kde2d(x, y, ...)
  ix <- findInterval(x, dens$x)
  iy <- findInterval(y, dens$y)
  ii <- cbind(ix, iy)
  return(dens$z[ii])
}

# Calculate the density for each point in your data
combined_filtered2_df$density <- get_density(combined_filtered2_df$nzero_value, combined_filtered2_df$Pearson, n = 100)

# Create the dotplot
dotplot <- ggplot(combined_filtered2_df, aes(x = nzero_value, y = Pearson, color = density)) +
    geom_point() +
    xlab("Number of Nonzero Feature Coefficients") +
    ylab("Pearson Values") +
    scale_color_viridis() +
    scale_x_continuous(breaks = seq(min(combined_filtered2_df$nzero_value), max(combined_filtered2_df$nzero_value), by = 1))  # Increase the number of breaks on the x-axis

# Save the dotplot
ggsave("dotplot2.png", dotplot, width = 20,)
