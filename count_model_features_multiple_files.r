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

# Function to process all files in a directory
process_directory <- function(directory_path) {
  # Get the list of .RData files in the directory
  file_list <- list.files(directory_path, pattern = "\\.RData$", full.names = TRUE)
   # Initialize an empty data frame to store the results
   result_df <- data.frame(filename = character(),
                           number_of_features = numeric(),
                           nzero_value = numeric(),
                           stringsAsFactors = FALSE)
  
   # Process each file and append the results to the data frame
   for (file_path in file_list) {
     result <- process_file(file_path)
     result_df <- rbind(result_df, result)
   }
  
  return(result_df)
}

# Define the input directory path
input_directory_path <- "C:/Users/johan/Desktop/regession_output_example/"

cat("imported directory path:\n", input_directory_path, "\n")

# Process the directory and get the results
output_df <- process_directory(input_directory_path)

# Print the resulting data frame
#print(output_df)

library(ggplot2)

# Create the density dotplot
ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
    geom_density_2d() +
    geom_point() +
    labs(x = "Number of Features", y = "nzero") +
    theme_minimal()

# Create the dotplot
ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
    geom_point() +
    labs(x = "Number of Features", y = "nzero") +
    theme_minimal()


