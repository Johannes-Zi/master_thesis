# Load required libraries
library(MASS)
library(ggplot2)
library(viridis)
library(glmnet)
library(ggExtra)


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
#input_directory_path <- "C:/Users/johan/Desktop/regession_output_example/"
input_directory_path <- "C:/Users/johan/Desktop/regression_output/regression_output/"


cat("imported directory path:\n", input_directory_path, "\n")

# Process the directory and get the results
output_df <- process_directory(input_directory_path, 20000, "\\Pearson.RData$")
output_df <- process_directory(input_directory_path, 20000, "\\Spearman.RData$")

# Print the resulting data frame
#print(output_df)

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

# Create the dotplot
ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
  geom_point() +
  labs(x = "Number of Features", y = "nzero") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, max(output_df$number_of_features)), ylim = c(0, max(output_df$nzero_value))) +
  geom_abline(intercept = 0, slope = 1, color = "red")

# Create the countplot
ggplot(output_df, aes(x = number_of_features)) +
  geom_bar() +
  labs(x = "Number of Features", y = "Count") +
  theme_minimal()

# # Create the dotplot
# ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
#   geom_point() +
#   labs(x = "Number of Features", y = "nzero") +
#   theme_minimal() +
#   coord_cartesian(xlim = c(0, max(output_df$number_of_features)), ylim = c(0, max(output_df$nzero_value))) +
#   geom_abline(intercept = 0, slope = 1, color = "red")

# # Create the marginal boxplot
# ggMarginal(x, type = "boxplot", fill = "transparent")

# # Display the plot
# print(y)


# Create the dotplot with counts
ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
  geom_point() +
  geom_count() +
  labs(x = "Number of Features", y = "nzero") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, max(output_df$number_of_features)), ylim = c(0, max(output_df$nzero_value))) +
  geom_abline(intercept = 0, slope = 1, color = "red")

  # Create the dotplot with counts
ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
  geom_point() +
  geom_count() +
  labs(x = "Number of Features", y = "nzero") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, max(output_df$number_of_features)), ylim = c(0, max(output_df$nzero_value))) +
  geom_abline(intercept = 0, slope = 1, color = "red")



  # Create the dotplot with counts
ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
  geom_point() +
  geom_count() +
  labs(x = "Number of Features", y = "nzero") +
  theme_minimal() +
  coord_cartesian(xlim = c(0, max(output_df$number_of_features)), ylim = c(0, max(output_df$nzero_value))) +
  geom_abline(intercept = 0, slope = 1, color = "red")


# # Define the function to get density of points in 2 dimensions
# get_density <- function(x, y, ...) {
#   dens <- MASS::kde2d(x, y, ...)
#   ix <- findInterval(x, dens$x)
#   iy <- findInterval(y, dens$y)
#   ii <- cbind(ix, iy)
#   return(dens$z[ii])
# }

# # Add a new column to your data frame for density
# output_df$density <- get_density(output_df$number_of_features, output_df$nzero_value, n = 100)

# # Create the plot
# ggplot(output_df, aes(x = number_of_features, y = nzero_value, color = density)) +
#   geom_point() +
#   labs(x = "Number of Features", y = "nzero") +
#   theme_minimal() +
#   coord_cartesian(xlim = c(0, max(output_df$number_of_features)), ylim = c(0, max(output_df$nzero_value))) +
#   geom_abline(intercept = 0, slope = 1, color = "red") +
#   scale_color_viridis()

# Create the plot
p <- ggplot(output_df, aes(x = number_of_features, y = nzero_value)) +
  geom_count(aes(color = ..n..)) +
  scale_color_viridis_c() +
  labs(x = "Number of Features", y = "nzero") +
  ggtitle(expression(paste(bold("Regression Spearman - "), "Number of non-zero Coefficients vs. Number of Features"))) +
  theme_bw() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_continuous(breaks = seq(0, max(output_df$number_of_features), by = 10)) +
  scale_y_continuous(breaks = seq(0, max(output_df$nzero_value), by = 10)) +
  coord_cartesian(xlim = c(0, max(output_df$number_of_features)), ylim = c(0, max(output_df$nzero_value))) +
  geom_abline(intercept = 0, slope = 1, color = "red")

  # Save the plot as a PNG file
ggsave("plot.png", plot = p, width = 12, height = 12, dpi = 600)
