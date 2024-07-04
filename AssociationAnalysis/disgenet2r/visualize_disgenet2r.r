




# Function to convert ratio character to numeric
convert_ratio_to_numeric <- function(ratio) {
  sapply(strsplit(ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# Load the DisGeNET results df from csv
combined_disgenet_results_df <- read.csv("DisGeNETresults_0.4pqcutoff_uptotop225correlations.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)
# Convert ratio columns to numeric
combined_disgenet_results_df$GeneRatio <- convert_ratio_to_numeric(combined_disgenet_results_df$GeneRatio)
str(combined_disgenet_results_df)

#'
#' Function that retrieves the combined_disgentes_results_df and a specific p.adjst threshold and calculates for each 
#' description the number of occurrences with a p.adjust below the threshold saved in the combined_disgenet_results_df 
#' as "DiseaseCount" and filters out the top N most abundant deseaes
extract_top_frequent_diseases_with_padjust_below_threshold <- function(combined_disgenet_results_df, p_adjust_threshold, max_represented_diseases) {
  # Filter out rows with p.adjust below the threshold
  combined_disgenet_results_df_filtered <- combined_disgenet_results_df[combined_disgenet_results_df$p.adjust < p_adjust_threshold, ]

  # Count the occurrences of each deseases
  description_counts <- table(combined_disgenet_results_df_filtered$Description)

  # Convert to a data frame
  description_counts_df <- as.data.frame(description_counts, stringsAsFactors = FALSE)
  names(description_counts_df) <- c("Description", "DiseaseCount")

  # Sort by count
  description_counts_df_sorted <- description_counts_df[order(-description_counts_df$DiseaseCount), ]

  # Reset the column names to standard
  rownames(description_counts_df_sorted) <- NULL

  # Apply the max_represented_diseases threshold
  description_counts_df_sorted_filtered <- description_counts_df_sorted[1:max_represented_diseases, ]

  # Merge with the original data frame to include the count of each description and throw out the rows with p.adjust above the threshold
  combined_disgenet_results_df_merged <- merge(combined_disgenet_results_df, description_counts_df_sorted_filtered, by = "Description")
  
  # Sort by Count (to get most abundant descriptions at the top)
  combined_disgenet_results_df_merged_sorted <- combined_disgenet_results_df_merged[order(-combined_disgenet_results_df_merged$DiseaseCount), ]

  # Return the sorted data frame
  return(combined_disgenet_results_df_merged_sorted)
}


#'
#' Function calculates a DeseaseCount, which is the number of occurrences of each description in the data frame. and 
#' then filters out the top N most abundant descriptions
extract_top_frequent_diseases_with_high_desease_representation <- function(data_frame, top_n = 20) {
  # Count the occurrences of each description
  description_counts <- table(data_frame$Description)
  
  # Convert to a data frame
  description_counts_df <- as.data.frame(description_counts, stringsAsFactors = FALSE)
  names(description_counts_df) <- c("Description", "DiseaseCount")
  
  # Sort by count
  description_counts_df_sorted <- description_counts_df[order(-description_counts_df$DiseaseCount), ]
  
  # Filter out the top N most abundant descriptions
  description_counts_df_sorted_filtered <- description_counts_df_sorted[1:top_n, ]
  
  # Reset the column names to standard
  rownames(description_counts_df_sorted_filtered) <- NULL
  
  # Merge with the original data frame to include the count of each description
  combined_disgenet_results_df_merged <- merge(data_frame, description_counts_df_sorted_filtered, by = "Description")

  # Sort by Count (to get most abundant descriptions at the top)
  combined_disgenet_results_df_merged_sorted <- combined_disgenet_results_df_merged[order(-combined_disgenet_results_df_merged$DiseaseCount), ]
  
  return(combined_disgenet_results_df_merged_sorted)
}


# Extract the top frequent diseases with p.adjust below the threshold
combined_disgenet_results_df_merged_1 <- extract_top_frequent_diseases_with_padjust_below_threshold(combined_disgenet_results_df, 0.3, 20)
# Filter out rows with p.adjust above 0.5 p.adjust threshold
combined_disgenet_results_df_merged_1_filtered <- combined_disgenet_results_df_merged_1[combined_disgenet_results_df_merged_1$p.adjust < 0.4, ]

# print(head(combined_disgenet_results_df_merged_1))
# print(nrow(combined_disgenet_results_df_merged_1))
# print(tail(combined_disgenet_results_df_merged_1))

# Extract the top 20 most abundant diseases
combined_disgenet_results_df_merged_2 <- extract_top_frequent_diseases_with_high_desease_representation(combined_disgenet_results_df, 20)
# Filter out rows with p.adjust above p.adjust threshold
combined_disgenet_results_df_merged_2_filtered <- combined_disgenet_results_df_merged_2[combined_disgenet_results_df_merged_2$p.adjust < 0.4, ]

# print(head(combined_disgenet_results_df_merged_2))
# print(tail(combined_disgenet_results_df_merged_2))
# print(nrow(combined_disgenet_results_df_merged_2))

library(ggplot2)

# Create the plot
ggplot(combined_disgenet_results_df_merged_1_filtered, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_size_continuous(name = "GeneRatio", range = c(2, 10)) +
  scale_color_gradient(name = "p.adjust", low = "#1d1db9", high = "#ce2323") +
  theme_minimal() +
  labs(title = "Filtered DisGeNET results",
       x = "Cluster", y = "Disease") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

