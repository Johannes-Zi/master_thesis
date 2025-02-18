library(ggplot2)
library(crayon)
library(dplyr)
library(reshape2)

#'
#' Import clinical metadata
import_clinical_meatadata <- function(file_path) {
  # Import clinical metadata
  clinical_metadata <- read.csv(file_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  
  # Rename column names
  colnames(clinical_metadata) <- gsub("age", "age", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("\\bAT\\b", "AT", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("AT_ET", "AT/ET", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("AVDO2", "AVDO2", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("BNP", "BNP", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("Cardiac.Index", "Cardiac Index", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("CVP", "CVP", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("DLCO", "DLCO", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("eGFR", "eGFR", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("heart.rate", "heart rate", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("height", "height", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("HZV_Fick", "HZV Fick", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("HZV_Thermodil", "HZV Thermodil", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("Kreatinin", "kreatinin", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("mPAP", "mPAP", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("NYHA", "NYHA", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("PA_diastolic", "PAdia", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("PA_systolic", "PAsys", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("Paradoxe_Septumbewegung", "PSM", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("PAWP", "PAWP", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("Perikarderguss", "PE", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("PVR", "PVR", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("RA_area", "RAS", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("Rrsys", "Rrsys", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("RVEDD", "RVEDD", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("\\bS\\b", "TAPSV", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("sPAP.excl..ZVD", "sPAP.est", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("TAPSE", "TAPSE", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("Tiffeneau.Index.FEV1.VC", "FEV1/FVC", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("ZVD", "CVP.est", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("VCI_diameter", "VCId", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("ven_SO2", "venSO2", colnames(clinical_metadata))
  colnames(clinical_metadata) <- gsub("\\bFEV1\\b", "FEV1", colnames(clinical_metadata))
  
  return(clinical_metadata)
}


#'
#' Get lower triangle of the correlation matrix
get_lower_tri <- function(cormat) {
  cormat[upper.tri(cormat)] <- NA
  return(cormat)
}


#'
#' Get upper triangle of the correlation matrix
get_upper_tri <- function(cormat) {
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}


#'
#' 
reorder_cormat <- function(cormat) {
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
}

#'
#' Correlate clinical parameters pairwise and visualize the correlation matrix as a heatmap
#' 
clinical_parameters_heatmap <- function(clinical_metadata) {

  # Chop of the first three columns, which are not relevant for the correlation analysis
  clinical_metadata <- clinical_metadata[, -c(1:3)]

  # Replace , with . in the column values
  clinical_metadata_numeric <- clinical_metadata %>%
  mutate(across(everything(), ~ as.numeric(gsub(",", ".", .))))

  # Check for missing values
  print(colSums(is.na(clinical_metadata_numeric)))

  # Replace missing values with the column mean
  clinical_metadata_numeric <- clinical_metadata_numeric %>%
  mutate(across(everything(), ~ ifelse(is.na(.), mean(., na.rm = TRUE), .)))
  #print(head(clinical_metadata_numeric))

  # Each clinical parameter column in the matrix can be seen as vector, we are
  # now computing the pairwise correlation between all vectors
  correlation_matrix <- cor(clinical_metadata_numeric, method = "pearson")
  #print(head(correlation_matrix))

  # Reorder the correlation matrix
  cormat <- reorder_cormat(correlation_matrix) 

  # Get the lower triangle of the correlation matrix
  upper_tri <- get_upper_tri(cormat)

  # Melt the correlation matrix to a long format
  melted_cormat  <- melt(upper_tri, na.rm = TRUE)

  # Create a ggheatmap
  ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value)) +
    geom_tile(color = "#848484") +
    scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white",
                         midpoint = 0, limit = c(-1, 1), space = "Lab",
                         name = "Pearson\nCorrelation    ") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, size = 22, hjust = 1),
      axis.text.y = element_text(angle = 45, vjust = 1, size = 22, hjust = 1, margin = margin(t = 0, r = 10, b = 0, l = 0)),
      panel.grid.major = element_line(color = "#dedede"),
      plot.background = element_rect(fill = "white", color = NA),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      legend.position = "bottom",
      ) +
    coord_fixed()

  # # Add correlation coefficients on the heatmap
  # ggheatmap <- ggheatmap + 
  # geom_text(aes(Var2, Var1, label = round(value, 2)), color = "#464646", 
  #           size = 3, na.rm = TRUE) +
  # theme(
  #     legend.justification = c(1, 0),
  #     legend.position = c(0.3, 0.85),
  #     legend.direction = "horizontal",
  #     legend.title = element_text(size = 14),
  #     legend.text = element_text(size = 13)
  # ) +
  # guides(fill = guide_colorbar(barwidth = 12, barheight = 2,
  #                              title.position = "top", title.hjust = 0.5))


  # Save the plot
  output_dir <- "AssociationAnalysis/clinicalParamsCorr/"
  output_file <- paste(output_dir, "clinical_parameters_correlation_heatmap.png", sep = "")
  ggsave(output_file, ggheatmap, width = 20, height = 20)
  # Save as svg
  output_file_svg <- paste(output_dir, "clinical_parameters_correlation_heatmap.svg", sep = "")
  ggsave(output_file_svg, ggheatmap, width = 20, height = 20)
}

# Import function to import the data
import_data <- function() {
    # Import clinical metadata
    clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
    clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)
    message(paste("Imported clinical metadata from: \n", clinical_metadata_path))
    
    # Create heatmap of clinical parameters
    clinical_parameters_heatmap(clinical_metadata)
}

import_data()
