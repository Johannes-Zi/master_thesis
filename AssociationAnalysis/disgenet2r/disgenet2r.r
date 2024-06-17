# if (!require("BiocManager", quietly = TRUE))
#     install.packages("BiocManager")

# BiocManager::install("clusterProfiler")
#BiocManager::install("org.Hs.eg.db")
#BiocManager::install("gprofiler2")
#BiocManager::install("enrichplot")

library(enrichplot)
library(gprofiler2)
library(clusterProfiler)
library(org.Hs.eg.db)
library(DOSE)
library(crayon)
library(ggplot2)

#install.packages("fastmap")


#'
#' Function that imports the correlattion results of BasicRegression.r and runs for the best Results a DisGeNET analysis for each clinical parameter
#' 

#'
#' Import clinical metadata
import_clinical_meatadata <- function(file_path) {
  # Import clinical metadata
  clinical_metadata <- read.csv(file_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)
  return(clinical_metadata)
}


#'
#' Import Regresion_results
#' 
import_regression_results <- function(toplevel_path, clinical_parameter_names){

    # Initialize list to store regression results
    regression_results_list <- list()

    # Iterate over all clinical parameters
    for (clinical_parameter in clinical_parameter_names) {
        # Create input path
        regression_results_path <- paste0(toplevel_path, clinical_parameter, "/", clinical_parameter, "_spearman_correlations.tsv")

        # Import regression results
        regression_results <- read.csv(regression_results_path, header = TRUE, sep = ";", stringsAsFactors = FALSE)

        # Append regression results to list
        regression_results_list[[clinical_parameter]] <- regression_results
    }


    # Return list with all regression results of each clinical parameter
    return(regression_results_list)
}


#'
#' Creates List with reduced datasets with the top regressions of each clinical parameter
#' 
create_reduced_regression_datasets <- function(regression_results_list){
    # Initialize list to store reduced regression datasets
    reduced_regression_dfs_list <- list()

    # Iterate over all clinical parameters
    for (clinical_parameter in names(regression_results_list)) {
        # Extract current regression results
        regression_results <- regression_results_list[[clinical_parameter]]

        # Transform the regression results into a df
        regression_results <- as.data.frame(regression_results)
        
        # Add colum with the absolute value of the correlation
        regression_results$absolute_spearman_corr <- abs(regression_results$spearman_corr)
        
        # Extract rows with absolute regression scores above 0.5
        reduced_regression_results <- regression_results[regression_results$absolute_spearman_corr > 0.5,]

        # Sort the reduced regression results by the absolute correlation
        reduced_regression_results <- reduced_regression_results[order(reduced_regression_results$absolute_spearman_corr, decreasing = TRUE),]

        # Print the number of reduced regression results
        #print(paste("Number of reduced regression results for ", clinical_parameter, ": ", nrow(reduced_regression_results)))

        # Create a reduced regression df and use only up to the top 225 absolute correlations
        reduced_regression_results <- reduced_regression_results[1:225,]

        # Filter out duplication in the gene column
        reduced_regression_results <- reduced_regression_results[!duplicated(reduced_regression_results$gene),]

        #print(paste("Number of reduced regression results after filtering out duplications for ", clinical_parameter, ": ", nrow(reduced_regression_results)))

        # Translate the Ensembl ids to NCBI gene names


        # Append reduced regression results df to list
        reduced_regression_dfs_list[[clinical_parameter]] <- reduced_regression_results
    }

    # Return list with reduced regression datasets
    return(reduced_regression_dfs_list)
}


#'
#'  Function that runs DisGeNET analysis for each clinical parameter
#'  
run_digenet_analysis <- function(reduced_regression_datasets_list){
    # Initialize list to store DisGeNET results
    disgenet_results_list <- list()

    # Initialize empty df to  store combined DisGeNET results
    combined_disgenet_results_df <- data.frame()

    # Extract clinical parameters
    clinical_parameters <- names(reduced_regression_datasets_list)

    # Iterate over all clinical parameters
    for (i in 1:length(reduced_regression_datasets_list)) {
        # Current clinical parameter name
        clinical_parameter <- clinical_parameters[i]

        message("\n##################")
        message(paste("current parameter:", cyan(clinical_parameter)))
        message("##################\n")

        # gene_names <- reduced_regression_datasets_list[[clinical_parameter]]$gene
        # gene_names_str <- paste(gene_names, collapse = ", ")
        # message(cyan("\nGene list:"), gene_names_str)

        # Current reduced regression results
        reduced_regression_results <- reduced_regression_datasets_list[[clinical_parameter]]

        # Extract gene names
        gene_names <- reduced_regression_results$gene
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
        current_results <- enrichDGN(gene = ncbi_gene_names, qvalueCutoff  = 0.4, pvalueCutoff  = 0.4)#,
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

# Import clinical metadata and extract column names
clinical_metadata_path <- "C:/Users/johan/OneDrive/dateien_cloud/Master/Semester_4/Masterarbeit/data/pulmanory_hypertension/clinical_data/version_26.05.24/isolated_metadata.csv"
clinical_metadata <- import_clinical_meatadata(clinical_metadata_path)

# Extract column names
column_names <- colnames(clinical_metadata[4:ncol(clinical_metadata)])

# Import regression results based on current cwd
regression_data_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/BasicRegression/combined_correlations/"
regression_results_list <- import_regression_results(regression_data_dir, column_names)

# Create reduced regression datasets
reduced_regression_datasets_list <- create_reduced_regression_datasets(regression_results_list)

# Run DisGeNET analysis
DisGeNETresults <- run_digenet_analysis(reduced_regression_datasets_list)

# Save the DisGeNET results df as csv
write.csv(DisGeNETresults$combined_disgenet_results_df, "DisGeNETresults_0.4pqcutoff_uptotop225correlations.csv", row.names = FALSE)










# Function to convert ratio character to numeric
convert_ratio_to_numeric <- function(ratio) {
  sapply(strsplit(ratio, "/"), function(x) as.numeric(x[1]) / as.numeric(x[2]))
}

# Load the DisGeNET results df from csv
combined_disgenet_results_df <- read.csv("DisGeNETresults_0.4pqcutoff_uptotop225correlations.csv", header = TRUE, sep = ",", stringsAsFactors = FALSE)

#disgenet_results_list <- DisGeNETresults$disgenet_results_list  # Complete results
#combined_disgenet_results_df <- DisGeNETresults$combined_disgenet_results_df    # Combined DisGeNET output dfs, prepared for visualization


print(head(combined_disgenet_results_df))
print(tail(combined_disgenet_results_df))
print(str(combined_disgenet_results_df))


# Convert ratio columns to numeric
combined_disgenet_results_df$GeneRatio <- convert_ratio_to_numeric(combined_disgenet_results_df$GeneRatio)
str(combined_disgenet_results_df)

# Count the occurrences of each description
description_counts <- table(combined_disgenet_results_df$Description)
print(str(description_counts))
# Convert to a data frame
description_counts_df <- as.data.frame(description_counts, stringsAsFactors = FALSE)
str(description_counts_df)
names(description_counts_df) <- c("Description", "DiseaseCount")
print(head(description_counts_df))
cat(str(description_counts_df))
# Sort by count
description_counts_df_sorted <- description_counts_df[order(-description_counts_df$DiseaseCount), ]

# Filter out the top N most abundant descriptions
description_counts_df_sorted_filtered <- description_counts_df_sorted[1:20, ]

# Reset the column names to standard
rownames(description_counts_df_sorted_filtered) <- NULL
print(description_counts_df_sorted_filtered)

# Merge with the original data frame to include the count of each description
combined_disgenet_results_df_merged <- merge(combined_disgenet_results_df, description_counts_df_sorted_filtered, by = "Description")
print(tail(combined_disgenet_results_df_merged))
print(str(combined_disgenet_results_df_merged))

# Sort by Count (to get most abundant descriptions at the top)
combined_disgenet_results_df_merged_sorted <- combined_disgenet_results_df_merged[order(-combined_disgenet_results_df_merged$DiseaseCount), ]
print(head(combined_disgenet_results_df_merged_sorted))
print(tail(combined_disgenet_results_df_merged_sorted))
print(str(combined_disgenet_results_df_merged_sorted))

# # Columns used in the ggplot
# columns_of_interest <- c("Cluster", "Description", "GeneRatio", "p.adjust")

# # Find rows with NA in any of the columns of interest
# rows_with_na <- apply(combined_disgenet_results_df_merged_sorted_filtered[columns_of_interest], 1, function(row) any(is.na(row)))


# print(rows_with_na)
# # Print rows with NA values
# na_entries <- combined_disgenet_results_df_merged_sorted_filtered[rows_with_na, ]
# print(na_entries)


library(ggplot2)

# Create the plot
ggplot(combined_disgenet_results_df_merged_sorted, aes(x = Cluster, y = Description)) +
  geom_point(aes(size = GeneRatio, color = p.adjust)) +
  scale_size_continuous(name = "GeneRatio", range = c(2, 10)) +
  scale_color_gradient(name = "p.adjust", low = "#1d1db9", high = "#ce2323") +
  theme_minimal() +
  labs(title = "DisGeNET results for top 20 most abundant deseases",
       x = "Cluster", y = "Deseas") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


























# # Access first DisGeNET result
# x1 <- disgenet_results_list[[2]]
# print(x1)



# barplot(x1, showCategory=20) 
# dotplot(x1, showCategory=30) + ggtitle("dotplot for ORA")
# cnetplot(x1, foldChange=geneList)
# cnetplot(x1, categorySize="pvalue", foldChange=geneList)
# cnetplot(x1, foldChange=geneList, circular = TRUE, colorEdge = TRUE)



# data(geneList)
# de <- names(geneList)[abs(geneList) > 2]
# edo <- enrichDGN(de)
# barplot(edo, showCategory=20)



# data(gcSample)
# str(gcSample) 


# ck <- compareCluster(geneCluster = gcSample, fun = enrichKEGG)
# str(ck)
# ck1 <- setReadable(ck, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
# str(ck1)
# head(ck1)

# dotplot(ck1)
