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
#'  Run GoEnrichment analysis for each clinical parameter
run_enrichment_analysis <- function(target_genes, output_dir, qvalue_cutoff, pvalue_cutoff, 
                                 query_strings){
  # Initialize list to store GoEnrich results
  GoEnrich_results_list <- list()

  # Initialize empty df to  store combined GoEnrich results
  combined_GoEnrich_results_df <- data.frame()

    # Extract gene names and translate them to varios other gene id formats
    if (TRUE) {
    # Extract gene names
    gene_names <- target_genes$gene
    message(paste("Number of genes represented in the correlation results: ", nrow(target_genes)))

    # Traslate the Ensemble ids to NCBI gene names
    translated_ENTREZ_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ENTREZID", OrgDb = "org.Hs.eg.db")))
    translated_ALIAS_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "ALIAS", OrgDb = "org.Hs.eg.db")))
    translated_GENENAME_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENENAME", OrgDb = "org.Hs.eg.db")))
    translated_SYMBOL_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")))
    translated_GENETYPE_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "GENETYPE", OrgDb = "org.Hs.eg.db")))
    translated_UNIPROT_ids <- suppressMessages(suppressWarnings(bitr(gene_names, fromType = "ENSEMBL", toType = "UNIPROT", OrgDb = "org.Hs.eg.db")))

    # Calculate failed mapping rate
    entrez_gene_names <- translated_ENTREZ_ids$ENTREZID
    failed_mappings <- length(gene_names) - length(entrez_gene_names)
    failed_mapping_percentage <- (failed_mappings / length(gene_names)) * 100
    message(paste("Failed to map ", failed_mappings, "gene names to NCBI gene names (", failed_mapping_percentage, "% )"))
    message(paste("Number of genes used for GoEnrich analysis (translated from ENSEMBL to ENTRZIDs): ", length(entrez_gene_names)))

    # Group the ALIAS ids by ENSEMBL ids - collapse multiple Alias ids to one string entry
    translated_ALIAS_ids_grouped <- translated_ALIAS_ids %>%
        group_by(ENSEMBL) %>%
        summarise(ALIAS = paste(ALIAS, collapse = "/")) %>%
        ungroup()
    }

    # Gather metadata of GoEnrich input and save as csv
    if (TRUE) {
      # Merge the metadata with the translated gene ids to a df based on the gene names
      # Create a data frame with all gene_names
      input_genes_metadata_df <- data.frame(ENSEMBL = gene_names)
      # Match and fill in the ids, use NA for unmatched
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_ENTREZ_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_ALIAS_ids_grouped, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_GENENAME_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_SYMBOL_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_GENETYPE_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)
      input_genes_metadata_df <- merge(input_genes_metadata_df, translated_UNIPROT_ids, by.x = "ENSEMBL", by.y = "ENSEMBL", all.x = TRUE)

      # Merge in the information of the input reduced correlation results
      input_genes_metadata_df <- merge(input_genes_metadata_df, target_genes, by.x = "ENSEMBL", by.y = "gene_id", all.x = TRUE)

      # Reorder the DataFrame columns
      cols <- colnames(input_genes_metadata_df)

      # Remove 'SYMBOL' to avoid duplication
      cols <- cols[cols != 'SYMBOL']

      # Insert 'SYMBOL' at the third position
      cols <- c(cols[1:2], 'SYMBOL', cols[3:length(cols)])

      # Reorder the DataFrame columns
      input_genes_metadata_df <- input_genes_metadata_df[, cols]

      # Save the metadata df as a csv file using file.path
      gene_names_output_path <- file.path(output_dir, "GoEnrich_input_genes_metadata.csv")
      write.csv(input_genes_metadata_df, gene_names_output_path, row.names = FALSE)

      # Print the metadata df column SYMBOL as a string with comma seperated values
      gene_names_print <- input_genes_metadata_df$UNIPROT
      gene_names_print <- paste(gene_names_print, collapse = ", ")
      message(paste("Gene names used for GoEnrich analysis: ", gene_names_print))

      # Export the gene names as a csv file
      gene_names_output_path <- file.path(output_dir, "GoEnrich_input_genes.csv")
      write.table(input_genes_metadata_df$UNIPROT, gene_names_output_path, row.names = FALSE)
      write.table(input_genes_metadata_df$UNIPROT, gene_names_output_path, row.names = FALSE, col.names = FALSE, quote = FALSE)

    }

    # Run GoEnrich analysis
    if (TRUE) {

      # Run GoEnrich analysis
      current_results <- enrichGO(gene = entrez_gene_names, OrgDb = org.Hs.eg.db, pvalueCutoff = pvalue_cutoff, qvalueCutoff = qvalue_cutoff)

      message("GoEnrich result insight:")
      print(current_results)
      flush.console()

      # Append GoEnrich results to list
      GoEnrich_results_list[["clinical_parameter"]] <- current_results

      # Extract result df
      current_result_df <- current_results@result
    }

    # Export unfiltered GoEnrich results as csv file
    if (TRUE) {
      # Change the rowname ID to GOid
      current_result_df <- current_result_df %>% 
        rename(
          GO_id = ID,
          GO_descreption = Description,
          GO_GeneRatio = GeneRatio,
          GO_BgRatio = BgRatio,
          GO_pvalue = pvalue,
          GO_p.adjust = p.adjust,
          GO_qvalue = qvalue,
          ENTREZIDs = geneID,
          GO_count = Count
        )

      # Save GoEnrich results as a csv file using file.path
      GoEnrich_output_path <- file.path(output_dir, "GO_enrichment_results.csv")
      write.csv(current_result_df, GoEnrich_output_path, row.names = FALSE)	
    }

    # Append GoEnrich results to combined results df
    if (TRUE) {
      # Remove rownames and use standard rownames
      rownames(current_result_df) <- NULL

      # Add cluster column and set to clinical parameter
      current_result_df$Cluster <- "clinical_parameter"

      # Move cluster column to first position
      current_result_df <- current_result_df[,c(ncol(current_result_df), 1:(ncol(current_result_df)-1))]

      # Append GoEnrich result df to combined df
      combined_GoEnrich_results_df <- rbind(combined_GoEnrich_results_df, current_result_df)
    }
  

  # Return list with GoEnrich results
  return(list(GoEnrich_results_list = GoEnrich_results_list, combined_GoEnrich_results_df = combined_GoEnrich_results_df))
}






#'
#' Main
if (TRUE) {

  # Run GoEnrich analysis
  if (TRUE) {
    output_dir <- "C:/Users/johan/VSCode_projects/bioinf_master/AssociationAnalysis/GoEnrich/runs_manual_input/v1/"
    qvalue_cutoff <- 0.1
    pvalue_cutoff <- 0.1
    query_strings <- c("pulmo", "arter", "hypertension")    # Keywords to create filtered GoEnrich results
    message(bold(cyan("\nRunning GoEnrich analysis")))
    message(underline("Output directory:\n"), output_dir)
    message(paste(underline("Querry strings: "), magenta(paste(query_strings, collapse = ", "))))
    message(paste(underline("Q-value cutoff: "), magenta(qvalue_cutoff)))
    message(paste(underline("P-value cutoff: "), magenta(pvalue_cutoff), "\n"))

    manual_gene_list <- c("ENSG00000133424,ENSG00000186716,ENSG00000069535,ENSG00000166313,ENSG00000182310,ENSG00000159958,
    ENSG00000160856,ENSG00000177455,ENSG00000184838,ENSG00000119411,ENSG00000189306,ENSG00000006756,ENSG00000132481,
    ENSG00000176293,ENSG00000182636,ENSG00000186204,ENSG00000186765,ENSG00000187912,ENSG00000103740,ENSG00000108786,
    ENSG00000116857,ENSG00000132970,ENSG00000165507,ENSG00000004478,ENSG00000006625,ENSG00000099617,ENSG00000114450,
    ENSG00000118420,ENSG00000134539,ENSG00000137265,ENSG00000138119,ENSG00000175707,ENSG00000204010,ENSG00000069188,
    ENSG00000070444,ENSG00000087237,ENSG00000107902,ENSG00000112182,ENSG00000121742,ENSG00000121807,ENSG00000122861,
    ENSG00000126217,ENSG00000130518,ENSG00000131398,ENSG00000133069,ENSG00000136854,ENSG00000138794,ENSG00000148498,
    ENSG00000164011,ENSG00000166428,ENSG00000168517,ENSG00000185052,ENSG00000188687,ENSG00000015676,ENSG00000113356,
    ENSG00000119688,ENSG00000128218,ENSG00000135919,ENSG00000136158,ENSG00000136573,ENSG00000156113,ENSG00000173272,
    ENSG00000183066,ENSG00000099864,ENSG00000105991,ENSG00000115112,ENSG00000115884,ENSG00000126266,ENSG00000149809,
    ENSG00000150594,ENSG00000163534,ENSG00000166780,ENSG00000169877,ENSG00000172164,ENSG00000173114,ENSG00000174233,
    ENSG00000179455,ENSG00000187091,ENSG00000070501,ENSG00000074527,ENSG00000076716,ENSG00000117480,ENSG00000128253,
    ENSG00000129422,ENSG00000135480,ENSG00000141294,ENSG00000142949,ENSG00000152413,ENSG00000160285,ENSG00000174307,
    ENSG00000175029,ENSG00000176463,ENSG00000178053,ENSG00000182489,ENSG00000182985,ENSG00000183018,ENSG00000197576,
    ENSG00000002933,ENSG00000006282,ENSG00000037280,ENSG00000074416,ENSG00000100448,ENSG00000107731,ENSG00000120738,
    ENSG00000124588,ENSG00000132773,ENSG00000134243,ENSG00000135476,ENSG00000139178,ENSG00000141905,ENSG00000148848,
    ENSG00000153774,ENSG00000163001,ENSG00000163600,ENSG00000170190,ENSG00000171502,ENSG00000188848")

    # Transform the manual gene list to a data frame
    manual_gene_list_df <- data.frame(gene_id = unlist(strsplit(manual_gene_list, ",")))

    GoEnrichresults <- run_enrichment_analysis(target_genes = manual_gene_list_df, output_dir = output_dir, qvalue_cutoff = qvalue_cutoff, pvalue_cutoff = pvalue_cutoff, query_strings = query_strings)
  }

  # Save the GoEnrich results df as csv
  if (FALSE) {	
    write.csv(GoEnrichresults$combined_GoEnrich_results_df, "GoEnrichresults_manual_input.csv", row.names = FALSE)
  }

}