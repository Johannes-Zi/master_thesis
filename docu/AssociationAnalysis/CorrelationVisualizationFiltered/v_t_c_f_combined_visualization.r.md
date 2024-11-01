# Idea
takes the filtered genes per parameter from [[visualize_top_correlations_filtered.r]] and sets them into context with the clinical parameters.

For the "number of gene occurrences" plot and the string analysis were performed with a reduced set of clinical parameters. 
The following parameters were excluded:
"age", "height", "weight", "NYHA", "ZVD", "heart.rate", "Paradoxe_Septumbewegung", "Perikarderguss"
This was done because either the parameters are not specific for the disease (e.g. age), or the distribution of the respective clinical values were not suitable to perform meaningful correlations (e.g. NYHA values 1,2,3,4). Interesting about this that the exclusion of those 8 Parameters reduced the gene set only by two genes, suggesting that those correlations were driven by a meaningful ATAC value distribution, less than meaningful clinical parameter values.
The SMW parameter was also removed, because there are only 7 data points represented, leading to many false positives - also very relevant point to mention in the thesis! the number of genes dropped from 123 to 117, but the clusters still remain.

The heatmap was created to se what clinical parameter share common genes in the filtered dataset. The map also displays how many of the total genes are represented in each clinical parameter.

![[number_of_entries_per_param.png]]

![[gene_frequency_distribution.png]]
![[string_db_graph.png]]
![[clinical_parameter_gene_associations_heatmap.png]]
### Heat-map Creation Process

1. **Reorder Percentage Matrix**:
   - Reorders the percentage matrix using hierarchical clustering to group similar clinical parameters together.

2. **Create Cross-Tabulation**:
   - Generates a matrix where rows represent clinical parameters and columns represent gene IDs, with cells indicating the presence (1) or absence (0) of a gene ID for a clinical parameter.

3. **Calculate Co-Occurrences**:
   - Performs matrix multiplication between the cross-tabulation matrix and its transpose to get the co-occurrence matrix.

4. **Calculate Percentages**:
   - Divides the co-occurrence matrix by the total number of gene IDs and multiplies by 100 to get the percentage matrix.

5. **Reorder the Percentage Matrix**:
   - Calls the `reorder_percentage_matrix` function to reorder the percentage matrix based on hierarchical clustering.

6. **Convert Matrix to Data Frame**:
   - Converts the reordered percentage matrix to a data frame for plotting.

7. **Reverse Y-Axis Order**:
   - Reorders the factor levels of the y-axis in reverse order to ensure the diagonal goes from the upper left to the lower right.

8. **Set Values Above Diagonal to NA**:
   - Sets the percentage values above the diagonal to NA to avoid redundancy.

9. **Create the Heatmap**:
   - Uses `ggplot2` to create the heatmap with appropriate aesthetics and themes, rotating x-axis labels by 45 degrees, centering the title, setting the plot background to white, and removing axis labels.
### String input
```
ENSG00000133424,ENSG00000186716,ENSG00000069535,ENSG00000166313,ENSG00000182310,ENSG00000159958,ENSG00000160856,ENSG00000177455,ENSG00000184838,ENSG00000119411,ENSG00000189306,ENSG00000006756,ENSG00000132481,ENSG00000176293,ENSG00000182636,ENSG00000186204,ENSG00000186765,ENSG00000187912,ENSG00000103740,ENSG00000108786,ENSG00000116857,ENSG00000132970,ENSG00000165507,ENSG00000004478,ENSG00000006625,ENSG00000099617,ENSG00000114450,ENSG00000118420,ENSG00000134539,ENSG00000137265,ENSG00000138119,ENSG00000175707,ENSG00000204010,ENSG00000069188,ENSG00000070444,ENSG00000087237,ENSG00000107902,ENSG00000112182,ENSG00000121742,ENSG00000121807,ENSG00000122861,ENSG00000126217,ENSG00000130518,ENSG00000131398,ENSG00000133069,ENSG00000136854,ENSG00000138794,ENSG00000148498,ENSG00000164011,ENSG00000166428,ENSG00000168517,ENSG00000185052,ENSG00000188687,ENSG00000015676,ENSG00000113356,ENSG00000119688,ENSG00000128218,ENSG00000135919,ENSG00000136158,ENSG00000136573,ENSG00000156113,ENSG00000173272,ENSG00000183066,ENSG00000099864,ENSG00000105991,ENSG00000115112,ENSG00000115884,ENSG00000126266,ENSG00000149809,ENSG00000150594,ENSG00000163534,ENSG00000166780,ENSG00000169877,ENSG00000172164,ENSG00000173114,ENSG00000174233,ENSG00000179455,ENSG00000187091,ENSG00000070501,ENSG00000074527,ENSG00000076716,ENSG00000117480,ENSG00000128253,ENSG00000129422,ENSG00000135480,ENSG00000141294,ENSG00000142949,ENSG00000152413,ENSG00000160285,ENSG00000174307,ENSG00000175029,ENSG00000176463,ENSG00000178053,ENSG00000182489,ENSG00000182985,ENSG00000183018,ENSG00000197576,ENSG00000002933,ENSG00000006282,ENSG00000037280,ENSG00000074416,ENSG00000100448,ENSG00000107731,ENSG00000120738,ENSG00000124588,ENSG00000132773,ENSG00000134243,ENSG00000135476,ENSG00000139178,ENSG00000141905,ENSG00000148848,ENSG00000153774,ENSG00000163001,ENSG00000163600,ENSG00000170190,ENSG00000171502,ENSG00000188848
```