## Used file-types
### Performance_Overview.txt
Generated in the context of the Stichit regression, representing a summary of  performances of all trained models.  [[Performance_Overview.txt | Detailed description.]]
### Elastic-net models (\*.RData) created during Stichit regression
[[ElNetCoeffCounts_elasticnet_modelStructure]]

## Plot types
Marcel told to focus only on the Pearson driven feature (segment) preselection at the end of the segmentation and the Pearson values of the cross validation (Spearman less relevant). Laura thinks this might be because we use a linear model and Pearson measures linear relationship. I think this might be maybe because we are searching for Segment that directly have an influence on the gene Expression - if there are non linear relationships, there might be other processes involved and thus the effect might be more difficult to understand based on the data we have - we are focusing on detecting simple DNA accessibility driven correlations.
### Pearson values multi-violin
[[PearsonValMulti]]
### MSE vaLues multi-violin
[[MSEValMulti]]
### MSE vaRiation multi-violin
[[MSEVarMulti]]
### MSE values x MSE variance paired-dotplot
[[MSEValXMSEVar]]
### Person values x Pearson variance paired-dotplot
[[PearsonValXVariance]]
### Number of features x number of nonzero feature coefficients
[[RegFeaturesXnzero]]
### Number of nonzero feature coefficients x Pearson values
[[NzeroXpearson]]

## Finished plots
[[RegPlotOverview.canvas|Overview]]