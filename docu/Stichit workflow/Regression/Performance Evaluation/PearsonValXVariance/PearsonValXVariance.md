## Description:
Plots the Pearson correlations that were computed for specific genes together with the variance between the correlations that were computed for the samples of the cross validation dataset. This plot can be only generated for standard cross validation. For the leaveOneOut cross validation can be no variance of the Pearson values be calculated when there is only one sample in the outer cross validation dataset and thus only one correlation calculated.

(plots in a range IQR +- 1.5 * IQR)
## Code
https://github.com/Johannes-Zi/master_thesis/blob/main/regression/performance_evaluation/PearsonValXVariance/PearsonValXVariance.r

## Runs performed
plots were only generated for the Pearson based feature (segment) preselection at the end of the segmentation process

| Feature segment preselection after segmentation | Standard CV      |
| ----------------------------------------------- | ---------------- |
| Pearson correlation based                       | [[fig_03051923]] |
