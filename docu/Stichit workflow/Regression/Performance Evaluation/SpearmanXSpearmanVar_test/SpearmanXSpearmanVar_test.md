==This plot was generated just out of interest==
## Description:
Plots the Spearman correlations that were computed for specific genes together with the variance between the correlations that were computed for the samples of the cross validation dataset. This plot can be only generated for standard cross validation. For the leaveOneOut cross validation can be no variance of the Spearman values be calculated when there is only one sample in the outer cross validation dataset and thus only one correlation calculated.

(plots in a range IQR +- 1.5 * IQR)

## Observations
* If the RNA Seq counts are pretty small, it can occur that during the normalization, upstream of the ElNet training, the values can be transformed to negative representations within the "normalized values training space" and this can lead to negative coefficients.
## Code
https://github.com/Johannes-Zi/master_thesis/blob/main/regression/performance_evaluation/SpearmanXSpearmanVar_test/SpearmanXSpearmanVar_test.r

## Runs performed
plots were only generated for the Pearson based feature (segment) preselection at the end of the segmentation process

| Feature segment preselection after segmentation | Standard CV      |
| ----------------------------------------------- | ---------------- |
| Pearson correlation based                       | [[fig_06051315]] |
