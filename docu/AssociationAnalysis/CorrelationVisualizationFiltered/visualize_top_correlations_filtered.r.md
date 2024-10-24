# Idea
Adapted code of [[visualize_top_correlations.r]], which includes an additional filter step that mimics the aspects that were considered for the target selection within the [[PotentialHitsOfInterest | manual curation]]. 

The following aspects are implemented in the filters:
* (arbitrary) more than two third of the data points have nonzero ATAC accessibility (otherwise the correlation are very unstable because there are many zero data points, what also subverts the concept of improved target selection by inclusion of the clinical parameters)
* exclude correlation that are driven by outliers
	* Both Spearman as well as Kendall correlation was strongly driven by outliers in the exception cases, thus a separate outlier assessment was necessary.
	* outlier selection by **IQR method** i is a robust, non-parametric approach that does not assume a normal distribution, making it suitable for datasets with unknown or non-normal distributions. It helps identify outliers based on percentiles, which is less influenced by extreme values compared to methods like Z-scores.
	* check if correlation is outlier driven by - own developed method - summary draft written by chatty - summary polished by me.
		* The null hypothesis of this permutation test is that the observed difference in Spearman correlations between the original dataset and the dataset without outliers is not significantly different from what would be expected by chance. In other words, the presence of the specific number of outliers does not have a significant effect on the correlation.
		1. **Observed Difference**: Calculate the observed difference in Spearman correlations between the original dataset and the dataset with the specified number of outliers removed.
		2. **Permutation Test**:
		    
		    - **Identify Outliers**: Identify the outliers in the original dataset. 
		    - **Random Selection**: Randomly select the same number of outliers as observed in the original dataset.
		    - **Create Reduced Dataset**: Remove the selected outliers from the original dataset to create a reduced dataset.
		    - **Recalculate Correlations**: Calculate the Spearman correlation for the reduced dataset and compute the difference from the correlation of the dataset without outliers.
		    - **Repeat**: Repeat the above steps a large number of times (1000 permutations) to build a distribution of correlation differences under the null hypothesis.
		3. **Calculate p-value**: Compare the observed difference in correlations to the distribution of differences obtained from the permutations. The p-value is calculated as the proportion of permuted differences that are as extreme or more extreme than the observed difference. This p-value indicates the probability of observing a difference in correlations as extreme as the observed difference under the null hypothesis.
	  This approach performs a permutation test by systematically removing outliers from the original dataset and recalculating the Spearman correlation. The observed difference in correlations between the original dataset and the dataset without radomly removed "outliers" is compared to a distribution of differences obtained from the permutations. This method was selected because it directly addresses the impact of outliers on the correlation, providing a robust assessment of the correlation's sensitivity to outliers. By considering multiple possible outlier sets, it ensures a thorough evaluation of the significance of the observed correlation difference. The null hypothesis is that the **removal of outliers has no significant effect on the correlation**. Under this null hypothesis, the correlations for the two datasets should not differ significantly if outliers were removed randomly (i.e., there is no systematic effect).
* exclude correlation below |0.4|
* exclude correlation in which the distributions of the disease and healthy data points significantly differ
	* The **Wilcoxon rank-sum test** (also known as the **Mann-Whitney U test**) is the appropriate choice when comparing two independent groups without assuming a normal distribution, making it ideal for non-parametric data. Since your data points are independent, and you are uncertain about the distribution, this test is robust for detecting differences in central tendency without the strict assumptions of parametric tests like the t-test.


- Add the additional filter step to the existing visualization code
- Merge the resulting plot directories and remove duplicates
- use resulting geneset for further analysis