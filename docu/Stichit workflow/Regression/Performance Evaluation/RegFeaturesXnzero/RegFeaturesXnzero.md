## Description:
Simply iterate over the .RData output files of the regression process. The files hold the elasticnet_model 's that were trained for the specific gene.
What the script does, is selecting for each gene the elasticnet_model which had the smallest lambda (lambda.min) and thus performed best at the cross validation step.
Then it is determined how many non-zero feature coefficients are in the min lambda model and how many features are there in total.
This is then done for all the gene specific models and plotted as 2d density plot

You can see here the number of non-zero feature coefficients vs the total number of features for each gene specific regression model.
## Code
https://github.com/Johannes-Zi/master_thesis/blob/main/regression/performance_evaluation/RegFeaturesXnzero/RegFeaturesXnzero.r

## Runs performed

| Feature segment preselection after segmentation | Standard CV        | LeaveOneOut CV     |
| ----------------------------------------------- | ------------------ | ------------------ |
| Spearman correlation based                      | [[fig_2504241636]] | [[fig_2504241638]] |
| Pearson correlation based                       | [[fig_2504241635]] | [[fig_2504241637]] |

## Initial Idea 
  * see tutorial https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
  * create a distribution for all elastic net models how many features they have and how many of them are actually zero
	  * maybe a dot plot that has on the one axis the number of features and on the other axis the percentage of the features that are zero
		  * if this is not working because there are two man y-dots with the exact location, there could be a density used or simply a bar chart created with the average percentage of zero for each of the abundances of feature counts
	  * use the respective models with the minimum lambda(like in the analysis of the exception cases) for the nonzero coeff counting counting (total number of features is in principle the same like in the segmentation output)