Als erstes eine df erzeugen, welcher die ols feature and nzero counts läd und dann die dazugehörigen elasticnet modelle -gene specific
dann die modelle mit nzero = anzahl der features und guten pvalues in den elasticnet anschauen und die coeffs der modelle + die Expressions im datensatz anschauen

erst mal mit der standard regression testen - war da glaube ich aufgefallen

## Description:

* There are OLS Models that have on zero coefficients, but have a good elastic-net performance - ==**we want to get an idea whats happening because it is against the expectation**==
	* Create scatterplot with predicted and actual expression values for 5-10 of the gene-specific OLS models
	* How to create the plot?
		* .RData represent the models (List - first entry is model (specific glmnet datatype))
			* the dataframe contains all models of all alpha and all lambdas that were generated/ evaluated during the training process
		* Use glmnet onboard shipped functions 
			- ==see tutorial== https://cran.r-project.org/web/packages/glmnet/glmnet.pdf
			* predict function
				* perform it for the 5-10 genes for all 26 samples
				* use min-lambda
					* There is a function/ parameter in glmnet to extract the model with the minimum lambda
					* The model with the minimum lambda represents the model that had the smallest cross validation MSE(predicted and actual RNASeq values) across all the different outer folds
						* the best model is evaluated bottom up - first between the inner folds, that are used for the parameter tuning via alpha and the the best performing model is used to predict the outer fold (the innerfold is obviously also fold specific thus the model is individual)
					* the model with the minimum lambda is not always the best - it can also be a crappy model that is very good at predicting a crappy outer fold test expression (eg. 0) - the more samples the lower the probability to get such crappy models
				* use as input the output of the segmentation: 
					* ATAC and RNASeq
						* samples as rows and columns are the regions(represent  optimized reduced intrinsic variation of the ATACSeq signal as segements)
						* last row are the RNASeq expression values
					* ==the segmentation output has to be transformed/scaled to work as input for the model!!==
						* the ATAC has to be and  the RNA can be transfomed(if you want want to compare the predicted directly) ( in the same way)
						* The transformation can be seen in the RScript around Line 191
						* The RScript functions can be directly applied on the segmentation dataframe
						* The dataframe has also be normalized with the same scaling parameters like used during the training - (mean and varaince - has to be retreived by running the traing again and printing it - is currently not stored in the output)
						* the transfomed dataframe can be then used as input for the model based prediction
						* Stichit scales, centers and log transforms the data
				* the predicted RNASeq expression values can be directly compared to the transformed RNASeq values, or can be transformed back to be compared to the original values from the segmentation. (to compare the log transformed ones is a bit nicer vor the visualization)
## Code
https://github.com/Johannes-Zi/master_thesis/blob/main/regression/performance_evaluation/NzeroXpearson/NzeroXpearson.r

## Runs performed
plots were only generated for the Pearson based feature (segment) preselection at the end of the segmentation process

| Feature segment preselection after segmentation | Standard CV      | LeaveOneOut CV   |
| ----------------------------------------------- | ---------------- | ---------------- |
| Pearson correlation based                       | [[fig_06051532]] | [[fig_06051540]] |

