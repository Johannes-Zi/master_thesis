Regression version with LeaveOneOut outer cross validation instead of performing the outer cross validation with a bigger subset of the actual dataset
	- this might be better because we are working with a smaller dataset and the training might be compromised to much by using a bigger set for the outer cv - but we need the outer cv for a reliable performance estimation

For a ==detailed workflow description== see the Notion documentation

# Description
This Version analyzes the gene-specific lambda.min models the were trained during the Regression 
##### paths:
```
# Linux
/projects/pulmonary_hypertension/work/analysis/clinical_data_stichit_run/regression/regression_output/

# Windows
\\OneDrive\\dateien_cloud\\Master\\Semester_4\\Masterarbeit\\data\\pulmanory_hypertension\\regression\\leaveOneOut_regression\\
```

# File overview
==For detailed in/output documentation see Notion docu==
### input
- \*\_Pearson.txt
- \*\_Spearman.txt
### output

* Performance_Overview.txt
*  \*\_Pearson.RData
- \*\_Spearman.RData
[[ElNetCoeffCounts_elasticnet_modelStructure]]
 