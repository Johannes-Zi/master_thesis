standardized Glmnet format - easy to access via the respective package functionalities.

```
# Example contetn of elasticnet_model

> str(elasticnet_model)
List of 6
 $ model      :List of 14
    ..$ lambda    : num [1:87] 8.72 7.94 7.24 6.59 6.01 ...  # lambda values used in the model
    ..$ cvm       : num [1:87] 0.999 0.999 0.998 0.995 0.991 ...  # cross-validated mean value
    ..$ cvsd      : num [1:87] 0.26 0.261 0.261 0.262 0.263 ...  # cross-validated standard deviation
    ..$ cvup      : num [1:87] 1.26 1.26 1.26 1.26 1.25 ...  # upper bound of the cross-validated values
    ..$ cvlo      : num [1:87] 0.74 0.739 0.736 0.733 0.728 ...  # lower bound of the cross-validated values
    ..$ nzero     : Named int [1:87] 0 1 1 3 5 5 6 6 7 8 ...  # number of non-zero coefficients for each lambda value
    .. ..- attr(*, "names")= chr [1:87] "s0" "s1" "s2" "s3" ...
    ..$ call      : language cv.glmnet(x = x_train, y = y_train, nfolds = as.numeric(argsL$innerCV),      keep = TRUE, alpha = x)  # function call used to create the model
    ..$ name      : Named chr "Mean-Squared Error"
    .. ..- attr(*, "names")= chr "mse"
    ..$ glmnet.fit:List of 12
    .. ..$ a0       : Named num [1:87] 4.34e-17 4.30e-17 4.25e-17 4.21e-17 4.10e-17 ...  # intercept values for each lambda value
    .. .. ..- attr(*, "names")= chr [1:87] "s0" "s1" "s2" "s3" ...
    .. ..$ beta     :Formal class 'dgCMatrix' [package "Matrix"] with 6 slots
    .. .. .. ..@ i       : int [1:650] 5 5 1 3 5 1 2 3 5 7 ...  # coefficient values for each lambda value
    .. .. .. ..@ p       : int [1:88] 0 0 1 2 5 10 15 21 27 34 ...
    .. .. .. ..@ Dim     : int [1:2] 8 87
    .. .. .. ..@ Dimnames:List of 2
    .. .. .. .. ..$ : chr [1:8] "chr6.41070974.41070993" "chr6.41071204.41071343" "chr6.41073354.41073453" "chr6.41073454.41073603" ...
    .. .. .. .. ..$ : chr [1:87] "s0" "s1" "s2" "s3" ...
    .. .. .. ..@ x       : num [1:650] 4.54e-03 9.43e-03 2.25e-03 9.02e-06 1.46e-02 ...
    .. .. .. ..@ factors : list()
    .. ..$ df       : int [1:87] 0 1 1 3 5 5 6 6 7 8 ...  # number of degrees of freedom for each lambda value
    .. ..$ dim      : int [1:2] 8 87
    .. ..$ lambda   : num [1:87] 8.72 7.94 7.24 6.59 6.01 ...  # lambda values used in the model
    .. ..$ dev.ratio: num [1:87] 0 0.00404 0.00833 0.01444 0.02879 ...  # deviance ratio
    .. ..$ nulldev  : num 20  # null deviance
    .. ..$ npasses  : int 523  # number of passes
    .. ..$ jerr     : int 0  # error code
    .. ..$ offset   : logi FALSE  # offset flag
    .. ..$ call     : language glmnet(x = x_train, y = y_train, alpha = x)  # function call used to create the model
    .. ..$ nobs     : int 21  # number of observations
    .. ..- attr(*, "class")= chr [1:2] "elnet" "glmnet"
 $ alpha      : num 0.05  # alpha parameter used in the model
 $ lambda     : num 1.63  # lambda parameter used in the model
 $ inner_folds:'data.frame':    21 obs. of  2 variables:
    ..$ sample_id: chr [1:21] "healthy-6818" "pah-6858" "ph-lung-6912" "healthy-6658" ...  # sample IDs
    ..$ fold_id  : int [1:21] 3 2 6 5 6 1 2 4 4 6 ...  # fold IDs
 $ sample_ids : chr [1:5] "pah-6872" "ph-lung-6885" "pah-6623" "ph-lung-6678" ...  # sample IDs used in the model
 $ outer_folds: int [1:5] 3 3 3 3 3  # number of outer cross-validation folds

# Copilot description for the elements within this structure:
# The elasticnet_model structure is a list with the following elements:

# model: A list containing various information about the model. It includes:
#   lambda: A numeric vector representing the lambda values used in the model.
#   cvm: A numeric vector representing the cross-validated mean value.
#   cvsd: A numeric vector representing the cross-validated standard deviation.
#   cvup: A numeric vector representing the upper bound of the cross-validated values.
#   cvlo: A numeric vector representing the lower bound of the cross-validated values.
#   nzero: A named integer vector representing the number of non-zero coefficients for each lambda value.
#   call: The function call used to create the model.
#   name: The name of the model.
#   glmnet.fit: A list containing information about the fitted model. It includes:
#     a0: A named numeric vector representing the intercept values for each lambda value.
#     beta: A sparse matrix representing the coefficient values for each lambda value.
#     alpha: A numeric value representing the alpha parameter used in the model.

# lambda: A numeric value representing the lambda parameter used in the model.

# inner_folds: A data frame containing information about the inner cross-validation folds. It includes:
#   sample_id: A character vector representing the sample IDs.
#   fold_id: An integer vector representing the fold IDs.

# sample_ids: A character vector representing the sample IDs used in the model.

# outer_folds: An integer vector representing the number of outer cross-validation folds.

# Each element provides specific information about the elastic net model, such as lambda values, cross-validated metrics, coefficient values, and sample IDs.
```