type.measure loss to use for cross-validation. Currently five options, not all available for all
models. The default is type.measure="deviance", which uses squared-error
for gaussian models (a.k.a type.measure="mse" there), deviance for logistic
and poisson regression, and partial-likelihood for the Cox model. type.measure="class"
applies to binomial and multinomial logistic regression only, and gives misclas-
sification error. type.measure="auc" is for two-class logistic regression only,
and gives area under the ROC curve. type.measure="mse" or type.measure="mae"
(mean absolute error) can be used by all models except the "cox"; they measure
rel’s concordance measure, only available for cox models.

==In Stichit is the standard option for loss measure used==