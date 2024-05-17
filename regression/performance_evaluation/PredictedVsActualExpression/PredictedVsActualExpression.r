library(glmnet)
library(ggplot2)

segmentation_path_1 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000125629_10/Segmentation_ENSG00000125629_10_Pearson.txt"
elnet_path_1 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000125629_10/Elasticnet_Regression_Model_Segmentation_ENSG00000125629_10_Pearson.RData"

segmentation_path_2 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000156232_10/Segmentation_ENSG00000156232_10_Pearson.txt"
elnet_path_2 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000156232_10/Elasticnet_Regression_Model_Segmentation_ENSG00000156232_10_Pearson.RData"

segmentation_path_3 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000249092_10/Segmentation_ENSG00000249092_10_Pearson.txt"
elnet_path_3 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000249092_10/Elasticnet_Regression_Model_Segmentation_ENSG00000249092_10_Pearson.RData"

segmentation_path_4 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000259056_10/Segmentation_ENSG00000259056_10_Pearson.txt"
elnet_path_4 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000259056_10/Elasticnet_Regression_Model_Segmentation_ENSG00000259056_10_Pearson.RData"

segmentation_path_5 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000211788_10/Segmentation_ENSG00000211788_10_Pearson.txt"
elnet_path_5 <- "C:/Users/johan/Desktop/PredictedVsActualExpression/ENSG00000211788_10/Elasticnet_Regression_Model_Segmentation_ENSG00000211788_10_Pearson.RData"


# Load elasticnet model
load(elnet_path_5)
# Load segmentation data
M<-read.table(segmentation_path_5,header=TRUE,sep="",row.names=1)
print(M)

# Normailize
M<-unique(M)
M<-data.frame(M)
M<-log2(M+1)
M<-data.frame(scale(M,center=TRUE, scale=TRUE))
print(M)

M_new <- M[, -ncol(M)]


lambda_min <- elasticnet_model$model$lambda.min
#coefficients <- coef(elasticnet_model$model, s = elasticnet_model$model$lambda.min)


predictions <- predict(elasticnet_model$model, newx = as.matrix(M_new), s = lambda_min)
print(predictions)

last_column_M = M[, ncol(M)]
df = data.frame(last_column_M, predictions)
print(df)


# Create a dotplot with a line of best fit
ggplot(df, aes(x = last_column_M, y = predictions)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE, color = "red") +
    xlab("Normalized actual expression values") +
    ylab("Normalized Predictions") +
    ggtitle("Predicted vs Actual expression values ENSG00000211788") +
    theme(plot.title = element_text(hjust = 0.5))

    ggsave("plot.png", plot = last_plot(), device = "png", path = getwd())