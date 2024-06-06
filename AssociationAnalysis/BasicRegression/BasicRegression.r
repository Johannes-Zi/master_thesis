# Create two vectors
x <- c(1, 2, 3, 4, 5)
y <- c(2, 3, 4, 5, 6)

# Compute Pearson correlation and p-value
pearson_test <- cor.test(x, y, method = "pearson")
print(paste("Pearson correlation:", pearson_test$estimate))
print(paste("Pearson p-value:", pearson_test$p.value))

# Compute Spearman correlation and p-value
spearman_test <- cor.test(x, y, method = "spearman")
print(paste("Spearman correlation:", spearman_test$estimate))
print(paste("Spearman p-value:", spearman_test$p.value))


