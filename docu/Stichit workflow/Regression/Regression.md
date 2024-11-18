Stichtit Regression step - based on pre-selected Segments of segmentation.
Performs OLS and [[ElasticnetBasics | elastic-net]] (using [[glmnetBasics | glmnet]])

## Performed runs
[[RegRuns]]
## Performance evaluation
[[RegPerformanceEvaluation]]
[[OnlyZeroCoeffOLS]]
## Related questions
[[QuestionCollectionStichitRegression]]

# Additional context infomration - currently not in thesis
In alter stichit version zusätzlich zu initialen ELnet:

2. ein OLs schritt, welcher die präselktierten top elemente der elastic nets als features verwendet und dann erneut ein modell trainiert. Nur diesmal halt nicht auf 80% sonder 100% der Daten.
