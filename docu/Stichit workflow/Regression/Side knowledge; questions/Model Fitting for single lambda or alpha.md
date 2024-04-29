Here's a step-by-step breakdown of what happens:

1. **Initialize the coefficients**: All coefficients are initialized, typically to zero.
    
2. **Calculate the cost**: The cost function is calculated. This includes the loss (how well the model fits the data) and the penalty term. The penalty term is a combination of L1 (Lasso) and L2 (Ridge) penalties, controlled by the `alpha` parameter. The `lambda` parameter controls the overall strength of the penalty.
    
3. **Update the coefficients**: The coefficients are updated in a way that reduces the cost. This involves calculating the gradient of the cost function with respect to each coefficient, and then adjusting the coefficient in the direction that reduces the cost. The size of the adjustment is determined by the learning rate. (==see Gradient descend==)
    
4. **Repeat steps 2 and 3**: Steps 2 and 3 are repeated until the cost stops decreasing significantly or a maximum number of iterations is reached. This process is also known as gradient descent.
    
5. **Store the model**: Once the algorithm has converged, the final set of coefficients is stored as the model for that `lambda`.