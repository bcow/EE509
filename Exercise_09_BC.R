load("data/Lab9.RData")

## Maximum Likelihood Poisson Regression

PR1 = glm(y ~ TDR, family=poisson(link="log"))

ic=c(0,0)
LnL = function(beta){
  -sum(dpois(y,exp(beta[1] + beta[2]*TDR),log=TRUE))
}
PR2 = nlm(LnL,ic)

### Lab Report Task 1

# 1.  Plot seedling densities as a function of TDR

plot(TDR,y, xlab="TDR", ylab="Seedling Counts")

# 2.  Fit the Poisson regression model using one of the methods above and turn in  the summary output (hint: for the second method you will need to define the initial condition vector ic) 

# Using approach 1
PR1 = glm(y ~ TDR, family=poisson(link="log"))
coef(PR1)

PR1.fit <- sort(fitted(PR1))
PR1.idx <- sort(fitted(PR1), index.return=TRUE)[[2]]

lines(TDR[PR1.idx], PR1.fit, lwd = 2)

sPR1 <- summary(PR1)
x <- seq(0,max(TDR),length=100)
ypred <- predict(PR1, list(x= x), type="resp" )

plot(x,ypred)

# 3.	Add regression lines to the plot
#     Hint 1: use “coef” to extract the regression coefficients from the GLM.
#     Hint 2: don't forget about the link function when plotting the lines


# 4.	Briefly describe how you would add model confidence and predictive intervals to these curves
# 5.	What would be an appropriate null model to compare to?  What metric would you use to compare the two models?
# 6.	Plot the calibration data of TDR vs. soil moisture.  Fit a Normal regression model to the calibration data, add the line to the plot, and report the summary table information for the fit