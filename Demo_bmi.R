rm(list=ls())
setwd("/Users/oscar/Desktop/causal_lasso_r")
library("car")
library("genlasso")
source("auxiliary.R")
## an application: the NHANES BMI dataset
library("ATE")
library("latex2exp")
data(nhanes_bmi)
#Treatment indicators:
z = nhanes_bmi$School_meal
#Response:
y = nhanes_bmi$BMI
#Covariates:
x = as.matrix(nhanes_bmi[, -c(1, 2)])
n = length(z)
x = scale(x)
lin = lm(y ~ z)
coef(lin)[2]
round(sqrt(hccm(lin, type = "hc2")[2, 2]), 3)


tau.m    = z
y1 = y[z==1]
y0 = y[z==0]

#proprnsity scores:
psore  = glm(z ~ x, family = binomial)$fitted.values
psore1 = psore[z==1]
psore0 = psore[z==0]
hist(psore)
hist(psore1)
hist(psore0)

###proxy for treatment effects:


for(i in 1:n)
{
  if(z[i]==1)
    tau.m[i] = y[i] - y0[which.min(abs(psore[i] - psore0))]
  
  if(z[i]==0)
    tau.m[i] = y1[which.min(abs(psore[i] - psore1))] - y[i]
}
hist(tau.m)


# Order propensity scores and proxys of individual treatment effects:

psore.order = psore[order(psore)]
tau.order = tau.m[order(psore)]

#Compute estimates:
beta_hat  =fused_lasso_bic(tau.order)

#Plot estimates
plot(tau.order ~ psore.order, type = "s", 
     col = "grey", lty = 2)
lines(beta_hat  ~ psore.order, type = "s")
abline(h = 0,col="red")
