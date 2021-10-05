
rm(list=ls())

setwd("/Users/oscar/Desktop/causal_lasso_r")
library("car")
library("Matching")
library("genlasso")
source("auxiliary.R")
source('abadie_model.r')
library(latex2exp)

## experimental data
data("lalonde")
head(lalonde)

y = lalonde$re78
z = lalonde$treat
x = as.matrix(lalonde[, c("age", "educ", "black",
                          "hisp", "married", "nodegr",
                          "re74", "re75")])
n = length(z)
x = scale(x)
lin = lm(y ~ z)
coef(lin)[2]
round(sqrt(hccm(lin, type = "hc2")[2, 2]), 3)


prgsore  = lm(y ~ x, weights = 1 - z)$fitted.values
prgsore1 = prgsore[z==1]
prgsore0 = prgsore[z==0]
hist(prgsore)
hist(prgsore1)
hist(prgsore0)

#Compute proxys for individual treatment effects:
tau.m    = z
y1 = y[z==1]
y0 = y[z==0]
for(i in 1:n)
{
  if(z[i]==1)
    tau.m[i] = y[i] - y0[which.min(abs(prgsore[i] - prgsore0))]
  
  if(z[i]==0)
    tau.m[i] = y1[which.min(abs(prgsore[i] - prgsore1))] - y[i]
}
hist(tau.m)


prgsore.order = prgsore[order(prgsore)]
tau.order = tau.m[order(prgsore)]



beta_hat  =fused_lasso_bic(tau.order)
plot(tau.order ~ prgsore.order, type = "s", 
     col = "grey", lty = 2)
lines(beta_hat  ~ prgsore.order, type = "s")
abline(h = 0,col="red")





##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
##################################################################################################################################################################################################################
## observational data




dat <- read.table("cps1re74.csv", header = TRUE)
# unemployed
dat$u74 <- as.numeric(dat$re74==0)
dat$u75 <- as.numeric(dat$re75==0)
y = dat$re78
z = dat$treat
x = as.matrix(dat[, c("age", "educ", "black",
                      "hispan", "married", "nodegree",
                      "re74", "re75", "u74", "u75")])

psore  = glm(z ~ x, family = binomial)$fitted.values


psore1 = psore[z==1]
psore0 = psore[z==0]
hist(psore)
hist(psore1)
hist(psore0)


tau1.m  = psore1
n1 = length(tau1.m)
y1 = y[z==1]
y0 = y[z==0]
for(i in 1:n1)
  tau1.m[i] = y1[i] - y0[which.min(abs(psore1[i] - psore0))]

hist(tau1.m)




psore1.order = psore1[order(psore1)]
tau1.order = tau1.m[order(psore1)]


beta_hat  =fused_lasso_bic(tau1.order)

plot(tau1.order ~ psore1.order, type = "s", 
     col = "grey", lty = 2,cex.lab=1.5,cex.main =2.4,cex.axis=1.4,xlab ="Propensity score",ylab= TeX("$\\tilde{\\tau}^{(1)}$"))
#axis(1, at=c(1,floor(n/2),n), labels=c(prgsore.order[1],prgsore.order[floor(n/2)],prgsore.order[n]),cex = 1.8)
lines(beta_hat  ~ psore1.order, type = "s")
abline(h = 0,col="red",lty=2)