## Nonparametric tests for independence, trivariate groupwise dependence
#### following Herwartz & Maxand (2018)

rm(list=ls())

library("energy")
library("copula")
library("SpatialNP")

set.seed(123)

## number of observations
n <- 100

## dependence factor 
rho <- 0.5


## simulation of trivariate sample following normal marginal distributions ("norm") and normal, t- or Clayton copula
mv <- mvdc(normalCopula(c(0,rho,rho),dim=3,dispstr="un"), c("norm", "norm", "norm"), list(list(mean = 0, sd = 1),list(mean = 0, sd = 1),list(mean = 0, sd = 1)))
x <- rMvdc(n, mv)

#mv <- mvdc(tCopula(c(0,rho,rho),dim=3,dispstr="un"), c("norm", "norm", "norm"), list(list(mean = 0, sd = 1),list(mean = 0, sd = 1),list(mean = 0, sd = 1)))
#x <- rMvdc(n, mv)

#mv <- mvdc(claytonCopula(c(0,rho,rho),dim=3,dispstr="un"), c("norm", "norm", "norm"), list(list(mean = 0, sd = 1),list(mean = 0, sd = 1),list(mean = 0, sd = 1)))
#x <- rMvdc(n, mv)


##Kojadinovic - multivariate version of indepTest
test <- multIndepTest(x,c(2,1))
test$fisher.pvalue

## Taskinen et al
srtest <- sr.indep.test(x[,1:2],x[,3], score="rank")
srtest$p.value

## Szekely/Rizzo

dcov <- indep.test(x[,1:2],x[,3], method = "dcov", R=1000)
dcov$p.value

## Wilks lambda test statistic
LR <- -n*log(det(cov(x))/(var(x[,3])*det(cov(x[,1:2]))))
