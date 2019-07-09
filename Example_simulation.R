
#######  NTP-MEM bivariate

source('functions.R')
library(mnormt)

oldw <- getOption("warn")
options(warn = -1)


d<-2
mu.xis <- 1
sigma2.xis <- 1 
alphas <- c(0,0)
betas <- c(1,1)
sigma2.us <- 0.5  
sigma2.es <- diag(c(0.25,0.25))
gammas <- 0.5

z=sample.MEM(n=100, alphas, betas, mu.xis, sigma2.xis, sigma2.us, sigma2.es, gammas)    # generate a random sample
theta=EM.NTP_MEM(z)   # estimation via EM algorithm
MI=MIapprox2(theta,z)  # information matrix
ep=sqrt(diag(solve(MI)))  # standard error of the parameters

