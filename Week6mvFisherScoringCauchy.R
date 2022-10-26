setwd("M:/aStatsMeth2022-23")

# Cauchy distribution

# unknown parameter
#theta[1]= alpha
#theta[2] = beta

# log likelihood funtion
 
logLik  <- function(theta, y)
{
   return(length(y)*log(theta[2])- sum(log(theta[2]^2+(y-theta[1])^2 ))  )
}

# Score vector
score <- function(theta, y)
{
 U1 =    2*sum( (y-theta[1])/(theta[2]^2+((y-theta[1]))^2 ))
 U2 =  n/theta[2] - 2*sum( theta[2]/(theta[2]^2+(y-theta[1])^2 ))

 return(c(U1, U2))
}

Fisher.info.scale <- function(theta, n)
{
 0.5*n/theta[2]^2
}

#=========================
# sim data 

#n=20 
#y = rt(n, 1)+2

#> y
# [1]  1.299800  4.761788  2.480328 -7.743641  2.403356  3.319482  3.057658
# [8] 10.591182  1.674542  2.017324  3.426539  1.536969  3.337224  1.857916
#[15] -8.308324  1.357775  2.680424  1.488331  1.693336  2.050363

#save(y, file = "Week7DataCauchy.load" )
# ===========================
# load data 

load(file = "Week6DataCauchy.load")
n = length(y)

# Find the MLE using Fisher's methof of scoring

theta0=c(0,1)
eps= 0.00001

theta = theta0
while(sum(score(theta, y)^2) > eps^2) theta = theta + score(theta, y)/Fisher.info.scale(theta, n)

# the MLE: 
theta
> theta
 2.0844634 0.7242023


theta.hat = theta

# Values of score vector, scale of the Fisher information matrix and the matrix
score(theta, y)
Fisher.info.scale(theta, n)

# Fisher information matrix  

Fisher.matrix = Fisher.info.scale(theta, n)*matrix(c(1,0,0,1), ncol=2, nrow=2)

> score(theta, y)
[1]  1.325544e-07 -6.330120e-06
> Fisher.info.scale(theta, n)
[1] 19.06691
> 

# estimated standard error of both parameters

ESE = 1/sqrt(Fisher.info.scale(theta, n))
ESE
 0.2290129

# Marginal 95\% CI 

# of alpha  
theta.hat[1]+c(-1,1)*1.96*ESE
> 1.635598 2.533329

# of beta
theta.hat[2]+c(-1,1)*1.96*ESE
> 0.275337 1.173067

# Reject H0 that theta = 1 at 5% significance level? no

#===========================
# now we run the Newton Raphson method using maxLik function in R.

maxLik(logLik, start = theta0, y=y)

> maxLik(logLik, start = theta0, y=y)
Maximum Likelihood estimation
Newton-Raphson maximisation, 9 iterations
Return code 1: gradient close to zero (gradtol)
Log-Likelihood: -21.0248 (2 free parameter(s))
Estimate(s): 2.084463 0.7242019

# the same values of the MLEs
 
