## Section 2.3
## Iterative estimation of the MLE: Single parameter case
## Example page 25 -- Truncated Poisson distribution



n <- 2423
ybar <- 3663 / n
U <- function(theta) 
{
  n * ybar / theta - n - n * exp(-theta) / (1 - exp(-theta)) 
}
J <- function(theta) 
{
  n * ybar / theta^2 - n * exp(-theta) / (1 - exp(-theta))^2 
}

theta0 <- ybar
U(theta0)
J(theta0)

#===============================================================
## Newton-Raphson
#===============================================================
theta1 <- theta0 + U(theta0) / J(theta0)
U(theta1)
J(theta1)

theta2 <- theta1 + U(theta1) / J(theta1)
U(theta2)
J(theta2)

#keep repeating this until the derivative is zero (or very close to zero).
#It happens at 5th iteration



##---------- Automating the process ------------------

newraph.TP <- function(n, sumy, eps = 1e-06){
  ybar <- sumy / n
  U <- function(theta)
    n * (ybar / theta) - n - n * (exp(-theta) / (1 - exp(-theta)))
  J <- function(theta)
    n * (ybar / theta^2) - n * (exp(-theta) / (1 - exp(-theta))^2)
  theta <- ybar
  while (abs(U(theta)-0)>eps){
    theta <- theta + U(theta) / J(theta)
  }
  list(theta=theta, U=U(theta), J=J(theta))
}

newraph.TP(n=2423, sumy=3663)
#0.8924961


#===============================================================
# Fisher's method of scoring for truncated Poisson distribution 
#===============================================================

fisherscoring.TP <- function(n, sumy, eps = 1e-06, theta0= sumy/n){
  ybar <- sumy / n
  U <- function(theta)
    n * (ybar / theta) - n - n * (exp(-theta) / (1 - exp(-theta)))
  I <- function(theta)
    n / (theta*(1 - exp(-theta))) - n * (exp(-theta) / (1 - exp(-theta))^2)
  theta <- theta0
  while (abs(U(theta))>eps)
  {
    theta <- theta + U(theta) / I(theta)
  }
  list(theta=theta, U=U(theta), I=I(theta))
}

fisherscoring.TP(n=2423, sumy=3663)


theta = seq(0.1, 2,0.01)
theta = seq(0.5, 1,0.01)
theta = seq(0.8, 1,0.01)
theta = seq(0.85, 0.95,0.001)
plot(theta, U(theta))
abline(h=0)

plot(theta, U(theta))
abline(h=0)



##------------- Using R inbuilt functions --------------

n <- 2423
ybar <- 3663 / n
Llik <- function(theta){
  n * ybar *log(theta) -n*theta -n*log(1 - exp(-theta))
 }

# optim function:
# optim(par, fn, method = c("Nelder-Mead", "BFGS", "CG", "L-BFGS-B", "SANN",
#                   "Brent"), control = list(), hessian = FALSE)

optim(ybar, Llik, method ="BFGS", control = list(fnscale=-1), 
      hessian = TRUE)
#0.8925626


# maxLik function:
# maxLik(logLik, grad = NULL, hess = NULL, start, method, constraints=NULL, ...)
#install.packages("maxLik")

library(maxLik)
maxLik(Llik, grad = NULL, hess = NULL, ybar, method="NR", constraints=NULL)
#0.8924961 

#Or try using the summary function for a nicer output:
MLE <- maxLik(Llik, grad = NULL, hess = NULL, ybar, method="NR", constraints=NULL)
summary(MLE)





