# Solution to Workshop 2 Question: Iterative MLE

# PART 1.
# Firstly reead in the values of y.
y <- c(1.09,-0.23,0.79,2.31,-0.81)

# Set theta0 to 0. You may want to vary this.
theta0 <- 0

# Create the function U which takes in parameter theta,
# and then calculates the value in the {} brackets.
# This is the score function.
U <- function(theta){2*sum((y-theta)/(1+(y-theta)^2))}

# Initialise the theta parameter to theta0
theta<- theta0

# Set the tolerance of the function to 10^{-6}.
eps <- 1/10^6

# Iterate through the function until the tolerance criteria is met.
# This is based off the fisher method of scoring which says:
# theta_{r+1} = theta_{r} + U(theta_{r})/I(theta_{r}).
# Here we use the result from Asymptotic variance part where
# we worked out that I(theta) = n/2= 5/2
while(abs(U(theta)) > eps) {theta <- theta + U(theta)*2/5 }

# Print out the value that theta becomes.
print(theta)

# PART 2.
# Make the function for the log-likelihood.
# This takes in theta and returns the value in {}.
loglik <- function(theta){-sum(log(1+(y-theta)^2))}

# Import the package maxLik. This takes care of the maximum likelihood
# for us, instead of having to hard code this. You may need to run the
# install command first before running the library.
install.packages("maxLik")
library(maxLik)

# Run the function on the loglikelihood.
# This takes in the function as the first argument,
# grad: Null means we use the numerical gradient.
# hess: Null means we use the conventional Hessian.
# method: "NR" is the Newton Raphson method.
# constraints: Null means we are performing unconstrained maximisation.
NRmax <- maxLik(loglik, grad = NULL, hess = NULL, -10, method="NR", constraints=NULL)

# Call the summary which gives us a bunch of useful information.
summary(NRmax)

# We can now create a plot should we wish that will help us
# visually verify this is the maximum of the system.

# Set some entry points from below the smallest y and above the largest.
# Have these points increase in increments of 0.01
x = seq(min(y)-1, max(y)+1, 0.01)

# Initilise a vector of length x with NAs (ready for the likelihoods).
loglik.x = rep(NA, length(x))

# For each value in the vector x, calculate the likelihood and
# then store it in the correct position in the loglik.x vector.
for (i in 1:length(x))
  loglik.x[i]= loglik(x[i])

# Set the optimal value as known from the NR MLE.
theta.hat =  0.6393

# Plot the values from the log-likelihood against the x values.
plot(x,  loglik.x, type = 'l', lwd=2, xlab = "theta", ylab="log likelihood")

# Add the optimal value point.
points(theta.hat, loglik(theta.hat), col=2, pch=15)

# Add a title to the plot.
title("Log likelihood and MLE (red point), Cauchy example")
