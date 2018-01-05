#
# R-code of the paper
#   van de Kassteele J, Hoogenveen RT, Engelfriet PM, van Baal PHM, Boshuizen HC (2012).
#   Estimating net transition probabilities from cross-sectional data with application to risk factors in chronic disease modeling.
#   Statistics in Medicine, 31(6), 533-543.
#
# Author: Jan van de Kassteele, RIVM
#

#
# Init
#

# Load packages
# (Install lpSolve package first)
library(mgcv)
library(Matrix)
library(lpSolve)

#
# Read data
#

# Read BMI data of Dutch males aged 0-85 in 2006/2007
bmi.data <- read.csv(file = "bmi data.csv")

# Show first records of bmi.data
# bmi.data has four variables: age and three BMI categories: normal, overweight, obese (counts)
head(bmi.data)

#
# Fit multinomial P-splines
#

# Number of classes
K <- 3

# Reshape bmi.data into person specific records (long format)
# We need this for the multinomial P-spline fit procedure
# 0 = normal, 1 = overweight, 2 = obese
bmi.data.long <- with(bmi.data, data.frame(
  age = rep(rep(age, times = K), times = c(normal, overweight, obese)),
  y = rep(0:(K - 1), times = c(sum(normal), sum(overweight), sum(obese)))))

# Fit multinomial P-splines
# See help(multinom) in the mgcv package
#
# Choose k large enough to allow enough flexibility, the penalty does the rest
# Because of the REML smoothness selection instead of BIC, the result is somewhat different from the fit in the paper
mod <- gam(
  formula = list(y ~ s(age, bs = "ps", k = 15), ~ s(age, bs = "ps", k = 15)),
  family = multinom(K = K - 1),
  data = bmi.data.long,
  method = "REML")

# Show summary of model fit
summary(mod)

#
# Predict multinomial P-splines
#

# Settings
age.pred <- 0:90
n.pred <- length(age.pred)

# Make prediction
# Returns a list of length two (fit and se.fit) containing an n.pred x K matrix
pred <- predict(mod,
  newdata = data.frame(age = age.pred),
  type = "response",
  se.fit = TRUE)

# Rename fit -> E.pi and se.fit -> se.pi
names(pred) <- c("E.pi", "se.pi")

# Add 95% CI
pred <- within(pred, {
  l.pi <- E.pi - qnorm(p = 0.975)*se.pi
  u.pi <- E.pi + qnorm(p = 0.975)*se.pi
})

# This is what we have created
str(pred)

#
# Optional: add simulations to pred object. Needed for 95% CI of transitions
#

# Number of simulations
n.sim <- 1000

# Extract spline coefficients and covariance matrix of spline coefficients
beta <- coef(mod)
V.beta <- vcov(mod)
n.beta <- length(beta)/(K-1) # Length of coefficient vector for each category

# Generate n.sim realisations for the coefficients
beta.sim <- rmvn(n = n.sim, mu = coef(mod), V = vcov(mod))

# Create model matrix for age.pred based on mod object
X <- predict(mod, newdata = data.frame(age = age.pred), type = "lpmatrix")

# Create empty n.pred x K x n.sim array in pred
pred$pi.sim <- array(0, dim = c(n.pred, K, n.sim))

# For each realisation of beta
for (i in 1:n.sim) {
  # Put beta.sim into block structure for easier computation of linear predictor
  beta.sim.block <- as.matrix(bdiag(split(x = beta.sim[i, ], f = rep(1:(K - 1), each = n.beta))))
  # Compute linear predictor eta. Zero is reference category
  eta <- cbind(0, X %*% beta.sim.block)
  # Compute pi and put it in pred object
  pred$pi.sim[, , i] <- exp(eta)/rowSums(exp(eta))
}

#
# Calculate net transitions
#

# Function for solving transportation problems using the simplex algorithm
transport.simplex <- function(supply, demand, cost, add.const.mat = NULL, add.const.dir = NULL, add.const.rhs = NULL) {
  # Input:
  # supply = K x 1 vector with supplies
  # demand = K x 1 vector with demands
  # cost   = K x K matrix with shipping costs
  # add.const.mat = matrix with additional constrains
  # add.const.dir = vector with additional directions
  # add.const.rhs = vector with additional right hand sides

  # Output:
  # trans = K x K matrix with transitions

  # Simplex
  sol <- lp(
    objective.in = as.vector(cost),
    const.mat = rbind(t(kronecker(diag(K), rep(1, K))), t(kronecker(rep(1, K), diag(K))), add.const.mat),
    const.dir = c(rep("=", 2*K), add.const.dir),
    const.rhs = c(supply, demand, add.const.rhs))
  trans <- matrix(sol$solution, K, K, byrow = TRUE)

  # Return output
  return(trans)
}

# Cost matrix
cost <- toeplitz(diffinv(1:(K-1)))

# Are simulations available?
exists.sim <- with(pred, exists("pi.sim"))

# Allocate K x K x n.pred-1 arrays in trans and trans.prob list objects for output
trans      <- list(E.trans      = array(NA, dim = c(K, K, n.pred - 1)))
trans.prob <- list(E.trans.prob = array(NA, dim = c(K, K, n.pred - 1)))

# If simulations are available, allocate K x K x n.pred-1 x n.sim arrays for them as well
if (exists.sim) {
  trans$trans.sim           <- array(NA, dim = c(K, K, n.pred - 1, n.sim))
  trans.prob$trans.prob.sim <- array(NA, dim = c(K, K, n.pred - 1, n.sim))
}

# Estimate transitions and transition probabilities
for (i in 1:(n.pred - 1)) {
  # Calculate net transitions and transition probabilities
  trans$E.trans[, , i] <- transport.simplex(
    supply = pred$E.pi[i, ],
    demand = pred$E.pi[i + 1, ],
    cost = cost)
  trans.prob$E.trans.prob[, , i] <- t(trans$E.trans[, , i]/pred$E.pi[i, ])
  # If simulations are available, do this n.sim times (takes a while if n.sim is large)
  if (exists.sim) {
    cat("age.pred =", i, "\n")
    for (j in 1:n.sim) {
      trans$trans.sim[, , i, j] <- transport.simplex(
        supply = pred$pi.sim[i, , j],
        demand = pred$pi.sim[i + 1, , j],
        cost = cost)
      trans.prob$trans.prob.sim[, , i, j] <- t(trans$trans.sim[, , i, j]/pred$pi.sim[i, , j])
    }
  }
}

# Calculate summary statistics, only if simulations are available
if (exists.sim) {
  trans <- within(trans, {
    E.trans <- apply(X = trans.sim, MARGIN = c(1, 2, 3), FUN = mean) # Replaces the previous E.trans
    l.trans <- apply(X = trans.sim, MARGIN = c(1, 2, 3), FUN = quantile, prob = 0.025)
    u.trans <- apply(X = trans.sim, MARGIN = c(1, 2, 3), FUN = quantile, prob = 0.975)
  })
  trans.prob <- within(trans.prob, {
    E.trans.prob <- apply(X = trans.prob.sim, MARGIN = c(1, 2, 3), FUN = mean) # Replaces the previous E.trans.prob
    l.trans.prob <- apply(X = trans.prob.sim, MARGIN = c(1, 2, 3), FUN = quantile, prob = 0.025)
    u.trans.prob <- apply(X = trans.prob.sim, MARGIN = c(1, 2, 3), FUN = quantile, prob = 0.975)
  })
}

# This is what we have created
str(trans)
str(trans.prob)

#
# Figures
#

# Names states (for titles and legend)
names.states <- c("normal weight", "overweight", "obese")

#
# Plot smoothed prevalences
#

# Set graphical parameters
par(mar = c(4.5, 4.5, 0.5, 0.5))

# Make the plot
plot.new()
plot.window(
  xlim = range(age.pred),
  ylim = c(0, 1))
# Smooth prevalences
with(pred, {
  for (i in 1:K) polygon(
    x = c(age.pred, rev(age.pred)),
    y = c(l.pi[, i], rev(u.pi[, i])),
    col = grey(0.95), border = NA)
  matlines(x = age.pred, y = E.pi, lty = 1, col = 1)
  matlines(x = age.pred, y = l.pi, lty = 3, col = grey(0.5))
  matlines(x = age.pred, y = u.pi, lty = 3, col = grey(0.5))
})
# Add data (as proportions)
with(bmi.data, matpoints(
  x = age,
  y = prop.table(as.matrix(bmi.data[, -1]), margin = 1),
  cex = 0.7, col = 1))
# Axes, titles and legend
axis(side = 1)
axis(side = 2)
box()
title(
  xlab = "Age",
  ylab = "Prevalence")
legend("topright",
  legend = names.states,
  pch = as.character(1:K),
  bty = "n")

#
# Plot transitions in K x K matrix
#

# Set graphical parameters
par(
  mar = c(3, 3, 1.5, 0.1),
  mfrow = c(K, K),
  mgp = c(1.8, 0.7, 0))

# Make the plot
for (i in 1:K) {
  for (j in 1:K) {
    plot.new()
    with(trans, {
      # If simulations are not available
      if (!exists.sim) {
        # Only plot the expected values
        plot.window(
          xlim = range(age.pred),
          ylim = range(E.trans[i, j, ]))
        lines(x = age.pred[-1], y = E.trans[i, j, ], lty = 1, col = 1)
      } else {
        # Plot expected values and confidence intervals
        plot.window(
          xlim = range(age.pred),
          ylim = range(c(l.trans[i, j, ], u.trans[i, j, ])))
        polygon(
          x = c(age.pred[-n.pred], rev(age.pred[-1])),
          y = c(l.trans[i, j, ], rev(u.trans[i, j, ])),
          col = grey(0.95), border = NA)
        lines(x = age.pred[-1], y = l.trans[i, j, ], lty = 3, col = grey(0.5))
        lines(x = age.pred[-1], y = u.trans[i, j, ], lty = 3, col = grey(0.5))
        lines(x = age.pred[-1], y = E.trans[i, j, ], lty = 1, col = 1)
      }
      # Axes and titles
      axis(1)
      axis(2)
      box()
      title(
        main = paste(names.states[i], "to", names.states[j]),
        xlab = "Age",
        ylab = "Transition")
    })
  }
}

#
# Plot transition probabilities in K x K matrix
#

# Set graphical parameters
par(
  mar = c(3, 3, 1.5, 0.1),
  mfrow = c(K, K),
  mgp = c(1.8, 0.7, 0))

# Make the plot
for (j in 1:K) {
  for (i in 1:K) {
    plot.new()
    with(trans.prob, {
      # If simulations are not available
      if (!exists.sim) {
        # Only plot the expected values
        plot.window(
          xlim = range(age.pred),
          ylim = range(E.trans.prob[i, j, ]))
        lines(x = age.pred[-1], y = E.trans.prob[i, j, ], lty = 1, col = 1)
      } else {
        # Plot expected values and confidence intervals
        plot.window(
          xlim = range(age.pred),
          ylim = range(c(l.trans.prob[i, j, ], u.trans.prob[i, j, ])))
        polygon(
          x = c(age.pred[-n.pred], rev(age.pred[-1])),
          y = c(l.trans.prob[i, j, ], rev(u.trans.prob[i, j, ])),
          col = grey(0.95), border = NA)
        lines(x = age.pred[-1], y = l.trans.prob[i, j, ], lty = 3, col = grey(0.5))
        lines(x = age.pred[-1], y = u.trans.prob[i, j, ], lty = 3, col = grey(0.5))
        lines(x = age.pred[-1], y = E.trans.prob[i, j, ], lty = 1, col = 1)
      }
      # Axes and titles
      axis(side = 1)
      axis(side = 2)
      box()
      title(
        main = paste(names.states[j], "to", names.states[i]),
        cex.main = 0.95,
        xlab = "Age",
        ylab = "Transition probability")
    })
  }
}
