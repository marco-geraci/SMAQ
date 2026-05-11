# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 5
# Marco Geraci
#
# This code reproduces analyses, figures, and simulations from:
# "Statistical Methods and Applications Using Quantiles"
#
# Copyright (C) 2026 Marco Geraci
#
# Licensed under the Creative Commons Attribution-NonCommercial 4.0
# International License (CC BY-NC 4.0).
#
# You are free to use, modify, and redistribute this material for
# non-commercial purposes, provided appropriate credit is given.
#
# Commercial use is not permitted without explicit written permission.
#
# License:
# https://creativecommons.org/licenses/by-nc/4.0/

###################################################
## Section 5.1
###################################################

## R Code 5.1

library(microbenchmark)
set.seed(564)
# n = 1000
y <- rnorm(1e3)
microbenchmark('quick' = sort(y, method = "quick"),
				'shell' = sort(y, method = "shell"),
				'radix' = sort(y, method = "radix"))
# n = 100000
y <- rnorm(1e5)
microbenchmark('quick' = sort(y, method = "quick"),
				'shell' = sort(y, method = "shell"),
				'radix' = sort(y, method = "radix"))

###################################################
## Section 5.2
###################################################

## R Code 5.2

library(lpSolve)
# generate data
set.seed(564)
n <- 5
print(y <- rnorm(n))
# define the LP problem
p <- 0.75
obj <- c(rep(p, n), rep(1 - p, n), 0)
S <- cbind(diag(2*n), 0)
A <- cbind(diag(n), -diag(n), 1)
b <- c(rep(0, 2*n), y)
constr <- c(rep(">=", 2*n), rep("=", n))
# solve
fit.lp <- lp("min", obj, rbind(S, A), constr, b)
fit.lp
# compare solution with standard R function
print(xi <- fit.lp$solution[2*n+1])
quantile(y, probs = p, names = FALSE)
# check linear constraints xi + u - v = y
u <- fit.lp$solution[1:n]
v <- fit.lp$solution[(n+1):(2*n)]
all.equal(xi + u - v, y)
# calculate loss function 'psi'
print(r <- u - v)
sum(r * (p - ifelse(r < 0, 1, 0)))

## R Code 5.3

constr <- c(rep(">=", n), rep("<=", n), rep("=", 1))
fitd.lp <- lp("max", y, rbind(diag(n),diag(n),
t(matrix(1, nrow = n))), constr, c(rep(0, n), rep(1, n), n*(1-p)))
fitd.lp
a <- fitd.lp$solution
d <- round(a - (1 - p), 3)
sum(y*d)
cbind(y, u, v, a, d, rank = rank(y))

###################################################
## Section 5.3
###################################################

## R Code 5.4

library(lqmm)
print(y)
# optimise - set verbose = TRUE to print out tracing information
p <- 0.75
fit.lqm <- lqm(y ~ 1, tau = p, control = lqmControl(verbose = FALSE))
coef(fit.lqm)

###################################################
## Section 5.4
###################################################

## R Code 5.5

library(lqmm)
print(y)
p <- 0.75
# define decreasing sequence of tuning constants
omicron <- c(sort(seq(0.1, 1, by = 0.1), decreasing = TRUE),
	0.05, 0.001, 0.0001)
# optimise for each value of the sequence
K <- length(omicron)
mu <- g <- rep(NA, K)
for(k in 1:K){
    fit.lqm <- lqm(y ~ 1, tau = p, control = lqmControl(method =
    "gs2", smooth = TRUE, omicron = omicron[k]))
    mu[k] <- coef(fit.lqm)
    # obtain quantile loss from asymmetric Laplace log-likelihood
    g[k] <- -(fit.lqm$logLik - n*log(p*(1 - p)/fit.lqm$scale))*
    fit.lqm$scale
}
c(estimate = mu[K], objective = g[K])

## R Code 5.6

print(y)
# IRLS algorithm
p <- 0.75
eps <- .Machine$double.eps^(2/3)
STOP <- FALSE; mu_old <- mean(y); k <- 1
while(!STOP & k < 100){
    # uncomment the line below to print out tracing information
    # cat("iteration", k, "; mu = ", mu_old, "\n")
    w <- 1/abs(y - mu_old + eps)
    mu_new <- (n*(2*p - 1) + sum(w*y))/sum(w)
    if(abs(mu_new - mu_old) < 1e-8)	STOP <- TRUE
    	else {mu_old <- mu_new; k <- k + 1}
}
c('mu' = mu_new, 'xi' = quantile(y, probs = p, names = FALSE))

