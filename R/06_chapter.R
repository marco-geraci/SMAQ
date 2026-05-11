# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 6
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
## Section 6.3
###################################################

## R Code 6.1

library(lqmm)
set.seed(123)
y <- rnorm(500)
fit <- mleAL(y)
fit
# quantile 0.4698466 of fitted AL
qal(fit$tau, mu = fit$m, sigma = fit$sigma, tau = fit$tau)
# sample quantile 0.4698466
quantile(y, probs = fit$tau, names = FALSE)
# quantile 0.75 of fitted normal
qnorm(0.75, mean = mean(y), sd = sd(y))
# sample quantile 0.75
quantile(y, probs = 0.75, names = FALSE)
# quantile 0.75 of fitted AL
qal(0.75, mu = fit$m, sigma = fit$sigma, tau = fit$tau)

## R Code 6.2

library(bbmle)
# maximise AL likelihood
loglik_al_constr <- function(mu, phi){
  sigma <- exp(phi)
  -sum(dal(x = y, mu = mu, sigma = sigma, tau = 0.75,
  log = TRUE))
}

fit <- mle2(loglik_al_constr, start = list(mu = 0, phi = 0))
thetahat <- as.numeric(coef(fit))
c(mu = thetahat[1], sigma = exp(thetahat[2]))

# similar to lqm
lqm(y ~ 1, tau = 0.75)

## R Code 6.3

loglik_al_constr <- function(mu, phi){
sigma <- exp(phi)
-sum(dal(x = y, mu = mu, sigma = sigma, tau = tau0,
log = TRUE))}

R <- 500
n <- 1000
mu <- 0
sigma <- 1
tau0 <- 0.75

set.seed(123)
nn <- list(NULL, c("Estimate", "SE"), c("Constr. AL",
"Sample quantile"))
res <- array(NA, dim = c(R, 2, 2), dimnames = nn)
CI <- array(NA, dim = c(R, 2, 2), dimnames = nn)

for(i in 1:R){
    y <- rnorm(n)
    fit <- mle2(loglik_al_constr, start = list(mu = 0, phi = 0))
    res[i,1,1] <- m <- coef(fit)[1]
    res[i,1,2] <- xip <- quantile(y, probs = tau0)
    s <- exp(coef(fit)[2])
    d <- density(y)
    res[i,2,1] <- se1 <- s/sqrt(n*tau0*(1-tau0))
    res[i,2,2] <- se2 <- sqrt(tau0*(1 - tau0))/
    (sqrt(n)*approxfun(d$x, d$y)(xip))
    CI[i,,1] <- m + c(-1,1)*qnorm(.975, 0, se1)
    CI[i,,2] <- xip + c(-1,1)*qnorm(.975, 0, se2)
}

# true population quantile
qnorm(tau0)
# average estimate and SE
apply(res, 3, colMeans)
# standard deviation of the estimates
apply(res[,1,], 2, sd)
# 95% confidence intervals
apply(CI - qnorm(tau0), 3, function(x) mean(
apply(x, 1, prod) < 0))

## R Code 6.4

library(evd)
R <- 100
n <- 1e5
k <- 100
p <- seq(0.99, 0.99999, length = k)
p0 <- 0.985; thresh <- qt(p0, df = 3)
range(n*(1-p))
xi <- array(NA, dim = c(k,R,2), dimnames =
list(quantile = 1:k,replication = 1:R, estimator =
c("empirical", "GP")))

set.seed(123)
for(i in 1:R){
    y <- rt(n, df = 3)
    Fn <- ecdf(y)
    p0i <- Fn(thresh)
    if(p[1] <= p0i) stop("Decrease threshold")
    fit <- fpot(y, threshold = thresh, model = "gpd")
    sigma <- fit$estimate['scale']
    gamma <- fit$estimate['shape']
    xi[,i,1] <- quantile(y, probs = p)
    xi[,i,2] <- thresh + sigma/gamma*(((1-p0i)/(1-p))^gamma - 1)
}

# calculate MSE
z <- apply(xi, 3, function(x, a) apply(sweep(x, 1, a,
FUN = "-")^2, 1, mean), a = qt(p, df = 3))

head(z)

###################################################
## Section 6.4
###################################################

## R Code 6.5

library(fitdistrplus)
library(EnvStats)

# select data
data("EPA.94b.tccb.df", package = "EnvStats")
dd <- subset(EPA.94b.tccb.df, EPA.94b.tccb.df$Area == "Reference")
y <- dd$TcCB
length(y)
# fit 3-parameter lognormal with quantile matching (QM)
pj <- c(0.05,.5,.95)
ss <- list(meanlog = mean(log(y)), sdlog = sd(log(y)), threshold = 0)
theta_qm <- qmedist(data = y, distr = "lnorm3",
probs = pj, start = ss)$estimate

# fit 3-parameter lognormal with maximum likelihood (ML)
theta_ml <- mledist(data = y, distr = "lnorm3",
start = ss)$estimate

# estimate quantiles over grid of p
p <- 1:49/50
Q_qm <- do.call(qlnorm3, c(list(p = p), as.list(theta_qm)))
Q_ml <- do.call(qlnorm3, c(list(p = p), as.list(theta_ml)))

# compare estimated quantiles with sample quantiles
xihat <- quantile(y, p)
sum(abs(xihat - Q_qm))
sum(abs(xihat - Q_ml))

