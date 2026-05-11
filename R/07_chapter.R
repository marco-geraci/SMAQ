# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 7
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
## Section 7.1
###################################################

## R Code 7.1

library(numDeriv)
library(bbmle)
# define exponential log-likelihood and quantiles as functions of theta
loglik_exp <- function(theta){
    -sum(dexp(y, rate = theta, log = TRUE))
}
fn <- function(theta, p){
    qexp(p = p, rate = theta)
}
# generate data
set.seed(123)
n <- 100
p <- 0.9
theta <- 2
y <- rexp(n, rate = theta)
# MLE of theta
fit <- mle2(loglik_exp, start = list(theta = 1/mean(y)))
# estimate quantile, SE and confidence interval
xi <- qexp(p = p, rate = coef(fit))
vtheta <- grad(fn, x = coef(fit), p = p)
xi_se <- sqrt(vtheta^2 * vcov(fit))
CI <- xi + c(-1,1)*qnorm(1-0.05/2, 0, xi_se)
c(xi = xi, se = xi_se, low = CI[1], up = CI[2])
# compare numerical and analytical calculation of standard error
abs(log(1-p))/(coef(fit)*sqrt(n))

###################################
# Reparameterise to ensure
# parameter estimates are positive
###################################

loglik_exp <- function(lograte){
    -sum(dexp(y, rate = exp(lograte), log = TRUE))
}
fn <- function(lograte, p){
    qexp(p = p, rate = exp(lograte))
}
fit <- mle2(loglik_exp, start = list(lograte = log(1/mean(y))))
xi <- qexp(p = p, rate = exp(coef(fit)))
vlograte <- grad(fn, x = coef(fit), p = p)
xi_se <- sqrt(vlograte^2 * vcov(fit))
CI <- xi + c(-1,1)*qnorm(1-0.05/2, 0, xi_se)
c(xi = xi, se = xi_se, low = CI[1], up = CI[2])

## R Code 7.2

# define gamma log-likelihood and quantiles as functions of theta
loglik_gamma <- function(shape, rate){
    -sum(dgamma(y, shape = shape, rate = rate, log = TRUE))
}
fn <- function(theta, p){
    qgamma(p = p, shape = theta[1], rate = theta[2])
}
# generate data
set.seed(123)
n <- 100
p <- 0.9
y <- rgamma(n, shape = 2, rate = 2)
ss <- list(shape = mean(y)^2/var(y), rate = mean(y)/var(y))
# MLE of theta
fit <- mle2(loglik_gamma, start = ss)
# estimate quantile, SE and confidence interval
xi <- qgamma(p = p, shape = coef(fit)[1], rate = coef(fit)[2])
vtheta <- matrix(grad(fn, x = coef(fit), p = p))
xi_se <- sqrt(as.numeric(crossprod(vtheta,
vcov(fit)) %*% vtheta))
CI <- xi + c(-1,1)*qnorm(1-0.05/2, 0, xi_se)
c(xi = xi, se = xi_se, low = CI[1], up = CI[2])

###################################
# Reparameterise to ensure
# parameter estimates are positive
###################################

loglik_gamma <- function(logshape, lograte){
    -sum(dgamma(y, shape = exp(logshape),
    rate = exp(lograte), log = TRUE))
}
fn <- function(theta, p){
    qgamma(p = p, shape = exp(theta[1]),
    rate = exp(theta[2]))
}
ss <- list(logshape = log(mean(y)^2/var(y)), lograte = log(mean(y)/var(y)))
fit <- mle2(loglik_gamma, start = ss)
xi <- qgamma(p = p, shape = exp(coef(fit)[1]),
rate = exp(coef(fit)[2]))
vlogtheta <- matrix(grad(fn, x = coef(fit), p = p))
xi_se <- sqrt(as.numeric(crossprod(vlogtheta,
vcov(fit)) %*% vlogtheta))
CI <- xi + c(-1,1)*qnorm(1-0.05/2, 0, xi_se)
c(xi = xi, se = xi_se, low = CI[1], up = CI[2])

## R Code 7.3

set.seed(123)
n <- 100
y <- rnorm(n)

# optimal smoothing with SJ
fhat1 <- density(y, bw = "SJ", kernel = "gaussian")
x1 <- fhat1$x
y1 <- fhat1$y

# undersmoothing (one fifth of the SJ bandwidth)
fhat2 <- density(y, bw = fhat1$bw/5, kernel = "gaussian")
x2 <- fhat2$x
y2 <- fhat2$y

# true density
x3 <- seq(-3, 3, by = 0.01)
y3 <- dnorm(x3)

# plot results
gg <- ggplot() + geom_line(aes(x = x1, y = y1), linetype = "dashed") +
geom_line(aes(x = x2, y = y2), linetype = "dotted") +
geom_line(aes(x = x3, y = y3), linetype = "solid") + xlab("y") +
ylab("Density")

gg

## R Code 7.4

library(Qtools)
set.seed(123)
n <- 100
y <- rbeta(n, shape1 = 0.5, shape2 = 0.5)
z <- logit(y)
r <- c(0.01, 1-0.01)

# optimal smoothing with SJ
fhat1 <- density(z, bw = "SJ", kernel = "gaussian",
from = logit(r[1]), to = logit(r[2]))
x1 <- invlogit(fhat1$x)
y1 <- fhat1$y/(x1*(1-x1))

# oversmoothing (20 times the SJ bandwidth)
fhat2 <- density(z, bw = fhat1$bw*20, kernel = "gaussian",
from = logit(r[1]), to = logit(r[2]))
x2 <- invlogit(fhat2$x)
y2 <- fhat2$y/(x2*(1-x2))

# true density
x3 <- seq(r[1], r[2], by = 0.01)
y3 <- dbeta(x3, shape1 = 0.5, shape2 = 0.5)

# plot results
gg <- ggplot() + geom_line(aes(x = x1, y = y1), linetype = "dashed") +
geom_line(aes(x = x2, y = y2), linetype = "dotted") +
geom_line(aes(x = x3, y = y3), linetype = "solid") + xlab("y") +
ylab("Density")

gg

## R Code 7.5

set.seed(123)
n <- 100
y <- rnorm(n)
p <- 0.5
alpha <- 0.05
# density estimation
fhat <- density(y, bw = "SJ", kernel = "gaussian")
# sample quantile
xi <- quantile(y, probs = p)
# density at quantile
f0 <- approx(fhat$x, fhat$y, xout = xi)$y
#standard error
se <- sqrt(p*(1-p)/n)/f0
# confidence interval
ci <- xi + outer(qnorm(1-alpha/2, 0, se), c(-1,1), FUN = "*")
# print results
c(xi = xi, low = ci[1], up = ci[2])

## R Code 7.5

set.seed(123)
n <- 100
y <- rbeta(n = n, shape1 = 0.5, shape2 = 0.5)
Qtools::qexact(x = y, probs = 0.1, level = 0.95)

## R Code 7.7

set.seed(123)
p <- 0.1
Y <- replicate(10000, rbeta(n = 100, shape1 = 0.5, shape2 = 0.5))
Q <- qbeta(p, shape1 = 0.5, shape2 = 0.5)
out <- t(apply(Y, 2, function(y, probs, level)
Qtools::qexact(y, probs, level), probs = p, level = 0.95))
# CIs that include true quantile Q have discordant signs of
# bounds centered around Q
I <- ((out[,2] - Q)*(out[,3] - Q)) < 0
# empirical coverage
mean(I)

## R Code 7.8

# sample from a population
y <- c(1.76, 5.98, 11.90, 12.66, 15.43)
n <- length(y)
print(xihat <- median(y))
# take B = 900 bootstrap samples with replacement from y
set.seed(321)
Y_star <- replicate(900, sample(y, n, replace = TRUE))

# calculate the median of each bootstrap sample
xi_star <- apply(Y_star, 2, median)

# first two samples
Y_star[,1:2]
# first two medians
xi_star[1:2]
# bootstrap estimate of the sampling distribution of the sample median
z <- sqrt(n)*(xi_star - xihat)
Gstar <- ecdf(z)

# 95% confidence interval for the population median (ordinary)
c(xihat - quantile(z, probs = 0.975, names = FALSE)/sqrt(n),
xihat - quantile(z, probs = 0.025, names = FALSE)/sqrt(n))
# 95% confidence interval for the population median (percentile method)
c(xihat + quantile(z, probs = 0.025, names = FALSE)/sqrt(n),
xihat + quantile(z, probs = 0.975, names = FALSE)/sqrt(n))
# 95% confidence interval for the population median (normal approximation)
c(xihat - qnorm(0.975)*sd(xi_star), xihat + qnorm(0.975)*sd(xi_star))

## R Code 7.9

library(Qtools)
library(discreteQ)

data("fars", package = "Qtools")
y <- fars$count
# sample size
length(y)
# summary statistics
summary(y)
# quantile levels
p <- c(0.25, .5, .75, .9, .95, .99)
# estimate quantiles and 95% CI based on 1000 bootstrap replications
set.seed(321)
system.time(Qhat <- discreteQ(y, bsrep = 1000, alpha = 0.05))
cbind(est = Qhat$lb.Q(p), lower = Qhat$Q(p), upper = Qhat$ub.Q(p))
# estimate mid-quantiles and 95% CIs
system.time({midQ <- midquantile(y, probs = p);
  midQci <- confint(midQ, level = 0.95)})
midQci
# recover quantiles from mid-quantiles
as.numeric(midq2q(midQ))

## R Code 7.10

library(reshape2)
library(fitdistrplus)
library(EnvStats)
data(EPA.94b.tccb.df)
dd <- subset(EPA.94b.tccb.df, EPA.94b.tccb.df$Area == "Reference")
y <- dd$TcCB
n <- length(y)

# fit lognormal with quantile matching
ss <- list(meanlog = mean(log(y)), sdlog = sd(log(y)),
threshold = min(y)-0.1)
theta_qm <- qmedist(data = y, distr = "lnorm3",
probs = c(0.05,.5,.95), start = ss)$estimate

# Q-Q plot ordinates
z <- sort(y)

# probabilities (shifted) corresponding to z
p <- (sort(rank(y, ties.method = "average")) - 0.5)/n

# Q-Q plot abscissas
xit <- do.call(qlnorm3, c(list(p = p), as.list(theta_qm)))

# reference line
que <- quantile(y, c(0.25, 0.75), names = FALSE)
qut <- do.call(qlnorm3, c(list(p = c(0.25, 0.75)),
as.list(theta_qm)))
slope <- diff(que)/diff(qut)
int <- que[1] - slope * qut[1]

# 95% confidence interval (asymptotic)
f <- do.call(dlnorm3, c(list(x = xit), as.list(theta_qm)))
se <- slope*sqrt(p*(1-p)/n)/f
zhat <- int + slope * xit
ci <- melt(zhat + outer(qnorm(0.975, 0, se), c(-1,1), FUN = "*"))
ci$xit <- rep(xit, 2)

# plot results
gg <- ggplot() + geom_point(aes(x = xit, y = z)) +
geom_line(data = ci, aes(x = xit, y = value, group = Var2),
linetype = "dashed") + geom_abline(intercept = int, slope = slope) +
xlab("Theoretical") + ylab("Empirical") + xlim(0, 2.45) + ylim(0, 2.45)

gg

###################################################
## Section 7.2
###################################################

## R Code 7.11

data(PlantGrowth)
y <- PlantGrowth$weight[PlantGrowth$group == "ctrl"]
n <- length(y)
p <- 0.5
print(xihat <- quantile(y, probs = p, names = FALSE))
xi0 <- 5

# Wald test
fhat <- density(y, bw = "SJ", kernel = "gaussian")
f <- approx(fhat$x, fhat$y, xout = xihat)$y
sigma <- sqrt(p*(1-p))/(sqrt(n)*f)
Tw <- (xihat - xi0)/sigma
2*(1 - pnorm(abs(Tw)))
# exact binomial
Tb <- sum(y <= xi0)
2*min(pbinom(Tb, n, p), 1 - pbinom(Tb, n, p) + dbinom(Tb, n, p))

## R Code 7.12

data(PlantGrowth)
y1 <- PlantGrowth$weight[PlantGrowth$group == "ctrl"]
y2 <- PlantGrowth$weight[PlantGrowth$group == "trt1"]
n1 <- length(y1)
n2 <- length(y2)
print(d <- median(y2) - median(y1))
# take B = 900 bootstrap samples with replacement from y1 and y2
set.seed(111)
B <- 900
Y1_star <- replicate(B, sample(y1, n1, replace = TRUE))
Y2_star <- replicate(B, sample(y2, n2, replace = TRUE))

# calculate the median of each bootstrap sample for y1 and y2
xi1_star <- apply(Y1_star, 2, median)
xi2_star <- apply(Y2_star, 2, median)

# take the difference
d_star <- xi2_star - xi1_star

# calculate p-value
Gstar <- (sum(d_star < 0) + 0.5*sum(d_star == 0))/B
print(pv <- 2*min(Gstar, 1 - Gstar))
# t-test
t.test(y1, y2, var.equal = FALSE)

## R Code 7.13

data(PlantGrowth)

p <- 0.5
y <- split(PlantGrowth$weight, PlantGrowth$group)
k <- length(y) - 1
print(xi <- sapply(y, quantile, probs = p, names = FALSE, type = 1))
# estimate quantile differences (treatments vs control)
print(delta <- kronecker(xi[-1], xi[1], "-"))
# estimate group-specific densities
fj <- sapply(y, function(z, p)
{fhat <- density(z, bw = "SJ", kernel = "gaussian")
return(approx(fhat$x, fhat$y, xout = quantile(z, probs = p))$y)}, p = p)

# calculate S
nj <- table(PlantGrowth$group)
n <- sum(nj)
mj <- nj*fj^2
S <- matrix(1, k, k)
diag(S) <- kronecker(mj[1], mj[-1], "+")/mj[-1]
S <- S*p*(1-p)/mj[1]

# test on equality of medians of the three groups
print(Tf <- as.numeric(t(delta)%*%solve(S)%*%delta)/k)
print(pv <- 1 - pf(Tf, df1 = k, df2 = n - (k + 1)))
# F-test (analysis of variance)
summary(aov(weight ~ group, PlantGrowth))

## R Code 7.14

data(PlantGrowth)
y1 <- PlantGrowth$weight[PlantGrowth$group == "ctrl"]
y2 <- PlantGrowth$weight[PlantGrowth$group == "trt1"]
n1 <- length(y1)
n2 <- length(y2)
n <- n1 + n2

p <- c(0.25, 0.75)
m <- length(p)
xi1 <- quantile(y1, probs = p, names = FALSE)
xi2 <- quantile(y2, probs = p, names = FALSE)

# estimate quantile differences (treatment vs control)
print(delta <- xi2 - xi1)
# estimate group-specific densities
fj <- sapply(list(y1, y2), function(z, p)
{fhat <- density(z, bw = "SJ", kernel = "gaussian")
return(approx(fhat$x, fhat$y, xout = quantile(z, probs = p))$y)},p = p)
rownames(fj) <- p

# calculate V
pp <- outer(p, p, FUN = pmin) - outer(p, p, FUN = "*")
ff1 <- outer(fj[,1], fj[,1], FUN = "*")
ff2 <- outer(fj[,2], fj[,2], FUN = "*")
print(V <- pp/ff1 * (n1*ff1 + n2*ff2)/(n2*ff2))
# define C
C <- matrix(c(-1,1), nrow = 1)

# test on equality of QTEs
Cd <- C%*%delta
CVC <- C%*%V%*%t(C)
print(Tf <- as.numeric(t(Cd)%*%solve(CVC)%*%(Cd)))
print(pv <- 1 - pf(Tf, df1 = m - 1, df2 = n * m - (m - 1)))

## Figure 7.5

library(quantreg)
set.seed(123)
n <- 200

dd1 <- data.frame(x = sort(rnorm(n, 0, 1)), y = sort(rnorm(n, 0, 1)))
dd2 <- data.frame(x = sort(rnorm(n, 0, 1)), y = sort(rnorm(n, 2, 1)))
dd3 <- data.frame(x = sort(rnorm(n, 0, 1)), y = sort(rnorm(n, 2, 3)))
dd4 <- data.frame(x = sort(rnorm(n, 0, 1)), y = sort(scale(rexp(n, 0.2))))

p <- (1:n)/(n+1)
x <- rep(0:1, each = n)
z <- with(dd1, c(x, y))
p1 <- anova(quantreg::rq(z ~ x, tau = p))$table$pvalue
z <- with(dd2, c(x, y))
p2 <- anova(quantreg::rq(z ~ x, tau = p))$table$pvalue
z <- with(dd3, c(x, y))
p3 <- anova(quantreg::rq(z ~ x, tau = p))$table$pvalue
z <- with(dd4, c(x, y))
p4 <- anova(quantreg::rq(z ~ x, tau = p))$table$pvalue

p1 <- if(p1 < 0.001) "< 0.001" else paste0("= ", fround(p1,3))
p2 <- if(p2 < 0.001) "< 0.001" else paste0("= ", fround(p2,3))
p3 <- if(p3 < 0.001) "< 0.001" else paste0("= ", fround(p3,3))
p4 <- if(p4 < 0.001) "< 0.001" else paste0("= ", fround(p4,3))
lab <- "       p "

sz <- 10
lsz <- 1

gg <- list()


gg[[1]] <- ggplot(dd1, aes(x = x, y = y)) + geom_point() + xlab("") + ylab(expression(y[2])) + theme(text = element_text(size=sz)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + ggtitle(paste0("(a)", lab, p1))

gg[[2]] <- ggplot(dd2, aes(x = x, y = y)) + geom_point() + xlab("") + ylab("") + theme(text = element_text(size=sz)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + ggtitle(paste0("(b)", lab, p2))

gg[[3]] <- ggplot(dd3, aes(x = x, y = y)) + geom_point() + xlab(expression(y[1])) + ylab(expression(y[2])) + theme(text = element_text(size=sz)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + ggtitle(paste0("(c)", lab, p3))

gg[[4]] <- ggplot(dd4, aes(x = x, y = y)) + geom_point() + xlab(expression(y[1])) + ylab("") + theme(text = element_text(size=sz)) + geom_abline(intercept = 0, slope = 1, linetype = "dashed") + ylim(c(-3,4)) + ggtitle(paste0("(d)", lab, p4))


marrangeGrob(gg, ncol = 2, nrow = 2, layout_matrix = matrix(1:4, 2, 2, byrow = TRUE), top = "")
