# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 13
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
## Section 13.1
###################################################

## R Code 13.1

library(lqmm)
set.seed(123)

M <- 4000
n <- 20
mu <- 0
omega_u <- 1
omega_e <- 1
p <- 0.75
id <- rep(1:M, each = n)

u <- rnorm(M, 0, omega_u)
e <- rnorm(M*n, 0, omega_e)
y <- mu + rep(u, each = n) + e

# true population-level quantile
mu + sqrt(omega_e^2 + omega_u^2)*qnorm(p)
# marginal sample quantile
quantile(y, probs = p, names = FALSE)
# true average of cluster-specific quantiles
# under symmetric random effects distribution
mu + omega_e*qnorm(p)
# average cluster-specific sample quantiles
print(xibar <- mean(tapply(y, id, quantile, probs = p)))
# lqmm estimate
fit <- lqmm(y ~ 1, random = ~1, group = id, tau = p,
data = data.frame(y = y, id = id), nK = 11)
coef(fit)

## R Code 13.2

library(nlme)
library(lqmm)

load(file = "../data/corals.Rdata")
# summary stats for staghorn corals dataset
summary(corals)
# calculate ICC and effective sample size with lme
fit.lme <- lme(Depth ~ 1, random = ~ 1|Species, data = corals)
# get random intercept variance, residual variance
omega <- as.numeric(VarCorr(fit.lme)[,1])
# ICC = sigma_u^2 / (sigma_u^2 + sigma_e^2)
print(ICC <- omega[1]/sum(omega))
# average cluster size
print(nbar <- mean(table(corals$Species)))
# total sample size
print(N <- nrow(corals))
# design effect
DEFF <- 1 + (nbar - 1)*ICC
# effective sample size
print(Ne <- N/DEFF)
# calculate ICC with lqmm at p = 0.9
ctrl <- lqmmControl(method = "df", LP_max_iter = 1000, LP_tol_ll = 1e-5)
fit.lqmm <- lqmm(fixed  = Depth ~ 1, random = ~ 1, group  = Species,
tau = 0.9, nK = 7, type = "normal", data = corals, control = ctrl)
print(psi_u <- VarCorr(fit.lqmm))
print(psi_e <- varAL(fit.lqmm$scale, tau = 0.9))
psi_u/(psi_u + psi_e)

## R Code 13.3

library(lqmm)

# set control parameters for estimation
ctrl <- lqmmControl(method = "df", LP_max_iter = 1000, LP_tol_ll = 1e-5)
# 2nd degree orthogonal polynomial on latitude
P <- poly(corals$Latitude, degree = 2)
corals$lat1 <- P[,1]; corals$lat2 <- P[,2]
# fit LQMM
fit <- lqmm(fixed  = Depth ~ lat1 + lat2, random = ~ 1,
group  = Species, tau = 0.9, nK = 7, type = "normal", data = corals,
control = ctrl)
fit

## R Code 13.4

print(u <- ranef(fit))
# fit a QR model for each species
xx <- by(corals, corals$Species, function(x)
  coef(rq(Depth ~ poly(Latitude, 2), data = x, tau = 0.9)))
xx <- do.call("rbind", xx)
# calculate deviation of each species intercept from the mean intercept
delta <- xx[,1] - mean(xx[,1])
# correlation between RQ deviations and LQMM predicted random effects
cor(delta, u[[1]])

## R Code 13.5

load(file = "../data/internet.Rdata")
summary(internet)

## R Code 13.6

# define z
internet$z <- 1*(internet$incomegr == "High-income countries")

# first approach

# z = 0
dd0 <- subset(internet, internet$z == 0)
fit0 <- nls(percentage ~ SSlogis(year, Asym, xmid, scal), data = dd0)

# z = 1
dd1 <- subset(internet, internet$z == 1)
fit1 <- nls(percentage ~ SSlogis(year, Asym, xmid, scal), data = dd1)

ss0 <- as.list(coef(fit0))
ss1 <- as.list(coef(fit1))
ss <- rep(NA, 6)
ss[c(1,3,5)] <- as.numeric(ss0)
ss[c(2,4,6)] <- as.numeric(ss1) - as.numeric(ss0)
names(ss) <- paste0("b", 1:6)
round(ss, 3)
ss_for_replication <- ss

# second approach

fit <- nls(percentage ~ SSlogis(year,
  Asym0 + Asym1*z,
  xmid0 + xmid1*z,
  scal0 + scal1*z), data = internet,
  start = list(Asym0 = 80, Asym1 = 10, xmid0 = 2005, xmid1 = 2,
  scal0 = 1, scal1 = 1))

ss <- coef(fit)
names(ss) <- paste0("b", 1:6)
round(ss, 3)

## R Code 13.7

ss <- ss_for_replication

#pak::pak("marco-geraci/nlqmm")
library(nlqmm)

# fit nlqmm for the median
nlqmm(model = percentage ~ Asym/(1 + exp((xmid - year)/scal)),
  fixed = Asym + xmid + scal ~ z, random = Asym + xmid ~ 1,
  group = iso3, covariance = "pdSymm", tau = 0.5, data = internet,
  start = as.list(ss))

###################################################
## Section 13.2
###################################################

## R Code 13.8

set.seed(123)
n <- 500
u <- runif(n)
y <- rep(0, n)

a0 <- function(u) qnorm(u)
a1 <- function(u) 0.3 + 0.6*u

for(i in 2:n){
  y[i] <- a0(u[i]) + a1(u[i]) * y[i-1]
}

dat <- data.frame(y = y[-1], ylag = y[-n])

p <- seq(0.1, 0.9, by = 0.1)
fit <- rq(y ~ ylag, tau = p, data = dat)

## Figure 13.8

p <- seq(0.1, 0.9, by = 0.1)
fit <- rq(y ~ ylag, tau = p, data = dat)
cf <- as.data.frame(t(coef(fit)))
cf$p <- p

p1 <- ggplot(data.frame(t = 1:n, y = y), aes(t, y)) +
  geom_line() +
  labs(x = "Time", y = expression(y[t]))

cf_long <- rbind(
  data.frame(p = cf$p, value = cf$`(Intercept)`, term = "Intercept"),
  data.frame(p = cf$p, value = cf$ylag, term = "Slope")
)

p2 <- ggplot(cf_long, aes(p, value)) +
  geom_line() +
  geom_point() +
  facet_wrap(~term, scales = "free_y") +
  labs(x = "p", y = "Estimate")

gridExtra::grid.arrange(p1, p2, ncol = 1)

## R Code 13.9

library(quantreg)
library(ggplot2)

data(UKgas)

y <- as.numeric(UKgas)
dat <- data.frame(
  y = y[-1],
  ylag = y[-length(y)]
)

# fit QAR(1)
fit <- rq(y ~ ylag, tau = c(0.1, 0.5, 0.9), data = dat)

# obtain standard errors
summary(fit, se = "nid")
