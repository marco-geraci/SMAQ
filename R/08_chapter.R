# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 8
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
## Section 8.3
###################################################

## R Code 8.1

library(lpSolve)
library(quantreg)
# generate data
set.seed(123)
n <- 500
x <- runif(n, 0, 3)
b <- c(1,3)
g <- c(1,2)
y <- cbind(1,x)%*%b + rnorm(n, mean = 0, sd = cbind(1,x)%*%g)
p <- 0.9
# define the LP problem
obj <- c(rep(p, n), rep(1 - p, n), rep(0, 2))
S <- cbind(diag(2*n), 0, 0)
A <- cbind(diag(n), -diag(n), 1, x)
b <- c(rep(0, 2*n), y)
constr <- c(rep(">=", 2*n), rep("=", n))
# solve
fit.lp <- lp("min", obj, rbind(S, A), constr, b)
fit.lp
# compare solution with function rq from quantreg package
print(bhat <- fit.lp$solution[-c(1:(2*n))])
fit.rq <- rq(y ~ x, tau = p, method = "br")
coef(fit.rq)
# check linear constraints xi + u - v = y
u <- fit.lp$solution[1:n]
v <- fit.lp$solution[(n+1):(2*n)]
all.equal(cbind(1,x)%*%bhat + u - v, y)
# calculate residuals and count 0s
r <- u - v
sum(mapply(identical, x = r, y = 0))
# calculate loss function 'psi'
sum(r * (p - ifelse(r < 0, 1, 0)))
fit.rq$rho

## R Code 8.2

library(quantreg)
library(lqmm)
set.seed(123)
n <- 200e3
x <- runif(n, 0, 3)
b <- c(1,3)
g <- c(1,2)
y <- cbind(1,x)%*%b + rnorm(n, mean = 0, sd = cbind(1,x)%*%g)
p <- 0.9

system.time(bhat.br <- rq(y ~ x, tau = p, method = "br")$coef)
system.time(bhat.fn <- rq(y ~ x, tau = p, method = "fn")$coef)
system.time(bhat.lqm <- lqm(y ~ x, tau = p)$theta)
# compare estimates
cbind(bhat.br, bhat.fn, bhat.lqm)

## R Code 8.3

library(quantreg)
library(conquer)
set.seed(123)
n <- 500e3
k <- 20
x <- matrix(runif(n*k, 0, 3), n, k)
b <- rep(1,k+1)
g <- rep(1,k+1)
y <- cbind(1,x)%*%b + rnorm(n, mean = 0, sd = cbind(1,x)%*%g)
p <- 0.9

system.time(bhat.fn <- rq(y ~ x, tau = p, method = "fn")$coef)
system.time(bhat.cq <- conquer(x, y, tau = p)$coef)
# true beta
beta <- b + qnorm(p)*g
# mean absolute relative error (average)
mean(abs(bhat.fn - beta)/beta)
mean(abs(bhat.cq - beta)/beta)

## R Code 8.4

library(quantreg)
data(faithful)
cutoff <- 3
x <- as.numeric(faithful$eruptions >= cutoff)
y <- faithful$waiting
p <- c(0.1, 0.5, 0.9)
fit <- rq(y ~ x, tau = p, method = "br")
fit

## R Code 8.5

# calculate J and xbar for generic x
n <- length(y)
X <- cbind(1, x)
J <- crossprod(X)/n
xbar <- c(1,mean(x))
# alternatively with binary x
# u <- mean(x == 1); J <- matrix(c(1,u,u,u), 2, 2); xbar <- c(1,u)

# Hall-Sheather bandwidth (hs = TRUE)
print(h <- quantreg::bandwidth.rq(p = p, n = n, hs = TRUE, alpha = 0.05))
# obtain dQ from two quantile models
bhat_low <- rq(y ~ x, tau = p - h, method = "br")$coef
bhat_up <- rq(y ~ x, tau = p + h, method = "br")$coef
dQ <- t(xbar) %*% (bhat_up - bhat_low)
as.numeric(dQ)
# difference quotient
sparsity <- dQ/(2*h)
as.numeric(sparsity)
# calculate variance-covariance
V <- outer(solve(J), as.numeric(p*(1-p)*sparsity^2/n), "*")
dimnames(V) <- list(c("intercept", "slope"),
c("intercept", "slope"), p)
V
# standard errors
apply(V, 3, function(x) sqrt(diag(x)))

## Code fragment

e <- lm(y ~ x)$residuals
1/dnorm(qnorm(p, 0, sd(e)), 0, sd(e))


## R Code 8.6

fit.s <- summary(fit, se = "iid", covariance = TRUE)
# extract variance-covariance matrix for each quantile
lapply(fit.s, function(x) x$cov)
# calculate 95% CI for p = 0.9
betahat <- coef(fit)[,3]
se <- sqrt(diag(fit.s[[3]]$cov))
betahat
cbind(betahat - qnorm(0.975)*se, betahat + qnorm(0.975)*se)

## R Code 8.7

load(file = "../data/fev.Rdata")
summary(fev)

## R Code 8.8

fit <- rq(fev ~ I(age - 10) + I(height - 61) + sex,
    tau = c(0.1, 1:3/4, 0.9), data = fev)
fit
anova(fit, joint = TRUE)
anova(fit, joint = FALSE)

## R Code 8.9

# NOT AVAILABLE DUE TO PROPRIETARY DATA

## R Code 8.10

# singletons only
dd <- subset(USbirths2001, USbirths2001$multiple == 0)
# some missing ga_group and tox
dd <- na.omit(dd[,c("bw", "ga_group", "female", "tox", "died", "cause")])
# about 3.3 million
nrow(dd)
# quantile regression at p = 0.1 (note method Frisch-Newton for large size)
fit.rq <- rq(bw ~ ga_group + female, data = dd, method = "fn", tau = 0.1)
summary(fit.rq, se = "nid")
# define sga as indicator of negative residuals
dd$sga <- as.numeric(residuals(fit.rq) < 0)
# logistic regression for mortality
fit.glm <- glm(died ~ ga_group + female + sga + tox,
   family = "binomial", data = dd)
round(exp(cbind(OR = coef(fit.glm), confint.default(fit.glm))), 3)

## R Code 8.11

fit2 <- rq(bw ~ cause, tau = c(0.05, 0.33, 0.66, 0.95), data = dd, method = "fn")
fit2
anova(fit2, joint = FALSE)

## R Code 8.12

# test location-shift hypothesis
ff <- fev ~ I(age - 10) + I(height - 61) + sex
kt <- KhmaladzeTest(ff, tau = 1:19/20, nullH = "location", data = fev)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)
# test location--scale-shift hypothesis
kt <- KhmaladzeTest(ff, tau = 1:19/20, nullH = "location-scale", data = fev)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)
# estimate the quantile process
fit.rqs <- rq(ff, tau = -1, data = fev) # use plot(fit.rqs) for plotting
# number of breakpoints
with(fit.rqs, ncol(sol))-1

## R Code 8.13

gg <- I(log(fev)) ~ I(age - 10) + I(height - 61) + sex
kt <- KhmaladzeTest(gg, tau = 1:19/20, nullH = "location", data = fev)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)

## Figure 8.9

dd1 <- data.frame(t(fit.rqs$sol[c(1,4),]))
dd2 <- data.frame(t(fit.rqs$sol[c(1,5),]))
dd3 <- data.frame(t(fit.rqs$sol[c(1,6),]))
dd4 <- data.frame(t(fit.rqs$sol[c(1,7),]))
names(dd1) <- names(dd2) <- names(dd3) <- names(dd4) <- c("x", "y")

lsz <- 1
sz <- 7

gg <- list()

xmin <- 0.9
ddr <- dd1[dd1$x > xmin,]

inset1 <- ggplot(dd1, aes(x = x, y = y)) +
  geom_line(linewidth = 0.5) +
  coord_cartesian(xlim = c(xmin, 1), ylim = c(3, 3.5), expand = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_bw(base_size = sz) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6),
    plot.margin = margin(1, 1, 1, 1),
    axis.title = element_blank()
  )

gg[[1]] <- ggplot(dd1, aes(x = x, y = y)) + geom_line(linewidth = lsz) + labs(x = "", y = "Regression quantile", title = "Intercept") + coord_cartesian(ylim = c(1.5, 3.7), clip = "off") + annotation_custom(ggplotGrob(inset1), xmin = 0.02, xmax = 0.65, ymin = 2.75, ymax = 3.75) + geom_line(data = ddr, linewidth = lsz * 1.3, colour = "grey")
gg[[2]] <- ggplot(dd2, aes(x = x, y = y)) + geom_line(linewidth = lsz) + xlab("Quantile index") + ylab("Regression quantile") + ggtitle("Age") + ylim(c(0,0.2))
gg[[3]] <- ggplot(dd3, aes(x = x, y = y)) + geom_line(linewidth = lsz) + xlab("") + ylab("") + ggtitle("Height") + ylim(c(0,0.13))
gg[[4]] <- ggplot(dd4, aes(x = x, y = y)) + geom_line(linewidth = lsz) + xlab("Quantile index") + ylab("") + ggtitle("Sex (male)")

marrangeGrob(gg, ncol = 2, nrow = 2, top = "")


## Figure 8.10

fit.rqs <- rq(ff, tau = -1, data = fev)
p <- fit.rqs$sol["tau",]
beta0_hat <- fit.rqs$sol["(Intercept)",]
sf <- rearrange(stepfun(p, c(-Inf, beta0_hat)))
beta0_hat <- sf(p)

xmin <- 0.9
dd1 <- data.frame(t(fit.rqs$sol[c(1,4),]))
names(dd1) <- c("x", "y")
dd1 <- dd1[dd1$x > xmin,]
dd2 <- data.frame(x = p, y = beta0_hat)
dd2 <- dd2[dd2$x > xmin,]

lsz <- 1

gg <- list()

gg[[1]] <- ggplot(dd1, aes(x = x, y = y)) + geom_line(linewidth = lsz) + xlab("Quantile index") + ylab("Regression quantile") + ggtitle("")
gg[[2]] <- ggplot(dd2, aes(x = x, y = y)) + geom_line(linewidth = lsz) + xlab("Quantile index") + ylab("") + ggtitle("")

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")

## Code fragment

p <- fit.rqs$sol["tau",]
Q0hat <- fit.rqs$sol["(Intercept)",]
sf <- rearrange(stepfun(p, c(-Inf, Q0hat)))
Q0hat <- sf(p)
plot(p, Q0hat, type = "l", xlim = c(0.9,1), ylim = c(3, 3.5))

## R Code 8.14

p <- c(0.1, 1:3/4, 0.9)
fit <- rq(fev ~ I(age - 10) + I(height - 61) + sex, tau = p, data = fev)
Qhat <- predict(fit)
colnames(Qhat) <- p
print(sel <- which(t(apply(Qhat, 1, diff)) < 0, arr.ind = TRUE))
fev[sel[1],]
Qhat[sel[1],]

## R Code 8.15

# standard quantile regression
fit.rq <- rq(dist ~ speed, tau = 1:19/20, data = cars)
Qhat <- predict(fit.rq)

# total number of quantile crossings across 19 vigintiles
sum(t(apply(Qhat, 1, diff)) < 0)
# test location-scale-shift hypothesis
kt <- KhmaladzeTest(dist ~ speed, tau = 1:19/20,
nullH = "location-scale", data = cars)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)
# restricted quantile regression
fit.rrq <- Qtools::rrq(dist ~ speed, tau = 1:19/20, data = cars)

## Figure 8.11

transalpha <- 0.5
lsz <- 1
gg <- list()

b <- coef(fit.rq)
k <- ncol(b)
gg[[1]] <- ggplot(cars, aes(x = speed, y = dist)) + geom_point(colour = "grey") + xlab("Speed (mph)") + ylab("Stopping distance (ft)") + ggtitle("")
for(i in 1:k) gg[[1]] <- gg[[1]] + geom_abline(intercept = b[1,i], slope = b[2,i], linetype = "dashed", linewidth = 0.3)
gg[[1]] <- gg[[1]] + geom_abline(intercept = b[1,c(5,10,15)], slope = b[2,c(5,10,15)], linetype = "solid", linewidth = lsz)

b <- coef(fit.rrq)
k <- ncol(b)
gg[[2]] <- ggplot(cars, aes(x = speed, y = dist)) + geom_point(colour = "grey") + xlab("Speed (mph)") + ylab("") + ggtitle("")
for(i in 1:k) gg[[2]] <- gg[[2]] + geom_abline(intercept = b[1,i], slope = b[2,i], linetype = "dashed", linewidth = 0.3)
gg[[2]] <- gg[[2]] + geom_abline(intercept = b[1,c(5,10,15)], slope = b[2,c(5,10,15)], linetype = "solid", linewidth = lsz)

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")


## R Code 8.16

summary(rq(dist ~ speed, data = cars, tau = 0.975), se = "nid")
h <- bandwidth.rq(0.975, nrow(cars), hs = TRUE)/2
Qhat <- predict(rq(dist ~ speed, data = cars, tau = c(0.975-h, 0.975+h)))
dQ <- apply(Qhat, 1, diff)
dQ[which(!(dQ > 0))]
cars[which(!(dQ > 0)),]

## R Code 8.17

set.seed(123)
nn <- 100
p <- 0.75
R <- 1000
y <- sort(rnorm(nn, 5, 2))
yr <- round(y)
ystar <- replicate(R, yr + runif(nn, -0.5, 0.5))

# original
quantile(y, probs = p, names = FALSE)
# rounded
quantile(yr, probs = p, names = FALSE)
# jittered
mean(apply(ystar, 2, quantile, probs = p))

