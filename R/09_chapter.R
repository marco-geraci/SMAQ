# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 9
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
## Section 9.1
###################################################

## R Code 9.1

library(quantreg)
data(Indometh)
# fit nonlinear median regression
fit.nlrq <- nlrq(conc ~ SSbiexp(time, A1, lrc1, A2, lrc2),
data = Indometh, tau = 0.5)
summary(fit.nlrq)
# fit nonlinear mean regression
fit.nls <- nls(conc ~ SSbiexp(time, A1, lrc1, A2, lrc2), data = Indometh)
summary(fit.nls)


# Figure 9.1

x <- Indometh$time
y <- Indometh$conc

p <- c(.1,.5,.9)
nt <- length(p)
fit <- list()
for(i in 1:nt){
    fit[[i]] <- nlrq(conc ~ SSbiexp(time, A1, lrc1, A2, lrc2),
    data = Indometh, tau = p[i])
}

nx <- 50
xx <- seq(min(x), max(x), length = nx)
yy <- c(
do.call(SSbiexp, args = c(list(input = xx), as.list(coef(fit[[1]])))),
do.call(SSbiexp, args = c(list(input = xx), as.list(coef(fit[[2]])))),
do.call(SSbiexp, args = c(list(input = xx), as.list(coef(fit[[3]]))))
)

dd <- data.frame(x = x, y = y)
qq <- data.frame(x = rep(xx, nt), y = as.numeric(yy), p = factor(rep(p, each = nx)))

psz <- 2
lsz <- 1
transalpha <- 0.4

ggplot() + geom_point(data = dd, aes(x = x, y = y),  shape = 16, size = psz, alpha = transalpha) +
geom_line(data = qq, aes(x = x, y = y, group = p, color = p), linetype = "solid", linewidth = lsz)  + scale_color_grey(start = 0.8, end = 0) + theme(legend.position = "none") + xlab("Time (hours after injection)") + ylab("Concentration (mcg/ml)")

## R Code 9.2

set.seed(123)

n <- 5000
x <- runif(n, 0, 4)
eta <- 2 + 0.5*x
y <- exp(eta)*rnorm(n,1,0.2)
all(y > 0)
z <- log(y)

# fit nonlinear model
system.time(fit.nlrq <- nlrq(y ~ exp(a + b*x), tau = 0.9,
start = list(a = 0, b = 0)))
summary(fit.nlrq, R = 200)

# fit log-linear model
system.time(fit.rq <- rq(z ~ x, tau = 0.9, method = "fn"))
summary(fit.rq, se = "boot", R = 200)

###################################################
## Section 9.2
###################################################

# Figure 9.2


library(grid)
library(gtable)

lambda <- c(-1, -0.5, 0, 0.5, 1)
lambda_pos <- c(0, 0.5, 1)

y.pos <- seq(0.05, 5, length.out = 500)
y.unit <- seq(0.001, 0.999, length.out = 500)

make_data <- function(y, fun, name, lambda) {
  do.call(rbind, lapply(lambda, function(l) {
    data.frame(y = y,
               h = fun(y, l),
               lambda = factor(l, levels = c(-1, -0.5, 0, 0.5, 1)),
               transformation = name)
  }))
}

dd.bc <- make_data(y.pos,  function(y, l) Qtools::bc(y, l), "Box--Cox", lambda)
dd.ao <- make_data(y.unit, function(y, l) Qtools::ao(y, l, symm = TRUE), "Symmetric AO", lambda_pos)
dd.mcja <- make_data(y.unit, function(y, l) Qtools::mcjI(y, l, symm = FALSE, dbounded = TRUE), "Asymmetric MCJ", lambda_pos)
dd.mcjs <- make_data(y.unit, function(y, l) Qtools::mcjI(y, l, symm = TRUE, dbounded = TRUE), "Symmetric MCJ", lambda_pos)

dd.bc <- subset(dd.bc, is.finite(h))
dd.ao <- subset(dd.ao, is.finite(h))
dd.mcja <- subset(dd.mcja, is.finite(h))
dd.mcjs <- subset(dd.mcjs, is.finite(h))

plot_fun <- function(d, title, show.legend = FALSE) {
  ggplot(d, aes(x = y, y = h, linetype = lambda)) +
    geom_line(linewidth = 0.6) + coord_cartesian(ylim = c(-6, 6)) + labs(x = "Y", y = "Transformed value", linetype = expression(lambda)) + ggtitle(title) +  theme(legend.position = if (show.legend) "top" else "none", legend.justification = "center", plot.title = element_text(hjust = 0.5))
}

gg <- list()

gg[[1]] <- ggplot(dd.bc, aes(x = y, y = h, linetype = lambda)) +
  geom_line(linewidth = lsz, colour = "black") +
  coord_cartesian(ylim = c(-6, 6)) +
  labs(x = "", y = expression(h(Y)), linetype = expression(lambda)) +
  ggtitle("Box-Cox") +
  theme(
  plot.title = element_text(hjust = 0.5),
  legend.position = "top",
  legend.justification = "center",
  legend.margin = margin(0,0,0,0),
  legend.box.margin = margin(50,0,-5,0),
  legend.key.height = unit(1, "cm"),
  legend.key.width = unit(3, "line"),
  legend.text = element_text(size = 15),
  legend.title = element_text(size = 15),
  legend.spacing.x = unit(0.7, 'cm')
) + scale_linetype_manual(
  labels = parse(text = c("-1", "-0.5", "0", "0.5", "1")),
  values = c("-1" = "dashed",
             "-0.5" = "dotdash",
             "0" = "solid",
             "0.5" = "dotted",
             "1" = "longdash"),
  drop = FALSE
)

gg[[2]] <- ggplot(dd.ao, aes(x = y, y = h, linetype = lambda)) +
  geom_line(linewidth = lsz, colour = "black") +
  coord_cartesian(ylim = c(-6, 6)) +
  labs(x = "", y = "",
       linetype = expression(lambda)) +
  ggtitle("Symmetric AO") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") + scale_linetype_manual(
  values = c("-1" = "dashed",
             "-0.5" = "dotdash",
             "0" = "solid",
             "0.5" = "dotted",
             "1" = "longdash"),
  drop = FALSE
)

gg[[3]] <- ggplot(dd.mcja, aes(x = y, y = h, linetype = lambda)) +
  geom_line(linewidth = lsz, colour = "black") +
  coord_cartesian(ylim = c(-6, 6)) +
  labs(x = "Y", y = "",
       linetype = expression(lambda)) +
  ggtitle("Asymmetric MCJ") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") + scale_linetype_manual(
  values = c("-1" = "dashed",
             "-0.5" = "dotdash",
             "0" = "solid",
             "0.5" = "dotted",
             "1" = "longdash"),
  drop = FALSE
)


gg[[4]] <- ggplot(dd.mcjs, aes(x = y, y = h, linetype = lambda)) +
  geom_line(linewidth = lsz, colour = "black") +
  coord_cartesian(ylim = c(-6, 6)) +
  labs(x = "Y", y = expression(h(Y)),
       linetype = expression(lambda)) +
  ggtitle("Symmetric MCJ") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "none") + scale_linetype_manual(
  values = c("-1" = "dashed",
             "-0.5" = "dotdash",
             "0" = "solid",
             "0.5" = "dotted",
             "1" = "longdash"),
  drop = FALSE
)

leg <- gtable_filter(ggplotGrob(gg[[1]]), "guide-box")

gg[[1]] <- gg[[1]] + theme(legend.position = "none")

g <- arrangeGrob(gg[[1]], gg[[2]], gg[[3]], gg[[4]], nrow = 2, ncol = 2)

grid.newpage()
grid.draw(
arrangeGrob(leg, g, ncol = 1,  heights = c(0.12, 1))
)

## R Code 9.3

library(Qtools)
library(quantreg)
data(airquality)

# fit linear model
fit.rq <- rq(Ozone ~ Solar.R, data = airquality, tau = 0.9)
summary(fit.rq, se = "nid")
# CUSUM lack-of-fit test
LOFTest(fit.rq, alpha = 0.05, B = 1000, seed = 987)

## R Code 9.4

# fit transformation model
fit.rqt <- tsrq(Ozone ~ Solar.R, data = airquality, tsf = "mcjI",
symm = TRUE, dbounded = FALSE, lambda = seq(1, 3, by = 0.005),
conditional = FALSE, tau = .9)
summary(fit.rqt, conditional = FALSE, se = "nid")
# CUSUM lack-of-fit test
LOFTest(fit.rqt, alpha = 0.05, B = 1000, seed = 987)

## Figure 9.3


p <- c(.1,.5,.9)
nt <- length(p)

fit <- tsrq(Ozone ~ Solar.R, data = airquality, tsf = "mcjI", symm = TRUE, dbounded = FALSE, lambda = seq(1, 3, by = 0.005),
conditional = FALSE, tau = p)

nx <- 200
xx <- seq(9, 334, length = nx)
qhat <- predict(fit, newdata = data.frame(Solar.R = xx), type = "response")
invisible(capture.output(dqhat <- predict(fit, newdata = data.frame(Solar.R = xx), type = "maref", namevec = "Solar.R")))

x <- airquality$Solar.R
y <- airquality$Ozone

dd <- data.frame(x = x, y = y)
qq1 <- data.frame(x = rep(xx, nt), y = as.numeric(qhat), p = factor(rep(p, each = nx)))
qq2 <- data.frame(x = rep(xx, nt), y = as.numeric(dqhat), p = factor(rep(p, each = nx)))

psz <- 2
lsz <- 1
transalpha <- 0.4

gg <- list()

gg[[1]] <- ggplot() + geom_point(data = dd, aes(x = x, y = y),  shape = 16, size = psz, alpha = transalpha) +
geom_line(data = qq1, aes(x = x, y = y, group = p, color = p), linetype = "solid", linewidth = lsz) + scale_color_grey(start = 0.8, end = 0) + theme(legend.position = "none") + xlab("Solar radiation (lang)") + ylab("Ozone (ppb)")

gg[[2]] <- ggplot() + geom_line(data = qq2, aes(x = x, y = y, group = p, color = p), linetype = "solid", linewidth = lsz) + scale_color_grey(start = 0.8, end = 0) + theme(legend.position = "none")  + xlab("Solar radiation (lang)") + ylab("QTE (ppb/lang)")

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")

###################################################
## Section 9.3
###################################################

## R Code 9.5

set.seed(123)
n <- 100
p <- 0.5
x <- sort(runif(n, 0, 2*pi))
y <- sin(x) + (1 + x)*rnorm(n, 0, 0.1)

# lambda set at an arbitrary value
fit.L2 <- Qtools::qsmspline(x = x, y = y, p = p, lambda = 0.1)
fit.L1 <- quantreg::rqss(y ~ qss(x, constraint= "N", lambda = 0.1), tau = p)

## Figure 9.4

psz <- 2
lsz <- 1
transalpha <- 0.4

dd <- data.frame(x = x, y = y)
qq <- data.frame(x = rep(x, 2), y = c(fit.L2[[5]], fit.L1$coefficients), penalty = factor(rep(c("L2", "L1"), each = n), levels = c("L1", "L2")))

g <- ggplot() + geom_point(data = dd, aes(x = x, y = y),  shape = 16, size = psz, alpha = transalpha) +
geom_line(data = qq, aes(x = x, y = y, linetype = penalty), linewidth = lsz, color = "black") + theme(legend.title = element_blank(), legend.position = "top", legend.spacing.x = unit(0.7, 'cm'), legend.key.width = unit(2, "line"), legend.text = element_text(size = 11)) + xlab("x") + ylab("y")
g


## R Code 9.6

# define grid - must be wide enough to include the minimum
lambda <- seq(0.1, 2, by = 0.1)
# selection of lambda for rqss based on SIC
sic <- rep(NA, length(lambda))
for(j in 1:length(lambda)){
    fit <- quantreg::rqss(y ~ qss(x, constraint= "N", lambda = lambda[j]),
    tau = p)
    sic[j] <- log(fit$fidelity/n) + 0.5/n*fit$edf*log(n)
}

# value of lambda that minimises SIC
lambda[which.min(sic)]

## R Code 9.7

set.seed(123)
n <- 100
p <- 0.5
x <- sort(runif(n, 0, 2*pi))
y <- sin(x) + (1 + x)*rnorm(n, 0, .1)
h_mean <- KernSmooth::dpill(x, y)
h <- h_mean * (p * (1 - p)/(dnorm(qnorm(p)))^2)^0.2
yhat <- rep(NA, n)
for(i in 1:n){
    x.c <- x[i] - x
    w <- dnorm(x.c/h)
    yhat[i] <- quantreg::rq(y ~ x.c, weights = w, tau = p)$coef[1]
}

## Figure 9.5

psz <- 2
lsz <- 1
transalpha <- 0.4

dd <- data.frame(x = x, y = y)
qq <- data.frame(x = x, y = yhat)

g <- ggplot() + geom_point(data = dd, aes(x = x, y = y),  shape = 16, size = psz, alpha = transalpha) +
geom_line(data = qq, aes(x = x, y = y), linetype = "solid", linewidth = lsz, color = "black") + xlab("x") + ylab("y")
g

###################################################
## Section 9.4
###################################################

## R Code 9.8

set.seed(123)
n <- 500
X <- cbind(1, x <- runif(n, 0, 5))
b <- c(0,1)
g <- c(0.5,0.5)
y <- X%*%b + exp(X%*%g)*rnorm(n)

ss <- list(a = 0.1, b = 0.1, sigma = 0.1, g = 0.1)
try(nlrq(y ~ a + b*x + sigma*exp(g*x), tau = 0.1, start = ss))
try(nlrq(y ~ a + b*x + sigma*exp(g*x), tau = 0.5, start = ss))
try(nlrq(y ~ a + b*x + sigma*exp(g*x), tau = 0.9, start = ss))

## R Code 9.9

library(nlme)
# fit a generalised least squares regression
# with exponential var function
fit <- gls(y ~ x, weights = varExp(form = ~ x))
# the fitted scale parameter sigma estimates
# exp(g[1]) from the data-generating model
gamma1 <- as.numeric(coef(fit$modelStruct))
print(thetahat <- c(coef(fit), gamma0 = log(fit$sigma),
gamma1 = gamma1))
# predict quantiles
probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
xx <- seq(0, 5, by = 0.05)
mu <- cbind(1, xx)%*%thetahat[c(1:2)]
sigma <- exp(cbind(1, xx)%*%thetahat[c(3:4)])
qhat <- sapply(probs, function(p) qnorm(p, mean = mu, sd = sigma))
# estimate QTE
qte <- gamma1*fit$sigma*exp(gamma1*xx)
qte <- coef(fit)[2] + outer(qte, qnorm(probs), "*")

## Figure 9.6

psz <- 2
lsz <- 1
transalpha <- 0.4

colnames(qhat) <- probs
df <- reshape2::melt(qhat)
df$x <- rep(xx, length(probs))
df$p <- factor(df$Var2)

colnames(qte) <- probs
qte <- reshape2::melt(qte)
qte$x <- rep(xx, length(probs))
qte$p <- factor(df$Var2)

dd <- data.frame(x = x, y = y)

gg <- ggplot() + geom_point(data = dd, aes(x = x, y = y),  shape = 16, size = psz, alpha = transalpha) + geom_line(data = df, aes(x = x, y = value, group = p, color = p), linewidth = lsz) + scale_color_grey(start = 0.8, end = 0) + theme(legend.title = element_blank(), legend.position = "none") + xlab("x") + ylab("y")

qq <- ggplot() + geom_line(data = qte, aes(x = x, y = value, group = p, color = p), linewidth = lsz) + theme(legend.title = element_blank(), legend.position = "none") + xlab("x") +  ylab("QTE") + scale_color_grey(start = 0.8, end = 0)

grid.arrange(gg, qq, heights = c(2, 1), ncol = 1)

## R Code 9.10

qgg <- function(p, mu, sigma, kappa){
	qp <- qgamma(p, shape = kappa, scale = 1)/kappa
	exp(mu)*qp^(sigma*sqrt(kappa))
}

dlgg <- function(x, mu, sigma, kappa, log = FALSE){
	z <- (x - mu)/sigma
	val1 <- (kappa - 0.5)*log(kappa) - log(sigma) - lgamma(kappa)
	val2 <- z*sqrt(kappa) - kappa*exp(z/sqrt(kappa))
	val <- val1 + val2
	if(!log) val <- exp(val)
	return(val)
}

loglik <- function(theta, x, y){
	p <- ncol(x)
	mu <- x %*% theta[1:p]
	sigma <- exp(x %*% theta[(p+1):(2*p)])
	kappa <- exp(x %*% theta[(2*p+1):(3*p)])
	val <- -sum(dlgg(y, mu, sigma, kappa, log = TRUE))
	return(val)
}

# obtain the MLE using the GG likelihood
data(ImmunogG, package = "Brq")
fit <- do.call(optim, args = list(par = rep(1, 6), fn = loglik, x = cbind(1,ImmunogG$Age), y = log(ImmunogG$IgG), method = "BFGS"))
-fit$value # loglik
fit$par # estimates
# predict quantiles using the GG distribution
xx <- seq(0, 6, by = 0.05)
muhat <- cbind(1, xx)%*%fit$par[c(1:2)]
sigmahat <- exp(cbind(1, xx)%*%fit$par[c(3:4)])
kappahat <- exp(cbind(1, xx)%*%fit$par[c(5:6)])
probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
qhat <- sapply(probs, function(p) qgg(p, mu = muhat, sigma = sigmahat, kappa = kappahat))


## R Code 9.11

# define quantile as a function of x, given the parameter theta and
# probability p
Q <- function(x, theta, p){
    mu <- sum(c(1,x)*theta[1:2])
    sigma <- exp(sum(c(1,x)*theta[3:4]))
    kappa <- exp(sum(c(1,x)*theta[5:6]))
    qgg(p, mu, sigma, kappa)
}
# estimate the QTE for the quantile, e.g., 0.1 along a sequence xx
library(numDeriv)
qte <- sapply(xx, function(x) grad(func = Q, x, theta = fit$par, p = 0.1))

## Figure 9.7

psz <- 2
lsz <- 1
transalpha <- 0.4

xx <- seq(0, 6, by = 0.05)
mu <- cbind(1, xx)%*%fit$par[c(1:2)]
sigma <- exp(cbind(1, xx)%*%fit$par[c(3:4)])
kappa <- exp(cbind(1, xx)%*%fit$par[c(5:6)])
probs <- c(0.10, 0.25, 0.50, 0.75, 0.90)
qhat <- sapply(probs, function(p) qgg(p, mu = mu, sigma = sigma, kappa = kappa))
colnames(qhat) <- probs
df <- reshape2::melt(qhat)
df$Age <- rep(xx, length(probs))
df$p <- factor(df$Var2)

gg <- ggplot() + geom_point(data = ImmunogG, aes(x = Age, y = IgG),  shape = 16, size = psz, alpha = transalpha) + geom_line(data = df, aes(x = Age, y = value, group = p, color = p), linewidth = lsz) + theme(legend.title = element_blank(), legend.position = "none") + xlab("Age (years)") + ylab("IgG (g/l)")  + scale_color_grey(start = 0.8, end = 0)

Q <- function(x, p, theta){
	mu <- sum(c(1,x)*theta[1:2])
	sigma <- exp(sum(c(1,x)*theta[3:4]))
	kappa <- exp(sum(c(1,x)*theta[5:6]))
	qgg(p, mu, sigma, kappa)
}

qte <- matrix(NA, length(xx), length(probs))
for(i in 1:length(probs)){
qte[,i] <- sapply(xx, function(x) numDeriv::grad(func = Q, x, theta = fit$par, p = probs[i]))
}
colnames(qte) <- probs
qte <- reshape2::melt(qte)
qte$Age <- rep(xx, length(probs))
qte$p <- factor(qte$Var2)

qq <- ggplot() + geom_line(data = qte, aes(x = Age, y = value, group = p, color = p), linewidth = lsz) + theme(legend.title = element_blank(), legend.position = "none") + xlab("Age (years)") +  ylab("QTE ((g/l)/years)")  + scale_color_grey(start = 0.8, end = 0)

grid.arrange(gg, qq, heights = c(2, 1), ncol = 1)

## R Code 9.12

set.seed(765)
R <- 500
p <- 0.75
b <- c(0, -2)
n <- 1000
res <- array(NA, dim = c(3,2,R,2))
nr <- array(NA, dim = c(R,2))
colnames(nr) <- dimnames(res)[[4]] <- c("unitquantreg","nlrq")

for(i in 1:R){
    x <- runif(n, -2, 2)
    y <- rbeta(n, shape1 = exp(x), shape2 = exp(2*x))
    y[y == 0] <- 0.001
    y[y == 1] <- 0.999
    dd <- data.frame(x = x, y = y)

    fit.uq <- try(unitquantreg::unitquantreg(formula = y ~ x,
    data = dd, family = "ulogistic", tau = p, link = "logit"),
    silent = TRUE)
    fit.nlrq <- try(nlrq(y ~ SSlogis(x, Asym, xmid, scal),
    data = dd, tau = p), silent = TRUE)

    if(!any(inherits(fit.uq, "try-error"), inherits(fit.nlrq, "try-error")))
    {res[,,i,1] <- summary(fit.uq)[[1]][,1:2]
    res[,,i,2] <- summary(fit.nlrq)$coefficients[,1:2]
    nr[i,1] <- mean(y - predict(fit.uq, type = "quantile") < 0)
    nr[i,2] <- mean(y - predict(fit.nlrq) < 0)}
}

## Table 9.1

tab1 <- rbind(
apply(res[2,2,,], 2, mean, na.rm = TRUE), # mean SE
apply(res[2,1,,], 2, sd, na.rm = TRUE), # SD of estimates
apply(nr, 2, mean, na.rm = TRUE) # prop. negative residuals
)

tab1 <- data.frame(tab1)
rownames(tab1) <- c("Mean SE", "SD of estimates", "PNR")
colnames(tab1) <- c("Est.", "SE")

addtorow <- list()
addtorow$pos <- list(-1,0,0,nrow(tab1))
addtorow$command <- c("\\toprule",
" & \\multicolumn{1}{c}{unitquantreg} & \\multicolumn{1}{c}{nlrq} \\\\",
"\\cmidrule(l){2-3}",
"\\bottomrule"
)

print.xtable(xtable(tab1, digits = c(0,4,4), align = c(rep("l",1),rep("r",2)), caption = c("Simulation results from R Code~\\ref{code:09_12}, showing the mean standard error (SE), standard deviation (SD) of the estimates, and proportion of negative residuals (PNR) for \\texttt{unitquantreg} and \\texttt{nlrq}.\\label{tab:09_01}")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.colnames = FALSE, hline.after = NULL, format.args = list(big.mark = ",", decimal.mark = "."))


## R Code 9.13

dlogis2 <- function(y, gamma, delta, log = FALSE){
	num <- delta*exp(gamma)*(y/(1-y))^delta
	den <- y*(1-y)*(1 + exp(gamma)*(y/(1-y))^delta)^2
	val <- num/den
	if(log) val <- log(val)
	return(val)
}

loglik_logis <- function(theta, x, y){
	k <- ncol(x)
	gamma <- x%*%theta[1:k]
	delta <- exp(theta[k+1])
	-sum(log(dlogis2(y, gamma, delta)))
}
# generate data from (unit) logistic with given p0
set.seed(765)
p0 <- 0.75
n <- 1000; delta <- 2; beta_p0 <- c(0, -2)
x <- runif(n, -2, 2)
eta <- cbind(1,x)%*%beta_p0
mu_p0 <- exp(eta)/(1 + exp(eta))
y <- unitquantreg::rulogistic(n, mu = mu_p0, theta = delta, tau = p0)
dd <- data.frame(x = x, y = y)

# fit (true) density with quantile-based parameterisation
fit.uq <- unitquantreg::unitquantreg(formula = y ~ x, data = dd,
family = "ulogistic", tau = p0, link = "logit")

# fit original density (i.e, not reparameterised)
fit.mle <- optim(par = c(0,0,0), loglik_logis, x = cbind(1,x), y = y,
hessian = TRUE)
deltahat <- exp(fit.mle$par[3])
# obtain the same estimates obtained with unitquantreg
betahat <- -(fit.mle$par[1:2] + c(log((1-p0)/p0), 0))/deltahat
# obtain standard errors with delta method
V_thetahat <- solve(fit.mle$hessian)
G <- matrix(0, nrow = 2, ncol = 3)
G[1,1] <- -1/deltahat; G[2,2] <- -1/deltahat; G[,3] <- -betahat
V_betahat <- G %*% V_thetahat %*% t(G)
se_betahat <- sqrt(diag(V_betahat))

## Table 9.2

tab2 <- data.frame(cbind(summary(fit.uq)[[1]][1:2,1:2], betahat, se_betahat))
rownames(tab2) <- c("Intercept", "Slope")
colnames(tab2) <- rep(c("Est.", "SE"), 2)

addtorow <- list()
addtorow$pos <- list(-1,0,0,0,nrow(tab2))
addtorow$command <- c("\\toprule",
" & \\multicolumn{2}{c}{unitquantreg} & \\multicolumn{2}{c}{MLE} \\\\",
"\\cmidrule(l){2-3}",
"\\cmidrule(l){4-5}",
"\\bottomrule"
)

print.xtable(xtable(tab2, digits = c(0,4,4,4,4), align = c(rep("l",1),rep("r",4)), caption = c("Estimates and standard errors (SE) of the $p_{0}$-specific regression coefficients in R Code \\ref{code:09_13} using \\texttt{unitquantreg} and standard maximum likelihood estimation (MLE).\\label{tab:09_02}")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.colnames = FALSE, hline.after = NULL, format.args = list(big.mark = ",", decimal.mark = "."))

## R Code 9.14

fit.cdfq <- cdfquantreg::cdfquantreg(y ~ x | 1, data = dd, family = "logit-logistic")
# account for parameterisation used in Smithson and Shou (2017)
d <- exp(-cdfquantreg:::coef.cdfqr(fit.cdfq)[3])
as.numeric(c(d*(-cdfquantreg:::coef.cdfqr(fit.cdfq)[1:2]), d))
# compare with MLE in previous R Code
c(fit.mle$par[1:2], exp(fit.mle$par[3]))

