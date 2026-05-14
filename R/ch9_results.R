# Statistical Methods and Applications Using Quantiles - Chapter 9

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(xtable)
library(Qtools)
library(nlme)
library(splines)
library(pspline)
library(ggplot2)

###################################################
# Simulation

set.seed(765)
R <- 500
tau <- 0.75
theta <- 2
b <- c(0, -2)
n <- 1000
res <- array(NA, dim = c(3,2,R,2))
nr <- array(NA, dim = c(R,2))
colnames(nr) <- dimnames(res)[[4]] <- c("mle","nlrq")

for(i in 1:R){
print(i)
x <- runif(n, -2, 2)
y <- rbeta(n, shape1 = exp(x), shape2 = exp(2*x))
y[y == 0] <- 0.001
y[y == 1] <- 0.999
dd <- data.frame(x = x, y = y)

fit.uq <- try(unitquantreg::unitquantreg(formula = y ~ x,
data = dd, family = "ulogistic", tau = tau, link = "logit"),
silent = TRUE)
fit.nlrq <- try(nlrq(y ~ SSlogis(x, Asym, xmid, scal),
data = dd, tau = tau), silent = TRUE)

if(!any(inherits(fit.uq, "try-error"), inherits(fit.nlrq, "try-error"))){
	res[,,i,1] <- summary(fit.uq)[[1]][,1:2]
	res[,,i,2] <- summary(fit.nlrq)$coefficients[,1:2]
	nr[i,1] <- mean(y - predict(fit.uq, type = "quantile") < 0)
	nr[i,2] <- mean(y - predict(fit.nlrq) < 0)
}
}

apply(res[2,2,,], 2, mean, na.rm = TRUE) # mean SE
apply(res[2,1,,], 2, sd, na.rm = TRUE) # SD of estimates
apply(nr, 2, mean, na.rm = TRUE) # prop. negative residuals


###################################################
# Save output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results"

save(list = c("nr","res"), file = file.path(path.out, "ch9_results.Rdata"), compress = "gzip")

