# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 1
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
### Section 1.2
###################################################

load(file = "../data/USbirths2001.Rdata")

p <- c(0.1,0.5,0.9)

sel <- USbirths2001$ga >= 37
fit.lm <- lm(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001, subset = sel)
fit.rq <- rq(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001, tau = p, method = "fn", subset = sel)
fit.rqs <- summary(fit.rq, se = "nid")

## Table 1.1

coefs.lm <- summary(fit.lm)$coefficients
coefs.rq <- lapply(fit.rqs, function(x) x$coefficients)

std.m <- coefs.lm[,2]
std.q <- sapply(coefs.rq, function(x) x[,2])

est.q <- cbind(sapply(coefs.rq, function(x) x[,1]),std.q)[,c(1,4,2,5,3,6)]
est.m <- cbind(coefs.lm[,1], std.m)

nn <- c("Intercept", "Sex (female)", "Multiple births", "Maternal age", "Smoking (yes)", "Drinking (yes)")
tab1 <- data.frame(est.m,est.q)
colnames(tab1) <- rep(c("Est.","SE"),4)
rownames(tab1) <- nn

print(xtable(tab1, digits = 1, align = c("l",rep("r",8)), caption = c("Birth weight mean and quantile effects for infants born at term, United States, 2001. Standard errors are reported in brackets.")), booktabs = TRUE, math.style.negative = TRUE)
