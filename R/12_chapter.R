# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 12
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
## Section 12.1
###################################################

## R Code 12.1

load(file = "../data/USbirths2001.Rdata")

# singletons born at term only
dd <- subset(USbirths2001, USbirths2001$ga >= 37 & USbirths2001$multiple == 0)
nrow(dd) # about 3 million

# fit logistic regression on LBW indicator
# conditional on sex (x = 1 if female)
y <- dd$bw
d <- 1*(y < 2500)
x <- dd$female
fit <- glm(d ~ x, family = binomial)

# conditional probabilities of birthweight less than 2500
xx <- cbind(1, c(0,1))
eta <- xx %*% coef(fit)
phat <- as.numeric(1/(1 + exp(-eta)))
names(phat) <- c("males", "females")
phat

# verify that the estimated probabilities correspond to
# the same empirical quantile (i.e., 2500)
round(quantile(y[x == 0], probs = phat[1], names = FALSE))
round(quantile(y[x == 1], probs = phat[2], names = FALSE))

###################################################
## Section 12.2
###################################################

## R Code 12.1

set.seed(123)
n <- 500
x <- sort(runif(n))
y <- (m <- 1 + 2*x) + (e <- rnorm(n))
length(unique(y)) == n
# > [1] TRUE
z <- sort(y)
range(pnorm(e))
# > [1] 0.00247881 0.99644571
taus <- seq(0.01, 0.99, by = 0.005)
F <- outer(m, z, function(mm, zz) pnorm(zz - mm))

# rq
t0 <- proc.time()['elapsed']
fit.rq <- rq(y ~ x, tau = taus)
Fhat_fun.rq <- predict(fit.rq, type = "Fhat", stepfun = TRUE)
Fhat.rq <- t(sapply(Fhat_fun.rq, function(F, z) F(z), z = z))
t1 <- proc.time()['elapsed']
print(t1 - t0)
# glm
t0 <- proc.time()['elapsed']
fit.glm <- lapply(z, function(zz, y, x)
  glm(1*(y <= zz) ~ x, family = binomial(link = "probit")),
  y = y, x = x)
Fhat.dr <- sapply(fit.glm, function(obj) obj$fitted)
Fhat_fun.dr <- apply(Fhat.dr, 1, function(u, z) stepfun(z, c(0,u)), z = z)
t1 <- proc.time()['elapsed']
print(t1 - t0)
mean(abs(Fhat.rq - F))
mean(abs(Fhat.dr - F))

## Figure 12.3

psz <- 2.5
lsz <- 1
transalpha <- 0.4
colg <- "grey70"

p <- 1:99/100
gr <- factor(rep(c("QR","DR"), each = length(z)), levels = c("QR", "DR"))

i <- 250
dd1 <- data.frame(x = c(z, z), y = c(Fhat_fun.rq[[i]](z), Fhat_fun.dr[[i]](z)), gr = gr)
dd2 <- data.frame(x = z, y = pnorm(z, mean = 1 + 2*x[i], sd = 1))

g <- ggplot(dd1, aes(x = x, y = y)) + geom_step(aes(linetype = gr)) + geom_line(data = dd2, colour = colg) + xlab("y") + ylab("Cumulative probability") + theme(legend.title = element_blank(), legend.position = "top", legend.spacing.x = unit(1.0, 'cm'), legend.key.width = unit(3, "line"), strip.background = element_rect(colour = colg, fill = colg), legend.text = element_text(size = 12))
g

###################################################
## Section 12.3
###################################################

## R Code 12.3

# fit logistic regression on LBW indicator
# conditional on multiple birth status (x = 0 if singleton)
nrow(USbirths2001) # about 3.5 million
y <- USbirths2001$bw
d <- 1*(y < 2500)
x <- 1*(USbirths2001$multiple == 1)

# LBW prevalence conditional on multiple status
print(pi0 <- mean(d[x == 0])) # (intercept conditional model)
print(pi1 <- mean(d[x == 1]))
# LBW prevalence difference (slope conditional model)
print(dp <- pi1 - pi0)
# probability linear model
fit <- lm(d ~ x)
fit
# LBW prevalence in the population (unconditional)
mean(d)
# expected percentage points change of LBW in the population
# for a change du in the singleton/multiple composition
du <- 0.1
(dp*du)*100

## R Code 12.4

y <- USbirths2001$bw
x <- 1*(USbirths2001$multiple == 1)
# multiple birth prevalence
print(px <- mean(x))
# LBW prevalence by sex
print(LBW0 <- mean(y[x == 0] <= 2500))
print(LBW1 <- mean(y[x == 1] <= 2500))
# LBW prevalence in the population
print(LBW <- mean(y <= 2500))
# define quantile levels
p <- round(c(LBW0, LBW), 3)
# fit unconditional and conditional QR
res <- matrix(NA, length(p), 2)
rownames(res) <- p
colnames(res) <- c("UQR", "QR")
for(j in 1:length(p)){
xi <- quantile(y, p[j])
fhat <- density(y)
w <- xi + (p[j] - 1*(y <= xi))/approxfun(fhat$x, fhat$y)(xi)
res[j,1] <- lm(w ~ x)$coef[2]
res[j,2] <- rq(y ~ x, tau = p[j], method = "fn")$coef[2]
}
round(res)
