# Statistical Methods and Applications Using Quantiles - Chapter 8

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(xtable)
library(Qtools)
library(nlme)
library(ggplot2)

###################################################
# Data

dd <- read.csv(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\_DONOTSHARE\\data_4.csv")

crohn <- subset(dd, dd$clinical.group %in% c("no disease", "control"))
crohn <- subset(crohn, !is.na(crohn$FC))
names(crohn)[names(crohn) == "clinical.group"] <- "group"
crohn <- crohn[,c("group","FC")]
crohn$group[crohn$group == "no disease"] <- "FDR"
crohn$group <- factor(crohn$group)
crohn <- crohn[order(crohn$group, crohn$FC),]
rownames(crohn) <- NULL
save(list = c("crohn"), file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\crohn.Rdata")

###################################################
# Analysis

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\crohn.Rdata")

fit.rq <- rq(FC ~ group, data = crohn)
fit.lm <- lm(FC ~ group, data = crohn)

summary(fit.rq, se = "nid")
round(summary(fit.lm)$coef, 4)


###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")

###################################################
# Analysis

# singletons only
dd <- subset(USbirths2001, USbirths2001$multiple == 0)
# some missing in ga and tox
dd <- na.omit(dd[,c("ga_group", "female", "tox", "died", "bw")])
# almost 3.5 million
nrow(dd)

fit <- rq(bw ~ ga_group + female, data = dd, method = "fn", tau = 0.1)
summary(fit, se = "nid")

# define sga as negative residuals
dd$sga <- as.numeric(residuals(fit) < 0)

fit.glm <- glm(died ~ female + ga_group + sga + tox, family = "binomial", data = dd)
round(exp(cbind(OR = coef(fit.glm), confint.default(fit.glm))), 3)


###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")

###################################################
# Analysis

# singletons only
dd <- subset(USbirths2001, USbirths2001$multiple == 0)

# sample 3.3 million infants who survived
dd <- dd[order(dd$cause),]
n <- sum(dd$cause == "ref")
set.seed(123); sel <- sample(1:n, 100e3)
k <- nrow(dd) - n
sel <- c(sel, (n+1):(n+k))
dd <- dd[sel,]

# compare bw
fit2 <- rq(bw ~ cause, tau = c(0.05, 0.33, 0.66, 0.95), data = dd, method = "fn")
fit.aov <- anova(fit2, joint = FALSE)

###################################################
# Save output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results"

save(list = c("fit2", "fit.aov", "dd"), file = file.path(path.out, "ch8_results.Rdata"), compress = "gzip")

###################################################
# Data

fev <- read.csv("C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\fev.csv", stringsAsFactors = TRUE)

save(list = c("fev"), file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\fev.Rdata")

###################################################
# Analysis
lsz <- 1
sz <- 7

kt <- KhmaladzeTest(fev ~ age + height + sex, tau = -1,
    nullH = "location", data = fev)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)
kt <- KhmaladzeTest(fev ~ age + height + sex, tau = -1,
    nullH = "location-scale", data = fev)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)

kt <- KhmaladzeTest(I(log(fev)) ~ age + height + sex, tau = -1,
    nullH = "location", data = fev)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)

anova(rq(log(fev) ~ age + height + sex, tau = 1:19/20, data = fev), joint = FALSE)

fit.rqs <- rq(fev ~ I(age - 10) + I(height - 61) + sex, tau = -1, data = fev)
ps <- fit.rqs$sol[1,]

qs <- c(c(1,0,0,0)%*%fit.rqs$sol[4:7,])

if(any(diff(qs) < 0)){
	sf <- rearrange(stepfun(ps, c(-Inf, qs)))
	qs <- sf(ps)
}


ps.wts <-(c(0,diff(ps))+c(diff(ps),0))/2
dens <- akj(qs, z = qs, p = ps.wts)$dens
dd1 <- data.frame(x = qs, y = dens)

qs <- c(c(1,0,0,1)%*%fit.rqs$sol[4:7,])
ps.wts <-(c(0,diff(ps))+c(diff(ps),0))/2
dens <- akj(qs, z = qs, p = ps.wts)$dens
dd2 <- data.frame(x = qs, y = dens)

qs <- c(c(1,19-10,0,0)%*%fit.rqs$sol[4:7,])
ps.wts <-(c(0,diff(ps))+c(diff(ps),0))/2
dens <- akj(qs, z = qs, p = ps.wts)$dens
dd3 <- data.frame(x = qs, y = dens)

qs <- c(c(1,0,74-61,0)%*%fit.rqs$sol[4:7,])
ps.wts <-(c(0,diff(ps))+c(diff(ps),0))/2
dens <- akj(qs, z = qs, p = ps.wts)$dens
dd4 <- data.frame(x = qs, y = dens)

gg <- list()
gg <- ggplot() + geom_line(data = dd1, aes(x = x, y = y), linetype = "solid", linewidth = lsz) +
geom_line(data = dd2, aes(x = x, y = y), linetype = "dashed", linewidth = lsz) +
geom_line(data = dd3, aes(x = x, y = y), linetype = "dotted", linewidth = lsz) +
geom_line(data = dd4, aes(x = x, y = y), linetype = "dotdash", linewidth = lsz) +
xlab("Birth weight (g)") +
ylab("Density")
print(gg)



###################################################
# Simulation to compare standard errors when estimating sparsity with Equation (8.15) and summary(..., se = "iid")

library(quantreg)

n <- 10000
R <- 500
p <- c(0.1, 0.5, 0.9)

res <- array(NA, dim = c(2,3,2,R))
for(i in 1:R){
x <- rbinom(n, 1, 0.5)
y <- 1 + 2*x + rnorm(n, 0, sd = 10)
fit <- rq(y ~ x, tau = p, method = "br")
xx <- cbind(1, x)
J <- matrix(rowMeans(apply(xx, 1, function(z) tcrossprod(z))), 2, 2)
xbar <- c(1,mean(x))

# Hall-Sheather bandwidth (hs = TRUE)
n <- length(y)
h <- quantreg::bandwidth.rq(p = p, n = n, hs = TRUE, alpha = 0.05)
# obtain dQ from two quantile models
bhat_low <- rq(y ~ x, tau = p - h, method = "br")$coef
bhat_up <- rq(y ~ x, tau = p + h, method = "br")$coef
dQ <- t(xbar) %*% (bhat_up - bhat_low)
# difference quotient
sparsity <- dQ/(2*h)
# calculate variance-covariance
V <- outer(solve(J), as.numeric(p*(1-p)*sparsity^2/n), "*")
res[,,1,i] <- apply(V, 3, function(x) sqrt(diag(x)))
res[,,2,i] <- sapply(summary(fit, se = "iid", covariance = TRUE), function(z) z$coef[,2])
}

apply(res[,1,1,], 1, mean)
apply(res[,1,2,], 1, mean)

apply(res[,2,1,], 1, mean)
apply(res[,2,2,], 1, mean)

apply(res[,3,1,], 1, mean)
apply(res[,3,2,], 1, mean)

###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\fev.Rdata")

###################################################
# Analysis

lsz <- 1

fit.rqs <- rq(fev ~ I(age - 10) + I(height - 61) + sex, tau = -1, data = fev)
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

gg <- list()

gg[[1]] <- ggplot(dd1, aes(x = x, y = y)) + geom_line(linewidth = 0.5) + xlab("") + ylab("") + xlim(xmin, 1) + ylim(3,3.5)

gg[[1]] <- ggplot(dd1, aes(x = x, y = y)) + geom_line(linewidth = lsz) + xlab("Quantile index") + ylab("Quantile effect") + ggtitle("")
gg[[2]] <- ggplot(dd2, aes(x = x, y = y)) + geom_line(linewidth = lsz) + xlab("Quantile index") + ylab("") + ggtitle("")

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")



###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\fev.Rdata")

###################################################
# Analysis

p <- c(0.1, 1:3/4, 0.9)
fit <- rq(fev ~ I(age - 10) + I(height - 61) + sex, tau = c(0.1, 1:3/4, 0.9), data = fev)
Qhat <- predict(fit)
colnames(Qhat) <- p
print(sel <- which(t(apply(Qhat, 1, diff)) < 0, arr.ind = TRUE))
fev[sel[1],]
Qhat[sel[1],]


###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\fev.Rdata")

###################################################
# Analysis

fit.rq <- rq(dist ~ speed, tau = 1:19/20, data = cars)
fit.rq
Qhat <- predict(fit.rq)
print(sel <- which(t(apply(Qhat, 1, diff)) < 0, arr.ind = TRUE))

kt <- KhmaladzeTest(dist ~ speed, tau = 1:19/20,
nullH = "location-scale", data = cars)
Qtools::KhmaladzeFormat(kt, epsilon = 0.05)
fit.rrq <- Qtools::rrq(dist ~ speed, tau = 1:19/20, data = cars)

transalpha <- 0.5
gg <- list()

b <- coef(fit.rq)
k <- ncol(b)
gg[[1]] <- ggplot(cars, aes(x = speed, y = dist)) + geom_point(alpha = transalpha) + xlab("Speed (mph)") + ylab("Stopping distance (ft)") + ggtitle("")
for(i in 1:k) gg[[1]] <- gg[[1]] + geom_abline(intercept = b[1,i], slope = b[2,i], linetype = "dashed")
gg[[1]] <- gg[[1]] + geom_abline(intercept = b[1,c(5,10,15)], slope = b[2,c(5,10,15)], linetype = "solid", linewidth = lsz)

b <- coef(fit.rrq)
k <- ncol(b)
gg[[2]] <- ggplot(cars, aes(x = speed, y = dist)) + geom_point(alpha = transalpha) + xlab("Speed (mph)") + ylab("Stopping distance (ft)") + ggtitle("")
for(i in 1:k) gg[[2]] <- gg[[2]] + geom_abline(intercept = b[1,i], slope = b[2,i], linetype = "dashed")
gg[[2]] <- gg[[2]] + geom_abline(intercept = b[1,c(5,10,15)], slope = b[2,c(5,10,15)], linetype = "solid", linewidth = lsz)

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")


