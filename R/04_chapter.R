# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 4
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
## Section 4.2
###################################################

## Figure 4.2

load(file = "../data/USbirths2001.Rdata")

ga <- USbirths2001$ga
ga <- ga[!is.na(ga)]
ga <- ga[ga > 24]

xx <- seq(min(ga), max(ga) + 3, length = 100)
dd <- data.frame(x = xx, y = ecdf(ga)(xx))

ggplot(dd, aes(x,y)) + geom_step(linetype = "solid", color = "black") + theme(legend.title = element_blank(), legend.position = "top", legend.spacing.x = unit(1.0, 'cm'), legend.key.width = unit(3, "line")) + xlab("Gestational age (weeks)") + ylab("ECDF")

###################################################
## Section 4.4
###################################################

## R Code 4.1

n <- c(10, 50, 100, 500, 1000, 5000, 10000)
R <- 5000
lambda <- 5

# CDF at jumps
round(ppois(0:5, lambda = lambda), 4)
# population median
print(xi <- qpois(0.5, lambda = lambda))
# CDF at the median
print(p <- ppois(xi, lambda = lambda))
# simulation for consistency
set.seed(456)
pn_jump <- rep(NA, length(n))
pn_step <- rep(NA, length(n))
for(j in 1:length(n)){
x <- replicate(R, rpois(n[j], lambda = lambda))
pn_jump[j] <- mean(apply(x, 2, quantile, probs = p, type = 1) == xi)
pn_step[j] <- mean(apply(x, 2, quantile, probs = 0.5, type = 1) == xi)
}

## Table 4.2

tab2 <- data.frame(rbind(pn_jump,pn_step))
colnames(tab2) <- n
rownames(tab2) <- c("Jump","Step")

addtorow <- list()
addtorow$pos <- list(-1,0,0,nrow(tab2))
addtorow$command <- c("\\toprule",
 "$n$ & \\multicolumn{1}{c}{$10$} & \\multicolumn{1}{c}{$50$} & \\multicolumn{1}{c}{$100$} & \\multicolumn{1}{c}{$500$} & \\multicolumn{1}{c}{$1000$} & \\multicolumn{1}{c}{$5000$} & \\multicolumn{1}{c}{$\\numprint{10000}$} \\\\",
"\\cmidrule(l){2-8}", "\\bottomrule"
)

print.xtable(xtable(tab2, digits = 2, align = c("l",rep("r",7)), caption = c("Consistency of the $p$th sample quantile when $p$ is at the jump or on the step.")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.colnames = FALSE, hline.after = NULL)

## R Code 4.2

library(Qtools)
set.seed(3874)
z <- rpois(1000, lambda = 5)
p <- ecdf(z)(sort(unique(z)))
Qhat <- quantile(z, probs = p, type = 1)
Qmid <- midquantile(z, probs = p)
tab <- data.frame(p, Qhat, Qmid$y)
colnames(tab) <- c("p","Q(p)","Qmid(p)")
all.equal(as.numeric(midq2q(Qmid)), as.numeric(Qhat))

## Table 4.3

tab3 <- data.frame(matrix(NA, 6, 8))
tab3[,1] <- rep(c("$p$", "$\\hat{Q}_{Y}(p)$", "$\\hat{\\mathcal{Q}}_{Y}(p)$"), 2)
tab3[1:3,2:8] <- t(tab)[,1:7]
tab3[4:6,2:7] <- t(tab)[,8:13]

dig <- matrix(3, 6, 9)
dig[c(2, 5),] <- 0

addtorow <- list()
addtorow$pos <- list(-1,3,nrow(tab3))
addtorow$command <- c("\\toprule","\\midrule","\\bottomrule")

print.xtable(xtable(tab3, digits = dig, align = c("l","l",rep("l",7)), caption = c("Quantiles and mid-quantiles of the Poisson data generated with R Code \\ref{code:04_02}.")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.rownames = FALSE, include.colnames = FALSE, hline.after = NULL, sanitize.text.function = identity)

###################################################
## Section 4.5
###################################################

## Simulation

# Normal

R <- 10000
types <- 1:10
nt <- length(types)
p <- 0:10/10
np <- length(p)
samp <- c(5, 10, 25, 50)
N <- length(samp)
res <- array(NA, dim = c(R, np, N, nt), dimnames = list(replication = 1:R, quantile = p, samplesize = samp, type = types))

set.seed(199)
for(n in 1:N){
	for(i in 1:R){
	x <- rnorm(samp[n], 0, 1)
		for(k in 1:(nt-1)){
			res[i,,n,k] <- quantile(x, type = k, probs = p)
		}
		res[i,,n,10] <- midquantile(x, probs = p)$y	
	}
}

true <- qnorm(p, 0, 1)
bias <- sds <- mse <- array(NA, dim = c(nt, np, N), dimnames = list(type = types, quantile = p, samplesize = samp))
for(n in 1:N){
bias[,,n] <- t((apply(res[,,n,], 3, function(x) apply(x, 2, mean)) - true))
sds[,,n] <- t(apply(res[,,n,], 3, function(x) apply(x, 2, sd)))
mse[,,n] <- bias[,,n]^2 + sds[,,n]^2
}

dd1 <- rbind(cbind(melt(bias[,-c(1,11),]), group = "Bias", distribution = "Normal"), cbind(melt(mse[,-c(1,11),]), group = "MSE", distribution = "Normal"))
dd1$type <- factor(dd1$type)
dd1$quantile <- factor(dd1$quantile)
dd1$samplesize <- factor(dd1$samplesize)
dd1$group <- factor(dd1$group)

# Poisson

R <- 10000
types <- 1:10
nt <- length(types)
p <- 0:10/10
np <- length(p)
samp <- c(5, 10, 25, 50)
N <- length(samp)
res2 <- array(NA, dim = c(R, np, N, nt), dimnames = list(replication = 1:R, quantile = p, samplesize = samp, type = types))

set.seed(199)
for(n in 1:N){
	for(i in 1:R){
	x <- rpois(samp[n], lambda = 1000)
		for(k in 1:(nt-1)){
			res2[i,,n,k] <- quantile(x, type = k, probs = p)
		}
		res2[i,,n,10] <- midquantile(x, probs = p)$y	
	}
}

true <- qpois(p, lambda = 1000)
bias <- sds <- mse <- array(NA, dim = c(nt, np, N), dimnames = list(type = types, quantile = p, samplesize = samp))
for(n in 1:N){
bias[,,n] <- t((apply(res2[,,n,], 3, function(x) apply(x, 2, mean)) - true))
sds[,,n] <- t(apply(res2[,,n,], 3, function(x) apply(x, 2, sd)))
mse[,,n] <- bias[,,n]^2 + sds[,,n]^2
}

dd2 <- rbind(cbind(melt(bias[,-c(1,11),]), group = "Bias", distribution = "Poisson"), cbind(melt(mse[,-c(1,11),]), group = "MSE", distribution = "Poisson"))
dd2$type <- factor(dd2$type)
dd2$quantile <- factor(dd2$quantile)
dd2$samplesize <- factor(dd2$samplesize)
dd2$group <- factor(dd2$group)


## Figures 4.5 and 4.6

f <- function(x) {
  r <- quantile(x, probs = c(0, 0.25, 0.5, 0.75, 1))
  names(r) <- c("ymin", "lower", "middle", "upper", "ymax")
  r
}

colg <- "grey70"
labs <- c('5' = "n = 5", '10' = "n = 10", '25' = "n = 25", '50' = "n = 50")

# normal
ggplot(dd1, aes(x = type, y = value)) + geom_boxplot(outlier.shape = NA) + stat_summary(fun.data = f, geom="boxplot") + facet_grid(group ~ samplesize, labeller = labeller(samplesize = labs), scales = "free") + theme(strip.background = element_rect(colour = colg, fill = colg), strip.text = element_text(face = "bold", size = 13)) + xlab("Quantile estimator type") + ylab("")

# Poisson
ggplot(dd2, aes(x = type, y = value)) + geom_boxplot(outlier.shape = NA) + stat_summary(fun.data = f, geom="boxplot") + facet_grid(group ~ samplesize, labeller = labeller(samplesize = labs), scales = "free") + theme(strip.background = element_rect(colour = colg, fill = colg), strip.text = element_text(face = "bold", size = 13)) + xlab("Quantile estimator type") + ylab("")
