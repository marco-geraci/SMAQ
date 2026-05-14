# Statistical Methods and Applications Using Quantiles - Chapter 4

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(xtable)
library(Qtools)
library(reshape2)

path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results"

###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")

###################################################
# Gestational age

ga <- USbirths2001$ga
ga <- ga[!is.na(ga)]
ga <- ga[ga > 24]

###################################################
# Simulation

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

###################################################
# Save output
save(list = c("ga", "dd1", "dd2"), file = file.path(path.out, "ch4_results.Rdata"))

