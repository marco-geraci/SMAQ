# Statistical Methods and Applications Using Quantiles - Chapter 7

###################################################
# Packages and functions

rm(list = ls())
library(ggplot2)
library(Qtools)
library(quantreg)
library(boot)
library(foreach)
library(doSNOW)
library(parallel)

boot.fn <- function(x, i, p) quantile(x[i], probs = p, type = 1)

density.fn <- function(x, bw, kernel, probs, from, to){
	xi <- quantile(x = x, probs = probs, names = FALSE)
	if(missing(from) & missing(to)){ 
	out <- density(x = x, bw = bw, kernel = kernel)
	}
	if(!missing(from) & missing(to)){ 
	out <- density(x = x, bw = bw, kernel = kernel, from = from)
	}
	if(missing(from) & !missing(to)){ 
	out <- density(x = x, bw = bw, kernel = kernel, to = to)
	}
	if(!missing(from) & !missing(to)){ 
	out <- density(x = x, bw = bw, kernel = kernel, from = from, to = to)
	}
	approx(x = out$x, y = out$y, xout = xi)$y
}

gn <- function(x, k, rate){
	factorial((2*k+1))/factorial(k)^2*pexp(x, rate = rate)^(k)*(1 - pexp(x, rate = rate))^(k)*dexp(x, rate = rate)
}

gn0 <- function(x, p, k, rate){
	p - integrate(f = gn, lower = 0, upper = x, k = k, rate = rate)$val
}

path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results"

###################################################
# Simulation section 7.1.2

set.seed(123)

R <- 10000
n <- c(100, 1000, 5000)
p <- c(0.1, 0.5)
alpha <- 0.05
bw <- c("nrd0", "SJ", "ucv")
Ks <- c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine", "optcosine")
dist <- c("norm", "chisq", "beta")

out_sim_01 <- expand.grid(bw = bw, K = Ks, p = p, dist = dist, n = n, stringsAsFactors = FALSE)
ns <- nrow(out_sim_01)

cl <- makeCluster(detectCores() - 1)
registerDoSNOW(cl)
out_sim_01$coverage <- foreach(s = 1:ns, .combine = "c", .export = c("out_sim_01", "density.fn")) %dopar% {

library(Qtools)

sel <- out_sim_01[s,]

theta <- switch(sel$dist,
		"norm" = c(mean = 0, sd = 1),
		"chisq" = c(df = 3),
		"beta" = c(shape1 = 0.5, shape2 = 0.5))

r <- switch(sel$dist,
		"chisq" = c(from = 0),
		"beta" = c(from = logit(0.001), to = logit(0.999)))

args <- c(list(bw = sel$bw, kernel = sel$K, probs = sel$p), as.list(r))

Y <- replicate(R, do.call(paste0("r", sel$dist), args = c(list(n = sel$n), as.list(theta))))
if(sel$dist == "beta") Y <- apply(Y, 2, Qtools::logit)

Q <- do.call(paste0("q", sel$dist), args = c(list(p = sel$p), as.list(theta)))
if(sel$dist == "beta") Q <- logit(Q)

xi <- apply(Y, 2, quantile, probs = sel$p)

fhat <- apply(Y, 2, function(y, args) do.call("density.fn", args = c(list(x = y), args)), args = args)

se <- sqrt(sel$p*(1-sel$p)/sel$n)*1/fhat
ci <- sweep(outer(qnorm(1-alpha/2, 0, se), c(-1,1), FUN = "*"), 1, xi - Q)
out <- mean((ci[,1]*ci[,2]) < 0, na.rm = TRUE)
return(out)
}
stopCluster(cl)

###################################################
# Simulation section 7.1.4

set.seed(321)
p <- 0.5
alpha <- 0.05
zalpha <- qnorm(c(alpha/2, 1-alpha/2), 0, 1)

N <- c(5, 11, 51, 101)
R <- 10000
b <- 1
xi <- qexp(p, rate = b)
fxi <- dexp(xi, rate = b)

out_sim_02 <- array(NA, dim = c(R, length(N), 7), dimnames = list(rep = 1:R, sample_size = N, type = c("ci_asymptotic", "ci_boot_norm", "ci_boot_perc", "width_asymptotic", "width_boot_norm", "width_boot_perc", "width_true")))

bias <- array(NA, dim = c(R, length(N), 2), dimnames = list(rep = 1:R, sample_size = N, type = c("sample_est", "boot_sample")))
SE <- array(NA, dim = c(R, length(N), 2), dimnames = list(rep = 1:R, sample_size = N, type = c("mean_sparsity", "sd")))

for(j in seq_along(N)){
n <- N[j]
k <- (n - 1)/2
Y <- replicate(R, rexp(n, rate = b))
xi_n <- apply(Y, 2, quantile, probs = p)
xi_n0 <- sqrt(n)*(xi_n - xi)
Qn <- quantile(xi_n0, probs = c(alpha/2, 1-alpha/2))

q1 <- uniroot(gn0, p = alpha/2, k = k, rate = b, interval = c(0, 4))$root
q2 <- uniroot(gn0, p = 1-alpha/2, k = k, rate = b, interval = c(0, 4))$root
Qsamp <- c(q1, q2) - xi # centered

for(i in 1:R){
	y <- Y[,i]
	B <- boot(y, statistic = boot.fn, R = 900, p = p)
	xi_n <- B$t0
	xi_star <- B$t
	xi_star0 <- sqrt(n)*(xi_star - xi_n)
	
	bias[i,j,1] <- xi_n
	bias[i,j,2] <- mean(xi_star)
	
	fhat <- density(y, bw = "SJ", kernel = "gaussian")
	f0 <- approx(fhat$x, fhat$y, xout = xi_n)
	SE[i,j,1] <- bspar <- sqrt(p*(1-p)/n)/f0$y
	SE[i,j,2] <- bsd <- sd(B$t)

	## asymptotic CI evaluated at true density
	#Qinf <- xi_n + zalpha*sqrt(p*(1-p))/fxi - xi # centered
	# asymptotic CI
	Qhatn <- xi_n + zalpha*bspar - xi # centered
	## normal CI (with or without bias correction)
	#Qnorm <- boot.ci(B, conf = 1-alpha, type = "norm")$normal[2:3] - xi  # centered, bias corrected
	Qnorm <- xi_n + zalpha*bsd - xi  # centered
	## percentile CI (backwards)
	Qperc <- boot.ci(B, conf = 1-alpha, type = "perc")$percent[4:5] - xi  # centered
	#Qstar <- quantile(xi_star0, probs = c(alpha/2, 1-alpha/2)) # approx = to boot.ci percentile method
	
	out_sim_02[i,j,1] <- prod(Qhatn) < 0
	out_sim_02[i,j,2] <- prod(Qnorm) < 0
	out_sim_02[i,j,3] <- prod(Qperc) < 0
	#out_sim_02[i,j,4] <- prod(Qsamp) < 0

	#yy <- seq(min(xi_star0), max(xi_star0), length = 300)	
	#Gn <- pnorm(yy, 0, sqrt(n)*bspar)
	#Gn <- pnorm(yy, 0, sqrt(p*(1-p))/fxi)
	out_sim_02[i,j,4] <- diff(Qhatn)
	out_sim_02[i,j,5] <- diff(Qnorm)
	out_sim_02[i,j,6] <- diff(Qperc)
	out_sim_02[i,j,7] <- diff(Qsamp)
}
}

apply(out_sim_02, 3, colMeans)
apply(bias, 3, colMeans) - xi
apply(SE, 3, colMeans)

###################################################
# Save output
save(list = c("out_sim_01", "out_sim_02", "bias", "SE"), file = file.path(path.out, "ch7_results.Rdata"))
