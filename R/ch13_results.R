# Statistical Methods and Applications Using Quantiles - Chapter 13

rm(list = ls())
library(quantreg)
library(Qtools)
library(aqmm)
library(nlqmm)
library(nlme)
library(reshape2)
library(tidyverse)
library(httr)
library(jsonlite)
library(ggplot2)
library(xtable)
library(ggrepel)

###############################################
###############################################

boot.data <- function(x, group, numeric = TRUE){

group.old <- group
	if(is.factor(group.old)){
		if(numeric){
			group <- as.numeric(levels(group))[group]
		} else {
			group <- as.character(levels(group))[group]
		}
	}
x <- x[order(group),]
group <- sort(group)
grp <- unique(group)
M <- length(grp)
n <- table(group)

id <- sort(sample(grp, M, replace = TRUE))
x$idnew <- group
dd <- do.call("rbind", lapply(id, function(i,x) subset(x, x$idnew == i), x = x))

nid <- rep(grp, n[id])
if(is.factor(group.old)){
	nid <- factor(nid, levels(group.old))
}
dd$newid <- nid
dd
}

###############################################
###############################################


###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\internet.Rdata")

tt <- table(internet$iso3)
table(tt)
sel <- names(tt)[tt > 29]
dds <- subset(internet, internet$iso3 %in% sel)
highlight <- c("CHN", "IND", "ITA", "USA")

###############################################
###############################################

tau <- c(.1, .5, .9)
nt <- length(tau)

internet$groupn <- 1*(internet$incomegr == "High-income countries")
dd1 <- subset(internet, internet$incomegr == "LMICs")
dd2 <- subset(internet, internet$incomegr == "High-income countries")

# get starting values
fit1 <- nls(percentage ~ SSlogis(year, Asym, xmid, scal), data = dd1)
summary(fit1)

fit2 <- nls(percentage ~ SSlogis(year, Asym, xmid, scal), data = dd2)
summary(fit2)

ss1 <- as.list(coef(fit1))
names(ss1) <- paste0("b", 1:3)
ss2 <- as.list(coef(fit2))
names(ss2) <- paste0("b", 1:3)
tmp <- as.numeric(ss2) - as.numeric(ss1)
ss3 <- rep(NA, 6)
ss3[c(1,3,5)] <- as.numeric(ss1)
ss3[c(2,4,6)] <- tmp
ss3 <- as.list(ss3)
names(ss3) <- paste0("b", 1:6) 

fit1.nlqmm <- fit2.nlqmm <- fit3.nlqmm <- list()
for(i in 1:nt){

# LMICs
fit1.nlqmm[[i]] <- nlqmm(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed = Asym + xmid + scal ~ 1, random = Asym + xmid ~ 1, group = iso3, data = dd1, start = ss1, tau = tau[i])

# High-income
fit2.nlqmm[[i]] <- nlqmm(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed =Asym + xmid + scal ~ 1, random = Asym  + xmid ~ 1, group = iso3, data = dd2, start = ss2, tau = tau[i])

# All (random Asym)
fit3.nlqmm[[i]] <- nlqmm(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed = Asym + xmid + scal ~ groupn, random = Asym ~ 1, group = iso3, data = dd, tau = tau[i], fit = TRUE, start = ss3)
}

fit1.nlme <- nlme(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed = Asym + xmid + scal ~ 1, random = pdDiag(Asym + xmid ~ 1), groups = ~ iso3, data = dd1, start = unlist(ss1))
fit2.nlme <- nlme(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed = Asym + xmid + scal ~ 1, random = pdDiag(Asym + xmid ~ 1), groups = ~ iso3, data = dd2, start = unlist(ss2))
fit3.nlme <- nlme(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed = Asym + xmid + scal ~ groupn, random = Asym ~ 1, groups = ~iso3, data = dd, start = unlist(ss3))


fit4.nlqmm <- list()

for(i in 1:nt){

# All (random Asym + xmid)
fit4.nlqmm[[i]] <- nlqmm(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed = Asym + xmid + scal ~ groupn, random = Asym + xmid ~ 1, covariance = "pdSymm", group = iso3, data = dd, tau = tau[i], fit = TRUE, start = ss3)

}

fit4.nlme <- nlme(model = percentage ~ Asym/(1+exp((xmid-year)/scal)), fixed = Asym + xmid + scal ~ groupn, random = pdSymm(form = Asym + xmid ~ 1), groups = ~iso3, data = dd, start = unlist(ss3))

save(list = c("fit1.nlqmm", "fit2.nlqmm", "fit3.nlqmm", "fit4.nlqmm", "fit1.nlme", "fit2.nlme", "fit3.nlme", "fit4.nlme"), file = "ch13_results.Rdata")

###############################################
###############################################

load(file = "ch13_results.Rdata")

###############################################
###############################################

# Bootstrap

set.seed(860)
M <- length(unique(dd$iso3))
npar <- length(c(fit4.nlqmm$theta_x, fit4.nlqmm$theta_z, fit4.nlqmm$sigma))
R <- 199
B <- array(NA, dim = c(R, npar))
ctrl <- lqmmControl(startQR = TRUE, method = "Nelder-Mead", verbose = FALSE, beta = 0.1)
dd$id <- factor(dd$iso3, levels(dd$iso3), labels = 1:M)
for(b in 1:R){
	print(b)
	ddb <- boot.data(dd, dd$id, numeric = TRUE)
	tmp <- try(update(fit4.nlqmm, data = ddb, control = ctrl, group = newid), silent = TRUE)
	if(class(tmp) != "try-error"){
		B[b,] <- c(tmp$theta_x, tmp$theta_z, tmp$sigma)
	}
}

dimnames(B) <- list(replicate = 1:R, c(paste0("b",1:6),paste0("csi",1:3),"sigma"))

save(list = c("fit1.nlqmm", "fit2.nlqmm", "fit3.nlqmm", "fit4.nlqmm", "fit1.nlme", "fit2.nlme", "fit3.nlme", "fit4.nlme", "B"), file = "ch13_results.Rdata")

###############################################
###############################################
