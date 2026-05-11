# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 11
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
## Section 11.1
###################################################

## R Code 11.1

p <- 0.1

# unconstrained
fit <- rq(fev ~ sex, tau = c(p, 1-p), data = fev)
coef(fit)
print(delta <- apply(coef(fit), 1, diff))

## R Code 11.2

# Asymptotic SE (nid)
xx <- crossprod(fit$x)
fit.s  <- summary(fit,  se = "nid", covariance = TRUE)

H1 <- fit.s[[1]]$Hinv
H2 <- fit.s[[2]]$Hinv

S12 <- (min(p,1-p) - p*(1-p)) * H1 %*% xx %*% t(H2)

Sigma <- rbind(
  cbind(fit.s[[1]]$cov,S12),
  cbind(t(S12),fit.s[[2]]$cov)
)

C <- cbind(-diag(2), diag(2))
V <- C %*% Sigma %*% t(C)
print(SEa <- sqrt(diag(V)))

# Bootstrap SE
library(boot)
delta_fun <- function(data, indices, p) {
  d <- data[indices, ]
  b1 <- coef(rq(fev ~ sex, tau = p, data = d))
  b2 <- coef(rq(fev ~ sex, tau = 1 - p, data = d))
  return(b1 - b2)
}
set.seed(123)
B <- boot::boot(data = fev, statistic = delta_fun, R = 1000, p = 0.1)
print(SEb <- apply(B$t, 2, sd))

## R Code 11.3

data(NHANES, package = "NHANES")
dds <- NHANES[,c("Age", "Gender", "BPSysAve")]
dds <- subset(dds, complete.cases(dds))
range(dds$Age)
# > [1]  8 80
dds$Age.c <- dds$Age - min(dds$Age)
p <- 0.05
fit <- qrr(BPSysAve ~ Age.c*Gender, data = dds, tau = c(p,1-p))
fit
exp(coef(fit))
# age effect among males
exp(coef(fit)['Age.c'] + coef(fit)['Age.c:Gendermale'])

set.seed(123)
summary(fit, se = "boot", R = 200)

## Figure 11.1

colg <- "grey70"
psz <- 2
transalpha <- 0.4

data(NHANES, package = "NHANES")
dds <- NHANES[,c("Age", "Gender", "BPSysAve")]
dds <- subset(dds, complete.cases(dds))
dds$Gender <- factor(dds$Gender, levels = c("female", "male"), labels = c("Female", "Male"))

p <- c(.05, .5, .95)
fit1 <- tsrq(BPSysAve ~ Age, tsf = "bc", data = dds, subset = dds$Gender == "Female", tau = p, lambda = seq(0, 1, by = 0.1))
fit2 <- tsrq(BPSysAve ~ Age, tsf = "bc", data = dds, subset = dds$Gender == "Male", tau = p, lambda = seq(0, 8, by = 1))

xx <- seq(8, 80, by = 0.5)
res <- array(NA, dim = c(length(xx), 3, 2), dimnames = list(NULL, p = p, Gender = c("Female", "Male")))
res[,,1] <- predict(fit1, newdata = data.frame(Age = xx))
res[,,2] <- predict(fit2, newdata = data.frame(Age = xx))

dd <- melt(res)
dd$Age <- rep(xx, 6)
dd$BPSysAve <- dd$value
dd$p <- factor(dd$p)

g <- ggplot() + geom_point(data = dds, aes(x = Age, y = BPSysAve), alpha = transalpha, colour = colg, size = psz) + geom_line(data = dd, aes(x = Age, y = BPSysAve, group = p), linewidth = 1.2) + facet_wrap(~Gender) + xlab("Age (years)") + ylab("Systolic blood pressure (mmHg)") + theme(strip.background = element_rect(colour = colg, fill = colg), strip.text = element_text(face = "bold", size = 13))
g

###################################################
## Section 11.2
###################################################

## R Code 11.4

fit1 <- qlss(formula = BPSysAve ~ Age, data = NHANES,
subset = NHANES$Gender == "female", tsf = "bc", probs = 0.05)
fit2 <- qlss(formula = BPSysAve ~ Age, data = NHANES,
subset = NHANES$Gender == "male", tsf = "bc", probs = 0.05)
fit1
fit2

## Figure 11.2

xx <- seq(8, 80, by = 0.5)
fitp1 <- predict(fit1, newdata = data.frame(Age = xx))
fitp2 <- predict(fit2, newdata = data.frame(Age = xx))

y <- c(fitp1$location$median, fitp2$location$median, fitp1$scale$IQR, fitp2$scale$IQR, fitp1$shape$skewness, fitp2$shape$skewness, fitp1$shape$shape, fitp2$shape$shape)
x <- rep(xx, 8)
stat <- factor(rep(c("Median","IQR","Skewness","Tail weight"), each = 2*length(xx)), levels = c("Median","IQR","Skewness","Tail weight"))
gender <- factor(rep(rep(c("Female", "Male"), each = length(xx)), 4), levels = c("Female", "Male"))

dd <- data.frame(x = x, y = y, stat = stat, gender = gender)

colg <- "grey70"

g <- ggplot(dd, aes(x = x, y = y)) + geom_line() + facet_grid(stat ~ gender, labeller = labeller(stat = c(Median = "Median (mmHg)", IQR = "IQR (mmHg)", Skewness = "Skewness",'Tail weight' = "Tail weight")), scales = "free") + xlab("Age (years)") + ylab(NULL) + theme( strip.background = element_blank(), strip.text = element_text(face = "bold", size = 11))
g

## Code fragment

set.seed(123)
x <- seq(8, 80, by = 0.5)
fit1p <- predict(fit1, newdata = data.frame(Age = x),
  interval = TRUE, level = 0.90, R = 500)
plot(fit1p, z = x, whichp = 0.05, interval = TRUE)

###################################################
## Section 11.3
###################################################

## Figure 11.3

library(dplyr)
library(tidyr)

data(NHANES, package = "NHANES")

dds <- NHANES[,c("Age","Gender","BPSysAve")]
dds <- subset(dds, complete.cases(dds))

dds$Gender <- factor(dds$Gender)

p_grid <- seq(.05,.95,.01)
ages <- c(20,40,60)

fit_f <- rq(BPSysAve ~ Age, tau = p_grid, data = subset(dds, Gender=="female"))
fit_m <- rq(BPSysAve ~ Age, tau = p_grid, data = subset(dds, Gender=="male"))

results <- list()

for(a in ages){

  qf <- as.numeric(predict(fit_f, newdata=data.frame(Age=a)))
  qm <- as.numeric(predict(fit_m, newdata=data.frame(Age=a)))

  # quantile difference
  delta <- as.numeric(qm - qf)

  # approximate comparison distribution
  upsilon <- approx(qm, p_grid, xout=qf, rule=2)$y

  results[[as.character(a)]] <- data.frame(
      p=p_grid,
      Female=qf,
      Male=qm,
      delta=delta,
      upsilon=upsilon,
      age=a
  )
}

res <- bind_rows(results)


plotdata <- res %>%
pivot_longer(c(Female,Male),
             names_to="Sex",
             values_to="sbp")

gg <- list()

colg <- "grey70"
lab_age <- function(x) paste("Age ", c(20, 40, 60))

gg[[1]] <- ggplot(plotdata, aes(p, sbp, linetype = Sex)) + geom_line(linewidth = 1.1, colour = "black") + scale_linetype_manual(values = c("solid", "dashed")) + facet_wrap(~age, labeller = labeller(age = lab_age)) + labs(x = "p", y = "Systolic blood pressure (mmHg)", linetype = "Sex") + theme_bw() + theme(legend.position = "top", strip.background = element_rect(colour = colg, fill = colg), strip.text = element_text(face = "bold", size = 13), legend.spacing.x = unit(0.7, 'cm'), legend.key.width = unit(2, "line"))

gg[[2]] <- ggplot(res, aes(p, delta, linetype = factor(age))) + geom_line(linewidth = 1.1, colour = "black") + scale_linetype_manual(values = c("solid", "dashed", "dotted")) + geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.6, colour = colg) + labs(x = "p", y = "Quantile difference", linetype = "Age") + theme(legend.position = "top", legend.spacing.x = unit(0.7, 'cm'), legend.key.width = unit(2, "line"))

gg[[3]] <- ggplot(res, aes(p, upsilon, linetype = factor(age))) + geom_line(linewidth = 1.1, colour = "black") + scale_linetype_manual(values = c("solid", "dashed", "dotted")) + geom_abline(intercept = 0, slope = 1, linetype = "solid", linewidth = 0.6, colour = colg) + labs(x = "p", y = "Comparison distribution", linetype = "Age") + theme(legend.position = "top", legend.spacing.x = unit(0.7, 'cm'), legend.key.width = unit(2, "line"))

lay <- rbind(c(1, 1), c(2, 3))
g1 <- arrangeGrob(grobs = gg[1:3], layout_matrix = lay, heights = c(2, 2))
marrangeGrob(list(g1), nrow = 1, ncol = 1, top = "")

## Figure 11.4

library(metR)
p_grid <- seq(.05,.95,.01)
age_grid <- seq(18,80,2)

heat <- expand.grid(p=p_grid, age=age_grid)
heat$upsilon <- NA

for(i in seq_along(age_grid)){
  a <- age_grid[i]
  qf <- predict(fit_f, newdata=data.frame(Age=a))
  qm <- predict(fit_m, newdata=data.frame(Age=a))
  upsilon <- approx(qm, p_grid, xout=qf, rule=2)$y
  heat$upsilon[heat$age==a] <- upsilon
}

heat$diff <- heat$upsilon - heat$p
heat <- subset(heat, heat$p > 0.05) # remove 0 values at p = 0.05

g <- ggplot(heat, aes(x = p, y = age, z = diff)) +
  geom_tile(aes(fill = diff)) +
  geom_contour(colour = "black", linewidth = 0.6) +
  geom_text_contour(
    colour = "black",
    skip = 0,
    nudge_x = -0.01,
    nudge_y = 3,
    rotate = FALSE
  ) +
  scale_fill_gradient2(
    low = grey(0.2),
    mid = "white",
    high = "red",
    midpoint = 0,
    name = expression(upsilon(p, a) - p)
  ) + labs(x = "p", y = "Age (years)"
 ) + theme_classic()

g

## Code fragment

p_grid <- seq(.05, .95, .01)
dds <- subset(NHANES, Gender == "female")
fit <- rq(BPSysAve ~ Age, tau = p_grid, data = dds)
# sparsity via difference quotient of the conditional quantile function
n <- length(fit$y)
h <- bandwidth.rq(p_grid, n, hs = TRUE)
if(any(p_grid + h > 1) | any(p_grid - h < 0)) warning("h is too large")
bhi <- mapply(function(p, x, y) rq.fit.fnb(x, y, p)$coef, p = p_grid + h,
  MoreArgs = list(x = fit$x, y = fit$y))
blo <- mapply(function(p, x, y) rq.fit.fnb(x, y, p)$coef, p = p_grid - h,
  MoreArgs = list(x = fit$x, y = fit$y))
x0 <- 40 # predict at age 40
dQ <- as.numeric(matrix(c(1,x0), nrow = 1) %*% (bhi - blo))
q <- pmax((dQ - 1e-11)/(2 * h), 1e-6)
# derivative of comparison function
deriv <- -coef(fit)[2,]/q

## Figure 11.5

dds <- NHANES[, c("Age", "Gender", "BPSysAve")]
dds <- subset(dds, complete.cases(dds))

p_grid <- seq(.05, .95, .01)
x0 <- 40

compute_derivative <- function(fit, x0, label){
  n <- length(fit$y)
  h <- bandwidth.rq(p_grid, n, hs = TRUE)
  if(any(p_grid + h > 1) | any(p_grid - h < 0)) warning("h is too large")
  bhi <- mapply(function(p, x, y) rq.fit.fnb(x, y, p)$coef, p = p_grid + h, MoreArgs = list(x = fit$x, y = fit$y))
  blo <- mapply(function(p, x, y) rq.fit.fnb(x, y, p)$coef, p = p_grid - h, MoreArgs = list(x = fit$x, y = fit$y))

  dQ <- as.numeric(matrix(c(1,x0), nrow = 1) %*% (bhi - blo))
  q <- pmax((dQ - 1e-11)/(2 * h), 1e-6)
  deriv <- -coef(fit)[2,]/q

  data.frame(
    p = p_grid,
	sparsity = q,
    derivative = deriv,
    group = label
  )
}

df_f <- compute_derivative(fit_f, x0 = x0, "Female")
df_m <- compute_derivative(fit_m, x0 = x0, "Male")

df_plot <- bind_rows(df_f, df_m)
names(df_plot)[4] <- "Sex"

ggplot(subset(df_plot, p >= 0.08 & p <= 0.92), aes(x = p, y = derivative, linetype = Sex)) + geom_smooth(se = FALSE, span = 0.2, linewidth = 1.2, color = "black") + ylim(c(-0.02,0)) + labs(x = "p", y = "Derivative of comparison distribution")

###################################################
## Section 11.4
###################################################

## R Code 11.5

load(file = "../data/ipums.Rdata")

pL <- seq(0.01, 0.10, by = 0.01)
pU <- seq(0.90, 0.99, by = 0.01)
dpL <- diff(c(0, pL))
dpU <- diff(c(pU, 1))

## Linear quantile model by sex
fit1 <- rq(incwage ~ sex, tau = c(0.1,0.9), data = ipums, method = "fn")
fitL <- rq(incwage ~ sex, tau = pL, data = ipums, method = "fn")
fitU <- rq(incwage ~ sex, tau = pU, data = ipums, method = "fn")

## Predicted conditional quantiles on original scale
res1 <- data.frame(sex = factor(c("Female","Male"),
levels = levels(ipums$sex)))
Qhat <- predict(fit1, newdata = res1)
res1$Q10 <- Qhat[,1]
res1$Q90 <- Qhat[,2]
res1$Quantile_ratio <- res1$Q90/res1$Q10
res1$SE_ratio <- sqrt(qr_var(fit1, res1, se = "nid"))

## Tail means using a Riemann approximation
Qhat <- predict(fitL, newdata = res1)
res1$lower_mean <- rowSums(Qhat*dpL)/sum(dpL)
Qhat <- predict(fitU, newdata = res1)
res1$upper_mean <- rowSums(Qhat*dpU)/sum(dpU)
res1$Tail_ratio <- res1$upper_mean/res1$lower_mean
res1$SE_Tratio <- sqrt(qtr_var(fitL, fitU, newdata = res1,
se = "iid", wL = dpL, wU = dpU))

## Code fragment

library(reldist)
tapply(ipums$incwage, ipums$sex, gini)

## R Code 11.6

## Linear quantile model by age, sex and their interaction
fit2 <- rq(incwage ~ age*sex, tau = c(0.1,0.9), data = ipums, method = "fn")
fitL <- rq(incwage ~ age*sex, tau = pL, data = ipums, method = "fn")
fitU <- rq(incwage ~ age*sex, tau = pU, data = ipums, method = "fn")

## Predicted conditional quantiles on original scale
res2 <- expand.grid(
  age = c(21, 30, 40, 50, 65),
  sex = factor(c("Female","Male"), levels = levels(ipums$sex))
)
Qhat <- predict(fit2, newdata = res2)
res2$Q10 <- Qhat[,1]
res2$Q90 <- Qhat[,2]
res2$Quantile_ratio <- res2$Q90/res2$Q10
res2$SE_ratio <- sqrt(qr_var(fit2, res2, se = "nid"))

## Tail means using a Riemann approximation
Qhat <- predict(fitL, newdata = res2)
res2$lower_mean <- rowSums(Qhat*dpL)/sum(dpL)
Qhat <- predict(fitU, newdata = res2)
res2$upper_mean <- rowSums(Qhat*dpU)/sum(dpU)
res2$Tail_ratio <- res2$upper_mean/res2$lower_mean
res2$SE_Tratio <- sqrt(qtr_var(fitL, fitU, newdata = res2,
se = "iid", wL = dpL, wU = dpU))

## Table 11.1

res1 <- data.frame(age = c("21-65", "21-65"), res1)

tab <- rbind(res1, res2)
tab <- tab[,c(2,1,3:10)]
tab[c(4:7,9:12),1] <- ""

addtorow <- list()
addtorow$pos <- list(-1,0,0,0,0,2,2,7,7,nrow(tab))
addtorow$command <- c("\\toprule",
" Sex & Age & \\multicolumn{1}{c}{Q(10)} & \\multicolumn{1}{c}{Q(90)} & \\multicolumn{1}{c}{Quantile} & \\multicolumn{1}{c}{SE} & \\multicolumn{1}{c}{$\\mu_{L}$} & \\multicolumn{1}{c}{$\\mu_{U}$} & \\multicolumn{1}{c}{Tail ratio} & \\multicolumn{1}{c}{SE}\\\\",
" &  &  &  & \\multicolumn{1}{c}{ratio} &  &  & &  & \\\\",
"\\cmidrule(l){3-6}","\\cmidrule(l){7-10}","\\cmidrule(l){3-6}","\\cmidrule(l){7-10}", "\\cmidrule(l){3-6}","\\cmidrule(l){7-10}","\\bottomrule"
)

print.xtable(xtable(tab, digits = c(0,0,0,0,0,2,2,0,0,2,2), align = c(rep("l",3),rep("r",8)), caption = c("Tail quantiles, 90:10 quantile ratio, lower- and upper-tail means ($\\mu_{L}$ and $\\mu_{U}$, respectively), and 90:10 tail ratio of income by sex and by sex and age, with standard error (SE) of the ratios, using the IPUMS income data, United States, 2024.")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.rownames = FALSE, include.colnames = FALSE, hline.after = NULL, format.args = list(big.mark = ",", decimal.mark = "."))
