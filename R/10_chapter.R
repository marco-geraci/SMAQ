# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 10
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
## Section 10.1
###################################################

## Figure 10.1

data("esterase", package = "Qtools")
bks <- (0:11)*100
ggplot(esterase, aes(x = Count)) + geom_histogram(aes(y = ..density..), breaks = bks) + scale_x_continuous(name = "Number of bindings", limits = c(0,1200)) + ylab("Density")

## Table 10.1

probs <- c(0.10, 0.50, 0.90)
y <- esterase$Count
x <- esterase$Esterase
logx <- log(esterase$Esterase, base = 2)
dd <- data.frame(x = x, logx = logx, y = y)

fit1 <- glm(y ~ x, data = dd, family = poisson())
fit2 <- glm(y ~ logx, data = dd, family = poisson())
fit3 <- glm.nb(y ~ logx, data = dd)

lambda <- fitted(fit1)
qhat1 <- sapply(probs, function(p) qpois(p, lambda = lambda))
nr1 <- apply(sweep(-qhat1, 1, y, "+") < 0, 2, mean)*100

lambda <- fitted(fit2)
qhat2 <- sapply(probs, function(p) qpois(p, lambda = lambda))
nr2 <- apply(sweep(-qhat2, 1, y, "+") < 0, 2, mean)*100

mu <- fitted(fit3)
qhat3 <- sapply(probs, function(p) qnbinom(p, size = fit3$theta, mu = mu))
nr3 <- apply(sweep(-qhat3, 1, y, "+") < 0, 2, mean)*100



nn <- c("Residual deviance (df)", "AIC", "Dispersion test p-value", "NB overdispersion", "PNR")

p1 <- round(AER::dispersiontest(fit1)[["p.value"]], 3)
p1 <- if(p1 < 0.001) "< 0.001"
res1 <- c(
paste0(round(deviance(fit1)), " (", fit1[['df.residual']], ")"),
round(AIC(fit1)),
p1,
"--",
paste(round(nr1), collapse = ", ")
)

p2 <- round(AER::dispersiontest(fit2)[["p.value"]], 3)
p2 <- if(p2 < 0.001) "< 0.001"
res2 <- c(
paste0(round(deviance(fit2)), " (", fit2[['df.residual']], ")"),
round(AIC(fit2)),
p2,
"--",
paste(round(nr2), collapse = ", ")
)

k <- as.numeric(fit3[['df.residual']])-1
res3 <- c(
paste0(round(deviance(fit3)), " (", k, ")"),
round(AIC(fit3)),
"--",
round(fit3$theta,1),
paste(round(nr3), collapse = ", ")
)

tab1 <- data.frame("Poisson" = res1, "Poisson (log-x)" = res2, "NB (log-x)" = res3)
rownames(tab1) <- nn

addtorow <- list()
# pos and command must have the same length
addtorow$pos <- list(-1,0,0,nrow(tab1))
addtorow$command <- c(
"\\toprule",
" Model & \\multicolumn{1}{c}{Poisson} & \\multicolumn{1}{c}{Poisson ($\\log_{2}x$)} & \\multicolumn{1}{c}{NB ($\\log_{2}x$)} \\\\",
"\\cmidrule(l){2-4}",
"\\bottomrule"
)

print.xtable(xtable(tab1, digits = c(0,0,0,2), align = c(rep("l",1),rep("r",3)), caption = c("Summary statistics for the generalised linear models fitted in the Esterase Assay dataset: residual deviance (degrees of freedom), Akaike information criterion (AIC), overdispersion test $p$-value, negative binomial (NB) overdispersion parameter, proportion of negative residuals (PNR, $\\%$) at three conditional quantiles ($p \\in\\{0.1,0.5,0.9\\}$).")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.colnames = FALSE, hline.after = NULL, format.args = list(big.mark = ",", decimal.mark = "."))

## Figure 10.2

library(reshape2)

pts1 <- data.frame(
  x = dd$x,
  y = dd$y,
  panel = "pois_x"
)

pts2 <- data.frame(
  x = dd$logx,
  y = dd$y,
  panel = "pois_logx"
)

pts3 <- data.frame(
  x = dd$logx,
  y = dd$y,
  panel = "nb_logx"
)

pts <- rbind(pts1, pts2, pts3)
pts$panel <- factor(
  pts$panel,
  levels = c("pois_x", "pois_logx", "nb_logx")
)

## fitted quantile curves: model 1
xx <- seq(min(dd$x), max(dd$x), length = 400)
lambda <- exp(cbind(1, xx) %*% coef(fit1))
qhat <- sapply(probs, function(p) qpois(p, lambda = lambda))

df1 <- melt(qhat)
df1$x <- rep(xx, length(probs))
df1$p <- factor(df1$Var2)
df1$panel <- "pois_x"

## fitted quantile curves: model 2
xx <- seq(min(dd$logx), max(dd$logx), length = 400)
lambda <- exp(cbind(1, xx) %*% coef(fit2))
qhat <- sapply(probs, function(p) qpois(p, lambda = lambda))

df2 <- melt(qhat)
df2$x <- rep(xx, length(probs))
df2$p <- factor(df2$Var2)
df2$panel <- "pois_logx"

## fitted quantile curves: model 3
xx <- seq(min(dd$logx), max(dd$logx), length = 400)
mu <- exp(cbind(1, xx) %*% coef(fit3))
qhat <- sapply(probs, function(p) qnbinom(p, size = fit3$theta, mu = mu))

df3 <- melt(qhat)
df3$x <- rep(xx, length(probs))
df3$p <- factor(df3$Var2)
df3$panel <- "nb_logx"

curves <- rbind(df1, df2, df3)
curves$panel <- factor(
  curves$panel,
  levels = c("pois_x", "pois_logx", "nb_logx")
)


psz <- 1
lsz <- 1
transalpha <- 0.4
colg <- "grey70"

g <- ggplot() +
  geom_point(
    data = pts,
    aes(x = x, y = y),
    shape = 16,
    size = psz,
    alpha = transalpha
  ) +
  geom_line(
    data = curves,
    aes(x = x, y = value, group = p, color = p),
    linewidth = lsz
  ) +
  facet_wrap(
    ~ panel,
    scales = "free_x",
    nrow = 1,
    labeller = as_labeller(c(
      pois_x    = "Poisson~(x)",
      pois_logx = "Poisson~(log[2]~x)",
      nb_logx   = "NB~(log[2]~x)"
    ), label_parsed)
  ) +
  coord_cartesian(ylim = c(0, 1350)) +
  scale_color_grey(start = 0.8, end = 0) +
  labs(x = "Esterase", y = "Number of bindings") +
  theme(legend.position = "none",
    strip.background = element_rect(colour = colg, fill = colg), strip.text = element_text(face = "bold", size = 13)
  )

g

## R Code 10.1

library(MASS)
data("esterase", package = "Qtools")
dd <- data.frame(logx = log(esterase$Esterase, base = 2),
y = esterase$Count)
fit.pois <- glm(y ~ logx, data = dd, family = poisson)
fit.pois
AER::dispersiontest(fit.pois)
summary(fit.nb <- glm.nb(y ~ logx, data = dd))

## R Code 10.2

# predicted quantiles of fitted negative binomial
probs <- c(0.10, 0.50, 0.90)
mu <- fitted(fit.nb)
qhat <- sapply(probs, function(p) qnbinom(p, size = fit.nb$theta, mu = mu))
# number of negative residuals
print(nr <- apply(sweep(-qhat, 1, dd$y, "+") < 0, 2, sum))
n <- nrow(esterase)
# p-value of test of the PNR at the median (probs[2])
binom.test(nr[2], n, p = probs[2])$p.value

###################################################
## Section 10.2
###################################################

## R Code 10.3

set.seed(456)
probs <- c(0.10, 0.50, 0.90)
fit.rq <- lapply(probs, function(p) rq.counts(y ~ logx, data = dd,
tau = p, M = 50))
names(fit.rq) <- probs
lapply(fit.rq, function(x) x$tTable)
# percentage of negative residuals
qhat <- sapply(fit.rq, function(x) predict(x, type = "response"))
apply(sweep(-qhat, 1, dd$y, "+") < 0, 2, mean)*100
# lack-of-fit test
lapply(fit.rq, LOFTest, alpha = 0.05, B = 1000, seed = 987)

## Figure 10.4

psz <- 2
lsz <- 1
transalpha <- 0.4

data("esterase", package = "Qtools")
y <- esterase$Count
x <- esterase$Esterase
logx <- log(esterase$Esterase, base = 2)
dd <- data.frame(x = x, logx = logx, y = y)
probs <- c(0.10, 0.50, 0.90)
set.seed(456)
fit4 <- lapply(probs, function(p) rq.counts(y ~ logx, data = dd, tau = p, M = 50))

xx <- seq(min(logx), max(logx), length = 400)
nd <- data.frame(logx = xx)
qhat <- sapply(fit4, function(x) predict(x, newdata = nd, type = "response"))
colnames(qhat) <- probs
df <- reshape2::melt(qhat)
df$x <- rep(xx, length(probs))
df$p <- factor(df$Var2)

g <- ggplot() + geom_point(data = dd, aes(x = logx, y = y),  shape = 16, size = psz, alpha = transalpha) + geom_line(data = df, aes(x = x, y = value, group = p, color = p), linewidth = lsz) + scale_color_grey(start = 0.8, end = 0) + theme(legend.title = element_blank(), legend.position = "none") + scale_x_continuous(expression(paste("Esterase (", log[2], ")"))) + ylab("Number of bindings") + ylim(c(0,1200))
g

## R Code 10.4

# find a model-implied p such that Q(p) is close to the NB mean
mu <- fitted(fit.nb)
pstar <- mean(pnbinom(floor(mu), size = fit.nb$theta, mu = mu))
pstar
# error when using mean as quantile pstar of NB
max(abs(mu - qnbinom(pstar, size = fit.nb$theta, mu = mu))/mu)
# fit quantile model for p = pstar
set.seed(123)
fit.rq_pstar <- rq.counts(y ~ logx, data = dd, tau = pstar)
# error when using the NB mean as an approximation to Q(pstar)
max(abs(mu - predict(fit.rq_pstar))/mu)
# compare standard errors of coefficients estimates
fit.rq_pstar$tTable[,1:2]
summary(fit.nb)$coef[,1:2]

###################################################
## Section 10.3
###################################################

## R Code 10.5

# generate data
set.seed(123)
n <- 200
x <- sample(0:3, n, replace = TRUE)
y <- floor(x) + floor(x + 1)*sample(1:6, n, replace = TRUE)
# unique values
print(z <- sort(unique(y)))
# fit a logistic regression for each threshold
X <- cbind(1,x)
b <- sapply(z, function(z) glm.fit(X, 1*(y <= z),
family = binomial())$coefficients)
# estimate CDF
F <- 1/(1 + exp(-X%*%b))
# estimate PMF
M <- t(diff(t(F)))
# estimate mid-CDF
G <- F - 0.5 * cbind(F[,1], M)

## R Code 10.6

library(Qtools)
set.seed(123)
n <- 1000
p <- c(1e-5, 0.15, 0.25, 0.5, 0.75, 0.85, 1-1e-5)
x <- runif(n, 0, 3)
y <- floor(1 + 2*x) + sample(1:6, n, replace = TRUE)
range(y)
fit <- midrq(y ~ x, tau = p, lambda = NULL, control = midrqControl(ecdf_est = "logit"))
r <- attr(fit$midFit$G, "range")
t(coef(fit))

## R Code 10.7

library(Qtools)
data("fars", package = "Qtools")
dd <- fars
# offset weekend/weekday
dd$offset <- log(ifelse(dd$weekend == 1, 2, 5))
# replace zero counts with small constant
dd$y <- ifelse(dd$count == 0, 0.1, dd$count)
# transform longitude and latitude
theta_c <- 37.5*pi/180
phi_c <- -96.5*pi/180
lat_rad <- dd$latitude*pi/180
lon_rad <- dd$longitude*pi/180
dd$lat_cos <- cos(lat_rad) - cos(theta_c)
dd$lat_sin <- sin(lat_rad) - sin(theta_c)
dd$lon_cos <- cos(lon_rad) - cos(phi_c)
dd$lon_sin <- sin(lon_rad) - sin(phi_c)

fit.midrq <- midrq(y ~ weekend + lat_cos + lat_sin + lon_cos + lon_sin, tau = c(0.75, 0.9), data = dd, offset = dd$offset, lambda = 0, control = midrqControl(ecdf_est = "logit"))

summary(fit.midrq, numerical = TRUE)
# mid-probability range
attr(fit.midrq$midFit$G, "range")
# weekend effect
100*(exp(coef(fit.midrq)["weekend",])-1)

## Figure 10.6

dd$wcount <- dd$count/exp(dd$offset)*7
dd$group <- factor(dd$weekend, levels = 0:1, labels = c("weekdays", "weekend days"))

gg <- list()

gg[[1]] <- ggplot(dd,aes(x = group, y = wcount, fill = group)) + geom_boxplot(outlier.alpha = 0.25) + scale_fill_manual(values = c("white", "grey70"), name = NULL,) + xlab(NULL) + ylab("Number of fatal crashes") + scale_y_sqrt("Fatal crashes per week", breaks = c(0, 10, 50, 100, 250, 500, 750), labels = function(x) ifelse(abs(x - eps) < 1e-12, "0", scales::label_number(accuracy = 1)(x))) + theme(legend.position = "none", panel.grid.minor = element_blank())

gg[[2]] <- ggplot(dd, aes(x = wcount, linetype = group)) + stat_ecdf(geom = "step", linewidth = 0.6) + ylab("ECDF") + labs(linetype = "") + scale_x_sqrt("Fatal crashes per week", breaks = c(0, 10, 50, 100, 250, 500, 750), labels = function(x) ifelse(abs(x - eps) < 1e-12, "0", scales::label_number(accuracy = 1)(x))) + theme(legend.position = "top", legend.title = element_text(size = 14), panel.grid.minor = element_blank())

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")

## Figure 10.7

library(usmap)
qhat <- predict(fit.midrq)
global_min <- min(qhat)
global_max <- max(qhat)

dd$pred <- qhat[,2]
dd1 <- subset(dd, dd$weekend == 0)
dd2 <- subset(dd, dd$weekend == 1)
center_point <- data.frame(lon = -96.5, lat = 37.5, label = "Geographic reference")
center_proj <- usmap_transform(center_point)

my_gradient <- scale_fill_gradient(low = "grey100", high = "black", name = "", na.value = "white", limits = c(global_min, global_max))

plot_usmap(regions = "counties", values = "pred", data = dd1, exclude = c("Alaska", "Hawaii")) + geom_sf(data = center_proj, color = "red", size = 3, shape = 17, fill = "red") + theme(panel.background=element_blank(), legend.position = "right") + my_gradient + theme_void()


## Figure 10.8

p <- 2:9/10
fit.p <- midrq(y ~ weekend + lat_cos + lat_sin + lon_cos + lon_sin, tau = p, data = dd, offset = dd$offset, lambda = 0, control = midrqControl(ecdf_est = "logit"))
se <- sapply(vcov(fit.p, numerical = TRUE), function(x) sqrt(diag(x)[2]))
b <- coef(fit.p)[2,]
plot(p, b, xlim = c(0,1), ylim = c(-1,0.4), type = "b")
lines(p, b - qnorm(0.975)*se, col = 2)
lines(p, b + qnorm(0.975)*se, col = 2)
wend <- data.frame(p = p, est = b, se = se, lower = b - qnorm(0.975)*se, upper = b + qnorm(0.975)*se)

psz <- 2
lsz <- 1

g <- ggplot(wend, aes(x = p, y = est)) + geom_line(linewidth = lsz) + geom_point(shape = 16, size = psz) + geom_line(aes(x = p, y = lower), linetype = "dashed", linewidth = lsz) + geom_line(aes(x = p, y = upper), linetype = "dashed", linewidth = lsz) + scale_x_continuous(name = "Quantile index", breaks = 0:10/10, limits = c(0,1)) + ylab("Regression quantile")
g