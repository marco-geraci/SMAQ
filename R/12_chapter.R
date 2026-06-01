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

## Figure 12.1

load(file = "../data/USbirths2001.Rdata")

# singletons born at term only
dd <- subset(USbirths2001, USbirths2001$ga >= 37 & USbirths2001$multiple == 0)
dd_male <- dd %>% filter(sex == "male")
dd_female <- dd %>% filter(sex == "female")

# random sample for figure
set.seed(123)
n <- 10000
USbirths2001sample <- dd[sample(seq_along(dd$bw), n),]
USbirths2001sample <- USbirths2001sample[,c("bw","sex")]

# LBW
Q <- 2500
USbirths2001sample$index <- seq_along(USbirths2001sample$bw)
USbirths2001sample$z <- 1*(USbirths2001sample$bw < Q)
dds10 <- subset(USbirths2001sample, z == 1)
dds_male <- USbirths2001sample %>% filter(sex == "male")
dds_female <- USbirths2001sample %>% filter(sex == "female")
dds10_male <- dds10 %>% filter(sex == "male")
dds10_female <- dds10 %>% filter(sex == "female")

dens_m <- density(dds_male$bw)
df_m <- data.frame(bw = dens_m$x, density = dens_m$y)
df_m_tail <- df_m %>% filter(bw < Q)
df_m_main <- df_m %>% filter(bw >= Q)
dens_f <- density(dds_female$bw)
df_f <- data.frame(bw = dens_f$x, density = dens_f$y)
df_f_tail <- df_f %>% filter(bw < Q)
df_f_main <- df_f %>% filter(bw >= Q)

pct <- mean(dd_male$bw < Q) * 100 # LBW prevalence for entire sample
pct_label <- paste0(sprintf("%.2f", pct), "%")
label_x <- 2000
label_y <- -max(df_m$density[df_m$bw < Q]) * 4.5

psz <- 2
lsz <- 1
transalpha <- 0.1
colg <- "grey70"
yl <- c(900,6000)

p_left <- ggplot() + geom_area(data = df_m_main, aes(x = bw, y = -density), fill = "grey70", alpha = 0.6) + geom_area(data = df_m_tail, aes(x = bw, y = -density), fill = "black", alpha = 0.8) + geom_vline(xintercept = Q, linetype = "dashed", linewidth = lsz) + coord_flip() + theme_void() + theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) + xlim(yl) + annotate("text", x = label_x, y = label_y, label = pct_label, color = "black", size = 3)

pct <- mean(dd_female$bw < Q) * 100 # LBW prevalence for entire sample
pct_label <- paste0(sprintf("%.2f", pct), "%")
label_x <- 2000
label_y <- max(df_f$density[df_f$bw < Q]) * 4.5

p_right <- ggplot() + geom_area(data = df_f_main, aes(x = bw, y = density), fill = "grey70", alpha = 0.6) + geom_area(data = df_f_tail, aes(x = bw, y = density), fill = "black", alpha = 0.8) + geom_vline(xintercept = Q, linetype = "dashed", linewidth = lsz) + coord_flip() + theme_void() + theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) + xlim(yl) + annotate("text", x = label_x, y = label_y, label = pct_label, color = "black", size = 3)

p_male <- ggplot(data = dds_male, aes(x = index, y = bw)) + geom_point(size = psz, alpha = transalpha, color = colg) + geom_point(data = dds10_male, size = psz, alpha = 0.6, color = "black") + geom_hline(yintercept = Q, linewidth = lsz, linetype = "dashed") + xlab("Index") + ylab("Birth weight (g)") + ylim(yl) + ggtitle("Males") + theme(axis.title.y = element_text(size = 9, vjust = 3.5))

p_female <- ggplot(data = dds_female, aes(x = index, y = bw)) + geom_point(size = psz, alpha = transalpha, color = colg) + geom_point(data = dds10_female, size = psz, alpha = 0.6, color = "black") + geom_hline(yintercept = Q, linewidth = lsz, linetype = "dashed") + xlab("Index") + ylab("") + scale_y_continuous(position = "right", limits = yl) + ylab("Birth weight (g)") + ggtitle("Females") + theme(axis.title.y = element_text(size = 9, vjust = 3.5))

g <- p_left + p_male + p_female + p_right + plot_layout(ncol = 4, widths = c(1, 5, 5, 1))




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

## R Code 12.2

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

nrow(USbirths2001) # about 3.5 million
# define variables
y <- USbirths2001$bw
d <- 1*(y < 2500)
x <- 1*(USbirths2001$multiple == 1)
# multiple birth prevalence in the population
print(px <- mean(x))
# LBW prevalence in the population (unconditional)
print(LBW <- mean(d))
# LBW prevalence conditional on multiple status
print(pi0 <- mean(d[x == 0])) # compare with lm intercept
print(pi1 <- mean(d[x == 1]))
# LBW prevalence difference between groups
print(dp <- pi1 - pi0) # compare with lm slope
# probability linear model
fit <- lm(d ~ x)
fit
# expected percentage points change of LBW in the population
# for a change du in the singleton/multiple composition
du <- 0.1
(dp*du)*100

## R Code 12.4

# use LBW prevalences (conditional and unconditional) as quantile levels
p <- round(c(pi0, LBW, pi1), 3)
# fit unconditional and conditional QR
res <- matrix(NA, length(p), 3)
res[,1] <- p
colnames(res) <- c("p", "UQR", "QR")
for(j in 1:length(p)){
xi <- quantile(y, p[j])
fhat <- density(y)
w <- xi + (p[j] - 1*(y <= xi))/approxfun(fhat$x, fhat$y)(xi)
res[j,2] <- round(lm(w ~ x)$coef[2])
res[j,3] <- round(rq(y ~ x, tau = p[j], method = "fn")$coef[2])
}
res