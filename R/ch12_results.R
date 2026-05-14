# Statistical Methods and Applications Using Quantiles - Chapter 12

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(Qtools)
library(ggplot2)
library(gridExtra)
library(patchwork)
library(tidyverse)
library(fastglm)

###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")

dd <- subset(USbirths2001, USbirths2001$ga >= 37 & USbirths2001$multiple == 0)

set.seed(123)
n <- 10000
USbirths2001sample <- dd[sample(seq_along(dd$bw), n),]
USbirths2001sample <- USbirths2001sample[,c("bw","sex")]

save(list = c("USbirths2001sample"), file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001sample.Rdata")

###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")

dd <- subset(USbirths2001, USbirths2001$ga >= 37 & USbirths2001$multiple == 0)
# Fit logistic regression on LBW indicator
# conditional on sex (x = 1 if female)
library(fastglm)
Q <- 2500
y <- dd$bw
z <- 1*(y < Q)
x <- 1*(dd$sex == "female")
fit <- fastglm(cbind(1,x), z, family = binomial(), method = 3)

# conditional probabilities of birthweight less than 2500
xx <- cbind(1, c(0,1))
eta <- xx %*% coef(fit)
phat <- as.numeric(1/(1 + exp(-eta)))
names(phat) <- levels(dd$sex)
phat

# verify conditional probabilities yield same
# conditional quantile (i.e., 2500)
round(quantile(y[x == 0], probs = phat[1], names = FALSE))
round(quantile(y[x == 1], probs = phat[2], names = FALSE))

# Unconditional LBW probability
ecdf(y)(Q)
# Unconditional LBW probability as weighted average of
# conditional probabilities
as.numeric(phat[1] + diff(phat)*mean(x))

###################################################
# Figure

p0_vals <- seq(0.01, 0.99, length.out = 200)
odds_ratios <- c(1.5, 2, 3, 5, 10)
df <- data.frame(p0 = rep(p0_vals, length(odds_ratios)), OR = rep(odds_ratios, each = length(p0_vals)), delta = NA)

for (OR in odds_ratios) {
  odds_0 <- p0_vals / (1 - p0_vals)
  odds_1 <- OR * odds_0
  delta <- (OR - 1) / (1 + 1/odds_0 + OR + odds_1)
  df$delta[df$OR == OR] <- rep(delta, length.out = length(p0_vals))
}

ggplot(df, aes(x = p0, y = delta, color = factor(OR))) + geom_line(linewidth = 1.2) + labs(x = expression(pi[0]), y = expression(pi[1]-pi[0]), color = "Odds ratio") + scale_color_grey(start = 1, end = 0.2) + theme_minimal() + theme(legend.position = "top")

###################################################
# Simulation

library(quantreg)
library(reshape2)
library(epitools)

n <- 1000
p <- c(0.1, 0.5, 0.75, 0.95, 0.99)
np <- length(p)
b <- c(0,-1)
g0 <- c(1, 5)
g1 <- c(0, 1)
R <- 1000
eff_labs <- c("Unit","Inflated")

true <- array(NA, dim = c(2, length(g0), length(g1), np), dimnames = list(model = c("QTE","RR"), size = eff_labs, var = c("Homoskedastic", "Heteroskedastic"), quantile = p))
for(j in 1:length(g1)){
for(h in 1:length(g0)){
print(eff_labs[h])
m <- as.numeric(cbind(1, c(0,1))%*%b)
s <- as.numeric(cbind(1, c(0,1))%*%c(g0[h],g1[j]))
print(diff(m)/max(s))
true["QTE",h,j,] <- b[2] + g1[j]*qnorm(p, 0, 1)
k <- qnorm(p, m[2], s[2])
pi0 <- pnorm(k, m[1], s[1])
true["RR",h,j,] <- exp(log(p/(1-p)) - log(pi0/(1-pi0)))
}}

res <- eff <- array(NA, dim = c(R, 2, length(g0), length(g1), np), dimnames = list(rep = 1:R, model = c("QTE", "RR"), size = eff_labs, var = c("Homoskedastic", "Heteroskedastic"), quantile = p))

for(j in 1:length(g1)){
se <- if(j == 1) "iid" else "nid"
for(h in 1:length(g0)){
	for(r in 1:np){
	m1 <- sum(c(1,1)*b)
	s1 <- sum(c(1,1)*c(g0[h],g1[j]))
	k <- qnorm(p[r], m1, s1)

for(i in 1:R){
set.seed(i)
x <- rbinom(n, 1, 0.5)
eta <- cbind(1, x)%*%b
s <- cbind(1, x)%*%c(g0[h],g1[j])
y <- eta + rnorm(n, 0, s)
y0 <- y[x == 0]
y1 <- y[x == 1]
#F0 <- ecdf(y0)
#F1 <- ecdf(y1)
z <- 1*(y <= k)

#fit.glm <- glm(z ~ x, family = binomial)
fit.rr <- try(riskratio(xtabs(~x+z)), silent = TRUE)
fit.rq <- rq(y ~ x, tau = p[r])
res[i,"QTE",h,j,r] <- summary(fit.rq, se = "boot")$coef[2,4]
#res[i,"OR",h,j,r] <- summary(fit.glm)$coef[2,4]
res[i,"RR",h,j,r] <- if(inherits(fit.rr, "try-error")) NA else fit.rr$p.value[2,3]
eff[i,"QTE",h,j,r] <- coef(fit.rq)[2]
#eff[i,"OR",h,j,r] <- exp(coef(fit.glm)[2])
eff[i,"RR",h,j,r] <- if(inherits(fit.rr, "try-error")) NA else fit.rr$measure[2,1]
}
}}}

###################################################
# Save output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results"
save(list = c("res", "eff", "true"), file = file.path(path.out, "ch12_results.Rdata"), compress = "gzip")

###################################################
# Table 12.1

<<tab:12_01, echo=FALSE, results = 'asis'>>=
load(file = file.path("_results/ch12_results.Rdata"))

p <- dimnames(res)$quantile
np <- length(p)

pv <- NULL
for(j in 1:2){
for(h in 1:2){
pv <- rbind(pv, apply(res[,,h,j,] < 0.05, 3, colMeans))
}}

tab <- matrix(NA, 8, np*2)
tab[,c(1,3,5,7,9)] <- round(matrix(melt(true)$value, ncol = np, byrow = FALSE), 2)
tab[,c(2,4,6,8,10)] <- pv
tab <- data.frame(rownames(pv), tab)

dig <- matrix(2, 8, np*2 + 2)

addtorow <- list()
addtorow$pos <- list(-1,0,0,0,0,0,4,4,4,nrow(tab))
addtorow$command <- c("\\toprule",
" Effect size & \\multicolumn{2}{c}{$p = 0.1$} & \\multicolumn{2}{c}{$p = 0.25$} \\ & \\multicolumn{2}{c}{$p = 0.5$} \\ & \\multicolumn{2}{c}{$p = 0.75$} \\ & \\multicolumn{2}{c}{$p = 0.9$} \\\\",
" & Value & Power & Value & Power & Value & Power & Value & Power & Value & Power \\\\",
"\\cmidrule(l){2-11}",
"& \\multicolumn{10}{c}{\\textit{Homoskedastic}}\\\\",
"\\cmidrule(l){2-11}",
"\\cmidrule(l){2-11}",
"& \\multicolumn{10}{c}{\\textit{Heteroskedastic}}\\\\",
"\\cmidrule(l){2-11}",
"\\bottomrule")
header <- function(x){
paste0('{\\bfseries ', x, '}')
}

print.xtable(xtable(tab, digits = dig, align = c("l",rep("r",11)), caption = c("Power.\\label{tab:12_01}")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.rownames = FALSE, include.colnames = FALSE, hline.after = NULL, sanitize.text.function = identity)
@

###################################################
# Figure

gg <- list()
psz <- 2.5
lsz <- 1
transalpha <- 0.4
colg <- "grey70"

p <- 1:99/100
gr <- factor(rep(c("log-RR","QTE"), each = length(p)), levels = c("QTE", "log-RR"))

b0 <- 0
b1 <- -1

#

g0 <- 1
g1 <- 0

cutoff <- qnorm(p, b0 + b1, g0 + g1)
R0 <- pnorm(cutoff, b0, g0)
R1 <- pnorm(cutoff, b0 + b1, g0 + g1)

deltap <- b1 + g1*qnorm(p)
pip <- pnorm(deltap/g0 + qnorm(p))
RR <- log(p/pip)

dd1 <- data.frame(p = c(p, p), y = c(RR, deltap), gr = gr)

#

g0 <- 5
g1 <- 0

cutoff <- qnorm(p, b0 + b1, g0 + g1)
R0 <- pnorm(cutoff, b0, g0)
R1 <- pnorm(cutoff, b0 + b1, g0 + g1)

deltap <- b1 + g1*qnorm(p)
pip <- pnorm(deltap/g0 + qnorm(p))
RR <- log(p/pip)

dd2 <- data.frame(p = c(p, p), y = c(RR, deltap), gr = gr)

#

g0 <- 1
g1 <- 1

cutoff <- qnorm(p, b0 + b1, g0 + g1)
R0 <- pnorm(cutoff, b0, g0)
R1 <- pnorm(cutoff, b0 + b1, g0 + g1)

deltap <- b1 + g1*qnorm(p)
pip <- pnorm(deltap/g0 + qnorm(p))
RR <- log(p/pip)

dd3 <- data.frame(p = c(p, p), y = c(RR, deltap), gr = gr)

#

g0 <- 5
g1 <- 1

cutoff <- qnorm(p, b0 + b1, g0 + g1)
R0 <- pnorm(cutoff, b0, g0)
R1 <- pnorm(cutoff, b0 + b1, g0 + g1)

deltap <- b1 + g1*qnorm(p)
pip <- pnorm(deltap/g0 + qnorm(p))
RR <- log(p/pip)

dd4 <- data.frame(p = c(p, p), y = c(RR, deltap), gr = gr)

#

yl <- range(c(dd1$y,dd2$y,dd3$y,dd4$y))

gg[[1]] <- ggplot(dd1, aes(x = p, y = y, group = gr)) + geom_line(aes(linetype = gr), linewidth = lsz) + xlab("p") + ylab("Effect size") + xlim(c(0,1)) + ylim(yl) + theme(legend.position = "none") + ggtitle("(a)")

gg[[3]] <- ggplot(dd2, aes(x = p, y = y, group = gr)) + geom_line(aes(linetype = gr), linewidth = lsz) + xlab("p") + ylab("") + xlim(c(0,1)) + ylim(yl) + theme(legend.position = "none") + ggtitle("(b)")

gg[[2]] <- ggplot(dd3, aes(x = p, y = y, group = gr)) + geom_line(aes(linetype = gr), linewidth = lsz) + xlab("p") + ylab("Effect size") + xlim(c(0,1)) + ylim(yl) + theme(legend.position = "none") + ggtitle("(c)")

gg[[4]] <- ggplot(dd4, aes(x = p, y = y, group = gr)) + geom_line(aes(linetype = gr), linewidth = lsz) + xlab("p") + ylab("") + xlim(c(0,1)) + ylim(yl) + theme(legend.position = "none") + ggtitle("(d)")

marrangeGrob(gg, ncol = 2, nrow = 2, top = "")

###################################################

y <- USbirths2001sample$bw
z <- 1*(y < 2500)
x <- 1*(USbirths2001sample$sex == "female")
cutoff <- seq(2000, 3000, by = 100)
p <- ecdf(y[x == 1])(cutoff)
deltap <- rq(y ~ x, tau = p, method = "fn")$coef[2,]
plot(p, deltap)

OR <- exp(sapply(cutoff, function(k, x, y) glm(I(y <= k) ~ x, family = binomial)$coef[2], y = y, x = x))

fit <- glm(z ~ x, family = binomial)

# conditional probabilities of birthweight less than 2500
xx <- cbind(1, c(0,1))
eta <- xx %*% coef(fit)
phat <- as.numeric(1/(1 + exp(-eta)))
names(phat) <- levels(USbirths2001sample$sex)
phat

# verify that the estimated probabilities correspond to
# the same empirical quantile (i.e., 2500)
round(quantile(y[x == 0], probs = phat[1], names = FALSE))
round(quantile(y[x == 1], probs = phat[2], names = FALSE))