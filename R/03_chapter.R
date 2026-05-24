# Statistical Methods and Applications Using Quantiles
# R code accompanying Chapter 3
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
### Section 3.1
###################################################

## R Code 3.1

qlss(fun = "qnorm", probs = c(0.1, 0.4))

## Table 3.1

load(file = "../data/EUSILC2015.Rdata")
K <- length(unique(eusilc2015$Country))
n <- length(eusilc2015$Country)

res1 <- tapply(eusilc2015$EDI, eusilc2015$Country, function(x) c(median(x)))
res2 <- tapply(eusilc2015$EDI, eusilc2015$Country, function(x) c(mean(x)))
res3 <- tapply(eusilc2015$EDI, eusilc2015$Country, function(x) c(ecdf(x)(mean(x))))
nn <- names(res1)

tab1 <- data.frame("Country code" = codes$Code[match(nn, codes$Country)], "Median" = res1, Mean = res2, "F" = res3)
rownames(tab1) <- nn

addtorow <- list()
addtorow$pos <- list(-1,0,0,nrow(tab1))
addtorow$command <- c("\\toprule",
" Country & Country code & \\multicolumn{1}{c}{Median} & \\multicolumn{1}{c}{Mean} & \\multicolumn{1}{c}{$F(\\mu)$} \\\\",
"\\cmidrule(l){3-5}","\\bottomrule"
)

print.xtable(xtable(tab1, digits = c(0,0,0,0,2), align = c(rep("l",2),rep("r",3)), caption = c("Median, mean and cumulative probability at the mean of the equivalised disposable income (EDI) in European countries, 2015.")), caption.placement="top", booktabs = TRUE, math.style.negative = TRUE, type="latex", comment=FALSE, add.to.row = addtorow, include.colnames = FALSE, hline.after = NULL, format.args = list(big.mark = ",", decimal.mark = "."))

## Figure 3.5

tmp <- tapply(eusilc2015$EDI, eusilc2015$Country, qlss, probs = c(0.1,0.4))
tmp <- lapply(tmp, function(x) unlist(x)[c(1:8)])
SK <- sapply(tmp, function(x) x[["shape.skewness.0.1"]])
PK <- sapply(tmp, function(x) x[["scale.IQR"]]/x[["scale.IPR.0.4"]])
TW <- sapply(tmp, function(x) x[["shape.shape.0.1"]])
K <- TW*PK

fit <- unlist(qlss("qnorm", probs = c(0.1,0.4)))
sk <- fit[["shape.skewness.0.1"]]
pk <- fit[["scale.IQR"]]/fit[["scale.IPR.0.4"]]
tw <- fit[["shape.shape.0.1"]]
k <- tw*pk

x <- names(SK)[order(SK)]
x <- codes$Code[match(x, codes$Country)]
x <- factor(x, levels = x)
dd0 <- data.frame(x = x, y = sort(SK))

x <- names(K)[order(K)]
x <- codes$Code[match(x, codes$Country)]
x <- factor(x, levels = x)
dd1 <- data.frame(x = x, y = sort(K))

x <- names(TW)[order(TW)]
x <- codes$Code[match(x, codes$Country)]
x <- factor(x, levels = x)
dd2 <- data.frame(x = x, y = sort(TW))

x <- names(PK)[order(PK)]
x <- codes$Code[match(x, codes$Country)]
x <- factor(x, levels = x)
dd3 <- data.frame(x = x, y = sort(PK))


sz <- 10
lsz <- 1
cl <- "black"

gg <- list()
gg[[1]] <- ggplot() + geom_bar(data = dd0, aes(x = x, y = y), stat = "identity", width = 0.5, color = "lightgrey", fill = "grey") + geom_hline(yintercept = sk, color = cl, linewidth = lsz, linetype = "dashed") + xlab("") + ylab("Skewness") + coord_cartesian(ylim = c(-1,1)) + theme(text = element_text(size=sz))
gg[[2]] <- ggplot() + geom_bar(data = dd1, aes(x = x, y = y), stat = "identity", width = 0.5, color = "lightgrey", fill = "grey") + geom_hline(yintercept = k, color = cl, linewidth = lsz, linetype = "dashed") + xlab("") + ylab("Kurtosis") + coord_cartesian(ylim = c(1, max(dd1$y))) + theme(text = element_text(size=sz))
gg[[3]] <- ggplot() + geom_bar(data = dd2, aes(x = x, y = y), stat = "identity", width = 0.5, color = "lightgrey", fill = "grey") + geom_hline(yintercept = tw, color = cl, linewidth = lsz, linetype = "dashed") + xlab("") + ylab("Tail weight") + coord_cartesian(ylim = c(1, max(dd2$y))) + theme(text = element_text(size=sz))
gg[[4]] <- ggplot() + geom_bar(data = dd3, aes(x = x, y = y), stat = "identity", width = 0.5, color = "lightgrey", fill = "grey") + geom_hline(yintercept = pk, color = cl, linewidth = lsz, linetype = "dashed") + xlab("Country") + ylab("Peakedness") + coord_cartesian(ylim = c(1, max(dd3$y))) + theme(text = element_text(size=sz))

marrangeGrob(gg, ncol = 1, nrow = 4, top = "")

## R Code 3.2

library(reldist)
tapply(eusilc2015$EDI, eusilc2015$Country, gini)

###################################################
## Section 3.2
###################################################

## Figure 3.6

load(file = "../data/sc2000wv.Rdata")

gg <- list()

gg[[1]] <- ggplot(sc2000wv, aes(x = rank, y = population)) +
  geom_line() +
  scale_x_continuous("Rank") +
  scale_y_log10(
    expression(paste("Population size (", log[10], ")")),
    labels = scales::label_number(big.mark = " ")
  ) +
  annotation_logticks(sides = "l") +
  theme(panel.grid.minor = element_blank())

gg[[2]] <- ggplot(sc2000wv, aes(x = population, y = after_stat(density))) +
  geom_histogram(bins = 10) +
  scale_x_continuous(
    expression(paste("Population size (", log[10], ")")),
    trans = "log10",
    labels = scales::label_number(big.mark = " ")
  ) +
  ylab("Density") + 
  annotation_logticks(sides = "b") + 
  theme(panel.grid.minor = element_blank(),
  axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)
  )

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")


## Figure 3.8

data(faithful)
cutoff <- 3
f <- faithful$eruptions >= cutoff

lsz <- 1
psz <- 1.6
transalpha <- 0.5

gg <- list()
gg[[1]] <- ggplot(faithful, aes(x = eruptions, y = waiting, group = f)) + geom_point(size = psz, alpha = transalpha, aes(shape = f)) + xlab("Duration (min)") + ylab("Waiting time (min)") + geom_vline(xintercept = cutoff, linetype = "dashed", color = "black", linewidth = lsz) +  theme(legend.position="none")
gg[[2]] <- ggplot(faithful, aes(x = waiting, linetype = f)) + stat_ecdf(geom = "step") + xlab("Waiting time (min)") + ylab("ECDF") +  theme(legend.position="none")

marrangeGrob(gg, ncol = 2, nrow = 1, top = "")

## R Code 3.3

data(faithful)
x <- as.numeric(faithful$eruptions >= cutoff)
y <- faithful$waiting
p <- c(0.1, 0.25, 0.5, 0.75, 0.9)
qhat <- tapply(y, x, quantile, probs = p)
muhat <- tapply(y, x, mean)
sdhat <- tapply(y, x, sd)
mat <- rbind(qhat[[1]], qhat[[2]], qhat[[2]] - qhat[[1]])
mat <- cbind(mat, mean = c(muhat, diff(muhat)), sd =
c(sdhat, sdhat[2]/sdhat[1]))
rownames(mat) <- c("x = 0", "x = 1", "Delta")
round(mat,1)


## R Code 3.4

ks.test(scale(y[x == 0]), "pnorm")
ks.test(scale(y[x == 1]), "pnorm")
