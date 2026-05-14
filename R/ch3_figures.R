# Statistical Methods and Applications Using Quantiles - Chapter 3

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(xtable)
library(diagram)
library(mvtnorm)
library(MASS)
library(Qtools)
library(sn)
library(pBrackets)
library(ggplot2)
library(grid)
library(sfsmisc)
library(lqmm)

###################################################
# Path output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\_figures\\"

###################################################
# Figure x - Shape 1
cx.lab <- 2.5
cx.ax <- 2
lw <- 2
cx.p <- 2
yl <- c(0, 300)
colg <- rgb(179, 179, 179, max = 255) # grey(0.7)

setEPS(width = 12.64, height = 9.32, family = "Times", reset = TRUE, pointsize = 12)

fig <- "geraci03_02"
#postscript(paste(path.out, fig, ".eps",sep=""), title = fig)

par(mar = c(7.1, 7.1, 2.1, 2.1), mgp = c(5, 2, 0))

plot(c(-3,3), c(0,.5), type = "n", axes = TRUE, xlab = "", ylab = "", xlim = c(-4, 4))

tmp <- function(x) dnorm(x, mean = 0, sd = 1)
curve(dnorm, -3, 3, add = TRUE, lwd = lw)

alpha <- c(-1, -3, -10)
delta <- alpha/sqrt(1+alpha^2)
omega <- (1 - (2/pi)*delta^2)^(-1/2)
phi <- omega^2 * (1 - (2/pi)*delta^2)

for( i in 3){

tmp <- function(x) dsn(x, xi = 0, omega = omega[i], alpha = alpha[i])
curve(tmp, -4, 2, lty = 2, add = T, col = colg, lwd = lw)
}

#dev.off()

n <- 10000
set.seed(89)
y <- c(rsn(n, xi = 0, omega = omega[1], alpha = alpha[1]),
	rsn(n, xi = 0, omega = omega[2], alpha = alpha[2]),
	rsn(n, xi = 0, omega = omega[3], alpha = alpha[3]))
f <- rep(1:3, each = n)
f <- factor(f, labels = c('a','b', 'c'))
x <- data.frame(f = f, y = c(y,y))
tapply(y, f, mean)


xb <- "y"
cx.lab <- 50
cx.ax <- 40
cx.leg <- 30

bp <- ggplot(x, aes(x=y, fill=f)) + geom_density(alpha=.3, kernel = "gaussian", bw = .5)
bp + theme(axis.title.x = element_text(size=cx.lab),
           axis.text.x  = element_text(size=cx.ax), axis.title.y = element_text(size=cx.lab), axis.text.y  = element_text(size=cx.ax), legend.title = element_text(size=cx.leg), legend.text = element_text(size = cx.leg), legend.position="top", legend.key.size = unit(20, "points")) + xlab("\ny") + ylab("Density\n") + scale_fill_discrete(name="none")




###################################################
# Figure 4 - FNS

cx.lab <- 1.7
cx.ax <- 2
lw <- 2
cx.p <- 2

qs <- c(0.5, 3.0, 8.0, 11.0, 13.0)

setEPS(width = 12.64, height = 3.32, family = "Times", reset = TRUE, pointsize = 12)

fig <- "geraci03_04"
postscript(paste(path.out, fig, ".eps",sep=""), title = fig)

par(mar = c(3.1, 3.1, 2.1, 2.1), mgp = c(5, 2, 0))

plot(range(qs), c(-1,1), type = "n", axes = F, ylab = "", xlab = "", xlim = c(0,16))
Arrows(-1, 0, 16, 0, lwd = lw, arr.width = .25, arr.length = .25, arr.type = "triangle")
segments(qs[1],-0.35,qs[1],0.67, lty = 2, lwd = lw)
segments(qs[2],-0.35,qs[2],0.35, lty = 2, lwd = lw)
segments(qs[3],-0.35,qs[3],0.67, lty = 2, lwd = lw)
segments(qs[4],-0.35,qs[4],0.35, lty = 2, lwd = lw)
segments(qs[5],-0.35,qs[5],0.67, lty = 2, lwd = lw)

text(qs[1], -0.5, expression(xi[0] == 0.5), cex = cx.lab)
text(qs[2], -0.5, expression(xi[0.25] == 3), cex = cx.lab)
text(qs[3], -0.5, expression(xi[0.5] == 8), cex = cx.lab)
text(qs[4], -0.5, expression(xi[0.75] == 11), cex = cx.lab)
text(qs[5], -0.5, expression(xi[1] == 13), cex = cx.lab)

#brackets(qs[2], 0.35, qs[4], 0.35, h = .25, curvature = 0.5, lwd = lw, lty = 1, ticks = NA, type = 4)
#brackets(qs[1], 0.35, qs[5], 0.35, h = .5, curvature = 0.5, lwd = lw, lty = 1, ticks = NA, type = 4)
y1 <- 0.38
y2 <- 0.70
delta <- 0.05
a <- 0.15

arrows(qs[2], y1, qs[3] - delta, y1, lwd = lw, code = 3, length=0.1)
arrows(qs[3] + delta, y1, qs[4], y1, lwd = lw, code = 3, length=0.1)

arrows(qs[1], y2, qs[3] - delta, y2, lwd = lw, code = 3, length=0.1)
arrows(qs[3] + delta, y2, qs[5], y2, lwd = lw, code = 3, length=0.1)

text((qs[3] + qs[2])/2, y1 + a, "25%", cex = cx.lab)
text((qs[4] + qs[3])/2, y1 + a, "25%", cex = cx.lab)
text((qs[3] + qs[1])/2, y2 + a, "50%", cex = cx.lab)
text((qs[5] + qs[3])/2, y2 + a, "50%", cex = cx.lab)

#points(4, 0, pch = 16, cex = cx.p)

dev.off()


<<03_04, echo=FALSE, fig.cap = "Five-number summary of recovery time in diseased patients in a hypothetical example.">>=
library(diagram)
cx.lab <- 1.3
cx.ax <- 1.3
lw <- 1.5

qs <- c(0.5, 3.0, 8.0, 11.0, 13.0)

par(mar = c(2.1, 1.1, 2.1, 1.1), oma = c(0,0,0,0))

plot(range(qs), c(-1,1), type = "n", axes = F, ylab = "", xlab = "", xlim = c(0,16))
Arrows(-1, 0, 16, 0, lwd = lw, arr.width = .25, arr.length = .25, arr.type = "triangle")
segments(qs[1],-0.35,qs[1],0.67, lty = 2, lwd = lw)
segments(qs[2],-0.35,qs[2],0.35, lty = 2, lwd = lw)
segments(qs[3],-0.35,qs[3],0.67, lty = 2, lwd = lw)
segments(qs[4],-0.35,qs[4],0.35, lty = 2, lwd = lw)
segments(qs[5],-0.35,qs[5],0.67, lty = 2, lwd = lw)

text(qs[1], -0.5, expression(xi[0] == 0.5), cex = cx.lab)
text(qs[2], -0.5, expression(xi[0.25] == 3), cex = cx.lab)
text(qs[3], -0.5, expression(xi[0.5] == 8), cex = cx.lab)
text(qs[4], -0.5, expression(xi[0.75] == 11), cex = cx.lab)
text(qs[5], -0.5, expression(xi[1] == 13), cex = cx.lab)

y1 <- 0.38
y2 <- 0.70
delta <- 0.05
a <- 0.15

arrows(qs[2], y1, qs[3] - delta, y1, lwd = lw, code = 3, length=0.1)
arrows(qs[3] + delta, y1, qs[4], y1, lwd = lw, code = 3, length=0.1)

arrows(qs[1], y2, qs[3] - delta, y2, lwd = lw, code = 3, length=0.1)
arrows(qs[3] + delta, y2, qs[5], y2, lwd = lw, code = 3, length=0.1)

text((qs[3] + qs[2])/2, y1 + a, "25%", cex = cx.lab)
text((qs[4] + qs[3])/2, y1 + a, "25%", cex = cx.lab)
text((qs[3] + qs[1])/2, y2 + a, "50%", cex = cx.lab)
text((qs[5] + qs[3])/2, y2 + a, "50%", cex = cx.lab)
@





