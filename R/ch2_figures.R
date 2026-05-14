# Statistical Methods and Applications Using Quantiles - Chapter 2

###################################################
# Packages and functions

rm(list = ls())
library(diagram)

###################################################
# Path output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\_figures\\"

###################################################
# Figure 4 - Quartiles

cx.lab <- 1.7
cx.ax <- 2
lw <- 2
cx.p <- 2

set.seed(789)
x <- rchisq(30, 3)
qs <- quantile(x, probs = 1:3/4)

setEPS(width = 12.64, height = 3.32, family = "Times", reset = TRUE, pointsize = 12)

fig <- "geraci02_04"
postscript(paste(path.out, fig, ".eps",sep=""), title = fig)

par(mar = c(3.1, 3.1, 2.1, 2.1), mgp = c(5, 2, 0))

plot(range(x),c(-1,1), type = "n", axes = F, ylab = "", xlab = "", xlim = c(0,14))
points(x, rep(0,length(x)), pch = 16, cex = cx.p)
Arrows(-1, 0, 13, 0, lwd = lw, arr.width = .25, arr.length = .25, arr.type = "triangle")
segments(qs[1],-0.35,qs[1],0.35, lty = 2, lwd = 2)
segments(qs[2],-0.35,qs[2],0.35, lty = 2, lwd = 2)
segments(qs[3],-0.35,qs[3],0.35, lty = 2, lwd = 2)
text(qs[1], -0.5, expression(xi[0.25]), cex = cx.lab)
text(qs[2], -0.5, expression(xi[0.5]), cex = cx.lab)
text(qs[3], -0.5, expression(xi[0.75]), cex = cx.lab)

text(qs[1] - 0.5, 0.5, "25%", cex = cx.lab)
text(qs[2] - 0.5, 0.5, "25%", cex = cx.lab)
text(qs[3] - 0.85, 0.5, "25%", cex = cx.lab)
text(qs[3] + 0.85, 0.5, "25%", cex = cx.lab)

dev.off()


<<02_04, echo=FALSE, fig.cap = "Quartiles divide the distribution in four parts weighting $25\\%$ each.">>=
library(diagram)
cx.lab <- 1.4
cx.ax <- 2
lw <- 2
cx.p <- 1.5
set.seed(789)
x <- rchisq(30, 3)
qs <- quantile(x, probs = 1:3/4)
par(mar = rep(1, 4), mai = c(0,1,1,1))
plot(range(x),c(-0.7,0.5), type = "n", axes = FALSE, ylab = "", xlab = "", xlim = c(0,14))
points(x, rep(0,length(x)), pch = 16, cex = cx.p)
Arrows(-1, 0, 13, 0, lwd = lw, arr.width = .25, arr.length = .25, arr.type = "triangle")
segments(qs[1],-0.35,qs[1],0.35, lty = 2, lwd = 2)
segments(qs[2],-0.35,qs[2],0.35, lty = 2, lwd = 2)
segments(qs[3],-0.35,qs[3],0.35, lty = 2, lwd = 2)
text(qs[1], -0.4, expression(xi[0.25]), cex = cx.lab)
text(qs[2], -0.4, expression(xi[0.5]), cex = cx.lab)
text(qs[3], -0.4, expression(xi[0.75]), cex = cx.lab)
text(qs[1] - 0.5, 0.2, "25%", cex = 1.3, srt = 90)
text(qs[2] - 0.5, 0.2, "25%", cex = 1.3, srt = 90)
text(qs[3] - 0.85, 0.2, "25%", cex = 1.3, srt = 90)
text(qs[3] + 0.85, 0.2, "25%", cex = 1.3, srt = 90)
@

###################################################
# Figure 7 - Diagram

setEPS(width = 12.64, height = 7.32, family = "Times", reset = TRUE, pointsize = 12)

cx.txt <- 2

names <- c("event", "real number", "probability")
M <- as.data.frame(matrix(nrow = 3, ncol = 3, byrow = TRUE, data = 0))
M[[2, 1]] <- "Y"
M[[3, 2]] <- "F[Y]"
M[[2, 3]] <- "Q[Y]"
C <- matrix(0,3,3)
C[2,3] <- 0.9

fig <- "geraci02_07"
postscript(file.path(path.out, "fig_02_07.eps",sep=""), title = fig)

par(mai = c(.5, .5, .5, .5), mar = c(0,1,1,1))
plotmat(M, pos = 3, curve = C, name = names, lwd = 1, arr.type = "triangle",
box.lwd = 2, box.type = "rect", box.prop = 1.0, dtext = 0.8,
box.cex = cx.txt, cex.txt = cx.txt)

dev.off()


<<02_07, echo=FALSE, fig.cap = "Diagram depicting the functional relationships between events, numbers and probabilities.">>=
library(diagram)
cx.lab <- 1.3
cx.ax <- 1.3
lw <- 1.5
cx.txt <- 1

windowsFont("Times")

names <- c("event", "real number", "probability")
M <- as.data.frame(matrix(nrow = 3, ncol = 3, byrow = TRUE, data = 0))
M[[2, 1]] <- "Y"
M[[3, 2]] <- "F[Y]"
M[[2, 3]] <- "Q[Y]"
C <- matrix(0,3,3)
C[2,3] <- 0.9

par(mar = c(2.1, 1.1, 2.1, 1.1), oma = c(0,0,0,0))
plotmat(M, pos = 3, curve = C, name = names, lwd = 1, arr.type = "triangle",
box.lwd = 2, box.type = "rect", box.prop = 1.0, dtext = 0.8,
box.cex = cx.txt, cex.txt = cx.txt)
@





















