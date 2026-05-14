# Statistical Methods and Applications Using Quantiles - Chapter 8

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(mvtnorm)
library(MASS)
library(Qtools)
library(ggplot2)
library(grid)
library(lqmm)

###################################################
# Path output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_figures\\"

###################################################
# Set ggplot global options

# lightgrey (#D3D3D3), gainsboro (#DCDCDC), whitesmoke (#F5F5F5)
theme_set(
  theme_minimal(base_size = 12) +
    theme(
      text = element_text(family = "serif"),
      plot.title = element_text(margin = margin(b = 10), vjust = 1),

      axis.title.x = element_text(margin = margin(t = 10)),
      axis.title.x.top = element_text(margin = margin(b = 10)),
      axis.text.x.top = element_text(margin = margin(b = 4)),

      axis.title.y = element_text(margin = margin(r = 10)),
      axis.title.y.right = element_text(margin = margin(l = 10)),
      axis.text.y.right = element_text(margin = margin(l = 4)),

      panel.grid.major = element_line(color = "gainsboro"),
      panel.grid.minor = element_line(color = "gainsboro")
    )
)


###################################################
# Figure 8
load(file.path("C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata"))

lsz <- 0.5
dd <- USbirths2001
n <- length(dd$bw)
y <- dd$bw
g <- dd$cause
z <- log(y)

fhat1 <- density(y[g == "ref"], from = 200, to = 7000)
x1 <- fhat1$x
y1 <- fhat1$y

fhat2 <- density(y[g == "sids"], from = 200, to = 7000)
x2 <- fhat2$x
y2 <- fhat2$y


g <- ggplot() + geom_line(aes(x = x1, y = y1), linetype = "solid", linewidth = lsz) +
geom_line(aes(x = x2, y = y2), linetype = "dashed", linewidth = lsz) +
xlab("Birth weight (g)") +
ylab("Density")


fig <- "fig_08_08"

png(paste(path.out, fig, ".png", sep = ""), width = 2100, height = 1050, units = "px", res = 300)
g
dev.off()


#pdf(paste(path.out, fig, ".pdf", sep = ""), width = 7, height = 3.5)
#g
#dev.off()



<<08_08, warning=FALSE, echo=FALSE, cache=TRUE, fig.cap="Estimated density of birth weights in infants who survived past one year of age (solid line) and infants who died of sudden infant death syndrome (dashed line) using the NCHS Births dataset.">>=
load(file.path("_datasets/USbirths2001.Rdata"))
lsz <- 0.5
dd <- USbirths2001
n <- length(dd$bw)
y <- dd$bw
g <- dd$cause
z <- log(y)

fhat1 <- density(y[g == "ref"], from = 200, to = 7000)
x1 <- fhat1$x
y1 <- fhat1$y

fhat2 <- density(y[g == "sids"], from = 200, to = 7000)
x2 <- fhat2$x
y2 <- fhat2$y


gg <- ggplot() + geom_line(aes(x = x1, y = y1), linetype = "solid", linewidth = lsz) +
geom_line(aes(x = x2, y = y2), linetype = "dashed", linewidth = lsz) +
xlab("Birth weight (g)") +
ylab("Density")
print(gg)
@
