# Statistical Methods and Applications Using Quantiles - Chapter 12

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(Qtools)
library(ggplot2)
library(patchwork)
library(tidyverse)
library(fastglm)

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
# Figure

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")
load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001sample.Rdata")

# singletons born at term only
dd <- subset(USbirths2001, USbirths2001$ga >= 37 & USbirths2001$multiple == 0)
nrow(dd) # about 3 million
dd_male <- dd %>% filter(sex == "male")
dd_female <- dd %>% filter(sex == "female")

psz <- 2
lsz <- 1
transalpha <- 0.4
colg <- "grey70"
yl <- c(900,6000)

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

pct <- mean(dd_male$bw < Q) * 100
pct_label <- paste0(sprintf("%.2f", pct), "%")
label_x <- 2000
label_y <- -max(df_m$density[df_m$bw < Q]) * 4.5

p_left <- ggplot() + geom_area(data = df_m_main, aes(x = bw, y = -density), fill = "grey70", alpha = 0.6) + geom_area(data = df_m_tail, aes(x = bw, y = -density), fill = "black", alpha = 0.8) + geom_vline(xintercept = Q, linetype = "dashed", linewidth = lsz) + coord_flip() + theme_void() + theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) + xlim(yl) + annotate("text", x = label_x, y = label_y, label = pct_label, color = "black", size = 3)

pct <- mean(dd_female$bw < Q) * 100
pct_label <- paste0(sprintf("%.2f", pct), "%")
label_x <- 2000
label_y <- max(df_f$density[df_f$bw < Q]) * 4.5

p_right <- ggplot() + geom_area(data = df_f_main, aes(x = bw, y = density), fill = "grey70", alpha = 0.6) + geom_area(data = df_f_tail, aes(x = bw, y = density), fill = "black", alpha = 0.8) + geom_vline(xintercept = Q, linetype = "dashed", linewidth = lsz) + coord_flip() + theme_void() + theme(plot.margin = margin(t = 20, r = 10, b = 20, l = 0)) + xlim(yl) + annotate("text", x = label_x, y = label_y, label = pct_label, color = "black", size = 3)

p_male <- ggplot(data = dds_male, aes(x = index, y = bw)) + geom_point(size = psz, alpha = transalpha, color = colg) + geom_point(data = dds10_male, size = psz, alpha = 0.6, color = "black") + geom_hline(yintercept = Q, linewidth = lsz, linetype = "dashed") + xlab("Index") + ylab("Birth weight (g)") + ylim(yl) + ggtitle("Males") + theme(axis.title.y = element_text(size = 9, vjust = 3.5))

p_female <- ggplot(data = dds_female, aes(x = index, y = bw)) + geom_point(size = psz, alpha = transalpha, color = colg) + geom_point(data = dds10_female, size = psz, alpha = 0.6, color = "black") + geom_hline(yintercept = Q, linewidth = lsz, linetype = "dashed") + xlab("Index") + ylab("") + scale_y_continuous(position = "right", limits = yl) + ylab("Birth weight (g)") + ggtitle("Females") + theme(axis.title.y = element_text(size = 9, vjust = 3.5))

g <- p_left + p_male + p_female + p_right + plot_layout(ncol = 4, widths = c(1, 5, 5, 1))

fig <- "fig_12_01"

png(paste(path.out, fig, ".png", sep = ""), width = 3000, height = 1500, units = "px", res = 300)
g
dev.off()



#pdf(paste(path.out, fig, ".pdf", sep = ""), width = 10, height = 5)
#g
#dev.off()

