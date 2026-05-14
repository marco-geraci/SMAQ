# Statistical Methods and Applications Using Quantiles - Chapter 10

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(Qtools)
library(lubridate)
library(ggplot2)
library(ggforce)
library(usmap)
library(tidyr)
library(dplyr, warn.conflicts = FALSE)

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

load("C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results\\ch10_results.Rdata")

qhat <- predict(fit.midrq)
global_min <- min(qhat)
global_max <- max(qhat)

dd$pred <- qhat[,2]
dd1 <- subset(dd, dd$weekend == 0)
dd2 <- subset(dd, dd$weekend == 1)
center_point <- data.frame(lon = -96.5, lat = 37.5, label = "Geographic reference")
center_proj <- usmap_transform(center_point)

my_gradient <- scale_fill_gradient(low = "grey100", high = "black", name = "", na.value = "white", limits = c(global_min, global_max))

fig <- "fig_10_07"

png(paste(path.out, fig, ".png", sep = ""), width = 2100, height = 1050, units = "px", res = 300)
plot_usmap(regions = "counties", values = "pred", data = dd1, exclude = c("Alaska", "Hawaii")) + geom_sf(data = center_proj, color = "red", size = 3, shape = 17, fill = "red") + theme(panel.background=element_blank(), legend.position = "right") + my_gradient + theme_void()
dev.off()

#pdf(paste(path.out, fig, ".pdf", sep = ""), width = 7, height = 3.5)
#plot_usmap(regions = "counties", values = "pred", data = dd1, exclude = c("Alaska", "Hawaii")) + geom_sf(data = center_proj, color = "red", size = 3, shape = 17, fill = "red") + theme(panel.background=element_blank(), legend.position = "right") + my_gradient + theme_void()
#dev.off()

#geom_sf_text(data = center_proj, aes(label = label), nudge_y = 50000) 
