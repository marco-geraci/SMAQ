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
# Data
n <- 2000

p <- c(1:3/4, 0.9)
R <- 100
res <- array(NA, dim = c(R, length(p), 2))

set.seed(123)
for(i in 1:R){
if((i %% 25) == 0)print(i)
x1 <- scale(rbeta(n, 5/3, 5/3))*5/6
x2 <- rbinom(n, 1, 0.2)
mu <- exp(1 + 1*x1 + 0*x2)
y <- sort(rpois(n, mu))
y[y %in% c(0)] <- 0.1
y[y %in% c(3,4)] <- 2
y[y %in% c(7:9)] <- 6
y[y %in% c(14:15)] <- 16
#y <- rnbinom(n, mu = mu, size = 2)
	for(j in 1:length(p)){
		fit1 <- rq.counts(y ~ x1 + x2, tau = p[j], M = 50)
		fit2 <- midrq(y ~ x1 + x2, tau = p[j], lambda = 0, control = midrqControl(ecdf_est = "logit"))
		res[i,j,1] <- fit1$tTable[3,5]
		est <- coef(fit2)
		se <- summary(fit2)$tTable[3,2]
		res[i,j,2] <- 2 * pt(-abs(est[3]/se), n - length(est))
	}
}

colMeans(res[,,1] < 0.1, na.rm = TRUE)
colMeans(res[,,2] < 0.1, na.rm = TRUE)


###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\FARS2020.Rdata")

dd <- dds %>% filter(year == 2020) %>% group_by(weekend, fips) %>% summarise(count = n(), .groups = "drop") %>% complete(weekend, fips, fill = list(count = 0))
dd <- merge(dd, centroids[,c("fips","longitude","latitude")], all.x = TRUE)
dd <- subset(dd, complete.cases(dd))
dd$weekend <- ifelse(dd$weekend == "weekend", 1, 0)

fars <- dd[,c("fips", "weekend", "longitude", "latitude", "count")]
save(list = c("fars"), file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\fars.Rdata")

#save(list = c("fars"), file = "C:\\Work\Research\\Projects\\I04 - R-Qtools\\Qtools_21\\Qtools\\data\\fars.rda")


###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\fars.Rdata")
dd <- fars

# offset weekend/weekday
dd$offset <- log(ifelse(dd$weekend == 1, 2, 5))
# replace zero counts with small constant
dd$y <- ifelse(dd$count == 0, 0.1, dd$count)
# trasform longitude and latitude
theta_c <- 37.5*pi/180
phi_c <- -96.5*pi/180
lat_rad <- dd$latitude*pi/180
lon_rad <- dd$longitude*pi/180
dd$lat_cos <- cos(lat_rad) - cos(theta_c)
dd$lat_sin <- sin(lat_rad) - sin(theta_c)
dd$lon_cos <- cos(lon_rad) - cos(phi_c)
dd$lon_sin <- sin(lon_rad) - sin(phi_c)

fit.midrq <- midrq(y ~ weekend + lat_cos + lat_sin + lon_cos + lon_sin,
tau = c(0.75, 0.9), data = dd, offset = dd$offset, lambda = 0,
control = midrqControl(ecdf_est = "logit"))

vcov(fit.midrq, numerical = TRUE, robust = TRUE)
summary(fit.midrq, numerical = TRUE, robust = TRUE)

p <- 2:9/10
fit.p <- midrq(y ~ weekend + lat_cos + lat_sin + lon_cos + lon_sin,
tau = p, data = dd, offset = dd$offset, lambda = 0,
control = midrqControl(ecdf_est = "logit"))
se <- sapply(vcov(fit.p, numerical = TRUE), function(x) sqrt(diag(x)[2]))
b <- coef(fit.p)[2,]
plot(p, b, xlim = c(0,1), ylim = c(-1,0.4), type = "b")
lines(p, b - qnorm(0.975)*se, col = 2)
lines(p, b + qnorm(0.975)*se, col = 2)
wend <- data.frame(p = p, est = b, se = se, lower = b - qnorm(0.975)*se, upper = b + qnorm(0.975)*se)

###################################################
# Save output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results"
save(list = c("fit.midrq", "dd", "wend"), file = file.path(path.out, "ch10_results.Rdata"), compress = "gzip")

###################################################
# Figure

k <- 10
x <- seq(min(dd$longitude), max(dd$longitude), length = k)
y <- seq(min(dd$latitude), max(dd$latitude), length = k)
df <- expand.grid(longitude = x, latitude = y)
df <- df %>% mutate(
    lon_rad = longitude*pi/180,
    lat_rad = latitude*pi/180,
    lon_cos = cos(lon_rad),
    lon_sin = sin(lon_rad),
    lat_cos = cos(lat_rad),
    lat_sin = sin(lat_rad)
)
df <- rbind(df, df)
df$weekend <- factor(rep(levels(dd$weekend), each = k^2), levels = levels(dd$weekend))
df$lonr <- cut(df$longitude, breaks = 10)
df$latr <- cut(df$latitude, breaks = 10)

qhat <- predict(fit.midrq, newdata = df)
global_min <- min(qhat)
global_max <- max(qhat)

df$pred <- qhat[,2]
dd1 <- subset(df, df$weekend == "weekday")
dd2 <- subset(df, df$weekend == "weekend")

ggplot(dd2, aes(x = latitude, y = pred)) + facet_wrap(~ lonr) + geom_line()
ggplot(dd2, aes(x = longitude, y = pred)) + facet_wrap(~ latr) + geom_line()

###################################################
# Figure

psz <- 2
lsz <- 1
transalpha <- 0.4

n <- 500
set.seed(123)
x1 <- scale(rbeta(n, 5/3, 5/3))*5/6
x2 <- rbinom(n, 1, 0.2)
mu <- exp(1 + 1*x1 + 0*x2)
y <- sort(rpois(n, mu))
y[y %in% c(0)] <- 0.1
y[y %in% c(3,4)] <- 2
y[y %in% c(7:9)] <- 6
y[y %in% c(14:15)] <- 16
w <- rep(0, max(y)+1)
w[as.numeric(names(table(y)))+1] <- table(y)
names(w) <- 0:max(y)


gg <- ggplot() + geom_segment(aes(x = 0:max(y), xend = 0:max(y), y = 0, yend = w), linewidth = lsz, color = "black") + xlab("y") + ylab("Frequency")

gg


