# Statistical Methods and Applications Using Quantiles - Chapter 1

###################################################
# Packages and functions

rm(list = ls())
library(quantreg)
library(xtable)
library(Qtools)

###################################################
# Data

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths.Rdata")

USbirths2001 <- subset(DATADy, DATADy$year == 2001)
USbirths2001$cod <- gsub("[[:space:]]", "", USbirths2001$cod)
USbirths2001$died <- factor(as.numeric(USbirths2001$cod != ""), levels = c(0,1))

range(USbirths2001$gest, na.rm = TRUE)
USbirths2001$ga_group <- cut(USbirths2001$gest, c(16.99, 27, 31, 36, 47))
USbirths2001$ga_group <- relevel(USbirths2001$ga_group, ref = "(36,47]")

USbirths2001$tox <- factor(USbirths2001$tox, levels = 0:3, c("ref", "alcohol", "smoking", "alcohol and smoking"))

# cause of death group
x1 <- ifelse(USbirths2001$cod == "089", 1, 0)
x2 <- ifelse(USbirths2001$cod == "135", 2, 0)
x3 <- ifelse(!USbirths2001$cod %in% c("", "135", "089"), 3, 0)
x <- x1+x2+x3
USbirths2001$cause <- factor(x, levels = 0:3, labels = c("ref", "lbw/pt", "sids", "other"))
USbirths2001$sex <- factor(USbirths2001$IsFemale, levels = 0:1, labels = c("male", "female"))

sel <- c("bw", "gest", "ga_group", "age_m", "sex", "IsFemale", "IsMult", "IsSmoker", "IsDrinker", "tox", "died", "cause")
USbirths2001 <- USbirths2001[,sel]
names(USbirths2001) <- c("bw", "ga", "ga_group", "age_mother", "sex", "female", "multiple", "smoking", "alcohol", "tox", "died", "cause")

summary(USbirths2001)
dim(USbirths2001) # [1] 3477802      12

save(list = c("USbirths2001"), file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")

###################################################
# Analysis

load(file = "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_datasets\\USbirths2001.Rdata")

p <- c(0.1,0.5,0.9)

run <- FALSE
if(run){

sel <- USbirths2001$ga >= 37
fit1.lm <- lm(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001, subset = sel)
fit1.rq <- rq(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001, tau = p, method = "fn", subset = sel)
fit1.rqs <- summary(fit1.rq, se = "nid")

sel <- USbirths2001$ga < 37
fit2.lm <- lm(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001, subset = sel)
fit2.rq <- rq(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001, tau = p, method = "fn", subset = sel)
fit2.rqs <- summary(fit2.rq, se = "nid")

fit3.lm <- lm(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001)
fit3.rq <- rq(bw ~ female + multiple + I(age_mother - 30) + smoking + alcohol, data = USbirths2001, tau = p, method = "fn")
fit3.rqs <- summary(fit3.rq, se = "nid")

bw <- USbirths2001$bw
bwm <- mean(bw)
bws <- sd(bw)
bw_term_single <- USbirths2001[USbirths2001$ga >= 37 & USbirths2001$multiple == 0,'bw']
bw_term_single <- as.numeric(na.omit(bw_term_single))
coefs.lm <- summary(fit1.lm)$coefficients
coefs.rq <- lapply(fit1.rqs, function(x) x$coefficients)
rdf <- sapply(fit1.rqs, function(x) x$rdf)
}

###################################################
# Save output
path.out <- "C:\\Work\\Research\\Projects\\Y03 - CRC\\Draft\\_results"

nn <- c(paste("fit",1:3,".lm",sep=""), paste("fit",1:3,".rq",sep=""), paste("fit",1:3,".rqs",sep=""))
save(list = nn, file = file.path(path.out, "ch1_results.Rdata"), compress = "gzip")

save(list = c("coefs.lm", "coefs.rq", "rdf", "bw_term_single"), file = file.path(path.out, "ch1_results_reduced.Rdata"))

###################################################
# Table 1 - Quantile effects at-term infants

load(file = file.path(path.out, "ch1_results_reduced.Rdata"))

std.m <- coefs.lm[,2]
std.q <- sapply(coefs.rq, function(x) x[,2])

est.q <- cbind(sapply(coefs.rq, function(x) x[,1]),std.q)[,c(1,4,2,5,3,6)]
est.m <- cbind(coefs.lm[,1], std.m)

nn <- c("Intercept", "Sex (female)", "Multiple births", "Maternal age", "Smoking (yes)", "Drinking (yes)")
tab1 <- data.frame(est.m,est.q)
colnames(tab1) <- rep(c("Est.","SE"),4)
rownames(tab1) <- nn

print(xtable(tab1, digits = 1, align = c("l",rep("r",8)), caption = c("Birth weight mean and quantile effects for infants born at term, United States, 2001. Standard errors are reported in brackets.\\label{tab:01_01}")), booktabs = TRUE, math.style.negative = TRUE)


###################################################
# Table 2 - Quantile effects preterm infants

load(file = file.path(path.out, "ch1_results.Rdata"))

std.m <- summary(fit2.lm)$coefficients[,2]
std.q <- sapply(fit2.rqs, function(x) x$coefficients[,2])

est.q <- cbind(coef(fit2.rq),std.q)[,c(1,4,2,5,3,6)]
est.m <- cbind(coef(fit2.lm), std.m)

tab2 <- data.frame(est.m,est.q)
colnames(tab2) <- rep(c("Est.","SE"),4)
rownames(tab2) <- nn

print(xtable(tab2, digits = 1, align = c("l",rep("r",8)), caption = c("Birth weight mean and quantile effects for infants born preterm, United States, 2001. Standard errors are reported in brackets.\\label{tab:01_02}")), booktabs = TRUE, math.style.negative = TRUE)

###################################################
# Table 3 - Quantile effects all infants

load(file = file.path(path.out, "ch1_results.Rdata"))

std.m <- summary(fit3.lm)$coefficients[,2]
std.q <- sapply(fit3.rqs, function(x) x$coefficients[,2])

est.q <- cbind(coef(fit3.rq),std.q)[,c(1,4,2,5,3,6)]
est.m <- cbind(coef(fit3.lm), std.m)

tab3 <- data.frame(est.m,est.q)
colnames(tab3) <- rep(c("Est.","SE"),4)
rownames(tab3) <- nn

print(xtable(tab3, digits = 1, align = c("l",rep("r",8)), caption = c("Birth weight mean and quantile effects for all infants, United States, 2001. Standard errors are reported in brackets.\\label{tab:01_03}")), booktabs = TRUE, math.style.negative = TRUE)

