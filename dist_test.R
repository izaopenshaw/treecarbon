# Distributional tests
df0 <- df[!is.na(df$C), ]
Y <- df$SOC
Y <- Y[!is.na(Y)]

# Plot multiple densities
library(tidyr)
df  %>%
  gather() %>%
  ggplot(aes(value)) + 
  facet_wrap(~ key, scales = "free") +
  geom_density(color = "blue", fill = "red")

# Plot Histogram
library("ggplot2")
ggplot(df, aes(x=TotalVolume,color=sex, fill=sex)) +
  geom_density(alpha=.5) +
  geom_histogram(color="black", fill="white",bins = 50) +
  facet_grid(sex ~ .)

hist(Y)
hist(Y, breaks=50)
hist(log(Y))
hist(log(Y), breaks=50)

# Try scaling & plotting: log, log10, 1/, scale() = Z value, scaling to a range: (x-x_min)/(x_max-x_min) set min/max values to exclude outliers

#======Fit distr
library(fitdistrplus)
descdist(Y, boot=100)

#======Test of Normality
qqnorm(Y, main='Carbon') ; qqline(Y)
shapiro.test(Y) # if Y~N pval > 0.05
ks.test(Y, 'pnorm') # if Y~N pval > 0.05

# Using the KS test with estimated parameters:
ks.test(Y, 'pnorm', mean=mean(Y), sd=sd(Y)) # Y~N: pvalue > 0.05
ks.test(Y, 'pnorm') # Y~N: pvalue > 0.05
ks.test(scale(Y), 'pnorm') # if Y~N pval > 0.05
ks.test(log(Y), 'pnorm') # if Y~N pval > 0.05
# log normal
ks.test(log(Y), 'plnorm') # if Y~N pval > 0.05  


library(moments) #install.packages("moments")
skewness(Y) # Z~N(0,1) skewness = 0
skewness(log(Y)) # Z~N(0,1) skewness = 0
skewness(log10(Y)) # Z~N(0,1) skewness = 0
kurtosis(Y) # Z~N(0,1) kurtosis = 3
kurtosis(log(Y)) # Z~N(0,1) kurtosis = 3

#======Gamma
library(MASS) #"beta", "cauchy", "chi-squared", "exponential", "gamma", "geometric", "log-normal", "lognormal", "logistic", "negative binomial", "normal", "Poisson", "t" and "weibull" are recognised, case being ignored.
gam <- fitdistr(Y, "gamma")
gam 
ks.test(Y, "pgamma", gam$estimate[1], gam$estimate[2]) # two-sided, exact

#======Weibull
wei <- fitdistr(Y, "weibull")
wei 
ks.test(Y, "pweibull", wei$estimate[1], wei$estimate[2]) # two-sided, exact

#======Exponential
ks.test(Y, "pexp")

expo <- fitdistr(Y, "exponential")
ks.test(Y, "pexp", expo$estimate[1], expo$estimate[2]) # two-sided, exact


#======Beta: Beta distribution can be used to analyze probabilistic experiments that have only two possible outcomes:
# success, with probability p; failure, with probability 1-p. (Bernoulli experiments)
#            PDF   Prob as a     Models the 
# Binomial  f(x)   parameter    no. of successes, x    X~Bin(n,p)
# Beta      f(p)   random var   prob.(p) of success    X~Beta(alpha, beta)

beta <- fitdistr(Y, "beta", start=list(shape1= 2, shape2=5), method="mle")
beta <- fitdistr(Y, "beta", start=list(shape1= 2, shape2=5), method="mme")
ks.test(Y, "pexp", expo$estimate[1], expo$estimate[2]) # two-sided, exact

x <- seq(round(min(Y),0), round(max(Y),0), length.out=1000)
dat <- data.frame(x=x, px=dexp(x, rate=as.numeric(expo$estimate)))
ggplot(dat, aes(x=x, y=px)) + geom_line()


# How to deal with Skew
# Time Series Object
library(forecast)
count <- ts(Y)
fit <- auto.arima(count, seasonal=FALSE)
fit

tsdisplay(residuals(fit), lag.max=30)

# Identify duplicates
duplicates <- Y[duplicated(Y)]

# Identify & remove Outliers
par(mfrow=c(1,2))
boxplot(Y)
Y1 <- subset(Y,!(Y > quantile(Y, probs=c(.01,.995)))[2] | Y < quantile(Y, probs=c(.01,.995))[1] )
boxplot(Y1)


install.packages("qualityTools")
library("qualityTools")
# calculates the most likely parameters for a given distribution
distribution(x, distribution = "weibull", start, ...)

#calculate McFadden's R-squared for model
with(summary(model), 1 - deviance/null.deviance)

