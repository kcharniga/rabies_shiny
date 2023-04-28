# Predicting the threshold for when PEP is recommended
# using final NASPHV survey results
# Aug 8, 2022
# Updated October 3, 2022 with feedback from NASPHV related to terrestrial rabies
# Updated 18 Nov 2022

library(dplyr)
library(jtools)
library(ggplot2)
library(visreg)

# Load in clean survey results data
dat <- readRDS("Results/Survey final results/nasphv_survey_final_res_18Nov2022.RDS")

# Drop missing (a few obs PR)
dat <-filter(dat, !(is.na(Median)))

# Try to predict whether PEP will be recommended using logistic regression first
dat$PEP_yes <- NA
dat$PEP_yes[dat$value == "Strongly agree"] <- 1
dat$PEP_yes[dat$value == "Agree"] <- 1
dat$PEP_yes[dat$value == "Disagree"] <- 0
dat$PEP_yes[dat$value == "Strongly disagree"] <- 0


# Try log transforming probabilities
dat$log_median <- log(dat$Median)

mylogit <- glm(PEP_yes ~ log_median, 
               data = dat, 
               family = "binomial")

summary(mylogit)

# ORs and 95% CIs
summ(mylogit, digits = 5, exp = TRUE)

# Visualize logistic regression output
effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            jitter = c(0.1, 0.05), point.alpha = 0.1, interval = TRUE) +
  ylab("Pr(PEP recommended = Yes)") +
  xlab("Median Pr(rabid|exposure)") 


predict(mylogit, data.frame(log_median=-5), type="response")

# Inflection point
# https://randomproblemsicomeacross.blogspot.com/2013/04/finding-inflection-point-in-logistic.html

# Equation ln(p/(1-p)) = a + bx

coef(mylogit)

p <- 0.5
x_int <- (log(p/(1-p)) - coef(mylogit)[1]) / coef(mylogit)[2]
x_int
exp(x_int) # 0.0004

library(scales)
scaleFUN <- function(x) {sprintf("%.1g", exp(x))}

effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            point.alpha = 0.1, interval = F, jitter = c(0, 0.1)) + 
  ylab("PEP yes") +
  xlab("Pr(rabid|exposure)") +
  geom_vline(xintercept = x_int, linetype = "dashed") +
  scale_x_continuous(breaks=c(-11.5,x_int,-4.5,-2,0), 
                     labels = scaleFUN) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.title = element_text(color = "black")) 




# CIs?
## grad the inverse link function
ilink <- family(mylogit)$linkinv
## add fit and se.fit on the **link** scale
dat <- bind_cols(dat, setNames(as_tibble(predict(mylogit, dat, se.fit = TRUE)[1:2]),
                               c('fit_link','se_link')))
## create the interval and backtransform
dat <- mutate(dat,
              fit_resp  = ilink(fit_link),
              right_upr = ilink(fit_link + (2 * se_link)),
              right_lwr = ilink(fit_link - (2 * se_link)))

# Check that this gives me the confidence intervals
effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            jitter = c(0.1, 0.05), point.alpha = 0.1, interval = TRUE) +
  ylab("Pr(PEP recommended = Yes)") +
  xlab("Median Pr(rabid|exposure)") +
  geom_ribbon(data = dat,
              aes(ymin = right_lwr, ymax = right_upr),
              alpha = 0.1, color = "blue")

summary(dat$Median)

# Use these predictions and R's spatial tools to get intersections for the CIs
library(sp)
library(rgeos)  ## for gIntersection()

# Upper limit
## Wrap your line up as a SpatialLines object
x <- dat$log_median
y <- dat$right_upr
SL1 <- SpatialLines(list(Lines(Line(cbind(x,y)), "A")))

## Create a horizontal SpatialLines object
SL2 <- SpatialLines(list(Lines(Line(cbind(range(x), 0.5)), "B")))

## Find their point(s) of intersection
coordinates(gIntersection(SL1, SL2))

#length(coordinates(gIntersection(SL1, SL2)))/2 # 46

x_int1 <- coordinates(gIntersection(SL1, SL2))[36,1] # This is trial and error (eye balling it)
exp(x_int1)

# Lower limit
## Wrap your line up as a SpatialLines object
x <- dat$log_median
y <- dat$right_lwr
SL1 <- SpatialLines(list(Lines(Line(cbind(x,y)), "A")))

## Create a horizontal SpatialLines object
SL2 <- SpatialLines(list(Lines(Line(cbind(range(x), 0.5)), "B")))

## Find their point(s) of intersection
coordinates(gIntersection(SL1, SL2))
length(coordinates(gIntersection(SL1, SL2)))/2

x_int2 <- coordinates(gIntersection(SL1, SL2))[1,1]
exp(x_int2)

library(scales)
scaleFUN <- function(x) {sprintf("%.1g", exp(x))}

effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            point.alpha = 0.1, interval = TRUE, jitter = c(0, 0.1)) + 
  ylab("PEP yes") +
  xlab("Pr(rabid|exposure)") +
  geom_vline(xintercept = x_int1, linetype = "dashed") +
  geom_vline(xintercept = x_int2, linetype = "dashed") +
  geom_vline(xintercept = x_int, linetype = "dashed") +
  #geom_hline(yintercept = 0.5) 
  scale_x_continuous(breaks=c(-11.5,-6.9,-4.5,-2,0), 
                     labels = scaleFUN) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.title = element_text(color = "black")) 
#labels=c(exp(-10), exp(-5), exp(-2), exp(0)))

# Threshold estimate
exp(x_int)

# 95% CI
exp(x_int1)
exp(x_int2)

# Predict the probability (p) 
probabilities <- predict(mylogit, type = "response")
predicted.classes <- ifelse(probabilities > 0.5, "pos", "neg")
head(predicted.classes)

# Check logistic regression assumptions

# Linearity assumption

# Select only numeric predictors
library(tidyverse)
mydata <- dat %>%
  dplyr::select(log_median) 
predictors <- colnames(mydata)
# Bind the logit and tidying the data for plot
mydata <- mydata %>%
  mutate(logit = log(probabilities/(1-probabilities))) %>%
  gather(key = "predictors", value = "predictor.value", -logit)

p1 <- ggplot(mydata, aes(logit, predictor.value))+
  geom_point(size = 0.5, alpha = 0.5) +
  geom_smooth(method = "loess") + 
  theme_bw() + 
  facet_wrap(~predictors, scales = "free_y")

# Check influential points
plot(mylogit, which = 4, id.n = 3)

# computes the standardized residuals (.std.resid) and the Cook’s distance (.cooksd)

# Extract model results
model.data <- broom::augment(mylogit) %>% 
  mutate(index = 1:n()) 

check <- filter(model.data, .cooksd > 4/nrow(dat))
check2 <- filter(model.data, abs(.std.resid) > 3)

model.data %>% top_n(3, .cooksd)

ggplot(model.data, aes(index, .std.resid)) + 
  geom_point(alpha = .5) +
  theme_bw()

model.data %>% 
  filter(abs(.std.resid) > 3)

# Put plots together to save

predicted <- predict(mylogit, type="response")

# Getting log odds values
log_odds = log(predicted / (1 - predicted))

par(mfrow = c(2, 1))
plot(x = log_odds, y = dat$log_median, xlab = "Logit", ylab = "Predictor value")
mtext(text = "A", col = "black", cex = 1.2, side = 3, adj = 0, line = 2, font = (face=2))
plot(mylogit, which = 4, id.n = 3)
#abline(h = 4/nrow(dat), lty = 2)
mtext(text = "B", col = "black", cex = 1.2, side = 3, adj = 0, line = 2, font = (face=2))

# save as 650 x 750 jpeg

################
# Remove influential points and run logistic regression
################
###
##
#

dat_no_outliers <- filter(model.data, .cooksd < 4/nrow(dat))

dat <- dat_no_outliers

mylogit <- glm(PEP_yes ~ log_median, 
               data = dat, 
               family = "binomial")

summary(mylogit)

# ORs and 95% CIs
summ(mylogit, digits = 5, exp = TRUE)

# Visualize logistic regression output
effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            jitter = c(0.1, 0.05), point.alpha = 0.1, interval = TRUE) +
  ylab("Pr(PEP recommended = Yes)") +
  xlab("Median Pr(rabid|exposure)") 


predict(mylogit, data.frame(log_median=-5), type="response")

# Inflection point
# https://randomproblemsicomeacross.blogspot.com/2013/04/finding-inflection-point-in-logistic.html

# Equation ln(p/(1-p)) = a + bx

coef(mylogit)

p <- 0.5
x_int <- (log(p/(1-p)) - coef(mylogit)[1]) / coef(mylogit)[2]
x_int
exp(x_int) # 0.0004

# CIs?
## grad the inverse link function
ilink <- family(mylogit)$linkinv
## add fit and se.fit on the **link** scale
dat <- bind_cols(dat, setNames(as_tibble(predict(mylogit, dat, se.fit = TRUE)[1:2]),
                               c('fit_link','se_link')))
## create the interval and backtransform
dat <- mutate(dat,
              fit_resp  = ilink(fit_link),
              right_upr = ilink(fit_link + (2 * se_link)),
              right_lwr = ilink(fit_link - (2 * se_link)))

# Check that this gives me the confidence intervals
effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            jitter = c(0.1, 0.05), point.alpha = 0.1, interval = TRUE) +
  ylab("Pr(PEP recommended = Yes)") +
  xlab("Median Pr(rabid|exposure)") +
  geom_ribbon(data = dat,
              aes(ymin = right_lwr, ymax = right_upr),
              alpha = 0.1, color = "blue")

summary(dat$Median)

# Use these predictions and R's spatial tools to get intersections for the CIs
library(sp)
library(rgeos)  ## for gIntersection()

# Upper limit
## Wrap your line up as a SpatialLines object
x <- dat$log_median
y <- dat$right_upr
SL1 <- SpatialLines(list(Lines(Line(cbind(x,y)), "A")))

## Create a horizontal SpatialLines object
SL2 <- SpatialLines(list(Lines(Line(cbind(range(x), 0.5)), "B")))

## Find their point(s) of intersection
coordinates(gIntersection(SL1, SL2))

#length(coordinates(gIntersection(SL1, SL2)))/2 # 46

x_int1 <- coordinates(gIntersection(SL1, SL2))[40,1] # This is trial and error (eye balling it)
exp(x_int1)

# Lower limit
## Wrap your line up as a SpatialLines object
x <- dat$log_median
y <- dat$right_lwr
SL1 <- SpatialLines(list(Lines(Line(cbind(x,y)), "A")))

## Create a horizontal SpatialLines object
SL2 <- SpatialLines(list(Lines(Line(cbind(range(x), 0.5)), "B")))

## Find their point(s) of intersection
coordinates(gIntersection(SL1, SL2))
length(coordinates(gIntersection(SL1, SL2)))/2

x_int2 <- coordinates(gIntersection(SL1, SL2))[18,1]
exp(x_int2)

library(scales)
scaleFUN <- function(x) {sprintf("%.1g", exp(x))}

effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            point.alpha = 0.1, interval = TRUE, jitter = c(0, 0.1)) + 
  ylab("PEP yes") +
  xlab("Pr(rabid|exposure)") +
  geom_vline(xintercept = x_int1, linetype = "dashed") +
  geom_vline(xintercept = x_int2, linetype = "dashed") +
  geom_vline(xintercept = x_int, linetype = "dashed") +
  #geom_hline(yintercept = 0.5) 
  scale_x_continuous(breaks=c(-11.5,-6.9,-4.5,-2,0), 
                     labels = scaleFUN) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme(axis.title = element_text(color = "black")) 
#labels=c(exp(-10), exp(-5), exp(-2), exp(0)))

# Threshold estimate
exp(x_int)

# 95% CI
exp(x_int1)
exp(x_int2)







############# Try to fit a logistic GAM model
# https://kaminsky.rocks/2020/04/notes-on-gams-in-r-with-a-binary-dependent-variables/

library(mgcv)

mod_gam <- gam(PEP_yes ~ s(log_median), 
               data = dat, 
               family = "binomial",
               method = "REML")

summary(mod_gam)

# Convert output to the probability scale and plot (partial effect with no intercept)
labs <- c(signif(exp(-12),1), 
          signif(exp(-8),1), 
          signif(exp(-4),1), 
          signif(exp(0),1))

plot(mod_gam, pages = 1, trans = plogis, main = "Logistic GAM", 
     xlab = "Median Pr(rabid|exposed)", xaxt = "n")
axis(1, at = c(-12, -8, -4, 0),
     labels = labs)
abline(h = 0.5, col = "blue")
abline(v = -4, col = "green")
abline(v = -4.45, col = "green")
abline(v = -3.6, col = "green")

exp(-4)
exp(-4.45)
exp(-3.6)

# Incorporate intercept
# I don't think this is necessary if there's only one variable in the model
# Each partial effect plot can be interpreted as showing the probability of the 
# outcome if all other variables were at their average value. At their own average 
# value, you get only the effect of the intercept.
# https://noamross.github.io/gams-in-r-course/chapter4
# plot(mod_gam, pages = 1, trans = plogis,
#      shift = coef(mod_gam)[1], seWithMean = TRUE,
#      main = "Logistic GAM", 
#      xlab = "Median Pr(rabid|exposed)", xaxt = "n" )
# axis(1, at = c(-12, -8, -4, 0),
#      labels = labs)
# abline(h = 0.5, col = "blue")
# abline(v = -6.3, col = "green")
# abline(v = -5.6, col = "green")
# abline(v = -7.35, col = "green")
# exp(-6.3)
# exp(-5.6)
# exp(-7.35)
# Model diagnostics, model converges
par(mfrow=c(2,2))
gam.check(mod_gam)

AIC(mod_gam) # GAM has lower AIC than logistic regression model

summary(mod_gam)$r.sq

visreg(mod_gam, gg = TRUE)

# Try changing basis dimensions, k
mod_gam2 <- gam(PEP_yes ~ s(log_median, k = 20), 
               data = dat, 
               family = "binomial",
               method = "REML")

summary(mod_gam2)
gam.check(mod_gam2)
AIC(mod_gam2)

plot(mod_gam2, pages = 1, trans = plogis, main = "Logistic GAM", 
     xlab = "Median Pr(rabid|exposed)", xaxt = "n")
axis(1, at = c(-12, -8, -4, 0),
     labels = labs)
abline(h = 0.5, col = "blue")

plot.gam(mod_gam2)

#############################
# Ordinal logistic regression

# Used to predict an ordinal dependent variable given one or
# more independent variables
# https://www.st-andrews.ac.uk/media/ceed/students/mathssupport/OrdinalexampleR.pdf
library(MASS)

dat$value_rank <- NA
dat$value_rank[dat$value == "Strongly agree"] <- 1
dat$value_rank[dat$value == "Agree"] <- 2
dat$value_rank[dat$value == "Disagree"] <- 3
dat$value_rank[dat$value == "Strongly disagree"] <- 4

dat$value_rank <- as.factor(dat$value_rank)

mod <- polr(value_rank ~ log_median,
            data=dat, Hess=T)

summary(mod)
mod$zeta

coef(mod)

# Visualize the thresholds
#https://stats.stackexchange.com/questions/89474/interpretation-of-ordinal-logistic-regression

mpg   <- seq(from=min(dat$log_median), to=max(dat$log_median), by=0.0293)
xbeta <- mpg*(coef(mod))
logistic_cdf <- function(x) {
  return( 1/(1+exp(-x) ) )
}

p1 <- logistic_cdf( mod$zeta[1] - xbeta )
p2 <- logistic_cdf( mod$zeta[2] - xbeta ) - logistic_cdf( mod$zeta[1] - xbeta )
p3 <- logistic_cdf( mod$zeta[3] - xbeta ) - logistic_cdf( mod$zeta[2] - xbeta )
p4 <- 1 - logistic_cdf( mod$zeta[3] - xbeta )

labs <- c(signif(exp(-12),1), 
          signif(exp(-8),1), 
          signif(exp(-4),1), 
          signif(exp(0),1))

plot(mpg, p1, type='l', ylab='Probability PEP recommended', ylim = c(0,1),
     xlab = "Median Pr(rabid|exposed)", main = "Ordinal logistic regression",
     xaxt = "n")
lines(mpg, p2, col='red')
lines(mpg, p3, col='blue')
lines(mpg, p4, col='green')
legend("topright", lty=1, col=c("black", "red", "blue", "green"),
       legend=c("Strongly agree", "Agree", "Disagree", "Strongly disagree"))
abline(v = -7.76, lty = "dashed")
axis(1, at = c(-12, -8, -4, 0),
     labels = labs) 


# What is the intersection of agree and disagree?
x <- mpg
y1 <- p2
y2 <- p3

# Plot the first curve
plot(x, y1, type = "l", xlab = "V1", ylab = "V2", col = "blue")
# Add the second
lines(x, y2, col = "red")

equivalent <- function(x, y, tol = 0.005) abs(x - y) < tol
xmin <- -8
xmax <- -6
intersection_indices <- which(equivalent(y1, y2) & x >= xmin & x <= xmax)
x[intersection_indices]
#> [1] 3.93 7.07
points(x[intersection_indices], y1[intersection_indices])





# Interpretation:standard interpretation of the ordered log-odds coefficient 
# is that for a one unit increase in the predictor, the response variable 
# level is expected to change by its respective regression coefficient in the 
# ordered log-odds scale.

# i.e. with 1% increase in median Pr(rabid|exposure), the log of odds of not
# recommending PEP decreases by 0.04?

# Calculate p values
coeffs <- coef(summary(mod))
p <- pnorm(abs(coeffs[, "t value"]), lower.tail = FALSE) * 2
cbind(coeffs, "p value" = round(p,4))

# Interpretation: test statistics and p-values, respectively for the 
# null hypothesis that an individual predictor’s regression coefficient 
# is zero given that the rest of the predictors are in the model

# ORs
exp(coef(mod))

# Predict probability that PEP will be recommended given Pr(rabid|exposure)
new_data <- data.frame("Median_100"= 1)

round(predict(mod, new_data, type = "p"), 3)

