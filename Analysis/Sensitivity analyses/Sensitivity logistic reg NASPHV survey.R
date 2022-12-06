# Predicting the threshold for when PEP is recommended
# using final NASPHV survey results
# Sensitivity where only using respondents with at least 9 years of experience
# Aug 16, 2022
# Updated October 4, 2022 following NASPHV feedback
# Updated 18 November 2022

library(dplyr)
library(jtools)
library(ggplot2)
library(visreg)

# Load in clean survey results data
dat <- readRDS("Results/Survey final results/nasphv_survey_sensitivity_18Nov2022.RDS")

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
exp(x_int) # 0.003

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

x_int1 <- coordinates(gIntersection(SL1, SL2))[26,1] # This is trial and error (eye balling it)
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

x_int2 <- coordinates(gIntersection(SL1, SL2))[5,1]
exp(x_int2)

library(scales)
scaleFUN <- function(x) {sprintf("%.1g", exp(x))}

effect_plot(mylogit, pred = log_median, plot.points = TRUE,
            point.alpha = 0.1, interval = TRUE) + #jitter = c(0, 0.05), 
  ylab("Pr(PEP recommended = Yes)") +
  xlab("Median Pr(rabid|exposure)") +
  geom_vline(xintercept = x_int1, linetype = "dashed") +
  geom_vline(xintercept = x_int2, linetype = "dashed") +
  #geom_hline(yintercept = 0.5) 
  scale_x_continuous(breaks=c(-11.5,-6.9,-4.5,-2,0), 
                     labels = scaleFUN)
#labels=c(exp(-10), exp(-5), exp(-2), exp(0)))

# Threshold estimate
exp(x_int)

# 95% CI
exp(x_int1)
exp(x_int2)

