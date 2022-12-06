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

###### Figure 1
# Plot colors
my_col_pal <- c("#ffd1d7", "#b28f81", "#76bae0", "#54483e")

ggplot(dat, aes(x = reorder(value, order), y = Median, color = region)) +
  geom_errorbar(aes(ymin = LL95CI, ymax = UL95CI), width = 0.15,
                position = position_dodge(width = 0.5)) +
  geom_point(size = 1.2, position = position_dodge(width = 0.5)) +
  scale_color_manual(values = my_col_pal) + # Use same colors as earlier plot
  theme_bw() +
  theme(axis.text.x = element_text(angle = 45, color = "black",
                                   vjust = 1, hjust = 1, size = 8),
        axis.title = element_text(color = "black", face = "bold",
                                  size = 12)) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(axis.title.x = element_blank()) +
  labs(y = "Pr(rabid|exposure)") +
  theme(axis.title.y = element_text(color = "black", size = 10, face = "bold")) +
  scale_y_continuous(limits = c(0,1),
                     breaks = c(0,0.2,0.4,0.6,0.8,1)) +
  labs(color = "") +
  ggtitle("Quantitative results for all surveyed states") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14))


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

######### Figure 2

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

# Predict the probability (p) of diabete positivity
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

