# Interpreting the probabilities in the risk assessment tool
# 18-2-22
# Updated Aug 8, 2022 with 2011-2020 data
# Updated Oct 3, 2022 according to feedback from NASPHV about terrestrial rabies
# Updated Nov 18, 2022 with corrected params for provoked

library(dplyr)
library(ggplot2)
library(patchwork)

# Pr(exposure)
exposure <- readRDS("R Shiny/final_df_for_tool_exposure_vaxed_Bayes_2011_2020_18Nov2022.RDS")

# Pr(death)
death <- readRDS("R Shiny/final_df_for_tool_death_vaxed_Bayes_2011_2020_18Nov2022.RDS")

# Exposure
summary(exposure$LL95CI)
summary(exposure$Median)
summary(exposure$UL95CI)

hist(exposure$Median, breaks = 30, col = "coral", main = "", xlab = "Pr(exposure)")

# Riskiest exposures - mostly unprovoked, ill raccoons and skunks as well as mongoose in PR
risky_exp <- filter(exposure, Median > 0.8)

coons <- filter(exposure, Animal == "Raccoon" &
                  Provoked == "Unprovoked" &
                  Healthy == "Ill or acting strangely")

# Least risky exposures - provoked, apparently healthy dogs and cats, followed closely by pocket pets
not_risky_exp <- filter(exposure, Median < 0.00001)


# Death
summary(death$LL95CI)
summary(death$Median)
summary(death$UL95CI)

# Plot both histograms together
par(mfrow=c(1,2))
hist(exposure$Median, breaks = 30, col = "coral", main = "", xlab = "Pr(exposure)")
hist(death$Median, breaks = 30, col = "lightblue", main = "", xlab = "Pr(death)")

# Box an whisker plots
death2 <- select(death, Median)
exposure2 <- select(exposure, Median)

death2$prob <- "Pr(death)"
exposure2$prob <- "Pr(exposure)"

dat <- rbind(death2, exposure2)

#par(mfrow=c(2,2))
#boxplot(Median ~ prob, dat, col = "#82C09A", ylab = "Probability", xlab = "")
#boxplot(Median ~ prob, dat, col = "#82C09A", ylim = c(0,0.08), ylab = "Probability", xlab = "")
#boxplot(Median ~ prob, dat, col = "#82C09A", ylim = c(0,0.005), ylab = "Probability", xlab = "")

# Still skewed and ugly
# Try violin plot on log scale

library(ggplot2)

p <- ggplot(data = dat, aes(x = prob, y = Median)) + 
  geom_violin(fill = '#82C09A', draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab(NULL) + 
  ylab('Probability')+
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
p # Don't save this one

# Need to transform data to deal with zeros
dat$Median_trans <- dat$Median + 0.0000000001

p <- ggplot(data = dat, aes(x = prob, y = log10(Median_trans))) + 
  geom_violin(fill = '#82C09A', draw_quantiles = c(0.25, 0.5, 0.75)) + 
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab(NULL) + 
  ylab('log10(probability)')+
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
p

# Don't use these plots because the data include a bunch of duplicates across states for national estimates

# Summary statistics for tool excluding duplicates for "all" states 
# (non-state-specific estimates)
library(dplyr)

# Pr(exposure)
exposure <- readRDS("R Shiny/df_for_tool_exposure_not_expanded_vaxed_Bayes_2011_2020_18Nov2022.RDS")

# Pr(death)
death <- readRDS("R Shiny/df_for_tool_death_not_expanded_vaxed_Bayes_2011_2020_18Nov2022.RDS")

# Exposure
summary(exposure$LL95CI)
summary(exposure$Median)
summary(exposure$UL95CI)

hist(exposure$Median, breaks = 30, col = "coral", main = "", xlab = "Pr(exposure)")

# Riskiest exposures - mostly unprovoked, ill raccoons and skunks as well as mongoose in PR
risky_exp <- filter(exposure, Median > 0.8)

# Least risky exposures - provoked, apparently healthy dogs and cats, followed closely by pocket pets
not_risky_exp <- filter(exposure, Median < 0.00001)


# Death
# First remove lick to intact skin
death <- filter(death, Severity != "Lick to intact skin")
summary(death$Median)

# Plot both histograms together
par(mfrow=c(1,2))
hist(exposure$Median, breaks = 30, col = "coral", main = "", xlab = "Pr(exposure)")
hist(death$Median, breaks = 30, col = "lightblue", main = "", xlab = "Pr(death)")

death2 <- select(death, Median)
exposure2 <- select(exposure, Median)

death2$prob <- "Pr(death)"
exposure2$prob <- "Pr(exposure)"

dat <- rbind(death2, exposure2)

# Try violin plot on log scale

library(ggplot2)

dat$prob[dat$prob == "Pr(death)"] <- "Pr(death|exposure)"
dat$prob[dat$prob == "Pr(exposure)"] <- "Pr(rabid|exposure)"
  

p <- ggplot(data = dat, aes(x = prob, y = Median)) + 
  geom_violin(fill = '#82C09A', draw_quantiles = c(0.25, 0.5, 0.75)) + 
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  xlab(NULL) + 
  ylab('Probability')+
  theme_bw() +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 14))
p

# NNT plots

p1 <- ggplot(death, aes(x = Median, y = 1/Median)) +
  geom_point(color = "#FF0F80") +
  theme_classic() +
  xlab("Pr(death|exposure)")+
  ylab("NNT") +
  scale_y_log10(labels = scales::number_format(big.mark = ",")) #+
  #geom_hline(yintercept = (32200*37)/3800, linetype = "dashed") # cost effectiveness threshold, 1/0.0032

p1

p2 <- ggplot(death, aes(x = Median, y = 3800*(1/Median))) + #Pieracci, E. et al. Vital Signs: Trends in Human Rabies Deaths and Exposures — United States, 1938–2018. MMWR Morb Mortal Wkly Rep 68, 524-528, doi:https://dx.doi.org/10.15585%2Fmmwr.mm6823e1 (2019).
  geom_point(color = "#FF0F80") +
  theme_classic() +
  xlab("Pr(death|exposure)")+
  ylab("PEP cost per death averted ($)") +
  scale_y_log10(labels = scales::number_format(big.mark = ",")) #+
  #geom_hline(yintercept = 32200*37, linetype = "dashed") # cost effectiveness threshold

p2

p1 + p2

# save 
