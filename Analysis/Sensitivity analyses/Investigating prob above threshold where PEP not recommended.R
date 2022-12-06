# Investigating scenarios in which prob of exposure was above the threshold and
# PEP was not recommended
# Using final NASPHV survey results
# Aug 15, 2022
# Updated October 4 with NASPHV feedback
# Updated 18 November 2022

library(dplyr)

# Load in clean survey results data
dat <- readRDS("Results/Survey final results/nasphv_survey_final_res_18Nov2022.RDS")

# Drop missing (a few obs PR)
dat <-filter(dat, !(is.na(Median)))

threshold <- 0.0006 # Upper 95% CI limit

check <- filter(dat, Median > threshold & 
                  value != "Strongly agree" &
                  value != "Agree")

summary(check$Median)

table(check$value)

sort(table(check$scenario))

sort(table(check$State))

OK <- filter(check, State == "OK")

# Strongly disagree (format table for paper)
SD <- filter(check, value == "Strongly disagree")

table(SD$Animal)
table(SD$State)

# Separate out national est from state specific
state_sp <- filter(SD, Animal != "Bobcat" &
                     Animal != "Mesocarnivores")

nat <- filter(SD, Animal == "Bobcat" |
                Animal == "Mesocarnivores")

nat2 <- nat
nat2$State <- NA
nat2$region <- NA

nat2 <- distinct(nat2)

all <- rbind(nat2, state_sp)

all <- all[order(all$Median, decreasing = T),]

all$LL95CI <- signif(all$LL95CI, digits = 2)
all$Median <- signif(all$Median, digits = 2)
all$UL95CI <- signif(all$UL95CI, digits = 2)


##########
# Disagree
D <- filter(check, value == "Disagree")

# Separate out national est from state specific
state_sp <- filter(D, Animal == "Cat" |
                     Animal == "Dog" |
                     Animal == "Fox")

nat <- filter(D, Animal != "Cat" &
                Animal != "Dog" &
                Animal != "Fox")

nat2 <- nat
nat2$State <- NA
nat2$region <- NA

nat2 <- distinct(nat2)

all <- rbind(nat2, state_sp)

all <- all[order(all$Median, decreasing = T),]


