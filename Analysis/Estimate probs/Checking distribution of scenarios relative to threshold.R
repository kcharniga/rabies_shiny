# Checking scenarios that fall above, below, and overlap the threshold
# April 18, 2023

library(dplyr)
library(ggplot2)
library(patchwork)

# Pr(exposure)
exposure <- readRDS("R Shiny/final_df_for_tool_exposure_vaxed_Bayes_2011_2020_18Nov2022.RDS")

threshold <- 0.0004

below <- filter(exposure, UL95CI < threshold)
above <- filter(exposure, LL95CI > threshold)

overlap <- nrow(exposure) - nrow(below) - nrow(above)

nrow(above)/nrow(exposure)
nrow(below)/nrow(exposure)
overlap/nrow(exposure)
