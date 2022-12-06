# Run function to create df for probability of exposure for BATS ONLY
# Aug 8, 2022
# Updated November 18, 2022 with corrected est. for provoked params

library(dplyr)
library(triangle)

set.seed(10) # Set random seed for reproducible analysis

#################################################
# Positivity rates for bat species by state
#################################################
pos_rates <- readRDS("Data/Cleaned/bat_species_prob_size.RDS")

# Nice bat names
pos_rates$Species <- stringr::str_to_sentence(pos_rates$Species)
pos_rates$CommonName <- stringr::str_to_sentence(pos_rates$CommonName)

df <- pos_rates

#########################################
# Create a df with all the combinations
#########################################

# Should be 30*4 = 120 combinations or rows
df$id <- 1:nrow(df)
my_list <- list()
for (i in 1:nrow(df)){
  bat_species <- filter(df, id == i)
  
    df_new <- data.frame(prob = rep(bat_species$prob,4), 
                         size = rep(bat_species$size,4),
                         common_name = rep(bat_species$CommonName,4),
                         species = rep(bat_species$Species,4),
                         sick = c("Yes","Yes","No","No"),
                         unprovoked = c("Yes","No","Yes","No"))
    
  # Add df to list
  my_list[[i]] <- df_new
  
}

# Combine dfs
new_df <- do.call("rbind", my_list)

###############################
# Use functions to calculate point est, LL, and UL of 95% CI
##############################################################

# Split data by ID for animal
ids_df <- select(df, CommonName, id) %>%
  rename(common_name = CommonName)

new_df <- left_join(new_df, ids_df, by = "common_name")

data_split <- split(new_df, ~ id)
x <- data_split[[1]]

# Use a function to calculate probabilities (point est and 95% CI)
source("Analysis/Functions/Functions to create df for BATS.R")

n_sims <- 1000000

# Probabilities for sick
prob_s_given_r <- 39/45
prob_s_given_not_r <- 721/3077

# Probabilities for healthy
prob_h_given_r <- 6/45
prob_h_given_not_r <- 2356/3077

# Probabilities for unprovoked
prob_unprov_given_r <-27/52
prob_unprov_given_not_r <- 111/444 # corrected
# Probabilities for provoked
prob_prov_given_r <- 25/52
prob_prov_given_not_r <- 333/444 # corrected

# Takes < 1 min
start_time <- Sys.time()

res <- lapply(data_split, function(x) calc_prob_new(size = x$size, 
                                                    prob = x$prob,
                                                    n_sims = n_sims,
                                                  
                                                    # Sick probabilities
                                                    prob_s_given_r = prob_s_given_r,
                                                    prob_s_given_not_r = prob_s_given_not_r,
                                                    # Healthy probabilities
                                                    prob_h_given_r = prob_h_given_r,
                                                    prob_h_given_not_r = prob_h_given_not_r,
                                                    # Unprovoked
                                                    prob_unprov_given_r = prob_unprov_given_r,
                                                    prob_unprov_given_not_r = prob_unprov_given_not_r,
                                                    # Provoked
                                                    prob_prov_given_r = prob_prov_given_r,
                                                    prob_prov_given_not_r = prob_prov_given_not_r))

end_time <- Sys.time()
end_time - start_time

# Combine dfs
final_df <- do.call("rbind", res)

# Add in other variables
final_df <- cbind(final_df, 
                  common_name = new_df$common_name,
                  species = new_df$species,
                  sick = new_df$sick,
                  unprovoked = new_df$unprovoked)

final_df <- data.frame(final_df)

# Re-order rows
final_df <- dplyr::select(final_df, common_name, species,
                          sick, unprovoked, lower_95CI, point_est, upper_95CI)

################################################################
df_final_tool2 <- final_df

df_final_tool2$lower_95CI <- as.numeric(df_final_tool2$lower_95CI)
df_final_tool2$point_est <- as.numeric(df_final_tool2$point_est)
df_final_tool2$upper_95CI <- as.numeric(df_final_tool2$upper_95CI)

# Clean up names a bit
df_final_tool2 <- rename(df_final_tool2,
                         Unprovoked = unprovoked,
                         Sick = sick,
                         LL_95CI = lower_95CI,
                         Median = point_est,
                         UL_95CI = upper_95CI)

# A bit more cleaning 
df_final_tool2$Unprovoked[df_final_tool2$Unprovoked == "Yes"] <- "Unprovoked"
df_final_tool2$Unprovoked[df_final_tool2$Unprovoked == "No"] <- "Provoked"

df_final_tool2$Sick[df_final_tool2$Sick == "No"] <- "Apparently healthy"
df_final_tool2$Sick[df_final_tool2$Sick == "Yes"] <- "Ill or acting strangely"

# Convert things to factors
df_final_tool2$common_name <- as.factor(df_final_tool2$common_name)
df_final_tool2$species <- as.factor(df_final_tool2$species)
df_final_tool2$Unprovoked <- as.factor(df_final_tool2$Unprovoked)
df_final_tool2$Sick <- as.factor(df_final_tool2$Sick)

# Change names
df_final_tool2 <- df_final_tool2 %>%
  rename(Healthy = Sick, Provoked = Unprovoked, CommonName = common_name,
         Species = species)

df_final_tool2 <- df_final_tool2 %>%
  rename(LL95CI = LL_95CI, UL95CI = UL_95CI)

# Reorder variables
df_final_tool2 <- df_final_tool2 %>%
  select(CommonName, Species, Provoked, Healthy, LL95CI, Median, UL95CI)

# Save final df for tool
#saveRDS(df_final_tool2, "final_df_for_tool_exposure_BATS_2011_2020_18Nov2022.RDS")
