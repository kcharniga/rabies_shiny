# Run function to create df for probability of death from rabies exposure for BATS
# Aug 9, 2022
# Updated November 18, 2022 with corrected est. for provoked params

library(dplyr)
library(triangle)

set.seed(10) # Set random seed for reproducible analysis

# Number of times to sample from distributions
n_sims <- 1000000

#################################################
# Positivity rates for animal categories by state
#################################################
pos_rates <- readRDS("Data/Cleaned/bat_species_prob_size.RDS")

# Nice bat names
pos_rates$Species <- stringr::str_to_sentence(pos_rates$Species)
pos_rates$CommonName <- stringr::str_to_sentence(pos_rates$CommonName)

df <- pos_rates

###################################
# Parameter table for location and severity of exposure
#######################################

loc_severity <- readRDS("Data/params_location_severity.RDS")

# Do NOT adjust single bite to head/neck
loc_severity$mode[loc_severity$location_cat == "Head/neck" & 
                    loc_severity$severity == "Single bite"] <- 0.535

loc_severity$max[loc_severity$location_cat == "Head/neck" & 
                   loc_severity$severity == "Single bite"] <- 0.79

########################################
# Calculating probabilities for p2
########################################

# Prepare a list containing p2 
# p2: location or severity of exposure (samples from triangular or binomial distribution or
# just uses point estimate)
my_list_of_p2 <- list()
for (i in 1:nrow(loc_severity)){
  
  if (loc_severity$distribution[i] == "Triangular"){
    p2 <- rtriangle(n = n_sims, 
                    a = loc_severity$min[i], 
                    b = loc_severity$max[i], 
                    c = loc_severity$mode[i])
  } 
  if (loc_severity$distribution[i] == "Binomial") {
    p2 <- rbeta(n = n_sims, 
                shape1 = 1 + 0, 
                shape2 = 1 + loc_severity$size[i] - 0) # Using conjugate prior
    
  } 
  if (loc_severity$distribution[i] == "Point estimate") {
    p2 <- loc_severity$point_est[i]
  }
  my_list_of_p2[[i]] <- p2
  
}

p2_list <- my_list_of_p2

#########################################
# Create a df with all the combinations
#########################################

# Should be 30 * 4 * 28 = 3360
df$id <- 1:nrow(df)
my_list <- list()
for (i in 1:nrow(df)){
  bat_species <- filter(df, id == i)
  
    new_df <- data.frame(prob = rep(bat_species$prob,112), 
                         size = rep(bat_species$size,112),
                         common_name = rep(bat_species$CommonName,112),
                         species = rep(bat_species$Species,112),
                         location_cat = rep(loc_severity$location_cat, 4),
                         severity = rep(loc_severity$severity,4),
                         sick = c(rep(c("Yes"), 28*2), 
                                  rep(c("No"), 28*2)),
                         # Need to make sure we get all the combinations
                         # The order here matches the order in the function
                         # which is:
                         # Sick & unprovoked; sick & provoked; healthy & unprovoked;
                         # Healthy & provoked
                         unprovoked = c(rep("Yes",28),
                                        rep("No",28),
                                        rep("Yes",28),
                                        rep("No",28)))
  # Add df to list
  my_list[[i]] <- new_df
  
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
source("Analysis/Functions/Functions to calc death BATS.R")

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

######################

# Takes 3 min
start_time <- Sys.time()

res <- lapply(data_split, function(x) calc_prob_death(size = x$size, 
                                                      prob = x$prob,
                                                      n_sims = n_sims,
                                                      p2_list = p2_list,
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

# Add in other variables using the df structure created earlier in the script
final_df <- cbind(final_df, 
                  common_name = new_df$common_name,
                  species = new_df$species,
                  location_cat = new_df$location_cat,
                  severity = new_df$severity,
                  sick = new_df$sick,
                  unprovoked = new_df$unprovoked)

final_df <- data.frame(final_df)

# Re-order rows
final_df <- dplyr::select(final_df, common_name, species, location_cat, severity, 
                          sick, unprovoked, lower_95CI, point_est, upper_95CI)

df_final_tool2 <- final_df

# Formatting
df_final_tool2$lower_95CI <- as.numeric(df_final_tool2$lower_95CI)
df_final_tool2$point_est <- as.numeric(df_final_tool2$point_est)
df_final_tool2$upper_95CI <- as.numeric(df_final_tool2$upper_95CI)

# Clean up names a bit
df_final_tool2 <- rename(df_final_tool2,
                         Unprovoked = unprovoked,
                         Sick = sick,
                         Lower_95CI = lower_95CI,
                         Point_est = point_est,
                         Upper_95CI = upper_95CI)

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
  rename(Healthy = Sick, Provoked = Unprovoked,
         LL95CI = Lower_95CI, UL95CI = Upper_95CI, Median = Point_est, 
         CommonName = common_name,
         Species = species)

# Reorder and fix names
df_final_tool2 <- df_final_tool2 %>%
  rename(Severity = severity, Location_cat = location_cat) %>%
  select(Location_cat, Severity, CommonName, Species, Provoked, Healthy, LL95CI, 
         Median, UL95CI)

# Correct single bite to head/neck
#saveRDS(df_final_tool2, "final_df_for_tool_death_BATS_2011_2020_18Nov2022.RDS")
