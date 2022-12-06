# Run function to create df for conditional probability of exposure
# 25-2-22
# Uses new method for adjusting positivity rates by healthy/provoked
# updated 30 Sept 2022 to stratify pos rates for certain animals based on whether state
# has terrestrial rabies
# Updated November 18, 2022 with corrected est. for provoked params


library(dplyr)
library(triangle)

set.seed(10) # Set random seed for reproducible analysis

#################################################
# Positivity rates for animal categories by state
#################################################
pos_rates <- readRDS("Data/Cleaned/prob_size_all_animals_2011_2020.RDS")

# Nice animal names
pos_rates$Animal <- stringr::str_to_title(pos_rates$Animal)
sort(unique(pos_rates$Animal))
pos_rates$Animal[pos_rates$Animal == "Arctic Fox"] <- "Arctic fox"
pos_rates$Animal[pos_rates$Animal == "Armadillo And Opossum"] <- "Armadillo and opossum"
pos_rates$Animal[pos_rates$Animal == "Large Carnivores"] <- "Large carnivores"
pos_rates$Animal[pos_rates$Animal == "Meso Carnivores"] <- "Mesocarnivores"
pos_rates$Animal[pos_rates$Animal == "Non-Native Wild"] <- "Non-native wild"
pos_rates$Animal[pos_rates$Animal == "Pocket Pets"] <- "Pocket pets"
pos_rates$Animal[pos_rates$Animal == "Rodents Of Unusual Size"] <- "Rodents of unusual size"
pos_rates$Animal[pos_rates$Animal == "Sheep And Goats"] <- "Sheep and goats"
pos_rates$Animal[pos_rates$Animal == "Small Rodents And Lagomorphs"] <- "Small rodents, lagomorphs, and Eulipotyphla"
pos_rates$Animal[pos_rates$Animal == "Unknown Or Not Provided"] <- "Unknown or not provided"
pos_rates$Animal[pos_rates$Animal == "Hoofed"] <- "Hoofed animals"

pos_rates$Animal[pos_rates$Animal == "Rodents Of Unusual Size Terr"] <- "Rodents of unusual size terr"
pos_rates$Animal[pos_rates$Animal == "Sheep And Goats Terr"] <- "Sheep and goats terr"
pos_rates$Animal[pos_rates$Animal == "Hoofed Terr"] <- "Hoofed animals terr"
pos_rates$Animal[pos_rates$Animal == "Cattle Terr"] <- "Cattle terr"
pos_rates$Animal[pos_rates$Animal == "Equine Terr"] <- "Equine terr"

df_for_function <- pos_rates

# Fix State for Mongoose and Arctic fox
df_for_function$State[df_for_function$Animal == "Mongoose"] <- "PR"
df_for_function$State[df_for_function$Animal == "Arctic fox"] <- "AK"

# Remove Nans (come back to these! Also expand out "All" to the different
# states later) ##########################
nans <- filter(df_for_function, is.na(prob))

df <- filter(df_for_function, !(is.na(prob)))

# Create indicator variable for whether RABV vaccine is licensed
# in the animal
df$vacc_cat <- 0
df$vacc_cat[df$Animal == "Cat"|
              df$Animal == "Dog"|
              df$Animal == "Ferret"|
              df$Animal == "Equine"|
              df$Animal == "Sheep and goats"|
              df$Animal == "Cattle"|
              df$Animal == "Equine terr"|
              df$Animal == "Sheep and goats terr"|
              df$Animal == "Cattle terr"] <- 1

licensed_vacc_cats <- filter(df, vacc_cat == 1)
table(licensed_vacc_cats$Animal)

#########################################
# Create a df with all the combinations
#########################################

df$id <- 1:nrow(df)
my_list <- list()
for (i in 1:nrow(df)){
  animal_state <- filter(df, id == i)
  
  if (animal_state$vacc_cat == 0){ # The order here matches the order in the function
  df_new <- data.frame(prob = rep(animal_state$prob,4), 
                       size = rep(animal_state$size,4),
                       State = rep(animal_state$State,4), 
                       Animal = rep(animal_state$Animal,4),
                       sick = c("Yes","Yes","No","No"),
                       unprovoked = c("Yes","No","Yes","No"),
                       vaxed = rep(0,4)) # NA
  } else{
    df_new <- data.frame(prob = rep(animal_state$prob,8), 
                         size = rep(animal_state$size,8),
                         State = rep(animal_state$State,8), 
                         Animal = rep(animal_state$Animal,8),
                         sick = rep(c("Yes","Yes","No","No"),2),
                         unprovoked = rep(c("Yes","No","Yes","No"),2),
                         vaxed = c(rep(1,4), rep(0,4))) # Yes; no/unknown
    
  }
  # Add df to list
  my_list[[i]] <- df_new
  
}

# Combine dfs
new_df <- do.call("rbind", my_list)

###############################
# Use functions to calculate point est, LL, and UL of 95% CI
##############################################################

# Split data by ID for state and animal
ids_df <- select(df, Animal, State, id)

new_df <- left_join(new_df, ids_df, by = c("Animal","State"))

data_split <- split(new_df, ~ id)
x <- data_split[[311]]

# Use a function to calculate probabilities (point est and 95% CI)
source("Analysis/Functions/New functions to create df with vax.R")

n_sims <- 1000000

# RABV vaccine efficacy in some animals
vacc_eff <- 1-0.95

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

# Takes 7 min?
start_time <- Sys.time()

res <- lapply(data_split, function(x) calc_prob_new(size = x$size, 
                                                    prob = x$prob,
                                                    vaxed = x$vaxed,
                                                    n_sims = n_sims,
                                                    vacc_eff = vacc_eff,
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

# Add in other variables (there is surely a better way to do this)
final_df <- cbind(final_df, 
                  State = new_df$State,
                  Animal = new_df$Animal,
                  sick = new_df$sick,
                  unprovoked = new_df$unprovoked,
                  vaxed = new_df$vaxed)

final_df <- data.frame(final_df)

# Re-order rows
final_df <- dplyr::select(final_df, Animal, State,
                          sick, unprovoked, vaxed, lower_95CI, point_est, upper_95CI)

df_final_tool2 <- final_df

######## Initially, skip to line 307, run through rest of the code and save an RDS file (the "not expanded one")
######## then come back to line 187 and run through the code, saving this dataset as final_df_for_tool_exposure_vaxed_Bayes_2011_2020_18Nov2022.RDS

########################################################
# Expand out "All" to include all combinations of states
# no terrestrial rabies
all_terr1 <- filter(final_df, 
                     Animal == "Cattle"|
                     Animal == "Rodents of unusual size"|
                     Animal == "Sheep and goats"|
                     Animal == "Equine"|
                     Animal == "Hoofed animals")

# terrestrial rabies
all_terr2 <- filter(final_df, 
                      Animal == "Cattle terr"|
                      Animal == "Rodents of unusual size terr"|
                      Animal == "Sheep and goats terr"|
                      Animal == "Equine terr"|
                      Animal == "Hoofed animals terr")

all_not_terr <- filter(final_df, State == "All" &
                         Animal != "Cattle"&
                         Animal != "Cattle terr"&
                         Animal != "Rodents of unusual size"&
                         Animal != "Rodents of unusual size terr"&
                         Animal != "Sheep and goats"&
                         Animal != "Sheep and goats terr"&
                         Animal != "Equine"&
                         Animal != "Equine terr"&
                         Animal != "Hoofed animals"&
                         Animal != "Hoofed animals terr" )

check <- filter(final_df, State == "All")
unique(check$Animal)
unique(all_not_terr$Animal)
unique(all_terr$Animal)

the_rest <- filter(final_df, State != "All")

# Check
table(the_rest$State)

states <- sort(unique(the_rest$State))

# National estimates not including those that will be separated based on terrestrial rabies status
all_not_terr$id <- 1:nrow(all_not_terr)


my_dfs <- list()
for (i in 1:nrow(all_not_terr)){
  # Filter by row
  row_of_interest <- all_not_terr[i,]
  # Add rows for all states
  df_all_states <- purrr::map_dfr(seq_len(length(states)), ~row_of_interest)
  # Add state names
  df_all_states$State <- states
  # Add df to list
  my_dfs[[i]] <- df_all_states
}

# Combine dfs
df_new_all <- do.call("rbind", my_dfs)
df_new_all <- select(df_new_all, -id)

##### Now take care of the terrestrial rabies estimates

# States without terrestrial rabies
states_no_terr <- c("WA","OR","NV","ID","IL","IN","MS","UT")
# states with terrestrial rabies
states_terr <- setdiff(states, states_no_terr)

all_terr1$id <- 1:nrow(all_terr1)
all_terr2$id <- 1:nrow(all_terr2)

# no terr rabies
my_dfs <- list()
for (i in 1:nrow(all_terr1)){
  # Filter by row
  row_of_interest <- all_terr1[i,]
  # Add rows for states
  df_all_states <- purrr::map_dfr(seq_len(length(states_no_terr)), ~row_of_interest)
  # Add state names
  df_all_states$State <- states_no_terr
  # Add df to list
  my_dfs[[i]] <- df_all_states
}

# Combine dfs
df_new_all2 <- do.call("rbind", my_dfs)
df_new_all2 <- select(df_new_all2, -id)

# terr rabies
my_dfs <- list()
for (i in 1:nrow(all_terr2)){
  # Filter by row
  row_of_interest <- all_terr2[i,]
  # Add rows for all states
  df_all_states <- purrr::map_dfr(seq_len(length(states_terr)), ~row_of_interest)
  # Add state names
  df_all_states$State <- states_terr
  # Add df to list
  my_dfs[[i]] <- df_all_states
}

# Combine dfs
df_new_all3 <- do.call("rbind", my_dfs)
df_new_all3 <- select(df_new_all3, -id)

df_new_all3$Animal[df_new_all3$Animal == "Cattle terr"] <- "Cattle"
df_new_all3$Animal[df_new_all3$Animal == "Equine terr"] <- "Equine"
df_new_all3$Animal[df_new_all3$Animal == "Sheep and goats terr"] <- "Sheep and goats"
df_new_all3$Animal[df_new_all3$Animal == "Rodents of unusual size terr"] <- "Rodents of unusual size"
df_new_all3$Animal[df_new_all3$Animal == "Hoofed animals terr"] <- "Hoofed animals"

table(df_new_all3$Animal)

# Combine the results for the expanded states with the 
# state-specific estimates
df_final_tool <- rbind(the_rest, df_new_all, df_new_all2, df_new_all3)

df_final_tool2 <- df_final_tool

################################################################

# Formatting
df_final_tool2$lower_95CI <- as.numeric(df_final_tool2$lower_95CI)
df_final_tool2$point_est <- as.numeric(df_final_tool2$point_est)
df_final_tool2$upper_95CI <- as.numeric(df_final_tool2$upper_95CI)

# Leave Nans out for now. Come back to them?

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
df_final_tool2$Animal <- as.factor(df_final_tool2$Animal)
df_final_tool2$State <- as.factor(df_final_tool2$State)
df_final_tool2$Unprovoked <- as.factor(df_final_tool2$Unprovoked)
df_final_tool2$Sick <- as.factor(df_final_tool2$Sick)

# Change names
df_final_tool2 <- df_final_tool2 %>%
  rename(Healthy = Sick, Provoked = Unprovoked)

df_final_tool2 <- df_final_tool2 %>%
  rename(LL95CI = LL_95CI, UL95CI = UL_95CI, Vaxed = vaxed)

# Reorder variables
df_final_tool2 <- df_final_tool2 %>%
  select(Animal, State, Provoked, Healthy, Vaxed, LL95CI, Median, UL95CI)

# Fixed vaxed
df_final_tool2$Vaxed[df_final_tool2$Vaxed == 1] <- "Yes"
df_final_tool2$Vaxed[df_final_tool2$Vaxed == 0 & 
                       df_final_tool2$Animal %in% licensed_vacc_cats$Animal] <- "No/unknown"
df_final_tool2$Vaxed[df_final_tool2$Vaxed == 0 &
                       !(df_final_tool2$Animal %in% licensed_vacc_cats$Animal)] <- "NA"

df_final_tool2 <- dplyr::rename(df_final_tool2, Vaccinated = Vaxed)

# Checking the df
dog_FL <- filter(df_final_tool2, Animal == "Dog" & State == "FL")

cat_CA <- filter(df_final_tool2, Animal == "Cat" & State == "CA")

cattle_OH <- filter(df_final_tool2, Animal == "Cattle" & State == "OH")

# Save final df for tool
#saveRDS(df_final_tool2, "final_df_for_tool_exposure_vaxed_Bayes_2011_2020_18Nov2022.RDS")

# save a version without the estimates expanded for "all" states
#saveRDS(df_final_tool2, "df_for_tool_exposure_not_expanded_vaxed_Bayes_2011_2020_18Nov2022.RDS")
