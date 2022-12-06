# Run function to create df for conditional probability of death from rabies exposure
# 25-2-22
# Uses new method for adjusting positivity rates by healthy/provoked
# Updated Aug 8, 2022 - we decided to NOT adjust single bite to head/neck
# updated 30 Sept 2022 to stratify pos rates for certain animals based on whether state
# has terrestrial rabies
# Updated November 18, 2022 with corrected est. for provoked params

library(dplyr)
library(triangle)

set.seed(10) # Set random seed for reproducible analysis

# Number of times to sample from distributions
n_sims <- 1000000

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

# Check
#hist(p2_list[[5]])

#########################################
# Create a df with all the combinations
#########################################

# Should be 328 * 28* 2 * 2 + + (4* nrow(check)*28) = 47488
df$id <- 1:nrow(df)
my_list <- list()
for (i in 1:nrow(df)){
  animal_state <- filter(df, id == i)
  
  if (animal_state$vacc_cat == 0){
  new_df <- data.frame(prob = rep(animal_state$prob,112), 
                       size = rep(animal_state$size,112),
                       State = rep(animal_state$State,112), 
                       Animal = rep(animal_state$Animal,112),
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
                                      rep("No",28)),
                       vaxed = rep(0,112)) # NA)
  } else{
    
    new_df <- data.frame(prob = rep(animal_state$prob,112*2), 
                         size = rep(animal_state$size,112*2),
                         State = rep(animal_state$State,112*2), 
                         Animal = rep(animal_state$Animal,112*2),
                         location_cat = rep(loc_severity$location_cat, 4*2),
                         severity = rep(loc_severity$severity,4*2),
                         sick = rep(c(rep(c("Yes"), 28*2), 
                                  rep(c("No"), 28*2)),
                                  2),
                         unprovoked = rep(c(rep("Yes",28),
                                        rep("No",28),
                                        rep("Yes",28),
                                        rep("No",28)),
                                        2),
                         vaxed = c(rep(1,112), rep(0,112))) # Yes; no/unknown
  }
  # Add df to list
  my_list[[i]] <- new_df
  
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
source("Analysis/Functions/New functions to calc death vaxed.R")

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

# RABV vaccine efficacy in some animals
vacc_eff <- 1-0.95

######################

# Takes 34 min.
start_time <- Sys.time()

res <- lapply(data_split, function(x) calc_prob_death(size = x$size, 
                                                      prob = x$prob,
                                                      vaxed = x$vaxed,
                                                      n_sims = n_sims,
                                                      vacc_eff = vacc_eff,
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
beepr::beep()

# Combine dfs
final_df <- do.call("rbind", res)

# Add in other variables using the df structure created earlier in the script (surely there is a better way to do this)
final_df <- cbind(final_df, 
                  State = new_df$State,
                  Animal = new_df$Animal,
                  location_cat = new_df$location_cat,
                  severity = new_df$severity,
                  sick = new_df$sick,
                  unprovoked = new_df$unprovoked,
                  vaxed = new_df$vaxed)

final_df <- data.frame(final_df)

# Re-order rows
final_df <- dplyr::select(final_df, Animal, State, location_cat, severity, 
                          sick, unprovoked, vaxed, lower_95CI, point_est, upper_95CI)

df_final_tool2 <- final_df

######## Initially, skip to line 379, run through rest of the code and save an RDS file (the "not expanded one")
######## then come back to line 263 and run through the code, saving this dataset as final_df_for_tool_death_vaxed_Bayes_2011_2020_18Nov2022.RDS

###############################################################
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

all <- filter(final_df, State == "All" &
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

the_rest <- filter(final_df, State != "All")

# Check
table(the_rest$State)

states <- sort(unique(the_rest$State))

all$id <- 1:nrow(all)

# Should have 1680*51 = 85680 rows in new "All" dataset
my_dfs <- list()
for (i in 1:nrow(all)){
  # Filter by row
  row_of_interest <- all[i,]
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

###########################################################################

# Formatting
df_final_tool2$lower_95CI <- as.numeric(df_final_tool2$lower_95CI)
df_final_tool2$point_est <- as.numeric(df_final_tool2$point_est)
df_final_tool2$upper_95CI <- as.numeric(df_final_tool2$upper_95CI)

# Leave Nans out for now. Come back to them?

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
df_final_tool2$Animal <- as.factor(df_final_tool2$Animal)
df_final_tool2$State <- as.factor(df_final_tool2$State)
df_final_tool2$Unprovoked <- as.factor(df_final_tool2$Unprovoked)
df_final_tool2$Sick <- as.factor(df_final_tool2$Sick)

# Change names
df_final_tool2 <- df_final_tool2 %>%
  rename(Healthy = Sick, Provoked = Unprovoked, Vaxed = vaxed, 
         LL95CI = Lower_95CI, UL95CI = Upper_95CI, Median = Point_est)

# Reorder and fix names
df_final_tool2 <- df_final_tool2 %>%
  rename(Severity = severity, Location_cat = location_cat) %>%
  select(Location_cat, Severity, Animal, State, Provoked, Healthy, Vaxed, LL95CI, 
         Median, UL95CI)

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

bats_MI <- filter(df_final_tool2, Animal == "Bat" & State == "MI")

# Save final df for tool
#saveRDS(df_final_tool2, "final_df_for_tool_death_vaxed_Bayes_2011_2020_18Nov2022.RDS")

# save a version without the estimates expanded for "all" states
#saveRDS(df_final_tool2, "df_for_tool_death_not_expanded_vaxed_Bayes_2011_2020_18Nov2022.RDS")

