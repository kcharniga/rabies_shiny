# Determine the threshold "p" in the logistic regression
# using final NASPHV survey results
# 10-fold cross validation and optimizing on F score

# Refs: https://www.theissaclee.com/post/logistic-regression-beta/
# https://towardsdatascience.com/accuracy-recall-precision-f-score-specificity-which-to-optimize-on-867d3f11124
# https://medium.com/swlh/cross-validation-techniques-to-assess-your-models-stability-3a4d55d90409#:~:text=For%20such%20data%2C%20a%20different%20cross-validation%20technique%20is,the%20same%20number%20of%20samples%20of%20each%20class.
# https://www.hpl.hp.com/techreports/2009/HPL-2009-359.pdf

library(ggplot2)
library(dplyr)
library(magrittr)
library(jtools)
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

# Data split into 5 (20% each)
dat$id <- 1:nrow(dat)

#Randomly shuffle the data
set.seed(10)
dat <- dat[sample(nrow(dat)),]

#Create 10 equally size folds
folds <- cut(seq(1,nrow(dat)),breaks=10,labels=FALSE)

my_list <- list()
#Perform 10 fold cross validation
for(i in 1:10){
  #Segement your data by fold using the which() function 
  testIndexes <- which(folds==i,arr.ind=TRUE)
  test <- dat[testIndexes, ]
  train <- dat[-testIndexes, ]
  
  # Estimate coefficients
  model <- glm(PEP_yes ~ log_median, 
               data = train, 
               family = binomial(link = "logit"))
  
  
  # Predictions using test data
  test2 <- select(test, PEP_yes, log_median)
  
  predicted_prob <- predict(model, test2, type = "response")
  
  # Make a function for confusion matrix and calculation of TP, TN, and accuracy for different p
  find_threshold <- function(p){
    
    # Make confusion matrix using p as threshold
    prediction <- as.integer(predicted_prob > p)
    confusion_mat <- addmargins(table(test$PEP_yes, prediction))
    
    # Labeling
    names(dimnames(confusion_mat)) <- c("True status", "Prediction")
    colnames(confusion_mat) <- c("No", "Yes", "Total")
    rownames(confusion_mat) <- c("No", "Yes", "Total")
    confusion_mat
    
    # False positives 
    FP <- confusion_mat[1,2]
    # True positives 
    TP <- confusion_mat[2,2]
    # False negatives
    FN <- (confusion_mat[2,1]) 
    
    return(list(FP = FP, TP = TP, FN = FN))
    
  }
  
  p <- seq(from = 0.35, to = 0.75, by = 0.05)
  
  res <- sapply(p, FUN = find_threshold)
  
  df <- as.data.frame(res)
  
  colnames(df) <- paste("p", p, sep = "_")
  
  # Add a variable for fold
  df$fold <- paste("fold_",i, sep = "")
  
  # Add results to list
  my_list[[i]] <- df
  
}

res <- do.call("rbind", my_list)

# Sum all TP, FP, and FN
FP <- res[c(1,4,7,10,13,16,19,22,25,28),]
TP <- res[c(2,5,8,11,14,17,20,23,26,29),]
FN <- res[c(3,6,9,12,15,18,21,24,27,30),]

# EQUATION  (2 * TP) / (2 * TP + FP + FN)
 # F scores for each decision boundary
(2*sum(unlist(TP$p_0.35))) / (2 * sum(unlist(TP$p_0.35)) + sum(unlist(FP$p_0.35)) + sum(unlist(FN$p_0.35)))
(2*sum(unlist(TP$p_0.4))) / (2 * sum(unlist(TP$p_0.4)) + sum(unlist(FP$p_0.4)) + sum(unlist(FN$p_0.4)))
(2*sum(unlist(TP$p_0.45))) / (2 * sum(unlist(TP$p_0.45)) + sum(unlist(FP$p_0.45)) + sum(unlist(FN$p_0.45)))
(2*sum(unlist(TP$p_0.5))) / (2 * sum(unlist(TP$p_0.5)) + sum(unlist(FP$p_0.5)) + sum(unlist(FN$p_0.5)))
(2*sum(unlist(TP$p_0.55))) / (2 * sum(unlist(TP$p_0.55)) + sum(unlist(FP$p_0.55)) + sum(unlist(FN$p_0.55)))
(2*sum(unlist(TP$p_0.6))) / (2 * sum(unlist(TP$p_0.6)) + sum(unlist(FP$p_0.6)) + sum(unlist(FN$p_0.6)))
(2*sum(unlist(TP$p_0.65))) / (2 * sum(unlist(TP$p_0.65)) + sum(unlist(FP$p_0.65)) + sum(unlist(FN$p_0.65)))
(2*sum(unlist(TP$p_0.7))) / (2 * sum(unlist(TP$p_0.7)) + sum(unlist(FP$p_0.7)) + sum(unlist(FN$p_0.7)))
(2*sum(unlist(TP$p_0.75))) / (2 * sum(unlist(TP$p_0.75)) + sum(unlist(FP$p_0.75)) + sum(unlist(FN$p_0.75)))


# Find the x values (risk thresholds) corresponding to these decision boundaries 
model <- glm(PEP_yes ~ log_median, 
             data = dat, 
             family = binomial(link = "logit"))

find_x_val <- function(p){
  x_int <- (log(p/(1-p)) - coef(model)[1]) / coef(model)[2]
  return(x_int)
}

x_vals <- sapply(p, FUN = find_x_val)

# Plot the thresholds on the all data
my_cols <- rainbow(9)
my_cols2 <- paletteer::paletteer_c("grDevices::rainbow", 9)

library(scales)
scaleFUN <- function(x) {sprintf("%.1g", exp(x))}

ggplot(dat, aes(x = log_median, y = PEP_yes)) + 
  #geom_point(alpha=.5) +
  geom_jitter(size = 1, width = 0, height = 0.1) +
  stat_smooth(method = "glm", se = F, method.args = list(family=binomial), 
              color = "black") +
  theme_bw() +
  geom_vline(xintercept = x_vals[1], linetype = "dashed", color = my_cols2[1], size = 1.05) +
  geom_vline(xintercept = x_vals[2], linetype = "dashed", color = my_cols2[2], size = 1.05) +
  geom_vline(xintercept = x_vals[3], linetype = "dashed", color = my_cols2[3], size = 1.05) +
  geom_vline(xintercept = x_vals[4], linetype = "dashed", color = my_cols2[4], size = 1.05) +
  geom_vline(xintercept = x_vals[5], linetype = "dashed", color = my_cols2[5], size = 1.05) +
  geom_vline(xintercept = x_vals[6], linetype = "dashed", color = my_cols2[6], size = 1.05) +
  geom_vline(xintercept = x_vals[7], linetype = "dashed", color = my_cols2[7], size = 1.05) +
  geom_vline(xintercept = x_vals[8], linetype = "dashed", color = my_cols2[8], size = 1.05) +
  geom_vline(xintercept = x_vals[9], linetype = "dashed", color = my_cols2[9], size = 1.05) +
  #geom_vline(xintercept = x_vals[10], linetype = "dashed", color = my_cols2[10], size = 1.05) +
  #geom_vline(xintercept = x_vals[11], linetype = "dashed", color = my_cols2[11], size = 1.05) +
  #geom_vline(xintercept = x_vals[12], linetype = "dashed", color = my_cols2[12], size = 1.05) +
  #geom_vline(xintercept = x_vals[13], linetype = "dashed", color = my_cols2[13], size = 1.05) +
  scale_x_continuous(breaks=c(-11.5,-6.9,-4.5,-2,0), 
                     labels = scaleFUN) +
  scale_y_continuous(labels = c(0, 0.25, 0.5, 0.75, 1),
                     breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  xlab("Pr(rabid|exposure)") +
  ylab("PEP yes") 

