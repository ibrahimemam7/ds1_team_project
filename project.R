# load relevant packages
library(funModeling)
library(tidyverse)
library(MASS)
library(mice)
library(car)
library(gtsummary)

#######################
## Declare Functions ##
#######################

#' having functions at the top of the file is important because they will be
#' available to use in the environment if the file is sourced (runs from top
#' to bottom)

#' function to determine the max number of degrees of freedom in the predictor
#' for a given response variable
max.predictors <- function(response_var, dataframe) {
  if(is.numeric(dataframe[[response_var]])) {
    
    # the variable is numeric
    tmp <- dataframe %>%
      summarise(n_observations = n())
    
    return(floor(tmp$n_observations/15))
      
  } else {
    
    # the variable is logical
    tmp <- dataframe %>% 
      summarise(
        true_count = sum(get(response_var)),
        false_count = sum(!get(response_var))
      )
    
    return(floor(min(tmp$true_count, tmp$false_count) / 15))
  
  }
}

###############
## Data Prep ##
###############

# import data into data frame
d <- read.csv("ds_project_data.csv")

# select relevant columns, keep all rows
d <- d[, c("Pittsburgh.Sleep.Quality.Index.Score",
           "Epworth.Sleepiness.Scale",
           "Berlin.Sleepiness.Scale",
           "Athens.Insomnia.Scale",
           "SF36.PCS",
           "SF36.MCS",
           "Age",
           "Gender",
           "BMI",
           "Time.from.transplant",
           "Liver.Diagnosis",
           "Recurrence.of.disease",
           "Rejection.graft.dysfunction",
           "Any.fibrosis",
           "Renal.Failure",
           "Depression",
           "Corticoid")]

# rename columns with long names for ease of coding
colnames(d)[1:4] <- c("PSQI", "ESS", "BSS", "AIS")

# check for duplicates, there are none
any(duplicated(d))

# check for any variables with high rate of missingness
status(d)

#' PSQI has ~32% missing data! This is above our chosen threshold of 30%
#' We will exclude PSQI from analysis.

# formatting columns as appropriate
d$Gender <- as.factor(d$Gender)
levels(d$Gender)

d$Liver.Diagnosis <- as.factor(d$Liver.Diagnosis)
levels(d$Liver.Diagnosis)

logical_columns <- c("BSS","Recurrence.of.disease", "Rejection.graft.dysfunction", "Any.fibrosis", "Renal.Failure", "Corticoid", "Depression")
d[, logical_columns] <- as.logical(c(d$BSS, d$Recurrence.of.disease, d$Rejection.graft.dysfunction, d$Any.fibrosis, d$Renal.Failure, d$Corticoid, d$Depression))

# remove invalid data and replace with NA
# ESS is score on a scale from 0-24, cannot have values higher than 24
sum(d$ESS > 24, na.rm = T)
d$ESS[d$ESS > 24] <- NA

# data for other scales is within expected ranges

###############################
## Exploratory Data Analysis ##
###############################

#-ESS, -AIS, -SF36.MCS, -SF36.PCS, -BSS

# create table with summary of baseline characteristics
d %>% 
  select(-PSQI, -ESS, -AIS, -SF36.MCS, -SF36.PCS, -BSS) %>% 
  tbl_summary(
    missing = "no",
    statistic = list(
      all_continuous() ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)")) %>%
  add_n() %>% 
  bold_labels()

# create histogram for ESS score
ggplot(data = d, mapping = aes(x = ESS)) +
  geom_histogram(binwidth = 1, width = 0.8, fill = "lightgrey", col = "black") +
  labs(
    x = "ESS Score",
    y = "Frequency",
    title = "Distribution of ESS Scores") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# create histogram for AIS score
ggplot(data = d, mapping = aes(x = AIS)) +
  geom_histogram(binwidth = 1, width = 0.8, fill = "lightgrey", col = "black") +
  labs(
    x = "AIS Score",
    y = "Frequency",
    title = "Distribution of AIS Scores") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# create histogram for SF36.MCS scores
ggplot(data = d, mapping = aes(x = SF36.MCS)) +
  geom_histogram(binwidth = 2, width = 0.8, fill = "lightgrey", col = "black") +
  labs(
    x = "SF36-MCS Score",
    y = "Frequency",
    title = "Distribution of SF36-MCS Scores") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

# create histogram for SF36.PCS scores
ggplot(data = d, mapping = aes(x = SF36.PCS)) +
  geom_histogram(binwidth = 2, width = 0.8, fill = "lightgrey", col = "black") +
  labs(
    x = "SF36-PCS Score",
    y = "Frequency",
    title = "Distribution of SF36-PCS Scores") +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5))

#################################
## Imputation for Missing Data ##
#################################

# set appropriate imputation method for each variable in the dataset 

#' "norm.nob" used for numerical variables, "logreg" used for binary variables.
#' Polyreg and polr methods are only mentioned because the vector must have a
#' length of 4. Additionally, no method will be used for variables
#' that have no missing values (no imputation necessary).
methods <- c("norm.nob", "logreg", "polyreg", "polr")

# create an imputed dataset with the methods specified above
d_imputed <- mice(d, defaultMethod = methods, print = F, seed = 7, m=1)

# check that correct method was used for each column
d_imputed$method

# extract the imputed dataset
d_complete <- complete(d_imputed, action = 1)

# check that no NAs remain
status(d_complete)

# correct imputed values that are outside of possible range for clinical scores
summary(d_complete)
d_complete$ESS[d_complete$ESS < 0] <- 0
d_complete$AIS[d_complete$AIS < 0] <- 0

#' it seems mice() converted BSS column from logical to numerical, so it should be
#' converted back to logical
d_complete$BSS <- as.logical(d_complete$BSS)

##############################
## Sleep Disturbance Models ##
##############################

############ ESS Model ################

# determine the max number of predictors for a model with ESS as a response variable
max.predictors("ESS", d_complete)

# create a liner model based on clinical literature
ess.mod <- lm(ESS ~ Gender + Time.from.transplant + BMI + Depression + Rejection.graft.dysfunction, data = d_complete)
summary(ess.mod)

# check if a simpler model is better
ess.mod.simple <- stepAIC(ess.mod, trace = F)
summary(ess.mod.simple)

# Check the AIC of the complex and simpler model
AIC(ess.mod, ess.mod.simple)

# compare the two models using deviance
deviance(ess.mod)
deviance(ess.mod.simple)

# compare the models via ANOVA since these models are nested
anova(ess.mod, ess.mod.simple)

#' AIC and anova indicate the simpler model is better, while deviance indicates the
#' complex model is better

# check for co-linearity
vif(ess.mod.simple)

# no values above 5, so there is no concern for co-linearity

################ AIS Model ##################

# check max df in predictors for model with AIS as response variable
max.predictors("AIS", d_complete)

# start with the predictors from literature
ais.mod <- lm(AIS ~ Gender + Time.from.transplant + BMI + Depression + Rejection.graft.dysfunction, data = d_complete)
summary(ais.mod)

# use AIC step back to see if the model can be simplified
ais.mod.simple <- stepAIC(ais.mod, trace = F)
summary(ais.mod.simple)

# compare the AIC of the two models
AIC(ais.mod, ais.mod.simple)

# compare the two models using deviance
deviance(ais.mod)
deviance(ais.mod.simple)

# compare the models via ANOVA since these models are nested
anova(ais.mod, ais.mod.simple)

#' the simpler model is better according to both ANOVA and AIC tests. Deviance
#' indicates that the more complex model is better.

# check for co-linearity
vif(ais.mod.simple)

# no values above 5, so there is no concern for co-linearity

################ BSS Model ##################

# check max df in predictors for model with BSS as response variable
max.predictors("BSS", d_complete)

# start with the predictors from literature
bss.mod <- glm(BSS ~ Gender + Time.from.transplant + BMI + Depression + Rejection.graft.dysfunction, data = d_complete, family = "binomial")
summary(bss.mod)

# use AIC step back to see if the model can be simplified
bss.mod.simple <- stepAIC(bss.mod, trace = F)
summary(bss.mod.simple)

# compare the AIC of the two models
AIC(bss.mod, bss.mod.simple)

# compare the deviance
deviance(bss.mod)
deviance(bss.mod.simple)

# use anova LRT to compare the two models
anova(bss.mod, bss.mod.simple, test = "Chisq")

# AIC and anova indicate that the simpler model is better, but deviance indicates the opposite

# check for co-linearity
vif(bss.mod.simple)

# no values above 5, so there is no concern for co-linearity