# load relevant packages
library(funModeling)
library(tidyverse)
library(MASS)
library(mice)

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

##########################################
## Multiple Imputation for Missing Data ##
##########################################

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

# determine the max number of predictors for a model with ESS as a response variable
max.predictors("BSS", d_complete)

# experiment with stepwise models
ess.full.mod <- lm(ESS~., data = d_complete)
ess.null.mod <- lm(ESS~1, data = d_complete)

ess.step.back <- stepAIC(ess.full.mod, trace = F)
summary(ess.step.back)

ess.step.forw <- stepAIC(ess.null.mod, direction = "forward", trace = F, scope = list(upper=ess.full.mod, lower=ess.null.mod))
summary(ess.step.forw)

# create a liner model based on clinical literature
ess.mod <- lm(ESS ~ Gender + Time.from.transplant + BMI + Depression, data = d_complete)
summary(ess.mod)

# create a simpler model without time from transplant and BMI, then compare
ess.mod1 <- lm(ESS ~ Gender + Depression, data = d_complete)
summary(ess.mod1)

AIC(ess.mod, ess.mod1)

# the simpler model (model 1) is better

# can try adding a predictor from the stepwise function and see if the model is improved
ess.mod2 <- lm(ESS ~ Gender + Depression + Rejection.graft.dysfunction, data = d_complete)
summary(ess.mod2)

AIC(ess.mod1, ess.mod2)

# the more complex model (model 2) is better

#' can try adding another predictor since there are still two degrees of freedom
#' in the predictors that can be used
ess.mod3 <- lm(ESS ~ Gender + Depression + Rejection.graft.dysfunction + Renal.Failure, data = d_complete)
summary(ess.mod3)

AIC(ess.mod2, ess.mod3)

# the simpler model (model 2) is better than model 3. Can try adding a different predictor
ess.mod4 <- lm(ESS ~ Gender + Depression + Rejection.graft.dysfunction + Recurrence.of.disease, data = d_complete)
summary(ess.mod4)

AIC(ess.mod2, ess.mod4)

# the simpler model (model 2) is better than model 4. Can try adding a different predictor
ess.mod5 <- lm(ESS ~ Gender + Depression + Rejection.graft.dysfunction + Liver.Diagnosis, data = d_complete)
summary(ess.mod5)

AIC(ess.mod2, ess.mod5)

# the simpler model (model 2) is better than model 5

# can try a more complex model
ess.mod6 <- lm(ESS ~ Gender + Depression + Rejection.graft.dysfunction + Recurrence.of.disease + Liver.Diagnosis, data = d_complete)
summary(ess.mod6)

AIC(ess.mod2, ess.mod6) 

# simpler model (model 2) is still better
