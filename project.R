# load relevant packages
library(funModeling)
library(tidyverse)
library(MASS)
library(mice)


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
