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

#' new columns to represent ESS and AIS as binary variable
#' an ESS score above 10 is indicative of excessive daytime sleepiness
#' an AIS score above 10 is considered clinical insomnia
d$ESS_binary <- d$ESS > 10
d$AIS_binary <- d$AIS > 10

##########################################
## Multiple Imputation for Missing Data ##
##########################################

# set appropriate imputation method for each variable in the dataset 

#' "norm.nob" used for numerical variables, "logreg" used for binary variables, and 
#' no method ("") used for PSQI since it will be excluded from analysis due to
#' excessive missing data. Additionally, no method will be used for variables
#' that have no missing values (no imputation necessary).
imp_methods <- c(
  "PSQI" = "",
  "ESS" = "norm.nob",
  "BSS" = "logreg",
  "AIS" = "norm.nob",
  "SF36.PCS" = "norm.nob",
  "SF36.MCS" = "norm.nob",
  "Age" = "norm.nob",
  "Gender" = "",
  "BMI" = "norm.nob",
  "Time.from.transplant" = "",
  "Liver.Diagnosis" = "",
  "Recurrence.of.disease" = "",
  "Rejection.graft.dysfunction" = "",
  "Any.fibrosis" = "",
  "Renal.Failure" = "",
  "Depression" = "",
  "Corticoid" = "",
  "ESS_binary" = "logreg",
  "AIS_binary" = "logreg"
)

# create 5 imputed data sets with the methods specified above
d_imputed <- mice(d, method = imp_methods, seed = 7, m = 5, print = FALSE)

# extract the first of 5 imputed data sets
imp1 <- complete(d_imputed, action = 1) 
status(imp1)
