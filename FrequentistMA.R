#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYESIAN META-ANALYSIS #################
#######################################################################################

################### Frequentist Meta Analysis #########################################
#----
# Clear global environment and load required packages
rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr) 

setwd("Y:/Paper/Drafts/pulsedtaVNS_pupil_meta/code")

# Source helper functions
source("HelperFunctions.R")

# Load "df.RDa"
load(file = "df.RDa")

# Filter sham and healthy samples
#df <- df %>% filter(Sample %in% c("healthy","high worriers"),
#                    Blindness %in% c("double blind","single blind"))

# Aggregate within study effect sizes
aggES <- agg(id     = Paper.No,
             es     = yi,
             var    = vi,
             data   = df,
             cor = .5,
             method = "BHHR")

# Merging aggregated ES with original dataframe 
MA <- merge(x = aggES, y = df, by.x = "id", by.y = "Paper.No")
MA <- unique(setDT(MA) [sort.list(id)], by = "id")
MA <- with(MA, MA[order(MA$es)])

##############################################

# I try to add this here as well (it takes out all the NaN values from es)
MA <- MA[complete.cases(MA$es), ]

#############################################


# Meta regression
m0<-myMareg(data = MA)[1:2]

data <- MA

# model <- mareg(formula = es~1,
#               var = var,
#               data = data,
#              slab = data$study)
model <- rma(yi, vi, data=data, method="FE", weighted=FALSE)

weights(model)

# Model summary
myCoef(m0$Modelsummary$coef)

# Model fit
myFit(m0$Modelsummary$fit)

# Model uniqueness
myUni(m0$`Model-Fit`$random)

