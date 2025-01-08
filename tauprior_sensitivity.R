# tau prior sensitivty analysis
#----
rm(list=ls())
dev.off()
if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, car, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, ggstatsplot, gtable, gridExtra, knitr, MAd, metafor, outliers, reactable, readr, readxl, 
               R.rsp, sjPlot, stringr, tidyr) 
source("HelperFunctions.R")
load(file = "df.RDa")
df <- df %>% filter(Stimulation.Type %in% input$stimtype,
                    Sample %in% input$Sample,
                    Blindness %in% input$blind,
                    Mean.Age >= input$age[1], Mean.Age <= input$age[2],
                    Percent.Females >= input$gender[1], Percent.Females <= input$gender[2],
                    Stimulation.side %in% input$stimside,
                    Part.of.the.ear.stimulated %in% input$partear,
                    Publication.Year >= input$pubyear[1], Publication.Year <= input$pubyear[2],
                    Control.Type %in% input$control,
                    Pupil.Dilation.Type %in% input$dilationtype,
                    Stimulation.duration.sec. >= input$stimduration[1], Stimulation.duration.sec. <= input$stimduration[2])

                    
aggES <- agg(id     = df$Paper.No,
             es     = yi,
             var    = vi,
             data   = df,
             cor = .5,
             method = "BHHR")
MA <- merge(x = aggES, y = df, by.x = "id", by.y = "Paper.No")
MA <- unique(setDT(MA) [sort.list(id)], by = "id")
MA <- with(MA, MA[order(MA$es)])

######################################

# adding my extra line to filetr for complete rows only here as well
MA <- MA[complete.cases(MA$es), ]

######################################

#----
# Our default analysis:
bma_default <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                        tau.prior = function(t) dhalfcauchy(t, scale = 0.5), 
                        mu.prior = c("mean" = 0, "sd" = 1.5))
summary(bma_default)

# Analysis with scale = 1
bma_scale1 <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                 tau.prior = function(t) dhalfcauchy(t, scale = 1), 
                 mu.prior = c("mean" = 0.64, "sd" = 1.5)) # here I used mean from a meta analysis on functional dyspepsia. If I also adapt sd I will get a flat line for the prior
summary(bma_scale1)
plot(bma_scale1, which=3, prior = TRUE)

# Analysis with scale = 2
bma_scale2 <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                        tau.prior = function(t) dhalfcauchy(t, scale = 2), 
                        mu.prior = c("mean" = 0, "sd" = 1.5))
summary(bma_scale2)
plot(bma_scale2, which=3, prior = TRUE)

# Analysis with scale = 0.1
bma_scale0.1 <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                        tau.prior = function(t) dhalfcauchy(t, scale = 0.1), 
                        mu.prior = c("mean" = 0, "sd" = 1.5))
summary(bma_scale0.1)
plot(bma_scale0.1, which=3, prior = TRUE)

