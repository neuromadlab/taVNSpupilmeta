### Main Effect Meta Analysis ####

# Clear global environment and load required packages
rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr) 

#setwd("D:/SynologyDrive/BON_general/SynologyDrive/Paper/Drafts/pulsedtaVNS_pupil_meta/code")
#setwd("D:/SynologyDrive/Paper/Drafts/pulsedtaVNS_pupil_meta/code")

# Source helper functions
source("HelperFunctions.R")

# Load "df.RDa"
load(file = "df.RDa")


df %>% 
  mutate(study.name = fct_cross(study,factor(Study.No),sep = ": Study ")) %>% 
  filter(Control.Type %in% c("earlobe sham","scapha sham"))-> df


df1 <- df %>% filter(Stimulation.Type %in% c("conventional")) 
df1 %>% 
  mutate(study.name = fct_cross(study,factor(Study.No),sep = ": Study ")) %>% 
  filter(Control.Type %in% c("earlobe sham","scapha sham"))-> df1

df2 <- df %>% filter(Stimulation.Type %in% c("pulsed"))
df2 %>% 
  mutate(study.name = fct_cross(study,factor(Study.No),sep = ": Study ")) %>% 
  filter(Control.Type %in% c("earlobe sham","scapha sham"))-> df2

# Aggregate within study effect sizes


# continuous
aggES1 <- agg(id     = study.name,
              es     = yi,
              var    = vi,
              data   = df1,
              cor = .5,
              method = "BHHR")

# pulsed
aggES2 <- agg(id     = study.name,
              es     = yi,
              var    = vi,
              data   = df2,
              cor = .5,
              method = "BHHR")


# Aggregate within study effect sizes
aggES <- agg(id     = study.name,
             es     = yi,
             var    = vi,
             data   = df,
             cor = .5,
             method = "BHHR")

# Merging aggregated ES with original dataframe 
MA <- merge(x = aggES, y = df, by.x = "id", by.y = "study.name")
MA <- unique(setDT(MA) [sort.list(id)], by = "id")
MA <- with(MA, MA[order(MA$es)])

# Generate bayesmeta-object "bma", which stores all relevant results
bma <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$id, 
                 tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                 mu.prior = c("mean" = 0, "sd" = 1.5))


# use this to then plot funnel plot with the metameta package (by D. Quintana)
res <- rma.uni(es, var, data=MA, slab = id) 


# Summary tables
summary(bma)



# continuous
MA1 <- merge(x = aggES1, y = df1, by.x = "id", by.y = "study.name")
MA1 <- unique(setDT(MA1) [sort.list(id)], by = "id")
MA1 <- with(MA1, MA1[order(MA1$es)])

# pulsed
MA2 <- merge(x = aggES2, y = df2, by.x = "id", by.y = "study.name")
MA2 <- unique(setDT(MA2) [sort.list(id)], by = "id")
MA2 <- with(MA2, MA2[order(MA2$es)])

# Generate bayesmeta-object "bma", which stores all relevant results
bma1 <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                  tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                  mu.prior = c("mean" = 0, "sd" = 1.5))

# Generate bayesmeta-object "bma", which stores all relevant results
bma2 <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                  tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                  mu.prior = c("mean" = 0, "sd" = 1.5))

########## 
#Meta-regressions
#pulsed vs. continuous
pulsed <- as.numeric(MA$Stimulation.Type == "pulsed")
continuous <- as.numeric(MA$Stimulation.Type == "conventional")

protocol <- cbind(pulsed,continuous)

bma.reg <- bmr(y = MA$es,sigma = sqrt(MA$var), X = protocol, labels = MA$id, 
               tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
               mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA$id)

bma.reg$qpredict(p=0.5, x = c(-1,1))
bma.reg$pred.interval(level =0.95, x = c(-1,1))

forestplot(bma.reg, X.mean=rbind("pulsed" = c(1,0),
                                 "continuous" = c(0,1),
                                 "difference" = c(1,-1)), slab  = bma.reg$labels)

#pulsed: Sensation-matched vs. not
MA %>% 
  filter(Stimulation.Type == "pulsed") -> MA_pulsed

sensation_matched <- as.numeric(MA_pulsed$Intensity.fixed.or.matched == "sensation-matched intensity for sham and tVNS")
not_matched <-  as.numeric(MA_pulsed$Intensity.fixed.or.matched != "sensation-matched intensity for sham and tVNS")

matching <- cbind(sensation_matched,not_matched)

bma.reg_match <- bmr(y = MA_pulsed$es,sigma = sqrt(MA_pulsed$var), X = matching, labels = MA_pulsed$id, 
               tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
               mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA_pulsed$id)

bma.reg_match$qpredict(p=0.5, x = c(-1,1))
bma.reg_match$pred.interval(level =0.95, x = c(-1,1))

bma.reg_match$labels <- as.character(bma.reg_match$labels)
forestplot(bma.reg_match, X.mean=rbind("sensation-matched" = c(1,0),
                                 "not matched" = c(0,1),
                                 "difference" = c(1,-1)))



##########################

# funnel plot with contour
library(metaviz)
library(metafor)
library(puniform)

forest_plot <- forestplot.bayesmeta(bma)

png("Forest.png", units="in", width=8, height=5, res=600)
forestplot.bayesmeta(bma) ## Forest plot
dev.off()

forest_plot <- forestplot.bayesmeta(bma)
forest_plot 

forest_plot1 <- forestplot.bayesmeta(bma1)
forest_plot1 

forest_plot2 <- forestplot.bayesmeta(bma2)
forest_plot2 

funnel_plot <- funnel(res)
funnel_plot

funnel_contour <- funnel(res, level=c(90, 95, 99), 
                         shade=c("cornsilk", "orange", "firebrick"), 
                         back="white",
                         label = "out", offset= -8,
                         refline=0, legend="topright",
                         col = case_when(MA$Stimulation.Type == "continuous"~ "plum",
                                         MA$Stimulation.Type == "pulsed"~ "mediumseagreen"),
                         cex = 1.3)



##############################################
funnel_contour

puni_star(yi = MA$es, vi = MA$var, 
          side = "right")


viz_sunset(res, 
           
           method = "RE",
           sig_level = 0.05, # changes sign levels for power (background color)
           power_stats = TRUE, # adds stats at the bottom under the power line
           sig_contours = TRUE, # YAY it's this one!
           text_size = 5,
           point_size = 2,
           xlab = "Effect",
           contours = TRUE,
           power_contours =  "continuous")


# Generate function "robustness2"
## Plots Bayes factor robustness plot for paper

SD=1.5 # let's try by defining SD manually. I took this value from bma$mu.prior


  narrow <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                      tau.prior = bma$tau.prior, #tauprior, 
                      mu.prior = c("mean" = 0, "sd" = (SD/2))) 
  default <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                       tau.prior = bma$tau.prior, #tauprior, 
                       mu.prior = c("mean" = 0, "sd" = SD))
  wide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = bma$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = SD+1))
  ultrawide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                         tau.prior = bma$tau.prior, #tauprior,  
                         mu.prior = c("mean" = 0, "sd" = SD+2))
  BFS <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
  SDS <- c(SD/2,SD,SD+1,SD+2)
  names <- c("Narrow", "Default","Wide","Ultrawide")
  ggplot(data = NULL, aes(SDS,BFS)) +
    geom_hline(yintercept = 0.1, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 0.33, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 1, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 3, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 0.03, color = "grey", size = 0.3, alpha = .6) + 
    #geom_hline(yintercept = 30, color = "grey", size = 0.3, alpha = .6) + 
    geom_point(aes(SDS,BFS), alpha = .75) + 
    geom_point(aes(1.5, default$bayesfactor[1,2]), alpha = .75) +
    geom_text(aes(label = round(BFS, digits = 2)), hjust = -0.75) +
    scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+1.5)) +
    #scale_y_continuous(limits = c(-1,BFS[4]+10)) +
    scale_y_continuous(limits = c(-0,BFS[4]+0.8)) +
    geom_line(aes(SDS,BFS), alpha = .4) +
    geom_text(aes(label = names),vjust = -1, size = 6) +
    labs(x = "Standard deviations", y = "BF01") +
    theme_classic(18) +
    #    ggtitle("A.") +
    #annotate("text", label = "Strong evidence for H0", x = SDS[4] + 2, y = 20, size = 3) + #10-30
    #annotate("text", label = "Very strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2, size = 3) + #30-100
    annotate("text", label = "Moderate evidence for H0", x = SDS[4] + 1, y = 6.5, size = 3) + #3-10
    annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 1, y = 1.1, size = 3) + #1-3
    annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 1, y = 0.4, size = 3) + #0.33-1
    annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 1, y = 0.15, size = 3) + #0.10-0.33
    annotate("text", label = "Strong evidence for H1", x = SDS[4] + 1, y = 0.05, size = 3) + #0.03-0.1
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))




# try a different approach
library(abtest)

# Bayesian A/B test with default settings
ab <- ab_test(data = MA)

# plot robustness check (i.e., prior sensitivity analysis)
plot_robustness(ab)
