### Main Effect Meta Analysis ####

# In this script, I'll create plots including all three analyses (all appers, continuous only, pulsed only)


# Clear global environment and load required packages
rm(list=ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr) 

#setwd("Y:/Paper/Drafts/pulsedtaVNS_pupil_meta/code")

# Source helper functions
source("HelperFunctions.R")

# Load "df.RDa"

load(file = "df.RDa")
df %>% 
  mutate(study.name = fct_cross(study,factor(Study.No),sep = ": Study ")) %>% 
  filter(Control.Type %in% c("earlobe sham","scapha sham"))-> df


df1 <- df %>% filter(Stimulation.Type %in% c("continuous"))
df2 <- df %>% filter(Stimulation.Type %in% c("pulsed"))

# Aggregate within study effect sizes
aggES <- agg(id     = study.name,
             es     = yi,
             var    = vi,
             data   = df,
             cor = .5,
             method = "BHHR")

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

# Merging aggregated ES with original dataframe 
MA <- merge(x = aggES, y = df, by.x = "id", by.y = "study.name")
MA <- unique(setDT(MA) [sort.list(id)], by = "id")
MA <- with(MA, MA[order(MA$es)])

# continuous
MA1 <- merge(x = aggES1, y = df1, by.x = "id", by.y = "study.name")
MA1 <- unique(setDT(MA1) [sort.list(id)], by = "id")
MA1 <- with(MA1, MA1[order(MA1$es)])

# pulsed
MA2 <- merge(x = aggES2, y = df2, by.x = "id", by.y = "study.name")
MA2 <- unique(setDT(MA2) [sort.list(id)], by = "id")
MA2 <- with(MA2, MA2[order(MA2$es)])


# Generate bayesmeta-object "bma", which stores all relevant results
bma <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$id, 
                 tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                 mu.prior = c("mean" = 0, "sd" = 1.5))

# continuous
bma1 <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                 tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                 mu.prior = c("mean" = 0, "sd" = 1.5))

# pulsed
bma2 <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                 tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                 mu.prior = c("mean" = 0, "sd" = 1.5))

# sue this to then plot funnel plot with the metameta package (by D. Quintana)
res <- rma.uni(es, var, data=MA, slab = study) 
res1 <- rma.uni(es, var, data=MA1, slab = study)
res2 <- rma.uni(es, var, data=MA2, slab = study)


# Summary tables
summary(bma)
summary(bma1)
summary(bma2)

###################################
# ROBUSTNESS PLOT
###################################

## Plots Bayes factor robustness plot for paper

SD=1.5 # let's try by defining SD manually. I took this value from bma$mu.prior

narrow <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$id, 
                    tau.prior = bma$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2))) 
default <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$id, 
                     tau.prior = bma$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
wide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$id, 
                  tau.prior = bma$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
ultrawide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$id, 
                       tau.prior = bma$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))
BFS <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
SDS <- c(SD/2,SD,SD+1,SD+2)

# continous
narrow1 <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                    tau.prior = bma$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2))) 
default1 <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                     tau.prior = bma$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
wide1 <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                  tau.prior = bma$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
ultrawide1 <- bayesmeta(y = MA1$es,sigma = sqrt(MA1$var), labels = MA1$id, 
                       tau.prior = bma$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))
BFS1 <- c(narrow1$bayesfactor[1,2], default1$bayesfactor[1,2], wide1$bayesfactor[1,2], ultrawide1$bayesfactor[1,2])
#SDS <- c(SD/2,SD,SD+1,SD+2)

# pulsed
narrow2 <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                    tau.prior = bma$tau.prior, #tauprior, 
                    mu.prior = c("mean" = 0, "sd" = (SD/2))) 
default2 <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                     tau.prior = bma$tau.prior, #tauprior, 
                     mu.prior = c("mean" = 0, "sd" = SD))
wide2 <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                  tau.prior = bma$tau.prior, #tauprior, 
                  mu.prior = c("mean" = 0, "sd" = SD+1))
ultrawide2 <- bayesmeta(y = MA2$es,sigma = sqrt(MA2$var), labels = MA2$id, 
                       tau.prior = bma$tau.prior, #tauprior,  
                       mu.prior = c("mean" = 0, "sd" = SD+2))
BFS2 <- c(narrow2$bayesfactor[1,2], default2$bayesfactor[1,2], wide2$bayesfactor[1,2], ultrawide2$bayesfactor[1,2])
#SDS <- c(SD/2,SD,SD+1,SD+2)

# plot

names <- c("Narrow", "Default","Wide","Ultrawide")
p1 <- ggplot(data = NULL, aes(SDS,BFS)) +
  geom_hline(yintercept = 0.1, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 0.33, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 1, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 3, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 0.03, color = "grey", size = 0.3, alpha = .6) + 
  #overall
  geom_point(aes(SDS,BFS), alpha = .75) + 
  geom_point(aes(1.5, default$bayesfactor[1,2]), alpha = .75) +
  geom_text(aes(label = round(BFS, digits = 2)), hjust = -0.75) +
  scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+1.5)) +
  #continuous
  geom_point(aes(SDS,BFS1), alpha = .75) + 
  geom_point(aes(1.5, default1$bayesfactor[1,2]), alpha = .75) +
  #geom_text(aes(label = round(BFS1, digits = 2)), hjust = -0.5) +
  #pulsed
  geom_point(aes(SDS,BFS2), alpha = .75) + 
  geom_point(aes(1.5, default2$bayesfactor[1,2]), alpha = .75) +
  #geom_text(aes(label = round(BFS2, digits = 2)), hjust = -0.75) +
  scale_y_continuous(limits = c(0,BFS[4]+12)) +
  geom_line(aes(SDS,BFS), alpha = .4) +
  geom_line(aes(SDS,BFS1), alpha = .4, color = "plum") + # continuous
  geom_line(aes(SDS,BFS2), alpha = .4, color = "mediumseagreen") + # pulsed
  geom_text(aes(label = names),vjust = -1, size = 6) +
  labs(x = "Standard deviations", y = "BF01") +
  theme_classic(18) +
  #    ggtitle("A.") +
  annotate("text", label = "Strong evidence for H0", x = SDS[4] + 2, y = 20, size = 3) + #10-30
  #annotate("text", label = "Very strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2, size = 3) + #30-100
  annotate("text", label = "Moderate evidence for H0", x = SDS[4] + 1, y = 6.5, size = 3) + #3-10
  annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 1, y = 1.1, size = 3) + #1-3
  annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 1, y = 0.4, size = 3) + #0.33-1
  annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 1, y = 0.15, size = 3) + #0.10-0.33
  #annotate("text", label = "Strong evidence for H1", x = SDS[4] + 1, y = 0.05, size = 3) + #0.03-0.1
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  #theme(legend.position = "top") +
  scale_fill_discrete(name = "Stimulation", labels = c("All", "Continuous", "Pulsed"))


# trying with log() (natural log)
names <- c("Narrow", "Default","Wide","Ultrawide")
p2 <- ggplot(data = NULL, aes(SDS,log(BFS))) +
  geom_hline(yintercept = -2.3, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = -1.1, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 0, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 1.10, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = -3.5, color = "grey", size = 0.3, alpha = .6) + 
  #overall
  geom_point(aes(SDS,log(BFS)), alpha = .75) + 
  geom_point(aes(1.5, log(default$bayesfactor[1,2])), alpha = .75) +
  geom_text(aes(label = round(log(BFS), digits = 2)), hjust = -0.5) +
  scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+1.5)) +
  #continuous
  geom_point(aes(SDS,log(BFS1)), alpha = .75) + 
  geom_point(aes(1.5, log(default1$bayesfactor[1,2])), alpha = .75) +
  geom_text(aes(label = round(BFS1, digits = 2)), vjust = -16, hjust = -0.5) +
  #pulsed
  geom_point(aes(SDS,log(BFS2)), alpha = .75) + 
  geom_point(aes(1.5, log(default2$bayesfactor[1,2])), alpha = .75) +
  geom_text(aes(label = round(log(BFS2), digits = 2)), vjust = 17, hjust = -0.5) +
  scale_y_continuous(limits = c(-2,log(BFS[4])+3)) +
  geom_line(aes(SDS,log(BFS)), alpha = .9, size = 1) +
  geom_line(aes(SDS,log(BFS1)), alpha = .9, color = "plum", size = 1) + # continuous
  geom_line(aes(SDS,log(BFS2)), alpha = .9, color = "mediumseagreen", size = 1) + # pulsed
  geom_text(aes(label = names),vjust = -1, size = 6) +
  labs(x = "Standard deviations", y = "ln(BF01)") +
  theme_classic(18) +
  #    ggtitle("A.") +
  annotate("text", label = "Strong evidence for H0", x = SDS[4] + 2, y = 3, size = 3) + #10-30
  #annotate("text", label = "Very strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2, size = 3) + #30-100
  annotate("text", label = "Moderate evidence for H0", x = SDS[4] + 1, y = 1.87, size = 3) + #3-10
  annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 1, y = 0.1, size = 3) + #1-3
  annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 1, y = -0.91, size = 3) + #0.33-1
  annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 1, y = -1.9, size = 3) + #0.10-0.33
  annotate("text", label = "Strong evidence for H1", x = SDS[4] + 1, y = -3, size = 3) + #0.03-0.1
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  scale_fill_discrete(name = "Stimulation", labels = c("All", "Continuous", "Pulsed"))



# temp - I need to use BF10 instead of BF01
names <- c("Narrow", "Default","Wide","Ultrawide")
p3 <- ggplot(data = NULL, aes(SDS,log(1/BFS))) + 
  geom_hline(yintercept = -1.1, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 0, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = 1.10, color = "grey", size = 0.3, alpha = .6) + 
  geom_hline(yintercept = -4.6, color = "grey", size = 0.3, alpha = .6) + 
  #overall
  geom_point(aes(SDS,log(1/BFS)), alpha = .75) + 
  geom_point(aes(1.5, log(1/(default$bayesfactor[1,2]))), alpha = .75) +
  geom_text(aes(label = round(log(1/BFS), digits = 2)), vjust = -0.5, hjust = -0.1, size = 3.5) +
  scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+1.5)) +
  #continuous
  geom_point(aes(SDS,log(1/BFS1)), alpha = .75) + 
  geom_point(aes(1.5, log(1/(default1$bayesfactor[1,2]))), alpha = .75) +
  geom_text(aes(label = round(log(1/BFS1), digits = 2)), vjust = 9, hjust = -0.1, size = 3.5) +
  #pulsed
  theme_cowplot(14) +
  geom_point(aes(SDS,log(1/BFS2)), alpha = .75) + 
  geom_point(aes(1.5, log(1/(default2$bayesfactor[1,2]))), alpha = .75) +
  geom_text(aes(label = round(log(1/BFS2), digits = 2)), vjust = -12, hjust = -0.1, size = 3.5) +
  scale_y_continuous(limits = c(-4.3,log(1/BFS2[4])+3)) +
  geom_line(aes(SDS,log(1/BFS)), alpha = .6, size=1) +
  geom_line(aes(SDS,log(1/BFS1)), alpha = .6, color = "plum", size=1) + # continuous
  geom_line(aes(SDS,log(1/BFS2)), alpha = .6, color = "mediumseagreen", size=1) + # pulsed
  #geom_text(aes(label = names),vjust = 2, size = 3) +
  #geom_text(aes(label = names),vjust = -18, size = 3) +
  geom_text(aes(label = names),vjust = 13, size = 3) +
  labs(x = "Standard deviations", y = "ln(BF10)") +
  labs(colour = "stimulation type") +

  #    ggtitle("A.") +
  annotate("text", label = "Strong evidence for H1", x = SDS[4] + 2, y = 3, size = 3) + #10-30
  #annotate("text", label = "Very strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2, size = 3) + #30-100
  #annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 1, y = 1.87, size = 3) + #3-10
  #annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 1, y = 0.1, size = 3) + #1-3
  #annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 1, y = -0.91, size = 3) + #0.33-1
  #annotate("text", label = "Moderate evidence for H0", x = SDS[4] + 1, y = -1.9, size = 3) + #0.10-0.33
  #annotate("text", label = "Strong evidence for H0", x = SDS[4] + 1, y = -3, size = 3) + #0.03-0.1
  theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
  #theme(legend.position = "top") +
  scale_fill_discrete(name = "Stimulation", labels = c("All", "Continuous", "Pulsed"))

ggplot2::ggsave("Robustness.png",p3, width = 5.5, height = 4, units = "in",dpi = 300,bg = "white")
