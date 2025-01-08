#######################################################################################
################### EFFECT OF taVNS ON HRV - A BAYSEIAN META-ANALYSIS #################
#######################################################################################

################### Helper functions #################################################
#----
#setwd("D:/SynologyDrive/BON_general/SynologyDrive/Paper/Drafts/pulsedtaVNS_pupil_meta/code")
#setwd("Y:/Paper/Drafts/pulsedtaVNS_pupil_meta/code")


# Load files
load(file = "df.RDa")
load(file = "screened.RDa")

# Generate function "priorposteriorlikelihood.ggplot"
## Plots prior, posterior and likelihood distribution
priorposteriorlikelihood.ggplot <- function(bma, lowerbound, upperbound) {
  effect <- seq(lowerbound, upperbound, length = 200)
  colors <- c("posterior" = "#D55E00", "prior" = "#009E73", "likelihood" = "#0072B2")
  devAskNewPage(ask = FALSE)
  ggplot(data = NULL, aes(x=effect,y=bma$dposterior(mu = effect))) + 
    geom_line(aes(x = effect, y = bma$likelihood(mu = effect), col = "likelihood")) +
    geom_line(aes(x = effect, y = bma$dposterior(mu = effect), col = "posterior")) + 
    geom_line(aes(x = effect, y = bma$dprior(mu = effect), col = "prior")) + 
    geom_vline(xintercept = 0, col = "gray") + 
    theme_minimal_hgrid(12) + 
    labs(x = "effect μ", y = "probability density", color = "legend") + 
    scale_color_manual(values = colors)
}

# Generate function "tauprior.ggplot"
## Plots tau prior distribution
tauprior.ggplot <- function(bma) {
  effect <- seq(0, 2, length = 200)
  devAskNewPage(ask = FALSE)
  ggplot(data = NULL, aes(x=effect,y=bma$dprior(tau = effect))) + 
    geom_line(aes(x = effect, y = bma$dprior(tau = effect))) + 
    geom_vline(xintercept = 0, col = "gray") + 
    theme_minimal_hgrid(12) + 
    labs(x = "heterogeneity τ", y = "probability density") 
}

# Generate function "robustness"
## Plots Bayes factor robustness plot for shiny app
robustness <- function(MA,SD, tauprior) {
  narrow <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                      tau.prior = tauprior, 
                      mu.prior = c("mean" = 0, "sd" = (SD/2)))
  default <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                       tau.prior = tauprior, 
                       mu.prior = c("mean" = 0, "sd" = 1.5))
  user <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = tauprior, 
                    mu.prior = c("mean" = 0, "sd" = SD))
  wide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = tauprior, 
                    mu.prior = c("mean" = 0, "sd" = SD+1))
  ultrawide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                         tau.prior = tauprior, 
                         mu.prior = c("mean" = 0, "sd" = SD+2))
  BFS <- c(narrow$bayesfactor[1,2], user$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
  SDS <- c(SD/2,SD,SD+1,SD+2)
  
  # this uses BF01
  # names <- c("Narrow", "User","Wide","Ultrawide")
  # deflabel <- data.frame(Ref = "Default", val = 1.5, stringsAsFactors = F)
  # ggplot(data = NULL, aes(SDS,BFS)) +
  #   geom_hline(yintercept = 1, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_hline(yintercept = 3, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_hline(yintercept = 10, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_hline(yintercept = 30, color = "grey", size = 0.3, alpha = .6) + 
  #   #geom_hline(yintercept = -3.4, color = "grey", size = 0.3, alpha = .6) + 
  #   #geom_hline(yintercept = -4.61, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_vline(xintercept = 1.5, color = "#D55E00", size = 0.3, alpha = .4) +
  #   geom_point(aes(SDS,BFS), alpha = .75) + 
  #   geom_text(aes(label = round(BFS, digits = 2)), hjust = -0.75) +
  #   geom_point(aes(1.5, default$bayesfactor[1,2]), alpha = .75) +
  #   scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+3)) +
  #   scale_y_continuous(limits = c(0,BFS[4]+10)) +
  #   geom_line(aes(SDS,BFS), alpha = .4) +
  #   geom_text(aes(label = names),vjust = -1) +
  #   geom_text(mapping = aes(x = val, y = 0, label = Ref, hjust = -.2, vjust = -1), data = deflabel, color = "#D55E00") +
  #   labs(x = "standard deviations", y = "BF01") +
  #   theme_cowplot(12) +
  #   annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 2, y = 2) + #1-3
  #   annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 2, y = 0) + #0-1 
  #   annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 2, y = -2) + 
  #   annotate("text", label = "Strong evidence for H1", x = SDS[4] + 2, y = -3) + 
  #   annotate("text", label = "Very strong evidence for H1", x = SDS[4] + 2, y = -4) + 
  #   annotate("text", label = "Decisive evidence for H1", x = SDS[4] + 2, y = -4.9)  
  
  # but I need to use BF10 instead of BF01
  names <- c("Narrow", "Default","Wide","Ultrawide")
  deflabel <- data.frame(Ref = "Default", val = 1.5, stringsAsFactors = F)
  ggplot(data = NULL, aes(SDS,(1/BFS))) +
    geom_hline(yintercept = 0.10, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 0.33, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 1, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 3, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 10, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 30, color = "grey", size = 0.3, alpha = .6) + 
    geom_vline(xintercept = 1.5, color = "#D55E00", size = 0.3, alpha = .4) +
    geom_point(aes(SDS,(1/BFS)), alpha = .75) + 
    geom_point(aes(1.5, (1/(default$bayesfactor[1,2]))), alpha = .75) +
    geom_text(aes(label = round((1/BFS), digits = 2)), vjust = -0.5, hjust = -0.75) +
    scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+1.5)) +
  
    scale_y_continuous(limits = c(-0.1,(1/BFS[1])+1)) +
    geom_line(aes(SDS,(1/BFS)), alpha = .4, size=1) +
    geom_text(aes(label = names),vjust = -1, size = 6) +
    labs(x = "Standard deviations", y = "BF10") +
    theme_classic(18) +
    #    ggtitle("A.") +
    annotate("text", label = "Strong evidence for H1", x = SDS[4] + 2, y = 13, size = 3) + #10-30
    #annotate("text", label = "Very strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2, size = 3) + #30-100
    annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 1, y = 4, size = 3) + #3-10
    annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 1, y = 2, size = 3) + #1-3
    annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 1, y = 0.66, size = 3) + #0.33-1
    annotate("text", label = "Moderate evidence for H0", x = SDS[4] + 1, y = 0.2, size = 3) + #0.10-0.33
    annotate("text", label = "Strong evidence for H0", x = SDS[4] + 1, y = -0.02, size = 3) + #0.03-0.1
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
    #theme(legend.position = "top") +
    scale_fill_discrete(name = "Stimulation", labels = c("All", "Continuous", "Pulsed"))
}

# Generate function "robustness2"
## Plots Bayes factor robustness plot for paper
robustness2 <- function(MA,SD, tauprior) {
  narrow <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                      tau.prior = tauprior, 
                      mu.prior = c("mean" = 0, "sd" = (SD/2)))
  default <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                       tau.prior = tauprior, 
                       mu.prior = c("mean" = 0, "sd" = SD))
  wide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                    tau.prior = tauprior, 
                    mu.prior = c("mean" = 0, "sd" = SD+1))
  ultrawide <- bayesmeta(y = MA$es,sigma = sqrt(MA$var), labels = MA$study, 
                         tau.prior = tauprior, 
                         mu.prior = c("mean" = 0, "sd" = SD+2))
  BFS <- c(narrow$bayesfactor[1,2], default$bayesfactor[1,2], wide$bayesfactor[1,2], ultrawide$bayesfactor[1,2])
  SDS <- c(SD/2,SD,SD+1,SD+2)
  # names <- c("Narrow", "Default","Wide","Ultrawide")
  # ggplot(data = NULL, aes(SDS,BFS)) +
  #   geom_hline(yintercept = 0.1, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_hline(yintercept = 0.33, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_hline(yintercept = 1, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_hline(yintercept = 3, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_hline(yintercept = 10, color = "grey", size = 0.3, alpha = .6) + 
  #   #geom_hline(yintercept = 30, color = "grey", size = 0.3, alpha = .6) + 
  #   geom_point(aes(SDS,BFS), alpha = .75) + 
  #   geom_point(aes(1.5, default$bayesfactor[1,2]), alpha = .75) +
  #   geom_text(aes(label = round(BFS, digits = 2)), hjust = -0.75) +
  #   scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+3)) +
  #   #scale_y_continuous(limits = c(-1,BFS[4]+10)) +
  #   scale_y_continuous(limits = c(-1,BFS[4]+3.2)) +
  #   geom_line(aes(SDS,BFS), alpha = .4) +
  #   geom_text(aes(label = names),vjust = -1, size = 6) +
  #   labs(x = "Standard deviations", y = "BF01") +
  #   theme_classic(18) +
  #   #    ggtitle("A.") +
  #   annotate("text", label = "Strong evidence for H0", x = SDS[4] + 2, y = 20, size = 3) + #10-30
  #   #annotate("text", label = "Very strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2, size = 3) + #30-100
  #   annotate("text", label = "Moderate evidence for H0", x = SDS[4] + 2, y = 6.5, size = 3) + #3-10
  #   annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 2, y = 2, size = 3) + #1-3
  #   annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 2, y = 0.4, size = 3) + #0.33-1
  #   annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 2, y = 0.15, size = 3) + #0.10-0.33
  #   #annotate("text", label = "Strong evidence for H1", x = SDS[4] + 2, y = -0.3, size = 3) + #0.03-0.1
  #   theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0)))
  #   #theme(axis.title.y = element_text(margin = margin(t = 0, r = 5, b = 0, l = 0)))
  
  # temp - I need to use BF10 instead of BF01
  names <- c("Narrow", "Default","Wide","Ultrawide")
  ggplot(data = NULL, aes(SDS,(1/BFS))) +
    geom_hline(yintercept = -2.3, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = -1.1, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 0, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = 1.10, color = "grey", size = 0.3, alpha = .6) + 
    geom_hline(yintercept = -3.5, color = "grey", size = 0.3, alpha = .6) + 
    #overall
    geom_point(aes(SDS,(1/BFS)), alpha = .75) + 
    geom_point(aes(1.5, (1/(default$bayesfactor[1,2]))), alpha = .75) +
    geom_text(aes(label = round((1/BFS), digits = 2)), vjust = -0.5, hjust = -0.75) +
    scale_x_continuous(breaks = c(SDS, 1.5), limits = c(SDS[1]-.5,SDS[4]+1.5)) +
    #continuous
    geom_point(aes(SDS,(1/BFS1)), alpha = .75) + 
    geom_point(aes(1.5,(1/(default1$bayesfactor[1,2]))), alpha = .75) +
    geom_text(aes(label = round((1/BFS1), digits = 2)), vjust = 24, hjust = -0.5) +
    #pulsed
    geom_point(aes(SDS,(1/BFS2)), alpha = .75) + 
    geom_point(aes(1.5,(1/(default2$bayesfactor[1,2]))), alpha = .75) +
    geom_text(aes(label = round((1/BFS2), digits = 2)), vjust = -11, hjust = -0.5) +
    scale_y_continuous(limits = c(-3.5,(1/BFS[4])+3)) +
    geom_line(aes(SDS,(1/BFS)), alpha = .4, size=1) +
    geom_line(aes(SDS,(1/BFS1)), alpha = .4, color = "plum", size=1) + # continuous
    geom_line(aes(SDS,(1/BFS2)), alpha = .4, color = "mediumseagreen", size=1) + # pulsed
    geom_text(aes(label = names),vjust = -1, size = 6) +
    geom_text(aes(label = names),vjust = -8, size = 6) +
    geom_text(aes(label = names),vjust = 14, size = 6) +
    labs(x = "Standard deviations", y = "ln(BF10)") +
    labs(colour = "stimulation type") +
    theme_classic(18) +
    #    ggtitle("A.") +
    annotate("text", label = "Strong evidence for H1", x = SDS[4] + 2, y = 3, size = 3) + #10-30
    #annotate("text", label = "Very strong evidence for H0", x = SDS[4] + 2, y = (BFS[4]+10+30)/2, size = 3) + #30-100
    annotate("text", label = "Moderate evidence for H1", x = SDS[4] + 1, y = 1.87, size = 3) + #3-10
    annotate("text", label = "Anecdotal evidence for H1", x = SDS[4] + 1, y = 0.1, size = 3) + #1-3
    annotate("text", label = "Anecdotal evidence for H0", x = SDS[4] + 1, y = -0.91, size = 3) + #0.33-1
    annotate("text", label = "Moderate evidence for H0", x = SDS[4] + 1, y = -1.9, size = 3) + #0.10-0.33
    annotate("text", label = "Strong evidence for H0", x = SDS[4] + 1, y = -3, size = 3) + #0.03-0.1
    theme(axis.title.y = element_text(margin = margin(t = 0, r = 20, b = 0, l = 0))) +
    #theme(legend.position = "top") +
    scale_fill_discrete(name = "Stimulation", labels = c("All", "Continuous", "Pulsed"))
}

# Generate function "is_outlier"
## Check for outliers
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}
#Generate meta regression
meta_reg <- function(MA,df,Predictor){
   if (Predictor %in% c("Stimulation type", "Risk of bias","Sensation-matching")){
    
    if (Predictor == "Stimulation type"){
      pulsed <- as.numeric(MA$Stimulation.Type == "pulsed")
      continuous <- as.numeric(MA$Stimulation.Type == "conventional")
      X <- cbind(pulsed,continuous)
    } 
    
    if (Predictor == "Risk of bias"){
      high_risk <- as.numeric(MA$risk_of_bias == "high risk")
      not_high <- as.numeric(MA$risk_of_bias != "high risk")
      X <- cbind(high_risk,not_high)
    } 
    
    if (Predictor == "Sensation-matching"){
      matched <- as.numeric(MA$Intensity.fixed.or.matched == "sensation-matched intensity for sham and tVNS")
      not_matched <-  as.numeric(MA$Intensity.fixed.or.matched != "sensation-matched intensity for sham and tVNS")
      X <- cbind(matched,not_matched)
    }  
   
    bma.reg <- bmr(y = MA$es,sigma = sqrt(MA$var), X = X, labels = MA$id, 
                   tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                   mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA$id)
    
    bma.reg$qpredict(p=0.5, x = c(-1,1))
    bma.reg$pred.interval(level =0.95, x = c(-1,1))
    bma.reg$labels <- as.character(bma.reg$labels, digits = 1)
    
    if (Predictor == "Stimulation type"){
      forestplot(bma.reg, X.mean=rbind("pulsed" = c(1,0),
                                       "conventional" = c(0,1),
                                       "difference" = c(1,-1)), slab  = bma.reg$labels, digits=1) 
    }
    if (Predictor == "Risk of bias"){
      forestplot(bma.reg, X.mean=rbind("high risk" = c(1,0),
                                       "low or some risk" = c(0,1),
                                       "difference" = c(1,-1)), slab  = bma.reg$labels, digits=1)
    }
    
    if (Predictor == "Sensation-matching"){
      forestplot(bma.reg, X.mean=rbind("matched" = c(1,0),
                                       "not-matched" = c(0,1),
                                       "difference" = c(1,-1)), slab  = bma.reg$labels, digits=1)
    }
   }
   if (Predictor %in% c("Task-related (rest vs. task)","Stimulation intensity (tVNS)","Stimulation intensity (tVNS-sham)","Stimulation length (pulsed stimulation only)","Stimulation frequency","Pupil size outcome (tonic vs. phasic, for conventional stimulation only)")){
     if (Predictor == "Task-related (rest vs. task)"){
       #Get variable to aggegrate depending on task as well
       df %>% 
         mutate(study.name = fct_cross(study,factor(Study.No),sep = ": S ")) %>% 
         mutate(study.name = fct_cross(study.name,factor(Task.related!="rest",labels = c("no","yes")),sep = ", Task: ")) -> df_sub_reg
         ## Aggregate effect sizes
         aggES_reg <- agg(id     = study.name,
                          es     = yi,
                          var    = vi,
                          data   = df_sub_reg,
                          cor = .5,
                          method = "BHHR")
         ## Merging aggregated ES with original dataframe 
         MA_reg <- merge(x = aggES_reg, y = df_sub_reg, by.x = "id", by.y = "study.name") 
         MA_reg <- unique(setDT(MA_reg) [sort.list(id)], by = "id")
         MA_reg <- with(MA_reg, MA_reg[order(MA_reg$es)])

       task <- as.numeric(MA_reg$Task.related != "rest")
       no_task <-  as.numeric(MA_reg$Task.related == "rest")
       
       tasks <- cbind(no_task,task)
       
       bma.reg <- bmr(y = MA_reg$es,sigma = sqrt(MA_reg$var), X = tasks, labels = MA_reg$id, 
                      tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                      mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA_reg$id)
       
       bma.reg$qpredict(p=0.5, x = c(-1,1))
       bma.reg$pred.interval(level =0.95, x = c(-1,1))
       
       bma.reg$labels <- as.character(bma.reg$labels)
       forestplot(bma.reg, X.mean=rbind("no task" = c(1,0),
                                        "task (cognitive or PLR)" = c(0,1),
                                        "difference" = c(1,-1)), slab  = as.character(bma.reg$labels), digits=1)
     }
     if (Predictor == "Stimulation intensity (tVNS)"){
       #Get variable to aggerate depending on task as well
       df %>% 
         mutate(study.name = fct_cross(study,factor(Study.No),sep = ": S ")) %>% 
         mutate(study.name = fct_cross(study.name,factor(Mean.Intensity.Stim.mA.),sep = ", mA: ")) -> df_sub_reg
       ## Aggregate effect sizes
       aggES_reg <- agg(id     = study.name,
                        es     = yi,
                        var    = vi,
                        data   = df_sub_reg,
                        cor = .5,
                        method = "BHHR")
       ## Merging aggregated ES with original dataframe 
       MA_reg <- merge(x = aggES_reg, y = df_sub_reg, by.x = "id", by.y = "study.name") 
       MA_reg <- unique(setDT(MA_reg) [sort.list(id)], by = "id")
       MA_reg <- with(MA_reg, MA_reg[order(MA_reg$Mean.Intensity.Stim.mA)])
       
       X <-cbind("int"=1,"mA"=MA_reg$Mean.Intensity.Stim.mA-mean(MA_reg$Mean.Intensity.Stim.mA))   
          
       bma.reg <- bmr(y = MA_reg$es,sigma = sqrt(MA_reg$var), X = X, labels = MA_reg$id, 
                      tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                      mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA_regA$id)
       
       bma.reg$labels <- as.character(bma.reg$labels)
       forestplot(bma.reg, X.mean=rbind("Average tVNS response" = c(1,0),
                                        "Change in tVNS response with 1 mA" = c(0,1),
                                        "tVNS response at max intensity" = c(1,max(X[,2])),
                                        "tVNS response at min intensity" = c(1,min(X[,2]))), slab  = as.character(bma.reg$labels), digits=2)
     }
     if (Predictor == "Stimulation intensity (tVNS-sham)"){
       #Get variable to aggerate depending on task as well
       df %>% 
         mutate(diff.intensity = round((Mean.Intensity.Stim.mA. - Mean.Intensity.Sham.mA.), digits =1)) %>% 
         mutate(study.name = fct_cross(study,factor(Study.No),sep = ": S ")) %>% 
         mutate(study.name = fct_cross(study.name,factor(diff.intensity),sep = ", mA: ")) -> df_sub_reg
       ## Aggregate effect sizes
       aggES_reg <- agg(id     = study.name,
                        es     = yi,
                        var    = vi,
                        data   = df_sub_reg,
                        cor = .5,
                        method = "BHHR")
       ## Merging aggregated ES with original dataframe 
       MA_reg <- merge(x = aggES_reg, y = df_sub_reg, by.x = "id", by.y = "study.name") 
       MA_reg <- unique(setDT(MA_reg) [sort.list(id)], by = "id")
       MA_reg <- with(MA_reg, MA_reg[order(MA_reg$diff.intensity)])
       
       X <-cbind("int"=1,"delta.mA"=MA_reg$diff.intensity-mean(MA_reg$diff.intensity))   
       
       bma.reg <- bmr(y = MA_reg$es,sigma = sqrt(MA_reg$var), X = X, labels = MA_reg$id, 
                      tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                      mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA_regA$id)
       
       bma.reg$labels <- as.character(bma.reg$labels)
       forestplot(bma.reg, X.mean=rbind("Average tVNS response" = c(1,0),
                                        "Change in tVNS response with 1 mA\ndiff tVNS-sham" = c(0,1),
                                        "tVNS response at max delta tVNS-sham" = c(1,max(X[,2])),
                                        "tVNS response at min delta tVNS-sham" = c(1,min(X[,2]))), slab  = as.character(bma.reg$labels), digits=1)
     }
     if (Predictor == "Stimulation frequency"){
       #Get variable to aggerate depending on task as well
       df %>% 
         mutate(study.name = fct_cross(study,factor(Study.No),sep = ": S ")) %>% 
         mutate(study.name = fct_cross(study.name,factor(`Impulse.frequency (Hz)`),sep = ", Hz: ")) -> df_sub_reg
       ## Aggregate effect sizes
       aggES_reg <- agg(id     = study.name,
                        es     = yi,
                        var    = vi,
                        data   = df_sub_reg,
                        cor = .5,
                        method = "BHHR")
       ## Merging aggregated ES with original dataframe 
       MA_reg <- merge(x = aggES_reg, y = df_sub_reg, by.x = "id", by.y = "study.name") 
       MA_reg <- unique(setDT(MA_reg) [sort.list(id)], by = "id")
       MA_reg <- with(MA_reg, MA_reg[order(MA_reg$`Impulse.frequency (Hz)`)])
       
       X <-cbind("int"=1,"Hz"=MA_reg$`Impulse.frequency (Hz)`-25)   
       
       bma.reg <- bmr(y = MA_reg$es,sigma = sqrt(MA_reg$var), X = X, labels = MA_reg$id, 
                      tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                      mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA_regA$id)
       
       bma.reg$labels <- as.character(bma.reg$labels)
       forestplot(bma.reg, X.mean=rbind("tVNS reponse at 25 Hz" = c(1,0),
                                        "Change in tVNS response with 5 Hz" = c(0,5),
                                        "tVNS response at max frequency" = c(1,max(X[,2])),
                                        "tVNS response at min frequency" = c(1,min(X[,2]))), slab  = as.character(bma.reg$labels), digits=1)
     }
     if (Predictor == "Stimulation length (pulsed stimulation only)"){
       #Get variable to aggerate depending on task as well
       df %>% 
         filter(Stimulation.Type == "pulsed") %>% 
         mutate(study.name = fct_cross(study,factor(Study.No),sep = ": S ")) %>% 
         mutate(study.name = fct_cross(study.name,factor(Stimulation.time),sep = ", s: ")) -> df_sub_reg
       ## Aggregate effect sizes
       aggES_reg <- agg(id     = study.name,
                        es     = yi,
                        var    = vi,
                        data   = df_sub_reg,
                        cor = .5,
                        method = "BHHR")
       ## Merging aggregated ES with original dataframe 
       MA_reg <- merge(x = aggES_reg, y = df_sub_reg, by.x = "id", by.y = "study.name") 
       MA_reg <- unique(setDT(MA_reg) [sort.list(id)], by = "id")
       MA_reg <- with(MA_reg, MA_reg[order(MA_reg$Stimulation.time)])
       
       X <-cbind("int"=1,"s"=MA_reg$Stimulation.time-mean(MA_reg$Stimulation.time))   
       
       bma.reg <- bmr(y = MA_reg$es,sigma = sqrt(MA_reg$var), X = X, labels = MA_reg$id, 
                      tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                      mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA_regA$id)
       
       bma.reg$labels <- as.character(bma.reg$labels)
       forestplot(bma.reg, X.mean=rbind("tVNS reponse for average pulse length" = c(1,0),
                                        "Change in tVNS response with 1s" = c(0,1),
                                        "tVNS response at max pulse length" = c(1,max(X[,2])),
                                        "tVNS response at min length" = c(1,min(X[,2]))), slab  = as.character(bma.reg$labels), digits=1)
     }
     if (Predictor == "Pupil size outcome (tonic vs. phasic, for conventional stimulation only)"){
     #Get variable to aggerate depending on task as well
     df %>% 
         filter(Stimulation.Type == "continuous") %>% 
         mutate(study.name = fct_cross(study,factor(Study.No),sep = ": S")) %>% 
         mutate(study.name = fct_cross(study.name,factor(Pupil.Dilation.Type),sep = ", Pupil: ")) -> df_sub_reg
     ## Aggregate effect sizes
     aggES_reg <- agg(id     = study.name,
                      es     = yi,
                      var    = vi,
                      data   = df_sub_reg,
                      cor = .5,
                      method = "BHHR")
     ## Merging aggregated ES with original dataframe 
     MA_reg <- merge(x = aggES_reg, y = df_sub_reg, by.x = "id", by.y = "study.name") 
     MA_reg <- unique(setDT(MA_reg) [sort.list(id)], by = "id")
     MA_reg <- with(MA_reg, MA_reg[order(MA_reg$es)])
     
     tonic <- as.numeric(MA_reg$Pupil.Dilation.Type == "tonic")
     phasic <-  as.numeric(MA_reg$Pupil.Dilation.Type == "phasic")
     
     tasks <- cbind(tonic,phasic)
     
     bma.reg <- bmr(y = MA_reg$es,sigma = sqrt(MA_reg$var), X = tasks, labels = MA_reg$id, 
                    tau.prior = function(t) dhalfcauchy(t, scale = 0.5), FE =T,
                    mu.prior = c("mean" = 0, "sd" = 1.5), slab = MA_reg$id)
     
     bma.reg$qpredict(p=0.5, x = c(-1,1))
     bma.reg$pred.interval(level =0.95, x = c(-1,1))
     
     bma.reg$labels <- as.character(bma.reg$labels)
     forestplot(bma.reg, X.mean=rbind("tonic" = c(1,0),
                                      "phasic" = c(0,1),
                                      "difference" = c(1,-1)), slab  = as.character(bma.reg$labels), digits=1)
   }
   }
}

################### Helper functions for frequentist analysis ########################

# Generate function "myMareg
## Meta analysis regression with summary, forest and funnel plot
myMareg <- function(formula = es~1,
                    var = var,
                    data,
                    addfit_forest = T) {
  model <- mareg(formula = formula,
  fixed = TRUE,
  var = var,
  data = data,
  slab = data$study)

 # model <- rma(yi, vi, data=data, method="REML", weighted=FALSE)
  
  result <- list(summary(model),
                 confint(model),
                 metafor::forest.rma(x = model, showweights = T, addfit = addfit_forest,
                                     order = "obs",
                                     xlim=c(-10,8)),
                 funnel(model, xlab = "Observed outcome"))
  listnames <- c("Modelsummary", "Model-Fit", "Forestplot", "Funnelplot")
  names(result) <- listnames
  return(result)
  
}

# Generate function "myCoef"
## Coefficients table
myCoef <- function(coef){
  coef<-round(coef,digits = 3)
  DT::datatable(coef,
                rownames=c("Intercept"),
                colnames=c("b","S.E.","z","lower CI","upper CI","p"))
}

# Generate function "myFit"
## Model fit
myFit <- function(fit){
  fit<-round(fit,digits = 3)
  DT::datatable(fit)
}

# Generate function "myUni"
## Model uniqueness
myUni<-function(rand){
  uni<-as.data.frame(rand)
  uni$estimate<-round(uni$estimate,digits = 2)
  uni$ci.lb<-round(uni$ci.lb,digits = 2)
  uni$ci.ub<-round(uni$ci.ub,digits = 2)
  DT::datatable(uni,
                colnames = c("estimate","lower CI","upper CI"))
}


