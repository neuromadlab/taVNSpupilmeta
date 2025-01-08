#######################################################################################
################### taVNS and pupil dilation - A BAYESIAN META-ANALYSIS #################
#######################################################################################

################### Shiny App v.2 21.12.2020 SERVER ###################################

# Define server logic
server <- function(input, output) {
  # Create parameters reactive for eye-tracking parameter columns in study overview
  parameters <- reactive({input$stimtype})
  # Create MA reactive for all outputs
  MA <- reactive({
    ## Create subset based on chosen inclusion criteria
    df_sub <- df %>% 
      mutate(study.name = fct_cross(study,factor(Study.No),sep = ": Study "))  %>% 
      filter(Stimulation.Type %in% input$stimtype,
                            Sample %in% input$Sample,
                            Blindness %in% input$blind,
                            Publication.Year >= input$pubyear[1], Publication.Year <= input$pubyear[2],
                            Stimulation.Type %in% input$stimtype,
                            Stimulation.side %in% input$stimside,
                            Part.of.the.ear.stimulated %in% input$partear,
                            taVNS.Device %in% input$device,
                            Intensity.fixed.or.matched %in% input$intensitytype,
                            `Impulse.frequency (Hz)` >= input$impulsefreq[1], `Impulse.frequency (Hz)` <= input$impulsefreq[2],
                            Eye.Side.Measured %in% input$eyeside,
                            Task.related %in% input$task,
                            Pupil.Dilation.Type %in% input$dilationtype,
                            Control.Type %in% input$control,
                            Ambient.Light.Yes.or.No %in% input$light,
                            Stimulation.time >= input$stimduration[1], Stimulation.time <= input$stimduration[2],
                            study %in%input$included,
                            risk_of_bias %in% input$rob)
    if (input$sex_yes == "yes") {
    df_sub <- df_sub %>%
      filter(Percent.Females >= input$gender[1], Percent.Females <= input$gender[2])}

    if (input$age_yes == "yes") {
      df_sub <- df_sub %>% filter(Mean.Age >= input$age[1], Mean.Age <= input$age[2])}

    
    ## Aggregate effect sizes
    aggES <- agg(id     = study.name,
                 es     = yi,
                 var    = vi,
                 data   = df_sub,
                 cor = .5,
                 method = "BHHR")
    ## Merging aggregated ES with original dataframe 
    MA <- merge(x = aggES, y = df_sub, by.x = "id", by.y = "study.name") 
    MA <- unique(setDT(MA) [sort.list(id)], by = "id")
    MA <- with(MA, MA[order(MA$es)])
    
    
    
   # MA <- MA[complete.cases(MA$es), ]
    
  })
  
  df_sub_reg <- reactive({
    ## Create subset based on chosen inclusion criteria
    df_sub <- df %>% 
      mutate(study.name = fct_cross(study,factor(Study.No),sep = ": Study "))  %>% 
      filter(Stimulation.Type %in% input$stimtype,
             Sample %in% input$Sample,
             Blindness %in% input$blind,
             Publication.Year >= input$pubyear[1], Publication.Year <= input$pubyear[2],
             Stimulation.Type %in% input$stimtype,
             Stimulation.side %in% input$stimside,
             Part.of.the.ear.stimulated %in% input$partear,
             taVNS.Device %in% input$device,
             Intensity.fixed.or.matched %in% input$intensitytype,
             `Impulse.frequency (Hz)` >= input$impulsefreq[1], `Impulse.frequency (Hz)` <= input$impulsefreq[2],
             Eye.Side.Measured %in% input$eyeside,
             Task.related %in% input$task,
             Pupil.Dilation.Type %in% input$dilationtype,
             Control.Type %in% input$control,
             Ambient.Light.Yes.or.No %in% input$light,
             Stimulation.time >= input$stimduration[1], Stimulation.time <= input$stimduration[2],
             study %in%input$included,
             risk_of_bias %in% input$rob)
    if (input$sex_yes == "yes") {
      df_sub <- df_sub %>%
        filter(Percent.Females >= input$gender[1], Percent.Females <= input$gender[2])}
    
    if (input$age_yes == "yes") {
      df_sub <- df_sub %>% filter(Mean.Age >= input$age[1], Mean.Age <= input$age[2])}
      df_sub_reg <- df_sub
  })
  
  # Create bma reactive needed for all outputs
  bma <- reactive({
    ## Generate bayesmeta-object "bma" depending on tau prior chosen
    if (input$tauprior == "Half cauchy") {
      bma <- bayesmeta(y = MA()$es,sigma = sqrt(MA()$var), labels = MA()$id, 
                       tau.prior = function(t) dhalfcauchy(t, scale = input$scaletau), 
                       mu.prior = c("mean" = input$mupriormean, "sd" = input$mupriorsd))
    } else if (input$tauprior == "Half student t") {
      bma <- bayesmeta(y = MA()$es,sigma = sqrt(MA()$var), labels = MA()$id, 
                       tau.prior = function(t) dhalfnormal(t, scale = input$scaletau), 
                       mu.prior = c("mean" = input$mupriormean, "sd" = input$mupriorsd))
    } else {
      bma <- bayesmeta(y = MA()$es,sigma = sqrt(MA()$var), labels = MA()$id, 
                       tau.prior = input$tauprior, 
                       mu.prior = c("mean" = input$mupriormean, "sd" = input$mupriorsd))
    }
  })
  # Study overview panel
  output$study <- DT::renderDataTable({
    MA <- as.data.frame(MA())
    MA <- MA %>% mutate_if(is.numeric, round, digits=3)
#    parameters <- base::subset(MA, select = c(parameters()))
    criteria <- base::subset(MA, select = c(study, Study.No ,es, var, Stimulation.Type, Sample, Total.N, Mean.Age, Percent.Females, Blindness, Stimulation.side, Part.of.the.ear.stimulated, Intensity.fixed.or.matched, Stimulation.time, Control.Type, Pupil.Dilation.Type, risk_of_bias, Task.related))
    colnames(criteria) <- c("Paper", "Study number", "Hedges' g", "Variance", "Stimulation Type", "Sample", "Total N", "Mean age","Percent Women", "Blindness", "Stimulation side", "Part of the ear stimulated", "Intensity calibration" ,"Stimulation Pulse Duration (sec.)", "Control Type","Pupil Dilation Outcome","Risk of Bias", "Task")
    MAclean <- criteria
    #manually change output type for skora and Bömmer
    MAclean %>% 
      mutate(Outcome.2 = ifelse((Paper %in% c("Skora L., 2024") & `Study number` == 2)|(Paper %in% c("Keute M., 2019", "D'Agostini M., 2022", "Borges U., 2021")),"phasic and tonic",as.character(`Pupil Dilation Outcome`))) %>% 
      #mutate(Outcome.2 = ifelse(Paper == "Skora L., 2024" & `Study number` == 2,"phasic and tonic",as.character(`Pupil Dilation Outcome`))) %>% 
      mutate(`Pupil Dilation Outcome` = as.factor(Outcome.2)) %>% 
      mutate(Task.2 = ifelse(Paper %in% c( "D'Agostini M., 2022"),"rest and task",as.character(Task))) %>% 
      mutate(Task = as.factor(Task.2)) %>% 
      select(-c(Outcome.2,Task.2)) -> MAclean
      
    
    DT::datatable(MAclean, extensions = "FixedColumns",
                  options = list(dom = 't', 
                                 pageLength = nrow(MAclean),
                                 scrollX = T, 
                                 fixedColumns = list(leftColumns = 2)))
  })
  
  
  ## Warning message if 3 or less studies are included
  output$warning <- renderPrint({
    MA <- as.data.frame(MA())
    if (nrow(MA) < 4) {print('WARNING: With the chosen inclusion criteria, 3 or less studies will be included in the analysis.')}
  })
  # Outliers panel
  output$boxplot <- renderPlot({
    MAo <- MA() %>% tibble::rownames_to_column(var = "outlier") %>% mutate(is_outlier=ifelse(is_outlier(es), es, as.numeric(NA)))
    MAo$study[which(is.na(MAo$is_outlier))] <- as.numeric(NA)
    ggplot(MAo, aes(x = factor(0), es)) +
      geom_boxplot(outlier.size = 3.5, outlier.colour = "#D55E00", outlier.shape = 18, fill = "lightgrey") +
      geom_text(aes(label=study),na.rm = T, nudge_y = 0.02, nudge_x = 0.05) +
      stat_boxplot(geom="errorbar", width = 0.05) +
      scale_x_discrete(breaks = NULL) +
      xlab(NULL) + ylab("Hedges' g") +
      theme_minimal_hgrid(12)
  }, width = 600, height = 600)
  
  # Forest Plot panel
  output$forest <- renderPlot({
    forestplot.bayesmeta(bma(), xlab = "Hedges' g")
  })
  # Funnel Plot panel
  output$funnel <- renderPlot({
    funnel.bayesmeta(bma(), main = "", FE= T)
  })
  # Statistics panel
  output$bf <- renderPrint ({
    bma()$bayesfactor[1,]
  })
  output$summary <- renderPrint({
    bma()$summary
  })
  output$ML <- renderPrint({
    bma()$ML
  })
  output$MAP <- renderPrint({
    bma()$MAP[1,]
  })
  # Full texts screened panel
  output$screened <- DT::renderDataTable({
    datatable(screened, escape = F, options = list(columnDefs = list(list(
      targets = 2,
      render = JS(
        "function(data, type, row, meta) {",
        "return type === 'display' && data != null && data.length > 30 ?",
        "'<span title=\"' + data + '\">' + data.substr(0, 30) + '...</span>' : data;",
        "}")
    ))),
    class = "display")
  })
  # Additional plots panel
  output$evupdate <- renderPlot({
    priorposteriorlikelihood.ggplot(bma(), lowerbound = 0 - (input$mupriormean + 1) * 1.5, upperbound = 0 + (input$mupriormean + 1) * 1.5)
  }, width = 800)
  output$joint <- renderPlot({
    plot.bayesmeta(bma(), which=2, main = "")
  }, width = 800)
  output$taupriorplot <- renderPlot({
    tauprior.ggplot(bma())
  }, width = 800)
  
  # Bayes factor robustness plot panel
  output$warning2 <- renderPrint({
    print("WARNING: Plot will not be computed, because an improper τ prior was chosen. Proper τ priors are 'Half student t' and 'Half cauchy'.")})
  output$robustplot <- renderPlot({
    if (input$robust == "Yes" &
        input$tauprior == "Half cauchy") {
      robustness(MA(),SD = input$mupriorsd, tauprior = function(t) dhalfcauchy(t, scale = input$scaletau))
    } else if (input$robust == "Yes" &
               input$tauprior == "Half student t") {
      robustness(MA(),SD = input$mupriorsd, tauprior = function(t) dhalfnormal(t, scale = input$scaletau))
    } 
  }, width = 800)
 # Meta-regression
output$warning3 <- renderPrint({
  print("WARNING: Meta-regression was not selected so will not be calculated at this moment.")})

output$metareg_stimtype <-renderPlot({
if (input$MetaReg == "Yes" &
   "Stimulation type" %in% input$Predictors){
  meta_reg(MA(),df_sub_reg(),"Stimulation type")
  }
}, width = 800)

output$metareg_rob <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Risk of bias" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Risk of bias")
  }
}, width = 800)

output$metareg_sensmatch <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Sensation-matching" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Sensation-matching")
  }
}, width = 800)

output$metareg_task <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Task-related (rest vs. task)" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Task-related (rest vs. task)")
  }
}, width = 800)

output$metareg_stimintensity <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Stimulation intensity (tVNS)" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Stimulation intensity (tVNS)")
  }
}, width = 800)

output$metareg_stimintensitydiff <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Stimulation intensity (tVNS-sham)" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Stimulation intensity (tVNS-sham)")
  }
}, width = 1200)

output$metareg_stimfreq <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Stimulation frequency" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Stimulation frequency")
  }
}, width = 800)

output$metareg_stimlength <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Stimulation length (pulsed stimulation only)" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Stimulation length (pulsed stimulation only)")
  }
}, width = 800)

output$metareg_pupilout <-renderPlot({
  if (input$MetaReg == "Yes" &
      "Pupil size outcome (tonic vs. phasic, for conventional stimulation only)" %in% input$Predictors){
    meta_reg(MA(),df_sub_reg(),"Pupil size outcome (tonic vs. phasic, for conventional stimulation only)")
  }
}, width = 800)
}