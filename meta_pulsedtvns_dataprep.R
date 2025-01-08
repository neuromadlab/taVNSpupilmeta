################### Data preparation ##################################################
#---- 
# Clear global environment and load required packages
rm(list = ls())

if (!require("pacman")) install.packages("pacman")
pacman::p_load(bayesmeta, compute.es, cowplot, data.table, dplyr, esc, forestplot, 
               ggplot2, knitr, MAd, metafor, reactable, readr, readxl, rmarkdown, 
               R.rsp, stringr, tidyr) 

#setwd("Y:/Paper/Drafts/pulsedtaVNS_pupil_meta/code")
#setwd("D:/SynologyDrive/Paper/Drafts/pulsedtaVNS_pupil_meta/code")
#setwd("/Users/lisa/Desktop/code_tVNS_pupil_meta")

# Import Master xlsx file, in which articles were coded
Master <- read_excel("master.xlsx")

# Generate column "study" serving as publication identifier (last name of 1st author + publication year)
firstauthor <- word(string = Master$Author, sep = ",") #Grab last name of the 1st author
study<-as.data.frame(cbind(firstauthor,Master$Publication.Year)) #Dataframe with 1st authors last name and publication year
study <- study %>% unite(study, sep = ", ") #Unite it into 1 column
Master <- cbind(study,Master)

# # Table showing reasons for exclusions
 excluded <- Master %>% filter(Master$Included == 0)
 excluded <- unique(setDT(excluded) [sort.list(study)], by = "DOI")
 reasons_exc <- table(excluded$Reason.for.exclusion)

# Generate data frame "df" consisting of all studies 
df <- Master%>% filter(Master$Included==1)

# Remove irrelevant columns
# df <- base::subset(df, select = -c(Reason.for.exclusion, Included, Pages, Publication.Title, 
#                                    Issue, Volume))

#df <- filter(df, SD_intervention != 0 |  SD_control != 0) # here i filter out everything that doesn't have a SD

# # Calculate effect sizes (Hedges' g) and store them in "es"
es <- as.data.frame(es<-escalc(m1i = df$mean_intervention,
                               sd1i =df$SD_intervention,
                               n1i =df$N_Intervention,
                               m2i = df$mean_control,
                               sd2i = df$SD_control,
                               n2i = df$N_Control,
                               measure = "SMD",
                               yi, sei, vi))

### Manual change!
# Add missing ES manually (Zhu, 2022; Villani, 2022)
#Zhu, 2022
#esc_t(2.22,0.31,49,49,49,"g")
#es[111,1:2] <- c(0.6241,0.0857) #check that we're using the correct row!
#Villani, 2022
esc_chisq(0.16,0.68,36,"g")
es[110,1:2] <- c(0.1307,0.1116)

df<-cbind(es,df)


# Generate factor variables
df$Design<-factor(df$Design,levels = c(0,1),labels = c("between","within"))
df$Sample<-factor(df$Sample,levels =c(0,1),labels=c("healthy","high worriers"))
df$Blindness<-factor(df$Blindness,levels=c(0,1),labels=c("double blind","single blind"))
df$Age.Category<-factor(df$Age.Category,levels=c(0),labels=c("adults"))
df$Country<-factor(df$Country,levels=c(0:7),labels=c("USA","China","UK","Netherlands","Belgium","Germany","Italy","Israel"))
df$Stimulation.Type<-factor(df$Stimulation.Type,levels=c(0,1),labels=c("conventional","pulsed"))
df$Stimulation.side<-factor(df$Stimulation.side,levels=c(0,1),labels=c("right","left"))
df$Part.of.the.ear.stimulated<-factor(df$Part.of.the.ear.stimulated,levels=c(0:4),labels=c("cymba concha","tragus","ear canal", "external acoustic meatus (EAM)","cymba concha and tragus"))
df$Control.Type<-factor(df$Control.Type,levels=c(0:4),labels=c("earlobe sham","no stimulation baseline","no stimulation sham","scapha sham","stimulation OFF phase"))
df$taVNS.Device<-factor(df$taVNS.Device,levels=c(0,1),labels=c("Cerbomed/tVNS technologies","other"))
df$Intensity.fixed.or.matched<-factor(df$Intensity.fixed.or.matched,levels=c(0:2),labels=c("sensation-matched intensity for sham and tVNS","fixed intensity for sham and tVNS","other"))
df$Task.related<-factor(df$Task.related,levels=c(0:2),labels=c("rest","task","pupillary light reflex"))
df$Pupil.Dilation.Type<-factor(df$Pupil.Dilation.Type,levels=c(0:1),labels=c("tonic","phasic"))
df$Eye.Tracker.Device<-factor(df$Eye.Tracker.Device,levels=c(0:4),labels=c("RED 500 eye tracker","SensoMotoric Instruments Eye-tracking glasses","Tobii T120 eye tracker","SR-Research Eyelink 1000","Sirius 3D Rotating Scheimpflug camera topographic system"))
df$Eye.Side.Measured<-factor(df$Eye.Side.Measured,levels=c(0:4),labels=c("right eye","left eye","both eyes","dominant eye","not specified"))
df$Colour.Fixation.Cross<-factor(df$Colour.Fixation.Cross,levels=c(0:6),labels=c("black","gray","white","islouminant colors","salmon","green","black/white"))
df$Colour.Orientation.Cross<-factor(df$Colour.Orientation.Cross,levels=c(0,1),labels=c("red","no cross"))
df$Colour.Background<-factor(df$Colour.Background,levels=c(0:3),labels=c("black","gray","isoluminant colors","blue"))
df$Ambient.Light.Yes.or.No<-factor(df$Ambient.Light.Yes.or.No,levels=c(0:2),labels=c("yes","no","no information"))
df$Correct.response.or.error<-factor(df$Correct.response.or.error,levels=c(0,1),labels=c("error","correct response"))
df$Direction_of_effect<-factor(df$Direction_of_effect,levels=c(0:2),labels=c("null","positive","negative"))
df$study<-as.factor(df$study)
df$risk_of_bias<-factor(df$`risk of bias Assessment`,labels=c("low risk","some concerns","high risk"))

########################################

### Manual change!
df$Publication.Year <- as.numeric(df$Publication.Year)

# #Adds a ' in front of numbers to prevent transformation into date format in Excel. Make sure to remove before using output!!!
#df$Mean.Age<-paste0("'",df$Mean.Age)
df$Percent.Females<-as.numeric(df$Percent.Females)
df$Mean.Age<-as.numeric(df$Mean.Age)
df$Mean.Intensity.Stim.mA.<-as.numeric(df$Mean.Intensity.Stim.mA.)
df$Mean.Intensity.Sham.mA.<-as.numeric(df$Mean.Intensity.Sham.mA.)
# df$Mean.Intensity.Stim.mA.<-paste0("'",df$Mean.Intensity.Stim.mA.)
# df$Mean.Intensity.Sham.mA.<-paste0("'",df$Mean.Intensity.Sham.mA.)
#df$Stimulation.duration.sec.<-paste0("'",df$Stimulation.duration.sec.)
df$Stimulation.duration.sec.<-as.numeric(df$Stimulation.duration.sec.)
df$mean_intervention<-as.numeric(df$mean_intervention)
df$mean_control<-as.numeric(df$mean_control)
df$mean_diff<-as.numeric(df$mean_diff)
df$SD_intervention<-as.numeric(df$SD_intervention)
df$SD_control<-as.numeric(df$SD_control)
df$SD_diff<-as.numeric(df$SD_diff)

# df$mean_intervention<-paste0("'",df$mean_intervention)
# df$mean_control<-paste0("'",df$mean_control)
# df$mean_diff<-paste0("'",df$mean_diff)
# df$SD_intervention<-paste0("'",df$SD_intervention)
# df$SD_control<-paste0("'",df$SD_control)
# df$SD_diff<-paste0("'",df$SD_diff)

df <-  df[complete.cases(df$yi), ] 
df <- df[complete.cases(df$vi), ] 

# manually remove D#agostini 2023 where applied intensity is = 0.2mA
df <- df %>%
  filter((Mean.Intensity.Stim.mA.!= c(0, 0.2)) %>% replace_na(TRUE))


# Add "https://doi.org/" infront of dois
Master$DOI <- paste0("https://doi.org/", Master$DOI)
Master$DOI <- paste0("<a href='", Master$DOI,"'target='_blank'>",Master$DOI,"</a>")

# Write "df" as .csv
write.csv2(df, "df.csv")

# Save "df" as .RDa
save(df, file = "df.RDa")

# Save overview of all screened studies
screened <- Master[,c(1,7,9,14:15)]
screened$Included <- factor(screened$Included,levels = c(0,1), labels = c("no","yes"))
screened <- unique(setDT(screened) [sort.list(study)], by = "DOI")
colnames(screened) <- c("Study", "Title", "DOI", "Design", "Sample")
save(screened, file = "screened.RDa")

