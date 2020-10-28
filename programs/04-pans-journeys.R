#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Feb 27, 2019                                          ##
##                                                             ##
#################################################################

remove(list=ls())

library(ggplot2)
library(ggpubr)

library(gridExtra)
library(scales)

timeline = TRUE
cohort = FALSE

source(here("programs", "fun-support.R"))
source(here("programs", "fun-plot.R"))



#--------------------------------------
# Step 2. Read clean dataFiles
#--------------------------------------

#---
#File with patient ID (de-ID), comorbidities, initial clinical presentation, etc.
patients.og <- read_csv(here("data", "patients.csv")) %>%
  column_to_rownames(var="id") 

#Select only a few columns
patients <- patients.og %>%
  select("age_onset", "age_1st_appt", "is_male", "NHW", "OCD", "foodprob", "anx", 
         "emotional", "mood", "agg", "sch", "reg", "sleep", "tics") %>%
  drop_na()

#---
#File with events
events.og <- read_csv(here("data", "medications.csv")) 

#Select only events related to patients in the previous file
events <- events.og %>%
  filter(id %in% rownames(patients)) %>%
  select("id", "medication", "start", "end")

#---
#File with clinical evaluations (outcomes)
outcomes.og <- read_csv(here("data","outcomes.csv"))

#Select only rows for the same patients listed in Profiles
outcomes <- outcomes.og %>%
  filter(id %in% rownames(patients)) %>%
  select("id", "gi_new", "cbiTotal", "daysSinceBirth")

#---
#Censoring
years = 2

#--------------------------------------
# Step 2. Clean data
#--------------------------------------


# Group by class of medication
medgroups <- lst(penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin"),
                 cephalosporin = c("cephalexin", "cefadroxil"), 
                 macrolide = c("azithromycin"),
                 nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin"),
                 corticosteroid.oral = c("prednisone", "maintenance prednisone", "decadron"),
                 corticosteroid.iv = c("solumedrol"),
                 immunoglobulins = c("ivig"),
                 dmard = c("rituximab", "methotrexate", "cellcept"))

medcolors= c("penicillin"="#66c2a5",
             "cephalosporin" = "#fc8d62",
             "macrolide" = "#8da0cb",
             "nsaid" = "#e7298a",
             "corticosteroid.oral" = "#a6d854",
             "corticosteroid.iv" = "#ffd92f",
             "antibody" = "green",
             "immunoglobulins" = "#e5c494",
             "dmard" = "#b3b3b3")


events <- events.og %>%
  cleanEvents(medgroups) %>%
  rightCensoring(years)





#--------------------------------------
# Step 3. Create timeline
#--------------------------------------

plist = list()



for(patient in rownames(patients)){
  
  #Get the events for one patient
  patIDs = which(events[,"id"] == patient)
  pat = events[patIDs, ]
  
  #profIndex = which(patients$id==patient)
  firstAppointment = ceiling(patients[patient,"age_1st_appt"]*365)
  firstOnset = ceiling(patients[patient, "age_onset"]*365)
  
  # add duration for all events and order
  pat$duration = pat[,"end"] - pat[,"start"]
  pat = pat[order(-pat$duration),]
  
  #------------------------
  # get the patient's scores
  #------------------------
  cliIDs  = which(outcomes[, "id"] == patient)
  cli = outcomes[cliIDs, ]
  cli = cli[complete.cases(cli),]
  
  
  model <- lm(gi_new~daysSinceBirth, cli)
  slope <- model$coefficients[2]
  
  mythreshold <- c("improving"="#67a9cf",
                   "stable"="gray80",
                   "worsening"="#ef8a62")
  gi_group = as.character(cut(slope, 
                              breaks=c(-1,-0.005,0.005,1), 
                              labels=names(mythreshold)))
  
  
  
  #Update values to adjust Clinical Values
  pat$start = pat$start + firstOnset
  pat$end = pat$end + firstOnset
  
  maxdays =max(pat$end, cli$daysSinceBirth)
  
  years = rev(paste("year", seq(1:ceiling(maxdays/365))))
  
  
  #------------------------
  # Plot
  #------------------------
  # predict <- predict(loess(gi_new~daysSinceBirth, data = cli))
  # cli$change <- cli$gi_new - predict
  # 
  #   ggplot(cli, aes(x=daysSinceBirth, y=gi_new, colour=change)) + 
  #   geom_line(aes(colour=gi_new), size=0.5) +
  #   geom_smooth(aes(colour=..y..), method="loess", size=1.5, se=FALSE) 

  
  gi_new <- ggplot(cli, aes(x=daysSinceBirth, y=gi_new)) +
    geom_smooth(aes(color=gi_group), method="lm", size=5, se=FALSE) +
    geom_smooth(color="gray70", method="loess", size=2, se=TRUE, alpha=0.1) +
    geom_line(size=0.5, colour="gray", linetype="dotted") +
    geom_point() + 
    geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
    geom_vline(xintercept = firstAppointment, linetype="dashed") +
    scale_color_manual(values=mythreshold,
                       limits = names(mythreshold),
                       name="Global Impairment Score") +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(maxdays/365)*365), 
                       breaks = seq(floor(firstOnset/365)*365,maxdays+365,365),
                       labels = paste("year", seq(floor(firstOnset/365)*365,maxdays+365,365)/365),
                       expand = c(0, 0.5)) +
    labs(subtitle="Global Impairment Score") +
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      legend.position="right", 
      #legend.position="none", 
      #legend.title = element_blank(),
      axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      #axis.ticks.x=element_blank(),
      axis.text.x=element_text(angle=90, vjust=0.5),
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      #Add a border
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  
  
  #------------------------
  # plot the patient timeline
  #------------------------

  
  
  #Plot
  gPat <- ggplot(pat) + 
    geom_segment(aes(x=start, xend=end, y=medication, yend=medication, colour=medication), 
                 size=8, lineend="butt") +
    scale_color_manual(values = medcolors,
                       limits = names(medcolors)) +
    geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
    geom_vline(xintercept = firstAppointment, linetype="dashed") +
    scale_y_discrete(limits = rev(names(medcolors))) +
    scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(maxdays/365)*365), 
                       breaks = seq(floor(firstOnset/365)*365,maxdays+365,365),
                       labels = paste("year", seq(floor(firstOnset/365)*365,maxdays+365,365)/365),
                       expand = c(0, 0.5)) +
    ggtitle(paste("Patient ", patient, sep="")) +
    labs(subtitle="Medication history") +
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      legend.position="right", 
      #legend.position="none", 
      legend.title = element_blank(),
      axis.title.x=element_blank(),
      #axis.text.x=element_text(angle=90, vjust=0.5),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      #Add a border
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  #-------
  gpanels <- ggarrange(gPat, gi_new, 
                       labels = c("A", "B"),
                       align = "v",
                       ncol = 1, nrow = 2, legend="right", common.legend = FALSE)
  
  var = paste("gPatient", patient, sep = "")
  
  #assign(var, gpanels)
  #plist[[patient]] = get(var)
  
  
  filepathname = paste(here("images", "journeys", 
                            paste0("patient", patient,".png", sep="")))
  ggexport(gpanels, filename=filepathname, height = 2000, width = 3000, res=300)
  
}

#Saving entire cohort
#gcohort <- ggarrange(plotlist=plist)

#ggexport(gcohort, filename="../../images/timeline-scores/cohort.png", height = 3000, width = 4000, res=300)


