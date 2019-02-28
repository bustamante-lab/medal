#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Feb 27, 2019                                          ##
##                                                             ##
#################################################################

remove(list=ls())

library(ggplot2)
library(gridExtra)
library(scales)

timeline = TRUE
cohort = FALSE

source("support-functions.R")


#--------------------------------------
# Step 1. Read file
#--------------------------------------


#---
#File with patient ID (de-ID), comorbidities, initial clinical presentation, etc.
profiles = read.csv("../../clinical/data-matrix-profiles-cluster.csv", stringsAsFactors = FALSE)

#Variables for stratification
strats <- c("is_male", "NHW", "OCD", "foodprob", "anx", 
            "emotional", "mood", "agg", "sch", "reg", "sleep", "tics")

#Select only a few columns
profiles = profiles[,c("id", "age_onset", "age_1st_appt", strats, "cluster")]

#Get unique patient IDs ordered
patients = sort(unique(profiles$id))

#---
#File with events
events = read.csv("../../medication/meddrug-12-2018.csv", stringsAsFactors = FALSE)

#Only select rows for the same patients listed in Profiles
rows = which(events$id %in% patients)

#Select only a few columns
events = events[rows,c("id", "medication", "start", "end")]

#---
#File with clinical evaluations (outcomes)
outcomes = read.csv("../../clinical/data-matrix-outcomes.csv", stringsAsFactors = FALSE)

#Only select rows for the same patients listed in Profiles
rows = which(outcomes$id %in% patients)

#Select only a few columns
outcomes = outcomes[rows,c("id", "gi_new", "daysSinceBirth")]



# Step 2. Clean Data ------------------------------


# Group by class of medication
medgroups = vector()
medgroups$penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin")
medgroups$cephalosporin = c("cephalexin", "cefadroxil")
medgroups$macrolide = c("azithromycin")
medgroups$nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin")
medgroups$hydrocortisone = c("prednisone", "maintenance prednisone", "decadron", "solumedrol")
medgroups$antibody = c("rituximab", "ivig")
medgroups$dmard = c("plaquenil", "methotrexate", "cellcept")


events = cleanEvents(events, medgroups)


#--------------------------------------
# Step 3. Create timeline
#--------------------------------------

plist = list()

for(patient in patients){
  
  #Get the events for one patient
  patIDs = which(events[,"id"] == patient)
  pat = events[patIDs, ]
  
  profIndex = which(profiles$id==patient)
  firstAppointment = ceiling(profiles[profIndex,"age_1st_appt"]*365)
  firstOnset = ceiling(profiles[profIndex, "age_onset"]*365)
  
  # add duration for all events and order
  pat$duration = pat[,"end"] - pat[,"start"]
  pat = pat[order(-pat$duration),]
  
  #------------------------
  # get the patient's scores
  #------------------------
  cliIDs  = which(outcomes[, "id"] == patient)
  cli = outcomes[cliIDs, ]
  cli = cli[complete.cases(cli),]
  
  gi_new <- ggplot(cli) + 
    scale_color_manual(values=c("gray25"), name="score") +
    #geom_smooth(aes(x=daysSinceBirth, y=gi_new, colour="gi_new"), method="loess", size=1, se=TRUE, fill="gray82", alpha=0.5) +
    geom_smooth(aes(x=daysSinceBirth, y=gi_new), method="lm", size=5, se=FALSE, color="gray40", alpha=0.3) +
    geom_line(aes(x=daysSinceBirth, y=gi_new), size=0.5, colour="gray", linetype="dotted") +
    geom_point(aes(x=daysSinceBirth, y=gi_new)) + 
    geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
    geom_vline(xintercept = firstAppointment, linetype="dashed") +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(6570/365)*365), 
                       breaks = seq(floor(firstOnset/365)*365,6570+365,365),
                       labels = paste("year", seq(floor(firstOnset/365)*365,6570+365,365)/365),
                       expand = c(0, 0.5)) +
    labs(subtitle="Global Impairment Score") +
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      legend.position="right", 
      #legend.position="none", 
      legend.title = element_blank(),
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
  #max = max(unique(pat2$year))
  max = 23
  years = rev(paste("year", seq(1:max)))
  grouping = names(medgroups)
  
  mycodes <- c("penicillin"="#1b9e77",
               "cephalosporin"="#d95f02",
               "macrolide"="#7570b3",
               "nsaid"="#e7298a",
               "hydrocortisone"="#66a61e",
               "antibody"="#e6ab02",
               "dmard"="#a6761d")
  
  #Update values to adjust Clinical Values
  pat$start = pat$start + firstOnset
  pat$end = pat$end + firstOnset
  
  
  #Plot
  gPat <- ggplot(pat) + 
    geom_segment(aes(x=start, xend=end, y=medication, yend=medication, colour=medication), 
                 size=8, lineend="butt") +
    scale_color_manual(values = mycodes,
                       limits = names(mycodes)) +
    geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
    geom_vline(xintercept = firstAppointment, linetype="dashed") +
    scale_y_discrete(limits = rev(names(mycodes))) +
    scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(6570/365)*365), 
                       breaks = seq(floor(firstOnset/365)*365,6570+365,365),
                       labels = paste("year", seq(floor(firstOnset/365)*365,6570+365,365)/365),
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
                       ncol = 1, nrow = 2, legend="right", common.legend = TRUE)
  
  var = paste("gPatient", patient, sep = "")
  
  #assign(var, gpanels)
  #plist[[patient]] = get(var)
  
  
  filepathname = paste("../../images/timeline-scores/patient",patient,".png", sep="")
  ggexport(gpanels, filename=filepathname, height = 3000, width = 4000, res=300)
  
}

#Saving entire cohort
#gcohort <- ggarrange(plotlist=plist)

#ggexport(gcohort, filename="../../images/timeline-scores/cohort.png", height = 3000, width = 4000, res=300)


