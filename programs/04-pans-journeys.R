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

source("fun-support.R")
source("fun-pans.R")



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
  
  mycodes <- c("penicillin"="#1b9e77",
               "cephalosporin"="#d95f02",
               "macrolide"="#7570b3",
               "nsaid"="#e7298a",
               "hydrocortisone"="#66a61e",
               "antibody"="#e6ab02",
               "dmard"="#a6761d")
  
  
  #Plot
  gPat <- ggplot(pat) + 
    geom_segment(aes(x=start, xend=end, y=medication, yend=medication, colour=medication), 
                 size=8, lineend="butt") +
    scale_color_manual(values = mycodes,
                       limits = names(mycodes)) +
    geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
    geom_vline(xintercept = firstAppointment, linetype="dashed") +
    scale_y_discrete(limits = rev(names(mycodes))) +
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
  
  
  filepathname = paste("../../images/timeline-scores/patient",patient,".png", sep="")
  ggexport(gpanels, filename=filepathname, height = 3000, width = 4000, res=300)
  
}

#Saving entire cohort
#gcohort <- ggarrange(plotlist=plist)

#ggexport(gcohort, filename="../../images/timeline-scores/cohort.png", height = 3000, width = 4000, res=300)


