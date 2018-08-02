###############################################################
#
# Project: Patient Journey
# Author: Arturo Lopez Pineda <arturolp@stanford.edu>
# Date: May 25, 2018
#
###############################################################

remove(list=ls())

library(ggplot2)
library(gridExtra)
library(scales)

timeline = TRUE
cohort = FALSE

#--------------------------------------
# Step 1. Read file
#--------------------------------------

events = read.csv("2-medication/medsEvents.csv", stringsAsFactors = FALSE)
dictionary = read.csv("2-medication/medsDictonary.csv", row.names = 1)
clinical = read.csv("1-clinical/clinicHistoryDeIdentified.csv", stringsAsFactors = FALSE)

# Column identifiers
patientIDcol = 1
umlsCol = 2
startCol = 4
endCol = 5

# Merge Prednisone (burst and maintenance)
events[which(events$medication == "Prednisone burst"), "medication"] = "Prednisone"
events[which(events$medication == "Maintenance prednisone"), "medication"] = "Prednisone"


#--------------------------------------
# Step 2. Create timeline
#--------------------------------------

plist = list()

# Temporary: remove all events without start day
events = events[which(!is.na(events$start)),]

# Obtain all the patient IDs
patients = unique(events[,patientIDcol])

for(patient in patients){
  
  #Get the events for one patient
  patIDs = which(events[,patientIDcol] == patient)
  pat = events[patIDs, ]
  firstAppointmentPat = pat$firstAppointment[1]
  
  #Obtain the grouping
  pat$group = dictionary[as.character(pat[,umlsCol]), "Grouping"]
  
  #For end events == NA, add same day as start + 1
  pat[is.na(pat[, endCol]), endCol] = pat[is.na(pat[, endCol]), startCol] + 1
  
  # add duration for all events and order
  pat$duration = pat[,endCol] - pat[,startCol]
  pat = pat[order(-pat$duration),]
  
  #------------------------
  # get the patient's scores
  #------------------------
  cliIDs  = which(clinical[, "id"] == patient)
  cli = clinical[cliIDs, ]
  
  firstAppointment = cli$daysSinceBirth[1] - cli$daysSinceFirstAppointment[1]
  firstOnset = cli$daysSinceBirth[1] - cli$daysPostOnset[1]
  
  
  gi_score  = cli[,c("gi_new", "tdTotal", "diTotal", "phiTotal", "ehiTotal", "sriTotal", "cbiTotal", "daysSinceBirth")]
  
  gi_score$tdTotal = rescale(gi_score$tdTotal, to=c(0,100), from=c(0,20))
  gi_score$diTotal = rescale(gi_score$diTotal, to=c(0,100), from=c(0,20))
  gi_score$phiTotal = rescale(gi_score$phiTotal, to=c(0,100), from=c(0,16))
  gi_score$ehiTotal = rescale(gi_score$ehiTotal, to=c(0,100), from=c(0,20))
  gi_score$sriTotal = rescale(gi_score$sriTotal, to=c(0,100), from=c(0,20))
  gi_score$cbiTotal = rescale(gi_score$cbiTotal, to=c(0,100), from=c(0,96))
  
  gi_table = cbind(score=gi_score[,"tdTotal"], daysSinceBirth=gi_score[,"daysSinceBirth"], type=rep("tdTotal", dim(gi_score)[1]))
  gi_table = rbind(gi_table, cbind(score=gi_score[,"diTotal"], daysSinceBirth=gi_score[,"daysSinceBirth"], type=rep("diTotal", dim(gi_score)[1])))
  gi_table = rbind(gi_table, cbind(score=gi_score[,"phiTotal"], daysSinceBirth=gi_score[,"daysSinceBirth"], type=rep("phiTotal", dim(gi_score)[1])))
  gi_table = rbind(gi_table, cbind(score=gi_score[,"ehiTotal"], daysSinceBirth=gi_score[,"daysSinceBirth"], type=rep("ehiTotal", dim(gi_score)[1])))
  gi_table = rbind(gi_table, cbind(score=gi_score[,"sriTotal"], daysSinceBirth=gi_score[,"daysSinceBirth"], type=rep("sriTotal", dim(gi_score)[1])))
  gi_table = rbind(gi_table, cbind(score=gi_score[,"cbiTotal"], daysSinceBirth=gi_score[,"daysSinceBirth"], type=rep("cbiTotal", dim(gi_score)[1])))
  gi_table = as.data.frame(gi_table)
  gi_table$score = as.numeric(as.character(gi_table$score))
  gi_table$daysSinceBirth = as.numeric(as.character(gi_table$daysSinceBirth))

  
  gi_new <- ggplot(gi_score) + 
    scale_color_manual(values=c("gray25"), name="score") +
    geom_smooth(aes(x=daysSinceBirth, y=gi_new, colour="gi_new"), method="loess", size=1, se=TRUE, fill="gray82", alpha=0.5) +
    geom_line(aes(x=daysSinceBirth, y=gi_new), size=0.5, colour="gray", linetype="dotted") +
    geom_point(aes(x=daysSinceBirth, y=gi_new)) + 
    geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
    geom_vline(xintercept = firstAppointment, linetype="dashed") +
    scale_y_continuous(limits = c(0, 100)) +
    scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(6570/365)*365), 
                       breaks=seq(floor(firstOnset/365)*365,6570+365,365),
                       labels=paste("year", seq(floor(firstOnset/365)*365,6570+365,365)/365),
                       expand = c(0, 0.5)) +
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      legend.position="right", 
      #legend.position="none", 
      legend.title = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.text.x=element_text(angle=90, vjust=0.5),
      #axis.title.y=element_blank(),
      #axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      # Force the plot into a rectangle aspect ratio
      #aspect.ratio = 0.3,
      #Add a border
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  gi_components <- ggplot(gi_table) +
    scale_color_manual(values=c("#2121D9", "#9999FF", "#D92121", "#21D921", "#BB994D", "#FF9326"), name="score") +
    geom_smooth(aes(x=daysSinceBirth, y=score, color=type), method="loess", size=0.7, se=FALSE) +
    geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
    geom_vline(xintercept = firstAppointment, linetype="dashed") +
    scale_y_continuous(limits = c(0, 100)) + 
    scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(6570/365)*365), 
                       breaks=seq(floor(firstOnset/365)*365,6570+365,365),
                       labels=paste("year", seq(floor(firstOnset/365)*365,6570+365,365)/365),
                       expand = c(0, 0.5)) +
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      legend.position="right", 
      #legend.position="none", 
      legend.title = element_blank(),
      axis.title.x=element_blank(),
      axis.text.x=element_blank(),
      axis.ticks.x=element_blank(),
      #axis.text.x=element_text(angle=90, vjust=0.5),
      #axis.title.y=element_blank(),
      #axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      # Force the plot into a rectangle aspect ratio
      #aspect.ratio = 0.3,
      #Add a border
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  #grid.arrange(g, gi_new, gi_components)
    
  
  #------------------------
  # plot the patient timeline
  #------------------------
  #max = max(unique(pat2$year))
  max = 23
  years = rev(paste("year", seq(1:max)))
  grouping = c("anti-inflamatory", "corticosteroid", "antibody", "antibiotic")
  
  mycodes <- c("NSAID", "Prednisone", "IVIG", "Rituximab", "Solumedrol", "Augmentin", "Amoxicillin", "Cefadroxil", "Clindamycin", "Azithromycin", "Cephalexin", "stop medication")
  
  mycolors <- c(
    #oral therapies
    "NSAID"="deepskyblue", "Prednisone"="deepskyblue3", 
    #IV therapies
    "IVIG"="violetred1", "Rituximab"="hotpink", "Solumedrol"="maroon3", 
    #Antibiotics 1
    "Augmentin"="gold1", "Amoxicillin"="orange",
    #Antibiotics 2
    "Cefadroxil"="limegreen", "Clindamycin"="yellowgreen",
    #Antibiotics 3
    "Azithromycin"="tomato", "Cephalexin"="firebrick2", 
    #Other
    "stop medication"="antiquewhite3"
    #slateblue
  )
  
  #Update values to adjust Clinical Values
  pat$firstAppointment = pat$firstAppointment + firstOnset
  pat$start = pat$start + firstOnset
  pat$end = pat$end + firstOnset
  
  
  #Plot
  if(timeline == TRUE){
    if(cohort == FALSE){
      g <- ggplot(pat) + 
        geom_segment(aes(x=start, xend=end, y=medication, yend=medication, colour=medication), 
                     size=8, lineend="butt") +
        scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(6570/365)*365), 
                           breaks=seq(floor(firstOnset/365)*365,6570+365,365),
                           labels=paste("year", seq(floor(firstOnset/365)*365,6570+365,365)/365),
                           expand = c(0, 0.5)) +
        scale_y_discrete(limits = rev(mycodes)) +
        ggtitle(paste("Patient ", patient, sep="")) +
        scale_color_manual(values=mycolors) +
        geom_vline(xintercept = firstOnset, linetype="dotted", color="red") +
        geom_vline(xintercept = firstAppointment, linetype="dashed") +
        theme(#Add a title
          plot.title = element_text(hjust = 0.5, size=15),
          #Remove elements
          legend.position="right", 
          #legend.position="none", 
          legend.title = element_blank(),
          axis.title.x=element_blank(),
          axis.text.x=element_text(angle=90, vjust=0.5),
          axis.title.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Force the plot into a rectangle aspect ratio
          #aspect.ratio = 0.3,
          #Add a border
          panel.border = element_rect(colour = "black", fill=NA, size=1)
        )
      
      
      #Plot the group
      g.box <- ggplotGrob(g)
      g.line <- ggplotGrob(gi_new)
      g.line2 <- ggplotGrob(gi_components)
      
      #Force the saame widht
      maxwidths <- grid::unit.pmax(g.box$widths[2:6], g.line$widths[2:6], g.line2$widths[2:6])
      g.box$widths[2:6] <- as.list(maxwidths)
      g.line$widths[2:6] <- as.list(maxwidths)
      g.line2$widths[2:6] <- as.list(maxwidths)
      
      #Plot the grid
      grid.arrange(g.box, g.line, g.line2, nrow=3, heights=c(0.6,0.2,0.2))
      grid.arrange(g.box, g.line, nrow=2, heights=c(0.7,0.3))
      
      
      
      
    }
    else if(cohort == TRUE){
      g <- ggplot(pat) + 
        geom_segment(aes(x=start, xend=end, y=medication, yend=medication, colour=medication), 
                     size=3, lineend="butt") +
        #scale_x_continuous(limits = c(0,9000), 
        #                   expand = c(0, 0.5)) +
        scale_x_continuous(limits = c(0,ceiling(max(pat$end)/365)*365), 
                           breaks=seq(0,max(pat$end)+365,365),
                           labels=seq(0,max(pat$end)+365,365)/365, 
                           expand = c(0, 0.5)) +
        scale_y_discrete(limits = rev(mycodes)) +
        ggtitle(paste("Patient ", patient, sep="")) +
        scale_color_manual(values=mycolors) +
        geom_vline(xintercept = firstAppointment, linetype="dashed") +
        theme(#Add a title
          plot.title = element_text(hjust = 0.5, size=14),
          #Remove elements
          legend.position="none", 
          legend.title = element_blank(),
          axis.title.x=element_blank(),
          #axis.text.x=element_blank(),
          #axis.text.x=element_text(angle=90, vjust=0.5),
          #axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.background = element_blank(),
          # Force the plot into a rectangle aspect ratio
          aspect.ratio = 0.5,
          #Add a border
          panel.border = element_rect(colour = "black", fill=NA, size=1)
        )
    }
  }
  assign(paste("patient", patient, sep = ""), g)
  plist[[patient]] <- g
  
  if(cohort == FALSE){
    patient.timeline = paste("images/timeline/patient_",patient,".png", sep="")
    #patient.timeline = paste("images/categories/patient_",patient,".png", sep="")
    ggsave(patient.timeline, width = 8, height = 6, dpi=300)
  }
  
}


if(cohort == TRUE){
  do.call("grid.arrange", c(plist[c(1:159, 161:169)], ncol = 19))
}



gi_new = clinical[which(clinical$id == patient), "gi_new"]
