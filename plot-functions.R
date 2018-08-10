###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda <arturolp@stanford.edu>
# Date: Aug 10, 2018
#
###############################################################
library(ggplot2)
library(gridExtra)
library(scales)


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

#TO DO:
getTimeline <- function(patient, firstAppointment, firstOnset){
  timeline = c()
  
  #Obtain the grouping
  pat$group = dictionary[as.character(pat[,umlsCol]), "Grouping"]
  
  #For end events == NA, add same day as start + 1
  pat[is.na(pat[, endCol]), endCol] = pat[is.na(pat[, endCol]), startCol] + 1
  
  # add duration for all events and order
  pat$duration = pat[,endCol] - pat[,startCol]
  pat = pat[order(-pat$duration),]
  
  #Update values to adjust Clinical Values
  pat$firstAppointment = pat$firstAppointment + firstOnset
  pat$start = pat$start + firstOnset
  pat$end = pat$end + firstOnset
  
  return(timeline)
}

#TO DO:
plotPatientTimeline <- function(patient, firstAppointment, firstOnset, mycolors=mycolors){
  timeline = getTimeline(patient, firstAppointment, firstOnset)
  
  max = max(unique(pat$year))
  years = rev(paste("year", seq(1:max)))
  
  #Plot
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
  
  return(g)
  
}


