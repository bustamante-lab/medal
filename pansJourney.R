###############################################################
#
# Project: Medal Algorithm
# Author: Arturo Lopez Pineda <arturolp@stanford.edu>
# Date: Apr 18, 2018
#
###############################################################

remove(list=ls())

library(ggplot2)
library(gridExtra)


#--------------------------------------
# Step 1. Read file
#--------------------------------------

events = read.csv("medication/medsEvents.csv", stringsAsFactors = FALSE)
dictionary = read.csv("medication/medsDictonary.csv", row.names = 1)

#Column identifiers
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
  
  #Obtain the colors and sizes
  #pat$color = dictionary[as.character(pat[,umlsCol]), "Color"]
  
  #Obtain start year
  pat$year.label = paste("year",ceiling(pat[, startCol]/365), sep=" ")
  pat$year = ceiling(pat[, startCol]/365)
  
  #Obtain end year
  pat$year.end.label = paste("year",ceiling(pat[, endCol]/365), sep=" ")
  pat$year.end = ceiling(pat[, endCol]/365)
  
  
  #order by time, adding start and end of segment
  pat2 = pat[order(pat[, startCol]),]
  
  #For end events == NA, add same day as start + 1
  pat2[is.na(pat2[, endCol]), endCol] = pat2[is.na(pat2[, endCol]), startCol] + 1
  

  # change all values to 1-year format
  for (i in 1:dim(pat2)[1]) {
    while (pat2[i, startCol] > 365) {
      pat2[i, startCol] = pat2[i, startCol] - 365
    }
    
    while (pat2[i, endCol] > 365) {
      pat2[i, endCol] = pat2[i, endCol] - 365
    }
  }
  
  # Add next year's event to current year for completion
    
  # add duration for all events and order
  pat2$duration = pat2[,endCol] - pat2[,startCol]
  pat2 = pat2[order(-pat2$duration),]
    
    #------------------------
    # plot the patient timeline
    #------------------------
    #max = max(unique(pat2$year))
  max = 23
    years = rev(paste("year", seq(1:max)))
    
    #mycolors <- c("lightseagreen", "plum2", "darksalmon", "papayawhip", "hotpink", "gold1", "deepskyblue", "tomato", "olivedrab1", "yellowgreen", "coral", "antiquewhite3")
    #mycodes <- c("NSAID", "Prednisone", "Solu-Medrol", "IVIg", "Rituximab", "Augmentin", "Azithromycin", "Cefadroxil", "Cephalexin", "Amoxicillin", "Clindamycin", "stop medication")
    
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
    
    
    g <- ggplot(pat2) + 
      geom_segment(aes(x=start, xend=end, y=year.label, yend=year.label, colour=medication), 
                   size=5, lineend="butt") + 
      #scale_colour_brewer(type="qual") +
      scale_y_discrete(limits = years) +
      ggtitle(paste("Patient ", patient, sep="")) +
      scale_color_manual(values=mycolors) +
      #scale_color_manual(values=mycolors, breaks=mycodes, labels=mycodes) +
      #labs(color = "UMLS") +
      #labs(fill = "UMLS") +
      theme(#Add a title
            #plot.title = element_text(hjust = 0.5, size=10),
            plot.title = element_text(hjust = 0.5, size=15),
            #Remove elements
            legend.position="bottom", 
            #legend.position="none", 
            legend.title = element_blank(),
            axis.title.x=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x=element_blank(),
            axis.title.y=element_blank(),
            #axis.text.y=element_blank(),
            axis.ticks.y=element_blank(),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.background = element_blank(),
            # Force the plot into a rectangle aspect ratio
            aspect.ratio = 1,
            #Add a border
            panel.border = element_rect(colour = "black", fill=NA, size=1)
          )
    
    assign(paste("patient", patient, sep = ""), g)
    plist[[patient]] <- g
    
    patient.timeline = paste("images/box/patient_",patient,".png", sep="")
    ggsave(patient.timeline, width = 10, height = 6, dpi=300)
    
}

#do.call("grid.arrange", c(plist[c(1:159, 161:164)], ncol = 15))
