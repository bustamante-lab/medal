###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda <arturolp@stanford.edu>
# Date: Aug 10, 2018
#
###############################################################
library(ggplot2)
#library(gridExtra)
#library(scales)
library(stringr)


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
plotPatientTimeline <- function(timeline, patient.label, firstAppointment, firstOnset, colors=mycolors){
  
  #Update values to adjust Clinical Values
  timeline$start = timeline$start + firstOnset
  timeline$end = timeline$end + firstOnset
  
  # Maximum years
  maxDay = max(timeline$end)
  
  #Plot
  g <- ggplot(timeline) + 
    geom_segment(aes(x=start, xend=end, y=medication, yend=medication, colour=medication), 
                 size=8, lineend="butt") +
    scale_x_continuous(limits = c(floor(firstOnset/365)*365,ceiling(maxDay/365)*365), 
                       breaks=seq(floor(firstOnset/365)*365,maxDay+365,365),
                       labels=paste("year", seq(floor(firstOnset/365)*365,maxDay+365,365)/365),
                       expand = c(0, 0.5)) +
    scale_y_discrete(limits = rev(names(colors))) +
    ggtitle(patient.label) +
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


unionPatients <- function(patient1, patient2){
  patientU = rbind(patient1, patient2)
  
  return(patientU)
}

#TO DO:
intersectPatients <- function(patient1, patient2){
  patientU = c()
  
  sequences = getSequences(patient1, patient2)
  
  patientID = paste(patient1$patientID[1], patient2$patientID[1], sep="-")
  
  medications = intersect(patient1$medication, patient2$medication)
  for(medication in medications){
    
    med1 = sequences[[medication]][1,]
    med2 = sequences[[medication]][2,]
    size = length(med1)
    
    
    #Assign a value where both are the same
    indSameLetter = intersect(which(med1==med2), which(med1 != '∅'))
    if(length(indSameLetter)>0){
      
      #Create an empty matrix
      combined = as.data.frame(matrix(rep('∅', size), nrow = 1), stringsAsFactors = FALSE)
      rownames(combined) = c(patientID)
      colnames(combined) = names(sequences[[medication]])
      
      #Assign medication character
      combined[, indSameLetter] = med1[,indSameLetter]
      
      #Reconstruct the format
      letterSequence = which(combined[1,] != '∅')
      
      current = letterSequence[1]
      start = as.numeric(str_replace(colnames(combined)[current], "day", ""))
      end = start
      for(i in 2:length(letterSequence)){
        ind = letterSequence[i]
        if(ind != (current+1)){
          
          #assign row to patientU
          end = as.numeric(str_replace(colnames(combined)[current], "day", ""))
          combinedMed = c(patientID, medication, start, end)
          patientU = rbind(patientU, combinedMed)
          
          #update start
          start = as.numeric(str_replace(colnames(combined)[ind], "day", ""))
        }
        current = ind
      }
      ind = letterSequence[length(letterSequence)]
      end = as.numeric(str_replace(colnames(combined)[ind], "day", ""))
      combinedMed = c(patientID, medication, start, end)
      patientU = rbind(patientU, combinedMed)
      
    }
    
  }
  
  colnames(patientU) = colnames(patient1)
  rownames(patientU) = NULL
  patientU = data.frame(patientU, stringsAsFactors=FALSE)
  patientU$start = as.numeric(as.character(patientU$start))
  patientU$end = as.numeric(as.character(patientU$end))
  
  
  return(patientU)
}

#TO DO:
averagePatients <- function(patient1, patient2){
  patientU = intersectPatients(patient1, patient2)
  
  #TO DO:
  medications = unique(c(patient1$medication, patient2$medication))
  for(medication in medications){
    numRows.p1 = length(which(patient1$medication == medication))
    numRows.p2 = length(which(patient2$medication == medication))
    
    print(paste(medication, numRows.p1, numRows.p2, sep=", "))
  }
  
  return(patientU)
}

