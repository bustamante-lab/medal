#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Feb 28, 2019                                          ##
##                                                             ##
#################################################################

library(ggplot2)
library(stringr)
library(reshape2)


mycolors <- c("penicillin"="#1b9e77",
              "cephalosporin"="#d95f02",
              "macrolide"="#7570b3",
              "nsaid"="#e7298a",
              "hydrocortisone"="#66a61e",
              "antibody"="#e6ab02",
              "dmard"="#a6761d")





plotTimeSeriesDrug <- function(cluster, events, profiles, medcolors=mycolors, medgroups){
  #Get the events for patients in a cluster
  patIDs = profiles[which(profiles[,"cluster"] == cluster), "id"]
  eveIDs = which(events[,"id"] %in% patIDs)
  pat = events[eveIDs, ]
  
  #Convert all days into months
  pat$start = round(pat$start/daysPerMonth)
  pat$end = round(pat$end/daysPerMonth)
  
  #Initialized matrix by medication
  drug <- matrix(rep(0, n*m), ncol=n, nrow=m)
  rownames(drug) = names(medgroups)
  colnames(drug) = seq(1:n)
  
  #Add one event if drug used in that month
  for(med in names(medgroups)){
    medIDs = which(pat[,"medication"]==med)
    clusterEvents = pat[medIDs,]
    for(i in 1:dim(clusterEvents)[1]){
      x = clusterEvents[i,"medication"]
      start = clusterEvents[i,"start"]
      end = clusterEvents[i,"end"]
      
      drug[x, start:end] = drug[x, start:end] + 1
    }
  }
  
  #Normalize by the number of patients in the cluster
  drug = drug/length(patIDs)
  
  
  #Reshape the matrix (Melt)
  melted_drug <- melt(drug)
  colnames(melted_drug) = c("med", "month", "value")
  melted_drug[which(melted_drug$value>1), "value"]=1 #avoiding duplicates
  
  #Remove all 0 values
  melted_drug = melted_drug[-which(melted_drug$value == 0),]
  
  # d <- data.frame(x=rep(1:20, 5), y=rnorm(100, 5, .2) + rep(1:5, each=20), z=rep(1:20, 5), grp=factor(rep(1:5, each=20)))
  # 
  # ggplot(d) +
  #   geom_path(aes(x, y, group=grp, alpha=z, color=grp), size=2)
  
  
  #Plot heatmap
  gPaths <- ggplot(melted_drug, aes(x=month, y=value)) +
    geom_path(aes(color = med), size=5, alpha=0.8) +
    geom_smooth(method="loess", color="gray30", se=FALSE, size=1.2, linetype = "dashed") +
    #Group by medication
    facet_wrap(.~med, ncol=1, scales="free_y") +
    #Reshape the scales
    scale_y_continuous(limits=c(0,1),
                       breaks = c(0,0.5,1),
                       labels = rev(c("all patients", "some", "none"))) +
    scale_x_continuous(breaks=c(seq(0, years, 1)*12),
                       labels=c(paste("year", seq(0,years, 1))))+
    #custom colours
    scale_color_manual(values=medcolors) +
    #Add labels
    labs(title=paste("Medication usage in cluster", cluster), x="Years of follow-up") +
    #set base size for all font elements
    theme_bw()+
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      legend.position="none", 
      #legend.position="none", 
      #legend.title = element_blank(),
      #axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      #axis.ticks.x=element_blank(),
      axis.text.x=element_text(angle=90, vjust=0.5),
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      #Facets
      strip.background = element_blank(),
      strip.text = element_text(size = 12, hjust=0)
    )
  
  return(gPaths)
}


plotCoOcurrenceTriangle <- function(cluster, events, profiles, medications){
  #Get the events for patients in a cluster
  patIDs = profiles[which(profiles[,"cluster"] == cluster), "id"]
  eveIDs = which(events[,"id"] %in% patIDs)
  pat = events[eveIDs, ]
  patients = unique(pat$id)
  
  #Convert all days into months
  pat$start = round(pat$start/daysPerMonth)
  pat$end = round(pat$end/daysPerMonth)
  
  #Create an empty matrix
  interactions = vector()
  
  #Fill in values
  for(col in 1:length(medications)){
    for(row in col:length(medications)){
      
      medsRow = medications[[row]]
      medsCol = medications[[col]]
      value = 0
      
      for(p in patients){
        indRow = which((events$id %in% p) & (events$medication %in% medsRow))
        indCol = which((events$id %in% p) & (events$medication %in% medsCol))
        if(length(indRow)>0 & length(indCol)>0){
          value = value + 1
        }
      }
      value = round(value/length(patients), digits=1)
      interactions = rbind(interactions, 
                           c(medications[col], medications[row], value))
    }
  }
  
  colnames(interactions) = c("A", "B", "value")
  interactions = as.data.frame(interactions, stringsAsFactors = TRUE)
  interactions$value <- as.numeric(as.character(interactions$value))

  
  
  # interactions = as.data.frame(matrix("",length(medications),length(medications)), 
  #                              stringsAsFactors = FALSE)
  # rownames(interactions) = medications
  # colnames(interactions) = medications
  # 
  # #Fill in values
  # for(col in 1:length(medications)){
  #   for(row in col:length(medications)){
  #     medsRow = medications[[row]]
  #     medsCol = medications[[col]]
  #     pats = 0
  #     for(p in patients){
  #       indRow = which((events$id %in% p) & (events$medication %in% medsRow))
  #       indCol = which((events$id %in% p) & (events$medication %in% medsCol))
  #       if(length(indRow)>0 & length(indCol)>0){
  #         #pats = pats + (length(indRow) + length(indCol))/2
  #         pats = pats + 1
  #       }
  #     }
  #     interactions[row,col]=round(pats/length(patients), digits=1)
  #   }
  # }
  
  #print(interactions)
  
  
  
  #Plot heatmap
  gTri <- ggplot(interactions, aes(x=A, y=B)) +
    geom_tile(aes(alpha=value), fill="steelblue", color="gray33") +
    geom_text(aes(label=value)) +
    scale_y_discrete(limits=rev(medications)) +
    scale_x_discrete(limits=medications) +
    ggtitle(paste("Co-occurrence in cluster", cluster)) +
    theme(#Add a title
      plot.title = element_text(hjust = 0, size=15),
      #Remove elements
      legend.position="none", 
      #legend.position="none", 
      #legend.title = element_blank(),
      axis.text=element_text(size=10),
      axis.title.x=element_blank(),
      axis.ticks.x=element_blank(),
      axis.text.x=element_text(angle=90, hjust=1),
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  return(gTri)
}


plotScores <- function(cluster, outcomes, profiles, maxDay){
  
  #Get the scores for patients in a cluster
  patIDs = profiles[which(profiles[,"cluster"] == cluster), "id"]
  outIDs = which(outcomes[,"id"] %in% patIDs)
  pat = outcomes[outIDs, ]
  patients = patIDs
  
  pat$id <- as.character(pat$id)
  
  #Plot heatmap
  gScore <- 
    ggplot(data = pat, aes(x = daysSinceBirth,  y = gi_new)) +
    #geom_point(aes(group=id, color=id), size=1)+
    geom_line(aes(group=id, color=id), size=0.5, linetype="dashed")+
    geom_smooth(aes(group=id, color=id), method = "loess", size=1, se=FALSE) +
    geom_smooth(method = "loess", size=2, se=FALSE, color="gray50") +
    scale_y_continuous(limits = c(0, 100), expand=c(0,0)) +
    scale_x_continuous(limits = c(0, maxDay) , 
                       breaks=seq(0,(maxDay+365),365),
                       labels=paste("year", seq(0,(maxDay+365),365)/365),
                       expand = c(0, 0)) +
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      #legend.position="right", 
      legend.position="none", 
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
  
  return(gScore)
}


#TO DO:
plotPatientTimeline <- function(timeline, patient.label, firstAppointment, firstOnset, colors=mycolors){
  
  #Update values to adjust Clinical Values
  timeline$start = timeline$start + firstOnset
  timeline$end = timeline$end + firstOnset
  
  # Maximum years
  maxDay = max(timeline$end)
  
  #Plot
  g <- ggplot(pat) + 
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
  
  return(g)
  
}


unionPatients <- function(patient1, patient2){
  patientID = paste(patient1$patientID[1], patient2$patientID[1], sep="-")
  
  patientU = rbind(patient1, patient2)
  patientU[,"patientID"] = patientID
  
  return(patientU)
}

#TO DO:
intersectPatients <- function(patient1, patient2){
  patientU = data.frame(patientID=character(), medication=character(), start=numeric(), end=numeric())
  
  sequences = getSequences(patient1, patient2)
  
  patientID = paste(patient1$patientID[1], patient2$patientID[1], sep="-")
  
  medications = intersect(patient1$medication, patient2$medication)
  
  for(medication in medications){
    
    med1 = sequences[[medication]][1,]
    med2 = sequences[[medication]][2,]
    size = length(med1)
    
    
    #Assign a value where both are the same
    indSameLetter = intersect(which(med1==med2), which(med1 != '∅'))
    if(length(indSameLetter)>1){
      
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
          combinedMed = data.frame("patientID"=patientID, "medication"=medication, "start"=start, "end"=end)
          patientU = rbind(patientU, combinedMed)
          
          #update start
          start = as.numeric(str_replace(colnames(combined)[ind], "day", ""))
        }
        current = ind
      }
      ind = letterSequence[length(letterSequence)]
      end = as.numeric(str_replace(colnames(combined)[ind], "day", ""))
      combinedMed = data.frame("patientID"=patientID, "medication"=medication, "start"=start, "end"=end)
      patientU = rbind(patientU, combinedMed)
      
    } else if(length(indSameLetter)==1){ #just one letter
      combinedMed = data.frame("patientID"=patientID, "medication"=medication, "start"=indSameLetter, "end"=indSameLetter)
      patientU = rbind(patientU, combinedMed)
    }
    
    
  }
  
  if(length(patientU)==0){
    patientU = data.frame(patientID=character(), medication=character(), start=numeric(), end=numeric())
    return(patientU)
  }
  colnames(patientU) = colnames(patient1)
  rownames(patientU) = NULL
  patientU = data.frame(patientU, stringsAsFactors=FALSE)
  patientU$medication = as.character(patientU$medication)
  patientU$start = as.numeric(as.character(patientU$start))
  patientU$end = as.numeric(as.character(patientU$end))
  patientU$patientID = as.character(patientU$patientID)
  
  
  return(patientU)
}


#TO DO:
averagePatients <- function(patient1, patient2){
  patientU = intersectPatients(patient1, patient2)
  #patientU = data.frame(patientID=character(), medication=character(), start=numeric(), end=numeric())
  
  
  patientID = paste(patient1$patientID[1], patient2$patientID[1], sep="-")
  
  #TO DO:
  medications = unique(c(as.character(patient1$medication), as.character(patient2$medication)))
  medications = medications[order(medications)]
  for(medication in medications){
    med.p1 = patient1[which(patient1$medication == medication),]
    med.p2 = patient2[which(patient2$medication == medication),]
    numRows.p1 = dim(med.p1)[1]
    numRows.p2 = dim(med.p2)[1]
    numRows.pU = ceiling((numRows.p1+numRows.p2)/2)
    
    # Create new combined rows
    for(i in 1:numRows.pU){
      if(i <= numRows.p1 && i <= numRows.p2){ #Both patients
        start = ceiling((med.p1[i,"start"]+med.p2[i,"start"])/2)
        end = ceiling((med.p1[i,"end"]+med.p2[i,"end"])/2)
      } else if(i <= numRows.p1 && i > numRows.p2){ #Only patient1
        dif = round((med.p1[i,"end"] - med.p1[i,"start"])/4)
        start = med.p1[i,"start"] + dif
        end = med.p1[i,"end"] - dif
      } else if(i > numRows.p1 && i <= numRows.p2){ #Only patient2
        dif = round((med.p2[i,"end"] - med.p2[i,"start"])/4)
        start = med.p2[i,"start"] + dif
        end = med.p2[i,"end"] - dif
      }
      
      combinedMed = data.frame("patientID"=patientID, "medication"=medication, "start"=start, "end"=end)
      patientU = rbind(patientU, combinedMed)
    }
    
    #print(paste(medication, numRows.p1, numRows.p2, numRows.pU, sep=", "))
  }
  
  patientU$medication = as.character(patientU$medication)
  patientU$patientID = as.character(patientU$patientID)
  
  
  return(patientU)
}


plotGI <- function(gi_score){
  
  if(typeof(gi_score) == "integer"){
    gi_score$id = paste("Patient", as.character(gi_score$id))
  }
  maxDay = max(gi_score$daysSinceFirstAppointment)
  
  gi_new <- ggplot(data = gi_score, aes(x = daysSinceFirstAppointment,  y = gi_new)) +
    geom_point(color="gray30", size=0.8) +
    geom_smooth(aes(group=id), color="gray50", method = "loess", size=0.5, se=FALSE, fill="gray82", alpha=0.5) +
    geom_smooth(method = "loess", size=2, se=TRUE, fill="gray82", alpha=0.5) +
    scale_y_continuous(limits = c(0, 100), expand=c(0,0)) +
    scale_x_continuous(limits = c(0, maxDay) , 
                       breaks=seq(0,(maxDay+365),365),
                       labels=paste("year", seq(0,(maxDay+365),365)/365),
                       expand = c(0, 0)) +
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      #legend.position="right", 
      legend.position="none", 
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
  
  return(gi_new)
}




