#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Feb 28, 2019                                          ##
##                                                             ##
#################################################################

# library(ggplot2)
# library(stringr)
# library(reshape2)
# library(scales)
# library(dendextend)
# #library(MASS)
# library("dplyr")
# library("colorspace")
# library("tidyr")



# medcolors= c("penicillin"="#66c2a5",
#              "cephalosporin" = "#fc8d62",
#              "macrolide" = "#8da0cb",
#              "nsaid" = "#e7298a",
#              "corticosteroid.oral" = "#a6d854",
#              "corticosteroid.iv" = "#ffd92f",
#              "antibody" = "green",
#              "immunoglobulins" = "#e5c494",
#              "dmard" = "#b3b3b3")

daysPerMonth = 30


plotDendrogram <- function(dend, k, title="Dendrogram (hierarchical clustering)"){
  
  ggd1 <- as.ggdend(dend)
  gClust <- ggplot(ggd1, horiz = FALSE) +
    ggtitle(title) +
    labs(x="Patients", y="Hierarchical clustering height") +
    theme_light(base_size = 14) +
    theme(legend.position="bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
          #axis.title = element_blank()
    )
  
  return(gClust)
  
  
}

plotMDS <- function(d, cluster, color.vector=mycolors, dim1, dim2, title="MDS"){
  
  #Merge cluster and data
  score = MASS::isoMDS(as.matrix(d), k=4)
  score = as.data.frame(score$points)
  #score = as.data.frame(prcomp(d, scale. = FALSE)$x)
  score$cluster = cluster
  
  #Get axis names
  nms <- names(score)
  xname <- nms[dim1]
  yname <- nms[dim2]
  
  # plot of observations
  gMDS1 <- ggplot(data = score, aes(x = !!ensym(xname), y = !!ensym(yname))) +
    geom_text_repel(aes(label = rownames(score), color=cluster),
                    #color = "grey30",
                    min.segment.length = unit(0.5, 'lines'),
                    segment.color = 'grey90', show.legend = FALSE) +
    geom_point(aes(colour = cluster), size=2) +
    stat_chull(aes(colour = cluster, fill = cluster), alpha = 0.1, geom = "polygon") +
    #stat_ellipse(aes(colour = cluster, fill=cluster), geom="polygon", alpha=0.1) +
    #geom_text(aes(label = rownames(score)), color = "grey30", size=2) +
    scale_color_manual(values=color.vector) +
    scale_fill_manual(values=color.vector) +
    labs(x=paste("MDS", dim1,sep=""), y=paste("MDS", dim2,sep="")) +
    ggtitle(title) +
    theme_light(base_size = 14) +
    theme(legend.position="bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank()
    )
  
  return(gMDS1)
  
}

plotTSNE <- function(score, color.vector, title="MDS"){
  
  # plot of observations
  gMDS1 <- ggplot(data = score, aes(x = V1, y = V2)) +
    geom_text_repel(aes(label = rownames(score), color=cluster),
                    #color = "grey30",
                    min.segment.length = unit(0.5, 'lines'),
                    segment.color = 'grey90', show.legend = FALSE) +
    geom_point(aes(colour = cluster), size=2) +
    stat_chull(aes(colour = cluster, fill = cluster), alpha = 0.1, geom = "polygon") +
    #stat_ellipse(aes(colour = cluster, fill=cluster), geom="polygon", alpha=0.1) +
    #geom_text(aes(label = rownames(scores)), color = "grey30", size=2) +
    scale_color_manual(values=color.vector) +
    scale_fill_manual(values=color.vector) +
    labs(x="TSNE 1", y="TSNE 2") +
    ggtitle(title) +
    theme_light(base_size = 14) +
    theme(legend.position="bottom",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.ticks = element_blank(),
          axis.text = element_blank())
}


plotTimeSeriesDrug <- function(cluster, events, profiles, medcolors=mycolors, medgroups, years){
  
  n = years*12
  m = length(names(medgroups))
  
  #Get the events for patients in a cluster
  patIDs = profiles[which(profiles[,"cluster"] == cluster), "id"]
  #eveIDs = which(events[,"id"] %in% patIDs)
  #pat = events[eveIDs, ]
  pat <- events %>%
    filter(id %in% patIDs) %>%
    as.data.frame()
  
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
    if(length(medIDs)>0){
      clusterEvents = pat[medIDs,]
      clusterEvents = rightCensoringMonths(clusterEvents, years)
      
      if(dim(clusterEvents)[1]>0){
        
        for(i in 1:dim(clusterEvents)[1]){
          x = clusterEvents[i,"medication"]
          start = clusterEvents[i,"start"]
          end = clusterEvents[i,"end"]
          
          drug[x, start:end] = drug[x, start:end] + 1
        }
      }
    }
  }
  
  #Normalize by the number of patients in the cluster
  drug = drug/length(patIDs)
  
  
  #Reshape the matrix (Melt)
  melted_drug <- reshape2::melt(drug)
  colnames(melted_drug) = c("med", "month", "value")
  melted_drug[which(melted_drug$value>1), "value"]=1 #avoiding duplicates
  
  #Set all 0 values to -1
  melted_drug[which(melted_drug$value == 0),"value"] = -1
  
  # d <- data.frame(x=rep(1:20, 5), y=rnorm(100, 5, .2) + rep(1:5, each=20), z=rep(1:20, 5), grp=factor(rep(1:5, each=20)))
  # 
  # ggplot(d) +
  #   geom_path(aes(x, y, group=grp, alpha=z, color=grp), size=2)
  
  
  #Plot heatmap
  gPaths <- ggplot(melted_drug, aes(x=month, y=value)) +
    geom_path(aes(color = med, group=1), size=5, lineend = "round") +
    geom_area(aes(fill = med, group=1), color=NA, alpha=0.3) +
    #geom_smooth(method="loess", color="gray30", se=FALSE, size=1.2, linetype = "dashed") +
    #Group by medication
    facet_wrap(.~med, ncol=1, scales="free_y") +
    #Reshape the scales
    scale_y_continuous(limits=c(0,1),
                       breaks = c(0,0.5,1),
                       labels = rev(c("all patients", "some", "none"))) +
    scale_x_continuous(limits=c(0,years*12),
                       breaks=c(seq(0, years, 1)*12),
                       labels=c(paste("year", seq(0,years, 1))))+
    #custom colours
    scale_color_manual(name="Medications", values=medcolors,
                      labels=str_replace(names(medcolors), "\\.", " ")) +
    scale_fill_manual(name="Medications", values=medcolors,
                      labels=str_replace(names(medcolors), "\\.", " ")) +
    #Add labels
    #labs(title=paste("Cluster", cluster), x="Years of follow-up") +
    labs(title=paste("Cluster", cluster)) +
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
      #axis.text.x=element_text(angle=90, vjust=0.5),
      axis.text.x=element_blank(),
      axis.title.y=element_blank(),
      axis.title.x=element_blank(),
      axis.ticks.y=element_blank(),
      axis.ticks.x=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      #Facets
      strip.background = element_blank(),
      strip.text = element_text(size = 12, hjust=0)
    )
  
  return(gPaths)
}


plotCoOcurrenceTriangle <- function(cluster, events, profiles, medications, years){
  #Get the events for patients in a cluster
  patIDs = profiles[which(profiles[,"cluster"] == cluster), "id"]
  eveIDs = which(events[,"id"] %in% patIDs)
  pat = events[eveIDs, ]
  patients = unique(pat$id)
  
  #Convert all days into months
  pat$start = round(pat$start/daysPerMonth)
  pat$end = round(pat$end/daysPerMonth)
  pat = rightCensoringMonths(pat, years)
  
  #Create an empty matrix
  interactions = vector()
  
  #Fill in values
  for(col in 1:length(medications)){
    for(row in col:length(medications)){
      
      medsRow = medications[[row]]
      medsCol = medications[[col]]
      value = 0
      
      for(p in patients){
        indRow = which((pat$id %in% p) & (pat$medication %in% medsRow))
        indCol = which((pat$id %in% p) & (pat$medication %in% medsCol))
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
  
  
  medlabs= c("penicillin"="pe.",
             "cephalosporin" = "ce.",
             "macrolide" = "ma.",
             "nsaid" = "ns.",
             "hydrocortisone" = "hc.",
             "antibody" = "ab.",
             "dmard" = "dm.")
  
  
  #Plot heatmap
  gTri <- ggplot(interactions, aes(x=A, y=B)) +
    #geom_tile(aes(alpha=value), fill="steelblue", color="gray33") +
    geom_tile(aes(fill=value), color="gray33") +
    geom_text(aes(label=value)) +
    #scale_y_discrete(limits=rev(medications)) +
    scale_y_discrete(limits=rev(medications), labels=medlabs) +
    #scale_x_discrete(limits=medications) +
    scale_x_discrete(limits=medications, labels=medlabs) +
    scale_fill_gradient2(low="white", mid="gray", midpoint=0.5, high="orangered") +
    ggtitle(paste("Cluster", cluster)) +
    theme_light(base_size = 14) +
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
      panel.background = element_blank()
      #panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  return(gTri)
}


plotScores <- function(cluster, outcomes, score, ylims=c(0,100,5),
                       profiles, years, mycolor="gray50", label="Score",
                       title=TRUE, xlabel=TRUE){
  maxDay = years*365
  
  #Get the scores for patients in a cluster
  patIDs = profiles[which(profiles[,"cluster"] == cluster), "id"]
  #outIDs = which(outcomes[,"id"] %in% patIDs)
  #pat = outcomes[outIDs, ]
  pat <- outcomes %>%
    filter(id %in% patIDs) %>%
    as.data.frame()
  #patients = patIDs
  
  #pat$id <- as.character(pat$id)
  
  
  #Normalize all scores to being at the same day (onset)
  for(i in 1:dim(pat)[1]){
    
    onset = ceiling(profiles[which(profiles$id==pat[i, "id"]), "age_onset"]*365)
    pat[i, "daysSinceOnset"] = pat[i, "daysSinceBirth"] - onset
    
  }
  
  pat <- pat %>%
    filter(daysSinceOnset <= maxDay)
  
  
  #Calculate direction of trend
   trend <- lm(formula=get(score)~daysSinceOnset, data=pat)
   
   slope = round(trend$coefficients[2], digits=2)
   
   if(trend$coefficients[2] >= 0){
     trend.color <- "#fc8d59" #orange
  } else{
     trend.color <- "#91bfdb" #blue
  }
  
  #Plot heatmap
  gScore <- 
    ggplot(data = pat, aes(x = daysSinceOnset,  y = get(score))) +
    #geom_point(aes(group=id, color=id), size=1)+
    #geom_line(aes(group=id, color=id), size=0.5, linetype="dashed")+
    #geom_smooth(aes(group=id), color="gray50", method = "lm", size=1, se=FALSE) +
    #geom_line(aes(group=id), stat="smooth", method = "lm",
    #          color = "gray50", size = 1, linetype ="solid", alpha = 0.4) +
    #geom_smooth(method = "loess", size=0.5, se=FALSE, color="gray50") +
    geom_smooth(method = "lm", size=2, se=FALSE, color=trend.color,
                na.rm = TRUE, formula=y~x) +
    #annotate("text", x=500, y=30, label= slope, size=15) +
    scale_y_continuous(limits = c(ylims[1], ylims[2]), 
                       breaks=seq(0, ylims[2], (ylims[2]-ylims[1])/ylims[3]), 
                       expand=c(0,0)) +
    scale_x_continuous(limits = c(0, maxDay) , 
                       breaks=seq(0,(maxDay+365),365),
                       labels=paste("year", seq(0,(maxDay+365),365)/365),
                       expand = c(0, 0)) +
    labs(y=label, 
      x="Time since onset") +
    theme(#Add a title
      #plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      #legend.position="right", 
      legend.position="none", 
      legend.title = element_blank(),
      axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      #axis.ticks.x=element_blank(),
      axis.text.x=element_text(angle=90, vjust=0.5),
      axis.title.y=element_blank(),
      #axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      # Force the plot into a rectangle aspect ratio
      #aspect.ratio = 0.3,
      #Add a border
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  if(title==TRUE){
    gScore = gScore +
      ggtitle(label)
      #ggtitle(paste("Cluster", cluster, sep=" "))
  }
  if(xlabel==FALSE){
    gScore = gScore +
      theme(axis.text.x = element_blank(),
            axis.ticks.x=element_blank())
  }
  
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


plotAgeBox <- function(profiles, color.vector){
  
  gAgeBox <- ggplot(profiles) +
    geom_boxplot(aes(x=as.character(cluster), y=age_onset, fill=as.character(cluster)), 
                 alpha=0.8, color="gray50", outlier.shape = 1, show.legend = FALSE) +
    geom_hline(aes(yintercept = 7.5, linetype="dashed"), color = "palevioletred", size=1.5) +
    scale_fill_manual(values=color.vector) +
    scale_y_continuous(limits=c(2,13), breaks=seq(2,13,1)) +
    scale_linetype_manual(name="", values="dashed", labels="cohort mean") +
    labs(title="Age", x="Cluster assignment", y="Age of onset") +
    coord_flip() +
    theme_light() +
    theme(
      legend.position=c(0.12, 0.1),
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  return(gAgeBox)
  
}


plotPercentage <- function(profiles, var, var.label=var, var.coding, color.vector, xlim=c(0,1)){
  
  profs <- profiles %>%
    mutate(var = recode(get(var), !!!var.coding)) %>%
    group_by(cluster, var) %>%
    summarise (n = n()) %>%
    mutate(freq = n / sum(n)) %>%
    ungroup %>%
    complete(cluster, var, fill = list(freq = 0))

  
  newcolors= c(rbind(darken(color.vector, 0.3), lighten(color.vector, 0.3)))
  names(newcolors) = interaction(as.character(profs$cluster), profs$var)
  
  gVar <- 
    ggplot(profs, aes(x=as.character(cluster), y=freq)) +
    geom_bar(aes(fill=interaction(as.character(cluster), var)),
             stat="identity", position="dodge2", width = 0.8, color="gray50", alpha=0.8) +
    geom_text(aes(label=scales::percent(freq, accuracy = 5L)),
              position=position_dodge2(width = 0.8), vjust=-0.3) +
    geom_text(aes(label=var), 
              position=position_dodge2(width = 0.8), vjust=-1.8) +
    scale_fill_manual(name="cluster", values=newcolors) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits=xlim, breaks=seq(0,1,0.25)) +
    labs(title=var.label, x="Cluster assignment", y="Percentage") +
    theme_light() +
    theme(
      legend.position="none",
      axis.title = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank())
  
  return(gVar)
  
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




