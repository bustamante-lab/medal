###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda <arturolp@stanford.edu>
# Date: Aug 10, 2018
#
###############################################################

remove(list=ls())

library(magrittr)
library(ggplot2)
library(dendextend)

source("medal-functions.R")
source("plot-functions.R")


# Step 1. Read file --------------------------------------------------

# data = as.data.frame(matrix(c("patient1","clindamycin", 1, 3,
#                               "patient1","clindamycin", 6, 8,
#                               "patient1","amoxicillin", 2, 3,
#                               "patient1","amoxicillin", 5, 6,
#                               "patient1","amoxicillin", 8, 9,
#                               "patient2","clindamycin", 3, 5,
#                               "patient2","clindamycin", 8, 10,
#                               "patient2","amoxicillin", 4, 5,
#                               "patient2","amoxicillin", 9, 10),
#                             nrow=9, ncol=4, byrow=TRUE))


events = read.csv("../medication/medsEvents.csv", stringsAsFactors = FALSE)
clinical = read.csv("../clinical/clinicHistoryDeIdentified.csv", stringsAsFactors = FALSE)

data=events[,c("id", "medication", "start", "end")]

colnames(data) = c("patientID", "medication", "start", "end")

# Step 2. Get patients and medications ------------------------------

patients = as.vector(sort(unique(data$patientID)))


# Clean data =======
#For single-day medication without end date
indexes = intersect(which(is.na(data$end)), which(!is.na(data$start)))
for(ind in indexes){
  data[ind, "end"] = data[ind, "start"]+1
}
#For single-day medication without start date
indexes = intersect(which(!is.na(data$end)), which(is.na(data$start)))
for(ind in indexes){
  data[ind, "start"] = data[ind, "end"]-1
}

#For negative start/end day
indexes = which(as.numeric(as.character(data$start)) < 1)
if(length(indexes) > 0 ){
  data[indexes, "start"] = 1
}
indexes = which(as.numeric(as.character(data$end)) < 1)
if(length(indexes) > 0 ){
  data[indexes, "end"] = 1
}

#Remove NA/NA
indexes = intersect(which(!is.na(data$start)), which(!is.na(data$end)))
data = data[indexes,]

# Merge Prednisone (burst and maintenance)
data[which(data$medication == "Prednisone burst"), "medication"] = "Prednisone"
data[which(data$medication == "Maintenance prednisone"), "medication"] = "Prednisone"



# Step 3. Create a distance matrix ----------------------

distMatrix = as.data.frame(matrix(rep(0, length(patients)*length(patients)), nrow = length(patients)), stringsAsFactors = FALSE)

for(i in 2:length(patients)){
  for(j in 1:(i-1)){
    
    p1=patients[i]
    p2=patients[j]
    
    pat1=data[which(data$patientID==p1),]
    pat2=data[which(data$patientID==p2),]
    
    distance =  medalDistance(pat1, pat2)
    
    print(paste("[",i,",",j,"] = ", distance, sep=""))
    distMatrix[i,j] = distance
    distMatrix[j,i] = distance
  }
}

#write.csv(distMatrix, "../distance-matrix-medal.csv")
distMatrix = read.csv("../distance-matrix-medal.csv", row.names = 1)

# Step 4. Compute a dendrogram ----------------------

k = 12



mybranch.colors

library(RColorBrewer)
n=8
color.vector = rep(brewer.pal(n, "Dark2"), ceiling(k/n))
color.vector = color.vector[1:k]


dend <- distMatrix %>% as.dist %>%
  hclust(method="complete") %>% as.dendrogram %>%
  set("branches_k_color", value = color.vector, k = k) %>% set("branches_lwd", 0.7) %>%
  set("labels_cex", 0.6) %>% set("labels_colors", value = color.vector, k = k) %>%
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5)
  #collapse_branch(tol = 11)  %>% hang.dendrogram(hang = 0)
ggd1 <- as.ggdend(dend)
ggplot(ggd1, horiz = FALSE) 

clusterCut = cutree(dend, k)

order = order.dendrogram(dend)



# Step 5. Plot summary patients ----------------------


#TO DO:
#Loop through the dendrogram inorder (using the tree structure)
#for(i in 1:length(dend)){
#  attr(dend[[i]], "members")
#  print(length(dend[[i]]))
#}

#Loop through the clusters
for(i in 1:k){
  patient.order = names(clusterCut[which(clusterCut==i)])
  patient.order = order[which(order %in% patient.order)]
  clusterName = paste("Cluster",i,"patients", paste(patient.order,collapse="-"), sep="-")
  print(clusterName)
  
  #For Patient 1
  patient1.ID = patient.order[1]
  patient1 = data[which(data$patientID==patient1.ID),]
  cli1 = clinical[which(clinical[, "id"] == patient1.ID), ]
  firstAppointment.p1 = cli1$daysSinceBirth[1] - cli1$daysSinceFirstAppointment[1]
  firstOnset.p1 = cli1$daysSinceBirth[1] - cli1$daysPostOnset[1]
  
  for(j in 1:(length(patient.order)-1)){
    
    #For Patient 2
    patient2.ID = patient.order[j+1]
    patient2 = data[which(data$patientID==patient2.ID),]
    cli2 = clinical[which(clinical[, "id"] == patient2.ID), ]
    firstAppointment.p2 = cli2$daysSinceBirth[1] - cli2$daysSinceFirstAppointment[1]
    firstOnset.p2 = cli2$daysSinceBirth[1] - cli2$daysPostOnset[1]
    
    #Update Values
    firstAppointment = floor((firstAppointment.p1 + firstAppointment.p2) /2)
    firstOnset = floor((firstOnset.p1 + firstOnset.p2) / 2)
    
    #calculate a composite timeline
    #timeline = intersectPatients(patient1, patient2)
    #timeline = unionPatients(patient1, patient2)
    timeline = averagePatients(patient1, patient2)
    
    #update patient1
    patient1 = timeline
    firstAppointment.p1 = firstAppointment
    firstOnset.p1 = firstOnset
    if(dim(patient1)[1] ==0){
      break()
    }
  }
  
  if(dim(patient1)[1] > 0){
    #Get Label
    if(nchar(patient1$patientID[1])>30){
      patient.label = paste("Cluster ", i, sep="")
      clusterName = paste("Cluster", i, sep="")
    } else{
      patient.label = paste("Patient ", timeline$patientID[1], sep="")
    }
    
    #Plot
    
    g <- plotPatientTimeline(patient1, patient.label, firstAppointment, firstOnset)
    #patient.timeline = paste("../images/cluster-union/",clusterName,".png", sep="")
    #patient.timeline = paste("../images/cluster-intersect/",clusterName,".png", sep="")
    patient.timeline = paste("../images/cluster-average/", clusterName,".png", sep="")
    ggsave(patient.timeline, width = 8, height = 6, dpi=300)
  }
}


# Step 6. Plot correlation triangles ----------------------

#oral anti-inflammatorics
A <- c("NSAID", "Prednisone")

#IV therapies
B <- c("IVIG", "Rituximab", "Solumedrol")

#Antibiotics Amoxicillin
C <- c("Augmentin", "Amoxicillin")

#Antibiotics Ceph+Lide
D <- c("Cefadroxil", "Clindamycin", "Azithromycin", "Cephalexin")

medications = list(A=A,B=B,C=C,D=D)

for(clust in 1:k){
  patients = names(clusterCut[which(clusterCut==clust)])
  patients = order[which(order %in% patients)]
  clusterName = paste("Cluster",clust,"patients", paste(patients,collapse="-"), sep="-")
  print(clusterName)
  
  #Get data for the entire cluster
  cluster = data[which(data$patientID %in% patients),]
  
  #Create an empty matrix
  interactions = as.data.frame(matrix("",4,4), stringsAsFactors = FALSE)
  rownames(interactions) = names(medications)
  colnames(interactions) = names(medications)
  
  #Fill in values
  for(col in 1:length(medications)){
    for(row in col:length(medications)){
      medsRow = medications[[row]]
      medsCol = medications[[col]]
      pats = 0
      for(p in patients){
        indRow = which((data$patientID %in% p) & (data$medication %in% medsRow))
        indCol = which((data$patientID %in% p) & (data$medication %in% medsCol))
        if(length(indRow)>0 & length(indCol)>0){
          #pats = pats + (length(indRow) + length(indCol))/2
          pats = pats + 1
        }
      }
      interactions[row,col]=round(pats/length(patients), digits=1)
    }
  }
  
  print(interactions)
}


# Step 6. Plot the impairment scores ----------------------

for(clust in 1:k){
 
  patients = names(clusterCut[which(clusterCut==clust)])
  patients = order[which(order %in% patients)]
  clusterName = paste("Cluster",clust, sep="-")
  print(clusterName)
  
  #Get data for the cluster: Global Impairment Score
  gi_score = clinical[which(clinical$id %in% patients),c("id", "gi_new", "daysPostOnset", "daysSinceBirth", "daysSinceFirstAppointment")]

  g <- plotGI(gi_score)
  patient.timeline = paste("../images/cluster-gi-score/", clusterName,".png", sep="")
  ggsave(patient.timeline, width = 8, height = 8, dpi=300)
}