###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda <arturolp@stanford.edu>
# Date: Aug 3, 2018
#
###############################################################

remove(list=ls())

source("medal-functions.R")

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

# Step 4. Compute a dendrogram ----------------------


