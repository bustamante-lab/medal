###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda <arturolp@stanford.edu>
# Date: June 28, 2018
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


events = read.csv("medication/medsEvents.csv", stringsAsFactors = FALSE)
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

#TODO: Update for all pairwise of patients
#p1="patient1"
#p2="patient2"
p1=1
p2=2
pat1=data[which(data$patientID==p1),]
pat2=data[which(data$patientID==p2),]

sequences = getSequences(pat1, pat2)


# Step 3. Medication Alignment Algorithm (Medal) ----------------------

medications = names(sequences)

for(medication in medications){
  
  if(dim(sequences[[medication]])[1] == 2){
    
    sequence1 = sequences[[medication]][1,]
    sequence2 = sequences[[medication]][2,]
    
    md = medalPairwise(sequence1, sequence2, TRUE)
    print(paste(medication, " [distance = ", md$distance, 
                ", size = ", length(md$alignment), "]",  sep =""))
  } 
}


