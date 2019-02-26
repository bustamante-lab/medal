#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Feb 26, 2019                                          ##
##                                                             ##
#################################################################


remove(list=ls())

library(magrittr)
library(ggplot2)
library(ggpubr)
library(dendextend)
library(factoextra)
library(NbClust)
library(cluster)

source("medal-functions.R")
source("support-functions.R")
source("plot-functions.R")


# Step 1. Read file --------------------------------------------------



#---
#File with patient ID (de-ID), comorbidities, initial clinical presentation, etc.
profiles = read.csv("../../clinical/data-matrix-profiles.csv", stringsAsFactors = FALSE)

#Variables for stratification
strats <- c("is_male", "NHW", "OCD", "foodprob", "anx", 
            "emotional", "mood", "agg", "sch", "reg", "sleep", "tics")

#Select only a few columns
profiles = profiles[,c("id", "age_onset", "age_1st_appt", strats)]

#Get unique patient IDs ordered
patients = sort(unique(profiles$id))
write.csv(patients, "patients.txt", row.names = FALSE, col.names = FALSE)

#---
#File with events
events = read.csv("../../medication/meddrug-12-2018.csv", stringsAsFactors = FALSE)

#Only select rows for the same patients listed in Profiles
rows = which(events$id %in% patients)

#Select only a few columns
events = events[rows,c("id", "medication", "start", "end")]

#---
#File with clinical evaluations (outcomes)
outcomes = read.csv("../../clinical/data-matrix-outcomes.csv", stringsAsFactors = FALSE)

#Only select rows for the same patients listed in Profiles
rows = which(outcomes$id %in% patients)

#Select only a few columns
outcomes = outcomes[rows,c("id", "gi_new", "daysSinceBirth")]



# Step 2. Clean Data ------------------------------


# Group by class of medication
medgroups = vector()
medgroups$penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin")
medgroups$cephalosporin = c("cephalexin", "cefadroxil")
medgroups$macrolide = c("azithromycin")
medgroups$nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin")
medgroups$hydrocortisone = c("prednisone", "maintenance prednisone", "decadron", "solumedrol")
medgroups$antibody = c("rituximab", "ivig")
medgroups$dmard = c("plaquenil", "methotrexate", "cellcept")


data = cleanEvents(events, medgroups)

write.csv(data, "../../clinical/data-matrix-clean.csv")



# Step 3. Create a distance matrix ----------------------


#-------
#Calling MEDAL in R

# distMatrix = as.data.frame(matrix(rep(0, length(patients)*length(patients)), nrow = length(patients)), stringsAsFactors = FALSE)
#
# for(i in 2:length(patients)){
#   for(j in 1:(i-1)){
#     
#     p1=patients[i]
#     p2=patients[j]
#     
#     pat1=data[which(data$patientID==p1),]
#     pat2=data[which(data$patientID==p2),]
#     
#     distance =  medalDistance(pat1, pat2)
#     
#     print(paste("[",i,",",j,"] = ", distance, sep=""))
#     distMatrix[i,j] = distance
#     distMatrix[j,i] = distance
#   }
# }


#-------
#Calling pyMEDAL

system('python3 ../pymedal/pymedal.py ../../clinical/data-matrix-clean.csv', wait=TRUE)


distMatrix = read.table("distance_mat.txt")
patientIDs = read.table("patientID.txt", sep=",")

pID = as.vector(unlist(patientIDs))
colnames(distMatrix) = pID
rownames(distMatrix) = pID

# Step 4.1 Choose number of clusters ----------------------

#d = scale(distMatrix, center=FALSE)
d = distMatrix

# Elbow method
elbow <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "wss") +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(title = "Elbow method")
# Silhouette method
silhouette <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "silhouette", 
                           print.summary = FALSE) +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(title = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
gapStat <- fviz_nbclust(x=d, diss=as.dist(d), hcut, nstart = 25, 
                        method = "gap_stat", nboot = 50, print.summary = FALSE,
                        maxSE=list(method="Tibs2001SEmax", SE.factor=1)) +
  geom_vline(xintercept = 4, linetype = 2) +
  labs(title = "Gap statistic method")

gpanels <- ggarrange(elbow, silhouette, gapStat,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3, legend="bottom", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure1-num-clusters.png", height = 3000, width = 2000, res=300)



# Step 4.1 Plot a dendrogram ----------------------

k = 4 #visual inspection of previous figure

mybranch.colors

library(RColorBrewer)
n=8
color.vector = rep(brewer.pal(n, "Dark2"), ceiling(k/n))
color.vector = color.vector[1:k]

# Create Dendrogram
dend <- distMatrix %>% as.dist %>%
  hclust(method="ward.D") %>% as.dendrogram %>%
  set("branches_k_color", value = color.vector, k = k) %>% set("branches_lwd", 0.7) %>%
  set("labels_cex", 0.6) %>% set("labels_colors", value = color.vector, k = k) %>%
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5)
ggd1 <- as.ggdend(dend)
gClust <- ggplot(ggd1, horiz = FALSE)

# Create PCA
k2 <- kmeans(d, centers = k, nstart = 25)
gPCA <- fviz_cluster(k2, data = d)

gpanels <- ggarrange(gClust, gPCA,
                     labels = c("A", "B"),
                     ncol = 2, nrow = 1, legend="bottom", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure2-dendro-pca.png", height = 2000, width = 3000, res=300)





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