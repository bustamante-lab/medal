#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Feb 27, 2019                                          ##
##                                                             ##
#################################################################


remove(list=ls())

#Plotting libraries
#library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

#Clustering libraries
library(factoextra)
library(NbClust)
#library(cluster)
library(aricode)
#library(tsne)
library(NMF) #for cluster purity and entropy

#Additional functions
#source("medal-functions.R")
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



# Step 2. Clean Data --------------------------------------------------


# Group by class of medication
medgroups = list()
medgroups$penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin")
medgroups$cephalosporin = c("cephalexin", "cefadroxil")
medgroups$macrolide = c("azithromycin")
medgroups$nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin")
medgroups$hydrocortisone = c("prednisone", "maintenance prednisone", "decadron", "solumedrol")
medgroups$antibody = c("rituximab", "ivig")
medgroups$dmard = c("plaquenil", "methotrexate", "cellcept")


data = cleanEvents(events, medgroups)

data = rightCensoring(data, 2)

write.csv(data, "../../clinical/data-matrix-clean.csv")



# Step 3. Create a distance matrix --------------------------------------------------


#-------
#Calling pyMEDAL

system('python3 ../pymedal/pymedal.py ../../clinical/data-matrix-clean.csv', wait=TRUE)


distMatrix = read.table("distance_mat.txt")
patientIDs = read.table("patientID.txt", sep=",")

pID = as.vector(unlist(patientIDs))
colnames(distMatrix) = pID
rownames(distMatrix) = pID

# Step 4 Choose number of clusters --------------------------------------------------

#d = scale(distMatrix, center=FALSE)
d = distMatrix

# Remove outlier (patient 37)
#e=d
#remInd = which(names(d)=="37")
#d = d[-remInd,]
#d = d[,-remInd]
k=4

# Elbow method
elbow <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "wss") +
  geom_vline(xintercept = k, linetype = "dashed", color="#5581B0", size=0.6) +
  labs(title = "Elbow method",
       y="Total within-clusters sum of squares")

# Silhouette method
silhouette <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "silhouette", 
                           print.summary = FALSE, barcolor = "white") +
  geom_vline(xintercept = k, linetype = "dashed", color="#5581B0", size=0.6) +
  labs(title = "Silhouette method")


# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
gapStat <- fviz_nbclust(x=d, diss=as.dist(d), hcut, nstart = 25, 
                        method = "gap_stat", nboot = 50, print.summary = FALSE,
                        maxSE=list(method="Tibs2001SEmax", SE.factor=1)) +
  #geom_vline(xintercept = 8, linetype = "dashed", color="#5581B0", size=0.6) +
  labs(title = "Gap statistic method")

gpanels <- ggarrange(elbow, silhouette, gapStat,
                     labels = c("A", "B", "C"),
                     ncol = 3, nrow = 1, legend="bottom", 
                     align="v", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure1-num-clusters.png", height = 1200, width = 4000, res=300)



# Step 5 Plots --------------------------------------------------------------

#Colors to be used
color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")
color.vector2 = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")
#color.vector3 = c("1"="mediumorchid3", "2"="darkturquoise", "3"="olivedrab3", "4"="orangered3")

#Plot Clustering strategies

dend = getDendrogram(d, k, color.vector)
gDend <- plotDendrogram(dend, k)
gMDSclus12 <- plotMDS(d, dend, color.vector, 1, 2, "MDS (hierarchical clustering)")
gMDSclus34 <- plotMDS(d, dend, color.vector, 3, 4, "MDS (hierarchical clustering)")

kmeans = getKMeansClusteringPCA(d, k)
#gMDSkmeans12 <- plotMDS(kmeans, color.vector, 1, 2, "MDS (k-means)")
#gMDSkmeans34 <- plotMDS(kmeans, color.vector, 3, 4, "MDS (k-means)")


# tsne1 = getHierarchicalClusteringTSNE(d, k)
# gTSNEclus <- plotTSNE(tsne1, color.vector, "TSNE (hierarchical clustering)")
# 
# tsne2 = getKMeansClusteringTSNE(d, k)
# gTSNEkmeans <- plotTSNE(tsne2, color.vector, "TSNE (k-means)")


# Combine plots and save
gpanels <- ggarrange(gDend, 
                     ggarrange(gMDSclus12, gMDSclus34,
                               labels = c("B", "C"),
                               align = "hv",
                               legend="bottom", common.legend = TRUE),
                     ncol = 1, nrow=2, 
                     labels = c("A"),
                     align = "h",
                     legend="bottom", common.legend = TRUE)
ggexport(gpanels, filename="../images/Figure2-dendro-mds.png", height = 4000, width = 4000, res=300)


#--------------------------------------------------

recoded = kmeans$cluster
recoded[which(recoded == 3)] = 1
recoded[which(recoded == 4)] = 2


#Calculating Normalized Mutual Information
NMI(recoded, cutree(dend,2), variant="sum")
NMI(kmeans$cluster, cutree(dend,4), variant="sum")
#https://course.ccs.neu.edu/cs6140sp15/7_locality_cluster/Assignment-6/NMI.pdf

#Calculating Cluster purity
purity(recoded, cutree(dend,2))
purity(kmeans$cluster, cutree(dend,4))
#https://www.rdocumentation.org/packages/NMF/versions/0.21.0/topics/purity

#Calculating Cluster entropy
entropy(recoded, cutree(dend,2))
entropy(kmeans$cluster, cutree(dend,4))
#https://www.rdocumentation.org/packages/NMF/versions/0.21.0/topics/purity


#Saving the cluster to profiles
assignment = cutree(dend,k)
index = which(profiles$id %in% names(assignment))
write.csv(cbind(profiles[index,], cluster=assignment),
          "../../clinical/data-matrix-profiles-cluster.csv")
