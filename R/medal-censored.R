#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Mar 14, 2019                                          ##
##                                                             ##
#################################################################


remove(list=ls())

library(magrittr)
library(ggplot2)
library(ggpubr)
library(ggrepel)

library(dendextend)
library(factoextra)
library(NbClust)
library(cluster)
library(aricode)

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
medgroups = list()
medgroups$penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin")
medgroups$cephalosporin = c("cephalexin", "cefadroxil")
medgroups$macrolide = c("azithromycin")
medgroups$nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin")
medgroups$hydrocortisone = c("prednisone", "maintenance prednisone", "decadron", "solumedrol")
medgroups$antibody = c("rituximab", "ivig")
medgroups$dmard = c("plaquenil", "methotrexate", "cellcept")


data = cleanEvents(events, medgroups)
data1 = rightCensoring(data, 1)
data2 = rightCensoring(data, 2)
data3 = rightCensoring(data, 3)
data4 = rightCensoring(data, 4)
data5 = rightCensoring(data, 5)
data6 = rightCensoring(data, 6)
data7 = rightCensoring(data, 7)

data1 = twoTailCensoring(data, 0, 1)
data2 = twoTailCensoring(data, 1, 2)
data3 = twoTailCensoring(data, 2, 3)
data4 = twoTailCensoring(data, 3, 4)
data5 = twoTailCensoring(data, 4, 5)
data6 = twoTailCensoring(data, 5, 6)
data7 = twoTailCensoring(data, 6, 7)

write.csv(data1, "../../clinical/data-matrix-clean-1year.csv")
write.csv(data2, "../../clinical/data-matrix-clean-2year.csv")
write.csv(data3, "../../clinical/data-matrix-clean-3year.csv")
write.csv(data4, "../../clinical/data-matrix-clean-4year.csv")
write.csv(data5, "../../clinical/data-matrix-clean-5year.csv")
write.csv(data6, "../../clinical/data-matrix-clean-6year.csv")
write.csv(data7, "../../clinical/data-matrix-clean-7year.csv")


# Step 3. Create a distance matrix ----------------------


#-------
#Calling pyMEDAL


year1 <- pyMedal(file ="../../clinical/data-matrix-clean-1year.csv")
year2 <- pyMedal(file ="../../clinical/data-matrix-clean-2year.csv")
year3 <- pyMedal(file ="../../clinical/data-matrix-clean-3year.csv")
year4 <- pyMedal(file ="../../clinical/data-matrix-clean-4year.csv")
year5 <- pyMedal(file ="../../clinical/data-matrix-clean-5year.csv")
year6 <- pyMedal(file ="../../clinical/data-matrix-clean-6year.csv")
year7 <- pyMedal(file ="../../clinical/data-matrix-clean-7year.csv")

# Step 4.1 Choose number of clusters ----------------------

# #d = scale(distMatrix, center=FALSE)
# d = distMatrix
# k = 4
# 
# #removing outlier 37 (index 6)
# k = 3
# d = d[-6,]
# d = d[,-6]
# 
# # Elbow method
# elbow <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "wss") +
#   geom_vline(xintercept = k, linetype = 2) +
#   labs(title = "Elbow method")
# # Silhouette method
# silhouette <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "silhouette", 
#                            print.summary = FALSE) +
#   geom_vline(xintercept = k, linetype = 2) +
#   labs(title = "Silhouette method")
# # Gap statistic
# # nboot = 50 to keep the function speedy. 
# # recommended value: nboot= 500 for your analysis.
# # Use verbose = FALSE to hide computing progression.
# set.seed(123)
# gapStat <- fviz_nbclust(x=d, diss=as.dist(d), hcut, nstart = 25, 
#                         method = "gap_stat", nboot = 50, print.summary = FALSE,
#                         maxSE=list(method="Tibs2001SEmax", SE.factor=1)) +
#   geom_vline(xintercept = k, linetype = 2) +
#   labs(title = "Gap statistic method")
# 
# gpanels <- ggarrange(elbow, silhouette, gapStat,
#                      labels = c("A", "B", "C"),
#                      ncol = 1, nrow = 3, legend="bottom", 
#                      align="v", common.legend = FALSE)
# ggexport(gpanels, filename="../images/Figure1-num-clusters-year3.png", height = 3000, width = 2000, res=300)



# Step 4.1 Plot a dendrogram ----------------------

#k = 4 # on visual inspection of previous figure (with 37)
#k = 3 #without outlier37

#removing outlier 37 (index 6)
#k = 3
#e = d
#d = d[-6,]
#d = d[,-6]
#color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#0080ff")
#color.vector2 = c("3"="#0080ff", "2"="#ca3542", "1"="#e8a631")

#for PCAs and dendrogram
color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")
color.vector2 = c("3"="#00a572", "4"="#0080ff", "2"="#ca3542", "1"="#e8a631")

#-------


#-------
# Create MDS
k = 3

e = year1
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
gMDS1 <- plotMDS(e, k, color.vector, "Year 1")


e = year2
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
gMDS2 <- plotMDS(e, k, color.vector, "Year 2")


e = year3
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
gMDS3 <- plotMDS(e, k, color.vector, "Year 3")

e = year4
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
gMDS4 <- plotMDS(e, k, color.vector, "Year 4")


e = year5
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
gMDS5 <- plotMDS(e, k, color.vector, "Year 5")

e = year6
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
gMDS6 <- plotMDS(e, k, color.vector, "Year 6")
#-------
gpanels <- ggarrange(gMDS1, gMDS2, gMDS3,gMDS4, gMDS5,
                     ncol = 3, nrow = 2, legend="bottom", common.legend = TRUE)
ggexport(gpanels, filename="../images/Figure2-mds-years-ind.png", height = 3000, width = 4000, res=300)



#K-means clustering
#k2 <- kmeans(d, centers = k, nstart = 25)

#Calculating Normalized Mutual Information
#NMI(k2$cluster, hclust.assignment, variant="sum")
#https://course.ccs.neu.edu/cs6140sp15/7_locality_cluster/Assignment-6/NMI.pdf


#Saving the cluster to profiles
#write.csv(cbind(profiles, cluster=hclust.assignment),
#          "../../clinical/data-matrix-profiles-cluster.csv")
