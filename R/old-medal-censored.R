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
data1 = rightCensoringExclusive(data, 1)
data2 = rightCensoringExclusive(data, 2)
data3 = rightCensoringExclusive(data, 3)
data4 = rightCensoringExclusive(data, 4)
data5 = rightCensoringExclusive(data, 5)
data6 = rightCensoringExclusive(data, 6)
data7 = rightCensoringExclusive(data, 7)

#data1 = twoTailCensoring(data, 0, 1)
#data2 = twoTailCensoring(data, 0, 2)
#data3 = twoTailCensoring(data, 0, 3)
#data4 = twoTailCensoring(data, 0, 4)
#data5 = twoTailCensoring(data, 0, 5)
#data6 = twoTailCensoring(data, 0, 6)
#data7 = twoTailCensoring(data, 0, 7)

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
color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#0080ff")
#color.vector2 = c("3"="#0080ff", "2"="#ca3542", "1"="#e8a631")

#for PCAs and dendrogram
#color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")
#color.vector2 = c("3"="#00a572", "4"="#0080ff", "2"="#ca3542", "1"="#e8a631")

color.vector3 = c("2"="#e8a631", "1"="#ca3542", "3"="#0080ff")


#-------


#-------
# Create MDS
k = 3

e = year1
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
clus = getHierarchicalClustering(e)
gMDS <- plotMDS(clus, color.vector, "Years 0-1")
prof = as.data.frame(cbind(id=rownames(clus), cluster=clus$cluster))
gTri1 <- plotCoOcurrenceTriangle(1, data1, prof, names(medgroups))
gTri2 <- plotCoOcurrenceTriangle(2, data1, prof, names(medgroups))
gTri3 <- plotCoOcurrenceTriangle(3, data1, prof, names(medgroups))
gpanels1 <- ggarrange(gMDS, gTri1, gTri2, gTri3,
                     ncol = 1, nrow = 4, legend="bottom", common.legend = TRUE)


e = year2
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
clus = getHierarchicalClustering(e)
gMDS <- plotMDS(clus, color.vector, "Years 0-2")
prof = as.data.frame(cbind(id=rownames(clus), cluster=clus$cluster))
gTri1 <- plotCoOcurrenceTriangle(1, data2, prof, names(medgroups))
gTri2 <- plotCoOcurrenceTriangle(2, data2, prof, names(medgroups))
gTri3 <- plotCoOcurrenceTriangle(3, data2, prof, names(medgroups))
gpanels2 <- ggarrange(gMDS, gTri1, gTri2, gTri3,
                      ncol = 1, nrow = 4, legend="bottom", common.legend = TRUE)

e = year3
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
clus = getHierarchicalClustering(e)
gMDS <- plotMDS(clus, color.vector, "Years 0-3")
prof = as.data.frame(cbind(id=rownames(clus), cluster=clus$cluster))
gTri1 <- plotCoOcurrenceTriangle(1, data3, prof, names(medgroups))
gTri2 <- plotCoOcurrenceTriangle(2, data3, prof, names(medgroups))
gTri3 <- plotCoOcurrenceTriangle(3, data3, prof, names(medgroups))
gpanels3 <- ggarrange(gMDS, gTri1, gTri2, gTri3,
                      ncol = 1, nrow = 4, legend="bottom", common.legend = TRUE)

e = year4
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
clus = getHierarchicalClustering(e)
gMDS <- plotMDS(clus, color.vector, "Years 0-4")
prof = as.data.frame(cbind(id=rownames(clus), cluster=clus$cluster))
gTri1 <- plotCoOcurrenceTriangle(1, data4, prof, names(medgroups))
gTri2 <- plotCoOcurrenceTriangle(2, data4, prof, names(medgroups))
gTri3 <- plotCoOcurrenceTriangle(3, data4, prof, names(medgroups))
gpanels4 <- ggarrange(gMDS, gTri1, gTri2, gTri3,
                      ncol = 1, nrow = 4, legend="bottom", common.legend = TRUE)


e = year5
ind = which(names(e)=="37")
e = e[-ind,]
e = e[,-ind]
clus = getHierarchicalClustering(e)
gMDS <- plotMDS(clus, color.vector, "Years 0-5")
prof = as.data.frame(cbind(id=rownames(clus), cluster=clus$cluster))
gTri1 <- plotCoOcurrenceTriangle(1, data5, prof, names(medgroups))
gTri2 <- plotCoOcurrenceTriangle(2, data5, prof, names(medgroups))
gTri3 <- plotCoOcurrenceTriangle(3, data5, prof, names(medgroups))
gpanels5 <- ggarrange(gMDS, gTri1, gTri2, gTri3,
                      ncol = 1, nrow = 4, legend="bottom", common.legend = TRUE)

#-------
gpanels <- ggarrange(gpanels1, gpanels2, gpanels3, gpanels4, gpanels5, 
                     ncol = 5, nrow = 1, legend="bottom", common.legend = TRUE)
#ggexport(gpanels, filename="../images/Figure2-mds-years-two.png", height = 1200, width = 5000, res=300)
#ggexport(gpanels, filename="../images/Figure2-mds-years-two.png", height = 3000, width = 4000, res=300)
ggexport(gpanels, filename="../images/Figure2-mds-years-exclusive.png", height = 3000, width = 4000, res=300)



#-------
# Inspect medication usage
#-------

medcolors= c("penicillin"="#66c2a5",
             "cephalosporin" = "#fc8d62",
             "macrolide" = "#8da0cb",
             "nsaid" = "#e7298a",
             "hydrocortisone" = "#a6d854",
             "antibody" = "#ffd92f",
             "dmard" = "#e5c494")


#Get all triangles
gTri1 <- plotCoOcurrenceTriangle(1, events, profiles, names(medgroups))
gTri2 <- plotCoOcurrenceTriangle(2, events, profiles, names(medgroups))
gTri3 <- plotCoOcurrenceTriangle(3, events, profiles, names(medgroups))

#Save plot
gpanels <- ggarrange(gTri1, gTri2, gTri3, 
                     labels = c("", "", ""),
                     ncol = 3, nrow = 1, legend="bottom", common.legend = TRUE)
ggexport(gpanels, filename="../images/Figure3-triangles.png", height = 2000, width = 5000, res=300)

