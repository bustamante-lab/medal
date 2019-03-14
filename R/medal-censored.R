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

write.csv(data, "../../clinical/data-matrix-clean-1year.csv")
write.csv(data, "../../clinical/data-matrix-clean-2year.csv")
write.csv(data, "../../clinical/data-matrix-clean-3year.csv")



# Step 3. Create a distance matrix ----------------------


#-------
#Calling pyMEDAL

system('python3 ../pymedal/pymedal.py ../../clinical/data-matrix-clean-1year.csv', wait=TRUE)


distMatrix = read.table("distance_mat.txt")
patientIDs = read.table("patientID.txt", sep=",")

pID = as.vector(unlist(patientIDs))
colnames(distMatrix) = pID
rownames(distMatrix) = pID

# Step 4.1 Choose number of clusters ----------------------

#d = scale(distMatrix, center=FALSE)
d = distMatrix
k=4

# Elbow method
elbow <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "wss") +
  geom_vline(xintercept = k, linetype = 2) +
  labs(title = "Elbow method")
# Silhouette method
silhouette <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "silhouette", 
                           print.summary = FALSE) +
  geom_vline(xintercept = k, linetype = 2) +
  labs(title = "Silhouette method")
# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
gapStat <- fviz_nbclust(x=d, diss=as.dist(d), hcut, nstart = 25, 
                        method = "gap_stat", nboot = 50, print.summary = FALSE,
                        maxSE=list(method="Tibs2001SEmax", SE.factor=1)) +
  geom_vline(xintercept = k, linetype = 2) +
  labs(title = "Gap statistic method")

gpanels <- ggarrange(elbow, silhouette, gapStat,
                     labels = c("A", "B", "C"),
                     ncol = 1, nrow = 3, legend="bottom", 
                     align="v", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure1-num-clusters-year1.png", height = 3000, width = 2000, res=300)



# Step 4.1 Plot a dendrogram ----------------------

k = 4 # on visual inspection of previous figure (with 37)
k=3 #without outlier37

#removing outlier 37 (index 6)
#k = 3
#e = d
#d = d[-6,]
#d = d[,-6]
color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#0080ff")
color.vector2 = c("3"="#0080ff", "2"="#ca3542", "1"="#e8a631")

#for PCAs and dendrogram
color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")
color.vector2 = c("3"="#00a572", "4"="#0080ff", "2"="#ca3542", "1"="#e8a631")

#-------
# Create Dendrogram
dend <- d %>% as.dist %>%
  hclust(method="ward.D") %>% as.dendrogram %>%
  set("branches_k_color", value = color.vector2, k = k) %>% set("branches_lwd", 0.7) %>%
  set("labels_cex", 0.6) %>% set("labels_colors", value = "grey30", k = k) %>%
  set("leaves_pch", 19) %>% set("leaves_cex", 0.5)
ggd1 <- as.ggdend(dend)
gClust <- ggplot(ggd1, horiz = FALSE)

#-------
# Create MDS
pca1 = prcomp(d, scale. = FALSE)
hclust.assignment = cutree(dend, k)
scores = as.data.frame(pca1$x)
scores = cbind(scores, cluster=as.character(hclust.assignment))

# plot of observations
gPCA12 <- ggplot(data = scores, aes(x = PC1, y = PC2)) +
  geom_text_repel(aes(label = rownames(scores)),
                  color = "grey30",
                  min.segment.length = unit(0.5, 'lines'),
                  segment.color = 'grey90') +
  geom_point(aes(colour = cluster), size=3) +
  stat_chull(aes(colour = cluster, fill = cluster), alpha = 0.1, geom = "polygon") +
  #stat_ellipse(aes(colour = cluster, fill=cluster), geom="polygon", alpha=0.1) +
  scale_color_manual(values=color.vector) +
  scale_fill_manual(values=color.vector) +
  labs(x="MDS1", y="MDS2") +
  ggtitle("MDS1 vs MDS2") +
  theme_light(base_size = 14) +
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())


gPCA23 <- ggplot(data = scores, aes(x = PC2, y = PC3)) +
  geom_text_repel(aes(label = rownames(scores)),
                  color = "grey30",
                  min.segment.length = unit(0.5, 'lines'),
                  segment.color = 'grey90') +
  geom_point(aes(colour = cluster), size=3) +
  stat_chull(aes(colour = cluster, fill = cluster), alpha = 0.1, geom = "polygon") +
  #stat_ellipse(aes(colour = cluster, fill=cluster), geom="polygon", alpha=0.1) +
  scale_color_manual(values=color.vector) +
  scale_fill_manual(values=color.vector) +
  labs(x="MDS2", y="MDS3") +
  ggtitle("MDS2 vs MDS3") +
  theme_light(base_size = 14) +
  theme(legend.position="bottom",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.ticks = element_blank(),
        axis.text = element_blank())

#-------
gpanels <- ggarrange(gClust, 
                     ggarrange(gPCA12, gPCA23,
                               ncol = 2, nrow=1, 
                               labels = c("B", "C"),
                               legend="bottom", common.legend = TRUE),
                     labels = c("A"),
                     ncol = 1, nrow = 2, legend="bottom", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure2-dendro-mds-rem37-year1.png", height = 3000, width = 4000, res=300)



#K-means clustering
k2 <- kmeans(d, centers = k, nstart = 25)

#Calculating Normalized Mutual Information
NMI(k2$cluster, hclust.assignment, variant="sum")
#https://course.ccs.neu.edu/cs6140sp15/7_locality_cluster/Assignment-6/NMI.pdf


#Saving the cluster to profiles
write.csv(cbind(profiles, cluster=hclust.assignment),
          "../../clinical/data-matrix-profiles-cluster.csv")
