#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Mar 6, 2019                                           ##
##                                                             ##
#################################################################

remove(list=ls())

library(ggplot2)
library(reshape2)

source("support-functions.R")


#--------------------------------------------------
# Step 1. Read file
#--------------------------------------------------


#---
#File with patient ID (de-ID), comorbidities, initial clinical presentation, etc.
profiles = read.csv("../../clinical/data-matrix-profiles-cluster.csv", stringsAsFactors = FALSE)

#Variables for stratification
strats <- c("is_male", "NHW", "OCD", "foodprob", "anx", 
            "emotional", "mood", "agg", "sch", "reg", "sleep", "tics")

#Select only a few columns
profiles = profiles[,c("id", "age_onset", "age_1st_appt", strats, "cluster")]

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


#--------------------------------------------------
# Step 2. Clean Data
#--------------------------------------------------


# Group by class of medication
medgroups = list()
medgroups$penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin")
medgroups$cephalosporin = c("cephalexin", "cefadroxil")
medgroups$macrolide = c("azithromycin")
medgroups$nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin")
medgroups$hydrocortisone = c("prednisone", "maintenance prednisone", "decadron", "solumedrol")
medgroups$antibody = c("rituximab", "ivig")
medgroups$dmard = c("plaquenil", "methotrexate", "cellcept")


events = cleanEvents(events, medgroups)

#--------------------------------------------------
# Step 3. Plot cluster summaries 
#--------------------------------------------------

clusters = unique(sort(profiles$cluster))
years = 7
daysPerMonth = 30
n = years*12 #7 years
m = length(names(medgroups))

#Loop through the clusters
for(cluster in clusters){
  
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
  
  #melted_drug$year = floor(melted_drug$month/12)
  
  
  #Plot heatmap
  ggplot(melted_drug, aes(x=month, y=med, fill=value)) + 
    #redrawing tiles to remove cross lines from legend
    geom_tile(colour="white",size=0.25, show.legend = FALSE)+
    #remove extra space
    scale_y_discrete(expand=c(0,0), breaks=names(medgroups)) +
    scale_x_continuous(expand=c(0,0),
                     breaks=c(seq(1:years)*12),
                     labels=c(paste("year", seq(1:years))))+
    #custom colours for cut levels and na values
    scale_fill_gradient2(low = "white", mid = "#c994c7" , high = "#dd1c77",
                         midpoint = 0.7,
                         limit = c(0,1), space = "Lab", 
                         name="Medication Usage") +
    #equal aspect ratio x and y axis
    #coord_fixed() +
    #remove axis labels, add title
    labs(title="Medication usage")+
    #set base size for all font elements
    theme(#Add a title
      plot.title = element_text(hjust = 0.5, size=15),
      #Remove elements
      legend.position="right", 
      #legend.position="none", 
      #legend.title = element_blank(),
      axis.title.x=element_blank(),
      #axis.text.x=element_blank(),
      #axis.ticks.x=element_blank(),
      axis.text.x=element_text(angle=90, vjust=0.5),
      axis.title.y=element_blank(),
      axis.ticks.y=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank(),
      #Add a border
      panel.border = element_rect(colour = "black", fill=NA, size=1)
    )
  
  
  
  
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