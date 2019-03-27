#################################################################
##                                                             ##
## Project: Medication Alignment Algorithm (Medal)             ##  
## Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)  ##
## Date: Mar 6, 2019                                           ##
##                                                             ##
#################################################################

remove(list=ls())

library(ggplot2)
library(ggpubr)

source("support-functions.R")
source("plot-functions.R")


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

medcolors= c("penicillin"="#66c2a5",
             "cephalosporin" = "#fc8d62",
             "macrolide" = "#8da0cb",
             "nsaid" = "#e7298a",
             "hydrocortisone" = "#a6d854",
             "antibody" = "#ffd92f",
             "dmard" = "#e5c494")



clusters = unique(sort(profiles$cluster))
years = 2

n = years*12 #7 years
m = length(names(medgroups))

#Removing patient 37 from the analysis
#clusters = clusters[-4]
#profiles[which(profiles[,"cluster"] == 3), "cluster"] = 5
#profiles[which(profiles[,"cluster"] == 4), "cluster"] = 3


#Get all paths
gPath1 <- plotTimeSeriesDrug(1, events, profiles, medcolors, medgroups, years)
gPath2 <- plotTimeSeriesDrug(2, events, profiles, medcolors, medgroups, years)
gPath3 <- plotTimeSeriesDrug(3, events, profiles, medcolors, medgroups, years)
gPath4 <- plotTimeSeriesDrug(4, events, profiles, medcolors, medgroups, years)


#Save plot
gpanels <- ggarrange(gPath1, gPath3, gPath2, gPath4,
                     labels = c("A", "B", "C", "D"),
                     ncol = 4, nrow = 1, legend="none", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure3-clusters.png", height = 4000, width = 5000, res=300)

  


#--------------------------------------------------
# Step 4. Plot co-occurence triangles 
#--------------------------------------------------

#Get all triangles
gTri1 <- plotCoOcurrenceTriangle(1, events, profiles, names(medgroups), years)
gTri2 <- plotCoOcurrenceTriangle(2, events, profiles, names(medgroups), years)
gTri3 <- plotCoOcurrenceTriangle(3, events, profiles, names(medgroups), years)
gTri4 <- plotCoOcurrenceTriangle(4, events, profiles, names(medgroups), years)

#Save plot
gpanels <- ggarrange(gTri1, gTri3, gTri2, gTri4,
                     labels = c("A", "B", "C", "D"),
                     ncol = 4, nrow = 1, legend="none", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure3-triangles.png", height = 2000, width = 6000, res=300)



#--------------------------------------------------
# Step 5. Plot the impairment scores
#--------------------------------------------------


  # #Get all impairment scores
  # gGIS1 <- plotScores(1, outcomes, profiles)
  # gGIS2 <- plotScores(2, outcomes, profiles)
  # gGIS3 <- plotScores(3, outcomes, profiles)
  # 
  # 
  # #Save plot
  # gpanels <- ggarrange(gGIS1, gGIS2, gGIS3, 
  #                      labels = c("", "", ""),
  #                      ncol = 3, nrow = 1, legend="none", common.legend = FALSE)
  # ggexport(gpanels, filename="../images/Figure3-scores.png", height = 2000, width = 5000, res=300)
  # 