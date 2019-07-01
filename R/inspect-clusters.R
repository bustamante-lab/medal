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

# #Get all triangles
# gTri1 <- plotCoOcurrenceTriangle(1, events, profiles, names(medgroups), years)
# gTri2 <- plotCoOcurrenceTriangle(2, events, profiles, names(medgroups), years)
# gTri3 <- plotCoOcurrenceTriangle(3, events, profiles, names(medgroups), years)
# gTri4 <- plotCoOcurrenceTriangle(4, events, profiles, names(medgroups), years)
# 
# #Save plot
# gpanels <- ggarrange(gTri1, gTri3, gTri2, gTri4,
#                      labels = c("A", "B", "C", "D"),
#                      ncol = 4, nrow = 1, legend="none", common.legend = FALSE)
# ggexport(gpanels, filename="../images/Figure3-triangles.png", height = 2000, width = 6000, res=300)



#--------------------------------------------------
# Step 5. Plot the impairment scores
#--------------------------------------------------
color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")


  #Get all impairment scores
   gGIS1 <- plotScores(1, outcomes, profiles, 2, color.vector[1])
   gGIS2 <- plotScores(2, outcomes, profiles, 2, color.vector[2])
   gGIS3 <- plotScores(3, outcomes, profiles, 2, color.vector[3])
   gGIS4 <- plotScores(4, outcomes, profiles, 2, color.vector[4])
   
   
  # #Save plot
   gpanels <- ggarrange(gGIS1, gGIS3, gGIS2, gGIS4,
                        labels = c("A", "B", "C", "D"),
                        ncol = 4, nrow = 1, legend="none", common.legend = FALSE)
   ggexport(gpanels, filename="../images/Figure4-scores.png", height = 1500, width = 5000, res=300)
  

#--------------------------------------------------
# Step 6. Plot demographics
#--------------------------------------------------

   gAgeBox <- plotAgeBox(profiles, color.vector)
   gSex <- plotPercentage(profiles, 
                          var="is_male", 
                          var.label = "Sex",
                          var.coding = c("0"="Female", "1"="Male"),
                          color.vector=color.vector)
   gEthnic <- plotPercentage(profiles, 
                          var="NHW", 
                          var.label = "Ethic group",
                          var.coding = c("0"="Other", "1"="White (non Hisp)"),
                          xlim=c(0,1.2),
                          color.vector=color.vector)
   gOCD <- plotPercentage(profiles, 
                             var="OCD", 
                             var.label = "Obsesive Compulsive Disorder",
                             var.coding = c("0"="No", "1"="Yes"),
                             xlim=c(0,1.1),
                             color.vector=color.vector)
   gFood <- plotPercentage(profiles, 
                          var="foodprob", 
                          var.label = "Food problems",
                          var.coding = c("0"="No", "1"="Yes"),
                          color.vector=color.vector)
   gAnx <- plotPercentage(profiles, 
                           var="anx", 
                           var.label = "Anxiety",
                           var.coding = c("0"="No", "1"="Yes"),
                           color.vector=color.vector)
   gEmotional <- plotPercentage(profiles, 
                          var="emotional", 
                          var.label = "Emotional",
                          var.coding = c("0"="No", "1"="Yes"),
                          color.vector=color.vector)
   gMood <- plotPercentage(profiles, 
                                var="mood", 
                                var.label = "Mood",
                                var.coding = c("0"="No", "1"="Yes"),
                                color.vector=color.vector)
   gAgg <- plotPercentage(profiles, 
                           var="agg", 
                           var.label = "Agg",
                           var.coding = c("0"="No", "1"="Yes"),
                           color.vector=color.vector)
   gSch <- plotPercentage(profiles, 
                          var="sch", 
                          var.label = "School problems",
                          var.coding = c("0"="No", "1"="Yes"),
                          color.vector=color.vector)
   gReg <- plotPercentage(profiles, 
                             var="reg", 
                             var.label = "Regression",
                             var.coding = c("0"="No", "1"="Yes"),
                             color.vector=color.vector)
   gSleep <- plotPercentage(profiles, 
                          var="sleep", 
                          var.label = "Sleep disorders",
                          var.coding = c("0"="No", "1"="Yes"),
                          color.vector=color.vector)
   gTics <- plotPercentage(profiles, 
                            var="tics", 
                            var.label = "Tics",
                            var.coding = c("0"="No", "1"="Yes"),
                            color.vector=color.vector)
   
   
   gpanels <- ggarrange(gAgeBox, gSex, gEthnic, gOCD, gFood, gAnx, 
                        gEmotional, gMood, gAgg, gSch, gReg, gSleep, gTics,
                        ncol = 3, nrow = 5, legend="none", common.legend = FALSE)
   ggexport(gpanels, filename="../images/Figure5-profiles.png", height = 4000, width = 5000, res=300)
   
   
   