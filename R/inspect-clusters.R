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
scores = c("gi_new", "cbiTotal", "ocd_score", "cybocs_total", 
           "aggresion_score" ,"cis_score", 
           "wpi_calculated", "sss_calculated",
           "arfi_total", "aggresion_score", "ygtss")
outcomes = outcomes[rows,c("id", scores,  "daysSinceBirth")]


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
   gGIS1 <- plotScores(1, outcomes, "gi_new", c(0,100,5), profiles, 2, color.vector[1], "Global Impairment", TRUE, FALSE)
   gGIS2 <- plotScores(2, outcomes, "gi_new", c(0,100,5), profiles, 2, color.vector[2], "Global Impairment", TRUE, FALSE)
   gGIS3 <- plotScores(3, outcomes, "gi_new", c(0,100,5), profiles, 2, color.vector[3], "Global Impairment", TRUE, FALSE)
   gGIS4 <- plotScores(4, outcomes, "gi_new", c(0,100,5), profiles, 2, color.vector[4], "Global Impairment", TRUE, FALSE)
   
   gcbi1 <- plotScores(1, outcomes, "cbiTotal", c(0,96,4), profiles, 2, color.vector[1], "Caregiver Burden", FALSE, TRUE)
   gcbi2 <- plotScores(2, outcomes, "cbiTotal", c(0,96,4), profiles, 2, color.vector[2], "Caregiver Burden", FALSE, TRUE)
   gcbi3 <- plotScores(3, outcomes, "cbiTotal", c(0,96,4), profiles, 2, color.vector[3], "Caregiver Burden", FALSE, TRUE)
   gcbi4 <- plotScores(4, outcomes, "cbiTotal", c(0,96,4), profiles, 2, color.vector[4], "Caregiver Burden", FALSE, TRUE)
   
   
   #Other scores
   gocd1 <- plotScores(1, outcomes, "ocd_score", c(0,15,3), profiles, 2, color.vector[1], "OCD Score", TRUE, FALSE)
   gocd2 <- plotScores(2, outcomes, "ocd_score", c(0,15,3), profiles, 2, color.vector[2], "OCD Score", TRUE, FALSE)
   gocd3 <- plotScores(3, outcomes, "ocd_score", c(0,15,3), profiles, 2, color.vector[3], "OCD Score", TRUE, FALSE)
   gocd4 <- plotScores(4, outcomes, "ocd_score", c(0,15,3), profiles, 2, color.vector[4], "OCD Score", TRUE, FALSE)
   
   gcyb1 <- plotScores(1, outcomes, "cybocs_total", c(0,40,4), profiles, 2, color.vector[1], "CYBOCS", FALSE, FALSE)
   gcyb2 <- plotScores(2, outcomes, "cybocs_total", c(0,40,4), profiles, 2, color.vector[2], "CYBOCS", FALSE, FALSE)
   gcyb3 <- plotScores(3, outcomes, "cybocs_total", c(0,40,4), profiles, 2, color.vector[3], "CYBOCS", FALSE, FALSE)
   gcyb4 <- plotScores(4, outcomes, "cybocs_total", c(0,40,4), profiles, 2, color.vector[4], "CYBOCS", FALSE, FALSE)
   
   gAgg1 <- plotScores(1, outcomes, "aggresion_score", c(0,100,5), profiles, 2, color.vector[1], "Aggression", FALSE, FALSE)
   gAgg2 <- plotScores(2, outcomes, "aggresion_score", c(0,100,5), profiles, 2, color.vector[2], "Aggression", FALSE, FALSE)
   gAgg3 <- plotScores(3, outcomes, "aggresion_score", c(0,100,5), profiles, 2, color.vector[3], "Aggression", FALSE, FALSE)
   gAgg4 <- plotScores(4, outcomes, "aggresion_score", c(0,100,5), profiles, 2, color.vector[4], "Aggression", FALSE, FALSE)
   
   gcis1 <- plotScores(1, outcomes, "cis_score", c(0,42,3), profiles, 2, color.vector[1], "Columbia Impairment", FALSE, TRUE)
   gcis2 <- plotScores(2, outcomes, "cis_score", c(0,42,3), profiles, 2, color.vector[2], "Columbia Impairment", FALSE, TRUE)
   gcis3 <- plotScores(3, outcomes, "cis_score", c(0,42,3), profiles, 2, color.vector[3], "Columbia Impairment", FALSE, TRUE)
   gcis4 <- plotScores(4, outcomes, "cis_score", c(0,42,3), profiles, 2, color.vector[4], "Columbia Impairment", FALSE, TRUE)
   
   garfi1 <- plotScores(1, outcomes, "arfi_total", c(0,15,3), profiles, 2, color.vector[1], "Avoidant Restrictive Food Intake", TRUE, FALSE)
   garfi2 <- plotScores(2, outcomes, "arfi_total", c(0,15,3), profiles, 2, color.vector[2], "Avoidant Restrictive Food Intake", TRUE, FALSE)
   garfi3 <- plotScores(3, outcomes, "arfi_total", c(0,15,3), profiles, 2, color.vector[3], "Avoidant Restrictive Food Intake", TRUE, FALSE)
   garfi4 <- plotScores(4, outcomes, "arfi_total", c(0,15,3), profiles, 2, color.vector[4], "Avoidant Restrictive Food Intake", TRUE, FALSE)
   
   gyg1 <- plotScores(1, outcomes, "ygtss", c(0,100,5), profiles, 2, color.vector[1], "Yale Global Tic Severity", FALSE, FALSE)
   gyg2 <- plotScores(2, outcomes, "ygtss", c(0,100,5), profiles, 2, color.vector[2], "Yale Global Tic Severity", FALSE, FALSE)
   gyg3 <- plotScores(3, outcomes, "ygtss", c(0,100,5), profiles, 2, color.vector[3], "Yale Global Tic Severity", FALSE, FALSE)
   gyg4 <- plotScores(4, outcomes, "ygtss", c(0,100,5), profiles, 2, color.vector[4], "Yale Global Tic Severity", FALSE, FALSE)
  
   gwpi1 <- plotScores(1, outcomes, "wpi_calculated", c(0,18,3), profiles, 2, color.vector[1], "Widespread Pain", FALSE, FALSE)
   gwpi2 <- plotScores(2, outcomes, "wpi_calculated", c(0,18,3), profiles, 2, color.vector[2], "Widespread Pain", FALSE, FALSE)
   gwpi3 <- plotScores(3, outcomes, "wpi_calculated", c(0,18,3), profiles, 2, color.vector[3], "Widespread Pain", FALSE, FALSE)
   gwpi4 <- plotScores(4, outcomes, "wpi_calculated", c(0,18,3), profiles, 2, color.vector[4], "Widespread Pain", FALSE, FALSE)
   
   gsss1 <- plotScores(1, outcomes, "sss_calculated", c(0,12,3), profiles, 2, color.vector[1], "Symptom Severity", FALSE, TRUE)
   gsss2 <- plotScores(2, outcomes, "sss_calculated", c(0,12,3), profiles, 2, color.vector[2], "Symptom Severity", FALSE, TRUE)
   gsss3 <- plotScores(3, outcomes, "sss_calculated", c(0,12,3), profiles, 2, color.vector[3], "Symptom Severity", FALSE, TRUE)
   gsss4 <- plotScores(4, outcomes, "sss_calculated", c(0,12,3), profiles, 2, color.vector[4], "Symptom Severity", FALSE, TRUE)
   
   
  # #Save plot
   gpanels <- ggarrange(gGIS1, gGIS3, gGIS2, gGIS4, 
                        gcbi1, gcbi3, gcbi2, gcbi4,
                        ncol = 4, nrow = 2, legend="none", common.legend = FALSE)
   ggexport(gpanels, filename="../images/Figure4-scores.png", height = 3000, width = 5000, res=300)
  
   
   gpanels <- ggarrange(gocd1, gocd3, gocd2, gocd4,
                        gcyb1, gcyb3, gcyb2, gcyb4,
                        gAgg1, gAgg3, gAgg2, gAgg4,
                        gcis1, gcis3, gcis2, gcis4,
                        ncol = 4, nrow = 4, legend="none", common.legend = FALSE)
   ggexport(gpanels, filename="../images/Figure4-other-scores1.png", height = 3000, width = 5000, res=300)
   
   
   gpanels <- ggarrange(garfi1, garfi3, garfi2, garfi4,
                        gyg1, gyg3, gyg2, gyg4,
                        gwpi1, gwpi3, gwpi2, gwpi4,
                        gsss1, gsss3, gsss2, gsss4,
                        ncol = 4, nrow = 4, legend="none", common.legend = FALSE)
   ggexport(gpanels, filename="../images/Figure4-other-scores2.png", height = 3000, width = 5000, res=300)
   
   

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
   
   
   