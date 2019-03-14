###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)
# Date: Feb 26, 2019
#
###############################################################


rightCensoring <- function(data, year){

  days = year*365
  dataC = data
  
  #Remove all rows with start date beyond right censorship (days)
  indexes = which(dataC$start>days)
  dataC = dataC[-indexes, ]
  
  #Remove all rows with start date beyond right censorship (days)
  indexes = which(dataC$end>days)
  dataC[indexes, "end"] = days
  
  return(dataC)
}


cleanEvents <- function(data, groups){
  
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
  
  
  # Group by class of medication
  data$medication <- tolower(data$medication)

  # for(medclass in names(groups)){
  #   print(toupper(medclass))
  #   for(med in groups[[medclass]]){
  #     print(med)
  #   }
  #   print("-----")
  # }

  #Group together all medications of the same class
  for(medclass in names(groups)){
    indexes = which(data$medication %in% groups[[medclass]])
    data[indexes, "medication"] = medclass
  }
  
  return(data)
}