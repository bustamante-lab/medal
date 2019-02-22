###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)
# Date: Feb 21, 2019
#
###############################################################


cleanForConsistency <- function(data){
  
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
  
  
  # Merge Prednisone (burst and maintenance)
  data[which(data$medication == "Prednisone burst"), "medication"] = "Prednisone"
  data[which(data$medication == "Maintenance prednisone"), "medication"] = "Prednisone"
  
}