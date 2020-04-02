###############################################################
#
# Project: Medication Alignment Algorithm (Medal)
# Author: Arturo Lopez Pineda (arturolp[at]stanford[dot]edu)
# Date: July 28, 2018
#
###############################################################
library(compiler)
matrixInitialization <- function(sequence1, sequence2){
  
  n = length(sequence1)
  m = length(sequence2)
  
  edit.matrix = matrix(rep(0, ((n+1)*(m+1))), ncol=n+1, nrow=m+1)
  arrow.matrix = matrix(rep("", ((n+1)*(m+1))), ncol=n+1, nrow=m+1)
  colnames(edit.matrix) = c("", sequence1)
  rownames(edit.matrix) = c("", sequence2)
  
  return(list(edit=edit.matrix, arrow=arrow.matrix))
}


matrixFill_uncompiled <- function(edit.matrix, arrow.matrix, arrow.labels, time=TRUE){
  
  n = dim(edit.matrix)[1] #last column 
  m = dim(edit.matrix)[2] #last row
  
  #Counts for printing time
  totalTime = n*m
  shortCount = 0
  completedTime = 0
  onePercent = ceiling(totalTime*0.001)
  
  #ptm <- proc.time()
  
  for(i in 2:n){
    
    for(j in 2:m){
      
      #cpu.time = ptm-proc.time()
      #print(paste("[",i,",",j,"] = ", cpu.time[2], sep=""))
      
      left = edit.matrix[i,j-1]
      diag = edit.matrix[i-1,j-1]
      top = edit.matrix[i-1,j]
      neighbors = c(left, diag, top)
      
      letterH1 = colnames(edit.matrix)[j-1]
      letterH2 = colnames(edit.matrix)[j]
      letterV1 = rownames(edit.matrix)[i-1]
      letterV2 = rownames(edit.matrix)[i]
      
      chars = sort(unique(c(colnames(edit.matrix), rownames(edit.matrix))))
      if(length(chars) > 2){
        char = chars[3]
      } else{
        char = chars[2]
      }
      
      #Start of medication
      startMedH = (letterH2 == char && letterH1 != char)
      startMedV = (letterV2 == char && letterV1 != char)
      
      #Continuing medication
      contMedH = (letterH2 == char && letterH1 == char)
      contMedV = (letterV2 == char && letterV1 == char)
      
      #End of medication
      endMedH = (letterH2 == '∅' && letterH1 != '∅') && letterH1 != "" ## NOTE
      endMedV = (letterV2 == '∅' && letterV1 != '∅')  && letterV1 != ""
      
      #Continuing gap
      contGapH = (letterH2 == '∅' && (letterH1 == '∅' || letterH1 == ""))
      contGapV = (letterV2 == '∅' && (letterV1 == '∅' || letterV1 == "" ))
      
      #Initialize values
      cellValue = -1
      direction = ""
      
      #New value ------------------
      
      #Same case scenario ------------------
      if((startMedH == TRUE && startMedV == TRUE) ||
         (contMedH == TRUE && contMedV == TRUE) ||
         (endMedH == TRUE && endMedV == TRUE) ||
         (contGapH == TRUE && contGapV == TRUE)) {
        cellValue = min(neighbors)
        direction = 2 #diag
      
      #Scenario switching ------------------
      } else if((startMedH == TRUE && contMedV == TRUE) ||
              (contMedH == TRUE && startMedV == TRUE) ||
              (contMedH == TRUE && endMedV == TRUE  ) ||
              (endMedH == TRUE && contMedV == TRUE  ) ||
              (endMedH == TRUE && contGapV == TRUE  ) ||
              (contGapH == TRUE && endMedV == TRUE  )) {
        cellValue = max(neighbors)
        direction = which(neighbors == max(neighbors))[1]
      
      #Opposite scenario ------------------
      } else if((startMedH == TRUE && endMedV == TRUE) ||
              (contMedH == TRUE && contGapV == TRUE) ||
              (endMedH == TRUE && startMedV == TRUE) ||
              (contGapH == TRUE && contMedV == TRUE)) {
        cellValue = max(neighbors) + 1
        direction = which(neighbors == max(neighbors))[1]
      
      #Gap involved ------------------
      } else if((startMedH == TRUE && contGapV == TRUE) ||
              (contGapH == TRUE && startMedV == TRUE)) {
        cellValue = min(neighbors) + 1
        direction = which(neighbors == min(neighbors))[1]
      }
      #Assign the new Value ------------------
      edit.matrix[i, j] = cellValue
      
      #Assign the arrow direction (L:left, D:diag, T:top)
      arrow.matrix[i,j] = arrow.labels[direction]
      
      #Print time ------------------
      if(time == TRUE){
        if(shortCount > onePercent){
          percentage = format(round(((completedTime/totalTime)*100), 2), nsmall = 2)
          print(paste(percentage, "%", sep=""))
          shortCount = 0
        }
        completedTime = completedTime + 1
        shortCount = shortCount + 1
      }
    } #end of j
  } #end of i
  
  write.table(edit.matrix, "azythromicin-p1-p2-matrix-armin.csv")
  return(list(edit=edit.matrix, arrow=arrow.matrix))
}

matrixFill <- cmpfun(matrixFill_uncompiled)

matrixTraceback <- function(edit.matrix){
  i = dim(edit.matrix)[1] #last column 
  j = dim(edit.matrix)[2] #last row
  
  alignment = c()
  nseq1 = ""
  nseq2 = ""
  
  while(i >= 1 && j >= 1){
    
    letterH2 = colnames(edit.matrix)[j]
    letterV2 = rownames(edit.matrix)[i]
    if(letterH2=="" && letterV2==""){
      break
    }
    
    step = c(i,j, edit.matrix[i,j])
    alignment = rbind(alignment, step)
    
    if(letterH2==letterV2){ #diagonal  
      i = i - 1
      j = j - 1
      nseq1 = c(letterH2, nseq1)
      nseq2 = c(letterV2, nseq2)
    } else{ #left or top? 
      
      if(i==1){
        direction = "left"
      }
      else if(j==1){
        direction = "top"
      }
      else if(edit.matrix[i,j-1] <= edit.matrix[i-1,j]){
        direction = "left"
      }
      else{
        direction = "top"
      }
      
      if(direction=="left"){ #left
        j = j - 1
        nseq1 = c(letterH2, nseq1)
        nseq2 = c("_", nseq2)
      }
      else{ #top
        i = i - 1
        nseq1 = c("_", nseq1)
        nseq2 = c(letterV2, nseq2)
      }
    }
    
  }
  colnames(alignment) = c("i", "j", "v")
  
  
  #alig = as.data.frame(rbind(unlist(strsplit(nseq1, " ")), unlist(strsplit(nseq2, " "))))
  alig = as.data.frame(rbind(nseq1, nseq2))
  rownames(alig) = c("patientH", "patientV")
  
  #orig = rbind(sequence1, sequence2)
  
  return(alig)
}

medalPairwise <- function(sequence1, sequence2, verbose=FALSE) {
  
  # First, Initialization Step ======================
  mi = matrixInitialization(sequence1, sequence2)
  edit.matrix = mi$edit
  arrow.matrix = mi$arrow
  #arrow.labels=c("left", "diag", "top")
  arrow.labels=c("←", "↖", "↑")
  
  if(verbose==TRUE){
    print("::Initialization [done]")
  }
  
  # Second, Matrix Fill Step ======================
  
  if(verbose==TRUE){
    mf = matrixFill(edit.matrix, arrow.matrix, arrow.labels, TRUE)
  }
  else{
    mf = matrixFill(edit.matrix, arrow.matrix, arrow.labels, FALSE)
  }
  edit.matrix = mf$edit
  arrow.matrix = mf$arrow
  
  if(verbose==TRUE){
    print("::Matrix Fill [done]")
  }
  
  
  # Third, Traceback Step ======================
  
  #alignment = matrixTraceback(edit.matrix, arrow.matrix, arrow.labels)
  alignment = matrixTraceback(edit.matrix)
  
  paste(unlist(alignment[1,]), collapse="")
  paste(unlist(alignment[2,]), collapse="")
  
  distance = length(which(alignment[1,]=="_"))
  paste("Edit-distance:", distance, sep=" ")
  
  if(verbose==TRUE){
    print("::Traceback [done]")
  }
  return(list(distance=distance, alignment=alignment))
}




getSequences <- function(pat1, pat2) {
  
  if(is.numeric(pat1$patientID)){
    pat1.label = paste("patient", unique(pat1$patientID), sep="")
  } else {
    pat1.label = unique(pat1$patientID)
  }
  if(is.numeric(pat2$patientID)){
    pat2.label = paste("patient", unique(pat2$patientID), sep="")
  } else {
    pat2.label = unique(pat2$patientID)
  }
  
  
  
  medications = as.vector(sort(unique(rbind(pat1,pat2)[,"medication"])))
  maxTime = max(as.numeric(as.character(rbind(pat1,pat2)[,"end"])), na.rm=TRUE)
  
  sequences = list()
  
  #Per Patient
  for(medication in medications){
    ind.pat1 = which(pat1$medication==medication)
    ind.pat2 = which(pat2$medication==medication)
    
    combined = rbind(pat1[ind.pat1,],pat2[ind.pat2,])
    maxTime = round(max(as.numeric(as.character(combined[,"end"]))))
    minTime = round(min(as.numeric(as.character(combined[,"start"]))))
    dif=maxTime-minTime + 1
    
    #If there is information for Patient 1
    if(length(ind.pat1)>0){
      
      #Create an empty matrix
      med1 = as.data.frame(matrix(rep('∅', dif), nrow = 1), stringsAsFactors = FALSE)
      rownames(med1) = c(pat1.label)
      colnames(med1) = paste("day", minTime:maxTime, sep="")
      
      #Fill in the values
      for(ind in ind.pat1){
        start = round(as.numeric(as.character(pat1[ind,"start"])))
        end = round(as.numeric(as.character(pat1[ind,"end"])))
        char = as.character(pat1[ind,"medication"])
        days = paste("day", start:end, sep="")
        med1[1, days] = toupper(substr(char,1,1))
      }
      
      
    } 
    #If there is information for Patient 2
    if(length(ind.pat2)>0){
      #Create an empty matrix
      med2 = as.data.frame(matrix(rep('∅', dif), nrow = 1), stringsAsFactors = FALSE)
      rownames(med2) = c(pat2.label)
      colnames(med2) = paste("day", minTime:maxTime, sep="")
      
      #Fill in the values
      for(ind in ind.pat2){
        start = round(as.numeric(as.character(pat2[ind,"start"])))
        end = round(as.numeric(as.character(pat2[ind,"end"])))
        char = as.character(pat2[ind,"medication"])
        days = paste("day", start:end, sep="")
        med2[1, days] = toupper(substr(char,1,1))
      }
    } 
    
    
    #Bind the medication per patients
    if(length(ind.pat1) > 0 && length(ind.pat2) == 0){
      sequences[[medication]]=med1
    } else if(length(ind.pat1) == 0 && length(ind.pat2) > 0){
      sequences[[medication]]=med2
    } else{
      sequences[[medication]]=rbind(med1, med2)
    }
    
  }
  
  return(sequences)
}

compactSequence1 <- function(sequence) {
  
  seq = c()
  count = 0
  for(i in 1:length(sequence)){
    if(sequence[i] != "∅"){
      seq = c(seq, sequence[i])
      count = 0
    }
    else{
      count = count + 1
      if(count <= 2){
        seq = c(seq, sequence[i])
      }
      
    }
  }
  
  return(seq)
  
}

compactSequence <- function(sequence) {
  
  seq = c()
  count.empty = 0
  count.med = 0
  for(i in 1:length(sequence)){
    if(sequence[i] != "∅"){
      count.empty = 0
      count.med = count.med + 1
      if(count.med <= 5){
        seq = c(seq, sequence[i])
      }
    }
    else{
      count.med = 0
      count.empty = count.empty + 1
      if(count.empty <= 2){
        seq = c(seq, sequence[i])
      }
      
    }
  }
  
  return(seq)
  
}

medalDistance <- function(pat1, pat2, verbose = FALSE){
  medal.distance = 0
  size.total = 0
  
  sequences = getSequences(pat1, pat2)
  
  medications = names(sequences)
  
  for(medication in medications){
    distance = 0
    size = 0
    pat1.empty = 0
    pat1.empty.orig = 0
    pat1.med = 0
    pat1.med.orig = 0
    pat2.empty = 0
    pat2.empty.orig = 0
    pat2.med = 0
    pat2.med.orig = 0
    
    if(dim(sequences[[medication]])[1] == 2){
      
      seq1 = unlist(sequences[[medication]][1,])
      #sequence1 = compactSequence(seq1)
      sequence1 = seq1
      seq2 = unlist(sequences[[medication]][2,])
      #sequence2 = compactSequence(seq2)
      sequence2 = seq2
      
      #Statistics
      ind = which(names(table(sequence1)) == "∅")
      pat1.empty = table(sequence1)[ind]
      ind = which(names(table(sequence1)) != "∅")
      pat1.med = table(sequence1)[ind]
      ind = which(names(table(seq1)) == "∅")
      pat1.empty.orig = table(seq1)[ind]
      ind = which(names(table(seq1)) != "∅")
      pat1.med.orig = table(seq1)[ind]
      
      
      ind = which(names(table(sequence2)) == "∅")
      pat2.empty = table(sequence2)[ind]
      ind = which(names(table(sequence2)) != "∅")
      pat2.med = table(sequence2)[ind]
      ind = which(names(table(seq2)) == "∅")
      pat2.empty.orig = table(seq2)[ind]
      ind = which(names(table(seq2)) != "∅")
      pat2.med.orig = table(seq2)[ind]
      
      md = medalPairwise(sequence1, sequence2, FALSE)
      distance = md$distance
      size = length(md$alignment)
      
    }
    
    else{
      if(rownames(sequences[[medication]]) == "patient1"){
        
        seq1 = unlist(sequences[[medication]][1,])
        sequence1 = compactSequence(seq1)
        
        ind = which(names(table(sequence1)) == "∅")
        pat1.empty = table(sequence1)[ind]
        ind = which(names(table(sequence1)) != "∅")
        pat1.med = table(sequence1)[ind]
        ind = which(names(table(seq1)) == "∅")
        pat1.empty.orig = table(seq1)[ind]
        ind = which(names(table(seq1)) != "∅")
        pat1.med.orig = table(seq1)[ind]
        
        distance = length(sequence1)
        size = distance
      }
      else{
        
        seq2 = unlist(sequences[[medication]][1,])
        sequence2 = compactSequence(seq2)
        
        ind = which(names(table(sequence2)) == "∅")
        pat2.empty = table(sequence2)[ind]
        ind = which(names(table(sequence2)) != "∅")
        pat2.med = table(sequence2)[ind]
        ind = which(names(table(seq2)) == "∅")
        pat2.empty.orig = table(seq2)[ind]
        ind = which(names(table(seq2)) != "∅")
        pat2.med.orig = table(seq2)[ind]
        
        distance = length(sequence2)
        size = distance
      }
    }
    if(verbose == TRUE){
    print(paste(medication, " [",
                "patient1 = ", pat1.empty, " (", pat1.empty.orig,  ") ", 
                pat1.med, " (", pat1.med.orig,")", ";  ",
                "patient2 = ", pat2.empty, " (", pat2.empty.orig,  ") ", 
                pat2.med,  " (", pat2.med.orig,")", ";  ",
                "distance = ", distance, ";  ",
                "size = ", size,
                "]", sep=""))
    }
    medal.distance = medal.distance + (distance*size)
  
    size.total = size.total + size
    
  }
  
  if(size.total > 0){
    medal.distance = medal.distance / size.total
  }
  return(medal.distance)
}



