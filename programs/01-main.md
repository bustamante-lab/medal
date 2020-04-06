Medication Alignment Algorithm (Medal)
================
Arturo Lopez Pineda
2020-04-05

yes

# Step 0. Load required libraries

``` r
remove(list=ls())

#General
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────── tidyverse 1.2.1 ──

    ## ✓ ggplot2 3.3.0     ✓ purrr   0.3.3
    ## ✓ tibble  2.1.3     ✓ dplyr   0.8.1
    ## ✓ tidyr   0.8.3     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(here)
```

    ## here() starts at /Users/lopezpia/Documents/GitHub/medal

``` r
#Plotting libraries
library(ggplot2)
library(ggpubr)
```

    ## Loading required package: magrittr

    ## 
    ## Attaching package: 'magrittr'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     set_names

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract

``` r
library(ggrepel)

#Clustering libraries
library(factoextra)
```

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
library(NbClust)
library(aricode)
library(tsne)
library(NMF) #for cluster purity and entropy
```

    ## Loading required package: pkgmaker

    ## Loading required package: registry

    ## Loading required package: rngtools

    ## Loading required package: cluster

    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''
    
    ## Warning: namespace 'Biobase' is not available and has been replaced
    ## by .GlobalEnv when processing object ''

    ## NMF - BioConductor layer [NO: missing Biobase] | Shared memory capabilities [NO: bigmemory] | Cores 11/12

    ##   To enable the Bioconductor layer, try: install.extras('
    ## NMF
    ## ') [with Bioconductor repository enabled]
    ##   To enable shared memory capabilities, try: install.extras('
    ## NMF
    ## ')

    ## 
    ## Attaching package: 'NMF'

    ## The following object is masked from 'package:aricode':
    ## 
    ##     entropy

``` r
#Additional functions
#source(here("programs", "fun-medal.R")) #PyMedal is used instead
source(here("programs", "fun-support.R"))
source(here("programs", "fun-plot.R"))
```

    ## 
    ## Attaching package: 'reshape2'

    ## The following object is masked from 'package:tidyr':
    ## 
    ##     smiths

    ## 
    ## Attaching package: 'scales'

    ## The following object is masked from 'package:purrr':
    ## 
    ##     discard

    ## The following object is masked from 'package:readr':
    ## 
    ##     col_factor

    ## 
    ## ---------------------
    ## Welcome to dendextend version 1.13.4
    ## Type citation('dendextend') for how to cite the package.
    ## 
    ## Type browseVignettes(package = 'dendextend') for the package vignette.
    ## The github page is: https://github.com/talgalili/dendextend/
    ## 
    ## Suggestions and bug-reports can be submitted at: https://github.com/talgalili/dendextend/issues
    ## Or contact: <tal.galili@gmail.com>
    ## 
    ##  To suppress this message use:  suppressPackageStartupMessages(library(dendextend))
    ## ---------------------

    ## 
    ## Attaching package: 'dendextend'

    ## The following object is masked from 'package:ggpubr':
    ## 
    ##     rotate

    ## The following object is masked from 'package:stats':
    ## 
    ##     cutree

    ## 
    ## Attaching package: 'MASS'

    ## The following object is masked from 'package:dplyr':
    ## 
    ##     select

 

# Step 1. Read file

``` r
#---
#File with patient ID (de-ID), comorbidities, initial clinical presentation, etc.
profiles = read.csv("../../data/patients.csv", stringsAsFactors = FALSE)

#Variables for stratification
strats <- c("is_male", "NHW", "OCD", "foodprob", "anx", 
            "emotional", "mood", "agg", "sch", "reg", "sleep", "tics")

#Select only a few columns
profiles = profiles[,c("id", "age_onset", "age_1st_appt", strats)]

#Get unique patient IDs ordered
patients = sort(unique(profiles$id))

#---
#File with events
events = read.csv("../../data/medications.csv", stringsAsFactors = FALSE)

#Only select rows for the same patients listed in Profiles
rows = which(events$id %in% patients)

#Select only a few columns
events = events[rows,c("id", "medication", "start", "end")]

#---
#File with clinical evaluations (outcomes)
outcomes = read.csv("../../data/outcomes.csv", stringsAsFactors = FALSE)

#Only select rows for the same patients listed in Profiles
rows = which(outcomes$id %in% patients)

#Select only a few columns
outcomes = outcomes[rows,c("id", "gi_new", "daysSinceBirth")]
```

 

# Step 2. Clean Data

``` r
# Group by class of medication
medgroups = list()
medgroups$penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin")
medgroups$cephalosporin = c("cephalexin", "cefadroxil")
medgroups$macrolide = c("azithromycin")
medgroups$nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin")
medgroups$corticosteroid.oral = c("prednisone", "maintenance prednisone", "decadron")
medgroups$corticosteroid.iv = c("solumedrol")
medgroups$immunoglobulins = c("ivig")
medgroups$dmard = c("rituximab", "methotrexate", "cellcept")


data = cleanEvents(events, medgroups)

data = rightCensoring(data, 2)

write.csv(data, "../../data/data-matrix-clean.csv")
```

 

# Step 3. Create a distance matrix

``` r
#-------
#Calling pyMEDAL

system('python3 ../pymedal/pymedal.py ../../data/data-matrix-clean.csv', wait=TRUE)


distMatrix = read.table("distance_mat.txt")
patientIDs = read.table("patientID.txt", sep=",")

pID = as.vector(unlist(patientIDs))
colnames(distMatrix) = pID
rownames(distMatrix) = pID
```

 

# Step 4 Choose number of clusters

``` r
#d = scale(distMatrix, center=FALSE)
d = distMatrix

# Remove outlier (patient 37)
#e=d
#remInd = which(names(d)=="37")
#d = d[-remInd,]
#d = d[,-remInd]
k=4

# Elbow method
elbow <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "wss") +
  geom_vline(xintercept = k, linetype = "dashed", color="#5581B0", size=0.6) +
  labs(title = "Elbow method",
       y="Total within-clusters sum of squares")

# Silhouette method
silhouette <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "silhouette", 
                           print.summary = FALSE, barcolor = "white") +
  geom_vline(xintercept = 2, linetype = "dashed", color="#5581B0", size=0.6) +
  geom_vline(xintercept = 4, linetype = "dashed", color="#5581B0", size=0.6) +
  labs(title = "Silhouette method")


# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
gapStat <- fviz_nbclust(x=d, diss=as.dist(d), hcut, nstart = 25, 
                        method = "gap_stat", nboot = 50, print.summary = FALSE,
                        maxSE=list(method="Tibs2001SEmax", SE.factor=1)) +
  labs(title = "Gap statistic method")

gpanels <- ggarrange(elbow, silhouette, gapStat,
                     labels = c("A", "B", "C"),
                     ncol = 3, nrow = 1, legend="bottom", 
                     align="v", common.legend = FALSE)
ggexport(gpanels, filename="../images/Figure1-num-clusters.png", height = 1200, width = 4000, res=300)
```

    ## file saved to ../images/Figure1-num-clusters.png

 

# Step 5. Plot clusters

``` r
#Colors to be used
color.vector = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")
color.vector2 = c("1"="#e8a631", "2"="#ca3542", "3"="#00a572", "4"="#0080ff")
#color.vector3 = c("1"="mediumorchid3", "2"="darkturquoise", "3"="olivedrab3", "4"="orangered3")

#Plot Clustering strategies

dend = getDendrogram(d, k, color.vector)
gDend <- plotDendrogram(dend, k)
gMDSclus12 <- plotMDS(d, as.character(cutree(dend, k)), color.vector, 1, 2, "MDS (hierarchical clustering)")
```

    ## initial  value 18.914825 
    ## final  value 18.914116 
    ## converged

``` r
gMDSclus34 <- plotMDS(d, as.character(cutree(dend, k)), color.vector, 3, 4, "MDS (hierarchical clustering)")
```

    ## initial  value 18.914825 
    ## final  value 18.914116 
    ## converged

``` r
kmeans = getKMeansClusteringPCA(d, k)
gMDSkmeans12 <- plotMDS(d, kmeans$cluster, color.vector, 1, 2, "MDS (k-means)")
```

    ## initial  value 18.914825 
    ## final  value 18.914116 
    ## converged

``` r
gMDSkmeans34 <- plotMDS(d, kmeans$cluster, color.vector, 3, 4, "MDS (k-means)")
```

    ## initial  value 18.914825 
    ## final  value 18.914116 
    ## converged

``` r
tsne1 = getHierarchicalClusteringTSNE(d, k, perplexity = 4)
gTSNEclus <- plotTSNE(tsne1, color.vector, "TSNE (hierarchical clustering)")
 
tsne2 = getKMeansClusteringTSNE(d, k, perplexity = 4)
gTSNEkmeans <- plotTSNE(tsne2, color.vector, "TSNE (k-means)")


# Combine plots and save
gpanels <- ggarrange(gDend, 
                     ggarrange(gMDSclus12, gTSNEclus,
                               labels = c("B", "C"),
                               align = "hv",
                               legend="bottom", common.legend = TRUE),
                     ncol = 1, nrow=2, 
                     labels = c("A"),
                     align = "h",
                     legend="bottom", common.legend = TRUE)
ggexport(gpanels, filename="../images/Figure2-dendro-mds.png", height = 4000, width = 4000, res=300)
```

    ## file saved to ../images/Figure2-dendro-mds.png

 

# Step 6. Comparison to K-means

``` r
recoded = kmeans$cluster
recoded[which(recoded == 3)] = 1
recoded[which(recoded == 4)] = 2


#Calculating Normalized Mutual Information
round(NMI(recoded, cutree(dend,2), variant="sum"), digits=2)
```

    ## [1] 0.16

``` r
round(NMI(kmeans$cluster, cutree(dend,4), variant="sum"), digits=2)
```

    ## [1] 0.37

``` r
#https://course.ccs.neu.edu/cs6140sp15/7_locality_cluster/Assignment-6/NMI.pdf

#Calculating Cluster purity
round(purity(recoded, cutree(dend,2)), digits=2)
```

    ## [1] 0.71

``` r
round(purity(kmeans$cluster, cutree(dend,4)), digits=2)
```

    ## [1] 0.6

``` r
#https://www.rdocumentation.org/packages/NMF/versions/0.21.0/topics/purity

#Calculating Cluster entropy
round(entropy(recoded, cutree(dend,2)), digits=2)
```

    ## [1] 0.81

``` r
round(entropy(kmeans$cluster, cutree(dend,4)), digits=2)
```

    ## [1] 0.6

``` r
#https://www.rdocumentation.org/packages/NMF/versions/0.21.0/topics/purity
```

 

# Step 7. Save cluster assignment

``` r
#Saving the cluster to profiles
assignment = cutree(dend,k)
index = which(profiles$id %in% names(assignment))
write.csv(cbind(profiles[index,], cluster=assignment),
          "../../data/data-matrix-profiles-cluster.csv")
```
