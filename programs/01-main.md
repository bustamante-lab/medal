Medication Alignment Algorithm (Medal)
================
Arturo Lopez Pineda
2020-05-14

# Step 0. Load required libraries

``` r
remove(list=ls())

#General
library(tidyverse)
```

    ## ── Attaching packages ────────────────────────────────────────────────────────────────────────── tidyverse 1.3.0 ──

    ## ✓ ggplot2 3.3.0     ✓ purrr   0.3.3
    ## ✓ tibble  3.0.0     ✓ dplyr   0.8.5
    ## ✓ tidyr   1.0.2     ✓ stringr 1.4.0
    ## ✓ readr   1.3.1     ✓ forcats 0.5.0

    ## Warning: package 'tibble' was built under R version 3.6.2

    ## ── Conflicts ───────────────────────────────────────────────────────────────────────────── tidyverse_conflicts() ──
    ## x dplyr::filter() masks stats::filter()
    ## x dplyr::lag()    masks stats::lag()

``` r
library(here)
```

    ## here() starts at /Users/lopezpia/Documents/GitHub/medal

``` r
library(dplyr)
library(knitr)

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
library(ggplotify)

#Clustering libraries
library(factoextra)
```

    ## Welcome! Want to learn more? See two factoextra-related books at https://goo.gl/ve3WBa

``` r
library(NbClust)
library(aricode)
library(Rtsne)
library(dendextend)
```

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

``` r
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
#Heatmap
library("pheatmap")

#ANOVA
library(lme4)
```

    ## Loading required package: Matrix

    ## 
    ## Attaching package: 'Matrix'

    ## The following objects are masked from 'package:tidyr':
    ## 
    ##     expand, pack, unpack

``` r
library(lmerTest)
```

    ## 
    ## Attaching package: 'lmerTest'

    ## The following object is masked from 'package:lme4':
    ## 
    ##     lmer

    ## The following object is masked from 'package:stats':
    ## 
    ##     step

``` r
#Additional functions
#source(here("programs", "fun-medal.R")) #PyMedal is used instead
source(here("programs", "fun-support.R"))
source(here("programs", "fun-plot.R"))
```

 

# Step 1. Read file

``` r
#---
#File with patient ID (de-ID), comorbidities, initial clinical presentation, etc.
patients.og <- read_csv(here("../data", "patients.csv")) %>%
  column_to_rownames(var="id") 
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double(),
    ##   clinical_history = col_logical(),
    ##   Imaging = col_logical(),
    ##   Medication_History = col_logical(),
    ##   HLA_genotyping = col_logical()
    ## )

    ## See spec(...) for full column specifications.

``` r
#Select only a few columns
patients <- patients.og %>%
  select("age_onset", "age_1st_appt", "is_male", "NHW", "OCD", "foodprob", "anx", 
            "emotional", "mood", "agg", "sch", "reg", "sleep", "tics") %>%
  drop_na()

#---
#File with events
events.og <- read_csv(here("../data", "medications.csv")) 
```

    ## Parsed with column specification:
    ## cols(
    ##   id = col_double(),
    ##   age_1st_appt = col_double(),
    ##   age_onset = col_double(),
    ##   medication = col_character(),
    ##   start = col_double(),
    ##   end = col_double(),
    ##   Duration = col_double(),
    ##   firstAppointment = col_double(),
    ##   Frequency = col_character(),
    ##   last_FU = col_double()
    ## )

``` r
#Select only events related to patients in the previous file
events <- events.og %>%
  filter(id %in% rownames(patients)) %>%
  select("id", "medication", "start", "end")

#---
#File with clinical evaluations (outcomes)
outcomes.og <- read_csv(here("../data","outcomes.csv"))
```

    ## Parsed with column specification:
    ## cols(
    ##   .default = col_double()
    ## )
    ## See spec(...) for full column specifications.

``` r
#Select only rows for the same patients listed in Profiles
outcomes <- outcomes.og %>%
  filter(id %in% rownames(patients)) %>%
  select("id", "gi_new", "cbiTotal", "daysSinceBirth")

#---
#Censoring
years = 2
```

 

# Step 2. Clean Data and Right Censor

``` r
# Group by class of medication
medgroups <- lst(penicillin = c("penicillin v", "penicillin g", "amoxicillin", "augmentin"),
                cephalosporin = c("cephalexin", "cefadroxil"), 
                macrolide = c("azithromycin"),
                nsaid = c("ibuprofen", "naproxen", "indomethacin", "sulindac", "aspirin"),
                corticosteroid.oral = c("prednisone", "maintenance prednisone", "decadron"),
                corticosteroid.iv = c("solumedrol"),
                immunoglobulins = c("ivig"),
                dmard = c("rituximab", "methotrexate", "cellcept"))

medcolors= c("penicillin"="#66c2a5",
             "cephalosporin" = "#fc8d62",
             "macrolide" = "#8da0cb",
             "nsaid" = "#e7298a",
             "corticosteroid.oral" = "#a6d854",
             "corticosteroid.iv" = "#ffd92f",
             "antibody" = "green",
             "immunoglobulins" = "#e5c494",
             "dmard" = "#b3b3b3")


events <- events.og %>%
  cleanEvents(medgroups) %>%
  rightCensoring(years)

write_csv(events, here("../data", "data-matrix-clean.csv"))
```

 

# Step 3. Create a distance matrix

``` r
#-------
#Calling pyMEDAL

pymedal <- here("pymedal", "pymedal.py")
input <- here("../data", "data-matrix-clean.csv")

system(paste('python3', pymedal, input), wait=TRUE)


#Reading the output
distMatrix = read.table("distance_mat.txt")
patientIDs = read.table("patientID.txt", sep=",")

pID = as.vector(unlist(patientIDs))
colnames(distMatrix) = pID
rownames(distMatrix) = pID
```

 

# Step 4. Choose number of clusters

``` r
d = distMatrix

# Elbow method
elbow <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "wss") +
  #geom_vline(xintercept = k, linetype = "dashed", color="#5581B0", size=0.6) +
  labs(title = "Elbow method",
       y="Total within-clusters sum of squares") +
  theme_light()

# Silhouette method
silhouette <- fviz_nbclust(x=d, diss=as.dist(d), hcut, method = "silhouette", 
                           print.summary = FALSE, barcolor = "white") +
  #geom_vline(xintercept = 2, linetype = "dashed", color="#5581B0", size=0.6) +
  #geom_vline(xintercept = 4, linetype = "dashed", color="#5581B0", size=0.6) +
  labs(title = "Silhouette method") +
  theme_light()


# Gap statistic
# nboot = 50 to keep the function speedy. 
# recommended value: nboot= 500 for your analysis.
# Use verbose = FALSE to hide computing progression.
set.seed(123)
gapStat <- fviz_nbclust(x=d, diss=as.dist(d), hcut, nstart = 25, 
                        method = "gap_stat", nboot = 50, print.summary = TRUE,
                        maxSE=list(method="Tibs2001SEmax", SE.factor=1)) +
  labs(title = "Gap statistic method") +
  theme_light()

clest <- fviz_nbclust(x=d, diss=as.dist(d), hcut, nstart = 25, 
                        method = "gap_stat", nboot = 50, print.summary = TRUE,
                        maxSE=list(method="globalSEmax", SE.factor=3)) +
  labs(title = "Clest method") +
  theme_light()

numKPanels <- ggarrange(elbow, silhouette, gapStat, clest,
                     labels = c("A", "B", "C", "D"),
                     ncol = 2, nrow = 2, legend="bottom", 
                     align="v", common.legend = FALSE)
ggexport(numKPanels, filename=here("images", "Figure4-b-num-clusters.png"), height = 2500, width = 3000, res=300)
```

    ## file saved to /Users/lopezpia/Documents/GitHub/medal/images/Figure4-b-num-clusters.png

``` r
numKPanels
```

![](01-main_files/figure-gfm/num-clusters-1.png)<!-- -->

 

# Step 5. Plot clusters

``` r
k=6
set.seed(123)

#Colors to be used
#colors.hclust <- c("1"="#a6cee3", "2"="#ff7f00", "3"="#b2df8a", 
#                   "4"="#6a3d9a", "5"="#fcbba1", "6"="#e31a1c")
#colors.kmeans <- c("1"="#fdbf6f", "2"="#1f78b4", "3"="#cab2d6", 
#                   "4"="#33a02c", "5"="#fb9a99", "6"="#b15928")

colors.hclust <- c("1"="#1b9e77", "2"="#7570b3", "3"="#e6ab02", 
                   "4"="#d95f02", "5"="#a6761d", "6"="#e7298a")
colors.kmeans <- c("1"="#e6ab02", "2"="#d95f02", "3"="#e7298a", 
                   "4"="#a6761d", "5"="#7570b3", "6"="#1b9e77")




#Plot Clustering strategies
dend = getDendrogram(d, k, colors.hclust)
gDend <- plotDendrogram(dend, k)
gMDSclus12 <- plotMDS(d, as.character(cutree(dend, k)), colors.hclust, 1, 2, "Hierarchical clustering (hclust)")
```

    ## initial  value 18.914825 
    ## final  value 18.914116 
    ## converged

``` r
#gMDSclus34 <- plotMDS(d, as.character(cutree(dend, k)), colors.hclust, 3, 4, "MDS (hierarchical clustering)")

kmeans = getKMeansClusteringPCA(d, k)
gMDSkmeans12 <- plotMDS(d, kmeans$cluster, colors.kmeans, 1, 2, "K-means")
```

    ## initial  value 18.914825 
    ## final  value 18.914116 
    ## converged

``` r
#gMDSkmeans34 <- plotMDS(d, kmeans$cluster, colors.kmeans, 3, 4, "MDS (k-means)")


#TSNE
tsne1 = as.data.frame(Rtsne(d, perplexity=4, is_distance=TRUE, initial_dims=5, theta=0)$Y)

tsne2 = getHierarchicalClusteringTSNE(d, k, tsne1)
gTSNEclus <- plotTSNE(tsne2, colors.hclust, "TSNE (hierarchical clustering)")
 
tsne3 = getKMeansClusteringTSNE(d, k, tsne1)
gTSNEkmeans <- plotTSNE(tsne3, colors.kmeans, "TSNE (k-means)")


#Combine plots and save
gpanels <- ggarrange(gDend,
                     ggarrange(gMDSclus12, gTSNEclus,
                               labels = c("B", "C"),
                               align = "hv",
                               legend="bottom", common.legend = TRUE),
                     ncol = 1, nrow=2,
                     labels = c("A"),
                     align = "h",
                     legend="bottom", common.legend = TRUE)
ggexport(gpanels, filename=here("images", "Figure5-a-dendro-mds.png"), height = 4000, width = 4000, res=300)
```

    ## file saved to /Users/lopezpia/Documents/GitHub/medal/images/Figure5-a-dendro-mds.png

``` r
#gpanels


gpanels <- ggarrange(ggarrange(gMDSkmeans12, gTSNEkmeans,
                               align = "h",
                               legend="bottom", common.legend = TRUE),
                     ggarrange(gMDSclus12, gTSNEclus,
                               align = "h",
                               legend="bottom", common.legend = TRUE),
                     ncol = 1, nrow=2,
                     labels = c("A", "B"),
                     align = "h",
                     legend="bottom", common.legend = FALSE)
ggexport(gpanels, filename=here("images", "Figure5-b-dendro-mds.png"), height = 4000, width = 4000, res=300)
```

    ## file saved to /Users/lopezpia/Documents/GitHub/medal/images/Figure5-b-dendro-mds.png

``` r
kmeans_annot <- as.data.frame(kmeans["cluster"]) %>%
  rename(kmeans=cluster)
hclust_annot <- as.data.frame(factor(cutree(dend,k))) %>%
  rename(hclust = colnames(.)[1])
ann_colors <- list(kmeans = colors.kmeans, hclust = colors.hclust)


gHeatMap <- as.ggplot(pheatmap(d, cutree_cols = 6, cutree_rows = 6,
                    scale="none",
                    show_colnames = TRUE,
                    show_rownames = TRUE,
                    clustering_distance_rows = "minkowski",
                    clustering_distance_cols = "minkowski",
                    clustering_method = "ward.D",
                    labels_col = paste0(rownames(d)),
                    labels_row = paste0(rownames(d)),
                    annotation_col = kmeans_annot,
                    annotation_row = hclust_annot,
                    annotation_colors = ann_colors,
                    legend = TRUE,
                    cluster_rows = TRUE,
                    cluster_cols = TRUE))
```

![](01-main_files/figure-gfm/plot-clusters-1.png)<!-- -->

``` r
# Combine plots and save
gpanels <- ggarrange(gHeatMap, 
                     ggarrange(gMDSclus12, gMDSkmeans12,
                               labels = c("B", "C"),
                               align = "hv",
                               ncol = 2, nrow=1, 
                               legend="bottom", common.legend = FALSE),
                     labels = c("A", ""),
                     align = "hv",
                     ncol = 1, nrow=2)

ggexport(gpanels, filename=here("images", "Figure5-heatmap-mds.png"), height = 4000, width = 3000, res=300)
```

    ## file saved to /Users/lopezpia/Documents/GitHub/medal/images/Figure5-heatmap-mds.png

``` r
gpanels
```

![](01-main_files/figure-gfm/plot-clusters-2.png)<!-- -->

 

# Step 6. Comparison to K-means

``` r
NMI <- c()
Entropy <- c()
Purity <- c()

for(i in 1:10){
  kmeans = getKMeansClusteringPCA(d, i)$cluster
  hclust = cutree(getDendrogram(d, i, colors.hclust), i)
  
  #Calculating Normalized Mutual Information
  #https://course.ccs.neu.edu/cs6140sp15/7_locality_cluster/Assignment-6/NMI.pdf
  NMI = c(NMI, round(NMI(kmeans, hclust, variant="sum"), digits=2))
  
  #Calculating Cluster purity
  #https://www.rdocumentation.org/packages/NMF/versions/0.21.0/topics/purity
  Purity = c(Purity, round(purity(kmeans, hclust), digits=2))
  
  #Calculating Cluster entropy
  #https://www.rdocumentation.org/packages/NMF/versions/0.21.0/topics/purity
  Entropy = c(Entropy, round(entropy(kmeans, hclust), digits=2))

}

#Print table
results <- as.data.frame(cbind(k=c(1:10), NMI, Purity, Entropy))

results %>%
  kable()
```

|  k |  NMI | Purity | Entropy |
| -: | ---: | -----: | ------: |
|  1 | 1.00 |   1.00 |     NaN |
|  2 | 0.86 |   0.98 |    0.13 |
|  3 | 0.49 |   0.74 |    0.47 |
|  4 | 0.50 |   0.67 |    0.47 |
|  5 | 0.81 |   0.81 |    0.18 |
|  6 | 0.90 |   0.93 |    0.09 |
|  7 | 0.93 |   0.95 |    0.06 |
|  8 | 0.84 |   0.83 |    0.17 |
|  9 | 0.90 |   0.90 |    0.09 |
| 10 | 0.91 |   0.93 |    0.08 |

``` r
res <- bind_rows(results %>% 
                   select(k, value=NMI) %>% 
                   mutate(Test="NMI"),
                 results %>% 
                   select(k, value=Purity) %>% 
                   mutate(Test="Purity"),
                 results %>% 
                   select(k, value=Entropy) %>% 
                   mutate(Test="Entropy"))
  

clus.test <- ggplot(res, aes(x=k, y=value, group=Test, color=Test)) +
  geom_vline(xintercept=6, color="gray30", linetype="dashed", size=1.5) +
  geom_point(size=2.5) +
  geom_line(size=1) +
  #geom_vline(xintercept=7, color="gray40", linetype="dashed") +
  scale_colour_manual(values=c("#7fc97f", "#beaed4", "#fdc086")) +
  scale_x_continuous(limits=c(1,10), breaks=c(0:10)) +
  ggtitle("Comparison between K-means and Hierarchical Clustering")+
  theme_light() +
  theme(legend.position = "bottom")

gpanel <- ggarrange(numKPanels, clus.test,
                    labels = c("", "E"),
                     align = "h",
                     ncol = 1, nrow=2)
```

    ## Warning: Removed 1 rows containing missing values (geom_point).

    ## Warning: Removed 1 row(s) containing missing values (geom_path).

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

``` r
ggexport(gpanel, filename=here("images", "Figure4-metrics.png"), height = 3600, width = 2000, res=300)
```

    ## file saved to /Users/lopezpia/Documents/GitHub/medal/images/Figure4-metrics.png

``` r
gpanel
```

![](01-main_files/figure-gfm/comparison-1.png)<!-- -->

 

# Step 7. Save cluster assignment

``` r
#Saving the cluster to profiles
assignment <- hclust_annot
assignment <- kmeans_annot 

out <- assignment %>%
  rename(cluster=kmeans) %>%
  rownames_to_column(var="id")

write.csv(out,here("dataExample", "cluster-assignment.csv"))

profiles <- patients %>%
  rownames_to_column(var="id") %>%
  filter(id %in% rownames(assignment)) %>%
  mutate(cluster=unlist(assignment))
```

 

# Step 8. Evaluate clusters

``` r
#Get all paths
gPath1 <- plotTimeSeriesDrug(1, events, profiles, medcolors, medgroups, years)
gPath2 <- plotTimeSeriesDrug(2, events, profiles, medcolors, medgroups, years)
gPath3 <- plotTimeSeriesDrug(3, events, profiles, medcolors, medgroups, years)
gPath4 <- plotTimeSeriesDrug(4, events, profiles, medcolors, medgroups, years)
gPath5 <- plotTimeSeriesDrug(5, events, profiles, medcolors, medgroups, years)
gPath6 <- plotTimeSeriesDrug(6, events, profiles, medcolors, medgroups, years)



#Get all impairment scores
gGIS1 <- plotScores(1, outcomes, "gi_new", c(0,100,5), profiles, 2, "gray50", "Global Impairment", TRUE, FALSE)
gGIS2 <- plotScores(2, outcomes, "gi_new", c(0,100,5), profiles, 2, "gray50", "Global Impairment", TRUE, FALSE)
gGIS3 <- plotScores(3, outcomes, "gi_new", c(0,100,5), profiles, 2, "gray50", "Global Impairment", TRUE, FALSE)
gGIS4 <- plotScores(4, outcomes, "gi_new", c(0,100,5), profiles, 2, "gray50", "Global Impairment", TRUE, FALSE)
gGIS5 <- plotScores(5, outcomes, "gi_new", c(0,100,5), profiles, 2, "gray50", "Global Impairment", TRUE, FALSE)
gGIS6 <- plotScores(6, outcomes, "gi_new", c(0,100,5), profiles, 2, "gray50", "Global Impairment", TRUE, FALSE)
   
gcbi1 <- plotScores(1, outcomes, "cbiTotal", c(0,96,4), profiles, 2, "gray50", "Caregiver Burden", TRUE, TRUE)
gcbi2 <- plotScores(2, outcomes, "cbiTotal", c(0,96,4), profiles, 2, "gray50", "Caregiver Burden", TRUE, TRUE)
gcbi3 <- plotScores(3, outcomes, "cbiTotal", c(0,96,4), profiles, 2, "gray50", "Caregiver Burden", TRUE, TRUE)
gcbi4 <- plotScores(4, outcomes, "cbiTotal", c(0,96,4), profiles, 2, "gray50", "Caregiver Burden", TRUE, TRUE)
gcbi5 <- plotScores(5, outcomes, "cbiTotal", c(0,96,4), profiles, 2, "gray50", "Caregiver Burden", TRUE, TRUE)
gcbi6 <- plotScores(6, outcomes, "cbiTotal", c(0,96,4), profiles, 2, "gray50", "Caregiver Burden", TRUE, TRUE)



#Save plot
gpanels <- ggarrange(gPath1, gPath2, gPath3, gPath4, gPath5, gPath6,
                     gGIS1, gGIS2, gGIS3, gGIS4, gGIS5, gGIS6,
                     gcbi1, gcbi2, gcbi3, gcbi4, gcbi5, gcbi6,
                     #labels = c("Clus1", "Clus2", "Clus3", "Clus4", "Clus5", "Clus6"),
                     heights = c(6,1,1),
                     align = "hv",
                     ncol = 6, nrow = 3, legend="none", common.legend = FALSE)
```

    ## Warning: Removed 97 rows containing missing values (position_stack).

    ## Warning: Removed 2 row(s) containing missing values (geom_path).

    ## Warning: Removed 97 rows containing missing values (position_stack).

    ## Warning: Removed 1 row(s) containing missing values (geom_path).

    ## Warning: Removed 164 rows containing missing values (position_stack).

    ## Warning: Removed 87 row(s) containing missing values (geom_path).

    ## Warning: Removed 30 rows containing missing values (position_stack).

    ## Warning: Removed 148 rows containing missing values (position_stack).

    ## Warning: Removed 7 row(s) containing missing values (geom_path).

    ## Warning: Removed 45 rows containing missing values (position_stack).

    ## Warning: Removed 10 row(s) containing missing values (geom_path).

    ## Warning: Graphs cannot be horizontally aligned unless the axis parameter is set.
    ## Placing graphs unaligned.

``` r
ggexport(gpanels, filename=here("images", "Figure6-clusters.png"),
         height = 4000, width = 5000, res=300)
```

    ## file saved to /Users/lopezpia/Documents/GitHub/medal/images/Figure6-clusters.png

``` r
gpanels
```

![](01-main_files/figure-gfm/evaluate-cluster-1.png)<!-- -->

 

# Step 9. ANOVA

``` r
#Get profiles for patients in each cluster
pat1 <- getClusterScores(profiles, 1, outcomes, years) %>% mutate(cluster="1")
pat2 <- getClusterScores(profiles, 2, outcomes, years) %>% mutate(cluster="2")
pat3 <- getClusterScores(profiles, 3, outcomes, years) %>% mutate(cluster="3")
pat4 <- getClusterScores(profiles, 4, outcomes, years) %>% mutate(cluster="4")
pat5 <- getClusterScores(profiles, 5, outcomes, years) %>% mutate(cluster="5")
pat6 <- getClusterScores(profiles, 6, outcomes, years) %>% mutate(cluster="6")

pat <- pat1 %>%
  bind_rows(pat2) %>%
  bind_rows(pat3) %>%
  bind_rows(pat4) %>%
  bind_rows(pat5) %>%
  bind_rows(pat6)

pat <- pat %>%
  mutate(one = case_when(cluster == 1 ~ 1, TRUE ~ 0)) %>%
  mutate(two = case_when(cluster == 2 ~ 1, TRUE ~ 0)) %>%
  mutate(three = case_when(cluster == 3 ~ 1, TRUE ~ 0)) %>%
  mutate(four = case_when(cluster == 4 ~ 1, TRUE ~ 0)) %>%
  mutate(five = case_when(cluster == 5 ~ 1, TRUE ~ 0)) %>%
  mutate(six = case_when(cluster == 6 ~ 1, TRUE ~ 0))

#Convert days to years
pat <- pat %>%
  mutate(years=daysSinceOnset/365)

#Plot the scores
gi_new <- ggplot(pat)+
  geom_smooth(aes(x=years, y=gi_new, group=cluster, col=cluster),
              method = "lm", size=2, se=FALSE, na.rm = TRUE, formula=y~x) +
  scale_color_manual(values = colors.kmeans) +
  theme_light()+
  ggtitle("Global Impairment Score by cluster")

cbiTotal <- ggplot(pat)+
  geom_smooth(aes(x=years, y=cbiTotal, group=cluster, col=cluster),
              method = "lm", size=2, se=FALSE, na.rm = TRUE, formula=y~x) +
  scale_color_manual(values = colors.kmeans) +
  theme_light()+
  ggtitle("Caregiver Burden by cluster")

gpanels <- ggarrange(gi_new, cbiTotal,
                     align = "hv",
                     ncol = 2, nrow = 1, legend="bottom", common.legend = TRUE)
ggexport(gpanels, filename=here("images", "Figure7-scores.png"),
         height = 1500, width = 3000, res=300)
```

    ## file saved to /Users/lopezpia/Documents/GitHub/medal/images/Figure7-scores.png

``` r
gpanels
```

![](01-main_files/figure-gfm/anova-1.png)<!-- -->

``` r
## Obtain summaries
summary(lmer(gi_new ~ cluster*years + (years | id), data = pat, REML = F))
```

    ## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
    ##   method [lmerModLmerTest]
    ## Formula: gi_new ~ cluster * years + (years | id)
    ##    Data: pat
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   4468.0   4536.1  -2218.0   4436.0      505 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.5921 -0.6057 -0.0361  0.5457  3.6545 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev. Corr 
    ##  id       (Intercept) 539.0    23.22         
    ##           years       129.7    11.39    -0.72
    ##  Residual             219.4    14.81         
    ## Number of obs: 521, groups:  id, 41
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error       df t value Pr(>|t|)   
    ## (Intercept)     33.1158    10.3135  42.0519   3.211  0.00254 **
    ## cluster2         5.1318    14.4075  40.1727   0.356  0.72356   
    ## cluster3        -4.1931    14.5985  42.2423  -0.287  0.77534   
    ## cluster4         6.4408    12.5250  40.7893   0.514  0.60986   
    ## cluster5         0.3012    14.5099  41.2740   0.021  0.98354   
    ## cluster6        24.8156    15.1526  40.8552   1.638  0.10916   
    ## years           -7.6932     5.9408  22.0965  -1.295  0.20869   
    ## cluster2:years  11.6572     8.2249  21.5192   1.417  0.17071   
    ## cluster3:years -13.6344    10.2222  42.3958  -1.334  0.18939   
    ## cluster4:years   0.2758     7.6187  25.0256   0.036  0.97141   
    ## cluster5:years -10.8159     8.6316  25.0612  -1.253  0.22175   
    ## cluster6:years  -4.8081     8.6326  21.9096  -0.557  0.58319   
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) clstr2 clstr3 clstr4 clstr5 clstr6 years  clst2: clst3:
    ## cluster2    -0.716                                                        
    ## cluster3    -0.706  0.506                                                 
    ## cluster4    -0.823  0.589  0.582                                          
    ## cluster5    -0.711  0.509  0.502  0.585                                   
    ## cluster6    -0.681  0.487  0.481  0.560  0.484                            
    ## years       -0.718  0.514  0.507  0.591  0.510  0.489                     
    ## clustr2:yrs  0.518 -0.719 -0.366 -0.427 -0.369 -0.353 -0.722              
    ## clustr3:yrs  0.417 -0.299 -0.655 -0.343 -0.297 -0.284 -0.581  0.420       
    ## clustr4:yrs  0.560 -0.401 -0.395 -0.698 -0.398 -0.381 -0.780  0.563  0.453
    ## clustr5:yrs  0.494 -0.354 -0.349 -0.407 -0.703 -0.336 -0.688  0.497  0.400
    ## clustr6:yrs  0.494 -0.354 -0.349 -0.407 -0.351 -0.717 -0.688  0.497  0.400
    ##             clst4: clst5:
    ## cluster2                 
    ## cluster3                 
    ## cluster4                 
    ## cluster5                 
    ## cluster6                 
    ## years                    
    ## clustr2:yrs              
    ## clustr3:yrs              
    ## clustr4:yrs              
    ## clustr5:yrs  0.537       
    ## clustr6:yrs  0.537  0.474

``` r
summary(lmer(gi_new ~ years + (years | id), 
             data = pat %>% filter(one==1), REML = F))$coefficients
```

    ##              Estimate Std. Error       df   t value   Pr(>|t|)
    ## (Intercept) 33.232018  10.309727 5.852491  3.223366 0.01870295
    ## years       -8.317455   5.876922 4.417811 -1.415274 0.22345300

``` r
summary(lmer(gi_new ~ years + (years | id), 
             data = pat %>% filter(two==1), REML = F))$coefficients
```

    ## boundary (singular) fit: see ?isSingular

    ##              Estimate Std. Error        df  t value     Pr(>|t|)
    ## (Intercept) 38.561533   4.751049  7.368225 8.116424 6.265316e-05
    ## years        3.801834   2.653336 74.634328 1.432850 1.560762e-01

``` r
summary(lmer(gi_new ~ years + (years | id), 
             data = pat %>% filter(three==1), REML = F))$coefficients
```

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## Model failed to converge with max|grad| = 0.00663392 (tol = 0.002, component 1)

    ##              Estimate Std. Error        df   t value     Pr(>|t|)
    ## (Intercept)  28.30487   5.861777  6.258895  4.828719 0.0025949492
    ## years       -19.74025   4.304396 20.003106 -4.586066 0.0001789895

``` r
summary(lmer(gi_new ~ years + (years | id), 
             data = pat %>% filter(four==1), REML = F))$coefficients
```

    ##              Estimate Std. Error        df   t value     Pr(>|t|)
    ## (Intercept) 39.733383   8.786994 11.721029  4.521840 0.0007417537
    ## years       -7.611023   5.904678  7.207885 -1.288982 0.2372285549

``` r
summary(lmer(gi_new ~ years + (years | id), 
             data = pat %>% filter(five==1), REML = F))$coefficients
```

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## Model failed to converge with max|grad| = 0.00342864 (tol = 0.002, component 1)

    ##              Estimate Std. Error       df   t value  Pr(>|t|)
    ## (Intercept)  32.77219  11.534359 5.485370  2.841267 0.0326025
    ## years       -16.66927   8.523312 5.147064 -1.955727 0.1062343

``` r
summary(lmer(gi_new ~ years + (years | id), 
             data = pat %>% filter(six==1), REML = F))$coefficients
```

    ## boundary (singular) fit: see ?isSingular

    ##              Estimate Std. Error        df   t value     Pr(>|t|)
    ## (Intercept)  58.86449  13.093115  5.013105  4.495836 0.0063838859
    ## years       -15.13776   3.249871 11.938430 -4.657957 0.0005603785

``` r
summary(lmer(cbiTotal ~ cluster*years + (years | id), data = pat, REML = F))
```

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## Model failed to converge with max|grad| = 0.00345686 (tol = 0.002, component 1)

    ## Linear mixed model fit by maximum likelihood . t-tests use Satterthwaite's
    ##   method [lmerModLmerTest]
    ## Formula: cbiTotal ~ cluster * years + (years | id)
    ##    Data: pat
    ## 
    ##      AIC      BIC   logLik deviance df.resid 
    ##   2017.8   2074.5   -992.9   1985.8      239 
    ## 
    ## Scaled residuals: 
    ##     Min      1Q  Median      3Q     Max 
    ## -2.1863 -0.5614 -0.0737  0.4785  3.7029 
    ## 
    ## Random effects:
    ##  Groups   Name        Variance Std.Dev. Corr
    ##  id       (Intercept) 222.39   14.913       
    ##           years        14.49    3.807   0.25
    ##  Residual              94.07    9.699       
    ## Number of obs: 255, groups:  id, 36
    ## 
    ## Fixed effects:
    ##                Estimate Std. Error      df t value Pr(>|t|)    
    ## (Intercept)      25.980      6.967  34.956   3.729 0.000679 ***
    ## cluster2          7.269      9.699  34.353   0.749 0.458692    
    ## cluster3         -1.371     11.267  40.685  -0.122 0.903726    
    ## cluster4         13.870      8.714  36.095   1.592 0.120195    
    ## cluster5         -9.131     10.711  33.136  -0.853 0.400032    
    ## cluster6         10.286     10.132  34.236   1.015 0.317126    
    ## years            -3.202      3.496   9.430  -0.916 0.382585    
    ## cluster2:years    4.218      4.995  12.143   0.845 0.414713    
    ## cluster3:years  -14.353     11.168 162.712  -1.285 0.200532    
    ## cluster4:years   -9.457      4.903  10.063  -1.929 0.082435 .  
    ## cluster5:years   -1.545      5.521  13.532  -0.280 0.783823    
    ## cluster6:years    1.094      4.783  10.779   0.229 0.823320    
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
    ## 
    ## Correlation of Fixed Effects:
    ##             (Intr) clstr2 clstr3 clstr4 clstr5 clstr6 years  clst2: clst3:
    ## cluster2    -0.718                                                        
    ## cluster3    -0.618  0.444                                                 
    ## cluster4    -0.799  0.574  0.494                                          
    ## cluster5    -0.650  0.467  0.402  0.520                                   
    ## cluster6    -0.688  0.494  0.425  0.550  0.447                            
    ## years       -0.239  0.172  0.148  0.191  0.156  0.165                     
    ## clustr2:yrs  0.167 -0.222 -0.104 -0.134 -0.109 -0.115 -0.700              
    ## clustr3:yrs  0.075 -0.054 -0.359 -0.060 -0.049 -0.052 -0.313  0.219       
    ## clustr4:yrs  0.171 -0.123 -0.105 -0.245 -0.111 -0.117 -0.713  0.499  0.223
    ## clustr5:yrs  0.152 -0.109 -0.094 -0.121 -0.210 -0.104 -0.633  0.443  0.198
    ## clustr6:yrs  0.175 -0.126 -0.108 -0.140 -0.114 -0.213 -0.731  0.512  0.229
    ##             clst4: clst5:
    ## cluster2                 
    ## cluster3                 
    ## cluster4                 
    ## cluster5                 
    ## cluster6                 
    ## years                    
    ## clustr2:yrs              
    ## clustr3:yrs              
    ## clustr4:yrs              
    ## clustr5:yrs  0.452       
    ## clustr6:yrs  0.521  0.463
    ## convergence code: 0
    ## Model failed to converge with max|grad| = 0.00345686 (tol = 0.002, component 1)

``` r
summary(lmer(cbiTotal ~ years + (years | id), 
             data = pat %>% filter(one==1), REML = F))$coefficients
```

    ## boundary (singular) fit: see ?isSingular

    ##              Estimate Std. Error       df    t value    Pr(>|t|)
    ## (Intercept) 25.792299   5.824493 5.886319  4.4282476 0.004639966
    ## years       -3.203289   4.165912 9.995665 -0.7689285 0.459718480

``` r
summary(lmer(cbiTotal ~ years + (years | id), 
             data = pat %>% filter(two==1), REML = F))$coefficients
```

    ##              Estimate Std. Error       df   t value    Pr(>|t|)
    ## (Intercept) 32.592183   6.416964 6.287147 5.0790659 0.001975285
    ## years        2.171493   5.049118 5.255525 0.4300738 0.684209342

``` r
summary(lmer(cbiTotal ~ years + (years | id), 
             data = pat %>% filter(three==1), REML = F))$coefficients
```

    ## boundary (singular) fit: see ?isSingular

    ##              Estimate Std. Error        df   t value    Pr(>|t|)
    ## (Intercept)  25.12319   6.658568  3.810596  3.773062 0.021346760
    ## years       -20.09366   5.545437 13.425377 -3.623459 0.002944212

``` r
summary(lmer(cbiTotal ~ years + (years | id), 
             data = pat %>% filter(four==1), REML = F))$coefficients
```

    ## boundary (singular) fit: see ?isSingular

    ##              Estimate Std. Error       df   t value     Pr(>|t|)
    ## (Intercept)  39.47174   5.781494 10.96931  6.827256 2.889161e-05
    ## years       -11.90322   3.132077 36.75828 -3.800425 5.257892e-04

``` r
summary(lmer(cbiTotal ~ years + (years | id), 
             data = pat %>% filter(five==1), REML = F))$coefficients
```

    ## boundary (singular) fit: see ?isSingular

    ## Warning: Model failed to converge with 1 negative eigenvalue: -4.0e+00

    ##              Estimate Std. Error df   t value     Pr(>|t|)
    ## (Intercept) 18.702558   3.130472 32  5.974358 1.169069e-06
    ## years       -5.788748   3.559149 32 -1.626442 1.136648e-01

``` r
summary(lmer(cbiTotal ~ years + (years | id), 
             data = pat %>% filter(six==1), REML = F))$coefficients
```

    ## Warning in checkConv(attr(opt, "derivs"), opt$par, ctrl = control$checkConv, :
    ## Model failed to converge with max|grad| = 0.00440242 (tol = 0.002, component 1)

    ##             Estimate Std. Error        df   t value   Pr(>|t|)
    ## (Intercept) 37.34323  10.787126  4.964901  3.461833 0.01821216
    ## years       -2.69536   1.937408 33.589539 -1.391219 0.17330593

``` r
summary(aov(gi_new ~ cluster * years + Error(id), data = pat))
```

    ## 
    ## Error: id
    ##         Df Sum Sq Mean Sq
    ## cluster  1   5171    5171
    ## 
    ## Error: Within
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## cluster         5  29967    5993  11.369 2.06e-10 ***
    ## years           1   4176    4176   7.921  0.00508 ** 
    ## cluster:years   5   8821    1764   3.347  0.00552 ** 
    ## Residuals     508 267802     527                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

``` r
summary(aov(cbiTotal ~ cluster * years + Error(id), data = pat))
```

    ## 
    ## Error: id
    ##         Df Sum Sq Mean Sq
    ## cluster  1   3617    3617
    ## 
    ## Error: Within
    ##                Df Sum Sq Mean Sq F value   Pr(>F)    
    ## cluster         5  14999  2999.8   9.674 1.95e-08 ***
    ## years           1    149   148.7   0.480 0.489305    
    ## cluster:years   5   7347  1469.4   4.739 0.000378 ***
    ## Residuals     242  75041   310.1                     
    ## ---
    ## Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1
