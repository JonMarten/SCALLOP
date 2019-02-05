# Mk. I annotation script

# Set WD and set up libraries
setwd("C:/Users/Jonathan/Documents/CEU/R_projects/SCALLOP/")
lapply(c("data.table", "dplyr", "stringr"), require, character.only = T)

# Read in INF1 results
clumpedPath <- "Z:/Factors/High_dimensional_genetics/Olink/INF1/INF1.clumped.tbl"
inf1 <- fread(clumpedPath, data.table=F)

# Reform columns to remove hyphens and separate data
prot <- str_split_fixed(inf1$Chromosome, ":", 2)
inf1$Chromosome <- as.numeric(prot[,2])
inf1$Protein <- prot[,1]
markerName <- str_split_fixed(inf1$MarkerName,"_",2)
inf1$psName <- markerName[,1]
