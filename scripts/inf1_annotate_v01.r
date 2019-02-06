# Mk. I annotation script. Lookup clumped hits from METAL in phenoscanner.

# Get phenoscanner
#install.packages("devtools")
#library(devtools)
#install_github("phenoscanner/phenoscanner")
#library(phenoscanner)

# Set WD and set up libraries
setwd("C:/Users/Jonathan/Documents/CEU/R_projects/SCALLOP/")
lapply(c("data.table", "dplyr", "stringr", "phenoscanner", "pheatmap"), require, character.only = T)

# Read in INF1 results
clumpedPath <- "Z:/Factors/High_dimensional_genetics/Olink/INF1/INF1.clumped.tbl"
inf1 <- fread(clumpedPath, data.table=F)
#inf1 <- inf1[1:10,] # Truncate for testing 

# Reform columns to remove hyphens and separate data
prot <- str_split_fixed(inf1$Chromosome, ":", 2)
markerName <- str_split_fixed(inf1$MarkerName,"_",2)
inf1 <- inf1 %>%
  mutate(Chromosome = as.numeric(prot[,2]),
         Protein = prot[,1],
         psName = markerName[,1]) %>%
  select(Protein, psName, Chromosome:N)

#Lookup in phenoscanner
infAnnotated <- data.frame()
allTraitsList <- character()
allGWASResults <- list()
#i = 1
for(i in 1:nrow(inf1)){
  infRow <- inf1[i,]
  gwas <- phenoscanner(snpquery = inf1$psName[i], catalogue = "GWAS", pvalue = 5e-8)
  if(nrow(gwas$results) > 0){
    gwas$results$traitName <- paste0(gwas$results$trait," (",gwas$results$pmid,")") 
  } else {
    gwas$results <- data.frame("traitName" = "None", stringsAsFactors = F)
  }
  
  infRow <- infRow %>%
    mutate(rsid = gwas$snps$rsid,
           consequence = gwas$snps$consequence,
           hgnc = gwas$snps$hgnc,
           gwas = paste(unlist(unique(gwas$results$traitName)), collapse = "; "))
  infAnnotated <- bind_rows(infAnnotated, infRow)
  allGWASResults[[i]] <- gwas$results
  allTraitsList <- c(allTraitsList, gwas$results$traitName)
  allTraitsList <- unique(allTraitsList)
  rm(gwas, infRow)
}

# Create matrix of associated traits
traitMat <- matrix(0,nrow = nrow(infAnnotated), 
                   ncol = length(allTraitsList), 
                   dimnames = list(paste0(infAnnotated$Protein, ":",infAnnotated$rsid),allTraitsList))

for(i in 1:nrow(infAnnotated)){
  if(infAnnotated$gwas[i] != "None"){
    colNums <- match(allGWASResults[[i]]$traitName, colnames(traitMat))
    traitMat[i,colNums] <- -log10(as.numeric(allGWASResults[[i]]$p)) 
  }    
}



# Plot heatmap
# png(filename=paste0("plots/urate_olink_correlations_All_corrplot_AnySigPheno_",timestamp,".png"),
#     type="cairo",
#     units="in",
#     width=6,
#     height=14,
#     pointsize=10,
#     res=96)
  pheatmap(traitMat[1:111,], color = colorRampPalette(c("white","firebrick2"))(200), cex = 0.2)
# dev.off()
