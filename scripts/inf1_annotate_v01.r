# Mk. I annotation script. Lookup clumped hits from METAL in phenoscanner.

# Get phenoscanner
#install.packages("devtools")
#library(devtools)
#install_github("phenoscanner/phenoscanner")
#library(phenoscanner)

# Setup
  setwd("C:/Users/Jonathan/Documents/CEU/R_projects/SCALLOP/")
  lapply(c("data.table", "dplyr", "stringr", "phenoscanner", "pheatmap"), require, character.only = T)

# Read in INF1 results
  clumpedPath <- "Z:/Factors/High_dimensional_genetics/Olink/INF1/INF1.clumped.tbl"
  inf1 <- fread(clumpedPath, data.table=F)

# Reform columns to remove hyphens and separate data
  prot <- str_split_fixed(inf1$Chromosome, ":", 2)
  markerName <- str_split_fixed(inf1$MarkerName,"_",2)
  inf1 <- inf1 %>%
    mutate(Chromosome = as.numeric(prot[,2]),
           Protein = prot[,1],
           psName = markerName[,1]) %>%
    select(Protein, psName, Chromosome:N)

# Submit SNPs to phenoscanner, 100 at a time
  for(i in 1:ceiling(nrow(inf1)/100)){
    startrow <- (i-1)*100 + 1
    endrow <- min(nrow(inf1), (i-1)*100 + 100)
    ps <- phenoscanner(snpquery = inf1$psName[startrow:endrow])
    
    if(!exists(psResults)){
      psResults <- ps
    } else {
      psResults$snps <- rbind(psResults$snps, ps$snps)
      psResults$results <- rbind(psResults$results, ps$results)
    }  
  }  
    
### Old code
    
#Lookup in phenoscanner, 1 at a time
  # Initialise empty data frames
    infAnnotated <- data.frame()
    allTraitsList <- character()
    allGWASResults <- list()
  # Loop over all SNPs
    for(i in 1:nrow(inf1)){
      infRow <- inf1[i,]
      gwas <- phenoscanner(snpquery = inf1$psName[i], 
                           catalogue = "GWAS", 
                           pvalue = (5e-8/nrow(infRow)))
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

traitMat[which(is.infinite(traitMat))] <- 300 # deal with infinite values from P = 0
traitMat <- traitMat[,-which(colnames(traitMat)=="None")] # remove SNPs with no GWAS hits

traitMatOut <- data.frame(traitMat)

fwrite(traitMatOut, file = "SCALLOP_INF1_phenoscanner_gwas_matrix.csv", sep = ",", row.names = T)

# Plot heatmap
 png(filename="SCALLOP_INF1_phenoscanner_gwas_heatmap.png",
     type="cairo",
     units="px",
     width=10000,
     height=10000,
     pointsize=10,
     res=400)
  pheatmap(traitMat, color = colorRampPalette(c("white","firebrick2"))(200), cex = 0.2)
 # dev.off()

# Heatmap only with phenos associated with >5 SNPs
counts <- apply(traitMat,MARGIN = 2, FUN = function(x){length(which(x!=0))})
traitMatFilt <- traitMat[,which(counts > 10)]
zeros <- apply(traitMatFilt,MARGIN = 1, FUN = sum) # remove SNPs with no phenos left
traitMatFilt <- traitMatFilt[which(zeros != 0),]
pdf("SCALLOP_INF1_phenoscanner_gwas_heatmap.pdf", width=10, height=10)
pheatmap(t(traitMatFilt), color = colorRampPalette(c("white","firebrick2"))(200), cex = 0.8)
dev.off()
