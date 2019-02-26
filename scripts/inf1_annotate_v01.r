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
  rm(psResults)
  for(i in 1:ceiling(nrow(inf1)/100)){
    startrow <- (i-1)*100 + 1
    endrow <- min(nrow(inf1), (i-1)*100 + 100)
    ps <- phenoscanner(snpquery = inf1$psName[startrow:endrow], pvalue = 0.05/nrow(inf1))
    # Convert PS output to sensible classes
    snpStrings <- c("snp","rsid","hg19_coordinates", "hg38_coordinates","a1","a2", "consequence","amino_acids","ensembl","hgnc")
    snpNums <<- c("pos_hg19","pos_hg38","afr","amr","eas","eur","sas","protein_position")
    snpNums <- names(ps$snps)[which(!names(ps$snps) %in% snpStrings)]
    
    # Replace "-" coding with NA to avoid warnings
      dash <- which(levels(ps$snps$protein_position)=="-")
      if(length(dash) > 0){
        levels(ps$snps$protein_position)[dash] <- NA
      }
      ps$snps <- ps$snps %>%
        mutate_at(snpStrings, as.character) %>%
        mutate_at(snpNums, function(x){as.numeric(as.character(x))})
    
    resStrings <- c("snp","rsid","hg19_coordinates", "hg38_coordinates","a1","a2", "trait","efo","study","pmid","ancestry","direction","unit","dataset")
    resNums <- names(ps$results)[which(!names(ps$results) %in% resStrings)]
  
    ps$results <- ps$results %>%
      mutate_at(resStrings, as.character) %>%
      mutate_at(resNums, function(x){as.numeric(as.character(x))})
  
  if(!exists("psResults")){
      psResults <- ps
    } else {
      psResults$snps <- rbind(psResults$snps, ps$snps)
      psResults$results <- rbind(psResults$results, ps$results)
    } 
    rm(ps, startrow, endrow)
  }  

# Merge phenoscanner output with inf
  psResults <- right_join(psResults$snps, psResults$results)
  names(psResults) <- paste0(names(psResults),".ps")
  psResults <- rename(psResults, psName = snp.ps)
  psResults <- full_join(inf1, psResults)
  psResults <- distinct(psResults) # remove duplicated rows, not sure why these appear and too lazy to figure it out right now.
    
# Look up proxies for SNPs with no results found  
#  missRow <- which(inf1$psName %in% psResults2$psName)
#  inf1Proxy <- inf1[missRow,]
#  psProxy <- phenoscanner(inf1Proxy$psName, proxies = "EUR", r2 = 0.8)
    
  psMelt <- psResults %>%
    select(Protein, trait.ps, p.ps) %>%
    mutate(logP = -log10(p.ps)) %>%
    select(-p.ps)
  psCast <- acast(psMelt, Protein ~ trait.ps)  
  
  
  pheatmap(psCast, color = colorRampPalette(c("white","firebrick2"))(200), cex = 0.2)
  
  psResults %>% 
  group_by(Protein) %>%
  summarise(n()) %>%
  data.frame()
  
psResults2 <- psResults %>%
  select(Protein,
         rsid = rsid.ps,
         Chromosome,
         Position,
         Allele1,
         Allele2,
         Freq1,
         Effect,
         StdErr,
         "P-value",
         Direction,
         N,
         a1.ps,
         a2.ps,
         eur.ps,
         consequence.ps,
         hgnc.ps,
         trait.ps,
         pmid.ps,
         ancestry.ps,
         beta.ps,
         se.ps,
         p.ps,
         n.ps,
         n_cases.ps,
         n_controls.ps,
         n_studies.ps,
         unit.ps)
         

fwrite(psResults2, file = "SCALLOP_INF1_phenoscanner_gwas_allResults.csv", na = "NA")
  
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
