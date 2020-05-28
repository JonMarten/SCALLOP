# Script to run phenoscanner on all SNPs that are Olink hits
# Set up and read in data -------------------------------------------------
# Setup
setwd("/home/jm2294/rds/rds-jmmh2-projects/olink_proteomics/scallop/jm2294")
lapply(c("data.table", "dplyr", "stringr", "phenoscanner", "pheatmap", "reshape2"), require, character.only = T)

# Read in INF1 results
inf1 <- fread("INF1.merge", data.table=F)
sigThresh <- 5e-8 / nrow(inf1)

# Reform columns to remove hyphens and separate data
#prot <- str_split_fixed(inf1$Chromosome, ":", 2)
markerName <- str_split_fixed(inf1$MarkerName,"_",2)
inf1 <- inf1 %>%
  mutate(psName = markerName[,1]) 

# Simplified run, running pleiotropic SNPs once only ----------------------
### Simplified query combining SNPs that are pQTLs for multiple proteins

infSmall <- inf1 %>% select(prot,psName) %>% 
  group_by(psName) %>% 
  summarise(numAssocs = n(), prots = paste(prot, collapse = "; ")) %>%
  data.frame()

for(i in 1:ceiling(nrow(infSmall)/100)){
  startrow <- (i-1)*100 + 1
  endrow <- min(nrow(infSmall), (i-1)*100 + 100)
  ps <- phenoscanner(snpquery = infSmall$psName[startrow:endrow], pvalue = sigThresh)
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
  
  if(!exists("psResultsSmall")){
    psResultsSmall <- ps
  } else {
    psResultsSmall$snps <- rbind(psResultsSmall$snps, ps$snps)
    psResultsSmall$results <- rbind(psResultsSmall$results, ps$results)
  } 
  rm(ps, startrow, endrow)
}  