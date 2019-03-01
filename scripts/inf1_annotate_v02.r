# Mk. II annotation script. Retains only the 'small' lookup from Mk. I, where each SNP is looked up only once regardless of the number of proteins it associates with.
# Set up and read in data -------------------------------------------------
# Setup
  setwd("C:/Users/Jonathan/Documents/CEU/R_projects/SCALLOP/")
  lapply(c("data.table", "dplyr", "stringr", "phenoscanner", "pheatmap", "reshape2"), require, character.only = T)
  

# Read in INF1 results
  clumpedPath <- "Z:/Factors/High_dimensional_genetics/Olink/INF1/INF1.clumped.tbl"
  inf1 <- fread(clumpedPath, data.table=F)
  sigThresh <- 5e-8 / nrow(inf1)

# Reform columns to remove hyphens and separate data
  prot <- str_split_fixed(inf1$Chromosome, ":", 2)
  markerName <- str_split_fixed(inf1$MarkerName,"_",2)
  inf1 <- inf1 %>%
    mutate(Chromosome = as.numeric(prot[,2]),
           Protein = prot[,1],
           psName = markerName[,1]) %>%
    select(Protein, psName, Chromosome:N)
  rm(prot, markerName, clumpedPath)

# Simplified run, running pleiotropic SNPs once only ----------------------
### Simplified query combining SNPs that are pQTLs for multiple proteins

infSmall <- inf1 %>% select(Protein,psName) %>% 
  group_by(psName) %>% 
  summarise(numAssocs = n(), prots = paste(Protein, collapse = "; ")) %>%
  data.frame()

rm(psResultsSmall)
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

# Cpnvert phenoscanner output to data frame instead of list
  psResultsSmall <- right_join(psResultsSmall$snps, psResultsSmall$results)
  
# Where a SNP is assocaited with a trait more than once, pick EA first, then largest study, or if no n is available, pick smallest P-value
  psResultsSmall <- psResultsSmall %>%
    mutate(snpTrait = paste0(snp, trait))
  
  dupes <- duplicated(psResultsSmall$snpTrait)
  dupePhenos <- unique(psResultsSmall$snpTrait[dupes])
  #dupeRows <- which(psResultsSmall$snpTrait %in% dupePhenos)
  
  dupesFiltered <- data.frame()
  for(i in 1:length(dupePhenos)){
    subset <- psResultsSmall %>% filter(snpTrait == dupePhenos[i])
    eaRows <- which(subset$ancestry == "European")
    if(length(eaRows) == 1){
      keeprow <- subset[eaRows,]
    } else {
      if (length(eaRows) > 1){
        subset <- subset[eaRows,]
        keeprow <- subset[which.max(na.exclude(subset$n)),]
        if(length(na.exclude(subset$n))>0) {
          keeprow <- subset[which.max(na.exclude(subset$n)),]
        } else {
          keeprow <- subset[which.min(na.exclude(subset$p)),]
        }
      } else {
        if(length(na.exclude(subset$n))>0) {
          keeprow <- subset[which.max(na.exclude(subset$n)),]
        } else {
          keeprow <- subset[which.min(na.exclude(subset$p)),]
        }
      }
    }
    dupesFiltered <- rbind(dupesFiltered, keeprow)
    rm(keeprow, subset, eaRows)
  }
  
psResultsSmall2 <- psResultsSmall #backup   
psResultsSmall <- psResultsSmall %>%
  filter(!snpTrait %in% dupePhenos)
psResultsSmall <- rbind(psResultsSmall, dupesFiltered)
    
# Merge phenoscanner output with infSmall  
  names(psResultsSmall) <- paste0(names(psResultsSmall),".ps")
  psResultsSmall <- rename(psResultsSmall, psName = snp.ps)
  psResultsSmall <- full_join(infSmall, psResultsSmall)
  psResultsSmall <- distinct(psResultsSmall) # remove duplicated rows, not sure why these appear and too lazy to figure it out right now.

# Look up proxies for SNPs with no results found  
#  missRow <- which(inf1$psName %in% psResults2$psName)
#  inf1Proxy <- inf1[missRow,]
#  psProxy <- phenoscanner(inf1Proxy$psName, proxies = "EUR", r2 = 0.8)

# Generate "RSID (protein)" identiers for each association
  noRS <- which(is.na(psResultsSmall$rsid.ps))
  psResultsSmall$rsid.ps[noRS] <- psResultsSmall$psName[noRS]
  psResultsSmall <- psResultsSmall %>%
    mutate(snpProt = paste0(rsid.ps," (gene: ",hgnc.ps,", Olink: ", prots,")"))

# List all SNPs and their associated traits  
  snpProtSummary <- psResultsSmall %>% 
    group_by(snpProt) %>%
    summarise(n(),
              traits = paste(unique(trait.ps), collapse = "; ")) %>% 
    data.frame
# List all traits and their proteins 
  traitSummary <- psResultsSmall %>% 
    group_by(trait.ps) %>%
    summarise(n(),
              Proteins = paste(unique(prots), collapse = "; ")) %>% 
    data.frame()

# Convert to a contingency table for easy pheatmapping
psMeltSmall <- psResultsSmall %>%
  select(snpProt, trait.ps, p.ps) %>%
  mutate(logP = -log10(p.ps),
         binP = ifelse(p.ps < sigThresh, 1, 0)) 

# Filter to only traits associated with >2 SNPs
psMeltSmallFilt <- psMeltSmall %>%
  filter(trait.ps %in% names(which(table(psMeltSmall$trait.ps) > 2))) %>%
  filter(!is.na(logP)) %>%  
  distinct()

  
psCastSmallFilt <- acast(psMeltSmallFilt, 
                         snpProt ~ trait.ps,
                         value.var = "binP",
                         fun.aggregate = function(x){x[1]}) # why max doesn't work here I have no idea, but it just returns Inf every time. 

psCastSmallFilt[which(is.na(psCastSmallFilt), arr.ind = T)] <- 0
psCastSmallFilt[which(is.infinite(psCastSmallFilt), arr.ind = T)] <- 1

pdf(file = "Y:/Projects/SCALLOP/INF1/SCALLOP_INF1_phenoscanner_gwas_small_heatmap_2plus_binary.pdf", height = 20, width = 20)
  pheatmap(t(psCastSmallFilt), color = colorRampPalette(c("white","firebrick2"))(200), cex = 1,border_color = "gray50")
dev.off()

write.csv(psResultsSmall, row.names=F, file = "Z:/Factors/High_dimensional_genetics/Olink/INF1/phenoscanner/SCALLOP_INF1_phenoscanner_gwas_annotations.csv")