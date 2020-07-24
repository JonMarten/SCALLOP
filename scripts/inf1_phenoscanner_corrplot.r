# Mk. III annotation script. Designed to run on CSD3 with final list of sentinel SNPs
# Set up and read in data -------------------------------------------------
# Setup
setwd("/home/jm2294/rds/rds-jmmh2-projects/olink_proteomics/scallop/jm2294")
lapply(c("data.table", "dplyr", "stringr", "phenoscanner", "pheatmap", "reshape2"), require, character.only = T)

# Read in INF1 conditionally independent results and merge to re-add betas
inf1 <- fread("INF1.merge", data.table=F)
fullRes <- fread("INF1.tbl", data.table = F)
pr <- str_split_fixed(fullRes$Chromosome, ":", 2)
fullRes <- fullRes %>%
  mutate(prot = pr[,1],
         Chr = as.numeric(pr[,2])) %>%
  mutate(SNPprot = paste0(prot, ":", MarkerName))
inf1 <- inf1 %>%
  mutate(SNPprot = paste0(prot, ":", MarkerName))
inf1 <- left_join(inf1, fullRes)

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

# Convert phenoscanner output to data frame instead of list
psResultsSmall <- right_join(psResultsSmall$snps, psResultsSmall$results)

# Where a SNP is assocaited with a trait more than once, pick EA first, then largest study, or if no n is available, pick smallest P-value. Remove any SNPs with no beta values in phenoscanner
psResultsSmall <- psResultsSmall %>%
  filter(!is.na(beta)) %>%
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

psResultsSmall <- psResultsSmall %>%
  filter(!snpTrait %in% dupePhenos)
psResultsSmall <- rbind(psResultsSmall, dupesFiltered)

# Merge phenoscanner output with inf1
inf1b <- inf1 %>% 
  select(psName, prots = prot, Allele1, Allele2, Effect, StdErr)

names(psResultsSmall) <- paste0(names(psResultsSmall),".ps")
psResultsSmall <- rename(psResultsSmall, psName = snp.ps)
psResultsSmall <- full_join(inf1b, psResultsSmall)
psResultsSmall <- distinct(psResultsSmall) # remove duplicated rows, not sure why these appear and too lazy to figure it out right now.

# check for beta direction match
psResultsSmall <- psResultsSmall %>%
  mutate(alleleMismatch = ifelse(toupper(a1.ps) == toupper(Allele1), 0, 1)) %>%
  mutate(betaFlip.ps = ifelse(alleleMismatch == 1, -beta.ps, beta.ps)) %>%
  mutate(betaMatch = ifelse(is.na(betaFlip.ps), 999, ifelse(sign(betaFlip.ps) == sign(Effect), 1, -1)))
  
# Look up proxies for SNPs with no results found  
#  missRow <- which(inf1$psName %in% psResults2$psName)
#  inf1Proxy <- inf1[missRow,]
#  psProxy <- phenoscanner(inf1Proxy$psName, proxies = "EUR", r2 = 0.8)

# Generate "RSID (protein)" identiers for each association
noRS <- which(is.na(psResultsSmall$rsid.ps))
psResultsSmall$rsid.ps[noRS] <- psResultsSmall$psName[noRS]
psResultsSmall <- psResultsSmall %>%
  mutate(snpProt = paste0(rsid.ps," (",hgnc.ps,"): ", prots))

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
fwrite(traitSummary, file = "trait_list.csv")

## group by EFO where possible
noEFO <- psResultsSmall %>% filter(efo.ps == "-") 
EFO <- psResultsSmall %>% filter(efo.ps != "-") 

noEfoSummary <- noEFO %>% 
  group_by(trait.ps) %>%
  summarise(n(),
            Traits = paste(unique(trait.ps), collapse = "; "),
            Proteins = paste(unique(prots), collapse = "; ")) %>% 
  data.frame()

efoSummary <- EFO %>% 
  group_by(efo.ps) %>%
  summarise(n(),
            Traits = paste(unique(trait.ps), collapse = "; "),
            Proteins = paste(unique(prots), collapse = "; ")) %>% 
  data.frame()
names(efoSummary) <- c("EFO","traits","proteins")
names(noEfoSummary) <- c("EFO","traits","proteins")
allSummary <- rbind(efoSummary, noEfoSummary)
fwrite(allSummary, file = "efo_list.csv")

#annotate with disease status column to pull out disease traits
  tr <- fread("efo_list_annotated.csv", data.table = F)
  psResultsSmall2 <- psResultsSmall %>%
    mutate(traitPlotName = ifelse(efo.ps == "-", 
                                  trait.ps, 
                                  tr$traitName[match(efo.ps, tr$EFO)])
           )
  diseases <- tr$traitName[tr$clinical == 1]
  inflam <- tr$traitName[tr$immune_mediated == 1]
  infdis <-  tr$traitName[tr$immune_mediated == 1 & tr$clinical == 1]
  
# Convert to a contingency table for easy pheatmapping
psMeltSmall <- psResultsSmall2 %>%
  select(snpProt, traitPlotName, p.ps, betaMatch) %>%
  mutate(logP = -log10(p.ps),
         binP = ifelse(p.ps < sigThresh, 1, 0)) %>%
  filter(traitPlotName %in% infdis)



# Filter to only traits associated with >0 SNPs
psMeltSmallFilt <- psMeltSmall %>%
  filter(traitPlotName %in% names(which(table(psMeltSmall$traitPlotName) > 0))) %>%
  filter(!is.na(logP)) %>%  
  distinct()

psCastSmallFilt <- acast(psMeltSmallFilt, 
                         snpProt ~ traitPlotName,
                         value.var = "betaMatch",
                         fun.aggregate = function(x){x[1]}) 

#psCastSmallFilt <- acast(psMeltSmallFilt, 
#                         snpProt ~ traitPlotName,
#                         value.var = "betaMatch",
#                         fun.aggregate = function(x){paste0(unique(x),collapse = ",")}) 

# why max doesn't work here I have no idea, but it just returns Inf every time. 

psCastSmallFilt[which(is.na(psCastSmallFilt), arr.ind = T)] <- 0
psCastSmallFilt[which(is.infinite(psCastSmallFilt), arr.ind = T)] <- 1

# To refine heatmap, output matrix as CSV and manually combine columns that tag multiple phenotypes already included, update plot names, then reimport. Would be nice to automate but probably too time consuming for too little gain
write.table(psCastSmallFilt, file = "heatmap_plot_matrix_toedit.csv", row.names = T, col.names = T, sep = ",")

plotdf <- fread("heatmap_plot_matrix_edited.csv", data.table = F)
plotmat <- as.matrix(plotdf[,-1])
rownames(plotmat) <- plotdf[,1]

pheatmap(mat = plotmat, 
         color = colorRampPalette(c("#4287f5","#ffffff","#e32222"))(3), 
         legend = T,
         main = "Olink pQTLs overlapping with QTLs for immune-mediated clinical outcomes",
         angle_col = "45", 
         filename = "SCALLOP_INF1_pQTL_immune-mediated_clinical_qtl_heatmap.png",
         width = 16,
         height = 10,
         treeheight_row = 100,
         treeheigh_col = 100,
         cellheight = 20,
         cellwidth = 20)

pheatmap(mat = plotmat, 
         color = colorRampPalette(c("#4287f5","#ffffff","#e32222"))(3), 
         legend = T,
         main = "Olink pQTLs overlapping with QTLs for immune-mediated clinical outcomes",
         angle_col = "45", 
         filename = "SCALLOP_INF1_pQTL_immune-mediated_clinical_qtl_heatmap_unclustered.png",
         width = 16,
         height = 10,
         cluster_rows = F,
         cluster_cols = F,
         cellheight = 20,
         cellwidth = 20)


write.csv(psResultsSmall, row.names=F, file = "SCALLOP_INF1_phenoscanner_gwas_annotations.csv")
# psResultsSmall<- fread("SCALLOP_INF1_phenoscanner_gwas_annotations.csv", data.table = F)

# Make filtered results file for paper
resTab <- psResultsSmall2 %>%
  filter(paste0(snpProt, traitPlotName) %in% paste0(psMeltSmallFilt$snpProt, psMeltSmallFilt$traitPlotName))
 # select("RSID" = rsid.ps,
 #        "chr" = chr.ps,
 #        "pos_hg19" = pos_hg19.ps,
 #        "A1" = Allele1,
 #        "A2" = Allele2,
 #        "protein" = prots,
 #        "Beta_olink" = Effect,
 #        "SE_olink" = StdErr,
 #        "Beta_phenoscanner" = betaFlip.ps,
 #        "Consequence_phenoscanner" = consequence.ps,
 #        
         
         

write.csv(resTab, quote = T, row.names = F, file = "SCALLOP_INF1_phenoscanner_gwas_filtered_associations.csv")