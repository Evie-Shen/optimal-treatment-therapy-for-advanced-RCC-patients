remove(list=ls())
library(tidyverse)
library(haven)
library(readxl)
library(stringr)
library(ggplot2)
library(pROC)
library(plyr)
library(tidyverse)
library(pspearman)
library(GGally)
library(ggfortify)
library(survival)
library(survminer)
library(grid)
library(gridExtra)

datafile <- "J:\\projects\\Research\\JAVELIN\\data\\41591_2020_1044_MOESM3_ESM.xlsx"

## read in each sheet of the datafile
clnc <- read_xlsx(datafile, sheet="S11_Clinical_data", skip = 1, col_names = TRUE, na=c("","NA")) 
ihc <- read_xlsx(datafile, sheet="S12_CD8_IHC_data", skip = 1, col_names = TRUE, na=c("","NA"))
immu <- read_xlsx(datafile, sheet="S14_ImmuneNet_decon_data", skip = 1, col_names = TRUE, na=c("","NA")) 
pathway <- read_xlsx(datafile, sheet="S16_Pathway_scores", skip = 1, col_names = TRUE, na=c("","NA")) 
wesSum <- read_xlsx(datafile, sheet="S18_WES_summary_stats", skip = 1, col_names = TRUE, na=c("","NA")) 
math <- read_xlsx(datafile, sheet="S19_MATH_scores", skip = 1, col_names = TRUE, na=c("","NA")) 
hla <- read_xlsx(datafile, sheet="S17_HLA_subtypes", skip = 1, col_names = TRUE, na=c("","NA")) 
rna <- read_xlsx(datafile, sheet="S13_Gene_expression_TPM", skip = 1, col_names = TRUE, na=c("","NA"))
mut <- read_xlsx(datafile, sheet="S21_Gene_mutation_status", skip = 1, col_names = TRUE, na=c("","NA"))
fcgr <- read_xlsx(datafile, sheet="S20_FCGR_mutation_data", skip = 1, col_names = TRUE, na=c("","NA"))

## add "M" in front of the gene name starting with number in rna
HUGO.num <- sort(rna$HUGO)[1:28]
rna$myrowname = row.names(rna)
rna_v2 <- rna %>% 
  mutate(HUGO2 = ifelse(HUGO %in% HUGO.num, paste0("M", HUGO, ".",myrowname), HUGO),
         HUGO = HUGO2) %>%
  select(-c("myrowname","HUGO2"))

## add "M" in front of the gene name starting with number in mut
HUGO.num <- sort(mut$HUGO)[1:27]
mut$myrowname = row.names(mut)
mut_v2 <- mut %>% 
  mutate(HUGO2 = ifelse(HUGO %in% HUGO.num, paste0("M", HUGO, ".",myrowname), HUGO),
         HUGO = HUGO2) %>%
  select(-c("myrowname","HUGO2"))

## remove mut genes with percentage of mutation "1" less than 4%
mut_percent <- cbind(mut_v2[,1], mut_num=apply(mut_v2[,-1],1,sum))
tmp <- mut_percent %>% filter(mut_num > 738*0.04) #2613 genes, >=30
tmp <- mut_percent %>% filter(mut_num >= 30) #2613 genes, >=30
mut_v3 <- mut_v2 %>% filter(HUGO %in% tmp$HUGO)

## remove rna genes with missing percentage greater than 95%
rna_percent <- cbind(rna_v2[,1], rna_num=apply(rna_v2[,-1],1,function(x)sum(x==0.01)))
tmp <- rna_percent %>% filter(rna_num < 727*0.95) #3992 genes, >=690
rna_v3 <- rna_v2 %>% filter(HUGO %in% tmp$HUGO)

names(immu)[1] <- "ID"
names(pathway)[1] <- "ID"
names(math)[1] <- "ID"

## merge all biomarkers with ID
dmain <- clnc %>% 
  mutate(TRT01N = ifelse(is.na(TRT01P), NA, ifelse(TRT01P=="Sunitinib", 0 ,1))) %>%
  full_join(ihc, by="ID") %>%
  full_join(immu, by="ID") %>%
  full_join(pathway, by="ID") %>%
  full_join(wesSum, by="ID") %>%
  full_join(math, by="ID") 

dset <- dmain %>%
  full_join(hla, by="ID")

## for column names with special character, replace them with "_"
colnames(dmain) = sapply(colnames(dmain), function(x) gsub("[*:-]", "_", x))
colnames(dmain) = sapply(colnames(dmain), function(x) gsub("\\(", "_", x))
colnames(dmain) = sapply(colnames(dmain), function(x) gsub("\\)", "_", x))

colnames(dset) = sapply(colnames(dset), function(x) gsub("[*:-]", "_", x))
colnames(dset) = sapply(colnames(dset), function(x) gsub("\\(", "_", x))
colnames(dset) = sapply(colnames(dset), function(x) gsub("\\)", "_", x))

rna_v2$HUGO = sapply(rna_v2$HUGO, function(x) gsub("[*:-]", "_", x))
rna_v2$HUGO = sapply(rna_v2$HUGO, function(x) gsub("\\(", "_", x))
rna_v2$HUGO = sapply(rna_v2$HUGO, function(x) gsub("\\)", "_", x))

rna_v3$HUGO = sapply(rna_v3$HUGO, function(x) gsub("[*:-]", "_", x))
rna_v3$HUGO = sapply(rna_v3$HUGO, function(x) gsub("\\(", "_", x))
rna_v3$HUGO = sapply(rna_v3$HUGO, function(x) gsub("\\)", "_", x))

mut_v2$HUGO = sapply(mut_v2$HUGO, function(x) gsub("[*:-]", "_", x))
mut_v2$HUGO = sapply(mut_v2$HUGO, function(x) gsub("\\(", "_", x))
mut_v2$HUGO = sapply(mut_v2$HUGO, function(x) gsub("\\)", "_", x))

mut_v3$HUGO = sapply(mut_v3$HUGO, function(x) gsub("[*:-]", "_", x))
mut_v3$HUGO = sapply(mut_v3$HUGO, function(x) gsub("\\(", "_", x))
mut_v3$HUGO = sapply(mut_v3$HUGO, function(x) gsub("\\)", "_", x))

save(dset, dmain, rna_v3, mut_v3, file="J:\\projects\\Research\\JAVELIN\\data\\anadat.RDa")

## generate a table to list all biomarker's name and their category
varTable <- rbind(
  cbind(Variable=names(clnc)[-1],Source="clnc"),
  cbind(Variable=names(ihc)[-1],Source="ihc"),
  cbind(Variable=names(immu)[-1],Source="immuNet"),
  cbind(Variable=names(pathway)[-1],Source="pathway"),
  cbind(Variable=names(wesSum)[-1],Source="wesSum"),
  cbind(Variable=names(math)[-1],Source="math"),
  cbind(Variable=names(hla)[-1],Source="hla"),
  cbind(Variable=rna_v2$HUGO,Source="rna"),
  cbind(Variable=mut_v2$HUGO,Source="mut"),
  cbind(Variable=names(fcgr)[-1],Source="fcgr")
  )
varTable <- as.data.frame(varTable) # 41464

## SEXN, PDL1N was not in the original data, add them into VarTable
varTable <- rbind(varTable, c("SEXN","clnc",NA), c("PDL1N","clnc",NA)) # 41466
varTable <- unique(varTable) # 41466

## for biomarker names with special character, replace them with "_"
varTable$Variable = sapply(varTable$Variable, function(x) gsub("[*:-]", "_", x))
varTable$Variable = sapply(varTable$Variable, function(x) gsub("\\(", "_", x))
varTable$Variable = sapply(varTable$Variable, function(x) gsub("\\)", "_", x))

varTable <- varTable %>%
  mutate(flag = ifelse((Variable %in% rna_v3$HUGO) & Source=="rna", "rna_v3",
                       ifelse((Variable %in% mut_v3$HUGO) & Source=="mut", "mut_v3", NA)))

varTable %>% group_by(Source) %>% dplyr::summarise(n=n())
varTable %>% group_by(flag) %>% dplyr::summarise(n=n())

save(varTable, file="J:\\projects\\Research\\JAVELIN\\data\\varTable.RDa")
