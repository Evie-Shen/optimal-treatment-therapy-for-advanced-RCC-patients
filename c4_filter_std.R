
rm(list=ls())
require(readxl)
library(plyr)
library(tidyverse)
library(ggfortify)
library(GGally)
library(gridExtra)
library(tidyr)
library(dplyr)
library(survival)
library(survminer)
require(gridExtra)
library(lmtest)

setwd("J:\\projects\\Research\\JAVELIN")

load("data\\anadat.RDa")
load("data\\varTable.RDa")

## the TCGA_cluster variable is not significant, so I removed it, it is a multicategory variable
dset <- dset %>% 
  mutate(SEXN = ifelse(is.na(SEX), NA, ifelse(SEX=="M", 0, 1)),
         PDL1N = ifelse(is.na(PDL1FL), NA, ifelse(PDL1FL=="N", 0, 1)),
         PFS_event = 1- as.numeric(PFS_P_CNSR)) %>%
  select(-c("PFS_P_CNSR","TRT01P","SEX","PDL1FL","TCGA_cluster"))

#### combine all pvalues together  #######

load("output\\singleP\\univariate_rna_bytrt.RDa")
load("output\\singleP\\univariate_rna_Int.RDa")
load("output\\singleP\\univariate_mut_bytrt.RDa")
load("output\\singleP\\univariate_mut_Int.RDa")
load("output\\singleP\\univariate_main_bytrt.RDa")
load("output\\singleP\\univariate_main_Int.RDa")

## for the by treatment pvalues, change from long to wide, the column name will add suffix of _0, _1
## then merge with int pvalue
main_pval <- univariate_main_bytrt %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(N, missPct, mean_value, sd_value, skewness_value, median_value, PFS_pvalue, PFS_beta)
  ) %>%
  inner_join(univariate_main_Int, by="biomarker") %>%
  inner_join(varTable %>% rename(biomarker=Variable), by="biomarker")

## for the by treatment pvalues, change from long to wide, the column name will add suffix of _0, _1
## then merge with int pvalue
rna_pval <- univariate_rna_bytrt %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(N, missPct, mean_value, sd_value, skewness_value, median_value, PFS_pvalue, PFS_beta)
  ) %>%
  inner_join(univariate_rna_Int, by="biomarker") %>%
  inner_join(varTable %>% filter(Source=="rna") %>% rename(biomarker=Variable), by="biomarker")

## for the by treatment pvalues, change from long to wide, the column name will add suffix of _0, _1
## then merge with int pvalue
mut_pval <- univariate_mut_bytrt %>%
  pivot_wider(
    names_from = Treatment,
    values_from = c(N, missPct, mean_value, sd_value, skewness_value, median_value, PFS_pvalue, PFS_beta)
  ) %>%
  inner_join(univariate_mut_Int, by="biomarker") %>%
  inner_join(varTable %>% filter(Source=="mut") %>% rename(biomarker=Variable), by="biomarker")

## combine results of all three biomarker tables (main, rna and mutation)
pval_all <- rbind(main_pval, rna_pval, mut_pval) #21834
## since mut and rna has the same gene name, to differentiate them, add "_mut" after mut's gene name
pval_all <- pval_all %>%
  mutate(biomarker_2 = ifelse(Source == "mut", paste0(biomarker,"_mut"), biomarker))
table(pval_all$Source) 

## filter out rna with missing percentage greater than 95% 
pval_all_2 <- pval_all %>%
  filter(missPct<=0.95 & missPct_0<=0.95 & missPct_1 <=0.95) # 21683
table(pval_all_2$Source) 


## select the biomarkers with p-values less than 0.05
pval_sig05 <- pval_all_2 %>%
  filter(PFS_pvalue_0 <= 0.05 | PFS_pvalue_1 <= 0.05 | PFS_int_Pval <= 0.05) #8331
table(pval_sig05$Source)

save(pval_all, pval_sig05, file = "output\\singleP\\pval_all_sig05.RDa")

###########################################################
### standardize data for significant BMx as well as age and sex

load("output\\singleP\\pval_all_sig05.RDa")
load("data\\rna_wide.RDa")
load("data\\mut_wide.RDa")

## from significant gene list, select the biomarkers not from rna and mut
main_gene_list <- pval_sig05 %>% 
  filter(!(Source %in% c("rna","mut"))) %>%
  select(biomarker)
## generate a dataset with the above selected biomarker as well as age and sex
dset_sub <- dset %>%
  select(ID, SEXN, AGE, main_gene_list$biomarker)

## from significant gene list, select the biomarkers from rna category and get their data
rna_gene_list <- pval_sig05 %>% 
  filter(Source == "rna") %>%
  select(biomarker)
rna_t = rna_v3 %>% 
  select(-HUGO) %>% 
  t() %>% 
  as.data.frame 
names(rna_t) <- rna_v3$HUGO
rna_t$ID = rownames(rna_t)

rna_sub <- rna_t %>%
  select(ID, rna_gene_list$biomarker)

## join rna data with the other biomarker in dset
## remove one biomarker CD8_INVASIVE_MARGIN_SURFACE_AREA since  
## more than half of population have missing values
dset_rna <- dset_sub %>%
  inner_join(rna_sub, by = "ID") %>%
  select(-"CD8_INVASIVE_MARGIN_SURFACE_AREA")

## found if the biomarker is multi-level (i.e. continuous) or bi-level
## if it is bi-level, i.e. 0 vs. 1, then we do not standardize it
## if it is multi-level, then we need standardize it
bmx_level <- apply(dset_rna, 2, function(x) return(sum(!is.na(unique(x)))))
bmx_level <- as.data.frame(bmx_level)
bmx_level$ID <- rownames(bmx_level)

## select the columns (i.e. biomarkers) which is multi-level
dset_rna_mlevel <- dset_rna %>%
  select(bmx_level[bmx_level$bmx_level>2,]$ID)

## select the columns (i.e. biomarkers) which is bi-level
dset_rna_2level <- dset_rna %>%
  select(ID, bmx_level[bmx_level$bmx_level==2,]$ID)

## for biomarker data with multi-level, standardize it using scale
dset_rna_mlevel_2 <- cbind(dset_rna_mlevel[,1],as.data.frame(scale(dset_rna_mlevel[,-1])))

## combine two-level biomarker data set with the standardized multi-level biomarker data
sig_data_std <-  dset_rna_2level %>%
  inner_join(dset_rna_mlevel_2, by="ID")

####### add mutational data

## from significant gene list, select the biomarkers from mut category and get their data
mut_gene_list <- pval_sig05 %>%
  filter(Source == "mut") %>%
  select(biomarker, biomarker_2)

mut_t = mut_v3 %>% 
  select(-HUGO) %>% 
  t() %>% 
  as.data.frame 
names(mut_t) <- mut_v3$HUGO
mut_t$ID = rownames(mut_t)

mut_sub <- mut_t %>%
  select(ID, mut_gene_list$biomarker) 
## add "_mut" to the column names since the mut gene and rna gene sometimes have the same name
colnames(mut_sub)[2:ncol(mut_sub)] <- sapply(colnames(mut_sub)[2:ncol(mut_sub)], function(x)paste0(x,"_mut"))

### add the mut_data with the previous combined data after standardization
sig_data_std_2 <-  sig_data_std %>%
  inner_join(mut_sub, by="ID") 

### combined the standadized biomarker data with the clinical outcome data
sig_data_std_3 <- dset %>%  ## total 691 subjects
  select(ID, PFS_P, PFS_event, TRT01N) %>%
  inner_join(sig_data_std_2, b="ID")


### calculate the biomarker missingness for each subjects
## total 691 subjects, 683 does not have any missing value
tmp <- apply(sig_data_std_3, 1, function(x)sum(is.na(x)))
missing_Sub <- as.data.frame(cbind(ID=sig_data_std_3$ID, missN=tmp))
missing_Sub$missN <- as.numeric(missing_Sub$missN)
table(missing_Sub$missN) 
# 0   4  24 
# 683   2   6 

### find the biomarkers with missing value to see if we can impute those biomarkers 
### some mut DNA has value 0 or 1, if we impute, we may imputed continuous value instead of 0 or 1
tmp <- apply(sig_data_std_3, 2, function(x)sum(is.na(x)))
missing_Var <- as.data.frame(cbind(Variable=names(sig_data_std_3), missN=tmp))
missing_Var$missN <- as.numeric(missing_Var$missN)
table(missing_Var$missN)
# 0    4   24 
# 8240   24   72 

####################################################

dset_std4ana <- sig_data_std_3 %>% 
  inner_join(missing_Sub %>% filter(missN==0), by="ID") %>%
  select(-"missN")

### save the standardize data
save(dset_std4ana, pval_sig05, file="output\\anaDset\\dset_std4ana.RDa")
