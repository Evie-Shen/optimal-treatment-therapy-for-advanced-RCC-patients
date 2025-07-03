
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

## the TCGA_cluster variable is not significant, so I removed it, it is a multicategory variable
dset <- dset %>% 
  mutate(SEXN = ifelse(is.na(SEX), NA, ifelse(SEX=="M", 0, 1)),
         PDL1N = ifelse(is.na(PDL1FL), NA, ifelse(PDL1FL=="N", 0, 1)),
         PFS_event = 1- as.numeric(PFS_P_CNSR)) %>%
  select(-c("PFS_P_CNSR","TRT01P","SEX","PDL1FL","TCGA_cluster"))

############## Function for UniVariate Analysis ########################
uni_analysis_byarm = function(dat) {
  output_list = list()
  ## get the association p-value for each arm separately
  for (trt in c(0,1)) {
    print(trt)
    
    ## get the data for one specific arm and split by biomarker
    tempdata_sub = dat %>% 
      filter(TRT01N==trt) %>% 
      group_split(biomarker)
    
    ## another function for regression analysis for a biomarker
    regression <- function(input_data) {
      biomarker <- input_data$biomarker[1]
      if (!is.na(biomarker)){
        
        input_data <- input_data[complete.cases(input_data),]
        
        n = nrow(input_data)
        m = sum(input_data[["value"]]==0.01)
        misspct = m/n
        mean_value = mean(input_data[["value"]], na.rm = T)
        sd_value = sd(input_data[["value"]], na.rm = T)
        skewness_value <- (n * sum((input_data[["value"]] - mean_value)^3)) / ((n - 1) * (n - 2) * sd_value^3)
        median_value = median(input_data[["value"]], na.rm=T)
        
        PFS_f = formula(paste0("Surv(","PFS_P",",","PFS_event",") ~ ", " value + AGE + SEXN"))
        PFS_fit = try(coxph(formula = PFS_f, data = input_data), silent = TRUE)
        if (class(PFS_fit)[1] == "try-error"){
          PFS_pvalue = NA_real_
        } else {
          PFS_pvalue = try(coef(summary(PFS_fit))["value","Pr(>|z|)"], silent = TRUE) %>% as.numeric()
          PFS_beta = try(coef(summary(PFS_fit))["value","coef"], silent = TRUE) %>% as.numeric()
        }
        
        ## put all results together
        data.frame("biomarker"=biomarker, 
                   "N" = n,
                   "missPct" = misspct,
                   "mean_value"=mean_value, 
                   "sd_value"=sd_value,
                   "skewness_value"=skewness_value,
                   "median_value"=median_value,
                   "PFS_pvalue"=PFS_pvalue,
                   "PFS_beta"=PFS_beta)
      }
    }
    
    ## split the data by biomarker
    regression_res_list = tempdata_sub %>% 
      purrr::map(~ regression(.)) 
    
    ## combine all biomarkers' results
    output_list[[trt+1]] = do.call(rbind.data.frame, regression_res_list) %>% mutate(Treatment = trt)
  }
  ## combine results from two treatment arms
  output_df = do.call(rbind.data.frame, output_list)
  return (output_df)
}

####################################################

dset_long = dset %>% 
  gather(key = "biomarker", value = "value", -ID, -PFS_P, -PFS_event, -AGE, -SEXN, -TRT01N) 

## call the function and get the p-value of main dataset 
univariate_main_bytrt = uni_analysis_byarm(dset_long)
save(univariate_main_bytrt, file="output\\singleP\\univariate_main_bytrt.RDa")


############## PFs rna Association ########################

rna <- rna_v3

## transpose the rna dataset
rna_t = rna %>% 
  select(-HUGO) %>% 
  t() %>% 
  as.data.frame 
names(rna_t) <- rna$HUGO
rna_t$ID = rownames(rna_t)

## combine the biomarker data with the clinical outcome
rna_wide = rna_t %>%
  left_join(x = dset %>% select(ID, PFS_P, PFS_event, TRT01N, AGE, SEXN), y = ., by = "ID") 

rna_long = rna_wide %>%
  gather(key = "biomarker", value = "value", -ID, -PFS_P, -PFS_event, -TRT01N, -AGE, -SEXN)

## call the function and get the p-value of rna genes 
univariate_rna_bytrt = uni_analysis_byarm(rna_long)
save(univariate_rna_bytrt, file="output\\singleP\\univariate_rna_bytrt.RDa")
rm(rna_long)
############## PFS DNA Association ########################

mut <- mut_v3

## transpose the mut dataset
mut_t = mut %>% 
  select(-HUGO) %>% 
  t() %>% 
  as.data.frame 
names(mut_t) <- mut$HUGO
mut_t$ID = rownames(mut_t)

## combine the biomarker data with the clinical outcome
mut_wide = mut_t %>%
  left_join(x = dset %>% select(ID, PFS_P, PFS_event, TRT01N, AGE, SEXN), y = ., by = "ID") 

mut_long = mut_wide %>%
  gather(key = "biomarker", value = "value", -ID, -PFS_P, -PFS_event, -TRT01N, -AGE, -SEXN)

## call the function and get the p-value of mutation genes 
univariate_mut_bytrt = uni_analysis_byarm(mut_long)
save(univariate_mut_bytrt, file="output\\singleP\\univariate_mut_bytrt.RDa")


