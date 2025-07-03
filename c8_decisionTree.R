
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
library(partykit)
library(glmnet)
library(randomForestSRC)

setwd("J:\\projects\\Research\\JAVELIN")

## load the data with 8331 biomarker and 683 (726) subjects
load("output\\anaDset\\dset_std4ana.RDa")
dset_all <- dset_std4ana

## load the data from previous code, including 
load("output\\pred_result\\final_subgroup_pred_RMST.RDa")

## we only need the selected 97 biomarker and the outcome (diff_1_0)
dset_4_tree <- rna_4_PRISM %>%
  select(-ID, -PFS_P, -PFS_event, -TRT01N, -mu_0, -mu_1, -pi_0, -pi_1)

##### Fit a decision tree #######

library(rpart)
fit1 <- rpart(diff_1_0 ~ ., data = dset_4_tree,
            minbucket = floor(dim(dset_4_tree)[1]*0.14), maxdepth = 4)
fit1

## plot your tree
library(rpart.plot)
windows()
rpart.plot(fit1)

#### 

mydset_HR <- rna_4_PRISM
mydset_HR$Subgrps <-  fit1$where
table(mydset_HR$Subgrps) 

## summarize the mean of delta RMST of each subgroup
mydset_HR %>%
  group_by(Subgrps) %>%
  summarise(mean_diff=mean(diff_1_0)) %>%
  select(Subgrps, mean_diff)

## save the data with subgroup variable
save(mydset_HR, file="output\\pred_result\\final_subgroup_pred_final_RMST_v10.RDa")
