
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
library(lmtest)
library(ranger)
library(glmnet)
library(randomForestSRC)

setwd("J:\\projects\\Research\\JAVELIN")

load("output\\anaDset\\dset_std4ana.RDa")
dset_all <- dset_std4ana

### load the prediction results for drug 0 
### load the prediction results for drug 1 
load("output\\pred_result\\arm0_train_pred_RMST.RDa")
load("output\\pred_result\\arm1_train_pred_RMST.RDa")

############

## pick up biomarkers with significant p-values in arm*biomarker interaction  effect
int_genelist <- pval_sig05 %>% 
  filter(PFS_int_Pval < 0.05) %>% 
  select(biomarker_2) %>% 
  filter(biomarker_2!="CD8_INVASIVE_MARGIN_SURFACE_AREA")    ## 3559 genes
## combine biomarkers with significant interaction p-values and 
## selected biomarkers from arm0 (36 biomarker) and arm1 (22 biomarkers)
## take unique set since some of them are overlap
genelist <- unique(c(int_genelist$biomarker_2, var_glmnet_trt0$biomarker_2, var_glmnet_trt1$biomarker_2)) #3578 genes

## put the selected biomarker list as well as outcome (y) variables together into one table
dset_sub <- dset_all_wmu_trt0 %>% 
  inner_join(dset_all_wmu_trt1 %>% select(ID, mu_1), by="ID") %>%
  select(ID, PFS_P, PFS_event, TRT01N, AGE, SEXN, mu_0, mu_1, genelist) %>%
  mutate(diff_1_0 = mu_1 - mu_0,
         pi_0 = 0.5139092,
         pi_1 = 0.4860908) 
# prop.table(table(dset_sub$TRT01N))

######################

x_train <- dset_sub %>% select(-ID, -PFS_P, -PFS_event, -TRT01N, -mu_0, -mu_1, -pi_0, -pi_1, -diff_1_0)
x_train = as.matrix(x_train)
y_train <- dset_sub %>% select(diff_1_0) %>% as.matrix
## set seeds so that the analysis results is repeatable
set.seed(42)

#### glmnet cross validation 

alpha_values <- seq(0.1, 1, by = 0.1) 
cv_results <- list()

for (alpha in alpha_values) {
  print(alpha)
  cv_fit <- cv.glmnet(x_train, y_train, family = "gaussian", nfolds = 10, alpha = alpha)
  cv_results[[paste0("alpha_", alpha)]] <- cv_fit
}

optimal_lambdas <- sapply(cv_results, function(cv_fit) cv_fit$lambda.1se)

## find the alpha which gave the minimum cv error
best_alpha_index <- which.min(sapply(cv_results, function(cv_fit) min(cv_fit$cvm)))
best_alpha <- alpha_values[best_alpha_index]
## find the lambda from the best alpha which gave the minimum cross validation error
best_lambda <- optimal_lambdas[best_alpha_index]
best_alpha
best_lambda

# best_alpha = 0.6
# best_lambda = 0.2523571 

#### glmnet filtering

## Fit the final model with the best alpha and its corresponding lambda
final_model <- glmnet(x_train, y_train, family = "gaussian", alpha = best_alpha, lambda = best_lambda)
coefs <- coef(final_model, s = "best_lambda")
## select the biomarkers with non-zeor coefficients, have biomarker name and their coefficients
glmnet_sel <- as.data.frame(cbind(rownames(coefs)[which(coefs != 0)],  coefs[which(coefs != 0)]))
names(glmnet_sel) <- c("biomarker_2","importance")
glmnet_sel %>% arrange(desc(importance)) # 98 including intercept

## merge the selected important BMx back to the p-value table
var_glmnet <- glmnet_sel %>% 
  inner_join(pval_sig05, by="biomarker_2") %>% 
  arrange(desc(importance)) %>%
  mutate(index = row_number()) %>%
  relocate(biomarker_2, importance, PFS_int_Pval, PFS_pvalue_1, PFS_pvalue_0, PFS_beta_1, PFS_beta_0, Source)
## var_glmnet may miss the age. if glmnet_sel has age, add the age variable
var_glmnet <- rbind.fill(var_glmnet, glmnet_sel %>% filter(biomarker_2=="AGE")) ## 97 variables, without age

## 97 variables selected

##### glmnet filtering

## put the selected 97 variables and outcome variables into one table
## this table will be used to fit a subgroup tree
rna_4_PRISM <- dset_sub %>% 
  select(ID, PFS_P, PFS_event, TRT01N, mu_0, mu_1, pi_0, pi_1, diff_1_0, var_glmnet$biomarker_2)

save(dset_sub, rna_4_PRISM, var_glmnet, file="output\\pred_result\\final_subgroup_pred_RMST.RDa")

