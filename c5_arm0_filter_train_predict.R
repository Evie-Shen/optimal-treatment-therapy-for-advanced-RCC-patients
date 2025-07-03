
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
library(randomForestSRC)
library(ranger)
library(glmnet)

setwd("J:\\projects\\Research\\JAVELIN")

load("output\\anaDset\\dset_std4ana.RDa")
dset_all <- dset_std4ana

### ARM 0 filter, traninng, prediction

## only pick up significant biomarkers from arm 0 analysis
genelist <- pval_sig05 %>% 
  filter(PFS_pvalue_0 < 0.05) %>% 
  select(biomarker_2) %>% 
  filter(biomarker_2!="CD8_INVASIVE_MARGIN_SURFACE_AREA")    ## 6670 genes
dset_sub <- dset_all %>% 
  select(ID, PFS_P, PFS_event, TRT01N, AGE, SEXN, genelist$biomarker_2) 

## subset data with only including patients from arm 0
dset_sub_arm <- dset_sub %>% filter(TRT01N==0) 

######################

x_train <- dset_sub_arm %>% select(-ID, -TRT01N,-PFS_P, -PFS_event)
x_train = as.matrix(x_train)
y_train <- dset_sub_arm %>% select(PFS_P, PFS_event) %>% rename("time"="PFS_P", "status"="PFS_event") %>% as.matrix
## set seeds so that the analysis results is repeatable
set.seed(42)

#### glmnet cv to get the best alpha and lambda who give the minimum cv error

alpha_values <- seq(0.1, 1, by = 0.1) 
cv_results <- list()

for (alpha in alpha_values) {
  print(alpha)
  cv_fit <- cv.glmnet(x_train, y_train, family = "cox", nfolds = 10, alpha = alpha)
  cv_results[[paste0("alpha_", alpha)]] <- cv_fit
}

optimal_lambdas <- sapply(cv_results, function(cv_fit) cv_fit$lambda.min)

## find the alpha which gave the minimum cv error
best_alpha_index <- which.min(sapply(cv_results, function(cv_fit) min(cv_fit$cvm)))
best_alpha <- alpha_values[best_alpha_index]
## find the lambda from the best alpha which gave the minimum cross validation error
best_lambda <- optimal_lambdas[best_alpha_index]
best_alpha
best_lambda

#### glmnet training

# best_alpha = 0.4
# best_lambda = 0.284805

## Fit the final model with the best alpha and its corresponding lambda
final_model <- glmnet(x_train, y_train, family = "cox", alpha = best_alpha, lambda = best_lambda) 
coefs <- coef(final_model, s = "best_lambda")
## select the biomarkers with non-zeor coefficients
glmnet_sel <- as.data.frame(cbind(rownames(coefs)[which(coefs != 0)],  coefs[which(coefs != 0)]))
names(glmnet_sel) <- c("biomarker_2","importance")
glmnet_sel %>% arrange(desc(importance))

## merge the selected important BMx back to the p-value table
var_glmnet <- glmnet_sel %>% 
  inner_join(pval_sig05, by="biomarker_2") %>% 
  arrange(desc(importance)) %>%
  mutate(index = row_number()) %>%
  relocate(index,importance, `PFS_pvalue_0`)
## var_glmnet may miss the age. if glmnet_sel has age, add the age variable
var_glmnet <- rbind.fill(var_glmnet, glmnet_sel %>% filter(biomarker_2=="AGE")) 


#### SRC training

## put the selected BMx as X and PFS as Y into a new dataset
rna_4_SRC <- dset_sub_arm %>% 
  select(PFS_P, PFS_event, var_glmnet$biomarker_2)

## using random Forest SRC to train a prediction model
o.grow <- rfsrc(Surv(PFS_P, PFS_event) ~ ., rna_4_SRC) 
print(o.grow)

##### SRC prediction

### using the predictive model from random Forest to do the prediction on all patients
### each patient will have a predicted survival curve

rna_4_pred <- dset_all %>% select(ID, var_glmnet$biomarker_2)
mu_pred <- predict.rfsrc(o.grow, rna_4_pred)
pred_matrix <- data.frame(mu_pred$time.interest, t(mu_pred$survival))
names(pred_matrix) <- c("rectime", rna_4_pred$ID)

## integral the area under the predicted curve to generate a mean survival time up to 15 months
RMST <- apply(pred_matrix[,-1], 2, 
              function(x){int_spline_0 <- splinefun(pred_matrix$rectime, x, method = "natural")
                return(integrate(f = int_spline_0, lower = 0, upper = 15 * 1)$value)
              })

dset_all$mu_0 = RMST
dset_all_wmu_trt0 <- dset_all
var_glmnet_trt0 <- var_glmnet
pred_matrix_trt0 <- pred_matrix %>% mutate(ARM = "ARM 0")

save(dset_all_wmu_trt0, var_glmnet_trt0, pred_matrix_trt0, file="output\\pred_result\\arm0_train_pred_RMST.RDa")

