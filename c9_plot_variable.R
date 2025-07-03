
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

load("output\\anaDset\\dset_std4ana.RDa")
dset_all <- dset_std4ana

load("output\\pred_result\\final_subgroup_pred_final_RMST_v10.RDa")
load("output\\pred_result\\arm0_train_pred_RMST.RDa")
load("output\\pred_result\\arm1_train_pred_RMST.RDa")

##################

## look at the arm1 vs. arm0 in each subgroup
mydset_HR_4plot <- mydset_HR %>%
  mutate(Subgroup = ifelse(Subgrps == 3, "Subgroup 1",
                           ifelse(Subgrps == 4, "Subgroup 2",
                                  ifelse(Subgrps == 6, "Subgroup 3", "Subgroup 4"))),
         Treatment = ifelse(TRT01N==0, "Sunitinib", "Avelumab+axitinib"))
  

p1 <- survfit(Surv(PFS_P, PFS_event) ~ Treatment, data=mydset_HR_4plot)
cols <-  c("#E41A1C", "#3060C0") # dark red vs. dark blue
names(cols) <- c("Sunitinib", "Avelumab+axitinib")
myp <- ggsurvplot_facet(p1, mydset_HR_4plot, facet.by = "Subgroup", 
                        conf.int = TRUE, size=1,
                 palette = cols, short.panel.labs=T, 
                 panel.labs.font=aes(size=12, face="bold"),
                 ylab="PFS Probability", xlab="Time (month)")
windows()
ggpar(myp, 
      font.main = c(12, "bold"),
      font.x = c(12, "bold"),
      font.y = c(12, "bold"),
      font.caption = c(12, "bold"), 
      font.legend = c(12, "bold"), 
      font.tickslab = c(12, "bold"))

#####################

## if assign as two group, good vs. bad
mydset_HR <- mydset_HR %>%
  mutate(two_grp = ifelse(Subgrps %in% c(6,7), "Good_group", "bad_grup"))

## calculate the HR of arm 1 vs. arm0 in each group
hr_results <- NULL
for (mygrp in unique(mydset_HR$two_grp)) {
  print(mygrp)
  tmp <- mydset_HR %>% filter(two_grp==mygrp) %>% 
    coxph(Surv(PFS_P, PFS_event) ~ TRT01N, data = .) %>% summary
  hr_results <- rbind(hr_results, c(mygrp, tmp$conf.int))
}
hr_results <- as.data.frame(hr_results)
colnames(hr_results) <- c("Group","HR","neg_HR","LB","UB")
hr_results %>% arrange(Group)

#####################

## for the important biomarker, split them to high vs. low two groups
mydset_HR <- mydset_HR %>%
  mutate(Treatment = ifelse(TRT01N==0, "Sunitinib", "Avelumab+axitinib")) %>%
  mutate(median_SYDE2 = median(SYDE2),
         median_SCIN = median(SCIN),
         median_B9991003_r5_Cell_cell_signaling = median(B9991003_rc5_Cell_cell_signaling),
         SYDE2_C = ifelse(SYDE2 >= median_SYDE2, "SYDE2 >= median", "SYDE2 < median"),
         SCIN_C = ifelse(SCIN >= median_SCIN, "SCIN >= median", "SCIN < median"),
         B9991003_rc5_Cell_cell_signaling_C = ifelse(B9991003_rc5_Cell_cell_signaling >= median_B9991003_r5_Cell_cell_signaling, "rc5_Cell_cell_signaling >= median", "rc5_Cell_cell_signaling < median"))


combo <- mydset_HR[mydset_HR$Treatment=="Avelumab+axitinib",]
chemo <- mydset_HR[mydset_HR$Treatment=="Sunitinib",]

# Visualize
#++++++++++++++++++++++++++++++++++++
# "#E41A1C" darkred,  "#F880A2" light red, 
# "#183060" dark blue, "#6798D7" light blue

############################################################
######## SYDE2 #######################################
############################################################

BMx_high = mydset_HR[mydset_HR$SYDE2_C=="SYDE2 >= median",]
BMx_low = mydset_HR[mydset_HR$SYDE2_C=="SYDE2 < median",]


### combo

cols <-  c("#6798D7", "#183060")
p1 <- survfit(Surv(PFS_P, PFS_event) ~ SYDE2_C, data=combo)
SP1 <- ggsurvplot(p1, data=combo, size=2,
                  palette = cols, pval=TRUE, pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("SYDE2 < median", "SYDE2 >= median"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="Avelumab+axitinib",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) 
SP11 <- SP1$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )

cols <-  c("#F880A2", "#E41A1C") 
p2 <- survfit(Surv(PFS_P, PFS_event) ~ SYDE2_C, data=chemo)
SP2 <- ggsurvplot(p2, data=chemo, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("SYDE2 < median", "SYDE2 >= median"), 
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="Sunitinib",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) 
SP22 <- SP2$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


cols <-  c("#183060", "#E41A1C") 
p3 <- survfit(Surv(PFS_P, PFS_event) ~ Treatment, data=BMx_high)
SP3 <- ggsurvplot(p3, data=BMx_high, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("Avelumab+axitinib", "Sunitinib"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(size=16), tables.y.text=F, 
                  title="SYDE2 >= median",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) 
SP33 <- SP3$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


cols <-  c("#6798D7", "#F880A2") 
p4 <- survfit(Surv(PFS_P, PFS_event) ~ Treatment, data=BMx_low)
SP4 <- ggsurvplot(p4, data=BMx_low, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("Avelumab+axitinib", "Sunitinib"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="SYDE2 < median",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold"))
SP44 <- SP4$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


p1 <- ggplot() + theme_void()
windows()
grid.arrange(SP44, SP33, p1, SP11, SP22, SP4$table, SP3$table, p1, SP1$table, SP2$table, 
             nrow = 2, ncol=5,
             widths = c(1, 1, 0.4, 1, 1), heights=c(5,1))


############################################################
######## SCIN #######################################
############################################################

BMx_high = mydset_HR[mydset_HR$SCIN_C=="SCIN >= median",]
BMx_low = mydset_HR[mydset_HR$SCIN_C=="SCIN < median",]


### combo

cols <-  c("#6798D7", "#183060") 
p1 <- survfit(Surv(PFS_P, PFS_event) ~ SCIN_C, data=combo)
SP1 <- ggsurvplot(p1, data=combo, size=2,
                  palette = cols, pval=TRUE, pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("SCIN < median", "SCIN >= median"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="Avelumab+axitinib",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) 
SP11 <- SP1$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )

cols <-  c("#F880A2", "#E41A1C") 
p2 <- survfit(Surv(PFS_P, PFS_event) ~ SCIN_C, data=chemo)
SP2 <- ggsurvplot(p2, data=chemo, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs =c("SCIN < median", "SCIN >= median"), 
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="Sunitinib",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) 
SP22 <- SP2$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


cols <-  c("#183060", "#E41A1C") 
p3 <- survfit(Surv(PFS_P, PFS_event) ~ Treatment, data=BMx_high)
SP3 <- ggsurvplot(p3, data=BMx_high, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("Avelumab+axitinib", "Sunitinib"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(size=16), tables.y.text=F, 
                  title="SCIN >= median",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) 
SP33 <- SP3$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


cols <-  c("#6798D7", "#F880A2")
p4 <- survfit(Surv(PFS_P, PFS_event) ~ Treatment, data=BMx_low)
SP4 <- ggsurvplot(p4, data=BMx_low, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("Avelumab+axitinib", "Sunitinib"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="SCIN < median",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold"))
SP44 <- SP4$plot +
  theme(
    plot.title = element_text(size = 22, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


p1 <- ggplot() + theme_void()
windows()
grid.arrange(SP44, SP33, p1, SP11, SP22, SP4$table, SP3$table, p1, SP1$table, SP2$table, 
             nrow = 2, ncol=5,
             widths = c(1, 1, 0.4, 1, 1), heights=c(5,1))



############################################################
######## B9991003_rc5_Cell_cell_signaling_C #######################################
############################################################

BMx_high = mydset_HR[mydset_HR$B9991003_rc5_Cell_cell_signaling_C=="rc5_Cell_cell_signaling >= median",]
BMx_low = mydset_HR[mydset_HR$B9991003_rc5_Cell_cell_signaling_C=="rc5_Cell_cell_signaling < median",]


### combo

cols <-  c("#6798D7", "#183060") 
p1 <- survfit(Surv(PFS_P, PFS_event) ~ B9991003_rc5_Cell_cell_signaling_C, data=combo)
SP1 <- ggsurvplot(p1, data=combo, size=2,
                  palette = cols, pval=TRUE, pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("Cell_cell_signaling < median", "Cell_cell_signaling >= median"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="Avelumab+axitinib",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) +
  guides(color = guide_legend(nrow = 2))
SP11 <- SP1$plot +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )

cols <-  c("#F880A2", "#E41A1C")
p2 <- survfit(Surv(PFS_P, PFS_event) ~ B9991003_rc5_Cell_cell_signaling_C, data=chemo)
SP2 <- ggsurvplot(p2, data=chemo, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs =  c("Cell_cell_signaling < median", "Cell_cell_signaling >= median"), 
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="Sunitinib",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) +
  guides(color = guide_legend(nrow = 2))
SP22 <- SP2$plot +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


cols <-  c("#183060", "#E41A1C") 
p3 <- survfit(Surv(PFS_P, PFS_event) ~ Treatment, data=BMx_high)
SP3 <- ggsurvplot(p3, data=BMx_high, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("Avelumab+axitinib", "Sunitinib"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(size=16), tables.y.text=F, 
                  title="B9991003_rc5_Cell_cell_signaling >= median",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold")) 
SP33 <- SP3$plot +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


cols <-  c("#6798D7", "#F880A2")
p4 <- survfit(Surv(PFS_P, PFS_event) ~ Treatment, data=BMx_low)
SP4 <- ggsurvplot(p4, data=BMx_low, size=2,
                  palette = cols, pval=TRUE,pval.method=TRUE,
                  pval.size = 6, pval.method.size = 6,
                  surv.median.line = "hv",
                  xlab="Time (month)", ylab="PFS Probability",
                  legend.labs = c("Avelumab+axitinib", "Sunitinib"),
                  risk.table = TRUE, tables.height = 0.15, risk.table.fontsize = 17 / .pt,
                  tables.theme = theme_cleantable(), tables.y.text=F, 
                  title="B9991003_rc5_Cell_cell_signaling < median",
                  gtheme = theme_bw(),
                  font.x = c(16, "bold"),
                  font.y = c(16, "bold"),
                  font.tickslab = c(12, "bold"))
SP44 <- SP4$plot +
  theme(
    plot.title = element_text(size = 16, face = "bold"),
    legend.title = element_blank(),
    legend.text = element_text(size = 14, face = "bold")
  )


p1 <- ggplot() + theme_void()
windows()
grid.arrange(SP44, SP33, p1, SP11, SP22, SP4$table, SP3$table, p1, SP1$table, SP2$table, 
             nrow = 2, ncol=5,
             widths = c(1, 1, 0.4, 1, 1), heights=c(5,1))
