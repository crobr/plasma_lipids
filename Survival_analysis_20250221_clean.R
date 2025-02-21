# Clean global environment
rm(list = ls())

#load libraries
library(survival)
library(survminer)
library(broom)
library(finalfit)
library(tidyverse)
library(tidytidbits)
library(survivalAnalysis)
library(ranger)
library(ggplot2)
library(dplyr)
library(ggfortify)
library(dcurves)
library(rms)
library(gtsummary)
library(formattable)

#Set working environment
working<-('users/your_file_location')
setwd(working)


##load original data
d<- read.csv("Clinical_MVA_noMDS_Final_20250210.csv") 

##KM analysis 
train<-subset(d, set %in% c("train"))
test<-subset(d, set %in% c("test"))
cut_reg_o<-surv_cutpoint(
  test,
  time = "time_OS",
  event = "status",
  progressbar = TRUE,
  variables=c("pred_CR", "pred_NR"))

cat_reg_o<- surv_categorize(cut_reg_o, variables= NULL)

fit_org<- survfit(Surv(time_OS, status) ~pred_CR+pred_NR, data=cat_reg_o)
ggsurvplot(fit_org,
           data=cat_reg_o, 
           pval=TRUE,
           # legend.title ="PC 19:2_22:6",
           legend.labs = c("Predicted low risk", "Predicted high risk" ),
           xlab = "Days",
           risk.table=TRUE, 
           conf.int=FALSE,
           tables.height = 0.2,
           font.y = c(14, face = "bold"), font.tickslab = c(14), 
           font.legend = c(14)) 


##Multivariable analysis 
cox1<- (coxph(Surv(time_OS, status) ~pred_NR+ELN_risk+AGE, data=test))
cs1<-concordance(cox1, timewt="n/G2")  # Uno's weighting
summary(cox1)

cox2<- (coxph(Surv(time_OS, status) ~ELN_risk+AGE, data=test))
cs2<-concordance(cox2, timewt="n/G2")  # Uno's weighting
summary(cox2)

cox3<- (coxph(Surv(time_OS, status) ~pred_NR+ELN_risk, data=test))
cs3<-concordance(cox3, timewt="n/G2")  # Uno's weighting
summary(cox3)

full_tbl<-tbl_regression(cox1,exponentiate = TRUE)

no_CR_tbl<-tbl_regression(cox2,exponentiate = TRUE) 

full_tbl$table_body$label<-c("Predicted NR", "ELN risk", "Age")
no_CR_tbl$table_body$label<-c("ELN risk", "Age")

CR_ELN_tbl<-tbl_regression(cox3,exponentiate = TRUE) 
CR_ELN_tbl$table_body$label<-c("Predicted NR", "ELN risk")

tbl_merge(
  list(full_tbl,no_CR_tbl),
  tab_spanner = c("**Clinical features with predicted NR,\n C-score: 0.789**", "**Clinical features only,\n  C-score: 0.711**")
)

####VALIDATION 
val<- read.csv("20250210_for_survival_Val.csv") #unprocessed


cut_reg_v<-surv_cutpoint(
  val,
  time = "time_OS",
  event = "status",
  progressbar = TRUE,
  variables=c("pred_CR","pred_NR"))
cat_reg_v<- surv_categorize(cut_reg_v, variables= NULL)

##KM analysis

fit_val<- survfit(Surv(time_OS, status) ~pred_CR+pred_NR, data=cat_reg_v)

ggsurvplot(fit_val,
           data=cat_reg_v, 
           pval=TRUE,
           legend.labs = c("Predicted low risk", "Predicted high risk" ),
           xlab = "Days",
           risk.table=TRUE, 
           conf.int=FALSE,
           tables.height = 0.2,
           font.y = c(14, face = "bold"), font.tickslab = c(14), 
           font.legend = c(14)) 

##Single lipid analysis 

list<- c("SM.d44.1.", "MGDG.20.4_22.5.",  "PC.17.0_18.1.", "SM.d16.2_27.0.",
         "AcCa.28.6.","AcHexZyE.28.6.", "SM.d35.2.", "PC.O.36.2.",
         "MGDG.18.2_22.6.", "PC.17.0_18.2.","PC.18.0_18.1.","DG.O.13.0_19.1.",
         "LPG.18.1.","ChE.18.1.",  "MePC.40.4.",  "DG.24.6.",
         "CerPE.d18.1_16.0.",  "SM.d43.3.","TG.10.0_18.2_18.2.")
##original 
data_org<- read.csv("Clinical_lip_noMDS_Final_Survival.csv") #unprocessed
cut_reg_lip<-surv_cutpoint(
  data_org,
  time = "time_OS",
  event = "status",
  progressbar = TRUE,
  variables=list)
cat_reg_lip<- surv_categorize(cut_reg_lip, variables= NULL)
cat_reg_lip<-data.frame(cat_reg_lip)

fit<- survfit(Surv(time_OS, status)~SM.d44.1., data=cat_reg_lip)
fit1<- survfit(Surv(time_OS, status)~MGDG.20.4_22.5., data=cat_reg_lip)
fit2<- survfit(Surv(time_OS, status)~PC.17.0_18.1., data=cat_reg_lip)
fit3<- survfit(Surv(time_OS, status)~SM.d16.2_27.0., data=cat_reg_lip)
fit4<- survfit(Surv(time_OS, status)~AcCa.28.6., data=cat_reg_lip) 
fit5<- survfit(Surv(time_OS, status)~DG.O.13.0_19.1., data=cat_reg_lip)
fit6<- survfit(Surv(time_OS, status)~AcHexZyE.28.6., data=cat_reg_lip)

ggsurvplot(fit6,
           data=cat_reg_lip, 
           pval=TRUE,
           legend.labs = c("AcHexZyE.28.6. high", "AcHexZyE.28.6. low" ),
           xlab = "Days",
           risk.table=TRUE, 
           conf.int=FALSE,
           tables.height = 0.2,
           font.y = c(14, face = "bold"), font.tickslab = c(14), 
           font.legend = c(14),) 
#Validation
#work computer
setwd('/Users/crobrien/OneDrive - UHN/Jones Lab/3. Lipids-Relapse AML/2. Plasma Sample Lipidomics Data/New folder/Final Data/Colarado validation')

data_val<-read.csv("Val_Col_lip_V2_adjusted.csv")
cut_reg_lip2<-surv_cutpoint(
  data_val,
  time = "time_OS",
  event = "status",
  progressbar = TRUE,
  variables=list)
cat_reg_lip2<- surv_categorize(cut_reg_lip2, variables= NULL)
cat_reg_lip2<-data.frame(cat_reg_lip2)

fit<- survfit(Surv(time_OS, status)~SM.d44.1., data=cat_reg_lip2)
fit1<- survfit(Surv(time_OS, status)~MGDG.20.4_22.5., data=cat_reg_lip2)
fit2<- survfit(Surv(time_OS, status)~PC.17.0_18.1., data=cat_reg_lip2)
fit3<- survfit(Surv(time_OS, status)~SM.d16.2_27.0., data=cat_reg_lip2)
fit4<- survfit(Surv(time_OS, status)~AcCa.28.6., data=cat_reg_lip2) 
fit5<- survfit(Surv(time_OS, status)~DG.O.13.0_19.1., data=cat_reg_lip2)
fit6<- survfit(Surv(time_OS, status)~AcHexZyE.28.6., data=cat_reg_lip2)

ggsurvplot(fit6,
           data=cat_reg_lip2, 
           pval=TRUE,
           legend.labs = c("AcHexZyE.28.6. high", "AcHexZyE.28.6. low" ),
           xlab = "Days",
           risk.table=TRUE, 
           conf.int=FALSE,
           tables.height = 0.2,
           font.y = c(14, face = "bold"), font.tickslab = c(14), 
           font.legend = c(14),) 

