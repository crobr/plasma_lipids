
# Clean global environment
rm(list = ls())

#Libraries
library(survival)
library(survminer)
library(M3C)

#work computer
setwd("users/your_file_location")


# df <- read.csv("Clinical_met_noMDS_Final.csv")
df<-read.csv("Clinical_lip_noMDS_Final.csv")

sfit<-survfit(Surv(time_OS, status)~1,data=df)

ggsurvplot(sfit)


##filter data based on variability 
lip_log2 <-t(subset (df, select=-c(1:2)))
filtered<-featurefilter(lip_log2, percentile=5, method = "var",topN=105)
idx<-data.frame(filtered$filtered_data)
rows<- rownames(idx)
top100<-df[rows]
