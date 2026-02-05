library(tidyverse)
library(dplyr)
library(glmmTMB)

data<-read.csv("~/Dropbox/ellie_projects/data/merged.working.csv")
str(data)
head(data)
# Experiment with “group by” in dplyr to do some calculation.
# Use mutate and summarise on your data. How are these different?
#   Commit and push to GitHub. Be sure to update your README!

data_2<- data%>%
  group_by(species)%>%
  count()

data<- data %>%
  mutate(avglogVPD = log10(avgVPD + 1)) #creates new columns

vpddat<-data %>%
  group_by(species)%>%
summarise(mean_vpd_spp=mean(avglogVPD,na.rm=TRUE)) 

#some spp have no VPD data

data<-data%>%
  filter(species!="COTO")%>%
  filter(species!="MYAU")%>%
  filter(species!="SUBSTRATE")%>%
  filter(species!="LANO")%>%
  filter(species!="MYGR")

data<-data%>%
  select(site,date,state,source,species,temp,season,lgdL,winter_year,avgVPD,avgTEMP,avgRH,avglogVPD)
sum(is.na(data$avglogVPD))

epfu_dat<-data %>%
  filter(species=="EPFU")

pesu_dat<-data %>%
  filter(species=="PESU")

mylu_dat<-data %>%
  filter(species=="MYLU")

mod1<-glmmTMB(avglogVPD~species + (1|site), data = data,family = gaussian(link = "identity"))
summary(mod1)

#try out pivot wider/longer

#try making a dataset where each row is a species, and each column is 

data_wide<-data%>%
  pivot_wider(names_from = state, 
              values_from = species) 

