rm(list = ls())

library(performance)
library(dplyr)
mylu.working.extrapolated<-read.csv("mylu.working.extrapolated")

mylu.working.extrapolated$phase.original = as.factor(mylu.working.extrapolated$phase.original)
mylu.working.extrapolated$season = as.factor(mylu.working.extrapolated$season)

######################
####### Part 1 #######
######################

### Hypothesis ###

#Vapor pressure deficit in early hibernation will significantly impact fungal loads in late hibernation

### examining VPD distribution ###

hist(mylu.working.extrapolated$avg_early_VPD)
hist(mylu.working.extrapolated$avgVPD)
hist(mylu.working.extrapolated$avglogVPD)
#none of this looks nice

### Model and interpretation ###

mod1<-lm(lgdL~avg_early_VPD, data = subset(mylu.working.extrapolated, season=="hiber_late"))
summary(mod1)

#Vapor pressure deficit from early hibernation has a statistically clear 
# negative impact on fungal loads in late hibernation. In other words, drier 
# conditions in early hibernation decrease fungal loads in late hibernation. 
# Note: I have not accounted for mortality here

### Model Diagnostics - doesn't look good ###

hist(resid(mod1)) #looks very left skewed. Remember this is the difference between
#your predicted and observed data, so you want this to be very normal looking

shapiro.test(resid(mod1))
#p val <0.5, this is not normal.

check_model(mod1) 

#POSTERIOR PREDICTIVE CHECK
# The observed data has a weird blip at lower fungal loads. This is probably related to 
# the left skew in residual histogram

#HOMOGENEITY OF VARIANCE
# does not quite seem horizontal

#NORMALITY OF RESIDUALS 
# Gets weird at the two edges

#Linearity
# this is not flat and debatable that its horizontal

mod1<-lm(lgdL~avg_early_VPD, data = subset(mylu.working.extrapolated, season=="hiber_late"))
plot(mod1)

# according to QQ plot, it looks very left skewed. 

#transformation:
# reflect then base 10 

### Plotting ###

library(ggplot2)
r=ggplot(data=subset(mylu.working.extrapolated, season=="hiber_late"), aes(x=avg_early_logVPD, y=lgdL))+ 
  geom_point()+
  stat_smooth(method = "lm")+
  theme_bw() + 
  theme(axis.title=element_text(size=20),axis.text=element_text(size=10),panel.grid = element_blank(), axis.line=element_line(),legend.position="top",legend.title=element_blank())
print(r)

#### Try with other VPD models ####

mod2<-lm(lgdL~avgVPD, data = subset(mylu.working.extrapolated, season=="hiber_earl"))
summary(mod2)

#Vapor pressure deficit has a statistically clear 
# positive impact on fungal loads in early hibernation. In other words, drier 
# conditions in early hibernation increase fungal loads in early hibernation. 

hist(resid(mod2)) 
shapiro.test(resid(mod2)) #not normal
plot(mod2)
#according to QQ this is light tailed
#Scale-location: residuals are high at low fitted values
# residuals vs leverage showing there might be a really weird outlier going on

#### what does data reflection look like ####

c<-max(mylu.working.extrapolated$avgVPD, na.rm = TRUE) + 0.01

mylu.working.extrapolated<-mylu.working.extrapolated%>%
  mutate(r.avgVPD=c-avgVPD)

mylu.working.extrapolated<-mylu.working.extrapolated%>%
  mutate(log.r.avgVPD=log10(r.avgVPD))

hist(mylu.working.extrapolated$log.r.avgVPD) #still looks ugly :(

######################
####### Part 2 #######
######################

