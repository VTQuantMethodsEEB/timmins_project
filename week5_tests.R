rm(list=ls())

library(dplyr)

data1<-read.csv("mylu.working.extrapolated")  
#I wanted to do this for comparing species, but the data file was too big to upload!!!

unique(data$state)
#this doesn't run

##Code corresponding to LE5 - Statistical Tests

##permutation tests####
  
#permutation test hypothesis:

# mean log VPD is different between NY MYLU bats and WI MYLU bats

library(ggplot2)

theme_set(theme_bw() +
            theme(
              plot.title = element_text(size = 30, hjust = 0.5),
              panel.grid = element_blank(),
              axis.title = element_text(size = 30),
              axis.text = element_text(size = 15),
              axis.line = element_line(),
              legend.text = element_text(size = 12),
              strip.text = element_text(size = 20, color = "black", family = "Arial"),
              strip.background = element_blank(),
            ))



##how to write your own permutation test##
set.seed(101)
#set.seed will set R's random number generator to start at the same place
#this ensures that when you, and I, and anyone else, does the test, we will all get the same results



#you always need to do something like this when you run a "for" loop
#you could also write res <- numeric(1000), which would give you a list of 1000 0's
#the important thing to have a vector already named "res"

data<-data1%>%
  filter(state %in% c("WI", "NY"))%>%
  filter(!is.na(avglogVPD))%>%
  select(state, avglogVPD)



#use jar example to illustrate sampling procedure
res <- NA ## set aside space for results
obs = mean(data$avglogVPD[data$state=="WI"]) - mean(data$avglogVPD[data$state=="NY"])

for (i in 1:1000) {
  
  VPDboot <- sample(data$avglogVPD)
  
  WIboot <- VPDboot[1:length(data$state[data$state=="WI"])] #this says assign the first six colonies to forest
  NYboot <- VPDboot[(length(data$state[data$state=="WI"])+1):length(data$state)] #this says assign the rest of the observations to field

  res[i] <- mean(WIboot) - mean(NYboot) #KL: this needed to be reversed
  
}
#what is our observed mean difference?
#obs <- mean(data$avglogVPD[data$state=="NY"]) -
#  mean(data$avglogVPD[data$state=="WI"])
obs #observed residuals _ KL this is not a residual

range(res)

hist(res,col="gray",las=1,main="",xlim = range(c(res, obs)))
abline(v=obs,col="red")

##so how do we get our p-value?
res[res>=obs]
length(res[res>=obs])
1/1000
#small p-value! 
mean(res>=obs)

#maybe need to do two-tailed test?

#using mean(permutations>=obs)) is a trick to calculate the proportion:
#the logical statement returns a logical (FALSE/TRUE) vector, which then gets converted to a 0/1 vector when you ask R to take the mean,
#so this is equivalent to counting the number of true values and dividing by the length


#what if we want a two-tailed?


#we could count the area in both tails
mean(abs(res)>=abs(obs))

#The difference in avgVPD between NY and WI is statistically significant (p ≈ 0.004)
#according to the two-tailed test

# classic test


#####Shapiro-Wilk Test#######
#Is our avglogVPD data normally distributed?##
#The null hypothesis is that the data are normally distributed
#P<0.05 indicates NOT normal

swt<-shapiro.test(data$avglogVPD)
swt

#this data is not normal!! :(

