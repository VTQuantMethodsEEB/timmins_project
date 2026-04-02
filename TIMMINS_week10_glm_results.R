rm(list = ls())
library(emmeans)
library(effects)
library(performance)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(glmmTMB)
library(DHARMa)
mylu.working.imputated <- read.csv("mylu.working.imputated")

#HYPOTHESIS

#from Alex's experiment and a priori knowledge, I hypothesize that VPD will decrease fungal loads
#and temperature should increase fungal loads

#I am examining this in the established phase only, because when parsing by
#phase, the null model is best in the invasion phase. I am also parsing to
#late hibernation to match the UV score model, and because it seems to
#be performing better than the early hibernation models

mylu.working.imputated <- mylu.working.imputated %>%
  mutate(adj_gdL = gdL^(1/4))

latehibergdLmod= glm(adj_gdL ~ avglogVPD*avgTEMP,
                         data=subset(mylu.working.imputated, ysw<9 & ysw>-1& 
                                       season=="hiber_late"&gd==1&phase.original=="established"),
                         family = Gamma(link="log"))
summary(latehibergdLmod)
plot(allEffects(latehibergdLmod))

simulationOutput <- simulateResiduals(fittedModel = latehibergdLmod, plot = F)
plot(simulationOutput) #this is still pretty bad

#Vapor pressure deficit significantly decreases fungal loads in late hibernation.
#Temperature significantly increases fungal loads in late hibernation.
#As temperature increases, the effect of vapor pressure deficit on fungal loads
#becomes 


mylu.working.avg <- mylu.working.imputated %>%
  filter(season == "hiber_late"&gd==1&phase.original=="established")%>%
  group_by(site, date,phase.original) %>%   
  summarise(
    avg_avglogVPD = mean(avglogVPD, na.rm = TRUE),
    avg_avgTEMP    = mean(avgTEMP, na.rm = TRUE),
    avg_adj_gdL     = mean(adj_gdL, na.rm = TRUE),
    n_bats        = n(),
    .groups = "drop")

nd4e <- mylu.working.imputated %>%
  filter(ysw >= 0 & ysw < 9, season == "hiber_late"&gd==1&phase.original=="established") %>%
  group_by(phase.original) %>%
  summarise(
    minVPD = min(avglogVPD, na.rm = TRUE),
    maxVPD = max(avglogVPD, na.rm = TRUE),
    minTEMP = min(avgTEMP, na.rm = TRUE),
    maxTEMP= max(avgTEMP, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  do(expand.grid(
    avglogVPD = seq(.$minVPD, .$maxVPD, length.out = 50),
    avgTEMP = seq(.$minTEMP, .$maxTEMP, length.out = 50),
    phase.original = .$phase.original,
    site = unique(mylu.working.imputated$site),
    season = "hiber_late"
  )) %>%ungroup()
nd4e$phat <- predict(latehibergdLmod, newdata = nd4e, type = "response", re.form = ~0)

nd4e <- nd4e %>%
  mutate(Temp_bin = cut(avgTEMP,breaks = quantile(avgTEMP, probs = seq(0, 1, length.out = 4), na.rm = TRUE),
      include.lowest = TRUE,labels = c("Cold", "Tepid", "Warm"))) %>%ungroup()

nd4e_bin <- nd4e %>%
  group_by(Temp_bin, avglogVPD) %>%
  summarise(phat = mean(phat), .groups = "drop")

late_fungalload_figure <- ggplot() +
  geom_point(
    data = mylu.working.avg,
    aes(x = avg_avglogVPD,y = avg_adj_gdL,color = site,size = n_bats),alpha = 0.6) +
  scale_color_viridis_d(option = "viridis") +
  scale_size_continuous(name = "Number of Bats Sampled") +
  guides(color = "none") +
  new_scale_color() +
  geom_line(data = nd4e_bin,aes(x = avglogVPD,y = phat,color = Temp_bin,
                                group = Temp_bin),linewidth = 2)+
  scale_color_viridis_d(option = "viridis", name = "Temperature") +
  guides(linetype = "none") +
   labs(y = expression("Transformed Fungal Load"),title = expression("Late Hibernation")) +
  scale_x_continuous(name = "log10 Vapor Pressure Deficit")+
  #scale_y_continuous(limits = c(0, 0.1))+
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

late_fungalload_figure




