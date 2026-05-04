rm(list = ls())
library(emmeans)
library(effects)
library(performance)
library(dplyr)
library(ggplot2)
library(ggnewscale)
library(glmmTMB)
library(DHARMa)
library(AICcmodavg)
library(bbmle)
library(patchwork)


recap <- read.csv("recap")

recap$phase.original = as.factor(recap$phase.original)
recap$phase.original = relevel(recap$phase.original, ref="invasion")

recap_sub<- recap%>%
  filter(!site %in% c("MAIDEN ROCK", "BAY CITY","NEDA MINE", "SOUTH LAKE MINE")) #remove large sites with low recap probability - doesn't indicate survival

mylu.working.imputated <- read.csv("mylu.working.imputated")

mylu.working.imputated$phase.original = as.factor(mylu.working.imputated$phase.original)
mylu.working.imputated$phase.original = relevel(mylu.working.imputated$phase.original, ref="invasion")


##########################################################################
################ FIGURES #################################################
##########################################################################

#set theme
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



########################################
########## Fungal loads #########
########################################


mylu.working.imputated <- mylu.working.imputated %>%
  mutate(adj_gdL = gdL^(1/4))

latehibergdLmod= glmmTMB(adj_gdL ~ avglogVPD * avgTEMP * phase.original  +  (1|site),
                         data=subset(mylu.working.imputated, ysw<9 & ysw>-1& season=="hiber_late"&gd ==1),family  = Gamma(link="log"))

simulationOutput <- simulateResiduals(fittedModel = latehibergdLmod, plot = F)
plot(simulationOutput) #suspicious activities here

summary(latehibergdLmod)
mylu.working.avg <- mylu.working.imputated %>%  
  filter(season == "hiber_late"&!is.na(phase.original))%>%
  group_by(site, date,phase.original) %>%   
  summarise(avg_avglogVPD = mean(avglogVPD, na.rm = TRUE),avg_gdL_adj = mean(adj_gdL, na.rm = TRUE),
    n_bats= n(),.groups = "drop")

temp_quants <- mylu.working.imputated %>%
  filter(ysw < 9, ysw > -1, season == "hiber_late") %>%
  group_by(phase.original) %>%
  summarise(q25 = quantile(avgTEMP, 0.25, na.rm = TRUE),q50 = quantile(avgTEMP, 0.50, na.rm = TRUE),
    q75 = quantile(avgTEMP, 0.75, na.rm = TRUE),
    .groups = "drop")

temp_vals <- data.frame(
  phase.original = c(rep("established", 3), rep("invasion", 3)),
  avgTEMP = c(
    3.976564, 5.539557, 7.235314,   # established (Q25, Q50, Q75)
    3.231299, 5.264804, 6.222857    # invasion (Q25, Q50, Q75)
  ),quantile = rep(c("Cold", "Tepid", "Warm"), 2))

nd_temp <- temp_vals %>%
  rowwise() %>%
  do(expand.grid(
    avglogVPD = seq(min(mylu.working.imputated$avglogVPD, na.rm = TRUE),max(mylu.working.imputated$avglogVPD, na.rm = TRUE),
      length.out = 50),
    avgTEMP = .$avgTEMP,phase.original = .$phase.original,quantile = .$quantile,
    site = unique(mylu.working.imputated$site),season = "hiber_late")) %>%
  ungroup()

nd_temp$phat <- predict(
  latehibergdLmod,newdata = nd_temp,type = "response",re.form = ~0)

nd_temp <- nd_temp %>%
  filter(phase.original != "invasion" |quantile == "Tepid")

late_fungalload_figure <- ggplot() +
  geom_point(data = mylu.working.avg,aes(x = avg_avglogVPD,y = avg_gdL_adj,color = site,size = n_bats),
    alpha = 0.6) +
  scale_color_viridis_d(option = "magma", begin = 0.3, end = 0.9) +
  scale_size_continuous(name = "Number of Bats Sampled") +
  guides(color = "none") +
  new_scale_color() +
  geom_line(
    data = nd_temp,
    aes(x = avglogVPD,y = phat,color = quantile,group = interaction(quantile, phase.original)),
    linewidth = 2) +
  scale_color_viridis_d(option = "inferno", begin = 0.3, end = 0.9, name = "Temperature") +
  facet_wrap(~ phase.original,scales = "free_x",labeller = labeller(
    phase.original = c("invasion" = "Invasion (years 0-3)","established" = "Established (years 4-8)"))) +
  labs(y = expression("Fungal Load (fourth-root transformed)"),x = expression("Vapor Pressure Deficit"),
       title = expression("Effect of Vapor Pressure Deficit on Fungal Loads (Late Hibernation)")) +
  scale_x_continuous(
    limits = range(mylu.working.avg$avg_avglogVPD, na.rm = TRUE),
    breaks = seq(
      floor(min(nd_temp$avglogVPD, na.rm = TRUE)),
      ceiling(max(nd_temp$avglogVPD, na.rm = TRUE)),
      by = 0.5),
    minor_breaks = seq(
      floor(min(nd_temp$avglogVPD, na.rm = TRUE)),
      ceiling(max(nd_temp$avglogVPD, na.rm = TRUE)),
      by = 0.1),
    labels = scales::math_format(10^.x)) +
  scale_y_continuous(limits = c(0, 1.1))

late_fungalload_figure

########################################
########## Tissue Invasion #########
########################################

latehiberuvmod=glmmTMB(uv ~ avglogVPD * lgdL * phase.original+ (1|site), 
                       data = subset(mylu.working.imputated, ysw >= 0 & ysw < 9 
                                     & season == "hiber_late"&gd==1), family = binomial(link="probit"))
summary(latehiberuvmod)

mylu.working.avg <- mylu.working.imputated %>%  #tissue invasion averages from a survey
  filter(season == "hiber_late"&!is.na(phase.original))%>%
  group_by(site, date,phase.original) %>%   
  summarise(avg_avglogVPD = mean(avglogVPD, na.rm = TRUE),avg_uv= mean(uv, na.rm = TRUE),
    avg_lgdL= mean(lgdL, na.rm = TRUE),n_bats= n(),
    .groups = "drop")

vpd_quants <- mylu.working.imputated %>%
  filter(ysw < 9 & ysw > -1) %>%
  group_by(phase.original) %>%
  summarise(q25 = quantile(avglogVPD, 0.25, na.rm = TRUE),q50 = quantile(avglogVPD, 0.50, na.rm = TRUE),q75 = quantile(avglogVPD, 0.75, na.rm = TRUE),
            .groups = "drop")

vpd_vals <- data.frame(
  phase.original = c(rep("established", 3), rep("invasion", 3)),
  avglogVPD = c(
    -2.115524, -1.885091, -1.645639,   # established
    -2.082432, -1.874709, -1.537849    # invasion
  ),quantile = rep(c("Wet", "Damp", "Dry"), 2))

nd_quant <- vpd_vals %>%
  rowwise() %>%
  do(expand.grid(
    lgdL = seq(min(mylu.working.imputated$lgdL, na.rm = TRUE),max(mylu.working.imputated$lgdL, na.rm = TRUE),
      length.out = 50),
    avglogVPD = .$avglogVPD,phase.original = .$phase.original,quantile= .$quantile,
    site= unique(mylu.working.imputated$site))) %>%
  ungroup() %>%
  mutate(phat = predict(latehiberuvmod, newdata = ., type = "response", re.form = ~0))

nd_quant <- nd_quant %>%
  mutate(line_style = ifelse(phase.original == "established", "dashed", "solid"))

# Plot
late_uv_figure <- ggplot() +
  # Points: observed averages
  geom_point(
    data = mylu.working.avg,
    aes(x = avg_lgdL,y = avg_uv,color = avg_avglogVPD,size = n_bats),
    alpha = 0.6) +
  scale_color_viridis_c(option = "viridis") +
  scale_size_continuous(name = "Number of Bats Sampled") +
  guides(color = "none") +
  new_scale_color() +
  geom_line(
    data = nd_quant,aes(x = lgdL,y = phat,color = quantile,group = interaction(quantile, phase.original),linetype = line_style),
    linewidth = 2) +
  scale_color_viridis_d(option = "viridis", name = "Early VPD") +
  guides(linetype = "none") +
  facet_wrap(~ phase.original, scales = "free_x",labeller = labeller(phase.original = c("invasion" = "Invasion (years 0-3)", "established" = "Established (years 4-8)"))) +
  scale_x_continuous(limits = range(mylu.working.avg$avg_lgdL, na.rm = TRUE),
                     breaks = seq(floor(min(nd_quant$lgdL, na.rm = TRUE)),
                                  ceiling(max(nd_quant$lgdL, na.rm = TRUE)),
                                  by = 1),
                     minor_breaks = seq(floor(min(nd_quant$lgdL, na.rm = TRUE)),
                                        ceiling(max(nd_quant$lgdL, na.rm = TRUE)),
                                        by = 0.1),
                     labels = scales::math_format(10^.x)) +
  labs(
    y = "Probability of Tissue Invasion",
    x = "Fungal Load",
    title = "Effect of Fungal Loads on Tissue Invasion (Late Hibernation)") 

late_uv_figure

########################################
########## recapture/survival #########
########################################

recapAftermod <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_fungal_load + phase.original + (1|site),
                          data = subset(recap_sub, ysw < 9 & ysw > -1), family = binomial)
summary(recapAftermod)

# Observed averages
recap.avg <- recap %>%
  filter(ysw < 9, !is.na(phase.original)) %>%
  group_by(site, early_date, phase.original) %>%
  summarise(across(c(early_fungal_load, early_RecapturedAfter, early_avglogVPD),
           ~mean(.x, na.rm = TRUE),.names = "avg_{.col}"),n_bats = n(),
    .groups = "drop")

# Prediction grid + predictions + binning
nd_recap_bin <- recap %>%
  filter(ysw < 9) %>%
  group_by(phase.original) %>%
  summarise(
    across(c(early_fungal_load, early_avglogVPD),list(min = ~min(.x, na.rm = TRUE),max = ~max(.x, na.rm = TRUE))),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  do(expand.grid(early_fungal_load = seq(.$early_fungal_load_min, .$early_fungal_load_max, length.out = 50),
    early_avglogVPD = seq(.$early_avglogVPD_min, .$early_avglogVPD_max, length.out = 50),
    phase.original= .$phase.original,site = unique(recap$site))) %>%
  ungroup() 


vpd_quants <- recap %>%
  filter(ysw < 9 & ysw > -1) %>%
  group_by(phase.original) %>%
  summarise(q25 = quantile(early_avglogVPD, 0.25, na.rm = TRUE),q50 = quantile(early_avglogVPD, 0.50, na.rm = TRUE),
    q75 = quantile(early_avglogVPD, 0.75, na.rm = TRUE),
    .groups = "drop")

vpd_vals <- data.frame(
  phase.original = c(rep("established", 3), rep("invasion", 3)),
  early_avglogVPD = c(
    -2.115524, -1.87226, -1.545598,   # established
    -2.056272, -1.85326, -1.492172    # invasion
  ),
  quantile = rep(c("Wet", "Damp", "Dry"), 2))


nd_recap_quant <- vpd_vals %>%
  rowwise() %>%
  do(expand.grid(
    early_fungal_load = seq(min(recap$early_fungal_load, na.rm = TRUE),max(recap$early_fungal_load, na.rm = TRUE),
      length.out = 50),
    early_avglogVPD = .$early_avglogVPD,phase.original= .$phase.original,quantile= .$quantile,
    site= unique(recap$site)
  )) %>%
  ungroup() %>%
  mutate( phat = predict(recapAftermod, newdata = ., type = "response", re.form = ~0))

# Plot
recap_figure <- ggplot() +
  geom_point(
    data = recap.avg %>% filter(phase.original == "established"),
    aes(avg_early_fungal_load, avg_early_RecapturedAfter,color = avg_early_avglogVPD, size = n_bats),
    alpha = 0.6) +
  scale_color_viridis_c(option = "viridis", guide = "none") +
  scale_size_continuous(name = "Number of Bats Sampled") +
  new_scale_color() +
  geom_line(data = nd_recap_quant %>% filter(phase.original == "established"),
    aes(early_fungal_load,phat,color = quantile,group = interaction(quantile, phase.original)),linetype = "solid",linewidth = 2)+
  scale_color_viridis_d(option = "viridis",name = "Early VPD") +
  facet_wrap(~phase.original, scales = "free_x",labeller = labeller(phase.original = c(
        invasion = "Invasion (years 0–3)",
        established = "Established (years 4–8)"))) +
  scale_x_continuous(limits = range(recap$early_fungal_load, na.rm = TRUE),
                     breaks = seq(floor(min(nd_recap_bin$early_fungal_load, na.rm = TRUE)),
                                  ceiling(max(nd_recap_bin$early_fungal_load, na.rm = TRUE)),
                                  by = 1),
                     minor_breaks = seq(floor(min(nd_recap_bin$early_fungal_load, na.rm = TRUE)),
                                        ceiling(max(nd_recap_bin$early_fungal_load, na.rm = TRUE)),
                                        by = 0.1),
                     labels = scales::math_format(10^.x))+
  labs(x = "Early Fungal Load",y = "Probability of Recapture",title = "Effect of Early Fungal Load on Recapture Probability") 

recap_figure 



#### site climate visualization


site_summary <- mylu.working.imputated %>%
  group_by(site,phase.original) %>%
  summarise(mean_TEMP = mean(avgTEMP, na.rm = TRUE),mean_roosttemp = mean(temp, na.rm = TRUE),
    mean_logVPD  = mean(avglogVPD, na.rm = TRUE),n_obs     = n(),
    .groups = "drop")

site_summary <- site_summary%>%
  filter(!(is.na(phase.original)))

ggplot(site_summary, aes(x = mean_logVPD, y = mean_TEMP)) +
  geom_point(aes(size = n_obs, color = site), alpha = 0.8) +
  geom_text(aes(label = site, color = site),vjust = -0.6,size = 4,show.legend = FALSE) +
  scale_size_continuous(name = "Number of Observations") +
  facet_wrap(~ phase.original) +   # 👈 add this line
  labs(x = "Mean logVPD",y = "Mean Temperature",title = "Mean Temperature vs Mean logVPD by Site") +
  theme(plot.title = element_text(face = "bold", hjust = 0.5),legend.position = "none")


##########################################################################
################ MODEL COMPARISONS #######################################
##########################################################################


########################################
########## Fungal loads #########
########################################

df <- mylu.working.imputated %>%
  subset(ysw < 9 & ysw > -1 & season == "hiber_late"&gd ==1) %>%
  tidyr::drop_na(adj_gdL, avglogVPD, avgTEMP, phase.original, site)

null = glmmTMB(adj_gdL~1 + (1|site),data=df,family = Gamma(link="log")); summary(null)
m = glmmTMB(adj_gdL ~ avglogVPD +  (1|site),data=df,family = Gamma(link="log"))
summary(m)
m2 = glmmTMB(adj_gdL ~ avgTEMP +  (1|site),data=df,family = Gamma(link="log"))
summary(m2)
m3 = glmmTMB(adj_gdL ~ avglogVPD + avgTEMP +  (1|site),data=df,family = Gamma(link="log"))
m4 = glmmTMB(adj_gdL ~ avglogVPD * avgTEMP +  (1|site),data=df,family = Gamma(link="log"))
m5= glmmTMB(adj_gdL ~ avglogVPD + phase.original +  (1|site),data=df,family = Gamma(link="log"))
m6= glmmTMB(adj_gdL ~ avglogVPD * phase.original +  (1|site),data=df,family = Gamma(link="log"))
m7= glmmTMB(adj_gdL ~ avglogVPD + avgTEMP + phase.original +  (1|site),data=df,family = Gamma(link="log"))
m8= glmmTMB(adj_gdL ~ avglogVPD * avgTEMP + phase.original +  (1|site),data=df,family = Gamma(link="log"))
m9= glmmTMB(adj_gdL ~ avglogVPD + avgTEMP * phase.original +  (1|site),data=df,family = Gamma(link="log"))
m10= glmmTMB(adj_gdL ~ avglogVPD * phase.original + avgTEMP +  (1|site),data=df,family = Gamma(link="log"))
m11= glmmTMB(adj_gdL ~ avglogVPD * avgTEMP * phase.original +  (1|site),data=df,family = Gamma(link="log"))
m12=glmmTMB(adj_gdL ~ avgTEMP + phase.original +  (1|site),data=df,family = Gamma(link="log"))
m13=glmmTMB(adj_gdL ~ avgTEMP * phase.original +  (1|site),data=df,family = Gamma(link="log"))
AICtab(null,m,m2,m3,m4,m5,m6,m7,m8,m9,m10,m11,m12,m13,
       sort=TRUE, weights=TRUE)

#model 11 best by 4.4 AIC


########################################
########## Tissue Invasion #########
########################################


df<-mylu.working.imputated %>%
  subset(ysw < 9 & ysw > -1 & season == "hiber_late"&gd ==1) %>%
  tidyr::drop_na(lgdL, avglogVPD, avgTEMP, phase.original, site,lgdL)

m00=glmmTMB(uv ~ avgTEMP + (1|site), data = df, family = binomial(link="logit"))
m0=glmmTMB(uv ~ avglogVPD + (1|site), data = df, family = binomial(link="logit"))
m=glmmTMB(uv ~ avgTEMP + lgdL + (1|site), data = df, family = binomial(link="logit"))
mt=glmmTMB(uv ~ avgTEMP * lgdL + (1|site), data = df, family = binomial(link="logit"))
mt2=glmmTMB(uv ~ avgTEMP * lgdL+phase.original + (1|site), data = df, family = binomial(link="logit"))
mt3=glmmTMB(uv ~ avgTEMP * lgdL*phase.original + (1|site), data = df, family = binomial(link="logit"))
mod=glmmTMB(uv ~ avglogVPD*avgTEMP + lgdL + (1|site), data = df, family = binomial(link="logit"))
m1=glmmTMB(uv ~ avglogVPD + lgdL + (1|site), data = df, family = binomial(link="logit"))
m2=glmmTMB(uv ~ avglogVPD + lgdL + avgTEMP + (1|site), data = df, family = binomial(link="logit"))
m3=glmmTMB(uv ~ avglogVPD + lgdL + phase.original + (1|site), data = df, family = binomial(link="logit"))
m4=glmmTMB(uv ~ avglogVPD + lgdL + avgTEMP + phase.original + (1|site), data = df, family = binomial(link="logit"))
m5 <- glmmTMB(uv ~ avglogVPD * lgdL + (1|site), data = df, family = binomial(link="logit"))
m6 <- glmmTMB(uv ~ avglogVPD + lgdL * avgTEMP + (1|site),  data = df, family = binomial(link="logit"))
m7 <- glmmTMB(uv ~ avglogVPD + lgdL * phase.original + (1|site), data = df, family = binomial(link="logit"))
m8 <- glmmTMB(uv ~ avglogVPD + lgdL * avgTEMP * phase.original + (1|site),  data = df, family = binomial(link="logit"))
m9 <- glmmTMB(uv ~ avglogVPD * lgdL + avgTEMP + (1|site),  data = df, family = binomial(link="logit"))
m10 <- glmmTMB(uv ~ avglogVPD * lgdL + phase.original + (1|site), data = df, family = binomial(link="logit"))
m11 <- glmmTMB(uv ~ avglogVPD * lgdL + avgTEMP + phase.original + (1|site), data = df, family = binomial(link="logit"))
m12=glmmTMB(uv ~ lgdL + (1|site), data = df, family = binomial(link="logit"))
m13=glmmTMB(uv ~ avglogVPD * lgdL * avgTEMP + (1|site), data = df, family = binomial(link="logit"))
m14=glmmTMB(uv ~ avglogVPD * lgdL * phase.original + (1|site), data = df, family = binomial(link="logit"))
m15=glmmTMB(uv ~ avglogVPD * lgdL * avgTEMP * phase.original + (1|site), data = df, family = binomial(link="logit"))

AICtab(m0,m00,m,mt, mt2, mt3,mod,m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13, m14, m15,sort = TRUE,delta = TRUE, weights = TRUE)

#m14 best by 4.4 AIC

########################################
########## Recapture #########
########################################


df<-recap %>%
  subset(ysw < 9 & ysw > -1 &(!(is.na(early_gdL)))) %>%
  tidyr::drop_na(early_fungal_load, early_avglogVPD, early_avgTEMP, phase.original, site)

null=glmmTMB(early_RecapturedAfter ~ 1 + (1|site),data=df,family=binomial(link="logit"));summary(null)
m0<- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + (1|site),
              data = df, family = binomial(link="logit")); summary(m0)
m01<- glmmTMB( early_RecapturedAfter ~ early_avgTEMP + (1|site),
               data = df, family = binomial(link="logit")); summary(m01)
m02<- glmmTMB( early_RecapturedAfter ~ phase.original + (1|site),
               data = df, family = binomial(link="logit")); summary(m02)
m03<- glmmTMB( early_RecapturedAfter ~ early_fungal_load + (1|site),
               data = df, family = binomial(link="logit")); summary(m03)

m1 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP + (1|site),
               data = df, family = binomial(link="logit"))

m2 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP + (1|site),
               data = df, family = binomial(link="logit"))

m3 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + phase.original + (1|site),
               data = df, family = binomial(link="logit"))

m4 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * phase.original + (1|site),
               data = df, family = binomial(link="logit"));summary(m4)

m5 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_fungal_load + (1|site),
               data = df, family = binomial(link="logit"))

m6 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_fungal_load + (1|site),
               data = df, family = binomial(link="logit"))

m7 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP + phase.original + (1|site),
               data = df, family = binomial(link="logit"))

m8 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP + phase.original + (1|site),
               data = df, family = binomial(link="logit"))

m9 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * phase.original + early_avgTEMP + (1|site),
               data = df, family = binomial(link="logit"))

m10 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP * phase.original + (1|site),
                data = df, family = binomial(link="logit"))

m11 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP * phase.original + (1|site),
                data = df, family = binomial(link="logit"))
m12 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m13 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m14 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_fungal_load + early_avgTEMP + (1|site),
                data = df, family = binomial(link="logit"))

m15 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP * early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m16 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP * early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))
m17 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + phase.original + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m18 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * phase.original + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m19 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_fungal_load + phase.original + (1|site),
                data = df, family = binomial(link="logit"))

m20 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + phase.original * early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m21 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * phase.original * early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"));summary(m21)
m22 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP + phase.original + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m23 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP + phase.original + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m24 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * phase.original + early_avgTEMP + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m25 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_fungal_load + early_avgTEMP + phase.original + (1|site),
                data = df, family = binomial(link="logit"))

m26 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP * phase.original + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m27 <- glmmTMB(early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP * early_fungal_load + phase.original + (1|site),
               data = df, family = binomial(link="logit"))

m28 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + phase.original * early_fungal_load + early_avgTEMP + (1|site),
                data = df, family = binomial(link="logit"))

m29 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP * phase.original + early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m30 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP * early_fungal_load + phase.original + (1|site),
                data = df, family = binomial(link="logit"))

m31 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * phase.original * early_fungal_load + early_avgTEMP + (1|site),
                data = df, family = binomial(link="logit"))

m32 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD + early_avgTEMP * phase.original * early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m33 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * early_avgTEMP + phase.original * early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m34 <- glmmTMB( early_RecapturedAfter ~ early_avglogVPD * phase.original + early_avgTEMP * early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m35 <- glmmTMB( early_RecapturedAfter ~ early_avgTEMP +  early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"))

m36 <- glmmTMB( early_RecapturedAfter ~ early_avgTEMP *  early_fungal_load + (1|site),
                data = df, family = binomial(link="logit"));summary(m36)

AICtab(null,m0, m01, m02, m03,m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11,m12, m13, m14, m15, m16,
       m17, m18, m19, m20, m21,m22, m23, m24, m25, m26, m27, m28,m29, m30, m31, m32, m33, m34,m35, m36,
       sort = TRUE,delta = TRUE,weights = TRUE)

#m19 and 22 are within 2 AIC

#I will select m19 because it is a simpler model, and has the lowest AIC score by 0.7
