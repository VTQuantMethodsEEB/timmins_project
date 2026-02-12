
rm(list=ls())

#mylu.working.extrapolated <- read.csv("~/Dropbox/ellie_projects/classes/Kate_class/timmins_project/mylu.working.extrapolated")

read.csv("mylu.working.extrapolated") 

library(tidyverse)
library(glmmTMB)
library(bbmle) 
library(effects)
library(ggnewscale)

mylu.working.extrapolated$phase.original = as.factor(mylu.working.extrapolated$phase.original)
mylu.working.extrapolated$phase.original = relevel(mylu.working.extrapolated$phase.original, ref="invasion")

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

#########################################
##### using VPD from early hiber on UV ####
#####################################


uvlatemod<-glmmTMB(uv ~ avg_early_logVPD * lgdL * phase.original + (1|site), data = subset(mylu.working.extrapolated, ysw >= 0 & ysw < 9 & season == "hiber_late"), family = binomial)
summary(uvlatemod)

mylu.working.avg <- mylu.working.extrapolated %>%  #tissue invasion averages from a survey
  filter(season == "hiber_late"&!is.na(phase.original))%>%
  group_by(site, date,phase.original) %>%   
  summarise(
    avg_avglogVPD = mean(avg_early_logVPD, na.rm = TRUE),
    avg_uv     = mean(uv, na.rm = TRUE),
    avg_lgdL     = mean(lgdL, na.rm = TRUE),
    n_bats        = n(),
    .groups = "drop")

nd4e <- mylu.working.extrapolated %>%
  filter(ysw >= 0 & ysw < 9, season == "hiber_late") %>%
  group_by(phase.original) %>%
  summarise(
    minVPD = min(avg_early_logVPD, na.rm = TRUE),
    maxVPD = max(avg_early_logVPD, na.rm = TRUE),
    minLgdL = min(lgdL, na.rm = TRUE),
    maxLgdL = max(lgdL, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  do(expand.grid(
    avg_early_logVPD = seq(.$minVPD, .$maxVPD, length.out = 50),
    lgdL  = seq(.$minLgdL, .$maxLgdL, length.out = 50),
    phase.original = .$phase.original,
    site = unique(mylu.working.extrapolated$site),
    season = "hiber_late"
  )) %>%
  ungroup()
nd4e$phat <- predict(uvlatemod, newdata = nd4e, type = "response", re.form = ~0)

# Create bins for avg_early_logVPD in nd4e
nd4e <- nd4e %>%
  mutate(
    VPD_bin = cut(avg_early_logVPD,
                  breaks = quantile(avg_early_logVPD, probs = seq(0,1, length.out = 4), na.rm = TRUE),
                  include.lowest = TRUE,
                  labels = c("Wet", "Damp", "Dry"))
  )

# Aggregate predictions per bin for smooth lines
nd4e_bin <- nd4e %>%
  group_by(VPD_bin, lgdL, phase.original) %>%
  summarise(phat = mean(phat), .groups = "drop") %>%
  mutate(line_style = ifelse(phase.original == "invasion", "longdash", "solid"))

# Plot
late_uv_figure <- ggplot() +
  # Points: observed averages
  geom_point(
    data = mylu.working.avg,
    aes(
      x = avg_lgdL,
      y = avg_uv,
      color = site,
      size = n_bats
    ),
    alpha = 0.6
  ) +
  scale_color_viridis_d(option = "viridis") +
  scale_size_continuous(name = "Number of Bats Sampled") +
  guides(color = "none") +
  new_scale_color() +
  # Lines: predicted phat by VPD bin
  geom_line(
    data = nd4e_bin,
    aes(
      x = lgdL,
      y = phat,
      color = VPD_bin,
      group = interaction(VPD_bin, phase.original),
      linetype = line_style
    ),
    linewidth = 2
  ) +
  scale_color_viridis_d(option = "viridis", name = "Early VPD") +
  guides(linetype = "none") +
  facet_wrap(~ phase.original, scales = "free_x",
             labeller = labeller(phase.original = c("invasion" = "Invasion (years 0-3)", "established" = "Established (years 4-8)"))) +
  labs(
    y = "Probability of Tissue Invasion",
    x = expression(log[10] ~ Fungal ~ Load),
    title = "Effect of Fungal Loads on Tissue Invasion (Late Hibernation)"
  ) 

late_uv_figure

#########################################
###### using VPD from late hiber on UV ####
#####################################

#scale_x_log10 for logging

mylu.working.extrapolated$phase.original = as.factor(mylu.working.extrapolated$phase.original)
mylu.working.extrapolated$phase.original = relevel(mylu.working.extrapolated$phase.original, ref="invasion")

uvlatemod_lateVPD<-glmmTMB(uv ~ avglogVPD * lgdL * phase.original + (1|site), data = subset(mylu.working.extrapolated, ysw >= 0 & ysw < 9 & season == "hiber_late"), family = binomial)
summary(uvlatemod_lateVPD)

mylu.working.avg <- mylu.working.extrapolated %>%  #tissue invasion averages from a survey
  filter(season == "hiber_late"&!is.na(phase.original))%>%
  group_by(site, date,phase.original) %>%   
  summarise(
    avg_avglogVPD = mean(avglogVPD, na.rm = TRUE),
    avg_uv     = mean(uv, na.rm = TRUE),
    avg_lgdL     = mean(lgdL, na.rm = TRUE),
    n_bats        = n(),
    .groups = "drop")

nd4e <- mylu.working.extrapolated %>%
  filter(ysw >= 0 & ysw < 9, season == "hiber_late") %>%
  group_by(phase.original) %>%
  summarise(
    minVPD = min(avglogVPD, na.rm = TRUE),
    maxVPD = max(avglogVPD, na.rm = TRUE),
    minLgdL = min(lgdL, na.rm = TRUE),
    maxLgdL = max(lgdL, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  do(expand.grid(
    avglogVPD = seq(.$minVPD, .$maxVPD, length.out = 50),
    lgdL  = seq(.$minLgdL, .$maxLgdL, length.out = 50),
    phase.original = .$phase.original,
    site = unique(mylu.working.extrapolated$site),
    season = "hiber_late"
  )) %>%
  ungroup()
nd4e$phat <- predict(uvlatemod_lateVPD, newdata = nd4e, type = "response", re.form = ~0)

# Create bins for avglogVPD in nd4e
nd4e <- nd4e %>%
  mutate(
    VPD_bin = cut(avglogVPD,
                  breaks = quantile(avglogVPD, probs = seq(0,1, length.out = 4), na.rm = TRUE),
                  include.lowest = TRUE,
                  labels = c("Wet", "Damp", "Dry"))
  )

# Aggregate predictions per bin for smooth lines
nd4e_bin <- nd4e %>%
  group_by(VPD_bin, lgdL, phase.original) %>%
  summarise(phat = mean(phat), .groups = "drop") %>%
  mutate(line_style = ifelse(phase.original == "established", "longdash", "solid"))

# Plot
late_uv_figure <- ggplot() +
  # Points: observed averages
  geom_point(
    data = mylu.working.avg,
    aes(
      x = avg_lgdL,
      y = avg_uv,
      color = site,
      size = n_bats
    ),
    alpha = 0.6
  ) +
  scale_color_viridis_d(option = "viridis") +
  scale_size_continuous(name = "Number of Bats Sampled") +
  guides(color = "none") +
  new_scale_color() +
  # Lines: predicted phat by VPD bin
  geom_line(
    data = nd4e_bin,
    aes(
      x = lgdL,
      y = phat,
      color = VPD_bin,
      group = interaction(VPD_bin, phase.original),
      linetype = line_style
    ),
    linewidth = 2
  ) +
  scale_color_viridis_d(option = "viridis", name = "Late VPD") +
  guides(linetype = "none") +
  facet_wrap(~ phase.original, scales = "free_x",
             labeller = labeller(phase.original = c("invasion" = "Invasion (years 0-3)", "established" = "Established (years 4-8)"))) +
  labs(
    y = "Probability of Tissue Invasion",
    x = expression(log[10] ~ Fungal ~ Load),
    title = "Effect of Fungal Loads on Tissue Invasion (Late Hibernation)"
  ) 

late_uv_figure


###################################################################
###### using VPD from early hiber on late fungal loads ############
###################################################################

mylu.working.extrapolated$phase.original = as.factor(mylu.working.extrapolated$phase.original)
mylu.working.extrapolated$phase.original = relevel(mylu.working.extrapolated$phase.original, ref="invasion")

latefungalloads<-glmmTMB(lgdL ~ avg_early_logVPD * avgTEMP * phase.original +  (1|site),data=subset(mylu.working.extrapolated, ysw<9 & ysw>-1& season=="hiber_late"),family = gaussian(link = "identity"))
summary(latefungalloads)
#VPD and temp both significant in invasion phase, marginally sig in established phase

#### with temp on x axis #####

mylu.working.avg <- mylu.working.extrapolated %>%  #tissue invasion averages from a survey
  filter(season == "hiber_late"&!is.na(phase.original))%>%
  group_by(site, date,phase.original) %>%   
  summarise(
    avg_avglogVPD = mean(avg_early_logVPD, na.rm = TRUE),
    avg_avgTEMP    = mean(avgTEMP, na.rm = TRUE),
    avg_lgdL     = mean(lgdL, na.rm = TRUE),
    n_bats        = n(),
    .groups = "drop")

nd4e <- mylu.working.extrapolated %>%
  filter(ysw >= 0 & ysw < 9, season == "hiber_late") %>%
  group_by(phase.original) %>%
  summarise(
    minVPD = min(avg_early_logVPD, na.rm = TRUE),
    maxVPD = max(avg_early_logVPD, na.rm = TRUE),
    minTEMP = min(avgTEMP, na.rm = TRUE),
    maxTEMP= max(avgTEMP, na.rm = TRUE)
  ) %>%
  rowwise() %>%
  do(expand.grid(
    avg_early_logVPD = seq(.$minVPD, .$maxVPD, length.out = 50),
   avgTEMP = seq(.$minTEMP, .$maxTEMP, length.out = 50),
    phase.original = .$phase.original,
    site = unique(mylu.working.extrapolated$site),
    season = "hiber_late"
  )) %>%
  ungroup()
nd4e$phat <- predict(latefungalloads, newdata = nd4e, type = "response", re.form = ~0)

nd4e <- nd4e %>%
  group_by(phase.original) %>%
  mutate(
    VPD_bin = cut(
      avg_early_logVPD,
      breaks = quantile(avg_early_logVPD, probs = seq(0, 1, length.out = 4), na.rm = TRUE),
      include.lowest = TRUE,
      labels = c("Wet", "Damp", "Dry")
    )
  ) %>%
  ungroup()
nd4e_bin <- nd4e %>%
  group_by(phase.original, VPD_bin, avgTEMP) %>%
  summarise(phat = mean(phat), .groups = "drop")

late_fungalload_figure <- ggplot() +
  geom_point(
    data = mylu.working.avg,
    aes(
      x = avg_avgTEMP,
      y = avg_lgdL,
      color = site,
      size = n_bats
    ),
    alpha = 0.6
  ) +
  scale_color_viridis_d(option = "viridis") +
  scale_size_continuous(name = "Number of Bats Sampled") +
  guides(color = "none") +
  new_scale_color() +
  geom_line(
    data = nd4e_bin,
    aes(
      x = avgTEMP,
      y = phat,
      color = VPD_bin,
      group = VPD_bin
    ),
    linewidth = 2
  )+
  scale_color_viridis_d(option = "viridis", name = "Early VPD") +
  guides(linetype = "none") +
  facet_wrap(~ phase.original, scales = "free_x",
             labeller = labeller(phase.original = c("invasion" = "Invasion (years 0-3)", "established" = "Established (years 4-8)"))) +
  labs(
    y = expression(log[10] ~ Fungal ~ Load),
    x = "Average Temperature",
    title = "Effect of Temperature on Fungal Loads (Late Hibernation)"
  ) 

late_fungalload_figure

