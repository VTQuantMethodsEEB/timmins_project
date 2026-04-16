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

mylu.working.imputated <- read.csv("mylu.working.imputated")

mylu.working.imputated$phase.original = as.factor(mylu.working.imputated$phase.original)
mylu.working.imputated$phase.original = relevel(mylu.working.imputated$phase.original, ref="invasion")

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

########## Tissue Invasion #########

latehiberuvmod=glmmTMB(uv ~ avglogVPD * lgdL * phase.original+ (1|site), 
                       data = subset(mylu.working.imputated, ysw >= 0 & ysw < 9 
                                     & season == "hiber_late"&gd==1), family = binomial(link="probit"))
summary(latehiberuvmod)

simulationOutput <- simulateResiduals(fittedModel = latehiberuvmod, plot = F)
plot(simulationOutput)


mylu.working.avg <- mylu.working.imputated %>%  #tissue invasion averages from a survey
  filter(season == "hiber_late"&!is.na(phase.original))%>%
  group_by(site, date,phase.original) %>%   
  summarise(
    avg_avglogVPD = mean(avglogVPD, na.rm = TRUE),
    avg_uv     = mean(uv, na.rm = TRUE),
    avg_lgdL     = mean(lgdL, na.rm = TRUE),
    n_bats        = n(),
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
    lgdL = seq(
      min(mylu.working.imputated$lgdL, na.rm = TRUE),
      max(mylu.working.imputated$lgdL, na.rm = TRUE),
      length.out = 50),
    avglogVPD = .$avglogVPD,
    phase.original  = .$phase.original,
    quantile        = .$quantile,
    site            = unique(mylu.working.imputated$site))) %>%
  ungroup() %>%
  mutate(phat = predict(latehiberuvmod, newdata = ., type = "response", re.form = ~0))

nd_quant <- nd_quant %>%
  mutate(
    line_style = ifelse(phase.original == "established", "dashed", "solid"))

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
    data = nd_quant,
    aes(x = lgdL,y = phat,color = quantile,
      group = interaction(quantile, phase.original),linetype = line_style),
    linewidth = 2) +
  scale_color_viridis_d(option = "viridis", name = "Early VPD") +
  guides(linetype = "none") +
  facet_wrap(~ phase.original, scales = "free_x",
             labeller = labeller(phase.original = c("invasion" = "Invasion (years 0-3)", "established" = "Established (years 4-8)"))) +
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


########## Fungal Loads #########

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
  summarise(
    avg_avglogVPD = mean(avglogVPD, na.rm = TRUE),
    avg_gdL_adj     = mean(adj_gdL, na.rm = TRUE),
    n_bats        = n(),
    .groups = "drop")

temp_quants <- mylu.working.imputated %>%
  filter(ysw < 9, ysw > -1, season == "hiber_late") %>%
  group_by(phase.original) %>%
  summarise(
    q25 = quantile(avgTEMP, 0.25, na.rm = TRUE),
    q50 = quantile(avgTEMP, 0.50, na.rm = TRUE),
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
    avglogVPD = seq(
      min(mylu.working.imputated$avglogVPD, na.rm = TRUE),
      max(mylu.working.imputated$avglogVPD, na.rm = TRUE),
      length.out = 50),
    avgTEMP = .$avgTEMP,
    phase.original = .$phase.original,
    quantile = .$quantile,
    site = unique(mylu.working.imputated$site),
    season = "hiber_late")) %>%
  ungroup()

nd_temp$phat <- predict(
  latehibergdLmod,newdata = nd_temp,type = "response",re.form = ~0)

nd_temp <- nd_temp %>%
  filter(phase.original != "invasion" |quantile == "Tepid")

late_fungalload_figure <- ggplot() +
  geom_point(
    data = mylu.working.avg,
    aes(x = avg_avglogVPD,y = avg_gdL_adj,color = site,size = n_bats),
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
