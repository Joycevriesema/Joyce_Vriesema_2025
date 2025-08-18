rm(list = ls())

# load libraries
library(tidyverse)    # includes dplyr, ggplot2 and tidyr
library(lme4)         # for mixed models
library(lmerTest)     # for tests of significance of mixed-effects models
library(MASS)         # for negative binomial models
library(emmeans)      # pairwise comparison
library(patchwork)    # combine multiple ggplot2 plots into one layout
library(multcomp)     # registers cld() method for emmGrid
library(multcompView) # generates the letters

# load data_water
data_water <- read.csv("data_water.csv") |>
  dplyr::mutate(transect_date_time = paste(transect_ID, substr(date_time_end, 1, 6), substr(date_time_start, 12, 16), sep = " "))|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"),
         habitat_detailed = case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Papyrus Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4") ~ "Papyrus Robana Mouth",
           transect_ID %in% c("pap_5", "pap_6") ~ "Papyrus Robana Far",
           transect_ID %in% c("tree_1", "tree_2") ~ "Tree Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Tree Robana Mid",
           transect_ID %in% c("tree_5", "tree_6") ~ "Tree Robana Mouth",
           transect_ID %in% c("tree_7", "tree_8") ~ "Tree Robana Far",
           TRUE ~ NA_character_),
         habitat_main = case_when(
           transect_ID %in% c("pap_1", "pap_2","pap_3", "pap_4","pap_5", "pap_6" ) ~ "Papyrus",
           transect_ID %in% c("tree_1", "tree_2","tree_3", "tree_4","tree_5", "tree_6","tree_7", "tree_8") ~ "Trees",
           TRUE ~ NA_character_),
         distance_to_river_mouth= case_when(
           transect_ID %in% c("pap_1", "pap_2", "pap_3", "pap_4", "tree_5", "tree_6") ~ "Mouth",
           transect_ID %in% c("tree_1", "tree_2", "tree_3", "tree_4") ~ "Mid",
           transect_ID %in% c("pap_5", "pap_6", "tree_7", "tree_8") ~ "Far"
         ),
         distance_to_river_mouth= factor(distance_to_river_mouth, levels= c("Mouth", "Mid", "Far")),
         habitat_detailed_2= case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4","tree_5", "tree_6") ~ "Robana Mouth",
           transect_ID %in% c("pap_5", "pap_6","tree_7", "tree_8") ~ "Robana Far",
           transect_ID %in% c("tree_1", "tree_2") ~ "Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Robana Mid",
           TRUE ~ NA_character_),
         river= case_when(
           transect_ID %in% c("pap_1", "pap_2", "tree_1", "tree_2") ~ "Mbalageti",
           transect_ID %in% c("pap_3", "pap_4", "pap_5", "pap_6", "tree_3", "tree_4","tree_5", "tree_6", "tree_7", "tree_8") ~ "Robana"
         ))
# from PCA script it is clear there are three outliers, which are rows 377, 622 and 740 t=from the water dataframe 
data_water <- data_water |>
  dplyr::slice(c(-377,-622, -740))

#### model Temperature ####
m1 <- lmer(temp ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m1)
# only river has a significant effect
summary(m1)
# significant higher temp at Robana *
# variance of transect_ID is very low

# try model without random effect of transect_ID
m1a <- lm(temp ~ habitat_main * river * distance_to_river_mouth, data = data_water)
anova (m1a)
# river has significant effect ***
# distance to river mouth significant *
# interaction habitat:river **
# interaction habitat:distance to river marginally
summary(m1a)
# Robana river significantly higher temp
# interaction between papyrus and Mbalageti --> lower temp
# Mbalageti mouth pap differs significantly from pap mouth Robana and pap far Robana

# likelihood ratio test to compare the linear mixed model and the fixed-effects model
anova(m1, m1a)
# AIC values of m1a is slightly lower, chisq is 0 and p=1 so no evidence to keep the random effect in the model

# pairwise comparison
emm <- emmeans(m1a, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm)
# significant difference between Mbalageti mouth papyrus and Robana mouth papyrus
# significant difference between Mbalageti mouth papyrus and Robana far papyrus

# create letters
cld_dist <- cld(emm, Letters = letters)
cld_dist
letters_df <- as.data.frame(cld_dist)
letters_df <- letters_df |>
  mutate(y_pos = 33) |>
  mutate(group_label = gsub(" ", "", .group)) |>
  ungroup() |>
  drop_na()

plot1 <-ggplot(data_water, aes(x= distance_to_river_mouth, y=temp, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  labs(x= "Distance to river mouth", y="Temperature (°C) ", title="a)")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  geom_text(
    data = letters_df,
    aes(x = distance_to_river_mouth, y = y_pos, label = group_label, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE
  )+
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10) 
  )
plot1


#### model for pH ####
m2 <- lmer(pH ~ habitat_main * river *distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m2)
# significant effect of habitat *
# significant effect of river *
# significant effect of distance to river mouth ***
summary(m2)
# low variance of transect ID
# Robana river significantly higher pH 
# mid distance not significantly higher than mouth
# far distance significantly lower pH than mouth

# model wihtout random effects
m2a <- lm(pH ~ habitat_main * river * distance_to_river_mouth, data = data_water)
anova(m2a)
# river significant effect **
# distance to river mouth both significant effect ***
summary(m2a)
# Robana significantly higher pH
# far distance significanlty lower than mouth

# compare models with likelihood ratio test
anova(m2,m2a)
# p=1, AIC of m2a is slightly lower,so no evidence to keep the random effect in the model

# pairwise comparison
emm2 <- emmeans(m2a, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm2)
# significant difference between Mbalageti mouth pap and Robana mouth far
# significant difference between Robana mouth pap and Robana far papyrus
# differencne between Robana mouth trees and Robana far trees
# Mbalageti Mid Trees - Robana Mid Trees 0.0038

# create letters
cld_dist2 <- cld(emm2, Letters = letters)
cld_dist2
letters_df2 <- as.data.frame(cld_dist2)
letters_df2 <- letters_df2 |>
  mutate(y_pos = 10) |>
  mutate(group_label = gsub(" ", "", .group)) |>
  ungroup() |>
  drop_na()

plot2 <-ggplot(data_water, aes(x= distance_to_river_mouth, y=pH, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth",y="pH ", title="b)")+
  geom_text(
    data = letters_df2,
    aes(x = distance_to_river_mouth, y = y_pos, label = group_label, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE
  )+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
plot2

#### model HDO ####
# papyrus
m3 <- lmer(HDO ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova (m3)
# significant effect of habitat *
# significant effect of river *
# significant effect distance to river mouth ***
summary(m3)
# significant higher HDO at Robana
# margianlly lower HDO at mid distance
# significantly lower HDO at far distance 
# variance of transect_ID is quite low

# model wihtout random effects
m3a <- lm(HDO ~ habitat_main * river * distance_to_river_mouth, data = data_water)
anova(m3a)
# distance to river mouth significant effect
summary(m3a)

# compare models 
anova(m3,m3a)
# m3a lower AIC and p=1 so no evidence to keep the random effect in the model

# pairwise comparison
emm3 <- emmeans(m3a, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm3)

plot3 <-ggplot(data_water, aes(x= distance_to_river_mouth, y=HDO, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth", y="HDO (mg/L)", title="c)")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
plot3

#### model HDO sat ####
m4 <- lmer(HDO_sat ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m4)
# significant effect of river *
# distance to river mouth significant effect **
summary(m4)
# quite some variation in transect_ID

# model wihtout random effects
m4a <- lm(HDO_sat ~ river + distance_to_river_mouth, data = data_water)
summary(m4a)

# compare models 
anova(m4,m4a)
# lower AIC for model m4 and p<0.05 so keep random effect in the model

# pairwise comparison
emm4 <- emmeans(m4, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm4)

plot4 <-ggplot(data_water, aes(x= distance_to_river_mouth, y=HDO_sat, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth", y="HDO saturation (%) ", title="d)")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
plot4

#### model Turbidity ####
# papyrus
m5 <- lmer(turb ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m5)
# river has a significant effect
summary(m5)
# significant higher turb at Robana *
# quite some variation in transect_ID

# check model wihtout random effects
m5a <- lm(turb ~ habitat_main * river * distance_to_river_mouth, data = data_water)
summary(m5a)
# Robana significantly higher turbidity
# mid distance not significantly lower turbidity than mouth
# far distance not significantly lower turbidity than mouth

# compare models 
anova(m5,m5a)
# p=1, difference in AIC is very small, so both models could be a good fit

# pairwise comparison
emm5 <- emmeans(m5a, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm5)

plot5 <-ggplot(data_water |> filter(turb<2000),
               aes(x= distance_to_river_mouth, y=turb, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth", y="Turbditity (NTU) ", title="e)")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
plot5
#### model Conductivity ####
m6 <- lmer(cond ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m6)
# significant effect of river *
# significant effect of distance to river mouth **
# interaction effect between habitat and distance to river mouth *
summary(m6)
# significant higher conductivity at Robana
# significant lower conductivity at mid distance compared to mouth
# no significant lower conductivity at far distance
# quite some variance between transect but more within transect

# model without random effects
m6a <- lm(cond ~ habitat_main * river * distance_to_river_mouth, data = data_water)
anova(m6a)
summary(m6a)
# Robana highly significant
# mid distance highly significant lower conductivity 
# far distance not significantly higher than mouth

# compare models 
anova(m6,m6a)
# AIC of the mixed model is little bit lower
# p <0.001 so keep the random effect

# pairwise comparison
emm6 <- emmeans(m6, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm6)

plot6 <-ggplot(data_water, aes(x= distance_to_river_mouth, y=cond, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth", y="Conductivity (µS/cm) ", title="f)")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
plot6


#### model for TDS ####
m7 <- lmer(TDS ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m7)
# significant effect of river *
# significant effect of distance to river mouth **
# interaction effect of habitat and distance to river mouth
summary(m7)
# significant higher TDS at Robana
# significant lower TDS at mid distance than mouth
# no significant higher TDS at far distance 

# model without random effect
m7a <- lm(TDS ~ habitat_main * river * distance_to_river_mouth, data = data_water)
anova(m7a)
# everything highly significant 
summary(m7a)
# compare models
anova(m7,m7a)
# AIC model m7 is lower, chisq is 81.8 and p<0.001, so keep the random effect in the model

# pairwise comparison
emm7 <- emmeans(m7, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm7)

# plot TDS
plot7 <-ggplot(data_water,aes( x= distance_to_river_mouth, y= TDS, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth", y="TDS (mg/L) ", title="g)")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
plot7

ggplot(data_water |>dplyr::filter(TDS>50), 
       aes( x= distance_to_river_mouth, y= TDS, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  scale_y_log10()+
  labs(y="TDS (mg/L) ", title="TDS")+
  
  plot1
#### model for ORP ####
m8 <- lmer(ORP ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m8)
# significant effect of river
# significant effect of distance to river mouth
summary(m8)
# quite some variance between transects
# significant lower ORP at Robana
# no significant lower ORP at mid distance
# significant higher ORP at far distance

# model without random effect
m8a <- lm(ORP ~ habitat_main * river * distance_to_river_mouth, data = data_water)
anova(m8a)
# both river and distance to river mouth, highly significant
summary(m8a)
# significant lower ORP at Robana river
# significant higher ORP at distance mid
# significant higher ORP at distance far
anova(m8,m8a)
# p=0.16  and AIC values the same so keep the random effect in the model

plot8 <-ggplot(data_water,aes( x= distance_to_river_mouth, y= ORP, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth", y="ORP (mV) ", title="h)")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))+
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(size = 13, face="bold"),
    axis.title.x = element_text(size = 13, face = "bold", margin = margin(t = 10)),
    axis.text.x = element_text(size = 11, face = "bold", color="black"), 
    axis.text = element_text(size = 11),
    legend.position = "right",
    legend.title = element_text(size = 13, face = "bold"),
    strip.text = element_text(size = 13, face = "bold"),
    panel.grid.major.x = element_line(color = "grey80", linetype = "solid"),
    panel.grid.major.y = element_line(color = "grey80", linetype = "solid"),
    panel.grid.minor = element_blank(),
    panel.spacing = unit(1.2, "lines"),
    panel.border = element_rect(color = "black", size = 1, fill = NA),
    plot.margin = margin(10, 10, 10, 10)
  )
plot8

### make combined plot ####
# combine the plots into 1
combined_plot <- (
  plot1 | plot2 | plot3
) /
  (
    plot4 | plot5 | plot6
  ) /
  (
    plot7 | plot8 | plot_spacer()
  ) +
  plot_layout(
    widths = c(1, 1, 1),       
    heights = c(1, 1, 1),     
    guides = "collect"         
  ) &
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "bottom"
  )

print(combined_plot)