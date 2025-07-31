rm(list = ls())

# load libraries
library(dplyr)
library(ggplot2)
library(lme4)
library(patchwork)
library(emmeans)
library(lmerTest)

# load data_water
data_water <- read.csv("data_water.csv") %>%
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


#### model Temperature ####
m1 <- lmer(temp ~ river + distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m1)
summary(m1)
# significant higher temp at Robana
# no significant higher temperature far than mouth
# no significant higher temperature at far distance compared to mouth
# variance of transect_ID is quite low

qqnorm(resid(m1)); qqline(resid(m1))
# residuals look good

# try model without random effect of transect_ID
m1a <- glm (temp ~ river + distance_to_river_mouth, data = data_water)
anova (m1a)
# river has an effect
# distance to river mouth has not
summary(m1a)
# Robana river significantly higher temp
# mid distance not significantly higher than mouth
# far not significantly higher than mouth

# likelihood ratio test to compare the linear mixed model and the fixed-effects model
anova(m1, m1a)
# AIC values of m1 is much lower, chisq is 8.1 and p<0.01 so this shows that transect_ID should be kept as an random factor in the model --> or is this maybe because i did not include habitat_main in the model

plot1 <-ggplot(data_water, aes(x= distance_to_river_mouth, y=temp, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth", y="Temperature (°C) ", title="a)")+
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
plot1


#### model for pH ####
m2 <- lmer(pH ~ river + distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m2)
# significant effect of river
# significant effect of distance to river mouth
summary(m2)
# low variance of transect ID
# Robana river significantly higher pH 
# mid distance not significantly higher than mouth
# far distance significantly lower pH than mouth

# model wihtout random effects
m2a <- glm(pH ~ river + distance_to_river_mouth, data = data_water)
anova(m2a)
# river and distance to river mouth both significant effect
summary(m2a)
# Robana significantly higher pH
# mid distance significantly higher than mouth
# far distance significanlty lower than mouth

# compare models with likelihood ratio test
anova(m2,m2a)
# p> 0.05 so the random effect can be removed from the model
plot2 <-ggplot(data_water, aes(x= distance_to_river_mouth, y=pH, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
  labs(x= "Distance to river mouth",y="pH ", title="b)")+
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
m3 <- lmer(HDO ~ river + distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova (m3)
# no significant effect of river
# significant effect distance to river mouth
summary(m3)
# no significant higher HDO at Robana
# no significant lower HDO at mid distance
# significantly lower HDO at far distance 
# variance of transect_ID is quite low

# model wihtout random effects
m3a <- glm(HDO ~ river + distance_to_river_mouth, data = data_water)
anova(m3a)
# river no significant effect
# distac=nce to river mouth significant effect
summary(m3a)
# Robana significantly higher HDO
# mid distance not significantly higher than mid
# far distance significantly lower than mouth

# compare models 
anova(m3,m3a)
# m3 lower AIC, p<0.05 so keep the random effect

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
m4 <- lmer(HDO_sat ~ river + distance_to_river_mouth + (1 | transect_ID), data = data_water)
anova(m4)
# no significant effect of river
# distance to river mouth significant effect
summary(m4)
# no signficant higher HDO at Robana
# no significant lower HDO at mid distance
# significant lower HDO at far distance
# quite some variation in transect_ID

# model wihtout random effects
m4a <- glm(HDO_sat ~ river + distance_to_river_mouth, data = data_water)
summary(m4a)
# Robana significantly higher HDO
# not significant lower HDO at mid distance
# significant lower HDO at far distance

# compare models 
anova(m4,m4a)
# lower AIC for model m4 and p<0.05 so keep random effect in the model

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
m5 <- lmer(turb ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_papyrus)
summary(m5)
# marginally significant higher turb at Robana
# no significant lower turbidity at far distance
# some variation in transect_ID

# model wihtout random effects
m5a <- glm(turb ~ river + distance_to_river_mouth, data = water_papyrus)
summary(m5a)
# Robana significantly higher turbidity
# far distance significantly lower turbidity

# compare models 
anova(m5,m5a)
# p= 0.15, difference in AIC= 0000.1, so both models could be a good fit

# tree
m5b <- lmer(turb ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_trees)
summary(m5b)
# Robana not significant higher turbidity
# mid distance has no significant lower turbidity
# far distance has no significant lower turbidity
# more variation within transects and not between transects

# model wihtout random effects
m5c <- glm(turb ~ river + distance_to_river_mouth, data = water_trees)
summary(m5c)
# Robana not significant higher turbidity
# mid distance has no significant lower turbidity
# far distance has no significant lower turbidity
anova(m5b,m5c)
# log likilihood is the same, p=1
# so simple model is a good choice

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
# papyrus
m6 <- lmer(cond ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_papyrus)
summary(m6)
# significant higher conductivity at Robana
# no significant lower conductivity at far distance
# some variance between transect but more within transect

# model wihtout random effects
m6a <- glm(cond ~ river + distance_to_river_mouth, data = water_papyrus)
summary(m6a)
# Robana significantly higher conductivity 
# far distance now suddenly significantly higher conductivity, is a bit suspicious

# compare models 
anova(m6,m6a)
# AIC of th emixed model is slighyl lower
# p= 0.07 so marginally significant that including the random effect improves the model
# keep the mixed model

# tree
m6b <- lmer(cond ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_trees)
summary(m6b)
# Robana not significant higher conductivity
# mid distance has significant lower conductivity
# far distance has no significant lower conductivity
# more variation within transects and not between transects

# model wihtout random effects
m6c <- glm(cond ~ river + distance_to_river_mouth, data = water_trees)
summary(m6c)
# Robana not significant higher conductivity
# mid distance has significant lower conductivity
# far distance has significant lower conductivity

anova(m6b,m6c)
# p <0.001 so significant effect shows a highly significant improvement when including the random effect.
# The AIC drops by 25 points, which is a large improvement (ΔAIC > 10 is usually considered strong evidence)

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
# papyrus 
m7 <- lmer(TDS ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_papyrus)
summary(m7)
# some variance between transects
# more variance within transects
# significant higher TDS at Robana
# no significant difference between mouth and far

# model without random effect
m7a <- glm(TDS ~ river + distance_to_river_mouth, data = water_papyrus)
summary(m7a)
# significant higher TDS at Robana river
# significant higher TDS at distance far
anova(m7,m7a)
# p= 0.08016, just above 0.05
# the mixed model has a slightly lower AIC, so better
# but the p value 0.08 suggests that adding the random effect only marginally improves the model 

# tree
m7b <- lmer(TDS ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_trees)
summary(m7b)
# significantly lower at mid-distance 
# no difference between mouth and far
# no difference between Mbalageti and Robana

# model without random effect
m7c <- glm(TDS ~ river + distance_to_river_mouth, data = water_trees)
summary(m7c)
anova(m7b,m7c)
# p value highly significant so mixed model is better

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
# papyrus 
m8 <- lmer(ORP ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_papyrus)
summary(m8)
# some variance between transects
# more variance within transects
# significant lower ORP at Robana
# no significant difference between mouth and far

# model without random effect
m8a <- glm(ORP ~ river + distance_to_river_mouth, data = water_papyrus)
summary(m8a)
# significant lower ORP at Robana river
# significant higher ORP at distance far
anova(m8,m8a)
# p= 0.6
# so simple model is a good fit

# tree
m8b <- lmer(ORP ~ river + distance_to_river_mouth + (1 | transect_ID), data = water_trees)
summary(m8b)
# Robana not significantly lower that Mbalageti
# not significantly lower at mid-distance 
# significantly higher at far distance 
# some variance between transects
# more variance within transects


# model without random effect
m8c <- glm(ORP ~ river + distance_to_river_mouth, data = water_trees)
summary(m8c)
# Robana significantly lower
# mid distance significantly lower
# far distance significantly higher 

anova(m8b,m8c)
# p=0.33
# so simple model would be a good fit

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
    widths = c(1, 1, 1),       # equal column widths
    heights = c(1, 1, 1),      # equal row heights
    guides = "collect"         # 👈 collect legend from all plots
  ) &
  theme(
    plot.margin = margin(5, 5, 5, 5),
    legend.position = "bottom"
  )

print(combined_plot)