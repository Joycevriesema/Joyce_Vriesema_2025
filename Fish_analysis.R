rm(list = ls())

# load libraries
library(dplyr)
library(tidyverse)
library(stringr)
library(ggplot2)
library(tidyr)
library(vegan)
library(lme4)
library(emmeans)
library(glmmTMB) # package for negative binomial model including random effects
library(MASS) # this packacge I used for a negative binomial model, because of overdispersion

fish_data <-read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=946923967&single=true&output=csv") %>%
  mutate(observation = str_trim(observation),
         observation = case_when(
           observation %in% c("O O", "O.", "OO") ~ "O", # remove spaces and dots so that al O's belong into the same class
           observation == "1" ~ "S1", # same for 1 into the S1 class
           TRUE ~ observation
         )) %>%
  mutate(observation = ifelse(observation == "OV", "O,V", observation)) %>%
  separate_rows(observation, sep = ",") |>
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
         ),
         fish_type= case_when(
           observation %in% c("H","V","O","B")~"Large predatory fish",
           observation %in% c("S1","S2","S3")~"Small schooling fish"
         ))|>
  group_by(transect_ID,river, habitat_detailed, habitat_detailed_2, distance_to_river_mouth, habitat_main, date, fish_type ) |>
  summarise(total_count= n())

ggplot(fish_data |>dplyr::filter(fish_type == "Small schooling fish"), 
       aes( x= distance_to_river_mouth, y= total_count, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  scale_y_log10()+
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
  )+
  labs(x= "Distance to river mouth", y="School count / 500 m ")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))

ggplot(fish_data |>dplyr::filter(fish_type == "Large predatory fish"), 
       aes( x= distance_to_river_mouth, y= total_count, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  scale_y_log10()+
  labs(y="School count / 500 m ", title="Large predatory fish")

#### model for small schools (papyrus and trees included) ####
fish_schools <- fish_data |>
  dplyr:: filter(fish_type == "Small schooling fish")
fish_schools$river <- as.factor(fish_schools$river)
fish_schools$transect_ID <- as.factor(fish_schools$transect_ID)  


# count data with fixed effects of river site and distance to river mouth and random effect of transect ID included so taking a poisson model
glmm_schools <- glmer(total_count ~ river + distance_to_river_mouth + (1 | transect_ID),
                         data = fish_schools,
                         family = poisson)
summary(glmm_schools)
# significantly lower school count at Robana river 
# significantly higher school count at mid distance compared to mouth
# not significantly higher school count at far distance compared to mouth
# some variance between transects

# test the poisson model for overdispersion
overdisp_fun <- function(glmm_schools) {
  rdf <- df.residual(glmm_schools)
  rp <- residuals(glmm_schools, type = "pearson")
  Pearson.chisq <- sum(rp^2)
  c(chisq = Pearson.chisq,
    ratio = Pearson.chisq / rdf,
    df = rdf,
    p = pchisq(Pearson.chisq, df = rdf, lower.tail = FALSE))
}

overdisp_fun(glmm_schools)
# dispersion ratio = 8.2 is quite high and p<0.001 so this poisson model is highly overdispersed meaning that the variance in my real data is much higher than this poisson model can fit

# try negative binomial model
nb_schools <- glmmTMB(
  total_count ~ river + distance_to_river_mouth + (1 | transect_ID),
  data = fish_schools,
  family = nbinom2
)
summary(nb_schools)
# significant lower school counts at Robana river 
# significant higher school counts at mid distance compared to mouth
# no significant higher school counts 
# variance of transect_ID is really low so exclude from the model


# negative binomial model without random effect
nb_schools2 <- glm.nb(total_count ~ river + distance_to_river_mouth, data = fish_schools)
summary(nb_schools2)
# Significant lower school count at Robana
# significant higher at mid distance
# not significant higher at far distance

# pairwise comparison
emm <- emmeans(nb_schools2, ~ river + distance_to_river_mouth)
pairs(emm)

# get the significan letters
cld_dist <- multcomp::cld(emm, Letters = letters, adjust = "tukey")
letters_df <- as.data.frame(cld_dist)
letters_df <- letters_df %>%
  group_by(river, distance_to_river_mouth) %>%
  mutate(
    y_pos = max(fish_schools$total_count) * 1.05 + (row_number() - 1) * 0.5
  ) %>%
  ungroup()

# make plot with significance letters
ggplot(fish_data |>dplyr::filter(fish_type == "Small schooling fish"), 
         aes( x= distance_to_river_mouth, y= total_count, fill= river))+
  geom_boxplot(position = position_dodge(width = 0.75))+
  facet_grid(~habitat_main)+
  scale_y_log10() +
  theme(text= element_text(size=14))+
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
  )+
  geom_text(
    data = letters_df,
    aes(x = distance_to_river_mouth, y = y_pos, label = .group, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE
  )+
  labs(x= "Distance to river mouth", y="School count / 500 m ")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))

 

#### old model for small schools in papyrus transects ####
school_pap <- fish_data %>% 
  filter(habitat_main == "Papyrus", fish_type == "Small schooling fish")
school_pap$river <- as.factor(school_pap$river)
school_pap$transect_ID <- as.factor(school_pap$transect_ID)

# poisson model
glmm_school_pap <- glmer(total_count ~ river + distance_to_river_mouth + (1 | transect_ID),
                         data = school_pap,
                         family = poisson)
summary(glmm_school_pap)
# Robana lower school counts than Mbalageti
# far from mouth significantly higher counts than mouth
# boundary singular fit indicates that including transect as random effect has no function in the model
# scaled residuals are large suggesting overdispersion

# model without transect_ID as random effect
glm_pap <- glm(total_count ~ river + distance_to_river_mouth, family = poisson, data = school_pap)
summary(glm_pap)
# Robana significantly lower than Mbalageti
# far from river signifcantly higher school counts
# checking for overdispersion --> residual variance/ degrees of freedom 210.45/ 21 = 10.02. normally a value close to 1 is acceptable because this would indicate that the standard variance of the mean is 1 and meets the assumption for poisson model. But since the value is 10.02 there is overdispersion
# so try negative bionomial model

# binomial model
nb_papyrus <- glm.nb(total_count ~ river + distance_to_river_mouth, data = school_pap)
summary(nb_papyrus)
# Robana signifcantly lower number of schools than Mbalageti p<0.001
# far from mouth higher p<0.05
# residual variance is now 24.468 on 21 degrees of freedom, so this model looks better

# make plot
ggplot(school_pap, aes( x= distance_to_river_mouth, y= total_count, fill= river))+
  geom_boxplot()+
  theme(text= element_text(size=14))+
  scale_y_log10()+
  labs(y="School count / 500 m ", title="Papyrus")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))


#### old model for small schools in tree transects ####
school_tree <- fish_data|>
  filter(habitat_main == "Trees", fish_type == "Small schooling fish")

glmm_school_tree <- glmer(total_count ~ river + distance_to_river_mouth + (1 | transect_ID),
                          data = school_tree,
                          family = poisson)
summary(glmm_school_tree)
# no difference between Robana and Mbalageti 
# distance mid is signifcantly higher than mouth
# no difference between far from mouth and mouth
# variance from transect_ID is quite low so maybe not include in model
# overdispersion ratio seems quite high 525/36

# so try negative binomial model
nb_trees <- glmer.nb(total_count ~ river + distance_to_river_mouth + (1 | transect_ID), data = school_tree)
summary(nb_trees)
# no significant difference between Robana and Mbalageti river
# significant higher number of fish as mid distance compared to mouth
# no significant difference between far and mouth
# variance of transect_ID is 0 so suggest to leave out of the model

nb_trees2 <- glm.nb(total_count ~ river + distance_to_river_mouth, data = school_tree)
summary(nb_trees2)
# residual deviance is 43.923 on 37 degrees of freedom, seems better fit
# no significant difference between Robana en Mbalageti
# significant effect of river mouth and mid
# no significant effect of river mouth and far

# pairwise comparison for distance to river mouth
emm <- emmeans(nb_trees, pairwise ~ river * distance_to_river_mouth)
emm

# make plot
ggplot(school_tree, aes(x= distance_to_river_mouth, y= total_count, fill= river))+
  geom_boxplot()+
  theme(text= element_text(size=14))+
  scale_y_log10()+
  labs(y="School count / 500 m ", title="Trees")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9")
  )