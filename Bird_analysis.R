rm(list = ls())

# load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(lme4)
library(MASS)
library(emmeans)


# load bird data and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_bird <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=0&single=true&output=csv") |>
  dplyr::filter (!date %in% c("5-Feb-2025", "8-Feb-2025"))|>
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
         #habitat_main = factor(habitat_main, levels = c("Papyrus", "Trees")),
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
         )) |>
  group_by(transect_ID,river, habitat_detailed, habitat_detailed_2, distance_to_river_mouth, habitat_main, species_name, date ) |>
  summarise(total_count= sum(count))|>
  filter(species_name=="Pied kingfisher")

# old exploratory graphs
ggplot(data_bird, aes(x= as.factor(date), y= total_count, fill= transect_ID))+
  geom_bar(stat="identity", position = "dodge")

ggplot(data_bird, aes(x= transect_ID, y= total_count, fill= transect_ID))+
  geom_boxplot()

# boxplot pied kingfishers
ggplot(data_bird, aes(x = distance_to_river_mouth, y = total_count, fill = river)) +
  geom_boxplot(
    color = "black",    
    outlier.shape = 21,    
    outlier.fill = "black",
    outlier.color = "black",
    outlier.size = 2
  ) +
  facet_grid(~ habitat_main) +
  scale_y_log10() +
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"
  ))+
  labs(
    x = "Distance to river mouth",
    y = "Count / 500 m",
    fill = "River"
  ) +
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

#### model analysis papyrus transects ####
# use Poisson regression model to test whether the total count of Pied kingfishers differs significantly between transects within each habitat type
# due to uneven number of the paired transect design it is only possible to compare within a habitat and not between habitats papyrus: n=6, trees: n=8

papyrus_data <- data_bird %>% filter(habitat_main == "Papyrus")
papyrus_data$river <- as.factor(papyrus_data$river)
papyrus_data$distance_to_river_mouth <- as.factor(papyrus_data$distance_to_river_mouth)
papyrus_data$transect_ID <- as.factor(papyrus_data$transect_ID)

glmm_papyrus <- glmer(total_count ~ river + distance_to_river_mouth + (1 | transect_ID),
                      data = papyrus_data,
                      family = poisson)
summary(glmm_papyrus)
# distance far from river mouth has a significantly higher number of pied kingfishers than mouth for papyrus
# p***
# overdispersion ratio seems quite high (deviance / residuals degrees of freedom)--> 373/35=10.66
# therefore try negative binomial model

# negative binomial model 
nb_papyrus <- glmer.nb(total_count ~ river + distance_to_river_mouth + (1 | transect_ID), data = papyrus_data)
summary(nb_papyrus)
# Robana not significant lower
# distance far singficant higher than mouth
# warning in model "boundary singular" --> variance fro the random effect of transect_ID is zero
# therefore remove this from the model, otherwise the model is overparameterized

# new model without transect_ID as random factor
nb_papyrus2 <- glm.nb(total_count ~ river + distance_to_river_mouth, data = papyrus_data)
summary(nb_papyrus2)
# now the residual variance is 38.843 and degrees of freedom 36 so ratio is better now
# distance far p***
# river Robana not significant

#### model analysis tree transects ####

# add visible tree volume
# load data of tree measurements
data_trees <- read.csv("data_bird_trees_only.csv")
data_trees$river <- as.factor(data_trees$river)
data_trees$distance_to_river_mouth <- as.factor(data_trees$distance_to_river_mouth)
data_trees$transect_ID <- as.factor(data_trees$transect_ID)

# use poisson model and include visible tree volume as parameter
glmm_trees <- glmer(total_count ~ river + distance_to_river_mouth + total_visible_tree_volume + (1 | transect_ID),
                    data = data_trees,
                    family = poisson) # gives a warning that some predictors are on different scales, consider rescaling
summary(glmm_trees)

# rescaling visivle tree volume
data_trees$total_visible_tree_volume_scaled <- scale(data_trees$total_visible_tree_volume)

glmm_trees <- glmer(total_count ~ river + distance_to_river_mouth + total_visible_tree_volume_scaled + (1 | transect_ID),
                    data = data_trees,
                    family = poisson)
summary(glmm_trees)
# overdispersion ratio seems quite high 448.9/49=9.16
# so try a negative binomial model
data_trees$distance_to_river_mouth <- relevel(data_trees$distance_to_river_mouth, ref = "Mouth")
nb_trees <- glmer.nb(total_count ~ river + distance_to_river_mouth + total_visible_tree_volume_scaled + (1 | transect_ID), data = data_trees)
summary(nb_trees)
# seems some variation among transects
# median of the scaled residuals is close to zero, so this seems a better model
# visible tree volume marginally signifcant 

# model without tree volume
nb_trees2 <- glmer.nb(total_count ~ river + distance_to_river_mouth + (1 | transect_ID), data = data_trees)
summary(nb_trees2)

# test which model is better
anova(nb_trees, nb_trees2, test = "LRT")
# AIC of the model with visible tree volume is little bit lower so could suggest this simpler model, although it is not significant

# pairwise comparison to see which combination of distance river mouth mid is different
emm <- emmeans(nb_trees, pairwise ~ river * distance_to_river_mouth)
emm
# Robana mouth signifcantly different from Mbalageti mid
# Robana far significantly different from Mbalageti mid


