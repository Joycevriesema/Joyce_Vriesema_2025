rm(list = ls())

# load libraries
library(tidyverse)  # includes dplyr, ggplot2 and tidyr
library(lme4)       # for mixed models
library(lmerTest)   # for tests of significance of mixed-effects models
library(MASS)       # for negative binomial models
library(multcomp)   # registers cld() method for emmGrid
library(multcompView) # generates the letters
library(emmeans)    # pairwise comparison

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
  group_by(transect_ID,river, habitat_detailed, habitat_detailed_2, distance_to_river_mouth, habitat_main, bird_species_ID, date ) |>
  summarise(total_count= sum(count))|>
  filter(bird_species_ID=="Pied kingfisher")

# save as csv file --> using for Structural Equation Modelling (SEM)
write.csv(data_bird, "data_bird.csv", row.names = FALSE)

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


#### bird counts corrected for visible tree volume along the shoreline
# join data_bird with data_trees, keep 94 observations of the bird df
data_trees <- read.csv("data_trees.csv")
data_bird_trees <- left_join(data_bird, data_trees, by = c("transect_ID", "habitat_detailed", "habitat_main", "habitat_detailed_2", "river", "distance_to_river_mouth"))

# calculate the number of birds per meter shoreline counted
# for papyrus this is the total number of bird counts / 500 m shoreline (the whole transect length)
# for tree this is the total number of bird counts / total visible circumfence
birds_meter_shoreline <- data_bird_trees |>
  mutate(
    birds_per_meter = case_when(
      habitat_main == "Papyrus" ~ total_count / 500,
      habitat_main == "Trees" ~ total_count / total_visible_circumference,
      TRUE ~ NA_real_  # catch-all for other habitat types if any
    ),
    birds_per_100m = birds_per_meter * 100
  )

# plot tree volume for the different habitat combinations with symbols
ggplot(birds_meter_shoreline|> dplyr:: filter(habitat_main == "Trees"), aes(x = total_visible_tree_volume, y = total_count)) +
  geom_point(aes(color = habitat_detailed, shape = habitat_detailed), size = 3) +
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
  labs(
    x = "Visible tree volume (m³)",
    y = "Total bird count",
    color = "Habitat type",
    shape = "Habitat type"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black"),
    legend.position = "right"
  )


# plot bird counts per 100 meter shoreline
ggplot(birds_meter_shoreline, aes(x = distance_to_river_mouth, y = birds_per_100m, fill = river)) +
  geom_boxplot(
    color = "black",     
    outlier.shape = 21,    
    outlier.fill = "black",
    outlier.color = "black",
    outlier.size = 2
  ) +
  facet_grid(~ habitat_main) +
  #scale_y_log10() +
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"
  ))+
  labs(
    x = "Distance to river mouth",
    y = "Count / 100 meter shoreline",
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

#### model analysis birds per meter shoreline ####
# make factors 
birds_meter_shoreline <- birds_meter_shoreline |>
  mutate(
    habitat_main = as.factor(habitat_main),
    river = as.factor(river),
    distance_to_river_mouth = factor(distance_to_river_mouth, levels = c("Mouth", "Mid", "Far")),  # set reference
    transect_ID = as.factor(transect_ID)
  )
m1 <- lmer(birds_per_100m ~ habitat_main + river + distance_to_river_mouth + (1 | transect_ID), data = birds_meter_shoreline)
anova(m1)
# habitat type had a significant effect on bird densities **
# river site does not have an affect on bird densities
# distance to river mouth has a significant effect on bird densities *
summary(m1)
# tree significant higher bird densities
# mid distance to river mouth has significant lower bird density than mouth
# far distance has not significantly lower bird densities than mouth

# pairwise comparisons for m1
emm <- emmeans(m1, ~ river+ distance_to_river_mouth)
pairs(emm)
# mouth significantly higher than mid
# mouth not significantly higher than far
# mid significantly lower than far

# get letters
cld_dist <- cld(emm, Letters = letters, adjust = "tukey")
letters_df <- as.data.frame(cld_dist)
letters_df <- letters_df |>
  mutate(
    y_pos = max(birds_meter_shoreline$total_count) * 1.05 + (row_number() - 1) * 0.5
  ) %>%
  ungroup()


plot_birds_meter_shoreline <- ggplot(birds_meter_shoreline, aes(x = distance_to_river_mouth, y = birds_per_100m, fill = river)) +
  geom_boxplot(position = position_dodge(width = 0.75),
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
    y = "Count / 100 meter shoreline",
    fill = "River"
  ) +
  geom_text(
    data = letters_df,
    aes(x = distance_to_river_mouth, y = y_pos, label = .group, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE
  )+
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
plot_birds_meter_shoreline
# save the plot
ggsave("plotbirds per meter shoreline.png", plot_birds_meter_shoreline, width = 12, height = 6, dpi = 300)

#### old model analysis papyrus transects  with count data ####
# use Poisson regression model to test whether the total count of Pied kingfishers differs significantly between transects within each habitat type
# due to uneven number of the paired transect design it is only possible to compare within a habitat and not between habitats papyrus: n=6, trees: n=8

papyrus_data <- data_bird |> filter(habitat_main == "Papyrus")
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

#### old model analysis tree transects with count data ####

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



