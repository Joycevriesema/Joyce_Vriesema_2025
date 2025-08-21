rm(list = ls())

# load libraries
library(tidyverse)    # includes dplyr, ggplot2 and tidyr
library(lme4)         # for mixed models
library(lmerTest)     # for tests of significance of mixed-effects models
library(MASS)         # for negative binomial models
library(multcomp)     # registers cld() method for emmGrid
library(multcompView) # generates the letters
library(emmeans)      # pairwise comparison

# load bird data and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_bird <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=0&single=true&output=csv") |>
  dplyr::filter (!date %in% c("5-Feb-2025", "8-Feb-2025"))|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"),
         habitat_detailed = case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Papyrus Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4") ~ "Papyrus Robana Mouth",
           transect_ID %in% c("pap_5", "pap_6") ~ "Papyrus Robana Mid",
           transect_ID %in% c("tree_1", "tree_2") ~ "Tree Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Tree Mbalageti Far",
           transect_ID %in% c("tree_5", "tree_6") ~ "Tree Robana Mouth",
           transect_ID %in% c("tree_7", "tree_8") ~ "Tree Robana Far",
           TRUE ~ NA_character_),
         habitat_main = case_when(
           transect_ID %in% c("pap_1", "pap_2","pap_3", "pap_4","pap_5", "pap_6" ) ~ "Papyrus",
           transect_ID %in% c("tree_1", "tree_2","tree_3", "tree_4","tree_5", "tree_6","tree_7", "tree_8") ~ "Trees",
           TRUE ~ NA_character_),
         distance_to_river_mouth = case_when(
           transect_ID %in% c("pap_1", "pap_2", "pap_3", "pap_4", "tree_5", "tree_6") ~ "Mouth",
           transect_ID %in% c("tree_1", "tree_2","pap_5", "pap_6" ) ~ "Mid",
           transect_ID %in% c("tree_7", "tree_8","tree_3", "tree_4") ~ "Far"
         ),
         distance_to_river_mouth= factor(distance_to_river_mouth, levels= c("Mouth", "Mid", "Far")),
         habitat_detailed2= case_when(
           transect_ID %in% c("pap_1", "pap_2") ~ "Mbalageti Mouth",
           transect_ID %in% c("pap_3", "pap_4","tree_5", "tree_6") ~ "Robana Mouth",
           transect_ID %in% c("tree_7", "tree_8") ~ "Robana Far",
           transect_ID %in% c("pap_5", "pap_6") ~ "Robana Mid",
           transect_ID %in% c("tree_1", "tree_2") ~ "Mbalageti Mid",
           transect_ID %in% c("tree_3", "tree_4") ~ "Mbalageti Far",
           TRUE ~ NA_character_),
         river= case_when(
           transect_ID %in% c("pap_1", "pap_2", "tree_1", "tree_2","tree_3", "tree_4") ~ "Mbalageti",
           transect_ID %in% c("pap_3", "pap_4", "pap_5", "pap_6","tree_5", "tree_6", "tree_7", "tree_8") ~ "Robana"))|>
  group_by(transect_ID,river, habitat_detailed, habitat_detailed2, distance_to_river_mouth, habitat_main, bird_species_ID, date ) |>
  summarise(total_count= sum(count))|>
  filter(bird_species_ID=="Pied kingfisher")

# merge data with meta data
# load data transect
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv")|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))

# join the two dataframes
data_bird <- data_bird |>
  left_join(data_transect, by = c("transect_ID", "date"))

# filter out the extra observations (run_ID contains "extra"). These observations were taken after 12 PM. Not every transects has been sampled evenly after 12 PM so exclude these observations from the dataset --> 10 observations will be excluded 
data_bird <- data_bird |>
  filter(!coalesce(str_detect(run_ID, regex("extra", ignore_case = TRUE)), FALSE))

# make start_time_fish time a factor with 6 levels from 7 until 12 hour
# make time a factor with 6 levels from 7 until 12 hour
data_bird <- data_bird |>
  mutate(
    # extract a clean time like "9:29", "10:03", or "10:03:15"
    start_time_token = str_extract(as.character(start_time_bird1),
                                   "\\b\\d{1,2}:\\d{2}(?::\\d{2})?\\b"),
    # add seconds if missing so hms() accepts it
    start_time_fixed = if_else(
      is.na(start_time_token), NA_character_,
      if_else(str_detect(start_time_token, "^\\d{1,2}:\\d{2}$"),
              paste0(start_time_token, ":00"),
              start_time_token)
    ),
    start_time_hms = hms(start_time_fixed),
    # hour factor with six levels (7–12). Values outside 7–12 will become NA.
    start_hour_bird = factor(hour(start_time_hms), levels = 7:12, ordered = TRUE)
  )



# save as csv file --> using for Structural Equation Modelling (SEM)
write.csv(data_bird, "data_bird.csv", row.names = FALSE)

data_bird$distance_to_river_mouth <- factor(data_bird$distance_to_river_mouth, levels= c("Mouth", "Mid", "Far"))

#### boxplot count pied kingfishers/transect ####
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


#### count pied kingfisher corrected for visible tree volume along the shoreline ####

# join data_bird with data_trees, keep 94 observations of the bird df
data_trees <- read.csv("data_trees.csv")
data_bird_trees <- data_bird |>
  left_join(
    data_trees,
    by = c("transect_ID", "habitat_detailed", "habitat_main",
           "habitat_detailed2", "river", "distance_to_river_mouth")
  )

# calculate the number of birds per meter shoreline counted
# for papyrus this is the total number of bird counts / 500 m shoreline (the whole transect length)
# for tree this is the total number of bird counts / total visible circumference
birds_meter_shoreline <- data_bird_trees |>
  mutate(
    birds_per_meter = case_when(
      habitat_main == "Papyrus" ~ total_count / 500,
      habitat_main == "Trees" ~ total_count / total_visible_circumference,
      TRUE ~ NA_real_  # catch-all for other habitat types if any
    ),
    birds_per_100m = birds_per_meter * 100
  )

#### count pied kingfisher for different visible tree volumes ####
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


#### boxplot and statistics count pied kingfishers/500m shoreline ####
# make factors
birds_meter_shoreline <- birds_meter_shoreline |> 
  mutate(
    habitat_main = as.factor(habitat_main),
    river = as.factor(river),
    distance_to_river_mouth = factor(distance_to_river_mouth, levels = c("Mouth", "Mid", "Far")),  # set reference
    transect_ID = as.factor(transect_ID),
    weather= as.factor(weather))



# make model
m1 <- lmer(birds_per_100m ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = birds_meter_shoreline)

anova(m1)
# habitat type had a significant effect on bird densities ***
# river does not have an affect on bird densities
# distance to river mouth has a significant effect on bird densities **
# interaction between habitat main and distance to river mouth **

summary(m1)
# tree significant higher bird densities
# mid distance to river mouth has significant lower bird density than mouth
# far distance has not significantly higher bird densities than mouth, but there is a trend

# pairwise comparisons for m1
emm <- emmeans(m1, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm)

# create letters
cld_dist <- cld(emm, Letters = letters)
cld_dist
letters_df <- as.data.frame(cld_dist)
letters_df <- letters_df |>
  mutate(y_pos = 80) |>
  mutate(group_label = gsub(" ", "", .group)) |>
  ungroup() |>
  drop_na()
  

# plot with letters
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
    aes(x = distance_to_river_mouth, y = y_pos, label = group_label, group= river),
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


# model with transect_ID, weather, time and date as random effects
birds_meter_shoreline <- birds_meter_shoreline |>
  mutate(
  date_f  = factor(date))

m2 <- lmer(
  birds_per_100m ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID) + (1 | date_f) + (1|weather) + (1|start_time_bird1),
  data = birds_meter_shoreline, REML = TRUE
)

anova(m2)
summary(m2)
# all random effects are really small in their variance 
# start_time has the smallest variance of all and transect the largest

# pairwise comparisons for m1
emm2 <- emmeans(m2, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
pairs(emm2)

# create letters
cld_dist2 <- cld(emm2, Letters = letters)
cld_dist2
letters_df2 <- as.data.frame(cld_dist2)
letters_df2 <- letters_df2 |>
  mutate(y_pos = 80) |>
  mutate(group_label = gsub(" ", "", .group)) |>
  ungroup() |>
  drop_na()

ggplot(birds_meter_shoreline, aes(x = distance_to_river_mouth, y = birds_per_100m, fill = river)) +
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
    data = letters_df2,
    aes(x = distance_to_river_mouth, y = y_pos, label = group_label, group= river),
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