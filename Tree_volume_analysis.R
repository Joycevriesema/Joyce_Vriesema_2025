rm(list = ls())

# load libraries
library(dplyr)
library(ggplot2)
library(tidyr)
library(ggrepel)
library(emmeans)
library(patchwork)

# load tree data
data_trees <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1205571162&single=true&output=csv") |>
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
         )) |>
  group_by(transect_ID,river, habitat_detailed, habitat_detailed_2, distance_to_river_mouth, habitat_main )

# calculate tree volume (as ellyptic cylinder) for every transect
data_trees$tree_volume <- pi * (data_trees$canopy_length / 2) * (data_trees$canopy_depth / 2) * data_trees$tree_heigth

# calculate visible tree volume = tree volume x visibility for every transect
data_trees <- data_trees |>
  mutate(visible_tree_volume = tree_volume * (tree_visibility / 100))

# calculate the total canopy length per transect
data_trees <- data_trees |>
  mutate(total_canopy_length = sum(canopy_length))

# sum tree volumes an visible tree volumes per transect, so that every transect has one value of visible tree volume for that transect
data_trees <- data_trees %>%
  #group_by(transect_ID) %>%
  summarise(total_tree_volume = sum(tree_volume),
            total_visible_tree_volume = sum(visible_tree_volume),
            total_canopy_length = sum(canopy_length))


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

# join data_bird with data_trees, keep 94 observations of the bird df
data_bird_trees <- left_join(data_bird, data_trees, by = c("transect_ID", "habitat_detailed", "habitat_main", "habitat_detailed_2", "river", "distance_to_river_mouth"))


# plot bird counts per tree volume
data_bird_trees_only <- data_bird_trees |> 
  filter(habitat_main == "Trees")
ggplot(data_bird_trees_only, aes(x = total_visible_tree_volume, y = total_count)) +
  geom_point(size = 2, color = "#0072B2")+
  geom_smooth(method = "lm", se = TRUE, color = "black", linetype = "solid") +
  #geom_text_repel(aes(label = transect_ID), size = 4, fontface = "bold") +
  labs(
    x = "Visible tree volume (m³)",
    y = "Total bird count"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    axis.title = element_text(face = "bold"),
    axis.text = element_text(color = "black")
  )

# plot with symbols for the different habitat combinations
ggplot(data_bird_trees_only, aes(x = total_visible_tree_volume, y = total_count)) +
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


write.csv(data_bird_trees_only, "data_bird_trees_only.csv", row.names = FALSE)


#### bird counts per meter shoreline ####
# to standardize bird counts along transect and to compare tree transects with papyrus transects
# for papyrus transects it is always the count per 500m transect
# but for tree it is actually the bird count per the total meter of tree diameter counted
# for example: diameter tree = 5 meter, and 50 trees --> 50 * 5= 250 meter shoreline counted, papyrus is 500 meter so for trees this would be 250 * 2= amount per 500 meter
birds_m_shoreline <- data_bird_trees |>
  mutate(
    birds_per_meter = case_when(
      habitat_main == "Papyrus" ~ total_count / 500,
      habitat_main == "Trees" ~ total_count / total_canopy_length,
      TRUE ~ NA_real_  # catch-all for other habitat types if any
    ),
    birds_per_100m = birds_per_meter * 100
  )

# plot
ggplot(birds_m_shoreline, aes(x = distance_to_river_mouth, y = birds_per_100m, fill = river)) +
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

