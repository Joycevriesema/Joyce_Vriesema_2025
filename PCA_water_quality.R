rm(list = ls())

# load libraries
library(dplyr)
library(vegan) # multivariate analysis of ecological community data 
library(psych) # for panel plots of multivariate datasets
library(tidyverse)
library(ggrepel)

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


# explore the correlations among the environmental factors in a panel pairs plot
psych::pairs.panels(pca_input,smooth=F,ci=T,ellipses=F,stars=T,method="pearson")
# very strong relation between conductivity and TDS, one on one correlation
# HDO and HDO_sat one on one correlation
# strong negative correlation between ORP and HDO & HDO_sat
# strong positive correlation between pH and HDO
# strong negative correlation between pH and ORP

# PCA with all data
pca_input <- data_water |>
  dplyr::select(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS)|>
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))

pca_result <- prcomp(pca_input,center=T, scale. = TRUE)
pca_result
summary(pca_result)

biplot(pca_result,xlab="PC1 48%",ylab="PC2 25%")

# points 376 and 697 suggest outliers
# point 376 has a very high turbidity, very low conductivity and low TDS, depth 0.46
# point 376 was recorded at 11:07:58, which coincides with the exact start time of the transect, seconds moved = 2 and meters moved = 2  suggesting the multiprobe may not yet have been fully submerged. The shallow depth (0.46m) supports this assumption. Partial submersion may also explain the unusually high HDO saturation and turbidity, as the sensors may have measured surface water or even been partially exposed to air.

# point 697: TDS relatively low, conductivity low, turbidity high, depth 0.47, suggest that the probe is held more close to the surface, point 696 is 0.81 meter and point 697 is 0.47 and 698 also 0.47 so maybe the probe was drawn a bit to the surface
# remove these two outliers from the data set, as these might nog be representative measurements
pca_input2 <- data_water |>
  dplyr::slice(c(-377,-698)) |>
  dplyr::select(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS)|>
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))

pca_result <- prcomp(pca_input2,center=T, scale. = TRUE)
pca_result
summary(pca_result)

biplot(pca_result, xlab = "PC1 (48%)", ylab = "PC2 (28%)")
# now the plot looks better, still there seems to be a cluster of so points to the right

# try PCA with averaged data
pca_input_avg <- data_water |>
  dplyr::slice(c(-377, -698)) |>
  dplyr::filter(if_all(c(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS),
                       ~ !is.na(.) & !is.infinite(.))) |>
  dplyr::group_by(transect_ID) |>
  dplyr::summarise(
    temp = mean(temp),
    pH = mean(pH),
    ORP = mean(ORP),
    turb = mean(turb),
    cond = mean(cond),
    HDO = mean(HDO),
    HDO_sat = mean(HDO_sat),
    TDS = mean(TDS),
    .groups = "drop"
  )

pca_result_avg <- prcomp(
  pca_input_avg |> dplyr::select(-transect_ID),
  center = TRUE,
  scale. = TRUE
)
pca_result_avg
summary(pca_result_avg)

biplot(pca_result_avg,xlab="PC1 53%",ylab="PC2 36%")
# so the raw data explains 48+25= 73% of all variation
# averaged data explains 53+36= 89% of all variation

# comparing the two biplots. the arrows of temp, turb, cond and TDS seem to change direction, so try o find out why this is

# Compute correlation matrices
cor_raw <- cor(pca_input2, use = "pairwise.complete.obs")
cor_avg <- cor(pca_input_avg %>% dplyr::select(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS), use = "pairwise.complete.obs")

# Convert matrices to long format data frames
cor_raw_df <- as.data.frame(as.table(cor_raw)) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation_raw = Freq)

cor_avg_df <- as.data.frame(as.table(cor_avg)) %>%
  rename(Var1 = Var1, Var2 = Var2, Correlation_avg = Freq)

# Join on Var1 and Var2
cor_compare <- left_join(cor_raw_df, cor_avg_df, by = c("Var1", "Var2")) %>%
  arrange(Var1, Var2)

# View table
print(cor_compare_unique)
# temp vs turb changes from 0.14 (weak positive) in raw to 0.51 (moderate positive) in averaged.
# ORP vs turb flips from 0.21 (positive) in raw to -0.15 (negative) in averaged.
# cond vs turb increases from 0.41 (moderate) in raw to 0.85 (very strong positive) in averaged.
# cond vs temp flips from nearly zero (-0.007) in raw to 0.56 (moderate positive) in averaged.
# pH vs turb flips from -0.16 (weak negative) in raw to 0.24 (weak positive) in averaged.
# TDS vs temp flips from -0.007 in raw to 0.56 (moderate positive) in averaged.
# TDS vs turb rises from 0.41 (moderate) in raw to 0.85 (very strong positive) in averaged
# cond vs pH flips from -0.19 (negative) in raw to +0.04 (near zero positive) in averaged.
# HDO vs turb shifts from -0.07 (near zero) in raw to 0.24 (weak positive) in averaged.
# HDO_sat vs turb shifts from -0.05 in raw to 0.27 in averaged.

# the averaging changed correlation and might hide variability 
# I want to use the clustering of the parameters for factors of water aeration and water clarity so it might be better to have enough variability

# make plot
# Step 1: Get row indices
rows_to_keep <- as.numeric(rownames(pca_input2))

# Step 2: Run PCA
pca_result <- prcomp(pca_input2, center = TRUE, scale. = TRUE)

# Step 3: Get habitat vector aligned to PCA rows
habitat_vector <- data_water$habitat_detailed_2[rows_to_keep]

# Step 4: Prepare scores dataframe
scores <- as.data.frame(pca_result$x)
scores$habitat <- habitat_vector

# Step 5: Extract loadings for arrows
loadings <- as.data.frame(pca_result$rotation[, 1:2])  # PC1 and PC2
loadings$varname <- rownames(loadings)
loadings$varname <- dplyr::recode(
  loadings$varname,
  "cond" = "conductivity",
  "turb" = "turbidity",
  "temp" = "temperature",
  "HDO_sat" = "HDO_saturation"
)
loadings$nudge_x <- 0
loadings$nudge_y <- 0

loadings$nudge_x[loadings$varname == "TDS"] <- 0.5
loadings$nudge_x[loadings$varname == "conductivity"] <- -0.3
loadings$nudge_x[loadings$varname == "turbidity"] <- -0.3
loadings$nudge_x[loadings$varname == "temperature"] <- 0.1
loadings$nudge_x[loadings$varname == "HDO"] <- 0.4
loadings$nudge_x[loadings$varname == "HDO_saturation"] <- 0.4
loadings$nudge_y[loadings$varname == "HDO_saturation"] <- -0.1
loadings$nudge_y[loadings$varname == "ORP"] <- -0.1
# Scale arrows for visibility (optional)
arrow_scale <- 4.5 # adjust this if arrows too small/big
loadings$PC1 <- loadings$PC1 * arrow_scale
loadings$PC2 <- loadings$PC2 * arrow_scale

# Step 6: Plot
library(ggplot2)

ggplot(scores, aes(x = PC1, y = PC2, color = habitat)) +
  geom_point(size = 1.5, alpha= 0.9) +
  # Add arrows
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.5) +
  # Add variable labels
  geom_text_repel(
    data = loadings,
    aes(x = PC1, y = PC2, label = varname),
    color = "black",
    size = 3.5,
    nudge_x = loadings$nudge_x,
    nudge_y = loadings$nudge_y,
    force_pull= 0.3,
    segment.color = "black",     # color of the connector line
    segment.size = 0.2,           # thickness of connector line
    box.padding = 0.4,            # space around text box to avoid overlap
    point.padding = 0.2,          # space around arrow tip point
    min.segment.length = 0        # show segment even if very short
  )+
  labs(
    title = "",
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100), "%)"),
    color = "Habitat type"
  ) +
  coord_fixed() +
  theme_minimal() +
  scale_color_manual(values = c("#332288", "#44AA99", "#DDCC77", "#AA4499", "#88CCEE"))



# extract PCA values and add as column to data water
# extract PCA scores
scores <- as.data.frame(pca_result$x)

# Step 5: Add metadata (e.g., habitat)
meta <- data_water[rows_to_keep, ]  # Match the same rows
pca_with_meta <- cbind(meta, scores)

# plot PC1 which represents water aeration
ggplot(pca_with_meta, aes(x = distance_to_river_mouth, y = PC1, fill= river)) +
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  labs(x= "Distance to river mouth", y= "Water aeration (PC1)")+
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

# plot PC2 which represents water clarity
ggplot(pca_with_meta, aes(x = distance_to_river_mouth, y = PC2, fill= river)) +
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  labs(x= "Distance to river mouth", y= "Water clarity (PC2)")+
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

# anova PC1
model_pc1 <- lm(PC1 ~ habitat_main + river + distance_to_river_mouth, data = pca_with_meta)
anova(model_pc1)
summary(model_pc1)

# anova PC2
model_pc2 <- lm(PC2 ~ habitat_main + river + distance_to_river_mouth, data = pca_with_meta)
anova(model_pc2)
summary(model_pc2)

# check assumptions