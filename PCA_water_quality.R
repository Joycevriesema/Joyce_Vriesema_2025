rm(list = ls())

# load libraries
library(tidyverse)  # includes dplyr, ggplot2
library(vegan)      # multivariate analysis of ecological community data   
library(psych)      # useful for panel plots of multivariate datasets
library(ggrepel)    # label placement for ggplot2 to avoid overlapping text
library(patchwork)  # combine multiple ggplot2 plots into one layout 

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


# explore the correlations among the factors in a panel pairs plot
psych::pairs.panels(data_water |> dplyr::select(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS),smooth=F,ci=T,ellipses=F,stars=T,method="pearson")
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

biplot(pca_result,xlab="PC1 49%",ylab="PC2 24%")

# points 376, 621 and 739 seem to be outliers, in the original dataframe of data_water these are rows 377, 622 and 740
# point 376 has a very high turbidity, very low conductivity and low TDS, depth 0.46
# point 376 was recorded at 11:07:58, which coincides with the exact start time of the transect, seconds moved = 2 and meters moved = 2  suggesting the multiprobe may not yet have been fully submerged. The shallow depth (0.46m) supports this assumption. Partial submersion may also explain the unusually high HDO saturation and turbidity, as the sensors may have measured surface water or even been partially exposed to air.
# there is some time issue with the Manta multimeter probe therefore it could be that the measurement already started but, while the transect was not yet started

# point 621: TDS relatively low, conductivity low, turbidity high, depth 0.47, suggest that the probe is held more close to the surface
# point 739 also has low conductivity and TDS.
# most values of conductivity in the df are above 100 and TDS above 70, these three data points deviate strongly
# remove these three outliers from the data set, as these might not be representative measurements

pca_input2 <- data_water |>
  dplyr::slice(c(-377,-622,-740)) |>
  dplyr::select(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS)|>
  filter(if_all(everything(), ~ !is.na(.) & !is.infinite(.)))

pca_result <- prcomp(pca_input2,center=T, scale. = TRUE)
pca_result
summary(pca_result)

biplot(pca_result, xlab = "PC1 (51%)", ylab = "PC2 (26%)")
# now the plot looks better, still there seems to be a cluster of points to the right

# try PCA with averaged data
pca_input_avg <- data_water |>
  dplyr::slice(c(-377,-622, -740)) |>
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

biplot(pca_result_avg,xlab="PC1 53%",ylab="PC2 37%")
# so the raw data explains 51+26= 76% of all variation
# averaged data explains 53+37= 90% of all variation

# comparing the two biplots, the arrows of temp, turb, cond and TDS seem to change direction
# try o find out why this is

# compute correlation matrices
cor_raw <- cor(pca_input2, use = "pairwise.complete.obs")
cor_avg <- cor(pca_input_avg |> dplyr::select(temp, pH, ORP, turb, cond, HDO, HDO_sat, TDS), use = "pairwise.complete.obs")

# convert matrices to long format data frames
cor_raw_df <- as.data.frame(as.table(cor_raw)) |>
  rename(Var1 = Var1, Var2 = Var2, Correlation_raw = Freq)

cor_avg_df <- as.data.frame(as.table(cor_avg)) |>
  rename(Var1 = Var1, Var2 = Var2, Correlation_avg = Freq)

# join on Var1 and Var2
cor_compare <- left_join(cor_raw_df, cor_avg_df, by = c("Var1", "Var2")) |>
  arrange(Var1, Var2)

# view table
print(cor_compare)
# temperature vs turbidity changes from 0.07 (weak positive) in raw to 0.50 (moderate positive) in averaged.
# ORP vs turbidity flips from 0.23 (positive) in raw to -0.08 (negative) in averaged.
# conductivity vs turbidity increases from 0.44 (moderate) in raw to 0.87 (very strong positive) in averaged.
# conductivity vs temp flips from (-0.11) in raw to 0.47 (moderate positive) in averaged.
# pH vs turbidity flips from -0.20 (weak negative) in raw to 0.15 (weak positive) in averaged.
# TDS vs temperature flips from -0.11 in raw to 0.47 (moderate positive) in averaged.
# TDS vs turbidity rises from 0.44 (moderate) in raw to 0.86 (very strong positive) in averaged
# HDO vs turbidity shifts from -0.11 (near zero) in raw to 0.2 (weak positive) in averaged.
# HDO_sat vs turbidity shifts from -0.09 in raw to 0.23 in averaged.

# the averaging pca changed correlation, or increased strength of correlations and might hide variability 
# I want to use the clustering of the parameters to give indication of chemical water quality and visible water quality, so it might be better to have enough variability in the data and choose for the raw data instead of averaging data

# make plot
# get row indices to align the meta data (categories) with the rows of the pca input
rows_to_keep <- as.numeric(rownames(pca_input2))

# run PCA
pca_result <- prcomp(pca_input2, center = TRUE, scale. = TRUE)

# add the meta data (habitat categories) with the rows used in the PCA
habitat_vector <- data_water$habitat_detailed_2[rows_to_keep]

# extract PCA scores into a dataframe and add the habitat categories
scores <- as.data.frame(pca_result$x)
scores$habitat <- habitat_vector

# extract loadings (arrows) for the first pricinipal components PC1 and PC2 
loadings <- as.data.frame(pca_result$rotation[, 1:2])  
# label the variable names
loadings$varname <- rownames(loadings)
loadings$varname <- dplyr::recode(
  loadings$varname,
  "cond" = "conductivity",
  "turb" = "turbidity",
  "temp" = "temperature",
  "HDO_sat" = "HDO_saturation"
)

# adjust the x and y positions of the labels
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

# scale arrows for visibility
arrow_scale <- 4.5 
loadings$PC1 <- loadings$PC1 * arrow_scale
loadings$PC2 <- loadings$PC2 * arrow_scale

# plot
plot1 <- ggplot(scores, aes(x = PC1, y = PC2, color = habitat)) +
  geom_point(size = 3.5, alpha= 0.9) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.5) +
  geom_text_repel(
    data = loadings,
    aes(x = PC1, y = PC2, label = varname),
    color = "black",
    size = 5.5,
    nudge_x = loadings$nudge_x,
    nudge_y = loadings$nudge_y,
    force_pull= 0.3,
    segment.color = "black",     
    segment.size = 0.2,          
    box.padding = 0.4,          
    point.padding = 0.2,      
    min.segment.length = 0
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
plot1

# from the PCA plot it seems that HDO saturation, HDO, pH and temperature are clustered over PC1 and ORP in opposite direction
# these factors give an indication about water aeration -->chemical water quality
# conductivity, TDS and turbidity cluster together over PC2, these factors give an indication about water clarity --> visible water quality
# seems some difference in aggregations of points between the two river types

# save the plot
ggsave("PCA.png", plot1 , width = 12, height = 6, dpi = 300)

# with these PCA results it is possible to plot the water parameters according to the PCA axis in one plot instead of every parameters separately
# therefore need to extract the PCA values and make a data frame of it

# extract PCA values and add as column to data water
scores <- as.data.frame(pca_result$x)

# add metadata of data water by matching rows
meta <- data_water[rows_to_keep, ]  
pca_with_meta <- cbind(meta, scores)


# test the effects of habitat, river and distance to river mouth on the PCA axis
# anova PC1 (chemical water quality)
model_pc1 <- lmer(PC1 ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID), data = pca_with_meta)
anova(model_pc1)
# habitat type significant **
# river significant effect *
# distance to river mouth highly significant ***
summary(model_pc1)
# papyrus significantly higher
# Robana significantly higher 
# far distance is significantly higher than mouth
# but zero variance of transect_ID, so remove from the model

# test without the random effect
model_pc1a <- lm(PC1 ~ habitat_main * river * distance_to_river_mouth, data = pca_with_meta)
anova(model_pc1a)
# habitat not significant
# river marginally significant 0.06
# distance to river mouth highly significant ***

# compare models
anova(model_pc1, model_pc1a)
# AIC lower for model pc1a and the p=1 so take out random effect

# pairwise comparisons for model_pc1a
emm <- emmeans(model_pc1a, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
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

# plot PC1 which represents water aeration (chemical water quality) + significance letters
plot1 <- ggplot(pca_with_meta, aes(x = distance_to_river_mouth, y = PC1, fill= river)) +
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  labs(x= "Distance to river mouth", y= "Chemical water quality (PC1)")+
  coord_cartesian(ylim = c(-5, 10))+
  scale_y_continuous(breaks = seq(-5, 10, by = 2.5))+
  geom_text(
    data = letters_df,
    aes(x = distance_to_river_mouth, y = Inf, label = group_label, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    vjust=1.8,
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
plot1

# anova PC2
model_pc2 <- lmer(PC2 ~ habitat_main * river * distance_to_river_mouth + (1 | transect_ID) , data = pca_with_meta)
anova(model_pc2)
# habitat not significant
# river significant effect **
# significant effect of distance to river mouth
# almost marginally signficant of interaction habitat:river and habitat:distance to river
summary(model_pc2)
# some variance between transects
# Robana significantly higher values 
# mid distance signficantly lower than mouth --> lower value for water clarity actually means a better water quality, cause lower values of turbditity, TDS and cond
# far distance significantly higher than mouth

# try model without random effect
model_pc2a <- lm(PC2 ~ habitat_main * river * distance_to_river_mouth, data = pca_with_meta)
anova(model_pc2a)
# every parameter is highly significant
summary(model_pc2a)

# compare models
anova (model_pc2, model_pc2a)
# AIC model pc2 is lower and the pvalue < 0.001 indicating to keep the random effect in the model

# pairwise comparisons for model_pc1a
emm2 <- emmeans(model_pc2, ~ river * distance_to_river_mouth * habitat_main, drop = TRUE)
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

# plot PC2 which represents water clarity (visible water quality)
plot2 <- ggplot(pca_with_meta, aes(x = distance_to_river_mouth, y = PC2, fill= river)) +
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  labs(x= "Distance to river mouth", y= "Visible water quality (PC2)")+
  coord_cartesian(ylim = c(-5, 10))+
  scale_y_continuous(breaks = seq(-5, 10, by = 2.5))+
  geom_text(
    data = letters_df2,
    aes(x = distance_to_river_mouth, y = Inf, label = group_label, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    vjust=1.8,
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

# make combined plot
plot2 <- plot2 + theme(strip.text.x = element_blank())

combined_plot <-
  (plot1 +
     theme(
       axis.title.x = element_blank(),
       axis.text.x  = element_blank(),
       axis.ticks.x = element_blank()
     )) /
  (plot2 + theme(strip.text.x = element_blank())) +
  plot_layout(guides = "collect", heights = c(1, 1)) &
  theme(
    legend.position = "right",
    axis.title.y = element_text(size = 10, face = "bold"),  # y-label size
    axis.text.y  = element_text(size = 12)                  # y tick size
  )

print(combined_plot)
ggsave("combined_plot.png", combined_plot, width = 12, height = 6, dpi = 300)


#### data preparation for Structural Equation Modelling SEM ####
# take average values for PC1 and PC2 for every transect + date combination so the number of observations is reduced to 94 and matches the number of observations in the different datasets
PCA_data <- pca_with_meta |>
  group_by(transect_run_ID,transect_ID,river, habitat_detailed, habitat_detailed_2, distance_to_river_mouth, habitat_main, date) |>
  dplyr::summarise(
    chemical_water_quality_PC1 = mean(PC1),
    visible_water_quality_PC2 = mean(PC2),
    .groups = "drop"
  )

# save as csv to use in SEM
write.csv(PCA_data, "PCA_data.csv", row.names = FALSE)



