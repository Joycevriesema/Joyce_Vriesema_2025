rm(list = ls())

# load libraries
library(tidyverse)  # includes dplyr, ggplot2
library(vegan)      # multivariate analysis of ecological community data   
library(psych)      # useful for panel plots of multivariate datasets
library(ggrepel)    # label placement for ggplot2 to avoid overlapping text
library(patchwork)  # combine multiple ggplot2 plots into one layout 
library(lme4)         # for mixed models
library(lmerTest)     # for tests of significance of mixed-effects models
library(multcomp)     # registers cld() method for emmGrid
library(multcompView) # generates the letters
library(emmeans)      # pairwise comparison
library(ellipse)      # drawing ellipses in plots
library(lubridate)

# load data_water
data_water <- read.csv("data_water.csv") |>
  dplyr::mutate(transect_date_time = paste(transect_ID, substr(date_time_end, 1, 6), substr(date_time_start, 12, 16), sep = " "))|>
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
         distance_to_river_mouth= case_when(
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
           transect_ID %in% c("pap_3", "pap_4", "pap_5", "pap_6","tree_5", "tree_6", "tree_7", "tree_8") ~ "Robana"))

# merge data with meta data
# load data transect
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv")|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))

data_water <- data_water |>
  left_join(data_transect, by = c("transect_ID", "date", "direction_fish", "run_ID", "transect_run_ID" ))


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

# add the meta data: habitat, river and distance with the rows used in the PCA
habitat_vector <- data_water$habitat_detailed2[rows_to_keep]
river_vector   <- data_water$river[rows_to_keep]
dist_vector     <- data_water$distance_to_river_mouth[rows_to_keep]

# extract PCA scores into a dataframe and add the meta data catgeories
scores <- as.data.frame(pca_result$x) |>
  mutate(
    river = factor(river_vector, levels = c("Mbalageti","Robana")),
    distance_to_river_mouth = factor(dist_vector, levels = c("Mouth","Mid","Far")),
    habitat = habitat_vector
  )

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

# draw ellipses per river x distance to river mouth
## Ellipses per river × distance_to_river_mouth (95%)
ellipse_df <- scores %>%
  group_by(river) %>%
  group_split() %>%
  lapply(function(df) {
    e <- as.data.frame(ellipse::ellipse(stats::cov(df[,c("PC1","PC2")], use="complete.obs"),
                                        centre = colMeans(df[,c("PC1","PC2")], na.rm = TRUE),
                                        level = 0.95, npoints = 200))
    colnames(e) <- c("x","y")
    e$river <- df$river[1]
    e$group <- e$river
    e
  }) %>% bind_rows()

# plot
plot1 <- ggplot(scores, aes(x = PC1, y = PC2,
                            color = river,
                            shape = distance_to_river_mouth)) +
  geom_point(size = 2, alpha = 0.9) +
  geom_segment(data = loadings,
               aes(x = 0, y = 0, xend = PC1, yend = PC2),
               inherit.aes = FALSE,
               arrow = arrow(length = unit(0.25, "cm")),
               color = "black", linewidth = 0.5, show.legend = FALSE) +
  geom_text_repel(data = loadings,
                  inherit.aes = FALSE,
                  aes(x = PC1, y = PC2, label = varname),
                  color = "black", size = 5.5,
                  nudge_x = loadings$nudge_x, nudge_y = loadings$nudge_y,
                  force_pull = 0.3, segment.color = "black",
                  segment.size = 0.2, box.padding = 0.4,
                  point.padding = 0.2, min.segment.length = 0) +
  labs(
    x = paste0("PC1 (", round(summary(pca_result)$importance[2, 1] * 100), "%)"),
    y = paste0("PC2 (", round(summary(pca_result)$importance[2, 2] * 100), "%)"),
    color = "River", shape = "Distance to river mouth"
  ) +
  coord_fixed() +
  theme_minimal() +
  scale_color_manual(values = c("Mbalageti" = "#0072B2",
                                "Robana"    = "#56B4E9")) +
  scale_shape_manual(values = c("Mouth" = 1,
                                "Mid"   = 2,
                                "Far"   = 0)) +
  guides(
    color = guide_legend(order = 1),
    shape = guide_legend(order = 2)
  )

plot1
geom_path(data = ellipse_df,
          aes(x = x, y = y,
              group = group),
          inherit.aes = FALSE,
          linewidth = 0.5,
          show.legend = TRUE) +

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

# filter out the extra observations because taken after 12 PM and not all transects have been sampled even,y after 12 pm 
pca_without_extra <- pca_with_meta |>
  filter(!coalesce(str_detect(run_ID, regex("extra", ignore_case = TRUE)), FALSE))

# make time a factor with 6 levels from 7 until 12 hour
pca_without_extra <- pca_without_extra|>
  mutate(
    time_parsed = if (inherits(time, "POSIXt")) time else
      lubridate::parse_date_time(as.character(time),
                                 orders = c("H:M:S","H:M","Y-m-d H:M:S","Y-m-d H:M")),
    hour_cont = lubridate::hour(time_parsed) +
      lubridate::minute(time_parsed)/60 +
      lubridate::second(time_parsed)/3600
  ) %>%
  filter(floor(hour_cont) %in% 7:12) %>%                             # keep 7–12 only
  mutate(hour = factor(floor(hour_cont), levels = 7:12, ordered = TRUE))


m2 <- lmer(PC1 ~ river * distance_to_river_mouth + (1|waves) + (1 | date) + (1 | transect_ID) + (1 | hour),data = pca_without_extra, REML = TRUE)

anova(m2)
# no significant effect of river, distance to river or interaction between river and distance to river mouth
summary(m2)
# variance date =2.36
# variance transect_ID = 0.12
# variance hour = 1.18
# variance waves = 0.78 
# variance of transect_ID and waves small so try different models

# model same as m2 but without transect_ID
m3 <- lmer(PC1 ~ river * distance_to_river_mouth + (1|waves) + (1 | date) + (1 | hour),data = pca_without_extra, REML = TRUE)
anova(m3) # nothing significant 
summary(m3)
# variance date = 2.27
# variance hour = 0.94
# variance waves = 0.19

# model same as m2 but without transect ID and waves as random effects
m4 <- lmer(PC1 ~ river * distance_to_river_mouth + (1 | date) + (1 | hour),data = pca_without_extra, REML = TRUE)

anova(m4)
# nothing signficant 
summary(m4)
# variance date = 1.90
# variance hour = 1.19

# compare model fit
AIC(m2,m3,m4) # model m2 has the lowest AIC value so seems the best model


#### model with all observations so after 12 pm, run_ID extra included ####
# make time a factor with 6 levels from 7 until 12 hour
pca_with_meta <- pca_with_meta|>
  mutate(
    time_parsed = if (inherits(time, "POSIXt")) time else
      lubridate::parse_date_time(as.character(time),
                                 orders = c("H:M:S","H:M","Y-m-d H:M:S","Y-m-d H:M")),
    hour_cont = lubridate::hour(time_parsed) +
      lubridate::minute(time_parsed)/60 +
      lubridate::second(time_parsed)/3600
  ) %>%
  filter(floor(hour_cont) %in% 7:12) %>%                             # keep 7–12 only
  mutate(hour = factor(floor(hour_cont), levels = 7:12, ordered = TRUE))

m1<- lmer(PC1 ~ river * distance_to_river_mouth + (1 | transect_ID) + (1|date)+ (1|hour) +(1|waves), data = pca_with_meta)
anova(m1)
# nothing significant
summary(m1)
# variance date = 2.48
# variance transect = 0.12
# variance hour = 0.88
# variance waves = 0.80

# model without transect_ID
m1a <- lmer(PC1 ~ river * distance_to_river_mouth + (1|waves) + (1 | date) + (1 | hour),data = pca_with_meta, REML = TRUE)
anova(m1a)
# river significant effect
# distance to river mouth
summary(m1a)
# variance date = 2.35
# variance hour = 0.70
# variance waves 0.86

# model without waves and transect_ID
m2a <- lmer(PC1 ~ river * distance_to_river_mouth + (1 | date) + (1 | hour),data = pca_with_meta, REML = TRUE)
anova(m2a)
# river and distance to river mouth significant 
summary(m2a)
# variance date 1.95
# variance hour 0.87

# compare models
AIC(m1,m1a,m2a) # lowest AIC value for m1

#### continue with m2 ####

# pairwise comparisons for m2
emm <- emmeans(m2, ~ river * distance_to_river_mouth, drop = TRUE)
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

# plot PC1 which represents water aeration + significance letters
plot1 <- ggplot(pca_without_extra, aes(x = distance_to_river_mouth, y = PC1, fill= river)) +
  geom_boxplot()+
  theme(text= element_text(size=14))+
  labs(x= "Distance to river mouth", y= "Water aeration (PC1)")+
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



#### analysis water turbidity ####
# work with the without_extra data
mpc2 <- lmer(PC2 ~ river * distance_to_river_mouth + (1|waves) + (1 | date) + (1 | transect_ID) + (1 | hour),data = pca_without_extra, REML = TRUE)
anova(mpc2)
# significant effect of river ***
# significant effect of interaction river and distance to river mouth
summary(mpc2)
# variance date = 0.91
# variance transect_ID = 0.21
# variance hour = 0.52
# variance waves = 0.15

# try model without random effect waves
mpc2a <- lmer(PC2 ~ river * distance_to_river_mouth + (1 | date) + (1 | transect_ID) + (1 | hour),data = pca_without_extra, REML = TRUE)
anova(mpc2a)
# only river significant effect
summary(mpc2a)
# variance date = 0.74
# variance transectID = 0.17
# variance hour = 0.46

# model without transect and waves as random effects
mpc2b <- lmer(PC2 ~ river * distance_to_river_mouth + (1 | date) + (1 | hour),data = pca_without_extra, REML = TRUE)
anova(mpc2b)
# river significant ***
# river x distance significant ***


# compare models
AIC(mpc2, mpc2a,mpc2b) # mpc2 best fit model

# pairwise comparisons for mpc2
emm2 <- emmeans(mpc2, ~ river * distance_to_river_mouth, drop = TRUE)
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

pca_with_meta$distance_to_river_mouth <- factor(pca_with_meta$distance_to_river_mouth, levels= c("Mouth", "Mid", "Far"))
# plot PC2 which represents water clarity (visible water quality)
plot2 <- ggplot(pca_with_meta, aes(x = distance_to_river_mouth, y = PC2, fill= river)) +
  geom_boxplot()+
  #facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  labs(x= "Distance to river mouth", y= "Water turbidity (PC2)")+
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



