rm(list = ls())

# load libraries
library(tidyverse)    # includes dplyr, ggplot2, tidyr
library(lme4)         # for mixed models
library(lmerTest)     # for tests of significance of mixed-effects models
library(MASS)         # for negative binomial models
library(emmeans)      # pairwise comparison
library(glmmTMB)      # for negative binomial model including random effects
library(multcomp)     # registers cld() method for emmGrid
library(multcompView) # generates the letters
library(car)          # for Anova test
library(lubridate)


# load fish data
fish_data <-read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=946923967&single=true&output=csv")|>
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
           transect_ID %in% c("pap_3", "pap_4", "pap_5", "pap_6","tree_5", "tree_6", "tree_7", "tree_8") ~ "Robana"),
         fish_type= case_when(
           observation %in% c("H","V","O","B")~"Large_predatory_fish",
           observation %in% c("S1","S2","S3")~"Small_schooling_fish"
         ))|>
  group_by(transect_ID,river, habitat_detailed, habitat_detailed2, distance_to_river_mouth, habitat_main, date, fish_type ) |>
  summarise(total_count= n()) |>
  pivot_wider(names_from = fish_type,
  values_from = total_count,
  names_prefix = "",
  values_fill = 0) 

# merge data with meta data
# load data transect
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv")|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))

fish_data <- fish_data |>
  left_join(data_transect, by = c("transect_ID", "date"))

# filter out the extra observations (run_ID contains "extra"). These observations were taken after 12 PM. Not every transects has been sampled evenly after 12 PM so exclude these observations from the dataset --> 10 observations will be excluded 
fish_data <- fish_data |>
  filter(!coalesce(str_detect(run_ID, regex("extra", ignore_case = TRUE)), FALSE))

# make start_time_fish time a factor with 6 levels from 7 until 12 hour
# make time a factor with 6 levels from 7 until 12 hour
fish_data <- fish_data |>
  mutate(
    # extract a clean time like "9:29", "10:03", or "10:03:15"
    start_time_token = str_extract(as.character(start_time_fish),
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
    start_hour_fish = factor(hour(start_time_hms), levels = 7:12, ordered = TRUE)
  )

# save as csv for later use in SEM
write.csv(fish_data, "fish_data.csv", row.names = FALSE)

# exploratory graph small schooling fish
ggplot(fish_data, aes( x= distance_to_river_mouth, y= Small_schooling_fish, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  #scale_y_log10()+
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

# exploratory graph large predatory fish
ggplot(fish_data, aes( x= distance_to_river_mouth, y= Large_predatory_fish, fill= river))+
  geom_boxplot()+
  facet_grid(~habitat_main)+
  theme(text= element_text(size=14))+
  scale_y_log10()+
  labs(y="School count / 500 m ", title="Large predatory fish")


#### model for small schools (papyrus and trees included) ####
fish_data$river <- as.factor(fish_data$river)
fish_data$transect_ID <- as.factor(fish_data$transect_ID) 
fish_data$habitat_main <- as.factor(fish_data$habitat_main)

# count data with fixed effects of river site and distance to river mouth and random effect of transect ID included so taking a poisson model and include interaction terms
glmm_schools <- glmer(Small_schooling_fish ~ river * distance_to_river_mouth + (1 | transect_ID) + (1|date)+ (1|start_hour_fish), data = fish_data,family = poisson)

summary(glmm_schools)
# significantly lower school count at Robana river 
# significantly higher school count at mid distance compared to mouth
# significantly higher school count at far distance compared to mouth
# little variance between transects

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
# dispersion ratio = 420 is quite high and p<0.001 so this poisson model is highly overdispersed meaning that the variance in my real data is much higher than this poisson model can fit

# try negative binomial model 
nb_schools <- glmmTMB(Small_schooling_fish ~ river * distance_to_river_mouth + (1 | transect_ID)+ (1|date)+ (1|start_hour_fish),
  data = fish_data,
  family = nbinom2
)
Anova(nb_schools)
# river significant ***
# distance to river significant ***
# interaction river x distance to river mouth significant *
summary(nb_schools)
# variance of all random effects are really small
# variance of transect ID and start_hour are very small

# negative binomial model without random effect of transect_ID
nb_schools2 <- glmmTMB(Small_schooling_fish ~ river * distance_to_river_mouth + (1 | start_hour_fish)+ (1|date),
                      data = fish_data,
                      family = nbinom2
)

Anova(nb_schools2)
# river significant ***
# distance to river significant ***
# interaction river x distance to river mouth significant *
summary (nb_schools2)
# variance start_hour = 0.08
# variance data = 1.01

# model with only date as random factor, transect_ID left out
nb_schools3 <- glmmTMB(Small_schooling_fish ~ river * distance_to_river_mouth + (1|date),
                       data = fish_data,
                       family = nbinom2
)

Anova(nb_schools3)
# river significant ***
# distance to river significant ***
# interaction river x distance to river mouth significant *
# compare models

summary(nb_schools3)

# compared models to decide which model is best
AIC (nb_schools, nb_schools2, nb_schools3)
anova(nb_schools, nb_schools2, nb_schools3)
# prefer the simpler model so nb_schools3

# pairwise comparison
emm <- emmeans(nb_schools2, ~ river * distance_to_river_mouth, drop = TRUE)
pairs(emm)

# create letters
cld_dist <- cld(emm, Letters = letters)
cld_dist
letters_df <- as.data.frame(cld_dist)
letters_df <- letters_df |>
  mutate(y_pos = 100) |>
  mutate(group_label = gsub(" ", "", .group)) |>
  ungroup() |>
  drop_na()


# make plot with significance letters
plot_fish <- ggplot(fish_data, aes( x= distance_to_river_mouth, y= Small_schooling_fish, fill= river))+
  geom_boxplot(position = position_dodge(width = 0.75))+
  #facet_grid(~habitat_main)+
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
    aes(x = distance_to_river_mouth, y = y_pos, label = group_label, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    vjust=1.1
  )+
  labs(x= "Distance to river mouth", y="School count / 500 m ")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))
plot_fish


emm2 <- emmeans(nb_schools3, ~ river * distance_to_river_mouth, drop = TRUE)
pairs(emm2)

# create letters
cld_dist2 <- cld(emm2, Letters = letters)
cld_dist2
letters_df2 <- as.data.frame(cld_dist2)
letters_df2 <- letters_df2 |>
  mutate(y_pos = 100) |>
  mutate(group_label = gsub(" ", "", .group)) |>
  ungroup() |>
  drop_na()


# make plot with significance letters
plot_fish2 <- ggplot(fish_data, aes( x= distance_to_river_mouth, y= Small_schooling_fish, fill= river))+
  geom_boxplot(position = position_dodge(width = 0.75))+
  #facet_grid(~habitat_main)+
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
    data = letters_df2,
    aes(x = distance_to_river_mouth, y = y_pos, label = group_label, group= river),
    color = "black",
    size = 4,
    fontface = "bold",
    position = position_dodge(width = 0.75),
    inherit.aes = FALSE,
    vjust=1.1
  )+
  labs(x= "Distance to river mouth", y="Dagaa school count / 500 m ")+
  scale_fill_manual(values = c(
    "Mbalageti" = "#0072B2",
    "Robana" = "#56B4E9"))
plot_fish2










# save the plot
ggsave("plot small fish schools.png", plot_fish, width = 12, height = 6, dpi = 300)


# model with random effects included of transect_ID,date,start_time, waves and water_green

nb_schools <- glmmTMB(Small_schooling_fish ~ river * distance_to_river_mouth + (1 | transect_ID) + (1|date) + (1|waves) + (1|water_green)+ (1|start_time_fish),
                      data = fish_data,
                      family = nbinom2
)
summary(nb_schools)

