rm(list = ls())

# load libraries
library(tidyverse)  # includes dplyr, ggplot2
library(lavaan)     # SEM modelling
library(lavaanPlot) # plot SEM model

#### data preparation ####
# combine datasets of fish, birds and water quality to have one big dataframe with 94 observations in total for every parameter
# load datasets: data_bird, fish_data and PCA_data
data_bird <- read.csv("data_bird_SEM.csv") |>
  pivot_wider(names_from = bird_species_ID,
              values_from = total_count,
              names_prefix = "",
              values_fill = 0) |>
  mutate(date = as.Date(date))
fish_data <- read.csv("fish_data.csv") |>
 mutate(date = as.Date(date))
PCA_data <- read.csv("PCA_data.csv") |>
  mutate(date = as.Date(date))


# fish_data and data_bird first need to be merged with corresponding transect_run_ID from data transect
# load data_transect and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv") |>
  dplyr::filter(date != "5-Feb-2025" & date != "8-Feb-2025") |>
  dplyr::select(transect_run_ID,run_ID, transect_ID, date) |>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))


# join dataframes
joined_data <- PCA_data |>
  left_join(
    fish_data |>
      dplyr::select(-one_of(setdiff(intersect(names(data_bird), names(fish_data)), "transect_run_ID"))),
    by = "transect_run_ID"
  ) %>%
  { 
    current_names <- names(.)
    left_join(
      .,
      data_bird |>
        dplyr::select(-one_of(setdiff(intersect(current_names, names(data_bird)), "transect_run_ID"))),
      by = "transect_run_ID"
    )
  }


# merge distance to river mouth dataframe with the joined dataframe
# river and habitat_main are still categorical variables, to use the lavaan package for the SEM, need to create dummy variables with values of 0 and 1
SEM_data <- joined_data |>
  mutate(
  habitat_dummy = ifelse(habitat_main == "Papyrus", 1, 0))    # 1 = papyrus, 0 = trees
  
  
  
#close_to_mouth= case_when (distance_to_river_mouth %in% "Mouth" ~ 1, 
                             distance_to_river_mouth %in% "Mid" ~ 0,
                              distance_to_river_mouth %in% "Far" ~ 0),
  close_to_human = case_when (distance_to_river_mouth %in% "Far" ~ 1,
                              distance_to_river_mouth %in% "Mid" ~ 0,
                              distance_to_river_mouth %in% "Mouth" ~ 0)
)


#### start SEM ####
# variables of interest
# river_dummy = river site with 1 = Mbalageti, 0 = Robana
# habitat_dummy = habitat_main with # 1 = papyrus, 0 = trees
# distance = distance to river mouth measured in meters
# chemical_water_quality_PC1 = PC1 values, variables included: HDO, HDO-saturation, temperature, pH and oxidation reduction potential
# visible_water_quality_PC2 = PC2 values, variables included: TDS, conductivity, turbidity 
# Small schooling fish = number of school counts, represents the number of dagaa school counts
# Large predatory fish = number of individual predatory fish
# Pied kingfisher = number of Pied kingfisher counts along the shoreline

# inspect linearity among the variables in a pairs panel plot
psych::pairs.panels(SEM_data %>% dplyr::select(habitat_dummy, chemical_water_quality_PC1, visible_water_quality_PC2, Small_schooling_fish, Large_predatory_fish, corrected_pied_count),
                    stars = T, ellipses = F)

# river has positive but weak correlation with chemical water quality
# river has negative correlation with visible water quality ***
# river positive correlation with schooling fish **
# river positive correlation with large predatory fish *
# river negative correlation with pied kingfisher ***
# distance to river mouth has negative correlation with chemical water quality * --> so higher distance from river mouth leads to lower values of the chemical water quality
# distance negative correlation with visible water quality
# visible water quality negative correlation with small schooling fish ** --> lower values of visible water quality leads to higher amount of school fish --> lower values of visible water quality indicate lower values for conductivity, turbidity and TDS so suggesting cleared water
# visible water quality has negative correlation with large predatory fish ** --> lower values of visible water quality indicate lower values for conductivity, turbidity and TDS so suggesting cleared water
# visible water has positive correlation with Pied kingfisher ***--> higher values for visible water quality, so higher value sof TDS, conductivity and turbidity --> less clear water would suggest higher number of pied kingfishers
# no correlation between chemical and visible water quality
# chemical water quality has positive correlation with small schooling fish *** --> higher values for pH, temp, HDO etc give higher number of fish schools
# chemical water quality has positive correlation with large predatory fish
# small schooling fish and large predatory fish positive correlation ***
# pied kingfisher has negative weaker correlation with small schooling fish *
# pied kingfisher has a negative correlation with river types ***
# pied kingfisher has a negative correlation with habitat *
# habitat has a negative correlation with small schooling fish **
# habitat has a negative correlation with large predatory fish ***

# standardize the continuous variables of interest that have different scale. This causes the multiple regression to yield standardized regression coefficients, independent of the scale of measurement of each. But for variables with values 0 and 1 don't standardize. The counts variables of fish and birds need to rescaled
SEM_data_std <- within(SEM_data, {
  chemical_water_quality_PC1 <- as.numeric(scale(chemical_water_quality_PC1))
  visible_water_quality_PC2  <- as.numeric(scale(visible_water_quality_PC2))
})|>
dplyr::mutate(
    dplyr::across(
      c(corrected_pied_count, Small_schooling_fish, Large_predatory_fish),
      ~ . / sd(., na.rm = TRUE) # dividing by their own SD
    )
  )

# run the multiple regression
multreg_std<-lm(corrected_pied_count ~ habitat_dummy + chemical_water_quality_PC1 + visible_water_quality_PC2 + Small_schooling_fish + Large_predatory_fish, data=SEM_data)
summary(multreg_std)
# the standardized coefficients give an estimate of the weight of each path in this diagram
# estimate values
# habitat -0.98
# chemical water quality 0.56
# visible water quality 5.35*
# small schooling fish -0.38 *
# large predatory fish 0.02
# These results suggest that visible water quality and habitat have the most important effect on number of pied kingfishers

# However the multiple regression tests how much variation is explained by a predictor that is not explained by the other predictors in the model. So when two predictors are strongly correlated, each may not explain variation in the response independently, why the regression of each with the response may be highly significant.
# path analysis (form of SEM) can analyse causal relations between the different predictors, instead of hypothesizing that they all directly affect the response variable

# make SEM model and fit with lavaan 
pied_model <- 'corrected_pied_count ~ habitat_dummy + visible_water_quality_PC2 + Small_schooling_fish
Small_schooling_fish ~ visible_water_quality_PC2 + chemical_water_quality_PC1'

pied_model.fit <- lavaan::sem(pied_model, estimator = "MLR", data= SEM_data_std)

summary(pied_model.fit, standardized=T, fit.measures=T,rsquare=T)
# CFI = 0.84, TLI= 0.45
# RMSEA= 0.19, SRMR= 0.07
# CFI/TLI are lower than conventional thresholds (>.90 or >.95), meaning fit isn’t great but modest
# SRMR = .074 is within acceptable range (<.08).
#Robust RMSEA = .186 is high (> .08).
# This suggests the model captures some structure but doesn’t fully explain the covariance in your data — which makes sense with only 84 cases and a small model.


# goodness and badness of fit measures of the model
# the goodness of fit measures CFI and TLI should be >0.9 for the model to be accepted, while the badness of fit measures RMSEA and SRMR should be be <0.1 for the model to be accepted
# if these parameters are not in these ranges, then your causal model is not supported by the data, and you test a different model. this can involve the inclusion of different predictors, or testing an alternative causal structure among the predictors.



# make the pathway diagram
lavaanPlot(pied_model.fit, coefs=T, stand=T, graph_options = list(rankdir= "LR"))

# estimating the indirect and total effect of PC2
# indirect effect of PC2 on pied count= -0.28 x -0.24 = 0.07
# total effect of PC2 on pied count is 0.28 + 0.07 = 0.35 so overal PC2 has a positive association but partially through negative effect on fish

# estimating the indirect and total effect of PC1



pied_model_indirect <- '
  # label paths
  corrected_pied_count ~ c1*visible_water_quality_PC2 + c2*Small_schooling_fish + habitat_dummy
  Small_schooling_fish ~ a1*visible_water_quality_PC2 + a2*chemical_water_quality_PC1

  # (in)direct + total effects
  indirect_PC2 := a1 * c2                  # PC2 -> fish -> pied
  direct_PC2   := c1                       # PC2 -> pied
  total_PC2    := c1 + (a1 * c2)           # direct + indirect

  indirect_PC1 := a2 * c2                  # PC1 -> fish -> pied
  total_PC1    := a2 * c2                  # (no direct PC1 -> pied path)
'

fit_mlr <- sem(pied_model_indirect, data = SEM_data_std, estimator = "MLR")
summary(fit_mlr, standardized = TRUE, fit.measures = TRUE)

modindices(fit_mlr, sort = TRUE, minimum.value = 5)



# new model without habitat
pied_modelx <- 'corrected_pied_count ~ visible_water_quality_PC2 + Small_schooling_fish
Small_schooling_fish ~ visible_water_quality_PC2 + chemical_water_quality_PC1'

pied_model.fit <- lavaan::sem(pied_modelx, estimator = "MLR", data= SEM_data_std)

summary(pied_model.fit, standardized=T, fit.measures=T,rsquare=T)

# make the pathway diagram
lavaanPlot(pied_model.fit, coefs=T, stand=T, graph_options = list(rankdir= "LR"))


#### SEM without tree counts in papyrus ####
# load datasets: data_bird, fish_data and PCA_data
data_bird2 <- read.csv("data_bird_SEM2.csv") |>
  pivot_wider(names_from = bird_species_ID,
              values_from = total_count,
              names_prefix = "",
              values_fill = 0) |>
  mutate(date = as.Date(date))
fish_data <- read.csv("fish_data.csv") |>
  mutate(date = as.Date(date))
PCA_data <- read.csv("PCA_data.csv") |>
  mutate(date = as.Date(date))


# fish_data and data_bird first need to be merged with corresponding transect_run_ID from data transect
# load data_transect and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv") |>
  dplyr::filter(date != "5-Feb-2025" & date != "8-Feb-2025") |>
  dplyr::select(transect_run_ID,run_ID, transect_ID, date) |>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))


# join dataframes
joined_data <- PCA_data |>
  left_join(
    fish_data |>
      dplyr::select(-one_of(setdiff(intersect(names(data_bird2), names(fish_data)), "transect_run_ID"))),
    by = "transect_run_ID"
  ) %>%
  { 
    current_names <- names(.)
    left_join(
      .,
      data_bird2 |>
        dplyr::select(-one_of(setdiff(intersect(current_names, names(data_bird2)), "transect_run_ID"))),
      by = "transect_run_ID"
    )
  }


# merge distance to river mouth dataframe with the joined dataframe
# river and habitat_main are still categorical variables, to use the lavaan package for the SEM, need to create dummy variables with values of 0 and 1
SEM_data <- joined_data |>
  mutate(
    habitat_dummy = ifelse(habitat_main == "Papyrus", 1, 0))    # 1 = papyrus, 0 = trees

SEM_data_std <- within(SEM_data, {
  chemical_water_quality_PC1 <- as.numeric(scale(chemical_water_quality_PC1))
  visible_water_quality_PC2  <- as.numeric(scale(visible_water_quality_PC2))
})|>
  dplyr::mutate(
    dplyr::across(
      c(corrected_pied_count, Small_schooling_fish, Large_predatory_fish),
      ~ . / sd(., na.rm = TRUE) # dividing by their own SD
    )
  )

# make SEM model and fit with lavaan 
pied_model <- 'corrected_pied_count ~ habitat_dummy + visible_water_quality_PC2 + Small_schooling_fish
Small_schooling_fish ~ visible_water_quality_PC2 + chemical_water_quality_PC1'

pied_model.fit <- lavaan::sem(pied_model, estimator = "MLR", data= SEM_data_std)

summary(pied_model.fit, standardized=T, fit.measures=T,rsquare=T)

lavaanPlot(pied_model.fit, coefs=T, stand=T, graph_options = list(rankdir= "LR"))






#### SEM with water turbidity and aeration being calculated inside the SEM ####
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
# from PCA script it is clear there are three outliers, which are rows 377, 622 and 740 t=from the water dataframe 
data_water <- data_water |>
  dplyr::slice(c(-377,-622, -740))

# merge with data transect and filter out the extr observations
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv")|>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))

data_water <- data_water |>
  left_join(
    data_transect,
    by = c("transect_ID", "date", "direction_fish", "run_ID", "transect_run_ID")
  ) |>
  filter(!str_starts(run_ID, "extra_"))

# make average values  
data_water_avg <- data_water |>
  group_by(transect_run_ID) |>
  summarise(across(c(TDS, turb, cond, pH, temp, ORP, HDO, HDO_sat), ~ mean(.x, na.rm = TRUE)))

# load datasets: data_bird, fish_data and PCA_data
data_bird <- read.csv("data_bird_SEM.csv") |>
  pivot_wider(names_from = bird_species_ID,
              values_from = total_count,
              names_prefix = "",
              values_fill = 0) |>
  mutate(date = as.Date(date))
fish_data <- read.csv("fish_data.csv") |>
  mutate(date = as.Date(date))

# merge with the fish and bird datasets
merged_data <- data_bird |>
  left_join(data_water_avg, by = "transect_run_ID") 

merged_data <- merged_data %>%
  left_join(
    fish_data %>% dplyr::select(-any_of(setdiff(names(merged_data), "transect_run_ID"))),
    by = "transect_run_ID"
  )|>
  mutate(habitat_dummy = ifelse(habitat_main == "Papyrus", 1, 0)) 


# standardize the data
merged_data_std <- merged_data %>%
  mutate(across(where(is.numeric) & !matches("habitat_dummy"), 
                ~ scale(.x)[, 1]))
# run the multiple regression
multreg_std<-lm(corrected_pied_count ~ habitat_dummy + turb + cond + TDS + temp + pH + ORP + HDO + HDO_sat + Small_schooling_fish + Large_predatory_fish, data=merged_data_std)
summary(multreg_std)

# leave out conductivity or TDS because 1:1 correlation give error in latent variables
# leave out HDO or HDO_sat because 1:1 correlation
new_pied_model <- '
  # latent variables
  water_turbidity =~ TDS + turb 
  water_aeration  =~ HDO + pH + temp

  # structural model
  corrected_pied_count ~ habitat_dummy + water_turbidity + Small_schooling_fish
  Small_schooling_fish ~ water_turbidity + water_aeration
'

model.fit <- lavaan::sem(new_pied_model, estimator = "MLR", data = merged_data_std)

summary(model.fit, fit.measures = TRUE, standardized = TRUE)
# CFI = 0.88, TLI= 0.81
# RMSEA= 0.14, SRMR= 0.09

# make the pathway diagram
lavaanPlot(model.fit, coefs=T, stand=T, graph_options = list(rankdir= "LR"))
# the low p value of chi-square indicates that this model has poor fit 
# number of parameters seems high and number of observations low
# so maybe better to use the pca scores, that reduces the numbers of parameters but it becomes a path model with only observed values instead of using latent variables, main focus is structural relationships (paths), not validating latent constructs
