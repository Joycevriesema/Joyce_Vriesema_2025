rm(list = ls())

# load libraries
library(tidyverse) # includes dplyr, ggplot2
library(lavaan)    # SEM modelling

#### data preparation ####
# combine datasets of fish, birds and water quality to have one big dataframe with 94 observations in total for every parameter
# load datasets: data_bird, fish_data and PCA_data
data_bird <- read.csv("data_bird.csv") |>
  pivot_wider(names_from = bird_species_ID,
              values_from = total_count,
              names_prefix = "",
              values_fill = 0) |>
  mutate(date = as.Date(date))
fish_data <- read.csv("fish_data.csv") |>
  pivot_wider(names_from = fish_type,
              values_from = total_count,
              names_prefix = "",
              values_fill = 0) |>
  mutate(date = as.Date(date))
PCA_data <- read.csv("PCA_data.csv") |>
  mutate(date = as.Date(date))

# fish_data and data_bird first need to be merged with corresponding transect_run_ID from data transect
# load data_transect and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv") %>%
  dplyr::filter(date != "5-Feb-2025" & date != "8-Feb-2025") |>
  dplyr::select(transect_run_ID,run_ID, transect_ID, date) |>
  mutate(date=as.Date(date, format= "%d-%b-%Y"))


# merge data_bird and data_transect
data_bird<- data_bird %>%
  dplyr::left_join(data_transect, by = c("transect_ID", "date"))

# merge fish_data and data_transect
fish_data <- fish_data %>%
  dplyr::left_join(data_transect, by = c("transect_ID", "date"))

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

# load the dataframe with distance to river mouth
distance_to_river <- read.csv("Distance_to_river_mouth.csv")|>
  rename(river= HubName,
         transect_ID= Transect,
         distance= HubDist) |>
  mutate(
    transect_ID = str_to_lower(transect_ID),
    transect_ID = str_replace(transect_ID, "([a-z]+)([0-9]+)$", "\\1_\\2")
  )

# merge distance to river mouth dataframe with the joined dataframe
# river and habitat_main are still categorical variables, to use the lavaan package for the SEM, need to create dummy variables with values of 0 and 1
SEM_data <- joined_data |>
  left_join(
    distance_to_river |>
      dplyr::select(transect_ID, river, distance),
    by = c("transect_ID", "river")
  )|>
  mutate(
  river_dummy   = ifelse(river == "Mbalageti", 1, 0),   # 1 = Mbalageti, 0 = Robana
  habitat_dummy = ifelse(habitat_main == "Papyrus", 1, 0)    # 1 = papyrus, 0 = trees
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
psych::pairs.panels(SEM_data %>% dplyr::select(river_dummy, habitat_dummy, distance, chemical_water_quality_PC1, visible_water_quality_PC2, `Small schooling fish`, `Large predatory fish`, `Pied kingfisher`),
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

# standardize the variables of interest. This causes the multiple regression to yield standardized regression coefficients, independent of the scale of measurement of each.
SEM_data_std <- SEM_data |>
  dplyr::select(distance, habitat_dummy, river_dummy, chemical_water_quality_PC1, visible_water_quality_PC2,
         `Small schooling fish`, `Large predatory fish`, `Pied kingfisher`) |>
  scale() |>
  as_tibble() |>
  rename(
    pied_kingfisher       = `Pied kingfisher`,
    small_schooling_fish  = `Small schooling fish`,
    large_predatory_fish  = `Large predatory fish`
  )

# run the multiple regression
multreg_std<-lm(pied_kingfisher ~ distance + habitat_dummy + river_dummy + chemical_water_quality_PC1 + visible_water_quality_PC2 + small_schooling_fish + large_predatory_fish, data=SEM_data_std)
summary(multreg_std)
# the standardized coefficients give an estimate of the weight of each path in this diagram
# estimate values
# distance 0.06
# habitat -0.22
# river -0.17
# chemical water quality 0.09
# visible water quality 0.27 *
# small schooling fish -0.25 *
# large predatory fish 0.16
# These results suggest that visible water quality and small schooling fish has the most important effect on number of pied kingfishers

# However the multiple regression tests how much variation is explained by a predictor that is not explained by the other predictors in the model. So when two predictors are strongly correlated, each may not explain variation in the response independently, why the regression of each with the response may be highly significant.
# path analysis (form of SEM) can analyse causal relations between the different predictors, instead of hypothesizing that they all directly affect the response variable

# make SEM model and fit with lavaan 
pied_model <- 'pied_kingfisher ~ distance + habitat_dummy + river_dummy + chemical_water_quality_PC1 + visible_water_quality_PC2 + small_schooling_fish + large_predatory_fish
small_schooling_fish ~ visible_water_quality_PC2 + chemical_water_quality_PC1 + large_predatory_fish
                  visible_water_quality_PC2 ~ river_dummy + distance
                  chemical_water_quality_PC1~ river_dummy + distance
                  large_predatory_fish ~ visible_water_quality_PC2 + chemical_water_quality_PC1
                  pied_kingfisher ~ habitat_dummy' 

pied_model.fit <- lavaan::sem(pied_model, data= SEM_data_std)

# goodness and badness of fit measures of the model
# the goodness of fit measures CFI and TLI should be \>0.9 for the model to be accepted, while the 'badness of fit' measures RMSEA and SRMR should be be \<0.1 for the model to be accepted
# if these parameters are not in these ranges, then your causal model is not supported by the data, and you test a different model. this can involve the inclusion of different predictors, or testing an alternative causal structure among the predictors.
