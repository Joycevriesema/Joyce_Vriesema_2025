rm(list = ls())

# load libraries
library(dplyr)
library(exiftoolr)
library(fuzzyjoin)
library(exifr)

# load data_water_photos and filter out old data (5-Feb-2025 & 8-Feb-2025) and very incomplete data (tree_5, tree_6, pap_6, tree_7 and tree_8 on 13-Feb-2025)
data_water_photos <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1825209523&single=true&output=csv") %>%
  dplyr::filter(date != "5-Feb-2025" & date != "8-Feb-2025") %>%
  dplyr::filter(!(transect_ID == "pap_6" & date == "13-Feb-2025") & 
                  !(transect_ID == "tree_7" & date == "13-Feb-2025") & 
                  !(transect_ID == "tree_8" & date == "13-Feb-2025") &
                  !(transect_ID == "tree_5" & date == "13-Feb-2025") &
                  !(transect_ID == "tree_6" & date == "13-Feb-2025"))

# load exif data, filter for jpg and JPG files and filter out old data (5-Feb-2025 & 8-Feb-2025)
files <- list.files("/Users/joycevriesema/Library/CloudStorage/GoogleDrive-j.vriesema.1@student.rug.nl/.shortcut-targets-by-id/1siXV0T2vKuXT_W5l-u1rm5aMLGcwJOad/2025_LakeFishBirdTransects_SLVC/Fish Videos and Water Photos", 
                    full.names = TRUE, 
                    ignore.case = TRUE, 
                    recursive = TRUE)
files <- files[grepl("\\.(jpe?g)$", files, ignore.case = TRUE)]

# split DateTimeOriginal into date and time and filter out old data (3-Feb-2025, 5-Feb-2025 & 8-Feb-2025)
exif_data <- read_exif(files, tags = c("FileName", "DateTimeOriginal")) %>%
  dplyr::mutate(
    FileName = tolower(sub("\\.(?i:jpe?g)$", "", FileName)),
    date = substr(DateTimeOriginal, 1, 10),
    time = substr(DateTimeOriginal, 12, 19)
  ) %>%
  dplyr::select(-DateTimeOriginal) %>%
  dplyr::filter(!date %in% c("2025:02:03", "2025:02:05", "2025:02:08")) %>%
  dplyr::select(FileName, time) 


# make the IDs from the sheet to match the same format as the FileName values from the exif data
data_water_photos <- data_water_photos %>%
  dplyr::mutate(
    photo_ID_start = tolower(sub("\\.(?i:jpe?g)$", "", photo_ID_start)),
    photo_ID_end   = tolower(sub("\\.(?i:jpe?g)$", "", photo_ID_end))
  )

# merge data by photo_ID_start and rename time into start_time
data_water_photos <- data_water_photos %>%
  dplyr::left_join(exif_data, by = c("photo_ID_start" = "FileName")) %>%
  dplyr::rename(start_time = time)

# merge data by photo_ID_end and rename time into end_time
data_water_photos <- data_water_photos %>%
  dplyr::left_join(exif_data, by = c("photo_ID_end" = "FileName")) %>%
  dplyr::rename(end_time = time)

# count seconds between start_time and end_time and remove "secs" from time_diff_sec
data_water_photos <- data_water_photos %>%
  dplyr::mutate(start_time = as.POSIXct(start_time, format = "%H:%M:%S"),
                end_time = as.POSIXct(end_time, format = "%H:%M:%S"),
                time_diff_sec = difftime(end_time, start_time, units = "secs"))%>%
  dplyr::mutate(time_diff_sec = as.numeric(time_diff_sec),
                start_time = substr(start_time, 12, 19),
                end_time = substr(end_time, 12, 19))

# calculate how many meters of the transect per second 
data_water_photos <- data_water_photos %>%
  dplyr::mutate(meters_per_sec = 500 / time_diff_sec)

# load data_water and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_water <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=15673116&single=true&output=csv") %>%
  dplyr::filter(date != "5-Feb-2025" & date != "8-Feb-2025") 

# create date_time variable for both datasets
data_water_photos <- data_water_photos %>%
  dplyr::mutate(date_time_start = paste(date, start_time),
                date_time_end = paste(date, end_time)) %>%
  dplyr::select(-date, -start_time, -end_time, -photo_ID_end, -photo_ID_start)

data_water <- data_water %>%
  dplyr::mutate(time = as.POSIXct(time, format = "%H:%M:%S"),
                time = substr(time, 12, 19),
                date_time = paste(date, time))

# add transect_ID to observations that fall between start and end time of a specific date 
data_water <- fuzzy_left_join(
  x = data_water,
  y = data_water_photos,
  by = c("date_time" = "date_time_start", "date_time" = "date_time_end"),
  match_fun = list(`>=`, `<=`))

# remove rows that do not below to a certain transect
data_water <- data_water %>%
  dplyr::filter(!is.na(transect_ID))

# load data_transect and filter out old data (5-Feb-2025 & 8-Feb-2025)
data_transect <- read.csv("https://docs.google.com/spreadsheets/d/e/2PACX-1vRCwiQGeumB9AuvRjnobaDJLq76NWyPQrvnPdvP58Qxv5SGMt4LMKjxMQMREGnYdoIkO1oCfTOcqp1Z/pub?gid=1366617186&single=true&output=csv") %>%
  dplyr::filter(date != "5-Feb-2025" & date != "8-Feb-2025")%>%
  dplyr::select(transect_run_ID,run_ID, transect_ID, date, direction_fish)

# merge data_water and data_transect
data_water <- data_water %>%
  dplyr::left_join(data_transect, by = c("transect_ID", "run_ID", "transect_run_ID", "date"))

# calculate meters_moved
data_water <- data_water %>%
  dplyr::mutate(date_time = as.POSIXct(date_time, format = "%d-%b-%Y %H:%M:%S"),
                date_time_start = as.POSIXct(date_time_start, format = "%d-%b-%Y %H:%M:%S"),
                seconds_moved = difftime(date_time, date_time_start, units = "secs")) %>%
  dplyr::mutate(seconds_moved = as.numeric(seconds_moved)) %>%
  dplyr::mutate(meters_moved = seconds_moved * meters_per_sec)

# subtract meters moved from 500 if direction_fish is "to coast"
data_water <- data_water %>%
  dplyr::mutate(meters_moved = ifelse(direction_fish == "to coast", 500 - meters_moved, meters_moved))

# save data_water
write.csv(data_water, "data_water.csv", row.names = F)