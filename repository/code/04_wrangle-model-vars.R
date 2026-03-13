rm(list=ls())

# load packages
library(dplyr)

# load data
site_info <- read.csv("data/raw/04_site-survey-info.csv")
weather_seasonal <- read.csv("data/cleaned/03_weather-vars-seasonal-cleaned.csv")
weather_annual <- read.csv("data/cleaned/03_weather-vars-annual-cleaned.csv")

# update site_info names to match weather
names(site_info)[1:3] <- c("Site", "Year", "Treatment")

# update site_info treatments to match bee div data
site_info <- site_info %>%
  mutate(Treatment=dplyr::recode(Treatment, "af"="Agroforestry", "c"="Control"))

# combine annually-resolved model vars data
vars_annual <- left_join(site_info, weather_annual)

# combine seasonally-resolved model vars data
site_info_seasonal <- rbind(site_info, site_info)
site_info_seasonal$Season <- "Late"
site_info_seasonal$Season[1:nrow(site_info)] <- "Early"

vars_seasonal <- left_join(site_info_seasonal, weather_seasonal)

# output
write.csv(vars_annual, "data/cleaned/04_model-vars-annual.csv", row.names=F)
write.csv(vars_seasonal, "data/cleaned/04_model-vars-seasonal.csv", row.names=F)
