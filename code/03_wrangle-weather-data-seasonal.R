rm(list=ls())

# load packages
library(dplyr)
library(tidyr)
library(stringr)
library(purrr)

# load data
midlands.precip <- read.csv("data/raw/03_Midlands-precip.csv")
anglia.precip <- read.csv("data/raw/03_East-Anglia-precip.csv")
midlands.temp <- read.csv("data/raw/03_Midlands-temp.csv")
anglia.temp <- read.csv("data/raw/03_East-Anglia-temp.csv")

# apply to sites
site1.precip <- midlands.precip
site1.precip$Site <- 1

site2.precip <- anglia.precip
site2.precip$Site <- 2

site3.precip <- anglia.precip
site3.precip$Site <- 3

site1.temp <- midlands.temp
site1.temp$Site <- 1

site2.temp <- anglia.temp
site2.temp$Site <- 2

site3.temp <- anglia.temp
site3.temp$Site <- 3

# combine sites
precip <- rbind(site1.precip, site2.precip, site3.precip )
temp <- rbind(site1.temp, site2.temp, site3.temp)

# set defaults
survey_years <- c(2010:2023)
sites <- c(1, 2, 3)
seasons <- c("Early", "Late")

# Month columns expected in temp/precip
month_cols <- c("jan","feb","mar","apr","may","jun","jul","aug","sep","oct","nov","dec")

##########################
# Basic checks 
##########################
stopifnot(all(c("year","Site") %in% names(temp)))
stopifnot(all(c("year","Site") %in% names(precip)))

missing_temp_months <- setdiff(month_cols, names(temp))
missing_prec_months <- setdiff(month_cols, names(precip))
if (length(missing_temp_months) > 0) stop("temp is missing month columns: ", paste(missing_temp_months, collapse=", "))
if (length(missing_prec_months) > 0) stop("precip is missing month columns: ", paste(missing_prec_months, collapse=", "))

# Ensure numeric months (Met data sometimes read as character)
temp[month_cols]   <- lapply(temp[month_cols], function(x) as.numeric(as.character(x)))
precip[month_cols] <- lapply(precip[month_cols], function(x) as.numeric(as.character(x)))

####################################################
# Compute seasonal and previous seasonal mean temp 
####################################################
get_season_temp <- function(df, site, yr, season) {
  # df is temp or precip
  row_this <- df %>% filter(Site == site, year == yr)
  row_prev <- df %>% filter(Site == site, year == yr-1)
  
  if (nrow(row_this) != 1) return(NULL)
  if (nrow(row_prev) != 1) return(NULL)
  
  if(season=="Early"){
    tibble(
      Temp = mean(c(row_this$mar, row_this$apr, row_this$may), na.rm = TRUE),
      Temp_prev  = mean(c(row_prev$dec, row_this$jan, row_this$feb), na.rm = TRUE)
    )} else if(season=="Late"){
      tibble(
        Temp = mean(c(row_this$jun, row_this$jul, row_this$aug), na.rm = TRUE),
        Temp_prev = mean(c(row_this$mar, row_this$apr, row_this$may), na.rm = TRUE)
      )}
}

temp_seasons <- expand_grid(Site = sites, Year = survey_years, Season=seasons) %>%
  mutate(tmp = pmap(list(Site, Year, Season),
                    ~ get_season_temp(temp, site = ..1, yr = ..2, season = ..3))) %>%
  unnest(tmp) 

####################################################
# compute seasonal and previous seasonal total raindays >1mm 
####################################################
get_season_precip <- function(df, site, yr, season) {
  # df is temp or precip
  row_this <- df %>% filter(Site == site, year == yr)
  row_prev <- df %>% filter(Site == site, year == yr-1)
  
  if (nrow(row_this) != 1) return(NULL)
  if (nrow(row_prev) != 1) return(NULL)
  
  if(season=="Early"){
    tibble(
      Precip = sum(c(row_this$mar, row_this$apr, row_this$may), na.rm = TRUE),
      Precip_prev  = sum(c(row_prev$dec, row_this$jan, row_this$feb), na.rm = TRUE)
    )} else if(season=="Late"){
      tibble(
        Precip = sum(c(row_this$jun, row_this$jul, row_this$aug), na.rm = TRUE),
        Precip_prev = sum(c(row_this$mar, row_this$apr, row_this$may), na.rm = TRUE)
      )}
}

precip_seasons <- expand_grid(Site = sites, Year = survey_years, Season=seasons) %>%
  mutate(tmp = pmap(list(Site, Year, Season),
                    ~ get_season_precip(precip, site = ..1, yr = ..2, season = ..3))) %>%
  unnest(tmp)

##########################
# combine
##########################
weather <- left_join(temp_seasons, precip_seasons)

####################################################
# calculate anomalies (difference from mean within 2011-2023 survey period)
####################################################
weather <- subset(weather, Year !=2010) %>%
  group_by(Site, Season) %>%
  mutate(
    Temp_anomaly        = Temp - mean(Temp, na.rm = TRUE),
    Precip_anomaly     = Precip - mean(Precip, na.rm = TRUE),
    Temp_prev_anomaly   = Temp_prev - mean(Temp_prev, na.rm = TRUE),
    Precip_prev_anomaly = Precip_prev - mean(Precip_prev, na.rm = TRUE)
  ) %>%
  ungroup()

# subset just to specific survey years
weather <- subset(weather, Year %in% c(2011, 2018, 2019, 2023))

##########################
# save
##########################
write.csv(weather, "data/cleaned/03_weather-vars-seasonal-cleaned.csv", row.names=F)
