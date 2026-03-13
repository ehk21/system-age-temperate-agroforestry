rm(list=ls())

# load packages ----
library(dplyr)
library(plotrix)

# load data ----
ek <- read.csv("data/raw/01_pantraps-raw-EK.csv")
ts <- read.csv("data/raw/01_pantraps-raw-TS.csv")
av <- read.csv("data/raw/01_pantraps-raw-AV.csv")

##########################################
# format and wrangle own pan trap data ----
##########################################

ek <- ek %>%
  rename(Spp=Species) %>% # change column header to match other datasets
  mutate(across(c(Treatment, Alley, Sample, Sex), toupper)) #capitalise

ek$Date <- as.Date(ek$Date, format = "%d/%m/%Y") #extract month and year
ek$Month <- format(ek$Date, "%B")
ek$Year <- format(ek$Date, "%Y")

##check everything entered correctly
unique(ek$Site) 
unique(ek$Treatment)

unique(ek$Alley) # need to remove margin pan traps (as these were removed from sample design after April) and correct spelling
ek <- subset(ek, !(Alley %in% c("M1", "M2", "M3")))
ek$Alley[ek$Alley=="S1"]<-"A1"
ek$Alley[ek$Alley=="S2"]<-"A2"

unique(ek$Sample) #correct S6.2 to S6. One sample without number recorded so will default to S1 as I don't think I will analyse to this resolution anyway.
ek$Sample[ek$Sample=="S6.2"] <- "S6"
ek$Sample[ek$Sample=="S-"] <- "S1"

unique(ek$Spp)
ek$Spp[ek$Spp=="Andrena chrysoscles" ] <- "Andrena chrysosceles"   #correct spelling
ek$Spp[ek$Spp=="Chrystotoxum bicinctum"] <- "Chrysotoxum bicinctum" 
ek$Spp[ek$Spp=="Andrena nigroanea" ] <- "Andrena nigroaenea"  
ek$Spp[ek$Spp=="Syrphus ribesi" ] <- "Syrphus ribesii"  

ek$Spp[ek$Spp %in% c("Sphaerophoria ", "Unidentified Spherophoria spp.")] <- "Sphaerophoria scripta" #assume female sphaerophoria and unidentified males are all scripta as all males were found to be scripta.
ek <- subset(ek, Spp != "Unidentifiable") #remove unidentifiable specimen from data
ek$Spp[ek$Spp %in% c("Unidenfitied Eupeodes genus.")] <- "Eupeodes corollae" # assign most common eupeodes spp at that site
ek$Spp[ek$Spp %in% c("Unidentified bacchini genus.") & ek$Site=="RH"] <- "Melanostoma mellinum" # assign most common Bacchinni spp at site
ek$Spp[ek$Spp %in% c("Unidentified bacchini genus.") & ek$Site=="WH"] <- "Platycheirus manicatus" # assign most common Bacchini spp at site
ek$Spp[ek$Spp == "Bombus terrestris/locorum"] <- "Bombus lucorum/terrestris" #correct spelling

unique(ek$TaxGroup)

# ## assign sex to bees but not hoverflies
# ## as functional differences can occur between sexes in the former but not generally the latter
# ek$Sex[ek$TaxGroup=="Syrphidae"]<-NA
# ek$Spp[ek$TaxGroup=="Anthophila"] <- paste0(ek$Spp[ek$TaxGroup=="Anthophila"]," (", ek$Sex[ek$TaxGroup=="Anthophila"],")")

## rename treatment variables
ek$Treatment[ek$Treatment == "AF"] <- "Agroforestry"
ek$Treatment[ek$Treatment == "C"] <- "Control"

## create genus and species column
ek$Genus <- sapply(strsplit(ek$Spp, " "), `[`, 1)
ek$Species <- sapply(strsplit(ek$Spp, " "), `[`, 2)

## Combine alley and sample codes to create unique sample per round
ek$Sample <- paste0(ek$Alley, ek$Sample)
ek <- subset(ek, select=c(Year, Month, Site, Treatment,  Sample, TaxGroup, Genus, Species, Spp, Sex))

## summarise abundance of each species per sample, month, site, treatment
ek <- ek %>% 
  group_by(Year, Month, Site, Treatment,  Sample, TaxGroup, Genus, Species, Spp, Sex) %>%
  summarise(Abundance=n())

##########################################
# format TS pan trap data ----
##########################################

##rename cols & replace variable names
ts <- ts %>% 
  rename(Spp=Species, Treatment=FarmingSystem) %>%
  mutate(Treatment=dplyr::recode(Treatment, "Arable"="Control"),
         Month=dplyr::recode(Month, "4"="April", "5"="May", "6"="June", "7"="July", "8"="August", "9"="September", "10"="October")
         )

##remove NAs
ts <- subset(ts, !is.na(ts$Spp)) 

##check each variable
unique(ts$Site) 
unique(ts$Treatment)
unique(ts$Month) 

unique(ts$Year) #convert to character to match other two dataframes
ts$Year <- as.character(ts$Year)

unique(ts$Spp)
ts <- subset(ts, Spp != "Apis mellifera") #remove honeybees
ts$Spp[ts$Spp=="Bombus lucorumterrestris"] <- "Bombus lucorum/terrestris"

unique(ts$Sex) #remove unknowns
ts$Sex[ts$Sex=="Unknown"] <- NA
ts$Sex[ts$TaxGroup=="Syrphidae"]<-NA

unique(ts$TaxGroup)

## create genus and species column
ts$Genus <- sapply(strsplit(ts$Spp, " "), `[`, 1)
ts$Species <- sapply(strsplit(ts$Spp, " "), `[`, 2)

# ## add sex to bee species names
# ts$Spp[ts$TaxGroup=="Anthophila"] <- paste0(ts$Spp[ts$TaxGroup=="Anthophila"]," (", ts$Sex[ts$TaxGroup=="Anthophila"],")")

# ## remove May 2018 because Staton et al. omitted due to a high number of missing samples in this month
# ts <- ts[!(ts$Year=="2018" & ts$Month=="May"),] # keeping for now as should be dealt with during rarefaction

# replace sample column name
names(ts)[which(grepl("SampleRef",names(ts)))] <- "Sample"

# remove unecessary columns
ts <- subset(ts, select=c(Year, Month, Site, Treatment,  Sample, TaxGroup, Genus, Species, Spp, Sex, Abundance))

##########################################
# format AV's pan trap data ----
##########################################

## check data
unique(av$Site) 
av <- subset(av, Site=="WH", select=c(Site, Treatment, Date, Repeat,Sample.station, Genus, Species, sp, Sex))
av$Site <- dplyr::recode(av$Site, "WH"=3)

unique(av$Treatment)
av <- av %>%
  mutate(Treatment = dplyr::recode(Treatment, "AF"="Agroforestry", "C"="Control"))

av$Date <- as.Date(av$Date, format = "%d/%m/%Y") #extract month and year
av$Month <- format(av$Date, "%B")
av$Year <- format(av$Date, "%Y")

names(av)[which(grepl("sp", names(av)))] <- "Spp"

unique(av$Genus)

unique(av$Species) 
av$Species[av$Species=="nigroanea"] <- "nigroaenea" # correct spelling

unique(av$Spp)
av$Spp[av$Spp=="Andrena nigroanea"] <- "Andrena nigroaenea" # correct spelling

av$Spp[av$Spp=="Lasioglossum semilaevis"] <- "Lasioglossum semilucens" #correct misnomer
av$Species[av$Species=="semilaevis"] <- "semilucens"

av$Spp[av$Spp=="Lasioglossum minutula"] <- "Lasioglossum minutissimum" #correct misnomer
av$Species[av$Spp=="Lasioglossum minutissimum"] <- "minutissimum"

unique(av$Sex) # convert to uppercase
av$Sex <- toupper(av$Sex)

# # combine species and sex
# av$Spp <- paste0(av$Spp, " (", av$Sex, ")")

av$TaxGroup <- "Anthophila"

# create sample column combining repeat(i.e. alley) and sample.station
av$Sample.station <- as.factor(av$Sample.station)
av$Sample <- paste0(av$Repeat, av$Sample.station)

av <- av %>% #summarise abundance of each species per month, site, treatment
  group_by(Year, Month, Site, Treatment,  Sample, TaxGroup, Genus, Species, Spp, Sex) %>%
  summarise(Abundance=n())

#####################
# combine dataframes ----
#####################
wrangled_df <- rbind(ek, ts, av)

#####################
# export wrangled data ----
#####################
write.csv(wrangled_df, "data/cleaned/01_pantraps-cleaned.csv", row.names=F)

