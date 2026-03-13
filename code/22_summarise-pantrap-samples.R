rm(list=ls())

# load packages
library(dplyr)

# load data
data <- read.csv("data/cleaned/01_pantraps-cleaned.csv")

# n.b. only summarising specimens here and not number of pan traps because pantraps where nothing was found are not included in the dataframe

# remove 2012 data as only 24 specimens caught across both treatments
data <- subset(data, Year != 2012)

# remove sep and october observations (only sampled in some years, likely different communities to june-august season)
data <- subset(data, !(Month %in% c("September", "October")))
late <- subset(data, (Month %in% c("September", "October")))

# remove May 2018 observations (discluded from Staton et al. due to high number of missing samples)
data <- subset(data, !(Year==2018 & Month=="May"))

# add season column
data$Season <- "Late"
data$Season[data$Month %in% c("March", "April", "May")] <- "Early"

# separate into taxonomic group
anthophila <- subset(data, TaxGroup=="Anthophila", select=c(Year, Season, Month, Site, Treatment, Sample, Spp, Genus, Abundance))
syrphidae <- subset(data, TaxGroup=="Syrphidae", select=c(Year, Season, Month, Site, Treatment, Sample, Spp, Genus, Abundance))
solitaries <- subset(data, TaxGroup=="Anthophila" & !(Genus=="Bombus"), select=c(Year, Season, Month, Site, Treatment, Sample, Spp, Genus, Abundance))

# subset syrphidae to just late season
syrphidae <- subset(syrphidae, Season=="Late")

# summarise
n.anthophila <- sum(anthophila$Abundance) #2459
n.syrphidae <- sum(syrphidae$Abundance) #1563
n.solitaries <- sum(solitaries$Abundance) #2402

syrphidae.spp <- length(unique(syrphidae$Spp)) #31
solitaries.spp <- length(unique(solitaries$Spp)) #39

syrphidae.genus <- length(unique(syrphidae$Genus)) #22
solitaries.genus <- length(unique(solitaries$Genus)) #7

syrphidae.af.spp <- length(unique(syrphidae$Spp[syrphidae$Treatment=="Agroforestry"])) #24
solitaries.af.spp <- length(unique(solitaries$Spp[solitaries$Treatment=="Agroforestry"])) #37
syrphidae.c.spp <- length(unique(syrphidae$Spp[syrphidae$Treatment=="Control"])) #25
solitaries.c.spp <- length(unique(solitaries$Spp[solitaries$Treatment=="Control"])) #30

syrphidae.af.genus <- length(unique(syrphidae$Genus[syrphidae$Treatment=="Agroforestry"])) #17
solitaries.af.genus <- length(unique(solitaries$Genus[solitaries$Treatment=="Agroforestry"])) #6
syrphidae.c.genus <- length(unique(syrphidae$Genus[syrphidae$Treatment=="Control"])) #18
solitaries.c.genus <- length(unique(solitaries$Genus[solitaries$Treatment=="Control"])) #6

