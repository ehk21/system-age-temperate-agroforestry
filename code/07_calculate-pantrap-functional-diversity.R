rm(list=ls())

# load packages
library(iNEXT.3D)
library(dplyr)
library(tidyr)
library(ggplot2)
library(tibble)
library(stringr)
library(cluster)

# load pan trap data
data <- read.csv("data/cleaned/01_pantraps-cleaned.csv")
bee_traits <- read.csv("data/cleaned/05_solitary-bee-traits-cleaned.csv")
syr_traits <- read.csv("data/cleaned/05_hoverfly-traits-cleaned.csv")

# remove 2012 data as only 24 specimens caught across both treatments
data <- subset(data, Year != 2012)

# remove sep and october observations
data <- subset(data, !(Month %in% c("September", "October")))

# remove May 2018 observations (discluded from Staton et al. due to high number of missing samples)
data <- subset(data, !(Year==2018 & Month=="May"))

# add season column
data$Season <- ifelse(data$Month %in% c("March", "April","May"), "Early", "Late")

##################
# solitary bees
##################

# subset pantrap data 
sol <- subset(data, TaxGroup=="Anthophila" & Genus != "Bombus", select=c(Year, Season, Month, Site, Treatment, Sample, Spp, Sex, Abundance))

# combine spp and sex
sol$Spp <- paste(sol$Spp, sol$Sex, sep=".")

# Combine identifiers into one community ID
sol <- sol %>%
  mutate(Community = paste(Year, Season, Site, Treatment, sep = "_"))

# Create community matrix (samples as rows, species as columns)
comm_matrix <- sol %>%
  group_by(Community, Spp) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Spp, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Community")

# ensure species names are consistent in traits & community df
bee_traits <- bee_traits %>%
  distinct(Species, .keep_all = TRUE) %>%
  filter(Species %in% colnames(comm_matrix)) %>%
  column_to_rownames("Species")

# subset to just traits of interest
bee_traits <- subset(bee_traits, select=c(Flight_period_.months., nesting_method, gens, ITD, Lecty))

bee_traits <- bee_traits %>%
  dplyr::select(Flight_period_.months., nesting_method, gens, ITD, Lecty) %>%
  mutate_if(is.character, as.factor)
  
# Calculate Gower distance (works with mixed data types)
gower_dist <- daisy(bee_traits, metric = "gower")
gower_matrix <- as.matrix(gower_dist)

# check sample sizes
sort(apply(comm_matrix, 1, function(z) sum(z > 0)))[1:10]
# 2019_Late_2_Control, 2019_Late_3_Control and 2019_Late_2_Agroforestry all have <5 observed species ... Need to remove or function won't run.
comm_matrix <- comm_matrix[-which(grepl("2019_Late_2|2019_Late_3", rownames(comm_matrix))),]

# run iNEXT.3D
fd_iNEXT <- iNEXT3D(
  data = as.data.frame(t(comm_matrix)),
  diversity = "FD",
  q = c(0,1,2),
  datatype = "abundance",
  FDdistM = gower_matrix,
  FDtype = "AUC",
  nboot = 50
)

# extract info
DataInfo <- fd_iNEXT[["FDInfo"]]

# calculate base coverage (i.e., min sample coverage of reference samples) using Chao & Jost, 2012
base.coverage <- min(DataInfo$"SC(n)"[DataInfo$"SC(n)">0.5])

# remove rows with sample coverage <0.5 (and paired sample points)
remove <- data.frame(Assemblage=DataInfo$Assemblage[DataInfo$"SC(n)"<0.5]) %>%
  separate(Assemblage, into = c("Year", "Season", "Site", "Treatment"), sep = "\\_") %>%
  mutate(Pairs = paste(Year, Season, Site, sep="."))

if(nrow(remove)>0){
  comm_matrix <- comm_matrix[-which(grepl(paste(remove$Pairs, collapse="|"), DataInfo$Assemblage)),]
}

# calculate diversity metrics at base coverage 
sol.fd.baseline <- estimate3D(
  data = as.data.frame(t(comm_matrix)),
  diversity = "FD",
  q = c(0,1,2),
  datatype = "abundance",
  base="coverage",
  # level=base.coverage, #0.66
  FDdistM = gower_matrix,
  FDtype = "AUC",  
  nboot = 50
)

# extract group columns
sol.fd.baseline <- sol.fd.baseline %>%
  separate(Assemblage, into = c("Year", "Season", "Site","Treatment"), sep = "\\_")

# convert 95% confidence interval to standard error (assuming normal distribution: 95% confidence interval is the mean +- 1.96*SE)
sol.fd.baseline <- sol.fd.baseline %>%
  mutate(SE = (qFD.UCL - qFD.LCL) / (2 * 1.96))

# add age column
sol.fd.baseline$Year <- as.numeric(sol.fd.baseline$Year)
origin_year <- c(`1` = 2014, `2` = 2015, `3` = 2009)
sol.fd.baseline$Age <- sol.fd.baseline$Year - origin_year[as.character(sol.fd.baseline$Site)]

# plot
ggplot(sol.fd.baseline, aes(x=Age, y=qFD, color=Treatment))+
  geom_point(aes(shape=as.factor(Site)), size=3) +
  geom_errorbar(aes(ymin=qFD-SE, ymax=qFD+SE), width=.4) +
  facet_grid(Season~Order.q, scales="free") +
  geom_smooth(method="lm", se=F)+
  theme_bw() +
  ggtitle("Functional diversity of solitary bees with system age") +
  ylab("Functional Diversity")

##################
# hoverflies
##################

# subset pantrap data 
syr <- subset(data, TaxGroup=="Syrphidae", select=c(Year, Season, Month, Site, Treatment, Sample, Spp, Sex, Abundance))

# replace 'merodontini strigatus with eumerus strigatus (to match traits df)
syr <- syr %>%
  mutate(Spp = dplyr::recode(Spp, "Merodontini strigatus"="Eumerus strigatus"))

# format traits dataframe
syr_traits <- syr_traits %>%
  mutate_if(is.character, as.factor)

# Combine identifiers into one community ID
syr <- syr %>%
  mutate(Community = paste(Year, Season, Site, Treatment, sep = "_"))

# Create community matrix (samples as rows, species as columns)
comm_matrix <- syr %>%
  group_by(Community, Spp) %>%
  summarise(Abundance = sum(Abundance), .groups = "drop") %>%
  pivot_wider(names_from = Spp, values_from = Abundance, values_fill = 0) %>%
  column_to_rownames("Community")

# Ensure species names are consistent
syr_traits <- syr_traits %>%
  distinct(Spp, .keep_all = TRUE) %>%
  filter(Spp %in% colnames(comm_matrix)) %>%
  column_to_rownames("Spp")

# Calculate Gower distance (works with mixed data types)
gower_dist <- daisy(syr_traits, metric = "gower")
gower_matrix <- as.matrix(gower_dist)

# remove unidentified specimens
rm <- which(colnames(comm_matrix)=="Unidentified bacchini genus.")
comm_matrix <- comm_matrix[,-rm]
comm_matrix <- as.matrix(comm_matrix)

# check sample sizes
sort(apply(comm_matrix, 1, function(z) sum(z > 0)))

# remove sites with <5 observed spp
comm_matrix <- comm_matrix[-which(grepl(c("2019_Early_1|2019_Early_2|2019_Early_3|2023_Early_1|2023_Early_2|2023_Early_3|2019_Late_1|2019_Late_2"), rownames(comm_matrix))),]

# run iNEXT.3D
fd_iNEXT <- iNEXT3D(
  data = as.data.frame(t(comm_matrix)),
  diversity = "FD",
  q = c(0,1,2),
  datatype = "abundance",
  FDdistM = gower_matrix,
  FDtype = "AUC",
  nboot = 50,
  conf=0.95
)

# extract info
DataInfo <- fd_iNEXT[["FDInfo"]]
# diversity_data <- as.data.frame(fd_iNEXT[["iNextEst"]][["coverage_based"]])

# calculate base coverage (i.e., min sample coverage of reference samples) using Chao & Jost, 2012
base.coverage <- min(DataInfo$"SC(n)"[DataInfo$"SC(n)">0.5])

# remove rows with sample coverage <0.5 (and paired sample points)
remove <- data.frame(Assemblage=DataInfo$Assemblage[DataInfo$"SC(n)"<0.5]) %>%
  separate(Assemblage, into = c("Year", "Season", "Site", "Treatment"), sep = "\\_") %>%
  mutate(Pairs = paste(Year, Season, Site, sep="."))

if(nrow(remove)>0){
  comm_matrix <- comm_matrix[-which(grepl(paste(remove$Pairs, collapse="|"), DataInfo$Assemblage)),]
}

# calculate diversity metrics at base coverage
syr.fd.baseline <- estimate3D(
  data = as.data.frame(t(comm_matrix)),
  diversity = "FD",
  q = c(0,1,2),
  datatype = "abundance",
  base="coverage",
  # level=base.coverage, #0.82
  FDdistM = gower_matrix,
  FDtype = "AUC",  
  nboot = 50
)

# extract group columns
syr.fd.baseline <- syr.fd.baseline %>%
  separate(Assemblage, into = c("Year", "Season", "Site", "Treatment"), sep = "\\_")

# convert 95% confidence interval to standard error (assuming normal distribution: 95% confidence interval is the mean +- 1.96*SE)
syr.fd.baseline <- syr.fd.baseline %>%
  mutate(SE = (qFD.UCL - qFD.LCL) / (2 * 1.96))

# add age column
syr.fd.baseline$Year <- as.numeric(syr.fd.baseline$Year)
origin_year <- c(`1` = 2014, `2` = 2015, `3` = 2009)
syr.fd.baseline$Age <- syr.fd.baseline$Year - origin_year[as.character(syr.fd.baseline$Site)]

# plot
ggplot(syr.fd.baseline, aes(x=Age, y=qFD, color=Treatment))+
  geom_point(aes(shape=as.factor(Site)), size=3) +
  geom_errorbar(aes(ymin=qFD-SE, ymax=qFD+SE), width=.4) +
  facet_grid(Order.q~Season, scales="free") +
  geom_smooth(method="lm", se=F)+
  theme_bw() +
  ggtitle("Functional diversity of hoverfly with system age") +
  ylab("Functional Diversity")

################
# save
################
sol.func <- subset(sol.fd.baseline, select=c(Year, Season, Site, Age, Treatment, Order.q, qFD, SE))
syr.func <- subset(syr.fd.baseline, select=c(Year, Season, Site, Age, Treatment, Order.q, qFD, SE))

write.csv(sol.func, "data/cleaned/07_solitary-func-div.csv", row.names=F)
write.csv(syr.func, "data/cleaned/07_syrphidae-func-div.csv", row.names=F)


