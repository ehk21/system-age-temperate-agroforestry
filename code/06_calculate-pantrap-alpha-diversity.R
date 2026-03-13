rm(list=ls())

# load packages
library(dplyr)
library(vegan)
library(tidyr)
library(ggplot2)
library(plotrix)
library(iNEXT)
library(gridExtra)

# load data
data <- read.csv("data/cleaned/01_pantraps-cleaned.csv")

# remove 2012 data as only 24 specimens caught across both treatments
data <- subset(data, Year != 2012)

# remove sep and october observations (only sampled in some years, likely different communities to june-august season)
data <- subset(data, !(Month %in% c("September", "October")))

# remove May 2018 observations (discluded from Staton et al. due to high number of missing samples)
data <- subset(data, !(Year==2018 & Month=="May"))

# add season column
data$Season <- "Late"
data$Season[data$Month %in% c("March", "April", "May")] <- "Early"

# separate into taxonomic group
syrphidae <- subset(data, TaxGroup=="Syrphidae", select=c(Year, Season, Month, Site, Treatment, Sample, Spp, Abundance))
solitaries <- subset(data, TaxGroup=="Anthophila" & !(Genus=="Bombus"), select=c(Year, Season, Month, Site, Treatment, Sample, Spp, Abundance))

# sum abundance for all sexes/castes and all samples per year, season, site and treatment
syrphidae <- syrphidae %>%
  group_by(Year, Season, Site, Treatment, Spp) %>%
  summarise(Abundance=sum(Abundance))

solitaries <- solitaries %>%
  group_by(Year, Season, Site, Treatment, Spp) %>%
  summarise(Abundance=sum(Abundance))

# pivot the dataframe to wide format
syrphidae.wide <- syrphidae %>%
  pivot_wider(names_from = Spp, values_from = Abundance, values_fill = 0)

solitaries.wide <- solitaries %>%
  pivot_wider(names_from = Spp, values_from = Abundance, values_fill = 0)

####################################
# solitary bees
####################################

# convert to matrix format
solitaries.mat <- as.matrix(solitaries.wide[5:ncol(solitaries.wide)])
rownames(solitaries.mat) <- paste(solitaries.wide$Year, solitaries.wide$Season, solitaries.wide$Site, solitaries.wide$Treatment, sep=".")

# check number of observed species per assemblage
sort(apply(solitaries.mat, 1, function(z) sum(z > 0)))

# run iNext
sol.div <- iNEXT(t(solitaries.mat), datatype="abundance", q=c(0,1,2))

# extract info
DataInfo <- sol.div[["DataInfo"]]
# diversity_data <- as.data.frame(sol.div[["iNextEst"]][["coverage_based"]])

# calculate base coverage (i.e., min sample coverage of reference samples) using Chao & Jost, 2012
base.coverage <- min(DataInfo$SC[DataInfo$SC>0.5])

# remove rows with sample coverage <0.5 (and paired sample points)
remove <- data.frame(Assemblage=DataInfo$Assemblage[DataInfo$SC<0.5]) %>%
  separate(Assemblage, into = c("Year", "Season", "Site", "Treatment"), sep = "\\.") %>%
  mutate(Pairs = paste(Year, Season, Site, sep="."))

if(nrow(remove>0)){
solitaries.mat <- solitaries.mat[-which(grepl(paste(remove$Pairs, collapse="|"), DataInfo$Assemblage)),]
}

# calculate diversity metrics at base coverage 
sol.baseline <- estimateD(
  t(solitaries.mat),
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  # level = base.coverage, # 0.7
  nboot = 50,
  conf = 0.95
)

# extract group columns
sol.baseline <- sol.baseline %>%
  separate(Assemblage, into = c("Year", "Season", "Site", "Treatment"), sep = "\\.")

# convert 95% confidence interval to standard error (assuming normal distribution: 95% confidence interval is the mean +- 1.96*SE)
sol.baseline <- sol.baseline %>%
  mutate(SE = (qD.UCL - qD.LCL) / (2 * 1.96))

# add age column (survey year - year of establishment)
sol.baseline$Year <- as.numeric(sol.baseline$Year)
origin_year <- c(`1` = 2014, `2` = 2015, `3` = 2009)
sol.baseline$Age <- sol.baseline$Year - origin_year[as.character(sol.baseline$Site)]

# rename Order.q column
sol.baseline <- sol.baseline %>%
  mutate(Order.q=dplyr::recode(Order.q,
                               "0"="Species Richness (q=0)",
                               "1"="Shannon Diversity (q=1)",
                               "2"="Simpson Diversity (q=2)"),
         Order.q = factor(Order.q,
                          levels=c("Species Richness (q=0)", 
                                 "Shannon Diversity (q=1)", 
                                 "Simpson Diversity (q=2)" )))

# plot
ggplot(sol.baseline, aes(x=Age, y=qD, color=Treatment))+
  geom_point(aes(shape=as.factor(Site)), size=3) +
  geom_errorbar(aes(ymin=qD-SE, ymax=qD+SE), width=.4) +
  facet_grid(Season~Order.q, scales="free") +
  geom_smooth(method="lm", se=F)+
  theme_bw() +
  ggtitle("Alpha diversity of solitary bees with system age") +
  ylab("Alpha Diversity")

#############
# hoverflies:
#############

# convert to matrix format
syrphidae.mat <- as.matrix(syrphidae.wide[5:ncol(syrphidae.wide)])
rownames(syrphidae.mat) <- paste(syrphidae.wide$Year, syrphidae.wide$Season, syrphidae.wide$Site, syrphidae.wide$Treatment, sep=".")

# check number of observed species per assemblage
sort(apply(syrphidae.mat, 1, function(z) sum(z > 0)))
# very few species observed in any early season assemblage. rarefaction/extrapolation unreliable. subset to just late season

syrphidae.mat <- syrphidae.mat[which(grepl("Late", rownames(syrphidae.mat))),]

# run iNext
syr.div <- iNEXT(t(syrphidae.mat), datatype="abundance", q=c(0,1,2))

# extract data
DataInfo <- syr.div[["DataInfo"]]
# diversity_data <- as.data.frame(syr.div[["iNextEst"]][["coverage_based"]])

# calculate base coverage (i.e., min sample coverage of reference samples) using Chao & Jost, 2012
base.coverage <- min(DataInfo$SC[DataInfo$SC>0.5])

# remove rows with sample coverage <0.5 (and paired sample points)
remove <- data.frame(Assemblage=DataInfo$Assemblage[DataInfo$SC<0.5]) %>%
  separate(Assemblage, into = c("Year", "Season", "Site", "Treatment"), sep = "\\.") %>%
  mutate(Pairs = paste(Year, Season, Site, sep="."))

if(nrow(remove)>0){
syrphidae.mat <- syrphidae.mat[-which(grepl(paste(remove$Pairs, collapse="|"), DataInfo$Assemblage)),]
}

# calculate diversity metrics at base coverage
syr.baseline <- estimateD(
  t(syrphidae.mat),
  q = c(0, 1, 2),
  datatype = "abundance",
  base = "coverage",
  # level = base.coverage, #0.7
  nboot = 50,
  conf = 0.95
)

# extract grouping columns
syr.baseline <- syr.baseline %>%
  separate(Assemblage, into = c("Year", "Season", "Site", "Treatment"), sep = "\\.")

# convert 95% confidence interval to standard error (assuming normal distribution: 95% confidence interval is the mean +- 1.96*SE)
syr.baseline <- syr.baseline %>%
  mutate(SE = (qD.UCL - qD.LCL) / (2 * 1.96))

# add age column
syr.baseline$Year <- as.numeric(syr.baseline$Year)
origin_year <- c(`1` = 2014, `2` = 2015, `3` = 2009)
syr.baseline$Age <- syr.baseline$Year - origin_year[as.character(syr.baseline$Site)]

# rename Order.q factors
syr.baseline <- syr.baseline %>%
  mutate(Order.q=dplyr::recode(Order.q,
                               "0"="Species Richness (q=0)",
                               "1"="Shannon Diversity (q=1)",
                               "2"="Simpson Diversity (q=2)"),
         Order.q = factor(Order.q,
                          levels=c("Species Richness (q=0)", 
                                   "Shannon Diversity (q=1)", 
                                   "Simpson Diversity (q=2)" )))

# plot
ggplot(subset(syr.baseline), aes(x=Age, y=qD, color=Treatment))+
  geom_point(aes(shape=Site), size=3) +
  geom_errorbar(aes(ymin=qD-SE, ymax=qD+SE), width=.4) +
  facet_grid(Order.q~Season, scales="free") +
  geom_smooth(method="lm", se=F)+
  theme_bw() +
  ggtitle("Alpha diversity of hoverflies with system age") +
  ylab("Alpha Diversity")

################
# save
################
sol.alpha.seasons <- subset(sol.baseline, select=c(Year, Season, Site, Age, Treatment, Order.q, qD, SE))
syr.alpha.seasons <- subset(syr.baseline, select=c(Year, Season, Site, Age, Treatment, Order.q, qD, SE))

write.csv(sol.alpha.seasons, "data/cleaned/06_solitary-alpha-div.csv", row.names=F)
write.csv(syr.alpha.seasons, "data/cleaned/06_syrphidae-alpha-div.csv", row.names=F)


