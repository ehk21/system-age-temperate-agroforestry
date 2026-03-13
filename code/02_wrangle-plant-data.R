rm(list=ls())

# load packages
library(dplyr)
library(tidyr)

# load data
av <- read.csv("data/raw/02_plants-raw-AV.csv")
ts <- read.csv("data/raw/02_plants-raw-TS.csv")
ek <- read.csv("data/raw/02_plants-raw-EK.csv")

###################
# wrangle AV data
###################

# convert to lengthwise
av.long <- av %>%
  pivot_longer(names_to="Taxon", values_to = "Percentage", cols=!(c("ID", "Date", "Site", "Treatment", "Field", "Transect", "Quadrat", "Distance.up.alley", "Distance.across.alley", "Crop.or.tree.strip")))

# format
av.long$Year <- 2011

av.long <- subset(av.long, select=c(Year, Site, Treatment, ID, Crop.or.tree.strip, Taxon, Percentage))
names(av.long)[which(names(av.long)=="Crop.or.tree.strip")]<- "Position"

av.long$Treatment <- dplyr::recode(av.long$Treatment, "AF"="Agroforestry", "C"="Control")
av.long$Site <- dplyr::recode(av.long$Site, "wh"=3)
av.long <- av.long[av.long$Percentage>0,]

# spellcheck species
av.long$Taxon <- gsub("\\.", " ", av.long$Taxon)

av.long$Taxon[av.long$Taxon=="Brassica sp"] <- "Brassica sp."
av.long$Taxon[av.long$Taxon=="Brasica napus"] <- "Brassica napus"
av.long$Taxon[av.long$Taxon=="Ploygonum aviculare"] <- "Polygonum aviculare"
av.long$Taxon[av.long$Taxon=="Fallopia convulvulus"] <- "Fallopia convolvulus"
av.long$Taxon[av.long$Taxon=="Cirsium sp"] <- "Cirsium sp."
av.long$Taxon[av.long$Taxon=="Galeopsis sp"] <- "Galeopsis sp."
av.long$Taxon[av.long$Taxon=="Capsella bursapastoris"] <- "Capsella bursa-pastoris"

# add common name column
common_names <- c(
  "Brassica sp." = "Mustard species",
  "Brassica napus" = "Oilseed rape",
  "Capsella bursa-pastoris" = "Shepherd's purse",
  "Centaurea nigra" = "Common knapweed",
  "Chenopodium album" = "Fat-hen",
  "Leucanthemum vulgare" = "Oxeye daisy",
  "Silene dioica" = "Red campion",
  "Solanum tuberosum" = "Potato",
  "Stellaria media" = "Common chickweed",
  "Trifolium repens" = "White clover",
  "Festuca pratensis" = "Meadow fescue",
  "Lamium album" = "White dead-nettle",
  "Rumex crispus" = "Curled dock",
  "Festuca rubra" = "Red fescue",
  "Poa trivialis" = "Rough meadow-grass",
  "Lotus corniculatus" = "Bird’s-foot trefoil",
  "Medicago lupulina" = "Black medick",
  "Polygonum aviculare" = "Knotgrass",
  "Trifolium hybridum" = "Alsike clover",
  "Alopecurus geniculatus" = "Marsh foxtail",
  "Trifolium pratense" = "Red clover",
  "Cirsium sp." = "Thistle species",
  "Fallopia convolvulus" = "Black-bindweed",
  "Cirsium arvense" = "Creeping thistle",
  "Galeopsis sp." = "Hemp-nettle species",
  "Viola arvensis" = "Field pansy",
  "Galium aparine" = "Cleavers"
)

av.long$Common.name <- common_names[av.long$Taxon]

# av.long <- av.long[!is.na(av.long$Common.name),]

av.long$Method <- "quadrat"

########################
# Wrangle TS data
########################

ts$FarmingSystem[ts$FarmingSystem=="Arable"] <- "Control"

names(ts)[which(names(ts)=="FarmingSystem")]<-"Treatment"
names(ts)[which(names(ts)=="SampleRef")]<-"ID"

ts$Position <- "c"
ts$Position[ts$Method=="us"]<-"t"

ts$ID <- paste(ts$Site, ts$ID, ts$Year, ts$Method, sep=".")

ts <- subset(ts, select=c(Year, Site, Treatment, ID, Position, Method, Taxon, Percentage))

# add common names
latin_to_common <- c(
  "Galium aparine" = "Cleavers",
  "Veronica persica" = "Persian speedwell",
  "Bryophyte" = "Moss sp.",  # Group, no common name
  "Poaceae" = "Grass sp.",    # Family
  "Picris echioides" = "Bristly oxtongue",
  "Tripleurospermum inodorum" = "Scentless mayweed",
  "Leucanthemum vulgare" = "Oxeye daisy",
  "Alopecurus myosuroides" = "Black-grass",
  "Epilobium" = "Willowherb (unspecified)",
  "Geranium dissectum" = "Cut-leaved crane's-bill",
  "Euphorbia" = "Spurge (unspecified)",
  "Dicot" = "Dicot sp.",  # Group
  "Triticum" = "Wheat",
  "Convolvulus arvensis" = "Field bindweed",
  "Stellaria media" = "Common chickweed",
  "Viola arvensis" = "Field pansy",
  "Papaver rhoeas" = "Common poppy",
  "Stachys officinalis" = "Betony",
  "Cirsium arvense" = "Creeping thistle",
  "Lamiaceae" = "Mint sp.",  # Family
  "Atriplex/Chenopodium" = "Goosefoot sp.",  # Mixed/ambiguous
  "Bromus sterilis" = "Barren brome",
  "Persicaria lapathifolia" = "Pale persicaria",
  "Myosotis arvensis" = "Field forget-me-not",
  "Poa trivialis" = "Rough meadow grass",
  "Phleum pratense" = "Timothy grass",
  "Trifolium" = "Clover (unspecified)",
  "Sonchus oleraceus" = "Smooth sow-thistle",
  "Trifolium pratense" = "Red clover",
  "Taraxacum agg." = "Dandelion (aggregate)",
  "Avena" = "Oat (unspecified)",
  "Rumex crispus" = "Curled dock",
  "Sonchus asper" = "Prickly sow-thistle",
  "Senecio vulgaris" = "Groundsel",
  "Ranunculus acris" = "Meadow buttercup",
  "Lactuca serriola" = "Prickly lettuce",
  "Sinapis arvensis" = "Charlock",
  "Lamium" = "Dead-nettle (unspecified)",
  "Silene latifolia" = "White campion",
  "Bromus hordeaceus" = "Soft brome",
  "Anthriscus sylvestris" = "Cow parsley",
  "Heracleum sphondylium" = "Hogweed",
  "Festuca rubra" = "Red fescue",
  "Lamium album" = "White dead-nettle",
  "Trifolium repens" = "White clover",
  "Rumex obtusifolius" = "Broad-leaved dock",
  "Urtica dioica" = "Stinging nettle",
  "Hordeum murinum" = "Wall barley",
  "Stachys sylvatica" = "Hedge woundwort",
  "Geranium molle" = "Dove's-foot crane's-bill",
  "Arrhenatherum elatius" = "False oat-grass",
  "Arctium lappa" = "Greater burdock",
  "Hordeum brachyantherum" = "Meadow barley",
  "Alopecurus pratensis" = "Meadow foxtail",
  "Dactylis glomerata" = "Cock's-foot",
  "Epilobium tetragonum" = "Square-stalked willowherb",
  "Plantago lanceolata" = "Ribwort plantain",
  "Epilobium hirsutum" = "Great willowherb",
  "Jacobaea vulgaris" = "Common ragwort",
  "Epilobium parviflorum" = "Hoary willowherb",
  "Rumex conglomeratus" = "Clustered dock",
  "Holcus lanatus" = "Yorkshire fog",
  "Artemisia vulgaris" = "Mugwort",
  "Medicago lupulina" = "Black medick",
  "Cirsium vulgare" = "Spear thistle",
  "Lolium perenne" = "Perennial ryegrass",
  "Lotus corniculatus" = "Bird’s-foot trefoil",
  "Achillea millefolium" = "Yarrow",
  "Trifolium dubium" = "Lesser trefoil",
  "Vicia tetrasperma" = "Smooth tare",
  "Centaurea nigra" = "Common knapweed",
  "Malva moschata" = "Musk mallow",
  "Cynosurus cristatus" = "Crested dog’s-tail",
  "Prunella vulgaris" = "Self-heal",
  "Galium verum" = "Lady's bedstraw"
)

ts$Common.name <- latin_to_common[ts$Taxon]

###################
# Wrangle 2023 data
###################

ek <- subset(ek, month %in% c("june","july"))

# add habitat column then remove margin transects
ek$transect <- substr(ek$quadrat, start=1, stop=1)

ek$habitat <- NA

ek$habitat[ek$site=="FE"&ek$treatment=="AF"&ek$transect %in% c("1", "3", "5", "7")]<- "treerow"
ek$habitat[ek$site=="FE"&ek$treatment=="AF"&ek$transect %in% c("2", "4", "6")]<- "crop"
ek$habitat[ek$site=="FE"&ek$treatment=="C"&ek$transect %in% c("3", "4", "5", "6")]<- "crop"
ek$habitat[ek$site=="FE"&ek$treatment=="AF"&ek$transect %in% c("8")]<- "clover ley"
ek$habitat[ek$site=="FE"&ek$treatment=="C"&ek$transect %in% c("1", "2")]<- "buckwheat"
ek$habitat[ek$site%in%c("WH","RH")&ek$treatment=="AF"&ek$transect %in% c("1", "3", "5")]<- "crop"
ek$habitat[ek$site%in%c("WH", "RH")&ek$treatment=="AF"&ek$transect %in% c("2", "4", "6")]<- "treerow"
ek$habitat[ek$site%in%c("WH","RH")&ek$treatment=="C"&!grepl("M", ek$transect)]<- "crop"
ek$habitat[grepl("M", ek$transect)] <- "margin"

ek <- ek[ek$habitat !="margin",]
ek$habitat[ek$habitat=="crop"]<-"c"
ek$habitat[ek$habitat=="treerow"]<-"t"
ek$habitat[ek$habitat=="buckwheat"]<-"c"
ek$habitat[ek$habitat=="clover ley"]<-"c"

# format
ek$ID <- paste(ek$site, ek$treatment, ek$quadrat, sep=".")

ek <- subset(ek, select=c(site, treatment, ID, habitat, spp, spp_perc))

ek$spp_perc[which(ek$spp_perc=="na")] <- 0
ek$spp_perc[which(ek$spp_perc == "<1")] <- 0.5
ek$spp_perc[which(ek$spp_perc == "")] <- 0

ek$site <- dplyr::recode(ek$site, "FE"=1, "RH"=2, "WH"=3)
ek$treatment <- dplyr::recode(ek$treatment, "AF"="Agroforestry", "C"="Control")

ek$Year <- 2023

names(ek)<-c("Site", "Treatment", "ID", "Position", "Common.name", "Percentage", "Year")

ek <- ek[ek$Common.name!="na",]

ek$Method <- "quadrat"

# correct
ek$Common.name[ek$Common.name=="red dead nettle"]<- "red dead-nettle"
ek$Common.name[ek$Common.name=="field speedwell"]<- "common field speedwell"
ek$Common.name[ek$Common.name=="brassica"]<- "mustard var."

# add latin names
common_to_latin <- c(
  "red dead-nettle" = "Lamium purpureum",
  "common field speedwell" = "Veronica persica",
  "common chickweed" = "Stellaria media",
  "common nettle" = "Urtica dioica",
  "groundsel" = "Senecio vulgaris",
  "scentless mayweed" = "Tripleurospermum inodorum",
  "oilseed rape" = "Brassica napus",
  "cow parsley" = "Anthriscus sylvestris",
  "field pansy" = "Viola arvensis",
  "shepherd's purse" = "Capsella bursa-pastoris",
  "meadow buttercup" = "Ranunculus acris",
  "forget-me-not" = "Myosotis arvensis",  # most common
  "dock" = "Rumex sp.",  # ambiguous
  "oxeye daisy" = "Leucanthemum vulgare",
  "cut-leaved cranesbill" = "Geranium dissectum",
  "common knapweed" = "Centaurea nigra",
  "selfheal" = "Prunella vulgaris",
  "creeping thistle" = "Cirsium arvense",
  "white clover" = "Trifolium repens",
  "field poppy" = "Papaver rhoeas",
  "knotgrass" = "Polygonum aviculare",
  "hogweed" = "Heracleum sphondylium",
  "goosegrass" = "Galium aparine",
  "mustard var." = "Sinapis arvensis",  # common assumption
  "bristly oxtongue" = "Picris echioides",
  "birdsfoot trefoil" = "Lotus corniculatus",
  "common ragwort" = "Jacobaea vulgaris",
  "fool's parsley" = "Aethusa cynapium",
  "scarlet pimpernel" = "Anagallis arvensis",
  "red clover" = "Trifolium pratense",
  "field speedwell" = "Veronica persica",  # same as common field speedwell
  "wild basil" = "Clinopodium vulgare",
  "white dead-nettle" = "Lamium album",
  "pale persicaria" = "Persicaria lapathifolia",
  "prickly sow thistle" = "Sonchus asper",
  "charlock" = "Sinapis arvensis",
  "musk mallow" = "Malva moschata",
  "narrow-leaved birdsfoot trefoil" = "Lotus tenuis",
  "black medick" = "Medicago lupulina",
  "brassica" = "Brassica sp.",  # vague
  "buckwheat" = "Fallopia convolvulus",
  "creeping buttercup" = "Ranunculus repens",
  "redshank" = "Persicaria maculosa"
)

ek$Taxon <- common_to_latin[ek$Common.name]

ek <- ek[1:nrow(ek)-1,]

#####################
# combine datasets
#####################
flora <- rbind(av.long, ts, ek)

#####################
# spellcheck
#####################
flora$Common.name <- tolower(flora$Common.name)
flora$Common.name[flora$Common.name=="prickly sow-thistle"] <- "prickly sow thistle"
flora$Common.name[flora$Common.name=="self-heal"] <- "selfheal"
flora$Common.name[flora$Common.name=="cut-leaved crane's-bill"] <- "cut-leaved cranesbill"
flora$Common.name[flora$Common.name=="rough meadow-grass"] <- "rough meadow grass"
flora$Common.name[flora$Common.name=="stinging nettle"] <- "common nettle"
flora$Common.name[flora$Common.name=="forget-me-not"] <- "field forget-me-not"
flora$Common.name[flora$Common.name=="mustard species"] <- "mustard var."
flora$Common.name[flora$Common.name=="clover (unspecified)"] <- "mustard var."
flora$Common.name[flora$Common.name=="dandelion (aggregate))"] <- "dandelion"

#####################
# save
#####################
write.csv(flora, "data/cleaned/02_plants-cleaned.csv", row.names = F)






