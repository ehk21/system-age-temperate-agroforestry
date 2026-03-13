rm(list=ls())

# load packages
library(ggplot2)
library(dplyr)
library(tidyr)
library(stringr)

# load pan trap data
data <- read.csv("data/cleaned/01_pantraps-cleaned.csv")

# load trait data
bee_traits_sg <- read.csv("data/raw/traits/05_bee-traits-Safeguard.csv")
bee_traits_sr <- read.csv("data/raw/traits/05_bee-traits-Roberts.csv")
syr_traits <- read.csv("data/raw/traits/05_hoverfly-traits-Staton-2022.csv")

###################
# solitary bees
###################

# add lecty column from Roberts database to Safeguard database
extra <- subset(bee_traits_sr, select=c(Species.name.in.database, Lecty..Cane.2020.))
names(extra) <- c("Species", "Lecty")
bee_traits <- left_join(bee_traits_sg, extra)

# subset pantrap data 
sol <- subset(data, TaxGroup=="Anthophila" & Genus != "Bombus", select=c(Year, Month, Site, Treatment, Sample, Spp, Sex, Abundance))

# subset traits data to only species found in pan traps
bee_traits <- subset(bee_traits, Species %in% unique(sol$Spp))

# check all spp present
length(bee_traits$Species)-length(unique(sol$Spp))
setdiff(unique(sol$Spp), bee_traits$Species)

# halictus tumulorum missing from traits data - will add data from other sources later.
# sol <- subset(sol, Spp %in% sol_traits$Species)

# subset traits to candidates for functional analyses and inspect
bee_traits_sub <- subset(bee_traits, 
                         select=c(Species, 
                                  nesting_method_excavator, 
                                  nesting_method_mason, 
                                  nesting_method_renter, 
                                  ITD_F_mm, 
                                  ITD_M_mm, 
                                  TongueLength_F_mm, 
                                  TongueLength_M_mm, 
                                  gens_.1, 
                                  gens_1, 
                                  gens_2, 
                                  gens_.2, 
                                  Flight_period_.months., 
                                  HairLength_F, 
                                  HairLength_M, 
                                  HairDensity_F, 
                                  HairDensity_M, 
                                  HairinessIndex_F,
                                  Lecty)) 

# remove hairiness metrics as lots of missing data
bee_traits_sub <- subset(bee_traits_sub, 
                         select=-c(HairLength_F, 
                                   HairLength_M, 
                                   HairDensity_F, 
                                   HairDensity_M, 
                                   HairinessIndex_F, 
                                   HairinessIndex_F)) 

# assess relationship between ITD and tongue length
ggplot(bee_traits_sub, aes(x=ITD_F_mm, y=TongueLength_F_mm)) +
  geom_point()+
  geom_smooth() +
  theme_bw()

ggplot(bee_traits_sub, aes(x=ITD_M_mm, y=TongueLength_M_mm)) +
  geom_point()+
  geom_smooth() +
  theme_bw()

# clearly colinear. Remove tongue length (diet specialism will be captured by lecty)
bee_traits_sub <- subset(bee_traits_sub, select=-c(TongueLength_F_mm, TongueLength_M_mm))

# missing data: ITD for males of 4 spp: lasioglossum parvulum, lasioglossum puncticolle, lasioglossum semilucens, hylaeus dilatatus. 

# will apply taxonomically-informed imputation (i.e. take these traits from similar species in the family, with the closest value of the female version of that trait)

impute <- function(traits.df, species, species.col, trait.col, comparitor.col){
    # find species in family with closest comparitor trait to the species with the missing trait
  fam <- traits.df[grepl(word(species, 1), traits.df[[species.col]]),]
  fam$diff <- abs(fam[[comparitor.col]]-fam[[comparitor.col]][ fam[[species.col]] == species ])
  closest <- fam[fam$diff>0,species.col][which.min(fam$diff[fam$diff>0])]
  closest.val <- fam[[comparitor.col]][fam[[species.col]]==closest]
  traits.df[[trait.col]][traits.df[[species.col]]==species] <- fam[[trait.col]][fam[[species.col]]==closest]
  
  print(paste("Species observed in family with closest value of", comparitor.col, "to", species, "is",closest, "(", closest.val, "vs", fam[[comparitor.col]][ fam[[species.col]] == species ], "). Value for", trait.col, "now imputated."))

  return(traits.df)
  }

for(i in which(is.na(bee_traits_sub$ITD_M_mm))){
  bee_traits_sub <- impute(traits.df=bee_traits_sub, 
                   species=bee_traits_sub$Species[i], 
                   species.col="Species",
                   trait.col="ITD_M_mm",
                   comparitor.col="ITD_F_mm")
  }

# Fill in missing lecty info from Falk field guide
bee_traits_sub$Lecty[bee_traits_sub$Species=="Nomada fabriciana"] <- "Polylectic"
bee_traits_sub$Lecty[bee_traits_sub$Species=="Nomada flavoguttata"] <- "Polylectic"
bee_traits_sub$Lecty[bee_traits_sub$Species=="Nomada ruficornis"] <- "Polylectic"
bee_traits_sub$Lecty[bee_traits_sub$Species=="Sphecodes gibbus"] <- "Mesolectic"
bee_traits_sub$Lecty[bee_traits_sub$Species=="Sphecodes niger"] <- "Mesolectic"

# convert broadly polylectic to polylectic
bee_traits_sub$Lecty[bee_traits_sub$Lecty=="Broadly polylectic"] <- "Polylectic"

# format traits dataframe
bee_traits_sub$Flight_period_.months. <- str_count(bee_traits_sub$Flight_period_.months., " ") + 1

# Process nesting_method
nesting_method_df <- bee_traits_sub %>%
  dplyr::select(Species, starts_with("nesting_method_")) %>%
  pivot_longer(cols = starts_with("nesting_method_"),
               names_to = "nesting_method",
               names_prefix = "nesting_method_",
               values_to = "present") %>%
  filter(present == "yes") %>%
  group_by(Species) %>%
  summarise(
    nesting_method = if_else(n() > 1, "multiple", first(nesting_method)),
    .groups = "drop"
  )

# Process gens
gens_df <- bee_traits_sub %>%
  dplyr::select(Species, starts_with("gens_")) %>%
  pivot_longer(cols = starts_with("gens_"),
               names_to = "gens",
               names_prefix = "gens_",
               values_to = "present") %>%
  filter(present == "yes") %>%
  mutate(
    gens_num = case_when(
      gens == ".1" ~ 0.5,
      gens == "1" ~ 1,
      gens == "2" ~ 2,
      gens == ".2" ~ 2.5
    )
  ) %>%
  group_by(Species) %>%
  summarise(
    gens = max(gens_num, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    gens = case_when(
      gens == 0.5 ~ "<1",
      gens == 1 ~ "1",
      gens == 2 ~ "2",
      gens == 2.5 ~ ">2"
    )
  )

# Combine back together
bee_traits_clean <- bee_traits_sub %>%
  dplyr::select(-starts_with("nesting_method_"), -starts_with("gens_")) %>%
  left_join(nesting_method_df, by = "Species") %>%
  left_join(gens_df, by = "Species")

# add halictus tumulorum (data from Roberts database and/or Falk ID)
tumulorum_traits <- subset(bee_traits_sr, Species.name.in.database=="Halictus tumulorum")

flight_cols <- which(colnames(tumulorum_traits) %in% c ("Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"))

tumulorum <- data.frame(Species="Halictus tumulorum",
                        ITD_F_mm=tumulorum_traits$ITD.Mean_f,
                        ITD_M_mm=NA, # will impute 
                        Flight_period_.months.=sum(tumulorum_traits[1,flight_cols], na.rm=T),
                        Lecty=tumulorum_traits$Lecty..Cane.2020.,
                        nesting_method=tumulorum_traits$Nesting.trait,
                        gens=tumulorum_traits$Voltinism)

tumulorum <- tumulorum %>%
  mutate(Lecty=dplyr::recode(Lecty, "Broadly polylectic"="Polylectic"),
         nesting_method=dplyr::recode(nesting_method, "Excavator: Ground"="excavator"),
         gens=dplyr::recode(gens, "Univoltine"="1"))

bee_traits_clean <-rbind(bee_traits_clean, tumulorum)

halictus <- subset(bee_traits_sg, Genus=="Halictus", select=c(Species, ITD_M_mm, ITD_F_mm))

bee_traits_clean <- impute(bee_traits_clean, "Halictus tumulorum", "Species", "ITD_M_mm", "ITD_F_mm")

# make lengthwise by sex
bee_traits_long <- bee_traits_clean %>%
  pivot_longer(
    cols = c(ITD_F_mm, ITD_M_mm, ),
    names_to = c(".value", "Sex"),
    names_pattern = "(ITD|TongueLength)_(F|M)_mm"
  )

bee_traits_long$gens <- factor(bee_traits_long$gens, levels=c("<1", "1", "2", ">2"))
bee_traits_long$Flight_period_.months. <- factor(bee_traits_long$Flight_period_.months.)
bee_traits_long$Lecty. <- factor(bee_traits_long$Lecty)

# Combine species and sex
bee_traits_long$Species <- paste(bee_traits_long$Species, bee_traits_long$Sex, sep=".")
sol$Spp <- paste(sol$Spp, sol$Sex, sep=".")

# subset again to spp sex combos found in sol
bee_traits_long <- subset(bee_traits_long, Species %in% unique(sol$Spp))

##############
# hoverflies 
##############

# subset pantrap data 
syr <- subset(data, TaxGroup=="Syrphidae", select=c(Year, Month, Site, Treatment, Sample, Spp, Sex, Abundance))

# replace 'merodontini strigatus with eumerus strigatus (to match traits df)
syr <- syr %>%
  mutate(Spp = dplyr::recode(Spp, "Merodontini strigatus"="Eumerus strigatus"))

# subset traits data to only hoverfly species found in pan traps
syr_traits <- subset(syr_traits, Spp %in% unique(syr$Spp))

# Remove overwintering phase
syr_traits <- subset(syr_traits, select=-c(OverwinteringPhase))

# assess relationship between voltinism & duration of development
plot(Voltinism ~ DurationOfDevelopment, data=syr_traits) # very colinear. remove duration of development
syr_traits <- subset(syr_traits, select=-c(DurationOfDevelopment))

###########
# output
###########
write.csv(bee_traits_long, "data/cleaned/05_solitary-bee-traits-cleaned.csv", row.names=F)
write.csv(syr_traits, "data/cleaned/05_hoverfly-traits-cleaned.csv", row.names=F)

