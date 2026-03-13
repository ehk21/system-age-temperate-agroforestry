rm(list=ls())

# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(iNEXT)
library(vegan)

# load data
flora <- read.csv("data/cleaned/02_plants-cleaned.csv")

# format
flora$Site <- factor(flora$Site)

# summarise sample sizes
flora.summ <- flora %>%
  group_by(Year, Site, Treatment, Position, Method) %>%
  summarise(n=length(unique(ID)))

# remove irrelevant species (non-insect pollinated or non-flowering)
to_exclude <- c(
  "potato",                         # Tuber, not insect-pollinated
  "meadow fescue",                  # Grass, wind-pollinated
  "red fescue",                     # Grass
  "rough meadow-grass",            # Grass
  "marsh foxtail",                 # Grass
  "grass sp.",                      # Ambiguous
  "black-grass",                    # Grass
  "wheat",                          # Self-/wind-pollinated crop
  "timothy grass",                  # Grass
  "oat (unspecified)",              # Self-/wind-pollinated crop
  "barren brome",                   # Grass
  "rough meadow grass",            # Duplicate of above, different format
  "soft brome",                     # Grass
  "wall barley",                    # Grass
  "false oat-grass",                # Grass
  "meadow barley",                  # Grass
  "meadow foxtail",                 # Grass
  "cock's-foot",                    # Grass
  "perennial ryegrass",            # Grass
  "crested dog’s-tail",            # Grass
  "ribwort plantain",              # Wind-pollinated
  "clustered dock",                # Wind-pollinated
  "broad-leaved dock",             # Wind-pollinated
  "curled dock",                   # Wind-pollinated
  "dock",                           # Ambiguous
  "knotgrass",                      # Wind/self-pollinated
  "moss sp.",                       # Non-vascular
  "spurge (unspecified)",           # Ambiguous; many Euphorbias not insect-pollinated
  "mint sp.",                       # Too vague
  "goosefoot sp.",                  # Wind/self-pollinated
  "clover (unspecified)",           # Too vague - red and white clover both found in same survey year so not needed
  "dead-nettle (unspecified)",      # Too vague
  "dicot sp.",                      # Ambiguous
  "black medick",                   # Can be self-pollinated
  "cleavers",                       # Galium aparine, self/wind-pollinated
  "goosegrass",                     # Synonym for cleavers
  "common nettle",                 # Wind-pollinated
  "hemp-nettle species"           # Ambiguous
)

flora <- flora[!flora$Common.name %in% to_exclude, ]

# combine 2023 data for each quadrat over all months
flora <- flora %>%
  group_by(Year, Site, Treatment, ID, Position, Taxon, Common.name, Method) %>%
  summarise(Percentage=sum(Percentage))

# remove zero values 
flora <- flora[flora$Percentage>0,]

# convert to presence-absence matrix format and add samples where no species were recorded 
flora <- flora %>%
  mutate(presence = 1) 

flora.wide <- subset(flora, select = c(-Common.name, -Percentage)) %>%
  pivot_wider(
    names_from = Taxon,
    values_from = presence,
    values_fill = list(presence = 0)
  )

samples.all <- read.csv("data/cleaned/11_plant-sample-names.csv")
samples <- unique(flora.wide$ID)
samples.other <- subset(samples.all, !(ID %in% samples))

flora.wide <- plyr::rbind.fill(flora.wide, samples.other)
flora.wide[is.na(flora.wide)] <- 0

#########################################
# analyse CROP data only (no treerow spp)
#########################################

# Filter to only crop samples
crop_flora <- flora.wide %>%
  filter(Position == "c") %>%
  dplyr::select(-Position, -Method)

meta_cols <- c("Year", "Site", "Treatment", "ID")
stopifnot(all(meta_cols %in% names(crop_flora)))

species_cols <- setdiff(names(crop_flora), meta_cols)

# ---- 1) Clean species columns to strict 0/1 ----
dat <- crop_flora %>%
  mutate(across(all_of(species_cols), ~ as.integer(.))) %>%
  mutate(across(all_of(species_cols), ~ ifelse(is.na(.), 0L, .))) %>%
  mutate(across(all_of(species_cols), ~ ifelse(. > 0, 1L, 0L)))

# ---- 2) Ensure one record per Year-Site-Treatment-ID ----
dat_id <- dat %>%
  group_by(Year, Site, Treatment, ID) %>%
  summarise(across(all_of(species_cols), ~ max(.x, na.rm = TRUE)), .groups = "drop")

# ---- 3) Build incidence-frequency vectors per assemblage ----
# For each Year-Site-Treatment:
#   T = number of unique IDs (sampling units)
#   freq(species) = number of IDs where species is present
inc_list <- dat_id %>%
  mutate(assemblage = paste(Year, Site, Treatment, sep = "_")) %>%
  group_by(assemblage) %>%
  group_map(~{
    spp <- .x %>% dplyr::select(all_of(species_cols))
    T <- as.integer(nrow(spp))
    freq <- as.integer(colSums(spp))  # incidence frequency per species
    c(T, freq)
  })

# Name list elements
assemblage_names <- dat_id %>%
  mutate(assemblage = paste(Year, Site, Treatment, sep = "_")) %>%
  distinct(assemblage) %>%
  pull(assemblage)

names(inc_list) <- assemblage_names

# remove assemblages where no plants were observed
rm <- c()
for(i in 1:length(inc_list)){
  assemb <- inc_list[[i]]
  seen <- sum(assemb[2:length(assemb)])
  print(seen)
  if(seen==0){rm<-c(rm,i)} 
}

inc_list_ok <- inc_list[-rm]

# ---- 4) Run iNEXT (Simpson diversity = q=2; incidence_freq datatype) ----
out <- iNEXT(
  x = inc_list_ok,
  q = 2,
  datatype = "incidence_freq",
  nboot = 50  
)

DataInfo<- out[["DataInfo"]]

simpson_rarefied <- estimateD(
  inc_list_ok, 
  q = 2,
  datatype = "incidence_freq",
  base = "coverage",
  nboot = 50,
  conf = 0.95
)

# coverage- rather than size-based rarefaction because sampling effort was very uneven (in terms of quadrat size, number of samples, etc.)
# so robust to varying quadrat numbers, size and survey frequency
# Simpson diversity (Q=2) was chosen because extrapolation is much more unbiased compared to Species Richness (Chao et al., 2014)
# As level=NULL, the function computes the diversity estimates for the minimum coverage value for samples extrapolated to double the size of the reference sample (Chao et al., 2014 (or 2012?))


simpson_rarefied <- simpson_rarefied %>%
  separate(Assemblage, into = c("Year", "Site", "Treatment"), sep = "\\_")

# add age column
simpson_rarefied$Year <- as.numeric(simpson_rarefied$Year)
simpson_rarefied$Age<-NA
simpson_rarefied$Age[simpson_rarefied$Site==3] <- simpson_rarefied$Year[simpson_rarefied$Site==3]-2009
simpson_rarefied$Age[simpson_rarefied$Site==2] <- simpson_rarefied$Year[simpson_rarefied$Site==2]-2015
simpson_rarefied$Age[simpson_rarefied$Site==1] <- simpson_rarefied$Year[simpson_rarefied$Site==1]-2014

# convert 95% confidence interval to standard error (assuming normal distribution: 95% confidence interval is the mean +- 1.96*SE)
simpson_rarefied <- simpson_rarefied %>%
  mutate(SE = (qD.UCL - qD.LCL) / (2 * 1.96))

# plot
pd <- position_dodge(width = 0.3)

ggplot(simpson_rarefied, aes(x=Age, y=qD, color=Treatment)) +
  geom_point(aes(shape=Site, 
                 fill=Treatment,
                 group=Treatment), 
             size = 3,
             alpha = 0.75,
             stroke = 0.9,
             colour = "black",
             position = pd) +
  geom_smooth(aes(linetype = Treatment, 
                  group = Treatment),
              method = "lm", 
              se = FALSE,
              linewidth = 0.8,
              colour = "black") +
  geom_errorbar(aes(ymin=qD-SE, 
                    ymax=qD+SE, 
                    group=Treatment), 
                width = .8,
                linewidth = 0.45,
                colour = "grey25",
                alpha = 0.9,
                position = pd) +
  theme_bw() +
  scale_shape_manual(values = c(22, 21, 24)) +   
  scale_fill_manual(values=c("grey20","white" )) +
  theme_bw() + 
  ylab("Hill-Chao diversity ('Simpson')") + 
  xlab("System age (years)") + 
  # ggtitle("Simpson Diversity of plants in the cropped area (not counting treerows)") +
  guides(
    fill = guide_legend(
      title = "Treatment",
      override.aes = list(
        shape = 21,
        colour = "black",
        size = 4,
        alpha = 1
      ))) +
  theme(legend.position = "bottom", 
        legend.title = element_text(size = 13), 
        legend.text = element_text(size = 11), 
        legend.box = "horizontal", 
        panel.grid.major = element_line(colour = "grey88", linewidth = 0.25), 
        panel.grid.minor = element_blank(), axis.title = element_text(size = 15), 
        axis.text = element_text(size = 12, colour = "black"), 
        strip.background = element_rect(fill="lightgrey", color="black"), 
        strip.placement = "outside", 
        strip.text = element_text(size = 13), 
        plot.title = element_text(size = 15, hjust = 0), 
        plot.margin = unit(c(.5, .5, .5, .5), "cm")) +
  scale_x_continuous(breaks=seq(min(simpson_rarefied$Age), max(simpson_rarefied$Age), by=2))


ggsave(filename="figs/11_plants-alpha-div.pdf", width=20, height=15, units="cm", dpi=300)

# ##############################
# # Create TREEROW species lists
# ##############################
# 
# # Filter to only crop samples
# treerow_flora <- flora %>%
#   filter(Position == "t")
# 
# treerows <- treerow_flora %>%
#   group_by(Year, Site, Treatment, Taxon, Common.name) %>%
#   summarise(Eh=sum(Percentage))
# 
# treerow_counts <- treerows %>%
#   group_by(Year, Site) %>%
#   summarise(N=length(unique((Taxon))))
# 
# wh_treerows <- subset(treerows, Site=="wh") %>%
#   mutate(Site_Year = paste(Site, Year, sep = "_")) %>%
#   mutate(presence = ifelse(!is.na(Eh), 1, 1)) 
# 
# wh_table <- subset(wh_treerows, select=c(Taxon, Common.name, Site_Year, presence)) %>%  # assume all rows represent presence
#   pivot_wider(names_from = Site_Year,
#               values_from = presence,
#               values_fill = list(presence = 0))  # fill absent values with 0
# 
# fe_treerows <- subset(treerows, Site=="fe") %>%
#   mutate(Site_Year = paste(Site, Year, sep = "_")) %>%
#   mutate(presence = ifelse(!is.na(Eh), 1, 1)) 
# 
# fe_table <- subset(fe_treerows, select=c(Taxon, Common.name, Site_Year, presence)) %>%  # assume all rows represent presence
#   pivot_wider(names_from = Site_Year,
#               values_from = presence,
#               values_fill = list(presence = 0))  # fill absent values with 0
# 
# rh_treerows <- subset(treerows, Site=="rh") %>%
#   mutate(Site_Year = paste(Site, Year, sep = "_")) %>%
#   mutate(presence = ifelse(!is.na(Eh), 1, 1)) 
# 
# rh_table <- subset(rh_treerows, select=c(Taxon, Common.name, Site_Year, presence)) %>%  # assume all rows represent presence
#   pivot_wider(names_from = Site_Year,
#               values_from = presence,
#               values_fill = list(presence = 0))  # fill absent values with 0



