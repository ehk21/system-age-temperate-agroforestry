rm(list=ls())

# load packages
library(dplyr)
library(vegan)
library(tidyr)
library(ggplot2)
library(plotrix)
library(iNEXT)
library(gridExtra)
library(DHARMa)

# load data
data <- read.csv("data/cleaned/01_pantraps-cleaned.csv")
syr_traits <- read.csv("data/cleaned/05_hoverfly-traits-cleaned.csv")
model.vars <- read.csv("data/cleaned/04_model-vars-seasonal.csv")

# Helper function for hoverfly predictions
get_ma_predictions_syr <- function(model, raw_data, age_mean, age_sd, order_q_label) {
  
  pred_grid <- expand.grid(
    Age_raw = seq(min(raw_data$Age.raw, na.rm=T), max(raw_data$Age.raw, na.rm=T), length.out = 50),
    Treatment = levels(raw_data$Treatment)
  )
  
  # Re-standardise Age using same arm::rescale transformation
  pred_grid$Age <- (pred_grid$Age_raw - age_mean) / (2 * age_sd)
  
  # Hold PC1.climate at mean (= 0 after standardisation)
  pred_grid$PC1.climate <- 0
  
  # Predict on log scale then back-transform
  pred_grid$log_qD <- predict(model, newdata = pred_grid)
  pred_grid$qD_pred <- exp(pred_grid$log_qD)
  pred_grid$Order.q <- order_q_label
  
  return(pred_grid)
}

#############
# format data:
#############

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

# remove species with aquatic larval stages
aquatic_species <- c(
  "Eristalis tenax",
  "Eristalis arbustorum",
  "Eristalis abusivus",
  "Eristalinus sepulchralis",
  "Myathropa florea",
  "Helophilus pendulus",
  "Helophilus hybridus",
  "Chrysogaster solstitialis",
  "Lejogaster metallina",
  "Neoascia tenur",
  "Neoascia podagrica",
  "Tropidia scita"
  )

syrphidae <- subset(syrphidae, !(Spp %in% aquatic_species))

# sum abundance for all sexes/castes and all samples per year, season, site and treatment
syrphidae <- syrphidae %>%
  group_by(Year, Season, Site, Treatment, Spp) %>%
  summarise(Abundance=sum(Abundance)) 

syrphidae <- subset(syrphidae, Season=="Late")

# pivot the dataframe to wide format
syrphidae.wide <- syrphidae %>%
  pivot_wider(names_from = Spp, values_from = Abundance, values_fill = 0)

#############
# calculate diversity:
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
                               "0"="Richness (q=0)",
                               "1"="Shannon (q=1)",
                               "2"="Simpson (q=2)"),
         Order.q = factor(Order.q,
                          levels=c("Richness (q=0)", 
                                   "Shannon (q=1)", 
                                   "Simpson (q=2)" )))


#############
# model and predict diversity:
#############

# combine diversity data with model vars
model.vars$Site <- factor(model.vars$Site)
syr.div <- left_join(syr.baseline, model.vars, by=c("Site", "Year", "Season", "Treatment"))

# factorise categorical variables
cols <- c("Site", "Treatment", "crop.stage", "Season", "Order.q")
syr.div[cols] <- lapply(syr.div[cols], factor)  ## as.factor() could also be used

# store raw age values (for future predictions)
syr.div$Age.raw <- syr.div$Age

# subset to each metric
syr.q0 <- subset(syr.div, Order.q=="Richness (q=0)")
syr.q1 <- subset(syr.div, Order.q=="Shannon (q=1)")
syr.q2 <- subset(syr.div, Order.q=="Simpson (q=2)")

# Store mean and SD of raw Age for re-standardisation
age.mean.raw.syr <- mean(syr.q0$Age, na.rm=T)
age.sd.raw.syr <- sd(syr.q0$Age, na.rm=T)

#############################################################
# inspect data
#############################################################
hist(syr.q0$qD)
hist(syr.q1$qD)
hist(syr.q2$qD)

shapiro.test(log(syr.q0$qD)) # normal when logged
shapiro.test(log(syr.q1$qD)) # normal when logged
shapiro.test(log(syr.q2$qD)) # normal when logged

#############################################################
# reduce climate variables using PCA to account for as many as poss whilst reducing the number of model terms
#############################################################
predictors.climate <- syr.q0 %>%
  dplyr::select(Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly)

# ggpairs(predictors.climate) +
#   theme_bw()

# Run PCA
pca.climate <- prcomp(predictors.climate, center = TRUE, scale. = TRUE)

# View PCA summary
summary(pca.climate) # axis 1 represents the majority of variance (58.7%). Axis 2 only represents 28.2% of variance.
pca.climate$rotation

# Extract first 2 PCs
pcs.climate <- as.data.frame(pca.climate$x[, 1:2])
names(pcs.climate) <- c("PC1.climate", "PC2.climate")

syr.q0 <- bind_cols(syr.q0, pcs.climate)
syr.q1 <- bind_cols(syr.q1, pcs.climate)
syr.q2 <- bind_cols(syr.q2, pcs.climate)

################################################################################
# Species richness (q=0)
################################################################################

# standardise all continuous predictors using arm::rescale
syr.q0 <- syr.q0 %>%
  mutate(across(c(Age, margin, boundary.area, Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly, PC1.climate, PC2.climate), ~ arm::rescale(.x)))

# assign weightings (accounting for log transformation)
syr.q0 <- syr.q0 %>%
  mutate(
    log_qD = log(qD),
    SE_log = SE / qD, # propagate according to big red book look up table
    w = 1 /(SE_log)
  ) %>%
  dplyr::filter(!is.infinite(w))

# specify global model
global_mod_q0 <- lm(log_qD ~ Age * Treatment + PC1.climate,
                    data = syr.q0,
                    weights = w,
                    na.action = "na.fail"
)

# Check the model assumptions
par(mfrow=c(2,2))
plot(global_mod_q0)
res <- simulateResiduals(fittedModel=global_mod_q0, plot=T)

# vif(global_mod_q0)
summary(global_mod_q0)

# Generate predictions
preds_q0_syr <- get_ma_predictions_syr(model=global_mod_q0, 
                                       raw_data=syr.q0, 
                                       age_mean=age.mean.raw.syr, 
                                       age_sd=age.sd.raw.syr, 
                                       order_q_label="Species Richness (q=0)")

################################################################################
# Shannon (q=1)
################################################################################

# standardise all continuous predictors using arm::rescale
syr.q1 <- syr.q1 %>%
  mutate(across(c(Age, margin, boundary.area, Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly, PC1.climate, PC2.climate), ~ arm::rescale(.x)))

# assign weightings (accounting for log transformation)
syr.q1 <- syr.q1 %>%
  mutate(
    log_qD = log(qD),
    SE_log = SE / qD, # propagate according to big red book look up table
    w = 1 / (SE_log)
  ) %>%
  dplyr::filter(!is.infinite(w))

# specify global model
global_mod_q1 <- lm(log_qD ~ Age * Treatment + PC1.climate,
                    data = syr.q1,
                    weights = w,
                    na.action = "na.fail"
)

# Check the model assumptions
par(mfrow=c(2,2))
plot(global_mod_q1)
res <- simulateResiduals(fittedModel=global_mod_q1, plot=T)

vif(global_mod_q1)
summary(global_mod_q1)

# Generate predictions
preds_q1_syr <- get_ma_predictions_syr(global_mod_q1, syr.q1, age.mean.raw.syr, age.sd.raw.syr, "Shannon Diversity (q=1)")

################################################################################
# Simpson(q=2)
################################################################################

# standardise all continuous predictors using arm::rescale
syr.q2 <- syr.q2 %>%
  mutate(across(c(Age, margin, boundary.area, Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly, PC1.climate, PC2.climate), ~ arm::rescale(.x)))

# assign weightings (accounting for log transformation)
syr.q2 <- syr.q2 %>%
  mutate(
    log_qD = log(qD),
    SE_log = SE / qD, # propagate according to big red book look up table
    w = 1 / (SE_log)
  ) %>%
  dplyr::filter(!is.infinite(w))

# specify global model
global_mod_q2 <- lm(log_qD ~ Age * Treatment + PC1.climate,
                    data = syr.q2,
                    weights = w,
                    na.action = "na.fail"
)

# Check the model assumptions
par(mfrow=c(2,2))
plot(global_mod_q2)
res <- simulateResiduals(fittedModel=global_mod_q2, plot=T)

vif(global_mod_q2)
summary(global_mod_q2)

# Generate predictions
preds_q2_syr <- get_ma_predictions_syr(global_mod_q2, syr.q2, age.mean.raw.syr, age.sd.raw.syr, "Simpson Diversity (q=2)")

#####################
# Combine and export predictions
#####################
pred <- bind_rows(preds_q0_syr, preds_q1_syr, preds_q2_syr) %>%
  mutate(Order.q=dplyr::recode(Order.q, 
                               "Species Richness (q=0)" = "q=0 (Richness)",
                               "Shannon Diversity (q=1)"="q=1 ('Shannon')", 
                               "Simpson Diversity (q=2)" = "q=2 ('Simpson')"))

#############
# plot
#############

syr.baseline <- syr.baseline %>%
  mutate(Order.q=dplyr::recode(Order.q, 
                               "Richness (q=0)" = "q=0 (Richness)",
                               "Shannon (q=1)"="q=1 ('Shannon')", 
                               "Simpson (q=2)" = "q=2 ('Simpson')")) 

pred$Order.q <- factor(pred$Order.q, levels = levels(syr.baseline$Order.q))

syr.baseline$Site <- factor(syr.baseline$Site)

pd <- position_dodge(width = 0.3)

syr.alpha.plot <-   
  ggplot(syr.baseline, aes(x = Age, y = qD, color=Treatment)) +
  geom_point(
    aes(shape = Site, 
        fill = Treatment, 
        group = Treatment),
    size = 3,
    alpha = 0.75,
    stroke = 0.9,
    colour = "black",
    position = pd
  ) +
  # error bars: keep them neutral + slightly lighter
  geom_errorbar(
    aes(ymin = qD - SE, ymax = qD + SE, group = Treatment),
    width = .8,
    linewidth = 0.45,
    colour = "grey25",
    alpha = 0.9,
    position = pd
  ) +
  
  # geom_smooth(
  #   aes(linetype = Treatment, group = Treatment),
  #   method = "lm", se = FALSE,
  #   linewidth = 0.8,
  #   colour = "black"
  # ) +
  geom_line(data=pred, aes(x=Age_raw, y=qD_pred, linetype=Treatment, group=Treatment),
            linewidth=.8,
            colour="black") +
  facet_grid(~ factor(Order.q,
                      levels = c("q=0 (Richness)", "q=1 ('Shannon')", "q=2 ('Simpson')"))) +
  scale_shape_manual(values = c(22, 21, 24)) +   
  scale_fill_manual(values=c("grey20","white" )) +
  theme_bw() + 
  ylab("Hoverfly alpha diversity (Hill–Chao)") + 
  xlab("System age (years)") + 
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
  scale_x_continuous(breaks=seq(min(syr.baseline$Age), max(syr.baseline$Age), by=2))


syr.alpha.plot

ggsave(plot=syr.alpha.plot, filename="figs/21_sensitivity-hoverflies-alpha-div-nonaquatic.pdf", width=20, height=10, units="cm", dpi=300)







