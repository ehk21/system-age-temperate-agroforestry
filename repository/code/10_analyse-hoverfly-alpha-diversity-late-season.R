rm(list=ls())

# load packages
library(arm) # standardise continuous independent variables to allow comparisons of the relative strength of parameter estimates.
library(MuMIn) # for multi-model inference
library(dplyr)
library(tidyr)
library(car) # for vif function 
library(DHARMa)
library(ggplot2)
library(GGally) # for correlation ggpairs plot
library(diagis) # for weighted SE
library(gridExtra)
library(fitdistrplus)

# load data
syr.div <- read.csv("data/cleaned/05_syrphidae-alpha-div.csv")
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

# combine diversity data with model vars
syr.div <- left_join(syr.div, model.vars)

# plot
ggplot(syr.div, aes(x=Age, y=qD, color=Treatment)) +
  geom_smooth(method="lm", se=F) +
  # geom_line(alpha=0.5,aes(color=factor(Site))) +
  geom_point(aes(shape=factor(Site)), alpha=.5, size=3) +
  geom_errorbar(aes(ymin=qD-SE, ymax=qD+SE), alpha=.5, width=.3) +
  theme_bw() +
  ylab("Diversity") +
  xlab("Age (Years)") +
  theme(axis.title.y=element_blank(),
        plot.margin=unit(c(.5,0,.5,.5), "cm")) +
  ggtitle("Hoverfly Alpha Diversity") +
  # facet_wrap(Order.q~Season, ncol=2)
  facet_wrap(Order.q~Season, ncol=2)


# factorise categorical variables
cols <- c("Site", "Treatment", "crop.stage", "Season", "Order.q")
syr.div[cols] <- lapply(syr.div[cols], factor)  ## as.factor() could also be used

# store raw age values (for future predictions)
syr.div$Age.raw <- syr.div$Age

# subset to each metric
syr.q0 <- subset(syr.div, Order.q=="Species Richness (q=0)")
syr.q1 <- subset(syr.div, Order.q=="Shannon Diversity (q=1)")
syr.q2 <- subset(syr.div, Order.q=="Simpson Diversity (q=2)")

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

par(mfrow=c(1,1))
descdist(syr.q0$qD)
descdist(syr.q1$qD)
descdist(syr.q2$qD)

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
  )

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

vif(global_mod_q0)
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
  )

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
  )

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
preds_all_syr <- bind_rows(preds_q0_syr, preds_q1_syr, preds_q2_syr)
preds_all_syr$Order.q <- factor(preds_all_syr$Order.q, levels = levels(syr.div$Order.q))

ggplot(preds_all_syr, aes(x=Age_raw, y=qD_pred, linetype=Treatment)) +
  geom_line() +
  facet_wrap(~Order.q) +
  theme_bw()

write.csv(preds_all_syr, "results/10_hoverfly-alpha-predictions.csv", row.names=F)


