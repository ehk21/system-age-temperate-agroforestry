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
library(equatiomatic)

# source functions
source("code/functions/08_get-predictions-function.R")
source("code/functions/08_extract-model-avg-results-function.R")

# load data
sol.div <- read.csv("data/cleaned/05_solitary-alpha-div.csv")
model.vars <- read.csv("data/cleaned/04_model-vars-seasonal.csv")

# combine diversity data with model vars
sol.div <- left_join(sol.div, model.vars)

# plot
ggplot(sol.div, aes(x=Age, y=qD, color=Treatment)) +
  geom_smooth(method="lm", se=F) +
  # geom_line(alpha=0.5,aes(color=factor(Site))) +
  geom_point(aes(shape=factor(Site)), alpha=.5, size=3) +
  geom_errorbar(aes(ymin=qD-SE, ymax=qD+SE), alpha=.5, width=.3) +
  theme_bw() +
  ylab("Diversity") +
  xlab("Age (Years)") +
  theme(axis.title.y=element_blank(),
        plot.margin=unit(c(.5,0,.5,.5), "cm")) +
  ggtitle("Solitary Bee Alpha Diversity") +
  # facet_wrap(Order.q~Season, ncol=2)
  facet_wrap(Order.q~Season, ncol=2)

# factorise categorical variables
cols <- c("Site", "Treatment", "crop.stage", "Season", "Order.q")
sol.div[cols] <- lapply(sol.div[cols], factor)  ## as.factor() could also be used

# # sensitivity: remove late season site 3 outlier
# sol.div <- subset(sol.div, !(Year==2018 & Season=="Late" & Site==3))

# store raw age values (for future predictions)
sol.div$Age.raw <- sol.div$Age

# subset to each metric
sol.q0 <- subset(sol.div, Order.q=="Species Richness (q=0)")
sol.q1 <- subset(sol.div, Order.q=="Shannon Diversity (q=1)")
sol.q2 <- subset(sol.div, Order.q=="Simpson Diversity (q=2)")

#############################################################
# inspect data 
#############################################################
# hist(sol.q0$qD)
# hist(sol.q1$qD)
# hist(sol.q2$qD)
# 
# shapiro.test(log(sol.q0$qD)) # normal when logged
# shapiro.test(log(sol.q1$qD)) # normal when logged
# shapiro.test(log(sol.q2$qD)) # normal when logged
# 
# par(mfrow=c(1,1))
# descdist(sol.q0$qD)
# descdist(sol.q1$qD)
# descdist(sol.q2$qD)
# #lognormal for all metrics (more robust to weightings than Gamma)

#############################################################
# reduce climate variables using PCA to account for as many as poss whilst reducing the number of model terms
#############################################################
predictors.climate <- sol.q0 %>%
  dplyr::select(Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly)

# ggpairs(predictors.climate) +
#   theme_bw()

# Run PCA
pca.climate <- prcomp(predictors.climate, center = TRUE, scale. = TRUE)

# # View PCA summary
# summary(pca.climate) # axis 1 represents the majority of variance (59.1%). Axis 2 only represents 23.2% of variance.
# pca.climate$rotation

# Extract first 2 PCs
pcs.climate <- as.data.frame(pca.climate$x[, 1:2])
names(pcs.climate) <- c("PC1.climate", "PC2.climate")

sol.q0 <- bind_cols(sol.q0, pcs.climate)
sol.q1 <- bind_cols(sol.q1, pcs.climate)
sol.q2 <- bind_cols(sol.q2, pcs.climate)

################################################################################
# Species richness (q=0)
################################################################################

# standardise all continuous predictors using arm::rescale
sol.q0 <- sol.q0 %>%
  mutate(across(c(Age, margin, boundary.area, Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly, PC1.climate, PC2.climate), ~ arm::rescale(.x)))

# assign weightings (accounting for log transformation)
sol.q0 <- sol.q0 %>%
  mutate(
    log_qD = log(qD),
    SE_log = SE / qD, # propagate according to big red book look up table
    w = 1 /(SE_log)
  )

# specify global model
global_mod_q0 <- lm(log_qD ~ Age + Treatment + Age:Treatment + Site + Season + Age:Season + crop.stage + PC1.climate,
                    data = sol.q0,
                    weights = w,
                    na.action = "na.fail"
                    )

vif(global_mod_q0)
summary(global_mod_q0)

# Check the model assumptions
par(mfrow=c(2,2))
plot(global_mod_q0)
res <- simulateResiduals(fittedModel=global_mod_q0, plot=T)
# accept

# # extract model equation for paper
# eq <- equatiomatic::extract_eq(global_mod_q0)

# create model sets
model_set_q0 <- dredge(global_mod_q0, rank = "AICc")
model_set_fixed_q0 <- dredge(global_mod_q0, rank = "AICc", fixed=c("Age", "Treatment", "Age:Treatment"))

# subset to models whose delta is 4 or less
model_subset_q0 <- subset(model_set_q0, delta <= 4)
model_subset_fixed_q0 <- subset(model_set_fixed_q0, delta <= 4)

# subset again to make sure the number of models is less than the number of site-year combinations
if(nrow(model_subset_q0)>=10){
  model_subset_q0 <- model_subset_q0[1:10,]
}
if(nrow(model_subset_fixed_q0)>=10){
  model_subset_fixed_q0 <- model_subset_fixed_q0[1:10,]
}

# Average top models with ΔAICc < 4
model_avg_q0 <- model.avg(model_subset_q0, fit=T)
model_avg_fixed_q0 <- model.avg(model_subset_fixed_q0, fit=T)

# look at summaries
summary(model_avg_q0)
summary(model_avg_fixed_q0)

# check confidence intervals
confint(model_avg_q0, full = F)
confint(model_avg_fixed_q0, full = F)

# check variable weights
sw(model_subset_q0)
sw(model_subset_fixed_q0)

# generate predictions
preds.q0<- get_ma_predictions(model_avg_fixed_q0, 
                              raw_data=sol.q0, 
                              age_mean=mean(sol.q0$Age.raw, na.rm=T),
                              age_sd=sd(sol.q0$Age.raw, na.rm=T),
                              order_q_label = "Species Richness (q=0)")

ggplot(preds.q0, aes(x=Age_raw, y=qD_pred, color=Treatment)) +
  geom_line() +
  facet_wrap(~Season) +
  theme_bw()

################################################################################
# Shannon (q=1)
################################################################################

# standardise all continuous predictors using arm::rescale
sol.q1 <- sol.q1 %>%
  mutate(across(c(Age, margin, boundary.area, Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly, PC1.climate, PC2.climate), ~ arm::rescale(.x)))

# assign weightings (accounting for log transformation)
sol.q1 <- sol.q1 %>%
  mutate(
    log_qD = log(qD),
    SE_log = SE / qD, # propagate according to big red book look up table
    w = 1 / (SE_log)
  )

# specify global model
global_mod_q1 <- lm(log_qD ~ Age + Treatment + Age:Treatment + Site + Season + Age:Season + crop.stage + PC1.climate,
                    data = sol.q1,
                    weights = w,
                    na.action = "na.fail"
)

vif(global_mod_q1)
summary(global_mod_q1)

# Check the model assumptions
res <- simulateResiduals(fittedModel=global_mod_q1, plot=T)

# create model sets
model_set_q1 <- dredge(global_mod_q1, rank = "AICc")
model_set_fixed_q1 <- dredge(global_mod_q1, rank = "AICc", fixed=c("Age", "Treatment", "Age:Treatment"))

# subset to models whose delta is 4 of less
model_subset_q1 <- subset(model_set_q1, delta <= 4)
model_subset_fixed_q1 <- subset(model_set_fixed_q1, delta <= 4)

# subset again to make sure the number of models is less than the number of site-year combinations
if(nrow(model_subset_q1)>=10){
  model_subset_q1 <- model_subset_q1[1:10,]
}
if(nrow(model_subset_fixed_q1)>=10){
  model_subset_fixed_q1 <- model_subset_fixed_q1[1:10,]
}

# Average top models with ΔAICc < 4
model_avg_q1 <- model.avg(model_subset_q1, fit=T)
model_avg_fixed_q1 <- model.avg(model_subset_fixed_q1, fit=T)

# look at summaries
summary(model_avg_q1)
summary(model_avg_fixed_q1)

# check confidence intervals
confint(model_avg_q1, full = F)
confint(model_avg_fixed_q1, full = F)

# check variable weights
sw(model_subset_q1)
sw(model_subset_fixed_q1)

# generate predictions
preds.q1<- get_ma_predictions(model_avg_fixed_q1, 
                              raw_data=sol.q1, 
                              age_mean=mean(sol.q1$Age.raw, na.rm=T),
                              age_sd=sd(sol.q1$Age.raw, na.rm=T),
                              order_q_label = "Shannon Diversity (q=1)")

ggplot(preds.q1, aes(x=Age_raw, y=qD_pred, color=Treatment)) +
  geom_line() +
  facet_wrap(~Season) +
  theme_bw()

################################################################################
# Simpson(q=2)
################################################################################

# standardise all continuous predictors using arm::rescale
sol.q2 <- sol.q2 %>%
  mutate(across(c(Age, margin, boundary.area, Temp_anomoly, Temp_prev_anomoly, Precip_anomoly, Precip_prev_anomoly, PC1.climate, PC2.climate), ~ arm::rescale(.x)))

# assign weightings (accounting for log transformation)
sol.q2 <- sol.q2 %>%
  mutate(
    log_qD = log(qD),
    SE_log = SE / qD, # propagate according to big red book look up table
    w = 1 / (SE_log)
  )

# specify global model
global_mod_q2 <- lm(log_qD ~ Age + Treatment + Age:Treatment + Site + Season + Age:Season + crop.stage + PC1.climate,
                    data = sol.q2,
                    weights = w,
                    na.action = "na.fail"
                    )

vif(global_mod_q2)
summary(global_mod_q2)

# Check the model assumptions
par(mfrow=c(2,2))
plot(global_mod_q2)
res <- simulateResiduals(fittedModel=global_mod_q2, plot=T)

# create model sets
model_set_q2 <- dredge(global_mod_q2, rank = "AICc")
model_set_fixed_q2 <- dredge(global_mod_q2, rank = "AICc", fixed=c("Age", "Treatment", "Age:Treatment"))

# subset to models whose delta is 4 of less
model_subset_q2 <- subset(model_set_q2, delta <= 4)
model_subset_fixed_q2 <- subset(model_set_fixed_q2, delta <= 4)

# subset again to make sure the number of models is less than the number of site-year combinations
if(nrow(model_subset_q2)>=10){
  model_subset_q2 <- model_subset_q2[1:10,]
}
if(nrow(model_subset_fixed_q2)>=10){
  model_subset_fixed_q2 <- model_subset_fixed_q2[1:10,]
}

# Average top models with ΔAICc < 4
model_avg_q2 <- model.avg(model_subset_q2, fit=T)
model_avg_fixed_q2 <- model.avg(model_subset_fixed_q2, fit=T)

# look at summaries
summary(model_avg_q2)
summary(model_avg_fixed_q2)

# check confidence intervals
confint(model_avg_q2, full = F)
confint(model_avg_fixed_q2, full = F)

# check variable weights
sw(model_subset_q2)
sw(model_subset_fixed_q2)

# generate predictions
preds.q2<- get_ma_predictions(model_avg_q2, 
                              raw_data=sol.q2, 
                              age_mean=mean(sol.q2$Age.raw, na.rm=T),
                              age_sd=sd(sol.q2$Age.raw, na.rm=T),
                              order_q_label = "Simpson Diversity (q=2)")

ggplot(preds.q2, aes(x=Age_raw, y=qD_pred, color=Treatment)) +
  geom_line() +
  facet_wrap(~Season) +
  theme_bw()

################################################################################
# combine predictions
################################################################################
preds.all <- bind_rows(preds.q0, preds.q1, preds.q2)
preds.all$Order.q <- factor(preds.all$Order.q, levels = levels(sol.div$Order.q))

################################################################################
# extract model results
################################################################################

# how variable names should appear in output tables and figures
pretty <- c(
  "Age" = "Age (yrs)",
  "TreatmentControl" = "Treatment: Control (vs ref)",
  "Age:TreatmentControl" = "Age × Treatment: Control (vs ref)",
  "SeasonLate" = "Season: Late (vs ref)",
  "PC1.climate" = "Climate PC1",
  "crop.stagecereal/ley" = "Crop Stage: Cereal/Ley (vs ref)",
  "crop.stageosr" = "Crop Stage: Oilseed Rape (vs ref)",
  "Site2" = "Site 2 (vs ref)",
  "Site3" = "Site 3 (vs ref)",
  "TreatmentControl:SeasonLate" = "Treatment: Control x Season: Late (vs ref)",
  "Age:SeasonLate" = "Age x Season: Late (vs ref)" 
)

models <- expand.grid(HillChao=c(0:2), Fixed=c("Fixed", "Not Fixed"))

models$Global <- paste0("global_mod_q", models$HillChao)

models$Avg <- ifelse(models$Fixed=="Fixed", paste0("model_avg_fixed_q",models$HillChao), paste0("model_avg_q",models$HillChao)) 


# run extraction across rows
res_list <- lapply(seq_len(nrow(models)), function(i) {

  out <- extract_model_avg(
    model_avg    = get(models$Avg[i]),
    global_model = get(models$Global[i]),   # your global fitted model
    hillchao     = models$HillChao[i],
    fixed        = models$Fixed[i],
    full         = F,
    pretty_names = pretty,   # or NULL
    present_char = "."
  )
  out
})

# bind into the two overall dataframes
variables <- do.call(rbind, lapply(res_list, `[[`, "vars"))
models_df <- do.call(rbind, lapply(res_list, `[[`, "models"))

# format vars dataframe into table for SM
var.table <- variables %>%
  mutate(
    CI = ifelse(
      is.na(Lwr) | is.na(Upr),
      NA_character_,
      paste0("[", round(Lwr, 3), ", ", round(Upr, 3), "]")
    ),
    Estimate = round(Estimate, 3),
    SE       = round(SE, 3),
    Weight   = round(Weight, 2)
  ) %>%
  dplyr::select(Fixed, Parameter, HillChao, Estimate, SE, CI, Weight) %>%
  pivot_wider(
    names_from  = HillChao,
    values_from = c(Estimate, SE, CI, Weight),
    names_glue  = "q{HillChao}_{.value}"
  ) %>%
  dplyr::select(Fixed, Parameter,
         q0_Estimate, q0_SE, q0_CI, q0_Weight,
         q1_Estimate, q1_SE, q1_CI, q1_Weight,
         q2_Estimate, q2_SE, q2_CI, q2_Weight) %>%
  arrange(Fixed)


################################################################################
# save model results and predictions
################################################################################
write.csv(variables, "results/08_bee-alpha-vars.csv", row.names=F)
write.csv(var.table, "results/08_bee-alpha-vars-table.csv", row.names=F)
write.csv(models_df, "results/08_bee-alpha-models.csv", row.names=F)
write.csv(preds.all, "results/08_bee-alpha-predictions.csv", row.names=F)

