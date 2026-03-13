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

# source functions
source("code/functions/08_get-predictions-function.R")
source("code/functions/08_extract-model-avg-results-function.R")

# load data
sol.div <- read.csv("data/cleaned/06_solitary-alpha-div.csv")
model.vars <- read.csv("data/cleaned/04_model-vars-seasonal.csv")

# remove late-season outlier
sol.div <- subset(sol.div, !(Site==3 & Season=="Late" & Year==2018))

# combine diversity data with model vars
sol.div <- left_join(sol.div, model.vars)

# # plot
# ggplot(sol.div, aes(x=Age, y=qD, color=Treatment)) +
#   geom_smooth(method="lm", se=F) +
#   # geom_line(alpha=0.5,aes(color=factor(Site))) +
#   geom_point(aes(shape=factor(Site)), alpha=.5, size=3) +
#   geom_errorbar(aes(ymin=qD-SE, ymax=qD+SE), alpha=.5, width=.3) +
#   theme_bw() +
#   ylab("Diversity") +
#   xlab("Age (Years)") +
#   theme(axis.title.y=element_blank(),
#         plot.margin=unit(c(.5,0,.5,.5), "cm")) +
#   ggtitle("Solitary Bee Alpha Diversity") +
#   # facet_wrap(Order.q~Season, ncol=2)
#   facet_wrap(Order.q~Season, ncol=2)

# factorise categorical variables
cols <- c("Site", "Treatment", "crop.stage", "Season", "Order.q")
sol.div[cols] <- lapply(sol.div[cols], factor)  ## as.factor() could also be used

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

# generate predictions
preds.q0<- get_ma_predictions(model_avg=model_avg_q0, 
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

# generate predictions
preds.q2<- get_ma_predictions(model_avg_fixed_q2, 
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
pred <- bind_rows(preds.q0, preds.q1, preds.q2)
pred$Order.q <- factor(pred$Order.q, levels = levels(sol.div$Order.q))

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

################################################################################
# save model results
################################################################################
# write.csv(variables, "results/19_sensitivity-bee-alpha-outlier-removed-vars.csv", row.names=F)
# write.csv(models_df, "results/19_sensitivity-bee-alpha-outlier-removed-models.csv", row.names=F)

#############
# plot
#############

# Base age plot
################

# format labels
sol.div <- sol.div %>%
  mutate(Order.q=dplyr::recode(Order.q, 
                               "Species Richness (q=0)" = "q=0 (Richness)",
                               "Shannon Diversity (q=1)"="q=1 ('Shannon')", 
                               "Simpson Diversity (q=2)" = "q=2 ('Simpson')"),
         Season=dplyr::recode(Season,
                              "Early"="Early season",
                              "Late"="Late season")) 
pred <- pred %>%
  mutate(Order.q=dplyr::recode(Order.q, 
                               "Species Richness (q=0)" = "q=0 (Richness)",
                               "Shannon Diversity (q=1)"="q=1 ('Shannon')", 
                               "Simpson Diversity (q=2)" = "q=2 ('Simpson')"),
         Season=dplyr::recode(Season,
                              "Early"="Early season",
                              "Late"="Late season"))


sol.div$Site <- factor(sol.div$Site)

pd <- position_dodge(width = 0.3)

pred <- pred %>%
  filter(!(Season == "Late season" & Age_raw>10))

sol.alpha.plot <-
  ggplot(sol.div, aes(x = Age, y = qD)) +
  geom_point(
    aes(shape = Site, fill = Treatment, group = Treatment),
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
  facet_grid(factor(Order.q,
                    levels = c("q=0 (Richness)", "q=1 ('Shannon')", "q=2 ('Simpson')")) ~ Season) +
  scale_shape_manual(values = c(22, 21, 24)) +   
  scale_fill_manual(values=c("grey20","white" )) +
  # coord_cartesian(ylim = c(0, 22.5), clip = "on") +
  theme_bw() + 
  ylab("Solitary bee alpha diversity (Hill–Chao)") + 
  xlab("System age (years)") + 
  ggtitle("A)") +
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
  scale_x_continuous(breaks=seq(min(sol.div$Age), max(sol.div$Age), by=2))


sol.alpha.plot

# ggsave(plot=sol.alpha.plot, filename="figs/19_sensitivity-bee-alpha-outlier-removed-div.pdf", width=20, height=20, units="cm", dpi=300)

# solitary bee alpha model results
##################################

# load data
bee_alpha <- read.csv("results/19_sensitivity-bee-alpha-outlier-removed-vars.csv")

# factorise categorical vars
cols <- c("Parameter", "term", "HillChao", "Fixed")
bee_alpha[cols] <- lapply(bee_alpha[cols], factor)  ## as.factor() could also be used

# reorder parameter levels
bee_alpha$Parameter <- factor(bee_alpha$Parameter, levels=c("Age (yrs)", 
                                                            "Treatment: Control (vs ref)",
                                                            "Age × Treatment: Control (vs ref)",
                                                            "Season: Late (vs ref)",
                                                            "Treatment: Control x Season: Late (vs ref)",
                                                            "Age x Season: Late (vs ref)",
                                                            "Climate PC1",
                                                            "Crop Stage: Cereal/Ley (vs ref)",
                                                            "Crop Stage: Oilseed Rape (vs ref)",
                                                            "Site 2 (vs ref)",
                                                            "Site 3 (vs ref)" 
))

bee_alpha$Fixed <- factor(bee_alpha$Fixed, levels=c("Not Fixed", "Fixed"))

levels(bee_alpha$Parameter)[c(7,8,9)] <- c("Weather PC1", "Crop: Cereal/Ley (vs ref)", "Crop: Oilseed Rape (vs ref)")

# bold fixed parameters
levels(bee_alpha$Parameter) # first three are fixed
for(i in c(1:3)){
  levels(bee_alpha$Parameter)[i] <- paste0("**", levels(bee_alpha$Parameter)[i], "**")
}

# signif labels 
bee_alpha$signif <- ""
for(i in 1:nrow(bee_alpha)){
  if(!is.na(bee_alpha$Weight[i]) & bee_alpha$Lwr[i]*bee_alpha$Upr[i]>0 & round(bee_alpha$Weight[i],2)>=0.6){
    bee_alpha$signif[i] <- "*"
  } }

# plot
alpha_plot <-
  ggplot(bee_alpha, aes(x = Parameter, y = Estimate, shape=HillChao)) +
  
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +
  
  geom_errorbar(aes(ymin = Lwr, ymax = Upr),
                position = position_dodge(width = 0.6), 
                width = 0.3,
                color="grey25",) +
  
  geom_point(position = position_dodge(width = 0.6), 
             size = 3, 
             color="grey25",
             fill="white",
             alpha=1) +
  
  facet_wrap(~Fixed, ncol=1, strip.position="right") +
  theme_bw() +
  labs(
    # title = "Solitary Bees: Alpha",
    y = "Estimate ± 95% CI",
    x = "Parameter",
    color = "Hill-Chao (q)",
    shape = "Hill-Chao (q)") +
  theme(
    legend.position = "bottom",
    axis.text.x=element_markdown( size=10,
                                  # angle=0,
                                  # vjust=1,
                                  # hjust=1
    ),
    axis.title=element_text(size=13),
    strip.text=element_text(size=13),
    legend.title=element_text(size=13),
    legend.text=element_text(size=11),
    legend.box = "horizontal",
    plot.title = element_text(size=15),
    # axis.title.y = element_blank(),
    # axis.text.y  = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  geom_text(aes(label = signif,
                y = Upr,
                x = Parameter),
            color = "grey25",
            vjust = .01,
            position = position_dodge(width = 0.6),
            size = 6, show.legend = FALSE) +
  geom_text(data = bee_alpha, inherit.aes = FALSE,
            aes(y = Lwr, x = Parameter, group = HillChao, label = format(round(Weight, 2), nsmall=2)),
            color="grey25",
            position = position_dodge(width = 0.6),
            size = 3, 
            angle=45,
            vjust = 1.35 ,
            hjust = 1,
            show.legend = FALSE) +
  coord_cartesian(ylim=c(-3.75,3.75)) +
  scale_color_grey() +
  scale_x_discrete(labels = function(x) gsub("\n", "<br>", stringr::str_wrap(x, width = 10)))  +
  scale_shape_manual(values = c(22, 21, 24)) +
  geom_vline(
    xintercept = seq(1.5, length(unique(bee_alpha$Parameter)) - 0.5, by = 1),
    color = "grey80",
    linewidth = 0.5
  ) +
  ggtitle("B)")


# ggsave(plot=alpha_plot, filename="figs/XX_sensitivity-solitaries-alpha-outlier-removed-beta-estimates.pdf", width=22, height=14, units="cm", dpi=300)

##############
# Combine plots
##############
grid.arrange(sol.alpha.plot, alpha_plot, ncol=1, heights=c(3,2))
combined.plot <- arrangeGrob(sol.alpha.plot, alpha_plot, ncol=1, heights=c(4.5,3))

ggsave(plot=combined.plot, filename="figs/19_sensitivity-solitaries-alpha-outlier-removed-combined.pdf", width=24, height=28, units="cm", dpi=300)
