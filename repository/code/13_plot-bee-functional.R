rm(list=ls())

library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(diagis)
library(ggh4x) # for customising facets

# load data
sol.func <- read.csv("data/cleaned/07_solitary-func-div.csv")
pred <- read.csv("results/09_bee-func-predictions.csv")

# format labels
sol.func <- sol.func %>%
  mutate(Order.q=dplyr::recode(Order.q, 
                               "0" = "q=0 (Richness)",
                               "1"="q=1 ('Shannon')", 
                               "2" = "q=2 ('Simpson')"),
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

sol.func$Site <- as.factor(sol.func$Site)

#############
# solitary bee functional diversity
#############
pd <- position_dodge(width = 0.3)

pred <- pred %>%
  filter(!(Season == "Late season" & Age_raw>10))

sol.func.plot <-
  ggplot(sol.func, aes(x = Age, y = qFD)) +
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
    aes(ymin = qFD - SE, ymax = qFD + SE, group = Treatment),
    width = .8,
    linewidth = 0.45,
    colour = "grey25",
    alpha = 0.9,
    position = pd
  ) +
  geom_line(data=pred, aes(x=Age_raw, y=qFD_pred, linetype=Treatment, group=Treatment),
            linewidth=.8,
            colour="black") +
    facet_grid2(factor(Order.q,
                       levels = c("q=0 (Richness)", "q=1 ('Shannon')", "q=2 ('Simpson')")) ~ Season) +
  scale_shape_manual(values = c(22, 21, 24)) +   
  scale_fill_manual(values=c("grey20","white" )) +
    # coord_cartesian(ylim = c(0, 22.5), clip = "on") +
  theme_bw() + 
  ylab("Solitary bee functional diversity (Hill–Chao)") + 
  xlab("System age (years)") + 
  ggtitle("A)") +
  guides(
    fill = guide_legend(
      title = "Treatment",
      override.aes = list(
        shape = 21,
        colour = "black",
        size = 4,
        func = 1
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
    scale_x_continuous(breaks=seq(min(sol.func$Age), max(sol.func$Age), by=2))
  

sol.func.plot

# ggsave(plot=sol.func.plot, filename="figs/13_solitaries-func-div.pdf", width=20, height=20, units="cm", dpi=300)

#############
# solitary bee functional model results
#############

# load data
bee_func <- read.csv("results/09_bee-func-vars.csv")

# factorise categorical vars
cols <- c("Parameter", "term", "HillChao", "Fixed")
bee_func[cols] <- lapply(bee_func[cols], factor)  ## as.factor() could also be used

# reorder parameter levels
bee_func$Parameter <- factor(bee_func$Parameter, levels=c("Age (yrs)", 
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

bee_func$Fixed <- factor(bee_func$Fixed, levels=c("Not Fixed", "Fixed"))

levels(bee_func$Parameter)[c(7,8,9)] <- c("Weather PC1", "Crop: Cereal/Ley (vs ref)", "Crop: Oilseed Rape (vs ref)")

# bold fixed parameters
levels(bee_func$Parameter) # first three are fixed
for(i in c(1:3)){
  levels(bee_func$Parameter)[i] <- paste0("**", levels(bee_func$Parameter)[i], "**")
}

# signif labels 
bee_func$signif <- ""
for(i in 1:nrow(bee_func)){
  if(!is.na(bee_func$Weight[i]) & bee_func$Lwr[i]*bee_func$Upr[i]>0 & round(bee_func$Weight[i],2)>=0.6){
    bee_func$signif[i] <- "*"
  } }

# plot
func_plot <-
  ggplot(bee_func, aes(x = Parameter, y = Estimate, shape=HillChao)) +
  
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
    # title = "Solitary Bees: func",
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
  geom_text(data = bee_func, inherit.aes = FALSE,
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
    xintercept = seq(1.5, length(unique(bee_func$Parameter)) - 0.5, by = 1),
    color = "grey80",
    linewidth = 0.5
  ) +
  ggtitle("B)")

# ggsave(plot=func_plot, filename="figs/13_solitaries-func-beta-estimates.pdf", width=22, height=14, units="cm", dpi=300)

##############
# Combine plots
##############
grid.arrange(sol.func.plot, func_plot, ncol=1, heights=c(3,2))
combined.plot <- arrangeGrob(sol.func.plot, func_plot, ncol=1, heights=c(4.5,3))

ggsave(plot=combined.plot, filename="figs/13_solitaries-func-combined.pdf", width=24, height=28, units="cm", dpi=300)
