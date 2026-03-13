rm(list=ls())

library(dplyr)
library(ggplot2)
library(viridis)
library(patchwork)
library(diagis)

# load data
syr.alpha <- read.csv("data/cleaned/05_syrphidae-alpha-div.csv")
pred <- read.csv("results/10_hoverfly-alpha-predictions.csv")

# format labels
syr.alpha$Order.q <- factor(syr.alpha$Order.q, levels=c("Species Richness (q=0)",
                                                        "Shannon Diversity (q=1)",
                                                        "Simpson Diversity (q=2)"))
pred$Order.q <- factor(pred$Order.q, levels=levels(syr.alpha$Order.q))
                                                        
syr.alpha <- syr.alpha %>%
  mutate(Order.q=dplyr::recode(Order.q, 
                               "Species Richness (q=0)" = "q=0 (Richness)",
                               "Shannon Diversity (q=1)"="q=1 ('Shannon')", 
                               "Simpson Diversity (q=2)" = "q=2 ('Simpson')")) 
pred <- pred %>%
  mutate(Order.q=dplyr::recode(Order.q, 
                               "Species Richness (q=0)" = "q=0 (Richness)",
                               "Shannon Diversity (q=1)"="q=1 ('Shannon')", 
                               "Simpson Diversity (q=2)" = "q=2 ('Simpson')")) 

syr.alpha$Site <- factor(syr.alpha$Site)

#############
# plot
#############

pd <- position_dodge(width = 0.3)

syr.alpha.plot <-   
  ggplot(syr.alpha, aes(x = Age, y = qD)) +
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
  geom_line(data=pred, aes(x=Age_raw, y=qD_pred, linetype=Treatment, group=Treatment),
            linewidth=.8,
            colour="black") +
  facet_grid(~ factor(Order.q,
                      levels = c("q=0 (Richness)", "q=1 ('Shannon')", "q=2 ('Simpson')"))) +
  scale_shape_manual(values = c(22, 21, 24)) +   
  scale_fill_manual(values=c("grey20","white" )) +
  # coord_cartesian(ylim = c(0, 22.5), clip = "on") +
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
  scale_x_continuous(breaks=seq(min(syr.alpha$Age), max(syr.alpha$Age), by=2))

syr.alpha.plot

ggsave(plot=syr.alpha.plot, filename="figs/14_hoverflies-alpha-div.pdf", width=20, height=10, units="cm", dpi=300)





