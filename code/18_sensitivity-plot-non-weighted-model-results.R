rm(list=ls())

# load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(GGally) # for correlation ggpairs plot
library(diagis) # for weighted SE
library(gridExtra)
library(stringr)
library(ggtext)

# load data
bee_alpha <- read.csv("results/16_sensitivity-bee-alpha-no-weightings-vars.csv")
bee_func <- read.csv("results/17_sensitivity-bee-func-no-weightings-vars.csv")

# factorise categorical vars
cols <- c("Parameter", "term", "HillChao", "Fixed")
bee_alpha[cols] <- lapply(bee_alpha[cols], factor)  ## as.factor() could also be used
bee_func[cols] <- lapply(bee_func[cols], factor)  ## as.factor() could also be used

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

bee_alpha$Fixed <- factor(bee_alpha$Fixed, levels=c("Not Fixed", "Fixed"))
levels(bee_alpha$Parameter)[c(7,8,9)] <- c("Weather PC1", "Crop: Cereal/Ley (vs ref)", "Crop: Oilseed Rape (vs ref)")

bee_func$Fixed <- factor(bee_func$Fixed, levels=c("Not Fixed", "Fixed"))
levels(bee_func$Parameter)[c(7,8,9)] <- c("Weather PC1", "Crop: Cereal/Ley (vs ref)", "Crop: Oilseed Rape (vs ref)")

# bold fixed parameters
levels(bee_alpha$Parameter) # first three are fixed
for(i in c(1:3)){
  levels(bee_alpha$Parameter)[i] <- paste0("**", levels(bee_alpha$Parameter)[i], "**")
}

levels(bee_func$Parameter) # first three are fixed
for(i in c(1:3)){
  levels(bee_func$Parameter)[i] <- paste0("**", levels(bee_func$Parameter)[i], "**")
}

# signif labels 
bee_alpha$signif <- ""
for(i in 1:nrow(bee_alpha)){
  if(!is.na(bee_alpha$Weight[i]) & bee_alpha$Lwr[i]*bee_alpha$Upr[i]>0 & round(bee_alpha$Weight[i],2)>=0.6){
    bee_alpha$signif[i] <- "*"
  } }

bee_func$signif <- ""
for(i in 1:nrow(bee_func)){
  if(!is.na(bee_func$Weight[i]) & bee_func$Lwr[i]*bee_func$Upr[i]>0 & round(bee_func$Weight[i],2)>=0.6){
    bee_func$signif[i] <- "*"
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
  ggtitle("A) Solitary bee alpha diversity")

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
    xintercept = seq(1.5, length(unique(bee_alpha$Parameter)) - 0.5, by = 1),
    color = "grey80",
    linewidth = 0.5
  ) +
  ggtitle("B) Solitary bee functional diversity")

##############
# Combine plots
##############
grid.arrange(alpha_plot, func_plot, ncol=1)
combined.plot <- arrangeGrob(alpha_plot, func_plot, ncol=1)

ggsave(plot=combined.plot, filename="figs/18_sensitivity-solitaries-no-weightings.pdf", width=24, height=25, units="cm", dpi=300)
