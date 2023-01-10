## ---------------------------
##
## Script name: Infaunal communities associated with mussel beds depending on their level of
##              euendolithic infestation and Perna lineages - DATA VISUALIZATION
##
## Purpose of script: Create beautiful graphical illustrations for the manuscript/thesis
##
## Author: Alexia Dievart
##
## Date Created: 2022-01-05
## Dates Updated: 2022-01-06
##
## Copyright (c) Alexia DIEVART 2023
## Email: alexia.dievart@hotmail.fr
##
## ---------------------------
##
## Notes: 
##   - Infestation levels : quadrats pooled for each site (n < 5)
##   - Perna lineage : quadrats pooled for Old Woman's River (n < 5)
## ---------------------------

###############################################################################################
# Section: Session setup ----------------------------------------------------------------------
###############################################################################################

# Set working directory
setwd("D:/Little Bull/Etudes/phD - Rhodes/1.1_QUADRATS/STATS")

# Load packages
library(pacman)
pacman::p_load(tidyverse,
               DHARMa,
               lme4,
               mgcv,
               gratia,
               ggplot2,
               ggtext,
               tidyr,
               dplyr,
               MuMIn,
               glue,
               glmmTMB,
               ggeffects,
               ggpubr,
               ggrepel,
               rstatix,
               multcomp,
               car)

# Set default ggplot theme
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA),
                  axis.text = element_text(colour = "black", size = 12),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 4), "mm"), size = 14),
                  legend.position = "none"))
my_colors <- c("darkgrey", "brown3")
site_order <- c("Mosselbaai", "Brenton-on-Sea", "Jeffreysbaai", "Old Woman's River", "Port Edward")
sites <- c("MB", "BoS", "JB", "OWR", "PEd")
inf <- c("Inf", "Non-inf")


###############################################################################################
# Section: Load data --------------------------------------------------------------------------
###############################################################################################

# Environmental variables - e.g. nb of live Perna, nb of byssal threads =======================
env <- read.csv("./RAW DATA/infauna_description.csv", dec = ",", header = T, sep = ";")
head(env)
dplyr::glimpse(env)
env <- env[!(is.na(env$live_perna)), ] # omit all missing quadrats

byssus <- read.csv("./RAW DATA/infauna_architecture.csv", dec = ",", header = T, sep = ";")
head(byssus)
byssus <- na.omit(byssus)

factor1 <- c("site", "infestation", "lineage", "quadrat")
env[,factor1] <- lapply(env[,factor1], factor)
byssus[,factor1] <- lapply(byssus[,factor1], factor)

View(env)
View(byssus)

# Community variables - e.g. abundance, biomass ==============================================
community <- read.csv("./RAW DATA/infauna_community.csv", dec = ",", header = T, sep = ";")

factor2 <- c("site", "bioregion", "infestation", "lineage", "quadrat", "species", "pic")
community[,factor2] <- lapply(community[,factor2], factor)
dplyr::glimpse(community)

abund_long <- community[,1:7]
abund_wide <- spread(abund_long, species, abundance)
abund_wide[is.na(abund_wide)] <- 0

biom_long <- community[,c(1:6,8)]
biom_wide <- spread(biom_long, species, biomass)
biom_wide[is.na(biom_wide)] <- 0

View(abund_wide)
View(biom_wide)

## Isolate Old Woman's River for analysis involving Perna lineage as a factor
byssus_owr <- byssus %>%
  filter(site == "Old Woman's River")



###############################################################################################
# Section: Visualize data ---------------------------------------------------------------------
###############################################################################################

# Mussel bed architecture =====================================================================

### Combination of multiple graphs

architecture <- ggarrange(barplot_livemussels2, barplot_deadmussels2, barplot_brokenshells2, boxplot_bt2,
          barplot_live_owr1, barplot_dead_owr1, barplot_broken_owr1, boxplot_bt_owr1,
          ncol = 4, nrow = 2, common.legend = F, legend = "none",
          labels = c("A", "B", "C", "D", "E", "F", "G", "H"), hjust = -2,
          font.label = list(size = 20, color = "grey30"))

annotate_figure(architecture, 
                left = text_grob("    Old Woman's River                                 Across all sites      ",
                                 rot = 90, size = 22, face = "bold"))



## Average number of live mussels per quadrat #################################################

### Across all sites (infested vs non-infested)

barplot_livemussels2 <- ggplot(data = live_mussels, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = mean_live, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_live - se_live, ymax = mean_live + se_live), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 44)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 44,
           label = c("A", "AB", "AB", "AB", "B")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of live mussels")
barplot_livemussels2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_live_owr1 <- ggplot(live_mussels_owr, aes(x = infestation, y = mean_live, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_live - se_live, ymax = mean_live + se_live), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 44)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of live mussels")
barplot_live_owr1

## Average number of dead mussels per quadrat #################################################

### Across all sites (infested vs non-infested)

barplot_deadmussels2 <- ggplot(data = dead_mussels, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = mean_dead, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_dead - se_dead, ymax = mean_dead + se_dead), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 7)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 7,
           label = c("C", "A", "B", "AB", "C")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of dead mussels")
barplot_deadmussels2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_dead_owr1 <- ggplot(dead_mussels_owr, aes(x = infestation, y = mean_dead, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_dead - se_dead, ymax = mean_dead + se_dead), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 7)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of dead mussels")
barplot_dead_owr1

## Average number of broken shells per quadrat ################################################

### Across all sites (infested vs non-infested)

barplot_brokenshells2 <- ggplot(data = broken_mussels, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = mean_broken, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_broken - se_broken, ymax = mean_broken + se_broken), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 34)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 34,
           label = c("A", "A", "A", "B", "A")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of broken shells")
barplot_brokenshells2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_broken_owr1 <- ggplot(broken_mussels_owr, aes(x = infestation, y = mean_broken, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_broken - se_broken, ymax = mean_broken + se_broken), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 34)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of broken shells")
barplot_broken_owr1

## Average number of byssal threads per mussels ###############################################

### Across all sites (infested vs non-infested)

boxplot_bt2 <- ggplot(byssus, aes(x = site, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(limits = c(0, 430)) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 430,
           label = c("AC", "B", "A", "B", "BC")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average number of byssal threads")
boxplot_bt2

### Old Woman's River (infested vs non-infested, eastern vs western)

boxplot_bt_owr1 <- ggplot(byssus_owr, aes(x = infestation, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 410)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of byssal threads")
boxplot_bt_owr1

# Community descriptors ================================================================

### Combination of multiple graphs

description <- ggarrange(barplot_abund_tot2, barplot_biom_tot2, barplot_rich2, barplot_shannon2, barplot_pielou2,
                         barplot_abund_tot_owr1, barplot_biom_tot_owr2, barplot_rich_owr1, barplot_shannon_owr1, barplot_pielou_owr1,
                          ncol = 5, nrow = 2, common.legend = F, legend = "none",
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J"), 
                          hjust = c(-2, -2, -1.2, -1.2, -1.2, -2, -1.5, -1.2, -2.5, -1.5),
                          font.label = list(size = 20, color = "grey30"))

annotate_figure(description, 
                left = text_grob("    Old Woman's River                                 Across all sites      ",
                                 rot = 90, size = 22, face = "bold"))

## Total abundance #####################################################################

### Across all sites (infested vs non-infested)

barplot_abund_tot2 <- ggplot(data = abund_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = ab_tot_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_se, ymax = ab_tot_mean + ab_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1050)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(0.85, 1.3, 1.95, 2.32, 2.95, 3.4, 3.85, 4.35, 4.9, 5.35),
           y = c(300, 300, 700, 1050, 700, 580, 500, 460, 220, 175),
           label = c("B", "B", "AB", "A", "AB", "AB", "B", "B", "B", "B"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total abundance")
barplot_abund_tot2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_abund_tot_owr1 <- ggplot(abund_tot_owr1, aes(x = infestation, y = ab_tot_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_se, ymax = ab_tot_mean + ab_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250),limits = c(0, 255)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total abundance")
barplot_abund_tot_owr1

## Total biomass ########################################################################

### Across all sites (infested vs non-infested)

barplot_biom_tot2 <- ggplot(data = biom_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = biom_tot_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_tot_mean - biom_tot_se, ymax = biom_tot_mean + biom_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 16100), breaks = c(0, 5000, 10000, 15000)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(0.95, 1.4, 1.85, 2.3, 3, 3.4, 4, 4.4, 5, 5.3),
           y = c(3700, 4300, 16100, 11600, 8400, 7800, 8500, 7900, 2400, 2300),
           label = c("AB", "AB", "A", "A", "AB", "AB", "AC", "AC", "BC", "B"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total biomass (mg)")
barplot_biom_tot2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_biom_tot_owr2 <- ggplot(biom_tot_owr1, aes(x = infestation, y = biom_tot_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_tot_mean - biom_tot_se, ymax = biom_tot_mean + biom_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total biomass (mg)")
barplot_biom_tot_owr2




## Species richness #####################################################################

### Across all sites (infested vs non-infested)

barplot_rich2 <- ggplot(data = richness1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = count_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_se, ymax = count_mean + count_se), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 55), breaks = c(0, 10, 20, 30, 40, 50)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 55,
           label = c("A", "AB", "AB", "B", "AB")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average species richness (S)")
barplot_rich2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_rich_owr1 <- ggplot(richness_owr1, aes(x = infestation, y = count_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_se, ymax = count_mean + count_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(limits = c(0, 55), breaks = c(0, 10, 20, 30, 40, 50)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average species richness (S)")
barplot_rich_owr1

## Species diversity - Shannon's Index H' ###############################################

### Across all sites (infested vs non-infested)

barplot_shannon2 <- ggplot(data = shannon, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = shannon_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_se, ymax = shannon_mean + shannon_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 4.2), breaks = c(1, 2, 3, 4)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_shannon_owr1 <- ggplot(shannon_owr, aes(x = infestation, y = shannon_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_se, ymax = shannon_mean + shannon_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(limits = c(0, 4.2), breaks = c(1, 2, 3, 4)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon_owr1

## Species diversity - Simpson's Index (λ) ###############################################

### Across all sites (infested vs non-infested)

barplot_simpson2 <- ggplot(data = simpson, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = simpson_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_se, ymax = simpson_mean + simpson_se), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.5, 1, 1.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Simpson's diversity (λ) ")
barplot_simpson2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_simpson_owr1 <- ggplot(simpson_owr, aes(x = infestation, y = simpson_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_se, ymax = simpson_mean + simpson_se), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.5, 1, 1.5)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(labels = inf) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Simpson's diversity (λ)")
barplot_simpson_owr1

## Species evenness - Pielou's Index (J) ###############################################

### Across all sites (infested vs non-infested)

barplot_pielou2 <- ggplot(data = pielou, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = pielou_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_se, ymax = pielou_mean + pielou_se), width = 0.2,
                position = position_dodge(0.9)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.5, 1, 1.5)) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 1.6,
           label = c("AB", "AB", "A", "AB", "B")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Pielou's evenness (J)")
barplot_pielou2

### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_pielou_owr1 <- ggplot(pielou_owr, aes(x = infestation, y = pielou_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_se, ymax = pielou_mean + pielou_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(limits = c(0, 1.6), breaks = c(0, 0.5, 1, 1.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Pielou's evenness (J)")
barplot_pielou_owr1

# Community analysis on abundance =================================================================

### Combination of multiple graphs

par(mfrow = c(2,3), oma=c(1,6,1,1))

### Across all sites (infested vs non-infested)

ordiplot(abund_bray, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(abund_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5)[env1$site],
       col = my_colors[env1$infestation])
ordihull(abund_bray, groups = env1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.4, 0.5, legend = c(levels(env1$site), levels(env1$infestation)),
       pch = c(16, 17, 18, 8, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "black", "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(-1.65, 1.2, "Stress = 0.180", bty = "n", cex = 1)

mtext("A", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(abund_disp_inf, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     xlab = "Dimension 1 (22.31 %)", ylab = "Dimension 2 (19.24%)",
     main = ""
)

mtext("B", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(abund_disp_site, hull = F, ellipse = T,
     pch = c(16:18, 8, 5), col = "black", label = F, cex = 2,
     xlab = "Dimension 1 (22.31%)", ylab = "Dimension 2 (19.24%)",
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     main = "")
legend(-0.57, 0.45, legend = levels(env1$site), pch = c(16, 17, 18, 8, 5),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)

mtext("C", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

### Old Woman's River (infested vs non-infested, eastern vs western)

ordiplot(abund_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(abund_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env1_owr$lineage],
       col = my_colors[env1_owr$infestation])
ordihull(abund_bray_owr, groups = env1_owr$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.4, -0.1, legend = c(levels(env1_owr$lineage), levels(env1_owr$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(-1.65, 1.2, "Stress = 0.116", bty = "n", cex = 1)

mtext("D", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

abund_owr_inf_plot <- plot(abund_disp_inf_owr, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21%)",
     main = ""
)

mtext("E", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

abund_owr_lineage_plot <- plot(abund_disp_lineage, hull = F, ellipse = T,
     pch = c(16:17), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21 %)",
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     main = "")

mtext("F", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

mtext("Across all sites", side = 2, line = 90, cex = 1.5, adj = 2.85, col = "black")
mtext("Old Woman's River", side = 2, line = 90, cex = 1.5, adj = 0.55, col = "black")




# Community analysis on biomass =================================================================

### Combination of multiple graphs

par(mfrow = c(2,3), oma=c(1,6,1,1))

### Across all sites (infested vs non-infested)

ordiplot(biom_bray, type = "none", xlim = c(-1.5, 1.5), ylim = c(-1.0, 1.0))
points(biom_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5)[env1$site],
       col = my_colors[env1$infestation])
ordihull(biom_bray, groups = env1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.65, 0.45, legend = c(levels(env1$site), levels(env1$infestation)),
       pch = c(16, 17, 18, 8, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "black", "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(-1.95, 1.4, "Stress = 0.180", bty = "n", cex = 1)

mtext("A", side = 3, line = -1, cex = 2, adj = -0.1, col = "grey30")

plot(biom_disp_inf, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5),
     xlab = "Dimension 1 (24.71 %)", ylab = "Dimension 2 (18.63 %)",
     main = ""
)

mtext("B", side = 3, line = -1, cex = 2, adj = -0.1, col = "grey30")

plot(biom_disp_site, hull = F, ellipse = T, label.cex = 2, label = F,
     pch = c(16:18, 8, 5), col = "black", cex = 2,
     xlab = "Dimension 1 (24.71 %)", ylab = "Dimension 2 (18.63 %)",
     xlim = c(-0.7, 0.7), ylim = c(-0.5, 0.5),
     main = "")
legend(-0.75, 0.0, legend = levels(env1$site), pch = c(16, 17, 18, 8, 5),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)

mtext("C", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

ordiplot(biom_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(biom_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env1_owr$lineage],
       col = my_colors[env1_owr$infestation])
ordihull(biom_bray_owr, groups = env1_owr$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.4, -0.1, legend = c(levels(env1_owr$lineage), levels(env1_owr$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 2)
legend(-1.67, 1.25, "Stress = 0.116", bty = "n", cex = 1)

mtext("D", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(biom_disp_inf_owr, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     xlab = "Dimension 1 (44.42 %)", ylab = "Dimension 2 (15.75 %)",
     main = ""
)

mtext("E", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(biom_disp_lineage, hull = F, ellipse = T,
     pch = c(16:17), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21 %)",
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     main = "")

mtext("F", side = 3, line = -1, cex = 2, adj = -0.12, col = "grey30")

mtext("Across all sites", side = 2, line = 90, cex = 1.5, adj = 2.9, col = "black")
mtext("Old Woman's River", side = 2, line = 90, cex = 1.5, adj = 0.50, col = "black")





# Cluster analysis =================================================================

par(mfrow = c(2,2), mar = c(3 + 1, 3 + 6, 1, 1))

## Cluster analysis - all sites on abundance
abund_clust <- hclust(abund_dist, "ward.D2")
levels(abund_wide$site) <- c("BoS", "JB", "MB", "OWR", "PEd")

abund_inf <- rev(levels(abund_wide$infestation))
abund_dend <- as.dendrogram(abund_clust)

labels_colors(abund_dend) <- my_colors[(abund_wide$infestation)[order.dendrogram(abund_dend)]]

labels(abund_dend) <- paste(abund_wide$site[order.dendrogram(abund_dend)],
                            " (", labels(abund_dend), ")", sep = "")
branches_attr_by_lists(abund_dend, abund_wide$site[order.dendrogram(abund_dend)],
                       attr = "lwd")
abund_dend <- set(abund_dend, "labels_cex", 1.1)
abund_dend <- hang.dendrogram(abund_dend, hang_height = 0.1)
plot(abund_dend, nodePar = list(cex = 0.007))

mtext("A", side = 3, line = -2, cex = 2, adj = -0.1, col = "grey30")

## Cluster analysis - all sites on biomass
biom_clust <- hclust(biom_dist, "ward.D2")
levels(biom_wide$site) <- c("BoS", "JB", "MB", "OWR", "PEd")

biom_inf <- rev(levels(biom_wide$infestation))
biom_dend <- as.dendrogram(biom_clust)

labels_colors(biom_dend) <- my_colors[(biom_wide$infestation)[order.dendrogram(biom_dend)]]

labels(biom_dend) <- paste(biom_wide$site[order.dendrogram(biom_dend)],
                           " (", labels(biom_dend), ")", sep = "")
biom_dend <- set(biom_dend, "labels_cex", 1.1)
biom_dend <- hang.dendrogram(biom_dend, hang_height = 0.005)
plot(biom_dend, nodePar = list(cex = 0.007))
legend(35, 2, legend = levels(abund_wide$infestation), fill = my_colors, bty = "n", cex = 1.5)

mtext("B", side = 3, line = -2, cex = 2, adj = -0.1, col = "grey30")

## Cluster analysis - OWR abundance
abund_clust_owr <- hclust(abund_dist_owr, "ward.D2")

abund_inf_owr <- rev(levels(abund_wide_owr$infestation))
abund_dend_owr <- as.dendrogram(abund_clust_owr)

labels_colors(abund_dend_owr) <- my_colors[(abund_wide_owr$infestation)[order.dendrogram(abund_dend_owr)]]

labels(abund_dend_owr) <- paste(abund_wide_owr$lineage[order.dendrogram(abund_dend_owr)],
                                " (", labels(abund_dend_owr), ")", sep = "")
abund_dend_owr <- set(abund_dend_owr, "labels_cex", 1.1)
abund_dend_owr <- hang.dendrogram(abund_dend_owr, hang_height = 0.005)
plot(abund_dend_owr, nodePar = list(cex = 0.007))

mtext("C", side = 3, line = -1.5, cex = 2, adj = -0.1, col = "grey30")

## Cluster analysis - OWR biomass
biom_clust_owr <- hclust(biom_dist_owr, "ward.D2")

biom_inf_owr <- rev(levels(biom_wide_owr$infestation))
biom_dend_owr <- as.dendrogram(biom_clust_owr)

labels_colors(biom_dend_owr) <- my_colors[(biom_wide_owr$infestation)[order.dendrogram(biom_dend_owr)]]

labels(biom_dend_owr) <- paste(biom_wide_owr$lineage[order.dendrogram(biom_dend_owr)],
                               " (", labels(biom_dend_owr), ")", sep = "")
biom_dend_owr <- set(biom_dend_owr, "labels_cex", 1.1)
biom_dend_owr <- hang.dendrogram(biom_dend_owr, hang_height = 0.005)
plot(biom_dend_owr, nodePar = list(cex = 0.007))

mtext("D", side = 3, line = -1.5, cex = 2, adj = -0.1, col = "grey30")

mtext("Across all sites", side = 2, line = 58, cex = 1.5, adj = 2.3, col = "black")
mtext("Old Woman's River", side = 2, line = 58, cex = 1.5, adj = 0.50, col = "black")

# Species that contribute the most to dissimilarity in abundance and biomass across all sites ====

contribution <- ggarrange(barplot_abund_cont, barplot_biom_cont,
                         ncol = 1, nrow = 2, common.legend = F, legend = "none",
                         labels = c("A", "B"), 
                         hjust = c(-2, -2.2),
                         font.label = list(size = 20, color = "grey30"))

annotate_figure(contribution, 
                left = text_grob("Across all sites",
                                 rot = 90, size = 22, face = "bold"))

## Abundance
View(abund_long)

abund_cont <- abund_long[abund_long$species %in% c("Acanthochitona garnoti",
                      "Actiniaria",
                      "Burnupena lagenaria",
                      "Helcion pruinosus",
                      "Mytilus galloprovincialis",
                      "Nereididae",
                      "Perna perna",
                      "Pseudonereis variegata",
                      "Tricolia neritina"),]
View(abund_cont)

abund_cont1 <- abund_cont %>%
    dplyr::group_by(site, infestation, species) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance),
    n = n(),
    se_ab  = sd_ab / sqrt(n)
  )
View(abund_cont1)

barplot_abund_cont <- ggplot(abund_cont1, aes(x = site, y = mean_ab, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - se_ab, ymax = mean_ab + se_ab), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = abund_cont_text, mapping = aes(x = x, y = y, label = label)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250, 300)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ species, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average abundance")
barplot_abund_cont

abund_cont_text <- data.frame(
  label = c("A", "AB", "B", "ABC", "C", 
            "A", "B", "A", "A", "A",
            "A", "B", "B", "B", "C",
            "A", "A", "A", "BC", "C",
            "A", "B", "AB", "C", "C",
            "AC", "AB", "BC", "B", "C",
            "AC", "A", "A", "B", "C",
            "AB", "A", "C", "AB", "B",
            "A", "B", "A", "A", "A"),
  species = c("Acanthochitona garnoti", "Acanthochitona garnoti", "Acanthochitona garnoti", "Acanthochitona garnoti", "Acanthochitona garnoti",
              "Actiniaria", "Actiniaria", "Actiniaria", "Actiniaria", "Actiniaria",
              "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria",
              "Helcion pruinosus",  "Helcion pruinosus",  "Helcion pruinosus",  "Helcion pruinosus",  "Helcion pruinosus",
              "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis",
              "Nereididae", "Nereididae", "Nereididae", "Nereididae", "Nereididae",
              "Perna perna", "Perna perna", "Perna perna", "Perna perna", "Perna perna",
              "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata",
              "Tricolia neritina", "Tricolia neritina", "Tricolia neritina", "Tricolia neritina", "Tricolia neritina"),
  infestation = c("infested"),
  x = c(1, 2, 3, 4, 5),
  y = c(20, 40, 30, 20, 10, 
        20, 80, 30, 20, 20,
        15, 50, 20, 20, 15,
        15, 10, 15, 20, 25,
        20, 185, 40, 10, 10,
        10, 20, 20, 30, 10,
        80, 330, 270, 20, 30,
        20, 40, 10, 20, 20,
        10, 40, 10, 10, 10)
)

## Biomass
View(biom_long)

biom_cont <- biom_long[biom_long$species %in% c("Acanthochitona garnoti",
                                                   "Burnupena lagenaria",
                                                   "Ischyromene huttoni",
                                                   "Mytilus galloprovincialis",
                                                   "Parvulastra exigua",
                                                   "Perna perna",
                                                   "Pseudonereis variegata"),]
View(biom_cont)

biom_cont1 <- biom_cont %>%
  dplyr::group_by(site, infestation, species) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    sd_biom = sd(biomass),
    n = n(),
    se_biom  = sd_biom / sqrt(n)
  )
View(biom_cont1)

barplot_biom_cont <- ggplot(biom_cont1, aes(x = site, y = mean_biom, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom - se_biom, ymax = mean_biom + se_biom), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  geom_text(data = biom_cont_text, mapping = aes(x = x, y = y, label = label)) +
  scale_y_continuous(breaks = c(0, 1000, 2000, 3000, 4000, 5000, 6000, 7000, 8000, 9000, 10000, 11000)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ species, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average biomass (mg)")
barplot_biom_cont

biom_cont_text <- data.frame(
  label = c("A", "A", "A", "A", "B",
            "A", "B", "B", "B", "A",
            "A", "B", "B", "B", "A",
            "A", "B", "AB", "C", "C",
            "A", "B", "B", "AB", "C",
            "A", "B", "B", "A", "A",
            "A", "A", "B", "A", "C"),
  species = c("Acanthochitona garnoti", "Acanthochitona garnoti", "Acanthochitona garnoti", "Acanthochitona garnoti", "Acanthochitona garnoti",
              "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria",
              "Ischyromene huttoni", "Ischyromene huttoni", "Ischyromene huttoni", "Ischyromene huttoni", "Ischyromene huttoni",
              "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis",
              "Parvulastra exigua", "Parvulastra exigua", "Parvulastra exigua", "Parvulastra exigua", "Parvulastra exigua",
              "Perna perna", "Perna perna", "Perna perna", "Perna perna", "Perna perna",
              "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata"),
  infestation = c("infested"),
  x = c(1, 2, 3, 4, 5),
  y = c(550, 1150, 800, 900, 300,
        300, 900, 500, 900, 400,
        500, 1500, 2200, 1800, 500,
        500, 2400, 800, 300, 300,
        1900, 500, 500, 900, 400,
        1100, 11100, 4300, 500, 500,
        500, 900, 300, 700, 500)
)
