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
## Dates Updated: 2022-01-06; 2023-02-05
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
                  axis.text = element_text(colour = "black", size = 14),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 4), "mm"), size = 14),
                  legend.position = "none"))
my_colors <- c("darkgrey", "brown3")
site_order <- c("Mosselbaai", "Brenton-on-Sea", "Jeffreysbaai", "Port Edward")
sites <- c("MB", "BR", "JB", "PEd")
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

architecture <- ggarrange(barplot_mussels2, boxplot_brokenshells2, boxplot_bt2,
                          barplot_mussels_owr, boxplot_broken_owr1, boxplot_bt_owr1,
          ncol = 3, nrow = 2, common.legend = F, legend = "none",
          labels = c("A", "B", "C", "D", "E", "F"), 
          hjust = c(-1, -1.2, -1.2, -1, -1.2, -1.2),
          font.label = list(size = 20, color = "grey30"))

annotate_figure(architecture, 
                left = text_grob("    Old Woman's River                                 Across all sites      ",
                                 rot = 90, size = 22, face = "bold"))



## Average number of mussels per quadrat #################################################

### Across all sites (infested vs non-infested)
mussels <- env1 %>%
  dplyr::group_by(site, infestation) %>%
  dplyr::summarise(
    mean_mussels = mean(tot_mussels),
    sd_mussels = sd(tot_mussels)
  )
head(mussels)

barplot_mussels2 <- ggplot(data = mussels, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = mean_mussels, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_mussels - sd_mussels, ymax = mean_mussels + sd_mussels), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 50)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 49,
           label = c("A", "AB", "AB", "B")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of mussels (live and dead)")
barplot_mussels2

### Old Woman's River (infested vs non-infested, eastern vs western)
mussels_owr <- env_owr %>%
  dplyr::group_by(lineage, infestation) %>%
  dplyr::summarise(
    mean_mussels = mean(tot_mussels),
    sd_mussels = sd(tot_mussels)
  )
head(mussels_owr)

barplot_mussels_owr <- ggplot(mussels_owr, aes(x = infestation, y = mean_mussels, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_mussels - sd_mussels, ymax = mean_mussels + sd_mussels), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 44)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of mussels (live and dead)")
barplot_mussels_owr

## Average number of broken shells per quadrat ################################################

### Across all sites (infested vs non-infested)
boxplot_brokenshells2 <- ggplot(data = env1, aes(x = site, y = nb_broken, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 10), breaks = seq(from = 0, to = 10, by = 2)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of broken shells")
boxplot_brokenshells2

### Old Woman's River (infested vs non-infested, eastern vs western)
boxplot_broken_owr1 <- ggplot(env_owr, aes(x = infestation, y = nb_broken, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(limits = c(0, 40)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of broken shells")
boxplot_broken_owr1

## Average number of byssal threads per mussels ###############################################

### Across all sites (infested vs non-infested)
boxplot_bt2 <- ggplot(byssus1, aes(x = site, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(limits = c(0, 430)) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 430,
           label = c("AC", "BC", "A", "C")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
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

description <- ggarrange(barplot_abund_tot2, barplot_biom_tot2, barplot_rich2, barplot_shannon2, barplot_simpson2, barplot_pielou2,
                         barplot_abund_tot_owr1, barplot_biom_tot_owr2, barplot_rich_owr1, barplot_shannon_owr1, barplot_simpson_owr1, barplot_pielou_owr1,
                          ncol = 6, nrow = 2, common.legend = F, legend = "none",
                          labels = c("A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L"), 
                          hjust = c(-1, -0.5, -0.5, -1, -0.5, -0.5,
                                    -1, -0.5, -2.5, -1, -0.5, -0.5),
                          font.label = list(size = 20, color = "grey30"))

annotate_figure(description, 
                left = text_grob("    Old Woman's River                                 Across all sites      ",
                                 rot = 90, size = 22, face = "bold"))

## Total abundance #####################################################################

### Across all sites (infested vs non-infested)
abund_tot <- abund_long1 %>%
  na.omit(abund_long) %>%
  dplyr::group_by(site, infestation, quadrat) %>%
  summarize(ab_total = sum(abundance))
View(abund_tot)

abund_tot1 <- abund_tot  %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    ab_tot_mean = mean(ab_total),
    ab_tot_sd = sd(ab_total)
  )
View(abund_tot1)

barplot_abund_tot2 <- ggplot(data = abund_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = ab_tot_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_sd, ymax = ab_tot_mean + ab_tot_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1050)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4.1),
           y = c(1050, 1050, 1050, 1050),
           label = c("A", "B", "B", "A"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total abundance")
barplot_abund_tot2

### Old Woman's River (infested vs non-infested, eastern vs western)
abund_tot_owr1 <- abund_tot_owr  %>%
  dplyr::group_by(infestation, lineage) %>%
  summarise(
    ab_tot_mean = mean(ab_total),
    ab_tot_sd = sd(ab_total)
  )
View(abund_tot_owr1)

barplot_abund_tot_owr1 <- ggplot(abund_tot_owr1, aes(x = infestation, y = ab_tot_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_sd, ymax = ab_tot_mean + ab_tot_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250),limits = c(0, 260)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total abundance")
barplot_abund_tot_owr1

## Total biomass ########################################################################

### Across all sites (infested vs non-infested)
biom_tot <- biom_long1 %>%
  na.omit(biom_long) %>%
  dplyr::group_by(site, infestation, quadrat) %>%
  summarize(biom_total = sum(biomass))
View(biom_tot)

biom_tot1 <- biom_tot  %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    biom_tot_mean = mean(biom_total),
    biom_mean_g = biom_tot_mean/1000,
    biom_tot_sd = sd(biom_total),
    biom_sd_g = biom_tot_sd/1000
  )
View(biom_tot1)

barplot_biom_tot2 <- ggplot(data = biom_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = biom_mean_g, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_mean_g - biom_sd_g, ymax = biom_mean_g + biom_sd_g), width = 0.2,
                position = position_dodge(0.9)) +
  annotate("text", 
           x = c(1, 2, 3.2, 4.1),
           y = c(20),
           label = c("A", "B", "AB", "A"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total biomass (in g)")
barplot_biom_tot2


### Old Woman's River (infested vs non-infested, eastern vs western)
biom_tot_owr1 <- biom_tot_owr  %>%
  dplyr::group_by(infestation, lineage) %>%
  summarise(
    biom_tot_mean = mean(biom_total),
    biom_mean_g = biom_tot_mean/1000,
    biom_tot_sd = sd(biom_total),
    biom_sd_g = biom_tot_sd/1000
  )
View(biom_tot_owr1)

barplot_biom_tot_owr2 <- ggplot(biom_tot_owr1, aes(x = infestation, y = biom_mean_g, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_mean_g - biom_sd_g, ymax = biom_mean_g + biom_sd_g), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total biomass (in g)")
barplot_biom_tot_owr2




## Species richness #####################################################################

### Across all sites (infested vs non-infested)
richness <- abund_long1 %>%
  dplyr::group_by(site, infestation, quadrat) %>%
  summarise(count = n_distinct(species))
View(richness)

richness1 <- richness %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
richness1

barplot_rich2 <- ggplot(data = richness1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = count_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 30), breaks = c(0, 10, 20, 30)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average species richness (S)")
barplot_rich2



### Old Woman's River (infested vs non-infested, eastern vs western)

barplot_rich_owr1 <- ggplot(richness_owr1, aes(x = infestation, y = count_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(limits = c(0, 40), breaks = c(0, 10, 20, 30, 40, 50)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average species richness (S)")
barplot_rich_owr1

## Species diversity - Shannon's Index H' ###############################################

### Across all sites (infested vs non-infested)
shannon1 <- abund_wide1 %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon)
  )
View(shannon1)

barplot_shannon2 <- ggplot(data = shannon1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = shannon_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(limits = c(0, 2.7), breaks = c(0, 0.5, 1, 1.5, 2, 2.5)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2.2, 3.2, 4.1),
           y = c(2.7),
           label = c("A", "AB", "AB", "B"), hjust = 1) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon2

### Old Woman's River (infested vs non-infested, eastern vs western)
barplot_shannon_owr1 <- ggplot(shannon_owr1, aes(x = infestation, y = shannon_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon_owr1

## Species diversity - Simpson's Index (λ) ###############################################

### Across all sites (infested vs non-infested)
simpson1 <- abund_wide1 %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson)
  )
View(simpson1)

barplot_simpson2 <- ggplot(data = simpson1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = simpson_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_sd, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1)) +
  annotate("text", 
           x = c(1, 2.2, 3.1, 4),
           y = c(1),
           label = c("A", "AB", "A", "B"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Simpson's diversity (λ) ")
barplot_simpson2

### Old Woman's River (infested vs non-infested, eastern vs western)
barplot_simpson_owr1 <- ggplot(simpson_owr1, aes(x = infestation, y = simpson_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_sd, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1, 1.5)) +
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
pielou1 <- abund_wide2 %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou)
  )
View(pielou1)

barplot_pielou2 <- ggplot(data = pielou1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = pielou_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1, 1)) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 1,
           label = c("A", "A", "A", "B")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Pielou's evenness (J)")
barplot_pielou2

### Old Woman's River (infested vs non-infested, eastern vs western)
barplot_pielou_owr1 <- ggplot(pielou_owr1, aes(x = infestation, y = pielou_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(limits = c(0, 1), breaks = c(0, 0.5, 1, 1.5)) +
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
points(abund_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 5)[env2$site],
       col = my_colors[env2$infestation])
ordihull(abund_bray, groups = env2$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.3, 0.3, legend = c(levels(env2$site), levels(env2$infestation)),
       pch = c(16, 17, 18, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-1.8, 1.35, "Stress = 0.123", bty = "n", cex = 1.5)

mtext("A", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(abund_disp_inf, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5),
     xlab = "Dimension 1 (31.82 %)", ylab = "Dimension 2 (23.93 %)",
     main = ""
)

mtext("B", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(abund_disp_site, hull = F, ellipse = T,
     pch = c(16:18, 8, 5), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (31.82 %)", ylab = "Dimension 2 (23.93 %)",
     xlim = c(-0.4, 0.4), ylim = c(-0.5, 0.5),
     main = "")

mtext("C", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

### Old Woman's River (infested vs non-infested, eastern vs western)

ordiplot(abund_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(abund_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env_owr1$lineage],
       col = my_colors[env_owr1$infestation])
ordihull(abund_bray_owr, groups = env_owr1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.35, -0.1, legend = c(levels(env_owr1$lineage), levels(env_owr1$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-1.8, 1.35, "Stress = 0.116", bty = "n", cex = 1.5)

mtext("D", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(abund_disp_inf_owr, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21 %)",
     main = ""
)

mtext("E", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

plot(abund_disp_lineage, hull = F, ellipse = T,
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

ordiplot(biom_bray, type = "none", xlim = c(-1.4, 1.4), ylim = c(-1.0, 1.0))
points(biom_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 5)[env2$site],
       col = my_colors[env2$infestation])
ordihull(biom_bray, groups = env2$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.55, 0.5, legend = c(levels(env2$site), levels(env2$infestation)),
       pch = c(16, 17, 18, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.15, 1.6, "Stress = 0.157", bty = "n", cex = 1.5)

mtext("A", side = 3, line = -1, cex = 2, adj = -0.1, col = "grey30")

plot(biom_disp_inf, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5),
     xlab = "Dimension 1 (36.50 %)", ylab = "Dimension 2 (21.77 %)",
     main = ""
)

mtext("B", side = 3, line = -1, cex = 2, adj = -0.1, col = "grey30")

plot(biom_disp_site, hull = F, ellipse = T,
     pch = c(16:18, 5), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (36.50 %)", ylab = "Dimension 2 (21.77 %)",
     xlim = c(-0.8, 0.8), ylim = c(-0.5, 0.5),
     main = "")

mtext("C", side = 3, line = -1, cex = 2, adj = -0.13, col = "grey30")

ordiplot(biom_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(biom_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env_owr1$lineage],
       col = my_colors[env_owr1$infestation])
ordihull(biom_bray_owr, groups = env_owr1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-1.35, -0.1, legend = c(levels(env_owr1$lineage), levels(env_owr1$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 0.8, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-1.85, 1.35, "Stress = 0.116", bty = "n", cex = 1.5)

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
levels(abund_wide1$site) <- c("BoS", "JB", "MB", "PEd")

abund_inf <- rev(levels(abund_wide1$infestation))
abund_dend <- as.dendrogram(abund_clust)

labels_colors(abund_dend) <- my_colors[(abund_wide1$infestation)[order.dendrogram(abund_dend)]]

labels(abund_dend) <- paste(abund_wide1$site[order.dendrogram(abund_dend)],
                            " (", labels(abund_dend), ")", sep = "")
branches_attr_by_lists(abund_dend, abund_wide1$site[order.dendrogram(abund_dend)],
                       attr = "lwd")
abund_dend <- set(abund_dend, "labels_cex", 1.5)
abund_dend <- hang.dendrogram(abund_dend, hang_height = 0.005)
plot(abund_dend, nodePar = list(cex = 0.007))
legend(35, 2, legend = levels(abund_wide1$infestation), fill = my_colors, bty = "n", cex = 1.5)

mtext("A", side = 3, line = -2, cex = 2, adj = -0.1, col = "grey30")

## Cluster analysis - all sites on biomass
biom_clust <- hclust(biom_dist, "ward.D2")
levels(env2$site) <- c("BoS", "JB", "MB", "PEd")

biom_inf <- rev(levels(env2$infestation))
biom_dend <- as.dendrogram(biom_clust)

labels_colors(biom_dend) <- my_colors[(env2$infestation)[order.dendrogram(biom_dend)]]

labels(biom_dend) <- paste(env2$site[order.dendrogram(biom_dend)],
                           " (", labels(biom_dend), ")", sep = "")
biom_dend <- set(biom_dend, "labels_cex", 1.5)
biom_dend <- hang.dendrogram(biom_dend, hang_height = 0.005)
plot(biom_dend, nodePar = list(cex = 0.007))
legend(35, 2, legend = levels(env2$infestation), fill = my_colors, bty = "n", cex = 1.5)

mtext("B", side = 3, line = -2, cex = 2, adj = -0.1, col = "grey30")

## Cluster analysis - OWR abundance
abund_clust_owr <- hclust(abund_dist_owr, "ward.D2")

abund_inf_owr <- rev(levels(env_owr1$infestation))
abund_dend_owr <- as.dendrogram(abund_clust_owr)

labels_colors(abund_dend_owr) <- my_colors[(env_owr1$infestation)[order.dendrogram(abund_dend_owr)]]

labels(abund_dend_owr) <- paste(env_owr1$lineage[order.dendrogram(abund_dend_owr)],
                                " (", labels(abund_dend_owr), ")", sep = "")
abund_dend_owr <- set(abund_dend_owr, "labels_cex", 1.5)
abund_dend_owr <- hang.dendrogram(abund_dend_owr, hang_height = 0.005)
plot(abund_dend_owr, nodePar = list(cex = 0.007))

mtext("C", side = 3, line = -1.5, cex = 2, adj = -0.1, col = "grey30")

## Cluster analysis - OWR biomass
biom_clust_owr <- hclust(biom_dist_owr, "ward.D2")

biom_inf_owr <- rev(levels(env_owr1$infestation))
biom_dend_owr <- as.dendrogram(biom_clust_owr)

labels_colors(biom_dend_owr) <- my_colors[(env_owr1$infestation)[order.dendrogram(biom_dend_owr)]]

labels(biom_dend_owr) <- paste(env_owr1$lineage[order.dendrogram(biom_dend_owr)],
                               " (", labels(biom_dend_owr), ")", sep = "")
biom_dend_owr <- set(biom_dend_owr, "labels_cex", 1.5)
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

abund_cont <- abund_long1[abund_long1$species %in% c("Talorchestia sp.",
                      "Actiniaria",
                      "Burnupena lagenaria",
                      "Mytilus galloprovincialis",
                      "Parisocladus perforatus",
                      "Perna perna",
                      "Parvulastra exigua",
                      "Pseudonereis variegata"),]
View(abund_cont)

abund_cont1 <- abund_cont %>%
    dplyr::group_by(site, infestation, species) %>%
  dplyr::summarise(
    mean_ab = mean(abundance),
    sd_ab = sd(abundance)
  )
View(abund_cont1)

barplot_abund_cont <- ggplot(abund_cont1, aes(x = site, y = mean_ab, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_ab - sd_ab, ymax = mean_ab + sd_ab), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_text(data = abund_cont_text, mapping = aes(x = x, y = y, label = label)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(limits = c(-5, 410), breaks = seq(from = 0, to = 400, by = 50)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ species, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average abundance")
barplot_abund_cont

abund_cont_text <- data.frame(
  label = c("B", "A", "B", "B", 
            "C", "A", "AB", "BC",
            "BC", "A", "AB", "C",
            "AB", "AC", "B", "C",
            "A", "B", "B", "B",
            "AB", "A", "A", "B",
            "AC", "C", "B", "AB",
            "A", "B", "BC", "AC"),
  species = c("Actiniaria", "Actiniaria", "Actiniaria", "Actiniaria",
              "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria",
              "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis",
              "Parisocladus perforatus", "Parisocladus perforatus", "Parisocladus perforatus", "Parisocladus perforatus",
              "Parvulastra exigua", "Parvulastra exigua", "Parvulastra exigua", "Parvulastra exigua",
              "Perna perna", "Perna perna", "Perna perna", "Perna perna",
              "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata",
              "Talorchestia sp.", "Talorchestia sp.", "Talorchestia sp.", "Talorchestia sp."),
  infestation = c("infested"),
  x = c(1, 2, 3, 4),
  y = c(20, 70, 30, 30, 
        10, 60, 30, 20,
        20, 240, 50, 10,
        70, 10, 80, 10,
        40, 20, 20, 20,
        90, 410, 320, 40,
        30, 50, 10, 30,
        10, 180, 50, 40)
)

## Biomass
View(biom_long)

biom_cont <- biom_long1[biom_long1$species %in% c("Actiniaria",
                                                "Burnupena lagenaria",
                                                "Mytilus galloprovincialis",
                                                "Ischyromene huttoni",
                                                "Perna perna",
                                                "Parvulastra exigua",
                                                "Pseudonereis variegata"),]
View(biom_cont)

biom_cont1 <- biom_cont %>%
  dplyr::group_by(site, infestation, species) %>%
  dplyr::summarise(
    mean_biom = mean(biomass),
    biom_g = mean_biom/1000,
    sd_biom = sd(biomass),
    sd_g = sd_biom/1000
  )
View(biom_cont1)

barplot_biom_cont <- ggplot(biom_cont1, aes(x = site, y = biom_g, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_g - sd_g, ymax = biom_g + sd_g), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  geom_text(data = biom_cont_text, mapping = aes(x = x, y = y, label = label)) +
  scale_y_continuous(limits = c(-3, 15), breaks = seq(from = 0, to = 15, by = 2.5)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ species, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average biomass (g)")
barplot_biom_cont

biom_cont_text <- data.frame(
  label = c("A", "B", "A", "A",
            "A", "B", "BC", "AC",
            "A", "ABC", "B", "C",
            "A", "B", "B", "A",
            "A", "B", "B", "B",
            "A", "AB", "AB", "B",
            "AC", "A", "B", "BC"),
  species = c("Actiniaria", "Actiniaria", "Actiniaria", "Actiniaria",
              "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria", "Burnupena lagenaria",
              "Ischyromene huttoni", "Ischyromene huttoni", "Ischyromene huttoni", "Ischyromene huttoni",
              "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis", "Mytilus galloprovincialis",
              "Parvulastra exigua", "Parvulastra exigua", "Parvulastra exigua", "Parvulastra exigua",
              "Perna perna", "Perna perna", "Perna perna", "Perna perna",
              "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata", "Pseudonereis variegata"
              ),
  infestation = c("infested"),
  x = c(1, 2, 3, 4),
  y = c(1, 2, 1, 1,
        0.5, 1.5, 1, 1,
        1, 2, 3, 1,
        1, 3.5, 1.5, 0.5,
        3, 1, 1, 1,
        2, 15, 5, 1,
        1, 1.5, 0.5, 1)
)
