## ---------------------------
##
## Script name: Infaunal communities associated with mussel beds depending on their level of
##              euendolithic infestation and Perna lineages
##
## Purpose of script: 
##    - Do infaunal communities differ between infested and non-infested mussel beds
##      across bioregions? 
##    - Do infaunal communities differ between infested and non-infested mussel beds and
##      eastern and western Perna lineages at Old Woman's River ?
##
## Author: Alexia Dievart
##
## Date Created: 2022-12-14
## Dates Updated: 2023-01-04, 2023-01-03
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
               car,
               MASS,
               emmeans,
               vegan)

# Set default ggplot theme
theme_set(theme_classic() +
            theme(panel.border = element_rect(colour = "black", fill = NA),
                  axis.text = element_text(colour = "black"),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  axis.title.x = element_text(margin = unit(c(2, 0, 0, 0), "mm")),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
                  legend.position = "none"))
my_colors <- c("darkgrey", "brown3")




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

## Create the contingency tables for community analyses
abund_wide1 <- abund_wide[,6:120]
abund_wide_owr <- abund_wide %>%
  filter(site == "Old Woman's River")
abund_wide1_owr <- abund_wide_owr[,6:120]
View(abund_wide1_owr)

biom_wide1 <- biom_wide[,6:120]
biom_wide_owr <- biom_wide %>%
  filter(site == "Old Woman's River")
biom_wide1_owr <- biom_wide_owr[,6:120]

## Clean the environmental dataset
View(env1)
env1 <- env[1:4]
env1$live_mussels <- env$live_perna + env$live_mytilus
env1$dead_mussels <- env$dead_perna + env$dead_mytilus
env1$nb_broken <- env$nb_broken
env1$mean_byssal <- env$mean_byssal
env1[,factor1] <- lapply(env1[,factor1], factor)

## Isolate Old Woman's River for analysis involving Perna lineage as a factor
byssus_owr <- byssus %>%
  filter(site == "Old Woman's River")

env_owr <- env %>%
  filter(site == "Old Woman's River")
env1_owr <- env1 %>%
  filter(site == "Old Woman's River")
View(env1_owr)



###############################################################################################
# Section: Explore data ----------------------------------------------------------------------
###############################################################################################

# Mussel bed architecture - statistical summary ===============================================

## Mean number of live mussels (Perna + Mytilus)
live_mussels <- env %>%
  dplyr::group_by(site, infestation) %>%
  dplyr::summarise(
    mean_live = mean(live_perna + live_mytilus),
    sd_live  = sd(live_perna + live_mytilus),
    n = n(),
    se_live  = sd_live / sqrt(n)
  )
head(live_mussels)

live_mussels_owr <- env %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage) %>%
  dplyr::summarise(
    mean_live = mean(live_perna + live_mytilus),
    sd_live  = sd(live_perna + live_mytilus),
    n = n(),
    se_live  = sd_live / sqrt(n)
  )
head(live_mussels_owr)

## Mean number of dead mussels (Perna + Mytilus)
dead_mussels <- env %>%
  dplyr::group_by(site, infestation) %>%
  dplyr::summarise(
    mean_dead = mean(dead_perna + dead_mytilus),
    sd_dead  = sd(dead_perna + dead_mytilus),
    n = n(),
    se_dead  = sd_dead / sqrt(n)
  )
head(dead_mussels)

dead_mussels_owr <- env %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage) %>%
  dplyr::summarise(
    mean_dead = mean(dead_perna + dead_mytilus),
    sd_dead  = sd(dead_perna + dead_mytilus),
    n = n(),
    se_dead  = sd_dead / sqrt(n)
  )
head(dead_mussels_owr)

## Mean number of broken mussels (Perna + Mytilus)
broken_mussels <- env %>%
  dplyr::group_by(site, infestation) %>%
  dplyr::summarise(
    mean_broken = mean(nb_broken),
    sd_broken  = sd(nb_broken),
    n = n(),
    se_broken  = sd_broken / sqrt(n)
  )
head(broken_mussels)

broken_mussels_owr <- env %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage) %>%
  dplyr::summarise(
    mean_broken = mean(nb_broken),
    sd_broken  = sd(nb_broken),
    n = n(),
    se_broken  = sd_broken / sqrt(n)
  )
head(broken_mussels_owr)




# Number of unique and common taxa per treatment =============================================

## Total number of unique species across all sites
abund_long %>% summarise(n_distinct(species)) 

## Number of unique species per infestation level across all sites
abund_long %>%
  group_by(infestation) %>%
  summarise(n_distinct(species))

## Number and identity of common species between infestation levels across all sites
Reduce(intersect, split(abund_long$species, abund_long$infestation))

## Number of unique species per infestation level and Perna lineage at Old Woman's River
abund_long %>%
  filter(site == "Old Woman's River") %>%
  group_by(infestation, lineage) %>%
  summarise(n_distinct(species))

## Number and identity of common species between infestation levels and Perna lineages at Old Woman's River
abund_owr <- abund_long %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage) %>%
  tidyr::unite(col = "id", infestation, lineage, sep = "_", remove = FALSE)
Reduce(intersect, split(abund_owr$species, abund_owr$id))




# Total abundance ===========================================================================
abund_tot <- abund_long %>%
  na.omit(abund_long) %>%
  dplyr::group_by(site, infestation, quadrat) %>%
  summarize(ab_total = sum(abundance))
View(abund_tot)

abund_tot1 <- abund_tot  %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    ab_tot_mean = mean(ab_total),
    ab_tot_sd = sd(ab_total),
    n = n(),
    ab_tot_se = ab_tot_mean / sqrt(n)
  )
View(abund_tot1)

abund_tot_owr <- abund_long %>%
  na.omit(abund_long) %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  summarize(ab_total = sum(abundance))
head(abund_tot_owr)

abund_tot_owr1 <- abund_tot_owr  %>%
  dplyr::group_by(infestation, lineage) %>%
  summarise(
    ab_tot_mean = mean(ab_total),
    ab_tot_sd = sd(ab_total),
    n = n(),
    ab_tot_se = ab_tot_mean / sqrt(n)
  )
View(abund_tot_owr1)





# Total biomass ===========================================================================
biom_tot <- biom_long %>%
  na.omit(biom_long) %>%
  dplyr::group_by(site, infestation, quadrat) %>%
  summarize(biom_total = sum(biomass))
View(biom_tot)

biom_tot1 <- biom_tot  %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    biom_tot_mean = mean(biom_total),
    biom_tot_sd = sd(biom_total),
    n = n(),
    biom_tot_se = biom_tot_mean / sqrt(n)
  )
View(biom_tot1)

biom_tot_owr <- biom_long %>%
  na.omit(biom_long) %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  summarize(biom_total = sum(biomass))
head(biom_tot_owr)

biom_tot_owr1 <- biom_tot_owr  %>%
  dplyr::group_by(infestation, lineage) %>%
  summarise(
    biom_tot_mean = mean(biom_total),
    biom_tot_sd = sd(biom_total),
    n = n(),
    biom_tot_se = biom_tot_mean / sqrt(n)
  )
View(biom_tot_owr1)



# Species richness ===========================================================================
richness <- abund_long %>%
  dplyr::group_by(site, infestation, quadrat) %>%
  summarise(count = n_distinct(species))
View(richness)

richness1 <- richness %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count),
    n = n(),
    count_se = count_mean / sqrt(n)
  )
richness1

richness_owr <- abund_long %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  summarise(count = n_distinct(species))
View(richness_owr)

richness_owr1 <- richness_owr %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count),
    n = n(),
    count_se = count_mean / sqrt(n)
  )
richness_owr1




# Species diversity - Shannon's Index (H') ========================================
library(vegan)
env$shannon <- diversity(abund_wide[,6:120], "shannon")
env1$shannon <- diversity(abund_wide[,6:120], "shannon")
View(env)

shannon <- env %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon),
    n = n(),
    shannon_se = shannon_mean / sqrt(n)
  )
View(shannon)

env_owr <- env %>%
  filter(site == "Old Woman's River")
View(env_owr)

shannon_owr <- env %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon),
    n = n(),
    shannon_se = shannon_mean / sqrt(n)
  )
View(shannon_owr)




# Species diversity - Simpson's Index (λ) ========================================
env$simpson <- diversity(abund_wide[,6:120], index="simpson")
env1$simpson <- diversity(abund_wide[,6:120], index="simpson")
View(env)

simpson <- env %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson),
    n = n(),
    simpson_se = simpson_mean / sqrt(n)
  )
View(simpson)

env_owr <- env %>%
  filter(site == "Old Woman's River")
View(env_owr)

simpson_owr <- env %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson),
    n = n(),
    simpson_se = simpson_mean / sqrt(n)
  )
View(simpson_owr)




# Species evenness - Pielou's Index (J) ========================================

# J = H'/ln(S)

richness2 <- abund_long %>%
  dplyr::group_by(site, infestation, lineage, quadrat) %>%
  summarise(count = n_distinct(species))
View(richness2)

env$richness <- richness2$count
env1$richness <- richness2$count

env$pielou <- env$shannon/log(env$richness)
env1$pielou <- env1$shannon/log(env1$richness)
View(env)

pielou <- env %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou),
    n = n(),
    pielou_se = pielou_mean / sqrt(n)
  )
View(pielou)

env_owr <- env %>%
  filter(site == "Old Woman's River")
View(env_owr)

pielou_owr <- env %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou),
    n = n(),
    pielou_se = pielou_mean / sqrt(n)
  )
View(pielou_owr)



###############################################################################################
# Section: Visualize data ---------------------------------------------------------------------
###############################################################################################

# Calculate and plot interesting variables for each infestation level and genetic lineages, 
# such as:
#   - Architectural complexity, with nb of live and dead mussels and broken shells
#   - Architectural complexity, with nb of byssal threads as a proxy
#   - Total abundance
#   - Species richness (S) - count = n_distinct(species)
#   - Species diversity - Shannon's Index (H')



# Mussel bed architecture - mean number of live mussels per quadrat ===================================

## Barplot of the average number of live mussels for each quadrat, site and infestation level
barplot_livemussels <- ggplot(env, aes(x = quadrat, y = (live_perna + live_mytilus), 
                                       color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Number of live mussels")
barplot_livemussels

barplot_livemussels1 <- ggplot(data = live_mussels, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = mean_live, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_live - se_live, ymax = mean_live + se_live), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 44,
           label = c("AB", "AB", "A", "AB", "B")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of live mussels")
barplot_livemussels1

## Identify outliers and extreme values for each quadrat, site and infestation level
env %>%
  group_by(site, infestation) %>%
  identify_outliers(live_perna)
### Identified extreme values were omitted in the analysis (i.e, OWR EI1)

## Barplot of the average number of live mussels per quadrat, infestation level and lineage at Old Woman's River
barplot_live_owr <- ggplot(env_owr, aes(x = quadrat, y = (live_perna + live_mytilus), 
                                       color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Number of live mussels")
barplot_live_owr

barplot_live_owr1 <- ggplot(live_mussels_owr, aes(x = infestation, y = mean_live, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_live - se_live, ymax = mean_live + se_live), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of live mussels")
barplot_live_owr1




# Mussel bed architecture - mean number of dead mussels per quadrat ===================================

## Barplot of the average number of dead mussels for each quadrat, site and infestation level
barplot_deadmussels <- ggplot(env, aes(x = quadrat, y = (dead_perna + dead_mytilus), 
                                       color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Number of dead mussels")
barplot_deadmussels

barplot_deadmussels1 <- ggplot(data = dead_mussels, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = mean_dead, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_dead - se_dead, ymax = mean_dead + se_dead), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 7.2,
           label = c("A", "B", "C", "AB", "C")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of dead mussels")
barplot_deadmussels1

## Identify outliers and extreme values for each quadrat, site and infestation level
env %>%
  group_by(site, infestation) %>%
  identify_outliers(dead_perna)
### Identified extreme values were omitted in the analysis (i.e, OWR WNI1)

## Barplot of the average number of dead mussels per quadrat, infestation level and lineage at Old Woman's River
barplot_dead_owr <- ggplot(env_owr, aes(x = quadrat, y = (dead_perna + dead_mytilus), 
                                        color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Number of dead mussels")
barplot_dead_owr

barplot_dead_owr1 <- ggplot(dead_mussels_owr, aes(x = infestation, y = mean_dead, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_dead - se_dead, ymax = mean_dead + se_dead), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of dead mussels")
barplot_dead_owr1




# Mussel bed architecture - mean number, sd and se of broken shells per quadrat ==============

## Barplot of the average number of broken shells for each quadrat, site and infestation level
barplot_brokenshells <- ggplot(env, aes(x = quadrat, y = nb_broken, 
                                       color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Number of broken shells")
barplot_brokenshells

barplot_brokenshells1 <- ggplot(data = broken_mussels, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = mean_broken, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_broken - se_broken, ymax = mean_broken + se_broken), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 30,
           label = c("A", "A", "A", "B", "A")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average number of broken shells")
barplot_brokenshells1

## Identify outliers and extreme values for each quadrat, site and infestation level
env %>%
  group_by(site, infestation) %>%
  identify_outliers(nb_broken)
### No identified outliers.

## Barplot of the average number of dead mussels per quadrat, infestation level and lineage at Old Woman's River
barplot_broken_owr <- ggplot(env_owr, aes(x = quadrat, y = nb_broken, 
                                        color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Number of broken shells")
barplot_broken_owr

barplot_broken_owr1 <- ggplot(broken_mussels_owr, aes(x = infestation, y = mean_broken, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_broken - se_broken, ymax = mean_broken + se_broken), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of broken shells")
barplot_broken_owr1




# Mussel bed architecture - mean number, sd and se of byssal threads per quadrat ==============

## Boxplot of the average number of byssal threads for each quadrat, site and infestation level
boxplot_bt <- ggplot(byssus, aes(x = quadrat, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Quadrat", y = "Average number of byssal threads")
boxplot_bt

boxplot_bt1 <- ggplot(byssus, aes(x = site, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 375,
           label = c("B", "A", "AC", "B", "BC")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average number of byssal threads")
boxplot_bt1

## Identify outliers and extreme values for each quadrat, site and infestation level
byssus %>%
  group_by(site, infestation, quadrat) %>%
  identify_outliers(nb_byssal)
### Identified extreme values were omitted in the analysis.

## Boxplot of the average number of byssal threads per quadrat, infestation level and Perna lineage at Old Woman's River
boxplot_bt_owr <- ggplot(byssus_owr, aes(x = quadrat, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Quadrat", y = "Average number of byssal threads")
boxplot_bt_owr

boxplot_bt_owr1 <- ggplot(byssus_owr, aes(x = infestation, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average number of byssal threads")
boxplot_bt_owr1




# Total abundance in infauna per quadrat ========================================================

## Barplot of the total abundance in infauna for each quadrat, site and infestation level
barplot_abund_tot <- ggplot(abund_tot, aes(x = quadrat, y = ab_total, 
                                        color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Total abundance")
barplot_abund_tot

barplot_abund_tot1 <- ggplot(data = abund_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = ab_tot_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_se, ymax = ab_tot_mean + ab_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  annotate("text", 
           x = c(0.85, 1.25, 1.85, 2.3, 2.8, 3.25, 3.8, 4.25, 4.8, 5.25),
           y = c(700, 1050, 700, 580, 300, 300, 500, 460, 220, 170),
           label = c("AB", "A", "AB", "AB", "B", "B", "B", "B", "B", "B"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total abundance")
barplot_abund_tot1

## Identify outliers and extreme values
abund_tot %>%
  group_by(site, infestation, quadrat) %>%
  identify_outliers(ab_total)
### No identified outliers.

## Barplot of the total abundance in infauna per quadrat, infestation level and lineage at Old Woman's River
barplot_abund_tot_owr <- ggplot(abund_tot_owr, aes(x = quadrat, y = ab_total, 
                                          color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Total abundance")
barplot_abund_tot_owr

barplot_abund_tot_owr1 <- ggplot(abund_tot_owr1, aes(x = infestation, y = ab_tot_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_se, ymax = ab_tot_mean + ab_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total abundance")
barplot_abund_tot_owr1




# Total biomass in infauna per quadrat ========================================================

## Barplot of the total abundance in infauna for each quadrat, site and infestation level
barplot_biom_tot <- ggplot(biom_tot, aes(x = quadrat, y = biom_total, 
                                           color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Total biomass")
barplot_biom_tot

barplot_biom_tot1 <- ggplot(data = biom_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = biom_tot_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_tot_mean - biom_tot_se, ymax = biom_tot_mean + biom_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  annotate("text", 
           x = c(0.85, 1.25, 1.85, 2.3, 2.8, 3.25, 3.8, 4.25, 4.8, 5.25),
           y = c(16000, 13000, 7000, 6000, 5000, 5100, 9000, 8500, 3000, 3000),
           label = c("A", "A", "AB", "AB", "AB", "AB", "AC", "AC", "BC", "B"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total biomass")
barplot_biom_tot1

## Identify outliers and extreme values
biom_tot %>%
  group_by(site, infestation, quadrat) %>%
  identify_outliers(biom_total)
### No identified outliers.

## Barplot of the total abundance in infauna per quadrat, infestation level and lineage at Old Woman's River
barplot_biom_tot_owr <- ggplot(biom_tot_owr, aes(x = quadrat, y = biom_total, 
                                                   color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Total biomass")
barplot_biom_tot_owr

barplot_biom_tot_owr1 <- ggplot(biom_tot_owr1, aes(x = infestation, y = biom_tot_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_tot_mean - biom_tot_se, ymax = biom_tot_mean + biom_tot_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total biomass")
barplot_biom_tot_owr1



# Species richness per quadrat ========================================================

## Barplot of the species richness for each quadrat, site and infestation level
barplot_rich <- ggplot(richness, aes(x = quadrat, y = count, 
                                           color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Species richness")
barplot_rich

barplot_rich1 <- ggplot(data = richness1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = count_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_se, ymax = count_mean + count_se), width = 0.2,
                position = position_dodge(0.9)) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 55,
           label = c("AB", "AB", "A", "B", "AB")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average species richness")
barplot_rich1

## Barplot of the species richness per quadrat, infestation level and lineage at Old Woman's River
barplot_rich_owr <- ggplot(richness_owr, aes(x = quadrat, y = count, 
                                                   color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Species richness")
barplot_rich_owr

barplot_rich_owr1 <- ggplot(richness_owr1, aes(x = infestation, y = count_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_se, ymax = count_mean + count_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average species richness")
barplot_rich_owr1




# Species diversity - Shannon's Index (H') ========================================

## Barplot of the species diversity for each quadrat, site and infestation level
barplot_shannon <- ggplot(env, aes(x = quadrat, y = shannon, 
                                     color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Shannon-Wiener diversity (H')")
barplot_shannon

barplot_shannon1 <- ggplot(data = shannon, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = shannon_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_se, ymax = shannon_mean + shannon_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon1

## Barplot of the species diversity per quadrat, infestation level and lineage at Old Woman's River
barplot_shannon_owr <- ggplot(env_owr, aes(x = quadrat, y = shannon, 
                                             color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Shannon-Wiener diversity (H')")
barplot_shannon_owr

barplot_shannon_owr1 <- ggplot(shannon_owr, aes(x = infestation, y = shannon_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_se, ymax = shannon_mean + shannon_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon_owr1




# Species diversity - Simpson's Index (λ)  ========================================

## Barplot of the species diversity for each quadrat, site and infestation level
barplot_simpson <- ggplot(env, aes(x = quadrat, y = simpson, 
                                   color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Simpson's diversity (λ) ")
barplot_simpson

barplot_simpson1 <- ggplot(data = simpson, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = simpson_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_se, ymax = simpson_mean + simpson_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Simpson's diversity (λ) ")
barplot_simpson1

## Barplot of the species diversity per quadrat, infestation level and lineage at Old Woman's River
barplot_simpson_owr <- ggplot(env_owr, aes(x = quadrat, y = simpson, 
                                           color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Simpson's diversity (λ) ")
barplot_simpson_owr

barplot_simpson_owr1 <- ggplot(simpson_owr, aes(x = infestation, y = simpson_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_se, ymax = simpson_mean + simpson_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Simpson's diversity (λ)")
barplot_simpson_owr1




# Species evenness - Pielou's Index (J)  ============================================

## Barplot of the species evenness for each quadrat, site and infestation level
barplot_pielou <- ggplot(env, aes(x = quadrat, y = pielou, 
                                   color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou

barplot_pielou1 <- ggplot(data = pielou, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = pielou_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_se, ymax = pielou_mean + pielou_se), width = 0.2,
                position = position_dodge(0.9)) +
  annotate("text", 
           x = c(1, 2, 3, 4, 5),
           y = 1.4,
           label = c("AB", "A", "AB", "AB", "B")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Pielou's evenness (J)")
barplot_pielou1

## Barplot of the species evenness per quadrat, infestation level and lineage at Old Woman's River
barplot_pielou_owr <- ggplot(env_owr, aes(x = quadrat, y = pielou, 
                                           color = infestation, fill = infestation)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou_owr

barplot_pielou_owr1 <- ggplot(pielou_owr, aes(x = infestation, y = pielou_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_se, ymax = pielou_mean + pielou_se), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Pielou's evenness (J)")
barplot_pielou_owr1





###################################################################################
# GRAPHICAL ANALYSIS - SUMMARY ----------------------------------------------------
###################################################################################

# Mussel bed architecture - mean number, sd and se of live mussels per quadrat ==============
##  - Few live mussels in Port Edward (i.e., PEd-NI-2, PEd-I-3 and PEd-NI-4)
##  - Few live mussels in eastern lineage (i.e., OWR-E-I-1 considered extreme)

# Mussel bed architecture - mean number, sd and se of dead mussels per quadrat ==============
##  - No dead mussels in Mosselbaai and Port Edward but high proportions of dead mussels in Jeffreysbaai
##    and Olw Woman's River (i.e., OWR-W-NI-1 considered extreme)
##  - More dead mussels in western lineage at Old Woman's River

# Mussel bed architecture - mean number, sd and se of byssal threads per quadrat ==============
##  - No visible difference between infestation levels across all sites
##  - No visible difference between infestation levels and Perna lineages at Old Woman's River





###############################################################################################
# Section: Statistical analyses ---------------------------------------------------------------
###############################################################################################

# Mussel bed architecture - mean number, sd and se of live mussels per quadrat ===============

## Infested vs Non-infested across all sites ##################################################

### Omit extreme values for statistical analysis
env_no <- env[!(row.names(env) %in% c("18")),]
View(env_no)

### Assumptions of normality and homogeneity of variances
ggqqplot(env_no, "live_perna", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_livemussels <- env_no %>%
  group_by(infestation, site) %>%
  shapiro_test(live_perna)
View(shap_livemussels)
#### Data is normally distributed

env_no %>% levene_test(live_perna ~ site * infestation)
#### Variance is not homogeneous

### ANOVA (using aov) with unbalanced design
live_mussels_aov <- aov(live_perna ~ site * infestation, data = env_no)
Anova(live_mussels_aov, type = "III")
#### There is a no significant effect of site or the infestation or the interaction on
#### the average number of live mussels.

live_mussels_aov1 <- aov(live_perna ~ site + infestation, data = env_no)
Anova(live_mussels_aov1, type = "II")
#### There is a significant effect of site on the average number of live mussels,
#### but no significant effect of infestation or the interaction.
#### ANOVA Type II is more powerful when interactions are not significant.

### Check ANOVA assumptions
plot(live_mussels_aov1, 1) # Homogeneity of variances
plot(live_mussels_aov1, 2) # Normality
#### Points 3, 4 and 28 were detected as outliers for homogeneity of variances
#### Points 2, 3 and 4 were detected as outliers for normality

### Post-hoc tests (using Tukey)
summary(glht(live_mussels_aov1, linfct = mcp(site = "Tukey")))
#### By comparing the average number of byssal threads between sites :
####    - Nb of live mussels in MB > PEd (and to a lesser extent to JB and OWR)

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Identify potential outliers
env_owr %>%
  group_by(infestation, lineage, quadrat) %>%
  identify_outliers(live_perna)

### Assumptions of normality and homogeneity of variances
ggqqplot(env_owr, "live_perna", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_livemussels_owr <- env_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(live_perna)
View(shap_livemussels_owr)
#### Data is normally distributed

env_owr %>% levene_test(live_perna ~ lineage * infestation)
#### Variance is homogeneous

### ANOVA (using aov) with unbalanced design
livemussels_owr_aov <- aov(live_perna ~ lineage * infestation, data = env_owr)
Anova(livemussels_owr_aov, type = "III")
#### There is no significant effect of the lineage or the infestation level on the number of byssal
#### threads, as well as no significant effect of the interaction.

livemussels_owr_aov1 <- aov(live_perna ~ lineage + infestation, data = env_owr)
Anova(livemussels_owr_aov1, type = "II")
#### There is a significant effect of the infestation level on the number of live mussels per quadrat,
#### but no significant effect of the lineage. 

### Check ANOVA assumptions
plot(livemussels_owr_aov1, 1) # Variances are more or less homogeneous
plot(livemussels_owr_aov1, 2) # Data is normally distributed
#### Points 2, 4 and 12 were detected as outliers for both homogeneity of variances and normality

### Post-hoc tests (using Tukey)
summary(glht(livemussels_owr_aov1, linfct = mcp(infestation = "Tukey")))
TukeyHSD((livemussels_owr_aov1))
emmeans(livemussels_owr_aov1, pairwise ~ lineage * infestation)
#### By comparing the average number of live mussels between infestation levels :
####    - Nb of live mussels in non-infested > infested, but p-value very close to 0.05




# Mussel bed architecture - mean number, sd and se of dead mussels per quadrat ===============

## Infested vs Non-infested across all sites ##################################################

### Omit extreme values for statistical analysis
env_no <- env[!(row.names(env) %in% c("23")),]
View(env_no)

### Assumptions of normality and homogeneity of variances
ggqqplot(env_no, "dead_perna", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

hist(env$dead_perna)
#### Data is normally distributed but follow a poisson/gamma distribution

env_no %>% levene_test(dead_perna ~ site * infestation)
#### Variance is not homogeneous

### ANOVA (using models)
deadmussels_lm <- lm(dead_perna ~ site * infestation, data = env_no)
summary(deadmussels_lm)
summary(deadmussels_lm)$r.squared
#### Linear model only explains 52 % of the variance.
DHARMa::simulateResiduals(fittedModel = deadmussels_lm, plot = T)
#### LM is not a good fit.

deadmussels_glm <- glm(dead_perna ~  site * infestation, data = env_no,
                      family = "poisson")
summary(deadmussels_glm)
DHARMa::simulateResiduals(fittedModel = deadmussels_glm, plot = T)
#### GLM with poisson is not a better fit, but quantile deviations detected.
#### Gamma distribution does not allow for the presence of null values. 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = deadmussels_glm, plot = F) 
plotResiduals(simulationOutput, env_no$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, env_no$site) # within-group deviations not significant

deadmussels_glm1 <- glm(dead_perna ~  site + infestation, data = env_no,
                       family = "poisson")
summary(deadmussels_glm1)
DHARMa::simulateResiduals(fittedModel = deadmussels_glm1, plot = T)
#### GLM with poisson without interaction is a better fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = deadmussels_glm1, plot = F) 
plotResiduals(simulationOutput, env_no$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, env_no$site) # within-group deviations not significant

### Post-hoc tests
summary(glht(deadmussels_glm1, linfct = mcp(site = "Tukey")))
#### By comparing the average number of dead mussels between sites :
####    - Nb of dead mussels in JB > BoS

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Identify potential outliers
env_owr %>%
  group_by(infestation, lineage, quadrat) %>%
  identify_outliers(dead_perna)

### Assumptions of normality and homogeneity of variances
ggqqplot(env_owr, "dead_perna", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_deadmussels_owr <- env_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(dead_perna)
View(shap_deadmussels_owr)
#### Data is normally distributed (except for E-I)

env_owr %>% levene_test(dead_perna ~ lineage * infestation)
#### Variance is homogeneous

### ANOVA (using aov) with unbalanced design
deadmussels_owr_aov <- aov(dead_perna ~ lineage * infestation, data = env_owr)
Anova(deadmussels_owr_aov, type = "III")
#### There is no significant effect of the lineage or the infestation level on the number of dead
#### mussels, as well as no significant effect of the interaction.

deadmussels_owr_aov1 <- aov(dead_perna ~ lineage + infestation, data = env_owr)
Anova(deadmussels_owr_aov, type = "II")
#### There is no significant effect of the lineage or the infestation level on the number of dead
#### mussels, as well as no significant effect of the interaction.

### Check ANOVA assumptions
plot(deadmussels_owr_aov1, 1) # Variances are not homogeneous
plot(deadmussels_owr_aov1, 2) # Data is not normally distributed
#### Points 1, 6 and 9 were detected as outliers for both homogeneity of variances and normality

### ANOVA (using models)
hist(env_owr$dead_perna)

deadmussels_owr_lm <- lm(dead_perna ~ lineage * infestation, data = env_owr)
summary(deadmussels_owr_lm)
summary(deadmussels_owr_lm)$r.squared
#### Linear model only explains 28 % of the variance.
DHARMa::simulateResiduals(fittedModel = deadmussels_owr_lm, plot = T)
#### LM is not a bad fit.

deadmussels_owr_glm <- glm(dead_perna ~  lineage * infestation, data = env_owr,
                      family = "poisson")
summary(deadmussels_owr_glm)
DHARMa::simulateResiduals(fittedModel = deadmussels_owr_glm, plot = T)
#### GLM with poisson is not a bad fit and has a better AIC.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = deadmussels_owr_glm, plot = F) 
plotResiduals(simulationOutput, env_owr$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, env_owr$lineage) # within-group deviations not significant




# Mussel bed architecture - mean number, sd and se of broken shells per quadrat ==============

## Infested vs Non-infested across all sites ##################################################

### Assumptions of normality and homogeneity of variances
ggqqplot(env, "nb_broken", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_broken <- env %>%
  group_by(infestation, site) %>%
  shapiro_test(nb_broken)
View(shap_broken)
#### Data is normally distributed

env %>% levene_test(nb_broken ~ site * infestation)
#### Variance is not homogeneous(close to p = 0.05)

### ANOVA (using aov) with unbalanced design
broken_aov <- aov(nb_broken ~ site * infestation, data = env)
Anova(broken_aov, type = "III")
#### There is no significant effect of site or infestation or interaction on the number of broken shells
#### per quadrat across all sites.

broken_aov1 <- aov(nb_broken ~ site + infestation, data = env)
Anova(broken_aov1, type = "II")

### Check ANOVA assumptions
plot(broken_aov1, 1) # Variances are not homogeneous
plot(broken_aov1, 2) # Data is not normally distributed
#### Points 18, 21 and 30 were detected as outliers for both homogeneity of variances and normality

#### Basic ANOVA is not the best for the analysis on broken shells but this will do for now.

### Post-hoc tests (using Tukey)
summary(glht(broken_aov1, linfct = mcp(site = "Tukey")))
TukeyHSD((broken_aov1))
#### By comparing the average number of broken shells per quadrat:
####  - Nb of broken shells in OWR > all sites.

### ANOVA (using models) - I tried several models but I can't seem to find one that's working for me
hist(env$nb_broken)

broken_lm <- lm(nb_broken ~ site * infestation, data = env)
summary(broken_lm)
summary(broken_lm)$r.squared
#### Linear model only explains 50 % of the variance.
DHARMa::simulateResiduals(fittedModel = broken_lm, plot = T)
#### LM is not a good fit.

broken_glm <- glm(nb_broken ~  site * infestation, data = env,
                       family = "poisson")
summary(broken_glm)
DHARMa::simulateResiduals(fittedModel = broken_glm, plot = T)
#### GLM with poisson is not a better fit.
#### Gamma distribution does not allow for the presence of null values. 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = broken_glm, plot = F) 
plotResiduals(simulationOutput, env$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, env$site) # within-group deviations not significant

broken_glm1 <- glm(nb_broken ~  site + infestation, data = env,
                  family = "poisson")
summary(broken_glm1)
DHARMa::simulateResiduals(fittedModel = broken_glm1, plot = T)
#### GLM with poisson and without the interaction is not a better fit.
#### Gamma distribution does not allow for the presence of null values. 

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Identify potential outliers
env_owr %>%
  group_by(infestation, lineage, quadrat) %>%
  identify_outliers(nb_broken)

### Assumptions of normality and homogeneity of variances
ggqqplot(env_owr, "nb_broken", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_broken_owr <- env_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(nb_broken)
View(shap_broken_owr)
#### Data is normally distributed (except for NI - E)

env_owr %>% levene_test(nb_broken ~ lineage * infestation)
#### Variance is homogeneous

### ANOVA (using aov) with unbalanced design
broken_owr_aov <- aov(nb_broken ~ lineage * infestation, data = env_owr)
Anova(broken_owr_aov, type = "III")
#### There is no significant effect of the lineage or the infestation level on the number of byssal
#### threads, as well as no significant effect of the interaction.

### Check ANOVA assumptions
plot(broken_owr_aov, 1) # Variances are barely homogeneous
plot(broken_owr_aov, 2) # Data is not normally distributed
#### Points 2, 10 and 16 were detected as outliers for both homogeneity of variances and normality

#### Let's use the ANOVA here

### Post-hoc tests (using Tukey)
TukeyHSD((broken_owr_aov))

### ANOVA (using models)
broken_owr_lm <- lm(nb_broken ~ lineage * infestation, data = env_owr)
summary(broken_owr_lm)
summary(broken_owr_lm)$r.squared
#### Linear model only explains 21 % of the variance.
DHARMa::simulateResiduals(fittedModel = broken_owr_lm, plot = T)
#### LM is not such a bad fit after all.
AIC(broken_owr_lm)
simulationOutput <- DHARMa::simulateResiduals(fittedModel = broken_owr_lm, plot = F) 
plotResiduals(simulationOutput, env_owr$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, env_owr$lineage) # within-group deviations not significant


broken_owr_glm1 <- glm(nb_broken ~  lineage * infestation, data = env_owr,
                      family = "poisson")
summary(broken_owr_glm1)
DHARMa::simulateResiduals(fittedModel = broken_owr_glm1, plot = T)
#### GLM with poisson is not a good fit.




# Mussel bed architecture - mean number, sd and se of byssal threads per quadrat ==============

## Infested vs Non-infested across all sites ##################################################

### Omit extreme values for statistical analysis
byssus_no <- byssus[!(row.names(byssus) %in% c("45","84")),]
View(byssus_no)

### Assumptions of normality and homogeneity of variances
ggqqplot(byssus_no, "nb_byssal", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_byssus <- byssus_no %>%
  group_by(infestation, site) %>%
  shapiro_test(nb_byssal)
View(shap_byssus)
#### Data is normally distributed

byssus_no %>% levene_test(nb_byssal ~ site * infestation)
#### Variance is homogeneous

### ANOVA (using aov) with unbalanced design
byssus_aov <- aov(nb_byssal ~ site * infestation, data = byssus_no)
Anova(byssus_aov, type = "III")
#### There is a significant effect of site on the average number of byssal threads (p < 0.05 **),
#### but no significant effect of infestation or the interaction.

byssus_aov1 <- aov(nb_byssal ~ site + infestation, data = byssus_no)
Anova(byssus_aov1, type = "II")

### Check ANOVA assumptions
plot(byssus_aov1, 1) # Homogeneity of variances
plot(byssus_aov1, 2) # Normality
#### Points 22, 130 and 26 were detected as outliers for both homogeneity of variances and normality

### Post-hoc tests (using Tukey)
summary(glht(byssus_aov1, linfct = mcp(site = "Tukey")))
TukeyHSD((byssus_aov1))
#### By comparing the average number of byssal threads between sites :
####    - Nb of byssal threads in JB > BoS, OWR and PEd
####    - Nb of byssal threads in MB > BoS, OWR

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Identify potential outliers
byssus_owr %>%
  group_by(infestation, lineage, quadrat) %>%
  identify_outliers(nb_byssal)
View(byssus_owr)

### Assumptions of normality and homogeneity of variances
ggqqplot(byssus_owr, "nb_byssal", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_byssus_owr <- byssus_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(nb_byssal)
View(shap_byssus_owr)
#### Data is not normally distributed for W[infested] and E[non infested]

byssus_owr %>% levene_test(nb_byssal ~ lineage * infestation)
#### Variance is homogeneous

### ANOVA (using aov) with unbalanced design
byssus_owr_aov <- aov(nb_byssal ~ lineage * infestation, data = byssus_owr)
Anova(byssus_owr_aov, type = "III")
#### There is no significant effect of the lineage or the infestation level on the number of byssal
#### threads, as well as no significant effect of the interaction.

### Check ANOVA assumptions
plot(byssus_owr_aov, 1) # Variances are not homogeneous
plot(byssus_owr_aov, 2) # Data is not normally distributed

### ANOVA (using models)
byssus_owr_lm <- lm(nb_byssal ~ lineage * infestation, data = byssus_owr)
summary(byssus_owr_lm)
summary(byssus_owr_lm)$r.squared
#### Linear model only explains 5 % of the variance.
DHARMa::simulateResiduals(fittedModel = byssus_owr_lm, plot = T)

byssus_owr_glm <- glm(nb_byssal ~  lineage * infestation, data = byssus_owr,
                      family = "poisson")
summary(byssus_owr_glm)
DHARMa::simulateResiduals(fittedModel = byssus_owr_glm, plot = T)
#### GLM with poisson is not a good fit.

byssus_owr_glm1 <- glm(nb_byssal ~  lineage * infestation, data = byssus_owr,
                      family = "Gamma")
summary(byssus_owr_glm1)
DHARMa::simulateResiduals(fittedModel = byssus_owr_glm1, plot = T)
### GLM with gamma is a good fit and has a better AIC.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = byssus_owr_glm1, plot = F) 
plotResiduals(simulationOutput, byssus_owr$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, byssus_owr$lineage) # within-group deviations not significant




# Total abundance in infauna ==================================================================

## Infested vs Non-infested across all sites ##################################################

## Assumptions of normality and homogeneity of variances
ggqqplot(abund_tot, "ab_total", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_abund_tot <- abund_tot %>%
  group_by(infestation, site) %>%
  shapiro_test(ab_total)
View(shap_abund_tot)
#### Data is normally distributed (except for Mosselbaai)

levene_test(abund_tot$ab_total ~ abund_tot$site * abund_tot$infestation)
#### Strangely, the test is not working when I don't specify the columns. 
#### Variance is barely homogeneous.

### ANOVA (using aov) with unbalanced design
abund_tot_aov <- aov(ab_total ~ site * infestation, data = abund_tot)
Anova(abund_tot_aov, type = "III")
#### There is a significant effect of site and infestation on the total abundance per quadrat,
#### but no effect of the interaction

abund_tot_aov1 <- aov(ab_total ~ site + infestation, data = abund_tot)
Anova(abund_tot_aov1, type = "II")
#### There is a significant effect of site on the total abundance per quadrat,
#### but no significant effect of infestation. 
#### ANOVA Type II is more powerful when interactions are not significant.

### Check ANOVA assumptions
plot(abund_tot_aov1, 1) # Variances are not homogeneous
plot(abund_tot_aov1, 2) # Data is not normaly distributed
#### Points 4, 6 and 7 were detected as outliers for normality and homogeneity of variances

### ANOVA (using models)
hist(abund_tot$ab_total)

abund_tot_lm <- lm(ab_total ~ site * infestation, data = abund_tot)
summary(abund_tot_lm)
summary(abund_tot_lm)$r.squared
#### Linear model only explains 70 % of the variance.
DHARMa::simulateResiduals(fittedModel = abund_tot_lm, plot = T)
#### LM is not a bad fit. 
AIC(abund_tot_lm)
Anova(abund_tot_lm)

simulationOutput <- DHARMa::simulateResiduals(fittedModel = abund_tot_lm, plot = F) 
plotResiduals(simulationOutput, abund_tot$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, abund_tot$site) # within-group deviations not significant

### Post-hoc tests (using Tukey)
emmeans(abund_tot_lm, pairwise ~ site * infestation)
#### By comparing the total abundance per site and infestation level:
####    - Total abundance in BoS NI > MB, OWR, PEd

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Identify potential outliers
abund_tot_owr %>%
  group_by(infestation, lineage, quadrat) %>%
  identify_outliers(ab_total)

### Assumptions of normality and homogeneity of variances
ggqqplot(abund_tot_owr, "ab_total", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_abund_tot_owr <- abund_tot_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(ab_total)
View(shap_abund_tot_owr)
#### Data is normally distributed

levene_test(abund_tot_owr$ab_total ~ abund_tot_owr$lineage * abund_tot_owr$infestation)
#### Strangely, this is not working again.

### ANOVA (using aov) with unbalanced design
abund_tot_owr_aov <- aov(ab_total ~ lineage * infestation, data = abund_tot_owr)
Anova(abund_tot_owr_aov, type = "III")
#### There is no significant effect of lineage or infestation or interaction on the total abundance.

abund_tot_owr_aov1 <- aov(ab_total ~ lineage + infestation, data = abund_tot_owr)
Anova(abund_tot_owr_aov1, type = "II")
#### There is no significant effect of lineage or infestation on total abundance.
#### ANOVA Type II is more powerful when interactions are not significant.

### Check ANOVA assumptions
plot(abund_tot_owr_aov1, 1) # Variances are not homogeneous
plot(abund_tot_owr_aov1, 2) # Normality

### ANOVA (using models)
abund_tot_owr_lm <- lm(ab_total ~ lineage * infestation, data = abund_tot_owr)
summary(abund_tot_owr_lm)
summary(abund_tot_owr_lm)$r.squared
#### Linear model only explains 10 % of the variance.
DHARMa::simulateResiduals(fittedModel = abund_tot_owr_lm, plot = T)
#### LM is not a bad fit though.





# Total biomass in infauna ==================================================================

## Infested vs Non-infested across all sites ##################################################

## Assumptions of normality and homogeneity of variances
ggqqplot(biom_tot, "biom_total", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_biom_tot <- biom_tot %>%
  group_by(infestation, site) %>%
  shapiro_test(biom_total)
View(shap_biom_tot)
#### Data is normally distributed (except for Mosselbaai)

levene_test(biom_tot$biom_total ~ biom_tot$site * biom_tot$infestation, data = biom_tot)
#### Strangely, the test is not working when I don't specify the columns. 
#### Variance is barely homogeneous.

### ANOVA (using aov) with unbalanced design
biom_tot_aov <- aov(biom_total ~ site * infestation, data = biom_tot)
Anova(biom_tot_aov, type = "III")
#### There is a significant effect of site and infestation on the total abundance per quadrat,
#### but no effect of the interaction

### Check ANOVA assumptions
plot(biom_tot_aov, 1) # Variances are not homogeneous
plot(biom_tot_aov, 2) # Data is not normaly distributed
#### Points 1, 2 and 3 were detected as outliers for normality and homogeneity of variances

### ANOVA (using models)
hist(biom_tot$biom_total)

biom_tot_lm <- lm(biom_total ~ site * infestation, data = biom_tot)
summary(biom_tot_lm)
summary(biom_tot_lm)$r.squared
#### Linear model only explains 50 % of the variance.
DHARMa::simulateResiduals(fittedModel = biom_tot_lm, plot = T)
#### LM is not a good fit as quantile deviations is detected.
AIC(biom_tot_lm)

simulationOutput <- DHARMa::simulateResiduals(fittedModel = biom_tot_lm, plot = F) 
plotResiduals(simulationOutput, biom_tot$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, biom_tot$site) # within-group deviations not significant

biom_tot_glm <- glm(biom_total ~  site * infestation, data = biom_tot,
                      family = "Gamma")
summary(biom_tot_glm)
DHARMa::simulateResiduals(fittedModel = biom_tot_glm, plot = T)
#### GLM is a much better fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = biom_tot_glm, plot = F) 
plotResiduals(simulationOutput, biom_tot$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, biom_tot$site) # within-group deviations not significant

biom_tot_glm1 <- glm(biom_total ~  site * infestation, data = biom_tot,
                    family = "Gamma"(link = "log"))
summary(biom_tot_glm1)
Anova(biom_tot_glm1)
DHARMa::simulateResiduals(fittedModel = biom_tot_glm1, plot = T)

### Post-hoc tests (using Tukey)
emmeans(biom_tot_glm1, pairwise ~ site * infestation)
#### By comparing the total biomass per site and infestation level:
####    - Total biomass not different between sites.

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Identify potential outliers
biom_tot_owr %>%
  group_by(infestation, lineage, quadrat) %>%
  identify_outliers(biom_total)

### Assumptions of normality and homogeneity of variances
ggqqplot(biom_tot_owr, "biom_total", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_biom_tot_owr <- biom_tot_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(biom_total)
View(shap_biom_tot_owr)
#### Data is normally distributed

levene_test(biom_tot_owr$biom_total ~ biom_tot_owr$lineage * biom_tot_owr$infestation, data = biom_tot_owr)
#### Variances are not homogeneous.

### ANOVA (using aov) with unbalanced design
biom_tot_owr_aov <- aov(biom_total ~ lineage * infestation, data = biom_tot_owr)
Anova(biom_tot_owr_aov, type = "III")
#### There is no significant effect of lineage or infestation or interaction on the total biomass

biom_tot_owr_aov1 <- aov(biom_total ~ lineage + infestation, data = biom_tot_owr)
Anova(biom_tot_owr_aov1, type = "II")
#### There is no significant effect of lineage or infestation on total biomass
#### ANOVA Type II is more powerful when interactions are not significant.

### Check ANOVA assumptions
plot(biom_tot_owr_aov, 1) # Variances are not homogeneous
plot(biom_tot_owr_aov, 2) # Data is not normally distributed

### ANOVA (using models)
hist(biom_tot_owr$biom_total)

biom_tot_owr_lm <- lm(biom_total ~ lineage * infestation, data = biom_tot_owr)
summary(biom_tot_owr_lm)
summary(biom_tot_owr_lm)$r.squared
#### Linear model only explains 10 % of the variance.
DHARMa::simulateResiduals(fittedModel = biom_tot_owr_lm, plot = T)
#### LM is not a bad fit though.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = biom_tot_owr_lm, plot = F) 
plotResiduals(simulationOutput, biom_tot_owr$infestation)  # within-group deviations not significant,
# but variance not homogeneous.
plotResiduals(simulationOutput, biom_tot_owr$site) # within-group deviations not significant




# Species richness ==================================================================

## Infested vs Non-infested across all sites ##################################################

## Assumptions of normality and homogeneity of variances
ggqqplot(richness, "count", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_rich <- richness %>%
  group_by(infestation, site) %>%
  shapiro_test(count)
View(shap_rich)
#### Data is normally distributed (except for OWR infested)

levene_test(count ~ richness$site * richness$infestation, data = richness)
#### Strangely, the test is not working at all.

### ANOVA (using aov) with unbalanced design
rich_aov <- aov(count ~ site * infestation, data = richness)
Anova(rich_aov, type = "III")
#### There is no significant effect of infestation on the species richness,
#### but there might be a significant effect of site (p-value close to 0.05)

### Check ANOVA assumptions
plot(rich_aov, 1) # Variances are homogeneous enough
plot(rich_aov, 2) # Data is normally distributed
#### Points 21, 23 and 19 were detected as outliers for normality and homogeneity of variances

### Post-hoc tests (using Tukey)
summary(glht(rich_aov, linfct = mcp(site = "Tukey")))
TukeyHSD((rich_aov1))
#### By comparing the species richness between sites:
####  - Species richness in OWR > MB

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Identify potential outliers
richness_owr %>%
  group_by(infestation, lineage, quadrat) %>%
  identify_outliers(count)

### Assumptions of normality and homogeneity of variances
ggqqplot(richness_owr, "count", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_rich_owr <- richness_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(count)
View(shap_rich_owr)
#### Data is normally distributed (except for Infested East)

levene_test(count ~ lineage * infestation, data = richness_owr)
#### Strangely, this is not working again.

### ANOVA (using aov) with unbalanced design
rich_owr_aov <- aov(count ~ lineage * infestation, data = richness_owr)
Anova(rich_owr_aov, type = "III")
#### There is no significant effect of lineage or infestation or interaction on the species richness.

rich_owr_aov1 <- aov(count ~ lineage + infestation, data = richness_owr)
Anova(rich_owr_aov1, type = "II")
#### There is no significant effect of lineage or infestation on total abundance.
#### ANOVA Type II is more powerful when interactions are not significant.

### Check ANOVA assumptions
plot(rich_owr_aov, 1) # Variances are not homogeneous
plot(rich_owr_aov, 2) # Normality




# Shannon-Wiener Diversity Index (H') ==================================================================

## Infested vs Non-infested across all sites ##################################################

## Assumptions of normality and homogeneity of variances
ggqqplot(env, "shannon", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_shannon <- env %>%
  group_by(infestation, site) %>%
  shapiro_test(shannon)
View(shap_shannon)
#### Data is normally distributed

levene_test(shannon ~ site * infestation, data = env)
#### Various are not homogeneous

### ANOVA (using aov) with unbalanced design
shannon_aov <- aov(shannon ~ site * infestation, data = env)
Anova(shannon_aov, type = "III")
#### There is no significant effect of infestation or site or interaction on species diversity
#### calculated with Shannon-Wiener diversity index (H')

### Check ANOVA assumptions
plot(shannon_aov, 1) # Variances are homogeneous
plot(shannon_aov, 2) # Data is normally distributed
#### Points 24, 28, 30 were detected as outliers for normality and homogeneity of variances

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Assumptions of normality and homogeneity of variances
ggqqplot(env_owr, "shannon", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_shannon_owr <- env_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(shannon)
View(shap_shannon_owr)
#### Data is normally distributed

levene_test(shannon ~ lineage * infestation, data = env_owr)
#### Variances are not homogeneous

### ANOVA (using aov) with unbalanced design
shannon_owr_aov <- aov(shannon ~ lineage * infestation, data = env_owr)
Anova(shannon_owr_aov, type = "III")
#### There is no significant effect of lineage or infestation or interaction on the species richness.

### Check ANOVA assumptions
plot(shannon_owr_aov, 1) # Variances are homogeneous
plot(shannon_owr_aov, 2) # Data is not normally distributed





# Simpson's Diversity Index (λ)  ==================================================================

## Infested vs Non-infested across all sites ##################################################

## Assumptions of normality and homogeneity of variances
ggqqplot(env, "simpson", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_simpson <- env %>%
  group_by(infestation, site) %>%
  shapiro_test(simpson)
View(shap_simpson)
#### Data is normally distributed

levene_test(simpson ~ site * infestation, data = env)
#### Various are not homogeneous

### ANOVA (using aov) with unbalanced design
simpson_aov <- aov(simpson ~ site * infestation, data = env)
Anova(simpson_aov, type = "III")
#### There is no significant effect of infestation or site or interaction on species diversity
#### calculated with Simpson's Diversity Index (λ)

### Check ANOVA assumptions
plot(simpson_aov, 1) # Variances are homogeneous
plot(simpson_aov, 2) # Data is normally distributed
#### Points 9, 10 and 28 were detected as outliers for normality and homogeneity of variances

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Assumptions of normality and homogeneity of variances
ggqqplot(env_owr, "simpson", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_simpson_owr <- env_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(simpson)
View(shap_simpson_owr)
#### Data is normally distributed

levene_test(simpson ~ lineage * infestation, data = env_owr)
#### Variances are not homogeneous

### ANOVA (using aov) with unbalanced design
simpson_owr_aov <- aov(simpson ~ lineage * infestation, data = env_owr)
Anova(simpson_owr_aov, type = "III")
#### There is no significant effect of lineage or infestation or interaction on the species diversity.

### Check ANOVA assumptions
plot(simpson_owr_aov, 1) # Variances are homogeneous enough
plot(simpson_owr_aov, 2) # Data is normally distributed
#### Points 2, 4 and 5 were detected as outliers for normality and homogeneity of variances





# Pielou's evenness (J)  ==================================================================

## Infested vs Non-infested across all sites ##################################################

## Assumptions of normality and homogeneity of variances
ggqqplot(env, "pielou", ggtheme = theme_bw()) +
  facet_grid(infestation ~ site, labeller = "label_both")

shap_pielou <- env %>%
  group_by(infestation, site) %>%
  shapiro_test(pielou)
View(shap_pielou)
#### Data is normally distributed

levene_test(pielou ~ site * infestation, data = env)
#### Various are homogeneous

### ANOVA (using aov) with unbalanced design
pielou_aov <- aov(pielou ~ site * infestation, data = env)
Anova(pielou_aov, type = "III")
#### There is a significant effect of site on the evenness index, but
#### no significant effect of infestation or interaction.

### Check ANOVA assumptions
plot(pielou_aov, 1) # Variances are homogeneous
plot(pielou_aov, 2) # Data is not normally distributed
#### Points 10, 26 and 28 were detected as outliers for normality and homogeneity of variances

hist(env$pielou) # Distribution could be considered bimodal.

### Post-hoc tests (using Tukey)
summary(glht(pielou_aov, linfct = mcp(site = "Tukey")))
TukeyHSD((pielou_aov))
#### By comparing the species evenness:
####    - Species evenness PEd > JB

## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Assumptions of normality and homogeneity of variances
ggqqplot(env_owr, "pielou", ggtheme = theme_bw()) +
  facet_grid(infestation ~ lineage, labeller = "label_both")

shap_pielou_owr <- env_owr %>%
  group_by(infestation, lineage) %>%
  shapiro_test(pielou)
View(shap_pielou_owr)
#### Data is normally distributed

levene_test(pielou ~ lineage * infestation, data = env_owr)
#### Variances are homogeneous

### ANOVA (using aov) with unbalanced design
pielou_owr_aov <- aov(pielou ~ lineage * infestation, data = env_owr)
Anova(pielou_owr_aov, type = "III")
#### There is no significant effect of lineage or infestation or interaction on the species diversity.

### Check ANOVA assumptions
plot(pielou_owr_aov, 1) # Variances are homogeneous enough
plot(pielou_owr_aov, 2) # Data is normally distributed
#### Points 2, 4 and 5 were detected as outliers for normality and homogeneity of variances




# Community analysis on abundance ========================================================

### Clean the dataset

colnames(abund_wide1) <- make.cepnames(colnames(abund_wide1))
# make.cepnames to shorten latin names into eight-character variables
View(abund_wide1)

## Infested vs Non-infested across all sites ##################################################

### Calculate the distance matrix
# Compute distance matrix using Bray-Curtis on transformed abundances
range(abund_wide1)
range(abund_wide1^0.5) 
range(abund_wide1^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
abund_dist <- vegdist(abund_wide1^0.25, method='bray')
abund_dist

### Calculate nMDS

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.005 --> Good fit

abund_bray <- metaMDS(abund_wide1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
abund_bray # Stress = 0.180 shows this nMDS is a fair fit to the data.
View(abund_bray)

### Exploring the results of the nMDS
plot(abund_bray, type = "t")
stressplot(abund_bray)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(abund_bray)
plot(abund_bray, type = "t", main = "goodness of fit")
points(abund_bray, display = "sites", cex = gof*100)

View(abund_wide1)

abund_fit <- envfit(abund_bray, env1, permu = 999) # Caution: envfit does not allow missing values
abund_fit
# The number of live mussels and dead mussels and the number of broken shells are the most
# strongly related to the first two ordination axes (r2 = 0.21, 0.28 and 0.20 respectively).
# Dead mussels are strongly related to the first ordination axis. 
# Live mussels and broken shells are strongly related to the second ordination axis.
install.packages("RVAideMemoire")
library(RVAideMemoire)
pairwise.factorfit(abund_bray, env1$site)

abund_fit_adj <- abund_fit
pvals_adj <- p.adjust(abund_fit$vectors$pvals, method = 'bonferroni')
abund_fit_adj$vectors$pvals <- pvals_adj
abund_fit_adj
# W/ a Bonferroni correction, the p-values are still the same.

plot(abund_bray, display = "sites")
plot(abund_fit_adj, p.max = 0.05) 
# Infestation is not significant. 

## Plot nMDS w/ infestation levels and sites

ordiplot(abund_bray, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(abund_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5)[env1$site],
       col = my_colors[env1$infestation])
ordihull(abund_bray, groups = env1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.35, 0.3, legend = c(levels(env1$site), levels(env1$infestation)),
       pch = c(16, 17, 18, 8, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "black", "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.70, 1.25, "Stress = 0.180", bty = "n", cex = 2)

## PERMANOVA
# Run PERMANOVA
abund_pmv <- adonis2(abund_wide1^0.25 ~ infestation * site, data = env1,
               permutations = 999, method = 'bray')
abund_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities DO NOT differ significantly between infestation levels (p = 0.82).
# But communities differ between sites (p = 0.001)

# Plot permuted F-values
densityplot(permustats(abund_pmv))

# Plot the average distance to median for infestation levels
abund_disp_inf <- betadisper(abund_dist, group = env1$infestation)
boxplot(abund_disp_inf)
# No difference in dispersion detected between infestation levels from the boxplot.
abund_disp_inf

# Is there a significant difference in dispersion between infestation levels ?
anova(abund_disp_inf) 
permutest(abund_disp_inf)
# No significant difference in dispersion between infestation levels
plot(abund_disp_inf, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     xlab = "Dimension 1 (22.31 %)", ylab = "Dimension 2 (19.24%)",
     main = ""
     )
round(100 * abund_disp_inf$eig / sum(abund_disp_inf$eig), 2)

# Plot the average distance to median for sites
abund_disp_site <- betadisper(abund_dist, group = env1$site)
boxplot(abund_disp_site)
abund_disp_site
# There is a clear difference in dispersion detected between sites from the boxplot.
anova(abund_disp_site)
permutest(abund_disp_site)
plot(abund_disp_site, hull = F, ellipse = T,
     pch = c(16:18, 8, 5), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (22.31%)", ylab = "Dimension 2 (19.24%)",
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     main = "")

## Cluster analysis
abund_clust <- hclust(abund_dist, "ward.D2")
levels(abund_wide$site) <- c("BoS", "JB", "MB", "OWR", "PEd")

abund_inf <- rev(levels(abund_wide$infestation))
abund_dend <- as.dendrogram(abund_clust)

labels_colors(abund_dend) <- my_colors[(abund_wide$infestation)[order.dendrogram(abund_dend)]]

labels(abund_dend) <- paste(abund_wide$site[order.dendrogram(abund_dend)],
                            " (", labels(abund_dend), ")", sep = "")
branches_attr_by_lists(abund_dend, abund_wide$site[order.dendrogram(abund_dend)],
                       attr = "lwd")
abund_dend <- set(abund_dend, "labels_cex", 1.5)
abund_dend <- hang.dendrogram(abund_dend, hang_height = 0.005)
plot(abund_dend, nodePar = list(cex = 0.007))
legend(35, 2, legend = levels(abund_wide$infestation), fill = my_colors, bty = "n", cex = 1.5)

## SIMPER on infestation
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)
abund_simp_inf <- simper(abund_wide1^0.25, group = env1$infestation)
abund_simp_inf
summary(abund_simp_inf)

## SIMPER on sites
abund_simp_site <- simper(abund_wide1^0.25, group = env1$site)
abund_simp_site
summary_simper <- summary(abund_simp_site)
View(summary_simper$`Old Woman's River_Port Edward`)




## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Calculate the distance matrix
# Compute distance matrix using Bray-Curtis on transformed abundances
range(abund_wide1_owr)
range(abund_wide1_owr^0.5) 
range(abund_wide1_owr^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
abund_dist_owr <- vegdist(abund_wide1_owr^0.25, method='bray')
abund_dist_owr

### Calculate nMDS

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.005 --> Good fit

abund_bray_owr <- metaMDS(abund_wide1_owr^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
abund_bray_owr # Stress = 0.116 shows this nMDS is a fair fit to the data.
View(abund_bray_owr)

### Exploring the results of the nMDS
plot(abund_bray_owr, type = "t")
stressplot(abund_bray_owr)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(abund_bray_owr)
plot(abund_bray_owr, type = "t", main = "goodness of fit")
points(abund_bray_owr, display = "sites", cex = gof*100)

View(abund_wide1)

abund_fit_owr <- envfit(abund_bray_owr, env1_owr, permu = 999) # Caution: envfit does not allow missing values
abund_fit_owr
# The mean number of byssal threads is the most strongly to the second ordination axis (r2 = 0.43)

abund_fit_adj_owr <- abund_fit_owr
pvals_adj_owr <- p.adjust(abund_fit_owr$vectors$pvals, method = 'bonferroni')
abund_fit_adj_owr$vectors$pvals <- pvals_adj_owr
abund_fit_adj_owr
# W/ a Bonferroni correction, this is not the case anymore.

plot(abund_bray_owr, display = "sites")
plot(abund_fit_adj_owr, p.max = 0.05) 
# Infestation and lineage are not significant. 

## Plot nMDS w/ infestation levels and lineages

ordiplot(abund_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(abund_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env1_owr$lineage],
       col = my_colors[env1_owr$infestation])
ordihull(abund_bray_owr, groups = env1_owr$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.1, -0.3, legend = c(levels(env1_owr$lineage), levels(env1_owr$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.45, 1.25, "Stress = 0.116", bty = "n", cex = 2)

## PERMANOVA
# Run PERMANOVA
abund_owr_pmv <- adonis2(abund_wide1_owr^0.25 ~ infestation * lineage, data = env1_owr,
                     permutations = 999, method = 'bray')
abund_owr_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities DO NOT differ significantly between infestation levels (p = 0.82).
# But communities differ between sites (p = 0.001)

# Plot permuted F-values
densityplot(permustats(abund_owr_pmv))

# Plot the average distance to median for infestation levels
abund_disp_inf_owr <- betadisper(abund_dist_owr, group = env1_owr$infestation)
boxplot(abund_disp_inf_owr)
# Difference in dispersion detected between infestation levels
abund_disp_inf_owr

# Is there a significant difference in dispersion between infestation levels ?
anova(abund_disp_inf_owr) 
permutest(abund_disp_inf_owr)
# No significant difference in dispersion between infestation levels
plot(abund_disp_inf_owr, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21%)",
     main = ""
)
round(100 * abund_disp_inf_owr$eig / sum(abund_disp_inf_owr$eig), 2)

# Plot the average distance to median for sites
abund_disp_lineage <- betadisper(abund_dist_owr, group = env1_owr$lineage)
boxplot(abund_disp_lineage)
abund_disp_lineage
# There is no clear difference in dispersion detected between lineages from the boxplot.
anova(abund_disp_lineage)
permutest(abund_disp_lineage)
plot(abund_disp_lineage, hull = F, ellipse = T,
     pch = c(16:17), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21 %)",
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     main = "")

## Cluster analysis
abund_clust_owr <- hclust(abund_dist_owr, "ward.D2")

abund_inf_owr <- rev(levels(abund_wide_owr$infestation))
abund_dend_owr <- as.dendrogram(abund_clust_owr)

labels_colors(abund_dend_owr) <- my_colors[(abund_wide_owr$infestation)[order.dendrogram(abund_dend_owr)]]

labels(abund_dend_owr) <- paste(abund_wide_owr$lineage[order.dendrogram(abund_dend_owr)],
                            " (", labels(abund_dend_owr), ")", sep = "")
abund_dend_owr <- set(abund_dend_owr, "labels_cex", 0.9)
abund_dend_owr <- hang.dendrogram(abund_dend_owr, hang_height = 0.005)
plot(abund_dend_owr, nodePar = list(cex = 0.007))
legend("topright", legend = levels(abund_wide$infestation), fill = my_colors)


## SIMPER on infestation
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)
abund_simp_inf_owr <- simper(abund_wide1_owr^0.25, group = env1_owr$infestation)
abund_simp_inf_owr
summary(abund_simp_inf_owr)

## SIMPER on sites
abund_simp_lineage <- simper(abund_wide1_owr^0.25, group = env1_owr$lineage)
abund_simp_lineage
summary(abund_simp_lineage)
View(summary_simper)




# Community analysis on biomass ========================================================

### Clean the dataset

colnames(biom_wide1) <- make.cepnames(colnames(biom_wide1))
# make.cepnames to shorten latin names into eight-character variables
View(biom_wide1)

## Infested vs Non-infested across all sites ##################################################

### Calculate the distance matrix
# Compute distance matrix using Bray-Curtis on transformed abundances
range(biom_wide1)
range(biom_wide1^0.5) 
range(biom_wide1^0.25) # Range about 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
biom_dist <- vegdist(biom_wide1^0.25, method='bray')
biom_dist

### Calculate nMDS

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.005 --> Good fit

biom_bray <- metaMDS(biom_wide1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
biom_bray # Stress = 0.190 shows this nMDS is a bit suspect fit to the data.
View(biom_bray)

### Exploring the results of the nMDS
plot(biom_bray, type = "t")
stressplot(biom_bray)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(biom_bray)
plot(biom_bray, type = "t", main = "goodness of fit")
points(biom_bray, display = "sites", cex = gof*100)

View(abund_wide1)

biom_fit <- envfit(biom_bray, env1, permu = 999) # Caution: envfit does not allow missing values
biom_fit
# The number of dead mussels and broken shells are the most strongly related to the first two
# ordination axes (r2 = 0.36 and 0.24). Both strongly related to the second ordination axis.

biom_fit_adj <- biom_fit
pvals_adj <- p.adjust(biom_fit$vectors$pvals, method = 'bonferroni')
biom_fit_adj$vectors$pvals <- pvals_adj
biom_fit_adj
# W/ a Bonferroni correction, the p-values are still the same.

plot(biom_bray, display = "sites")
plot(biom_fit_adj, p.max = 0.05) 
# Infestation is not significant. 

## Plot nMDS w/ infestation levels and sites

ordiplot(biom_bray, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(biom_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 8, 5)[env1$site],
       col = my_colors[env1$infestation])
ordihull(biom_bray, groups = env1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.1, 0.3, legend = c(levels(env1$site), levels(env1$infestation)),
       pch = c(16, 17, 18, 8, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "black", "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.45, 1.25, "Stress = 0.180", bty = "n", cex = 2)

## PERMANOVA
# Run PERMANOVA
biom_pmv <- adonis2(biom_wide1^0.25 ~ infestation * site, data = env1,
                     permutations = 999, method = 'bray')
biom_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities DO NOT differ significantly between infestation levels (p = 0.82).
# But communities differ between sites (p = 0.001)

# Plot permuted F-values
densityplot(permustats(biom_pmv))

# Plot the average distance to median for infestation levels
biom_disp_inf <- betadisper(biom_dist, group = env1$infestation)
boxplot(biom_disp_inf)
# Difference in dispersion detected between infestation levels for infaunal biomass
biom_disp_inf

# Is there a significant difference in dispersion between infestation levels ?
anova(biom_disp_inf) 
permutest(biom_disp_inf)
# No significant difference in dispersion between infestation levels
plot(biom_disp_inf, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5),
     xlab = "Dimension 1 (24.71 %)", ylab = "Dimension 2 (18.63 %)",
     main = ""
)
round(100 * biom_disp_inf$eig / sum(biom_disp_inf$eig), 2)

# Plot the average distance to median for sites
biom_disp_site <- betadisper(biom_dist, group = env1$site)
boxplot(biom_disp_site)
biom_disp_site
# There is a clear difference in dispersion detected between sites from the boxplot.
anova(biom_disp_site)
permutest(biom_disp_site)
plot(biom_disp_site, hull = F, ellipse = T,
     pch = c(16:18, 8, 5), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (24.71 %)", ylab = "Dimension 2 (18.63 %)",
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5),
     main = "")


## Cluster analysis - all sites on biomass
biom_clust <- hclust(biom_dist, "ward.D2")
levels(biom_wide$site) <- c("BoS", "JB", "MB", "OWR", "PEd")

biom_inf <- rev(levels(biom_wide$infestation))
biom_dend <- as.dendrogram(biom_clust)

labels_colors(biom_dend) <- my_colors[(biom_wide$infestation)[order.dendrogram(biom_dend)]]

labels(biom_dend) <- paste(biom_wide$site[order.dendrogram(biom_dend)],
                            " (", labels(biom_dend), ")", sep = "")
biom_dend <- set(biom_dend, "labels_cex", 1.5)
biom_dend <- hang.dendrogram(biom_dend, hang_height = 0.005)
plot(biom_dend, nodePar = list(cex = 0.007))
legend(35, 2, legend = levels(abund_wide$infestation), fill = my_colors, bty = "n", cex = 1.5)


## SIMPER on infestation
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)
biom_simp_inf <- simper(biom_wide1^0.25, group = env1$infestation)
biom_simp_inf
summary(biom_simp_inf)

## SIMPER on sites
biom_simp_site <- simper(biom_wide1^0.25, group = env1$site)
biom_simp_site
summary_biom_simper <- summary(biom_simp_site)
View(summary_biom_simper$`Old Woman's River_Port Edward`)




## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

### Calculate the distance matrix
# Compute distance matrix using Bray-Curtis on transformed abundances
range(biom_wide1_owr)
range(biom_wide1_owr^0.5) 
range(biom_wide1_owr^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
biom_dist_owr <- vegdist(biom_wide1_owr^0.25, method='bray')
biom_dist_owr

### Calculate nMDS

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.005 --> Good fit

biom_bray_owr <- metaMDS(biom_wide1_owr^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
biom_bray_owr # Stress = 0.116 shows this nMDS is a fair fit to the data.

### Exploring the results of the nMDS
plot(biom_bray_owr, type = "t")
stressplot(biom_bray_owr)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(biom_bray_owr)
plot(biom_bray_owr, type = "t", main = "goodness of fit")
points(biom_bray_owr, display = "sites", cex = gof*100)

biom_fit_owr <- envfit(biom_bray_owr, env1_owr, permu = 999) # Caution: envfit does not allow missing values
biom_fit_owr
# There is no environmental variables that strongly with ordination axes.

biom_fit_adj_owr <- biom_fit_owr
pvals_adj_owr <- p.adjust(biom_fit_owr$vectors$pvals, method = 'bonferroni')
biom_fit_adj_owr$vectors$pvals <- pvals_adj_owr
biom_fit_adj_owr
# W/ a Bonferroni correction, this is the same. 

plot(biom_bray_owr, display = "sites")
plot(biom_fit_adj_owr, p.max = 0.05) 
# Infestation and lineage are not significant. 

## Plot nMDS w/ infestation levels and lineages

ordiplot(biom_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(biom_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env1_owr$lineage],
       col = my_colors[env1_owr$infestation])
ordihull(biom_bray_owr, groups = env1_owr$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.1, -0.3, legend = c(levels(env1_owr$lineage), levels(env1_owr$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.45, 1.25, "Stress = 0.116", bty = "n", cex = 2)

## PERMANOVA
# Run PERMANOVA
biom_owr_pmv <- adonis2(biom_wide1_owr^0.25 ~ infestation * lineage, data = env1_owr,
                         permutations = 999, method = 'bray')
biom_owr_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities DO NOT differ significantly between infestation levels (p = 0.82).
# But communities differ between sites (p = 0.001)

# Plot permuted F-values
densityplot(permustats(biom_owr_pmv))

# Plot the average distance to median for infestation levels
biom_disp_inf_owr <- betadisper(biom_dist_owr, group = env1_owr$infestation)
boxplot(biom_disp_inf_owr)
# Difference in dispersion detected between infestation levels
biom_disp_inf_owr

# Is there a significant difference in dispersion between infestation levels ?
anova(biom_disp_inf_owr) 
permutest(biom_disp_inf_owr)
# No significant difference in dispersion between infestation levels
plot(biom_disp_inf_owr, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     xlab = "Dimension 1 (44.42 %)", ylab = "Dimension 2 (15.75 %)",
     main = ""
)
round(100 * biom_disp_inf_owr$eig / sum(biom_disp_inf_owr$eig), 2)

# Plot the average distance to median for sites
biom_disp_lineage <- betadisper(biom_dist_owr, group = env1_owr$lineage)
boxplot(biom_disp_lineage)
biom_disp_lineage
# There is no potential difference in dispersion detected between lineages from the boxplot.
anova(biom_disp_lineage)
permutest(biom_disp_lineage)
plot(biom_disp_lineage, hull = F, ellipse = T,
     pch = c(16:17), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21 %)",
     xlim = c(-0.5, 0.5), ylim = c(-0.4, 0.4),
     main = "")

## Cluster analysis - OWR biomass
biom_clust_owr <- hclust(biom_dist_owr, "ward.D2")

biom_inf_owr <- rev(levels(biom_wide_owr$infestation))
biom_dend_owr <- as.dendrogram(biom_clust_owr)

labels_colors(biom_dend_owr) <- my_colors[(biom_wide_owr$infestation)[order.dendrogram(biom_dend_owr)]]

labels(biom_dend_owr) <- paste(biom_wide_owr$lineage[order.dendrogram(biom_dend_owr)],
                                " (", labels(biom_dend_owr), ")", sep = "")
biom_dend_owr <- set(biom_dend_owr, "labels_cex", 1.5)
biom_dend_owr <- hang.dendrogram(biom_dend_owr, hang_height = 0.005)
plot(biom_dend_owr, nodePar = list(cex = 0.007))


## SIMPER on infestation
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)
biom_simp_inf_owr <- simper(biom_wide1_owr^0.25, group = env1_owr$infestation)
biom_simp_inf_owr
summary(biom_simp_inf_owr)

## SIMPER on sites
biom_simp_lineage <- simper(biom_wide1_owr^0.25, group = env1_owr$lineage)
biom_simp_lineage
summary(biom_simp_lineage)




###################################################################################
# STATISTICAL ANALYSIS - SUMMARY ----------------------------------------------------
###################################################################################

# Mussel bed architecture - mean number, sd and se of live mussels per quadrat ==============
##  - No statistical difference between infestation levels across sites,
##    but significant statistical difference between sites
##  - Barely significant statistical difference between infestation levels at Old Woman's River,
##    but no statistical difference between Perna lineages

# Mussel bed architecture - mean number, sd and se of byssal threads per quadrat ==============
##  - No statistical difference between infestation levels across sites,
##    but significant statistical difference between sites
##  - No statistical difference between infestation levels and Perna lineages at Old Woman's River

