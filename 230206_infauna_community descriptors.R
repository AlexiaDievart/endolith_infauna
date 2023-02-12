## ---------------------------
##
## Script name: Infaunal communities (descriptive) associated with manipulated mussel beds depending on their level
##              of euendolithic infestation and their transplant site
##
## Purpose of script: 
##    - Is there a difference between the total abundance of infauna between each quadrat,
##      each euendolithic infestation treatment and each site ?
##    - Is there a difference between the species richness (S) between each quadrat,
##      each euendolithic infestation treatment and each site ?
##    - Is there a difference between the species diversity, calculated w/ the Shannon-Wiener index (H') 
##      between each quadrat, each euendolithic infestation treatment and each site ?
##    - Is there a difference between the species diversity, calculated w/ the Simpson index (λ) 
##      between each quadrat, each euendolithic infestation treatment and each site ?
##    - Is there a difference between the species evenness, calculated w/ the Pielou's index (J) 
##      between each quadrat, each euendolithic infestation treatment and each site ?
##
## Author: Alexia Dievart
##
## Date Created: 2023-02-05
## Dates Updated: 
##
## Copyright (c) Alexia DIEVART 2023
## Email: alexia.dievart@hotmail.fr
##
## ---------------------------
##
## Notes: 
##   - Quadrats cannot be specified as a random factor (n < 5).
##   - Analyses conducted ACROSS ALL SITES excluded Old Woman's River, as the effect of site and Perna lineages
##     are intricately confounded. 
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

# Environmental variables - total number of mussels and number of broken shells ===============
env <- read.csv("./RAW DATA/infauna_description.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(env)

factor1 <- c("site", "infestation", "lineage", "quadrat")
env[,factor1] <- lapply(env[,factor1], factor)

env <- env %>% drop_na(live_perna)

View(env)

env1 <- env[, c(1:4, 9:10)]
env1 <- env1[env1$site !="Old Woman's River",]
View(env1)

# Community variables - e.g. abundance, biomass ===============================================

# Load dataset ################################################################################
community <- read.csv("./RAW DATA/infauna_community.csv", dec = ",", header = T, sep = ";")

factor2 <- c("site", "bioregion", "infestation", "lineage", "quadrat", "higher_group", "species", "pic")
community[,factor2] <- lapply(community[,factor2], factor)
dplyr::glimpse(community)

View(community)

# Isolate ABUNDANCE ###########################################################################
abund_long <- community[,c(1:5,7:8)]
abund_long1 <- abund_long[abund_long$site != "Old Woman's River",]
View(abund_long)

abund_wide <- spread(abund_long, species, abundance)
abund_wide[is.na(abund_wide)] <- 0
View(abund_wide)

abund_wide1 <- abund_wide[abund_wide$site != "Old Woman's River",]
View(abund_wide1)

abund_long_owr <- abund_long[abund_long$site == "Old Woman's River",]
View(abund_long_owr)

abund_wide_owr <- spread(abund_long_owr, species, abundance)
abund_wide_owr[is.na(abund_wide_owr)] <- 0
View(abund_wide_owr)




# Isolate BIOMASS #############################################################################
biom_long <- community[,c(1:5,7, 9)]
biom_long1 <- biom_long[biom_long$site != "Old Woman's River",]
View(biom_long1)

biom_wide <- spread(biom_long, species, biomass)
biom_wide[is.na(biom_wide)] <- 0
View(biom_wide)

biom_long_owr <- biom_long[biom_long$site == "Old Woman's River",]
View(biom_long_owr)

biom_wide_owr <- spread(biom_long_owr, species, biomass)
biom_wide_owr[is.na(biom_wide_owr)] <- 0
View(biom_wide_owr)




###############################################################################################
# Section: DESCRIPTION of infaunal communities ------------------------------------------------
###############################################################################################

###############################################################################################
# Unique species and common species ===========================================================
###############################################################################################

# Total number of species across all sites (including OWR) ####################################
abund_long %>% summarise(n_distinct(species))

# Across all sites (excluding OWR) ############################################################

# Number of unique species per infestation level across all sites (excluding OWR)
abund_long1 %>%
  group_by(infestation) %>%
  summarise(n_distinct(species))
# n(infested) = 80 and n(non-infested) = 80

# Number and identity of common species between infestation levels across all sites (excluding OWR)
Reduce(intersect, split(abund_long1$species, abund_long1$infestation))

# At Old Woman's River #########################################################################

# Number of unique species per infestation level and Perna lineage at Old Woman's River
abund_long %>%
  filter(site == "Old Woman's River") %>%
  group_by(infestation, lineage) %>%
  summarise(n_distinct(species))

# Number and identity of common species between infestation levels and Perna lineages at Old Woman's River
abund_owr <- abund_long %>%
  filter(site == "Old Woman's River") %>%
  dplyr::group_by(infestation, lineage) %>%
  tidyr::unite(col = "id", infestation, lineage, sep = "_", remove = FALSE)
Reduce(intersect, split(abund_owr$species, abund_owr$id))




###############################################################################################
# Total abundance across all sites  ===========================================================
###############################################################################################

# Curate dataset ##############################################################################
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

# Barplot of the total abundance of infauna across all sites (excluding OWR) ##################

# Barplot of the total abundance of infauna per quadrat, site and infestation levels
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

# Barplot of the total abundance of infauna per site and infestation levels
barplot_abund_tot1 <- ggplot(data = abund_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = ab_tot_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_sd, ymax = ab_tot_mean + ab_tot_sd), width = 0.2,
                position = position_dodge(0.9)) +
#  annotate("text", 
#           x = c(0.85, 1.25, 1.85, 2.3, 2.8, 3.25, 3.8, 4.25, 4.8, 5.25),
#           y = c(700, 1050, 700, 580, 300, 300, 500, 460, 220, 170),
#           label = c("AB", "A", "AB", "AB", "B", "B", "B", "B", "B", "B"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total abundance")
barplot_abund_tot1

# Statistical analyses on the total abundance of infauna #########################################

# ANOVA w/ unbalanced design (i.e., not equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_totab <- ggboxplot(abund_tot, x = "site", y = "ab_total", color = "infestation",
                                palette = my_colors)
boxplot_totab

interaction_totab <- ggline(abund_tot, x = "site", y = "ab_total", color = "infestation",
                                  add = c("mean_se", "dotplot"),
                                  palette = my_colors)
interaction_totab

hist(abund_tot$ab_total)

# Type III ANOVA w/ interaction
anova_totab <- aov(ab_total ~ infestation * site, data = abund_tot)
Anova(anova_totab, type = "III") # Interaction is not significant

plot(anova_totab, 1)
plot(anova_totab, 2)
anova_totab_residuals <- residuals(object = anova_totab)
shapiro.test(x = anova_totab_residuals) # Residuals are NOT normally distributed (p = 0.02459)

# Type II ANOVA without interaction --> to use because more powerful
anova_totab1 <- aov(ab_total ~ infestation + site, data = abund_tot)
Anova(anova_totab1, type = "II")

plot(anova_totab1, 1)
leveneTest(ab_total ~ infestation * site, data = abund_tot) # Variances are homogeneous

plot(anova_totab1, 2)
anova_totab1_residuals <- residuals(object = anova_totab1)
shapiro.test(x = anova_totab1_residuals) # Residuals are NOT normally distributed (p = 0.0415)

# ANOVA using models

abund_tot_lm <- lm(ab_total ~ site * infestation, data = abund_tot)
summary(abund_tot_lm)
summary(abund_tot_lm)$r.squared
#### Linear model only explains 70 % of the variance.
DHARMa::simulateResiduals(fittedModel = abund_tot_lm, plot = T)
#### LM is not a bad fit. 
AIC(abund_tot_lm)
Anova(abund_tot_lm)

abund_tot_glm <- glm(ab_total ~ site * infestation, data = abund_tot,
                    family = "Gamma")
summary(abund_tot_glm)
DHARMa::simulateResiduals(fittedModel = abund_tot_glm, plot = T)
#### GLM is a much better fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = abund_tot_glm, plot = F) 
plotResiduals(simulationOutput, abund_tot$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, abund_tot$site) # within-group deviations not significant
Anova(abund_tot_glm)

abund_tot_glm1 <- glm(ab_total ~ infestation + site, data = abund_tot,
                     family = "Gamma")
summary(abund_tot_glm1)
DHARMa::simulateResiduals(fittedModel = abund_tot_glm1, plot = T)
#### GLM is a much better fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = abund_tot_glm1, plot = F) 
plotResiduals(simulationOutput, abund_tot$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, abund_tot$site) # within-group deviations not significant
Anova(abund_tot_glm1)

# Pairwise comparisons
summary(glht(abund_tot_glm1, linfct = mcp(site = "Tukey")))



###############################################################################################
# Total biomass across all sites  ============================================================
###############################################################################################

# Curate dataset ##############################################################################
biom_tot <- biom_long1 %>%
  na.omit(biom_long) %>%
  dplyr::group_by(site, infestation, quadrat) %>%
  summarize(biom_total = sum(biomass))
View(biom_tot)

biom_tot1 <- biom_tot  %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    biom_tot_mean = mean(biom_total),
    biom_tot_sd = sd(biom_total)
  )
View(biom_tot1)

# Barplot of the total biomass of infauna across all sites (excluding OWR) ####################

# Barplot of the total abundance of infauna for each quadrat, site and infestation level
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

# Barplot of the total abundance of infauna for site and infestation level
barplot_biom_tot1 <- ggplot(data = biom_tot1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = biom_tot_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_tot_mean - biom_tot_sd, ymax = biom_tot_mean + biom_tot_sd), width = 0.2,
                position = position_dodge(0.9)) +
#  annotate("text", 
#           x = c(0.85, 1.25, 1.85, 2.3, 2.8, 3.25, 3.8, 4.25, 4.8, 5.25),
#           y = c(16000, 13000, 7000, 6000, 5000, 5100, 9000, 8500, 3000, 3000),
#           label = c("A", "A", "AB", "AB", "AB", "AB", "AC", "AC", "BC", "B"), hjust = 1) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average total biomass (mg)")
barplot_biom_tot1

# Statistical analyses on the total biomass of infauna #########################################

# ANOVA w/ unbalanced design (i.e., not equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_totbiom <- ggboxplot(biom_tot, x = "site", y = "biom_total", color = "infestation",
                           palette = my_colors)
boxplot_totbiom

interaction_totbiom <- ggline(biom_tot, x = "site", y = "biom_total", color = "infestation",
                            add = c("mean_se", "dotplot"),
                            palette = my_colors)
interaction_totbiom

hist(biom_tot$biom_total)

# Type III ANOVA w/ interaction
biom_tot_aov <- aov(biom_total ~ site * infestation, data = biom_tot)
Anova(biom_tot_aov, type = "III") # No significant interaction

# Type II ANOVA without interaction --> to use because more powerful
biom_tot_aov1 <- aov(biom_total ~ infestation + site, data = biom_tot)
Anova(biom_tot_aov1, type = "II")

plot(biom_tot_aov1, 1)
leveneTest(biom_total ~ site * infestation, data = biom_tot) # Variances are homogeneous

plot(biom_tot_aov1, 2)
anova_totbiom1_residuals <- residuals(object = biom_tot_aov1)
shapiro.test(x = anova_totbiom1_residuals) # Residuals are NOT normally distributed (p = 0.001)

# ANOVA (using models)

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
Anova(biom_tot_glm)

biom_tot_glm1 <- glm(biom_total ~  site + infestation, data = biom_tot,
                    family = "Gamma")
summary(biom_tot_glm1)
DHARMa::simulateResiduals(fittedModel = biom_tot_glm1, plot = T)
#### GLM is a much better fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = biom_tot_glm1, plot = F) 
plotResiduals(simulationOutput, biom_tot$infestation)  # within-group deviations not significant
plotResiduals(simulationOutput, biom_tot$site) # within-group deviations not significant
Anova(biom_tot_glm1)

AIC(biom_tot_glm); AIC(biom_tot_glm1)

### Post-hoc tests (using Tukey)
summary(glht(biom_tot_glm1, linfct = mcp(site = "Tukey")))




###############################################################################################
# Species richness across all sites  ==========================================================
###############################################################################################

# Curate dataset ##############################################################################
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

# Barplot of the species richness of infauna across all sites (excluding OWR) ##################

# Barplot of species richness of infauna per quadrat, per site and per infestation
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

# Barplot of species richness of infauna per site and per infestation
barplot_rich1 <- ggplot(data = richness1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = count_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9)) +
#  annotate("text", 
#           x = c(1, 2, 3, 4, 5),
#           y = 55,
#           label = c("AB", "AB", "A", "B", "AB")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average species richness")
barplot_rich1

# Statistical analyses on the species richness of infauna #########################################

# ANOVA w/ unbalanced design (i.e., not equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_rich <- ggboxplot(richness, x = "site", y = "count", color = "infestation",
                             palette = my_colors)
boxplot_rich

interaction_rich <- ggline(richness, x = "site", y = "count", color = "infestation",
                              add = c("mean_se", "dotplot"),
                              palette = my_colors)
interaction_rich

hist(richness$count) # Normal looking distribution

# Type III ANOVA with unbalanced design
rich_aov <- aov(count ~ infestation * site, data = richness)
Anova(rich_aov, type = "III") # No significant interaction

# Type II ANOVA with unbalanced design
rich_aov1 <- aov(count ~ infestation + site, data = richness)
Anova(rich_aov1, type = "II")

plot(rich_aov1, 1)
leveneTest(count ~ infestation * site, data = richness) # Variances are homogeneous

plot(rich_aov1, 2)
anova_rich_residuals <- residuals(object = rich_aov1)
shapiro.test(x = anova_rich_residuals) # Residuals are normally distributed

# Pairwise comparisons
summary(glht(rich_aov1, linfct = mcp(site = "Tukey")))





###############################################################################################
# Species diversity w/ Shannon diversity (H') across all sites  ===============================
###############################################################################################

# Curate dataset ##############################################################################
library(vegan)

abund_wide2 <- abund_wide1
View(env1)
View(abund_wide2)

shannon <- diversity(abund_wide1[,6:120], "shannon")
View(shannon)
abund_wide2 <- abund_wide2 %>%
  add_column(shannon, .after = "quadrat")

shannon1 <- abund_wide1 %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon)
  )
View(shannon1)

# Barplot of the Shannon index of infauna across all sites (excluding OWR) ##################

# Barplot of Shannon index of infauna per quadrat, per site and per infestation
barplot_shannon <- ggplot(abund_wide1, aes(x = quadrat, y = shannon, 
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

# Barplot of Shannon index of infauna per site and per infestation
barplot_shannon1 <- ggplot(data = shannon1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = shannon_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon1

# Statistical analyses on the Shannon diversity of infauna #########################################

# ANOVA w/ unbalanced design (i.e., not equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_shannon <- ggboxplot(abund_wide1, x = "site", y = "shannon", color = "infestation",
                          palette = my_colors)
boxplot_shannon

interaction_shannon <- ggline(abund_wide1, x = "site", y = "shannon", color = "infestation",
                           add = c("mean_se", "dotplot"),
                           palette = my_colors)
interaction_shannon

hist(abund_wide1$shannon)

# Type III ANOVA with unbalanced design
shannon_aov <- aov(shannon ~ infestation * site, data = abund_wide1)
Anova(shannon_aov, type = "III") # No significant interaction

# Type II ANOVA with unbalanced design
shannon_aov1 <- aov(shannon ~ infestation + site, data = abund_wide1)
Anova(shannon_aov1, type = "II")

plot(shannon_aov1, 1)
leveneTest(shannon ~ infestation * site, data = abund_wide1) # Variances are homogeneous

plot(shannon_aov1, 2)
anova_shannon_residuals <- residuals(object = shannon_aov1)
shapiro.test(x = anova_shannon_residuals) # Residuals are normally distributed

# Pairwise comparisons
summary(glht(shannon_aov1, linfct = mcp(site = "Tukey")))




###############################################################################################
# Species diversity w/ Simpson diversity (λ) across all sites  ===============================
###############################################################################################

# Curate dataset ##############################################################################
View(abund_wide2)

simpson <- diversity(abund_wide1[,6:120], index="simpson")
abund_wide2 <- abund_wide2 %>%
  add_column(simpson, .after = "shannon")

simpson1 <- abund_wide1 %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson)
  )
View(simpson1)

# Barplot of the Simpson index of infauna across all sites (excluding OWR) ####################

# Barplot of Simpson index of infauna per quadrat, per site and per infestation
barplot_simpson <- ggplot(abund_wide1, aes(x = quadrat, y = simpson, 
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

# Barplot of Simpson index per site and per infestation
barplot_simpson1 <- ggplot(data = simpson1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = simpson_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_sd, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Simpson's diversity (λ) ")
barplot_simpson1

# Statistical analyses on the Simpson diversity of infauna #########################################

# ANOVA w/ unbalanced design (i.e., not equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_simpson <- ggboxplot(abund_wide1, x = "site", y = "simpson", color = "infestation",
                             palette = my_colors)
boxplot_simpson

interaction_simpson <- ggline(abund_wide1, x = "site", y = "simpson", color = "infestation",
                              add = c("mean_se", "dotplot"),
                              palette = my_colors)
interaction_simpson

hist(abund_wide1$simpson)

# Type III ANOVA with unbalanced design
simpson_aov <- aov(simpson ~ site * infestation, data = abund_wide1)
Anova(simpson_aov, type = "III") # No significant interaction

# Type II ANOVA with unbalanced design
simpson_aov1 <- aov(simpson ~ infestation + site, data = abund_wide1)
Anova(simpson_aov1, type = "II") # No significant interaction

plot(simpson_aov1, 1)
leveneTest(simpson ~ infestation * site, data = abund_wide1) # Variances are homogeneous

plot(simpson_aov1, 2)
anova_simpson_residuals <- residuals(object = simpson_aov1)
shapiro.test(x = anova_simpson_residuals) # Residuals are normally distributed

# Pairwise comparisons
summary(glht(simpson_aov1, linfct = mcp(site = "Tukey")))





###############################################################################################
# Pielou's evenness (J) across all sites  =====================================================
###############################################################################################

# Curate dataset ##############################################################################
View(abund_wide2)
View(richness)

# J = H'/ln(S)
abund_wide2 <- abund_wide2 %>%
  add_column(richness = richness$count, .after = "simpson")
pielou <- abund_wide2$shannon/log(abund_wide2$richness)
abund_wide2 <- abund_wide2 %>%
  add_column(pielou, .after = "richness")

pielou1 <- abund_wide2 %>%
  dplyr::group_by(site, infestation) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou)
  )
View(pielou1)

# Barplot of the Pielou's index of infauna across all sites (excluding OWR) ####################

# Barplot of the species evenness for each quadrat, site and infestation level
barplot_pielou <- ggplot(abund_wide2, aes(x = quadrat, y = pielou, 
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

# Barplot of the species evenness for each site and infestation level
barplot_pielou1 <- ggplot(data = pielou1, aes(x = site, fill = infestation)) + 
  geom_bar(aes(y = pielou_mean, color = infestation),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
#  annotate("text", 
#           x = c(1, 2, 3, 4, 5),
#           y = 1.4,
#           label = c("AB", "A", "AB", "AB", "B")) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = my_colors) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Infestation", fill = "Infestation",
       x = "Site", y = "Average Pielou's evenness (J)")
barplot_pielou1

# Statistical analyses on the Pielou's diversity of infauna #########################################

# ANOVA w/ unbalanced design (i.e., not equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_pielou <- ggboxplot(abund_wide2, x = "site", y = "pielou", color = "infestation",
                             palette = my_colors)
boxplot_pielou

interaction_pielou <- ggline(abund_wide2, x = "site", y = "pielou", color = "infestation",
                              add = c("mean_se", "dotplot"),
                              palette = my_colors)
interaction_pielou

hist(abund_wide2$pielou) # Normal-looking distribution

# Type III ANOVA with unbalanced design
pielou_aov <- aov(pielou ~ site * infestation, data = abund_wide2)
Anova(pielou_aov, type = "III") # No significant interaction

# Type II ANOVA with unbalanced design
pielou_aov1 <- aov(pielou ~ infestation + site, data = abund_wide2)
Anova(pielou_aov1, type = "II")

plot(pielou_aov1, 1)
leveneTest(pielou ~ infestation * site, data = abund_wide2) # Variances are homogeneous

plot(pielou_aov1, 2)
anova_pielou_residuals <- residuals(object = pielou_aov1)
shapiro.test(x = anova_pielou_residuals) # Residuals are normally distributed

# Pairwise comparisons
summary(glht(pielou_aov1, linfct = mcp(site = "Tukey")))


###############################################################################################
# Total abundance at Old Woman's River  =======================================================
###############################################################################################

# Curate dataset ####
abund_tot_owr <- abund_long_owr %>%
  na.omit(abund_long_owr) %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  summarize(ab_total = sum(abundance))
head(abund_tot_owr)

abund_tot_owr1 <- abund_tot_owr  %>%
  dplyr::group_by(infestation, lineage) %>%
  summarise(
    ab_tot_mean = mean(ab_total),
    ab_tot_sd = sd(ab_total)
  )
View(abund_tot_owr1)

# Barplot of the total abundance of infauna at OWR ####################

# Barplot of the total abundance in infauna per quadrat, infestation level and lineage at Old Woman's River
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

# Barplot of the total abundance in infauna per infestation level and lineage at Old Woman's River
barplot_abund_tot_owr1 <- ggplot(abund_tot_owr1, aes(x = infestation, y = ab_tot_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_sd, ymax = ab_tot_mean + ab_tot_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total abundance")
barplot_abund_tot_owr1

# Statistical analyses on the total abundance of infauna at OWR ###################################

# ANOVA w/ balanced design (i.e., equal sample sizes within levels of independent groupings)
table(abund_tot_owr$lineage, abund_tot_owr$infestation)

# Visualize data
boxplot_abundtot_owr <- ggboxplot(abund_tot_owr, x = "lineage", y = "ab_total", color = "infestation",
                                    palette = my_colors)
boxplot_abundtot_owr

interaction_abundtot_owr <- ggline(abund_tot_owr, x = "lineage", y = "ab_total", color = "infestation",
                                      add = c("mean_se", "dotplot"),
                                      palette = my_colors)
interaction_abundtot_owr

# ANOVA w/ interaction
anova_abundtot_owr <- aov(ab_total ~ infestation * lineage, data = abund_tot_owr)
Anova(anova_abundtot_owr) # Interaction is not significant

# ANOVA w/out interaction
anova_abundtot_owr1 <- aov(ab_total ~ infestation + lineage, data = abund_tot_owr)
Anova(anova_abundtot_owr1)

plot(anova_abundtot_owr1, 1)
leveneTest(ab_total ~ infestation * lineage, data = abund_tot_owr) # Variances are homogeneous

plot(anova_abundtot_owr1, 2)
anova_abundtotowr_residuals <- residuals(object = anova_abundtot_owr1)
shapiro.test(x = anova_abundtotowr_residuals) # Residuals are normally distributed




###############################################################################################
# Total biomass at OWR  =======================================================================
###############################################################################################

# Curate dataset ####
View(biom_long_owr)

biom_tot_owr <- biom_long_owr %>%
  na.omit(biom_long) %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  summarize(biom_total = sum(biomass))
head(biom_tot_owr)

biom_tot_owr1 <- biom_tot_owr  %>%
  dplyr::group_by(infestation, lineage) %>%
  summarise(
    biom_tot_mean = mean(biom_total),
    biom_tot_sd = sd(biom_total)
    )
View(biom_tot_owr1)

# Barplot of the total biomass of infauna at OWR ####

# Barplot of the total biomass of infauna for each quadrat, site and infestation level
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

# Barplot of the total biomass of infauna for each site and infestation level
barplot_biom_tot_owr1 <- ggplot(biom_tot_owr1, aes(x = infestation, y = biom_tot_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = biom_tot_mean - biom_tot_sd, ymax = biom_tot_mean + biom_tot_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average total biomass")
barplot_biom_tot_owr1

# Statistical analyses on the total biomass of infauna at OWR ###################################

# ANOVA w/ balanced design (i.e., equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_biomtot_owr <- ggboxplot(biom_tot_owr, x = "lineage", y = "biom_total", color = "infestation",
                                  palette = my_colors)
boxplot_biomtot_owr

interaction_biomtot_owr <- ggline(biom_tot_owr, x = "lineage", y = "biom_total", color = "infestation",
                                   add = c("mean_se", "dotplot"),
                                   palette = my_colors)
interaction_biomtot_owr

# ANOVA w/ interaction
anova_biomtot_owr <- aov(biom_total ~ infestation * lineage, data = biom_tot_owr)
Anova(anova_biomtot_owr) # Interaction is not significant

# ANOVA w/out interaction
anova_biomtot_owr1 <- aov(biom_total ~ infestation + lineage, data = biom_tot_owr)
Anova(anova_biomtot_owr1)

plot(anova_biomtot_owr1, 1)
leveneTest(biom_total ~ infestation * lineage, data = biom_tot_owr) # Variances are NOT homogeneous

plot(anova_biomtot_owr1, 2)
anova_biomtotowr_residuals <- residuals(object = anova_biomtot_owr1)
shapiro.test(x = anova_biomtotowr_residuals) # Residuals are normally distributed


###############################################################################################
# Species richness across at OWR  =============================================================
###############################################################################################

# Curate dataset ##############################################################################
View(abund_long_owr)

richness_owr <- abund_long_owr %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  summarise(count = n_distinct(species))
View(richness_owr)

richness_owr1 <- richness_owr %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
richness_owr1

# Barplot of the species richness at OWR ####

# Barplot of the species richness per quadrat, infestation level and lineage at Old Woman's River
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

# Barplot of the species richness per infestation level and lineage at Old Woman's River
barplot_rich_owr1 <- ggplot(richness_owr1, aes(x = infestation, y = count_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average species richness")
barplot_rich_owr1

# Statistical analyses on the species richness of infauna at OWR ###################################

# ANOVA w/ balanced design (i.e., equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_rich_owr <- ggboxplot(richness_owr, x = "lineage", y = "count", color = "infestation",
                                 palette = my_colors)
boxplot_rich_owr

interaction_rich_owr <- ggline(richness_owr, x = "lineage", y = "count", color = "infestation",
                                  add = c("mean_se", "dotplot"),
                                  palette = my_colors)
interaction_rich_owr

# ANOVA w/ interaction
anova_rich_owr <- aov(count ~ infestation * lineage, data = richness_owr)
Anova(anova_rich_owr) # Interaction is not significant

# ANOVA w/out interaction
anova_rich_owr1 <- aov(count ~ infestation + lineage, data = richness_owr)
Anova(anova_rich_owr1)

plot(anova_rich_owr1, 1)
leveneTest(count ~ infestation * lineage, data = richness_owr) # Variances are homogeneous

plot(anova_rich_owr1, 2)
anova_richtowr_residuals <- residuals(object = anova_rich_owr1)
shapiro.test(x = anova_richtowr_residuals) # Residuals are normally distributed


###############################################################################################
# Species diversity w/ Shannon diversity (H') at OWR  =========================================
###############################################################################################

# Curate dataset ##############################################################################
abund_wide_owr1 <- abund_wide_owr
View(abund_wide_owr1)

shannon_owr <- diversity(abund_wide_owr[,6:75], "shannon")
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(shannon = shannon_owr, .after = "quadrat")

shannon_owr1 <- abund_wide_owr1 %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon)
  )
View(shannon_owr1)

# Barplot of the species diversity at OWR ####

# Barplot of the species diversity for each quadrat, lineage and infestation level
barplot_shannon_owr <- ggplot(abund_wide_owr1, aes(x = quadrat, y = shannon, 
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

# Barplot of the species diversity per infestation and lineage
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

# Statistical analyses on the Shannon index of infauna at OWR ###################################

# ANOVA w/ balanced design (i.e., equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_shannon_owr <- ggboxplot(abund_wide_owr1, x = "lineage", y = "shannon", color = "infestation",
                              palette = my_colors)
boxplot_shannon_owr

interaction_shannon_owr <- ggline(abund_wide_owr1, x = "lineage", y = "shannon", color = "infestation",
                               add = c("mean_se", "dotplot"),
                               palette = my_colors)
interaction_shannon_owr

# ANOVA w/ interaction
anova_shannon_owr <- aov(shannon ~ infestation * lineage, data = abund_wide_owr1)
Anova(anova_shannon_owr) # Interaction is not significant

# ANOVA w/out interaction
anova_shannon_owr1 <- aov(shannon ~ infestation + lineage, data = abund_wide_owr1)
Anova(anova_shannon_owr1)

plot(anova_shannon_owr1, 1)
leveneTest(shannon ~ infestation * lineage, data = abund_wide_owr1) # Variances are homogeneous

plot(anova_shannon_owr1, 2)
anova_shannonowr_residuals <- residuals(object = anova_shannon_owr1)
shapiro.test(x = anova_shannonowr_residuals) # Residuals are normally distributed


###############################################################################################
# Species diversity w/ Simpson diversity (λ) at OWR  =========================================
###############################################################################################

# Curate dataset ##############################################################################
View(abund_wide_owr1)

simpson_owr <- diversity(abund_wide_owr[,6:75], "simpson")
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(simpson = simpson_owr, .after = "shannon")

simpson_owr1 <- abund_wide_owr1 %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson)
  )
View(simpson_owr1)

# Barplot of the Simpson's index at OWR ####

# Barplot of the species diversity per quadrat, infestation level and lineage at Old Woman's River
barplot_simpson_owr <- ggplot(abund_wide_owr1, aes(x = quadrat, y = simpson, 
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

# Barplot of the species diversity per infestation level and lineage at Old Woman's River
barplot_simpson_owr1 <- ggplot(simpson_owr1, aes(x = infestation, y = simpson_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_sd, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Simpson's diversity (λ)")
barplot_simpson_owr1

# Statistical analyses on the Simpson index of infauna at OWR ###################################

# ANOVA w/ balanced design (i.e., equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_simpson_owr <- ggboxplot(abund_wide_owr1, x = "lineage", y = "simpson", color = "infestation",
                                 palette = my_colors)
boxplot_simpson_owr

interaction_simpson_owr <- ggline(abund_wide_owr1, x = "lineage", y = "simpson", color = "infestation",
                                  add = c("mean_se", "dotplot"),
                                  palette = my_colors)
interaction_simpson_owr

# ANOVA w/ interaction
anova_simpson_owr <- aov(simpson ~ infestation * lineage, data = abund_wide_owr1)
Anova(anova_simpson_owr) # Interaction is not significant

# ANOVA w/out interaction
anova_simpson_owr1 <- aov(simpson ~ infestation + lineage, data = abund_wide_owr1)
Anova(anova_simpson_owr1)

plot(anova_simpson_owr1, 1)
leveneTest(simpson ~ infestation * lineage, data = abund_wide_owr1) # Variances are homogeneous

plot(anova_simpson_owr1, 2)
anova_simpsonowr_residuals <- residuals(object = anova_simpson_owr1)
shapiro.test(x = anova_simpsonowr_residuals) # Residuals are normally distributed


###############################################################################################
# Pielou's evenness (J) at OWR  ===============================================================
###############################################################################################

# Curate dataset ##############################################################################
View(abund_wide_owr1)
View(richness_owr)

# J = H'/ln(S)
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(richness = richness_owr$count, .after = "simpson")
pielou_owr <- abund_wide_owr1$shannon/log(abund_wide_owr1$richness)
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(pielou = pielou_owr, .after = "richness")

pielou_owr1 <- abund_wide_owr1 %>%
  dplyr::group_by(lineage, infestation) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou)
  )
View(pielou_owr1)

# Barplot of the Pielou's index of infauna at OWR ####################

# Barplot of the species evenness per quadrat, infestation level and lineage at Old Woman's River
barplot_pielou_owr <- ggplot(abund_wide_owr1, aes(x = quadrat, y = pielou, 
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

# Barplot of the species evenness per infestation level and lineage at Old Woman's River
barplot_pielou_owr1 <- ggplot(pielou_owr1, aes(x = infestation, y = pielou_mean, fill = infestation)) + 
  geom_bar(aes(color = infestation), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  scale_color_manual(values = my_colors) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Infestation", y = "Average Pielou's evenness (J)")
barplot_pielou_owr1

# Statistical analyses on the Pielou index of infauna at OWR ###################################

# ANOVA w/ balanced design (i.e., equal sample sizes within levels of independent groupings)

# Visualize data
boxplot_pielou_owr <- ggboxplot(abund_wide_owr1, x = "lineage", y = "pielou", color = "infestation",
                                 palette = my_colors)
boxplot_pielou_owr

interaction_pielou_owr <- ggline(abund_wide_owr1, x = "lineage", y = "pielou", color = "infestation",
                                  add = c("mean_se", "dotplot"),
                                  palette = my_colors)
interaction_pielou_owr

# ANOVA w/ interaction
anova_pielou_owr <- aov(pielou ~ infestation * lineage, data = abund_wide_owr1)
Anova(anova_pielou_owr) # Interaction is not significant

# ANOVA w/out interaction
anova_pielou_owr1 <- aov(pielou ~ infestation + lineage, data = abund_wide_owr1)
Anova(anova_pielou_owr1)

plot(anova_pielou_owr1, 1)
leveneTest(pielou ~ infestation * lineage, data = abund_wide_owr1) # Variances are homogeneous

plot(anova_pielou_owr1, 2)
anova_pielouowr_residuals <- residuals(object = anova_pielou_owr1)
shapiro.test(x = anova_pielouowr_residuals) # Residuals are normally distributed


