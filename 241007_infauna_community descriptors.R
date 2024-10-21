## ---
##
## Script name: General descriptors of macrofaunal communities associated with manipulated mussel beds depending on their level
##              of euendolithic corrosion, either across a biogeographical gradient or between genetic lineages
##
## Purpose of script: 
##    - EUENDOLITHS x BIOGEOGRAPHY : Compare the following variables between corrosion status and among transplant sites
##      (excluding Old Woman's River) 
##        - Unique vs Common species
##        - Total macrofaunal abundance
##        - Total macrofaunal biomass
##        - Species richness (S)
##        - Species diversity : Shannon-Wiener (H')
##        - Species diversity : Simpson
##        - Pielou's evenness (J)
##    - EUENDOLITHS x PERNA LINEAGES : Compare the above variables between corrosion status and Perna lineages at Old Woman's River
##
## Find in the manuscript :
##    - EUENDOLITHS x BIOGEOGRAPHY
##        - Figure 2b : Barplot of species richness between corroded and non-corroded mussel beds for each site (except Old Woman's River)
##        - Figure 2c : Barplot of Shannon diversity index between corroded and non-corroded mussel beds for each site (except Old Woman's River)
##        - Figure 2d : Barplot of Pielou's evenness between corroded and non-corroded mussel beds for each site (except Old Woman's River)
##        - Table S6 : Results of the series of ANOVAs on the Gamma GLMs without interactions
##            - Table S6a : Total macrofaunal abundance --> Significant effect of site
##            - Table S6b : Total macrofaunal biomass --> Significant effect of site
##        - Table S7 : Pairwise comparisons (when significant effect detected) on
##            - Table S7a : Total macrofaunal abundance --> BR = JB > MB = PEd
##            - Table S7b :  Total macrofaunal biomass --> BR > MB = PEd and JB > PEd
##        - Table S8 : Results of the series of ANOVAs on
##            - Table S8a : Species richness (S)
##            - Table S8b : Shannon-Wiener index (H') --> Significant effect of site
##            - Table S8c : Simpson's diversity index (λ) --> Significant effect of site
##            - Table S8d : Pielou's evenness (J) --> Significant effect of site
##        - Table S9 : Pairwise comparisons (when significant effect detected) on
##            - Table S9a : Shannon-Wiener index (H') --> PEd > MB
##            - Table S9b : Simpson's diversity index (λ) --> PEd > MB + JB
##            - Table S9c : Pielou's evenness (J) --> PEd > all
##    - EUENDOLITHS x PERNA LINEAGES
##        - Table S10 : Results of the series of ANOVAs on 
##            - Table S10a : Total macrofaunal abundance 
##            - Table S10b : Total macrofaunal biomass
##            - Table S10c : Species richness (S)
##            - Table S10d : Shannon-Wiener index (H')
##            - Table S10e : Simpson's diversity index (λ)
##            - Table S10f : Pielou's evenness (J)
##
## Author: Alexia Dievart
##
## Date Created: 2024-02-01
## Dates Updated: 2024-10-07; 2024-10-16
##
## Copyright (c) Alexia DIEVART 2024
## Email: alexia.dievart@hotmail.fr
##
## ---
##
## Notes: 
##   - Quadrats cannot be specified as a random factor (n < 5).
##   - Analyses conducted ACROSS ALL SITES excluded Old Woman's River, as the effect of site and Perna lineages
##     are intricately confounded. 
## ---

##
# SECTION: Session setup ----
##

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

# Create vectors for graphical visualization
col_corr <- c("darkgrey", "brown3")
site_order <- c("Mosselbaai", "Brenton-on-Sea", "Jeffreysbaai", "Port Edward")
sites <- c("MB", "BR", "JB", "PEd")
inf <- c("C", "NC")





##
# SECTION: Load data ----
##

# Environmental variables - architectural complexity ====
arch <- read.csv("./RAW DATA/infauna_architectural complexity.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(arch)

factor1 <- c("site", "corrosion", "lineage", "quadrat")
arch[,factor1] <- lapply(arch[,factor1], factor) # Transform relevant columns in factor
arch <- arch %>% drop_na(live_perna) # Drop the missing quadrats at collection

View(arch)



# Environmental variables - number of byssal threads ====
byssus <- read.csv("./RAW DATA/infauna_byssal.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(byssus)

byssus[,factor1] <- lapply(byssus[,factor1], factor)
byssus <- na.omit(byssus) # Drop the quadrats for which byssal threads were not counted
View(byssus)




# Community variables - e.g. abundance, biomass ====
community <- read.csv("./RAW DATA/infauna_community.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(community)

factor2 <- c("site", "bioregion", "corrosion", "lineage", "quadrat", "higher_group", "species", "code_sp")
community[,factor2] <- lapply(community[,factor2], factor)

# Drop the lines without count (abundance = NA)
community <- community %>% drop_na(abundance)
View(community)





##
# SECTION: Macrofaunal communities ACROSS ALL SITES ----
##

# Curate data sets ====

# Unique and common species ####
abund_long <- community[,c(1:5,7:9)] #  Select only the relevant variables in 'community' data set
abund_long1 <- abund_long[abund_long$site != "Old Woman's River",] # Exclude OWR from the data set
View(abund_long1)

# Total number of taxa across all sites (including OWR) = 115
abund_long %>% summarise(n_distinct(species))
# Including 18 taxa identified to higher taxonomic levels and 7 taxa identified to the genus level.

# Number of unique taxa per corrosion level (excluding OWR)
abund_long1 %>%
  group_by(corrosion) %>%
  summarise(n_distinct(species))
# There are 80 unique taxa in corroded mussel beds, and 80 unique taxa in non-corroded mussel beds.

# Number and identity of common taxa between corrosion levels (excluding OWR)
Reduce(intersect, split(abund_long1$species, abund_long1$corrosion))
# There are 60 common taxa between corroded and non-corroded mussel beds.

# I would love to add a bit of script here that would give me a list of the unique species present in the
# corroded and a list of unique species present in the non-corroded mussel beds, without the common species.



# Total macrofaunal abundance ####
abund_tot <- abund_long1 %>%
  dplyr::group_by(site, corrosion, quadrat) %>%
  summarize(tot_abund = sum(abundance))
View(abund_tot)

abund_tot1 <- abund_tot  %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    mean_abund = mean(tot_abund),
    sd_abund = sd(tot_abund))
View(abund_tot1)



# Total macrofaunal biomass ####
biom_long <- community[,c(1:5,7:8, 10)] #  Select only the relevant variables in 'community' data set
biom_long1 <- biom_long[biom_long$site != "Old Woman's River",] # Exclude OWR from the data set
View(biom_long1)

biom_tot <- biom_long1 %>%
  na.omit(biom_long) %>%
  dplyr::group_by(site, corrosion, quadrat) %>%
  summarize(tot_biom = sum(biomass))
View(biom_tot)

biom_tot1 <- biom_tot %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    mean_biom = mean(tot_biom),
    sd_biom = sd(tot_biom)
  )
View(biom_tot1)



# Species richness (S) ####
richness <- abund_long1 %>%
  dplyr::group_by(site, corrosion, quadrat) %>%
  summarise(count = n_distinct(species))
View(richness)

richness1 <- richness %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count))
View(richness1)



# Species diversity w/ Shannon diversity (H') ####
library(vegan)

abund_long2 <- abund_long1[,-6] # Transform the 'abund_long1' dataset into a wide format
abund_wide <- spread(abund_long2, code_sp, abundance)
abund_wide[is.na(abund_wide)] <- 0
View(abund_wide)

abund_wide1 <- abund_wide # Create a new wide-format dataset to add the community indexes
View(arch); View(abund_wide1)

shannon <- diversity(abund_wide[,6:104], "shannon") # Calculate Shannon's diversity index (H')
abund_wide1 <- abund_wide1 %>%
  add_column(shannon, .after = "quadrat")
View(abund_wide1)

shannon1 <- abund_wide1 %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon)
  )
View(shannon1)



# Species diversity w/ Simpson diversity (λ) ####
simpson <- diversity(abund_wide[,6:104], index = "simpson")
abund_wide1 <- abund_wide1 %>%
  add_column(simpson, .after = "shannon")
View(abund_wide1)

simpson1 <- abund_wide1 %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson)
  )
View(simpson1)



# Pielou's evenness (J) across all sites (excluding OWR) ####
View(abund_wide1) ; View(richness)

abund_wide1 <- abund_wide1 %>%
  add_column(richness = richness$count, .after = "simpson")
pielou <- abund_wide1$shannon/log(abund_wide1$richness) # Formula : J = H'/ln(S)
abund_wide1 <- abund_wide1 %>%
  add_column(pielou, .after = "richness")
View(abund_wide1)

pielou1 <- abund_wide1 %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou)
  )
View(pielou1)





# Visualize data sets ====

# GRAPH : Total macrofaunal abundance ####

# Barplot of the total macrofaunal abundance per quadrat, site and corrosion levels
barplot_abund_tot <- ggplot(abund_tot, aes(x = quadrat, y = tot_abund, 
                                           color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Total abundance")
barplot_abund_tot
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).
# Quadrats that could introduce biases : PEd-E-NI-2 and PEd-E-NI-4 (loss of mussels due to wave action).

# Barplot of the total macrofaunal abundance per site and corrosion levels
barplot_abund_tot1 <- ggplot(data = abund_tot1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = mean_abund, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average total abundance")
barplot_abund_tot1



# GRAPH : Total macrofaunal biomass ####

# Barplot of the total macrofaunal biomass per quadrat, site and corrosion levels
barplot_biom_tot <- ggplot(biom_tot, aes(x = quadrat, y = tot_biom, 
                                         color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Total biomass")
barplot_biom_tot
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).
# Quadrats that could introduce biases : PEd-E-I-1 and PEd-E-NI-4 (loss of mussels due to wave action)

# Barplot of the total macrofaunal biomass per site and corrosion level
barplot_biom_tot1 <- ggplot(data = biom_tot1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = mean_biom, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom - sd_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average total biomass (mg)")
barplot_biom_tot1



# GRAPH : Species richness (S)

# Barplot of macrofaunal species richness per quadrat, site and corrosion level
barplot_rich <- ggplot(richness, aes(x = quadrat, y = count, 
                                     color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Species richness")
barplot_rich
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).

# Barplot of macrofaunal species richness per site and corrosion levels
barplot_rich1 <- ggplot(data = richness1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = count_mean, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5, 4.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average species richness")
barplot_rich1



# GRAPH : Shannon-Wiener diversity index (H') ####

# Barplot of macrofaunal Shannon index per quadrat, site and corrosion level
barplot_shannon <- ggplot(abund_wide1, aes(x = quadrat, y = shannon, 
                                           color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Shannon-Wiener diversity (H')")
barplot_shannon
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).

# Barplot of macrofaunal Shannon index per site and corrosion levels
barplot_shannon1 <- ggplot(data = shannon1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = shannon_mean, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon1



# GRAPH : Simpson's diversity index (λ) ####

# Barplot of macrofaunal Simpson index per quadrat, site and corrosion level
barplot_simpson <- ggplot(abund_wide1, aes(x = quadrat, y = simpson, 
                                           color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Simpson's diversity (λ) ")
barplot_simpson
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).

# Barplot of macrofaunal Simpson index per site and corrosion levels
barplot_simpson1 <- ggplot(data = simpson1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = simpson_mean, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_sd, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average Simpson's diversity (λ) ")
barplot_simpson1



# GRAPH : Pielou's evenness (J) ####

# Barplot of the species evenness per quadrat, site and corrosion level
barplot_pielou <- ggplot(abund_wide1, aes(x = quadrat, y = pielou, 
                                          color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, nrow = 2) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).

# Barplot of the species evenness per site and corrosion levels
barplot_pielou1 <- ggplot(data = pielou1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = pielou_mean, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average Pielou's evenness (J)")
barplot_pielou1





##
# STATS : Total macrofaunal abundance ====
##

# Assumptions ####
table(abund_tot$site, abund_tot$corrosion)
# Data does not have equal sample sizes within levels of independent groupings.

# Visualize data
hist(abund_tot$tot_abund)

boxplot_tot_abund <- ggboxplot(abund_tot, x = "site", y = "tot_abund", color = "corrosion",
                               palette = col_corr)
boxplot_tot_abund

interaction_tot_abund <- ggline(abund_tot, x = "site", y = "tot_abund", color = "corrosion",
                                add = c("mean_se", "dotplot"),
                                palette = col_corr)
interaction_tot_abund


# Classic ANOVAs with outliers ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_tot_abund <- aov(tot_abund ~ corrosion * site, data = abund_tot)
Anova(anova_tot_abund, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_tot_abund, 1) # Points 4, 6 and 7 detected as outliers.
leveneTest(tot_abund ~ corrosion * site, data = abund_tot) # Variances are homogeneous

plot(anova_tot_abund, 2) # Points fall approximately along the line
anova_tot_abund_residuals <- residuals(object = anova_tot_abund) 
shapiro.test(x = anova_tot_abund_residuals) # Residuals are NOT normally distributed (p = 0.01398).



# Type II ANOVA without interaction --> NOT SELECTED
anova_tot_abund1 <- aov(tot_abund ~ corrosion + site, data = abund_tot)
Anova(anova_tot_abund1, type = "II") # There is a significant effect of site (p < 0.05) but NO effect of corrosion (p = 0.5338).

plot(anova_tot_abund1, 1) # Points 4, 6 and 7 detected as outliers.
leveneTest(tot_abund ~ corrosion * site, data = abund_tot)

plot(anova_tot_abund1, 2) # Points 4, 6 and 7 detected as outliers.
anova_tot_abund1_residuals <- residuals(object = anova_tot_abund1)
shapiro.test(x = anova_tot_abund1_residuals) # Residuals are NOT normally distributed (p = 0.03627).


# ANOVAs on models ####

# ANOVA w/ interaction using LM --> NOT SELECTED
lm_tot_abund <- lm(tot_abund ~ site * corrosion, data = abund_tot)
Anova(lm_tot_abund, type = "III") # Interaction is not significant.

summary(lm_tot_abund)
summary(lm_tot_abund)$r.squared # Linear model only explains 70 % of the variance.

DHARMa::simulateResiduals(fittedModel = lm_tot_abund, plot = T) # LM with interaction is not a bad fit. 

AIC(lm_tot_abund) # AIC(lm) = 365.1747



# ANOVA without interaction using LM --> NOT SELECTED
lm1_tot_abund <- lm(tot_abund ~ site + corrosion, data = abund_tot)
Anova(lm1_tot_abund, type = "III") # There is a significant effect of site (p < 0.05) but NO effect of corrosion (p = 0.5338). 

summary(lm1_tot_abund)
summary(lm1_tot_abund)$r.squared # Linear model only explains 64 % of the variance.

DHARMa::simulateResiduals(fittedModel = lm1_tot_abund, plot = T) # LM without interaction is not a good fit. 

AIC(lm1_tot_abund) # AIC(lm1) = 365.9085



# ANOVA w/ interaction using GLM (family = "Gamma") --> NOT SELECTED
glm_tot_abund <- glm(tot_abund ~ site * corrosion, data = abund_tot,
                     family = "Gamma")
Anova(glm_tot_abund, type = "III") # Interaction is not significant.

summary(glm_tot_abund)

DHARMa::simulateResiduals(fittedModel = glm_tot_abund, plot = T) # GLM is a good fit, except for the KS test. 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = glm_tot_abund, plot = F) 
plotResiduals(simulationOutput, abund_tot$corrosion)  # No significant within-group deviation.

AIC(glm_tot_abund) # AIC(glm) = 346.9155



# ANOVA w/out interaction using GLM (family = "Gamma") --> Table S6a
glm_tot_abund1 <- glm(tot_abund ~ corrosion + site, data = abund_tot,
                      family = "Gamma")
Anova(glm_tot_abund1) # There is a significant effect of site (p < 0.05) but NO effect of corrosion (p = 0.4988).

summary(glm_tot_abund1)

DHARMa::simulateResiduals(fittedModel = glm_tot_abund1, plot = T) # GLM w/out interaction is a great fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = glm_tot_abund1, plot = F) 
plotResiduals(simulationOutput, abund_tot$corrosion)  # within-group deviations not significant

AIC(glm_tot_abund1) # AIC(glm) = 344.2506



# Pairwise comparisons of total macrofaunal abundance between sites --> Table S7a
summary(glht(glm_tot_abund1, linfct = mcp(site = "Tukey")))





##
# STATS : Total macrofaunal biomass ====
##

# Assumptions #### 

# Frequency table
table(biom_tot$site, biom_tot$corrosion)
# Data does not have equal sample sizes within levels of independent groupings

# Visualize data
hist(biom_tot$tot_biom)

boxplot_tot_biom <- ggboxplot(biom_tot, x = "site", y = "tot_biom", color = "corrosion",
                              palette = col_corr)
boxplot_tot_biom

interaction_tot_biom <- ggline(biom_tot, x = "site", y = "tot_biom", color = "corrosion",
                               add = c("mean_se", "dotplot"),
                               palette = col_corr)
interaction_tot_biom



# Classic ANOVAs with outliers ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_tot_biom <- aov(tot_biom ~ corrosion * site, data = biom_tot)
Anova(anova_tot_biom, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_tot_biom, 1) # Points 1, 2 and 3 detected as outliers.
leveneTest(tot_biom ~ corrosion * site, data = biom_tot) # Variances are homogeneous (p = 0.487).

plot(anova_tot_biom, 2) # Points 1, 2 and 3 detected as outliers.
anova_tot_biom_residuals <- residuals(object = anova_tot_biom)
shapiro.test(x = anova_tot_biom_residuals) # Residuals are NOT normally distributed (p = 0.0117). 



# Type II ANOVA without interaction --> NOT SELECTED
anova_tot_biom1 <- aov(tot_biom ~ corrosion + site, data = biom_tot)
Anova(anova_tot_biom1, type = "II") # Significant effect of site (p < 0.05) but NO effect of corrosion (p = 0.619056).

plot(anova_tot_biom1, 1) # Points 2, 3, and 6 detected as outliers.
leveneTest(tot_biom ~ corrosion * site, data = biom_tot)

plot(anova_tot_biom1, 2) # Points 2, 3 and 4 detected as outliers.
anova_tot_biom1_residuals <- residuals(object = anova_tot_biom1)
shapiro.test(x = anova_tot_biom1_residuals) # Residuals are NOT normally distributed (p = 0.001)



# ANOVAs on models ####

# ANOVA w/ interaction using LM --> NOT SELECTED
lm_tot_biom <- lm(tot_biom ~ site * corrosion, data = biom_tot)
Anova(lm_tot_biom, type = "III") # Interaction is not significant.

summary(lm_tot_biom)
summary(lm_tot_biom)$r.squared # Linear model only explains 50 % of the variance.

DHARMa::simulateResiduals(fittedModel = lm_tot_biom, plot = T) # LM is not a good fit.

AIC(lm_tot_biom) # AIC(lm) = 543.7281



# ANOVA w/out interaction using LM --> NOT SELECTED
lm_tot_biom1 <- lm(tot_biom ~ site + corrosion, data = biom_tot)
Anova(lm_tot_biom1, type = "II") # There is a significant effect of site (p = 0.0015) but NO effect of corrosion (p = 0.619056).

summary(lm_tot_biom1)
summary(lm_tot_biom1)$r.squared # Linear model only explains 48 % of the variance.

DHARMa::simulateResiduals(fittedModel = lm_tot_biom1, plot = T) # Although, LM is a good fit.

AIC(lm_tot_biom1) # AIC(lm) = 538.7832



# ANOVA w/ interaction using GLM (family = "Gamma") --> NOT SELECTED
glm_tot_biom <- glm(tot_biom ~ site * corrosion, data = biom_tot,
                    family = "Gamma")
Anova(glm_tot_biom, type = "III") # Interaction is not significant.

summary(glm_tot_biom)

DHARMa::simulateResiduals(fittedModel = glm_tot_biom, plot = T) # GLM is a good fit. 
simulationOutput <- DHARMa::simulateResiduals(fittedModel = glm_tot_biom, plot = F) 
plotResiduals(simulationOutput, biom_tot$corrosion)  # No significant within-group deviation.

AIC(glm_tot_biom) # AIC(glm) = 512.8858



# ANOVA w/out interaction using GLM (family = "Gamma") --> Table S6b
glm_tot_biom1 <- glm(tot_biom ~  corrosion + site, data = biom_tot,
                     family = "Gamma")
Anova(glm_tot_biom1, type = "II") # There is a significant effect of site (p < 0.05) but NO effect of corrosion (p = 0.6215).

summary(glm_tot_biom1)

DHARMa::simulateResiduals(fittedModel = glm_tot_biom1, plot = T) # GLM is a very good fit.
simulationOutput <- DHARMa::simulateResiduals(fittedModel = glm_tot_biom1, plot = F) 
plotResiduals(simulationOutput, biom_tot$corrosion)  # Within-group deviations not significant

AIC(glm_tot_biom1) # AIC(glm1) = 507.0966



# Pairwise comparisons of total macrofaunal biomass between sites --> Table S7b
summary(glht(glm_tot_biom1, linfct = mcp(site = "Tukey")))





##
# STATS : Species richness (S) ====
##

# Assumptions ####

# Frequency table
table(richness$site, richness$corrosion)
# Data does not have equal sample sizes within levels of independent groupings.

# Visualize data
hist(richness$count) # Normal-looking distribution

boxplot_rich <- ggboxplot(richness, x = "site", y = "count", color = "corrosion",
                          palette = col_corr)
boxplot_rich

interaction_rich <- ggline(richness, x = "site", y = "count", color = "corrosion",
                           add = c("mean_se", "dotplot"),
                           palette = col_corr)
interaction_rich




# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_rich <- aov(count ~ corrosion * site, data = richness)
Anova(anova_rich, type = "III") # Interaction is not significant.

plot(anova_rich, 1) # Points 13, 19 and 21 detected as outliers.
leveneTest(count ~ corrosion * site, data = richness) # Variances are homogeneous (p = 0.8264).

plot(anova_rich, 2) # Points fall approximately along the line.
anova_rich_residuals <- residuals(object = anova_rich) 
shapiro.test(x = anova_rich_residuals) # Residuals are normally distributed (p = 0.9074).



# Type II ANOVA without interaction --> Table S8a
anova_rich1 <- aov(count ~ corrosion + site, data = richness)
Anova(anova_rich1, type = "II") # There is NO significant effect of site (p = 0.09736) or corrosion (p = 0.87476).

plot(anova_rich1, 1) # Points 13, 19 and 21 detected as outliers.
leveneTest(count ~ corrosion * site, data = richness)

plot(anova_rich1, 2)  # Points fall approximately along the line.
anova_rich1_residuals <- residuals(object = anova_rich1)
shapiro.test(x = anova_rich1_residuals) # Residuals are normally distributed (p = 0.868).



# Pairwise comparisons of species richness between sites
summary(glht(anova_rich1, linfct = mcp(site = "Tukey")))





##
# STATS : Species diversity w/ Shannon diversity (H') ====
##

# Assumptions ####

# Frequency table
table(abund_wide1$site, abund_wide1$corrosion)
# Data does not have equal sample sizes within levels of independent groupings.

# Visualize data
hist(abund_wide1$shannon)

boxplot_shannon <- ggboxplot(abund_wide1, x = "site", y = "shannon", color = "corrosion",
                             palette = col_corr)
boxplot_shannon

interaction_shannon <- ggline(abund_wide1, x = "site", y = "shannon", color = "corrosion",
                              add = c("mean_se", "dotplot"),
                              palette = col_corr)
interaction_shannon



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_shannon <- aov(shannon ~ corrosion * site, data = abund_wide1)
Anova(anova_shannon, type = "III") # Interaction is not significant. 

plot(anova_shannon, 1) # Points 3, 9 and 21 detected as outliers.
leveneTest(shannon ~ corrosion * site, data = abund_wide1) # Variances are homogeneous (p = 0.4076).

plot(anova_shannon, 2) # Points fall approximately along the line.
anova_shannon_residuals <- residuals(object = anova_shannon) 
shapiro.test(x = anova_shannon_residuals) # Residuals are normally distributed (p = 0.747).



# Type II ANOVA without interaction --> Table S8b
anova_shannon1 <- aov(shannon ~ corrosion + site, data = abund_wide1)
Anova(anova_shannon1, type = "II") # Significant effect of site (p = 0.01978) but NO effect of corrosion (p = 0.46455).

plot(anova_shannon1, 1) # Points 8, 19 and 21 detected as outliers.
leveneTest(shannon ~ corrosion * site, data = abund_wide1)

plot(anova_shannon1, 2)  # Points fall approximately along the line.
anova_shannon1_residuals <- residuals(object = anova_shannon1)
shapiro.test(x = anova_shannon1_residuals) # Residuals are normally distributed (p = 0.9007).



# Pairwise comparisons of Shannon diversity between sites --> Table S9a
summary(glht(anova_shannon1, linfct = mcp(site = "Tukey")))






##
# STATS : Simpson diversity (λ) ====
##

# Assumptions ####

# Frequency table
table(abund_wide1$site, abund_wide1$corrosion)
# Data does not have equal sample sizes within levels of independent groupings.

# Visualize data
hist(abund_wide1$simpson) # Data distribution is left skewed.

boxplot_simpson <- ggboxplot(abund_wide1, x = "site", y = "simpson", color = "corrosion",
                             palette = col_corr)
boxplot_simpson

interaction_simpson <- ggline(abund_wide1, x = "site", y = "simpson", color = "corrosion",
                              add = c("mean_se", "dotplot"),
                              palette = col_corr)
interaction_simpson



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_simpson <- aov(simpson ~ corrosion * site, data = abund_wide1)
Anova(anova_simpson, type = "III") # Interaction is NOT significant.

plot(anova_simpson, 1) # Points 3, 8 and 9 detected as outliers.
leveneTest(simpson ~ corrosion * site, data = abund_wide1) # Variances are homogeneous (p = 0.4405).

plot(anova_simpson, 2) # Points fall approximately along the line.
anova_simpson_residuals <- residuals(object = anova_simpson) 
shapiro.test(x = anova_simpson_residuals) # Residuals are normally distributed (p = 0.9921).

# Type II ANOVA without interaction --> Table S8c
anova_simpson1 <- aov(simpson ~ corrosion + site, data = abund_wide1)
Anova(anova_simpson1, type = "II") # Significant effect of site (p = 0.01387) but NO effect of corrosion (p = 0.18014).

plot(anova_simpson1, 1) # Points 8, 19 and 21 detected as outliers.
leveneTest(simpson ~ corrosion * site, data = abund_wide1)

plot(anova_simpson1, 2)  # Points fall approximately along the line.
anova_simpson1_residuals <- residuals(object = anova_simpson1)
shapiro.test(x = anova_simpson1_residuals) # Residuals are normally distributed (p = 0.989). 



# Pairwise comparisons of Simpson's diversity between sites --> Table S9b
summary(glht(anova_simpson1, linfct = mcp(site = "Tukey")))





##
# STATS : Pielou's evenness (J) ====
##

# Assumptions ####

# Frequency table
table(abund_wide1$site, abund_wide1$corrosion)
# Data does not have equal sample sizes within levels of independent groupings

# Visualize data
hist(abund_wide1$pielou) # Normal-looking distribution

boxplot_pielou <- ggboxplot(abund_wide1, x = "site", y = "pielou", color = "corrosion",
                            palette = col_corr)
boxplot_pielou

interaction_pielou <- ggline(abund_wide1, x = "site", y = "pielou", color = "corrosion",
                             add = c("mean_se", "dotplot"),
                             palette = col_corr)
interaction_pielou



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_pielou <- aov(pielou ~ corrosion * site, data = abund_wide1)
Anova(anova_pielou, type = "III") # Interaction is NOT significant.

plot(anova_pielou, 1) # Points 3, 4 and 9 detected as outliers.
leveneTest(pielou ~ corrosion * site, data = abund_wide1) # Variances are homogeneous (p = 0.4394).

plot(anova_pielou, 2) # Points fall approximately along the line.
anova_pielou_residuals <- residuals(object = anova_pielou) 
shapiro.test(x = anova_pielou_residuals) # Residuals are normally distributed (p = 0.4629).

# Type II ANOVA without interaction --> Table S8d
anova_pielou1 <- aov(pielou ~ corrosion + site, data = abund_wide1)
Anova(anova_pielou1, type = "II") # Significant effect of site (p = 0.001044) but NO effect of corrosion (p = 0.259934).

plot(anova_pielou1, 1) # Points 4, 12 and 25 detected as outliers.
leveneTest(pielou ~ corrosion * site, data = abund_wide1)

plot(anova_pielou1, 2)  # Points fall approximately along the line.
anova_pielou1_residuals <- residuals(object = anova_pielou1)
shapiro.test(x = anova_pielou1_residuals) # Residuals are normally distributed (p = 0.1279).



# Pairwise comparisons of Pielou's evenness between sites
summary(glht(anova_pielou1, linfct = mcp(site = "Tukey")))








##
# SECTION: Macrofaunal communities AT OWR ----
##

# Curate data sets ====

# Unique and common species ####

# Create a data set for Old Woman's River
abund_long_owr <- abund_long[abund_long$site == "Old Woman's River",]
View(abund_long_owr)

# Number of unique taxa per corrosion level and Perna lineage
abund_long_owr %>%
  group_by(corrosion, lineage) %>%
  summarise(n_distinct(species))
# Results : I-E = 39 (+ Opisthobranchia), I-W = 50, NI-E = 53 and NI-W = 43.

# Number and identity of common taxa between corrosion levels
abund_long_owr <- abund_long_owr %>%
  dplyr::group_by(corrosion, lineage) %>%
  tidyr::unite(col = "id", corrosion, lineage, sep = "_", remove = FALSE)
Reduce(intersect, split(abund_long_owr$species, abund_long_owr$id))



# Total macrofaunal abundance ####
abund_long_owr <- abund_long[abund_long$site == "Old Woman's River",]
View(abund_long_owr)

abund_tot_owr <- abund_long_owr %>%
  na.omit(abund_long_owr) %>%
  dplyr::group_by(corrosion, lineage, quadrat) %>%
  summarize(tot_abund = sum(abundance))
head(abund_tot_owr)

abund_tot_owr1 <- abund_tot_owr  %>%
  dplyr::group_by(corrosion, lineage) %>%
  summarise(
    mean_abund = mean(tot_abund),
    sd_abund = sd(tot_abund)
  )
View(abund_tot_owr1)



# Total macrofaunal biomass ####
biom_long_owr <- biom_long[biom_long$site == "Old Woman's River",]
View(biom_long_owr)

biom_tot_owr <- biom_long_owr %>%
  na.omit(biom_long) %>%
  dplyr::group_by(corrosion, lineage, quadrat) %>%
  summarize(tot_biom = sum(biomass))
head(biom_tot_owr)

biom_tot_owr1 <- biom_tot_owr  %>%
  dplyr::group_by(corrosion, lineage) %>%
  summarise(
    mean_biom = mean(tot_biom),
    sd_biom = sd(tot_biom)
  )
View(biom_tot_owr1)



# Species richness (S) ###
richness_owr <- abund_long_owr %>%
  dplyr::group_by(corrosion, lineage, quadrat) %>%
  summarise(count = n_distinct(species))
View(richness_owr)

richness_owr1 <- richness_owr %>%
  dplyr::group_by(lineage, corrosion) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
richness_owr1



# Species diversity w/ Shannon diversity (H') ####
abund_long_owr1 <- abund_long_owr[,-7] # Transform the 'abund_long1' dataset into a wide format
abund_wide_owr <- spread(abund_long_owr1, species, abundance)
abund_wide_owr[is.na(abund_wide_owr)] <- 0
View(abund_wide_owr)

abund_wide_owr1 <- abund_wide_owr # Create a new wide-format dataset to add the community indexes

shannon_owr <- diversity(abund_wide_owr[,6:75], "shannon") # Calculate Shannon's diversity index (H')
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(shannon = shannon_owr, .after = "quadrat")

shannon_owr1 <- abund_wide_owr1 %>%
  dplyr::group_by(lineage, corrosion) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon)
  )
View(shannon_owr1)



# Species diversity w/ Simpson diversity (λ) ####
simpson_owr <- diversity(abund_wide_owr[,6:75], "simpson")
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(simpson = simpson_owr, .after = "shannon")

simpson_owr1 <- abund_wide_owr1 %>%
  dplyr::group_by(lineage, corrosion) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson)
  )
View(simpson_owr1)



# Pielou's evenness (J) ####
View(abund_wide_owr1) ; View(richness_owr)

abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(richness = richness_owr$count, .after = "simpson")
pielou_owr <- abund_wide_owr1$shannon/log(abund_wide_owr1$richness) # Formula : J = H'/ln(S)
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(pielou = pielou_owr, .after = "richness")

pielou_owr1 <- abund_wide_owr1 %>%
  dplyr::group_by(lineage, corrosion) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou)
  )
View(pielou_owr1)





# Visualize data sets ====

# GRAPH : Total macrofaunal abundance ####

# Barplot of the total macrofaunal abundance per quadrat, corrosion level and lineage at OWR
barplot_abund_tot_owr <- ggplot(abund_tot_owr, aes(x = quadrat, y = tot_abund, 
                                                   color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Total abundance")
barplot_abund_tot_owr
# Quadrats that could introduce biases : OWR-E-I-4 (very high total abundance)

# Barplot of the total macrofaunal abundance per corrosion level and lineage at OWR
barplot_abund_tot_owr1 <- ggplot(abund_tot_owr1, aes(x = corrosion, y = mean_abund, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average total abundance")
barplot_abund_tot_owr1



# GRAPH : Total macrofaunal biomass ####

# Barplot of the total macrofaunal biomass per quadrat, corrosion level and lineage at OWR
barplot_biom_tot_owr <- ggplot(biom_tot_owr, aes(x = quadrat, y = tot_biom, 
                                                 color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Total biomass")
barplot_biom_tot_owr

# Barplot of the total macrofaunal biomass per corrosion level and lineage at OWR
barplot_biom_tot_owr1 <- ggplot(biom_tot_owr1, aes(x = corrosion, y = mean_biom, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_biom - sd_biom, ymax = mean_biom + sd_biom), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average total biomass")
barplot_biom_tot_owr1



# GRAPH : Species richness (S) ####

# Barplot of the species richness per quadrat, corrosion level and lineage at OWR
barplot_rich_owr <- ggplot(richness_owr, aes(x = quadrat, y = count, 
                                             color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Species richness")
barplot_rich_owr

# Barplot of the species richness per corrosion level and lineage at Old Woman's River
barplot_rich_owr1 <- ggplot(richness_owr1, aes(x = corrosion, y = count_mean, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = count_mean - count_sd, ymax = count_mean + count_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average species richness")
barplot_rich_owr1



# GRAPH : Shannon-Wiener diversity index (H') ####

# Barplot of the species diversity for each quadrat, lineage and corrosion level at OWR
barplot_shannon_owr <- ggplot(abund_wide_owr1, aes(x = quadrat, y = shannon, 
                                                   color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Shannon-Wiener diversity (H')")
barplot_shannon_owr

# Barplot of the species diversity per corrosion and lineage at OWR
barplot_shannon_owr1 <- ggplot(shannon_owr1, aes(x = corrosion, y = shannon_mean, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average Shannon-Wiener diversity (H')")
barplot_shannon_owr1



# GRAPH : Simpson's diversity index (λ) ####

# Barplot of the species diversity per quadrat, corrosion level and lineage at OWR
barplot_simpson_owr <- ggplot(abund_wide_owr1, aes(x = quadrat, y = simpson, 
                                                   color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Simpson's diversity (λ) ")
barplot_simpson_owr

# Barplot of the species diversity per corrosion level and lineage at OWR
barplot_simpson_owr1 <- ggplot(simpson_owr1, aes(x = corrosion, y = simpson_mean, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = simpson_mean - simpson_sd, ymax = simpson_mean + simpson_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average Simpson's diversity (λ)")
barplot_simpson_owr1



# GRAPH : Pielou's evenness (J) ####

# Barplot of the species evenness per quadrat, corrosion level and lineage at OWR
barplot_pielou_owr <- ggplot(abund_wide_owr1, aes(x = quadrat, y = pielou, 
                                                  color = corrosion, fill = corrosion)) + 
  geom_bar(stat = "identity", position = position_dodge()) +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, nrow = 1) +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Quadrat", y = "Pielou's evenness (J)")
barplot_pielou_owr

# Barplot of the species evenness per infestation level and lineage at OWR
barplot_pielou_owr1 <- ggplot(pielou_owr1, aes(x = corrosion, y = pielou_mean, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "none", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average Pielou's evenness (J)")
barplot_pielou_owr1





##
# STATS : Total macrofaunal abundance ====
##

# Assumptions ####

# Frequency table
table(abund_tot_owr$lineage, abund_tot_owr$corrosion)
# Data has equal sample sizes within levels of independent groupings.

# Visualize data
hist(abund_tot_owr$tot_abund)

boxplot_abundtot_owr <- ggboxplot(abund_tot_owr, x = "lineage", y = "tot_abund", color = "corrosion",
                                  palette = col_corr)
boxplot_abundtot_owr

interaction_abundtot_owr <- ggline(abund_tot_owr, x = "lineage", y = "tot_abund", color = "corrosion",
                                   add = c("mean_se", "dotplot"),
                                   palette = col_corr)
interaction_abundtot_owr

# Assumption of normality
shapiro.test(abund_tot_owr$tot_abund) # Total abundance is NOT linear (p = 0.01644).
ggdensity(abund_tot_owr$tot_abund) # Data distribution looks right skewed.
ggqqplot(abund_tot_owr$tot_abund) # Points fall nicely along the line.



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_tot_abund_owr <- aov(tot_abund ~ corrosion * lineage, data = abund_tot_owr)
Anova(anova_tot_abund_owr, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_tot_abund_owr, 1) # Points 1, 4 and 6 detected as outliers.
leveneTest(tot_abund ~ corrosion * lineage, data = abund_tot_owr) # Variances are homogeneous (p = 0.172).

plot(anova_tot_abund_owr, 2) # Points 1, 4 and 6 detected as outliers
anova_tot_abund_owr_residuals <- residuals(object = anova_tot_abund_owr)
shapiro.test(x = anova_tot_abund_owr_residuals) # Residuals are normally distributed (p = 0.583).

# Type II ANOVA without interaction --> Table S10a
anova_tot_abund_owr1 <- aov(tot_abund ~ corrosion + lineage, data = abund_tot_owr)
Anova(anova_tot_abund_owr1, type = "II") # No significant effects of corrosion levels (p = 0.6105) or lineages (p = 0.4668).

plot(anova_tot_abund_owr1, 1) # Points 1, 4 and 12 detected as outliers.
leveneTest(tot_abund ~ corrosion * lineage, data = abund_tot_owr)

plot(anova_tot_abund_owr1, 2) # Points 1, 4 and 12 detected as outliers.
anova_tot_abund_owr1_residuals <- residuals(object = anova_tot_abund_owr1)
shapiro.test(x = anova_tot_abund_owr1_residuals) # Residuals are normally distributed (p = 0.2802).





##
# STATS : Total macrofaunal biomass ====
##

# Assumptions ####

# Visualize data
hist(biom_tot_owr$tot_biom)
hist(log(biom_tot_owr$tot_biom))

boxplot_biomtot_owr <- ggboxplot(biom_tot_owr, x = "lineage", y = "tot_biom", color = "corrosion",
                                 palette = col_corr)
boxplot_biomtot_owr

interaction_biomtot_owr <- ggline(biom_tot_owr, x = "lineage", y = "tot_biom", color = "corrosion",
                                  add = c("mean_se", "dotplot"),
                                  palette = col_corr)
interaction_biomtot_owr



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_tot_biom_owr <- aov(tot_biom ~ corrosion * lineage, data = biom_tot_owr)
Anova(anova_tot_biom_owr, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_tot_biom_owr, 1) # Points 1, 4 and 14 detected as outliers.
leveneTest(tot_biom ~ corrosion * lineage, data = biom_tot_owr) # (!) Variances are NOT homogeneous (p = 0.037).

plot(anova_tot_biom_owr, 2) # Points 1, 4 and 14 detected as outliers
anova_tot_biom_owr_residuals <- residuals(object = anova_tot_biom_owr)
shapiro.test(x = anova_tot_biom_owr_residuals) # Residuals are normally distributed (p = 0.08594).

# Type II ANOVA w/out interaction --> NOT SELECTED
anova_tot_biom_owr1 <- aov(tot_biom ~ corrosion + lineage, data = biom_tot_owr)
Anova(anova_tot_biom_owr1, type = "II") # No significant effect of corrosion levels (p = 0.73) or lineages (p = 0.6176).

plot(anova_tot_biom_owr1, 1) # Points 1, 4 and 12 detected as outliers.
leveneTest(tot_biom ~ corrosion * lineage, data = biom_tot_owr) # Variances are NOT homogeneous (p = 0.037)

plot(anova_tot_biom_owr1, 2) # Points 1, 4 and 12 detected as outliers.
anova_tot_biom_owr1_residuals <- residuals(object = anova_tot_biom_owr1)
shapiro.test(x = anova_tot_biom_owr1_residuals) # Residuals are normally distributed (p = 0.09078).

# Type II ANOVA w/out interaction on log-transformed data --> Table S10b
anova_tot_biom_owr2 <- aov(log(tot_biom) ~ corrosion + lineage, data = biom_tot_owr)
Anova(anova_tot_biom_owr2, type = "II") # No significant effect of infestation levels (p = 0.9591) or lineages (p = 0.5835).

plot(anova_tot_biom_owr2, 1) # Points 1, 4 and 15 detected as outliers.
leveneTest(log(tot_biom) ~ corrosion * lineage, data = biom_tot_owr) # Variances are homogeneous (p = 0.07973).

plot(anova_tot_biom_owr2, 2) # Points 1, 4 and 15 detected as outliers.
anova_tot_biom_owr2_residuals <- residuals(object = anova_tot_biom_owr2)
shapiro.test(x = anova_tot_biom_owr2_residuals) # Residuals are normally distributed (p = 0.1482).

AIC(anova_tot_biom_owr1) ; AIC(anova_tot_biom_owr2)





##
# STATS : Species richness ====
##

# Assumptions ####

# Frequency table
table(abund_tot_owr$lineage, abund_tot_owr$corrosion)
# Data has equal sample sizes within levels of independent groupings.

# Visualize data
hist(richness_owr$count)

boxplot_rich_owr <- ggboxplot(richness_owr, x = "lineage", y = "count", color = "corrosion",
                              palette = col_corr)
boxplot_rich_owr

interaction_rich_owr <- ggline(richness_owr, x = "lineage", y = "count", color = "corrosion",
                               add = c("mean_se", "dotplot"),
                               palette = col_corr)
interaction_rich_owr



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_rich_owr <- aov(count ~ corrosion * lineage, data = richness_owr)
Anova(anova_rich_owr, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_rich_owr, 1) # Points 5, 6 and 9 detected as outliers.
leveneTest(count ~ corrosion * lineage, data = richness_owr) # Variances are homogeneous (p = 0.5042).

plot(anova_rich_owr, 2) # Points 5, 6 and 9 detected as outliers
anova_rich_owr_residuals <- residuals(object = anova_rich_owr)
shapiro.test(x = anova_rich_owr_residuals) # Residuals are normally distributed (p = 0.6937).



# Type II ANOVA w/out interaction --> Table S10c
anova_rich_owr1 <- aov(count ~ corrosion + lineage, data = richness_owr)
Anova(anova_rich_owr1, type = "II") # No significant effect of corrosion levels (p = 0.1263) or lineages (p = 0.6447).

plot(anova_rich_owr1, 1) # Points 5, 6 and 9 detected as outliers.
leveneTest(count ~ corrosion * lineage, data = richness_owr)

plot(anova_rich_owr1, 2) # Points 5, 6 and 9 detected as outliers
anova_rich_owr1_residuals <- residuals(object = anova_rich_owr1)
shapiro.test(x = anova_rich_owr1_residuals) # Residuals are normally distributed (p = 0.4113).





##
# STATS : Species diversity w/ Shannon diversity (H') ====
##

# Assumptions ####

# Frequency table
table(abund_tot_owr$lineage, abund_tot_owr$corrosion)
# Data has equal sample sizes within levels of independent groupings.

# Visualize data
hist(abund_wide_owr1$shannon)

boxplot_shannon_owr <- ggboxplot(abund_wide_owr1, x = "lineage", y = "shannon", color = "corrosion",
                                 palette = col_corr)
boxplot_shannon_owr

interaction_shannon_owr <- ggline(abund_wide_owr1, x = "lineage", y = "shannon", color = "corrosion",
                                  add = c("mean_se", "dotplot"),
                                  palette = col_corr)
interaction_shannon_owr



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_shannon_owr <- aov(shannon ~ corrosion * lineage, data = abund_wide_owr1)
Anova(anova_shannon_owr, type = "III") # Interaction is not significant, so type II ANOVA is more powerful. 

plot(anova_shannon_owr, 1) # Points 2, 5 and 6 detected as outliers.
leveneTest(shannon ~ corrosion * lineage, data = abund_wide_owr1) # Variances are homogeneous (p = 0.8312).

plot(anova_shannon_owr, 2) # Points 2, 5 and 6 detected as outliers.
anova_shannon_owr_residuals <- residuals(object = anova_shannon_owr)
shapiro.test(x = anova_shannon_owr_residuals) # Residuals are normally distributed (p = 0.9793).



# Type II ANOVA w/out interaction --> Table S10d
anova_shannon_owr1 <- aov(shannon ~ corrosion + lineage, data = abund_wide_owr1)
Anova(anova_shannon_owr1, type = "II") # No significant effect of corrosion levels (p = 0.1208) or lineages (p = 0.3016).

plot(anova_shannon_owr1, 1) # Points 2, 4 and 6 detected as outliers.
leveneTest(shannon ~ corrosion * lineage, data = abund_wide_owr1)

plot(anova_shannon_owr1, 2) # Points 2, 4 and 6 detected as outliers
anova_shannon_owr1_residuals <- residuals(object = anova_shannon_owr1)
shapiro.test(x = anova_shannon_owr1_residuals) # Residuals are normally distributed (p = 0.5612).





##
# STATS : Species diversity w/ Simpson diversity (λ) ====
##

# Assumptions ####

# Frequency table
table(abund_tot_owr$lineage, abund_tot_owr$corrosion)
# Data has equal sample sizes within levels of independent groupings.

# Visualize data
hist(abund_wide_owr1$simpson) # Data distribution is left skewed.
shapiro.test(x = simpson) # Data is not normally distributed (p = 0.01573).

boxplot_simpson_owr <- ggboxplot(abund_wide_owr1, x = "lineage", y = "simpson", color = "corrosion",
                                 palette = col_corr)
boxplot_simpson_owr

interaction_simpson_owr <- ggline(abund_wide_owr1, x = "lineage", y = "simpson", color = "corrosion",
                                  add = c("mean_se", "dotplot"),
                                  palette = col_corr)
interaction_simpson_owr



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_simpson_owr <- aov(simpson ~ corrosion * lineage, data = abund_wide_owr1)
Anova(anova_simpson_owr, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_simpson_owr, 1) # Points 2, 4 and 5 detected as outliers.
leveneTest(simpson ~ corrosion * lineage, data = abund_wide_owr1) # Variances are homogeneous (p = 0.6669).

plot(anova_simpson_owr, 2) # Points 2, 4 and 5 detected as outliers
anova_simpson_owr_residuals <- residuals(object = anova_simpson_owr)
shapiro.test(x = anova_simpson_owr_residuals) # Residuals are normally distributed (p = 0.6943).



# Type II ANOVA w/out interaction --> Table S10e
anova_simpson_owr1 <- aov(simpson ~ corrosion + lineage, data = abund_wide_owr1)
Anova(anova_simpson_owr1, type = "II") # No significant effect of corrosion levels (p = 0.1593) or lineages (p = 0.386).

plot(anova_simpson_owr1, 1) # Points 2, 4 and 5 detected as outliers.
leveneTest(simpson ~ infestation * lineage, data = abund_wide_owr1)

plot(anova_simpson_owr1, 2) # Points 2, 4 and 5 detected as outliers
anova_simpson_owr1_residuals <- residuals(object = anova_simpson_owr1)
shapiro.test(x = anova_simpson_owr1_residuals) # Residuals are normally distributed (p = 0.9246).





##
# STATS : Pielou's evenness (J) ====
##

# Assumptions ####

# Frequency table
table(abund_tot_owr$lineage, abund_tot_owr$corrosion)
# Data has equal sample sizes within levels of independent groupings

# Visualize data
hist(abund_wide_owr1$pielou)

boxplot_pielou_owr <- ggboxplot(abund_wide_owr1, x = "lineage", y = "pielou", color = "corrosion",
                                palette = col_corr)
boxplot_pielou_owr

interaction_pielou_owr <- ggline(abund_wide_owr1, x = "lineage", y = "pielou", color = "corrosion",
                                 add = c("mean_se", "dotplot"),
                                 palette = col_corr)
interaction_pielou_owr



# Classic ANOVAs ####

# Type III ANOVA w/ interaction --> NOT SELECTED
anova_pielou_owr <- aov(pielou ~ corrosion * lineage, data = abund_wide_owr1)
Anova(anova_pielou_owr, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_pielou_owr, 1) # Points 2, 4 and 5 detected as outliers.
leveneTest(pielou ~ corrosion * lineage, data = abund_wide_owr1) # Variances are homogeneous (p = 0.6906).

plot(anova_pielou_owr, 2) # Points 2, 4 and 5 detected as outliers
anova_pielou_owr_residuals <- residuals(object = anova_pielou_owr)
shapiro.test(x = anova_pielou_owr_residuals) # Residuals are normally distributed (p = 0.8791)



# Type II ANOVA w/out interaction --> Table S10f
anova_pielou_owr1 <- aov(pielou ~ corrosion + lineage, data = abund_wide_owr1)
Anova(anova_pielou_owr1, type = "II") # No significant effect of infestation levels (p = 0.1852) or lineages (p = 0.2498).

plot(anova_pielou_owr1, 1) # Points 2, 4 and 12 detected as outliers.
leveneTest(pielou ~ corrosion * lineage, data = abund_wide_owr1)

plot(anova_pielou_owr1, 2) # Points 2, 4 and 12 detected as outliers
anova_pielou_owr1_residuals <- residuals(object = anova_pielou_owr1)
shapiro.test(x = anova_pielou_owr1_residuals) # Residuals are normally distributed (p = 0.9538).




