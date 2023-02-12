## ---------------------------
##
## Script name: Architectural complexity associated with manipulated mussel beds depending on their level
##              of euendolithic infestation and their transplant site
##
## Purpose of script: 
##    - Is there a difference between the total number of mussels between each quadrat,
##      each euendolithic infestation treatment and each site ?
##    - Is there a difference between the number of broken shells between each quadrat,
##      each euendolithic infestation treatment and each site ?
##    - Is there a difference between the average number of byssal threads per mussel between each quadrat,
##      each euendolithic infestation treatment and each site ?
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

# Environmental variables - number of byssal threads ==========================================
byssus <- read.csv("./RAW DATA/infauna_architecture.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(byssus)
byssus <- na.omit(byssus)

byssus[,factor1] <- lapply(byssus[,factor1], factor)

View(byssus)





###############################################################################################
# Section: Architectural complexity -----------------------------------------------------------
###############################################################################################

###############################################################################################
# Architectural complexity ACROSS ALL SITES ===================================================
###############################################################################################

# Curate dataset ##############################################################################

# Select only the relevant variables and exclude Old Woman's River and NA
env1 <- env[, c(1:4, 9:10)]
env1 <- env1[env1$site !="Old Woman's River",]
View(env1)

byssus1 <- byssus[byssus$site !="Old Woman's River",]
View(byssus1)

# Create a unique identifier for each quadrat
env1 <- env1 %>%
  dplyr::group_by(site, infestation, lineage, quadrat) %>%
  tidyr::unite(col = "id", site, infestation, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

byssus1 <- byssus1 %>%
  dplyr::group_by(site, infestation, lineage, quadrat) %>%
  tidyr::unite(col = "id", site, infestation, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

# Transform the dataset from wide to long format
env2 <- env1 %>% pivot_longer(cols = c('nb_broken', 'tot_mussels'),
                    names_to = 'mussels',
                    values_to = 'nb')
View(env2)





# Barplot of the total number of mussels and number of broken shells for each quadrat #########
barplot_mussels <- ggplot(env2, aes(x = quadrat, y = nb, color = mussels, fill = mussels)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ site * infestation * lineage, ncol = 4) +
  theme(legend.position = "bottom")
barplot_mussels

# Boxplot of the number of byssal threads per mussel for each quadrat #########################
boxplot_bt <- ggplot(byssus1, aes(x = quadrat, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, ncol = 4) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Quadrat", y = "Average number of byssal threads")
boxplot_bt

boxplot_bt1 <- ggplot(byssus1, aes(x = site, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 375,
           label = c("B", "A", "AC", "B")) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average number of byssal threads")
boxplot_bt1





# Statistical analyses on the total number of mussels #########################################

# ANOVA w/ unbalanced design (i.e., not equal sample sizes within levels of independent groupings)
View(env1)
table(env1$site, env1$infestation) # Frequency table

# Visualize data
boxplot_totmussels <- ggboxplot(env1, x = "site", y = "tot_mussels", color = "infestation",
                                palette = my_colors)
boxplot_totmussels

interaction_totmusssels <- ggline(env1, x = "site", y = "tot_mussels", color = "infestation",
                                  add = c("mean_se", "dotplot"),
                                  palette = my_colors)
interaction_totmusssels

# Type III ANOVA w/ interaction
anova_totmussels <- aov(tot_mussels ~ infestation * site, data = env1)
Anova(anova_totmussels, type = "III") # Interaction is not significant, so type II is more powerful

plot(anova_totmussels, 1) # Points 2, 3 and 4 detected as outliers. No difference in results when excluded from ANOVA.
leveneTest(tot_mussels ~ infestation * site, data = env1) # Variances are homogeneous (but better when excluding outlier quadrats)

plot(anova_totmussels, 2) # Points fall approximately along the line
anova_totmussels_residuals <- residuals(object = anova_totmussels) 
shapiro.test(x = anova_totmussels_residuals) # Residuals are normally distributed

# Type II ANOVA without interaction --> to use because more powerful (Table S1)
anova_totmussels1 <- aov(tot_mussels ~ infestation + site, data = env1)
Anova(anova_totmussels1, type = "II") # Significant effect of site but not of euendolithic infestation.

plot(anova_totmussels1, 1)
plot(anova_totmussels1, 2)
anova_totmussels1_residuals <- residuals(object = anova_totmussels1)
shapiro.test(x = anova_totmussels1_residuals) # Residuals are normally distributed

# Re-run the analyses without the following quadrats: PEd-INF-3, PEd-NINF-2 and PEd-NINF-4
# env1_no <- env1[-c(2, 4, 7), ]
# View(env1_no)

# boxplot_totmussels_no <- ggboxplot(env1_no, x = "site", y = "tot_mussels", color = "infestation",
#                                palette = my_colors)
# boxplot_totmussels_no

# interaction_totmusssels_no <- ggline(env1_no, x = "site", y = "tot_mussels", color = "infestation",
#                                  add = c("mean_se", "dotplot"),
#                                  palette = my_colors)
# interaction_totmusssels_no

# anova_totmussels_no <- aov(tot_mussels ~ infestation * site, data = env1_no)
# Anova(anova_totmussels_no, type = "III")
# Same results than ANOVA including the outlier quadrats.

# Pairwise comparisons between sites --> to use (Table S2)
summary(glht(anova_totmussels1, linfct = mcp(site = "Tukey")))





# Statistical analyses on the number of broken shells #########################################

# Visualize data
boxplot_broken <- ggboxplot(env1, x = "site", y = "nb_broken", color = "infestation",
                                palette = my_colors)
boxplot_broken

interaction_broken <- ggline(env1, x = "site", y = "nb_broken", color = "infestation",
                                  add = c("mean_se", "dotplot"),
                                  palette = my_colors)
interaction_broken

# Type III ANOVA w/ interaction
anova_broken <- aov(nb_broken ~ infestation * site, data = env1)
Anova(anova_broken, type = "III") # Interaction is not significant, so type II is more powerful

# Type II ANOVA without interaction --> to use because more powerful (Table S3)
anova_broken1 <- aov(nb_broken ~ infestation + site, data = env1)
Anova(anova_broken1, type = "II") # Significant effect of site but not of euendolithic infestation.

plot(anova_broken1, 1) # Points 6, 22 and 27 detected in outliers. 
leveneTest(nb_broken ~ infestation * site, data = env1) # Variances are homogeneous

plot(anova_broken1, 2)
anova_broken1_residuals <- residuals(object = anova_broken1)
shapiro.test(x = anova_broken1_residuals) # Residuals are NOT normally distributed

# Fortunately, an ANOVA is not very sensitive to moderate deviations from normality; simulation studies, 
# using a variety of non-normal distributions, have shown that the false positive rate is not affected 
# very much by this violation of the assumption (Glass et al. 1972, Harwell et al. 1992, Lix et al. 1996).

# Non-parametric equivalent to two-way ANOVA

# Align-and-rank data for nonparametric factorial ANOVA
# install.packages("ARTool")
# library(ARTool)

# anova_broken2 <- art(nb_broken ~ infestation * site, data = env1)
# anova(anova_broken2)
# Same results with non-parametric ANOVA.

# Pairwise comparisons between sites --> to use
summary(glht(anova_broken1, linfct = mcp(site = "Tukey")))





# Statistical analyses on the number of byssal threads per mussel #########################################

# Type III ANOVA (using aov) with unbalanced design
byssus_aov <- aov(nb_byssal ~ site * infestation, data = byssus1)
Anova(byssus_aov, type = "III") # No significant interaction

# Type II ANOVA w/ unbalanced design --> to use because more powerful (Table S4)
byssus_aov1 <- aov(nb_byssal ~ infestation + site, data = byssus1) 
Anova(byssus_aov1, type = "II") # Effect of site on the number of byssal threads

# Check ANOVA assumptions
plot(byssus_aov1, 1) # Points 19, 73 and 77 detected as outliers
leveneTest(nb_byssal ~ site * infestation, data = byssus1) # Variances are homogeneous

plot(byssus_aov1, 2) # Points fall more or less on the line
byssus_aov1_residuals <- residuals(object = byssus_aov1)
shapiro.test(x = byssus_aov1_residuals) # Residuals are normally distributed

# Pairwise comparisons between sites --> to use (Table S5)
summary(glht(byssus_aov1, linfct = mcp(site = "Tukey")))





###############################################################################################
# Architectural complexity ACROSS AT OWR ======================================================
###############################################################################################

# Curate dataset ##############################################################################

# Select only the relevant variables and exclude Old Woman's River and NA
env_owr <- env[, c(1:4, 9:10)]
env_owr <- env_owr[env_owr$site =="Old Woman's River",]
View(env_owr)

byssus_owr <- byssus[byssus$site =="Old Woman's River",]
View(byssus_owr)

# Create a unique identifier for each quadrat
env_owr <- env_owr %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  tidyr::unite(col = "id", infestation, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

byssus_owr <- byssus_owr %>%
  dplyr::group_by(infestation, lineage, quadrat) %>%
  tidyr::unite(col = "id", infestation, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

# Transform the dataset from wide to long format
env_owr1 <- env_owr %>% pivot_longer(cols = c('nb_broken', 'tot_mussels'),
                              names_to = 'mussels',
                              values_to = 'nb')
View(env_owr1)





# Barplot of the total number of mussels and number of broken shells for each quadrat #########
barplot_mussels_owr <- ggplot(env_owr1, aes(x = quadrat, y = nb, color = mussels, fill = mussels)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ infestation * lineage, ncol = 4) +
  theme(legend.position = "bottom")
barplot_mussels_owr

# Boxplot of the number of byssal threads per mussel for each quadrat #########################
boxplot_bt_owr <- ggplot(byssus_owr, aes(x = quadrat, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, ncol = 4) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Quadrat", y = "Average number of byssal threads")
boxplot_bt_owr

boxplot_bt_owr1 <- ggplot(byssus_owr, aes(x = infestation, y = nb_byssal, fill = infestation)) + 
  geom_boxplot() +
  facet_wrap(~ lineage, ncol = 2) +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = my_colors) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Infestation",
       x = "Site", y = "Average number of byssal threads")
boxplot_bt_owr1





# Statistical analyses on the total number of mussels #########################################

# ANOVA w/ balanced design (i.e., equal sample sizes within levels of independent groupings)
View(env_owr)
table(env_owr$lineage, env_owr$infestation) # Frequency table

# Visualize data
boxplot_totmussels_owr <- ggboxplot(env_owr, x = "lineage", y = "tot_mussels", color = "infestation",
                                palette = my_colors)
boxplot_totmussels_owr

interaction_totmusssels_owr <- ggline(env_owr, x = "lineage", y = "tot_mussels", color = "infestation",
                                  add = c("mean_se", "dotplot"),
                                  palette = my_colors)
interaction_totmusssels_owr

# ANOVA w/ interaction
anova_totmussels_owr <- aov(tot_mussels ~ infestation * lineage, data = env_owr)
Anova(anova_totmussels_owr) # Interaction is not significant

# ANOVA w/out interaction
anova_totmussels_owr1 <- aov(tot_mussels ~ infestation + lineage, data = env_owr)
Anova(anova_totmussels_owr1) 

plot(anova_totmussels_owr1, 1)
plot(anova_totmussels_owr1, 2)
anova_totmussels_owr1_residuals <- residuals(object = anova_totmussels_owr1)
shapiro.test(x = anova_totmussels_owr1_residuals) # Residuals are normally distributed

# Pairwise comparisons between sites --> to use (Table S2)
summary(glht(anova_totmussels_owr1, linfct = mcp(lineage = "Tukey")))





# Statistical analyses on the number of broken shells #########################################

# Visualize data
boxplot_broken_owr <- ggboxplot(env_owr, x = "lineage", y = "nb_broken", color = "infestation",
                            palette = my_colors)
boxplot_broken_owr

interaction_broken_owr <- ggline(env_owr, x = "lineage", y = "nb_broken", color = "infestation",
                             add = c("mean_se", "dotplot"),
                             palette = my_colors)
interaction_broken_owr

# Type I ANOVA w/ interaction
anova_broken_owr <- aov(nb_broken ~ infestation * lineage, data = env_owr)
Anova(anova_broken_owr) # Interaction is not significant, so type II is more powerful

# Type I ANOVA without interaction 
anova_broken_owr1 <- aov(nb_broken ~ infestation + lineage, data = env_owr)
Anova(anova_broken_owr1)

plot(anova_broken_owr1, 1) # Points 2, 10 and 16 detected in outliers. 
leveneTest(nb_broken ~ infestation * site, data = env_owr) # Variances are homogeneous

plot(anova_broken_owr1, 2)
anova_broken_owr1_residuals <- residuals(object = anova_broken_owr1)
shapiro.test(x = anova_broken_owr1_residuals) # Residuals are normally distributed






# Statistical analyses on the number of byssal threads per mussel #########################################

# Type III ANOVA (using aov) with balanced design
byssus_aov_owr <- aov(nb_byssal ~ infestation * lineage, data = byssus_owr)
Anova(byssus_aov_owr) # No significant interaction

byssus_aov_owr1 <- aov(nb_byssal ~ infestation + lineage, data = byssus_owr) 
Anova(byssus_aov_owr1) 

# Check ANOVA assumptions
plot(byssus_aov_owr1, 1) # Points 3, 11 and 42 detected as outliers
leveneTest(nb_byssal ~ lineage * infestation, data = byssus_owr) # Variances are homogeneous

plot(byssus_aov_owr1, 2) # Points fall more or less on the line
byssus_aov_owr1_residuals <- residuals(object = byssus_aov_owr1)
shapiro.test(x = byssus_aov_owr1_residuals) # Residuals are NOT normally distributed


