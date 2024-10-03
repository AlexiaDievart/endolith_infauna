## ---
##
## Script name: Architectural within-bed complexity associated with manipulated mussel beds depending on their level
##              of euendolithic corrosion
##
## Purpose of script: 
##    - EUENDOLITHS x BIOGEOGRAPHY : Compare the following variables between corrosion status and among transplant sites
##      (excluding Old Woman's River) 
##        - Total number of mussels (live + dead, Perna + Mytilus)
##        - Number of live mussels (Perna + Mytilus)
##        - Number of dead mussels (Perna + Mytilus)
##        - Number of broken mussels (Perna + Mytilus)
##        - Average number of byssal threads/mussel
##    - EUENDOLITHS x PERNA LINEAGES : Compare the above variables between corrosion status and Perna lineages at Old Woman's River
##
## Find in the manuscript :
##    - EUENDOLITHS x BIOGEOGRAPHY
##        - Figure 2a : Boxplot of the average number of byssal threads per mussel between corroded and non-corroded mussel beds
##                      for each site (except Old Woman's River)
##        - Table S3 : Results of the series of two-way ANOVAs on architectural complexity variables
##            - Table S3a : Total number of mussels --> Significant effect of site
##            - Table S3b : Number of live mussels --> Significant effect of site
##            - Table S3c : Number of dead mussels
##            - Table S3d : Number of shell fragments
##            - Table S3e : Average number of byssal threads per mussel --> Significant effect of site
##        - Table S4 : Pairwise comparisons (when significant effect detected) on architectural complexity variables
##            - Table S4a : Total number of mussels --> PEd < MB
##            - Table S4b : Number of live mussels --> PEd < MB
##            - Table S4c : Average number of byssal threads per mussel --> JB > BR = PEd and MB > BR
##    - EUENDOLITHS x PERNA LINEAGES
##        - Figure 2e : Boxplot of the average number of byssal threads per mussel between corroded and non-corroded mussel beds,
##                      and eastern vs westerm mussel beds
##        - Table S3 : Results of the series of two-way ANOVAs on architectural complexity variables
##            - Table S3f : Total number of mussels
##            - Table S3g : Number of live mussels
##            - Table S3h : Number of dead mussels --> Significant effect of lineage
##            - Table S3i : Number of shell fragments
##            - Table S3j : Average number of byssal threads per mussel --> Significant effect of lineage
##        - Table S4 : Pairwise comparisons (when significant effect detected) on architectural complexity variables
##            - Table S4d : Number of dead mussels --> W > E
##            - Table S4e : Average number of byssal threads per mussel --> W > E
##
## Author: Alexia Dievart
##
## Date Created: 2023-02-05
## Dates Updated: 2024-09-06 ; 2024-09-23 ; 2024-09-30 ; 2024-10-03
##
## Copyright (c) Alexia DIEVART 2023
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
# Section: Session setup ----
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
# Section: Load data ----
##

# Environmental variables - total number of mussels and number of broken shells ====
arch <- read.csv("./RAW DATA/infauna_habitat.csv", dec = ",", header = T, sep = ";")
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





##
# SECTION: Architectural complexity ACROSS ALL SITES ----
##

# GOAL : Is there a significant difference between corroded and non-corroded mussel beds across all sites
# (excluding Old Woman's River) in terms of total number of mussels, number of live mussels, dead mussels 
# and broken shells per quadrat, and average number of byssal threads / mussel ?

# Curate data set ====

# Exclude Old Woman's River
arch1 <- arch[arch$site != "Old Woman's River",]
View(arch1)

byssus1 <- byssus[byssus$site != "Old Woman's River",]
View(byssus1)

# Calculate habitat architectural complexity parameters (i.e., total number of mussels, number of live and dead mussels and broken shells,
# average number of byssal threads/mussel) for each quadrat
arch2 <- arch1 %>%
  dplyr::group_by(site, corrosion, lineage, quadrat) %>%
  dplyr::summarize(
    total = tot_mussels, 
    live = live_perna + live_mytilus,
    dead = dead_perna + dead_mytilus,
    broken = nb_broken,
    mean_byssal = mean_byssal
    )
View(arch2)

# Create a unique identifier for each quadrat in each data set
arch2 <- arch2 %>%
  dplyr::group_by(site, corrosion, lineage, quadrat) %>%
  tidyr::unite(col = "id", site, corrosion, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

byssus1 <- byssus1 %>%
  dplyr::group_by(site, corrosion, lineage, quadrat) %>%
  tidyr::unite(col = "id", site, corrosion, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

# Transform the data set from wide to long format for graphical purposes
arch3 <- arch2[,1:9] %>% pivot_longer(cols = c('live', 'dead', 'broken'),
                              names_to = 'mussels',
                              values_to = 'nb')
View(arch3)





# Visualize data set ====

# Barplot of the number of live, dead and broken mussels in each quadrat
barplot_mussels <- ggplot(arch3, aes(x = quadrat, y = nb, color = mussels, fill = mussels)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ site * corrosion * lineage, ncol = 4) +
  theme(legend.position = "bottom")
barplot_mussels

# Boxplot of the number of live, dead and broken mussels for each corrosion status and site
boxplot_mussels <- ggplot(arch3, aes(x = site, y = nb, fill = corrosion)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ mussels) +
  scale_fill_manual(values = col_corr) +
  scale_x_discrete(labels = c("BR", "JB", "MB", "PEd")) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Site", y = "Mussel bed structure")
boxplot_mussels
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).
# Quadrats that could introduce biases : PEd-E-I-3, PEd-E-NI-2 and PEd-E-NI-4 (loss of mussels due to wave action).

# Boxplot of the average number of byssal threads per mussel for each quadrat
boxplot_bt <- ggplot(byssus1, aes(x = quadrat, y = nb_byssal, fill = corrosion)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ site, ncol = 4) +
  scale_fill_manual(values = col_corr) +
  scale_x_discrete(limits = site_order, labels = sites) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Quadrat", y = "Average number of byssal threads")
boxplot_bt

# Boxplot of the average number of byssal threads per mussel for each corrosion status and site --> Figure 2a
boxplot_bt1 <- ggplot(byssus1, aes(x = site, y = nb_byssal, fill = corrosion)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 430,
           label = c("AC", "BC", "A", "C"),
           size = 6) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Site", y = "Average number of byssal threads",
       title = "Byssal threads")
boxplot_bt1
# Missing quadrats : BR-W-I-3 (not stored properly), JB-W-I-4 and JB-W-NI-4 (lost and not replaced on the field),
#                    PEd-E-I-4 (lost and not replaced on the field).
# Quadrats with no byssal thread counts : PEd-E-I-1 and PEd-E-NI-1 (byssal threads not counted).





##
# STATS : Total number of mussels ACROSS ALL SITES ====
##

# Assumptions ####
View(arch2)
table(arch2$site, arch2$corrosion) # Frequency table
# Because there are not equal sample sizes within levels of independent groupings : ANOVA with unbalanced design.

# Visualize data
boxplot_totmussels <- ggboxplot(arch2, x = "site", y = "total", color = "corrosion",
                                palette = col_corr)
boxplot_totmussels

interaction_totmusssels <- ggline(arch2, x = "site", y = "total", color = "corrosion",
                                  add = c("mean_se", "dotplot"),
                                  palette = col_corr)
interaction_totmusssels

# Assumption of normality

shapiro.test(arch2$total) # (!) Total number of mussels (live + dead) is NOT linear (p = 0.01207).
ggdensity(arch2$total) # Data distribution is right skewed.
ggqqplot(arch2$total) # Points fall approx. along the line.
# Outliers identified but NOT omitted from the present analysis : PEd-NC-E-2 ; PEd-NC-E4 ; PEd-C-E-3.



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_totmussels <- aov(total ~ corrosion * site, data = arch2)
Anova(anova_totmussels, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_totmussels, 1) # 3 points detected as potential outliers (Points 26, 27 and 28).
leveneTest(total ~ corrosion * site, data = arch2) # Variances are homogeneous (p = 0.06699).

plot(anova_totmussels, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_totmussels_residuals <- residuals(object = anova_totmussels) 
shapiro.test(x = anova_totmussels_residuals) # Residuals are normally distributed (p = 0.1702).



# Type II ANOVA without interaction (with outliers) --> Table S3a
anova_totmussels1 <- aov(total ~ corrosion + site, data = arch2)
Anova(anova_totmussels1, type = "II") # There is a significant effect of site (p = 0.01915) BUT NO effect of corrosion (p = 0.86068).

plot(anova_totmussels1, 1) # 3 points detected as potential outliers (Points 26, 27 and 28).

plot(anova_totmussels1, 2) # Points fall approx. along the line.
anova_totmussels1_residuals <- residuals(object = anova_totmussels1)
shapiro.test(x = anova_totmussels1_residuals) # Residuals are normally distributed (p = 0.3269). 



# Pairwise comparisons between sites --> Table S4a
summary(glht(anova_totmussels1, linfct = mcp(site = "Tukey")))
# Results : Total number of mussels in MB > PEd (p = 0.0157)
# Discussion : In PEd, 3 quadrats have a very low number of whole mussels (live and dead), probably due
#              to the mussels slipping through the mesh out of the quadrat due to wave action.





##
# STATS : Nb of live mussels ACROSS ALL SITES ====
##

# Assumptions #####

# Visualize data set
boxplot_live <- ggboxplot(arch2, x = "site", y = "live", color = "corrosion",
                            palette = col_corr)
boxplot_live

interaction_live <- ggline(arch2, x = "site", y = "live", color = "corrosion",
                             add = c("mean_se", "dotplot"),
                             palette = col_corr)
interaction_live


# Assumption of normality

shapiro.test(arch2$live) # (!) Number of live mussels (Mytilus + Perna) is NOT linear (p = 0.04475).
ggdensity(arch2$live) # Data distribution is right skewed.
ggqqplot(arch2$live) # Points fall approx. along the line.
# Outliers identified but NOT omitted from the present analysis : PEd-NC-E-2 ; PEd-NC-E4 ; PEd-C-E-3.



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_livemussels <- aov(live ~ corrosion * site, data = arch2)
Anova(anova_livemussels, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_livemussels, 1) # 3 points detected as potential outliers (Points 26, 27 and 28).
leveneTest(live ~ corrosion * site, data = arch2) # Variances are homogeneous (p = 0.1299).

plot(anova_livemussels, 2) # 3 points detected as potential outliers. Points DO NOT fall along the line.
anova_livemussels_residuals <- residuals(object = anova_livemussels) 
shapiro.test(x = anova_livemussels_residuals) # Residuals are normally distributed (p = 0.08841).



# Type II ANOVA without interaction --> Table S3b
anova_livemussels1 <- aov(live ~ corrosion + site, data = arch2)
Anova(anova_livemussels1, type = "II")  # There is a significant effect of site (p = 0.01465) BUT NO effect of corrosion (p = 0.92638).

plot(anova_livemussels1, 1) # 3 points detected as potential outliers (Points 26, 27 and 28).

plot(anova_livemussels1, 2) # Points fall approx. along the line.
anova_livemussels1_residuals <- residuals(object = anova_livemussels1)
shapiro.test(x = anova_livemussels1_residuals) # Residuals are normally distributed (p = 0.1679). 



# Pairwise comparisons between sites --> Table S4b
summary(glht(anova_livemussels1, linfct = mcp(site = "Tukey")))
# Results : Number of live mussels in MB > PEd (p = 0.0170).
# Discussion : In PEd, 3 quadrats have a very low number of live mussels, probably due
#              to the mussels slipping through the mesh out of the quadrat due to wave action.





##
# STATS : Nb of dead mussels ACROSS ALL SITES ====
##

# Assumptions ####

# Visualize data set
boxplot_dead <- ggboxplot(arch2, x = "site", y = "dead", color = "corrosion",
                          palette = col_corr)
boxplot_dead

interaction_dead <- ggline(arch2, x = "site", y = "dead", color = "corrosion",
                           add = c("mean_se", "dotplot"),
                           palette = col_corr)
interaction_dead

# Assumption of normality

shapiro.test(arch2$dead) # (!) Number of dead mussels (Mytilus + Perna) is NOT linear (p < 0.05).
ggdensity(arch2$dead) # Data distribution is left skewed and count a lot of '0'.
ggqqplot(arch2$dead) # Points DO NOT fall approx. along the line.
# For this analysis, we will perform an ANOVA on a GLM model. 

# Models with outliers ####

# Poisson GLM
model_deadmussels <- glm(dead ~ corrosion * site, 
                         family = "poisson",
                         data = arch2)
summary(model_deadmussels) # AIC = 74.413

# Check for over/underdispersion in the model
E2 <- resid(model_deadmussels, type = "pearson")
N <- nrow(arch2)
p <- length(coef(model_deadmussels))
sum(E2^2) / (N - p) 
# This model produces a bit of overdispersion (result = 1.198718).

# Negative Binomial GLM
model_deadmussels1 <- glm.nb(dead ~ corrosion * site,
                         data = arch2)
summary(model_deadmussels1) # AIC = 72.422

# Check for over/underdispersion in the model
E2 <- resid(model_deadmussels1, type = "pearson")
N <- nrow(arch2)
p <- length(coef(model_deadmussels1)) + 1 # '+1' is for the variance parameter in negative binomial models.
sum(E2^2) / (N - p) 
# This model produces a lot of underdispersion (result = 0.5672704).

# Zero-inflated Poisson GLM
library(pscl)

model_deadmussels2 <- zeroinfl(dead ~ corrosion * site,
                               dist = "poisson",
                               data = arch2)
summary(model_deadmussels2)

model_deadmussels2a <- zeroinfl(dead ~ corrosion + site,
                               dist = "poisson",
                               data = arch2)
summary(model_deadmussels2a)


# Check for over/underdispersion in the model
E2 <- resid(model_deadmussels2, type = "pearson")
N <- nrow(arch2)
p <- length(coef(model_deadmussels2))
sum(E2^2) / (N - p) 
# This is a little bit overdispersed, but really close to 1 (result = 1.074154).

# Zero-Inflated Negative Binomial GLM
model_deadmussels3 <- zeroinfl(dead ~ corrosion * site,
                               dist = "negbin",
                               data = arch2)
summary(model_deadmussels3)
# This last GLM w/out interaction produces NA, which makes sense as there is 0 dead mussels in PEd and MB.

# Check for over/underdispersion in the model
E2 <- resid(model_deadmussels3, type = "pearson")
N <- nrow(arch2)
p <- length(coef(model_deadmussels3))
sum(E2^2) / (N - p) 
# Exact same dispersion as for the Zero-Inflated Poisson GLM (result = 1.074151).

# Compare the two most probable models
library(lmtest)

lrtest(model_deadmussels2a, model_deadmussels3) # There is no difference between the zero-inflated models (p = 0.9033).

AIC(model_deadmussels2a, model_deadmussels3) # Zero-Inflated Poisson model w/out interaction has a lower AIC.

library(stats)
Anova(model_deadmussels2a, type = "II", test = "F") # --> Table S3c

# Results : Number of dead mussels in PEd = MB (0) < JB and BR but besides that, there is no difference between site or corrosion status.

# Discussion : No whole dead mussels were recorded in quadrats from PEd and MB, while quadrats from JB and MB contained a few (< 10),
# which could contribute to the architectural complexity of the mussel bed.
# For example, Desis formidabilis is a spider that use dead mussel shells to protect itself during high tide.
# Question : Is there more Desis formidabilis in quadrats with more whole dead mussels ?





##
# STATS : Nb of broken mussels ACROSS ALL SITES ====
##

# Assumptions ####

# Visualize data
boxplot_broken <- ggboxplot(arch2, x = "site", y = "broken", color = "corrosion",
                          palette = col_corr)
boxplot_broken

interaction_broken <- ggline(arch2, x = "site", y = "broken", color = "corrosion",
                           add = c("mean_se", "dotplot"),
                           palette = col_corr)
interaction_broken

# Assumption of normality

shapiro.test(arch2$broken) #(!) Number of broken mussels (Mytilus + Perna) is NOT linear (p = 0.00168).
shapiro.test(sqrt(arch2$broken)) # Sqrt-transformed number of broken mussels is linear (p = 0.07)

ggdensity(sqrt(arch2$broken)) # Transformed data distribution looks linear.
ggqqplot(sqrt(arch2$broken)) # Transformed points fall approx. along the line.



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_brokenmussels <- aov(sqrt(broken) ~ corrosion * site, data = arch2)
Anova(anova_brokenmussels, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_brokenmussels, 1) # 3 points detected as potential outliers (Points 5, 12 and 23).
leveneTest(sqrt(broken) ~ corrosion * site, data = arch2) # Variances are homogeneous (p = 0.8668).

plot(anova_brokenmussels, 2) # 3 points detected as potential outliers. Points fall approximately along the line.
anova_brokenmussels_residuals <- residuals(object = anova_brokenmussels) 
shapiro.test(x = anova_brokenmussels_residuals) # Residuals are normally distributed (p = 0.6839).



# Type II ANOVA without interaction --> Table S3d
anova_brokenmussels1 <- aov(sqrt(broken) ~ corrosion + site, data = arch2)
Anova(anova_brokenmussels1, type = "II")  # There is NO significant effect of corrosion or site.

plot(anova_brokenmussels1, 1) # 3 points detected as potential outliers (Points 5, 12 and 23).

plot(anova_brokenmussels1, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_brokenmussels1_residuals <- residuals(object = anova_brokenmussels1)
shapiro.test(x = anova_brokenmussels1_residuals) # Residuals are normally distributed (p = 0.767). 





##
# STATS : Average nb of byssal threads ACROSS ALL SITES ====
##

# Assumptions #####
View(byssus1)

boxplot_byssus <- ggboxplot(byssus1, x = "site", y = "nb_byssal", color = "corrosion",
                            palette = col_corr)
boxplot_byssus

interaction_byssus <- ggline(byssus1, x = "site", y = "nb_byssal", color = "corrosion",
                             add = c("mean_se", "dotplot"),
                             palette = col_corr)
interaction_byssus

# Assumptions of normality

shapiro.test(byssus1$nb_byssal) # (!) Average number of byssal threads per mussel is NOT linear (p = 0.01199).
ggdensity(byssus1$nb_byssal) # Data distribution is left skewed.
ggqqplot(byssus1$nb_byssal) # Points fall approx. along the line.



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_byssus <- aov(nb_byssal ~ corrosion * site, data = byssus1)
Anova(anova_byssus, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_byssus, 1) # 3 points detected as potential outliers (Points 19, 73 and 77).
leveneTest(nb_byssal ~ corrosion * site, data = byssus1) # Variances are homogeneous (p = 0.2721).

plot(anova_byssus, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_byssus_residuals <- residuals(object = anova_byssus) 
shapiro.test(x = anova_byssus_residuals) # Residuals are normally distributed (p = 0.5166).



# Type II ANOVA without interaction --> Table S3e
anova_byssus1 <- aov(nb_byssal ~ corrosion + site, data = byssus1)
Anova(anova_byssus1, type = "II")  # There is a significant effect of site (p < 0.05) BUT NO significant effect of corrosion (p = 0.7868).

plot(anova_byssus1, 1) # 3 points detected as potential outliers (Points 19, 73 and 77).

plot(anova_byssus1, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_byssus1_residuals <- residuals(object = anova_byssus1)
shapiro.test(x = anova_byssus1_residuals) # Residuals are normally distributed (p = 0.0938). 



# Pairwise comparisons between sites --> Table S4c
summary(glht(anova_byssus1, linfct = mcp(site = "Tukey")))
# Results : JB > BR = PEd and MB > BR.





##
# SECTION: Architectural complexity at OWR ------------------------------------------
##

# GOAL : Is there a significant difference between corroded and non-corroded mussel beds, and eastern and western Perna perna lineages,
# at Old Woman's River, in terms of total number of mussels, number of live mussels, dead mussels and broken shells per quadrat,
# and average number of byssal threads / mussel ?

# Curate data set ====

# Select only Old Woman's River and the relevant variables in each data set
arch_owr <- arch[arch$site =="Old Woman's River",]
View(arch_owr)

byssus_owr <- byssus[byssus$site == "Old Woman's River",]
View(byssus_owr)

# Calculate habitat architectural complexity parameters (i.e., total number of mussels, number of live and dead mussels and broken shells, 
# average number of byssal threads/mussel) for each quadrat
arch_owr1 <- arch_owr %>%
  dplyr::group_by(site, corrosion, lineage, quadrat) %>%
  dplyr::summarize(
    total = tot_mussels, 
    live = live_perna + live_mytilus,
    dead = dead_perna + dead_mytilus,
    broken = nb_broken,
    mean_byssal = mean_byssal
  )
View(arch_owr1)

# Create a unique identifier for each quadrat in each data set
arch_owr1 <- arch_owr1 %>%
  dplyr::group_by(site, corrosion, lineage, quadrat) %>%
  tidyr::unite(col = "id", site, corrosion, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

byssus_owr <- byssus_owr %>%
  dplyr::group_by(site, corrosion, lineage, quadrat) %>%
  tidyr::unite(col = "id", site, corrosion, lineage, quadrat, sep = "_", remove = FALSE) %>%
  dplyr::mutate(id = as.factor(id))

# Transform the data set from wide to long format for graphical purposes
arch_owr2 <- arch_owr1[,1:9] %>% pivot_longer(cols = c('live', 'dead', 'broken'),
                                      names_to = 'mussels',
                                      values_to = 'nb')
View(arch_owr2)





# Visualize data set ====

# Barplot of the number of live, dead and broken mussels in each quadrat
barplot_mussels_owr <- ggplot(arch_owr2, aes(x = quadrat, y = nb, color = mussels, fill = mussels)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  facet_wrap(~ corrosion * lineage) +
  theme(legend.position = "bottom")
barplot_mussels_owr

# Boxplot of the number of live, dead and broken mussels for each corrosion status and site
boxplot_mussels_owr <- ggplot(arch_owr2, aes(x = site, y = nb, fill = corrosion)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ mussels) +
  scale_fill_manual(values = col_corr) +
  scale_x_discrete(labels = c("BR", "JB", "MB", "PEd")) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Site", y = "Mussel bed structure")
boxplot_mussels_owr

# Boxplot of the average number of byssal threads per mussel for each quadrat
boxplot_bt_owr <- ggplot(byssus_owr, aes(x = quadrat, y = nb_byssal, fill = corrosion)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  facet_wrap(~ lineage, ncol = 4) +
  scale_fill_manual(values = col_corr) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Quadrat", y = "Average number of byssal threads")
boxplot_bt_owr

# Boxplot of the average number of byssal threads per mussel for each corrosion status and site --> Figure 2e
ann_text_bt_owr1 <- data.frame(label = c("A", "B"), lineage = c("EAST", "WEST"), corrosion = "corroded")

boxplot_bt_owr1 <- ggplot(byssus_owr, aes(x = corrosion, y = nb_byssal, fill = corrosion)) + 
  geom_boxplot() +
  theme(plot.title = element_text(hjust = 0.5)) +
  scale_fill_manual(values = col_corr) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(limits = c(0, 410)) +
  geom_text(data = ann_text_bt_owr1, label = ann_text_bt_owr1$label, y = 400,
            size = 6, x = 1.5) +
  facet_wrap(~ lineage, nrow = 1) +
  theme(legend.position = "bottom", legend.title = element_text(face = "bold")) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average number of byssal threads",
       title = "Byssal threads")
boxplot_bt_owr1





##
# STATS : Total number of mussels at OWR ============================================
##

# Assumptions ####

View(arch_owr1)
table(arch_owr1$lineage, arch_owr1$corrosion) # Frequency table
# Because equal sample sizes within levels of independent groupings : ANOVA with balanced design.

# Visualize data set
boxplot_totmussels_owr <- ggboxplot(arch_owr1, x = "lineage", y = "total", color = "corrosion",
                                    palette = col_corr)
boxplot_totmussels_owr

interaction_totmusssels_owr <- ggline(arch_owr1, x = "lineage", y = "total", color = "corrosion",
                                      add = c("mean_se", "dotplot"),
                                      palette = col_corr)
interaction_totmusssels_owr

# Assumption of normality

shapiro.test(arch_owr1$total) # Total number of mussels (live + dead) is linear (p = 0.7832).
ggdensity(arch_owr1$total) # Data distribution looks linear.
ggqqplot(arch_owr1$total) # Points fall perfectly along the line.

View(arch_owr1)



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_totmussels_owr <- aov(total ~ corrosion * lineage, data = arch_owr1)
Anova(anova_totmussels_owr, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_totmussels_owr, 1) # 3 points detected as potential outliers (Points 2, 15 and 16).
leveneTest(total ~ corrosion * site, data = arch_owr1) # Variances are homogeneous (p = 0.4431).

plot(anova_totmussels_owr, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_totmussels_owr_residuals <- residuals(object = anova_totmussels_owr) 
shapiro.test(x = anova_totmussels_owr_residuals) # Residuals are normally distributed (p = 0.7707).



# Type II ANOVA without interaction (with outliers) --> Table S3f
anova_totmussels_owr1 <- aov(total ~ corrosion + lineage, data = arch_owr1)
Anova(anova_totmussels_owr1, type = "II") # There is NO significant effect of lineage (p = 0.09113) or corrosion (p = 0.11653).

plot(anova_totmussels_owr1, 1) # 3 points detected as potential outliers (Points 26, 27 and 28).

plot(anova_totmussels_owr1, 2) # Points fall approx. along the line.
anova_totmussels_owr1_residuals <- residuals(object = anova_totmussels_owr1)
shapiro.test(x = anova_totmussels_owr1_residuals) # Residuals are normally distributed (p = 0.7435). 





##
# STATS : Number of live mussels at OWR ============================================
##

# Assumptions ####

View(arch_owr1)
table(arch_owr1$lineage, arch_owr1$corrosion) # Frequency table

# Visualize data set
boxplot_livemussels_owr <- ggboxplot(arch_owr1, x = "lineage", y = "live", color = "corrosion",
                                    palette = col_corr)
boxplot_livemussels_owr

interaction_livemusssels_owr <- ggline(arch_owr1, x = "lineage", y = "live", color = "corrosion",
                                      add = c("mean_se", "dotplot"),
                                      palette = col_corr)
interaction_livemusssels_owr

# Assumption of normality

shapiro.test(arch_owr1$live) # Number of live mussels (Perna + Mytilus) is linear (p = 0.5383)
ggdensity(arch_owr1$live) # Data distribution looks linear.
ggqqplot(arch_owr1$live) # Points fall perfectly along the line.

View(arch_owr1)



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_livemussels_owr <- aov(live ~ corrosion * lineage, data = arch_owr1)
Anova(anova_livemussels_owr, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_livemussels_owr, 1) # 3 points detected as potential outliers (Points 1, 2 and 16).
leveneTest(live ~ corrosion * site, data = arch_owr1) # Variances are homogeneous (p = 0.4173).

plot(anova_livemussels_owr, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_livemussels_owr_residuals <- residuals(object = anova_livemussels_owr) 
shapiro.test(x = anova_livemussels_owr_residuals) # Residuals are normally distributed (p = 0.9061).



# Type II ANOVA without interaction (with outliers) --> Table S3g
anova_livemussels_owr1 <- aov(live ~ corrosion + lineage, data = arch_owr1)
Anova(anova_livemussels_owr1, type = "II") # There is NO significant effect of lineage (p = 0.17666) or corrosion (p = 0.05526).

plot(anova_livemussels_owr1, 1) # 3 points detected as potential outliers (Points 1, 2 and 16).

plot(anova_livemussels_owr1, 2) # Residuals fall approx. along the line.
anova_livemussels_owr1_residuals <- residuals(object = anova_livemussels_owr1)
shapiro.test(x = anova_livemussels_owr1_residuals) # Residuals are normally distributed (p = 0.8007). 





##
# STATS : Number of dead mussels at OWR ============================================
##

# Assumptions ####

# Visualize data set
boxplot_deadmussels_owr <- ggboxplot(arch_owr1, x = "lineage", y = "dead", color = "corrosion",
                                     palette = col_corr)
boxplot_deadmussels_owr

interaction_deadmusssels_owr <- ggline(arch_owr1, x = "lineage", y = "dead", color = "corrosion",
                                       add = c("mean_se", "dotplot"),
                                       palette = col_corr)
interaction_deadmusssels_owr

# Assumption of normality

shapiro.test(arch_owr1$dead) # (!) Number of dead mussels (Perna + Mytilus) is NOT linear (p = 0.04769)
ggdensity(arch_owr1$dead) # Data distribution is slightly left skewed.
ggqqplot(arch_owr1$dead) # Points fall approx. along the line.

View(arch_owr1)



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_deadmussels_owr <- aov(dead ~ corrosion * lineage, data = arch_owr1)
Anova(anova_deadmussels_owr, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_deadmussels_owr, 1) # 3 points detected as potential outliers (Points 3, 7 and 13).
leveneTest(dead ~ corrosion * site, data = arch_owr1) # Variances are homogeneous (p = 0.5441).

plot(anova_deadmussels_owr, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_deadmussels_owr_residuals <- residuals(object = anova_deadmussels_owr)
shapiro.test(x = anova_deadmussels_owr_residuals) # Residuals are normally distributed (p = 0.6644).



# Type II ANOVA without interaction (with outliers) --> Table S3h
anova_deadmussels_owr1 <- aov(dead ~ corrosion + lineage, data = arch_owr1)
Anova(anova_deadmussels_owr1, type = "II") # There is a significant effect of lineage (p = 0.01939) BUT NO effect of corrosion (p = 0.08778).

plot(anova_deadmussels_owr1, 1) # 3 points detected as potential outliers (Points 26, 27 and 28).

plot(anova_deadmussels_owr1, 2) # Residuals fall approx. along the line.
anova_deadmussels_owr1_residuals <- residuals(object = anova_deadmussels_owr1)
shapiro.test(x = anova_deadmussels_owr1_residuals) # Residuals are normally distributed (p = 0.2509). 



# Pairwise comparisons between Perna perna lineages --> Table S4d
summary(glht(anova_deadmussels_owr1, linfct = mcp(lineage = "Tukey")))



##
# STATS : Number of broken shells at OWR ============================================
##

# Assumptions ####

# Visualize data set
boxplot_broken_owr <- ggboxplot(arch_owr1, x = "lineage", y = "broken", color = "corrosion",
                                     palette = col_corr)
boxplot_broken_owr

interaction_broken_owr <- ggline(arch_owr1, x = "lineage", y = "broken", color = "corrosion",
                                       add = c("mean_se", "dotplot"),
                                       palette = col_corr)
interaction_broken_owr

# Assumptions of normality

shapiro.test(arch_owr1$broken) # Number of broken shells is linear (p = 0.053)
ggdensity(arch_owr1$broken) # Data distribution is left skewed and nearly bimodal.
ggqqplot(arch_owr1$broken) # Points fall approx. along the line.

View(arch_owr1)



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_brokenmussels_owr <- aov(broken ~ corrosion * lineage, data = arch_owr1)
Anova(anova_brokenmussels_owr, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_brokenmussels_owr, 1) # 3 points detected as potential outliers (Points 1, 9 and 15).
leveneTest(broken ~ corrosion * site, data = arch_owr1) # Variances are homogeneous (p = 0.473).

plot(anova_brokenmussels_owr, 2) # 3 points detected as potential outliers. Points fall approximately along the line.
anova_brokenmussels_owr_residuals <- residuals(object = anova_brokenmussels_owr) 
shapiro.test(x = anova_brokenmussels_owr_residuals) # Residuals are normally distributed (p = 0.0912).



# Type II ANOVA without interaction (with outliers) --> Table S3i
anova_brokenmussels_owr1 <- aov(broken ~ corrosion + lineage, data = arch_owr1)
Anova(anova_brokenmussels_owr1, type = "II") # There is NO significant effect of lineage (p = 0.1709) or corrosion (p = 0.2687).

plot(anova_brokenmussels_owr1, 1) # 3 points detected as potential outliers (Points 26, 27 and 28).

plot(anova_brokenmussels_owr1, 2) # Points fall approx. along the line.
anova_brokenmussels_owr1_residuals <- residuals(object = anova_brokenmussels_owr1)
shapiro.test(x = anova_brokenmussels_owr1_residuals) # Residuals are normally distributed (p = 0.1441). 





##
# STATS : Average nb of byssal threads at OWR ####
##

# Assumptions ####

# Visualize data set
boxplot_byssus_owr <- ggboxplot(byssus_owr, x = "lineage", y = "nb_byssal", color = "corrosion",
                                palette = col_corr)
boxplot_byssus_owr

interaction_byssus_owr <- ggline(byssus_owr, x = "lineage", y = "nb_byssal", color = "corrosion",
                                 add = c("mean_se", "dotplot"),
                                 palette = col_corr)
interaction_byssus_owr

# Assumption of normality

shapiro.test(byssus_owr$nb_byssal) # (!) The average number of byssal threads per mussel is NOT linear (p < 0.05).
ggdensity(byssus_owr$nb_byssal) # Data distribution is left skewed. Suitable transformations include square root, cube root and log.
ggqqplot(byssus_owr$nb_byssal) # Points fall approx. along the line.

shapiro.test(log(byssus_owr$nb_byssal)) # Log-transformed average number of byssal threads per mussel is linear (p = 0.7048).
ggdensity(log(byssus_owr$nb_byssal)) # Log-transformed data distribution looks linear.

View(arch_owr1)



# ANOVAs with outliers ####

# Type III ANOVA w/ interaction (with outliers)
anova_byssus_owr <- aov(log(nb_byssal) ~ corrosion * lineage, data = byssus_owr)
Anova(anova_byssus_owr, type = "III") # Interaction is not significant, so type II is more powerful.

plot(anova_byssus_owr, 1) # 3 points detected as potential outliers (Points 3, 11 and 42).
leveneTest(log(nb_byssal) ~ corrosion * lineage, data = byssus_owr) # Variances are homogeneous (p = 0.357).

plot(anova_byssus_owr, 2) # 3 points detected as potential outliers. Points fall approx. along the line.
anova_byssus_owr_residuals <- residuals(object = anova_byssus_owr)
shapiro.test(x = anova_byssus_owr_residuals) # Residuals are normally distributed (p = 0.1889).



# Type II ANOVA without interaction (with outliers)
anova_byssus_owr1 <- aov(log(nb_byssal) ~ corrosion + lineage, data = byssus_owr)
Anova(anova_byssus_owr1, type = "II") # There is NO significant effect of corrosion (p = 0.18862) or lineage (p = 0.05935).

plot(anova_byssus_owr1, 1) # 3 points detected as potential outliers (Points 3, 11 and 42).

plot(anova_byssus_owr1, 2) # Points fall approx. along the line.
anova_byssus_owr1_residuals <- residuals(object = anova_byssus_owr1)
shapiro.test(x = anova_byssus_owr1_residuals) # Residuals are normally distributed (p = 0.2127). 





# ANOVAs without outliers ####

# Omit identified outliers
byssus_owr_no <- byssus_owr[-c(3, 11, 42), ]
View(byssus_owr_no)
# Outliers OWR-W-I-3, OWR-E-I-2 and OWR-E-NI-2 were omitted from the analysis because of their higher number of byssal threads.

# Frequency table
table(byssus_owr_no$corrosion, byssus_owr_no$lineage)
# Data does not have equal sample sizes within levels of independent groupings.

# Visualize data set
boxplot_byssus_owr_no <- ggboxplot(byssus_owr_no, x = "lineage", y = "nb_byssal", color = "corrosion",
                                   palette = col_corr)
boxplot_byssus_owr_no

interaction_byssus_owr_no <- ggline(byssus_owr_no, x = "lineage", y = "nb_byssal", color = "corrosion",
                                    add = c("mean_se", "dotplot"),
                                    palette = col_corr)
interaction_byssus_owr_no



# Type III ANOVA w/ interaction (without outliers)
anova_byssus_owr_no <- aov(log(nb_byssal) ~ corrosion * lineage, data = byssus_owr_no)
Anova(anova_byssus_owr_no, type = "III") # Interaction is not significant, so type II ANOVA is more powerful.

plot(anova_byssus_owr_no, 1) # Points 3, 30 and 38 detected as outliers.
leveneTest(log(nb_byssal) ~ lineage * corrosion, data = byssus_owr_no) # Variances are homogeneous (p = 0.5837).

plot(anova_byssus_owr_no, 2) # Points fall approximately along the line.
anova_byssus_owr_no_residuals <- residuals(object = anova_byssus_owr_no) 
shapiro.test(x = anova_byssus_owr_no_residuals) # Residuals are normally distributed (p = 0.665).



# Type II ANOVA w/out interaction (without outliers) --> Table S3j
anova_byssus_owr_no1 <- aov(log(nb_byssal) ~ corrosion + lineage, data = byssus_owr_no)
Anova(anova_byssus_owr_no1, type = "II") # There is a significant effect of lineage (p < 0.05) BUT NOT corrosion (p = 0.058882).

plot(anova_byssus_owr_no1, 1) # Points 3, 30 and 38 detected as outliers.
leveneTest(log(nb_byssal) ~ corrosion * lineage, data = byssus_owr_no) # Variances are homogeneous (p = 0.5837).

plot(anova_byssus_owr_no1, 2) # Points fall approximately along the line.
anova_byssus_owr_no1_residuals <- residuals(object = anova_byssus_owr_no1) 
shapiro.test(x = anova_byssus_owr_no1_residuals) # Residuals are normally distributed (p = 0.6506).



# Pairwise comparisons between Perna lineages --> Table S4e
summary(glht(anova_byssus_owr_no1, linfct = mcp(lineage = "Tukey")))