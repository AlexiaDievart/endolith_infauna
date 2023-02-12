## ---------------------------
##
## Script name: Infaunal communities (community analysis) associated with manipulated mussel beds
##              depending on their level of euendolithic infestation and their transplant site
##
## Purpose of script: 
##
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
env <- read.csv("./RAW DATA/infauna_description2.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(env)

factor1 <- c("site", "infestation", "lineage", "quadrat")
env[,factor1] <- lapply(env[,factor1], factor)

env <- env %>% drop_na(nb_broken)

View(env)

env1 <- env[, c(1, 2, 4:7)]
env2 <- env1 %>%
  filter(site != "Old Woman's River")
factor3 <- c("site", "infestation", "quadrat")
env2[,factor3] <- lapply(env2[,factor3], factor)

env_owr <- env %>%
  filter(site == "Old Woman's River")
env_owr1 <- env_owr[, c(1:7)]
env_owr1[,factor1] <- lapply(env_owr1[,factor1], factor)
View(env_owr1)


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
View(abund_long1)

abund_wide <- spread(abund_long, species, abundance)
abund_wide[is.na(abund_wide)] <- 0
View(abund_wide)

abund_wide1 <- abund_wide[abund_wide$site != "Old Woman's River",]
factor4 <- c("site", "bioregion", "infestation", "lineage", "quadrat")
abund_wide1[,factor4] <- lapply(abund_wide1[,factor4], factor)

View(abund_wide1)

abund_long_owr <- abund_long[abund_long$site == "Old Woman's River",]
View(abund_long_owr)

abund_wide_owr <- spread(abund_long_owr, species, abundance)
abund_wide_owr[is.na(abund_wide_owr)] <- 0
View(abund_wide_owr)

# Isolate BIOMASS #############################################################################
biom_long <- community[,c(1:3, 5,7, 9)]
biom_long1 <- biom_long[biom_long$site != "Old Woman's River",]
View(biom_long1)

biom_wide <- spread(biom_long1, species, biomass)
biom_wide[is.na(biom_wide)] <- 0
View(biom_wide)

biom_long_owr <- biom_long[biom_long$site == "Old Woman's River",]
View(biom_long_owr)

biom_wide_owr <- spread(biom_long_owr, species, biomass)
biom_wide_owr[is.na(biom_wide_owr)] <- 0
View(biom_wide_owr)



###############################################################################################
# Section: Community analysis on ABUNDANCE across all sites -----------------------------------
###############################################################################################

# Curate the dataset ####
View(abund_wide1)

# Isolate abundances from environmental variables
abund_wide2 <- abund_wide1[,c(6:120)]
View(abund_wide2)
colnames(abund_wide2) <- make.cepnames(colnames(abund_wide2))

# Calculate the distance matrix ####
# Compute distance matrix using Bray-Curtis on transformed abundances
range(abund_wide2)
range(abund_wide2^0.5) 
range(abund_wide2^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
abund_dist <- vegdist(abund_wide2^0.25, method='bray')
abund_dist

# Calculate nMDS ####

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.005 --> Good fit

abund_bray <- metaMDS(abund_wide2^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
abund_bray # Stress = 0.123 shows this nMDS is a fair fit to the data.

# Exploring the results of the nMDS ####
plot(abund_bray, type = "t")
stressplot(abund_bray)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(abund_bray)
plot(abund_bray, type = "t", main = "goodness of fit")
points(abund_bray, display = "sites", cex = gof*100)

abund_fit <- envfit(abund_bray, env1, permu = 999) # Caution: envfit does not allow missing values
abund_fit
# The total number of mussels is strongly related to the first ordination axes (r2 = 0.36)
# Site is a significant predictor of infaunal community on abundance

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

# Plot nMDS w/ infestation levels and sites ####
ordiplot(abund_bray, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(abund_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 5)[env2$site],
       col = my_colors[env2$infestation])
ordihull(abund_bray, groups = env2$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.35, 0.3, legend = c(levels(env2$site), levels(env2$infestation)),
       pch = c(16, 17, 18, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.70, 1.25, "Stress = 0.123", bty = "n", cex = 2)

# PERMANOVA ####

# Run PERMANOVA
abund_pmv <- adonis2(abund_wide2^0.25 ~ infestation * site, data = env2,
                     permutations = 999, method = 'bray')
abund_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities DO NOT differ significantly between infestation levels (p = 0.366).
# But communities differ between sites (p = 0.001)

# Plot permuted F-values
densityplot(permustats(abund_pmv))

# Plot the average distance to median for infestation levels
abund_disp_inf <- betadisper(abund_dist, group = env2$infestation)
boxplot(abund_disp_inf)
# No difference in dispersion detected between infestation levels from the boxplot.

# Is there a significant difference in dispersion between infestation levels ?
anova(abund_disp_inf) 
permutest(abund_disp_inf)
# No significant difference in dispersion between infestation levels
plot(abund_disp_inf, hull = F, ellipse = T, segments = T, seg.col = my_colors, seg.lwd = 1,
     pch = c(16, 16), col = my_colors, label = T, cex = 2, lwd = 2,
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5),
     xlab = "Dimension 1 (31.82 %)", ylab = "Dimension 2 (23.93 %)",
     main = ""
)
round(100 * abund_disp_inf$eig / sum(abund_disp_inf$eig), 2)

# Plot the average distance to median for sites
abund_disp_site <- betadisper(abund_dist, group = env2$site)
boxplot(abund_disp_site)
# There is a clear difference in dispersion detected between sites from the boxplot.

anova(abund_disp_site)
permutest(abund_disp_site)
plot(abund_disp_site, hull = F, ellipse = T,
     pch = c(16:18, 8, 5), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (31.82 %)", ylab = "Dimension 2 (23.93 %)",
     xlim = c(-0.4, 0.4), ylim = c(-0.5, 0.5),
     main = "")
round(100 * abund_disp_site$eig / sum(abund_disp_site$eig), 2)

## Cluster analysis ####
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

## SIMPER on infestation
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)
abund_simp_inf <- simper(abund_wide2^0.25, group = env2$infestation)
abund_simp_inf
summary(abund_simp_inf)

## SIMPER on sites
abund_simp_site <- simper(abund_wide2^0.25, group = env2$site)
abund_simp_site
summary_simper <- summary(abund_simp_site)
View(summary_simper$`Brenton-on-Sea_Mosselbaai`)
View(summary_simper$`Brenton-on-Sea_Jeffreysbaai`)
View(summary_simper$`Brenton-on-Sea_Port Edward`)
View(summary_simper$Jeffreysbaai_Mosselbaai)
View(summary_simper$`Jeffreysbaai_Port Edward`)
View(summary_simper$`Mosselbaai_Port Edward`)



###############################################################################################
# Section: Community analysis on ABUNDANCE at OWR ---------------------------------------------
###############################################################################################

abund_wide_owr1 <- abund_wide_owr[,c(6:75)]
View(abund_wide_owr1)
View(env_owr1)

# Calculate the distance matrix #####
# Compute distance matrix using Bray-Curtis on transformed abundances
range(abund_wide_owr1)
range(abund_wide_owr1^0.5) 
range(abund_wide_owr1^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
abund_dist_owr <- vegdist(abund_wide_owr1^0.25, method='bray')
abund_dist_owr

# Calculate nMDS ####

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.005 --> Good fit

abund_bray_owr <- metaMDS(abund_wide_owr1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
abund_bray_owr # Stress = 0.116 shows this nMDS is a fair fit to the data.

### Exploring the results of the nMDS
plot(abund_bray_owr, type = "t")
stressplot(abund_bray_owr)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(abund_bray_owr)
plot(abund_bray_owr, type = "t", main = "goodness of fit")
points(abund_bray_owr, display = "sites", cex = gof*100)

abund_fit_owr <- envfit(abund_bray_owr, env_owr1, permu = 999) # Caution: envfit does not allow missing values
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

# Plot nMDS w/ infestation levels and lineages ####

ordiplot(abund_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(abund_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env_owr1$lineage],
       col = my_colors[env_owr1$infestation])
ordihull(abund_bray_owr, groups = env_owr1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.1, -0.3, legend = c(levels(env_owr1$lineage), levels(env_owr1$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.45, 1.25, "Stress = 0.116", bty = "n", cex = 2)

# PERMANOVA ####
# Run PERMANOVA
abund_owr_pmv <- adonis2(abund_wide_owr1^0.25 ~ infestation * lineage, data = env_owr1,
                         permutations = 999, method = 'bray')
abund_owr_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels

# Plot permuted F-values
densityplot(permustats(abund_owr_pmv))

# Plot the average distance to median for infestation levels
abund_disp_inf_owr <- betadisper(abund_dist_owr, group = env_owr1$infestation)
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
     xlab = "Dimension 1 (42.64 %)", ylab = "Dimension 2 (13.21 %)",
     main = ""
)
round(100 * abund_disp_inf_owr$eig / sum(abund_disp_inf_owr$eig), 2)

# Plot the average distance to median for sites
abund_disp_lineage <- betadisper(abund_dist_owr, group = env_owr1$lineage)
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


## Cluster analysis ####
abund_clust_owr <- hclust(abund_dist_owr, "ward.D2")

abund_inf_owr <- rev(levels(env_owr1$infestation))
abund_dend_owr <- as.dendrogram(abund_clust_owr)

labels_colors(abund_dend_owr) <- my_colors[(env_owr1$infestation)[order.dendrogram(abund_dend_owr)]]

labels(abund_dend_owr) <- paste(env_owr1$lineage[order.dendrogram(abund_dend_owr)],
                                " (", labels(abund_dend_owr), ")", sep = "")
abund_dend_owr <- set(abund_dend_owr, "labels_cex", 0.9)
abund_dend_owr <- hang.dendrogram(abund_dend_owr, hang_height = 0.005)
plot(abund_dend_owr, nodePar = list(cex = 0.007))
legend("topright", legend = levels(env_owr1$infestation), fill = my_colors)

## SIMPER on infestation
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)
abund_simp_inf_owr <- simper(abund_wide_owr1^0.25, group = env_owr1$infestation)
abund_simp_inf_owr
summary(abund_simp_inf_owr)

## SIMPER on sites
abund_simp_lineage <- simper(abund_wide_owr1^0.25, group = env_owr1$lineage)
abund_simp_lineage
summary(abund_simp_lineage)
View(summary_simper)





###############################################################################################
# Section: Community analysis on BIOMASS across all sites -------------------------------------
###############################################################################################

# Curate dataset ####
View(biom_wide1)
biom_wide1 <- biom_wide[,c(5:104)]
View(env2)

# Calculate the distance matrix ####
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
biom_bray # Stress = 0.157 shows this nMDS is a bit suspect fit to the data.
View(biom_bray)

### Exploring the results of the nMDS
plot(biom_bray, type = "t")
stressplot(biom_bray)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(biom_bray)
plot(biom_bray, type = "t", main = "goodness of fit")
points(biom_bray, display = "sites", cex = gof*100)

biom_fit <- envfit(biom_bray, env2, permu = 999) # Caution: envfit does not allow missing values
biom_fit
# The total number of mussels is strongly related to the first and second axis (r2 = 0.35).
# Site is a strong predictor of communities

biom_fit_adj <- biom_fit
pvals_adj <- p.adjust(biom_fit$vectors$pvals, method = 'bonferroni')
biom_fit_adj$vectors$pvals <- pvals_adj
biom_fit_adj
# W/ a Bonferroni correction, the p-values are still the same.

plot(biom_bray, display = "sites")
plot(biom_fit_adj, p.max = 0.05) 
# Infestation is not significant. 

# Plot nMDS w/ infestation levels and sites ####

ordiplot(biom_bray, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0))
points(biom_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 5)[env2$site],
       col = my_colors[env2$infestation])
ordihull(biom_bray, groups = env2$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.1, 0.3, legend = c(levels(env2$site), levels(env2$infestation)),
       pch = c(16, 17, 18, 5, 15, 15), 
       col = c("black", 'black', "black", "black", "black", "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.45, 1.25, "Stress = 0.157", bty = "n", cex = 2)

## PERMANOVA ####
# Run PERMANOVA
biom_pmv <- adonis2(biom_wide1^0.25 ~ infestation * site, data = env2,
                    permutations = 999, method = 'bray')
biom_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities DO NOT differ significantly between infestation levels (p = 0.82).
# But communities differ between sites (p = 0.001)

# Plot permuted F-values
densityplot(permustats(biom_pmv))

# Plot the average distance to median for infestation levels
biom_disp_inf <- betadisper(biom_dist, group = env2$infestation)
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
     xlab = "Dimension 1 (36.50 %)", ylab = "Dimension 2 (21.77 %)",
     main = ""
)
round(100 * biom_disp_inf$eig / sum(biom_disp_inf$eig), 2)

# Plot the average distance to median for sites
biom_disp_site <- betadisper(biom_dist, group = env2$site)
boxplot(biom_disp_site)
biom_disp_site
# There is a clear difference in dispersion detected between sites from the boxplot.
anova(biom_disp_site)
permutest(biom_disp_site)
plot(biom_disp_site, hull = F, ellipse = T,
     pch = c(16:18, 8, 5), col = "black", label = T, cex = 2,
     xlab = "Dimension 1 (36.50 %)", ylab = "Dimension 2 (21.77 %)",
     xlim = c(-0.5, 0.5), ylim = c(-0.5, 0.5),
     main = "")

# Cluster analysis - all sites on biomass ####
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

## SIMPER on infestation
# SIMPER shows you which species are responsible for the difference between observed groups
# How to read the results of SIMPER analysis ?
# contr: contribution to dissimilarity between infested and non infested
# sd: standard deviation of contribution (is the species response consistent?)
# ratio: ratio between 'contr' and 'sd' (high ratio = high, consistent contribution)
# av.: average abundance per group
# cumsum: cumulative contribution (rule of thumb = species till 70% are investigated)
biom_simp_inf <- simper(biom_wide1^0.25, group = env2$infestation)
biom_simp_inf
summary(biom_simp_inf)

## SIMPER on sites
biom_simp_site <- simper(biom_wide1^0.25, group = env2$site)
biom_simp_site
summary_biom_simper <- summary(biom_simp_site)
View(summary_biom_simper$BoS_MB)
View(summary_biom_simper$BoS_JB)
View(summary_biom_simper$BoS_PEd)
View(summary_biom_simper$JB_MB)
View(summary_biom_simper$JB_PEd)
View(summary_biom_simper$MB_PEd)



## Infested vs Non-infested / Western vs Eastern at Old Woman's River #####################

View(env_owr1)
View(biom_wide_owr1)
biom_wide_owr1 <- biom_wide_owr[, c(6:75)]

### Calculate the distance matrix
# Compute distance matrix using Bray-Curtis on transformed abundances
range(biom_wide_owr1)
range(biom_wide_owr1^0.5) 
range(biom_wide_owr1^0.25) # Range < 10
# If working on community data, down-weigth the abundant species to a range between 0 and 10
biom_dist_owr <- vegdist(biom_wide_owr1^0.25, method='bray')
biom_dist_owr

### Calculate nMDS

# How to interpret a nMDS ?
# If stress value approaching 0.3 --> Ordination is arbitrary
# If stress value around or above 0.2 --> Suspect
# If stress value equal to or below 0.1 --> Fair
# If stress valye equal to or below 0.005 --> Good fit

biom_bray_owr <- metaMDS(biom_wide_owr1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
biom_bray_owr # Stress = 0.116 shows this nMDS is a fair fit to the data.

### Exploring the results of the nMDS
plot(biom_bray_owr, type = "t")
stressplot(biom_bray_owr)
# No evidence of large scatter around the line suggested that our nMDS is not a bad fit.
gof <- goodness(biom_bray_owr)
plot(biom_bray_owr, type = "t", main = "goodness of fit")
points(biom_bray_owr, display = "sites", cex = gof*100)

biom_fit_owr <- envfit(biom_bray_owr, env_owr1, permu = 999) # Caution: envfit does not allow missing values
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
points(biom_bray_owr, display = "sites", cex = 3, pch = c(16, 17)[env_owr1$lineage],
       col = my_colors[env_owr1$infestation])
ordihull(biom_bray_owr, groups = env_owr1$infestation, lty = "dotted", lwd = 3,
         col = my_colors)
legend(-2.1, -0.3, legend = c(levels(env_owr1$lineage), levels(env_owr1$infestation)),
       pch = c(16, 17, 15, 15), 
       col = c("black", 'black', "darkgrey", "brown3"),
       cex = 1.2, box.lty = 0, bg = "transparent", pt.cex = 3)
legend(-2.45, 1.25, "Stress = 0.116", bty = "n", cex = 2)

## PERMANOVA
# Run PERMANOVA
biom_owr_pmv <- adonis2(biom_wide_owr1^0.25 ~ infestation * lineage, data = env_owr1,
                        permutations = 999, method = 'bray')
biom_owr_pmv
# In front of the tilde --> dependent variable = community matrix (w/ a square-root transformation)
# After the tilde --> independent variable = infestation levels
# Communities DO NOT differ significantly between infestation levels (p = 0.82).
# But communities differ between sites (p = 0.001)

# Plot permuted F-values
densityplot(permustats(biom_owr_pmv))

# Plot the average distance to median for infestation levels
biom_disp_inf_owr <- betadisper(biom_dist_owr, group = env_owr1$infestation)
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
biom_disp_lineage <- betadisper(biom_dist_owr, group = env_owr1$lineage)
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

biom_inf_owr <- rev(levels(env_owr1$infestation))
biom_dend_owr <- as.dendrogram(biom_clust_owr)

labels_colors(biom_dend_owr) <- my_colors[(env_owr1$infestation)[order.dendrogram(biom_dend_owr)]]

labels(biom_dend_owr) <- paste(env_owr1$lineage[order.dendrogram(biom_dend_owr)],
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