## ---
##
## Script name: Data visualization of macrofaunal communities associated with mussel beds depending on their level of euendolithic
##              corrosion either along a biogeographical gradient or between Perna perna lineages
##
## Purpose of script: 
##    - EUENDOLITHS x BIOGEOGRAPHY : Compare the macrofaunal communities between corrosion status and among transplant sites
##      (excluding Old Woman's River) 
##        - Specific abundance
##        - Specific biomass
##    - EUENDOLITHS x PERNA LINEAGES : Compare the above variables between corrosion status and Perna lineages at Old Woman's River
##
## Find in the manuscript :
##
## Author: Alexia Dievart
##
## Date Created: 2022-01-05
## Dates Updated: 2022-01-06; 2023-02-05; 2024-02-12; 2024-04-30; 2024-10-21
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
               vegan,
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
                  plot.title = element_text(colour = "grey30", size = 16, hjust = 0, face = "bold"),
                  axis.text = element_text(colour = "black", size = 16),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 4), "mm"), size = 14),
                  strip.text.x = element_text(colour = "black", size = 16, hjust = 0.5, face = "bold"),
                  legend.position = "none"))

# Create vectors for graphical visualization
col_corr <- c("darkgrey", "brown3")
site_order <- c("Mosselbaai", "Brenton-on-Sea", "Jeffreysbaai", "Port Edward")
sites <- c("MB", "BR", "JB", "PEd")
inf <- c("C", "NC")
my_pch <- c(16, 17, 18, 15)




##
# SECTION: Load data ----
##

# Environmental variables ====
env <- read.csv("./RAW DATA/infauna_architectural complexity.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(env)

factor1 <- c("site", "corrosion", "lineage", "quadrat")
env[,factor1] <- lapply(env[,factor1], factor)

env1 <- env %>% drop_na(nb_broken) # Omit NA values
env1[is.na(env1)] <- 0
View(env1)

# Isolate desired variables in a new data set - i.e., tot. nb of mussels, live and dead mussels, byssal threads
arch <- env1[, c(1, 2, 4, 9:10, 11)]
add_arch <- env1 %>%
  summarise(
    live_mussels = live_perna + live_mytilus,
    dead_mussels = dead_perna + dead_mytilus
  )
arch <- arch %>%
  add_column(add_arch, .after = "tot_mussels")
View(arch)



# ACROSS ALL SITES: Environmental variables ####
arch1 <- arch %>%
  filter(site != "Old Woman's River")
glimpse(arch1)

arch1$site <- as.character(arch1$site) # Reset the number of levels in the factor 'site' to 4, without Old Woman's River
arch1$site <- as.factor(arch1$site)



# AT OWR: Environmental variables ####
arch_owr <- env1[, c(1:4, 9:11)]
add_arch_owr <- env1 %>%
  summarise(
    live_mussels = live_perna + live_mytilus,
    dead_mussels = dead_perna + dead_mytilus
  )
arch_owr <- arch_owr %>%
  add_column(add_arch_owr, .after = "tot_mussels")
View(arch_owr)

arch_owr <- arch_owr %>%
  filter(site == "Old Woman's River")
View(arch_owr)





# Community variables - e.g. abundance, biomass ====
community <- read.csv("./RAW DATA/infauna_community.csv", dec = ",", header = T, sep = ";")

factor2 <- c("site", "bioregion", "corrosion", "lineage", "quadrat", "higher_group", "species", "code_sp", "pic")
community[,factor2] <- lapply(community[,factor2], factor)
dplyr::glimpse(community)

View(community)



# ACROSS ALL SITES: Abundance ####
abund_long <- community[,c(1:5, 8:9)]

abund_long1 <- abund_long[abund_long$site != "Old Woman's River",]
View(abund_long1)

abund_wide <- spread(abund_long1, code_sp, abundance)
abund_wide[is.na(abund_wide)] <- 0
View(abund_wide)

abund_wide1 <- abund_wide[,c(6:105)] # Isolate abundances from environmental variables
View(abund_wide1)



# ACROSS ALL SITES: Biomass ####
biom_long <- community[,c(1:5, 8, 10)]

biom_long1 <- biom_long[biom_long$site != "Old Woman's River",]
View(biom_long1)

biom_wide <- spread(biom_long1, code_sp, biomass)
biom_wide[is.na(biom_wide)] <- 0
View(biom_wide)

biom_wide1 <- biom_wide[,c(6:105)]
View(biom_wide1)



# AT OWR: Abundance ####
abund_long_owr <- abund_long[abund_long$site == "Old Woman's River",]
View(abund_long_owr)

abund_wide_owr <- spread(abund_long_owr, code_sp, abundance)
abund_wide_owr[is.na(abund_wide_owr)] <- 0
View(abund_wide_owr)

abund_wide_owr1 <- abund_wide_owr[,c(6:75)]



# AT OWR: Biomass ####
biom_long_owr <- biom_long[biom_long$site == "Old Woman's River",]
View(biom_long_owr)

biom_wide_owr <- spread(biom_long_owr, code_sp, biomass)
biom_wide_owr[is.na(biom_wide_owr)] <- 0
View(biom_wide_owr)

biom_wide_owr1 <- biom_wide_owr[, c(6:75)]
View(biom_wide_owr1)





##
# SECTION: Visualize data ----
##

# FIGURE 3 ====
layout.matrix <- matrix(c(1, 4, 2, 5, 3, 6), nrow = 2, ncol = 3)

layout.matrix

layout(mat = layout.matrix,
       heights = c(1, 1), # Heights of the two rows
       widths = c(3, 3, 2)) # Widths of the two columns

layout.show(6)

par(mfrow = c(2,3), oma = c(1,6,1,1))



# FIGURE 3a - ACROSS ALL SITES : nMDS on abundances ####

# Load Bray-Curtis dissimilarity matrix on abundances across all sites
abund_bray <- metaMDS(abund_wide1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
abund_bray # Stress = 0.122 shows this nMDS is a fair fit to the data.

abund_fit <- envfit(abund_bray, arch1, permu = 999)

abund_fit_adj <- abund_fit
pvals_adj <- p.adjust(abund_fit$vectors$pvals, method = 'bonferroni')
abund_fit_adj$vectors$pvals <- pvals_adj
abund_fit_adj



# PLOT: nMDS on abundances between sites and infestation levels
par(mar = c(5, 5, 3, 1))
plot(abund_bray, type = "none", xlim = c(-0.7, 1), ylim = c(-0.8, 0.8),
     cex.lab = 1.5, cex.axis = 1.5)
title(main = "Abundances", cex.main = 2, adj = 0)
points(abund_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 15)[arch1$site],
       col = col_corr[arch1$corrosion])
plot(abund_fit_adj, p.max = 0.05, col = 'navy', cex = 1.2)

legend(-1.3, 0.95, "Stress = 0.122", bty = "n", cex = 1.5)

mtext("A", side = 3, line = 0.5, cex = 2, adj = -0.15, col = "grey30")






# FIGURE 3b - ACROSS ALL SITES : nMDS on biomasses ####

# Load Bray-Curtis dissimilarity matrix on biomasses across all sites
biom_bray <- metaMDS(biom_wide1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
biom_bray # Stress = 0.157 shows this nMDS is a bit suspect fit to the data.

biom_fit <- envfit(biom_bray, arch1, permu = 999)

biom_fit_adj <- biom_fit
pvals_adj <- p.adjust(biom_fit$vectors$pvals, method = 'bonferroni')
biom_fit_adj$vectors$pvals <- pvals_adj
biom_fit_adj



# PLOT: nMDS on biomasses between sites and infestation levels
par(mar = c(5, 5, 3, 1))
ordiplot(biom_bray, type = "none", xlim = c(-1.0, 1.4), ylim = c(-1.2, 1.2),
         cex.lab = 1.5, cex.axis = 1.5)
title(main = "Biomasses", cex.main = 2, adj = 0)
points(biom_bray, display = "sites", cex = 3, pch = c(16, 17, 18, 15)[arch1$site],
       col = col_corr[arch1$corrosion])
plot(biom_fit_adj, p.max = 0.05, col = 'navy', cex = 1.2) 
legend(-1.95, 1.4, "Stress = 0.157", bty = "n", cex = 1.5)

mtext("B", side = 3, line = 0.5, cex = 2, adj = -0.15, col = "grey30")





# FIGURE 3 (top row) - Legend across all sites ####
par(mar = c(0, 0, 0, 0))
plot.new()
legend(0.2, 0.9, legend = levels(arch1$site),
       pch = c(16, 17, 18, 15),
       col = c("black", 'black', "black", "black"),
       cex = 1.7, box.lty = 0, bg = "transparent", pt.cex = 2.8,
       title = as.expression(bquote(bold("Sampling sites"))),
       title.adj = 0.03)
legend(0.2, 0.45, legend = c("corroded", "non-corroded"),
       pch = c(20, 20),
       col = c("darkgrey", "brown3"),
       cex = 1.7, box.lty = 0, bg = "transparent", pt.cex = 4,
       title = as.expression(bquote(bold("Corrosion level"))),
       title.adj = 0.03)





# FIGURE 3c - AT OWR : nMDS on abundances ####

# Load Bray-Curtis dissimilarity matrix on abundances at OWR
abund_bray_owr <- metaMDS(abund_wide_owr1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
abund_bray_owr # Stress = 0.117 shows this nMDS is a fair fit to the data.

abund_fit_owr <- envfit(abund_bray_owr, arch_owr, permu = 999)

abund_fit_adj_owr <- abund_fit_owr
pvals_adj_owr <- p.adjust(abund_fit_owr$vectors$pvals, method = 'bonferroni')
abund_fit_adj_owr$vectors$pvals <- pvals_adj_owr
abund_fit_adj_owr



# PLOT: nMDS on abundances between Perna lineages and infestation levels
par(mar = c(5, 5, 3, 1))
ordiplot(abund_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0),
         cex.lab = 1.5, cex.axis = 1.5)
title(main = "Abundances", cex.main = 2, adj = 0)
points(abund_bray_owr, display = "sites", cex = 3, pch = c('+', 'o')[arch_owr$lineage],
       col = col_corr[arch_owr$corrosion])
ordihull(abund_bray_owr, groups = arch_owr$corrosion, lty = "dotted", lwd = 3,
         col = col_corr)
plot(abund_fit_adj_owr, p = 0.05, col = 'navy', cex = 1.2)
legend(-1.76,1.2, "Stress = 0.116", bty = "n", cex = 1.5)

mtext("C", side = 3, line = 0.5, cex = 2, adj = -0.15, col = "grey30")





# FIGURE 3f - AT OWR : nMDS on biomasses ####

# Load Bray-Curtis dissimilarity matrix on abundances at OWR
biom_bray_owr <- metaMDS(biom_wide_owr1^0.25, distance = "bray", trace = F, trymax = 100)
# w/ trymax = number of random starts AND trace = F tells R to no output all of these random starts
biom_bray_owr # Stress = 0.116 shows this nMDS is a fair fit to the data

biom_fit_owr <- envfit(biom_bray_owr, arch_owr, permu = 999)

biom_fit_adj_owr <- biom_fit_owr
pvals_adj_owr <- p.adjust(biom_fit_owr$vectors$pvals, method = 'bonferroni')
biom_fit_adj_owr$vectors$pvals <- pvals_adj_owr
biom_fit_adj_owr



# PLOT: nMDS on biomasses between Perna lineages and infestation levels
par(mar = c(5, 5, 3, 1))
ordiplot(biom_bray_owr, type = "none", xlim = c(-1.0, 1.0), ylim = c(-1.0, 1.0),
         cex.lab = 1.5, cex.axis = 1.5)
title(main = "Biomasses", cex.main = 2, adj = 0)
ordihull(biom_bray_owr, groups = arch_owr$corrosion, lty = "dotted", lwd = 3,
         col = col_corr)
points(biom_bray_owr, display = "sites", cex = 3, pch = c('+', 'o')[arch_owr$lineage],
       col = col_corr[arch_owr$corrosion])
plot(biom_fit_adj_owr, p = 0.05, col = 'navy', cex = 1.2)
legend(-1.8, 1.2, "Stress = 0.116", bty = "n", cex = 1.5)

mtext("D", side = 3, line = 0.5, cex = 2, adj = -0.15, col = "grey30")





# FIGURE 3 (bottom row) - Legend at OWR ####
par(mar = c(0, 0, 0, 0))
plot.new()
legend(0.18, 0.9, legend = c("Eastern lineage", "Western lineage"),
       pch = c('+', 'o'),
       col = c("black", 'black', "black", "black"),
       cex = 1.7, box.lty = 0, bg = "transparent", pt.cex = 2.8,
       title = as.expression(bquote(bolditalic(.("Perna perna")~bold(.("lineage"))))),
       title.adj = 0.08)
legend(0.2, 0.6, legend = c("corroded", "non-corroded"),
       pch = c(20, 20),
       col = c("darkgrey", "brown3"),
       cex = 1.7, box.lty = 0, bg = "transparent", pt.cex = 4,
       title = as.expression(bquote(bold("Corrosion level"))),
       title.adj = 0.05)

mtext("Across all sites", side = 2, line = 78, cex = 2, adj = 2.7, col = "black")
mtext("Old Woman's River", side = 2, line = 78, cex = 2, adj = 0.65, col = "black")




