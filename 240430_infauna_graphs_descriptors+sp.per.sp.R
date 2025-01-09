## ---
##
## Script name: Data visualization of within-bed architectural complexity and general macrofaunal community descriptors
##              depending on their level of euendolithic corrosion, either along a biogeographical gradient or between genetic lineages
##
## Purpose of script: 
##
## Find in the manuscript :
##
## Author: Alexia Dievart
##
## Date Created: 2022-01-05
## Dates Updated: 2022-01-06; 2023-02-05; 2023-12-12; 2024-02-05; 2024-04-30; 2024-10-21
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
                  axis.text = element_text(colour = "black", size = 14),
                  axis.text.x = element_text(angle = 0, hjust = 0.5),
                  axis.title = element_text(size = 14),
                  axis.title.x = element_blank(),
                  axis.title.y = element_text(margin = unit(c(0, 4, 0, 0), "mm")),
                  plot.title = element_text(face = "bold", size = 14),
                  strip.text = element_text(size = 14),
                  legend.position = "none"))

# Create vectors for graphical visualization
col_corr <- c("darkgrey", "brown3")
site_order <- c("Mosselbaai", "Brenton-on-Sea", "Jeffreysbaai", "Port Edward")
sites <- c("MB", "BR", "JB", "PEd")
inf <- c("C", "NC")





##
# Section: Load data ----
##

# Environmental variables - number of byssal threads ====
byssus <- read.csv("./RAW DATA/infauna_byssal.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(byssus)

factor1 <- c("site", "corrosion", "lineage", "quadrat")
byssus[,factor1] <- lapply(byssus[,factor1], factor)
byssus <- na.omit(byssus) # Drop the quadrats for which byssal threads were not counted
View(byssus)

# ACROSS ALL SITES: Number of byssal threads
byssus1 <- byssus[byssus$site != "Old Woman's River",]
View(byssus1)

# AT OWR : Number of byssal threads
byssus_owr <- byssus %>%
  filter(site == "Old Woman's River")
View(byssus_owr)



# Community variables - e.g. abundance, biomass ====
community <- read.csv("./RAW DATA/infauna_community.csv", dec = ",", header = T, sep = ";")
dplyr::glimpse(community)

factor2 <- c("site", "bioregion", "corrosion", "lineage", "quadrat", "higher_group", "species", "code_sp")
community[,factor2] <- lapply(community[,factor2], factor)

# Drop the lines without count (abundance = NA)
community <- community %>% drop_na(abundance)
View(community)



# ACROSS ALL SITES : Community variables
community1 <- community[community$site != "Old Woman's River",]
View(community1)

# ACROSS ALL SITES : Total macrofaunal abundance
abund_long1 <- community1[,c(1:5,8:9)] # Long format
View(abund_long1)

abund_wide <- spread(abund_long1, code_sp, abundance) # Wide format
abund_wide[is.na(abund_wide)] <- 0
View(abund_wide)

abund_tot <- abund_long1 %>%
  na.omit(abund_long) %>%
  dplyr::group_by(site, corrosion, quadrat) %>%
  summarize(ab_total = sum(abundance))
View(abund_tot)

abund_tot1 <- abund_tot  %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    ab_tot_mean = mean(ab_total),
    ab_tot_sd = sd(ab_total)
  )
View(abund_tot1)



# ACROSS ALL SITES : Total macrofaunal biomass
biom_long1 <- community1[,c(1:5,8,10)] # Long format
View(biom_long1)

biom_wide <- spread(biom_long1, code_sp, biomass) # Wide format
biom_wide[is.na(biom_wide)] <- 0
View(biom_wide)

biom_tot <- biom_long1 %>%
  na.omit(biom_long) %>%
  dplyr::group_by(site, corrosion, quadrat) %>%
  summarize(biom_total = sum(biomass))
View(biom_tot) # (!) Biomass is expressed in mg.

biom_tot1 <- biom_tot  %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    biom_tot_mean = mean(biom_total),
    biom_tot_sd = sd(biom_total)
  )
View(biom_tot1)



# ACROSS ALL SITES : Species richness
richness <- abund_long1 %>%
  dplyr::group_by(site, corrosion, quadrat) %>%
  summarise(count = n_distinct(code_sp))
View(richness)

richness1 <- richness %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
View(richness1)



# ACROSS ALL SITES : Species diversity w/ Shannon diversity (H')
library(vegan)

abund_wide1 <- abund_wide
View(abund_wide1)

shannon <- diversity(abund_wide[,6:105], "shannon")
abund_wide1 <- abund_wide1 %>%
  add_column(shannon, .after = "quadrat")

shannon1 <- abund_wide1 %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    shannon_mean = mean(shannon),
    shannon_sd = sd(shannon)
  )
View(shannon1)



# ACROSS ALL SITES : Species diversity w/ Simpson diversity (λ)
View(abund_wide1)

simpson <- diversity(abund_wide[,6:105], index = "simpson")
abund_wide1 <- abund_wide1 %>%
  add_column(simpson, .after = "shannon")

simpson1 <- abund_wide1 %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson)
  )
View(simpson1)



# ACROSS ALL SITES : Pielou's evenness (J)
View(abund_wide1) ; View(richness)

abund_wide1 <- abund_wide1 %>%
  add_column(richness = richness$count, .after = "simpson")
pielou <- abund_wide1$shannon/log(abund_wide1$richness) # Formula : J = H'/ln(S)
abund_wide1 <- abund_wide1 %>%
  add_column(pielou, .after = "richness")

pielou1 <- abund_wide1 %>%
  dplyr::group_by(site, corrosion) %>%
  summarise(
    pielou_mean = mean(pielou),
    pielou_sd = sd(pielou)
  )
View(pielou1)



# ACROSS ALL SITES: Contribution
cumsum_abund <- read.csv("./RAW DATA/infauna_cumsum_sites_abund.csv", dec = ",", header = T, sep = ";")
View(cumsum_abund)

cumsum_abund_long <- gather(cumsum_abund, site, contr, BR.vs.JB:MB.vs.PEd, factor_key = T)
View(cumsum_abund_long)

cumsum_biom <- read.csv("./RAW DATA/infauna_cumsum_sites_biom.csv", dec = ",", header = T, sep = ";")
View(cumsum_biom)

cumsum_biom_long <- gather(cumsum_biom, site, contr, BR.vs.JB:MB.vs.PEd, factor_key = T)
View(cumsum_biom_long)




# AT OWR : Community variables
community_owr <- community %>%
  filter(site == "Old Woman's River")
View(community_owr)

# AT OWR : Total infaunal abundance
abund_long_owr <- community_owr[,c(1:5,8,9)] # Long format
View(abund_long_owr)

abund_wide_owr <- spread(abund_long_owr, code_sp, abundance) # Wide format
abund_wide_owr[is.na(abund_wide_owr)] <- 0
View(abund_wide_owr)

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



# AT OWR : Total infaunal biomass
biom_long_owr <- community_owr[,c(1:5,8,10)] # Long format
View(biom_long_owr)

biom_wide_owr <- spread(biom_long_owr, code_sp, biomass)
biom_wide_owr[is.na(biom_wide_owr)] <- 0
View(biom_wide_owr)

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



# AT OWR : Species richness (S)
richness_owr <- abund_long_owr %>%
  dplyr::group_by(corrosion, lineage, quadrat) %>%
  summarise(count = n_distinct(code_sp))
View(richness_owr)

richness_owr1 <- richness_owr %>%
  dplyr::group_by(lineage, corrosion) %>%
  summarise(
    count_mean = mean(count),
    count_sd = sd(count)
  )
richness_owr1



# AT OWR : Species diversity w/ Shannon diversity (H')
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



# AT OWR : Species diversity w/ Simpson diversity (λ)
simpson_owr <- diversity(abund_wide_owr[,6:75], "simpson")
abund_wide_owr1 <- abund_wide_owr1 %>%
  add_column(simpson = simpson_owr, .after = "shannon")

simpson_owr1 <- abund_wide_owr1 %>%
  dplyr::group_by(lineage, corrosion) %>%
  summarise(
    simpson_mean = mean(simpson),
    simpson_sd = sd(simpson)
  )
simpson_owr1



# AT OWR : Pielou's evenness (J)
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
pielou_owr1





##
# SECTION: Visualize data ----
##

# FIGURE 2 ====

description <- ggarrange(boxplot_bt2, barplot_abund_tot2, barplot_shannon2, barplot_pielou2,
                         boxplot_bt_owr1, barplot_abund_tot_owr1, barplot_shannon_owr1, barplot_pielou_owr1,
                         ncol = 4, nrow = 2, common.legend = F, legend = "none",
                         labels = c("A", "B", "C", "D", "E", "F", "G", "H"),
                         vjust = 1.15,
                         hjust = c(-1.5, -1.5, -1, -1, -1.5, -1.5, -1, -1),
                         font.label = list(size = 20, color = "grey30"))
description

annotate_figure(description, 
                left = text_grob("  Old Woman's River                                Across all sites       ",
                                 rot = 90, size = 18, face = "bold"))


# FIGURE 2a - ACROSS ALL SITES : Number of byssal threads ####
boxplot_bt2 <- ggplot(byssus1, aes(x = site, y = nb_byssal, fill = corrosion)) + 
  geom_boxplot() +
  scale_fill_manual(values = col_corr) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(limits = c(0, 450), expand = c(0,0)) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 430,
           label = c("AC", "BC", "A", "C"),
           size = 6) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  labs(fill = "Corrosion",
       x = "Site", y = "Average number of byssal threads",
       title = "Byssal threads")
boxplot_bt2



# FIGURE 2b - ACROSS ALL SITES : Total macrofaunal abundance ####
barplot_abund_tot2 <- ggplot(data = abund_tot1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = ab_tot_mean, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = ab_tot_mean - ab_tot_sd, ymax = ab_tot_mean + ab_tot_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 1100), expand = c(0,0)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = c(1050, 1050, 1050, 1050),
           label = c("A", "B", "B", "A"),
           size = 6) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average total macrofaunal abundance",
       title = "Total abundance")
barplot_abund_tot2



# FIGURE 2c - ACROSS ALL SITES : Shannon-Wiener index (H') ####
barplot_shannon2 <- ggplot(data = shannon1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = shannon_mean, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_y_continuous(limits = c(0, 3.15), breaks = c(0, 1, 2, 3), expand = c(0,0)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 3,
           label = c("A", "AB", "AB", "B"),
           size = 6) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average Shannon-Wiener diversity (H')",
       title = "Shannon's index (H')")
barplot_shannon2



# FIGURE 2d - ACROSS ALL SITES : Pielou's index (J) ####
barplot_pielou2 <- ggplot(data = pielou1, aes(x = site, fill = corrosion)) + 
  geom_bar(aes(y = pielou_mean, color = corrosion),
           stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_x_discrete(limits = site_order, labels = sites) +
  scale_y_continuous(limits = c(0, 1.05), breaks = c(0, 0.5, 1, 1), expand = c(0,0)) +
  annotate("text", 
           x = c(1, 2, 3, 4),
           y = 1,
           label = c("A", "A", "A", "B"),
           size = 6) +
  geom_vline(xintercept = c(1.5, 2.5, 3.5), linetype = "dotted") +
  scale_color_manual(values = col_corr) +
  scale_fill_manual(values = col_corr) +
  labs(color = "Corrosion", fill = "Corrosion",
       x = "Site", y = "Average Pielou's evenness (J)",
       title = "Pielou's index (J)")
barplot_pielou2



# FIGURE 2e - AT OWR : Number of byssal threads ####
ann_text_bt_owr1 <- data.frame(label = c("A", "B"), lineage = c("EAST", "WEST"), corrosion = "corroded")

boxplot_bt_owr1 <- ggplot(byssus_owr, aes(x = corrosion, y = nb_byssal, fill = corrosion)) + 
  geom_boxplot() +
  scale_fill_manual(values = col_corr) +
  scale_y_continuous(limits = c(0, 420), expand = c(0,0)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  geom_text(data = ann_text_bt_owr1, label = ann_text_bt_owr1$label, y = 400,
            size = 6, x = 1.5) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average number of byssal threads",
       title = "Byssal threads")
boxplot_bt_owr1



# FIGURE 2f - AT OWR : Total macrofaunal abundance  ####
barplot_abund_tot_owr1 <- ggplot(abund_tot_owr1, aes(x = corrosion, y = mean_abund, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean_abund - sd_abund, ymax = mean_abund + sd_abund), width = 0.2,
                position = position_dodge(0.9)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  scale_y_continuous(breaks = c(0, 50, 100, 150, 200, 250),limits = c(0, 270), expand = c(0,0)) +
  scale_x_discrete(labels = inf) +
  facet_wrap(~ lineage, nrow = 1) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average total infaunal abundance",
       title = "Total abundance")
barplot_abund_tot_owr1


# FIGURE 2g - AT OWR : Shannon-Wiener index (H') ####
barplot_shannon_owr1 <- ggplot(shannon_owr1, aes(x = corrosion, y = shannon_mean, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = shannon_mean - shannon_sd, ymax = shannon_mean + shannon_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(expand = c(0,0), limits = c(0,3.2)) +
  facet_wrap(~ lineage, nrow = 1) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average Shannon-Wiener diversity (H')",
       title = "Shannon's index (H')")
barplot_shannon_owr1



# FIGURE 2h - AT OWR : Pielou's index (J) ####
barplot_pielou_owr1 <- ggplot(pielou_owr1, aes(x = corrosion, y = pielou_mean, fill = corrosion)) + 
  geom_bar(aes(color = corrosion), stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = pielou_mean - pielou_sd, ymax = pielou_mean + pielou_sd), width = 0.2,
                position = position_dodge(0.9)) +
  scale_x_discrete(labels = inf) +
  scale_y_continuous(limits = c(0, 1.1), breaks = c(0, 0.5, 1, 1.5), expand = c(0,0)) +
  scale_fill_manual(values = col_corr) +
  scale_color_manual(values = col_corr) +
  facet_wrap(~ lineage, nrow = 1) +
  labs(fill = "Corrosion",
       x = "Corrosion", y = "Average Pielou's evenness (J)",
       title = "Pielou's index (J)")
barplot_pielou_owr1





# FIGURE 4 ====
contribution <- ggarrange(barplot_abund_cont, barplot_biom_cont,
                          ncol = 2, nrow = 1, common.legend = T, legend = "right",
                          labels = c("A", "B"), 
                          hjust = c(-0.5, -1.5),
                          font.label = list(size = 20, color = "grey30"))

annotate_figure(contribution, 
                left = text_grob("Across all sites",
                                 rot = 90, size = 22, face = "bold"))


View(cumsum_abund_long)

col_sp <- c("ACTINI" = "pink", "BURLAG" = "ivory4", "ISCHUT" = "mediumorchid3", "MYTGAL" = "dodgerblue2", "NUCDUB" = "palegreen",
            "PAREXI" = "salmon2", "PARPER" = "antiquewhite3", "PERPER" = "indianred3", "PSEPOD" = "gold1", "SPHPOL" = "cyan3",
            "PSEPOD" = "lightgray", "Others" = "white")
label_sp <- c("ACTINI" = expression(italic("Actiniaria")), "BURLAG" = expression(italic("Burnupena lagenaria")), 
              "ISCHUT" = expression(italic("Ischyromene huttoni")), "MYTGAL" = expression(italic("Mytilus galloprovincialis")),
              "NUCDUB" = expression(italic("Nucella dubia")), "PAREXI" = expression(italic("Parvulastra exigua")), 
              "PARPER" = expression(italic("Parisocladus perforatus")), "PERPER" = expression(italic("Perna perna")),
              "PSEPOD" = expression(italic("Pseudonereis podocirra")), "SPHPOL" = expression(italic("Sphaeramene polytylotos")), 
              "Others")
View(label_sp)

barplot_abund_cont <- ggplot(cumsum_abund_long, aes(x = site, y = contr, fill = factor(species, levels = c(setdiff(species, "Others"), "Others")))) +
  geom_col(position = "fill", col = "black") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_discrete(labels = c("BR vs JB", "BR vs MB", "BR vs PEd", "JB vs MB", "JB vs PEd", "MB vs PEd")) +
  scale_fill_manual(values = col_sp, labels = label_sp) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold", size = 20),
        legend.text = element_text(hjust = 0, size = 16),
        axis.text.x = element_text(size = 16),
        axis.title.y = element_text(size = 16)) +
  labs(fill = "Species",
       x = "Site by site comparison", y = "Species contribution to similarity",
       title = "Abundances")
barplot_abund_cont

barplot_biom_cont <- ggplot(cumsum_biom_long, aes(x = site, y = contr, fill = factor(species, levels = c(setdiff(species, "Others"), "Others")))) +
  geom_col(position = "fill", col = "black") +
  scale_y_continuous(labels = scales::percent, expand = c(0,0)) +
  scale_x_discrete(labels = c("BR vs JB", "BR vs MB", "BR vs PEd", "JB vs MB", "JB vs PEd", "MB vs PEd")) +
  scale_fill_manual(values = col_sp, labels = label_sp) +
  theme(legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(hjust = 0, size = 16),
        axis.text.x = element_text(size = 16)) +
  labs(fill = "Species",
       x = "Site by site comparison", y = "\n",
       title = "Biomasses")
barplot_biom_cont

