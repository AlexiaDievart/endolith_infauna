# Repository description

Datasets and R scripts concerning the scientific publication _Biogeography, not intraspecific trait variation, determines macrofaunal communities associated with mussel beds_ by Dievart _et al._ (2024) [in revision].
  
## General information

**Citation**: The correct citation will be added once the paper will be accepted for publication.

**Principal investigator**: [Alexia M. A. DIEVART](https://scholar.google.com/citations?user=1CQgX5kAAAAJ&hl=fr&oi=ao), Coastal Research Group, Rhodes University (South Africa)

**Co-investigators**:
* [Christopher D. McQUAID](https://scholar.google.com/citations?user=uNl9g6wAAAAJ&hl=fr&oi=ao)
* [Gerardo I. ZARDI](https://scholar.google.com/citations?user=s8019k0AAAAJ&hl=fr&oi=ao)
* [Katy R. NICASTRO](https://scholar.google.com/citations?user=UUOXLPcAAAAJ&hl=fr&oi=ao)
* [Pierre W. FRONEMAN](https://scholar.google.com/citations?user=G5tEQu4AAAAJ&hl=fr&oi=ao)

**Corresponding investigator**: Alexia M. A. DIEVART, alexia.dievart@hotmail.fr

## Abstract
The abstract will be added once the paper will be accepted for publication.

## Data collection and curation

### Goals
The goals of this main investigation were to :
* Determine if euendolithic corrosion and its indirect effects on mussels would influence the within-bed architectural complexity (i.e., number of byssal threads per mussel) across bioregions or between _P. perna_ genetic lineages
* Determine if euendolithic corrosion would influence macrofaunal communities associated with mussel beds across bioregions, as euendolithic corrosion has been demonstrated to improve mussel bed microclimate and alleviate temperature and desiccation stress to associated macrofauna.
* Determine if euendolithic corrosion interacts with _P. perna_ genetic lineages to influence macrofaunal communities associated with mussel beds within the contact zone. 

### Material & Methods

To assess the effects of euendolithic corrosion on infaunal communities associated with both *Perna perna* lineages, experimental mussel beds were deployed on six transplant sites along the south and east coasts of South Africa (Figure 1a):
* Mosselbaai (34°10'58.6"S, 22°09'29.2"E) - western lineage
* Brenton-on-Sea (34°04'31.7"S 23°01'28.1"E) - western lineage
* Jeffreysbaai (34°01'33.2"S 24°55'50.4"E) - western lineage
* Old Woman’s River (33°28'55.7"S 27°09'08.2"E) - both lineages
* Port Edward (31°03'23.5"S, 30°13'40.6"E) - eastern lineage
  

The first manipulative experiment deployed in the distribution area of pure genetic lineages (Barker, 2021; Cunha et al., 2014; Zardi et al., 2007b) consisted of two treatments: (a) 100% non-corroded *Perna perna* individuals (category A-B), or (b) 100% corroded *P. perna individuals* (category C-D), as defined in Kaehler (1999) (Figure 1b,d). 

The second manipulative (common-garden) experiment conducted within the overlapping area (i.e., Old Woman's River) consisted in four treatments: (a) 100% non-corroded eastern *P. perna*, (b) 100% non-corroded western *P. perna*, (c) 100% corroded eastern *P. perna*, and (d) 100% corroded western *P. perna* (Figure 1b,e). 

![Fig_1_Material and methods_V5](https://github.com/user-attachments/assets/d9cc365d-4c39-486b-b126-edf20ba1df60)

An additional experiment was conducted in a laboratory setting to test if the presence of the mesh on the experimental mussel beds we deployed on the rocky shores could have masked the beneficial effects of euendolithic corrosion by increasing shading of the mussel shell, and consisted in two treatments: (a) 100% mesh-covered non-corroded _P. perna_ or (b) 100% mesh-covered corroded _P. perna_ (Annex S2).

## Data and file overview

### Data files

* **infauna_byssal.csv** includes the exact number of byssal threads for 3 randomly selected mussels for each quadrat, each corrosion level, each *Perna* lineage (only used in data analyses for Old Woman's River) and each site (including the common garden experiment at Old Woman's River). 
* **infauna_architectural complexity.csv** includes the number of live and dead mussels - with details of the mussel species in question - and the number of shell fragments, as well as the average number of byssal threads for each quadrat (calculated from 'infauna_byssal.csv'), each corrosion level, each *Perna* lineage (only used in data analyses for Old Woman's River) and each site (including the common garden experiment at Old Woman's River). 
* **infauna_community.csv** includes the abundance (count) and biomass (in mg) for each infaunal species for each quadrat, each corrosion level, each *Perna* lineage (only used in data analyses for Old Woman's River) and each site (including the common garden experiment at Old Woman's River) in a long format.
* **infauna_infrared_control.csv** (Annex S2) presents the shell temperatures of 5 randomly selected mussels for mesh-covered corroded and non-corroded experimental mussel beds (n = 3) recorded every 5 minutes for 90 minutes on 2 different dates. This was an additional experiment to ensure our experimental set up did not influence our results.


### R scripts

* **240930_infauna_architectural complexity.R** includes the statistical analyses performed to determine (a) the effects of euendolithic corrosion and biogeography on within-bed architectural complexity (e.g., total number of mussels, number of live mussels, average number of byssal threads per mussel) across all sites (except Old Woman's River), and (2) the effects of euendolithic corrosion and _P. perna_ lineages on within-bed architectural complexity at Old Woman's River.
* **241007_infauna_community descriptors.R** includes the statistical analyses performed to determine (a) the effects of euendolithic corrosion and biogeography on general macrofaunal community descriptors (e.g., total abundance, total biomass, species richness, etc.) across all sites (except Old Woman's River), and (2) the effects of euendolithic corrosion and _P. perna_ lineages on general macrofaunal community descriptors at Old Woman's River.
* **240430_infauna_graphs_descriptors+sp.per.sp** (Figures 2 and 4) includes the graphical visualization of the analyses applied on macrofaunal community descriptors (e.g., average number of byssal threads per mussel, total abundance, Shannon's diversity index, etc.) and on the dissimilarity in macrofaunal communities between sites (i.e., on abundances and biomasses). 
* **241021_infauna_graphs_nMDS.R** (Figure 3) includes the graphical visualization of the results of the multivariate analyses applied on macrofaunal communities associated with mussel beds.
* **240920_infauna_infrared_control.R** (Annex S2) includes the statistical analyses presented in Annex S2, where we tried to determine if the mesh that secured the experimental mussel beds to the rocky shores could have been an experimental design flaw that masked the beneficial effects of euendolithic corrosion. It was not ! 

## Expectations from the statistical analyses

### Across all sites (excluding Old Woman's River)

* Compare the architectural complexity between infested and non-infested mussel beds in terms of the number of live and dead mussels, the number of broken mussels and the average number of byssal threads for each quadrat
* Compare the "community descriptors" between infested and non-infested mussel beds in terms of average total abundance and biomass per quadrat, species richness, species diversity (Shannon-Wiener and Simpson's indexes) and species evenness
* Compare the infaunal communities between infested and non-infested mussel beds in terms of abundance and biomass
* Identify the species that contribute the most to dissimilarities between communities if any

### At Old Woman's River

* Compare the architectural complexity between infested and non-infested mussel beds AND eastern and western *Perna* lineages in terms of the number of live and dead mussels (*Perna* and *Mytilus* separately), the number of broken mussels and the average number of byssal threads for each quadrat
* Compare the "community descriptors" between infested and non-infested mussel beds AND eastern and western *Perna* lineages in terms of average total abundance and biomass per quadrat, species richness, species diversity (Shannon-Wiener and Simpson's indexes) and species evenness
* Compare the infaunal communities between infested and non-infested mussel beds AND eastern and western *Perna* lineages in terms of abundance and biomass
* Identify the species that contribute the most to dissimilarities between communities if any

### Annex S2 : Influence of the experimental set up on the thermal mitigation offered by photoautotrophic euendoliths between corroded and non-corroded mussel beds
* Compare the shell temperatures of mesh-covered corroded and non-corroded mussels in a setting that replicated our experimental design actually deployed in the field.

## Summary of the results

### EUENDOLITHS x BIOGEOGRAPHY across all sites (excluding Old Woman's River)

#### Within-bed architectural complexity
Within-bed architectural complexity did not differ significantly between corroded and non-corroded mussel beds at any site (Table S3). However, the total number of mussels (F(3,23) = 4.039, p = 0.019; Table S3a), the number of live mussels (F(3,23) = 4.333, p = 0.015; Table S3b) and the number of byssal threads per mussel (F(3,73) = 9.034, p < 0.0001; Table S3e) differed significantly among sites (Table S3a,b,e, Figure 2a). Across all sites (excluding Old Woman’s River), the total number of mussels and of live mussels was significantly lower at Port Edward than at Mossel Bay (p < 0.05; Table S4a,b). Meanwhile, the average number of byssal threads per mussel was significantly higher at Jeffreys Bay than at both Brenton-on-Sea and Port Edward, and significantly higher at Mossel Bay than at Brenton-on-Sea (p < 0.05 in all cases; Table S4d, Figure 2a). 

#### Macrofaunal communities
A total of 115 invertebrate taxa were recorded across all treatments and sites (Table S5), including 18 taxa identified to phylum, order, class, or family levels, and seven taxa identified to the genus level. Across all sites (excluding Old Woman’s River), 80 taxa were associated with euendolith-corroded mussel beds and 80 taxa with the non-corroded mussel beds, with 60 taxa in common. 

Across all sites (excluding Old Woman’s River), total macrofaunal abundance (X2(3,28) = 72.276, p < 0.0001; Table S6a) and biomass (X2(3,28) = 35.414, p < 0.0001; Table S6b), as well as species diversity, calculated using Shannon-Wiener index (F(3,23) = 3.905, p = 0.022; Table S8b) or Simpson’s index (F(3,23) = 4.351, p = 0.014; Table S8c) and species evenness (F(3,23) = 7.442, p = 0.001; Table S8d), differed significantly among sites, but not between levels of euendolithic corrosion (Tables S6 and S8b,c,d, Figure 2b,c,d). Species richness did not differ significantly among sites or between corrosion levels (Table S8a). At Mossel Bay and Port Edward, total macrofaunal abundance and biomass were significantly lower than at Brenton-on-Sea and Jeffreys Bay (Table S7a, Figure 2b), and at Brenton-on-Sea (Table S7b) respectively. Conversely, both species diversity indexes, Shannon-Wiener (H’) and Simpson’s (λ), were significantly higher at Port Edward than at Mossel Bay (Table S9a,b, Figure 2c), while Pielou’s evenness was significantly higher at Port Edward than all other sites (Table S9c, Figure 2d).

![Fig_2_Structural + Community descriptors_V4](https://github.com/user-attachments/assets/7bf884a2-a31a-4691-9dba-9f18e70b9f3d)

Macrofaunal abundances (F(3,27) = 11.345, p = 0.001; Table S11a) and biomasses (F(3,27) = 12.685, p = 0.001; Table S11b) significantly differed among sites, Old Woman’s River excluded (PERMANOVA, Table S11a,b, Figure 3a,b), but not between levels of euendolithic corrosion. Although significant differences were found in multivariate dispersion between sites for both species-specific macrofaunal abundance (F(3,24) = 8.162, p < 0.001) and biomass (F(3,24) = 7.001, p = 0.002), the nMDS plots showed clear group separations by sites (Figure 3a,b), indicating that the significant effect detected was at least partially, a result of differences among centroids and not only due to a difference in the spread of the macrofaunal data between groups. In addition, hierarchical clustering analyses mirrored the PERMANOVA results (see R scripts).

![Fig_3_nMDS and hierarchical clustering_V5](https://github.com/user-attachments/assets/018af990-49be-48c7-94a2-32d8d1a1c18b)

SIMPER analyses showed that the number of species collectively accounting for ≈ 50% of the dissimilarity among sites was higher for abundance than for biomass (Tables S13, S14). Amongst others, the following species collectively accounted for 50% of the dissimilarity between sites (excluding Old Woman’s River) (Tables S13, S14, Figure 4): the isopods Parisocladus perforatus (H. Milne Edwards, 1840), Sphaeramene polytylotos Barnard, 1914 and Ischyromene huttoni (G. Thomson, 1879), the juvenile bivalves Perna perna and Mytilus galloprovincialis Lamarck, 1819, the unidentified Actiniaria species (probably juvenile Bunodactis reynaudi (Milne Edwards, 1857)), the whelks Burnupena lagenaria (Lamarck, 1822) and Nucella dubia (Krauss, 1848), the starfish Parvulastra exigua (Lamarck, 1816), and the worm Pseudonereis podocirra (Schmarda, 1861).

![Fig_4_Contribution of species to abund and biom_V5](https://github.com/user-attachments/assets/02725ec2-0180-4cae-9d76-6c30679c0f17)

### EUENDOLITHS x PERNA LINEAGES at Old Woman's River

#### Within-bed architectural complexity
At Old Woman’s River, both the number of dead mussels (F(1,13) = 7.110, p = 0.019; Table S3h) and the average number of byssal threads per mussel (F(1,42) = 7.669, p = 0.008; Table S3j) differed significantly between P. perna lineages (Table S3h,j, Figure 2e). The number of dead mussels and the average number of byssal threads per mussel were significantly higher in western P. perna mussel beds than eastern mussel beds (p < 0.05; Table S4d,e, Figure 2e).

#### Macrofaunal communities
At Old Woman’s River, 40 taxa were associated with eastern corroded mussel beds, 50 with western corroded mussel beds, 53 with eastern non-corroded mussel beds, and 43 with western non-corroded mussel beds, with 25 taxa in common.

At Old Woman’s River, total macrofaunal abundance and biomass, species richness, species diversity, and species evenness, were not significantly different between corrosion levels or P. perna lineages (Table S10, Figure 2f,g,h).

There were no significant differences in macrofaunal communities between eastern and western P. perna lineages, nor between levels of euendolithic corrosion, at Old Woman’s River (PERMANOVA, Table S11, Figure 3).
