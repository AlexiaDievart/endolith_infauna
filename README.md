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

![Fig_1_Material and methods_V4](https://github.com/user-attachments/assets/1580c0ca-154e-4e6d-8bcb-f1605b17f64d)

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
Within-bed architectural complexity did not differ significantly between infested corroded and non-infested corroded mussel beds at any site (Table S3). However, the total number of mussels, the number of live mussels and the number of byssal threads per mussel differed significantly among sites (Table S3a,b,e, Figure 2a). Across all sites (excluding Old Woman’s River), the total number of mussels and of live mussels was significantly lower at Port Edward than at Mossel Bay (p < 0.05; Table S4a,b).  Meanwhile, the average number of byssal threads per mussel was significantly higher at Jeffreys Bay than at both Brenton-on-Sea and Port Edward , and significantly higher at Mossel Bay than at Brenton-on-Sea (p  <  0.05 in all cases; Table S4da, Figure 2a). 

### EUENDOLITHS x PERNA LINEAGES at Old Woman's River

#### Within-bed architectural complexity
At Old Woman’s River, both the number of dead mussels and the average number of byssal threads per mussel differed significantly between P. perna lineages (Table S3h,j, Figure 2e). The number of dead mussels and the average number of byssal threads per mussel were significantly higher in western _P. perna_ mussel beds than eastern mussel beds (Table S4d,e, Figure 2e).
