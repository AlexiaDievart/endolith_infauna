# Repository description

Datasets and R scripts concerning the scientific publication _Biogeography, not intraspecific trait variation, determines macrofaunal communities associated with mussel beds_ by Dievart _et al._ (2024) [in revision].
  
## General information

**Citation**: 

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

### Geographic location

To assess the effects of euendolithic corrosion on infaunal communities associated with both *Perna perna* lineages, experimental mussel beds were deployed on six transplant sites along the south and east coasts of South Africa (Figure 1a):
* Mosselbaai (34°10'58.6"S, 22°09'29.2"E) - western lineage
* Brenton-on-Sea (34°04'31.7"S 23°01'28.1"E) - western lineage
* Jeffreysbaai (34°01'33.2"S 24°55'50.4"E) - western lineage
* Old Woman’s River (33°28'55.7"S 27°09'08.2"E) - both lineages
* Port Edward (31°03'23.5"S, 30°13'40.6"E) - eastern lineage
  

The first manipulative experiment deployed in the distribution area of pure genetic lineages (Barker, 2021; Cunha et al., 2014; Zardi et al., 2007b) consisted of two treatments: (a) 100% non-infested *Perna perna* individuals (category A-B), or (b) 100% infested *P. perna individuals* (category C-D), as defined in Kaehler (1999) (Figure 1b,d). 

The second manipulative (common-garden) experiment conducted within the overlapping area (i.e., Old Woman's River) consisted in four treatments: (a) 100% non-infested eastern *P. perna*, (b) 100% non-infested western *P. perna*, (c) 100% infested eastern *P. perna*, and (d) 100% infested western *P. perna* (Figure 1b,e). 

![Fig_1_Material and methods_V4](https://github.com/user-attachments/assets/1580c0ca-154e-4e6d-8bcb-f1605b17f64d)

## Data and file overview

### Data files

* **infauna_byssal.csv** includes the number of byssal threads for 3 randomly selected mussels for each quadrat, each corrosion level, each *Perna* lineage (only important for Old Woman's River) and each site. 
* **infauna_description.csv** includes the number of live and dead mussels, the number of broken mussels and the average number of byssal threads for each quadrat, each infestation level, each *Perna* lineage (only important for Old Woman's River) and each site. 
* **infauna_community** includes the abundance (count) and biomass (in mg) for each infaunal species for each quadrat, each infestation level, each *Perna* lineage (only important for Old Woman's River) and each site, in a long format.
* **infauna_infrared_control** (Annex S2) presents the shell temperatures of 5 randomly selected mussels for mesh-covered corroded and non-corroded experimental mussel beds (n = 3) recorded every 5 minutes for 90 minutes on 2 different dates. 


### R scripts

* **230109_infauna_analyses.R** includes the all the statistical analyses. Enjoy, it is the longest script I ever wrote and it is not the cleanest for sure.
* **231010_infauna_graphs.R** summarized all the graphs made to illustrate the chapter.
* **240920_infauna_infrared_control.R** includes the statistical analyses presented in Annex S2, where we tried to determine if the mesh that secured the experimental mussel beds to the rocky shores could have been an experimental design flaw that masked the beneficial effects of euendolithic corrosion. It was not ! 

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
