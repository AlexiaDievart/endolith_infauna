# Repository description

R scripts and datasets to test the influence of euendolithic infestation on infaunal communities associated with Perna perna mussel beds along the South African coastline

## General information

**Citation**: 

**Principal investigator**: [Alexia M. A. DIEVART](https://scholar.google.com/citations?user=1CQgX5kAAAAJ&hl=fr&oi=ao), Coastal Research Group, Rhodes University (South Africa)

**Co-investigators**:
* [Christopher D. McQUAID](https://scholar.google.com/citations?user=uNl9g6wAAAAJ&hl=fr&oi=ao)
* [Gerardo I. ZARDI](https://scholar.google.com/citations?user=s8019k0AAAAJ&hl=fr&oi=ao)
* [Katy R. NICASTRO](https://scholar.google.com/citations?user=UUOXLPcAAAAJ&hl=fr&oi=ao)
* [Pierre W. FRONEMAN](https://scholar.google.com/citations?user=G5tEQu4AAAAJ&hl=fr&oi=ao)

**Corresponding investigator**: Alexia M. A. DIEVART, alexia.dievart@hotmail.fr

## Data collection and curation

### Geographic location

To assess the effects of euendolithic infestation on infaunal communities associated with both *Perna perna* lineages, artificial mussel beds were deployed on six transplant sites along the south and east coasts (Figure 1):
* Mosselbaai (34°10'58.6"S, 22°09'29.2"E) - western lineage
* Brenton-on-Sea (34°04'31.7"S 23°01'28.1"E) - western lineage
* Jeffreysbaai (34°01'33.2"S 24°55'50.4"E) - western lineage
* Old Woman’s River (33°28'55.7"S 27°09'08.2"E) - both lineages
* Port Edward (31°03'23.5"S, 30°13'40.6"E) - eastern lineage

![Figure_Map of transplant sites](https://user-images.githubusercontent.com/87645412/211604337-d043f6f5-e19e-4d0e-b850-91a66f685ae3.png)

The manipulative experiment deployed in the distribution area of pure genetic lineages (Barker, 2021; Cunha et al., 2014; Zardi et al., 2007b) consisted of two treatments: (a) 100% non-infested *Perna perna* individuals (category A-B), or (b) 100% infested *P. perna individuals* (category C-D), as defined in Kaehler (1999) (Figure 2). 
The common-garden experiment conducted within the overlapping area (i.e., Old Woman's River) consisted in four treatments: (a) 100% non-infested eastern *P. perna*, (b) 100% non-infested western *P. perna*, (c) 100% infested eastern *P. perna*, and (d) 100% infested western *P. perna* (Figure 2). 

![Figure_Experimental design](https://user-images.githubusercontent.com/87645412/211604855-0c65c130-0211-4916-88d1-6e5aab75f9a5.jpg)

## Data and file overview

### Data files

* **infauna_architecture.csv** includes the number of byssal threads for 3 mussels for each quadrat, each infestation level, each *Perna* lineage (only important for Old Woman's River) and each site.
* **infauna_description.csv** includes the number of live and dead mussels, the number of broken mussels and the average number of byssal threads for each quadrat, each infestation level, each *Perna* lineage (only important for Old Woman's River) and each site. 
* **infauna_community** includes the abundance (count) and biomass (in mg) for each infaunal species for each quadrat, each infestation level, each *Perna* lineage (only important for Old Woman's River) and each site, in a long format.


### R scripts

* **230109_infauna_analyses.R** includes the all the statistical analyses. Enjoy, it is the longest script I ever wrote and it is not the cleanest for sure.
* **231010_infauna_graphs.R** summarized all the graphs made to illustrate the chapter. 

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
