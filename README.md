# Morphometric and Molecular Delineation of Melanaphis sacchari (Zehntner, 1897) and Melanaphis sorghi (Theobald, 1904)
## ***Nibouche S, Costet L, Medina RF, Holt JR, Sadeyen J, Zoogones AS, Moretti S, Halbert SE, Blackman RL

### PLoS One 2020, DOI:x

This is code to replicate the analyses of the morphometric data reported in our manuscript 'Morphometric and Molecular Delineation of Melanaphis sacchari and Melanaphis sorghi'

## Abstract

*Melanaphis sacchari* (Zehntner) and *Melanaphis sorghi* (Theobald) are major worldwide crop pests causing direct feeding damage on sorghum and transmitting viruses to sugarcane. In the USA, Mexico and Caribbean, an invasive genotype of sugarcane aphid has become a major pest of sorghum since 2013. The taxonomic status of *M. sacchari* varies among authors. Some consider that the species *M. sorghi* is a synonym of *M. sacchari*, despite no formal study demonstrating this synonymy. In this study, based on the comparison of samples collected from the whole distribution area, we use both morphometric and molecular data to confirm that *M. sacchari* and *M. sorghi* are distinct species. The invasive genotype recently introduced in the Americas is found to be *M. sorghi*.

## Data analysis
### Data
Morphometric data are available at [DOI:10.18167/DVN1/PDPDS4](http://dx.doi.org/10.18167/DVN1/PDPDS4)

### Load packages
```
library(ade4)
library(adegenet)
```
### dataset import and selection of 12 traits
```
morpho <- read.table ("morphometry_apterous_mean_and_ratio.csv", head = T, sep = ",", dec =".", fileEncoding = "latin1")
morpho2<-morpho[,c(1:3,8,26:32,34:36)]
```
