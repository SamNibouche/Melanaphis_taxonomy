# Morphometric and Molecular Delineation of _Melanaphis sacchari_ (Zehntner, 1897) and _Melanaphis sorghi_ (Theobald, 1904)
## ***Nibouche S, Costet L, Medina RF, Holt JR, Sadeyen J, Zoogones AS, Moretti S, Halbert SE, Blackman RL***

### PLoS One 2020, DOI:x

This is code to replicate the analyses of the morphometric data reported in our manuscript 'Morphometric and Molecular Delineation of _Melanaphis sacchari_ and _Melanaphis sorghi_'

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
### DAPC
#### filtering sacchari and sorghi specimens, excluding MLL-B
```
filter<-(morpho2$species=="sorghi"|morpho2$species=="sacchari")&morpho2$MLL!="B"
morpho3<-morpho2[filter,]
```
#### elimination of specimens with NAs, except paratype NHM-1915-81
```
morphoDAPC<-rbind(na.omit(morpho3), morpho3[morpho3$specimen=="NHM-1915-81",])
```
#### data standardization
```
morphoDAPCstd<-data.frame(morphoDAPC[,c(1:3)],scale(morphoDAPC[,c(4:14)], center=TRUE, scale=TRUE))
```
#### replacing three missing values by zeros in paratype NHM-1915-81
```
morphoDAPCstd[morphoDAPCstd$specimen=="NHM-1915-81",is.na(morphoDAPCstd[morphoDAPCstd$specimen=="NHM-1915-81",])]<-0
```
#### sorting data by species
```
morphoDAPCstd<-morphoDAPCstd[order(morphoDAPCstd$species),]
```

#### DAPC
```
grp<-find.clusters(morphoDAPCstd[,4:14], max.n.clust = 4, stat="BIC",choose.n.clust = FALSE, 
                   method = "kmeans", criterion = "diffNgroup", n.pca =11)
dapc1<-dapc(morphoDAPCstd[,4:14],grp$grp,var.loadings=TRUE,n.pca=5,n.da=1)
```
#### DAPC using alternative interactive options
```
#grp<-find.clusters(morphoDAPCstd[,4:14], max.n.clust = 4, stat="BIC",choose.n.clust = TRUE,
#                   method = "kmeans", criterion = "diffNgroup")
# dapc1<-dapc(morphoDAPCstd[,4:14],grp$grp,var.loadings=TRUE,n.da=1)
```
#### Figure 3
```
loadingplot(dapc1$var.contr,srt=90,cex.lab=0.8,adj = -0.2,xlab="Morphological traits",main="")
```
#### Figure 1
```
compoplot(dapc1, posi="bottom",
          txt.leg=paste("Cluster", 1:2), lab="",
          xlab="specimens", col=c("blue","orange"))
```
#### Table 2 confusion matrix : MLL vs. morphology
```
assignations<-data.frame(morphoDAPCstd[,c(1:3)],dapc1$posterior)
assignations$assign.morpho[assignations$X1>=0.8]<-"sorghi"
assignations$assign.morpho[assignations$X2>=0.8]<-"sacchari"
assignations$assign.morpho[assignations$X1<0.8&assignations$X2<0.8]<-"undetermined"
table(assignations$MLL,assignations$assign.morpho)
```
#### Figure 2
```
data.graph<-morpho3
data.graph$species<-as.factor(as.vector(data.graph$species))
par(mfrow=c(2,2))
boxplot(data.graph[,"pt_cauda"]~data.graph$species,col=c("orange","blue"),xlab="",ylab="pt:cauda")
boxplot(data.graph[,"HindTibia_pt"]~data.graph$species,col=c("orange","blue"),xlab="",ylab="HindTibia:pt")
boxplot(data.graph[,"pt_siph"]~data.graph$species,col=c("orange","blue"),xlab="",ylab="pt:siph")
```
