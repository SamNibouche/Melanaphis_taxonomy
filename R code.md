# Data analysis
### Data
Morphometric data are available at [DOI:10.18167/DVN1/PDPDS4](http://dx.doi.org/10.18167/DVN1/PDPDS4)
### Load packages
```
library(ade4)
library(adegenet)
library(MASS)
library(ggplot2)
```
### Import the dataset and select the 12 morphological traits
```
morpho <- read.table ("morphometry_apterous_raw.csv", head = T, sep = ",", dec =".", fileEncoding = "latin1")
```
### Computation of means of left and right sides, and ratio
```
morpho$pt<-rowMeans(morpho[,c("pt_r","pt_l")], na.rm=TRUE)
morpho$HindTibia<-rowMeans(morpho[,c("HindTibia_r","HindTibia_l")], na.rm=TRUE)
morpho$Ant<-rowMeans(morpho[,c("Ant_r","Ant_l")], na.rm=TRUE)
morpho$htII<-rowMeans(morpho[,c("htII_r","htII_l")], na.rm=TRUE)
morpho$VIb<-rowMeans(morpho[,c("VIb_r","VIb_l")], na.rm=TRUE)
morpho$siph<-rowMeans(morpho[,c("siph_r","siph_l")], na.rm=TRUE)
morpho$siphBW<-rowMeans(morpho[,c("siphBW_r","siphBW_l")], na.rm=TRUE)
morpho$pt_cauda<-morpho$pt/morpho$cauda
morpho$HindTibia_pt<-morpho$HindTibia/morpho$pt
morpho$Ant_BL<-morpho$Ant/morpho$BL
morpho$urs_htII<-morpho$urs/morpho$htII
morpho$pt_VIb<-morpho$pt/morpho$VIb
morpho$pt_siph<-morpho$pt/morpho$siph
morpho$cauda_urs<-morpho$cauda/morpho$urs
morpho$siph_BL<-morpho$siph/morpho$BL
morpho$siph_siphBW<-morpho$siph/morpho$siphBW
morpho$siph_cauda<-morpho$siph/morpho$cauda

morpho2<-morpho[,c(1:3,40,52:61)]
```
### DAPC
#### filter sacchari and sorghi specimens, excluding MLL-B
```
filter<-(morpho2$species=="sorghi"|morpho2$species=="sacchari")&morpho2$MLL!="B"
morpho3<-morpho2[filter,]
```
#### eliminate specimens with NAs, except Theobald's paratype NHM-1915-81
```
morphoDAPC<-rbind(na.omit(morpho3), morpho3[morpho3$specimen=="NHM-1915-81",])
```
#### standardize data 
```
morphoDAPCstd<-data.frame(morphoDAPC[,c(1:3)],scale(morphoDAPC[,c(4:14)], center=TRUE, scale=TRUE))
```
#### replace three missing values by zeros in paratype NHM-1915-81
```
morphoDAPCstd[morphoDAPCstd$specimen=="NHM-1915-81",is.na(morphoDAPCstd[morphoDAPCstd$specimen=="NHM-1915-81",])]<-0
```
#### sort data by species
```
morphoDAPCstd<-morphoDAPCstd[order(morphoDAPCstd$species),]
```
#### DAPC
```
grp<-find.clusters(morphoDAPCstd[,4:14], max.n.clust = 4, stat="BIC",choose.n.clust = FALSE, method = "kmeans", criterion = "diffNgroup", n.pca =11)
dapc1<-dapc(morphoDAPCstd[,4:14],grp$grp,var.loadings=TRUE,n.pca=5,n.da=1)

# DAPC using alternative interactive options
#grp<-find.clusters(morphoDAPCstd[,4:14], max.n.clust = 4, stat="BIC",choose.n.clust = TRUE, method = "kmeans", criterion = "diffNgroup")
# dapc1<-dapc(morphoDAPCstd[,4:14],grp$grp,var.loadings=TRUE,n.da=1)
```
#### Figures
```
# Figure 1
compoplot(dapc1, posi="bottom",
          txt.leg=paste("Cluster", 1:2), lab="",
          xlab="specimens", col=c("blue","orange"))

# Figure 2
data.graph<-morpho3
data.graph$species<-as.factor(as.vector(data.graph$species))
par(mfrow=c(2,2))
boxplot(data.graph[,"pt_cauda"]~data.graph$species,col=c("orange","blue"),xlab="",ylab="pt:cauda")
boxplot(data.graph[,"HindTibia_pt"]~data.graph$species,col=c("orange","blue"),xlab="",ylab="HindTibia:pt")
boxplot(data.graph[,"pt_siph"]~data.graph$species,col=c("orange","blue"),xlab="",ylab="pt:siph")

# Figure 3
loadingplot(dapc1$var.contr,srt=90,cex.lab=0.8,adj = -0.2,xlab="Morphological traits",main="")
```
#### Table 2 confusion matrix : MLL vs. morphology
```
assignations<-data.frame(morphoDAPCstd[,c(1:3)],dapc1$posterior)
assignations$assign.morpho[assignations$X1>=0.8]<-"sorghi"
assignations$assign.morpho[assignations$X2>=0.8]<-"sacchari"
assignations$assign.morpho[assignations$X1<0.8&assignations$X2<0.8]<-"undetermined"
table(assignations$MLL,assignations$assign.morpho)
```
### LDA
#### eliminate specimens with NAs
```
morphoLDA<-na.omit(morpho2)
```
#### standardize data 
```
morphoLDAstd<-data.frame(morphoLDA[,c(1:3)],scale(morphoLDA[,c(4:14)], center=TRUE, scale=TRUE))
```
#### train dataset without _M. sorghi_
```
train<-morphoLDAstd[morphoLDAstd$species!="sorghi",]
```
#### test dataset with _M.sorghi_ only
```
test<-morphoLDAstd[morphoLDAstd$species=="sorghi",]
```
#### LDA
```
cl<-as.factor(as.vector(train$species))
cl2<-as.factor(as.vector(test$species))

z<-lda(train[,c(4:14)],cl)
w<-predict(z,test[,c(4:14)])

res.lda<-as.matrix(train[,c(4:14)])%*%z$scaling
result<-rbind(res.lda,w$x)
colours<-as.factor(c(as.vector(cl),as.vector(cl2)))

col1<-as.vector(colours)
result2<-data.frame(result,col1)
result2$pred[result2$col1=="sorghi"]<-"test"
result2$pred[result2$col1!="sorghi"]<-"train"
```
#### Figure 4
```
p<-ggplot(result2, aes(x=LD1,y=LD2, colour=col1,shape=pred)) +
  geom_point(size=4) +
  stat_ellipse(level=0.68,type="norm")+
  scale_colour_manual(values = c("green","orange","blue","purple"), 
                      name="Species",
                      labels = c("M. indosacchari", "M. sacchari", "M. sorghi", "M. sorini")
                      )+
  scale_shape_manual(values = c(17,19),
                     name="")+
  theme(legend.position = "bottom",
        legend.background = element_blank(),
        legend.key = element_blank(),
        legend.text=element_text(face="italic", size=12),
        plot.background =  element_blank(), 
        panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), 
        panel.border = element_blank(), 
        panel.background = element_rect(fill = "white", colour = "black"),
        axis.line = element_line(size=.4),
        axis.line.y.left=element_line(size=.4)
  )
p
```

