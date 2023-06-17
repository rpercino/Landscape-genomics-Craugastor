#load library
library(ape)
library(pegas)
library(seqinr)
library(ggplot2)
library(ggmap)
library(rgdal)
library(colorspace)
library(adegenet)


## DAPC for denovo30 ###

#### denovo30_0.4_0.01_25X_64samp.recode.vcf ######
df<-read.vcf("/Volumes/LaCie/outfiles_vcf_oct19/denovo30_26X_0.5_0.01_.vcf.recode.vcf", to=1801) #64samples
pop<-read.csv("/Users/RutiPer/Documents/Craugastor_Analysis_Genomics_Ipyrad/data_denovo8/populations_vcf_denovo8.csv", sep="\t", header = TRUE)
pop<-pop[-c(1, 2, 3, 4, 5, 6, 7, 15, 24, 25, 26, 27, 31, 39, 46, 52, 54, 57, 75, 76, 77, 78, 87, 88, 89, 90, 91, 92, 93, 94), ]

df<-read.vcf("/Users/RutiPer/Documents/Craugastor_Analysis_Genomics_Ipyrad/data_denovo8/denovo8_1.recode.vcf", to=30816)
pop<-read.csv("/Users/RutiPer/Documents/Craugastor_Analysis_Genomics_Ipyrad/data_denovo8/populations_vcf_denovo8.csv", sep="\t", header = TRUE)

#
names <- as.vector(rownames(pop))
x <- df2genind(df, sep="/", NA.char="./.", ncode=NULL, ind.names=names, loc.names=NULL, pop=pop$pop, ploidy=2)

#identify the ideal number of clusters in the data
grp<-find.clusters(x, max.n.clust=10)

#perform the dapc
dapc1<-dapc(x, grp$grp)

#find the optimal number of PCs to retain by finding the optimal a-score
aScore <- optim.a.score(dapc1)

#OR: find the optimal number of PCs to retain by cross-validation
mat <- tab(x, NA.method="mean")
xval <- xvalDapc(mat, grp$grp, n.pca.max = 100, training.set = 0.9,
                 result = "groupMean", center = TRUE, scale = FALSE,
                 n.pca = NULL, n.rep = 30, xval.plot = TRUE)
xval[2:6]

#redo DAPC using the optimized number of PCs
dapc2 <- dapc(x, grp$grp, n.pca=10)

#Generate plots
myCol<- c("darkseagreen3","orange", "steelblue3")
myCol<- c("darkseagreen3", "steelblue3", "orange")
myCol<- c("yellow", "orange", "darkseagreen", "steelblue3")
scatter(dapc2, scree.da=TRUE, bg="white", pch=20, cell=3, 
        cstar=1, col=myCol, solid=1, cex=2,clab=0, leg=FALSE,
        txt.leg=paste("Cluster",1:4), posi.da="bottomright")

###
scatter(dapc1)
scatter(dapc2, posi.da="topleft", col=myCol)
scatter(dapc2, posi.da="bottomright", pch=19, clab=0, cstar=1, col=myCol)
assignplot(dapc2)
compoplot(dapc2, posi="bottomleft",
          txt.leg=paste("Cluster", 1:3), lab="",
          ncol=1, xlab="individuals", col=myCol)

dapc2
