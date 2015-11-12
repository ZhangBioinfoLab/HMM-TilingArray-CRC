#-------------------INFO----------------------
# author: Junjiang Lin
# date: 2015-05-06
# email: junjianglin@cs.toronto.com

#
# This r script is used to detect outliers and quantitatively visualize them.


#######################   Initialization   ################

#libraries to load
library(cluster)


cat("\n","Initializing the environment and loading data...","\n")


###BATCH1###

#######################    Preprocessing   ###############
input = "batch1.RData"
output = "batch1"
dat = get(load(input))


## Calculating IACs for all pairs of samples and examining the distribution of IACs in the dataset:
IAC=cor(dat,method="pearson") 


# Draw IAC histogram
pdf(paste(output,"_IAC_hist.pdf",sep=""))
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
dev.off()

## Performing hierachical clustering (average linkage) using 1-IAC as a distance metric:
cluster1 = hclust(as.dist(1-IAC),method="average")
pdf(paste(output,"_IAC_cluster.pdf",sep=""))
plot(cluster1,cex=0.7,labels=dimnames(dat)[[2]])
dev.off()


## Another way to visualize outliers is to calculate the mean IAC for each array and examine this distribution:
meanIAC=apply(IAC,2,mean) 
sdCorr=sd(meanIAC) 
numbersd=(meanIAC-mean(meanIAC))/sdCorr 
pdf(paste(output,"_meanIAC.pdf",sep=""))
plot(numbersd)
dev.off()

pdf(paste(output,"_meanIAC_abline.pdf",sep=""))
plot(numbersd)
abline(h=-1.3)
dev.off()

##identify and remove outliers
sdout=-1.3
outliers=dimnames(dat)[[2]][numbersd<sdout]
dat2=dat[,numbersd>sdout] 


#names of outliers in first round: 
#"CRC148_CTR" "CRC076_CRC" "CRC157_CTR" "CRC180_CTR" "CRC161_CTR"
#"CRC136_CTR" "CRC187_CTR" "CRC078_CRC" "CRC164_CTR" "CRC163_CTR"
#"CRC114_CTR" "CRC059_CRC" "CRC173_CTR" "CRC108_CTR" "CRC075_CRC"

#######################    After FIRST round Preprocessing   ###############
IAC2=cor(dat2,method="pearson") 

pdf(paste(output,"_IAC_hist_after.pdf",sep=""))
hist(IAC2,sub=paste("Mean=",format(mean(IAC2[upper.tri(IAC2)]),digits=3)))
dev.off()

cluster2=hclust(as.dist(1-IAC2),method="average") 
pdf(paste(output,"_IAC_cluster_after.pdf",sep=""))
plot(cluster2,cex=0.7,labels=dimnames(dat2)[[2]])
dev.off()

meanIAC2=apply(IAC2,2,mean) 
sdCorr2=sd(meanIAC2) 
numbersd2=(meanIAC2-mean(meanIAC2))/sdCorr2 
pdf(paste(output,"_meanIAC_after.pdf",sep=""))
plot(numbersd2)
dev.off()

pdf(paste(output,"_meanIAC_abline_after.pdf",sep=""))
plot(numbersd2)
abline(h=-3)
dev.off()

sdout2=-3 
outliers2=dimnames(dat2)[[2]][numbersd2<sdout2] 
dat3=dat2[,numbersd2>sdout2] 

#names of outliers in second round:
#"CRC091_CRC" "CRC004_CRC" "CRC049_CRC"

#######################    After Second round Preprocessing   ###############
IAC3=cor(dat3,method="pearson") 

pdf(paste(output,"_IAC_hist_after2.pdf",sep=""))
hist(IAC3,sub=paste("Mean=",format(mean(IAC3[upper.tri(IAC3)]),digits=3)))
dev.off()

cluster3=hclust(as.dist(1-IAC3),method="average") 
pdf(paste(output,"_IAC_cluster_after2.pdf",sep=""))
plot(cluster3,cex=0.7,labels=dimnames(dat3)[[2]])
dev.off()

meanIAC3=apply(IAC3,2,mean) 
sdCorr3=sd(meanIAC3) 
numbersd3=(meanIAC3-mean(meanIAC3))/sdCorr3 
pdf(paste(output,"_meanIAC_after2.pdf",sep=""))
plot(numbersd3)
dev.off()

pdf(paste(output,"_meanIAC_abline_after2.pdf",sep=""))
plot(numbersd3)
abline(h=-3)
dev.off()

save(dat3,file="batch1_noOutliers.RData")
#Summaries:   18 outliers were removed and 175 samples were remained for further analysis
#all outliers:
#"CRC148_CTR" "CRC076_CRC" "CRC157_CTR" "CRC180_CTR" "CRC161_CTR"
#"CRC136_CTR" "CRC187_CTR" "CRC078_CRC" "CRC164_CTR" "CRC163_CTR"
#"CRC114_CTR" "CRC059_CRC" "CRC173_CTR" "CRC108_CTR" "CRC075_CRC"
#"CRC091_CRC" "CRC004_CRC" "CRC049_CRC"

###BATCH2###
#######################    Preprocessing   ###############
input = "batch2.RData"
output = "batch2"
dat = get(load(input))


## Calculating IACs for all pairs of samples and examining the distribution of IACs in the dataset:
IAC=cor(dat,method="pearson") 


# Draw IAC histogram
pdf(paste(output,"_IAC_hist.pdf",sep=""))
hist(IAC,sub=paste("Mean=",format(mean(IAC[upper.tri(IAC)]),digits=3)))
dev.off()

## Performing hierachical clustering (average linkage) using 1-IAC as a distance metric:
cluster1 = hclust(as.dist(1-IAC),method="average")
pdf(paste(output,"_IAC_cluster.pdf",sep=""))
plot(cluster1,cex=0.7,labels=dimnames(dat)[[2]])
dev.off()


## Another way to visualize outliers is to calculate the mean IAC for each array and examine this distribution:
meanIAC=apply(IAC,2,mean) 
sdCorr=sd(meanIAC) 
numbersd=(meanIAC-mean(meanIAC))/sdCorr 
pdf(paste(output,"_meanIAC.pdf",sep=""))
plot(numbersd)
dev.off()

pdf(paste(output,"_meanIAC_abline.pdf",sep=""))
plot(numbersd)
abline(h=-3.3)
dev.off()

##identify and remove outliers
sdout=-3.3
outliers=dimnames(dat)[[2]][numbersd<sdout]
dat2=dat[,numbersd>sdout] 


#names of outliers in first round: 
#"CRC298_CRC" "CRC340_CRC" "CRC224_CTR" "CRC301_CRC"

#######################    After FIRST round Preprocessing   ###############
IAC2=cor(dat2,method="pearson") 

pdf(paste(output,"_IAC_hist_after.pdf",sep=""))
hist(IAC2,sub=paste("Mean=",format(mean(IAC2[upper.tri(IAC2)]),digits=3)))
dev.off()

cluster2=hclust(as.dist(1-IAC2),method="average") 
pdf(paste(output,"_IAC_cluster_after.pdf",sep=""))
plot(cluster2,cex=0.7,labels=dimnames(dat2)[[2]])
dev.off()

meanIAC2=apply(IAC2,2,mean) 
sdCorr2=sd(meanIAC2) 
numbersd2=(meanIAC2-mean(meanIAC2))/sdCorr2 
pdf(paste(output,"_meanIAC_after.pdf",sep=""))
plot(numbersd2)
dev.off()

pdf(paste(output,"_meanIAC_abline_after.pdf",sep=""))
plot(numbersd2)
abline(h=-1.5)
dev.off()

sdout2=-1.5
outliers2=dimnames(dat2)[[2]][numbersd2<sdout2] 
dat3=dat2[,numbersd2>sdout2] 

#names of outliers in second round:
#"CRC345_CRC" "CRC241_CTR" "CRC204_CTR" "CRC258_CTR" "CRC306_CRC"
#"CRC309_CRC" "CRC331_CRC" "CRC305_CRC" "CRC281_CTR" "CRC324_CRC"
#"CRC270_CTR" "CRC293_CTR" "CRC261_CTR" "CRC358_CRC" "CRC310_CRC"
#"CRC308_CRC" "CRC299_CRC" "CRC307_CRC"

#######################    After SECOND round Preprocessing   ###############
IAC3=cor(dat3,method="pearson") 

pdf(paste(output,"_IAC_hist_after2.pdf",sep=""))
hist(IAC3,sub=paste("Mean=",format(mean(IAC3[upper.tri(IAC3)]),digits=3)))
dev.off()

cluster3=hclust(as.dist(1-IAC3),method="average") 
pdf(paste(output,"_IAC_cluster_after2.pdf",sep=""))
plot(cluster3,cex=0.7,labels=dimnames(dat3)[[2]])
dev.off()

meanIAC3=apply(IAC3,2,mean) 
sdCorr3=sd(meanIAC3) 
numbersd3=(meanIAC3-mean(meanIAC3))/sdCorr3 
pdf(paste(output,"_meanIAC_after2.pdf",sep=""))
plot(numbersd3)
dev.off()

pdf(paste(output,"_meanIAC_abline_after2.pdf",sep=""))
plot(numbersd3)
abline(h=-2.3)
dev.off()

sdout3=-2.3
outliers3=dimnames(dat3)[[2]][numbersd3<sdout3] 
dat4=dat3[,numbersd3>sdout3] 

#names of outliers in third round:
#"CRC223_CTR" "CRC363_CRC"

#######################    After THIRD round Preprocessing   ###############
IAC4=cor(dat4,method="pearson") 

pdf(paste(output,"_IAC_hist_after3.pdf",sep=""))
hist(IAC4,sub=paste("Mean=",format(mean(IAC4[upper.tri(IAC4)]),digits=3)))
dev.off()

cluster4=hclust(as.dist(1-IAC4),method="average") 
pdf(paste(output,"_IAC_cluster_after3.pdf",sep=""))
plot(cluster4,cex=0.7,labels=dimnames(dat4)[[2]])
dev.off()

meanIAC4=apply(IAC4,2,mean) 
sdCorr4=sd(meanIAC4) 
numbersd4=(meanIAC4-mean(meanIAC4))/sdCorr4
pdf(paste(output,"_meanIAC_after3.pdf",sep=""))
plot(numbersd4)
dev.off()

save(dat4,file="batch2_noOutliers.RData")
## 24 outliers were removed and 176 samples in batch2 were remained for further analysis
# all outliers
#"CRC298_CRC" "CRC340_CRC" "CRC224_CTR" "CRC301_CRC"
#"CRC345_CRC" "CRC241_CTR" "CRC204_CTR" "CRC258_CTR" "CRC306_CRC"
#"CRC309_CRC" "CRC331_CRC" "CRC305_CRC" "CRC281_CTR" "CRC324_CRC"
#"CRC270_CTR" "CRC293_CTR" "CRC261_CTR" "CRC358_CRC" "CRC310_CRC"
#"CRC308_CRC" "CRC299_CRC" "CRC307_CRC"
#"CRC223_CTR" "CRC363_CRC"

