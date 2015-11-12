#-------------------INFO----------------------
# author: Junjiang Lin
# date: 2015-05-06
# email: junjianglin@cs.toronto.com

#
# This r script is used to get tiling array raw data with sample name and group info
#

#######################   Initialization   ################
rm(list=ls(all=TRUE))
library(AffyTiling)
library(limma)


batch1 <- "/home/zhanglab1/zl/No_Backup/Data/Petronis_Colorectal/CRC-1st Batch"

batch2 <- "/home/zhanglab1/zl/No_Backup/Data/Petronis_Colorectal/CRC-2nd Batch"

bpmap <- "/home/zhanglab1/yueli/Projects/colon_cancer_final/bpmap/Hs35b_P04R_v01-3_NCBIv36.bpmap"

batch1CEL <- grep("leftOutData", list.files(batch1, "CEL$", recursive=TRUE, full.names=TRUE), invert=TRUE, value=TRUE)

batch2CEL <- list.files(batch2, "CEL$", recursive=TRUE, full.names=TRUE)


######### get sample keys
batch1_sampleinfo <- read.delim("/home/zhanglab1/yueli/Projects/colon_cancer_final/data/batch1_sample_keys.txt", check.names=F)

batch2_sampleinfo <- read.delim("/home/zhanglab1/yueli/Projects/colon_cancer_final/data/batch2_sample_keys.txt", check.names=F)


sampleinfo_b1 <- batch1_sampleinfo[,c("FileName","Sample","Group")]
sampleinfo_b2 <- batch2_sampleinfo[,c("FileName","Sample","Group")]

sampleinfo_b1$Group <- sub("HEA", "CTR", sampleinfo_b1$Group)
sampleinfo_b2$Group <- sub("HEA", "CTR", sampleinfo_b2$Group)

sampleinfo_b1$celpaths <- batch1CEL[match(sampleinfo_b1$FileName, basename(batch1CEL))]
sampleinfo_b2$celpaths <- batch2CEL[match(sampleinfo_b2$FileName, basename(batch2CEL))]


readCEL <- function(x) {
    
        einter <- AnalyzeTilingCelFiles(CEL_filenames = x,
            BPMAP_filename = bpmap, readOnlyNCBI=FALSE, ReturnRawIntensities=TRUE)
    
        df <- matrix( as.numeric( einter[,-c(1:3)] ), ncol= (ncol(einter) - 3) )
    
        rownames(df) <- einter[,1]
    
        df  
}

data_b1 <- do.call("cbind",lapply(sampleinfo_b1$celpaths, readCEL))
colnames(data_b1) <- sprintf("%s_%s", sampleinfo_b1$Sample, sampleinfo_b1$Group)

message(sprintf("Finished processing batch1 sets!"))

data_b2 <- do.call("cbind",lapply(sampleinfo_b2$celpaths, readCEL))
colnames(data_b2) <- sprintf("%s_%s", sampleinfo_b2$Sample, sampleinfo_b2$Group)
message(sprintf("Finished processing batch2 sets!"))


# probeInfo from match_REprobes.R
load("/home/zhanglab1/yueli/Projects/colon_cancer_final/bpmap/probeInfo.RData")


# use only NCBI probes
data_b1 <- data_b1[rownames(data_b1) %in% probeInfo$UniqueID,]

data_b2 <- data_b2[rownames(data_b2) %in% probeInfo$UniqueID,]

save(data_b1,file="/home/zhanglab/junjianglin/project_data/cir_cancer/thesis_project/rawData/batch1.RData")
save(data_b2,file="/home/zhanglab/junjianglin/project_data/cir_cancer/thesis_project/rawData/batch2.RData")







