#-------------------INFO----------------------
# author: Junjiang Lin
# date: 2015-05-12
# email: junjianglin@cs.toronto.com

#
# This r script is used to merge two batchs of data and generate a number of ramdom sampling split

# Input List:
# 
# 1.A start seed
# 2.Primary output name
# 3.Percentage of training data
# 4.the number of random sampling
#
# Output List:
# 1.Randomly splitted train and test data

#######################   Initialization   ################

library(optparse)
rm(list=ls(all=TRUE))

option_list = list(
        make_option(c("-o","--output"),dest = "output",
                    help="the primary path of output file"),
        make_option(c("-s","--seed"),dest = "seed",default=1,
        			help="seed, [DEFAULT]=1"),
        make_option(c("-p","--trainP"),dest = "trainP",default=0.7,
        			help="the percentage of training data,[DEFAULT]=0.7"),
        make_option(c("-d","--data"),dest = "dataOption",default="both",
        			help="the data you want to load, batch1,batch2,or both, [Default]=both"),
        make_option(c("-n","--noRMA"),dest = "noRMA",action="store_true",default=F,
        			help="whether to save noRMA data as well, [DEFAULT]= False")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$output)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}


library(AffyTiling)

cat("\n","Initializing the environment and loading data...","\n")
output = opt$output
seed = opt$seed
trainP = opt$trainP
noRMA = opt$noRMA
dataOption = opt$dataOption
#######################	  0.Preprocessing and Merge   ########################
batch1 <- "/home/zhanglab1/zl/No_Backup/Data/Petronis_Colorectal/CRC-1st Batch"

batch2 <- "/home/zhanglab1/zl/No_Backup/Data/Petronis_Colorectal/CRC-2nd Batch"

bpmap <- "/home/zhanglab1/yueli/Projects/colon_cancer_final/bpmap/Hs35b_P04R_v01-3_NCBIv36.bpmap"

batch1CEL <- grep("leftOutData", list.files(batch1, "CEL$", recursive=TRUE, full.names=TRUE), invert=TRUE, value=TRUE)

batch2CEL <- list.files(batch2, "CEL$", recursive=TRUE, full.names=TRUE)




batch1_sampleinfo <- read.delim("/home/zhanglab1/yueli/Projects/colon_cancer_final/data/batch1_sample_keys.txt", check.names=F)
outliers1 = c("CRC148","CRC076","CRC157","CRC180","CRC161","CRC136","CRC187","CRC078","CRC164","CRC163",
		      "CRC114","CRC059","CRC173","CRC108","CRC075","CRC091","CRC004","CRC049")
outlier_idx1 = match(outliers1,batch1_sampleinfo$Sample)
batch1_sampleinfo = batch1_sampleinfo[-outlier_idx1,]

batch2_sampleinfo <- read.delim("/home/zhanglab1/yueli/Projects/colon_cancer_final/data/batch2_sample_keys.txt", check.names=F)
outliers2 = c("CRC298","CRC340","CRC224","CRC301","CRC345","CRC241","CRC204","CRC258","CRC306","CRC309","CRC331","CRC305", 
	         "CRC281","CRC324","CRC270","CRC293","CRC261","CRC358","CRC310","CRC308","CRC299","CRC307","CRC223","CRC363")
outlier_idx2 = match(outliers2,batch2_sampleinfo$Sample)
batch2_sampleinfo = batch2_sampleinfo[-outlier_idx2,]

if(dataOption == 'both'){
	celfiles.all <- c(batch1CEL, batch2CEL)
	sampleinfo <- rbind(batch1_sampleinfo[,c("FileName","Sample","Group")], batch2_sampleinfo[,c("FileName","Sample","Group")])
}else if(dataOption == 'batch1'){
	celfiles.all <- c(batch1CEL)
	sampleinfo <- batch1_sampleinfo[,c("FileName","Sample","Group")]
}else if(dataOption == 'batch2'){
	celfiles.all <- c(batch2CEL)
	sampleinfo <- batch2_sampleinfo[,c("FileName","Sample","Group")]
}else{
	message("No such data option")
	q(save = "no")
}

sampleinfo$Group <- sub("HEA", "CTR", sampleinfo$Group)

sampleinfo$celpaths <- celfiles.all[match(sampleinfo$FileName, basename(celfiles.all))]

#######################    1.random partition   ###############

readCEL_noRMA <- function(x) {
	
		einter <- AnalyzeTilingCelFiles(CEL_filenames = x,
			BPMAP_filename = bpmap,ReturnRawIntensities=TRUE)
	
		df <- matrix( as.numeric( einter[,-c(1:3)] ), ncol= (ncol(einter) - 3) )
	
		rownames(df) <- einter[,1]
	
		df	
}

readCEL_RMA <- function(x) {
	
		einter <- AnalyzeTilingCelFiles(CEL_filenames = x,
			BPMAP_filename = bpmap,)
	
		df <- matrix( as.numeric( einter[,-c(1:3)] ), ncol= (ncol(einter) - 3) )
	
		rownames(df) <- einter[,1]
	
		df	
}

cases.all <- subset(sampleinfo, Group=="CRC")
ctrls.all <- subset(sampleinfo, Group=="CTR")


set.seed(seed); cases.1 <- sample(1:nrow(cases.all), round(nrow(cases.all) * trainP))
set.seed(seed); ctrls.1 <- sample(1:nrow(ctrls.all), round(nrow(ctrls.all) * trainP))
cases.2 <- c(1:nrow(cases.all))[!c(1:nrow(cases.all)) %in% cases.1]
ctrls.2 <- c(1:nrow(ctrls.all))[!c(1:nrow(ctrls.all)) %in% ctrls.1]

trainset <- rbind(cases.all[cases.1,], ctrls.all[ctrls.1,])
testset <- rbind(cases.all[cases.2,], ctrls.all[ctrls.2,])

train.matrix <- readCEL_RMA(trainset$celpaths)

colnames(train.matrix) <- sprintf("%s_%s", trainset$Sample, trainset$Group)


message(paste("Finished processing training sets ",seed,sep=""))


test.matrix <- readCEL_RMA(testset$celpaths)

colnames(test.matrix) <- sprintf("%s_%s", testset$Sample, testset$Group)

message(paste("Finished processing training sets ",seed,sep=""))


save(train.matrix,file=paste(output,"/train",seed,".RData",sep=""))
save(test.matrix,file=paste(output,"/test",seed,".RData",sep=""))


# Save the last dataset without RMA for comparison
if(noRMA){
	train.matrix <- readCEL_noRMA(trainset$celpaths)

	colnames(train.matrix) <- sprintf("%s_%s", trainset$Sample, trainset$Group)


	message(paste("Finished processing training sets ",seed,sep=""))


	test.matrix <- readCEL_noRMA(testset$celpaths)

	colnames(test.matrix) <- sprintf("%s_%s", testset$Sample, testset$Group)

	message(paste("Finished processing training sets ",seed,sep=""))

	save(train.matrix,file=paste(output,"/train",seed,"_noRMA.RData",sep=""))
	save(test.matrix,file=paste(output,"/test",seed,"._noRMA.RData",sep=""))	
}














