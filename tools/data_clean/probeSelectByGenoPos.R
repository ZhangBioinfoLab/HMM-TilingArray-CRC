#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-31
# email: junjianglin@cs.toronto.com

#
# This r script is used to select a range of probes within some threshold of a list of specific probes
# For example, the script can be used to select a neighbourhood of 200 bp around Restriction Enzyme site
#
# Input List:
# 1. Microarray data in RData or txt format
# 2. ProbeInfo data in RData format
# 3. range threshold(for example, 100 means neighbourhood 200bp, upstream 100 and downstream 100)

#
# Output List:
# 1. Microarry data with probes in designated regions(RData and txt)
# 2. within range idx

library(optparse)
rm(list=ls(all=TRUE))

option_list = list(
        make_option(c("-i","--input"),dest = "input",
                    help="the microarry data(ExpressionSet RData) or matrix in txt"),
        make_option(c("-p","--probe"),dest = "probeInfo",
                    help="the probeInfo RData contains information about probe"),
        make_option(c("-r","--range"),dest = "range",default=100,type="integer",
                    help="the range threshold for choosing probes,[default] = 100"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output file")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$input) || is.null(opt$output)|| is.null(opt$probeInfo)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}
              
#libraries to load
library(parallel)
library(tools)
library(Biobase)

cat("\n","Initializing the environment and loading data...","\n")
input = opt$input
probeInfo = opt$probeInfo
output = opt$output
range = opt$range


# load or read training data
if(file_ext(input) == 'txt'){
    dataSet = read.delim(file = input)
} else if(file_ext(input) == 'RData'){
    dataSet = exprs(get(load(input)))
}
probeInfo = get(load(probeInfo))

##################### define useful functions ######################

addAnnotation = function(raw.matrix){
        raw.dataframe = data.frame(raw.matrix)
        name_list = strsplit(rownames(raw.matrix),'[;-]')
        df_name = data.frame(do.call(rbind,name_list))
        annotation = data.frame(chr = df_name$X2,pos = as.numeric(as.character(df_name$X3)))
        return(annotation)
}

as.numeric.factor = function(x) {as.numeric(levels(x))[x]}
getProbesWithinRERange_helper = function(re_idx,re_start,train_pos,range){
        withinRange = numeric()
        k = re_idx - 1
        while(k > 0){
                if( (re_start - train_pos[k]) <= range){
                        withinRange = c(withinRange,k)
                        k = k-1
                }else{
                        break
                }
        }
        k = re_idx + 1
        while (TRUE){
                if((train_pos[k] - re_start) <= range){
                        withinRange = c(withinRange,k)
                        k = k+1
                }else{
                        break
                }
        }
        return(withinRange)
}

getProbesWithinRERange = function(train.matrix = train.matrix, probeInfo = probeInfo, range = 100){
	train_anno = addAnnotation(train.matrix)
	re_idx = which(probeInfo$hasREsite == T)
    re_starts = as.numeric.factor(probeInfo[probeInfo$hasREsite == T,]$Start)
    withinRange = numeric()
	withinRange = c(withinRange,re_idx)
    withinList = mcmapply(getProbesWithinRERange_helper,
                        re_idx,re_starts,MoreArgs = list(train_anno$pos,range))
    withinList = c(withinList,recursive=T)
    withinRange = c(withinRange,withinList)
	withinRange = sort(unique(withinRange))
	return(withinRange)
}

cat("\n","Start calculating withinRange ...","\n")
withinRange = getProbesWithinRERange(dataSet,probeInfo,range)
dataMatrix = dataSet[withinRange,]
cat("\n","saving data","\n")
save(dataMatrix,file=paste(output,"_RE",range,"bp.RData",sep=""))
write.table(withinRange,file=paste(output,"_RE",range,"bp_idx.txt",sep=""),quote=F,row.names=F,col.names=F)
write.table(dataMatrix,file=paste(output,"_RE",range,"bp.txt",sep=""),quote=F,sep="\t")



	
	
