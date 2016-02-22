#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-31
# email: junjianglin@cs.toronto.com

#
# This r script is used to select a range of probes within some threshold of a list of specific probes
# For example, the script can be used to select a neighbourhood of 5 probes around Restriction Enzyme site
#
# Input List:
# 1. Microarray data in RData or txt format
# 2. ProbeInfo data in RData format
# 3. range threshold(for example, 5 means neighbourhood 10 probes, upstream 5 and downstream 5)

#
# Output List:
# 1. Microarry data with probes in designated regions(RData and txt)
# 2. within range index

library(optparse)
rm(list=ls(all=TRUE))

option_list = list(
        make_option(c("-i","--input"),dest = "input",
                    help="the microarry data(ExpressionSet RData) or matrix in txt"),
        make_option(c("-p","--probe"),dest = "probeInfo",
                    help="the probeInfo RData contains information about probe"),
        make_option(c("-r","--range"),dest = "range",default=5,type="integer",
                    help="the range threshold for choosing probes,[default] = 5"),
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

getProbesWithinRERange = function(train.matrix = train.matrix, probeInfo = probeInfo, range = 5){
	re_idx = which(probeInfo$hasREsite == T)
    withinRange = numeric()
    for(i in re_idx){
        withinRange = c(withinRange,(i-range):(i+range))
    }
	withinRange = sort(unique(withinRange))
	return(withinRange)
}

cat("\n","Start calculating withinRange ...","\n")
withinRange = getProbesWithinRERange(dataSet,probeInfo,range)
dataMatrix = dataSet[withinRange,]
cat("\n","saving data","\n")
save(dataMatrix,file=paste(output,"_RE",range,".RData",sep=""))
write.table(withinRange,file=paste(output,"_RE",range,"_idx.txt",sep=""),quote=F,row.names=F,col.names=F)
write.table(dataMatrix,file=paste(output,"_RE",range,".txt",sep=""),quote=F,sep="\t")






