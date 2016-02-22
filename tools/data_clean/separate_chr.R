#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-20
# email: junjianglin@cs.toronto.com

#
# This r script is used to separate the microarry expression data 
#
# Input List:
# 1. Original expression data set 
#
# Output List:
# 1. separated microarray data in RData format

#######################   Initialization   ################
library(optparse)
rm(list=ls(all=TRUE))

option_list = list(
        make_option(c("-i","--input"),dest = "input",
                    help="the sample and feature data"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output prefix file")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$input) || is.null(opt$output)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}


#libraries to load
library(tools)

cat("\n","Initializing the environment and loading data...","\n")
input = opt$input
output = opt$output
if(file_ext(input) == 'txt'){
    dataSet = read.delim(file = input)
} else if(file_ext(input) == 'RData'){
    dataSet = get(load(input))
}

##################### define useful functions ######################
addAnnotation = function(feature_list){
        name_list = strsplit(feature_list,'[;-]')
        df_name = data.frame(do.call(rbind,name_list))
        annotation = data.frame(chr = df_name$X2, pos = as.numeric(as.character(df_name$X3)))
        return(annotation)
}

##################### adding annotation and separate ######################
cat("\n","start adding annotation and separate the data by chromosomes...","\n")
data = list(exprs = dataSet,anno=addAnnotation(rownames(dataSet)))
data_splitted = split(as.data.frame(data$exprs),as.factor(data$anno$chr))
for( i in 1:length(data_splitted)){
    name = names(data_splitted)[i]
    exprs = as.matrix(data_splitted[[i]])
    save(exprs,file=paste(output,"_",name,".RData",sep=""))
}





