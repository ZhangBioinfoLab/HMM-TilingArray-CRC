#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-05-4
# email: junjianglin@cs.toronto.com

#
# This r script is used to get labels info from the training and testing data with 
# sample names in format like "CRC038_CRC"
#
# Input List:
# 1.Data matrix after the training and testing partitioning (RData)
# 
#
# Output List:
# 1.The label data in txt format with output name specified 


#######################   Initialization   ################
library(optparse)
rm(list=ls(all=TRUE))


option_list = list(
        make_option(c("-i","--input"),dest = "input",
                    help="the matrix RData"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output file")

        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$input) || is.null(opt$output)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}


#libraries to load

cat("\n","Initializing the environment and loading data...","\n")
dataMatrix = get(load(opt$input))
output = opt$output
#######################   0. Preprocessing   ###############
cat("\n","preprocessing sample data...","\n")
sample_names = colnames(dataMatrix)
labels = do.call('rbind',strsplit(sample_names,"_"))[,2]

cat("\n","saving label data...","\n")
write.table(labels,file=paste(output,'_labels.txt',sep=''),quote=F,col.names=F,row.names=F)



























