#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-20
# email: junjianglin@cs.toronto.com

#
# This r script is used to make a training and testing set based on given samples using package caret
#
# Input List:
# 1. Original all samples data set with features but not labels(RData or txt)
# 2. Sample dimension, 0 means samples are in row, 1 means samples are in column
# 3. Labels of the sample in txt format
# 4. percentage of training data
#
# Output List:
# 1. training data in RData and Txt
# 2. training data labels in Txt
# 3. testing data  in RData and Txt
# 4. testing data labels in Txt

#######################   Initialization   ################
library(optparse)
rm(list=ls(all=TRUE))
set.seed(1)

option_list = list(
        make_option(c("-i","--input"),dest = "input",
                    help="the sample and feature data"),
        make_option(c("-d","--dim"),dest = "dim",type="integer",default=1,
        			help="0 means samples are in row, 1 means samples are in column,\n[default] = 1, which means samples are in the column"),
        make_option(c("-l","--label"),dest = "label",
        			help="label file"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output prefix file"),
        make_option(c("-p","--percentage"),dest = "percentage",type="double",default=0.7,
                    help="percentage of training data,\n[default] = 0.7")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$input) || is.null(opt$label) || is.null(opt$output)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}


#libraries to load
library(Biobase)
library(caret)
library(tools)

cat("\n","Initializing the environment and loading data...","\n")
input = opt$input
output = opt$output
dim = opt$dim
labels = read.table(opt$label)[,1]
p = opt$percentage
if(file_ext(input) == 'txt'){
    dataFrame = read.delim(file = input)
} else if(file_ext(input) == 'RData'){
    dataFrame = get(load(input))
}

if(dim == 1){
    dataFrame = as.data.frame(t(as.matrix(dataFrame)))
}

#######################   making training and testing   ################
cat("\n","start making training and testing data...","\n")
inTrain = createDataPartition(y=labels,p=p)[[1]]
training = dataFrame[inTrain,]
testing = dataFrame[-inTrain,]
training_labels = data.frame(label=labels[inTrain])
rownames(training_labels) = rownames(training)
testing_labels = data.frame(label=labels[-inTrain])
rownames(testing_labels) = rownames(testing)
if(dim == 1){
    training = t(training)
    testing = t(testing)
}

save(training,file=paste(output,"_training",".RData",sep=""))
save(testing,file=paste(output,"_testing",".RData",sep=""))
write.table(training,file=paste(output,"_training",".txt",sep=""),quote=F,sep="\t")
write.table(testing,file=paste(output,"_testing",".txt",sep=""),quote=F,sep="\t")
write.table(training_labels,file=paste(output,"_trainingLabels",".txt",sep=""),quote=F,col.names=T,row.names=T,sep="\t")
write.table(testing_labels,file=paste(output,"_testingLabels",".txt",sep=""),quote=F,col.names=T,row.names=T,sep="\t")


