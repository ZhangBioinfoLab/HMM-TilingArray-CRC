#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-06-01
# email: junjianglin@cs.toronto.com

#
# Randomly select some number of features from the standard microarray data
# This is for verifying null hypothesis is not correct
#
# Input List:
# 1. Training data: microarray data for training ( Rdata, or txt format)
# 2. Testing data: microarray data 
# 3. The number of features you want to select

#
# Output List:
# 1. training data in RData(biobase) and Txt
# 2. testing data  in RData(biobase) and Txt


#------------------------
rm(list=ls(all=TRUE))

#######################   Initialization   ################
library(optparse)

option_list = list(
        make_option(c("-i","--training"),dest = "training",
                    help="the traing set of tilling array(after normalization),RData or txt"),
        make_option(c("-t","--testing"),dest = "testing",
                    help="the testing set of tilling array(after normalization),RData or txt"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output file"),
        make_option(c("-s","--seed"),type="integer",dest = "seed",default=1,
                    help = "the seed for reproducible research,
                    [default] = 1"),
        make_option(c("-n","--number"),dest="number",default=2000,
                    help = "the number of features to retain. [DEFAULT]=2000")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$training) || is.null(opt$output)|| is.null(opt$testing) ){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}

#libraries to load
library(tools)


cat("\n","Initializing the environment and loading data...","\n")
seed = opt$seed
training = opt$training
testing = opt$testing
output = opt$output
number = opt$number
# load or read training data
if(file_ext(training) == 'txt'){
    train.matrix = as.matrix(read.delim(file = training))
} else if(file_ext(training) == 'RData'){
    tmp = get(load(training))
    if(class(tmp) == "ExpressionSet"){
        train.matrix = exprs(tmp)
    }else{
        train.matrix = tmp
    }
}

# load or read testing data
if(file_ext(testing) == 'txt'){
    test.matrix = as.matrix(read.delim(file = testing))
} else if(file_ext(testing) == 'RData'){
    tmp = get(load(testing))
    if(class(tmp) == "ExpressionSet"){
        test.matrix = exprs(tmp)
    }else{
        test.matrix = tmp
    }   
}

#######################   Randomly select index   ################
nProbe = dim(train.matrix)[1]
select_number = number
if(select_number > nProbe){
    select_number = nProbe
}

set.seed(seed)
select_idx = sample(nProbe,select_number)
trainSet = t(train.matrix[select_idx,])
testSet = t(test.matrix[select_idx,])

cat("\n","Start saving training and testing dataset for machine learning ...","\n")
save(trainSet,file=paste(output,"_trainRand",".RData",sep=""))
save(testSet,file=paste(output,"_testRand",".RData",sep=""))
write.table(trainSet,file=paste(output,"_trainRand",".txt",sep=""),quote=F,sep="\t")
write.table(testSet,file=paste(output,"_testRand",".txt",sep=""),quote=F,sep="\t")

