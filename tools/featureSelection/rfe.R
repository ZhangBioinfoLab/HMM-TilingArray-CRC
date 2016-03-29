#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-20
# email: junjianglin@cs.toronto.com

#
# This r script is used to choose the appropriate feature size by using recursive feature elimination
#
# Input List:
# 1. Training data( RData or txt)
# 2. Training labels(txt)
# 3. Testing data(RData or txt)
# 4. Testing labels(txt)
#
# Output List:
# 1. classification auc summary
# 2. classification roc plot(only when choose all classifiers)



#######################   Initialization   ################
library(optparse)
rm(list=ls(all=TRUE))


option_list = list(
        make_option(c("-i","--training"),dest = "training",
                    help="the training set RData"),
        make_option(c("-t","--testing"),dest = "testing",
                    help="the testing set RData"),
        make_option(c("-I","--trainingLabel"),dest = "trainingLabel",
                    help="label file for training data"),
        make_option(c("-T","--testingLabel"),dest = "testingLabel",
                    help="label file for testing data"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output file"),
        make_option(c("-s","--seed"),dest = "seed",default=1,type="integer",
                    help="the seed, [Default] = 1")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$training) || is.null(opt$output)|| is.null(opt$testing)|| is.null(opt$testingLabel)|| is.null(opt$trainingLabel)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}
              
#libraries to load
library(caret)
library(scales)
library(tools)
library(ggplot2)


cat("\n","Initializing the environment and loading data...","\n")
seed = opt$seed
set.seed(seed)
training = opt$training
testing = opt$testing
trainingLabel = opt$trainingLabel
testingLabel = opt$testingLabel
output = opt$output


# load or read training data
if(file_ext(training) == 'txt'){
    trainSet = t(read.delim(file = training))
} else if(file_ext(training) == 'RData'){
    trainSet = get(load(training))
}

# load or read testing data
if(file_ext(testing) == 'txt'){
    testSet = t(read.delim(file = testing))
} else if(file_ext(testing) == 'RData'){
    testSet = get(load(testing))
}

trainLabel = as.factor(read.table(file=trainingLabel)[,1])
testLabel = as.factor(read.table(file=testingLabel)[,1])


#######################   Preprocessing  ################
# Changing all features' name to valid R names
cat("\n","preprocessing the data","\n")
#if(ncol(trainSet) == 1){
#        cat("\n","No training features","\n")
#        out = data.frame(name = opt$input,numOfFeatures=0)
#        write.table(out,file=opt$summary,quote=F,append=T,col.names=F,row.names=F)
#        q(save="no")
#}

names(trainSet) = make.names(names(trainSet))
names(testSet) = make.names(names(testSet))

#control should be in a higher level
case = "CRC"
control = "CTR"
trainLabel = factor(trainLabel,levels=c(control,case))
testLabel = factor(testLabel,levels=c(control,case))

numOfFeatures = ncol(trainSet)

#######################  Recursive Feature Elimination  ################

#centering and scaling
normalization <- preProcess(trainSet)
trainSet <- predict(normalization, trainSet)
trainSet <- as.data.frame(trainSet)
subsets = c(1,seq(10,numOfFeatures,10))

ctrl <- rfeControl(functions = rfFuncs,
                   method = "repeatedcv",
                   repeats = 3,
                   verbose = FALSE,
                   returnResamp = "all")

rfProfile <- rfe(trainSet, trainLabel,
                 sizes = subsets,
                 rfeControl = ctrl)

plot1 <- plot(rfProfile, type = c("g", "o"))

#plot2 <- xyplot(rfProfile,
#                type = c("g", "p", "smooth"),
#                ylab = "RMSE CV Estimates")

#print(plot1, split=c(1,1,1,2), more=TRUE)
#print(plot2, split=c(1,2,1,2))
print(plot1)

pred_vars = predictors(rfProfile)


