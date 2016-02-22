#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-20
# email: junjianglin@cs.toronto.com

#
# This r script is used to streamline the building of classifier
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


# This r script is used to do classification with 
# 1.svm, 
# 2.neural net,
# 3.random forest, 
# 4.recursive partitioning, 
# 5.stochastic gradient boosting, 
# 6.c5.0,
# 7.boosted logistic regression
# 8.knn


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
                    help="the seed, [Default] = 1"),
        make_option(c("-c","--classifier"),dest = 'classifier',default="all",
                    help="choosing classifers from svm,nnet,rf,rp,gbm,c50,logit,knn.(pls follow the order)[Default]=all")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$training) || is.null(opt$output)|| is.null(opt$testing)|| is.null(opt$testingLabel)|| is.null(opt$trainingLabel)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}
              
#libraries to load
library(caret)
library(ROCR)
library(scales)
library(tools)
library(ggplot2)
library(pROC)
library(doMC)
registerDoMC(cores=3)

cat("\n","Initializing the environment and loading data...","\n")
seed = opt$seed
set.seed(seed)
training = opt$training
testing = opt$testing
trainingLabel = opt$trainingLabel
testingLabel = opt$testingLabel
output = opt$output
classifier = opt$classifier
classifier = tolower(strsplit(classifier,",")[[1]])

# load or read training data
if(file_ext(training) == 'txt'){
    trainSet = read.delim(file = training)
} else if(file_ext(training) == 'RData'){
    trainSet = get(load(training))
}

# load or read testing data
if(file_ext(testing) == 'txt'){
    testSet = read.delim(file = testing)
} else if(file_ext(testing) == 'RData'){
    testSet = get(load(testing))
}

trainLabel = as.factor(as.character(read.table(file=trainingLabel)[,1]))
testLabel = as.factor(as.character(read.table(file=testingLabel)[,1]))


#######################   Outscourcing Functions  #############
# ROC
roceval <- function(myscores, labels_true, model) {
        
        pred <- prediction(myscores, labels_true)
        
        perf <- performance( pred, "tpr", "fpr" )
        
        auc <- unlist(slot(performance( pred, "auc" ), "y.values"))
        
        fpr <- unlist(slot(perf, "x.values"))
        
        tpr <- unlist(slot(perf, "y.values"))
        
        cutoffval <- unlist(slot(perf, "alpha.values"))	
        
        rocdf <- data.frame(x= fpr, y=tpr, auc=auc, 
                            
                            cutoff=cutoffval, 
                            
                            evalTypeAUC=sprintf("%s (%s)", model, percent(auc)), 
                            
                            model=model, curveType="ROC")
        
        return(rocdf)
}


# PRC
prceval <- function(myscores, labels_true, model) {
        
        pred <- prediction(myscores, labels_true)
        
        perf <- performance( pred, "prec", "rec" )
        
        rec <- unlist(slot(perf, "x.values"))
        
        prec <- unlist(slot(perf, "y.values"))
        
        cutoffval <- unlist(slot(perf, "alpha.values"))	
        
        prec[is.nan(prec)] <- 0
        
        prec[length(prec)] <- 0
        
        rec[is.nan(rec)] <- 0
        
        auc <- integrate(approxfun(cbind(rec, prec)), lower=0, upper=1,
                         subdivisions=1e5, stop.on.error = FALSE)$value					
        
        prcdf <- data.frame(x=rec, y=prec, auc=auc,
                            evalTypeAUC=sprintf("%s (%s)", model, 
                                                percent(auc)), model=model, curveType="PRC")
        
        return(prcdf)
}

#plotting
evalLinePlot <- function(mydf, curve, mytitle=NA) {
        
        if(curve=="ROC") {
                x.title <- "False positive rate"	
                y.title <- "True positive rate"
        } else {
                x.title <- "Recall"	
                y.title <- "Precision"		
        }
        
        gg <- ggplot(mydf, aes(x=x, y=y, color=evalTypeAUC)) + 
                
                theme_bw() +
                
                geom_line(aes(linetype=evalTypeAUC), size=0.7) +
                
                scale_x_continuous(x.title, labels=percent) +
                
                scale_y_continuous(y.title, labels=percent) + 		
                
                theme(axis.text.x = element_text(colour = "black"), 
                      
                      axis.text.y = element_text(colour = "black"),
                      
                      axis.title.x = element_text(colour = "black"),				
                      
                      legend.title= element_blank(),	
                      
                      legend.position=c(0.8,0.2),
                      
                      plot.title = element_text(colour = "black"),
                      
                      legend.background = element_blank())
        
        if(!is.na(mytitle)) gg <- gg + ggtitle(mytitle)
        
        if(curve=="ROC") gg + 
                geom_abline(intercept = 0, slope = 1, colour="grey", linetype=2) else gg	
}

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

#impute missing data if na exists
if(any(is.na(trainSet)) | any(is.na(testSet))){
    library(zoo)
    trainSet = na.spline(trainSet)
    testSet = na.spline(testSet)
}

#control should be in a higher level
if("CRC" %in% levels(trainLabel)){
  case = "CRC"
  control = "CTR"
  trainLabel = factor(trainLabel,levels=c(control,case))
  testLabel = factor(testLabel,levels=c(control,case))
}

roc_all = c()
resample_list = c()
sink(file=paste(output,"_ModelSummary",".txt",sep=""))
sink()
summary_list = c()

#get rid of highly correlated features( cor > 0.8 )

 idx_cor = findCorrelation(cor(trainSet),cutoff=0.8)
 if(length(idx_cor)>0){
   trainSet = trainSet[,-idx_cor]
   testSet = testSet[,-idx_cor]
 }


#some presetting for recursive feature elimination
numOfFeatures = ncol(trainSet)
if(numOfFeatures<200){
    subsets = seq(numOfFeatures)
  }else{
    subsets = c(1,seq(10,numOfFeatures,10))
  }
normalization <- preProcess(trainSet)
trainSet<- predict(normalization, trainSet)
trainSet<- as.data.frame(trainSet)

normalization2 = preProcess(testSet)
testSet = predict(normalization2,testSet)
testSet = as.data.frame(testSet)

#some general functions

#######################   1.SVM(RBF) with cv  ###############
if ('svm' %in% classifier || 'all' %in% classifier){
  cat("\n","SVM with repeated cv","\n")
  set.seed(seed)
  ## recursive feature elimination with svm
  #caretFuncs$summary <- twoClassSummary
  ctrl <- rfeControl(functions = caretFuncs,
                     method = "repeatedcv",
                     repeats = 2,
                     verbose = T,
                     returnResamp="final"
                     )


  svmProfileROC <- rfe(trainSet, trainLabel,
                    sizes = subsets,
                    rfeControl = ctrl,
                    method = "svmLinear",
                    #metric='ROC',
                    #trControl=trainControl(classProbs=TRUE)
                    )

  plot1 <- plot(svmProfileROC, type = c("g", "o"))
  pdf(file=paste(output,"_SVM_RLE",".pdf",sep=""))
  print(plot1)
  dev.off()
  pred_vars = predictors(svmProfileROC)
  cat("\n",sprintf("%s features were retained for SVM", length(pred_vars)),"\n")
  trainSet_svm = trainSet[,pred_vars]
  testSet_svm = testSet[,pred_vars]

  ## building svm classifier
  svmctrl = trainControl(method = 'repeatedcv',
                         repeats=3,
                         summaryFunction = twoClassSummary,
                         classProbs = TRUE,
                         returnResamp = "final")
  #tuning parameters
  svmgrid = expand.grid(C = exp(seq(log(0.01),log(10),length=8)),
                        sigma = exp(seq(log(0.01),log(1),length=5)))

  svm_tr = caret::train(trainSet_svm,trainLabel,
                method = 'svmRadial',
                trControl = svmctrl,
                tuneGrid = svmgrid,
                #tuneLength = 10,
                metric = 'ROC',
                preProc = c("center","scale"),
                verbose = F)

  plot <- plot(svm_tr,metric='ROC')
  pdf(file=paste(output,"_SVM_Hyper",".pdf",sep=""))
  print(plot)
  dev.off()
  svm_pred = predict(svm_tr,newdata=testSet_svm)
  svm_pred_prob = predict(svm_tr,newdata = testSet_svm, type='prob')
  svm_table = caret::confusionMatrix(data = svm_pred, testLabel)
  svm_roc = roceval(svm_pred_prob[[case]], as.numeric(testLabel), "SVM")
  svm_gg =  evalLinePlot(svm_roc, "ROC", "Test ROC for SVM")  
  roc_all = c(roc_all,list(svm_roc))
  resample_list = c(resample_list,list(SVM=svm_tr))
  #writing model summary 

  roc_obj = roc(testLabel,svm_pred_prob[[control]],auc=T)
  svm_all = coords(roc_obj,"best")
  svm_all$model = 'svm'
  svm_all$auc = as.numeric(roc_obj$auc)
  svm_all$size = length(pred_vars)
  summary_list = c(summary_list,list(svm_all))
  #ggsave(sprintf("%s/svm_9.png", ggout), gg_svm, width=6.4, height=4.8)
  sink(file=paste(output,"_ModelSummary",".txt",sep=""),append=T)
  cat("\n","Here is the final model for SVM","\n")
  cat("opt size of variables is ",svmProfileROC$optsize,"\n")
  print(svmProfileROC$fit)
  print(svm_tr)
  print(svm_tr$finalModel)
  print(svm_all)
  sink()
}


######################    2.Neural Network(single layer nnet) with repeatedcv  ####################
if("nnet" %in% classifier || 'all' %in% classifier){
  cat("\n","Neural Network(single layer nnet) with repeatedcv","\n")
  set.seed(seed)
  nnetctrl = trainControl(method = 'repeatedcv',
                        repeats = 2,
                        classProbs = T,
                        summaryFunction = twoClassSummary)

  nnetgrid = expand.grid(size = seq(50,700,length=5),
                        decay = c(1e-4,1e-3,1e-2,1e-1))
  nnet_tr = caret::train(trainSet,trainLabel,
               method = 'mlpWeightDecay',
               trControl = nnetctrl,
               tuneGrid = nnetgrid,
               metric = 'ROC',
               preProc = c("center","scale"))

  nnet_pred = predict(nnet_tr,newdata=testSet)
  nnet_pred_prob = predict(nnet_tr,newdata=testSet, type = 'prob')
  nnet_table = caret::confusionMatrix(data = nnet_pred, testLabel)
  nnet_roc = roceval(nnet_pred_prob[[case]], as.numeric(testLabel), "nnet")
  nnet_gg =  evalLinePlot(nnet_roc, "ROC", "Test ROC for Neural Network")  
  roc_all = c(roc_all,list(nnet_roc))
  resample_list = c(resample_list,list(NNET=nnet_tr))
  #ggsave(sprintf("%s/nnet_9.png", ggout), gg_nnet, width=6.4, height=4.8)
}


#######################   3.General Random Forest with repeated cv  #################
if("rf" %in% classifier || 'all' %in% classifier){

  cat("\n","general random forest with repeatedcv","\n")
  set.seed(seed)

  ##  recursive feature elimination
  rfFuncs$summary <- twoClassSummary
  ctrl <- rfeControl(functions = rfFuncs,
                     method = "repeatedcv",
                     #repeats = 1,
                     repeats = 2,
                     number = 10,
                     verbose = T)

  rfProfileROC <- rfe(trainSet, trainLabel,
                    sizes = subsets,
                    rfeControl = ctrl,
                    metric='ROC')

  plot1 <- plot(rfProfileROC, type = c("g", "o"))
  pdf(file=paste(output,"_RF_RLE",".pdf",sep=""))
  print(plot1)
  dev.off()
  pred_vars = predictors(rfProfileROC)
  cat("\n",sprintf("%s features were retained for random forest", length(pred_vars)),"\n")
  trainSet_rf = trainSet[,pred_vars]
  testSet_rf = testSet[,pred_vars]
  # classifier training

  rfctrl = trainControl(method = 'repeatedcv',
                        repeats = 3,
                        classProbs = T,
                        summaryFunction = twoClassSummary,
                        returnResamp = "final")
  rf_tr = caret::train(trainSet_rf,trainLabel,
               method = 'rf',
               tuneLength = 20,
               trControl = rfctrl,
               metric = 'ROC',
               preProc = c("center","scale"),
               verbose = F,
               ntree = 2000)
  plot2 <- plot(rf_tr,metric='ROC')
  pdf(file=paste(output,"_RF_Hyper",".pdf",sep=""))
  print(plot2)
  dev.off()
  # rf.profileROC <- rfe(trainSet_rf, trainLabel, 
  #                            sizes=subsets,
  #                            rfeControl=ctrl,
  #                            method="rf",
  #                            preProc = c("center","scale"),
  #                            tuneLength = 1,
  #                            metric = "ROC",
  #                            trControl = rfctrl)

  rf_pred = predict(rf_tr,newdata=testSet_rf)
  rf_pred_prob = predict(rf_tr,newdata=testSet_rf,type = "prob")
  rf_table = caret::confusionMatrix(data = rf_pred, testLabel)
  rf_roc = roceval(rf_pred_prob[[case]], as.numeric(testLabel), "RF")
  rf_gg =  evalLinePlot(rf_roc, "ROC", "Test ROC for Random Forest")
  roc_all = c(roc_all,list(rf_roc))  
  resample_list = c(resample_list,list(RF=rf_tr))
  #ggsave(sprintf("%s/rf_9.png", ggout), gg_rf, width=6.4, height=4.8)

  roc_obj = roc(testLabel,rf_pred_prob[[control]],auc=T)
  rf_all = coords(roc_obj,"best")
  rf_all$model = 'rf'
  rf_all$auc = as.numeric(roc_obj$auc)
  rf_all$size = length(pred_vars)
  summary_list = c(summary_list,list(rf_all))
  sink(file=paste(output,"_ModelSummary",".txt",sep=""),append=T)
  cat("\n","Here is the final model for RF","\n")
  cat("opt size of variables is ",rfProfileROC$optsize,"\n")
  print(rfProfileROC$fit)
  print(rf_tr)
  print(rf_tr$finalModel)
  print(rf_all)
  sink()
}


######################   4.Recursive Partitioning with repeatedcv ###########
if("rp" %in% classifier || 'all' %in% classifier){
  cat("\n","recursive partitioning with repeatedcv","\n")
  set.seed(seed)
  rpctrl = trainControl(method = 'repeatedcv',
                        repeats = 3,
                        classProbs = T,
                        summaryFunction = twoClassSummary)

  rp_tr = caret::train(trainSet,trainLabel,
               method = 'rpart',
               tuneLength = 15,
               trControl = rpctrl,
               metric = 'ROC',
               preProc = c("center","scale"))

  rp_pred = predict(rp_tr,newdata=testSet)
  rp_pred_prob = predict(rp_tr,newdata=testSet,type = "prob")
  rp_table = caret::confusionMatrix(data = rp_pred, testLabel)
  rp_roc = roceval(rp_pred_prob[[case]], as.numeric(testLabel), "RP")
  rp_gg =  evalLinePlot(rf_roc, "ROC", "Test ROC for Recursive Partitioning")  
  roc_all = c(roc_all,list(rp_roc))
  resample_list = c(resample_list,list(RP=rp_tr))
  #ggsave(sprintf("%s/rp_9.png", ggout), gg_rf, width=6.4, height=4.8)
}


##################### 5.Stochastic Gradient Boosting with repeatedcv #################
if("gbm" %in% classifier || 'all' %in% classifier){
  cat("\n","stochastic gradient boosting with repeatedcv","\n")
  set.seed(seed)
  gbmctrl = trainControl(method = 'repeatedcv',
                        repeats = 3,
                        classProbs = T,
                        summaryFunction = twoClassSummary,
                        returnResamp = "final")

  gbmGrid = expand.grid(interaction.depth = c(1, 5, 9),
                          n.trees = (1:40)*50,
                          shrinkage = c(0.001,0.01,0.1),
                          n.minobsinnode = 10)
  gbm_tr = caret::train(trainSet,trainLabel,
                method = 'gbm',
                tuneGrid = gbmGrid,
                verbose = F,
                trControl = gbmctrl,
                metric = 'ROC',
                preProc = c("center","scale"))
  plot <- plot(gbm_tr,metric='ROC')
  pdf(file=paste(output,"_GBM_Hyper",".pdf",sep=""))
  print(plot)
  dev.off()
  gbm_pred = predict(gbm_tr,newdata=testSet)
  gbm_pred_prob = predict(gbm_tr,newdata=testSet,type = "prob")
  gbm_table = caret::confusionMatrix(data = gbm_pred, testLabel)
  gbm_roc = roceval(gbm_pred_prob[[case]], as.numeric(testLabel), "gbm")
  gbm_gg =  evalLinePlot(gbm_roc, "ROC", "Test ROC for Stochastic Gradient Boosting")  
  roc_all = c(roc_all,list(gbm_roc))
  resample_list = c(resample_list,list(GBM=gbm_tr))

  roc_obj = roc(testLabel,gbm_pred_prob[[control]],auc=T)
  gbm_all = coords(roc_obj,"best")
  gbm_all$model = 'gbm'
  gbm_all$auc = as.numeric(roc_obj$auc)
  gbm_all$size = sum(relative.influence(gbm_tr$finalModel,length(gbm_tr$finalModel$train.error))>0)
  feat_imp = relative.influence(gbm_tr$finalModel,length(gbm_tr$finalModel$train.error))
  feat_imp_sorted = sort(feat_imp[feat_imp>0],decreasing=T)
  write.table(feat_imp_sorted,file=paste(output,"gbm_weight",".txt",sep=""),sep="\t",quote=F,col.names=F)
  summary_list = c(summary_list,list(gbm_all))
  sink(file=paste(output,"_ModelSummary",".txt",sep=""),append=T)
  cat("\n","Here is the final model for GBM","\n")
  print(gbm_tr)
  print(gbm_tr$finalModel)
  print(gbm_all)
  sink()
}


##################### 6. C5.0 with repeatedcv #################
if("c50" %in% classifier || 'all' %in% classifier){
  cat("\n","C5.0 with repeatedcv","\n")
  set.seed(seed)
  c50ctrl = trainControl(method = 'repeatedcv',
                         repeats = 3,
                         classProbs = T,
                         summaryFunction = twoClassSummary)
  c50Grid = expand.grid(model="tree",
                        trials=c(1:100),
                        winnow=FALSE)
  c50_tr = caret::train(trainSet,trainLabel,
                 method = 'C5.0',
                 tuneGrid = c50Grid,
                 trControl = c50ctrl,
                 metric = 'ROC',
                 preProc = c("center","scale"))

  c50_pred = predict(c50_tr,newdata=testSet)
  c50_pred_prob = predict(c50_tr,newdata=testSet,type = "prob")
  c50_table = caret::confusionMatrix(data = c50_pred, testLabel)
  c50_roc = roceval(c50_pred_prob[[case]], as.numeric(testLabel), "C5.0")
  c50_gg =  evalLinePlot(c50_roc, "ROC", "Test ROC for C5.0")  
  roc_all = c(roc_all,list(c50_roc))
  resample_list = c(resample_list,list(C50=c50_tr))
}


##################### 7.boosted logistic regression with repeatedcv #################
if("logit" %in% classifier || 'all' %in% classifier){
  cat("\n","boosted logistic regression with repeatedcv","\n")
  set.seed(seed)
  logitctrl = trainControl(method = 'repeatedcv',
                         repeats = 3,
                         classProbs = T,
                         summaryFunction = twoClassSummary)

  logit_tr = caret::train(trainSet,trainLabel,
                 method = 'LogitBoost',
                 tuneLength = 10,
                 trControl = logitctrl,
                 metric = 'ROC',
                 preProc = c("center","scale"))

  logit_pred = predict(logit_tr,newdata=testSet)
  logit_pred_prob = predict(logit_tr,newdata=testSet,type = "prob")
  logit_table = caret::confusionMatrix(data = logit_pred, testLabel)
  logit_roc = roceval(logit_pred_prob[[case]], as.numeric(testLabel), "logitBoost")
  logit_gg =  evalLinePlot(logit_roc, "ROC", "Test ROC for boosted logistic regression") 
  roc_all = c(roc_all,list(logit_roc))
  resample_list = c(resample_list,list(LOGIT=logit_tr))
}

##################### 8.knn with repeatedcv #################
if("knn" %in% classifier || 'all' %in% classifier){
  cat("\n","knn with repeatedcv","\n")
  set.seed(seed)
  knnctrl = trainControl(method = 'repeatedcv',
                         repeats = 3,
                         classProbs = T,
                         summaryFunction = twoClassSummary)

  knn_tr = caret::train(trainSet,trainLabel,
                 method = 'knn',
                 tuneLength = 20,
                 trControl = knnctrl,
                 metric = 'ROC',
                 preProc = c("center","scale"))

  knn_pred = predict(knn_tr,newdata=testSet)
  knn_pred_prob = predict(knn_tr,newdata=testSet,type = "prob")
  knn_table = caret::confusionMatrix(data = knn_pred, testLabel)
  knn_roc = roceval(knn_pred_prob[[case]], as.numeric(testLabel), "logitBoost")
  knn_gg =  evalLinePlot(knn_roc, "ROC", "Test ROC for knn") 
  roc_all = c(roc_all,list(knn_roc))
  resample_list = c(resample_list,list(KNN=knn_tr))
}

# #####################  Simple Voting  ################
# predDF = data.frame(svm_pred,
#                     nnet_pred,
#                     rf_pred,
#                     rp_pred,
#                     gbm_pred,
#                     c50_pred,
#                     logit_pred)

# comb_pred = factor(apply(predDF,1,function(row) ifelse(mean(row==case) > 0.5,case,control)),
#                    levels=(c(control,case)))
# comb_table = caret::confusionMatrix(data = comb_pred, testLabel)
# predDF$comb_pred = comb_pred

######################   Summary   ##################

cat("\n","generating summary","\n")
# ROC Curve 
#combineRoc = rbind(rf_roc,svm_roc,nnet_roc,rp_roc,gbm_roc,c50_roc,logit_roc)
generateSummary = function(roc_all,resample_list,summary_list,output){

  ##generating roc plots
  combineRoc = c()
  for(roc in roc_all){
    combineRoc = rbind(combineRoc,roc)
  }
  all_gg = evalLinePlot(combineRoc, "ROC", "Test ROC for All Models")  
  ggsave(all_gg,file=paste(output,"_MLROC",".pdf",sep=""))


  ##summary result
  summary = do.call(rbind,summary_list)
  nsum = ncol(summary)
  summary = summary[,c(4,5,(1:nsum)[-c(4,5)])]  #reorder
  write.table(summary,file=paste(output,"_MLAucSummary",".txt",sep=""),quote=F,col.names=T,row.names=F,sep="\t",append=T)

  ##resampling result
  resamps <- resamples(resample_list)
  sink(file=paste(output,"_ModelSummary",".txt",sep=""),append=T)
  cat("\n","Here is the between model information","\n")
  print(summary(resamps))
  sink()
  #resampling visualization 1
  pdf(file=paste(output,"_betweenModel1",".pdf",sep=""))
  bwplot = bwplot(resamps)
  print(bwplot)
  dev.off()
  #resampling visualization 2
  pdf(file=paste(output,"_betweenModel2",".pdf",sep=""))
  dotplot1 = dotplot(resamps, metric = "ROC")
  print(dotplot1)
  dev.off()
}

generateSummary(roc_all,resample_list,summary_list,output)
