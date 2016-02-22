library(randomForest)
library(boot)
library(e1071) # svm
library(ggplot2)
library(gdata)
library(gplots)
library(scales)
library(ROCR)
library(gridExtra)
library(reshape)
library(RColorBrewer)
#library(multicore)
#library(rtracklayer)
library(GenomicRanges)
library(limma)
library(glmnet)

rm(list=ls(all=TRUE))

#localdir <- "/Users/mike/Desktop/cirDNA_colon_cancer/colon_cancer_final"
localdir = "/Users/Lin/Box Sync/Research/Projects_git/cirDNA_project/output"

serverdir <- "/home/zhanglab1/yueli/Projects/colon_cancer_final"

if(file.exists(localdir)) mydir <- localdir else mydir <- serverdir

#source(sprintf("%s/scripts/classification/eval_functions.R", mydir))
source(sprintf("%s/eval_fun.R", mydir))

################ Input ################
#load(sprintf("%s/results/all_samples_pta/pta.RData", mydir))
load(sprintf("%s/train_test_set_9.rdata", mydir))

#trainset <- ptadata$train
#testset <- ptadata$test
trainset = trainset_9
testset = testset_9

################ Output ################
evalOut <- sprintf("%s/results/classification/CRC_classification.eps", mydir)

results.RData <- sprintf("%s/results/classification/CRC_classification.RData", mydir)

limmaOut <- sprintf("%s/results/classification/limma_PTA.csv", mydir)

localdir <- sprintf("%s/results/classification", mydir)

#dpdir <- "~/Dropbox/CRC_MS/YL/"

################ Apply Limma filter ################
design <- model.matrix(~0+factor(as.numeric(!trainset$class)+1))

colnames(design) <- c("CRC", "CTR")

x <- trainset

x$class <- NULL

x <- t(x)

#linear model fit
fit <- lmFit(x,design)

#contrast
contrast.matrix <- makeContrasts(CRC-CTR,levels=design)

fit2 <- contrasts.fit(fit, contrast.matrix)

fit2 <- eBayes(fit2)

limma.tab <- topTable(fit2, number=nrow(fit2))

sig.probesID <- subset(limma.tab, P.Value < 0.01)$ID

trainset <- trainset[,c(1, which(colnames(trainset) %in% sig.probesID))]

testset <- testset[,c(1, which(colnames(testset) %in% sig.probesID))]

trainset$class <- factor(trainset$class)

testset$class <- factor(testset$class)

################ RF ################

myRF <- randomForest(class ~ ., data=trainset, ntree=2000,
	xtest=testset[,-1], ytest=testset[,1],
	importance=FALSE, keep.forest=FALSE)


################ SVM ################
mysvm <- svm(class ~ ., data=trainset)
pred.tr.svm <- predict(mysvm, trainset[,-1], decision.values = TRUE)
pred.te.svm <- predict(mysvm, testset[,-1], decision.values = TRUE)


################ GLMNET ################
lrfit <- glmnet(as.matrix(trainset[,-1]), trainset$class, family="binomial")
pred.tr.glmnet <- predict(lrfit, as.matrix(trainset[,-1]), s=0.01, type = "response")
pred.te.glmnet <- predict(lrfit, as.matrix(testset[,-1]), s=0.01, type = "response")


# save outputs
save(trainset, testset, limma.tab, myRF, mysvm, pred.tr.svm, pred.te.svm, pred.tr.glmnet, pred.te.svm, file=results.RData)

write.csv(limma.tab, file=limmaOut, row.names=F)

################ Plot ################
trainLegend <- sprintf("Train (CRC:%s; CTR:%s)", table(trainset$class)[2], table(trainset$class)[1])

testLegend <- sprintf("Test (CRC:%s; CTR:%s)", table(testset$class)[2], table(testset$class)[1])

evalCombined.roc <- rbind(
	roceval(myRF$votes[,2], trainset$class, trainLegend, "RF"),	
	roceval(attr(pred.tr.svm, "decision.values")[,1], trainset$class, trainLegend, "SVM"),
	roceval(pred.tr.glmnet, trainset$class, trainLegend, "LR"),
	
	roceval(myRF$test$votes[,2], testset$class, testLegend, "RF"),	
	roceval(attr(pred.te.svm, "decision.values")[,1], testset$class, testLegend, "SVM"),
	roceval(pred.te.glmnet, testset$class, testLegend, "LR")
)

evalCombined.prc <- rbind(
	prceval(myRF$votes[,2], trainset$class, trainLegend, "RF"),	
	prceval(attr(pred.tr.svm, "decision.values")[,1], trainset$class, trainLegend, "SVM"),
	prceval(pred.tr.glmnet, trainset$class, trainLegend, "LR"),
	
	prceval(myRF$test$votes[,2], testset$class, testLegend, "RF"),	
	prceval(attr(pred.te.svm, "decision.values")[,1], testset$class, testLegend, "SVM"),
	prceval(pred.te.glmnet, testset$class, testLegend, "LR")
)

ggline.roc <- evalLinePlot(evalCombined.roc, "ROC", "ROC performance of PTA regions")

ggline.prc <- evalLinePlot(evalCombined.prc, "PRC", "Precision-Recall performance of PTA regions")


cairo_ps(evalOut, width=9.6, height=8.2); grid.arrange(ggline.roc, ggline.prc); dev.off()


#system(sprintf("rsync -avz --exclude screenlog.0 %s %s", localdir, dpdir))

# sync scripts from local to bc
# rsync -avz /Users/mike/Desktop/cirDNA_colon_cancer/colon_cancer_final/scripts/classification yueli@bc.ccbr.utoronto.ca:/home/zhanglab1/yueli/Projects/colon_cancer_final/scripts/



