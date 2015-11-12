#-------------------INFO----------------------
#
# author: Junjiang Lin
# date: 2015-03-19
# email: junjianglin@cs.toronto.com

#
# This r script is used to remove batch effects by using the fsva function in sva package
#
# Input List:
# 1.ExpressionSet data from Microarray having the phenotype and intensity info (RData)
# 2.A list of variables the user want see before and after batch removal, 
#   the name(Case-sensitive) of variables must be in the phenotype data of ExpressionSet
#
# Output List:
# 1. The corrected data, Both matrix(txt) and ExpressionSet(RData) format
# 2. Visualization showing the before and after batch removal(PCA will be used to visualize)

#######################   Initialization   ################
library(optparse)
rm(list=ls(all=TRUE))
set.seed(1)

option_list = list(
        make_option(c("-i","--input"),dest = "training",
                    help="the normarlized train RData"),
        make_option(c("-t","--test"),dest = "testing",
                    help="the normarlized test RData"),
        make_option(c("-o","--output"),dest = "output",
                    help="the primary name of output file"),
        make_option(c("-s","--sample"),dest = "samples",
                    help="the sample covariates files(separate by common, for example, \"batch1_info,batch2_info\")"),
        make_option(c("-r","--report"),dest = "report",action="store_true",default=TRUE,
                    help="whether to draw 2-D plots and test report to show the batch effects are removed,default[True]"),
        make_option(c("-m","--method"),dest = "method",default="Combat",
                    help="what method to remove batch effect,pls choose from Combat, SVA, Limma default[Combat]"),
        make_option(c("-n","--vars_num"),dest="vars_num",type='integer',default=100000,
                    help="the number of most variable features the SVA will use to compute the surrogate variable,only valid when use SVA.
                    default[100000]")
        )
opt = parse_args(OptionParser(option_list=option_list))
if(is.null(opt$training) || is.null(opt$output) || is.null(opt$samples)){
        print_help(OptionParser(option_list=option_list))
        q(save = "no")
}


#libraries to load
library(Biobase)
library(limma)
library(ggplot2)
library(pamr)
library(plyr)
library(sva)


cat("\n","Initializing the environment and loading data...","\n")
dataMatrix = get(load(opt$training))
dataMatrix_test = get(load(opt$testing))
sampleinfo = opt$samples
output = opt$output
report = opt$report
vars_n = opt$vars_num
method = opt$method
#######################   0. Preprocessing   ###############
cat("\n","preprocessing sample info data...","\n")

sampleinfo_list = strsplit(sampleinfo,",")[[1]]  # the list of sample info 


getDataSet = function(sampleinfo_list,dataMatrix){
    batch_number = length(sampleinfo_list)
    if(batch_number==1){
    sampleinfo = read.delim(sampleinfo_list[1],check.names=F)[,c('FileName','Sample','sub-batch','Array Name','Group','Age at Sample Collection (yrs)',
                                                    'Donor Sex','Expiration date','Scanning date','Lot#')]
    sampleinfo['Batch'] = 1
    vars_list = c("Subbatch","Lot","Age","Gender","Group")
    }else if(batch_number==2){
        batch1_info = read.delim(sampleinfo_list[1],check.names=F)[,c('FileName','Sample','sub-batch','Array Name','Group','Age at Sample Collection (yrs)',
                                                        'Donor Sex','Expiration date','Scanning date','Lot#')]
        batch1_info['Batch'] = 1
        batch2_info = read.delim(sampleinfo_list[2],check.names=F)[,c('FileName','Sample','sub-batch','Array Name','Group','Age at Sample Collection (yrs)',
                                                        'Donor Sex','Expiration date','Scanning date','Lot#')]
        batch2_info['Batch'] = 2
        sampleinfo = rbind(batch1_info,batch2_info)
        vars_list = c("Batch","Subbatch","Lot","Age","Gender","Group")
    }

    sampleinfo = rename(sampleinfo, c("sub-batch"="Subbatch", "Array Name"="ArrayName","Age at Sample Collection (yrs)"="Age","Donor Sex"="Gender",
                          "Expiration date"="ExpirationDate","Scanning date"="ScanningDate","Lot#"="Lot"))

    sampleinfo$Batch = factor(sampleinfo$Batch)
    sampleinfo$Subbatch = factor(sampleinfo$Subbatch)
    sampleinfo$Lot = factor(sampleinfo$Lot)
    levels(sampleinfo$Group) = list(CRC="CRC",CTR="CTR",CTR="HEA")
    data_samples = do.call(rbind,strsplit(colnames(dataMatrix),"[_]"))
    sampleinfo = sampleinfo[which(sampleinfo$Sample %in% data_samples[,1]),]
    rownames(sampleinfo) = paste(sampleinfo$Sample,'_',sampleinfo$Group,sep="")
    sampleinfo = sampleinfo[colnames(dataMatrix),]

    cat("\n","Encapsulation of ExpressionSet...","\n")
    pData = sampleinfo
    exprs = dataMatrix
    return(list(pData,exprs,vars_list,batch_number))
}
trainData = getDataSet(sampleinfo_list,dataMatrix)
pData_train = trainData[[1]]
eData_train = trainData[[2]]
vars_list = trainData[[3]]
batch_number = trainData[[4]]
testData = getDataSet(sampleinfo_list,dataMatrix_test)
pData_test = testData[[1]]
eData_test = testData[[2]]

#######################   1. Using PCA to visualize the original data  ###############
if(report == TRUE){
    getPCAScore = function(eData){
        n = nrow(eData)
        set.seed(1)
        rand_n = sample(n,n/3)
        pca_expr = prcomp(t(eData)[,rand_n],center=T)
        pca_score = data.frame(pca_expr$x[,1:5])  #only use first five for verification
        return(pca_score)
    }
    cat("\n","Using PCA to processing the original data...","\n")
    print(system.time({
        pca_score_train = getPCAScore(eData_train)
        cat("\n","PCA for training data done...","\n")}))
    print(system.time({
        pca_score_test = getPCAScore(eData_test)
        cat("\n","PCA for testing data done...","\n")}))
    
    
    # visualize different batch effect

    # for(i in vars_list){
    #         q = ggplot(pca_score,aes(x=pca_score$PC1,y=pca_score$PC2,color=as.factor(dataSet[[i]])))+
    #                 geom_point()+labs(title=paste("classification on",i,"before SVA",sep=" "),
    #                                   x="PC1",y="PC2")+
    #                 scale_color_manual(i,values=levels(as.factor(as.numeric(dataSet[[i]]))))
    #         ggsave(q,file=paste(output,"_beforeSVA","_",i,".pdf",sep=""))
    # }



########################   1.1 Test the correlation between all variables and original data  ###############
    cat("\n","correlation testing before correction...","\n")
    getTestReport = function(vars_list,pData,pca_score){
        design_formula = as.formula(paste("~ ",paste(paste(vars_list,collapse="+"))))
        design = model.matrix(design_formula,data=pData)
        fit = lmFit(t(pca_score),design)
        test_result = decideTests(fit,adjust.method="BH",p.value=0.05)
        return(test_result)
    }
    test_result_train = getTestReport(vars_list,pData_train,pca_score_train)
    test_result_test = getTestReport(vars_list,pData_test,pca_score_test)
    write.table(test_result_train,file=paste(output,"_","testReportTrain_before",method,".txt",sep=""),quote=F,sep="\t")
    write.table(test_result_train,file=paste(output,"_","testReportTest_before",method,".txt",sep=""),quote=F,sep="\t")
}
#######################   2. Using some methods to remove batch and other effects  ###############


#remove Batch Effect using sva and fsva function
rmBatch_sva = function(pData,eData_train,eData_test,vars_n){
    vars_num = dim(eData_train)[1]
    if(vars_num>vars_n){
        vars_num=vars_n
    }
    cat("\n","the number of variances used would be ",vars_num,"\n")
    #if(batch_number==1){
        #mod = model.matrix(~Subbatch+Lot+Age+Gender+Group,data=pData)
        #mod0 = model.matrix(~Subbatch+Lot+Age+Gender+1,data=pData)
        print(system.time({
            mod = model.matrix(~Group,data=pData)
            mod0 = model.matrix(~1,data=pData)
            cat("\n","build up model matrix","\n")}))
        print(system.time({
            n.sv = num.sv(eData_train,mod,vfilter=vars_num)
            cat("\n","done getting the number of surrogate variables","\n")}))
        print(system.time({
            svaobj = sva(eData_train,mod,mod0,vfilter=vars_num,n.sv=n.sv)
            cat("\n","done getting the sva obj","\n")}))
        print(system.time({
            fsvaobj = fsva(eData_train,mod,svaobj,eData_test,method="fast")
            cat("\n","done fsva part","\n")}))
        
        return(fsvaobj)
}
rmBatch_combat = function(pData,eData,batch_number){
    if(batch_number==1){
        modcombat = model.matrix(~Lot+Gender+1, data=pData)
        print(system.time({
            combat_edata = ComBat(dat=eData, batch=pData$Subbatch,mod=modcombat)
            cat("\n","done getting the combat data","\n")}))
        return(combat_edata)
    }else if(batch_number==2){
        modcombat = model.matrix(~Subbatch+Lot+Gender+1,data=pData)
        print(system.time({
            combat_edata = ComBat(dat=eData,batch=pData$Batch,mod=modcombat)
            cat("\n","done getting the combat data","\n")}))
        return(combat_edata)
    }

}
rmBatch_limma = function(pData,eData,batch_number){
    if(batch_number==1){
        mod = model.matrix(~Age+Gender+1, data=pData)
        print(system.time({
            limma_edata = removeBatchEffect(eData, batch=pData$Subbatch,batch2=pData$Lot,covariates=mod)
            cat("\n","done getting the limma adjusted data","\n")}))
        return(limma_edata)
    }else if(batch_number==2){
        mod = model.matrix(~Age+Lot+Gender+1,data=pData)
        print(system.time({
            limma_edata = removeBatchEffect(eData,batch=pData$Batch,batch2 = pData$Subbatch,mod=mod)
            cat("\n","done getting the limma adjusted data","\n")}))
        return(limma_edata)
    }
}
if(method=="SVA"){
    cat("\n","Using sva to remove batch effects...","\n")
    fsvaobj = rmBatch_sva(pData_train,eData_train,eData_test,vars_n)
    eData_train = fsvaobj$db
    eData_test = fsvaobj$new
    save(eData_train,file=paste(output,"_SVA_train.RData",sep=""))
    save(eData_test,file=paste(output,"_SVA_test.RData",sep=""))
}else if(method=="Combat"){
    cat("\n","Using combat to remove batch effects in training data...","\n")
    eData_train = rmBatch_combat(pData_train,eData_train,batch_number)
    cat("\n","Using combat to remove batch effects in testing data...","\n")
    eData_test = rmBatch_combat(pData_test,eData_test,batch_number)
    save(eData_train,file=paste(output,"_Combat_train.RData",sep=""))
    save(eData_test,file=paste(output,"_Combat_test.RData",sep=""))
}else if(method=="Limma"){
    cat("\n","Using limma to remove batch effects in training data...","\n")
    eData_train = rmBatch_limma(pData_train,eData_train,batch_number)
    cat("\n","Using limma to remove batch effects in testing data...","\n")
    eData_test = rmBatch_limma(pData_test,eData_test,batch_number)
    save(eData_train,file=paste(output,"_Limma_train.RData",sep=""))
    save(eData_test,file=paste(output,"_Limma_test.RData",sep=""))
}



#######################   3. Using PCA to visualize the corrected data  ###############
if(report==TRUE){
    cat("\n","Using PCA to processing the corrected data...","\n")
    pca_score_train = getPCAScore(eData_train)
    pca_score_test = getPCAScore(eData_test)

    # visualize different batch effect
    
    # for(i in vars_list){
    #         q = ggplot(pca_score,aes(x=pca_score$PC1,
    #                                          y=pca_score$PC2,color=as.factor(dataSet[[i]])))+
    #                 geom_point()+labs(title=paste("classification on",i,"after SVA",sep=" "),
    #                                   x="PC1",y="PC2")+
    #                 scale_color_manual(i,values=levels(as.factor(as.numeric(dataSet[[i]]))))
    #         ggsave(q,file=paste(output,"_afterSVA","_",i,".pdf",sep=""))
    # }

########################   3.1 Test the correlation between all variables and corrected data  ###############
    cat("\n","correlation testing after correction...","\n")
    test_result_train = getTestReport(vars_list,pData_train,pca_score_train)
    test_result_test = getTestReport(vars_list,pData_test,pca_score_test)
    write.table(test_result_train,file=paste(output,"_","testReportTrain_after",method,".txt",sep=""),quote=F,sep="\t")
    write.table(test_result_test,file=paste(output,"_","testReportTest_after",method,".txt",sep=""),quote=F,sep="\t")

}

########################   4 Save the corrected data in RData and matrix format  ###############

#big matrix can not be written, so here we split the matrix
# if(report==TRUE){
#     n = nrow(exprs(dataSet))
#     n1 = as.integer(n/3)
#     n2 = as.integer(n*2/3) 
#     write.table(exprs(dataSet)[1:n1,],file=paste(output,"_","LimmaBatch_matrix.txt",sep=""),quote=F,sep="\t")
#     write.table(exprs(dataSet)[(n1+1):n2,],file=paste(output,"_","LimmaBatch_matrix.txt",sep=""),quote=F,sep="\t",col.names=F,append=T)
#     write.table(exprs(dataSet)[(n2+1):n,],file=paste(output,"_","LimmaBatch_matrix.txt",sep=""),quote=F,sep="\t",col.names=F,append=T)
# }



